//! Memory-mapped, popcount-sorted fingerprint store.
//!
//! # Layout
//!
//! ```text
//! offset 0    header      64 bytes
//!             ids         count * u64
//!             popcounts   count * u32   (padded to an 8-byte boundary)
//!             payload     count * n_words * u64
//! ```
//!
//! Records are sorted by popcount ascending, so the candidate band from
//! [`crate::candidate_band`] is a contiguous slice that binary search locates in
//! `O(log n)`. Everything outside it is provably below the threshold and is never paged in.
//!
//! Every section starts on an 8-byte boundary, and an mmap base is page-aligned, so the
//! `u64` views into the payload are always correctly aligned. The constructors assert it
//! rather than assuming it.
//!
//! The payload is written in native byte order: an index is a local cache, not an
//! interchange format. The header records the writer's byte order and [`Index::open`]
//! refuses a file produced on the other kind of machine instead of returning silent
//! nonsense.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::mem::align_of;
use std::path::Path;

use memmap2::Mmap;

use crate::{candidate_band, popcount, tanimoto, Hit};

/// Magic bytes at the head of every index file.
pub const MAGIC: [u8; 8] = *b"FPSEARCH";
/// Format version. Bumped whenever the layout changes in a way older readers cannot handle.
pub const FORMAT_VERSION: u32 = 1;
/// Byte-order sentinel, written natively and compared on read.
const BYTE_ORDER_SENTINEL: u32 = 0x0102_0304;
/// Size of the fixed header, in bytes. A multiple of 8 so the sections after it stay aligned.
const HEADER_BYTES: usize = 64;

/// The similarity metric an index was built for.
///
/// The popcount bound in [`crate::max_possible_tanimoto`] is only valid for the standard
/// binary Tanimoto. Recording the metric means a future count-based variant cannot silently
/// reuse a bound that would prune true hits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum Metric {
    /// Standard binary Tanimoto: `|a & b| / |a | b|`.
    TanimotoBinary = 0,
}

impl Metric {
    fn from_u32(raw: u32) -> Option<Self> {
        match raw {
            0 => Some(Metric::TanimotoBinary),
            _ => None,
        }
    }
}

/// Anything that stops an index being read or written.
#[derive(Debug)]
pub enum IndexError {
    Io(std::io::Error),
    /// The file does not start with [`MAGIC`].
    NotAnIndex,
    /// Written by a different format version.
    UnsupportedVersion {
        found: u32,
        supported: u32,
    },
    /// Written on a machine with the opposite byte order.
    ByteOrderMismatch,
    /// The metric field holds a value this build does not know.
    UnknownMetric(u32),
    /// The file is shorter than its own header says it should be.
    Truncated {
        expected: usize,
        actual: usize,
    },
    /// A query fingerprint of a different fold width than the index.
    ///
    /// Folding to a different width makes different substructures share a bit, so the two
    /// are not comparable. Rejected rather than compared.
    WidthMismatch {
        index_bits: u32,
        query_bits: u32,
    },
    /// A record whose fingerprint is not the index's width.
    RecordWidthMismatch {
        expected: usize,
        actual: usize,
    },
    /// A fingerprint with no bits set. Its Tanimoto against anything is undefined.
    EmptyFingerprint {
        id: u64,
    },
}

impl std::fmt::Display for IndexError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IndexError::Io(e) => write!(f, "io error: {e}"),
            IndexError::NotAnIndex => write!(f, "not an fpsearch index (bad magic)"),
            IndexError::UnsupportedVersion { found, supported } => write!(
                f,
                "index format version {found} is not supported (this build reads {supported})"
            ),
            IndexError::ByteOrderMismatch => {
                write!(
                    f,
                    "index was written on a machine with the opposite byte order"
                )
            }
            IndexError::UnknownMetric(raw) => write!(f, "unknown metric discriminant {raw}"),
            IndexError::Truncated { expected, actual } => {
                write!(
                    f,
                    "index is truncated: header implies {expected} bytes, file is {actual}"
                )
            }
            IndexError::WidthMismatch {
                index_bits,
                query_bits,
            } => write!(
                f,
                "query is {query_bits}-bit but the index was built at {index_bits}-bit; \
                 folding widths are not comparable"
            ),
            IndexError::RecordWidthMismatch { expected, actual } => {
                write!(f, "record is {actual} words, index width is {expected}")
            }
            IndexError::EmptyFingerprint { id } => {
                write!(
                    f,
                    "fingerprint {id} has no bits set; its similarity is undefined"
                )
            }
        }
    }
}

impl std::error::Error for IndexError {}

impl From<std::io::Error> for IndexError {
    fn from(e: std::io::Error) -> Self {
        IndexError::Io(e)
    }
}

/// Collects fingerprints, sorts them by popcount, and writes an index file.
///
/// Sorting happens once here so that every subsequent query gets a contiguous candidate
/// band for free.
pub struct IndexBuilder {
    n_words: usize,
    n_bits: u32,
    records: Vec<(u64, u32, Vec<u64>)>,
}

impl IndexBuilder {
    /// Start an index for fingerprints of `n_bits` bits.
    ///
    /// `n_bits` is the fold width, which is what makes two fingerprints comparable; the
    /// word count is derived from it.
    pub fn new(n_bits: u32) -> Self {
        let n_words = (n_bits as usize).div_ceil(64);
        IndexBuilder {
            n_words,
            n_bits,
            records: Vec::new(),
        }
    }

    /// Number of fingerprints added so far.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether nothing has been added yet.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Add one fingerprint.
    ///
    /// An all-zero fingerprint is refused: it has no defined similarity to anything, and
    /// admitting it would push that error into every future query instead of this call.
    pub fn push(&mut self, id: u64, fingerprint: &[u64]) -> Result<(), IndexError> {
        if fingerprint.len() != self.n_words {
            return Err(IndexError::RecordWidthMismatch {
                expected: self.n_words,
                actual: fingerprint.len(),
            });
        }
        let bits = popcount(fingerprint);
        if bits == 0 {
            return Err(IndexError::EmptyFingerprint { id });
        }
        self.records.push((id, bits, fingerprint.to_vec()));
        Ok(())
    }

    /// Sort by popcount and write the index to `path`.
    pub fn write(mut self, path: impl AsRef<Path>) -> Result<(), IndexError> {
        self.records.sort_by_key(|(id, bits, _)| (*bits, *id));

        let count = self.records.len() as u64;
        let file = File::create(path)?;
        let mut out = BufWriter::new(file);

        let mut header = [0u8; HEADER_BYTES];
        header[0..8].copy_from_slice(&MAGIC);
        header[8..12].copy_from_slice(&FORMAT_VERSION.to_ne_bytes());
        header[12..16].copy_from_slice(&BYTE_ORDER_SENTINEL.to_ne_bytes());
        header[16..20].copy_from_slice(&(self.n_words as u32).to_ne_bytes());
        header[20..24].copy_from_slice(&self.n_bits.to_ne_bytes());
        header[24..28].copy_from_slice(&(Metric::TanimotoBinary as u32).to_ne_bytes());
        header[28..32].copy_from_slice(&0u32.to_ne_bytes()); // reserved
        header[32..40].copy_from_slice(&count.to_ne_bytes());
        out.write_all(&header)?;

        for (id, _, _) in &self.records {
            out.write_all(&id.to_ne_bytes())?;
        }
        for (_, bits, _) in &self.records {
            out.write_all(&bits.to_ne_bytes())?;
        }
        // Pad the popcount section so the payload starts on an 8-byte boundary.
        let popcount_bytes = self.records.len() * 4;
        let padding = (8 - (popcount_bytes % 8)) % 8;
        out.write_all(&vec![0u8; padding])?;

        for (_, _, fingerprint) in &self.records {
            for word in fingerprint {
                out.write_all(&word.to_ne_bytes())?;
            }
        }
        out.flush()?;
        Ok(())
    }
}

/// How much of the database a query actually touched.
///
/// The point of the popcount bound is that `examined` is far below `total`. Reporting a
/// query time without this number leaves the reader unable to tell whether the bound did
/// anything.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SearchStats {
    /// Fingerprints in the index.
    pub total: usize,
    /// Fingerprints whose popcount fell inside the candidate band.
    pub examined: usize,
    /// Fingerprints skipped without a single word comparison.
    pub pruned: usize,
}

impl SearchStats {
    /// Fraction of the database skipped, between 0.0 and 1.0.
    pub fn pruned_fraction(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        self.pruned as f64 / self.total as f64
    }
}

/// A memory-mapped fingerprint index.
pub struct Index {
    mmap: Mmap,
    n_words: usize,
    n_bits: u32,
    metric: Metric,
    count: usize,
    ids_offset: usize,
    popcounts_offset: usize,
    payload_offset: usize,
}

impl Index {
    /// Open an index file.
    ///
    /// The file is mapped, not read: a 256 GB index costs no resident memory here, and the
    /// popcount band means a query pages in only the region it can possibly match.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, IndexError> {
        let file = File::open(path)?;
        // SAFETY: the map is read-only and the file is not modified while mapped. A
        // concurrent truncation would be a torn read, which is the standard caveat on mmap.
        let mmap = unsafe { Mmap::map(&file)? };

        if mmap.len() < HEADER_BYTES || mmap[0..8] != MAGIC {
            return Err(IndexError::NotAnIndex);
        }
        let version = u32::from_ne_bytes(mmap[8..12].try_into().unwrap());
        if version != FORMAT_VERSION {
            return Err(IndexError::UnsupportedVersion {
                found: version,
                supported: FORMAT_VERSION,
            });
        }
        let sentinel = u32::from_ne_bytes(mmap[12..16].try_into().unwrap());
        if sentinel != BYTE_ORDER_SENTINEL {
            return Err(IndexError::ByteOrderMismatch);
        }
        let n_words = u32::from_ne_bytes(mmap[16..20].try_into().unwrap()) as usize;
        let n_bits = u32::from_ne_bytes(mmap[20..24].try_into().unwrap());
        let raw_metric = u32::from_ne_bytes(mmap[24..28].try_into().unwrap());
        let metric = Metric::from_u32(raw_metric).ok_or(IndexError::UnknownMetric(raw_metric))?;
        let count = u64::from_ne_bytes(mmap[32..40].try_into().unwrap()) as usize;

        let ids_offset = HEADER_BYTES;
        let popcounts_offset = ids_offset + count * 8;
        let popcount_bytes = count * 4;
        let padding = (8 - (popcount_bytes % 8)) % 8;
        let payload_offset = popcounts_offset + popcount_bytes + padding;
        let expected = payload_offset + count * n_words * 8;
        if mmap.len() < expected {
            return Err(IndexError::Truncated {
                expected,
                actual: mmap.len(),
            });
        }

        let index = Index {
            mmap,
            n_words,
            n_bits,
            metric,
            count,
            ids_offset,
            popcounts_offset,
            payload_offset,
        };
        // Every section offset is a multiple of 8 and an mmap base is page-aligned, so this
        // holds by construction. Asserted because the `u64` views below depend on it.
        assert_eq!(
            index.mmap.as_ptr() as usize % align_of::<u64>(),
            0,
            "mmap base is not 8-byte aligned"
        );
        Ok(index)
    }

    /// Number of fingerprints in the index.
    pub fn len(&self) -> usize {
        self.count
    }

    /// Whether the index holds no fingerprints.
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Fold width in bits. A query built at a different width is not comparable.
    pub fn n_bits(&self) -> u32 {
        self.n_bits
    }

    /// Fingerprint width in 64-bit words.
    pub fn n_words(&self) -> usize {
        self.n_words
    }

    /// The metric this index was built for.
    pub fn metric(&self) -> Metric {
        self.metric
    }

    /// Identifiers, in popcount order.
    pub fn ids(&self) -> &[u64] {
        let bytes = &self.mmap[self.ids_offset..self.ids_offset + self.count * 8];
        debug_assert_eq!(bytes.as_ptr() as usize % align_of::<u64>(), 0);
        // SAFETY: the slice is in-bounds, 8-byte aligned (asserted in `open`), and `u64`
        // has no invalid bit patterns.
        unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u64, self.count) }
    }

    /// Popcounts, ascending. This ordering is what makes the candidate band contiguous.
    pub fn popcounts(&self) -> &[u32] {
        let bytes = &self.mmap[self.popcounts_offset..self.popcounts_offset + self.count * 4];
        debug_assert_eq!(bytes.as_ptr() as usize % align_of::<u32>(), 0);
        // SAFETY: as above, for `u32`.
        unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u32, self.count) }
    }

    /// The `i`th fingerprint, in popcount order.
    pub fn fingerprint(&self, i: usize) -> &[u64] {
        assert!(
            i < self.count,
            "index {i} out of range for {} records",
            self.count
        );
        let start = self.payload_offset + i * self.n_words * 8;
        let bytes = &self.mmap[start..start + self.n_words * 8];
        // SAFETY: as above.
        unsafe { std::slice::from_raw_parts(bytes.as_ptr() as *const u64, self.n_words) }
    }

    /// First position whose popcount is at least `bits`.
    fn lower_bound(&self, bits: u32) -> usize {
        self.popcounts().partition_point(|&p| p < bits)
    }

    /// First position whose popcount is greater than `bits`.
    fn upper_bound(&self, bits: u32) -> usize {
        self.popcounts().partition_point(|&p| p <= bits)
    }

    /// Top-`k` most similar fingerprints scoring at least `threshold`.
    ///
    /// Returns hits in descending score order, with the statistics describing how much of
    /// the database the popcount bound let the query skip.
    ///
    /// The heap is preallocated to `k + 1` and reused, so the hot loop does not allocate.
    pub fn search(
        &self,
        query: &[u64],
        threshold: f64,
        k: usize,
    ) -> Result<(Vec<Hit>, SearchStats), IndexError> {
        let query_bits = (query.len() * 64) as u32;
        if query.len() != self.n_words {
            return Err(IndexError::WidthMismatch {
                index_bits: self.n_bits,
                query_bits,
            });
        }
        let query_popcount = popcount(query);
        if query_popcount == 0 {
            return Err(IndexError::EmptyFingerprint { id: u64::MAX });
        }

        let (low, high) = candidate_band(query_popcount, threshold);
        let start = self.lower_bound(low);
        let end = self.upper_bound(high);
        let examined = end.saturating_sub(start);
        let stats = SearchStats {
            total: self.count,
            examined,
            pruned: self.count - examined,
        };
        let mut hits = self.scan_range(query, threshold, k, start, end);
        hits.sort_by(|a, b| b.cmp(a));
        Ok((hits, stats))
    }

    /// Every hit at or above `threshold`, without a `k` ceiling.
    pub fn search_all(
        &self,
        query: &[u64],
        threshold: f64,
    ) -> Result<(Vec<Hit>, SearchStats), IndexError> {
        self.search(query, threshold, self.count)
    }

    /// Scan `[start, end)` and return the best `k` hits in that range.
    ///
    /// Shared by the serial and parallel paths so the two cannot drift apart.
    fn scan_range(
        &self,
        query: &[u64],
        threshold: f64,
        k: usize,
        start: usize,
        end: usize,
    ) -> Vec<Hit> {
        if k == 0 || start >= end {
            return Vec::new();
        }
        let mut heap: std::collections::BinaryHeap<std::cmp::Reverse<Hit>> =
            std::collections::BinaryHeap::with_capacity(k + 1);
        for (offset, &id) in self.ids()[start..end].iter().enumerate() {
            let score = tanimoto(query, self.fingerprint(start + offset))
                .expect("width checked by caller; empty records refused at build time");
            if score < threshold {
                continue;
            }
            if heap.len() < k {
                heap.push(std::cmp::Reverse(Hit { id, score }));
            } else if let Some(std::cmp::Reverse(weakest)) = heap.peek() {
                if score > weakest.score {
                    heap.pop();
                    heap.push(std::cmp::Reverse(Hit { id, score }));
                }
            }
        }
        heap.into_iter().map(|std::cmp::Reverse(hit)| hit).collect()
    }

    /// Top-`k` search with the candidate band split across `threads` worker threads.
    ///
    /// The scan is embarrassingly parallel — every candidate is scored independently — and
    /// the serial version is memory-bandwidth bound well below what the machine can supply,
    /// so this is close to a linear win. Each worker keeps its own heap and the results are
    /// merged at the end, which costs `threads * k` and needs no shared state during the
    /// scan.
    ///
    /// Uses scoped threads from the standard library; there is no thread-pool dependency.
    pub fn search_parallel(
        &self,
        query: &[u64],
        threshold: f64,
        k: usize,
        threads: usize,
    ) -> Result<(Vec<Hit>, SearchStats), IndexError> {
        let query_bits = (query.len() * 64) as u32;
        if query.len() != self.n_words {
            return Err(IndexError::WidthMismatch {
                index_bits: self.n_bits,
                query_bits,
            });
        }
        let query_popcount = popcount(query);
        if query_popcount == 0 {
            return Err(IndexError::EmptyFingerprint { id: u64::MAX });
        }

        let (low, high) = candidate_band(query_popcount, threshold);
        let start = self.lower_bound(low);
        let end = self.upper_bound(high);
        let examined = end.saturating_sub(start);
        let stats = SearchStats {
            total: self.count,
            examined,
            pruned: self.count - examined,
        };
        let threads = threads.max(1);
        if k == 0 || examined == 0 {
            return Ok((Vec::new(), stats));
        }
        if threads == 1 || examined < threads * 1024 {
            // Below this size the thread spawn costs more than the scan saves.
            let mut hits = self.scan_range(query, threshold, k, start, end);
            hits.sort_by(|a, b| b.cmp(a));
            return Ok((hits, stats));
        }

        let chunk = examined.div_ceil(threads);
        let mut merged: Vec<Hit> = std::thread::scope(|scope| {
            let handles: Vec<_> = (0..threads)
                .map(|t| {
                    let lo = start + t * chunk;
                    let hi = (lo + chunk).min(end);
                    scope.spawn(move || self.scan_range(query, threshold, k, lo.min(end), hi))
                })
                .collect();
            handles
                .into_iter()
                .flat_map(|h| h.join().expect("search worker panicked"))
                .collect()
        });

        merged.sort_by(|a, b| b.cmp(a));
        merged.truncate(k);
        Ok((merged, stats))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fingerprint_from_bits(bits: &[usize], n_words: usize) -> Vec<u64> {
        let mut fp = vec![0u64; n_words];
        for &bit in bits {
            fp[bit / 64] |= 1u64 << (bit % 64);
        }
        fp
    }

    fn temp_path(name: &str) -> std::path::PathBuf {
        let mut path = std::env::temp_dir();
        path.push(format!("fpsearch-test-{}-{}.idx", name, std::process::id()));
        path
    }

    #[test]
    fn round_trips_and_sorts_by_popcount() {
        let path = temp_path("roundtrip");
        let mut builder = IndexBuilder::new(128);
        builder
            .push(10, &fingerprint_from_bits(&[0, 1, 2, 3], 2))
            .unwrap();
        builder.push(11, &fingerprint_from_bits(&[0], 2)).unwrap();
        builder
            .push(12, &fingerprint_from_bits(&[0, 1], 2))
            .unwrap();
        builder.write(&path).unwrap();

        let index = Index::open(&path).unwrap();
        assert_eq!(index.len(), 3);
        assert_eq!(index.n_bits(), 128);
        assert_eq!(index.metric(), Metric::TanimotoBinary);
        // Ascending popcount, which is what the band relies on.
        assert_eq!(index.popcounts(), &[1, 2, 4]);
        assert_eq!(index.ids(), &[11, 12, 10]);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn search_agrees_with_a_brute_force_scan() {
        // The whole point of the band is that it changes cost, not answers. This asserts
        // the pruned search returns exactly what an exhaustive scan would.
        let path = temp_path("bruteforce");
        let n_words = 4;
        let mut builder = IndexBuilder::new(256);
        let mut raw = Vec::new();
        for id in 0..400u64 {
            let bits: Vec<usize> = (0..256)
                .filter(|b| (id.wrapping_mul(2_654_435_761).wrapping_add(*b as u64)) % 7 == 0)
                .collect();
            if bits.is_empty() {
                continue;
            }
            let fp = fingerprint_from_bits(&bits, n_words);
            builder.push(id, &fp).unwrap();
            raw.push((id, fp));
        }
        builder.write(&path).unwrap();
        let index = Index::open(&path).unwrap();

        let query = raw[3].1.clone();
        for &threshold in &[0.0_f64, 0.3, 0.5, 0.7, 0.9] {
            let (hits, stats) = index.search_all(&query, threshold).unwrap();

            let mut expected: Vec<Hit> = raw
                .iter()
                .filter_map(|(id, fp)| {
                    let score = tanimoto(&query, fp).unwrap();
                    (score >= threshold).then_some(Hit { id: *id, score })
                })
                .collect();
            expected.sort_by(|a, b| b.cmp(a));

            assert_eq!(hits.len(), expected.len(), "threshold {threshold}");
            for (got, want) in hits.iter().zip(expected.iter()) {
                assert_eq!(got.id, want.id, "threshold {threshold}");
                assert!((got.score - want.score).abs() < 1e-12);
            }
            assert!(stats.examined <= stats.total);
            assert_eq!(stats.examined + stats.pruned, stats.total);
        }
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn top_k_is_the_k_best() {
        let path = temp_path("topk");
        let mut builder = IndexBuilder::new(128);
        for id in 1..=50u64 {
            let bits: Vec<usize> = (0..(id as usize % 40 + 1)).collect();
            builder.push(id, &fingerprint_from_bits(&bits, 2)).unwrap();
        }
        builder.write(&path).unwrap();
        let index = Index::open(&path).unwrap();

        let query = fingerprint_from_bits(&(0..20).collect::<Vec<_>>(), 2);
        let (top, _) = index.search(&query, 0.0, 5).unwrap();
        let (all, _) = index.search_all(&query, 0.0).unwrap();
        assert_eq!(top.len(), 5);
        for (a, b) in top.iter().zip(all.iter().take(5)) {
            assert!((a.score - b.score).abs() < 1e-12);
        }
        // Descending order.
        for pair in top.windows(2) {
            assert!(pair[0].score >= pair[1].score);
        }
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn the_band_actually_prunes() {
        // If this stops being true the index has degenerated into a linear scan.
        let path = temp_path("prunes");
        let mut builder = IndexBuilder::new(256);
        for id in 1..=500u64 {
            let n = (id as usize % 200) + 1;
            builder
                .push(id, &fingerprint_from_bits(&(0..n).collect::<Vec<_>>(), 4))
                .unwrap();
        }
        builder.write(&path).unwrap();
        let index = Index::open(&path).unwrap();

        let query = fingerprint_from_bits(&(0..100).collect::<Vec<_>>(), 4);
        let (_, stats) = index.search_all(&query, 0.9).unwrap();
        assert!(
            stats.pruned_fraction() > 0.5,
            "expected the 0.9 band to skip most of the database, skipped {:.2}",
            stats.pruned_fraction()
        );
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn parallel_search_returns_exactly_what_the_serial_one_does() {
        // Splitting the band across threads must change the cost and nothing else. A
        // per-thread heap that merged incorrectly would still return plausible hits.
        let path = temp_path("parallel");
        let n_words = 4;
        let mut builder = IndexBuilder::new(256);
        for id in 1..=5000u64 {
            let bits: Vec<usize> = (0..256)
                .filter(|b| (id.wrapping_mul(2_654_435_761).wrapping_add(*b as u64)) % 5 < 2)
                .collect();
            if bits.is_empty() {
                continue;
            }
            builder
                .push(id, &fingerprint_from_bits(&bits, n_words))
                .unwrap();
        }
        builder.write(&path).unwrap();
        let index = Index::open(&path).unwrap();

        let query = index.fingerprint(17).to_vec();
        for &threads in &[1usize, 2, 4, 8] {
            for &k in &[1usize, 10, 50] {
                let (serial, s1) = index.search(&query, 0.4, k).unwrap();
                let (parallel, s2) = index.search_parallel(&query, 0.4, k, threads).unwrap();
                assert_eq!(s1, s2, "stats differ at {threads} threads");
                assert_eq!(serial.len(), parallel.len(), "{threads} threads, k={k}");
                for (a, b) in serial.iter().zip(parallel.iter()) {
                    assert!(
                        (a.score - b.score).abs() < 1e-12,
                        "{threads} threads, k={k}: {a:?} vs {b:?}"
                    );
                }
            }
        }
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn a_query_of_another_width_is_refused_not_compared() {
        let path = temp_path("width");
        let mut builder = IndexBuilder::new(128);
        builder.push(1, &fingerprint_from_bits(&[0, 5], 2)).unwrap();
        builder.write(&path).unwrap();
        let index = Index::open(&path).unwrap();

        let wrong_width = fingerprint_from_bits(&[0, 5], 4); // 256-bit
        match index.search(&wrong_width, 0.5, 10) {
            Err(IndexError::WidthMismatch {
                index_bits,
                query_bits,
            }) => {
                assert_eq!(index_bits, 128);
                assert_eq!(query_bits, 256);
            }
            other => panic!("expected WidthMismatch, got {other:?}"),
        }
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn an_empty_fingerprint_is_refused_at_build_time() {
        let mut builder = IndexBuilder::new(128);
        match builder.push(7, &[0u64; 2]) {
            Err(IndexError::EmptyFingerprint { id }) => assert_eq!(id, 7),
            other => panic!("expected EmptyFingerprint, got {other:?}"),
        }
    }

    #[test]
    fn a_non_index_file_is_rejected() {
        let path = temp_path("garbage");
        std::fs::write(&path, b"this is not an index at all").unwrap();
        assert!(matches!(Index::open(&path), Err(IndexError::NotAnIndex)));
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn a_truncated_index_is_rejected() {
        let path = temp_path("truncated");
        let mut builder = IndexBuilder::new(128);
        for id in 1..=20u64 {
            builder
                .push(id, &fingerprint_from_bits(&[0, id as usize % 100], 2))
                .unwrap();
        }
        builder.write(&path).unwrap();
        let full = std::fs::read(&path).unwrap();
        std::fs::write(&path, &full[..full.len() - 8]).unwrap();
        assert!(matches!(
            Index::open(&path),
            Err(IndexError::Truncated { .. })
        ));
        std::fs::remove_file(&path).ok();
    }
}
