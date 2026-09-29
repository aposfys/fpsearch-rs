//! Tanimoto similarity over packed binary fingerprints.
//!
//! Fingerprints are stored as slices of `u64` words. Intersection and union counts use
//! `count_ones`, which lowers to a hardware popcount on aarch64, and on x86-64 when the
//! build enables `popcnt` (for example with `-C target-cpu=native`).
//!
//! The bound in [`max_possible_tanimoto`] depends only on the two popcounts, so a database
//! sorted by popcount lets a thresholded query skip whole regions without comparing a
//! single fingerprint.

use std::cmp::Ordering;

pub mod fps;
pub mod index;

pub use index::{Index, IndexBuilder, IndexError, Metric, SearchStats};

/// Errors that must not be silently absorbed into a similarity score.
#[derive(Debug, PartialEq, Eq)]
pub enum FpError {
    /// Fingerprints of different widths are not comparable.
    WidthMismatch { left: usize, right: usize },
    /// Tanimoto of two all-zero fingerprints is undefined, not 1.0.
    BothEmpty,
}

/// Number of set bits in a packed fingerprint.
#[inline]
pub fn popcount(fp: &[u64]) -> u32 {
    fp.iter().map(|word| word.count_ones()).sum()
}

/// Number of bits set in both fingerprints.
#[inline]
pub fn intersection(a: &[u64], b: &[u64]) -> u32 {
    a.iter()
        .zip(b.iter())
        .map(|(left, right)| (left & right).count_ones())
        .sum()
}

/// Tanimoto coefficient of two equal-width fingerprints.
///
/// Two empty fingerprints are an error rather than a perfect match: returning 1.0 there
/// would make every malformed record the best hit for every query.
pub fn tanimoto(a: &[u64], b: &[u64]) -> Result<f64, FpError> {
    if a.len() != b.len() {
        return Err(FpError::WidthMismatch {
            left: a.len(),
            right: b.len(),
        });
    }
    let both = intersection(a, b);
    let total = popcount(a) + popcount(b);
    if total == 0 {
        return Err(FpError::BothEmpty);
    }
    Ok(f64::from(both) / f64::from(total - both))
}

/// Upper bound on the Tanimoto of two fingerprints, from their popcounts alone.
///
/// Valid only for the standard binary Tanimoto. Applying it to a count-based or weighted
/// variant prunes true hits, which is why the index records the metric it was built for.
#[inline]
pub fn max_possible_tanimoto(popcount_a: u32, popcount_b: u32) -> f64 {
    let (smaller, larger) = if popcount_a <= popcount_b {
        (popcount_a, popcount_b)
    } else {
        (popcount_b, popcount_a)
    };
    if larger == 0 {
        return 0.0;
    }
    f64::from(smaller) / f64::from(larger)
}

/// Popcount range worth scanning for a query, given a similarity threshold.
///
/// Returns the inclusive `[low, high]` popcount band. Everything outside it is provably
/// below `threshold` and is never touched.
pub fn candidate_band(query_popcount: u32, threshold: f64) -> (u32, u32) {
    if threshold <= 0.0 {
        return (0, u32::MAX);
    }
    let query = f64::from(query_popcount);
    let mut low = (query * threshold).ceil() as u32;
    let mut high = (query / threshold).floor() as u32;
    // The float product and quotient can round one popcount inside the true band (7 / 0.07
    // is just under 100), so widen until the next popcount out cannot reach the threshold.
    // The bound uses the same division as the score, so the two agree exactly.
    while low > 0 && max_possible_tanimoto(low - 1, query_popcount) >= threshold {
        low -= 1;
    }
    while high < u32::MAX && max_possible_tanimoto(high + 1, query_popcount) >= threshold {
        high += 1;
    }
    (low, high)
}

/// One hit from a search.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hit {
    pub id: u64,
    pub score: f64,
}

impl Eq for Hit {}

impl Ord for Hit {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score
            .partial_cmp(&other.score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| other.id.cmp(&self.id))
    }
}

impl PartialOrd for Hit {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_fingerprints_score_one() {
        let fp = [0b1011u64, 0u64];
        assert_eq!(tanimoto(&fp, &fp).unwrap(), 1.0);
    }

    #[test]
    fn disjoint_fingerprints_score_zero() {
        assert_eq!(tanimoto(&[0b1010u64], &[0b0101u64]).unwrap(), 0.0);
    }

    #[test]
    fn matches_hand_computation() {
        // a = {0,1,2}, b = {1,2,3}: intersection 2, union 4.
        let score = tanimoto(&[0b0111u64], &[0b1110u64]).unwrap();
        assert!((score - 0.5).abs() < 1e-12);
    }

    #[test]
    fn two_empty_fingerprints_are_an_error_not_a_perfect_match() {
        assert_eq!(tanimoto(&[0u64], &[0u64]), Err(FpError::BothEmpty));
    }

    #[test]
    fn different_widths_are_rejected() {
        assert_eq!(
            tanimoto(&[0b1u64], &[0b1u64, 0u64]),
            Err(FpError::WidthMismatch { left: 1, right: 2 })
        );
    }

    #[test]
    fn the_bound_is_never_below_the_true_score() {
        // The pruning bound is only safe if it never underestimates. Exhaustive over a
        // small width, because an unsafe bound silently loses true hits.
        for a in 0u64..256 {
            for b in 0u64..256 {
                if a == 0 && b == 0 {
                    continue;
                }
                let truth = tanimoto(&[a], &[b]).unwrap();
                let bound = max_possible_tanimoto(popcount(&[a]), popcount(&[b]));
                assert!(bound >= truth - 1e-12, "bound {bound} < truth {truth}");
            }
        }
    }

    #[test]
    fn candidate_band_contains_every_qualifying_popcount() {
        let threshold = 0.7;
        let query = 40u32;
        let (low, high) = candidate_band(query, threshold);
        for other in 0u32..200 {
            if max_possible_tanimoto(query, other) >= threshold {
                assert!(
                    other >= low && other <= high,
                    "band excluded popcount {other}"
                );
            }
        }
    }

    #[test]
    fn candidate_band_survives_floating_point_rounding() {
        // 7 / 0.07 evaluates to just under 100, so a band computed from the float quotient
        // alone would drop a popcount-100 record that scores exactly 0.07.
        let (low, high) = candidate_band(7, 0.07);
        assert!(low <= 7 && high >= 100, "band [{low}, {high}] excludes 100");

        for step in 1..=100u32 {
            let threshold = f64::from(step) / 100.0;
            for query in 1u32..=256 {
                let (low, high) = candidate_band(query, threshold);
                for other in 0u32..=1024 {
                    if max_possible_tanimoto(query, other) >= threshold {
                        assert!(
                            other >= low && other <= high,
                            "t={threshold} q={query}: band [{low}, {high}] excluded {other}"
                        );
                    }
                }
            }
        }
    }
}
