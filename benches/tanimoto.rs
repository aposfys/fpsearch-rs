//! Kernel benchmarks. The numbers in the README come from here, on stated hardware.
//!
//! Two things are measured separately on purpose: the raw kernel throughput, which is a
//! memory-bandwidth question, and the effect of the popcount bound, which is an algorithmic
//! one. Reporting only the combined figure hides which of the two is actually doing the
//! work.

use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use fpsearch::{candidate_band, max_possible_tanimoto, popcount, tanimoto};

const WORDS: usize = 32; // 2048 bits
const DATABASE: usize = 100_000;

fn pseudo_random_fingerprint(seed: u64, words: usize) -> Vec<u64> {
    // A tiny xorshift, so the benchmark has no dependency beyond criterion and the
    // generated data is identical on every machine that runs it.
    let mut state = seed | 1;
    (0..words)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        })
        .collect()
}

fn bench_kernel(c: &mut Criterion) {
    let query = pseudo_random_fingerprint(42, WORDS);
    let database: Vec<Vec<u64>> = (0..DATABASE)
        .map(|i| pseudo_random_fingerprint(i as u64 + 1, WORDS))
        .collect();

    let mut group = c.benchmark_group("tanimoto");
    group.throughput(Throughput::Elements(DATABASE as u64));

    group.bench_function("scan_all", |b| {
        b.iter(|| {
            let mut best = 0.0f64;
            for fp in &database {
                let score = tanimoto(black_box(&query), black_box(fp)).unwrap();
                if score > best {
                    best = score;
                }
            }
            best
        })
    });

    let threshold = 0.7;
    let query_popcount = popcount(&query);
    let (low, high) = candidate_band(query_popcount, threshold);

    group.bench_function("scan_with_popcount_bound", |b| {
        b.iter(|| {
            let mut best = 0.0f64;
            for fp in &database {
                let other = popcount(fp);
                if other < low || other > high {
                    continue;
                }
                if max_possible_tanimoto(query_popcount, other) < best {
                    continue;
                }
                let score = tanimoto(black_box(&query), black_box(fp)).unwrap();
                if score > best {
                    best = score;
                }
            }
            best
        })
    });

    group.finish();
}

criterion_group!(benches, bench_kernel);
criterion_main!(benches);
