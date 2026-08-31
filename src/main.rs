//! `fpsearch` — build an index, query it, and measure it.
//!
//! Argument parsing is hand-rolled because the CLI has four subcommands and a handful of
//! flags; a parser crate would be more dependency than the surface justifies.

use std::process::ExitCode;
use std::time::Instant;

use fpsearch::{fps, Index, IndexBuilder};

const USAGE: &str = "\
fpsearch — Tanimoto similarity search over binary fingerprints

USAGE:
    fpsearch build <input.fps> <output.idx>
    fpsearch query <index.idx> <query-hex> [--threshold T] [--top-k K] [--threads N]
    fpsearch bench <index.idx> [--queries N] [--threshold T] [--top-k K] [--threads N]
    fpsearch info  <index.idx>

INPUT FORMAT (.fps)
    One record per line, `id<TAB>hex`. Lines starting with # are ignored.

OPTIONS
    --threshold T   Minimum Tanimoto to report        (default 0.7)
    --top-k K       Maximum hits to return            (default 10)
    --queries N     Queries to time, drawn from the index itself (default 1000)
    --threads N     Worker threads for the scan     (default: available cores)
";

fn flag_value(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .cloned()
}

fn parse_flag<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> Result<T, String> {
    match flag_value(args, name) {
        None => Ok(default),
        Some(raw) => raw
            .parse()
            .map_err(|_| format!("{name} expects a value, got {raw:?}")),
    }
}

fn build(input: &str, output: &str) -> Result<(), String> {
    let text = std::fs::read_to_string(input).map_err(|e| format!("reading {input}: {e}"))?;
    let started = Instant::now();
    let (records, n_bits) = fps::parse(&text).map_err(|e| format!("{input}: {e}"))?;
    if records.is_empty() {
        return Err(format!("{input} contains no fingerprints"));
    }

    let mut builder = IndexBuilder::new(n_bits);
    let mut refused = 0usize;
    for record in &records {
        // An all-zero fingerprint has no defined similarity to anything. Counting the
        // refusals is the point: a silent skip would leave the index quietly smaller than
        // the input and nothing would say so.
        match builder.push(record.id, &record.fingerprint) {
            Ok(()) => {}
            Err(fpsearch::IndexError::EmptyFingerprint { .. }) => refused += 1,
            Err(e) => return Err(format!("{input}: {e}")),
        }
    }
    let kept = builder.len();
    builder
        .write(output)
        .map_err(|e| format!("writing {output}: {e}"))?;
    let elapsed = started.elapsed();

    println!("built {output}");
    println!("  fingerprints  {kept}");
    println!(
        "  width         {n_bits} bits ({} words)",
        (n_bits as usize).div_ceil(64)
    );
    if refused > 0 {
        println!("  refused       {refused} (no bits set; similarity undefined)");
    }
    println!("  elapsed       {:.2?}", elapsed);
    Ok(())
}

fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

fn query(
    index_path: &str,
    hex: &str,
    threshold: f64,
    k: usize,
    threads: usize,
) -> Result<(), String> {
    let index = Index::open(index_path).map_err(|e| format!("{index_path}: {e}"))?;
    let fingerprint = fps::from_hex(hex.trim(), 1).map_err(|e| format!("query: {e}"))?;

    let started = Instant::now();
    let (hits, stats) = index
        .search_parallel(&fingerprint, threshold, k, threads)
        .map_err(|e| format!("search: {e}"))?;
    let elapsed = started.elapsed();

    for hit in &hits {
        println!("{}\t{:.4}", hit.id, hit.score);
    }
    // To stderr, so piping the hits to another tool does not pick up the commentary.
    eprintln!(
        "{} hits in {:.2?} — examined {} of {} ({:.1}% skipped by the popcount bound)",
        hits.len(),
        elapsed,
        stats.examined,
        stats.total,
        stats.pruned_fraction() * 100.0
    );
    Ok(())
}

fn bench(
    index_path: &str,
    n_queries: usize,
    threshold: f64,
    k: usize,
    threads: usize,
) -> Result<(), String> {
    let index = Index::open(index_path).map_err(|e| format!("{index_path}: {e}"))?;
    if index.is_empty() {
        return Err("index is empty".into());
    }

    // Queries are drawn from the index itself, spread evenly across the popcount ordering
    // so the timing is not dominated by one end of the distribution.
    let step = (index.len() / n_queries.max(1)).max(1);
    let queries: Vec<Vec<u64>> = (0..index.len())
        .step_by(step)
        .take(n_queries)
        .map(|i| index.fingerprint(i).to_vec())
        .collect();

    // One untimed pass, so the measurement reflects steady state rather than the cost of
    // faulting the mapping in for the first time.
    for q in &queries {
        let _ = index.search_parallel(q, threshold, k, threads);
    }

    let mut timings = Vec::with_capacity(queries.len());
    let mut examined_total = 0usize;
    for q in &queries {
        let started = Instant::now();
        let (_, stats) = index
            .search_parallel(q, threshold, k, threads)
            .map_err(|e| format!("search: {e}"))?;
        timings.push(started.elapsed().as_secs_f64() * 1e6);
        examined_total += stats.examined;
    }
    timings.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mean = timings.iter().sum::<f64>() / timings.len() as f64;
    let median = timings[timings.len() / 2];
    let p95 = timings[((timings.len() as f64 * 0.95) as usize).min(timings.len() - 1)];
    let examined_mean = examined_total as f64 / queries.len() as f64;

    println!(
        "index         {} fingerprints, {} bits",
        index.len(),
        index.n_bits()
    );
    println!(
        "queries       {} at threshold {threshold}, top-{k}, {threads} thread(s)",
        timings.len()
    );
    println!("median        {median:.1} µs");
    println!("mean          {mean:.1} µs");
    println!("p95           {p95:.1} µs");
    println!(
        "examined      {examined_mean:.0} of {} per query ({:.1}% skipped)",
        index.len(),
        (1.0 - examined_mean / index.len() as f64) * 100.0
    );
    println!(
        "throughput    {:.2} M comparisons/s over examined candidates",
        examined_mean / mean
    );
    Ok(())
}

fn info(index_path: &str) -> Result<(), String> {
    let index = Index::open(index_path).map_err(|e| format!("{index_path}: {e}"))?;
    let popcounts = index.popcounts();
    println!("path          {index_path}");
    println!("fingerprints  {}", index.len());
    println!(
        "width         {} bits ({} words)",
        index.n_bits(),
        index.n_words()
    );
    println!("metric        {:?}", index.metric());
    if !popcounts.is_empty() {
        let sum: u64 = popcounts.iter().map(|&p| p as u64).sum();
        println!(
            "popcount      min {} median {} max {} mean {:.1}",
            popcounts[0],
            popcounts[popcounts.len() / 2],
            popcounts[popcounts.len() - 1],
            sum as f64 / popcounts.len() as f64
        );
    }
    Ok(())
}

fn run(args: &[String]) -> Result<(), String> {
    let threshold: f64 = parse_flag(args, "--threshold", 0.7)?;
    let k: usize = parse_flag(args, "--top-k", 10)?;
    let threads: usize = parse_flag(args, "--threads", default_threads())?;

    match args.first().map(String::as_str) {
        Some("build") => {
            let input = args.get(1).ok_or("build needs an input .fps file")?;
            let output = args.get(2).ok_or("build needs an output .idx path")?;
            build(input, output)
        }
        Some("query") => {
            let index = args.get(1).ok_or("query needs an index path")?;
            let hex = args.get(2).ok_or("query needs a fingerprint in hex")?;
            query(index, hex, threshold, k, threads)
        }
        Some("bench") => {
            let index = args.get(1).ok_or("bench needs an index path")?;
            let n: usize = parse_flag(args, "--queries", 1000)?;
            bench(index, n, threshold, k, threads)
        }
        Some("info") => info(args.get(1).ok_or("info needs an index path")?),
        Some("--help") | Some("-h") | None => {
            print!("{USAGE}");
            Ok(())
        }
        Some(other) => Err(format!("unknown subcommand {other:?}\n\n{USAGE}")),
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("fpsearch: {message}");
            ExitCode::FAILURE
        }
    }
}
