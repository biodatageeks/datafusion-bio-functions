use std::fmt::Write;
use std::time::Instant;

/// Simulates the bio-formats CIGAR formatting cost.
/// Formats 19M ops (matching the real BAM) from (len, op_char) pairs to strings.
fn main() {
    // Simulate ~19M reads, mostly simple CIGARs (e.g. "100M", "76M")
    // with some complex ones
    let simple_ops: Vec<(u32, char)> = vec![(100, 'M')];
    let complex_ops: Vec<(u32, char)> = vec![(30, 'M'), (2, 'D'), (50, 'M'), (1, 'I'), (20, 'M')];

    let num_reads = 19_300_000u64;
    let mut buf = String::with_capacity(64);

    let start = Instant::now();
    for i in 0..num_reads {
        buf.clear();
        let ops = if i % 10 == 0 {
            &complex_ops
        } else {
            &simple_ops
        };
        for &(len, op) in ops {
            let _ = write!(buf, "{}{}", len, op);
        }
        // Simulate StringBuilder append (just measure the formatting)
        std::hint::black_box(&buf);
    }
    let elapsed = start.elapsed();

    println!("Reads formatted:    {}", num_reads);
    println!("Format time:        {:.3}s", elapsed.as_secs_f64());
    println!(
        "Per-read:           {:.0}ns",
        elapsed.as_nanos() as f64 / num_reads as f64
    );
}
