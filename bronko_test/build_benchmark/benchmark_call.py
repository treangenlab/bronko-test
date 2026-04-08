#!/usr/bin/env python3

import subprocess
import csv
import re
from pathlib import Path

SIZES = [
    500, 1000, 2500, 5000, 10000, 25000,
    50000, 100000, 250000, 500000,
    1000000, 2500000
]

## NEED TO CHANGE THE 1000X depending on depth
depth = 30000

DB_DIR = Path("databases/random_genomes")
READ_DIR = Path(f"data/random_genome_seq_data/{depth}x") 
OUT_BASE = Path("outputs/random_genome_call")

BRONKO_CMD = ["bronko", "call"]
OUTPUT_CSV = f"benchmark_results_call_{depth}.csv"


def run_benchmark(size):
    db_path = DB_DIR / f"{size}.bkdb"
    read_path = READ_DIR / f"{size}.fastq"
    out_dir = OUT_BASE / str(size)
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = (
        ["/usr/bin/time", "-v"]
        + BRONKO_CMD
        + ["-d", str(db_path)]
        + ["-r", str(read_path)]
        + ["-o", str(out_dir)]
        + ["-t", "10"]
    )

    result = subprocess.run(
        cmd,
        stderr=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        text=True
    )

    stderr = result.stderr

    time_match = re.search(r"Elapsed \(wall clock\) time.*: (.*)", stderr)
    mem_match = re.search(r"Maximum resident set size.*: (\d+)", stderr)

    elapsed = time_match.group(1) if time_match else "NA"
    memory_kb = mem_match.group(1) if mem_match else "NA"

    return elapsed, memory_kb


def main():
    rows = []

    for size in SIZES:
        print(f"Running size: {size}")
        elapsed, memory = run_benchmark(size)
        print("  ", elapsed, memory)
        rows.append([size, elapsed, memory])

    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["genome_size", "elapsed_time", "max_rss_kb"])
        writer.writerows(rows)

    print(f"\nResults written to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()