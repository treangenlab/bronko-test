#!/usr/bin/env python3

import subprocess
import csv
import re
from pathlib import Path

SIZES = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]

## NEED TO CHANGE THE 1000X depending on depth
depth = 30000

DB_DIR = Path("databases/sars")
READ_DIR = Path(f"data/sars_genome_seq_data/") 
OUT_BASE = Path("outputs/sars")

BRONKO_CMD = ["bronko", "call"]
OUTPUT_CSV = f"benchmark_results_call_sars_{depth}.csv"


def run_benchmark(size):
    db_path = DB_DIR / f"{size}.bkdb"
    read_path = READ_DIR / f"{depth}.fastq"
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
        writer.writerow(["number_files", "elapsed_time", "max_rss_kb"])
        writer.writerows(rows)

    print(f"\nResults written to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()