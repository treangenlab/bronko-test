import subprocess
import csv
import re
from pathlib import Path


FASTA_DIR = Path("data/random_genomes")     
OUTPUT_DIR = Path("databases/random_genomes")      
BRONKO_CMD = ["bronko", "build"]   
OUTPUT_CSV = "benchmark_results.csv"

def run_benchmark(fasta_path):
    output_prefix = Path(OUTPUT_DIR) / fasta_path.stem

    cmd = (
        ["/usr/bin/time", "-v"]
        + BRONKO_CMD
        + ["-g", str(fasta_path)]
        + ["-o", str(output_prefix)]
        + ["-t", str(10)]
    )
    
    result = subprocess.run(
        cmd,
        stderr=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        text=True
    )

    stderr = result.stderr

    # Extract metrics
    time_match = re.search(r"Elapsed \(wall clock\) time.*: (.*)", stderr)
    mem_match = re.search(r"Maximum resident set size.*: (\d+)", stderr)

    elapsed = time_match.group(1) if time_match else "NA"
    memory_kb = mem_match.group(1) if mem_match else "NA"

    return elapsed, memory_kb


def main():
    fasta_files = sorted(FASTA_DIR.glob("*.fasta"))

    rows = []

    for fasta in fasta_files:
        print(f"Running: {fasta.name}")
        elapsed, memory = run_benchmark(fasta)
        rows.append([fasta.name, elapsed, memory])
        print(elapsed, memory)

    # Write CSV
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "elapsed_time", "max_rss_kb"])
        writer.writerows(rows)

    print(f"\nResults written to {OUTPUT_CSV}")
    

if __name__ == "__main__":
    main()