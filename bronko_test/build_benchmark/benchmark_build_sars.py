import subprocess
import csv
import re
from pathlib import Path


FASTA_DIR = Path("data/sars_genomes")     
OUTPUT_DIR = Path("databases/sars")      
BRONKO_CMD = ["bronko", "build"]   
OUTPUT_CSV = "benchmark_results_sars.csv"

def run_benchmark(fasta_paths, num):
    
    output_prefix = Path(OUTPUT_DIR) / str(num)

    cmd = (
        ["/usr/bin/time", "-v"]
        + BRONKO_CMD
        + ["-g"] + fasta_paths[:num]
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
    
    num_genomes = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]

    for num in num_genomes:
        print(f"Running: {num}")
        elapsed, memory = run_benchmark(fasta_files, num)
        rows.append([str(num), elapsed, memory])
        print(elapsed, memory)

    # Write CSV
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["number_files", "elapsed_time", "max_rss_kb"])
        writer.writerows(rows)

    print(f"\nResults written to {OUTPUT_CSV}")
    

if __name__ == "__main__":
    main()