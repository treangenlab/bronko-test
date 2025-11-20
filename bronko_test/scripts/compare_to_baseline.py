import subprocess
import argparse

def build_index(reference_fasta, index_prefix):
    """Build Bowtie2 index from a reference genome FASTA file."""
    command = [
        "bowtie2-build",
        reference_fasta,
        index_prefix
    ]
    try:
        subprocess.run(command, check=True)
        print(f"Index built successfully with prefix {index_prefix}")
    except subprocess.CalledProcessError as e:
        print(f"Error building index: {e}")

def run_bowtie2(reference, reads1, reads2, output_sam, threads=4):
    """Run Bowtie2 to align paired-end reads to a reference genome."""
    command = [
        "bowtie2",
        "-x", reference,
        "-1", reads1,
        "-2", reads2,
        "-S", output_sam,
        "--threads", str(threads)
    ]
    
    try:
        subprocess.run(command, check=True)
        print(f"Alignment complete. Output saved to {output_sam}")
    except subprocess.CalledProcessError as e:
        print(f"Error running Bowtie2: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build Bowtie2 index and align paired-end reads to a reference genome")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome FASTA file")
    parser.add_argument("-1", "--reads1", required=True, help="Path to forward reads (R1)")
    parser.add_argument("-2", "--reads2", required=True, help="Path to reverse reads (R2)")
    parser.add_argument("-o", "--output", required=True, help="Output SAM file")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")
    
    args = parser.parse_args()
    
    prefix = args.reference.strip().split('.fasta')[0]
    
    # Build Bowtie2 index
    build_index(args.reference, prefix)
    
    # Align reads
    run_bowtie2(prefix, args.reads1, args.reads2, args.output, args.threads)