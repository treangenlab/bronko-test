import os
import subprocess
import pysam
from Bio import SeqIO
from collections import defaultdict
import sys
import csv
from src_py.build import build
from src_py.screen import screen
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd
import time
import math
import pandas as pd
import numpy as np
from scipy.spatial import distance
import pysam
from cyvcf2 import VCF
from collections import namedtuple
from pathlib import Path
import psutil
import threading

def load_tsv(file_path):
    """
    Load the TSV file and return the full DataFrame.
    """
    df = pd.read_csv(file_path, sep='\t')
    return df

def count_nonzero_bases(row, min=0):
    """Count how many of A/C/G/T are > 0"""
    return sum([row[base] > min for base in ["A", "C", "G", "T"]])

def count_af_above_threshold(row, min_af=0.05, min_reads=0):
    """Counts the number of ACGT above the allele threshold"""
    row = row.copy()
    depth = row[["A", "C", "G", "T"]].sum()
    min_depth = min_reads / min_af
    if depth > min_depth:
        row[["A", "C", "G", "T"]] = row[["A", "C", "G", "T"]] / depth
        return sum([row[base] > min_af for base in ["A","C","G","T"]])
    else:
        return 0

def analyze_diversity(df1, df2, output_folder, min_af=0.05, end_filter=True, k=15):
    """
    Compare base diversity between two DataFrames row by row.
    Returns lists of positions and their marker type:
    - ('purple', '.') → both have diversity
    - ('blue', 'x')   → only df1
    - ('red', 'x')    → only df2
    """
    positions = df1["global_index"]
    markers = []
    print(k)
    with open(f'{output_folder}/variants.tsv', "w") as f, open(f'{output_folder}/counts.tsv', "w") as f_full:
        f.write('reference\tindex\tglobal_index\tref\tpresent\tA_iv\tC_iv\tG_iv\tT_iv\tA_bt\tC_bt\tG_bt\tT_bt\n')
        f_full.write('reference\tindex\tglobal_index\tref\tpresent\tA_iv\tC_iv\tG_iv\tT_iv\tA_bt\tC_bt\tG_bt\tT_bt\n')
        for i in range(len(positions)):
            
            d1 = count_af_above_threshold(df1.iloc[i], min_af)
            d2 = count_af_above_threshold(df2.iloc[i], min_af)
            row1 = df1.iloc[i]
            row2 = df2.iloc[i]
            max_seq_length = max(df1[df1["reference"] == row1["reference"]]["index"])
            length_filter = k if end_filter else 0
            
            
            
            if (row1["index"]>length_filter and row1["index"] < max_seq_length-length_filter+1):
                if d1 >= 2 and d2 >= 2:
                    markers.append(('purple', '|', row1["global_index"]))
                    f.write(f"{row1["reference"]}\t{row1["index"]}\t{row1["global_index"]}\t{row1["ref"]}\tBoth\t{row1["A"]}\t{row1["C"]}\t{row1["G"]}\t{row1["T"]}\t{row2["A"]}\t{row2["C"]}\t{row2["G"]}\t{row2["T"]}\n")
                    f_full.write(f"{row1["reference"]}\t{row1["index"]}\t{row1["global_index"]}\t{row1["ref"]}\tBoth\t{row1["A"]}\t{row1["C"]}\t{row1["G"]}\t{row1["T"]}\t{row2["A"]}\t{row2["C"]}\t{row2["G"]}\t{row2["T"]}\n")
                elif d1 >= 2:
                    markers.append(('blue', '|', row1["global_index"]))
                    f.write(f"{row1["reference"]}\t{row1["index"]}\t{row1["global_index"]}\t{row1["ref"]}\tIV\t{row1["A"]}\t{row1["C"]}\t{row1["G"]}\t{row1["T"]}\t{row2["A"]}\t{row2["C"]}\t{row2["G"]}\t{row2["T"]}\n")
                    f_full.write(f"{row1["reference"]}\t{row1["index"]}\t{row1["global_index"]}\t{row1["ref"]}\tIV\t{row1["A"]}\t{row1["C"]}\t{row1["G"]}\t{row1["T"]}\t{row2["A"]}\t{row2["C"]}\t{row2["G"]}\t{row2["T"]}\n")
                elif d2 >= 2:
                    markers.append(('red', '|', row1["global_index"]))
                    f.write(f"{row1["reference"]}\t{row1["index"]}\t{row1["global_index"]}\t{row1["ref"]}\tBT\t{row1["A"]}\t{row1["C"]}\t{row1["G"]}\t{row1["T"]}\t{row2["A"]}\t{row2["C"]}\t{row2["G"]}\t{row2["T"]}\n")
                    f_full.write(f"{row1["reference"]}\t{row1["index"]}\t{row1["global_index"]}\t{row1["ref"]}\tBT\t{row1["A"]}\t{row1["C"]}\t{row1["G"]}\t{row1["T"]}\t{row2["A"]}\t{row2["C"]}\t{row2["G"]}\t{row2["T"]}\n")
                else:
                    f_full.write(f"{row1["reference"]}\t{row1["index"]}\t{row1["global_index"]}\t{row1["ref"]}\tNeither\t{row1["A"]}\t{row1["C"]}\t{row1["G"]}\t{row1["T"]}\t{row2["A"]}\t{row2["C"]}\t{row2["G"]}\t{row2["T"]}\n")

    return markers


def concatenate_segments(df):
    """
    Convert a dataframe with multiple segments into a unified position space.
    Returns the updated dataframe and a dictionary of segment offsets.
    """
    df_sorted = df.sort_values(by=['reference', 'index'])
    offsets = {}
    cumulative = 0
    new_positions = []

    for segment, group in df_sorted.groupby('reference'):
        offsets[segment] = cumulative
        new_positions.extend(group['index'] + cumulative)
        cumulative += group['index'].max()

    df_sorted = df_sorted.copy()
    df_sorted['global_index'] = new_positions
    return df_sorted, offsets


def plot_depth(file1, file2, out, label1="IV tool", label2="Bowtie2", k=15, min_af=0.05, end_filter=True):
    df1 = load_tsv(file1) ## my tool
    df2 = load_tsv(file2) ## other tool for comparison (Bowtie2 for now)
    
    df1, offsets = concatenate_segments(df1)
    df2, _ = concatenate_segments(df2)

    positions = df1["global_index"]
    depth1 = df1[["A", "C", "G", "T"]].sum(axis=1)
    depth2 = df2[["A", "C", "G", "T"]].sum(axis=1)

    markers = analyze_diversity(df1, df2, k=k, min_af=min_af, output_folder=out, end_filter=end_filter)

    fig, axs = plt.subplots(3, 1, figsize=(14, 10), 
                            gridspec_kw={'height_ratios': [4, 1, 2]}, 
                            sharex=True)

    ax1, ax2, ax3 = axs

    # Plot depth
    ax1.plot(positions, depth1, label=label1, color='blue', alpha=0.7)
    ax1.plot(positions, depth2, label=label2, color='red', alpha=0.7)
    ax1.set_ylabel("Coverage")
    ax1.set_title("Read Depth and Base Diversity")
    ax1.legend()
    ax1.grid(axis='y')

    ivtool_iv = 0
    bt_iv = 0
    both_iv = 0

    # Diversity markers
    for color, marker, pos in markers:
        if color == 'purple':
            y = 3
            both_iv += 1
        elif color == 'blue':
            y = 2
            ivtool_iv += 1
        elif color == 'red':
            y = 1
            bt_iv += 1
        ax2.plot(pos, y, marker=marker, color=color, markersize=10)

    ax2.set_yticks([1, 2, 3])
    ax2.set_yticklabels([f"{label2} only (n={bt_iv})", f"{label1} only (n={ivtool_iv})", f"Both (n={both_iv})"])
    ax2.set_ylabel("Diversity")
    ax2.set_ylim([0.5, 3.5])
    ax2.grid(True, axis='x', linestyle=':', alpha=0.3)
    ax2.xaxis.set_minor_locator(MultipleLocator(100))


    # Minor allele frequencies (purple only)
    diff_positions = []
    freq_diffs = []
    diff_colors = []
    diff_shapes = []
    
    for color, marker, pos in markers:
        row1 = df1[df1['global_index'] == pos][["A", "C", "G", "T"]].values[0]
        row2 = df2[df2['global_index'] == pos][["A", "C", "G", "T"]].values[0]

        sum1 = row1.sum()
        sum2 = row2.sum()

        if sum1 == 0 or sum2 == 0:
            continue ## TODO: maybe change this to plot in gray or something

        freq1 = row1 / sum1
        freq2 = row2 / sum2

        # Use second-highest (minor) allele from each sample
        minor_idx1 = np.argsort(freq1)[-2]
        minor_idx2 = np.argsort(freq2)[-2]

        # Match alleles: use same minor if they agree, otherwise treat separately
        if minor_idx1 == minor_idx2:
            diff = freq1[minor_idx1] - freq2[minor_idx2]
            shape = 'o'
        else:
            # Take average diff for simplicity, or choose more consistent logic if needed
            diff1 = freq1[minor_idx1] - freq2[minor_idx1]
            diff2 = freq1[minor_idx2] - freq2[minor_idx2]
            diff = np.nanmean([diff1, diff2])  # skip if both are nan
            shape = '^'

        diff_positions.append(pos)
        freq_diffs.append(diff)
        diff_colors.append(color)
        diff_shapes.append(shape)

    # Plot frequency differences
    for x, y, c, m in zip(diff_positions, freq_diffs, diff_colors, diff_shapes):
        # if c=="purple":
        ax3.plot(x, y, color=c, marker=m)

    ax3.axhline(0, linestyle='--', color='gray', alpha=0.5)
    ax3.set_ylabel(f"Minor Allele Frequency Difference\n({label1} − {label2})")
    ax3.set_xlabel("Reference Position")
    padding = 0.05  # 5% padding
    ymin = min(freq_diffs)
    ymax = max(freq_diffs)
    yrange = ymax - ymin
    ax3.set_ylim([ymin - padding * yrange, ymax + padding * yrange])
    
        
    for segment, offset in offsets.items():
        ax3.axvline(x=offset, color='black', linestyle='--', alpha=0.8)
        ax2.axvline(x=offset, color='black', linestyle='--', alpha=0.8)
        ax1.axvline(x=offset, color='black', linestyle='--', alpha=0.8)

        mid = offset + df1[df1['reference'] == segment]['index'].max() / 2
        ax3.text(mid, ax3.get_ylim()[1], segment, ha='center', va='bottom', fontsize=8)
    
    ax3.set_title("Minor Allele Frequency Differences at Variant Sites")
    ax3.grid(axis='y')

    plt.tight_layout()
    plt.savefig(f'{out}/comparison.png')
    return ivtool_iv, bt_iv, both_iv

def run_cmd_shell_mem(cmd, label):
    """Run a command with timing and peak memory usage tracking."""
    print(f"\n▶️  Running: {label}")
    print(f"   Command: {' '.join(cmd)}")

    start = time.time()
    process = subprocess.Popen(cmd)
    ps = psutil.Process(process.pid)

    peak_mem = 0
    while process.poll() is None:
        try:
            mem = ps.memory_info().rss / (1024**2)  # MB
            peak_mem = max(peak_mem, mem)
        except psutil.NoSuchProcess:
            break
        time.sleep(0.1)

    process.wait()
    end = time.time()
    total_time = end-start

    print(f"✅  {label} completed in {end - start:.2f} seconds | Peak memory: {peak_mem:.2f} MB")
    return total_time, peak_mem

def build_bowtie2_index(ref_fasta, index_prefix, threads = 10):
    """
    Build Bowtie2 index from reference FASTA.
    """
    cmd = ["bowtie2-build", ref_fasta, index_prefix, "--threads", str(threads)]
    t, m = run_cmd_shell_mem(cmd, "Bowtie2 Build Index")
    return t, m

def align_with_bowtie2(fastq_file, index_prefix, output_sam, threads = 10):
    """
    Align reads using Bowtie2 and output a SAM file.
    """
    cmd = ["bowtie2", "-x", str(index_prefix), "-U", str(fastq_file), "-S", str(output_sam), "-p", str(threads)]
    t, m = run_cmd_shell_mem(cmd, "Bowtie2 Align to Index SE")
    return t, m


def align_paired_with_bowtie2(fastq_r1, fastq_r2, index_prefix, output_sam, threads=10):
    """
    Align paired-end reads using Bowtie2 and output a SAM file.
    """
    cmd = ["bowtie2", "-x", str(index_prefix), "-1", str(fastq_r1), "-2", str(fastq_r2), "-S", str(output_sam), "-p", str(threads)]
    t, m = run_cmd_shell_mem(cmd, "Bowtie2 Build Index PE")
    return t, m



def parse_alignment(sam_file, ref_fasta):
    """
    Parse SAM file and count bases per position, separated by strand.
    Returns a dict of segment -> list of dicts with counts per position.
    """
    samfile = pysam.AlignmentFile(sam_file, "r")

    # Load reference sequences
    ref_seqs = {record.id: str(record.seq).upper() for record in SeqIO.parse(ref_fasta, "fasta")}

    # Initialize counts
    base_counts = {
        segment: [defaultdict(int) for _ in range(len(seq))]
        for segment, seq in ref_seqs.items()
    }

    mapped_reads = 0
    unmapped_reads = 0

    for read in samfile.fetch(until_eof=True):
        if read.is_unmapped:
            unmapped_reads += 1
            continue

        mapped_reads += 1
        segment = samfile.get_reference_name(read.reference_id)
        is_reverse = read.is_reverse

        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            if qpos is not None and rpos is not None:
                base = read.query_sequence[qpos].upper()
                if base in "ACGT":
                    key = base.lower() if is_reverse else base
                    base_counts[segment][rpos][key] += 1

    return base_counts, ref_seqs, mapped_reads, unmapped_reads

def write_counts_to_tsv(base_counts, ref_seqs, output_tsv):
    """
    Write counts to TSV with the format:
    reference, index, ref, A, C, G, T, a, c, g, t
    """
    with open(output_tsv, "w", newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["reference", "index", "ref", "A", "C", "G", "T", "a", "c", "g", "t"])

        for segment, counts_list in base_counts.items():
            ref_seq = ref_seqs[segment]
            for i, counts in enumerate(counts_list):
                ref_base = ref_seq[i] if i < len(ref_seq) else "N"
                row = [
                    segment,
                    i + 1,
                    ref_base,
                    counts.get("A", 0),
                    counts.get("C", 0),
                    counts.get("G", 0),
                    counts.get("T", 0),
                    counts.get("a", 0),
                    counts.get("c", 0),
                    counts.get("g", 0),
                    counts.get("t", 0),
                ]
                writer.writerow(row)
                

def run_cmd(cmd, label):
    """Wrapper for a subprocess to run an individual command with some additional output and timing information"""
    print(f"\n▶️  Running: {label}")
    print(f"   Command: {' '.join(cmd)}")
    start = time.time()
    subprocess.run(cmd, check=True)
    end = time.time()
    print(f"✅  {label} completed in {end - start:.2f} seconds")
    
def run_cmd_mem(cmd, label):
    """Run a command with timing and peak memory usage tracking."""
    print(f"\n▶️  Running: {label}")
    print(f"   Command: {' '.join(cmd)}")

    start = time.time()
    process = subprocess.Popen(cmd)
    ps = psutil.Process(process.pid)

    peak_mem = 0
    while process.poll() is None:
        try:
            mem = ps.memory_info().rss / (1024**2)  # MB
            peak_mem = max(peak_mem, mem)
        except psutil.NoSuchProcess:
            break
        time.sleep(0.1)

    process.wait()
    end = time.time()
    total_time = end-start

    print(f"✅  {label} completed in {end - start:.2f} seconds | Peak memory: {peak_mem:.2f} MB")
    return total_time, peak_mem
    
def run_cmd_with_pipe(cmd1, cmd2, label):
    """Wrapper for a subprocess to run a command piped into the other with additional output and timing information"""
    print(f"\n▶️  Running: {label}")
    print(f"   Command: {' '.join(cmd1)} | {' '.join(cmd2)}")
    start = time.time()
    
    proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    proc2 = subprocess.run(cmd2, stdin=proc1.stdout, check=True)
    
    proc1.stdout.close()
    proc1.wait()

    end = time.time()
    print(f"✅  {label} completed in {end - start:.2f} seconds")

    
def run_cmd_with_pipe_and_mem(cmd1, cmd2, label):
    """Run cmd1 | cmd2 with timing and peak memory tracking for both (no hang)."""
    print(f"\n▶️  Running: {label}")
    print(f"   Command: {' '.join(cmd1)} | {' '.join(cmd2)}")

    start = time.time()
    proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    ps1 = psutil.Process(proc1.pid)
    peak_mem = {"value": 0}

    def monitor_memory():
        ps2 = None
        while proc1.poll() is None or (ps2 and ps2.is_running()):
            try:
                if ps1.is_running():
                    mem = ps1.memory_info().rss / (1024**2)
                    peak_mem["value"] = max(peak_mem["value"], mem)
            except psutil.NoSuchProcess:
                pass
            if ps2 and ps2.is_running():
                try:
                    mem = ps2.memory_info().rss / (1024**2)
                    peak_mem["value"] = max(peak_mem["value"], mem)
                except psutil.NoSuchProcess:
                    pass
            time.sleep(0.1)

    # Start memory monitoring thread
    monitor_thread = threading.Thread(target=monitor_memory, daemon=True)
    monitor_thread.start()

    # Run the second command safely (subprocess.run handles the pipe properly)
    proc2 = subprocess.run(cmd2, stdin=proc1.stdout, check=True)
    proc1.stdout.close()
    proc1.wait()

    monitor_thread.join(timeout=0.2)

    end = time.time()
    total_time = end - start
    print(f"✅  {label} completed in {total_time:.2f} seconds | Peak memory: {peak_mem['value']:.2f} MB")
    return total_time, peak_mem["value"]
    
def run_ivar_pipeline(sam_file, ref_fasta, output_prefix, threads=10, min_af=0.01):
    """
    Run iVar iSNV detection pipeline using samtools mpileup piped to ivar variants.

    Args:
        sorted_bam: path to input sorted BAM file
        ref_fasta: path to reference FASTA file
        output_prefix: prefix for iVar output files

    Returns:
        Path to final iVar output TSV file
    """
    bam = f"{output_prefix}.bam"
    sorted_bam = f"{output_prefix}.sorted.bam"
    ivar_output = f"{output_prefix}_ivar.tsv"

    times = []
    mems = []
    # Step 1: Index reference FASTA (required by samtools)
    t, m = run_cmd_mem(["samtools", "faidx", str(ref_fasta)], "Index Reference FASTA")
    times.append(t)
    mems.append(m)
    
    # Step 2: Convert SAM to BAM
    t, m = run_cmd_mem(["samtools", "view", "-bS", sam_file, "-o", bam, "-@", str(threads)], "SAM to BAM")
    times.append(t)
    mems.append(m)

    # Step 3: Sort BAM
    t, m = run_cmd_mem(["samtools", "sort", bam, "-o", sorted_bam, "-@", str(threads)], "Sort BAM")
    times.append(t)
    mems.append(m)

    # Step 4: Index sorted BAM
    t, m = run_cmd_mem(["samtools", "index", sorted_bam, "-@", str(threads)], "Index Sorted BAM")
    times.append(t)
    mems.append(m)
    samtools_mem = max(mems)
    samtools_time = sum(times)
    
    # Step 5: Run mpileup + ivar
    samtools_cmd = [
        "samtools", "mpileup",
        "-aa", "-A", "-d", "0", "-B", "-Q", "0",
        "--reference", str(ref_fasta),
        str(sorted_bam)
    ]

    ivar_cmd = [
        "ivar", "variants",
        "-p", f"{output_prefix}_ivar", 
        "-q", "20",
        "-t", str(min_af),
        "-r", str(ref_fasta)
    ]

    # run_cmd_with_pipe(samtools_cmd, ivar_cmd, "iVar Variant Calling (with mpileup included)")
    t, m = run_cmd_with_pipe_and_mem(samtools_cmd, ivar_cmd, "iVar Variant Calling (with mpileup included)")
    times.append(t)
    mems.append(m)
    
    max_mem = m
    total_time = sum(times)

    print(f"\n🎉 Finished! Final iVar output: {ivar_output}")
    return ivar_output, total_time, max_mem, samtools_time, samtools_mem

def run_lofreq_pipeline(sam_file, ref_fasta, output_prefix, threads=4):
    """
    Run LoFreq iSNV pipeline from SAM to VCF (without indel realignment).
    
    Steps:
        1. Convert SAM to BAM
        2. Sort BAM
        3. Index BAM
        4. Call variants with LoFreq
    
    Args:
        sam_file: path to input SAM file
        ref_fasta: path to reference FASTA file
        output_prefix: prefix for BAM and VCF files
        threads: number of threads to use in parallel call
    
    Returns:
        Path to final VCF file
    """

    bam = f"{output_prefix}.bam"
    sorted_bam = f"{output_prefix}.sorted.bam"
    vcf = f"{output_prefix}.vcf"
    
    times = []
    mems = []
    
    # Step 1: Index the fasta file
    t, m = run_cmd_mem(["samtools", "faidx", ref_fasta, "-@", str(threads)], "Index Reference FASTA")
    times.append(t)
    mems.append(m)
    
    # Step 2: Convert SAM to BAM
    t, m = run_cmd_mem(["samtools", "view", "-bS", sam_file, "-o", bam, "-@", str(threads)], "SAM to BAM")
    times.append(t)
    mems.append(m)
    
    # Step 3: Sort BAM
    t, m = run_cmd_mem(["samtools", "sort", bam, "-o", sorted_bam, "-@", str(threads)], "Sort BAM")
    times.append(t)
    mems.append(m)
    
    # Step 4: Index sorted BAM
    t, m = run_cmd_mem(["samtools", "index", sorted_bam, "-@", str(threads)], "Index Sorted BAM")
    times.append(t)
    mems.append(m)
    samtools_mem = max(mems)
    samtools_time = sum(times)
    
    # Step 5: LoFreq variant calling
    if os.path.exists(vcf):
        run_cmd(["rm", vcf], "Removing existing VCF file")
        
    t, m = run_cmd_mem([
        "lofreq", "call-parallel", "--pp-threads", str(threads),
        "-f", ref_fasta,
        "-o", vcf,
        sorted_bam
    ], "LoFreq Variant Calling")
    times.append(t)
    mems.append(m)
    
    total_time = sum(times)
    max_mem = m

    print(f"\n🎉 Finished! Final VCF: {vcf}")
    return vcf, total_time, max_mem, samtools_time, samtools_mem

            
Variant = namedtuple("Variant", ["chrom", "pos", "ref", "alt"])

def get_segment_lengths(vcf: VCF):
    """
    Extract contig lengths from VCF header (predicted VCF only).
    Returns dict: {chrom: length}
    """
    contig_lengths = {}
    for line in vcf.raw_header.splitlines():
        if line.startswith("##contig="):
            # Example: ##contig=<ID=PV249936.1,length=7555>
            parts = line.strip().lstrip("##contig=<").rstrip(">").split(",")
            kv = dict(p.split("=") for p in parts)
            contig_lengths[kv["ID"]] = int(kv["length"])
    return contig_lengths

def parse_vcf_variants(vcf_path, af_cutoff=0.01, edge_exclude_k=0, segment_lengths=None):
    """Parse a VCF file into (minor, major) variant sets with optional AF cutoff and edge exclusion."""
    vcf = VCF(vcf_path)
    variants_minor = set()
    variants_major = set()

    if segment_lengths is None:
        segment_lengths = get_segment_lengths(vcf)

    for rec in vcf:
        chrom = rec.CHROM
        pos = rec.POS
        alts = rec.ALT or []

        for i, alt in enumerate(alts):
            if alt is None:
                continue

            af = rec.INFO.get("AF")
            af_val = (
                af[i] if isinstance(af, (list, tuple)) and i < len(af)
                else af[0] if isinstance(af, (list, tuple))
                else af if af is not None
                else 1.0
            )

            if af_val < af_cutoff:
                continue

            if edge_exclude_k > 0 and chrom in segment_lengths:
                if pos <= edge_exclude_k or pos >= (segment_lengths[chrom] - edge_exclude_k + 1):
                    continue

            variant = Variant(chrom, pos, rec.REF, alt)
            if af_val <= 0.5:
                variants_minor.add(variant)
            else:
                variants_major.add(variant)

    return variants_minor, variants_major


def parse_ivar_variants(tsv_path, af_cutoff=0.01, edge_exclude_k=0, segment_lengths=None):
    """Parse an iVar TSV file into (minor, major) variant sets."""
    variants_minor = set()
    variants_major = set()

    with open(tsv_path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                chrom = row["REGION"]
                pos = int(row["POS"])
                ref = row["REF"]
                alt = row["ALT"]
                af = float(row["ALT_FREQ"])
            except (KeyError, ValueError, TypeError):
                continue

            if af < af_cutoff:
                continue

            if '+' in alt or '-' in alt:
                continue

            if edge_exclude_k > 0 and segment_lengths and chrom in segment_lengths:
                seg_len = segment_lengths[chrom]
                if pos <= edge_exclude_k or pos >= (seg_len - edge_exclude_k + 1):
                    continue

            variant = Variant(chrom, pos, ref, alt)
            if af <= 0.5:
                variants_minor.add(variant)
            else:
                variants_major.add(variant)

    return variants_minor, variants_major

def load_counts_tsv(tsv_file):
    """
    Loads a TSV file with reference base and counts.
    Returns a dict of (reference, index) -> row_dict
    """
    data = {}
    with open(tsv_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            key = (row["reference"], int(row["index"]))
            data[key] = row
    return data

def build_variant_table_rows(
    bronko_counts_file,
    bowtie2_counts_file,
    bronko_calls,
    lofreq_calls,
    ivar_calls,
    output_file=None  # Optional: write TSV
):
    """
    Build rows of variant info from base count tables and sets of variant calls.
    Returns a list of dicts (each row), and optionally writes to a TSV.
    """
    bronko_counts = load_counts_tsv(bronko_counts_file)
    bowtie2_counts = load_counts_tsv(bowtie2_counts_file)

    all_positions = set((v.chrom, v.pos, v.ref, v.alt) for v in bronko_calls) \
                  | set((v.chrom, v.pos, v.ref, v.alt) for v in lofreq_calls) \
                  | set((v.chrom, v.pos, v.ref, v.alt) for v in ivar_calls)
                  

    rows = []

    for chrom, idx, ref, alt in sorted(all_positions):
        tools = []
        if any((v.chrom, v.pos, v.ref, v.alt) == (chrom, idx, ref, alt) for v in bronko_calls): tools.append("bronko")
        if any((v.chrom, v.pos, v.ref, v.alt) == (chrom, idx, ref, alt) for v in lofreq_calls): tools.append("lofreq")
        if any((v.chrom, v.pos, v.ref, v.alt) == (chrom, idx, ref, alt) for v in ivar_calls): tools.append("ivar")

        tools_str = ",".join(sorted(tools))

        ref_base = (
            bowtie2_counts.get((chrom, idx), {}).get("ref")
            or bronko_counts.get((chrom, idx), {}).get("ref")
            or "N"
        )

        def get_counts(d):
            row = d.get((chrom, idx), {})
            return {
                "A": row.get("A", "0"), "C": row.get("C", "0"),
                "G": row.get("G", "0"), "T": row.get("T", "0"),
                "a": row.get("a", "0"), "c": row.get("c", "0"),
                "g": row.get("g", "0"), "t": row.get("t", "0"),
            }

        bt2 = get_counts(bowtie2_counts)
        bronko = get_counts(bronko_counts)

        row = {
            "reference": chrom,
            "index": idx,
            "ref": ref_base,
            "alt": alt,
            "tools": tools_str,
            "A_bt2": bt2["A"], "C_bt2": bt2["C"], "G_bt2": bt2["G"], "T_bt2": bt2["T"],
            "a_bt2": bt2["a"], "c_bt2": bt2["c"], "g_bt2": bt2["g"], "t_bt2": bt2["t"],
            "A_bronko": bronko["A"], "C_bronko": bronko["C"], "G_bronko": bronko["G"], "T_bronko": bronko["T"],
            "a_bronko": bronko["a"], "c_bronko": bronko["c"], "g_bronko": bronko["g"], "t_bronko": bronko["t"]
        }

        rows.append(row)


    if len(rows) > 0 and output_file:
        fieldnames = list(rows[0].keys())
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    return rows

def print_stats(pred_variants, truth_variants, title='Minor Variant Comparison Summary'):
    """Prints summary statistics of two variants with a title

    Args:
        pred_variants (_type_): _description_
        truth_variants (_type_): _description_

    Returns:
        _type_: _description_
    """
    tp = pred_variants & truth_variants
    fp = pred_variants - truth_variants
    fn = truth_variants - pred_variants

    precision = len(tp) / (len(tp) + len(fp)) if (tp or fp) else 1.0
    recall = len(tp) / (len(tp) + len(fn)) if (tp or fn) else 1.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) else 1.0
    
    print(f"\n📊 {title}")
    print(f"True Positives (TP): {len(tp)}")
    print(f"False Positives (FP): {len(fp)}")
    print(f"False Negatives (FN): {len(fn)}")
    print(f"Precision:  {precision:.4f}")
    print(f"Recall:     {recall:.4f}")
    print(f"F1 Score:   {f1:.4f}")

    return {
        "TP": tp, "FP": fp, "FN": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1
    }

def compare_variant_sets(bronko_vcf, lofreq_vcf, ivar_tsv,
                         af_cutoff=0.01,
                         edge_exclude_k=0):
    """
    Compares predicted VCF against ground truth (LoFreq) VCF.
    Filters predicted variants near edges (if segment lengths available).
    """
    segment_lengths = get_segment_lengths(VCF(bronko_vcf))

    bronko_minor, bronko_major     = parse_vcf_variants(bronko_vcf, af_cutoff, edge_exclude_k, segment_lengths)
    lofreq_minor, lofreq_major = parse_vcf_variants(lofreq_vcf, af_cutoff, edge_exclude_k, segment_lengths)
    ivar_minor, ivar_major     = parse_ivar_variants(ivar_tsv, af_cutoff, edge_exclude_k, segment_lengths)

    # Shared and unique variants (minor only, could parameterize)
    shared_all     = bronko_minor & lofreq_minor & ivar_minor
    shared_lofreq  = bronko_minor & lofreq_minor
    shared_ivar    = bronko_minor & ivar_minor
    unique_bronko    = bronko_minor - (lofreq_minor | ivar_minor)
    unique_lofreq  = lofreq_minor - (bronko_minor | ivar_minor)
    unique_ivar    = ivar_minor - (bronko_minor | lofreq_minor)

    shared_all_mj     = bronko_major & lofreq_major & ivar_major
    shared_lofreq_mj  = bronko_major & lofreq_major
    shared_ivar_mj    = bronko_major & ivar_major
    unique_bronko_mj    = bronko_major - (lofreq_major | ivar_major)
    unique_lofreq_mj  = lofreq_major - (bronko_major | ivar_major)
    unique_ivar_mj    = ivar_major - (bronko_major | lofreq_major)
    
    print("\n🔍 Minor Variant Overlap Summary")
    print(f"Shared by all three:          {len(shared_all)}")
    print(f"Shared by your tool + LoFreq: {len(shared_lofreq)}")
    print(f"Shared by your tool + iVar:   {len(shared_ivar)}")
    print(f"Unique to your tool:          {len(unique_bronko)}")
    print(f"Unique to LoFreq:             {len(unique_lofreq)}")
    print(f"Unique to iVar:               {len(unique_ivar)}")
    
    print("\n🔍 Major Variant Overlap Summary")
    print(f"Shared by all three:          {len(shared_all_mj)}")
    print(f"Shared by your tool + LoFreq: {len(shared_lofreq_mj)}")
    print(f"Shared by your tool + iVar:   {len(shared_ivar_mj)}")
    print(f"Unique to your tool:          {len(unique_bronko_mj)}")
    print(f"Unique to LoFreq:             {len(unique_lofreq_mj)}")
    print(f"Unique to iVar:               {len(unique_ivar_mj)}")

    # Pairwise stats
    stats_vs_lofreq_mj = print_stats(bronko_major, lofreq_major, title="Bronko vs. LoFreq (Major)")
    stats_vs_ivar_mj   = print_stats(bronko_major, ivar_major, title="Bronko vs. iVar (Major)")
    stats_vs_lofreq = print_stats(bronko_minor, lofreq_minor, title="Bronko vs. LoFreq (Minor)")
    stats_vs_ivar   = print_stats(bronko_minor, ivar_minor, title="Bronko vs. iVar (Minor)")

    return {
        "minor": {
            "bronko": bronko_minor,
            "lofreq": lofreq_minor,
            "ivar": ivar_minor,
            "shared_all": shared_all,
            "shared_lofreq": shared_lofreq,
            "shared_ivar": shared_ivar,
            "unique_bronko": unique_bronko,
            "unique_lofreq": unique_lofreq,
            "unique_ivar": unique_ivar,
            "bronko_vs_lofreq": stats_vs_lofreq,
            "bronko_vs_ivar": stats_vs_ivar
        },
        "major": {
            "bronko": bronko_major,
            "lofreq": lofreq_major,
            "ivar": ivar_major,
            "shared_all": shared_all_mj,
            "shared_lofreq": shared_lofreq_mj,
            "shared_ivar": shared_ivar_mj,
            "unique_bronko": unique_bronko_mj,
            "unique_lofreq": unique_lofreq_mj,
            "unique_ivar": unique_ivar_mj,
            "bronko_vs_lofreq": stats_vs_lofreq_mj,
            "bronko_vs_ivar": stats_vs_ivar_mj
        }
    }
    
def clean_sample_id(filename):
    """
    Removes known FASTQ extensions from a file path to extract the sample ID.
    """
    base = os.path.basename(filename)
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if base.endswith(ext):
            return base[:-len(ext)]
    return os.path.splitext(base)[0]
    
def bench(fastq_files, r1s, r2s, ref_fasta, output_folder, k=19, min_af=0.03, re_run=True, threads=10):
        
    ## collect all the samples
    all_samples = []
    if fastq_files:
        for se in fastq_files:
            sample_id = clean_sample_id(se)
            all_samples.append((sample_id, se, None))

    if r1s and r2s:
        for r1, r2 in zip(r1s, r2s):
            sample_id = clean_sample_id(r1)
            all_samples.append((sample_id, r1, r2))
            
    os.makedirs(output_folder, exist_ok=True)

    overview = f"{output_folder}/overview.tsv"
    minor_counts = f"{output_folder}/minor_variants.tsv"
    major_counts = f"{output_folder}/major_variants.tsv"
    
    with open(overview, "w") as f, open(major_counts, "w") as f_mj, open(minor_counts, "w") as f_mn:
        
        overview_headers = [
            "Sample",
            "BT2_Time", "BT2_LoFreq_Time", "BT2_iVar_Time", "Samtools_Time", "Bronko_Time",
            "BT2_Mem", "BT2_LoFreq_Mem", "BT2_iVar_Mem", "Samtools_Mem", "Bronko_Mem",
            # LoFreq comparison (Minor)
            "Minor_TP_LoFreq", "Minor_FP_LoFreq", "Minor_FN_LoFreq",
            "Minor_Precision_LoFreq", "Minor_Recall_LoFreq", "Minor_F1_LoFreq",
            # LoFreq comparison (Major)
            "Major_TP_LoFreq", "Major_FP_LoFreq", "Major_FN_LoFreq",
            "Major_Precision_LoFreq", "Major_Recall_LoFreq", "Major_F1_LoFreq",
            # iVar comparison (Minor)
            "Minor_TP_iVar", "Minor_FP_iVar", "Minor_FN_iVar",
            "Minor_Precision_iVar", "Minor_Recall_iVar", "Minor_F1_iVar",
            # iVar comparison (Major)
            "Major_TP_iVar", "Major_FP_iVar", "Major_FN_iVar",
            "Major_Precision_iVar", "Major_Recall_iVar", "Major_F1_iVar",
            # Shared/unique (Minor)
            "Minor_Shared_All", "Minor_Shared_LoFreq", "Minor_Shared_iVar",
            "Minor_Unique_Bronko", "Minor_Unique_LoFreq", "Minor_Unique_iVar",
            # Shared/unique (Major)
            "Major_Shared_All", "Major_Shared_LoFreq", "Major_Shared_iVar",
            "Major_Unique_Bronko", "Major_Unique_LoFreq", "Major_Unique_iVar"
        ]
        
        counts_headers = [
            "Sample", "reference", "index", "ref", "alt", "tools", "A_bt2", "C_bt2", "G_bt2", "T_bt2", "a_bt2", "c_bt2", "g_bt2", "t_bt2",
            "A_bronko", "C_bronko", "G_bronko", "T_bronko", "a_bronko", "c_bronko", "g_bronko", "t_bronko"
        ]
        
        f.write("\t".join(overview_headers) + "\n")
        f_mj.write("\t".join(counts_headers) + "\n")
        f_mn.write("\t".join(counts_headers) + "\n")
        
        for sample_id, r1, r2 in all_samples:
            
            bench_output_folder = f"{output_folder}/{sample_id}_k{k}"

            result, time_info, mem_info, minor_info, major_info = bench_individual(r1, ref_fasta, bench_output_folder, r1_file=r1 if r2 else None, r2_file=r2, k=k, min_af=min_af, re_run_bench=re_run, threads=threads)

            # Extract stats
            minor_lofreq = result["minor"]["bronko_vs_lofreq"]
            major_lofreq = result["major"]["bronko_vs_lofreq"]
            minor_ivar = result["minor"]["bronko_vs_ivar"]
            major_ivar = result["major"]["bronko_vs_ivar"]

            row = [
                sample_id,
                time_info["bowtie2"],
                time_info["bowtie2+lofreq"],
                time_info["bowtie2+ivar"],
                time_info["samtools"],
                time_info["bronko"],
                
                mem_info["bowtie2"],
                mem_info["bowtie2+lofreq"],
                mem_info["bowtie2+ivar"],
                mem_info["samtools"],
                mem_info["bronko"],

                # Minor (LoFreq)
                len(minor_lofreq["TP"]),
                len(minor_lofreq["FP"]),
                len(minor_lofreq["FN"]),
                round(minor_lofreq["precision"], 4),
                round(minor_lofreq["recall"], 4),
                round(minor_lofreq["f1"], 4),

                # Major (LoFreq)
                len(major_lofreq["TP"]),
                len(major_lofreq["FP"]),
                len(major_lofreq["FN"]),
                round(major_lofreq["precision"], 4),
                round(major_lofreq["recall"], 4),
                round(major_lofreq["f1"], 4),

                # Minor (iVar)
                len(minor_ivar["TP"]),
                len(minor_ivar["FP"]),
                len(minor_ivar["FN"]),
                round(minor_ivar["precision"], 4),
                round(minor_ivar["recall"], 4),
                round(minor_ivar["f1"], 4),

                # Major (iVar)
                len(major_ivar["TP"]),
                len(major_ivar["FP"]),
                len(major_ivar["FN"]),
                round(major_ivar["precision"], 4),
                round(major_ivar["recall"], 4),
                round(major_ivar["f1"], 4),

                # Shared/unique minor
                len(result["minor"]["shared_all"]),
                len(result["minor"]["shared_lofreq"]),
                len(result["minor"]["shared_ivar"]),
                len(result["minor"]["unique_bronko"]),
                len(result["minor"]["unique_lofreq"]),
                len(result["minor"]["unique_ivar"]),

                # Shared/unique major
                len(result["major"]["shared_all"]),
                len(result["major"]["shared_lofreq"]),
                len(result["major"]["shared_ivar"]),
                len(result["major"]["unique_bronko"]),
                len(result["major"]["unique_lofreq"]),
                len(result["major"]["unique_ivar"]),
            ]

            f.write("\t".join(map(str, row)) + "\n")
            
            for row in major_info:
                values = [sample_id] + [str(row.get(col, "")) for col in counts_headers[1:]]
                f_mj.write("\t".join(values) + "\n")

            for row in minor_info:
                values = [sample_id] + [str(row.get(col, "")) for col in counts_headers[1:]]
                f_mn.write("\t".join(values) + "\n")
         
        
            
def bench_individual(fastq_file, ref_fasta, output_folder, r1_file=None, r2_file=None, k=15, min_af=0.01, re_run_bench=False, threads=30):
    
    ## Bowtie2 based benchmark
    print('RUNNING BOWTIE2 BENCHMARK')
    bt_start_time = time.time()
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    
    index_prefix = f"{output_folder}/bt_ref_index"
    sam_file = f"{output_folder}/bt_bench.sam"
    if r1_file and r2_file:
        base_name = clean_sample_id(r1_file)
    else:
        base_name = clean_sample_id(fastq_file)

    if re_run_bench:
        ## build bowtie2 index and run bowtie2
        bowtie_build_time, bowtie_build_mem = build_bowtie2_index(ref_fasta, index_prefix, threads = threads)
        if r1_file and r2_file:
            bowtie_align_time, bowtie_align_mem = align_paired_with_bowtie2(r1_file, r2_file, index_prefix, sam_file, threads = threads)     
        else:
            bowtie_align_time, bowtie_align_mem = align_with_bowtie2(fastq_file, index_prefix, sam_file, threads = threads)
        bt_align_end = time.time()
        print(f'Finished Bowtie Alignment in {bt_align_end-bt_start_time}')
        max_bowtie_mem = max(bowtie_build_mem, bowtie_align_mem)
        
    ## run lofreq
    if re_run_bench:
        lo_freq_start = time.time()
        vcf_lofreq, lofreq_time, lofreq_mem, samtools_time, samtools_mem = run_lofreq_pipeline(
            sam_file=sam_file,
            ref_fasta=ref_fasta,
            output_prefix=f"{output_folder}/{base_name}_lofreq",
            threads=threads
        )
        lo_freq_finish = time.time()
        print(f'Finished bowtie2 + lofreq in {lo_freq_finish - bt_start_time}')
    else:
        vcf_lofreq = f"{output_folder}/{base_name}_lofreq.vcf"
    
    ## run ivar
    if re_run_bench:
        ivar_start = time.time()
        ivar_out, ivar_time, ivar_mem, samtools_time, samtools_mem = run_ivar_pipeline(
            sam_file=sam_file,
            ref_fasta=ref_fasta,
            output_prefix=f"{output_folder}/{base_name}_ivar",
            threads=threads, 
            min_af=min_af)
        ivar_finish = time.time()
        print(f'Finished bowtie2 + ivar in {bt_align_end-bt_start_time+ivar_finish-ivar_start}')
    else:
        ivar_out = f"{output_folder}/{base_name}_ivar_ivar.tsv"
        
    if re_run_bench:
        print('\nStarting bowtie2 parsing')
        bowtie_parse_start = time.time()
        base_counts, ref_seqs, mapped, unmapped = parse_alignment(sam_file, ref_fasta)

        write_counts_to_tsv(base_counts, ref_seqs, output_tsv=f"{output_folder}/bt_results.tsv")
        print("\nNucleotide Counts at Each Reference Position:")
        print(f"\nTotal Reads Mapped: {mapped}")
        print(f"Total Reads Unmapped: {unmapped}")
        bt_end_time = time.time()
    
        print(f'Parsing Bowtie2 Samfile Finished in {round(bt_end_time-bowtie_parse_start,2)}s')
        print(f'Both Bowtie2 Pipelines Finished in {round(bt_end_time-bt_start_time,2)}s')


    ## Run bronko
    print("\n\n Running bronko:")
    iv_start_time = time.time()
    print(k)
    # build(ref_fasta, k, out=output_folder)
    # test_final_time = time.time()
    # print(f"bronko build finished in {test_final_time - test_time}")
    # screen(fastq_file, f"{output_folder}/iv_kmer.pkl", k, output_folder, min_count)
    if r1_file and r2_file:
        bronko_cmd = [
            "cargo",
            "run",
            "--",
            "call",
            "-g", ref_fasta,
            "-1", r1_file,
            "-2", r2_file,
            "--threads", str(threads),
            "--output", output_folder,
            "-k", str(k),
            "--pileup",
            "--min-af", str(min_af),
            "--keep-kmer-info"
        ]
        bronko_time, bronko_mem = run_cmd_mem(bronko_cmd, "Running bronko")
    else:
        bronko_cmd = [
            "bronko",
            "call", 
            "-g", ref_fasta,
            "-r", fastq_file,
            "--threads", str(threads),
            "--output", output_folder,
            "-k", str(k),
            "--pileup",
            "--min-af", str(min_af),
            "--keep-kmer-info"
        ]
        bronko_time, bronko_mem = run_cmd_mem(bronko_cmd, "Running bronko")
    iv_end_time = time.time()
    print(f'bronko finished in {round(iv_end_time-iv_start_time,2)}s')

    # # OLD TEMPORARY VARIANT CALLING DIRECTLY FROM READ COUNTS        
    # iv, bt_iv, overlap = plot_depth(f'{output_folder}/bronko.tsv', f"{output_folder}/bt_results.tsv", output_folder, min_af=min_af, end_filter=True, k=k)
    # total = iv + bt_iv + overlap
        
    # print('\nComparison Using Raw Filters:')
    # print(f'Unique to IV tool: {iv}/{total}={round(iv/total,2)}')
    # print(f'Unique to BT: {bt_iv}/{total}={round(bt_iv/total,2)}')
    # print(f'Shared between both: {overlap}/{total}={round(overlap/total,2)}')
    
    # print(f'\nSensitivity: {overlap}/{overlap+bt_iv}={round(overlap/(overlap+bt_iv),2)}')
    # print(f'Precision: {overlap}/{overlap+iv}={round(overlap/(overlap+iv),2)}')
    # print(f'Accuracy: {overlap}/{overlap+iv+bt_iv}={round(overlap/(overlap+iv+bt_iv),2)}')
    
    result_comparison = compare_variant_sets(f'{output_folder}/{base_name}.vcf', vcf_lofreq, ivar_out, af_cutoff=min_af, edge_exclude_k=k)
    
    minor_rows = build_variant_table_rows(
        bronko_counts_file=f"{output_folder}/{base_name}.tsv",
        bowtie2_counts_file=f"{output_folder}/bt_results.tsv",
        bronko_calls=result_comparison['minor']['bronko'],
        lofreq_calls=result_comparison['minor']['lofreq'],
        ivar_calls=result_comparison['minor']['ivar'],
        output_file=f"{output_folder}/minor_variant_comparison.tsv"
    )

    major_rows = build_variant_table_rows(
        bronko_counts_file=f"{output_folder}/{base_name}.tsv",
        bowtie2_counts_file=f"{output_folder}/bt_results.tsv",
        bronko_calls=result_comparison['major']['bronko'],
        lofreq_calls=result_comparison['major']['lofreq'],
        ivar_calls=result_comparison['major']['ivar'],
        output_file=f"{output_folder}/major_variant_comparison.tsv"
    )
    # major_allele_dists, minor_allele_dists = mahalanobis_distance_between_major_minor(f'{output_folder}/iv_output.tsv', f"{output_folder}/bt_results.tsv")
    # print(major_allele_dists, minor_allele_dists)
    
    ## Print out the final time comparisons
    if re_run_bench:
        print('\nFinal Time Comparison:')
        print(f'Only Bowtie2: {round(bt_align_end-bt_start_time, 2)}')
        print(f'Bowtie2 (with direct sam parsing):{round(bt_end_time-bt_start_time - (lo_freq_finish-lo_freq_start),2)}s')
        print(f'Bowtie2 + lofreq: {round(lo_freq_finish - bt_start_time,2)}')
        print(f'Bowtie2 + ivar: {round(bt_align_end-bt_start_time+ivar_finish-ivar_start, 2)}')
        print(f'bronko:{round(iv_end_time-iv_start_time,2)}s')
    
    ## return time info and results
    if re_run_bench:
        time_info = {"bowtie2": round(bt_align_end-bt_start_time, 2),
                    "samtools": round(samtools_time, 2),
                    "bowtie2+lofreq": round(lo_freq_finish - bt_start_time,2), 
                    "bowtie2+ivar": round(bt_align_end-bt_start_time+ivar_finish-ivar_start, 2),
                    "bronko": round(iv_end_time-iv_start_time,2)}
        mem_info = {"bowtie2": round(max_bowtie_mem, 2),
                    "samtools": round(samtools_mem, 2),
                    "bowtie2+lofreq": round(lofreq_mem,2), 
                    "bowtie2+ivar": round(ivar_mem, 2),
                    "bronko": round(bronko_mem,2)}
    else:
        time_info = {"bowtie2": None,
                    "samtools": None, 
                    "bowtie2+lofreq": None, 
                    "bowtie2+ivar": None,
                    "bronko": round(iv_end_time-iv_start_time,2)}
        mem_info = {"bowtie2": None, 
                    "bowtie2+lofreq": None, 
                    "samtools": None,
                    "bowtie2+ivar": None,
                    "bronko": round(bronko_mem,2)}
    
    return result_comparison, time_info, mem_info, minor_rows, major_rows

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python pipeline.py <reads.fastq> <reference.fasta> <output_folder>")
        sys.exit(1)
    bench(sys.argv[1], sys.argv[2], sys.argv[3])
