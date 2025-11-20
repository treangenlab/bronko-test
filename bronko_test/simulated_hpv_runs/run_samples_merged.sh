#!/bin/bash

# Parameters
BASE_DIR="/home/Users/rdd4/GENOMICON-Seq-OLD"
READ_DIR="$BASE_DIR/bronko_final_data/HPV_baseline"
REF="$BASE_DIR/input_data_ampliseq/HPV16.fa"
OUT_BASE="/home/Users/rdd4/bronko_test/simulated_hpv_runs/HPV_baseline"
# ADAPTER_FASTA="/home/Users/rdd4/GENOMICON-Seq/input_data_ampliseq/mapping_reference/all_primers.fasta"
THREADS=10

MIN_AF=0.001
K=21

# Loop over sample numbers 1 to 10
for i in {1..10}; do
    SAMPLE="sample_$i"
    SAMPLE_DIR="$READ_DIR/$SAMPLE/generated_reads"
    MERGED_DIR="$READ_DIR/$SAMPLE/merged_reads"
    mkdir -p "$MERGED_DIR"

    echo "Merging reads for $SAMPLE..."

    # Merge each replicate's R1 and R2 across PCR runs and run fastp
    for rep in 2; do
        for read in 1 2; do
            out_file="$MERGED_DIR/rep${rep}_R${read}.fastq.gz"
            cat "$SAMPLE_DIR/tech_replicate_${rep}_PCR_"*_R${read}.fastq.gz > "$out_file"
        done
    done

    echo "Running iv.py bench for $SAMPLE..."
    python /home/Users/rdd4/bronko/iv.py bench \
        --r1 "$MERGED_DIR/rep2_R1.fastq.gz" \
        --r2 "$MERGED_DIR/rep2_R2.fastq.gz" \
        -fa "$REF" \
        -o "$OUT_BASE/$SAMPLE" \
        --min-af "$MIN_AF" \
        -k "$K" \
        --threads "$THREADS" \
        --no-rerun

    # echo "Running iv.py bench for $SAMPLE..."
    # python /home/Users/rdd4/bronko/iv.py bench \
    #     --r1 "$MERGED_DIR/rep1_R1.fastq.gz" "$MERGED_DIR/rep2_R1.fastq.gz" \
    #     --r2 "$MERGED_DIR/rep1_R2.fastq.gz" "$MERGED_DIR/rep2_R2.fastq.gz" \
    #     -fa "$REF" \
    #     -o "$OUT_BASE/$SAMPLE" \
    #     --min-af "$MIN_AF" \
    #     -k "$K" \
    #     --threads "$THREADS"
done