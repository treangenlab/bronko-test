#!/bin/bash
#
# Merge PCR replicates and run `iv.py bench` for every simulated dataset.
#
# By default this loops over every mr*/ dataset in $DATA_ROOT and every sample
# inside it. Optionally restrict to specific datasets by passing them as args:
#
#   ./run_samples_merged.sh                          # all mr* datasets
#   ./run_samples_merged.sh mr1e-4_c500_cy15         # one dataset
#   ./run_samples_merged.sh mr1e-5_c*                # a subset (shell glob)
#   SAMPLES="1 2 3" ./run_samples_merged.sh          # only samples 1-3
#   SKIP_EXISTING=0 ./run_samples_merged.sh          # redo finished samples
#   IVAR_ALL=1 ./run_samples_merged.sh               # count iVar's PASS=FALSE rows too
#
set -uo pipefail

# Parameters
BASE_DIR="/dodo/rdd4/GENOMICON-Seq"
DATA_ROOT="$BASE_DIR/bronko_simulated_data"
REF="$BASE_DIR/input_data_ampliseq/HPV16.fa"
OUT_ROOT="/home/Users/rdd4/bronko-test/bronko_test/simulated_hpv_runs_expanded"
IV_PY="/home/Users/rdd4/bronko-test/bronko_test/iv.py"
# ADAPTER_FASTA="/home/Users/rdd4/GENOMICON-Seq/input_data_ampliseq/mapping_reference/all_primers.fasta"
THREADS=20

MIN_AF=0.001
K=21

SAMPLES="${SAMPLES:-1 2 3 4 5 6 7 8 9 10}"
REPS="${REPS:-1}"
# bench compares against iVar's PASS=TRUE calls only; set IVAR_ALL=1 to include
# the rows iVar itself rejected.
IVAR_ALL="${IVAR_ALL:-0}"
# Skip a sample whose overview.tsv already exists (so the run is resumable).
SKIP_EXISTING="${SKIP_EXISTING:-0}"

FAIL_LOG="$OUT_ROOT/failed_runs.txt"

IVAR_ARGS=()
if [ "$IVAR_ALL" = "1" ]; then
    IVAR_ARGS+=(--ivar-all)
fi

# Datasets: use the ones given on the command line, otherwise every mr* folder.
if [ "$#" -gt 0 ]; then
    DATASETS=("$@")
    # tolerate trailing slashes / full paths from tab-completion
    DATASETS=("${DATASETS[@]%/}")
    DATASETS=("${DATASETS[@]##*/}")
else
    mapfile -t DATASETS < <(cd "$DATA_ROOT" && ls -d mr*/ 2>/dev/null | sed 's#/$##' | sort)
fi

if [ "${#DATASETS[@]}" -eq 0 ]; then
    echo "No datasets found in $DATA_ROOT" >&2
    exit 1
fi

echo "Datasets (${#DATASETS[@]}): ${DATASETS[*]}"
echo "Samples: $SAMPLES"
echo

n_ok=0
n_fail=0
n_skip=0

for DATASET in "${DATASETS[@]}"; do
    READ_DIR="$DATA_ROOT/$DATASET"
    OUT_BASE="$OUT_ROOT/$DATASET"

    if [ ! -d "$READ_DIR" ]; then
        echo "!! Missing dataset dir: $READ_DIR -- skipping"
        continue
    fi

    echo "================================================================"
    echo "Dataset: $DATASET"
    echo "================================================================"

    for i in $SAMPLES; do
        SAMPLE="sample_$i"
        SAMPLE_DIR="$READ_DIR/$SAMPLE/generated_reads"
        MERGED_DIR="$READ_DIR/$SAMPLE/merged_reads"

        if [ ! -d "$SAMPLE_DIR" ]; then
            echo "!! Missing reads: $SAMPLE_DIR -- skipping"
            continue
        fi

        # bench() writes the overview header before it runs the sample, so a crashed
        # run can leave a header-only file behind -- require a data row, not just bytes.
        n_rows=$(wc -l < "$OUT_BASE/$SAMPLE/overview.tsv" 2>/dev/null || echo 0)
        if [ "$SKIP_EXISTING" = "1" ] && [ "$n_rows" -ge 2 ]; then
            echo "== $DATASET/$SAMPLE already done -- skipping"
            n_skip=$((n_skip + 1))
            continue
        fi

        mkdir -p "$MERGED_DIR"

        echo "Merging reads for $DATASET/$SAMPLE..."

        # Merge each replicate's R1 and R2 across PCR runs
        merge_ok=1
        for rep in $REPS; do
            for read in 1 2; do
                out_file="$MERGED_DIR/rep${rep}_R${read}.fastq.gz"
                cat "$SAMPLE_DIR/tech_replicate_${rep}_PCR_"*_R${read}.fastq.gz > "$out_file" || merge_ok=0
            done
        done

        if [ "$merge_ok" -ne 1 ]; then
            echo "!! Merge failed for $DATASET/$SAMPLE"
            echo "$DATASET/$SAMPLE merge" >> "$FAIL_LOG"
            n_fail=$((n_fail + 1))
            continue
        fi

        echo "Running iv.py bench for $DATASET/$SAMPLE..."
        python "$IV_PY" bench \
            --r1 "$MERGED_DIR/rep1_R1.fastq.gz" \
            --r2 "$MERGED_DIR/rep1_R2.fastq.gz" \
            -fa "$REF" \
            -o "$OUT_BASE/$SAMPLE" \
            --min-af "$MIN_AF" \
            -k "$K" \
            --threads "$THREADS" \
            --no-rerun \
            "${IVAR_ARGS[@]}"

        if [ $? -eq 0 ]; then
            n_ok=$((n_ok + 1))
        else
            echo "!! bench failed for $DATASET/$SAMPLE"
            echo "$DATASET/$SAMPLE bench" >> "$FAIL_LOG"
            n_fail=$((n_fail + 1))
        fi

        # echo "Running iv.py bench for $SAMPLE (both replicates)..."
        # python "$IV_PY" bench \
        #     --r1 "$MERGED_DIR/rep1_R1.fastq.gz" "$MERGED_DIR/rep2_R1.fastq.gz" \
        #     --r2 "$MERGED_DIR/rep1_R2.fastq.gz" "$MERGED_DIR/rep2_R2.fastq.gz" \
        #     -fa "$REF" \
        #     -o "$OUT_BASE/$SAMPLE" \
        #     --min-af "$MIN_AF" \
        #     -k "$K" \
        #     --threads "$THREADS"
    done
done

echo
echo "Done. ok=$n_ok failed=$n_fail skipped=$n_skip"
if [ "$n_fail" -gt 0 ]; then
    echo "Failures logged in $FAIL_LOG"
fi
