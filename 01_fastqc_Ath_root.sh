#!/bin/bash
set -euo pipefail

# CONFIG
METADATA="/prj/pflaphy-ironph/Transcriptome_data/metadata_all_samples.csv"
FASTQC_DIR="/prj/pflaphy-ironph/QC/fastqc"
TARGET_DIR="$FASTQC_DIR/Ath_root"
LOG="$FASTQC_DIR/fastqc_ath_root.log"

mkdir -p "$TARGET_DIR"
: > "$LOG"

echo "Start extracting A. thaliana root FastQC: $(date)" | tee -a "$LOG"

tail -n +2 "$METADATA" | while IFS=',' read -r sample_id _ _ species _ _ _ _ tissue _ _ _ file_path_R1 file_path_R2; do

    # clean fields
    species=$(echo "$species" | sed 's/^ *//;s/ *$//;s/^"//;s/"$//')
    tissue=$(echo "$tissue" | sed 's/^ *//;s/ *$//;s/^"//;s/"$//')
    file1=$(echo "$file_path_R1" | sed 's/^ *//;s/ *$//;s/^"//;s/"$//')
    file2=$(echo "$file_path_R2" | sed 's/^ *//;s/ *$//;s/^"//;s/"$//')
    
    if [[ "$species" != "Arabidopsis_thaliana" || "$tissue" != "root" ]]; then
        continue
    fi

    for fp in "$file1" "$file2"; do
        [[ -z "$fp" ]] && continue

        fname=$(basename "$fp")
        sample_base="${fname%%.*}"

       html=$(find "$FASTQC_DIR" -name "${sample_base}_fastqc.html" | head -n 1)
       zipf=$(find "$FASTQC_DIR" -name "${sample_base}_fastqc.zip" | head -n 1)

        if [[ -f "$html" ]]; then
            cp -n "$html" "$TARGET_DIR/"
            echo "COPIED: $html" | tee -a "$LOG"
        fi

        if [[ -f "$zipf" ]]; then
            cp -n "$zipf" "$TARGET_DIR/"
            echo "COPIED: $zipf" | tee -a "$LOG"
        fi

    done

done

echo "Done: $(date)" | tee -a "$LOG"
