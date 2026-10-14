#!/bin/bash

# Export paths to activate java
export PATH=/home/soft/jdk-23.0.2/bin:$PATH
export LD_LIBRARY_PATH=/home/soft/jdk-23.0.2/lib:$LD_LIBRARY_PATH

# Configuration
TRIMMOMATIC_JAR="/home/soft/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERS="/home/soft/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
OUTPUT_DIR="/HK_annotation_base/trimmed2"
THREADS=16
MEMORY="-Xmx8G"

# Trimmomatic parameters (paired outputs only)
PARAMS="ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process all R1/R2 pairs
for R1_FILE in *_R1_001.fastq.gz; do
    PREFIX="${R1_FILE%_R1_001.fastq.gz}"
    R2_FILE="${PREFIX}_R2_001.fastq.gz"

    # Define output files (paired only)
    OUTPUT_R1="${OUTPUT_DIR}/${PREFIX}_trimmed_R1_paired.fq.gz"
    OUTPUT_R2="${OUTPUT_DIR}/${PREFIX}_trimmed_R2_paired.fq.gz"

    echo "Processing $PREFIX..."
    
    # Single-line Trimmomatic command (no backslashes needed)
    java $MEMORY -jar "$TRIMMOMATIC_JAR" PE -threads $THREADS -phred33 "$R1_FILE" "$R2_FILE" "$OUTPUT_R1" /dev/null "$OUTPUT_R2" /dev/null $PARAMS

    echo "Saved: $(basename "$OUTPUT_R1") and $(basename "$OUTPUT_R2")"
done

echo "Done! Paired outputs in: $OUTPUT_DIR"
