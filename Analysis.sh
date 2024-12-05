#!/bin/bash

# Create Conda environments
mamba create -n nanoseq -c bioconda nanopolish minimap2 samtools bedtools
mamba deactivate
mamba create -n m6a -c bioconda m6anet
mamba deactivate

# Set file paths
ENSEMBL_DIR="/home/onggf/data/ensembl"
FASTA_FILE="${ENSEMBL_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
GTF_FILE="${ENSEMBL_DIR}/Homo_sapiens.GRCh38.113.gtf"
BAM_FILE="${ENSEMBL_DIR}/sorted_sample.bam"
GENOME_FILE="${ENSEMBL_DIR}/transcripts.hg38.fa"

# Download Ensembl files (if not done yet)
if [ ! -f "$FASTA_FILE" ]; then
    wget -P $ENSEMBL_DIR https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
    gunzip ${FASTA_FILE}.gz
fi

if [ ! -f "$GTF_FILE" ]; then
    wget -P $ENSEMBL_DIR https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
    gunzip ${GTF_FILE}.gz
fi

# Process GTF file for transcripts and generate bed file
awk -F"\t" '{OFS="\t"; if($3 == "transcript"){print $1,$4-1,$5,$9,$6,$7}}' $GTF_FILE > transcripts.hg38.bed

# Extract FASTA sequences for transcripts using bedtools
bedtools getfasta -fi $FASTA_FILE -bed transcripts.hg38.bed -s > transcripts.hg38.fa

# Index the genome
minimap2 -x map-ont -d transcriptome.human.mmi transcripts.hg38.fa

# Mapping reads to transcriptome
MINIMAP_OUTPUT="sample.sam"
minimap2 -t 64 -ax splice -uf -k14 transcriptome.human.mmi /path/to/sample.fastq.gz > $MINIMAP_OUTPUT

# Convert SAM to BAM, sort, and index with samtools
samtools view -@ 32 -b -h -o sample.bam -T $FASTA_FILE $MINIMAP_OUTPUT
samtools sort -@ 32 sample.bam -o sorted_sample.bam
samtools index -@ 32 sorted_sample.bam

# Perform Nanopolish event alignment
NANO_OUTPUT_DIR="path/to/data/02.nanopolish"
mkdir -p $NANO_OUTPUT_DIR

nanopolish index -d path/to/fast5_pass/ \
--sequencing-summary path/to/sequencing_summary_run.txt \
path/to/data/sample.fastq.gz

nanopolish eventalign --reads /path/to/sample.fastq.gz \
--bam path/to/data/01.processed/sorted_sample.bam \
--genome path/to/transcripts.hg38.fa \
--scale-events \
--signal-index \
--summary path/to/data/02.m6anet/summary.txt \
--threads 50 > path/to/data/02.m6anet/eventalign.txt

# M6ANet processing and inference
M6ANET_DIR="path/to/data/03.m6anet"
mkdir -p $M6ANET_DIR

# Data preparation for M6ANet
m6anet dataprep --eventalign $NANO_OUTPUT_DIR/eventalign.txt --out_dir $M6ANET_DIR/output --n_processes 50

# M6ANet inference
m6anet inference --input_dir $M6ANET_DIR/output --out_dir /path/to/m6anet_output --n_processes 32 --num_iterations 1000

echo "Nanopore analysis pipeline completed!"
