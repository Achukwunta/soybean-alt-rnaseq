#!/bin/bash

: '
This workflow is used to generate read counts for differential gene expression.
Control samples: healthy soybean with 4 replications
Treatment samples: soybean infected with alternaria with 4 replications

Step 1: Quality check for raw reads (FastQC, MultiQC)
Step 2: Trim off bad quality reads and adapter sequences (Trimmomatic)
Step 3: Map reads to the reference genome (HISAT2)
Step 4: Sort the bam files (Samtools)
Step 5: Generate quantification metrics (featureCounts)

'

# Set duration for running the script
SECONDS=0

# Change the currect directory to a suitable working directory
cd /mnt/e/glycine_max

# create directories
mkdir -p raw_reads/fastq_report
mkdir -p raw_reads/multiqc_report

# Step 1: Quality check on raw reads

# Run fastqc for quality check on raw reads
fastqc raw_reads/*.fastq.gz -o raw_reads/fastq_report
# Aggregate fastqc raw files using multiqc to generate a single report for all the reads
multiqc -d -s raw_reads/fastq_report/*_fastqc.zip -o raw_reads/multiqc_report

mkdir -p trimmed_reads/fastq_report
mkdir -p trimmed_reads/multiqc_report
mkdir -p trimmed_reads/paired
mkdir -p trimmed_reads/unpaired

raw_reads_dir="raw_reads/"
paired_trimmed_reads="trimmed_reads/paired/"
unpaired_trimmed_reads="trimmed_reads/unpaired/"
threads=16

# Step 2: Trim off bad quality reads and adapter sequences

# Trimmomatic: Trim poor quality reads and adapter sequence
for R1_read in ${raw_reads_dir}*R1.fastq.gz; do
    # Get corresponding R2 file#
    base_trim=$(basename "${R1_read%_R1.fastq.gz}")
    R2_read="${raw_reads_dir}${base_trim}_R2.fastq.gz"

    # Construct output file names for paired and unpaired reads
    paired_R1="${paired_trimmed_reads}${base_trim}.trimmed.paired.R1.fastq.gz"
    unpaired_R1="${unpaired_trimmed_reads}${base_trim}.trimmed.unpaired.R1.fastq.gz"
    paired_R2="${paired_trimmed_reads}${base_trim}.trimmed.paired.R2.fastq.gz"
    unpaired_R2="${unpaired_trimmed_reads}${base_trim}.trimmed.unpaired.R2.fastq.gz"

    # Run Trimmomatic with adjusted parameters
    trimmomatic PE -threads $threads \
        "$R1_read" "$R2_read" \
        "$paired_R1" "$unpaired_R1" \
        "$paired_R2" "$unpaired_R2" \
        ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:2:True \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:5:18 MINLEN:30

done

echo "Trimmomatic finished running!"

# Run quality on trimmed-reads using FastQC to ensure you're working with good quality reads
fastqc trimmed_reads/paired/*.fastq.gz -o trimmed_reads/fastq_report
fastqc trimmed_reads/unpaired/*.fastq.gz -o trimmed_reads/fastq_report


# Aggregate trimmed reads using multiqc to generate a single report for read quality
multiqc -d -s trimmed_reads/fastq_report/*_fastqc.zip -o trimmed_reads/multiqc_report


# Step 3: Align reads to reference genome


# Prepare reads for alignment to reference genome
# Unzip trimmed reads (optional)
gzip -d trimmed_reads/paired/*.fastq.gz
gzip -d trimmed_reads/unpaired/*.fastq.gz


#mkdir -p ref_genome_sb/genome_index

# Index reference genome
#hisat2-build -p 16 \
#    ref_genome_sb/Gmax_880_v6.0.fa \
#    ref_genome_sb/genome_index

#echo "HISAT2 finished indexing!"


# Create directories and sub-directories to contain mapped reads
mkdir -p mapped_reads/paired
mkdir -p mapped_reads/unpaired

mapped_paired="mapped_reads/paired/"
mapped_unpaired="mapped_reads/unpaired/"
genome_index="ref_genome_sb/genome_index"

# Map reads to reference genome

# Part A: map paired-end reads
echo "HISAT2 Alignment: Mapping paired-end reads to the reference genome..."
for trimmed_paired_read in ${paired_trimmed_reads}*.fastq; do

    # Define base name for R1 and R2 trimmed paired reads
    base_paired="${trimmed_paired_read%.trimmed.paired.R1.fastq}"
    base_paired="${base_paired%.trimmed.paired.R2.fastq}"
    base_paired=$(basename "$base_paired")

    # Construct input file names for paired reads
    R1_trimmed_paired="${paired_trimmed_reads}${base_paired}.trimmed.paired.R1.fastq"
    R2_trimmed_paired="${paired_trimmed_reads}${base_paired}.trimmed.paired.R2.fastq"

# Run HISAT2 for paired-end reads and save alignment summary to a log file
    hisat2 -p $threads \
        -x $genome_index \
        -1 $R1_trimmed_paired \
        -2 $R2_trimmed_paired \
        -S ${mapped_paired}${base_paired}.sam \
        2> ${mapped_paired}${base_paired}.hisat2_paired.log

    # Convert SAM to BAM and sort
    samtools view -bS ${mapped_paired}${base_paired}.sam | samtools sort -o ${mapped_paired}${base_paired}.sorted.bam

    # Remove intermediate SAM file
    rm ${mapped_paired}${base_paired}.sam

done
echo "HISAT2 finished mapping paired-end reads!"
echo "Samtools: finished sorting paired-end reads!"



# Part B: map unpaired (single-end) reads
echo "HISAT2 Alignment: Mapping unpaired reads (single-end reads) to the reference genome..."
for trimmed_unpaired_read_R1 in "${unpaired_trimmed_reads}"*.trimmed.unpaired.R1.fastq; do

    # Get the corresponding R2 read
    trimmed_unpaired_read_R2="${trimmed_unpaired_read_R1%.R1.fastq}.R2.fastq"

    # Define base name for trimmed unpaired reads
    base_unpaired=$(basename "${trimmed_unpaired_read_R1%.trimmed.unpaired.R1.fastq}")

# Run HISAT2 for unpaired reads and save alignment summary to a log file
    hisat2 -p $threads \
        -x $genome_index \
        -U "$trimmed_unpaired_read_R1","$trimmed_unpaired_read_R2" \
        -S "${mapped_unpaired}${base_unpaired}.sam" \
        2> ${mapped_unpaired}${base_unpaired}.hisat2_unpaired.log


    # Convert SAM to BAM and sort
    samtools view -bS ${mapped_unpaired}${base_unpaired}.sam | samtools sort -o ${mapped_unpaired}${base_unpaired}.sorted.bam
    # Remove intermediate SAM file
    rm ${mapped_unpaired}${base_unpaired}.sam

done
echo "HISAT2 finished mapping (unpaired)!"
echo "Samtools: finished sorting single-end reads!"

# Visualize mapped reads (paired-end and single-end) using Multiqc
destination_dir="mapped_reads/"
cp "$mapped_paired"/*.hisat2_paired.log "$destination_dir"
cp "$mapped_unpaired"/*.hisat2_unpaired.log "$destination_dir"

# Visualize the alignment reports for both paired-end and single-end reads
multiqc -d -s mapped_reads/*.log -o mapped_reads/multiqc_report


# Step 4: Generate read counts

mkdir -p mapped_reads/counts

# Paired sorted bam files
featureCounts -p -t exon -g gene_id -a ref_genome_sb/gmax_gene_v6.1.gtf\
 -o mapped_reads/counts/counts_paired.txt -T 8 mapped_reads/paired/*sorted.bam

# Unpaired sorted bam files
featureCounts -t exon -g gene_id -a ref_genome_sb/gmax_gene_v6.1.gtf\
 -o mapped_reads/counts/counts_unpaired.txt -T 8 mapped_reads/unpaired/*sorted.bam

echo "FeatureCounts completed successfully!!"
echo "reads ready for downstream analysis!!!"

