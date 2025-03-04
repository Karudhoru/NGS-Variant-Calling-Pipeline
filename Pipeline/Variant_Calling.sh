#!/bin/bash
set -e

LOGFILE="variant_calling.log"
exec > >(tee -a "$LOGFILE") 2>&1

# Define a logging function with a timestamp.
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log "Starting variant calling pipeline..."

# Prefetch and fastq-dump steps.
log "Prefetching ERR13985875..."
prefetch ERR13985875 --progress
log "Prefetch complete."

log "Running fastq-dump to split FASTQ files..."
fastq-dump --split-files ERR13985875/
log "fastq-dump complete."

# Run FastQC on raw FASTQ files.
log "Running FastQC on raw FASTQ files..."
fastqc ERR13985875_1.fastq ERR13985875_2.fastq
log "FastQC on raw FASTQ files complete."

# Create directory for trimmed reads and run Trimmomatic.
log "Creating directory for trimmed reads and running Trimmomatic..."
mkdir  Trimmed
cd Trimmed
trimmomatic PE ../ERR13985875_1.fastq ../ERR13985875_2.fastq \
    ERR13985875_1_paired.fastq ERR13985875_1_unpaired.fastq \
    ERR13985875_2_paired.fastq ERR13985875_2_unpaired.fastq \
    SLIDINGWINDOW:4:20 MINLEN:50
log "Trimmomatic processing complete."

# Run FastQC on trimmed files.
log "Running FastQC on trimmed FASTQ files..."
fastqc ERR13985875_1_paired.fastq ERR13985875_2_paired.fastq ERR13985875_1_unpaired.fastq ERR13985875_2_unpaired.fastq
log "FastQC on trimmed FASTQ files complete."
cd ..

# Alignment with Bowtie2.
log "Build an index for the reference genome"
bowtie2-build ../Reference_Genome/Human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz ../Reference_Genome/Human/index
log "Index successfully built on the reference genome"

log "Changing to Trimmed directory for alignment..."
cd Trimmed
log "Running Bowtie2 alignment..."
bowtie2 --no-unal -p 8 -x ../../Reference_Genome/Human/index \
    -1 ERR13985875_1_paired.fastq -2 ERR13985875_2_paired.fastq -S alignment.sam
log "Bowtie2 alignment complete."
mv alignment.sam ../
cd ..

# Convert SAM to BAM, sort, and index using SAMtools.
log "Converting SAM to BAM..."
samtools view -Sb -o alignment.bam alignment.sam
log "Sorting BAM file..."
samtools sort -O bam -o sorted.bam alignment.bam
log "Indexing sorted BAM file..."
samtools index sorted.bam

# Add or replace read groups with GATK.
log "Adding or replacing read groups with GATK..."
gatk AddOrReplaceReadGroups \
    -I sorted.bam \
    -O sorted_rg.bam \
    --RGID 1 \
    --RGLB lib1 \
    --RGPL ILLUMINA \
    --RGPU unit1 \
    --RGSM SampleName
log "Read groups added/replaced."

# Mark duplicates using GATK MarkDuplicates.
log "Marking duplicates with GATK..."
gatk MarkDuplicates \
    -I sorted_rg.bam \
    -R ../Reference_Genome/Human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
    -M metrics.txt \
    -O unique_reads.bam
log "Duplicates marked."

# Rescore the bases using GATK BaseRecalibrator.
log "Performing BaseRecalibration with GATK BaseRecalibrator..."
gatk BaseRecalibrator \
    -R ../Reference_Genome/Human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
    -I unique_reads.bam \
    --known-sites ../Reference_Genome/Human/known_sites.vcf \
    -O recal_data.table
log "BaseRecalibration complete."

log "Applying BQSR with GATK ApplyBQSR..."
gatk ApplyBQSR \
    -R ../Reference_Genome/Human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
    -I unique_reads.bam \
    --bqsr-recal-file recal_data.table \
    -O recalibrated.bam
log "BQSR applied."

# Call variants using GATK HaplotypeCaller.
log "Calling variants with GATK HaplotypeCaller..."
gatk HaplotypeCaller \
    -R ../Reference_Genome/Human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
    -I recalibrated.bam \
    -O output.vcf
log "Variant calling complete."

# Filter variants using GATK VariantFiltration.
log "Filtering variants with GATK VariantFiltration..."
gatk VariantFiltration \
    -R ../Reference_Genome/Human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
    -V output.vcf \
    -O filtered_output.vcf \
    --filter-expression "QD < 2.0" \
    --filter-name "LowQD" \
    --filter-expression "FS > 60.0" \
    --filter-name "HighFS"
log "Variant filtration complete."

# Annotate variants using GATK VariantAnnotator.
log "Annotating variants with GATK VariantAnnotator..."
gatk VariantAnnotator \
    -R ../Reference_Genome/Human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
    -V filtered_output.vcf \
    -I recalibrated.bam \
    -O output.annotated.vcf \
    -A Coverage \
    -A QualByDepth \
    -A MappingQualityRankSumTest
log "Variant annotation complete."


log "Variant calling pipeline finished successfully."
