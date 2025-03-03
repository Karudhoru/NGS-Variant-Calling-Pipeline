
```bash
#!/bin/bash
prefetch ERR13985875 --progress
fastq-dump --split-files ERR13985875/

# Run FastQC on raw FASTQ files
fastqc ERR13985875_1.fastq ERR13985875_2.fastq

# Create directory for trimmed reads and run Trimmomatic
mkdir -p Trimmed
cd Trimmed
trimmomatic PE ../ERR13985875_1.fastq ../ERR13985875_2.fastq \
    ERR13985875_1_paired.fastq ERR13985875_1_unpaired.fastq \
    ERR13985875_2_paired.fastq ERR13985875_2_unpaired.fastq \
    SLIDINGWINDOW:4:20 MINLEN:50
# Run FastQC on trimmed files
fastqc ERR13985875_1_paired.fastq ERR13985875_2_paired.fastq ERR13985875_1_unpaired.fastq ERR13985875_2_unpaired.fastq
cd ..

# Change to trimmed directory, perform alignment with Bowtie2, then move SAM file
cd Trimmed
bowtie2 --no-unal -p 8 -x ../../Reference_Genome/Human/index \
    -1 ERR13985875_1_paired.fastq -2 ERR13985875_2_paired.fastq -S alignment.sam
mv alignment.sam ../
cd ..

# Convert SAM to BAM, sort, and index using SAMtools
samtools view -Sb -o alignment.bam alignment.sam
samtools sort -O bam -o sorted.bam alignment.bam
samtools index sorted.bam

# Add or replace read groups with GATK
gatk AddOrReplaceReadGroups \
    -I sorted.bam \
    -O sorted_rg.bam \
    --RGID 1 \
    --RGLB lib1 \
    --RGPL ILLUMINA \
    --RGPU unit1 \
    --RGSM SampleName

# Mark duplicates using GATK MarkDuplicates
gatk MarkDuplicates \
    -I sorted_rg.bam \
    -R ../Reference_Genome/Human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
    -M metrics.txt \
    -O unique_reads.bam

# Rescore the bases using GATK BaseRecalibrator 
gatk BaseRecalibrator \
    -R ../Reference_Genome/Human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
    -I unique_reads.bam \
    --known-sites ../Reference_Genome/Human/known_sites.vcf \
    -O recal_data.table

gatk ApplyBQSR \
    -R ../Reference_Genome/Human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
    -I unique_reads.bam \
    --bqsr-recal-file recal_data.table \
    -O recalibrated.bam

# Call variants using GATK HaplotypeCaller
gatk HaplotypeCaller \
    -R ../Reference_Genome/Human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
    -I recalibrated.bam \
    -O output.vcf

# Filter variants using GATK VariantFiltration
gatk VariantFiltration \
    -R ../Reference_Genome/Human/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
    -V output.vcf \
    -O filtered_output.vcf \
    --filter-expression "QD < 2.0" \
    --filter-name "LowQD" \
    --filter-expression "FS > 60.0" \
    --filter-name "HighFS"


