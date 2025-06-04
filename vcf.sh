set -e

# Location to download genome from
genome_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz"

# Load required software
module load BWA/0.7.18-GCCcore-13.3.0
module load SAMtools/1.18-GCC-12.3.0
module load BCFtools/1.18-GCC-12.3.0

# Make required folders
mkdir data/genomes mkdir results/sam results/bam results/vcf results/bcf

# Download E. coli genome
wget -O data/genomes/ecoli_rel606.fna.gz  $genome_url
gunzip data/genomes/ecoli_rel606.fna.gz


# Index the genome
bwa index data/genomes/ecoli_rel606.fna


# Loop over reads and do alignment and variant calling
for fwd in data/trimmed_fastq/*_1.paired.fastq.gz
do
   sample=$(basename $fwd _1.paired.fastq.gz )

   # Run alignment
   echo "Aligning $sample"
   rev=data/trimmed_fastq_small/${sample}_2.paired.fastq.gz
   bwa mem data/genomes/ecoli_rel606.fna $fwd $rev > results/sam/$sample.sam 

   # Convert to BAM and sort 
   samtools view -S -b results/sam/$sample.sam > results/bam/$sample.bam
   samtools sort -o results/bam/$sample.sorted.bam results/bam/$sample.bam

   # Variant calling
   echo "Calling variants in $sample"
   bcftools mpileup -O b -o results/bcf/$sample.bcf -f data/genomes/ecoli_rel606.fna results/bam/$sample.sorted.bam
   bcftools call --ploidy 1 -m  -v -o results/vcf/$sample.vcf  results/bcf/$sample.bcf

done
