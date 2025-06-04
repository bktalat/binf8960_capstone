#!/bin/bash
#SBATCH --job-name=QC
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --error=QC_%j.err
#SBATCH --out=QC_%j.out
#SBATCH --mem=32gb
#SBATCH --time=72:00:00

# The goal is to QC, Trim and align reads to a reference genome

# Loading modules from the cluster before the analysis

ml FastQC/0.11.9-Java-11
ml Trimmomatic/0.39-Java-13
ml BWA/0.7.18-GCCcore-13.3.0
ml SAMtools/1.18-GCC-12.3.0
ml BCFtools/1.18-GCC-12.3.0
ml MultiQC/1.14-foss-2022a

#Data folder
data_folder="/home/bt51523/data"

# The genome and reads were downloaded from online databases, in this case NCBI.
# The first step is to evaluate the quality of the reads with FastQC

# going in the data folder with cd
#cd ~/data/
#fastqc -t 8 --nogroup --noextract *.fastq.gz
#multiqc -o ~/data/ ~/data/

# Trimming the reads using Trimmomatic, file containing the Nextera adapters sequences is.
# Trim raw reads with trimmomatic

echo "Running Trimmomatic"

mkdir ./data/trimmed_fastq/

TRIMMOMATIC="java -jar /apps/eb/Trimmomatic/0.39-Java-13/trimmomatic-0.39.jar" 
for fwd in ./data/*_1.fastq.gz
do
  sample=$(basename $fwd _1.fastq.gz)
  echo "    Trimming sample $sample" 
  $TRIMMOMATIC PE ./data/${sample}_1.fastq.gz ./data/${sample}_2.fastq.gz  \
      ./data/trimmed_fastq/${sample}_1.paired.fastq.gz ./data/trimmed_fastq/${sample}_1.unpaired.fastq.gz \
      ./data/trimmed_fastq/${sample}_2.paired.fastq.gz ./data/trimmed_fastq/${sample}_2.unpaired.fastq.gz \
      ILLUMINACLIP:./data/NexteraPE-PE.fa:2:30:10:5:True SLIDINGWINDOW:4:20

done


# Run FASTQC on trimmed reads
echo "Running FASTQC on trimmed reads"
module load FastQC/0.11.9-Java-11

mkdir ./results/
mkdir ./results/fastqc_trimmed_results
fastqc -o ./results/fastqc_trimmed_results ./data/trimmed_fastq/*.paired.fastq.gz  # Only paired reads

# Run MultiQC to compile the trimmed FASTQC results

echo "Runing MultiQC on trimmed reads"
module load MultiQC/1.14-foss-2022a
multiqc -o results/fastqc_trimmed_results results/fastqc_trimmed_results

# Indexing the genome for alignment
bwa index $data_folder/ecoli_ref.fna

# Align the samples to the reference genome
for fwd in $data_folder/trimmed_fastq/*_PE_1.fq
do
	sample=$(basename $fwd _PE_1.fq)
	echo "Aligning $sample"
	rev=$data_folder/trimmed/${sample}_PE_2.fq
	bwa mem $data_folder/ecoli_ref.fna $fwd $rev > ~/results/$sample.sam

	#convert to BAM and sort the files
	samtools view -S -b ~/results/$sample.sam > ~/results/$sample.bam
	samtools sort -o ~/results/*.sorted.bam ~/ecoli/results/*.bam

	#variant calling
	echo "Variant calling in file $sample"
	bcftools mpileup -O b -o ~/results/$sample.bcf -f $data_folder/ecoli_ref.fna ~/results/$sample.bam
	bcftools call --ploidy -1 -m -v -o ~/results/$sample.vcf ~/results/$sample.bcf

done
