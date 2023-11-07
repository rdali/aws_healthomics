
################################################################################
# This script runs a small HiC analysis from Fastq file to .hic file. It:
# 1- Downloads hg38
# 2- Indexes genome with bowtie2
# 3- Digests genome with MboI
# 4- Trims fastq reads with Trim Galore
# 5- Aligns reads to genome using HiCUP
# 6- Sorts bam file
# 7- Produces .hic file
# 8- Creating report with MultiQc
#
# This script was tested on a Mac and requires several libraries whose
# installation was included in script_dependencies.txt.
# The script is expecting a folder called "fastq" that has 2 fastq files
# with the naming ${id}_{1,2}.fastq.gz
# To use this script:
# 1- install required tools
# 2- download fastq files in a directory named fastq
# 3- Edit the "Set Variable" section to point to your tools
################################################################################

## Set Variables:
dir=`pwd`
id=SRR1658570_1Msample
bowtie_path=~/tools/bowtie2-2.5.2-sra-linux-x86_64/bowtie2
R_path=/usr/bin/R
juicer_path=~/tools/juicer_tools_1.22.01.jar


## create folder structure:
echo "... ... ... Creating Folder Structure ... ... ..."
mkdir -p genome/index
mkdir -p genome/digest
mkdir bams

## 1- Download genome file hg38 from UCSC:
echo "... ... ... Downloading hg38 Genome File ... ... ..."
curl --output genome/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

## 2- Create bowtie2 index
echo "... ... ... Indexing hg38 Genome file with Bowtie2 ... ... ..."
bowtie2-build genome/hg38.fa.gz genome/index/hg38


## 3- Create Digest file:
echo "... ... ... Creating MboI digest of hg38 Genome ... ... ..."
hicup_digester --genome hg38 --re1 ^GATC,MboI genome/hg38.fa.gz --outdir genome/digest/
mv genome/digest/Digest_hg38_MboI_*.txt genome/digest/hg38_MboI.txt

## 4- Trim adapters and bad sequences with trim galore:
echo "... ... ... Trimming Fastq Reads ... ... ..."
trim_galore --paired --phred33 -o ${dir}/fastq/ ${dir}/fastq/${id}_1.fastq.gz ${dir}/fastq/${id}_2.fastq.gz

## 5- Align reads to hg38 with HICUP:
echo "... ... ... Aligning reads to hg38 with HiCUP ... ... ..."

echo "
Outdir: ${dir}/bams/
Threads: 5
Quiet:0
Keep:0
Zip:1
Bowtie2: ${bowtie_path}
R: ${R_path}
Index: ${dir}/genome/index/hg38
Digest: ${dir}/genome/digest/hg38_MboI.txt
Format: Sanger
Longest: 800
Shortest: 0
${dir}/fastq/${id}_1_val_1.fq.gz
${dir}/fastq/${id}_2_val_2.fq.gz" > hicup_${id}.conf

hicup -c ${dir}/hicup_${id}.conf


## 6- Sort bam:
echo "... ... ... Sorting Bam ... ... ..."
samtools sort -n bams/${id}_1_val_1.${id}_2_val_2.hicup.bam -o bams/${id}.sorted.bam -O BAM


## 7- Create .hicfile:
echo "... ... ... Creating .hicFile ... ... ..."
utils/createHiCfile.sh -b bams/${id}.sorted.bam -n ${id} -s -g hg38 -j ${juicer_path}
## you can visualize .hic file on https://aidenlab.org/juicebox/


## 8- Create report with mutiqc:
echo "... ... ... Creating report with MultiQc ... ... ..."
multiqc .


