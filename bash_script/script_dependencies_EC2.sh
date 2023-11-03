## This is a user script for an EC2 Linux instance to set up library dependencies to run the hic bash script:
## This is writen by LogicWorks & tested Nov 2023
## Versions might change, but the process is unlikly to change
## tested on Amazon Linux 2023 AMI 2023.2.2 on EC2 t2.xlarge



## Amazon Linux EC2:

## create a tools folder:
mkdir ~/tools
cd ~/tools
## add tools to PATH: export PATH=${PATH}:~/tools

## install pip: tested pip 21.3.1 (python 3.9)
sudo yum -y install python-pip

## install java: OpenJDK Runtime Environment Corretto-11.0.21.9.1 (build 11.0.21+9-LTS)
sudo yum -y install java-11-amazon-corretto

## install perl: tested v5.32.1
sudo yum -y install perl

## install cutadapt: tested 4.5
pip install cutadapt

## install fastqc: v0.12.1
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
chmod 755 ./FastQC/fastqc
## add to path: export PATH=${PATH}:~/tools/FastQC

## install trim galore: tested 0.6.10
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
## add to path: export PATH=${PATH}:~/tools/TrimGalore-0.6.10

## intall HICUP: tested v0.9.2
wget -O hicup.tar.gz https://github.com/StevenWingett/HiCUP/archive/refs/tags/v0.9.2.tar.gz
tar -xvzf hicup.tar.gz
## add to path: export PATH=${PATH}:~/tools/HiCUP-0.9.2

## install samtools: tested v1.18
sudo yum -y group install "Development Tools"
sudo yum -y install ncurses-devel bzip2-devel xz-devel
cd /tmp
wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
tar xvjf samtools-1.18.tar.bz2
cd samtools-1.18
./configure --prefix=/usr/local
make
sudo make install
cd ~/tools

## install bowtie2: tested 2.5.2
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.2/bowtie2-2.5.2-sra-linux-x86_64.zip
unzip bowtie2-2.5.2-sra-linux-x86_64.zip
## add to path: export PATH=${PATH}:~/tools/bowtie2-2.5.2-sra-linux-x86_64

## install MultiQc for reports: tested 1.17
pip install multiqc

## download juicer_tools: tested 1.22.01
wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar

## intall R:  tested 4.3.1
sudo yum -y install R
# gets intalled in /usr/bin/R


## Add tools to PATH:
echo "export PATH=${PATH}:~/tools:~/tools/FastQC:~/tools/TrimGalore-0.6.10:~/tools/HiCUP-0.9.2:~/tools/bowtie2-2.5.2-sra-linux-x86_64" >> ~/.bashrc

