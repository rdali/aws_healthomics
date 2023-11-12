# Running HiC on AWS

### Introduction

This repo contains code to analyze a [HiC experiment](https://www.cell.com/cms/10.1016/j.cell.2014.11.021/attachment/368a9ec6-81b4-430a-a8e9-28906e2cb77a/mmc5.mp4) to understand the difference between:  
1- bash script.   
2- Nextflow.    
3- AWS Healthomics.

It contains sample HiC sequences of 1M reads extracted from [SRR1658570](https://ddbj.nig.ac.jp/resource/sra-run/SRR1658570) from [Rao et al.](https://www.cell.com/fulltext/S0092-8674(14)01497-4). The DNA was digested with MboI.

This mini analysis consists of 8 steps, namely:

1- Download hg38 human genome.  
2- Index genome with bowtie2.  
3- Digest genome with MboI.  
4- Trim fastq reads with Trim Galore.  
5- Align reads to genome using HiCUP.  
6- Sort bam file.  
7- Produce .hic file.  
8- Create report with MultiQc.  

You can visualize the resulting `SRR1658570_1Msample.hic` file by uploading it to the [Aiden Lab platform](https://aidenlab.org/juicebox/).

**Disclaimer:** This script is made to show the difference in tools while making it a complete example without genome advance setup. In a real world environment, steps that are run at difference cadence should be separated; genome setup tasks are run once per reference and would be pulled out into their separate workflow. But for the purposes of this example, setting up the genome before running the alignment would be easier for users who have less experience setting up required files. 

## Running bash script analysis on AWS EC2

You will need to provision an EC2 instance with sufficient compute, RAM and storage. Genomics experiments tend to be demanding. Even this mini test run cannot work on t2.micros. It requires at least 8 GB of RAM and ~ 20 GB of storage.

This script was successfully tested on `c6a.xlarge` (4 vCPU, 8 GB RAM) with an attached EBS store of over 25 GB. Smaller instances resulted in errors.

1- Once to you have your instance provisioned, ssh into the instance.

**Optional:** color your prompt for better visibility:
```
export PS1='\[\033[36m\]\u\[\033[m\]@\[\033[32m\]\h:\[\033[33;1m\]\W\[\033[m\]\$ '
```

2- Clone the code repo:
```
curl -L -O https://github.com/rdali/aws_healthomics/archive/refs/heads/main.zip
unzip main.zip
cd aws_healthomics-main/
```

3- Install pipeline dependencies:  
```
bash ./bash_script/script_dependencies_EC2.sh
```

4- Source bashrc to append libraries to PATH:  
```source ~/.bashrc```

5- Run the processing script:  
```
bash ./bash_script/hic_script.sh
```

This will take ~1-2 hours to run depending on the vcpus chosen, the longest step being the genome indexing. To run the command in the background and avoid losing your session:
```
setsid nohup bash ./bash_script/hic_script.sh &
```

To monitor the process, you can use `ps x`

Once the processing is successfully completed, you should find various results files, including a .hic file: `./hic/SRR1658570_1Msample.hic`. You can upload the file at the following [site](https://aidenlab.org/juicebox/). You should see a heatmap of contact intensities across the human genome.

## Running Nextflow workflow on AWS EC2

You will need to provision an EC2 instance with sufficient compute, RAM and storage. Genomics experiments tend to be demanding. Even this mini test run cannot work on t2.micros. It requires at least 8 GB of RAM and ~ 35 GB of storage.

This script was successfully tested on `t2.2xlarge` (8 vCPU, 32 GB RAM) with an attached EBS store of over 40 GB. Smaller instances resulted in errors.

1- Once to you have your instance provisioned, ssh into the instance.

**Optional:** color your prompt for better visibility:
```
export PS1='\[\033[36m\]\u\[\033[m\]@\[\033[32m\]\h:\[\033[33;1m\]\W\[\033[m\]\$ '
```

2- Clone the code repo:
```
curl -L -O https://github.com/rdali/aws_healthomics/archive/refs/heads/main.zip
unzip main.zip
```

3- Create a tools folder:  
```
mkdir -p ~/tools
cd ~/tools
```

4- Install Java:  
```
sudo yum -y install java-11-amazon-corretto
```

5- Install conda dependency:   
```
sudo yum install -y libgcrypt-devel
```

6- Install Nextflow:  
```
curl -s https://get.nextflow.io | bash
```

7- Install Conda:
```
wget https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh
bash Anaconda3-2023.09-0-Linux-x86_64.sh
```

8- Add tools to PATH and source PATH:
```
echo "export PATH=${PATH}:~/tools:~/anaconda3/bin/" >> ~/.bashrc
source ~/.bashrc
```

9- Update conda:  
```
conda update conda
```

10- Install dependencies:
```
cd ~/aws_healthomics-main
mkdir logs
bash ./bash_script/script_dependencies_EC2.sh
```

11- Run the Nextflow pipeline:
```
nextflow run nextflow_workflow/hic_nextflow.nf -profile conda
```

**Note:** To get Nextflow reports, use the following:

```
nextflow run nextflow_workflow/hic_nextflow.nf -profile conda -with-report -with-trace -with-timeline -with-dag
```

To monitor the process, you can use `ps x`

Once the processing is successfully completed, you should find various results files, including a .hic file: `./hic/SRR1658570_1Msample.hic`. You can upload the file at the following [site](https://aidenlab.org/juicebox/). You should see a heatmap of contact intensities across the human genome.

On a `t2.2xlarge`, the run stats were:

```
Completed at: 08-Nov-2023 23:13:23
Duration    : 2h 40m 22s
CPU hours   : 2.7
Succeeded   : 8
Used GBs of storage: 26 GB
```

## Running AWS Healthomics

