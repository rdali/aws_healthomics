/*
========================================================================================
   Info
========================================================================================

This is an example HiC nextflow script. It is truncated and simplified to demonstrate
basic concepts, specifically, the difference betwen a bash script and a nextflow
workflow.

This workflow runs the following steps:
# 1- Downloads hg38
# 2- Indexes genome with bowtie2
# 3- Digests genome with MboI
# 4- Trims fastq reads with Trim Galore
# 5- Aligns reads to genome using HiCUP
# 6- Sorts bam file
# 7- Produces .hic file
# 8- Creating report with MultiQc

*/



/*
========================================================================================
   Params
========================================================================================
*/


// Script parameters

// env params
params.dir = "."

// tool paths
params.bowtie_path = "~/tools/bowtie2-2.5.2-sra-linux-x86_64/bowtie2"
params.r_path = "/usr/bin/R"
params.juicer_path = "~/tools/juicer_tools_1.22.01.jar"

// genome params
params.genomeUrl = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
params.genome = "hg38"
params.genomeFa = "hg38.fa.gz"
params.genomeFaPath = "${params.dir}/genome/${params.genomeFa}"
params.index = "${params.dir}/genome/index/${params.genome}"
params.enzyme = "MboI"
params.enzymeMotif = "^GATC"
params.digestfile = "${params.genome}_${params.enzyme}.txt"
params.digestPath = "${params.dir}/genome/digest/${params.digestfile}"


// sample params
params.id = "SRR1658570_1Msample"
params.readsFa = "${params.dir}/fastq/${params.id}_{1,2}.fastq.gz"

// process params
// TRIM_READS
params.trimR1 = "${params.id}_1_val_1.fq.gz"
params.trimR2 = "${params.id}_2_val_2.fq.gz"

// BAM_SORT
params.sortedBam = "${params.id}.sorted.bam"

// CREATE_HICFILE
params.create_hic_script = "${params.dir}/scripts/createHiCfile.sh"
params.hicfile = "${params.id}.hic"

// MULTIQC
params.logPath = "${params.dir}/logs"

/*
========================================================================================
   Processes
========================================================================================
*/

process GET_GENOME {
   // process downloads genome

   publishDir "${params.dir}/genome", mode:'copy'

   output:
   path "${params.genomeFa}", emit: genomeFa 

   script:
   """
   curl --output ${params.genomeFa} ${params.genomeUrl}
   """

}

process INDEX_GENOME {
   // process indexes genome with bowtie2

   publishDir "${params.dir}/genome/index", mode:'copy'

   input:
   file genome

   output:
   path "${params.genome}.1.bt2", emit: genomeIndex
   path "${params.genome}.rev.1.bt2", emit: genomeIndex_outrev

   script:
   """
   bowtie2-build ${genome} ${params.genome}
   """

}

process DIGEST_GENOME {
   // process digests genome with provided enzyme

   publishDir "${params.dir}/genome/digest", mode:'copy'

   input:
   file genomeFa

   output:
   path "${params.digestfile}", emit: genomeDigest 

   script:
   """
   hicup_digester --genome ${params.genome} --re1 ${params.enzymeMotif},${params.enzyme} ${genomeFa}
   mv Digest_${params.genome}_${params.enzyme}_*.txt ${params.digestfile}
   """

}

process TRIM_READS {
   // process trims the fastq reads with trim galore

	publishDir "${params.dir}/fastq", mode:'copy'
   publishDir "${params.logPath}", pattern: "*_trimming_report.txt", mode:'copy'

  input:
	tuple val( sample_id ), path( reads )

  output:
    path "${params.trimR1}", emit: fq1
    path "${params.trimR2}", emit: fq2
    path "*"

  script:
	"""
	trim_galore --paired --phred33 ${reads}
	"""
}


process HICUP_ALIGN {
   // process aligns reads to genome

   publishDir "${params.dir}/bams", mode:'copy'
   publishDir "${params.logPath}", pattern: "hicup_*summary*.txt"


  input:
   path fq1
   path fq2
   path genomeDigest
   path genomeIndex

  output:
   path "${params.id}_1_val_1.${params.id}_2_val_2.hicup.bam", emit: bam
   path "*"

  script:
   """
   echo "
   Outdir: .
   Threads: 5
   Quiet:0
   Keep:0
   Zip:1
   Bowtie2: ${params.bowtie_path}
   R: ${params.r_path}
   Index: ${params.index}
   Digest: ${params.digestPath}
   Format: Sanger
   Longest: 800
   Shortest: 0
   ${fq1}
   ${fq2}" > hicup_${params.id}.conf

   hicup -c ./hicup_${params.id}.conf
   """
}

process BAM_SORT {
   // process sorts bam

   publishDir "${params.dir}/bams", mode:'copy'


  input:
   path bam

  output:
    path "${params.sortedBam}", emit: bamSorted

  script:
   """
   samtools sort -n ${bam} -o ${params.sortedBam} -O BAM
   """
}


process CREATE_HICFILE {
   // process creates .hic file

   publishDir "${params.dir}/hic", mode:'copy'

  input:
   path bam_sorted
   path create_hic_script

  output:
    path "${params.hicfile}", emit: hicfile

  script:
   """
   ./${create_hic_script} -b ${bam_sorted} -n ${params.id} -s -g ${params.genome} -j ${params.juicer_path}
   """
}

process MULTIQC {
   // rums multiqc report
   publishDir "${params.dir}/multiqc", mode:'copy'

   input:
   path '*'

   output:
   path 'multiqc_report.html'

   script:
   """
   multiqc ${params.logPath}
   """
}

/*
========================================================================================
   Create Channels
========================================================================================
*/

read_pairs_ch = Channel.fromFilePairs(params.readsFa, checkIfExists: true)
hic_file_ch = Channel.fromPath(params.create_hic_script, checkIfExists: true)

/*
========================================================================================
   Main Workflow
========================================================================================
*/

log.info """\
   HiC  P I P E L I N E
   ===================================
   Sample: ${params.id}
   Genome: ${params.genome}
   Enzyme: ${params.enzyme}
   Reads : ${params.readsFa}
   """
   .stripIndent(true)



workflow {
   GET_GENOME()
   INDEX_GENOME(GET_GENOME.out.genomeFa)
   DIGEST_GENOME(GET_GENOME.out.genomeFa)
	TRIM_READS(read_pairs_ch)
   HICUP_ALIGN(TRIM_READS.out.fq1, TRIM_READS.out.fq2, DIGEST_GENOME.out.genomeDigest, INDEX_GENOME.out.genomeIndex)
   BAM_SORT(HICUP_ALIGN.out.bam)
   CREATE_HICFILE(BAM_SORT.out.bamSorted, hic_file_ch)
   MULTIQC(HICUP_ALIGN.out.bam)
}




workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
  Script Complete
========================================================================================
*/
