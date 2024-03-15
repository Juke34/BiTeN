#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// Import
import static groovy.io.FileType.FILES
import java.nio.file.*

//*************************************************
// STEP 0 - parameters
//*************************************************

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 * See https://www.nextflow.io/docs/latest/config.html#configuration
 */

// Input/output params
params.help = false
params.reads_folder = "/path/to/reads/foder/"
params.genome = "/path/to/genome.fa"
params.outdir = "results"
params.pattern_reads = "*.fastq.gz" // Extension used to detect reads in folder


// Read feature params
params.single_end = true // Boolean to see if we have a single end or paired end data set
params.stranded = false // Boolean to see if we have a single or stranded data set

// Extra parameter provided by the user to the tool
params.bowtie2_options = ''


//*************************************************
// STEP 1 - LOG INFO
//*************************************************

// ------------ First an header printed in all cases -----------------
log.info """
IRD
.-./`) .-------.     ______
\\ .-.')|  _ _   \\   |    _ `''.
/ `-' \\| ( ' )  |   | _ | ) _  \\
 `-'`\"`|(_ o _) /   |( ''_'  ) |
 .---. | (_,_).' __ | . (_) `. |
 |   | |  |\\ \\  |  ||(_    ._) '
 |   | |  | \\ `'   /|  (_.\\.' /
 |   | |  |  \\    / |       .'
 '---' ''-'   `'-'  '-----'`


BiTeN - Bioinformatics Template in Nextflow
BiTeN is a template for developing pipeline in Nextflow
============================================================
"""

// ------------ A help printed only when --help is called -----------------

if (params.help) { exit 0, helpMSG() }

// Help Message
def helpMSG() {
    log.info """
    ********* HELP *********

        Usage example:
    nextflow run -profile docker main.nf --genome test/hpv16.fa --reads test

    --help                      prints the help section

        Input Reads:
    --reads                     path to the directory containing the reads 
    --pattern_reads             pattern to match the read files. In the case of single end data it would looks like: "*.fastq.gz"
                                                                 In the case of paired end data it would looks like: "*_{R1,R2}_001.fastq.gz" or "*_{1,2}.fastq.gz"        
    --single_end                Boolean to inform if we have a single end or paired end data. (default: ${params.single_end})
    --stranded                  Boolean to inform if we have a single or stranded data. (default: ${params.stranded})

        Input Genome:
    --genome                    path to the genome file in fasta format

        Alignment
    --bowtie2_options           Parameter to tune the bowtie2 aligner behaviour.  (default: ${params.bowtie2_options})
    """
}

// ------------ When --help is called we never go further. If no help asked, let's report to the user the parameters taken into account by the pipeline -----------------

log.info """
General Parameters
    genome                     : ${params.genome}
    reads                      : ${params.reads}
    reads pattern              : ${params.pattern_reads}
    single_end                 : ${params.single_end}
    outdir                     : ${params.outdir}
  
Alignment Parameters
 bowtie2 parameters
     bowtie2_options            : ${params.bowtie2_options}
 
 """

//*************************************************
// STEP 2A - Include needed modules
//*************************************************

include { bowtie2_index; bowtie2 } from "$baseDir/modules/bowtie2.nf"
// When using the same process several times like here with  fastqc you must provide a specific name
// by call using this structure "fastqc as fastqc_raw" where the process fastqc will be available here with the name fastqc_raw
include { fastqc as fastqc_raw; fastqc as fastqc_ali } from "$baseDir/modules/fastqc.nf"
include { samtools_sam2bam; samtools_sort  } from "$baseDir/modules/samtools.nf"

//*************************************************
// STEP 2B - Include needed subworkflows if outside of this file. See Sub-workflow paragraph
//*************************************************

//*************************************************
// STEP 3 - Deal with parameters
//*************************************************

// check profile
if (
    workflow.profile.contains('singularity') ||
    workflow.profile.contains('docker')
  ) { "executer selected" }
else { exit 1, "No executer selected: -profile docker/singularity"}

// check input (file or folder?)
def list_files = []
def pattern_reads =  "${params.pattern_reads}"
File input_reads = new File(params.reads)
if(input_reads.exists()){
    if ( input_reads.isDirectory()) {
       
        input_reads.eachFileRecurse(FILES){
            if (it.name =~ ~/\\*.fastq.gz/){
                list_files.add(it)
            }
        }
        samples_number = list_files.size()
        log.info "The ${params.reads} input folder contains ${samples_number} file(s) with pattern ${params.pattern_reads}, let's analyze that..."
        pattern_reads="${input_reads}/${params.pattern_reads}"
    }
    else {
        exit 1, "The input ${params.reads} is a file! A folder is expected\n"
    }
} else {
    exit 1, "The input ${params.reads} does not exists!\n"
}

//*************************************************
// Main Workflow - 
//*************************************************
// It can connect several sub workflows
// Here we have only one called ALIGN. If we do not want any subworkflow at all you will have to remove the "ALIGN(reads,genome)" line
// and then move all the code from ALIGN here excepted:
//workflow ALIGN {
//
//    take:
//        reads
//        genome
//
//    main:
//}
//*************************************************


workflow {

    main:
        Channel.fromFilePairs("${pattern_reads}", size: params.single_end ? 1 : 2, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find reads matching ${params.reads}!\n" }
            .set {reads}
        Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" }
            .set {genome}
        ALIGN(reads,genome)
}

//*************************************************
// Sub-Workflow
//*************************************************
// Sub-Workflow align 
// For clarity you may decide to move this part into a folder name subworflows in a file called e.g. align.nf
// To make it accessible from here you will have to import the subworklow as follow:
// include { ALIGN } from "${baseDir}/subworkflows/ALIGN.nf"
// A subworkflow behaves like a process, in the case your main workflow needs to get access to a result 
// emited by the sub-subworklow, you must use the emit: statement at the end of the sub-subworklow.
//*************************************************

workflow ALIGN {

    take:
        reads
        genome

    main:

        // ------------------- QC -----------------
        fastqc_raw(reads)
        
        // ------------------- BOWTIE2 -----------------
        bowtie2_index(genome) // index
        bowtie2(reads, bowtie2_index.out.collect(), genome) // align

        // ------------------- SAMTOOLS -----------------
        samtools_sam2bam(bowtie2.out.tuple_sample_sam)
        // sort
        samtools_sort(samtools_sam2bam.out.tuple_sample_bam)
        
}


//*************************************************
// Information to report at the end of the pipeline
//*************************************************

workflow.onComplete {
    log.info ( workflow.success ? "\nBiTeN pipeline complete!\n" : "Oops .. something went wrong\n" )

    log.info """
    BiTeN Pipeline execution summary
    --------------------------------------
    Completed at : ${workflow.complete}
    UUID         : ${workflow.sessionId}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Exit Status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """

    // Move pipeline execution information files into the result folder. 
    // It is safe to keep a copy in the result to keep track to what has been done to generate the results.
    int num = 0;
    String save = "${params.pipeline_report}"
    File file = new File("${params.outdir}", save);
    while(file.exists()) {
        save = "${params.pipeline_report}" + "_" + (num++) ;
        file = new File("${params.outdir}", save);
    }
    Files.move(new File("${params.pipeline_report}").toPath(), file.toPath(), StandardCopyOption.REPLACE_EXISTING);
}

