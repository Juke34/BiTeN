/*
Here are described all processes related to FastQC
FastQC is a quality control tool for high throughput sequence data.
See https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
*/

// This process run fastqc with default parameters
process fastqc { 
    label 'fastqc'
    tag "$sample_id" 
    publishDir "${params.outdir}/FastQC", mode: 'copy'  

    input: 
    tuple val(sample_id), path(reads)

    output: 
    path ("fastqc_${sample_id}_logs")

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -t ${task.cpus} -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """

}
