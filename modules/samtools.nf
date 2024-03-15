/*
Here are described all processes related to Samtools
Samtools is a suite of programs for interacting with high-throughput sequencing data.
See http://www.htslib.org/doc/samtools.html
*/

// Process to convert sam files into bam files
process samtools_sam2bam { 
    label 'samtools'  
    tag "$sample" 
    publishDir "${params.outdir}/Hisat2_alignments", mode: 'copy' 

    input:
        tuple val(sample), path(sam)

    output:
        tuple val(sample), path ("*.bam"), emit: tuple_sample_bam

    script:
        """
            samtools view -@ ${task.cpus} ${sam} -b -o ${sam.baseName}.bam 
        """
}

// Process to sort bam files
process samtools_sort {
    label 'samtools'
    tag "$sample" 
    publishDir "${params.outdir}/Hisat2_alignments_sorted", mode: 'copy'  

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path ("*_sorted.bam"), emit: tuple_sample_sortedbam

    script:
        """
            samtools sort -@ ${task.cpus} ${bam} -o ${bam.baseName}_sorted.bam 
        """
}