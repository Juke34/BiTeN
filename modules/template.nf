/*
Here are described all processes related to XXX
Explain what does the tool
Link to tool documentation
*/

// Explain the process
process process_name { // A process is the basic processing primitive to execute a user script. The process tag is mandatory followed by the name of the process. See https://www.nextflow.io/docs/latest/process.html#processes
    label 'mylabel'    // allows the annotation of processes with mnemonic identifier of your choice. Can be use to link resource (tool/cpu) from config files. See https://www.nextflow.io/docs/latest/process.html#label
    tag "$sample"      // to associate each process execution with a custom label, so that it will be easier to identify them in the log file or in the trace execution report. See https://www.nextflow.io/docs/latest/process.html#tag
    publishDir "${params.outdir}/MychoosenName", mode: 'copy'  // to publish the process output files to a specified folder. See https://www.nextflow.io/docs/latest/process.html#publishdir

    input:             // The input block allows you to define the input channels of a process, similar to function arguments. See https://www.nextflow.io/docs/latest/process.html#inputs
        tuple val(sample), path(bam) // The tuple qualifier allows you to group multiple values into a single input definition. See https://www.nextflow.io/docs/latest/process.html#input-type-tuple

    output:             // The output block allows you to define the output channels of a process, similar to function outputs. https://www.nextflow.io/docs/latest/process.html#outputs
        tuple val(sample), path ("*_matchOutputFileOrFolder"), emit: tuple_sample_sortedbam // The tuple qualifier allows you to output multiple values in a single channel. See https://www.nextflow.io/docs/latest/process.html#output-type-tuple
                                                                                            // The emit option can be added to the process output definition to assign a name identifier. This name can be used to reference the channel within the caller scope. See https://www.nextflow.io/docs/edge/dsl2.html#process-named-output 
    script:             // The script block defines, as a string expression, the script that is executed by the process. See https://www.nextflow.io/docs/edge/process.html#script
        """
            toolname [ -option ] [ arguments ]
        """
} // End of the process