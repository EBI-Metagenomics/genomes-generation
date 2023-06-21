/*
 * fastp
*/

process FASTP {
    tag "${name}"

    publishDir "${params.outdir}/intermediate_steps/qc/fastp", mode: 'copy', pattern: "*html"
    publishDir "${params.outdir}/intermediate_steps/qc/fastp", mode: 'copy', pattern: "*json"
    publishDir "${params.outdir}/intermediate_steps/qc", mode: 'copy', pattern: "*fastq*"

    container 'quay.io/biocontainers/fastp:0.23.1--h79da9fb_0'
    label 'fastp'

    input:
    tuple val(name), path(reads)

    output:
    tuple val(name), path("${name}_fastp*.fastq.gz"), emit: output_reads
    path "*_fastp.json", emit: json
    path "*_fastp.html", emit: html

    script:
    input_reads = reads.collect()
    if ( input_reads.size() == 1 ) {
        input_reads = "--in1 ${input_reads[0]}";
        output_reads = "--out1 ${name}_fastp.fastq.gz";
    }
    else if ( input_reads.size() == 2 ) {
        input_reads = "--in1 ${input_reads[0]} --in2 ${input_reads[1]} --detect_adapter_for_pe";
        output_reads = "--out1 ${name}_fastp_1.fastq.gz --out2 ${name}_fastp_2.fastq.gz";
    }
    else {
        print('Incorrect number of reads')
        exit 1
    }

    """
    fastp -w ${task.cpus} \
    ${input_reads} \
    ${output_reads} \
    --json ${name}_fastp.json \
    --html ${name}_fastp.html \
    -l 10 -x 10 -q 15 -u 10
    """
}
