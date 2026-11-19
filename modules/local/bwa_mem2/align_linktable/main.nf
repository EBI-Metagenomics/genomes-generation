process ALIGNMENT_LINKTABLE {

    /*
    This module aligns reads to the reference using specified arguments and produces output in the form of FASTQ.GZ, BAM and BAM.BAI files.
    It also runs:
     - eukcc binlinks
    */

    label 'process_medium'

    tag "${meta.id} align to ${ref_fasta}"

    container 'quay.io/microbiome-informatics/bwa_eukcc:2.3_2.3.1'

    input:
    tuple val(meta), path(reads), path(ref_fasta), path(bins), path(depth)

    output:
    tuple val(meta), path("*.links.csv"), emit: links_table
    tuple val(meta), path("*.idxstats") , emit: idxstats
    path "versions.yml"                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_args = task.ext.alignment_args

    """
    mkdir -p output

    echo " ---> index fasta"
    bwa-mem2 index ${ref_fasta}

    echo " ---> mapping files to assembly"
    bwa-mem2 mem -M \
      -t ${task.cpus} \
      ${ref_fasta} \
      ${reads} | \
    samtools view -@ ${task.cpus} ${samtools_args} - | \
    samtools sort -@ ${task.cpus} -O bam - -o output/${meta.id}_sorted.bam

    echo "samtools index sorted bam"
    samtools index -@ ${task.cpus} output/${meta.id}_sorted.bam

    echo " ---> samtools idxstats sorted bam"
    samtools idxstats --threads ${task.cpus} output/${meta.id}_sorted.bam > ${prefix}.idxstats

    echo "linktable"
    mkdir -p bins
    cp ${bins}/* bins/
    rm -f bins/*unbinned*

    export BINS=\$(ls bins | wc -l)
    if [ \$BINS -eq 0 ]; then
        echo "Bins directory is empty"
        touch ${meta.id}.links.csv
    else
        binlinks --ANI 99 \
        --within 1500 \
        --out ${meta.id}.links.csv \
        --debug \
        --binfolder bins \
        --alignment output/${meta.id}_sorted.bam
    fi

    rm -rf output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2> /dev/null)
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
        EukCC: \$(eukcc -v | grep -o '[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+')
    END_VERSIONS
    """
}
