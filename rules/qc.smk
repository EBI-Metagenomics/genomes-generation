def get_reads(wc):
    """
    function to get read pairs from sample.csv that
    the user provided
    """
    for row in rows:
        if row['project'] ==wc.project and row['run'] == wc.sample:
            return({"fq1": row['fastq_1'], "fq2": row["fastq_2"]})

rule fastp:
    input:
        unpack(get_reads)
    output:
        fq1 = "{data}/reads/trimmed/{project}/{sample}/{sample}_1.fastq.gz",
        fq2 = "{data}/reads/trimmed/{project}/{sample}/{sample}_2.fastq.gz",
        html = "{data}/reads/trimmed/{project}/{sample}/{sample}.html",
        json = "{data}/reads/trimmed/{project}/{sample}/{sample}.json",
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: ((attempt-1) * 10000) + 4000
    singularity:
        "docker://nanozoo/fastp"
    shell:
        """
        fastp -w {threads} \
            --detect_adapter_for_pe \
            --in1 {input.fq1} \
            --in2 {input.fq2} \
            --html "{output.html}" \
            --json "{output.json}" \
            --out1 "{output.fq1}" \
            --out2 "{output.fq2}"
        """

rule bmtagger:
    input:
        fq1 = "{data}/reads/trimmed/{project}/{sample}/{sample}_1.fastq.gz",
        fq2 = "{data}/reads/trimmed/{project}/{sample}/{sample}_2.fastq.gz",
    output:
        fq1z = "{data}/reads/qc/{project}/{sample}/{sample}_1.fastq.gz",
        fq2z = "{data}/reads/qc/{project}/{sample}/{sample}_2.fastq.gz",
        fq1 = temp("{data}/reads/qc/{project}/{sample}/{sample}_1.fastq"),
        fq2 = temp("{data}/reads/qc/{project}/{sample}/{sample}_2.fastq"),
    params:
        out = "{data}/reads/qc/{project}/{sample}/",
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: ((attempt-1) * 10000) + 20000,
    singularity:
        #"docker://nanozoo/metawrap"
        "/hps/research/finn/saary/database/singularity/hexmek-container-metawrap.img"
    shell:
        """
        R1=$(realpath {output.fq1})
        R2=$(realpath {output.fq2})
        D=$(pwd)
        gunzip -c {input.fq1} > {output.fq1} & 
        gunzip -c {input.fq2} > {output.fq2} &
        wait
        cd {params.out}
        metawrap read_qc \
         -1 $R1 \
         -2 $R2 \
         --skip-trimming \
         --skip-pre-qc-report \
         --skip-post-qc-report \
         -t {threads} \
         -x hg38 \
         -o .
         gzip -c final_pure_reads_1.fastq > {wildcards.sample}_1.fastq.gz &
         gzip -c final_pure_reads_2.fastq > {wildcards.sample}_2.fastq.gz &
         wait

         rm final_pure_reads_* host*
         rm -r bmtagger_tmp || true

         # go back and remove trimmed files as they are not needed anymore
         cd $D
         rm {input.fq1}
         rm {input.fq2}
         """

