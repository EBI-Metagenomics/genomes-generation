rule metaSPAdes:
    input:
        fq1 = "{data}/reads/qc/{project}/{sample}/{sample}_1.fastq.gz",
        fq2 = "{data}/reads/qc/{project}/{sample}/{sample}_2.fastq.gz",
    output:
      scaffolds = "{data}/assembly/spades/{project}/{sample}/spades_output/scaffolds.fasta",
      ver = "{data}/assembly/spades/{project}/{sample}/spades_output/version.txt"
    params:
        out = "{data}/assembly/spades/{project}/{sample}/spades_output"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: ((attempt-1) * 100000) + 301000,
        mem_gb = lambda wildcards, attempt: ((attempt-1) * 100) + 300
    singularity:
        "docker://microbiomeinformatics/spades:3.15.3"
    shell:
        """
        if [[ -f {params.out}/spades.log ]]; then
          # lets continue where we stopped last
          spades.py --continue -o {params.out} 
        else
          spades.py  --version > {output.ver}
          spades.py -1 {input.fq1} \
              -2 {input.fq2} \
              --meta -o {params.out} \
              -t {threads} -m  {resources.mem_gb}
            fi
        rm -r {params.out}/K[0-9]*
        rm -r {params.out}/corrected
        """
