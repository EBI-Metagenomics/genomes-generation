rule cmseq_cov:
  input:
      fasta = "{data}/stats/{project}/{sample}/bins/{bin}.fa",
      fq1 = "{data}/reads/qc/{project}/{sample}/{sample}_1.fastq.gz",
      fq2 = "{data}/reads/qc/{project}/{sample}/{sample}_2.fastq.gz",
  output:
      cov = "{data}/stats/{project}/{sample}/coverage/cov_{bin}.csv",
      bam = temp("{data}/stats/{project}/{sample}/bins/{bin}.bam"),
      bai = temp("{data}/stats/{project}/{sample}/bins/{bin}.bam.bai"),
  threads: 4
  params:
      mincov = 1
  singularity:
      "docker://saardock/cmseq:latest"
  resources:
      mem_mb = lambda wildcards, attempt: ((attempt-1) * 10000) + 10000
  shell:
    """
      minimap2 -ax sr -t {threads} {input.fasta} {input.fq1} {input.fq2}  | 
              samtools view -q 20 -Sb - | \
              samtools sort -@ {threads} -O bam - -o {output.bam}
              samtools index {output.bam}

      # cmseq
      breadth_depth.py \
          --combine \
          --mincov {params.mincov} \
          {output.bam} > {output.cov}
    """
