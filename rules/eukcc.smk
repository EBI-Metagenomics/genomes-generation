rule align:
    input:
        scaffolds = "{data}/assembly/{project}/{assembly}",
        fq1 = "{data}/reads/{project}/{sample}_1.fastq.gz",
        fq2 = "{data}/reads/{project}/{sample}_2.fastq.gz",
    output:
        bam = "{data}/bam/{project}/{sample}/reads.bam",
        bai = "{data}/bam/{project}/{sample}/reads.bam.bai",
    threads: 4
    singularity:
        "docker://microbiomeinformatics/eukcc:release-eukcc2"
    resources:
        mem_mb = lambda wildcards, attempt: ((attempt-1) * 10000) + 10000
    shell:
        """
        minimap2 -ax sr -t {threads} {input.scaffolds} {input.fq1} {input.fq2}  | 
                samtools view -q 20 -Sb - | \
                samtools sort -@ {threads} -O bam - -o {output.bam}
                samtools index {output.bam}
        """

rule linktable:
    input:
        bam = "{data}/bam/{project}/{sample}/reads.bam",
        bai = "{data}/bam/{project}/{sample}/reads.bam.bai",
        bindir = "{data}/binning/{project}/{sample}/{binner}_bins",
    output:
        links = "{data}/eukcc/{project}/{sample}/{binner}/links.csv",
    params:
        wd = "{data}/eukcc/{project}/{sample}/{binner}/",
    threads: 1
    singularity:
        "docker://microbiomeinformatics/eukcc:release-eukcc2"
    resources:
        mem_mb = lambda wildcards, attempt: ((attempt-1) * 10000) + 2000
    shell:
        """
        # if no bins
        if [[ $( ls {input.bindir}/*.fa | wc -l ) -lt 1 ]]; then
          touch {output.links}
        else
          python3 scripts/binlinks.py  --ANI 99 --within 1500 \
              --out {output.links} {input.bindir} {input.bam}
        fi
        """

rule get_EukCC_db:
    output:
        db = directory("{data}/dbs/eukcc/eukcc2_db_ver_1.1")
    params:
        tar = "{data}/dbs/eukcc/eukcc2_db_ver_1.1.tar.gz",
        d = "{data}/dbs/eukcc"
    singularity:
        "docker://microbiomeinformatics/eukcc:release-eukcc2"
    resources:
        mem_mb = lambda wildcards, attempt: ((attempt-1) * 8000) + 4000
    shell:
        """
        wget -O {params.tar} http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz
        cd {params.d}
        tar -xzvf eukcc2_db_ver_1.1.tar.gz
        """

def eukcc_db(wc):
    if "eukcc2db" in config.keys():
        return config['eukcc2db']
    else:
        return expand("{data}/dbs/eukcc/eukcc2_db_ver_1.1",
                      data = wc.data)

rule EukCC:
    input:
        db = eukcc_db,
        bindir = "{data}/binning/{project}/{sample}/{binner}_bins",
        links = "{data}/eukcc/{project}/{sample}/{binner}/links.csv",
    output:
        csv = "{data}/eukcc/{project}/{sample}/{binner}/eukcc.csv",
        bins = directory("{data}/eukcc/{project}/{sample}/{binner}/merged_bins"),
    params:
        out = "{data}/eukcc/{project}/{sample}/{binner}/",
    threads: 4
    singularity:
        "docker://microbiomeinformatics/eukcc:release-eukcc2"
    resources:
        mem_mb = lambda wildcards, attempt: ((attempt-1) * 8000) + 28000
    shell:
        """
        rm -r {params.out}workdir || true
        # if no bins
        if [[ $( ls {input.bindir} | wc -l ) -lt 1 ]]; then
          touch {output.csv}
          mkdir -p {output.bins}
          exit 0
        fi

        eukcc --debug folder \
            --improve_percent 10 \
            --n_combine 1 \
            --threads {threads} \
            --improve_ratio     5 \
            --links {input.links} \
            --min_links 100 \
            --suffix .fa \
            --db {input.db} \
            --out {params.out} \
            --prefix {wildcards.project}_{wildcards.sample}_{wildcards.binner}_merged. \
            {input.bindir}
        rm -r {params.out}/refine_workdir || true
        """

rule EukCC_MAG:
    input:
        fa = "{data}/MAGs/fa/{MAG}.fa",
        db = "{data}/dbs/eukcc/eukcc2_db_ver_1.1",
    output:
        csv = "{data}/qc/eukcc/{MAG}/eukcc.csv", tsv = temp("{data}/qc/eukcc/{MAG}/eukcc.tsv")
    params:
        out = "{data}/qc/eukcc/{MAG}/"
    threads: 4
    singularity:
        "docker://microbiomeinformatics/eukcc:release-eukcc2"
    resources:
        mem_mb = lambda wildcards, attempt: ((attempt-1) * 8000) + 28000
    shell:
        """
        rm -r {params.out}workdir || true
        eukcc --debug single \
            --threads {threads} \
            --db {input.db} \
            --out {params.out} \
            {input.fa}
        cp {output.tsv} {output.csv}
        """
