rule binning:
    input:
        fq1 = "{data}/reads/{project}/{sample}_1.fastq.gz",
        fq2 = "{data}/reads/{project}/{sample}_2.fastq.gz",
        scaffolds = "{data}/assembly/{project}/{assembly}"
    output:
        fq1 = temp("{data}/reads/qc/{project}/{sample}/{sample}_1.fastq"),
        fq2 = temp("{data}/reads/qc/{project}/{sample}/{sample}_2.fastq"),
        concoct = directory("{data}/binning/{project}/{sample}/concoct_bins"),
        metabat2 = directory("{data}/binning/{project}/{sample}/metabat2_bins"),
        maxbin2 = directory("{data}/binning/{project}/{sample}/maxbin2_bins"),
        unbinned = directory("{data}/binning/{project}/{sample}/unbinned"),
    params:
        out = "{data}/binning/{project}/{sample}/"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: ((attempt-1) * 20000) + 20000,
        mem_GB = lambda wildcards, attempt: ((attempt-1) * 20) + 20
    #singularity:
        #"docker://nanozoo/metawrap"
    shell:
        """
        gunzip -c {input.fq1} > {output.fq1}
        gunzip -c {input.fq2} > {output.fq2}

        metawrap binning \
            -o {params.out} \
            -t {threads} \
            -l {config[mincontiglength]} \
            -m {resources.mem_GB} \
            -a {input.scaffolds} \
            --concoct --metabat2 --maxbin2 {output.fq1} {output.fq2}

        # rename bins with name
        cd {params.out}concoct_bins
        for f in *.fa ; do mv -- "$f" "{wildcards.project}_{wildcards.sample}_concoct_$f" || true ; done
        cd {params.out}metabat2_bins
        for f in *.fa ; do mv -- "$f" "{wildcards.project}_{wildcards.sample}_metabat2_$f" || true ; done
        cd {params.out}maxbin2_bins
        for f in *.fa ; do mv -- "$f" "{wildcards.project}_{wildcards.sample}_maxbin2_$f" || true ; done

        # remove unbinned bins
        cd {params.out}
        mkdir unbinned
        mv *_bins/*unbinned.fa unbinned/ || true

        # remove temp files
        cd {params.out}
        rm -r work_files
        """
