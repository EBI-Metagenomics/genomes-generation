checkpoint binrefine:
    input:
        concoct = "{data}/binning/{project}/{sample}/concoct_bins",
        metabat2 = "{data}/binning/{project}/{sample}/metabat2_bins",
        #maxbin2 = "{data}/binning/{project}/{sample}/maxbin2_bins",
    output:
        out = directory("{data}/binrefine/{project}/{sample}/metawrap/metawrap_bins")
    params:
        out = "{data}/binrefine/{project}/{sample}/metawrap/",
        d = "{data}/binrefine/{project}/{sample}/metawrap/metawrap_50_10_bins"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: ((attempt-1) * 10000) + 40000,
        mem_GB = lambda wildcards, attempt: ((attempt-1) * 10) + 40
    #singularity:
        #"docker://nanozoo/metawrap"
    shell:
        """
        #rm -r {params.out}
        #mkdir -p {params.out}
        metawrap bin_refinement -o {params.out} \
            -A {input.concoct} -B {input.metabat2} \
            -m {resources.mem_GB} \
            -t {threads}\
            -c 50 -x 10 \
            --quick  || true

        mv {params.d} {output.out}
        cd {params.out}metawrap_bins
        for f in *.fa ; do mv -- "$f" "{wildcards.project}_{wildcards.sample}_metawrap_$f" || true ; done
        """

