checkpoint dRep_sample:
    input:
      concoct = "{data}/binning/{project}/{sample}/concoct_bins",
      metabat2 ="{data}/binning/{project}/{sample}/metabat2_bins",
      concoct_merged = "{data}/eukcc/{project}/{sample}/concoct/merged_bins",
      metabat2_merged = "{data}/eukcc/{project}/{sample}/metabat2/merged_bins",
      concoct_eukcc = "{data}/eukcc/{project}/{sample}/concoct/eukcc.csv",
      metabat2_eukcc = "{data}/eukcc/{project}/{sample}/metabat2/eukcc.csv",
    output:
        quality = "{data}/binrefine/{project}/{sample}/drep99/quality.csv",
        genomes = "{data}/binrefine/{project}/{sample}/drep99/genomes.txt",
        bins = directory("{data}/binrefine/{project}/{sample}/drep99/dereplicated_genomes"),
        input_bins = directory("{data}/binrefine/{project}/{sample}/drep99/input_genomes"),
    params:
        out = "{data}/binrefine/{project}/{sample}/drep99/"
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: ((attempt-1) * 20000) + 30000
    singularity:
        "docker://quay.io/biocontainers/drep:3.2.0--py_0"
    shell:
        """
        mkdir -p {params.out}
        rm -r  {params.out}
        mkdir -p {params.out}
        BD=$(pwd)
        # prepare drep quality file
        echo "genome,completeness,contamination" > {output.quality}
        cat {input.concoct_eukcc}  {input.metabat2_eukcc}|\
            awk '{{if($2 - 5*$3 >=50){{print $0}}}}' |\
            sort -k 2,3 -n |\
            tr '\\t' ',' >> {output.quality}
        if [[ $(wc -l <{output.quality}) -lt 2  ]];
        then
            touch {output.quality}
            touch {output.genomes}
            mkdir -p {output.input_bins} || true
            mkdir -p {output.bins} || true
            exit 0
        fi

        cat {output.quality}
        mkdir {output.input_bins}
        cd {output.input_bins}
        ln -s {input.concoct}/*.fa . || true
        ln -s {input.concoct_merged}/*.fa . || true
        ln -s {input.metabat2}/*.fa . || true
        ln -s {input.metabat2_merged}/*.fa . || true
        ls
        ls > ../genomes.txt
        cd ..
        cat quality.csv | grep -v "genome,completeness,contaminatio" |\
            cut -f 1 -d , > .f
        grep -f .f genomes.txt | sed -e 's/^/input_genomes\//'  > .g
        mv .g genomes.txt
        rm .f
        cd $BD

        mkdir -p {params.out}
        cd {params.out}

        cat genomes.txt
        if [[ $(wc -l <genomes.txt) -lt 2  ]];
        then
            echo "No genomes"
            cp -r input_genomes dereplicated_genomes
        else
          dRep dereplicate -p {threads} \
              . -g genomes.txt \
              -pa 0.80 -sa 0.99 -nc 0.40 \
              -cm larger \
              --genomeInfo quality.csv\
              -comp 49 -con 21
        fi
        """

checkpoint dRep_MAGs:
    input:
        fa = "{data}/MAGs/raw",
        genomes = "{data}/MAGs/genomes.txt",
        quality = "{data}/MAGs/drep_quality.csv",
    output:
        genomes = temp("{data}/MAGs/drep/genomes.txt"),
        mag_list = "{data}/MAGs/drep/MAGs.txt",
        mags = directory("{data}/MAGs/drep/dereplicated_genomes"),
    params:
        out = "{data}/MAGs/drep"
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: ((attempt-1) * 20000) + 30000
    singularity:
        "docker://quay.io/biocontainers/drep:3.2.0--py_0"
    shell:
        """
        mkdir -p {params.out}
        rm -r  {params.out}
        mkdir -p {params.out}
        BD=$(pwd)
        cat {input.genomes} | sed -e 's/^/..\/raw\//'  > {output.genomes}

        mkdir -p {params.out}
        cd {params.out}

        N=$(wc -l genomes.txt) 
        if [[ $(wc -l <genomes.txt) -lt 2  ]];
        then
            cp -r {input.fa} dereplicated_genomes
        else
          dRep dereplicate -p {threads} \
              . -g genomes.txt \
              -pa 0.80 -sa 0.95 -nc 0.30 \
              -cm larger \
              --genomeInfo {input.quality}\
              -comp 49 -con 21
        fi
        ls dereplicated_genomes > MAGs.txt
        """
