def busco_db(wc):
    if "buscodb" in config.keys():
        return config['buscodb']
    else:
        return expand("{data}/dbs/busco",
                      data = wc.data)

rule BUSCO_DB:
    output:
        directory("{data}/dbs/busco")
    params:
        d =  "{data}/dbs/busco/"
    shell:
        """
        mkdir -p {params.d}
        cd {params.d}
        wget -r -np -R "index.html*"  https://busco-data.ezlab.org/v5/data/
        mv busco-data.ezlab.org/v5/data/* .
        rm -r busco-data.ezlab.org
        #extract all gzipped files
        cd lineages
        find . -type f -iname "*tar.gz" -exec rm "{{}}" \;
        cd ../placement_files
        find . -type f -iname "*tar.gz" -exec rm "{{}}" \;
        cd ../information
        find . -type f -iname "*tar.gz" -exec rm "{{}}" \;
        """

rule BUSCO_MAG:
    input:
      db = busco_db,
      fa = "{data}/MAGs/fa/{MAG}.fa",
    output:
        summary = "{data}/qc/busco/{MAG}/short_summary.specific.txt"
    params:
        d = "{data}/qc/busco/{MAG}/"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: ((attempt-1) * 6000) + 28000
    singularity:
        "docker://ezlabgva/busco:v5.2.2_cv1"
    shell:
        """
        FNA=$(realpath {input.fa})
        DB=$(realpath {input.db})
        mkdir -p {params.d}
        if [ ! -z "$(ls -A {params.d})" ]; then
          echo "Folder not empty. Deleting"
          rm -rf {params.d}
        fi
        mkdir -p {params.d}
        cd {params.d}
        busco  --offline \
              -i $FNA \
              -m genome \
              -o out \
              --auto-lineage-euk \
              --download_path $DB \
              -c {threads} || touch out/short_summary.specific_failed.out.txt

        cp out/short_summary.specific*.out.txt short_summary.specific.txt
        rm -r out
        """

