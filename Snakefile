import csv
import shutil


configfile: "config/config.yaml"
config['data'] = os.path.abspath(config['data'])

include: "rules/qc.smk"
include: "rules/spades.smk"
include: "rules/binning.smk"
include: "rules/eukcc.smk"
include: "rules/busco.smk"
include: "rules/drep.smk"
include: "rules/binrefine.smk"
include: "rules/cmseq.smk"

# loading samples from sample table
# very simple can be replaced as you see fit
rows = []
with open(config['sample_table']) as fin:
    for row in csv.DictReader(fin, delimiter=config['delimiter']):
        rows.append(row)

clean_reads = expand(expand("{data}/reads/qc/{{project}}/{{sample}}/{{sample}}_1.fastq.gz",
    data = config['data']), zip,
    project = [x["project"] for x in rows],
    sample = [x["run"] for x in rows])

assembly_files = expand(expand("{data}/assembly/spades/{{project}}/{{sample}}/spades_output/scaffolds.fasta",
    data = config['data']), zip,
    project = [x["project"] for x in rows],
    sample = [x["run"] for x in rows])

bin_files = expand(expand("{data}/binning/{{project}}/{{sample}}/concoct_bins",
    data = config['data']), zip,
    project = [x["project"] for x in rows],
    sample = [x["run"] for x in rows])

proka_bins = expand(expand("{data}/binrefine/{{project}}/{{sample}}/metawrap/metawrap_bins",
    data = config['data']), zip,
    project = [x["project"] for x in rows],
    sample = [x["run"] for x in rows])

coverage_files =  expand(expand(
    "{data}/stats/{{project}}/{{sample}}/coverage.csv",
    data = config['data']), zip,
    project = [x["project"] for x in rows],
    sample = [x["run"] for x in rows])

final_summary = expand("{data}/MAGs/qc.csv",
              data = config['data'])

rule all:
    input:
      reads = clean_reads,
      assembly = assembly_files,
      bins =  bin_files,
      prokaryotes = proka_bins,
      coverage = coverage_files,
      #summary = final_summary

def eMAGs_per_sample(wc):
    mag_folders = []
    return expand(expand("{data}/binrefine/{{project}}/{{sample}}/drep99/dereplicated_genomes",
          data = wc.data), zip,
          project = [x["project"] for x in rows],
          sample = [x["run"] for x in rows])

def qual_per_sample(wc):
    mag_folders = []
    return expand(expand("{data}/binrefine/{{project}}/{{sample}}/drep99/quality.csv",
          data = wc.data), zip,
          project = [x["project"] for x in rows],
          sample = [x["run"] for x in rows])

def MAGs_per_sample(wc):
    # euk MAGs are in the folder
    mag_folders = []
    mag_folders.extend(expand(expand("{data}/binrefine/{{project}}/{{sample}}/drep99/dereplicated_genomes",
          data = wc.data), zip,
          project = [x["project"] for x in rows],
          sample = [x["run"] for x in rows]))
    mag_folders.extend(expand(expand("{data}/binrefine/{{project}}/{{sample}}/metawrap/metawrap_bins",
          data = wc.data), zip,
          project = [x["project"] for x in rows],
          sample = [x["run"] for x in rows]))
    return mag_folders

def both_mags(wc):
    from glob import glob
    fas = []
    folder_euks = checkpoints.dRep_sample.get(data = wc.data, project = wc.project, sample = wc.sample).output.bins
    folder_bacs = checkpoints.binrefine.get(data = wc.data, project = wc.project, sample = wc.sample).output.out
    for p in [folder_euks, folder_bacs]:
        fas.extend(glob(os.path.join(p,"*.fa")))
    return fas

def coverage_per_run(wc):
    from glob import glob
    folder = checkpoints.get_both_mags.get(data = wc.data, project = wc.project, sample = wc.sample).output.bins
    fasta = [os.path.basename(x).rsplit(".", 1)[0] for x in glob(os.path.join(folder, "*.fa"))]
    return expand("{data}/stats/{project}/{sample}/coverage/cov_{bin}.csv",
            data = wc.data,
            project = wc.project,
            sample = wc.sample,
            bin=fasta)

checkpoint get_both_mags:
    input:
        both_mags
    output:
        bins = directory("{data}/stats/{project}/{sample}/bins")
    shell:
        """
        mkdir -p {output}
        cp {input} {output}
        """

rule coverage:
    input:
        coverage_per_run
    output:
        "{data}/stats/{project}/{sample}/coverage.csv"
    shell:
        """
        touch {output}
        """


rule aggregate_MAG:
    input:
        dirs = eMAGs_per_sample,
        quality =  qual_per_sample
    output:
        fa = directory("{data}/MAGs/raw"),
        genomes = "{data}/MAGs/genomes.txt",
        quality = "{data}/MAGs/drep_quality.csv",
    run:
        import shutil
        import glob
        import os
        os.makedirs(output.fa, exist_ok=True)
        with open(output.genomes, "w") as fout:
          for folder in input.dirs:
            for fa in glob.glob("{}/*.fa".format(folder)):
                shutil.copy(fa, output.fa)
                fout.write("{}\n".format(os.path.basename(fa)))
        with open(output.quality, "w") as fout:
            cout = csv.DictWriter(fout, fieldnames=['genome','completeness', 'contamination'])
            cout.writeheader()
            for f in input.quality:
                with open(f) as fin:
                    for row in csv.DictReader(fin):
                        cout.writerow(row)

rule rename_mag:
    input:
        "{data}/MAGs/raw/{MAG}.fa"
    output:
        "{data}/MAGs/fa/{MAG}.fa"
    shell:
        """
        python3 scripts/rename_multifasta_prefix.py --clean \
            -p {wildcards.MAG} \
            -f {input} > {output}
        """

def eukcc_files(wc):
    mags = []
    with checkpoints.dRep_MAGs.get(data = wc.data).output.mag_list.open() as f:
        for line in f:
            mags.append(line.strip().replace(".fa",""))
    f= expand("{data}/qc/eukcc/{MAG}/eukcc.csv",
          data = wc.data,
          MAG = mags)
    return(f)

def busco_files(wc):
    mags = []
    with checkpoints.dRep_MAGs.get(data = wc.data).output.mag_list.open() as f:
        for line in f:
            mags.append(line.strip().replace(".fa",""))
    return expand("{data}/qc/busco/{MAG}/short_summary.specific.txt",
          data = wc.data,
          MAG = mags)

def mag_files(wc):
    mags = []
    with checkpoints.dRep_MAGs.get(data = wc.data).output.mag_list.open() as f:
        for line in f:
            mags.append(line.strip().replace(".fa",""))
    return expand("{data}/MAGs/fa/{MAG}.fa",
          data = wc.data,
          MAG = mags)

rule QC_table:
    input:
        eukcc = eukcc_files,
        busco = busco_files,
        raw = "{data}/MAGs/raw",
        mag_list = "{data}/MAGs/drep/MAGs.txt",
        fna =  mag_files,
    output:
        csv = "{data}/MAGs/qc.csv"
    shell:
        """
        python3 scripts/create_qc_table.py \
            --mags {input.raw} \
            --qc_dir {wildcards.data}/qc \
            --output {output.csv}
        """
