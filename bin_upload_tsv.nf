params.metabat_depths = "euks_depth.txt.gz"
params.genomes = "genomes"
params.assembly_software = "assembly_software.txt"
params.stats = "eukcc_final_qc.csv"
params.metagenome = ""
params.biomes = ""

include { COVERAGE_RECYCLER as COVERAGE_RECYCLER_EUK } from './modules/local/coverage_recycler/main'
include { BAT                                        } from './modules/local/cat/bat/main'
include { BAT_TAXONOMY_WRITER                        } from './modules/local/bat_taxonomy_writer/main'
include { PREPARE_UPLOAD_FILES                       } from './subworkflows/local/prepare_upload'


workflow BIN_UPLOAD_TSV {
    
    meta_genomes = [[id: "bins"], params.genomes]
    cov_genome_ch = Channel.of(meta_genomes)

    COVERAGE_RECYCLER_EUK(
        cov_genome_ch,
        params.metabat_depths
    )
   
    genomes_ch = Channel.fromPath("${params.genomes}/*.fa")

    BAT( genomes_ch, params.cat_db_folder, params.cat_taxonomy_db )

    BAT_TAXONOMY_WRITER( BAT.out.bat_names.collect() )

    taxonomy_euks  = BAT_TAXONOMY_WRITER.out.all_bin2classification
    coverage_euks  = COVERAGE_RECYCLER_EUK.out.mag_coverage.map{ meta, coverage_file -> coverage_file }.collect()
    genomes_list = Channel.fromPath("${params.genomes}/*.fa").collect()

    // RNAs is set to false as default for euks currently. skipping

    PREPARE_UPLOAD_FILES(
        genomes_list,
        [],
        params.assembly_software,
        params.stats,
        [],
        coverage_euks,
        [],
        [],
        taxonomy_euks,
        []
    )

}
