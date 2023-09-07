/*
    ~~~~~~~~~~~~~~~~~~
     Input validation
    ~~~~~~~~~~~~~~~~~~
*/

include { RENAME_MAXBIN } from '../modules/local/rename_maxbin/main'
assemblies =            channel.fromPath("/Users/kates/Desktop/EBI/pipelines/eukrecover/nextflow/workflows/*.fasta", checkIfExists: true)
/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow GGP {
    groupAssemblies = { fasta_file ->
            def cluster = fasta_file.toString().tokenize("/")[-1].tokenize(".")[0]
                def meta = [:]
                meta.id = cluster
                meta.library_layout = "paired"
                meta.single_end = false
            return tuple(meta, fasta_file)
        }
    tuple_assemblies = assemblies.map(groupAssemblies)
    tuple_assemblies.view()
    RENAME_MAXBIN(tuple_assemblies)
}
