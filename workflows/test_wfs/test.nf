/*
    ~~~~~~~~~~~~~~~~~~
     Input validation
    ~~~~~~~~~~~~~~~~~~
*/
accession = channel.value("SRR6311991")
fasta = channel.fromPath("/Users/kates/Desktop/EBI/pipelines/eukrecover/tests/resources/ERR.fasta")
cm_model = channel.fromPath("/Users/kates/Desktop/EBI/pipelines/eukrecover/tests/resources/empty_file")
/*
    ~~~~~~~~~~~~~~~~~~
     DBs
    ~~~~~~~~~~~~~~~~~~
*/

/*
    ~~~~~~~~~~~~~~~~~~
     Steps
    ~~~~~~~~~~~~~~~~~~
*/
include { DETECT_RRNA } from '../../modules/detect_rrna'

/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/



workflow {
    DETECT_RRNA(accession, fasta, cm_model )
}
