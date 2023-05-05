PROKKA(
            GUNC.out.cluster_gunc_result.filter({
                it[2].name.contains('_complete.txt')
            }).map({ cluster_name, cluster_fasta, cluster_gunc ->
                return tuple(cluster_name, cluster_fasta)
            })
        )