#!/usr/bin/env python

import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Propagate taxonomy from dRep cluster reps to all cluster members."
    )
    parser.add_argument(
        "-c", "--cdb", type=Path, required=True,
        help="Path to dRep Cdb.csv file (maps genomes to clusters)"
    )
    parser.add_argument(
        "-w", "--wdb", type=Path, required=True,
        help="Path to dRep Wdb.csv file (maps clusters to representatives)"
    )
    parser.add_argument(
        "-t", "--taxonomy", type=Path, required=True,
        help="Path to taxonomy file (TSV output from gtdb_to_ncbi_majority_vote.py)"
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True,
        help="Output path for the resulting TSV file"
    )
    return parser.parse_args()


def strip_fasta_ext(name):
    # peel compression suffix if present
    for comp in (".gz", ".zip"):
        if name.endswith(comp):
            name = name[: -len(comp)]
            break
    # peel fasta-like suffixes (case-insensitive)
    fasta_exts = (".fa", ".fna", ".fasta", ".fas", ".fsa")
    low = name.lower()
    for ext in fasta_exts:
        if low.endswith(ext):
            name = name[: -len(ext)]
            break
    return name


def propagate_taxonomy(cdb_path: Path, wdb_path: Path, taxonomy_path: Path, output_path: Path) -> None:
    """
    Propagate taxonomy from cluster representatives to all cluster members.

    Args:
        cdb_path (Path): Path to the Cdb.csv file from dRep (maps genomes to clusters)
        wdb_path (Path): Path to the Wdb.csv file from dRep (maps clusters to representatives)
        taxonomy_path (Path): Path to the representative taxonomy TSV
        output_path (Path): Path to write the output TSV
    """
    # Load files
    cdb = pd.read_csv(cdb_path)
    wdb = pd.read_csv(wdb_path)
    taxonomy = pd.read_csv(taxonomy_path, sep="\t")

    # Rename column in wdb to match taxonomy file
    rep2cluster = wdb[["cluster", "genome"]].rename(columns={"genome": "Genome ID"})
    
    # Strip fasta extensions from Genome IDs to match taxonomy file
    rep2cluster["Genome ID"] = rep2cluster["Genome ID"].apply(strip_fasta_ext)

    # Merge reps with taxonomy to associate each cluster with taxonomy
    cluster2tax = rep2cluster.merge(
        taxonomy,
        on="Genome ID",
        how="left"
    )

    # Merge all cluster members with cluster2tax
    merged = cdb.merge(
        cluster2tax.drop(columns=["Genome ID"]),
        left_on="secondary_cluster",
        right_on="cluster",
        how="left"
    )

    # Rename column listing cluster members to "Genome ID"
    merged = merged.rename(columns={"genome": "Genome ID"})
    merged["Genome ID"] = merged["Genome ID"].apply(strip_fasta_ext)

    # Output: name columns identical to taxonomy file
    colnames = taxonomy.columns.tolist()
    result = merged[colnames]

    result.to_csv(output_path, sep="\t", index=False)


def main():
    args = parse_args()
    propagate_taxonomy(args.cdb, args.wdb, args.taxonomy, args.output)


if __name__ == "__main__":
    main()