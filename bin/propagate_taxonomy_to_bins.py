#!/usr/bin/env python

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Propagate taxonomy from dRep cluster reps to all cluster members.")
    parser.add_argument(
        "-c", "--cdb", type=Path, required=True,
        help="Path to dRep Cdb.csv file"
    )
    parser.add_argument(
        "-t", "--taxonomy", type=Path, required=True,
        help="Path to taxonomy file (TSV with representative and lineage)"
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True,
        help="Output path for the resulting TSV file"
    )
    return parser.parse_args()


def propagate_taxonomy(cdb_path: Path, taxonomy_path: Path, output_path: Path) -> None:
    """
    Propagate taxonomy from cluster representatives to all cluster members.

    Args:
        cdb_path (Path): Path to the Cdb.csv file from dRep
        taxonomy_path (Path): Path to the representative taxonomy TSV
        output_path (Path): Path to write the output TSV
    """
    # Load Cdb.csv
    cdb = pd.read_csv(cdb_path)

    # Ensure required columns exist
    required_cols = {'genome', 'cluster', 'cluster_rep'}
    if not required_cols.issubset(cdb.columns):
        raise ValueError(f"Cdb file must contain columns: {required_cols}")

    # Load taxonomy file
    taxonomy = pd.read_csv(taxonomy_path, sep='\t', names=["representative", "lineage"], header=None)

    # Extract representatives and their cluster
    reps = cdb[cdb['cluster_rep']][['cluster', 'genome']].rename(columns={'genome': 'representative'})

    # Merge representatives with their taxonomy
    rep_tax = reps.merge(taxonomy, on='representative', how='left')

    if rep_tax['lineage'].isnull().any():
        missing = rep_tax[rep_tax['lineage'].isnull()]
        raise ValueError(f"Missing taxonomy for representatives:\n{missing}")

    # Merge full Cdb with taxonomy via cluster
    merged = cdb.merge(rep_tax[['cluster', 'lineage']], on='cluster', how='left')

    # Output: genome (member), lineage
    result = merged[['genome', 'lineage']].rename(columns={'genome': 'cluster_member'})

    # Write output
    result.to_csv(output_path, sep='\t', index=False)


def main():
    args = parse_args()
    propagate_taxonomy(args.cdb, args.taxonomy, args.output)


if __name__ == "__main__":
    main()