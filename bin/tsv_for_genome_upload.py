#!/usr/bin/env python3
# coding=utf-8

import argparse
import os
import sys
import pandas as pd
import logging

LIMIT_RRNA = 80
LIMIT_TRNA = 18

DEFAULT_BINNING_SOFTWARE = "MGnify-genomes-generation-pipeline_v1.0.0"
DEFAULT_BINNING_SOFTWARE_PARAMS = "default"

EUK_PATH = os.path.join("genomes_drep", "eukaryotes", "genomes")
PROK_PATH = os.path.join("genomes_drep", "prokaryotes", "genomes")

STATS_SOFTWARE = {
    "eukaryotes": "EukCC_v2.1.0",
    "prokaryotes": "CheckM2_v1.0.1"
}

COLUMNS = {
    "genome_name": "genome_name",
    "accessions": "accessions",
    "assembly_software": "assembly_software",
    "binning_software": "binning_software",
    "binning_parameters": "binning_parameters",
    "stats_generation_software": "stats_generation_software",
    "completeness": "completeness",
    "contamination": "contamination",
    "metagenome": "metagenome",
    "co-assembly": "co-assembly",
    "rRNA_presence": "rRNA_presence",
    "taxonomy_lineage": "NCBI_lineage",
    "genome_coverage": "genome_coverage",
    "genome_path": "genome_path",
    "environment_biome": "broad_environment",
    "environment_feature": "local_environment",
    "environment_material": "environmental_medium"
}

metagenomes = ["activated carbon metagenome", "activated sludge metagenome",
    "aerosol metagenome", "air metagenome", "algae metagenome", "alkali sediment metagenome",
    "amphibian metagenome", "anaerobic digester metagenome", "anchialine metagenome",
    "annelid metagenome", "ant fungus garden metagenome", "ant metagenome",
    "aquaculture metagenome", "aquatic eukaryotic metagenome", "aquatic metagenome",
    "aquatic viral metagenome", "aquifer metagenome", "ballast water metagenome",
    "bat gut metagenome", "bat metagenome", "beach sand metagenome", "beetle metagenome",
    "bentonite metagenome", "bioanode metagenome", "biocathode metagenome",
    "biofilm metagenome", "biofilter metagenome", "biofloc metagenome",
    "biogas fermenter metagenome", "bioleaching metagenome", "bioreactor metagenome",
    "bioreactor sludge metagenome", "bioretention column metagenome", "biosolids metagenome",
    "bird metagenome", "blood metagenome", "bog metagenome", "book metagenome",
    "bovine gut metagenome", "bovine metagenome", "brine metagenome", "canine metagenome",
    "cave metagenome", "cetacean metagenome", "chemical production metagenome",
    "chicken gut metagenome", "ciliate metagenome", "clay metagenome", "clinical metagenome",
    "cloud metagenome", "coal metagenome", "cold seep metagenome", "cold spring metagenome",
    "compost metagenome", "concrete metagenome", "coral metagenome", "coral reef metagenome",
    "cow dung metagenome", "crab metagenome", "crude oil metagenome",
    "Crustacea gut metagenome", "crustacean metagenome", "ctenophore metagenome",
    "decomposition metagenome", "desalination cell metagenome", "dietary supplements metagenome",
    "dinoflagellate metagenome", "drinking water metagenome", "dust metagenome",
    "ear metagenome", "echinoderm metagenome", "egg metagenome", "electrolysis cell metagenome",
    "endophyte metagenome", "epibiont metagenome", "estuary metagenome", "eukaryotic metagenome",
    "eukaryotic plankton metagenome", "eye metagenome", "factory metagenome", "feces metagenome",
    "feline metagenome", "fermentation metagenome", "fertilizer metagenome",
    "fish gut metagenome", "fishing equipment metagenome", "fish metagenome",
    "floral nectar metagenome", "flotsam metagenome", "flower metagenome",
    "food contamination metagenome", "food fermentation metagenome", "food metagenome",
    "food production metagenome", "fossil metagenome", "freshwater metagenome",
    "freshwater sediment metagenome", "frog metagenome", "fuel tank metagenome",
    "fungus metagenome", "gas well metagenome", "gill metagenome", "glacier lake metagenome",
    "glacier metagenome", "gonad metagenome", "grain metagenome", "granuloma metagenome",
    "groundwater metagenome", "gut metagenome", "halite metagenome",
    "herbal medicine metagenome", "honeybee metagenome", "honey metagenome", "horse metagenome",
    "hospital metagenome", "hot springs metagenome", "human bile metagenome",
    "human blood metagenome", "human brain metagenome", "human eye metagenome",
    "human feces metagenome", "human gut metagenome", "human hair metagenome",
    "human lung metagenome", "human metagenome", "human milk metagenome",
    "human nasopharyngeal metagenome", "human oral metagenome",
    "human reproductive system metagenome", "human saliva metagenome",
    "human semen metagenome", "human skeleton metagenome", "human skin metagenome",
    "human sputum metagenome", "human tracheal metagenome", "human urinary tract metagenome",
    "human vaginal metagenome", "human viral metagenome", "HVAC metagenome",
    "hydrocarbon metagenome", "hydrothermal vent metagenome", "hydrozoan metagenome",
    "hypersaline lake metagenome", "hyphosphere metagenome", "hypolithon metagenome",
    "ice metagenome", "indoor metagenome", "industrial waste metagenome",
    "insect gut metagenome", "insect metagenome", "insect nest metagenome",
    "internal organ metagenome", "interstitial water metagenome", "invertebrate gut metagenome",
    "invertebrate metagenome", "jellyfish metagenome", "karst metagenome", "koala metagenome",
    "lagoon metagenome", "lake water metagenome", "landfill metagenome", "leaf litter metagenome",
    "leaf metagenome", "lichen crust metagenome", "lichen metagenome", "liver metagenome",
    "lung metagenome", "macroalgae metagenome", "mangrove metagenome", "manure metagenome",
    "marine metagenome", "marine plankton metagenome", "marine sediment metagenome",
    "marsh metagenome", "marsupial metagenome", "medical device metagenome", "metagenome",
    "microbial eukaryotic metagenome", "microbial fuel cell metagenome",
    "microbial mat metagenome", "microeukaryotic metagenome", "milk metagenome",
    "mine drainage metagenome", "mine metagenome", "mine tailings metagenome",
    "mite metagenome", "mixed culture metagenome", "mollusc metagenome", "money metagenome",
    "moonmilk metagenome", "mosquito metagenome", "moss metagenome", "mouse gut metagenome",
    "mouse metagenome", "mouse skin metagenome", "mud metagenome", "museum specimen metagenome",
    "musk metagenome", "nematode metagenome", "neuston metagenome", "nutrient bag metagenome",
    "oasis metagenome", "oil field metagenome", "oil metagenome",
    "oil production facility metagenome", "oil sands metagenome", "oral metagenome",
    "oral-nasopharyngeal metagenome", "oral viral metagenome", "outdoor metagenome",
    "ovine metagenome", "oyster metagenome", "painting metagenome", "paper pulp metagenome",
    "parasite metagenome", "parchment metagenome", "peat metagenome", "periphyton metagenome",
    "permafrost metagenome", "photosynthetic picoeukaryotic metagenome", "phycosphere metagenome",
    "phyllosphere metagenome", "phytotelma metagenome", "pig gut metagenome", "pig metagenome",
    "pipeline metagenome", "pitcher plant inquiline metagenome", "placenta metagenome",
    "plant metagenome", "plastic metagenome", "plastisphere metagenome", "pollen metagenome",
    "pond metagenome", "poultry litter metagenome", "power plant metagenome", "primate metagenome",
    "probiotic metagenome", "protist metagenome", "psyllid metagenome", "rat gut metagenome",
    "rat metagenome", "reproductive system metagenome", "respiratory tract metagenome",
    "retting metagenome", "rhizoplane metagenome", "rhizosphere metagenome",
    "rice paddy metagenome", "riverine metagenome", "rock metagenome",
    "rock porewater metagenome", "rodent metagenome", "root associated fungus metagenome",
    "root metagenome", "runoff metagenome", "saline spring metagenome", "saltern metagenome",
    "salt lake metagenome", "salt marsh metagenome", "salt mine metagenome",
    "salt pan metagenome", "sand metagenome", "scorpion gut metagenome",
    "sea anemone metagenome", "seagrass metagenome", "sea squirt metagenome",
    "sea urchin metagenome", "seawater metagenome", "sediment metagenome", "seed metagenome",
    "semen metagenome", "shale gas metagenome", "sheep gut metagenome", "sheep metagenome",
    "shoot metagenome", "shrew metagenome", "shrimp gut metagenome", "silage metagenome",
    "skin metagenome", "slag metagenome", "sludge metagenome", "snake metagenome",
    "snow metagenome", "soda lake metagenome", "soda lime metagenome", "soil crust metagenome",
    "soil metagenome", "solid waste metagenome", "spider metagenome", "sponge metagenome",
    "starfish metagenome", "steel metagenome", "stomach metagenome", "stromatolite metagenome",
    "subsurface metagenome", "surface metagenome", "symbiont metagenome", "synthetic metagenome",
    "tannin metagenome", "tar pit metagenome", "termitarium metagenome",
    "termite fungus garden metagenome", "termite gut metagenome", "termite metagenome",
    "terrestrial metagenome", "tick metagenome", "tidal flat metagenome", "tin mine metagenome",
    "tobacco metagenome", "tomb wall metagenome", "tree metagenome",
    "upper respiratory tract metagenome", "urban metagenome", "urinary tract metagenome",
    "urine metagenome", "urogenital metagenome", "vaginal metagenome", "viral metagenome",
    "volcano metagenome", "wallaby gut metagenome", "wasp metagenome", "wastewater metagenome",
    "wetland metagenome", "whale fall metagenome", "whole organism metagenome", "wine metagenome",
    "Winogradsky column metagenome", "wood decay metagenome", "zebrafish metagenome"]

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Allows to create a tsv of metadata for MAG upload to ENA.")

    parser.add_argument('--debug', action='store_true', help="logging.DEBUG output")
    parser.add_argument('-o', '--output', type=str, required=False, help="Output files", default="final_table_for_uploader.tsv")

    parser.add_argument('-mp', '--mags-proks', type=str, required=False, help="Prokaryotic MAGs list of files", nargs='*')
    parser.add_argument('-me', '--mags-euks', type=str, required=False, help="Eukaryotic MAGs list of files",
                        nargs='*')
    parser.add_argument('-a', '--assembly-software-file', type=str, required=True,
                        help="File with assembler for each assembly")
    parser.add_argument('--coassemblies', type=str, default=None,
                        help="enables processing of MAGs generated from co-assemblies i.e. they won't follow the "
                             "conventional ERR[0-9] nomenclature. Takes as input a tsv with header "
                             "'name\tassembler\tcomma-separated runs'")
    parser.add_argument('-b', '--binning-software', type=str, required=False, default=DEFAULT_BINNING_SOFTWARE,
                        help="Binning software that was used for binning")
    parser.add_argument('-p', '--binning-params', type=str, required=False, default=DEFAULT_BINNING_SOFTWARE_PARAMS,
                        help="Binning parameters that were used")
    parser.add_argument('-se', '--stats-euks', type=str, required=False, help="path to eukaryotic stats output")
    parser.add_argument('-sp', '--stats-proks', type=str, required=False, help="path to prokaryotic stats output")
    parser.add_argument('--metagenome', type=str, required=True, choices=metagenomes,
                        help="choose the most appropriate metagenome "
                             "from https://www.ebi.ac.uk/ena/browser/view/408169?show=tax-tree")
    parser.add_argument('-ce', '--coverage-euks', type=str, required=False,
                        help="path to eukaryotic coverage files", nargs="*")
    parser.add_argument('-cp', '--coverage-proks', type=str, required=False, help="path to prokaryotic coverage files",
                        nargs="*")
    parser.add_argument('-rna', '--rna-outs', type=str, required=False, help="path to tRNA and rRNA .out files",
                        nargs="*")
    parser.add_argument('-te', '--tax-euks', type=str, required=False, help="path to eukaryotic taxonomy")
    parser.add_argument('-tp', '--tax-proks', type=str, required=False, help="path to prokaryotic taxonomy")
    parser.add_argument('--biomes', type=str, required=True, help="comma-separated environment parameters "
                                                                             "(biome,feature,material)")
    parser.add_argument('--absolute-path', type=str, required=True, help="Absolute path to result folder of pipeline")

    args = parser.parse_args()

    if args.coassemblies is not None and not os.path.exists(args.coassemblies):
        print("Co-assembly description file does not exist")
        sys.exit(1)

    if not os.path.exists(args.assembly_software_file):
        print(f'Assembly software file {args.assembly_software_file} does not exist')
        sys.exit(1)

    if len(args.biomes.split(',')) != 3:
        print(f'Environment variables must be 3: biome, feature, and material. Got {args.biomes}')
        sys.exit(1)

    return args


class MAGupload:
    def __init__(self):
        args = parse_args()
        logLevel = logging.INFO

        if args.debug:
            logLevel = logging.DEBUG
        logging.basicConfig(
            format="%(asctime)s %(message)s",
            datefmt="%d-%m-%Y %H:%M:%S: ",
            level=logLevel,
        )

        self.euk_mag = args.mags_euks
        self.prok_mag = args.mags_proks

        self.coassembly = args.coassemblies
        if self.coassembly:
            logging.error("GGP doesn't support co-assemblies")

        self.assembly_software_file = args.assembly_software_file
        self.binning_software = args.binning_software
        self.binning_params = args.binning_params
        self.stats_euks = args.stats_euks if args.stats_euks else None
        self.stats_proks = args.stats_proks if args.stats_proks else None
        self.metagenome = args.metagenome
        self.biomes = args.biomes.split(',')
        self.coverage_euks = args.coverage_euks if args.coverage_euks else None
        self.coverage_proks = args.coverage_proks if args.coverage_proks else None
        self.rna = args.rna_outs if args.rna_outs else None
        self.tax_euks = args.tax_euks if args.tax_euks else None
        self.tax_proks = args.tax_proks if args.tax_proks else None
        self.absolute_path = args.absolute_path
        self.output_file = args.output

    def process_mags(self):
        # genomes
        genomes_list, stats_software, genome_paths = self.get_genomes_info()
        self.output_table = pd.DataFrame({COLUMNS["genome_name"]: genomes_list})
        self.output_table.set_index(COLUMNS["genome_name"], inplace=True)
        self.output_table[COLUMNS["genome_path"]] = genome_paths
        # run accessions
        accessions = [item.split('_')[0] for item in genomes_list]
        self.output_table[COLUMNS["accessions"]] = accessions

        # assembly software
        self.output_table[COLUMNS["assembly_software"]] = self.get_assembly_software(genomes_list)

        # binning
        self.output_table[COLUMNS["binning_software"]] = [self.binning_software for _ in range(len(genomes_list))]
        self.output_table[COLUMNS["binning_parameters"]] = [self.binning_params for _ in range(len(genomes_list))]

        # stats
        self.output_table[COLUMNS["stats_generation_software"]] = stats_software
        stats = self.get_stats(self.stats_euks)
        stats.update(self.get_stats(self.stats_proks))
        self.output_table[COLUMNS["completeness"]] = [stats[item][0] for item in genomes_list]
        self.output_table[COLUMNS["contamination"]] = [stats[item][1] for item in genomes_list]

        # coverage
        coverage = self.get_coverage(self.coverage_proks)
        coverage.update(self.get_coverage(self.coverage_euks))

        self.output_table[COLUMNS["genome_coverage"]] = [coverage[item] for item in genomes_list]

        # study info
        self.output_table[COLUMNS["metagenome"]] = [self.metagenome for _ in range(len(genomes_list))]
        self.output_table[COLUMNS["co-assembly"]] = ["False" for _ in range(len(genomes_list))]
        self.output_table[COLUMNS["environment_biome"]] = [self.biomes[0] for _ in range(len(genomes_list))]
        self.output_table[COLUMNS["environment_feature"]] = [self.biomes[1] for _ in range(len(genomes_list))]
        self.output_table[COLUMNS["environment_material"]] = [self.biomes[2] for _ in range(len(genomes_list))]
        self.output_table[COLUMNS["rRNA_presence"]] = self.get_rna(genomes_list)
        taxonomy, unclassified = self.get_taxonomy(genomes_list)
        self.output_table[COLUMNS["taxonomy_lineage"]] = taxonomy
        # Remove ".fa" from each value in the specified column
        self.output_table.index = self.output_table.index.str.replace(".fa", "")

        # Remove unclassified mags
        if unclassified:
            print(f"You have {len(unclassified)} unclassified genomes")
            self.output_table = self.output_table.drop(unclassified)
        # output to file
        self.output_table.to_csv(self.output_file, sep='\t', index=True, header=True)

    def get_genomes_info(self):
        genomes, stats_software, paths = [[] for _ in range(3)]
        if self.euk_mag:
            genomes.extend([os.path.basename(i).replace('.gz', '') for i in self.euk_mag])
            stats_software.extend([STATS_SOFTWARE["eukaryotes"] for _ in range(len(self.euk_mag))])
            paths.extend([os.path.join(self.absolute_path, EUK_PATH, os.path.basename(i)) for i in self.euk_mag])
        if self.prok_mag:
            genomes.extend([os.path.basename(i).replace('.gz', '') for i in self.prok_mag])
            stats_software.extend([STATS_SOFTWARE["prokaryotes"] for _ in range(len(self.prok_mag))])
            paths.extend([os.path.join(self.absolute_path, PROK_PATH, os.path.basename(i)) for i in self.prok_mag])
        return genomes, stats_software, paths

    def get_assembly_software(self, genomes):
        assembly_software = {}
        software_list = []
        with open(self.assembly_software_file, 'r') as file_in:
            for line in file_in:
                line = line.strip().split('\t')
                assembly_software[line[0]] = line[1]
        for genome in genomes:
            name = genome.split('_')[0]
            software_list.append(assembly_software[name])
        return software_list

    def get_stats(self, input_file):
        stats = {}
        if not input_file:
            return stats
        with open(input_file, 'r') as file_in:
            for line in file_in:
                if "completeness" in line:
                    continue
                line = line.strip().split(',')
                comp = round(float(line[1]), 2)
                if comp == 100.0:
                    comp = 100
                cont = round(float(line[2]), 2)
                stats[line[0]] = [comp, cont]
        return stats

    def get_coverage(self, files):
        coverage = {}
        if not files:
            return coverage
        for coverage_file in files:
            with open(coverage_file, 'r') as file_in:
                for line in file_in:
                    line = line.strip().split('\t')
                    coverage[line[0]] = round(float(line[1]), 2)
        return coverage

    def get_rna(self, genomes):
        rrna, trna = {}, {}
        if not self.rna:
            return ['False' for _ in genomes]
        for filename in self.rna:
            if filename.endswith('_rRNAs.out'):
                genome = filename.split('_rRNAs.out')[0] + '.fa'
                rrna[genome] = self.check_rna(filename, LIMIT_RRNA)
            if filename.endswith('_tRNA_20aa.out'):
                genome = filename.split('_tRNA_20aa.out')[0] + '.fa'
                trna[genome] = self.check_rna(filename, LIMIT_TRNA)
        final_decision = []
        for genome in genomes:
            if genome not in trna or genome not in rrna:
                final_decision.append('False')
            else:
                if rrna[genome] and trna[genome]:
                    final_decision.append('True')
                else:
                    final_decision.append('False')
        return final_decision

    def check_rna(self, filename, limit):
        with open(filename, 'r') as file_in:
            for line in file_in:
                cur_rna_count = float(line.strip().split('\t')[1])
                if cur_rna_count < limit:
                    return False
        return True

    def process_tax_file(self, filename, type):
        lineage = {}
        unclassified = []
        with open(filename, 'r') as file_in:
            for line in file_in:
                if 'classification' in line:
                    continue
                else:
                    line = line.strip().split('\t')
                    if type == 'euks':
                        lineage[line[0]] = line[3]
                        if line[3] == 'Unclassified':
                            unclassified.append(line[0].replace('.fa', ''))
                    else:
                        lineage[line[0] + '.fa'] = line[2]
                        if line[2] == 'Unclassified':
                            unclassified.append(line[0])
        return lineage, unclassified

    def get_taxonomy(self, genomes):
        lineage = {}
        unclassified = []
        if self.tax_euks:
            lineage_euks, unclassified_euks = self.process_tax_file(self.tax_euks, type='euks')
            lineage.update(lineage_euks)
            unclassified.extend(unclassified_euks)
        if self.tax_proks:
            lineage_proks, unclassified_proks = self.process_tax_file(self.tax_proks, type='proks')
            lineage.update(lineage_proks)
            unclassified.extend(unclassified_proks)
        final_tax = [lineage[genome] for genome in genomes]
        return final_tax, unclassified


def main():
    MAGupload().process_mags()
    logging.info('Completed')


if __name__ == "__main__":
    main()

