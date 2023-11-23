#!/usr/bin/env python3
# coding=utf-8

import argparse
import os
import sys
import pandas as pd
import logging

DEFAULT_BINNING_SOFTWARE = "mgbinrefinder_v1.0.0"
DEFAULT_BINNING_SOFTWARE_PARAMS = "default"

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
    parser.add_argument('-o', '--output', type=str, required=False, help="Script output folder", default="for_uploader")

    parser.add_argument('-m', '--mag-folder', type=str, required=True, help="MAGs storage folder")
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
    parser.add_argument('-se', '--stats-euks', type=str, required=True, help="path to eukaryotic stats output")
    parser.add_argument('-sp', '--stats-proks', type=str, required=True, help="path to prokaryotic output")

    args = parser.parse_args()

    if not os.path.exists(args.mag_folder):
        print(f'MAG folder {args.mag_folder} does not exist')
        sys.exit(1)

    if args.coassemblies is not None and not os.path.exists(args.coassemblies):
        print("Co-assembly description file does not exist")
        sys.exit(1)

    if not os.path.exists(args.assembly_software_file):
        print(f'Assembly software file {args.assembly_software_file} does not exist')
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

        self.mag_folder = args.mag_folder
        self.mag_output_dir = args.output
        if not os.path.exists(self.mag_output_dir):
            os.makedirs(self.mag_output_dir)

        self.coassembly = args.coassemblies
        if self.coassembly:
            logging.error("GGP doesn't support co-assemblies")

        self.assembly_software_file = args.assembly_software_file
        self.binning_software = args.binning_software
        self.binning_params = args.binning_params
        self.stats_euks = args.stats_euks
        self.stats_proks = args.stats_proks

    def process_mags(self):
        # names
        genomes_list, stats_software = self.get_genomes_info()
        self.output_table = pd.DataFrame({COLUMNS["genome_name"]: genomes_list})
        self.output_table.set_index(COLUMNS["genome_name"], inplace=True)

        # run accessions
        accessions = [item.split('_')[0] for item in genomes_list]
        self.output_table[COLUMNS["accessions"]] = accessions

        # assembly software
        software = self.get_assembly_software()
        self.output_table[COLUMNS["assembly_software"]] = [software[item] for item in genomes_list]

        # binning
        self.output_table[COLUMNS["binning_software"]] = [self.binning_software for _ in range(len(genomes_list))]
        self.output_table[COLUMNS["binning_parameters"]] = [self.binning_params for _ in range(len(genomes_list))]

        # stats
        self.output_table[COLUMNS["stats_generation_software"]] = stats_software
        stats = self.get_stats(self.stats_euks)
        stats.update(self.get_stats(self.stats_proks))
        self.output_table[COLUMNS["completeness"]] = [stats[item][0] for item in genomes_list]
        self.output_table[COLUMNS["contamination"]] = [stats[item][1] for item in genomes_list]

        print(self.output_table)
        #self.output_table.to_csv("output", sep='\t', index=True, header=True)

    def get_genomes_info(self):
        genomes = []
        stats_software = []
        for subfolder in os.listdir(self.mag_folder):
            cur_genomes = os.listdir(os.path.join(self.mag_folder, subfolder))
            genomes.extend(cur_genomes)
            stats_software.extend([STATS_SOFTWARE[subfolder] for _ in range(len(cur_genomes))])
        return genomes, stats_software

    def get_assembly_software(self):
        assembly_software = {}
        with open(self.assembly_software_file, 'r') as file_in:
            for line in file_in:
                line = line.strip().split('\t')
                assembly_software[line[0]] = line[1]
        return assembly_software

    def get_stats(self, input_file):
        stats = {}
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
"""

        self.euks = True if self.args.euks else False
        self.rnaDir = self.args.RNA_folder
        self.taxInfo = self.args.taxInfo
        self.binningSW = self.args.binning_sw.replace('_', ' ')
        self.binningParams = self.args.binning_pars
        self.completenessSW = self.args.completeness_sw
        self.metagenome = self.args.metagenome
        self.envVars = self.args.biomes
        self.stats = self.args.stats
        self.force = True if self.args.force else False
        self.covDir = self.args.cov_folder

        self.taxNCBI = "PATH/TO/NCBI_TAXDUMP_PREFERENTIALLY_"
        self.ar53_metadata = "PATH/TO/ar53_metadata_r214.tsv"
        self.bac120_metadata = "PATH/TO/bac120_metadata_r214.tsv"
"""

def main():
    ena_uploader = MAGupload().process_mags()
    logging.info('Completed')


if __name__ == "__main__":
    main()

