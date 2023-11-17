import argparse
import os
import sys
import glob
import re
import json
import csv
import pandas as pd
from ena_portal_api.ena_handler import EnaApiHandler
from tqdm import tqdm
import xml.dom.minidom as minidom
import requests
import gtdb_to_ncbi_majority_vote_v2

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

handler = EnaApiHandler()
run_identifier = re.compile("([E|S|D]R[R|S]\d{6,})")

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

def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description="Allows to create a tsv of metadata for MAG upload to ENA.")
    
    parser.add_argument('-m', '--mag_folder', type=str, required=True, help="MAGs storage folder")

    parser.add_argument('--euks', action='store_true', help="computes xml for eukaryotic MAGs")
    parser.add_argument('--stats', type=str, default=None, help="path to checkM/EukCC output")
    parser.add_argument('--binning_sw', type=str, required=True, help="binning software (format: software_vX.X)")
    parser.add_argument('--binning_pars', type=str, required=True, help="parameters used with binning software ",
        default="MaxBin2, MetaBat2, Concoct with default parameter of the metaWRAP pipeline. Bin refinement " +
        "module used from metaWRAP with default parameters")
    parser.add_argument('--completeness_sw', type=str, required=True, help="software used to determine MAG " +
        "completeness (checkM/EukCC)")
    parser.add_argument('--RNA_folder', type=str, help="folder where RNA detection output " +
        "is stored")
    parser.add_argument('--taxInfo', type=str, required=True, help="either GTDB output folder, or eukaryotes' " + 
        "taxonomy file")
    parser.add_argument('--metagenome', type=str, required=True, help="choose the most appropriate metagenome " +
        "from https://www.ebi.ac.uk/ena/browser/view/408169?show=tax-tree")
    parser.add_argument('--biomes', type=str, required=True, help="comma-separated environment parameters " +
        "(biome,feature,material)")
    parser.add_argument('--cov_folder', type=str, required=True, help="folder for coverage depth")
    parser.add_argument('--coassemblies', type=str, default=None, help="enables processing " +
        "of MAGs generated from co-assemblies i.e. they won't follow the conventional ERR[0-9]"  +
        "nomenclature. Takes as input a tsv with header 'name\tassembler\tcomma-separated runs'")
    
    parser.add_argument('--force', action='store_true', help="forces reset of sample xml's backups")
    
    args = parser.parse_args(argv)

    if not os.path.exists(args.mag_folder):
        print('MAG folder "{}" does not exist'.format(args.mag_folder))
        sys.exit(1)

    if not args.euks and not os.path.exists(args.RNA_folder):
        print("The specified RNA_folder path does not exist. If RNA detection " + 
            "has not been computed, please refer to https://www.ebi.ac.uk/seqdb/" +
            "confluence/display/MetaGen/MAG+upload+workflow+in+Codon")
        sys.exit(1)

    if args.coassemblies is not None and not os.path.exists(args.coassemblies):
        print("Co-assembly description file does not exist")
        sys.exit(1)

    if not os.path.exists(args.taxInfo):
        print("The specified taxonomic annotation file does not exist")
        sys.exit(1)
    else:
        isFolder = os.path.isdir(args.taxInfo)
        if isFolder:
            if args.euks:
                print("Eukaryotic upload requires a file describing taxonomy (not a folder)")
                sys.exit(1)
            else:
                if not os.path.isfile(os.path.join(args.taxInfo, "gtdbtk.bac120.summary.tsv")):
                    print("Specified taxInfo folder does not contain any GTDB annotation file.")
                    sys.exit(1)
        elif not isFolder and not args.euks:
            print("Prokaryotic upload requires GTDB output folder (not a file)")
            sys.exit(1)

    if args.stats is not None and not os.path.exists(args.stats):
        print("The specified stats path does not exist")
        sys.exit(1)

    if not "_v" in args.binning_sw:
        print("Incorrect binning software format [software_vX.X]")
        sys.exit(1)

    if not re.search(",*,", args.biomes):
        print("Unrecognised environment variables [format: 'biome, feature, material']")
        sys.exit(1)
            
    return args

def extract_MAG_info(magDir, coassembly, binner="", binParams=""):
    fastaList = glob.glob(os.path.join(magDir, "*.gz"))

    # ASSUMING THIS IS STILL THE MAG NAME STRUCTURE!!
    allMAGpaths = [f for f in fastaList if re.search(r'_?bin.[0-9]+.fa(sta)?.gz', f)]

    # extract all coassembly name-to-runs mapping from file if file != None
    coassDict = {}
    if coassembly:
        with open (coassembly, 'r') as f:
            lines = f.readlines()[1:]
            for line in lines:
                mapping = line.strip().split('\t')
                coassName = mapping[0]
                coassembler = mapping[1]
                coassRuns = mapping[2].split(',')
                coassDict[coassName] = {}
                coassDict[coassName]["runs"] = coassRuns
                coassDict[coassName]["assembler"] = coassembler

    coassMAGPaths = [path for coassName in coassDict for path in allMAGpaths if coassName in os.path.basename(path)]

    # all the following information inherits run, bin number, and co-assembly=True/False from the MAG name
    MAGinfoDict, MAGruns = {}, []
    missingRuns, multipleRuns = [], []
    for p in allMAGpaths:
        isCoassembly = False
        fileName = os.path.basename(p)
        bin = (re.findall("bin.[0-9]+", fileName)[-1]).replace("bin.", "") # if we need it, the bin number
        
        # MAGprefix is an identifier I came up with for unique genome name
        # It's the run for single-assembly MAGs and the co-assembly name otherwise

        if p in coassMAGPaths:
            MAGprefix = fileName.split("bin")[0].rstrip('_')
            run = coassDict[MAGprefix]["runs"]
            isCoassembly = True
        else:
            run = run_identifier.findall(p)
            if len(run) == 0:
                missingRuns.append(os.path.basename(p))
                continue
            elif len(run) > 1:
                multipleRuns.append(os.path.basename(p))
                continue
            MAGprefix = run[0]
            MAGruns.append(MAGprefix)
    
        MAGinfoDict[fileName] = {
            "genome_name" : fileName,
            "accessions" : run,
            "bin" : bin,
            "genome_path" : p,
            "co-assembly" : isCoassembly,
            "MAGprefix" : MAGprefix,
            "binning_software" : binner,
            "binning_parameters" : binParams
        }

        if isCoassembly and len(MAGinfoDict[fileName]["MAGprefix"]) > 14:
            raise ValueError("Co-assemblies names have to be shorter than 15 characters.")

        if p in coassMAGPaths:
            MAGinfoDict[fileName]["assembly_sw"] = coassDict[MAGprefix]["assembler"]

    if missingRuns:
        raise ValueError("O  or more MAG names don't contain any run, nor they appear " + 
            "as a co-assembly generated MAG. Check folder content, co-assembly file " +
            "list, or rename MAGs respecting the format runName_bin.x.fa.\nMAGs " +
            "affected: {}.\nFilenames to edit: {}".format(len(missingRuns), missingRuns))
    if multipleRuns:
        raise ValueError("One or more MAG names contain more than one run. If they were " +
            "generated from a co-assembly, add them in the co-assembly description " +
            "file, and use the '--coassemblies' option. Otherwise, check your MAG " +
            "names and rename them after the format runName_bin.x.fa.\nFilenames " +
            "to edit: {}".format(multipleRuns))
    if not MAGruns and not coassDict:
        raise FileNotFoundError("No MAGs available for upload. Check whether filenames include " +
            "the name of the run they were generated from, and that all files are in .gz format.")

    return MAGinfoDict

def extract_assembler(assemblyAccession):
    manifestXml = minidom.parseString(requests.get("https://www.ebi.ac.uk" +
        "/ena/browser/api/xml/" + assemblyAccession).text)

    program = manifestXml.getElementsByTagName("PROGRAM")
    assembler = program[0].childNodes[0].nodeValue
     
    return assembler

def extract_assembler_from_ena(MAGinfoDict):
    # retrieving original accessions from runs included in MAG names
    allRuns = []
    for MAG in MAGinfoDict:
        allRuns.extend(MAGinfoDict[MAG]["MAGruns"])

    runsSet, studySet, run_to_assembler = set(allRuns), set(), {}
 
    for r in runsSet:
        run_info = handler.get_run(r, fields="secondary_study_accession")
        studySet.add(run_info["secondary_study_accession"])
                    
        if not studySet:
            raise ValueError("No study corresponding to runs found.")

    for s in studySet:
        ena_info = handler.get_study_runs(s)
        if ena_info == []:
            raise IOError("No runs found on ENA for project {}.".format(s))
        for run, item in enumerate(ena_info):
            runAccession = ena_info[run]["run_accession"]
            if runAccession in runsSet:
                try:
                    sampleAccession = ena_info[run]["sample_accession"]
                    assembly = handler.get_assembly_from_sample(sampleAccession)
                    analysisAccession = assembly["analysis_accession"]
                    assembler = extract_assembler(analysisAccession)
                except:
                    assembler = "not provided"

                run_to_assembler[runAccession] = assembler

    return run_to_assembler

def determine_biomes(lineageTab, metagenome):
    lineageTab = [elem for elem in lineageTab.split(',')]

    if len(lineageTab) != 3:
        print("Environment variables must be 3: biome, feature, and material.")
        sys.exit(1)

    m = check_metagenome_int(metagenome)
    if not m:
        if not check_metagenome_string(metagenome):
            print("Unrecognised metagenome")
            sys.exit(1)
        else:
            lineageTab.append(metagenome)

    return lineageTab

def check_metagenome_int(rsp):
    try:
        if abs(int(rsp)) < len(metagenomes):
            return metagenomes[int(rsp)]
        else:
            return False
    except:
        return False

def check_metagenome_string(rsp):
    return rsp in metagenomes

def addStatToDict(statsDict, name, comp, cont):
    completeness = round(float(comp), 2)
    if completeness == 100.0:
        completeness = 100
    contamination = round(float(cont), 2)
    statsDict[name + ".gz"] = [completeness, contamination]

def extract_MAG_stats(MAGNames, MAGinfoDict, stats):
    compContDict = {}
    MAGName, completeness, contamination = "", 0, 0
    uniquePrefixes = set()
    for MAG in MAGinfoDict:
        uniquePrefixes.add(MAGinfoDict[MAG]["MAGprefix"])

    with open(stats, 'r') as f:
        try:
            delim = '\t'
            file = csv.reader(f, delimiter=delim)
            firstLine = file.__next__()
            if len(firstLine) < 3:
                raise IndexError
        except IndexError:
            delim=','

        f.seek(0)
        file = csv.reader(f, delimiter=delim)
        for line in file:
            MAGName = line[0]
            if ".fa" not in MAGName:
                MAGName += ".fa"
            if MAGName in MAGNames:
                completeness = line[1]
                contamination = line[2]
                addStatToDict(compContDict, MAGName, completeness, contamination)
    
    if len(compContDict) != len(MAGinfoDict):
        missingStats = [MAG for MAG in MAGNames if MAG+".gz" not in compContDict]
        raise IOError("Stats missing for MAGs {}.".format(','.join(m for m in missingStats)))

    return compContDict

def compute_RNA_presence(MAGNames, rnaDir, upDir, euks):
    tRNADict, rRNADict = {}, {}
    if not euks:
        # RNA detection
        # This process can take a while, therefore a backup file is generated if anything fails midway through
        tRNAFiles = glob.glob(os.path.join(rnaDir, "*tRNA_20aa.out"))
        rRNAFiles = glob.glob(os.path.join(rnaDir, "*rRNAs.out"))
        backupRrnaFile = os.path.join(upDir, "rRNAbackup.json")
        backupTrnaFile = os.path.join(upDir, "tRNAbackup.json")

        tRNAcounter, rRNAcounter = 0, 0
        # tRNA
        if not os.path.exists(backupTrnaFile):
            with open(backupTrnaFile, 'w') as file:
                pass
        with open(backupTrnaFile, "r+") as file:
            try:
                tRNAbackupDict = json.loads(file.read())
                tRNADict = dict(tRNAbackupDict)
                tqdm.write("\tA tRNA backup file has been found")
            except json.decoder.JSONDecodeError:
                tRNAbackupDict = {}
            for f in tRNAFiles:
                if not os.path.exists(f):
                    raise IOError("\t\t Failing to retrieve some or all tRNA detection files " +
                        "(*tRNA_20aa.out).Please relaunch the procedure described at https://" +
                        "www.ebi.ac.uk/seqdb/confluence/pages/viewpage.action?spaceKey=" +
                        "MetaGen&title=MAGs+analysis+workflow")
                
                # need for numerical indices, as no header exists
                MAGname = os.path.basename(f.replace("_tRNA_20aa.out", ".fa.gz"))
                if MAGname not in tRNAbackupDict:
                    try:
                        tempCSV = pd.read_csv(f, sep='\t', header=None)
                        tRNApresent = True
                        if int(tempCSV.iloc[0][1]) < 18:
                            tRNApresent = False
                    except:
                        tRNApresent = False

                    tRNADict[MAGname] = tRNApresent
                    
                    tRNAcounter += 1
                    if (tRNAcounter%10 == 0) or (len(tRNAFiles) - len(tRNAbackupDict) == tRNAcounter):
                        file.seek(0)
                        file.write(json.dumps(tRNADict))
                        file.truncate()

                tRNADict = {**tRNADict, **tRNAbackupDict}


        # rRNA
        if not os.path.exists(backupRrnaFile):
            with open(backupRrnaFile, 'w') as file:
                pass
        with open(backupRrnaFile, "r+") as file:
            try:
                rRNAbackupDict = json.loads(file.read())
                rRNADict = dict(rRNAbackupDict)
                tqdm.write("\tA rRNA backup file has been found")
            except json.decoder.JSONDecodeError:
                rRNAbackupDict = {}
            for f in rRNAFiles:
                if not os.path.exists(f):
                    raise IOError("Failing to retrieve some or all rRNA detection files " +
                        "(*rRNAs.out).Please relaunch the procedure described at https://" + 
                        "www.ebi.ac.uk/seqdb/confluence/pages/viewpage.action?spaceKey=" +
                        "MetaGen&title=MAGs+analysis+workflow")
                
                # need for numerical indices, as no header exists
                MAGname = os.path.basename(f.replace("_rRNAs.out", ".fa.gz"))
                if MAGname not in rRNAbackupDict:
                    try:
                        tempCSV = pd.read_csv(f, sep='\t', header=None)
                        rRNApresent = True
                        rRNAvalues = tempCSV.iloc[:,2].values
                        if 0.00 in rRNAvalues or not pd.Series(rRNAvalues > 80).all():
                            rRNApresent = False
                    except:
                        rRNApresent = False

                    rRNADict[MAGname] = rRNApresent

                    rRNAcounter += 1
                    if (rRNAcounter%10 == 0) or (len(rRNAFiles) - len(rRNAbackupDict) == rRNAcounter):
                        file.seek(0)
                        file.write(json.dumps(rRNADict))
                        file.truncate()

                rRNADict = {**rRNADict, **rRNAbackupDict}

    RNADict = {}
    for m in MAGNames:
        gzName = m + ".gz"
        
        if gzName not in tRNADict:
            tRNADict[gzName] = False
        if gzName not in rRNADict:
            rRNADict[gzName] = False
        
        if tRNADict[gzName] and rRNADict[gzName]:
            RNADict[gzName] = "yes"
        else:
            RNADict[gzName] = "no"
    
    return RNADict

def extract_MAG_coverage(coverageDir):
    covPaths = glob.glob(os.path.join(coverageDir, "*_coverage", "*MAGcoverage.txt"))
    covDict = {}
    for cov_path in covPaths:
        with open(cov_path, 'r') as c:
            for line in c:
                name = line.split('\t')[0] + ".gz"
                coverage = line.split('\t')[1].strip()
                covDict[name] = coverage
    
    return covDict

def extract_Bacteria_lineage(levels, taxDictNCBI):
    for i in reversed(range(len(levels))):
        name = levels[i].split("__")[1]
        if not name == "":
            if i == 6:
                newDef = name
            elif i == 0:
                newDef = "uncultured {}".format(name)
            elif i == 1 or i == 2 or i == 3 or i == 4:
                newDef = "uncultured {} bacterium".format(name)
            elif i == 5:
                newDef = "uncultured {} sp.".format(name)
            try:
                tax_lineage = taxDictNCBI[newDef][1]
                break
            except KeyError:
                if (i == 0 or i == 1) and name.lower().endswith("bacteria"):
                    try:
                        newDef = "uncultured {}".format(name.lower().replace("bacteria", "bacterium"))
                        tax_lineage = taxDictNCBI[newDef][1]
                        break
                    except KeyError:
                        pass
                elif i == 2:
                    if name.lower() == "Deltaproteobacteria":
                        newDef = "uncultured delta proteobacterium"
                        try:
                            tax_lineage = taxDictNCBI[newDef][1]
                        except KeyError:
                            pass
                else:
                    pass

    return tax_lineage

def extract_Archaea_lineage(levels, taxDictNCBI):
    for i in reversed(range(len(levels))):
        name = levels[i].split("__")[1]
        if not name == "":
            if i == 6:
                newDef = name
            elif i == 0:
                newDef = "uncultured archaeon"
            elif i == 1:
                if "Euryarchaeota" in name:
                    newDef = "uncultured euryarchaeote"
                elif "Candidatus" in name:
                    newDef = "{} archaeon".format(name)
                else:
                    newDef = "uncultured {} archaeon".format(name)
            if i == 2 or i == 3 or i == 4:
                newDef = "uncultured {} archaeon".format(name)
            if i == 5:
                newDef = "uncultured {} sp.".format(name)
            try:
                tax_lineage = taxDictNCBI[newDef][1]
                break
            except KeyError:
                if "Candidatus" in newDef:
                    if i == 1:
                        try:
                            newDef = newDef.replace("Candidatus ", '')
                            tax_lineage = taxDictNCBI[newDef][1]
                            break
                        except KeyError:
                            pass
                    elif i == 4:
                        try:
                            newDef = newDef.replace("uncultured ", '')
                            tax_lineage = taxDictNCBI[newDef][1]
                            break
                        except KeyError:
                            pass

    return tax_lineage

def compute_MAG_taxonomy(taxInfo, euks, taxNCBI, archaea, bac120):
    taxDictMAG = {}

    if euks:
        '''
            Expected input from BAT for eukaryotic taxonomy:
            #bin   classification  reason  lineage lineage scores
            ERR598983_bin.99.fa     classified      based on 11317/13161 ORFs \
            1;131567;2759;33090;3041;1035538;13792;41873;38832;38833;564608* \
            1.00;1.00;1.00;1.00;1.00;1.00;1.00;0.99;0.99;0.94;0.94
        '''
        
        with open(taxInfo, 'r') as eukaryote:
            next(eukaryote)
            for line in eukaryote:
                info = line.split('\t')
                binName = info[0]
                lineage = info[3].split(';')
                taxDictMAG[binName] = lineage
    else:
        '''
            Expected input from GTDB for prokaryotic taxonomy:
            #bin   classification  extraFields...
            ERR598983_bin.99.fa \
            d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Levilactobacillus;s__
        '''
        
        taxDictNCBI = {}
        with open(taxNCBI, 'r') as t:
            for line in t:
                features = line.split('|')
                taxID = features[0].strip()
                taxName = features[1].strip()
                taxLineage = features[2].strip()
                taxDictNCBI[taxName] = [taxID, taxLineage]
        
        gtdbDir = glob.glob(taxInfo)[0]
        tax = gtdb_to_ncbi_majority_vote_v2.Translate()
        taxDictGTDB = tax.run(gtdbDir, archaea, bac120, "gtdbtk")
        tax_NCBI_lineage = ""
        
        for MAGname in taxDictGTDB:
            levels = taxDictGTDB[MAGname].split(';')
            if "Archaea" in taxDictGTDB[MAGname]:
                tax_NCBI_lineage = extract_Archaea_lineage(levels, taxDictNCBI)
            elif "Bacteria" in taxDictGTDB[MAGname]:
                tax_NCBI_lineage = extract_Bacteria_lineage(levels, taxDictNCBI)
            fastaName = str(MAGname) + ".fa.gz"
            taxClassification = ("Taxonomy was originally computed"
                "with GTDBtk, which assigned the following taxonomic "
                "annotation: {}".format(taxDictGTDB[MAGname]))
            taxDictMAG[fastaName] = [tax_NCBI_lineage, taxClassification]

    return taxDictMAG

def create_MAG_dictionary(ena_uploader):
    MAGinfoDict = extract_MAG_info(ena_uploader.magDir, ena_uploader.coassembly, 
                                    ena_uploader.bins, ena_uploader.binningSW, ena_uploader.binningParams)
    assemblerDict = extract_assembler_from_ena(MAGinfoDict)
    biomes = determine_biomes(ena_uploader.envVars, ena_uploader.metagenome)

    simpleNames = [MAG[:-3] for MAG in MAGinfoDict.keys()] # remove .gz extension
    MAGstats = extract_MAG_stats(simpleNames, MAGinfoDict, ena_uploader.stats)
    RNA_presence = compute_RNA_presence(simpleNames, ena_uploader.rnaDir, ena_uploader.MAG_UPLOAD_DIR, ena_uploader.euks)
    taxDict = compute_MAG_taxonomy(ena_uploader.taxInfo, ena_uploader.euks, ena_uploader.taxNCBI) # every key has both NCBI lineage and GTDB lineage
    coverageDepth = extract_MAG_coverage(ena_uploader.covDir)
    
    for MAGname in MAGinfoDict:
        MAGinfoDict[MAGname]["genome_name"] = "GENOMENAMEOFSOMESORT"
        MAGinfoDict[MAGname]["assembly_sw"] = assemblerDict[MAGinfoDict["accessions"]] # unless it's coassembly - then we already defined it in the first function 
        MAGinfoDict[MAGname]["binning_sw"] = "WE'LLFINDAWAY"
        MAGinfoDict[MAGname]["stats_generation_sw"] = ena_uploader.completenessSW
        MAGinfoDict[MAGname]["completeness"] = str(MAGstats[MAGname][0])
        MAGinfoDict[MAGname]["contamination"] = str(MAGstats[MAGname][1])        
        MAGinfoDict[MAGname]["environmentBiome"] = biomes[0]
        MAGinfoDict[MAGname]["environmentFeature"] = biomes[1]
        MAGinfoDict[MAGname]["environmentMaterial"] = biomes[2]
        MAGinfoDict[MAGname]["metagenome"] = biomes[3]
        MAGinfoDict[MAGname]["NCBI_lineage"] = taxDict[MAGname][0]
        MAGinfoDict[MAGname]["genome_coverage"] = coverageDepth[MAGname]
        MAGinfoDict[MAGname]["RNA_presence"] = RNA_presence[MAGname]
        
    return MAGinfoDict

def write_MAG_tsv(MAGs, uploadDir):
    with open(os.path.join(uploadDir, 'genomes_metadata.tsv'), 'wb') as f:
        dict_writer = csv.DictWriter(f, MAGs.keys(), delimiter='\t')
        dict_writer.writeheader()
        dict_writer.writerows(MAGs)

def root():
    ena_uploader = MAGupload()
    
    # tsv creation
    MAGs = create_MAG_dictionary(ena_uploader)
    write_MAG_tsv(MAGs, ena_uploader.MAG_UPLOAD_DIR)

class MAGupload:
    def __init__(self, argv=sys.argv[1:]):
        self.args = parse_args(argv)
        self.database = self.args.database
        self.magDir = self.args.mag_folder.rstrip('\/')
        self.dir = self.magDir.rsplit('/', 1)[0]
        self.coassembly = self.args.coassemblies if self.args.coassemblies else ""
        
        self.MAG_UPLOAD_DIR = self.create_MAG_upload_dir(self.dir, self.bins)
        
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
    
    def create_MAG_upload_dir(self, d, bins):
        upload_dir = os.path.join(d, 'MAG_upload')
        if bins:
            upload_dir = upload_dir.replace("MAG", "bin")
        os.makedirs(upload_dir, exist_ok=True)
        
        return upload_dir

# TODO remove once introduced in mi_automation
if __name__ == "__main__":
    root()
    tqdm.write('Completed')
