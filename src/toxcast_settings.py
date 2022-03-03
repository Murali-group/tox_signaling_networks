# This file contains various settings for running and testing different interactomes, weighting parameters, and so on.

#import toxcast_utils as t_utils

VERSION = ''

INTERACTOMES = {}

# may specify the path to the evidence file used for this version, otherwise the default is used
EVIDENCE_FILES = {}

# zscore penalty versions. Adds a cost to the super-source -> rec and from the tf -> super-target based on the z-score of the assay
ZSCORE_PENALTY = []

# bulk add these versions
for version in ["2018_01-toxcast-d2d-p1_5-u1_25"]:
    INTERACTOMES[version] = "./inputs/2018_01-toxcast-net/2019-02-18-human-ppi-d2d-dir-weighted-cap0_99.txt"
    ZSCORE_PENALTY.append(version)
    EVIDENCE_FILES[version] = "./inputs/2018_01-toxcast-net/2018_01pathlinker-nokegg.tsv"

for version in ["2018_01-toxcast-p1_5"]:
    INTERACTOMES[version] = "./inputs/2018_01-toxcast-net/2019-02-18-human-ppi-weighted-cap0_99.txt"
    ZSCORE_PENALTY.append(version)
    EVIDENCE_FILES[version] = "./inputs/2018_01-toxcast-net/2018_01pathlinker-nokegg.tsv"

ALLOWEDVERSIONS = sorted(INTERACTOMES.keys())

# -log of this number is added to each edge before running cyclinker.
# this effectively penalizes each path by the # of edges in the path * -log(edge_penalty)
EDGE_PENALTIES = {
    "2018_01-toxcast-p1_5": 1.5,
    "2018_01-toxcast-d2d-p1_5-u1_25": 1.5,
}

UNDIR_PENALTIES = {
    "2018_01-toxcast-d2d-p1_5-u1_25": 1.25,
}

## DATADIR is the path to the data/ directory checked into SVN.  
DATADIR = '/data/jeff-law/data/svn-data'
PATHLINKERDATADIR = '/data/jeff-law/projects/2015-03-pathlinker/data/pathway-specific-interactomes'
# populate the interactomes paths
#for version in INTERACTOMES:
#    if "%s" in INTERACTOMES[version]:
#        if version in ["netpath-pathlinker-signaling-children-reg", "kegg-pathlinker-signaling-children-reg"]:
#            INTERACTOMES[version] = INTERACTOMES[version] % PATHLINKERDATADIR
#        else:
#            INTERACTOMES[version] = INTERACTOMES[version] % DATADIR

# these are specified by input options
INPUTSPREFIX = ''
RESULTSPREFIX = ''
REC_TFS_FILE = "%s/rec-tfs/%s-rec-tfs.txt"
# Instead of create a different interactome for each version, I will simply post-process the output of cyclinker
#INTERACTOME_FILES = "%s/interactome/%s-interactome.txt"
# follows the convention inputs/version/version-interactome.txt where version is the version name
SPLIT_REC_TFS_INTERACTOME = '%s/%s-interactome.txt'
CHEMICAL_MAP = ''
#chemDSStoName, chemNametoDSS = t_utils.getChemicalNameMaps()
# this is currently the 'scope' we are using to generate the random networks
# 'permute-dir-undir' randomly swaps edges separately in the directed and undirected graphs
DEFAULT_SCOPE = "permute-dir-undir"


def set_version(version):
    global VERSION, RESULTSPREFIX, INPUTSPREFIX, INTERACTOME
    global EDGE_PENALTY, REC_TFS_PENALTY, SPLIT_REC_TFS_INTERACTOME

    VERSION = version
    print("Using version %s" % (VERSION))

    INPUTSPREFIX = "inputs/versions/%s/" % VERSION
    RESULTSPREFIX = "outputs/%s/weighted/" % VERSION
    INTERACTOME = INTERACTOMES[VERSION]
    # Also create a new 
    if VERSION in UNDIR_PENALTIES:
        INTERACTOME = "%s/%s-interactome-undir-penalty.txt" % (INPUTSPREFIX, VERSION)

    # also setup some other variables for running cyclinker for each version
    if VERSION in ZSCORE_PENALTY:
        REC_TFS_PENALTY = True
    else:
        REC_TFS_PENALTY = False
    if VERSION in EDGE_PENALTIES:
        EDGE_PENALTY = EDGE_PENALTIES[VERSION]
    else:
        EDGE_PENALTY = None

    return INPUTSPREFIX, RESULTSPREFIX, INTERACTOME
