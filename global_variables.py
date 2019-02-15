import os

################################################################################
#################### Chromosome and contig names ###############################
################################################################################

# Chromosome and scaffold names for later use
CHR_GRCh38 = ["chromosome."+str(x) for x in range(1,23)] \
           + ["chromosome."+str(x) for x in ["MT","X","Y"]]

EGYPTREF_SCAFFOLDS = ["fragScaff_scaffold_"+str(x)+"_pilon" for x in range(0,41)] \
                + ["original_scaffold_"+str(x)+"_pilon" for x in range(41,145)]

EGYPTREFV2_SCAFFOLDS = ["fragScaff_scaffold_"+str(x)+"_pilon" for x in range(0,226)] \
                + ["original_scaffold_"+str(x)+"_pilon" for x in range(226,1728)]
# Note: The EGYPTREFV2_SCAFFOLDS for which no repeats have been detected are:
# original_scaffold_{1078,499,447,1298,778,956,1014,583,471,1349,303,1632,1186,1643,399,535,1662,1067,1724,1572,1701,985,719,1711,1101,318,1122,731}_pilon

CEGYPT_CONTIGS = ["Contig"+str(x) for x in range(0,360)]

# Read in contig names from pre-generated file which reads in the fasta headers
CEGYPTV2_CONTIGS = []
with open("data/file.contigsetv2.seqnames.txt","r") as f_in:
    for line in f_in:
        # Remove ">" at start and "|arrow" at end
        CEGYPTV2_CONTIGS.append(line.split("|")[0][1:])

CHR_YORUBA = [x for x in CHR_GRCh38 if not x in ["chromosome.MT","chromosome.Y"]]

YORUBA_SCAFFOLDS = []
if os.path.exists("seq_YORUBA/yoruba_scaffold_to_genbank.txt"):
    with open("seq_YORUBA/yoruba_scaffold_to_genbank.txt") as f_in:
        for line in f_in:
            s = line.split("\t")
            if not "HS_" in line:
                YORUBA_SCAFFOLDS.append("chromosome."+s[0])
            else:
                YORUBA_SCAFFOLDS.append(s[0])

AK1_SCAFFOLDS = []
if os.path.exists("seq_AK1/ak1_scaffold_to_genbank.txt"):
    with open("seq_AK1/ak1_scaffold_to_genbank.txt") as f_in:
        for line in f_in:
            AK1_SCAFFOLDS.append(line.split("\t")[0])

# For plotting etc. we sometimes want the longest SCAFFOLDS, since these are not
# named according to size, here are the longest ones
LONGEST_AK1_SCAFFOLDS = [
    "Scaffold0147","Scaffold0001","Scaffold00019","Scaffold0008", \
    "Scaffold0151","Scaffold0148","Scaffold00033","Scaffold0002", \
    "Scaffold00022","Scaffold0152","Scaffold0150","Scaffold00034", \
    "Scaffold00063","Scaffold00068","Scaffold0007","Scaffold00025", \
    "Scaffold00066","Scaffold00032","Scaffold0010","Scaffold0011", \
    "Scaffold00030_pilon","Scaffold00012","Scaffold0009","Scaffold00011", \
    "Scaffold00067","Scaffold00027","Scaffold0142","Scaffold0012", \
    "Scaffold0056","Scaffold0013"
]

LONGEST_EGYPTREFV2_SCAFFOLDS = ["fragScaff_scaffold_"+str(x)+"_pilon" for x in \
    [100,170,6,123,149,184,89,195,205,163,201,76,155,29,68,137,80,61,154,147, \
     116,212,196,158,9,26,186,194,98]] + ["original_scaffold_1041_pilon"]


################################################################################
####################### PacBio-related variables ###############################
################################################################################

# There are 5 PacBio libraries from the same individual, each sequences in 
# various sequencing runs  
PACBIO_SAMPLES = ["r54171","r54172","r54212","r54214","r54217"]

# The naming convention for folders is the sample name, _, then the seqrun ID, 
# The naming convention for files is the same, but for some reason the "r" of
# the samples has been replaced by "m"; also sum files are in subdirectories
# Since there seems no apparent system to the file naming, I here just map the
# samples to the corresponding Pacbio filenames (without ending, but the file
# basename is always the same)
PACBIO_SAMPLES_TO_SEQRUN_PATH = { \
    "r54171": ["r54171_180507_074037/m54171_180507_074037", \
               "r54171_180508_081816/m54171_180508_081816", \
               "r54171_180509_085337/m54171_180509_085337", \
               "r54171_180509_190202/m54171_180509_190202", \
               "r54171_180510_051157/m54171_180510_051157", \
               "r54171_180511_073925/m54171_180511_073925", \
               "r54171_180511_174954/m54171_180511_174954", \
               "r54171_180512_040316/m54171_180512_040316", \
               "r54171_180512_141733/m54171_180512_141733", \
               "r54171_180513_003153/m54171_180513_003153", \
               "r54171_180514_191117/m54171_180514_191117", \
               "r54171_180515_052445/m54171_180515_052445", \
               "r54171_180515_153940/m54171_180515_153940"],\
    "r54172": ["r54172_20180226_063627/1_A08/m54172_180226_064443", \
               "r54172_20180227_060945/1_A08/m54172_180227_061743", \
               "r54172_20180227_060945/2_B08/m54172_180227_162339", \
               "r54172_20180227_060945/3_C08/m54172_180228_023312", \
               "r54172_20180301_065149/2_B08/m54172_180301_170719"], \
    "r54212": ["r54212_20180207_084734/1_A05/m54212_180207_085743"], \
    "r54214": ["r54214_20180225_094705/1_A08/m54214_180225_095639", \
               "r54214_20180226_063218/1_A08/m54214_180226_064236", \
               "r54214_20180226_063218/2_B08/m54214_180226_164754", \
               "r54214_20180227_074241/1_A08/m54214_180227_075436", \
               "r54214_20180227_074241/2_B08/m54214_180227_180004", \
               "r54214_20180228_083736/1_A05/m54214_180228_084706", \
               "r54214_20180301_092943/1_A08/m54214_180301_094052", \
               "r54214_20180301_092943/2_B08/m54214_180301_194631", \
               "r54214_20180301_092943/3_C08/m54214_180302_055606", \
               "r54214_20180303_091311/1_A08/m54214_180303_092301", \
               "r54214_20180304_073054/1_A05/m54214_180304_074025", \
               "r54214_20180304_073054/2_B05/m54214_180304_174558", \
               "r54214_20180304_073054/3_C05/m54214_180305_035534", \
               "r54214_20180304_073054/4_D05/m54214_180305_140511", \
               "r54214_20180304_073054/5_E05/m54214_180306_001437", \
               "r54214_20180304_073054/6_F05/m54214_180306_102433", \
               "r54214_20180304_073054/7_G05/m54214_180306_203421", \
               "r54214_20180304_073054/8_H05/m54214_180307_064357", \
               "r54214_20180308_072240/1_A01/m54214_180308_073253", \
               "r54214_20180308_072240/2_B01/m54214_180308_173821", \
               "r54214_20180309_085608/1_A01/m54214_180309_090535", \
               "r54214_20180309_085608/2_B01/m54214_180309_191107", \
               "r54214_20180309_085608/3_C01/m54214_180310_052041", \
               "r54214_20180309_085608/4_D01/m54214_180310_153039", \
               "r54214_20180309_085608/5_E01/m54214_180311_014012", \
               "r54214_20180309_085608/6_F01/m54214_180311_114949", \
               "r54214_20180312_065341/1_A08/m54214_180312_071349", \
               "r54214_20180313_083026/1_A08/m54214_180313_083936", \
               "r54214_20180314_082924/1_A05/m54214_180314_083852"], \
    "r54217": ["r54217_20180205_093834/1_A01/m54217_180205_095019"]
}


################################################################################
####################### 10X-related variables ##################################
################################################################################

ILLUMINA_10X_LIBS = [
    "NDHX00201-AK654_L4",
    "NDHX00201-AK654_L5",
    "NDHX00201-AK654_L6",
    "NDHX00201-AK654_L7",
    "NDHX00201-AK655_L4",
    "NDHX00201-AK655_L5",
    "NDHX00201-AK655_L6",
    "NDHX00201-AK655_L7",
    "NDHX00201-AK656_L4",
    "NDHX00201-AK656_L5",
    "NDHX00201-AK656_L6",
    "NDHX00201-AK656_L7",
    "NDHX00201-AK657_L4",
    "NDHX00201-AK657_L5",
    "NDHX00201-AK657_L6",
    "NDHX00201-AK657_L7"
]