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