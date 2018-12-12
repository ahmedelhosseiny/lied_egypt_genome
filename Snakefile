#kate:syntax python;

#######################################
### Analyzing an Egyptian genome
#######################################

from Bio import SeqIO
import os


# Rules to be executed for new assemblies: compute_content_and_assembly_numbers,
# repeatmasker_summary_table_egyptrefv2, align_assemblies_with_mummer_all, 

################################################################################
############### Writing some general statistics to file ########################
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
######### Preprocessing of the first Novogene assembly called EGYPTREF #########
################################################################################

# Writing the scaffolds of the Egyptian genome to separate fasta files because
# processing the whole assembly often takes too much time
rule write_scaffold_fastas:
    input: "data/pilon.fasta"
    output: expand("seq_EGYPTREF/Homo_sapiens.EGYPTREF.dna.{scaffold}.fa", \
                   scaffold=EGYPTREF_SCAFFOLDS)
    run:
        with open(input[0], "r") as f_in:
            i = 0
            for record in SeqIO.parse(f_in,"fasta"):            
                with open(output[i], "w") as f_out:
                    SeqIO.write(record, f_out, "fasta")
                    i += 1

# Writing the contigs of the Egyptian genome to separate fasta files because
# processing the whole assembly often takes too much time
rule write_contig_fastas:
    input: "data/file.contigset.fasta"
    output: expand("seq_CEGYPTREF/Homo_sapiens.CEGYPTREF.dna.{contig}.fa", \
                   contig=CEGYPT_CONTIGS)
    run:
        with open(input[0], "r") as f_in:
            for record in SeqIO.parse(f_in,"fasta"):            
                # Remove the trailing "|arrow" since pipe symbols can cause
                # problems and the "|arrow" is not needed
                record.id = record.id[:-6]
                record.name = ''
                record.description = ''
                record.seq = record.seq.upper()
                out_fname = "seq_CEGYPTREF/Homo_sapiens.CEGYPTREF.dna." + \
                               record.id+".fa"
                with open(out_fname, "w") as f_out:
                    SeqIO.write(record, f_out, "fasta")

# Making a file with all contigs; this is the same as the 
# data/file.contigset.fastafile provided by Novogene, but the "|arrow" in the
# sequence names removed
# Note: The original contig file has upper and lower-case letters, don't know
# why! Perhaps repeatmasking was done with them? Anyway, I convert the sequences
# to upper case, because repeatmasking is done later and should be as for the 
# other assemblies
rule combine_contigs_to_primary_assembly:
    input: expand("seq_CEGYPTREF/Homo_sapiens.CEGYPTREF.dna.{contig}.fa", \
                   contig=CEGYPT_CONTIGS)
    output: "seq_CEGYPTREF/Homo_sapiens.CEGYPTREF.dna.primary_assembly.fa"
    shell: "cat {input} > {output}"

# Copy the assembled sequence
rule cp_and_rename_assembly:
    input: "data/pilon.fasta"
    output: "seq_EGYPTREF/Homo_sapiens.EGYPTREF.dna.primary_assembly.fa"
    shell: "cp {input} {output}"


################################################################################
####### Preprocessing of the second Novogene assembly called EGYPTREFV2 ########
################################################################################

# Writing the scaffolds of the Egyptian genome to separate fasta files because
# processing the whole assembly often takes too much time
rule write_scaffold_fastas_v2:
    input: "data/pilon_v2.fasta"
    output: expand("seq_EGYPTREFV2/Homo_sapiens.EGYPTREFV2.dna.{scaffold}.fa", \
                   scaffold=EGYPTREFV2_SCAFFOLDS)
    run:
        with open(input[0], "r") as f_in:
            i = 0
            for record in SeqIO.parse(f_in,"fasta"):            
                with open(output[i], "w") as f_out:
                    SeqIO.write(record, f_out, "fasta")
                    i += 1

# Copy the assembled sequence
rule cp_and_rename_assembly_v2:
    input: "data/pilon_v2.fasta"
    output: "seq_EGYPTREFV2/Homo_sapiens.EGYPTREFV2.dna.primary_assembly.fa"
    shell: "cp {input} {output}"


################################################################################
############### Computing some statistics for all assemblies ###################
################################################################################

# Just getting the header lines of the individual sequences in the fasta
rule scaffold_names:
    input: "seq_{assembly}/{fname}.fa"
    output: "results/{assembly}/scaffold_names_{fname}.txt"
    shell: "cat {input} | grep '>' > {output}"
    
# Quantifying the sequence content individually for all scaffolds
rule sequence_content:
    input: "seq_{assembly}/{fname}.fa"
    output: "results/{assembly}/num_bases_{fname}.txt"
    script: "scripts/sequence_content.py"

# Quantifying the sequence content over all scaffolds
rule sequence_content_overall:
    input: "seq_{assembly}/{fname}.fa"
    output: "results/{assembly}/num_all_{fname}.txt"
    script: "scripts/sequence_content_overall.py"
            
# Compute N50 and other related values as statistic for the assembly
rule compute_assembly_stats:
    input: "results/{assembly}/num_bases_{fname}.txt"
    output: "results/{assembly}/assembly_stats_{fname}.txt"
    script: "scripts/compute_assembly_stats.py"

# Computing all info numbers:
rule compute_content_and_assembly_numbers:
    input: expand( \
           "results/{assembly}/{task}_Homo_sapiens.{assembly}.dna.primary_assembly.txt", \
           assembly = ["GRCh38","EGYPTREF","AK1","YORUBA","CEGYPTREF","EGYPTREFV2"], \
           task = ["scaffold_names","num_bases","num_all","assembly_stats"])


################################################################################
############### Finding mammalian core genes as QC using busco #################
################################################################################

# Downloading the Busco lineage information
rule download_linage:
    output: temp("busco_lineage/mammalia_odb9.tar.gz")
    shell: "wget -P busco_lineage https://busco.ezlab.org/datasets/mammalia_odb9.tar.gz"
    
# ... and extracting it; the output files are just two of the many files in this
# archive
rule extract_lineage:
    input: "busco_lineage/mammalia_odb9.tar.gz"
    output: "busco_lineage/mammalia_odb9/lengths_cutoff",
            "busco_lineage/mammalia_odb9/scores_cutoff"
    shell: "tar --directory busco_lineage -xvzf {input}"

# Running Busco on a genome file
# --force: Deleting results folder; start new run
# --tmp: Likely /tmp is too small, so make a new tmp folder on scratch (also 
#  this can be accessed much quicker)
# --blast_single_core: There is a (known!) bug, that blast sometimes fails in
# multi-cpu mode. I also observe this for GRCh38, with exactly the corresponding
# error message; therefore, this is run with a single core.
# Note: According to Busco documentation, 3.1Gbp genome assessment with 12 CPUs 
# takes 6 days and 15 hours
# I use a separate environment for busco, because, as of now, its newest version
# cannot be used together with the repeatmasker and installing it together 
# would result in downgrading of augustus, blast, boost and busco to older 
# versions.
rule run_busco:
    input: "busco_lineage/mammalia_odb9/lengths_cutoff",
           "seq_{assembly}/Homo_sapiens.{assembly}.dna.{chr_or_type}.fa"
    output: "busco_{assembly}/run_busco_{assembly}_{chr_or_type}/short_summary_busco_{assembly}_{chr_or_type}.txt",
            "busco_{assembly}/run_busco_{assembly}_{chr_or_type}/full_table_busco_{assembly}_{chr_or_type}.tsv",
    threads: 12
    conda: "envs/busco.yaml"
    shell:  "workdir=$PWD; cd /scratch; " + \
            "rm -rf /scratch/run_busco_{wildcards.assembly}_{wildcards.chr_or_type}; " + \
            "rm -rf /scratch/tmp_busco_{wildcards.assembly}_{wildcards.chr_or_type}; " + \
            "mkdir /scratch/tmp_busco_{wildcards.assembly}_{wildcards.chr_or_type}; " + \
            "cd /scratch; " + \
            "run_busco --in $workdir/{input[1]} " + \
            "--out busco_{wildcards.assembly}_{wildcards.chr_or_type} " + \
            "--lineage_path $workdir/busco_lineage/mammalia_odb9 " + \
            "--mode genome " + \
            "--force " + \
            "--cpu 12 " + \
#           "--blast_single_core " + \
            "--tmp /scratch/tmp_busco_{wildcards.assembly}_{wildcards.chr_or_type}; " + \
            "rm -rf /scratch/tmp_busco_{wildcards.assembly}_{wildcards.chr_or_type}; " + \
            "mkdir -p busco_{wildcards.assembly}; " + \
            "cd $workdir; "
            "rm -rf busco_{wildcards.assembly}/run_busco_{wildcards.assembly}_{wildcards.chr_or_type}; " + \
            "rsync -avz /scratch/run_busco_{wildcards.assembly}_{wildcards.chr_or_type} busco_{wildcards.assembly}/; " + \
            "rm -rf /scratch/run_busco_{wildcards.assembly}_{wildcards.chr_or_type}; "

# Running busco on the entire primary assembly...
rule run_busco_primary_assembly:
    input: "busco_EGYPTREF/run_busco_EGYPTREF_primary_assembly/short_summary_busco_EGYPTREF_primary_assembly.txt",
           "busco_GRCh38/run_busco_GRCh38_primary_assembly/short_summary_busco_GRCh38_primary_assembly.txt"

# ... and running busco chromosome or scaffold-wise
rule run_busco_chromosomewise:
    input: expand("busco_EGYPTREF/run_busco_EGYPTREF_{scaffolds}/short_summary_busco_EGYPTREF_{scaffolds}.txt", \
                  scaffolds=EGYPTREF_SCAFFOLDS),
           expand("busco_GRCh38/run_busco_GRCh38_{chrom}/short_summary_busco_GRCh38_{chrom}.txt", \
                  chrom=CHR_GRCh38)

# Make a comparison table for the busco analysis for EGYPTREF
rule summary_busco_egyptref:
    input: expand("busco_EGYPTREF/run_busco_EGYPTREF_{scaffolds}/full_table_busco_EGYPTREF_{scaffolds}.tsv", \
                  scaffolds=EGYPTREF_SCAFFOLDS)
    output: "busco_EGYPTREF/busco_summary.txt"
    script: "scripts/busco_summary.py"

# Make a comparison table for the busco analysis for GRCh38
rule summary_busco_grch38:
    input: expand("busco_GRCh38/run_busco_GRCh38_{chrom}/full_table_busco_GRCh38_{chrom}.tsv", \
                  chrom=CHR_GRCh38)
    output: "busco_GRCh38/busco_summary.txt"
    script: "scripts/busco_summary.py"


################################################################################
###################### Getting reference sequences #############################
################################################################################

# Downloading all GRCh38 sequence data available from Ensembl (release 93,
# but note, that on sequence level, the release shouldn't make a difference)
rule download_GRCh38:
    output: "seq_GRCh38/Homo_sapiens.GRCh38.{dna_type}.{chr_or_type}.fa.gz"
    run: 
        # Remove target dir to obtain file name for download
        base = output[0].split("/")[1]
        shell("wget -P seq_GRCh38 " + \
        "ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/{base}")
              
# Download README
rule download_GRCh38_readme:
    output: "seq_GRCh38/README"
    shell: "wget -P seq_GRCh38 " + \
           "ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/README"

# Downloading all GRCh38 sequence files available under the ENSEMBLE release 93
# FTP address 
CHR_OR_TYPE = ["chromosome."+str(x) for x in range(1,23)] \
       + ["chromosome."+str(x) for x in ["MT","X","Y"]] \
       + ["nonchromosomal","primary_assembly","toplevel","alt"]
rule download_GRCh38_all:
    input: expand("seq_GRCh38/"+ \
                  "Homo_sapiens.GRCh38.{dna_type}.{chr_or_type}.fa.gz", \
            dna_type=["dna","dna_rm","dna_sm"],chr_or_type=CHR_OR_TYPE),
           "seq_GRCh38/README"

# Uncompressing fasta files, needed e.g. for Busco analysis
# -d decompress; -k keep archive; -c to stdout
rule uncompress_fasta:
    input: "seq_GRCh38/{fname}.fa.gz"
    output: "seq_GRCh38/{fname}.fa"
    resources: io=1
    shell: "gzip -cdk {input} > {output}"


################################################################################
#### Getting the genome assembly (only chromosomes) of a Yoruba individual #####
################################################################################

# Getting the chromosome sequences from Genbank for 1000G individual NA19240,
# which is a Yoruba female
# Getting the assembly report
rule get_yoruba_assembly_report:
    output: "seq_YORUBA/GCA_001524155.4_NA19240_prelim_3.0_assembly_report.txt"
    shell: "wget -P seq_YORUBA ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/524/155/GCA_001524155.4_NA19240_prelim_3.0/GCA_001524155.4_NA19240_prelim_3.0_assembly_report.txt"

# Getting the scaffold names and the corresponding Genbank_ids from the report
rule get_yoruba_scaffold_names:
    input: "seq_YORUBA/GCA_001524155.4_NA19240_prelim_3.0_assembly_report.txt"
    output: "seq_YORUBA/yoruba_scaffold_to_genbank.txt"
    shell: "cat {input} | grep -v '#' | cut -f 1,5 > {output}"

YORUBA_SCAFFOLD_TO_GENBANK = {}
if os.path.exists("seq_YORUBA/yoruba_scaffold_to_genbank.txt"):
    with open("seq_YORUBA/yoruba_scaffold_to_genbank.txt") as f_in:
        for line in f_in:
            s = line.strip("\n").split("\t")
            if not "HS_" in line:
                YORUBA_SCAFFOLD_TO_GENBANK["chromosome."+s[0]] = s[1]
            else:
                YORUBA_SCAFFOLD_TO_GENBANK[s[0]] = s[1]

rule get_yoruba:
    output: "seq_YORUBA/Homo_sapiens.YORUBA.dna.{id}.fa"
    params: genbank_id=lambda wildcards: YORUBA_SCAFFOLD_TO_GENBANK[wildcards.id]
    script: "scripts/get_genbank_seqs.py"

rule get_yoruba_all:
    input: expand("seq_YORUBA/Homo_sapiens.YORUBA.dna.{scaffold}.fa", \
                  scaffold=YORUBA_SCAFFOLDS)

rule yoruba_primary_assembly:
    input: expand("seq_YORUBA/Homo_sapiens.YORUBA.dna.{scaffolds}.fa", \
                  scaffolds=YORUBA_SCAFFOLDS)
    output: "seq_YORUBA/Homo_sapiens.YORUBA.dna.primary_assembly.fa"
    shell: "cat {input} > {output}"


################################################################################
#### Getting the genome assembly (all scaffolds) of a Korean individual ########
################################################################################

# Getting the assembly report
rule get ak1_assembly_report:
    output: "seq_AK1/GCA_001750385.2_AK1_v2_assembly_report.txt"
    shell: "wget -P seq_AK1 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/750/385/GCA_001750385.2_AK1_v2/GCA_001750385.2_AK1_v2_assembly_report.txt"

# Getting the scaffold names and the corresponding Genbank_ids from the report
rule get_ak1_scaffold_names:
    input: "seq_AK1/GCA_001750385.2_AK1_v2_assembly_report.txt"
    output: "seq_AK1/ak1_scaffold_to_genbank.txt"
    shell: "cat {input} | grep -v '#' | cut -f 1,5 > {output}"

AK1_SCAFFOLD_TO_GENBANK = {}
if os.path.exists("seq_AK1/ak1_scaffold_to_genbank.txt"):
    with open("seq_AK1/ak1_scaffold_to_genbank.txt") as f_in:
        for line in f_in:
            AK1_SCAFFOLD_TO_GENBANK[line.split("\t")[0]]=line.strip().split("\t")[-1]

rule get_ak1:
    input: "seq_AK1/ak1_scaffold_to_genbank.txt"
    output: "seq_AK1/Homo_sapiens.AK1.dna.{scaffold}.fa"
    params: genbank_id=lambda wildcards: AK1_SCAFFOLD_TO_GENBANK[wildcards.scaffold]
    script: "scripts/get_genbank_seqs.py"

rule get_ak1_all:
    input: expand("seq_AK1/Homo_sapiens.AK1.dna.{scaffold}.fa", \
                  scaffold=AK1_SCAFFOLDS)

# Construct one file with all AK1 sequences (called "primary_assembly")
rule ak1_primary_assembly:
    input: expand("seq_AK1/Homo_sapiens.AK1.dna.{scaffold}.fa", \
                  scaffold=AK1_SCAFFOLDS)
    output: "seq_AK1/Homo_sapiens.AK1.dna.primary_assembly.fa"
    shell: "cat {input} > {output}"


################################################################################
######################### Repeat masking with repeatmasker #####################
################################################################################

# Running repeatmasker on the Egyptian genome assembly
# I use a separate environment for repeatmasker, because, as of now, it cannot 
# be used together with the newest busco version and installing it together 
# would result in downgrading of augustus, blast, boost and busco to older 
# versions.
# -s  Slow search; 0-5% more sensitive, 2-3 times slower than default
# -q  Quick search; 5-10% less sensitive, 2-5 times faster than default
# -qq Rush job; about 10% less sensitive, 4->10 times faster than default
# -html Creates an additional output file in xhtml format
# -gff Creates an additional Gene Feature Finding format output
# Note: Result file 
# "repeatmasked_{assembly}/Homo_sapiens.{assembly}.dna.{chr_or_type}.fa.cat.gz"
# is not in the output file list, because depending on the size, either this
# file or the uncompressed file
# "repeatmasked_{assembly}/Homo_sapiens.{assembly}.dna.{chr_or_type}.fa.cat"
# will be generated.
# Temporarily outcommented output files (in case they will also be zipped):
# "repeatmasked_{assembly}/Homo_sapiens.{assembly}.dna.{chr_or_type}.fa.out",
# "repeatmasked_{assembly}/Homo_sapiens.{assembly}.dna.{chr_or_type}.fa.out.gff",
# "repeatmasked_{assembly}/Homo_sapiens.{assembly}.dna.{chr_or_type}.fa.out.html",
rule run_repeatmasker:
    input: "seq_{assembly}/Homo_sapiens.{assembly}.dna.{chr_or_type}.fa"
    output: "repeatmasked_{assembly}/Homo_sapiens.{assembly}.dna.{chr_or_type}.fa.masked",
            "repeatmasked_{assembly}/Homo_sapiens.{assembly}.dna.{chr_or_type}.fa.tbl"
    threads: 12
    conda: "envs/repeatmasker.yaml"
    shell: "workdir=$PWD; cd /scratch; " + \
           "rm -rf /scratch/repeatmasked_{wildcards.assembly}_{wildcards.chr_or_type}; " + \
           "mkdir -p /scratch/repeatmasked_{wildcards.assembly}_{wildcards.chr_or_type}; " + \
           "RepeatMasker -species human " + \
           "             -dir /scratch/repeatmasked_{wildcards.assembly}_{wildcards.chr_or_type} " + \
           "             -pa 12 " + \
           "             -xsmall " + \
           "             -q " + \
           "             -html " + \
           "             -gff $workdir/{input}; " + \
           "cd $workdir; "
           "rsync -avz /scratch/repeatmasked_{wildcards.assembly}_{wildcards.chr_or_type}/ repeatmasked_{wildcards.assembly}/; " + \
           "rm -rf /scratch/repeatmasked_{wildcards.assembly}_{wildcards.chr_or_type}; "

# Running repeatmasker on the primary assembly ...
rule run_repeatmasker_primary_assembly:
    input: "repeatmasked_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.tbl",
           "repeatmasked_EGYPTREF/Homo_sapiens.EGYPTREF.dna.primary_assembly.fa.tbl"

# ... and on the individual scaffolds
# Note: For AK1 and YORUBA, there are the following scaffolds, for which 
# repeatmasker does not generate the output files, because "No repetitive 
# sequences were detected". For these, I generated the fa.masked file manually
# by copying the input file; further, I manually 
# AK1: Scaffold1437, Scaffold0873, Scaffold2568, Scaffold2697, Scaffold2245
# Yoruba: HS_LKPB_CHRUN_SCAFFOLD_{256,672,296,980,165,472,1137,1147,1547,1411,
#         1355,1638,1014}
rule run_repeatmasker_chromosomewise:
    input: expand("repeatmasked_GRCh38/Homo_sapiens.GRCh38.dna.{x}.fa.tbl", \
                  x=CHR_GRCh38),
           expand("repeatmasked_EGYPTREF/Homo_sapiens.EGYPTREF.dna.{x}.fa.tbl", \
                  x=EGYPTREF_SCAFFOLDS),
           expand("repeatmasked_AK1/Homo_sapiens.AK1.dna.{x}.fa.tbl", \
                  x=AK1_SCAFFOLDS),
           expand("repeatmasked_YORUBA/Homo_sapiens.YORUBA.dna.{x}.fa.tbl", \
                  x=YORUBA_SCAFFOLDS)

# Summarising the chromosome-wise repeatmasker summary files for Egyptref
rule repeatmasker_summary_table_egyptref:
    input: expand("repeatmasked_EGYPTREF/Homo_sapiens.EGYPTREF.dna.{x}.fa.tbl", \
                  x=EGYPTREF_SCAFFOLDS)
    output: "repeatmasked_EGYPTREF/summary.txt"
    script: "scripts/repeatmasker_summary.py"

# Summarising the contig-wise repeatmasker summary files for Egyptref
rule repeatmasker_summary_table_cegyptref:
    input: expand("repeatmasked_CEGYPTREF/Homo_sapiens.CEGYPTREF.dna.{x}.fa.tbl", \
                  x=CEGYPT_CONTIGS)
    output: "repeatmasked_CEGYPTREF/summary.txt"
    script: "scripts/repeatmasker_summary.py"

# Summarising the scaffold-wise repeatmasker summary files for Egyptrefv2
rule repeatmasker_summary_table_egyptrefv2:
    input: expand("repeatmasked_EGYPTREFV2/Homo_sapiens.EGYPTREFV2.dna.{x}.fa.tbl", \
                  x=EGYPTREFV2_SCAFFOLDS)
    output: "repeatmasked_EGYPTREFV2/summary.txt"
    script: "scripts/repeatmasker_summary.py"

# Summarising the chromosome-wise repeatmasker summary files for GRCh38
rule repeatmasker_summary_table_grch38:
    input: expand("repeatmasked_GRCh38/Homo_sapiens.GRCh38.dna.{x}.fa.tbl", \
                  x=CHR_GRCh38)
    output: "repeatmasked_GRCh38/summary.txt"
    script: "scripts/repeatmasker_summary.py"
    
# Summarising the chromosome-wise repeatmasker summary files for AK1
rule repeatmasker_summary_table_ak1:
    input: expand("repeatmasked_AK1/Homo_sapiens.AK1.dna.{x}.fa.tbl", \
                  x=AK1_SCAFFOLDS)
    output: "repeatmasked_AK1/summary.txt"
    script: "scripts/repeatmasker_summary.py"

# Summarising the chromosome-wise repeatmasker summary files for Yoruba
rule repeatmasker_summary_table_yoruba:
    input: expand("repeatmasked_YORUBA/Homo_sapiens.YORUBA.dna.{x}.fa.tbl", \
                  x=YORUBA_SCAFFOLDS)
    output: "repeatmasked_YORUBA/summary.txt"
    script: "scripts/repeatmasker_summary.py"

# Making a repeatmasker stat table over all chromosomes, one line for EGYPTREF,
# one line for GRCh38
rule comparison_repeatmasker:
    input: expand("repeatmasked_{assembly}/summary.txt", \
                  assembly=["EGYPTREF","EGYPTREFV2","CEGYPTREF","AK1","YORUBA","GRCh38"])
    output: "results/repeatmasker_comparison.txt"
    script: "scripts/repeatmasker_comparison.py"


################################################################################
########### Reference to assembly genome alignment with lastz ##################
################################################################################

# Computing genome alignments using lastz
# [unmask] Attaching this to the chromosome filename instructs lastz to ignore 
# masking information and treat repeats the same as any other part of the 
# chromosome -> We do NOT want this, alignments will be crappy with it!!
# Parameters used for quick and dirty, alignment (lastz manual), taking minutes
# --notransition Don't allow any match positions in seeds to be satisified by 
#                transitions (lowers seeding sensitivity and reduces runtime)
# --nogapped Eliminates the computation of gapped alignments
# --step 20 Lowers seeding senitivity reducing runtime and memory (factor 3.3)
# Parameters from the Korean reference genome AK1 (Seo et al. 2016)
# --gapped Perform gapped extension of HSPs after first reducing them to anchor 
#          points
# --gap=600,150 Gap open and gap extension penalty
# --hspthresh=4500 Set the score threshold for the x-drop extension method; HSPs
#                  scoring lower are discarded.
# --seed 12of19 Seeds require a 19bp word with matches in 12 specific positions
# --notransition Don't allow any match positions in seeds to be satisified by 
#                transitions
# --ydrop=15000 Set the threshold for terminating gapped extension; this
#               restricts the endpoints of each local alignment by 
#               limiting the local region around each anchor in which 
#               extension is performed
# --chain Perform chaining of HSPs with no penalties
# Parameters from another Korean reference genome, KOREF (Cho et al. 2016)
# --step 19 Offset between the starting positins of successive target words 
#           considered for potential seeds
# --hspthresh 3000 Set the score threshold for the x-drop extension method; HSPs
#                  scoring lower are discared.
# --gappedthresh 3000 Set the threshold for gapped extension; alignments scoring
#                     lower than score are discarded.
# --seed 12of19 Seeds require a 19bp word with matches in 12 specific positions
# --minScore 3000 ? kenttools?
# --linearGap medium ? kenttools?
rule align_with_lastz:
    input: "repeatmasked_GRCh38/Homo_sapiens.GRCh38.dna.{chr}.fa.masked",
           "repeatmasked_{assembly}/Homo_sapiens.{assembly}.dna.{scaffold}.fa.masked"
    output: "align_lastz_GRCh38_vs_{assembly}/{chr}_vs_{scaffold}.maf",
            "align_lastz_GRCh38_vs_{assembly}/dotplots/{chr}_vs_{scaffold}.rdotplot"
    conda: "envs/lastz.yaml"
    shell: "lastz {input[0]} {input[1]} " + \
                                  "--gapped " + \
                                  "--gap=600,150 " + \
                                  "--hspthresh=4500 " + \
                                  "--seed=12of19 " + \
                                  "--notransition " + \
                                  "--ydrop=15000 " + \
                                  "--chain " + \
                                  "--format=maf " + \
                                  "--rdotplot={output[1]} " + \
                                  ">{output[0]}"

# Plot the dotplot output of lastz
#rule individual_lastz_dotplot:
#    input: "align_lastz_GRCh38_vs_EGYPTREF/dotplots/{chr}_vs_{scaffold}.rdotplot"
#    output: "align_lastz_GRCh38_vs_EGYPTREF/dotplots/{chr}_vs_{scaffold}.pdf"
#    script: "scripts/dotplot.R"

# Plotting for one scaffold the dotplot versus all chromosomes
rule dotplots_scaffold_vs_chromosomes:
    input: "results/GRCh38/num_bases_Homo_sapiens.GRCh38.dna.primary_assembly.txt",
           "results/{assembly}/num_bases_Homo_sapiens.{assembly}.dna.primary_assembly.txt",
           expand("align_lastz_GRCh38_vs_{{assembly}}/dotplots/{chr}_vs_{{scaffold}}.rdotplot", \
                  chr=CHR_GRCh38)
    output: "align_lastz_GRCh38_vs_{assembly}/dotplots/{scaffold}.pdf"
    script: "scripts/scaffold_vs_grch38.R"            

# Plotting the dotplots for all scaffolds
rule dotplots_scaffold_vs_chromosomes_all:
    input: expand("align_lastz_GRCh38_vs_EGYPTREF/dotplots/{scaffold}.pdf", \
                  scaffold=EGYPTREF_SCAFFOLDS),
#            expand("align_lastz_GRCh38_vs_YORUBA/dotplots/{scaffold}.pdf", \
#                  scaffold=YORUBA_SCAFFOLDS[:23]),
            expand("align_lastz_GRCh38_vs_CEGYPTREF/dotplots/{contig}.pdf", \
                  contig=CEGYPT_CONTIGS)#,
#            expand("align_lastz_GRCh38_vs_AK1/dotplots/{scaffold}.pdf", \
#                   scaffold=LONGEST_AK1_SCAFFOLDS),
#            expand("align_lastz_GRCh38_vs_EGYPTREFV2/dotplots/{scaffold}.pdf", \
#                   scaffold=LONGEST_EGYPTREFV2_SCAFFOLDS)

# All versus all comparisons of reference and Egyptian genome
rule align_all_vs_all:
    input: expand("align_lastz_GRCh38_vs_EGYPTREF/{chr}_vs_{scaffold}.maf", \
                  chr=CHR_GRCh38, scaffold=EGYPTREF_SCAFFOLDS)

# Computing the GRCh38 recovery rate using the mafTools package 
# (as in Cho et al.). Using mafTools program mafPairCoverage, it is necessary
# to first combine all scaffold maf files for a chromosome, and then run 
# mafTransitiveClosure
rule combine_maf_files_for_recovery:
    input: expand("align_lastz_GRCh38_vs_EGYPTREF/{{chr}}_vs_{scaffold}.maf", \
                   scaffold=EGYPTREF_SCAFFOLDS)
    output: "align_lastz_GRCh38_vs_EGYPTREF/recovery/{chr}_alignments.maf"
    run: 
        shell("cat {input[0]} > {output}")
        for filename in input[1:]:
            # Append to large file; some file only have comments, no alignments
            # therefore we need to add & true because other wise the exit code
            # would indicate an error
            shell("cat {filename} | grep -v '#' >> {output} & true")

rule transitive_closure:
    input: "align_lastz_GRCh38_vs_EGYPTREF/recovery/{chr}_alignments.maf"
    output: "align_lastz_GRCh38_vs_EGYPTREF/recovery/{chr}.transclos"
    params: chr_number=lambda wildcards: wildcards.chr.split(".")[1]
    shell: "./ext_tools/mafTools/bin/mafTransitiveClosure " + \
           "--maf {input} > {output}"

rule maftools_coverage:
    input: "align_lastz_GRCh38_vs_EGYPTREF/recovery/{chr}.transclos"
    output: "align_lastz_GRCh38_vs_EGYPTREF/recovery/{chr}.coverage"
    params: chr_number=lambda wildcards: wildcards.chr.split(".")[1]
    shell: "./ext_tools/mafTools/bin/mafPairCoverage " + \
           "--maf {input} --seq1 {params.chr_number} --seq2 \* > {output}"

rule recovery:
    input: expand("align_lastz_GRCh38_vs_EGYPTREF/recovery/{chr}.coverage", \
                   chr=CHR_GRCh38)
    output: "align_lastz_GRCh38_vs_EGYPTREF/recovery/recovery.txt"
    run:
        pass


################################################################################
########### Reference to assembly genome alignment with mummer #################
################################################################################

# Genome alignments using mummer4
rule align_with_mummer:
    input: "repeatmasked_GRCh38/Homo_sapiens.GRCh38.dna.{chr}.fa.masked",
           "repeatmasked_EGYPTREF/Homo_sapiens.EGYPTREF.dna.{scaffold}.fa.masked"
    output: "align_mummer_GRCh38_vs_EGYPTREF/{chr}_vs_{scaffold}.delta"
    conda: "envs/mummer.yaml"
    shell: "nucmer " + \
           "-p align_mummer_GRCh38_vs_EGYPTREF/{wildcards.chr}_vs_{wildcards.scaffold} " + \
           "{input[0]} {input[1]}"

rule plot_mummer:
    input: "align_mummer_GRCh38_vs_EGYPTREF/{chr}_vs_{scaffold}.filter"
    output: "align_mummer_GRCh38_vs_EGYPTREF/dotplots/{chr}_vs_{scaffold}.gp",
            "align_mummer_GRCh38_vs_EGYPTREF/dotplots/{chr}_vs_{scaffold}.rplot",
            "align_mummer_GRCh38_vs_EGYPTREF/dotplots/{chr}_vs_{scaffold}.fplot",
            "align_mummer_GRCh38_vs_EGYPTREF/dotplots/{chr}_vs_{scaffold}.ps"
    conda: "envs/mummer.yaml"
    shell: "mummerplot " + \
           "--postscript " + \
           "-p align_mummer_GRCh38_vs_EGYPTREF/dotplots/{wildcards.chr}_vs_{wildcards.scaffold} " + \
           "{input[0]}; " + \
           "gnuplot {output[0]}"

# All versus all comparisons of reference and Egyptian genome
rule align_all_vs_all_mummer:
    input: expand("align_mummer_GRCh38_vs_EGYPTREF/{chr}_vs_{scaffold}.delta", \
                  chr=CHR_GRCh38, scaffold=EGYPTREF_SCAFFOLDS)

# Comparing the entire GRCh38 assembly with the entire EGYPTREF assembly
# --mum:         Use anchor matches that are unique in both the reference and 
#                query (false)
# --threads=NUM: Use NUM threads (# of cores)
# Path to conda environment:
# /data/lied_egypt_genome/lied_egypt_genome/.snakemake/conda/5fab0d7a
rule align_assemblies_with_mummer:
    input: ref="seq_{a1}/Homo_sapiens.{a1}.dna.primary_assembly.fa",
           query="seq_{a2}/Homo_sapiens.{a2}.dna.primary_assembly.fa"
    output: "align_mummer_{a1}_vs_{a2}/assemblies/{a1}_vs_{a2}.delta"
    conda: "envs/mummer.yaml"
    shell: "nucmer " + \
           "--mum " + \
           "--threads=24 "
           "-p align_mummer_{wildcards.a1}_vs_{wildcards.a2}/assemblies/{wildcards.a1}_vs_{wildcards.a2} " + \
           "{input[0]} {input[1]}"

rule align_assemblies_with_mummer_all:
    input: expand("align_mummer_GRCh38_vs_{a}/assemblies/GRCh38_vs_{a}.delta", \
                   a=["EGYPTREF","CEGYPTREF","EGYPTREFV2","AK1","YORUBA","GRCh38"])

# All versus all dotplots of reference and Egyptian genome
rule all_vs_all_dotplots_mummer:
    input: expand("align_mummer_GRCh38_vs_EGYPTREF/dotplots/{chr}_vs_{scaffold}.gp", \
                  chr=CHR_GRCh38, scaffold=EGYPTREF_SCAFFOLDS)

# Plotting the dotplots for all scaffolds
rule mummer_dotplots_scaffold_vs_chromosomes_all:
    input: expand("align_lastz_GRCh38_vs_EGYPTREF/dotplots/{scaffold}.pdf", \
                  scaffold=EGYPTREF_SCAFFOLDS)

# Filtering the mummer alignments: Query sequences can be mapped to reference 
# sequences with -q, this allows the user to exclude chance and repeat 
# alignments, leaving only the best alignments between the two data sets (i.e.
# use the -q option for mapping query contigs to their best reference location)
# -u: float; Set the minimum alignment uniqueness, i.e. percent of the alignment 
#     matching to unique reference AND query sequence [0, 100], default 0
# -l: int; Set the minimum alignment length, default 0
# -i: float; Set the minimum alignment identity [0, 100], default 0
# -1: 1-to-1 alignment allowing for rearrangements
# -r: Maps each position of each reference to its best hit in the query, 
#     allowing for query overlaps (intersection of -r and -q alignments)
rule delta_filter_mummer:
    input: "align_mummer_GRCh38_vs_EGYPTREF/{chr}_vs_{scaffold}.delta"
    output: "align_mummer_GRCh38_vs_EGYPTREF/{chr}_vs_{scaffold}.filter"
    conda: "envs/mummer.yaml"
    shell: "delta-filter -l 10000 -u 0 -q {input} > {output}"

# Running the tool nucdiff to compare two assemblies based on alignment with 
# mummer, which is also performed by the nucdiff tool
rule run_nucdiff:
    input: ref="seq_{a1}/Homo_sapiens.{a1}.dna.{chr}.fa", \
           query="seq_{a2}/Homo_sapiens.{a2}.dna.primary_assembly.fa"
    output: "nucdiff_{a1}_vs_{a2}/results/{a1}_vs_{a2}_{chr}_ref_snps.gff", \
            "nucdiff_{a1}_vs_{a2}/results/{a1}_vs_{a2}_{chr}_ref_struct.gff", \
            "nucdiff_{a1}_vs_{a2}/results/{a1}_vs_{a2}_{chr}_ref_blocks.gff", \
            "nucdiff_{a1}_vs_{a2}/results/{a1}_vs_{a2}_{chr}_ref_snps.vcf", \
            "nucdiff_{a1}_vs_{a2}/results/{a1}_vs_{a2}_{chr}_query_snps.gff", \
            "nucdiff_{a1}_vs_{a2}/results/{a1}_vs_{a2}_{chr}_query_struct.gff", \
            "nucdiff_{a1}_vs_{a2}/results/{a1}_vs_{a2}_{chr}_query_blocks.gff", \
            "nucdiff_{a1}_vs_{a2}/results/{a1}_vs_{a2}_{chr}_query_snps.vcf", \
            "nucdiff_{a1}_vs_{a2}/results/{a1}_vs_{a2}_{chr}_stat.out"
    params: outdir=lambda wildcards: "nucdiff_"+wildcards.a1+"_vs_"+wildcards.a2
    conda: "envs/nucdiff.yaml"
    shell: "nucdiff {input.ref} {input.query} {params.outdir} " + \
           "{wildcards.a1}_vs_{wildcards.a2}_{wildcards.chr} " + \
           "--vcf yes " + \
           "--filter_opt '-l 1000 -i 99' "
           "--proc 24"

rule run_nucdiff_all:
    input: expand("nucdiff_GRCh38_vs_EGYPTREF/results/GRCh38_vs_EGYPTREF_{chr}_stat.out", \
                  chr=CHR_GRCh38)

rule run_nucdiff_for_assembly:
    input: ref="seq_{a1}/Homo_sapiens.{a1}.dna.primary_assembly.fa", \
           query="seq_{a2}/Homo_sapiens.{a2}.dna.primary_assembly.fa", \
           delta="align_mummer_{a1}_vs_{a2}/assemblies/{a1}_vs_{a2}.delta"
    output: "nucdiff_{a1}_vs_{a2}/assemblies/results/{a1}_vs_{a2}_ref_snps.gff", \
            "nucdiff_{a1}_vs_{a2}/assemblies/results/{a1}_vs_{a2}_ref_struct.gff", \
            "nucdiff_{a1}_vs_{a2}/assemblies/results/{a1}_vs_{a2}_ref_blocks.gff", \
            "nucdiff_{a1}_vs_{a2}/assemblies/results/{a1}_vs_{a2}_ref_snps.vcf", \
            "nucdiff_{a1}_vs_{a2}/assemblies/results/{a1}_vs_{a2}_query_snps.gff", \
            "nucdiff_{a1}_vs_{a2}/assemblies/results/{a1}_vs_{a2}_query_struct.gff", \
            "nucdiff_{a1}_vs_{a2}/assemblies/results/{a1}_vs_{a2}_query_blocks.gff", \
            "nucdiff_{a1}_vs_{a2}/assemblies/results/{a1}_vs_{a2}_query_snps.vcf", \
            "nucdiff_{a1}_vs_{a2}/assemblies/results/{a1}_vs_{a2}_stat.out"
    params: outdir=lambda wildcards: "nucdiff_"+wildcards.a1+"_vs_"+wildcards.a2
    conda: "envs/nucdiff.yaml"
    shell: "nucdiff {input.ref} {input.query} {params.outdir} " + \
           "{wildcards.a1}_vs_{wildcards.a2} " + \
           "--vcf yes " + \
           "--filter_opt '-l 1000 -i 99' " + \
           "--delta_file {input.delta} " + \
           "--proc 24"

################################################################################
######## Processing Illumina PE data for the assembled individual ##############
################################################################################

# The Illumina library sample names
ILLUMINA_SAMPLES = ["NDES00177","NDES00178","NDES00179","NDES00180","NDES00181"]
ILLUMINA_SAMPLES_TO_LANES = {
    "NDES00177": [4,5,6,7],
    "NDES00178": [1,4,5,6,7],
    "NDES00179": [4,5,6,7],
    "NDES00180": [1,4,5,6,7],
    "NDES00181": [4,5,6,7]
}
ILLUMINA_LIBS = []
for sample in ILLUMINA_SAMPLES:
    ILLUMINA_LIBS += [sample+"_L"+str(x) for x in \
                      ILLUMINA_SAMPLES_TO_LANES[sample]]

# Mapping the Illumina PE data to the scaffolds
# -a STR: Algorithm for constructing BWT index. Chosen option: 
#         bwtsw: Algorithm implemented in BWT-SW. This method works with the 
#         whole human genome.
# -p STR: Prefix of the output database [same as db filename] 
rule bwa_index:
    input: "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.fa"
    output: "bwa_index/Homo_sapiens.{assembly}.dna.primary_assembly.amb",
            "bwa_index/Homo_sapiens.{assembly}.dna.primary_assembly.ann",
            "bwa_index/Homo_sapiens.{assembly}.dna.primary_assembly.bwt",
            "bwa_index/Homo_sapiens.{assembly}.dna.primary_assembly.pac",
            "bwa_index/Homo_sapiens.{assembly}.dna.primary_assembly.sa"
    conda: "envs/bwa.yaml"
    shell: "bwa index -a bwtsw " + \
                     "-p bwa_index/Homo_sapiens." + \
                     "{wildcards.assembly}.dna.primary_assembly " + \
                     "{input}"

rule bwa_index_all:
    input: expand("bwa_index/Homo_sapiens.{assembly}.dna.primary_assembly.sa", \
                  assembly=["EGYPTREF","GRCh38"])

rule bwa_mem:
    input: index = "bwa_index/Homo_sapiens.{assembly}.dna.primary_assembly.sa",
           fastq_r1 = "data/02.DES/{lib}_1.fq.gz",
           fastq_r2 = "data/02.DES/{lib}_2.fq.gz"
    output: "map_bwa_{assembly}/{lib}.bam"
    shell: "bwa mem -t 48 " + \
           "bwa_index/Homo_sapiens.{wildcards.assembly}.dna.primary_assembly "+\
           "{input.fastq_r1} {input.fastq_r2} " + \
           " | samtools sort -@48 -o {output} -"

rule bwa_mem_all:
    input: expand("map_bwa_{assembly}/{lib}.bam", \
                  assembly=["EGYPTREF","GRCh38"], lib=ILLUMINA_LIBS)

# Making some stats about the overall number of bases etc.
rule illumina_bamstats:
    input: "map_bwa_{assembly}/{lib}.bam"
    output: "map_bwa_{assembly}/{lib}.bamstats"
    shell: "samtools stats {input} > {output}"

rule illumina_bamstats_all:
    input: expand("map_bwa_{assembly}/{lib}.bamstats", \
                  assembly=["EGYPTREF","GRCh38"], lib=ILLUMINA_LIBS)

# For SNP calling and other things that are done for the Illumnina PE data
rule symlink_illumina_wgs_dir:
    output: directory("data/02.DES")
    shell: "ln -s /data/lied_egypt_genome/raw/P101HW18010820-01_human_2018.08.29/00.data/02.DES {output}"

# Some QC: Here, fastqc for all Illumina PE WGS files
rule run_fastqc:
    input: "data/02.DES/{lib}_{read}.fq.gz"
    output: html="illumina_qc/fastqc/{lib}_{read}_fastqc.html",
            zip="illumina_qc/fastqc/{lib}_{read}_fastqc.zip"
    conda: "envs/fastqc.yaml"
    shell: "fastqc --outdir illumina_qc/fastqc/ {input[0]}"

rule run_fastqc_all:
    input: expand("illumina_qc/fastqc/{lib}_{read}_fastqc.html", \
                  lib=ILLUMINA_LIBS, read=["1","2"])


################################################################################
################### SNP Calling for 10 Egyptian individuals ####################
#################### (one the reference genome individual) #####################
################################################################################

# If possible, the variant calling (vc) tasks are performed with 24 threads and
# with 80Gb of memory, such that 2 tasks can be run on a node in parallel;
# The reason is that the 80Gb of memory are needed and thus only two tasks per
# node can be run anyway.

# These are the sample IDs
# The EGYPTREF Illumina sequencing data is actually not in the same folder, but
# in order to do the SNP calling with this sample the same way as with the 
# other 9 samples, I made a folder EGYPTREF there in which I generated symlinks
# which link to the fastq files within the 02.DES folder of the data from 
# EGYPTREF and I symlinked the individual fastq files and used a naming 
# convention that follows the same pattern than that for the 9 Egyptian samples
EGYPT_SAMPLES = ["EGYPTREF","LU18","LU19","LU2","LU22","LU23","LU9","PD114", \
                 "PD115","PD82"] #,"TEST"

# These are additional IDs after the sample IDs, given by Novogene, e.g 
# H75TCDMXX is the ID of the sequencer, L1 is the first lane
EGYPT_SAMPLES_TO_PREPLANES = {
    "EGYPTREF": ["NDES00177_L4","NDES00177_L5","NDES00177_L6","NDES00177_L7", \
                 "NDES00178_L1","NDES00178_L4","NDES00178_L5","NDES00178_L6", \
                 "NDES00178_L7", \
                 "NDES00179_L4","NDES00179_L5","NDES00179_L6","NDES00179_L7", \
                 "NDES00180_L1","NDES00180_L4","NDES00180_L5","NDES00180_L6", \
                 "NDES00180_L7", \
                 "NDES00181_L4", "NDES00181_L5","NDES00181_L6", "NDES00181_L7"],
    "LU18":     ["NDHG02363_H75HVDMXX_L1", "NDHG02363_H75TCDMXX_L1", \
                 "NDHG02363_H75HVDMXX_L2", "NDHG02363_H75TCDMXX_L2", \
                 "NDHG02363_H75FVDMXX_L1", "NDHG02363_H75FVDMXX_L2"],
    "LU19":     ["NDHG02358_H7777DMXX_L1", "NDHG02358_H7777DMXX_L2"],
    "LU2":      ["NDHG02365_H75FVDMXX_L1", "NDHG02365_H75FVDMXX_L1"],
    "LU22":     ["NDHG02364_H75LLDMXX_L1", "NDHG02364_H75LLDMXX_L2"],
    "LU23":     ["NDHG02366_H75FVDMXX_L1", "NDHG02366_H75FVDMXX_L2"],
    "LU9":      ["NDHG02362_H772LDMXX_L1", "NDHG02362_H772LDMXX_L2"],
    "PD114":    ["NDHG02360_H772LDMXX_L1", "NDHG02360_H772LDMXX_L2"],
    "PD115":    ["NDHG02361_H772LDMXX_L1", "NDHG02361_H772LDMXX_L2"],
    "PD82":     ["NDHG02359_H772LDMXX_L1", "NDHG02359_H772LDMXX_L2"]#,
#    "TEST":     ["PROTOCOL_SEQUENCER_L1", "PROTOCOL_SEQUENCER_L2"] # the last is for testing purposes
}

# Symlinking the raw data directory
rule symlink_data_for_variant_detection:
    output: directory("data/raw_data")
    shell: "ln -s /data/lied_egypt_genome/raw_data {output}"

# Getting the latest dbsnp version for GRCh38, this is version 151; I am 
# getting the VCF file deposited under GATK, which is very slightly larger than
# the file under VCF, but I didn't check the precise difference and there is 
# no note in the READMEs.
rule get_known_snps_from_dbsnp:
    output: "dbsnp_GRCh38/All_20180418.vcf.gz"
    shell: "wget -P dbsnp_GRCh38 " + \
           "ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/GATK/All_20180418.vcf.gz"

# ... and getting its index
rule get_index_of_known_snps_from_dbsnp:
    output: "dbsnp_GRCh38/All_20180418.vcf.gz.tbi"
    shell: "wget -P dbsnp_GRCh38 " + \
           "ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/GATK/All_20180418.vcf.gz.tbi"

# Replacing trailing "chr" and then zipping...
rule unzip_and_reformat_dbsnp:
    input: "dbsnp_GRCh38/All_20180418.vcf.gz"
    output: "dbsnp_GRCh38/dbsnp.vcf.gz"
    shell: "zcat {input} | " + \
           "awk '{{gsub(/^chr/,\"\"); print}}' | " + \
           "bgzip > {output}"

# ... and indexing again
rule tabix_dbsnp:
    input: "dbsnp_GRCh38/dbsnp.vcf.gz"
    output: "dbsnp_GRCh38/dbsnp.vcf.gz.tbi"
    shell: "tabix -p vcf {input}"

### 1. map reads to genome
# Mapping to reference/assembly using bwa
rule vc_bwa_mem:
    input: index = "bwa_index/Homo_sapiens.{assembly}.dna.primary_assembly.sa",
           fastq_r1 = "data/raw_data/{sample}/{sample}_{infolane}_1.fq.gz",
           fastq_r2 = "data/raw_data/{sample}/{sample}_{infolane}_2.fq.gz"
    output: temp("variants_{assembly}/{sample}_{infolane}.sam")
    wildcard_constraints: sample="[A-Z,0-9]+", infolane="[A-Z,0-9,_]+"
    conda: "envs/variant_calling.yaml"
    shell: "bwa mem -t 24 " + \
           "bwa_index/Homo_sapiens.{wildcards.assembly}.dna.primary_assembly " + \
           "{input.fastq_r1} {input.fastq_r2} > {output}"

### 2. CleanSam
# Cleans the provided SAM/BAM, soft-clipping beyond-end-of-reference alignments 
# and setting MAPQ to 0 for unmapped reads
# java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar share/picard-2.18.9-0/picard.jar
rule vc_clean_sam:
    input: "variants_{assembly}/{sample}_{infolane}.sam"
    output: temp("variants_{assembly}/{sample}_{infolane}.cleaned.sam")
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx80g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \ 
           "CleanSam " + \
           "I={input} " + \
           "O={output}"

### 3. Sort Sam -> output: bam + idx
# Sorting by coordinates, making an index and outputting as bam
rule vc_sort_and_index_sam:
    input: "variants_{assembly}/{sample}_{infolane}.cleaned.sam"
    output: temp("variants_{assembly}/{sample}_{infolane}.cleaned.bam")
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx80g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "SortSam " + \
           "I={input} " + 
           "O={output} " + \
           "SORT_ORDER=coordinate " + \
           "CREATE_INDEX=true"

### 4. Fix Mate Pair Information
# verify mate-pair information between mates and fix if needed
rule vc_fix_mates:
    input: "variants_{assembly}/{sample}_{infolane}.cleaned.bam"
    output: temp("variants_{assembly}/{sample}_{infolane}.fixed.bam")
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx80g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "FixMateInformation " + \
           "I={input} " + \
           "O={output} " + \
           "SORT_ORDER=coordinate " + \
           "CREATE_INDEX=true"

### 5. Mark Duplicates
rule vc_mark_duplicates:
    input: "variants_{assembly}/{sample}_{infolane}.fixed.bam"
    output: temp("variants_{assembly}/{sample}_{infolane}.rmdup.bam"),
            "variants_{assembly}/{sample}_{infolane}.rmdup.txt"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx80g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "MarkDuplicates " + \
           "I={input} " + \
           "O={output[0]} " + \
           "M={output[1]} " + \
           "REMOVE_DUPLICATES=true "+ \
           "ASSUME_SORTED=coordinate " + \
           "CREATE_INDEX=true"

### 6. merge *.bam files
# Merging all bam Files for a sample
rule vc_merge_bams_per_sample:
    input: lambda wildcards: \
           expand("variants_{assembly}/{sample}_{infolane}.rmdup.bam", \
           assembly = wildcards.assembly, \
           sample=wildcards.sample, \
           infolane = EGYPT_SAMPLES_TO_PREPLANES[wildcards.sample])
    output: temp("variants_{assembly}/{sample}.merged.bam")
    params:
        picard_in=lambda wildcards, input: "I="+" I=".join(input)
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx80g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "MergeSamFiles " + \
           "{params.picard_in} " + \
           "O={output} " + \
           "SORT_ORDER=coordinate " + \
           "CREATE_INDEX=true " + \
           "USE_THREADING=24"

### 7. Collect Alignment Summary Metrics
rule vc_alignment_metrics:
    input: "variants_{assembly}/{sample}.merged.bam"
    output: "variants_{assembly}/{sample}.stats.txt"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx80g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "CollectAlignmentSummaryMetrics " + \
           "I={input} " + \
           "O={output}"

### 8. Replace Read Groups
rule vc_replace_read_groups:
    input: "variants_{assembly}/{sample}.merged.bam"
    output: temp("variants_{assembly}/{sample}.merged.rg.bam")
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx80g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "AddOrReplaceReadGroups " + \
           "I={input} " + \
           "O={output} " + \
           "RGID={wildcards.sample} " + \
           "RGPL=illumina " + \
           "RGLB={wildcards.sample} " + \
           "RGPU=unit1 " + \
           "RGSM={wildcards.sample} " + \
           "CREATE_INDEX=true"

### 9. Realign

# Therefore, generate a sequence dictionary for use with picard tools first
rule vc_seq_dict:
    input: "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.fa"
    output: "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.dict"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx80g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "CreateSequenceDictionary " + \
           "R={input} " + \
           "O={output}"

rule vc_reorder:
    input: "variants_{assembly}/{sample}.merged.rg.bam",
           "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.fa",
           "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.dict"
    output: temp("variants_{assembly}/{sample}.merged.rg.ordered.bam")
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx80g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "ReorderSam " + \
           "I={input[0]} " + \
           "O={output} " + \
           "R={input[1]} " + \
           "CREATE_INDEX=true"

### Formally, 10. RealignerTargetCreator and 11. IndelRealigner
# Note: I do not perform this anymore, since according to the GATK website: 
# "Note that indel realignment is no longer necessary for variant discovery if 
# you plan to use a variant caller that performs a haplotype assembly step, such
# as HaplotypeCaller or MuTect2. "

### 12. Base Quality Recalibration
# --knownSites / -knownSites: A database of known polymorphic sites. This 
#                             algorithm treats every reference mismatch as an 
#                             indication of error. However, real genetic 
#                             variation is expected to mismatch the reference, 
#                             so it is critical that a database of known 
#                             polymorphic sites (e.g. dbSNP) is given to the 
#                             tool in order to mask out those sites.
# --covariate / -cov: One or more covariates to be used in the recalibration. 
#                     Can be specified multiple times. Note that the ReadGroup 
#                     and QualityScore covariates are required and do not need 
#                     to be specified. Also, unless --no_standard_covs is 
#                     specified, the Cycle and Context covariates are standard 
#                     and are included by default. Use the --list argument to 
#                     see the available covariates. 

# Therefore, the fasta file needs to be indexed
rule vc_index_fasta:
    input: "seq_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output: "seq_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
    shell: "samtools faidx {input}"

rule vc_base_recalibrator:
    input: "variants_{assembly}/{sample}.merged.rg.ordered.bam",
           "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.fa"#,
#           "dbsnp_{assembly}/All_20180418.vcf"
#           "dbsnp_{assembly}/All_20180418.vcf.tbi"
    output: "variants_{assembly}/{sample}.recal.csv"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx80g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/opt/gatk-3.8/GenomeAnalysisTK.jar " + \
           "-T BaseRecalibrator " + \
           "-R {input[1]} " + \
           "-I {input[0]} " + \
           "-cov ReadGroupCovariate " + \
           "-cov QualityScoreCovariate " + \
           "-cov CycleCovariate " + \
           "-cov ContextCovariate " + \
           "-o {output} " + \
#           "-knownSites {input[2]} " + \
           "-nct 24"

### 13. Print Reads
rule vc_print_reads:
    input: "variants_{assembly}/{sample}.merged.rg.ordered.bam",
           "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.fa",
           "variants_{assembly}/{sample}.recal.csv"
    output: "variants_{assembly}/{sample}.final.bam"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx80g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/opt/gatk-3.8/GenomeAnalysisTK.jar " + \
           "-T PrintReads " + \
           "-R {input[1]} " + \
           "-I {input[0]} " + \
           "-o {output[0]} " + \
           "-BQSR {input[2]} " + \
           "-nct 24"

# Making some simple stats about the number of mapped reads etc.
rule vc_flagstat:
    input: "variants_{assembly}/{sample}.final.bam"
    output: "variants_{assembly}/{sample}.final.flagstat.txt"
    shell: "samtools flagstat {input} > {output}"

# Making some stats about the overall number of bases etc.
rule vc_bamstats:
    input: "variants_{assembly}/{sample}.final.bam"
    output: "variants_{assembly}/{sample}.final.bamstats.txt"
    shell: "samtools stats {input} > {output}"

rule vc_flagstat_all:
    input: expand("variants_{assembly}/{sample}.final.flagstat.txt", \
                  assembly=["GRCh38"], sample=EGYPT_SAMPLES),
           expand("variants_{assembly}/{sample}.final.bamstats.txt", \
                  assembly=["GRCh38"], sample=EGYPT_SAMPLES)

### 14. variant calling with GATK-HC
# use GATK Haplotypecaller with runtime-optimized settings
# --genotyping_mode / -gt_mode: Specifies how to determine the alternate alleles
#                               to use for genotyping (DISCOVERY: The genotyper 
#                               will choose the most likely alternate allele
# --dbsnp / -D: dbSNP file rsIDs from this file are used to populate the ID 
#               column of the output. Also, the DB INFO flag will be set when 
#               appropriate. dbSNP is not used in any way for the calculations
#               themselves. 
# --emitRefConfidence / -ERC: Mode for emitting reference confidence scores
#                             Records whether the trimming intervals are going 
#                             to be used to emit reference confidence
#                             GVCF: Reference model emitted with condensed 
#                             non-variant blocks, i.e. the GVCF format. 
# -variant_index_type LINEAR -variant_index_parameter 128000 IW added because
# the GATK program told me so and otherwise would exit with error.
# java -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=4 -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar gatk 
rule vc_snp_calling_with_gatk_hc:
    input: "variants_{assembly}/{sample}.final.bam",
           "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.fa"#,
#           "dbsnp_{assembly}/All_20180418.vcf"#,
#           "dbsnp_{assembly}/All_20180418.vcf.tbi"
    output: "variants_{assembly}/{sample}.vcf"
    wildcard_constraints: sample="[A-Z,0-9]+"
    conda: "envs/variant_calling.yaml"
    shell: "java -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=4 -Xmx80g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/opt/gatk-3.8/GenomeAnalysisTK.jar " + \
           "-T HaplotypeCaller " + \
           "-R {input[1]} " + \
           "-I {input[0]} " + \
#           "--dbsnp {input[2]} " + \
           "--genotyping_mode DISCOVERY " + \
           "-o {output} " + \
           "-ERC GVCF " + \
           "-variant_index_type LINEAR " + \
           "-variant_index_parameter 128000 " + \
           "-nct 24"

### 15. Perform joint genotyping on gVCF files produced by HaplotypeCaller

# Therefore, split GATK HC files chromosome-wise
rule vc_split_gatk_hc_chromosomewise:
    input: "variants_{assembly}/{sample}.vcf"
    output: "variants_{assembly}/{sample}.chromosome.{chr}.vcf"
    shell: "head -n 1000 {input} | grep '#' > {output}; " + \
           "cat {input} | grep -P '^{wildcards.chr}\t' >> {output}; "

rule vc_split_gatk_hc_chromosomewise_all:
    input: expand("variants_{assembly}/{sample}.{chr}.vcf", \
                  assembly=["GRCh38"], \
                  sample=EGYPT_SAMPLES, \
                  chr=CHR_GRCh38)

# --intervals / -L: One or more genomic intervals over which to operate
rule vc_joint_genotyping:
    input: vcfs=expand("variants_{{assembly}}/{sample}.chromosome.{{chr}}.vcf", \
                       sample=EGYPT_SAMPLES),
           ref="seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.fa",
#           dbsnp="dbsnp_{assembly}/dbsnp.vcf.gz"
    output: "variants_{assembly}/egyptians.chromosome.{chr}.vcf"
    conda: "envs/variant_calling.yaml"
    params: variant_files=lambda wildcards, input: " --variant " + \
                                                   " --variant ".join(input.vcfs)
    shell: "rm -rf /scratch/tmp_gatk_{wildcards.chr}; " + \
           "mkdir /scratch/tmp_gatk_{wildcards.chr}; " + \
           "java -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=4 -Xmx80g " + \
           "-Djava.io.tmpdir=/scratch/tmp_gatk_{wildcards.chr} " + \
           "-jar .snakemake/conda/d590255f/opt/gatk-3.8/GenomeAnalysisTK.jar " + \
           "-T GenotypeGVCFs " + \
           "-R {input.ref} " + \
           "{params.variant_files} " + \
#           "--dbsnp {input.dbsnp} " + \
           "--num_threads 24 " + \
           "--intervals {wildcards.chr} " + \
           "-o {output}; " + \
           "rm -r /scratch/tmp_gatk_{wildcards.chr}; "

rule vc_compress_vcf:
    input: "variants_{assembly}/egyptians.{chr}.vcf"
    output: "variants_{assembly}/egyptians.{chr}.vcf.gz"
    shell: "cat {input} | bgzip > {output} "

rule vc_joint_genotyping_all:
    input: expand("variants_{assembly}/egyptians.{chr}.vcf.gz", \
                  assembly=["GRCh38"],chr=CHR_GRCh38)

# Combining the chromosome-wise vcfs
rule vc_combine_chromosomewise_vcf:
    input: expand("variants_{{assembly}}/egyptians.{chr}.vcf.gz", chr=CHR_GRCh38)
    output: "variants_{assembly}/egyptians.vcf.gz"
    conda: "envs/genotype_pcs.yaml"
    shell: "vcf-concat --pad-missing {input} | bgzip > {output}"

# and indexing them
rule vc_index_egyptians:
    input: "variants_{assembly}/egyptians.vcf.gz"
    output: "variants_{assembly}/egyptians.vcf.gz.tbi"
    conda: "envs/genotype_pcs.yaml"    
    shell: "tabix -p vcf {input}"

# Doing the variant calling for all 10 samples
# and collecting alignment summary stats
rule vc_snp_calling_with_gatk_hc_all:
    input: expand("variants_GRCh38/{sample}.vcf", sample=EGYPT_SAMPLES),
           expand("variants_GRCh38/{sample}.stats.txt", sample=EGYPT_SAMPLES)

# Extracting unmapped reads of the reference Egyptian for Axel, who wants to
# run SOAP-denovo to compute a de-novo short-read assembly based on them
rule get_unmapped_reads:
    input: "variants_{assembly}/{sample}.final.bam"
    output: "variants_{assembly}/{sample}.final.unmapped.bam"
    shell: "samtools view -b -f 4 {input} > {output}"

# Extracting reads that map to the mitochondrium
# -F 4 extract mapped reads
# -u output uncompressed sam
# -b output bam
# What's done: 1) Writing header, 2) Getting MT aligned reads, 3) sam to bam
# conversion 4) Removing intermediate sam file
rule mt_reads:
    input: "variants_{assembly}/{sample}.final.bam"
    output: "variants_{assembly}/{sample}.final.MT.bam"
    shell: "samtools view -H {input} > {output}.sam; " + \
           "samtools view {input} | grep -P \"\tMT\t\" >> {output}.sam; " + \
           "samtools view -hb {output}.sam > {output}; " + \
           "rm {output}.sam"


################################################################################
##### Assembly assessment and correction (Things related to PacBio data) #######
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

rule symlink_pacbio:
    output: directory("data/01.pacbio")
    shell: "ln -s /data/lied_egypt_genome/raw/P101HW18010820-01_human_2018.08.29/00.data/01.pacbio {output}"

# Making some stats about the overall number of bases etc.
rule pb_bamstats:
    input: "data/01.pacbio/{pb_files}.subreads.bam"
    output: "pacbio/{pb_files}.bamstats"
    shell: "samtools stats {input} > {output}"

rule pb_bamstats_all:
    input: expand("pacbio/{pb_files}.bamstats", \
            pb_files = [item for subl in PACBIO_SAMPLES_TO_SEQRUN_PATH.values() \
                       for item in subl])

# Conversion of bam to fasta sequences
rule pb_bam2fasta:
    input: "data/01.pacbio/{pb_files}.subreads.bam"
    output: temp("pacbio/{pb_files}.fa")
    shell: "samtools fasta -t {input} > {output}"

# Constructing one large fasta file
rule pb_combined_fasta:
    input: expand("pacbio/{pb_files}.fa", \
            pb_files = [item for subl in PACBIO_SAMPLES_TO_SEQRUN_PATH.values() \
                       for item in subl])
    output: "pacbio/pb_EGYPTREF.fa"
    shell: "cat {input} > {output}"

# Get and install assembler
rule install_wtdbg2_assembler:
    output: "wtdbg2/wtdbg2",
            "wtdbg2/wtpoa-cns"
    shell: "git clone https://github.com/ruanjue/wtdbg2; " + \
           "cd wtdbg2 && make "

# Run assembler wtdbg2
# WTDBG: De novo assembler for long noisy sequences
# Author: Jue Ruan <ruanjue@gmail.com>
# Version: 2.2 (20181111)
# Usage: wtdbg2 [options] -i <reads.fa> -o <prefix> [reads.fa ...]
# Options:
#  -i <string> Long reads sequences file (REQUIRED; there can be multiple -i), []
#  -o <string> Prefix of output files (REQUIRED), []
#  -t <int>    Number of threads, 0 for all cores, [4]
#  -f          Force to overwrite output files
#  -x <string> Presets, comma delimited, []
#             rsII/rs: -p 21 -S 4 -s 0.05 -L 5000
#           sequel/sq: -p 0 -k 15 -AS 2 -s 0.05 -L 5000
#        nanopore/ont: -p 19 -AS 2 -s 0.05 -L 5000
#       corrected/ccs: -p 0 -k 19 -AS 4 -s 0.5 -L 5000
#              Example: '-e 3 -x ont -S 1' in parsing order, -e will be infered by seq depth, -S will be 1
#  -g <number> Approximate genome size (k/m/g suffix allowed) [0]
#  -X <float>  Choose the best <float> depth from input reads(effective with -g) [50]
#  -L <int>    Choose the longest subread and drop reads shorter than <int> [0]
#  -k <int>    Kmer fsize, 0 <= k <= 25, [0]
#  -p <int>    Kmer psize, 0 <= p <= 25, [21]
#              k + p <= 25, seed is <k-mer>+<p-homopolymer-compressed>
#  -K <float>  Filter high frequency kmers, maybe repetitive, [1000.05]
#              >= 1000 and indexing >= (1 - 0.05) * total_kmers_count
#  -E <int>    Min kmer frequency, [2]
#  -S <float>  Subsampling kmers, 1/(<-S>) kmers are indexed, [4.00]
#              -S is very useful in saving memeory and speeding up
#              please note that subsampling kmers will have less matched length
#  -l <float>  Min length of alignment, [2048]
#  -m <float>  Min matched length by kmer matching, [200]
#  -A          Keep contained reads during alignment
#  -s <float>  Min similarity, calculated by kmer matched length / aligned length, [0.05]
#  -e <int>    Min read depth of a valid edge, [3]
#  -q          Quiet
#  -v          Verbose (can be multiple)
#  -V          Print version information and then exit
#  --help      Show more options
rule assembl_with_wtdbg2:
    input: fasta="pacbio/pb_EGYPTREF.fa",
           binary="wtdbg2/wtdbg2"
    output: "pacbio/assembly_wtdbg2/EGYPTREF_wtdbg2.ctg.lay.gz"
    conda: "envs/wtdbg.yaml"
    params: out_base=lambda wildcards, output: output[0][:-11]
    shell: "./{input.binary} -i {input.fasta} " + \
                           " -t 48 " + \
                           " -g 3g " + \
                           " -L 10000 " + \
#                           " -S 6 " + \
                           " -o {params.out_base} "

# WTPOA-CNS: Consensuser for wtdbg using PO-MSA
# Author: Jue Ruan <ruanjue@gmail.com>
# Version: 1.2
# Usage: wtpoa-cns [options]
# Options:
#  -t <int>    Number of threads, [4]
#  -d <string> Reference sequences for SAM input, will invoke sorted-SAM input mode and auto set '-j 100 -W 0 -w 1000'
#  -u          Only process reference regions present in/between SAM alignments
#  -p <string> Similar with -d, but translate SAM into wtdbg layout file
#  -i <string> Input file(s) *.ctg.lay from wtdbg, +, [STDIN]
#              Or sorted SAM files when having -d
#  -o <string> Output files, [STDOUT]
#  -f          Force overwrite
#  -j <int>    Expected max length of node, or say the overlap length of two adjacent units in layout file, [1500] bp
#  -M <int>    Match score, [2]
#  -X <int>    Mismatch score, [-5]
#  -I <int>    Insertion score, [-2]
#  -D <int>    Deletion score, [-4]
#  -B <int>    Bandwidth, [96]
#  -W <int>    Window size in the middle of the first read for fast align remaining reads, [200]
#              If $W is negative, will disable fast align, but use the abs($W) as Band align score cutoff
#  -w <int>    Min size of aligned size in window, [$W * 0.5]
#              In sorted-SAM input mode, -w is the sliding window size [2000]
#  -A          Abort TriPOA when any read cannot be fast aligned, then try POA
#  -R <int>    Realignment bandwidth, 0: disable, [16]
#  -C <int>    Min count of bases to call a consensus base, [3]
#  -F <float>  Min frequency of non-gap bases to call a consensus base, [0.5]
#  -N <int>    Max number of reads in PO-MSA [20]
#              Keep in mind that I am not going to generate high accurate consensus sequences here
#  -v          Verbose
#  -V          Print version information and then exit
rule consensus_with_wtdbg2:
    input: contig_layout="pacbio/assembly_wtdbg2/EGYPTREF_wtdbg2.ctg.lay.gz",
           binary="wtdbg2/wtpoa-cns"
    output: "pacbio/assembly_wtdbg2/EGYPTREF_wtdbg2.ctg.lay.fa"
    conda: "envs/wtdbg.yaml"
    shell: "./{input.binary} -i {input.contig_layout} " + \
                           " -t 48 " + \
                           " -o {output[0]}"


################################################################################
##### Assembly assessment and correction (Things related to 10X data) ##########
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

rule symlink_10X:
    output: directory("data/03.10X")
    shell: "ln -s /data/lied_egypt_genome/raw/P101HW18010820-01_human_2018.08.29/00.data/03.10X {output}"

# Some QC: Here, fastqc for all Illumina PE WGS files
rule run_10x_fastqc:
    input: "data/03.10X/{lib}_{read}.fq.gz"
    output: html="illumina_10x_qc/fastqc/{lib}_{read}_fastqc.html",
            zip="illumina_10x_qc/fastqc/{lib}_{read}_fastqc.zip"
    conda: "envs/fastqc.yaml"
    shell: "fastqc --outdir illumina_10x_qc/fastqc/ {input[0]}"

rule run_10x_fastqc_all:
    input: expand("illumina_10x_qc/fastqc/{lib}_{read}_fastqc.html", \
                  lib=ILLUMINA_10X_LIBS, read=["1","2"])


################################################################################
##### Population stratification analysis using Eigenstrat (e.g. PC plots) ######
################################################################################

# Downloading 1000 genomes data
rule download_1000g_genotypes:
    output: "1000_genomes/ALL.chr{chr}_GRCh38.genotypes.20170504.vcf.gz"
    shell: "wget -P 1000_genomes/ http://ftp.1000genomes.ebi.ac.uk/vol1/" + \
                                  "ftp/release/20130502/supporting/" + \
                                  "GRCh38_positions/" + \
                                  "ALL.chr{wildcards.chr}_GRCh38.genotypes.20170504.vcf.gz"

# Downloading 1000 genomes data (index)
rule download_1000g_genotypes_index:
    output: "1000_genomes/ALL.chr{chr}_GRCh38.genotypes.20170504.vcf.gz.tbi"
    shell: "wget -P 1000_genomes/ http://ftp.1000genomes.ebi.ac.uk/vol1/" + \
                                  "ftp/release/20130502/supporting/" + \
                                  "GRCh38_positions/" + \
                                  "ALL.chr{wildcards.chr}_GRCh38.genotypes.20170504.vcf.gz.tbi"

# Downloading 1000 genomes data (Readme)
rule download_1000g_genotypes_readme:
    output: "1000_genomes/README_GRCh38_liftover_20170504.txt"
    shell: "wget -P 1000_genomes/ http://ftp.1000genomes.ebi.ac.uk/vol1/" + \
                                  "ftp/release/20130502/supporting/" + \
                                  "GRCh38_positions/" + \
                                  "README_GRCh38_liftover_20170504.txt"

# Get the ped file which contains the population of the samples (and more info)
rule download_1000g_genotypes_ped:
    output: "1000_genomes/integrated_call_samples_v2.20130502.ALL.ped"
    shell: "wget -P 1000_genomes/ http://ftp.1000genomes.ebi.ac.uk/vol1/" + \
                                  "ftp/release/20130502/" + \
                                  "integrated_call_samples_v2.20130502.ALL.ped"

rule download_1000g_genotypes_all:
    input: expand("1000_genomes/ALL.chr{chr}_GRCh38.genotypes.20170504.vcf.gz", \
                   chr=[str(x) for x in range(1,23)]+["X","Y"]), \
           expand("1000_genomes/ALL.chr{chr}_GRCh38.genotypes.20170504.vcf.gz.tbi", \
                   chr=[str(x) for x in range(1,23)]+["X","Y"]), \
           "1000_genomes/README_GRCh38_liftover_20170504.txt", \
           "1000_genomes/integrated_call_samples_v2.20130502.ALL.ped"

# Selecting 1000G individuals for inclusion in Stratification analysis
# We select individuals belonging to populations of interest and which
# have phase3 genotypes avaliable (this is column 14) 
# ACB: African Caribbeans in Barbados
# ASW: Americans of African Ancestry in SW USA 	
# CEU: Utah Residents (CEPH) with Northern and Western European Ancestry
# ESN: Esan in Nigeria
# FIN: Finnish in Finland
# GBR: British in England and Scotland
# GWD: Gambian in Western Divisions in the Gambia
# IBS: Iberian Population in Spain
# LWK: Luhya in Webuye, Kenya
# MSL: Mende in Sierra Leone
# TSI: Toscani in Italia
# YRI: Yoruba in Ibadan, Nigeria
POPULATIONS_1000G = [
"ACB","ASW","CEU","ESN","FIN","GBR","GWD","IBS","LWK","MSL","TSI","YRI"
]
POPULATIONS_AFR = ["ACB","ASW","ESN","GWD","LWK","MSL","YRI"]
POPULATIONS_EUR = ["CEU","FIN","GBR","IBS","TSI"]
rule gp_select_pop_from_1000g:
    input: "1000_genomes/integrated_call_samples_v2.20130502.ALL.ped"
    output: "genotype_pcs/keep_indiv.txt"
    run:
        with open(input[0],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                s = line.split("\t")
                if s[6] in POPULATIONS_1000G and s[13] == "1":
                    f_out.write(s[1]+"\n")

# Make an annotation file with sample names and population for plotting
rule gp_make_pop_annotation:
    input: "1000_genomes/integrated_call_samples_v2.20130502.ALL.ped"
    output: "genotype_pcs/annotation_EGYPT_AFR_EUR_GRCh38.txt"
    run:
        with open(input[0],"r") as f_in, open(output[0],"w") as f_out:
            # First write the egyptians
            for egyptian in EGYPT_SAMPLES:
                # PD114, PD115, PD82 are from Upper Egypt
                if egyptian[:2] == "PD":
                    f_out.write(egyptian+"\tEGU\tEGY\n")
                # LU18, LU19, LU2, LU22. LU23, LU9, and Egyptref are from Delta
                else:
                    f_out.write(egyptian+"\tEGD\tEGY\n")
            # Then the 1000G samples
            for line in f_in:
                s = line.split("\t")
                pop = s[6]
                if pop in POPULATIONS_1000G and s[13] == "1":
                    if pop in POPULATIONS_AFR:
                        f_out.write(s[1]+"\t"+pop+"\tAFR\n")
                    elif pop in POPULATIONS_EUR:
                        f_out.write(s[1]+"\t"+pop+"\tEUR\n")

# Selecting from the VCF files those individuals that are to be used
# Keeping only variants with at last 5% MAF
# Keeping only variants not violating Hardy-Weinberg-Equilibrium
# Keeping only bi-allelic variants (min-allele = max-allele = 2)
rule gp_select_1000g_individual_genotypes:
    input: "1000_genomes/ALL.chr{chr}_GRCh38.genotypes.20170504.vcf.gz",
           "genotype_pcs/keep_indiv.txt"
    output: "genotype_pcs/AFR_EUR.chr{chr}_GRCh38.vcf.gz"
    params: log_base=lambda wildcards, output: output[0][:-7]
    conda: "envs/genotype_pcs.yaml"
    shell: "vcftools --gzvcf {input[0]} " + \
                    "--keep {input[1]} " + \
                    "--min-alleles 2 " + \
                    "--max-alleles 2 " + \
                    "--maf 0.05 " + \
                    "--hwe 0.000001 " + \
                    "--recode-INFO-all " + \
                    "--recode " + \
                    "--out {params.log_base} " + \
                    "--stdout | bgzip > {output[0]}"

# Compressing and indexing of files to be used with vcf-merge
rule gp_index_1000g:
    input: "genotype_pcs/AFR_EUR.chr{chr}_GRCh38.vcf.gz"
    output: "genotype_pcs/AFR_EUR.chr{chr}_GRCh38.vcf.gz.tbi"
    conda: "envs/genotype_pcs.yaml"
    shell: "tabix -p vcf {input}"

# Getting the list of SNPs for genotype PCs from the 1000 Genomes samples
rule gp_get_1000g_snps:
    input: "genotype_pcs/AFR_EUR.chr{chr}_GRCh38.vcf.gz"
    output: "genotype_pcs/snps_chr{chr}.txt"
    shell: "zcat {input} | grep -v '#' | cut -f 1,2 > {output}"

# Here, we select from the SNPs called for the egyptians those, which are also
# kept from the 1000 genomes samples, i.e. 5% MAF, HWE, bi-allelic
rule gp_select_matching_egyptian_snps:
    input: "variants_GRCh38/egyptians.chromosome.{chr}.vcf.gz",
           "genotype_pcs/snps_chr{chr}.txt"
    output: "genotype_pcs/egyptians.chromosome.{chr}.vcf.gz"
    params: log_base=lambda wildcards, output: output[0][:-7]
    conda: "envs/genotype_pcs.yaml"
    shell: "vcftools --gzvcf {input[0]} " + \
                    "--positions {input[1]} " + \
                    "--recode-INFO-all " + \
                    "--recode " + \
                    "--out {params.log_base} " + \
                    "--stdout | bgzip > {output[0]}"

# Compressing and indexing of files to be used with vcf-merge
rule gp_index_egyptians:
    input: "genotype_pcs/egyptians.chromosome.{chr}.vcf.gz"
    output: "genotype_pcs/egyptians.chromosome.{chr}.vcf.gz.tbi"
    conda: "envs/genotype_pcs.yaml"
    shell: "tabix -p vcf {input}"

# Merging the vcf-files of 1000 genomes with our SNP calls for the egyptians
rule gp_merge_1000g_with_egyptians:
    input: "genotype_pcs/egyptians.chromosome.{chr}.vcf.gz",
           "genotype_pcs/egyptians.chromosome.{chr}.vcf.gz.tbi",
           "genotype_pcs/AFR_EUR.chr{chr}_GRCh38.vcf.gz",
           "genotype_pcs/AFR_EUR.chr{chr}_GRCh38.vcf.gz.tbi"
    output: "genotype_pcs/EGYPT_AFR_EUR.chr{chr}_GRCh38.vcf.gz"
    conda: "envs/genotype_pcs.yaml"
    shell: "vcf-merge {input[0]} {input[2]} | bgzip > {output[0]}"

rule gp_merge_1000g_with_egyptians_all:
    input: expand("genotype_pcs/EGYPT_AFR_EUR.chr{chr}_GRCh38.vcf.gz", \
                   chr=[str(x) for x in range(1,23)]+["X","Y"])

# Concatenate the vcf file from several chromosomes
# --pad-missing: Write '.' in place of missing columns. Useful for joining chrY 
# with the rest.
rule gp_concatenate_chr_vcfs:
    input: expand("genotype_pcs/EGYPT_AFR_EUR.chr{chr}_GRCh38.vcf.gz", \
                   chr=[str(x) for x in range(1,23)]+["X","Y"])
    output: "genotype_pcs/EGYPT_AFR_EUR_GRCh38.vcf.gz"
    conda: "envs/genotype_pcs.yaml"
    shell: "vcf-concat --pad-missing {input} | bgzip > {output}"

# VCF file for the egyptians only (all chromosomes) here gets annotated with 
# the rsids from dbsnp
rule gp_cp_egyptian_vcfs_and_annotate_rsids:
    input: vcf="variants_GRCh38/egyptians.vcf.gz",
           vcf_index="variants_GRCh38/egyptians.vcf.gz.tbi",
           dbsnp="dbsnp_GRCh38/dbsnp.vcf.gz",
           dbsnp_index="dbsnp_GRCh38/dbsnp.vcf.gz.tbi"
    output: "genotype_pcs/egyptians.vcf.gz"
    conda: "envs/genotype_pcs.yaml"
    shell: "bcftools annotate --annotations {input.dbsnp} " + \
                             "--columns ID " + \
                             "--output {output} " + \
                             "--output-type z " + \
                             "{input.vcf} "

# Filtering egyptian only variants
# Additional to maf and number of alleles, we also exclude all SNPs with 
# missing data here, since we have only 10 individuals
# --max-missing <float>: Exclude sites on the basis of the proportion of missing
#                        data (defined to be between 0 and 1, where 0 allows 
#                        sites that are completely missing and 1 indicates no 
#                        missing data allowed).
# Further, we select the set of SNPs after LD pruning using 1000 genomes data
# The reason is that LD pruning on only the 10 individuals will not select
# truly independent SNPs because the 10 individuals are not enough.
rule gp_filter_for_egyptian_only_pcs:
    input: "genotype_pcs/egyptians.vcf.gz"
    output: "genotype_pcs/EGYPT_GRCh38.vcf.gz"
    params: log_base=lambda wildcards, output: output[0][:-7]
    conda: "envs/genotype_pcs.yaml"
    shell: "vcftools --gzvcf {input[0]} " + \
                    "--min-alleles 2 " + \
                    "--max-alleles 2 " + \
                    "--maf 0.05 " + \
                    "--hwe 0.000001 " + \
                    "--max-missing 1 " + \
                    "--recode-INFO-all " + \
                    "--recode " + \
                    "--out {params.log_base} " + \
                    "--stdout | bgzip > {output[0]}"

# Converting vcf files to plink binary format (bed/bim/fam) for preparing for
# Eigenstrat analysis
rule gp_vcf_to_plink:
    input: "genotype_pcs/{set}_GRCh38.vcf.gz"
    output: "genotype_pcs/plink/{set}_GRCh38.bed",
            "genotype_pcs/plink/{set}_GRCh38.bim",
            "genotype_pcs/plink/{set}_GRCh38.fam"
    params: out_base=lambda wildcards, output: output[0][:-4]
    conda: "envs/genotype_pcs.yaml"
    shell: "plink2 --vcf {input} " + \
                  "--make-bed " + \
                  "--out {params.out_base}"

# Removal of regions of high LD and/or known inversions from Abraham 2014, 
# i.e. Fellay 2009:
# chr6:25 Mb–33.5 Mb, (see also Wang 2009)
# chr5:44 Mb–51.5 Mb, chr8:8 Mb–12 Mb, chr11:45 Mb–57 Mb
# Therefore, make lists of SNPs in the respective regions to be removed,
# Then: Concatenate all the SNPs to be removed
rule gp_find_snps_from_high_ld_regions:
    input: "genotype_pcs/plink/{set}_GRCh38.bed", 
           "genotype_pcs/plink/{set}_GRCh38.bim",
           "genotype_pcs/plink/{set}_GRCh38.fam"
    output: "genotype_pcs/plink/{set}_GRCh38_6_25-33.5.snplist",
            "genotype_pcs/plink/{set}_GRCh38_5_44-51.5.snplist",
            "genotype_pcs/plink/{set}_GRCh38_8_8-12.snplist",
            "genotype_pcs/plink/{set}_GRCh38_11_45-57.snplist",
            "genotype_pcs/plink/{set}_GRCh38_exclusion.snplist" 
    params: in_base = lambda wildcards, input: input[0][:-4],
            chr6_base = lambda wildcards, output: output[0][:-8],
            chr5_base = lambda wildcards, output: output[1][:-8],
            chr8_base = lambda wildcards, output: output[2][:-8],
            chr11_base = lambda wildcards, output: output[3][:-8]
    conda: "envs/genotype_pcs.yaml"
    shell: 
        "plink2 --bfile {params.in_base} " + \ 
               "--chr 6 " + \
               "--from-mb 25 " + \
               "--to-mb 33.5 " + \
               "--write-snplist " + \
               "--out {params.chr6_base}; " + \
        "plink2 --bfile {params.in_base} " + \ 
               "--chr 5 " + \
               "--from-mb 44 " + \
               "--to-mb 51.5 " + \
               "--write-snplist " + \
               "--out {params.chr5_base}; " + \
        "plink2 --bfile {params.in_base} " + \ 
               "--chr 8 " + \
               "--from-mb 8 " + \
               "--to-mb 12 " + \
               "--write-snplist " + \
               "--out {params.chr8_base}; " + \
        "plink2 --bfile {params.in_base} " + \ 
               "--chr 11 " + \
               "--from-mb 45 " + \
               "--to-mb 57 " + \
               "--write-snplist " + \
               "--out {params.chr11_base}; " + \
        "cat {output[0]} {output[1]} {output[2]} {output[3]}  > {output[4]} "

# Now exclude the SNPs from these regions
rule gp_exclude_snps_from_high_ld_regions:
    input: "genotype_pcs/plink/{set}_GRCh38.bed", 
           "genotype_pcs/plink/{set}_GRCh38.bim",
           "genotype_pcs/plink/{set}_GRCh38.fam",
           "genotype_pcs/plink/{set}_GRCh38_exclusion.snplist"
    output: "genotype_pcs/plink/{set}_GRCh38_wo_ldregions.bed", 
            "genotype_pcs/plink/{set}_GRCh38_wo_ldregions.bim",
            "genotype_pcs/plink/{set}_GRCh38_wo_ldregions.fam"
    params: in_base = lambda wildcards, input: input[0][:-4],
            out_base = lambda wildcards, output: output[0][:-4]
    conda: "envs/genotype_pcs.yaml"
    shell: "plink2 --bfile {params.in_base} " + \
                  "--exclude {input[3]} " + \
                  "--make-bed " + \
                  "--out {params.out_base} "

# LD prune the PLINK files; therefore, first make a list of SNPs in LD (and not 
# in LD)(i.e. to be removed or not)
# Parameters for indep-pairwise: [window size]<kb> [step size (variant ct)] 
# [VIF threshold]
# Explanation Plink website): the command above that specifies 50 5 0.5 would 
# a) consider a window of 50 SNPs, 
# b) calculate LD between each pair of SNPs in the window, 
# c) remove one of a pair of SNPs if the LD is greater than 0.5, 
# d) shift the window 5 SNPs forward and repeat the procedure
# Abraham 2014 used: 1000 10 0.02
# Anderson 2010 used: 50 5 0.2
# Wang 2009 used: 100 ? 0.2
# Fellay 2009 used: 1500 150 0.2 
# Watch out: the LD pruned SNPs (to be kept, to be pruned out) are always taken
# from the entire set of Egyptian, European and African samples we have; there
# we use as LD pruning parameter: 1000 10 0.2
rule gp_find_ld_pruned_snps:
    input: "genotype_pcs/plink/EGYPT_AFR_EUR_GRCh38_wo_ldregions.bed", 
           "genotype_pcs/plink/EGYPT_AFR_EUR_GRCh38_wo_ldregions.bim",
           "genotype_pcs/plink/EGYPT_AFR_EUR_GRCh38_wo_ldregions.fam"
    output: "genotype_pcs/plink/EGYPT_AFR_EUR_GRCh38_wo_ldregions.prune.in", 
            "genotype_pcs/plink/EGYPT_AFR_EUR_GRCh38_wo_ldregions.prune.out"
    params: in_base = lambda wildcards, input: input[0][:-4]
    conda: "envs/genotype_pcs.yaml"
    shell: "plink2 --bfile {params.in_base} " + \
                  "--indep-pairwise 1000 10 0.2 " + \
                  "--out {params.in_base} "

# Now exclude the pruned SNPs
# We use the same set of pruned SNPs, those from the entire Egypt/Eur/Afr 
# samples also for other subsets, because the LD pruning will not work well
# for small numbers of samples
# "--exclude {input[3]} " + \ in case the file is 
# genotype_pcs/plink/EGYPT_AFR_EUR_GRCh38_wo_ldregions.prune.out
# For the merged 1000g and egyptian data there should be no difference,
# but for the egyptian only data there is, because we only want to keep IDs also
# in the 1000 G data (and no rare variants, even if they are not rare in 
# egyptians)
rule gp_exclude_ld_pruned_snps:
    input: "genotype_pcs/plink/{set}_GRCh38_wo_ldregions.bed", 
            "genotype_pcs/plink/{set}_GRCh38_wo_ldregions.bim",
            "genotype_pcs/plink/{set}_GRCh38_wo_ldregions.fam",
            "genotype_pcs/plink/EGYPT_AFR_EUR_GRCh38_wo_ldregions.prune.in"
    output: "genotype_pcs/plink/{set}_GRCh38_wo_ldregions_pruned.bed", 
            "genotype_pcs/plink/{set}_GRCh38_wo_ldregions_pruned.bim",
            "genotype_pcs/plink/{set}_GRCh38_wo_ldregions_pruned.fam"
    params: in_base = lambda wildcards, input: input[0][:-4],
            out_base = lambda wildcards, output: output[0][:-4]
    conda: "envs/genotype_pcs.yaml"
    shell: "plink2 --bfile {params.in_base} " + \
                  "--extract {input[3]} " + \
                  "--make-bed " + \
                  "--out {params.out_base}"

# Conversion from bed/bim/fam to ped/map
rule gp_convert_to_ped_map:
    input: "genotype_pcs/plink/{set}_GRCh38_wo_ldregions_pruned.bed", 
           "genotype_pcs/plink/{set}_GRCh38_wo_ldregions_pruned.bim",
           "genotype_pcs/plink/{set}_GRCh38_wo_ldregions_pruned.fam"
    output: "genotype_pcs/plink/{set}_GRCh38_wo_ldregions_pruned.ped", 
            "genotype_pcs/plink/{set}_GRCh38_wo_ldregions_pruned.map"
    params: in_base = lambda wildcards, input: input[0][:-4]
    conda: "envs/genotype_pcs.yaml"
    shell: "plink2 --bfile {params.in_base} " + \
                  "--recode " + \
                  "--out {params.in_base} "

# Write the parameter file needed by the Eigensoft convertf program
rule gp_eigentstrat_parameter_file:
    input: "genotype_pcs/plink/{set}_GRCh38_wo_ldregions_pruned.ped", 
           "genotype_pcs/plink/{set}_GRCh38_wo_ldregions_pruned.map"
    output: "genotype_pcs/eigenstrat/{set}_GRCh38.ped2eigenstrat.params",
    params: gout="genotype_pcs/eigenstrat/{set}_GRCh38.eigenstratgeno",
            sout="genotype_pcs/eigenstrat/{set}_GRCh38.snp",
            iout="genotype_pcs/eigenstrat/{set}_GRCh38.ind"
    run: 
        with open(output[0],"w") as f_out:
            f_out.write("genotypename:    "+input[0]+"\n")
            f_out.write("snpname:         "+input[1]+"\n") 
            f_out.write("indivname:       "+input[0]+"\n")
            f_out.write("outputformat:    EIGENSTRAT\n")
            f_out.write("genotypeoutname: "+params.gout+"\n")
            f_out.write("snpoutname:      "+params.sout+"\n")
            f_out.write("indivoutname:    "+params.iout+"\n")
            f_out.write("familynames:     NO\n")

# This is the actual conversion from ped format to the eigenstrat input format
rule gp_ped_to_eigentstrat:
    input: "genotype_pcs/plink/{set}_GRCh38_wo_ldregions_pruned.ped", 
           "genotype_pcs/plink/{set}_GRCh38_wo_ldregions_pruned.map",
           "genotype_pcs/eigenstrat/{set}_GRCh38.ped2eigenstrat.params"
    output: "genotype_pcs/eigenstrat/{set}_GRCh38.eigenstratgeno",
            "genotype_pcs/eigenstrat/{set}_GRCh38.snp",
            "genotype_pcs/eigenstrat/{set}_GRCh38.ind"
    conda: "envs/genotype_pcs.yaml"
    shell: "convertf -p {input[2]}"

# Running Eigensofts smartpca module which computes the population PCs
# The smartpca parameters:
# -i example.geno  : genotype file in any format (see ../CONVERTF/README)
# -a example.snp   : snp file in any format (see ../CONVERTF/README)
# -b example.ind   : indiv file in any format (see ../CONVERTF/README)
# -k k             : (Default is 10) number of principal components to output
# -o example.pca   : output file of principal components.  Individuals removed
#                    as outliers will have all values set to 0.0 in this file.
# -p example.plot  : prefix of output plot files of top 2 principal components.
#                    (labeling individuals according to labels in indiv file)
# -e example.eval  : output file of all eigenvalues
# -l example.log   : output logfile
# -m maxiter       : (Default is 5) maximum number of outlier removal iterations.
#                    To turn off outlier removal, set -m 0.
# -t topk          : (Default is 10) number of principal components along which 
#                    to remove outliers during each outlier removal iteration.
# -s sigma         : (Default is 6.0) number of standard deviations which an
#                    individual must exceed, along one of topk top principal
# 		             components, in order to be removed as an outlier.
# Watch out: The clusters smartpca.perl scripts is used, but calling the conda
# smartpca. ploteig seems to not be included in conda, and again is used by 
# from the cluster. This is not working properly though and not producing the 
# PDF plot. Probably it would be sufficient to use the smartpca binary.
rule gp_eigensoft_smartpca:
    input: "genotype_pcs/eigenstrat/{set}_GRCh38.eigenstratgeno",
           "genotype_pcs/eigenstrat/{set}_GRCh38.snp",
           "genotype_pcs/eigenstrat/{set}_GRCh38.ind"
    output: "genotype_pcs/eigenstrat/{set}_GRCh38.pca",
            "genotype_pcs/eigenstrat/{set}_GRCh38.eval",
            "genotype_pcs/eigenstrat/{set}_GRCh38.log",
            "genotype_pcs/eigenstrat/{set}_GRCh38.pca.evec",
    params: out_base = lambda wildcards, output: output[0][:-4]
    conda: "envs/genotype_pcs.yaml"
    shell: "smartpca.perl -i {input[0]} " + \
                         "-a {input[1]} " + \
                         "-b {input[2]} " + \
                         "-o {output[0]} " + \
                         "-p {params.out_base}.plot " + \
                         "-e {output[1]} " + \
                         "-l {output[2]} " + \
                         "-m 0; "

# Computing the Tracy-Widom statistics to evaluate the statistical 
# significance of each principal component identified by pca
rule gp_tracy_widom_pval:
    input: "genotype_pcs/eigenstrat/{set}_GRCh38.eval",
           "data/misc/twtable"
    output: "genotype_pcs/eigenstrat/{set}_GRCh38.tw"
    conda: "envs/genotype_pcs.yaml"
    shell: "twstats -i {input[0]} " + \
                   "-t {input[1]} " + \
                   "-o {output[0]} "

# Plotting the PCs
rule gp_plot_gt_pcs:
    input: "genotype_pcs/eigenstrat/EGYPT_AFR_EUR_GRCh38.pca.evec",
           "genotype_pcs/annotation_EGYPT_AFR_EUR_GRCh38.txt"
    output: "genotype_pcs/figures/EGYPT_AFR_EUR_GRCh38_pca_1vs2.pdf",
            "genotype_pcs/figures/EGYPT_AFR_EUR_GRCh38_pca_1vs3.pdf",
            "genotype_pcs/figures/EGYPT_AFR_EUR_GRCh38_pca_1vs4.pdf",
            "genotype_pcs/figures/EGYPT_AFR_EUR_GRCh38_pca_2vs3.pdf",
            "genotype_pcs/figures/EGYPT_AFR_EUR_GRCh38_pca_2vs4.pdf",
            "genotype_pcs/figures/EGYPT_AFR_EUR_GRCh38_pca_3vs4.pdf",
            "genotype_pcs/figures/EGYPT_AFR_EUR_GRCh38_scree_plot.pdf",
            "genotype_pcs/figures/EGYPT_AFR_EUR_GRCh38_pca_3d.pdf"
    params: out_path = "genotype_pcs/figures/"
    conda: "envs/genotype_pcs.yaml"
    script: "scripts/plot_gt_pcs.R"

rule gp_plot_gt_pcs_egyptians:
    input: "genotype_pcs/eigenstrat/EGYPT_GRCh38.pca.evec",
           "genotype_pcs/annotation_EGYPT_AFR_EUR_GRCh38.txt"
    output: "genotype_pcs/figures/EGYPT_GRCh38_pca_1vs2.pdf",
            "genotype_pcs/figures/EGYPT_GRCh38_pca_1vs3.pdf",
            "genotype_pcs/figures/EGYPT_GRCh38_pca_1vs4.pdf",
            "genotype_pcs/figures/EGYPT_GRCh38_pca_2vs3.pdf",
            "genotype_pcs/figures/EGYPT_GRCh38_pca_2vs4.pdf",
            "genotype_pcs/figures/EGYPT_GRCh38_pca_3vs4.pdf",
            "genotype_pcs/figures/EGYPT_GRCh38_scree_plot.pdf",
            "genotype_pcs/figures/EGYPT_GRCh38_pca_3d.pdf"
    params: out_path = "genotype_pcs/figures/"
    conda: "envs/genotype_pcs.yaml"
    script: "scripts/plot_gt_pcs_egyptians.R"


################################################################################
##### Computing depth for Novogene file of GRCh38 reads unmapped to assembly ###
################################################################################

# Symlink file
rule ng_symlink_unmapped_reads:
    output: "novogene/all.new1.bam"
    shell: "ln -s /data/lied_egypt_genome/raw/" + \
           "P101HW18010820-01_human_2018.10.11_align/all.new1.bam " + \
           "{output}"

# Sort the reads by coordinate, needed for computing coverage
rule ng_sort_bam:
    input: "novogene/all.new1.bam"
    output: "novogene/all.new1.sorted.bam"
    shell: "samtools sort --threads 24 " + \
                         "--output-fmt BAM " + \
                         "{input} > {output}"

# Getting coverage
rule ng_coverage:
    input: "novogene/all.new1.sorted.bam"
    output: "novogene/coverage.txt"
    shell: "samtools depth {input} > {output}"

rule ng_occ_vs_dp:
    input: "novogene/coverage.txt"
    output: "novogene/coverage_occurrence_vs_depth.txt"
    shell: "cat {input} | cut -f 3 | sort | uniq -c > {output}"

rule ng_get_bases_not_in_assembly:
    input: "novogene/coverage_occurrence_vs_depth.txt"
    output: "novogene/uncovered_bases.txt"
    run:
        with open(input[0],"r") as f_in, open(output[0],"w") as f_out:
            sum = 0
            for line in f_in:
                s = line.strip(" ").split(" ")
                if s[1] == "0":
                    continue
                sum += int(s[0])
            f_out.write("UNCOVERED BASES: "+str(sum)+"\n")

rule ng_plot_uncovered_bases_distribution:
    input: "novogene/coverage_occurrence_vs_depth.txt"
    output: "novogene/dist_uncovered_dp.pdf"
    script: "scripts/dist_uncovered_dp.R"


################################################################################
## Extracting various information within specified genes (i.e. gene-centric) ###
################################################################################

# These are the genes of interest (ABCC7=CFTR, HD=HDDC3?, Factor V=F5)
GENES = ["CFTR","HDDC3","DMD","BRCA1","BRCA2","TP53","EGFR","APP","PSEN1","F5", \
         "CARD11","LAMA4","MRC1","USH2A"]

# Therefore, obtain a recent Ensembl annotation file first
rule get_ensembl_gene_annotation_gtf:
    output: temp("annotations/Homo_sapiens.GRCh38.94.gtf.gz")
    shell: "wget -P annotations " + \
           "ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/" + \
           "Homo_sapiens.GRCh38.94.gtf.gz "

rule unzip_ensembl_gene_annotation_gtf:
    input: "annotations/Homo_sapiens.GRCh38.94.gtf.gz"
    output: "annotations/Homo_sapiens.GRCh38.94.gtf"
    shell: "gzip -d {input}"

rule gc_get_gene_annotation:
    input: "annotations/Homo_sapiens.GRCh38.94.gtf"
    output: "gene_centric/{gene}/{gene}.gtf"
    shell: "cat {input} | grep '#' > {output}; " + \
           "cat {input} | grep 'gene_name \"{wildcards.gene}\";' >> {output}"

# How many bases left and right of gene boundaries to consider 
WINDOW = {
    "CFTR": [100000,100000],
    "HDDC3": [100000,100000],
    "DMD": [100000,100000],
    "BRCA1": [100000,100000],
    "BRCA2": [100000,100000],
    "TP53": [100000,100000],
    "EGFR": [100000,100000],
    "APP": [100000,100000],
    "PSEN1": [100000,100000],
    "F5": [100000,100000],
    "CARD11": [100000,100000],
    "LAMA4": [100000,100000],
    "MRC1": [100000,100000],
    "USH2A": [100000,100000]
}
rule gc_get_start_end_position:
    input: "gene_centric/{gene}/{gene}.gtf"
    output: "gene_centric/{gene}/{gene}.bed"
    run:
        with open(input[0],"r") as f_in, open(output[0],"w") as f_out:
            f_out.write("# Custom bed file for region around gene\n")
            for line in f_in:
                if line[0] == "#":
                    continue
                s = line.split("\t")
                if s[2] == "gene":
                    chr = s[0]
                    start = str(int(s[3]) - WINDOW[wildcards.gene][0])
                    end = str(int(s[4]) + WINDOW[wildcards.gene][1])
                    strand = s[6]
                    f_out.write("\t".join([chr,start,end,'.','.',strand])+"\n")

rule gc_get_overlapping_genes:
    input: "annotations/Homo_sapiens.GRCh38.94.gtf",
           "gene_centric/{gene}/{gene}.bed"
    output: "gene_centric/{gene}/{gene}_overlapping.gtf"
    run:
        with open(input[1],"r") as f_in:
            for line in f_in:
                if line[0] == "#":
                    continue
                s = line.split('\t')
                [q_chrom,q_start,q_end] = s[:3]
        with open(input[0],"r") as f_in, open(output[0],"w") as f_out:
                for line in f_in:
                    if line[0] == "#":
                        continue
                    s = line.split("\t")
                    chrom,start,end = s[:3]
                    if chrom == q_chrom:
                        if q_start<start<q_end or q_start<end<q_end:
                            f_out.write(line)

# Selecting the Illumina reads of the Egyptians mapped to this gene 
# -b: output bam
# -L FILE: only include reads overlapping this BED FILE
rule gc_get_mapped_egyptref_reads:
    input: bam="variants_GRCh38/{sample}.final.bam",
           bed="gene_centric/{gene}/{gene}.bed"
    output: "gene_centric/{gene}/{sample}.{gene}.bam",
            "gene_centric/{gene}/{sample}.{gene}.bam.bai",
    shell: "samtools view -b -L {input.bed} {input.bam} > {output[0]}; " + \
           "samtools index {output[0]} "

rule gc_get_mapped_egyptref_reads_all:
    input: expand("gene_centric/{gene}/{sample}.{gene}.bam", \
                   sample=EGYPT_SAMPLES, \
                   gene=GENES),
          expand("gene_centric/{gene}/{gene}_overlapping.gtf", gene=GENES)

# Get the SNP calls of the Egyptians for the specified genes
rule gc_get_variants:
    input: vcf="variants_GRCh38/egyptians.vcf.gz",
           bed="gene_centric/{gene}/{gene}.bed"
    output: "gene_centric/{gene}/{gene}_egyptians.vcf.gz"
    shell: "vcftools --gzvcf {input.vcf} " + \
                    "--bed {input.bed} " + \
                    "--recode " + \
                    "--recode-INFO-all " + \
                    "--stdout " + \
                    "| bgzip > {output}"

rule gc_get_dbsnp_variants:
    input: vcf="dbsnp_GRCh38/dbsnp.vcf.gz",
           bed="gene_centric/{gene}/{gene}.bed"
    output: "gene_centric/{gene}/{gene}_dbsnp.vcf.gz"
    conda: "envs/genotype_pcs.yaml"
    shell: "vcftools --gzvcf {input.vcf} " + \
                    "--bed {input.bed} " + \
                    "--recode " + \
                    "--recode-INFO-all " + \
                    "--stdout " + \
                    "| bgzip > {output}"

# Combining the chromosome-wise vcfs and removing trailing "chr"
rule gc_combine_chromosomewise_1000g_vcf:
    input: expand("1000_genomes/ALL.chr{chr}_GRCh38.genotypes.20170504.vcf.gz", \
                   chr=[str(x) for x in range(1,23)]+["X","Y"])
    output: "1000_genomes/ALL.GRCh38.genotypes.20170504.vcf.gz"
    conda: "envs/genotype_pcs.yaml"
    shell: "vcf-concat --pad-missing {input} | bgzip > {output}"

# needed?
#rule index_1000g_vcf:
#    input: "1000_genomes/ALL.GRCh38.genotypes.20170504.vcf.gz"
#    output: "1000_genomes/ALL.GRCh38.genotypes.20170504.vcf.gz.tbi"
#    shell: "tabix -p vcf {input}"

rule gc_get_1000g_variants:
    input: vcf="1000_genomes/ALL.GRCh38.genotypes.20170504.vcf.gz",
           bed="gene_centric/{gene}/{gene}.bed"
    output: "gene_centric/{gene}/{gene}_1000g.vcf.gz"
    conda: "envs/genotype_pcs.yaml"
    shell: "vcftools --gzvcf {input.vcf} " + \
                    "--bed {input.bed} " + \
                    "--recode " + \
                    "--recode-INFO-all " + \
                    "--stdout " + \
                    "| bgzip > {output}"

rule gc_get_variants_all:
    input: expand("gene_centric/{gene}/{gene}_egyptians.vcf.gz",gene=GENES),
           expand("gene_centric/{gene}/{gene}_dbsnp.vcf.gz",gene=GENES),
           expand("gene_centric/{gene}/{gene}_1000g.vcf.gz",gene=GENES)

# Preprocessing and filtering the Annovar annotated files from Axel
# Here, first filter genes that are within the specified gene, are exonic and 
# variants with CADD score >=20
rule gc_filter_annotated_files:
    input: "data/annovar/{gene}.{sample}.avinput.hg38_multianno.txt"
    output: "gene_centric/{gene}/{gene}.{sample}_annovar_filtered.txt"
    run:
        with open (input[0],"r") as f_in, open (output[0],"w") as f_out:
            for line in f_in:
                s = line.split("\t")
                if line[:3] == "Chr":
                    f_out.write(line)
                    continue
                gene = s[6]
                cadd_phred = s[69]
                if cadd_phred in ["NA","."]:
                    continue
                function = s[5]
                if gene == wildcards.gene and \
                   float(cadd_phred) >= 20 and \
                   function == "exonic":
                   f_out.write(line)

rule gc_concat_filter_annotated_files:
    input: expand("gene_centric/{gene}/{gene}.{sample}_annovar_filtered.txt", \
                   gene=GENES, sample=EGYPT_SAMPLES)
    output: "gene_centric/annovar_filtered.txt"
    run:
        out_lines = []
        i = 0
        for filename in input:
            with open (filename,"r") as f_in:
                for line in f_in:
                    if line[:3] == "Chr":
                        header = "Sample\t"+line
                        continue
                    sample = filename.split(".")[1].split("_")[0]
                    s = line.split("\t")
                    chrom = s[0]
                    pos = s[1]
                    out_lines.append([sample]+s)
            # Sort by sample tertiary
            out_lines = sorted(out_lines, key=lambda x: x[0])
            # Sort py position secondary
            out_lines = sorted(out_lines,key=lambda x: x[2])
            # Sort by chromosome primarily
            out_lines = sorted(out_lines,key=lambda x: x[1])
            with open (output[0],"w") as f_out:
                f_out.write(header)
                for elem in out_lines:
                    f_out.write("\t".join(elem))

# Preprocessing and filtering the Annovar annotated files from Axel, here all
# variants have been annotated.
# Here, first filter variants that are exonic and rare according to no (="NA")
# or according to low gnomAD_exome_ALL
rule av_filter_egyptian_variants:
    input: "data/annovar/sample.{sample}.avinput.hg38_multianno.txt"
    output: "annovar/sample.{sample}_annovar_egyptian.txt"
    run:
        with open (input[0],"r") as f_in, open (output[0],"w") as f_out:
            for line in f_in:
                if line[:3] == "Chr":
                    f_out.write(line)
                    continue
                s = line.split("\t")
                function = s[5]
                gnomAD_exome = s[11]
                if function == "exonic":
                    if gnomAD_exome in ["NA","."]:
                        f_out.write(line)
                    elif float(gnomAD_exome) <= 0.001:
                        f_out.write(line)

# This is without applying any filtering
rule av_nofilter_variants:
    input: "data/annovar/sample.{sample}.avinput.hg38_multianno.txt"
    output: "annovar/sample.{sample}_annovar_nofilter.txt"
    shell: "cp {input} {output}"

rule av_nofilter_variants_all:
    input: expand("annovar/sample.{sample}_annovar_nofilter.txt", \
                  sample=EGYPT_SAMPLES)

# Filter variants with CADD score greater 20
rule av_filter_deleterious:
    input: "annovar/sample.{sample}_annovar_egyptian.txt"
    output: "annovar/sample.{sample}_annovar_egyptiandeleterious.txt" 
    run:
        with open (input[0],"r") as f_in, open (output[0],"w") as f_out:
            for line in f_in:
                if line[:3] == "Chr":
                    f_out.write(line)
                    continue
                s = line.split("\t")
                cadd_phred = s[69]
                if cadd_phred in ["NA","."]:
                    continue
                if float(cadd_phred) >= 20:
                        f_out.write(line)

# Concatenate the file for all samples, add a 'sample' column and sort by 
# position
rule av_concat_files:
    input: expand("annovar/sample.{sample}_annovar_{{filter}}.txt", \
                   sample=EGYPT_SAMPLES)
    output: "annovar/annovar_{filter}.txt"
    run:
        out_lines = []
        i = 0
        for filename in input:
            with open (filename,"r") as f_in:
                for line in f_in:
                    if line[:3] == "Chr":
                        header = "Sample\t"+line
                        continue
                    sample = filename.split(".")[1].split("_")[0]
                    s = line.split("\t")
                    chrom = s[0]
                    pos = s[1]
                    out_lines.append([sample]+s)
            # Sort by sample tertiary
            out_lines = sorted(out_lines, key=lambda x: x[0])
            # Sort py position secondary
            out_lines = sorted(out_lines,key=lambda x: int(x[2]))
            # Sort by chromosome primarily
            out_lines = sorted(out_lines,key=lambda x: x[1])
            with open (output[0],"w") as f_out:
                f_out.write(header)
                for elem in out_lines:
                    f_out.write("\t".join(elem))

# Filter (deleterious) variants to keep only those in more than one Egyptian
# Sort the list by chromosome (-k 2 -g -s, .g: numeric, -s: stable to keep 
# sorting by position) and then numerically by position (-k 3 -g)
rule av_deleterius_recurrent:
    input: "annovar/annovar_{filter}.txt"
    output: "annovar/annovar_{filter}_recurrent.txt"
    shell: "head -n 1 {input} > {output}; " + \
           "cat {input} | " + \
           "cut -f 2-100 | " + \
           "sort | " + \
           "uniq -c | " + \
           "grep -v '      1 ' | " + \
           "sed -e 's/      //g' | " + \
           "sed -e 's/     //g' | " + \
           "sed -e 's/ /\t/1' | " + \
           "sort -k 3 -g -t \"\t\" | " + \
           "sort -k 2 -g -t \"\t\" -s >> {output}"

rule av_recurrent_all:
    input: expand("annovar/annovar_{filter}_recurrent.txt", \
                  filter=["egyptian","egyptiandeleterious"])
