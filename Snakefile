#kate:syntax python;

#######################################
### Analyzing an Egyptian genome
#######################################

from Bio import SeqIO
import os


# Chromosome and scaffold names for later use
CHR_GRCh38 = ["chromosome."+str(x) for x in range(1,23)] \
           + ["chromosome."+str(x) for x in ["MT","X","Y"]]
EGYPT_SCAFFOLDS = ["fragScaff_scaffold_"+str(x)+"_pilon" for x in range(0,41)] \
                + ["original_scaffold_"+str(x)+"_pilon" for x in range(41,145)]

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
           "results/GRCh38/{task}_Homo_sapiens.GRCh38.dna.primary_assembly.txt", \
           task = ["scaffold_names","num_bases","num_all","assembly_stats"]),
           expand( \
           "results/EGYPTREF/{task}_Homo_sapiens.EGYPTREF.dna.primary_assembly.txt", \
           task = ["scaffold_names","num_bases","num_all","assembly_stats"])
                       
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
                  scaffolds=EGYPT_SCAFFOLDS),
           expand("busco_GRCh38/run_busco_GRCh38_{chrom}/short_summary_busco_GRCh38_{chrom}.txt", \
                  chrom=CHR_GRCh38)

# Make a comparison table for the busco analysis for EGYPTREF
rule summary_busco_egyptref:
    input: expand("busco_EGYPTREF/run_busco_EGYPTREF_{scaffolds}/full_table_busco_EGYPTREF_{scaffolds}.tsv", \
                  scaffolds=EGYPT_SCAFFOLDS)
    output: "busco_EGYPTREF/busco_summary.txt"
    script: "scripts/busco_summary.py"

# Make a comparison table for the busco analysis for GRCh38
rule summary_busco_grch38:
    input: expand("busco_GRCh38/run_busco_GRCh38_{chrom}/full_table_busco_GRCh38_{chrom}.tsv", \
                  chrom=CHR_GRCh38)
    output: "busco_GRCh38/busco_summary.txt"
    script: "scripts/busco_summary.py"

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

# Copy the assembled sequence
rule cp_and_rename_assembly:
    input: "data/pilon.fasta"
    output: "seq_EGYPTREF/Homo_sapiens.EGYPTREF.dna.primary_assembly.fa"
    shell: "cp {input} {output}"

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
rule run_repeatmasker_chromosomewise:
    input: expand("repeatmasked_GRCh38/Homo_sapiens.GRCh38.dna.{x}.fa.tbl", \
                  x=CHR_GRCh38),
           expand("repeatmasked_EGYPTREF/Homo_sapiens.EGYPTREF.dna.{x}.fa.tbl", \
                  x=EGYPT_SCAFFOLDS)

# Summarising the chromosome-wise repeatmasker summary files for Egyptref
rule repeatmasker_summary_table_egyptref:
    input: expand("repeatmasked_EGYPTREF/Homo_sapiens.EGYPTREF.dna.{x}.fa.tbl", \
                  x=EGYPT_SCAFFOLDS)
    output: "repeatmasked_EGYPTREF/summary.txt"
    script: "scripts/repeatmasker_summary.py"

# Summarising the chromosome-wise repeatmasker summary files for GRCh38
rule repeatmasker_summary_table_grch38:
    input: expand("repeatmasked_GRCh38/Homo_sapiens.GRCh38.dna.{x}.fa.tbl", \
                  x=CHR_GRCh38)
    output: "repeatmasked_GRCh38/summary.txt"
    script: "scripts/repeatmasker_summary.py"

# Making a repeatmasker stat table over all chromosomes, one line for EGYPTREF,
# one line for GRCh38
rule comparison_repeatmasker:
    input: expand("repeatmasked_{assembly}/summary.txt", \
                  assembly=["EGYPTREF","GRCh38"])
    output: "results/repeatmasker_comparison.txt"
    script: "scripts/repeatmasker_comparison.py"

# Writing the scaffolds of the Egyptian genome to separate fasta files because
# processing the whole assembly often takes too much time
rule write_scaffold_fastas:
    input: "data/pilon.fasta"
    output: expand("seq_EGYPTREF/Homo_sapiens.EGYPTREF.dna.{scaffold}.fa", \
                   scaffold=EGYPT_SCAFFOLDS)
    run:
        with open(input[0], "r") as f_in:
            i = 0
            for record in SeqIO.parse(f_in,"fasta"):            
                with open(output[i], "w") as f_out:
                    SeqIO.write(record, f_out, "fasta")
                    i += 1
            
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
           "repeatmasked_EGYPTREF/Homo_sapiens.EGYPTREF.dna.{scaffold}.fa.masked"
    output: "align_lastz_GRCh38_vs_EGYPTREF/{chr}_vs_{scaffold}.maf",
            "align_lastz_GRCh38_vs_EGYPTREF/dotplots/{chr}_vs_{scaffold}.rdotplot"
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
rule individual_lastz_dotplot:
    input: "align_lastz_GRCh38_vs_EGYPTREF/dotplots/{chr}_vs_{scaffold}.rdotplot"
    output: "align_lastz_GRCh38_vs_EGYPTREF/dotplots/{chr}_vs_{scaffold}.pdf"
    script: "scripts/dotplot.R"

# Plotting for one scaffold the dotplot versus all chromosomes
rule dotplots_scaffold_vs_chromosomes:
    input: expand("align_lastz_GRCh38_vs_EGYPTREF/dotplots/{chr}_vs_{{scaffold}}.rdotplot", \
                  chr=CHR_GRCh38)
    output: "align_lastz_GRCh38_vs_EGYPTREF/dotplots/{scaffold}.pdf"
    script: "scripts/scaffold_vs_grch38.R"            

# Plotting the dotplots for all scaffolds
rule dotplots_scaffold_vs_chromosomes_all:
    input: expand("align_lastz_GRCh38_vs_EGYPTREF/dotplots/{scaffold}.pdf", \
                  scaffold=EGYPT_SCAFFOLDS)

# All versus all comparisons of reference and Egyptian genome
rule align_all_vs_all:
    input: expand("align_lastz_GRCh38_vs_EGYPTREF/{chr}_vs_{scaffold}.maf", \
                  chr=CHR_GRCh38, scaffold=EGYPT_SCAFFOLDS)

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
                  chr=CHR_GRCh38, scaffold=EGYPT_SCAFFOLDS)

# All versus all dotplots of reference and Egyptian genome
rule all_vs_all_dotplots_mummer:
    input: expand("align_mummer_GRCh38_vs_EGYPTREF/dotplots/{chr}_vs_{scaffold}.gp", \
                  chr=CHR_GRCh38, scaffold=EGYPT_SCAFFOLDS)

# Plotting for one scaffold the dotplot versus all chromosomes
#rule mummer_dotplots_scaffold_vs_chromosomes:
#    input: expand("align_mummer_GRCh38_vs_EGYPTREF/dotplots/{chr}_vs_{{scaffold}}.ps", \
#                  chr=CHR_GRCh38)
#    output: "align_mummer_GRCh38_vs_EGYPTREF/dotplots/{scaffold}.ps",
#            "align_mummer_GRCh38_vs_EGYPTREF/dotplots/{scaffold}.pdf"
#    run: 
        # Putting 5 plots together horizontal; needs to be done 5 times
#        for h in range(5):
#            files_horizontal_merge = " ".join(input[h*5:h*2*h*5])
#            shell("convert "+files_horizontal_merge+" +append "+output[0]+str(h))
        # Putting the 5 horizontal plots together vertically
#        files_vertical_merge = " ".join([output[0]+str(h) for h in range(5)])
#        shell("convert "+files_vertical_merge+" -append "+output[0])
#        shell("convert {output[0]} {output[1]}")

# Plotting the dotplots for all scaffolds
rule mummer_dotplots_scaffold_vs_chromosomes_all:
    input: expand("align_lastz_GRCh38_vs_EGYPTREF/dotplots/{scaffold}.pdf", \
                  scaffold=EGYPT_SCAFFOLDS)

# Filtering the mummer alignments: Query sequences can be mapped to reference 
# sequences with -q, this allows the user to exclude chance and repeat 
# alignments, leaving only the best alignments between the two data sets (i.e.
# use the -q option for mapping query contigs to their best reference location)
# -u: float; Set the minimum alignment uniqueness, i.e. percent of the alignment 
#     matching to unique reference AND query sequence [0, 100], default 0
# -l: int; Set the minimum alignment length, default 0
rule delta_filter_mummer:
    input: "align_mummer_GRCh38_vs_EGYPTREF/{chr}_vs_{scaffold}.delta"
    output: "align_mummer_GRCh38_vs_EGYPTREF/{chr}_vs_{scaffold}.filter"
    conda: "envs/mummer.yaml"
    shell: "delta-filter -l 10000 -u 0 -q {input} > {output}"

# Computing the GRCh38 recovery rate using the mafTools package 
# (as in Cho et al.). Using mafTools program mafPairCoverage, it is necessary
# to first combine all chromosome/scaffold maf files, and then run 
# mafTransitiveClosure
rule combine_maf_files_for_recovery:
    input: expand("align_lastz_GRCh38_vs_EGYPTREF/{chr}_vs_{scaffold}.maf", \
                  chr=CHR_GRCh38, scaffold=EGYPT_SCAFFOLDS)
    output: "align_lastz_GRCh38_vs_EGYPTREF/recovery/all_alignments.maf"
    run: 
        shell("cat {input[0]} > {output}")
        for filename in input[1:]:
            # Append to large file; some file only have comments, no alignments
            # therefore we need to add & true because other wise the exit code
            # would indicate an error
            shell("cat {filename} | grep -v '#' >> {output} & true")

rule transitive_closure:
    input: "align_lastz_GRCh38_vs_EGYPTREF/recovery/all_alignments.maf"
    output: "align_lastz_GRCh38_vs_EGYPTREF/recovery/all_transclos.maf"
    shell: "./ext_tools/mafTools/bin/mafPairCoverage " + \
           "--maf --seq1 1 --seq2 * {input} > {output}"

rule maftools_coverage:
    input: "align_lastz_GRCh38_vs_EGYPTREF/recovery/all_transclos.maf"
    output: "align_lastz_GRCh38_vs_EGYPTREF/recovery/all.coverage"
    shell: "./ext_tools/mafTools/bin/mafPairCoverage " + \
           "--maf {input} > {output}"

rule recovery:
    input: "align_lastz_GRCh38_vs_EGYPTREF/recovery/all.coverage"
    output: "align_lastz_GRCh38_vs_EGYPTREF/recovery/recovery.txt"
    run:
        pass


# Processing Illumina PE data

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

# For SNP calling and other things that are done for the Illumnin PE data
rule symlink_illumina_wgs_dir:
    output: directory("data/02.DES")
    shell: "ln -s /data/lied_egypt_genome/raw/02.DES {output}"

# Some QC: Here, fastqc for all Illumina PE WGS files
rule run_fastqc:
    input: "data/02.DES/{lib}_{read}.fq.gz"
    output: html="illumina_qc/fastqc/{lib}_{read}_fastqc.html",
            zip="illumina_qc/fastqc/{lib}_{read}_fastqc.zip"
    conda: "envs/fastqc.yaml"
    shell: "fastqc --outdir illumina_qc/fastqc/ {input[0]}"

rule run_fastqc_all:
    input: expand("illumina_qc/fastqc/{lib}_{read}_fastqc.html", lib=ILLUMINA_LIBS, \
                                                          read=["1","2"])

################################################################################
################### SNP Calling for 9 Egyptian individuals #####################
################################################################################

# If possible, the variant calling (vc) tasks are performed with 8 threads and
# with 35Gb of memory, such that 5 tasks can be run on a node in parallel

# These are the sample IDs
EGYPT_SAMPLES = ["LU18","LU19","LU2","LU22","LU23","LU9","PD114","PD115","PD82","TEST"]

# These are additional IDs after the sample IDs, given by Novogene, e.g 
# H75TCDMXX is the ID of the sequencer, L1 is the first lane
EGYPT_SAMPLES_TO_PREPLANES = {
    "LU18":  ["NDHG02363_H75HVDMXX_L1", "NDHG02363_H75TCDMXX_L1", \
              "NDHG02363_H75HVDMXX_L2", "NDHG02363_H75TCDMXX_L2", \
              "NDHG02363_H75FVDMXX_L1", "NDHG02363_H75FVDMXX_L2"],
    "LU19":  ["NDHG02358_H7777DMXX_L1", "NDHG02358_H7777DMXX_L2"],
    "LU2":   ["NDHG02365_H75FVDMXX_L1", "NDHG02365_H75FVDMXX_L1"],
    "LU22":  ["NDHG02364_H75LLDMXX_L1", "NDHG02364_H75LLDMXX_L2"],
    "LU23":  ["NDHG02366_H75FVDMXX_L1", "NDHG02366_H75FVDMXX_L2"],
    "LU9":   ["NDHG02362_H772LDMXX_L1", "NDHG02362_H772LDMXX_L2"],
    "PD114": ["NDHG02360_H772LDMXX_L1", "NDHG02360_H772LDMXX_L2"],
    "PD115": ["NDHG02361_H772LDMXX_L1", "NDHG02361_H772LDMXX_L2"],
    "PD82":  ["NDHG02359_H772LDMXX_L1", "NDHG02359_H772LDMXX_L2"],
    "TEST":  ["PROTOCOL_SEQUENCER_L1", "PROTOCOL_SEQUENCER_L2"] # the last is for testing purposes
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

### 1. map reads to genome
# Mapping to reference/assembly using bwa
rule vc_bwa_mem:
    input: index = "bwa_index/Homo_sapiens.{assembly}.dna.primary_assembly.sa",
           fastq_r1 = "data/raw_data/{sample}/{sample}_{infolane}_1.fq.gz",
           fastq_r2 = "data/raw_data/{sample}/{sample}_{infolane}_2.fq.gz"
    output: "variants_{assembly}/{sample}_{infolane}.sam"
    wildcard_constraints: sample="[A-Z,0-9]+", infolane="[A-Z,0-9,_]+"
    conda: "envs/variant_calling.yaml"
    shell: "bwa mem -t 8 " + \
           "bwa_index/Homo_sapiens.{wildcards.assembly}.dna.primary_assembly " + \
           "{input.fastq_r1} {input.fastq_r2} > {output}"

### 2. CleanSam
# Cleans the provided SAM/BAM, soft-clipping beyond-end-of-reference alignments 
# and setting MAPQ to 0 for unmapped reads
# java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar share/picard-2.18.9-0/picard.jar
rule vc_clean_sam:
    input: "variants_{assembly}/{sample}_{infolane}.sam"
    output: "variants_{assembly}/{sample}_{infolane}.cleaned.sam"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \ 
           "CleanSam " + \
           "I={input} " + \
           "O={output}"

### 3. Sort Sam -> output: bam + idx
# Sorting by coordinates, making an index and outputting as bam
rule vc_sort_and_index_sam:
    input: "variants_{assembly}/{sample}_{infolane}.cleaned.sam"
    output: "variants_{assembly}/{sample}_{infolane}.cleaned.bam"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "SortSam " + \
           "I={input} " + 
           "O={output} " + \
           "SORT_ORDER=coordinate " + \
           "CREATE_INDEX=true"

### 4. Fix Mate Pair Information
# verify mate-pair information between mates and fix if needed
rule vc_fix_mates:
    input: "variants_{assembly}/{sample}_{infolane}.cleaned.bam"
    output: "variants_{assembly}/{sample}_{infolane}.fixed.bam"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "FixMateInformation " + \
           "I={input} " + \
           "O={output} " + \
           "SORT_ORDER=coordinate " + \
           "CREATE_INDEX=true"

### 5. Mark Duplicates
rule vc_mark_duplicates:
    input: "variants_{assembly}/{sample}_{infolane}.fixed.bam"
    output: "variants_{assembly}/{sample}_{infolane}.rmdup.bam",
            "variants_{assembly}/{sample}_{infolane}.rmdup.txt"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
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
    output: "variants_{assembly}/{sample}.merged.bam"
    params:
        picard_in=lambda wildcards, input: "I="+" I=".join(input)
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "MergeSamFiles " + \
           "{params.picard_in} " + \
           "O={output} " + \
           "SORT_ORDER=coordinate " + \
           "CREATE_INDEX=true " + \
           "USE_THREADING=8"

### 7. Collect Alignment Summary Metrics
rule vc_alignment_metrics:
    input: "variants_{assembly}/{sample}.merged.bam"
    output: "variants_{assembly}/{sample}.stats.txt"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "CollectAlignmentSummaryMetrics " + \
           "I={input} " + \
           "O={output}"

### 8. Replace Read Groups
rule vc_replace_read_groups:
    input: "variants_{assembly}/{sample}.merged.bam"
    output: "variants_{assembly}/{sample}.merged.rg.bam"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
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
    shell: "java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "CreateSequenceDictionary " + \
           "R={input} " + \
           "O={output}"

rule vc_realign:
    input: "variants_{assembly}/{sample}.merged.rg.bam",
           "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.fa",
           "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.dict"
    output: "variants_{assembly}/{sample}.merged.rg.ordered.bam"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/share/picard-2.18.9-0/picard.jar " + \
           "ReorderSam " + \
           "I={input[0]} " + \
           "O={output} " + \
           "R={input[1]} " + \
           "CREATE_INDEX=true"

### 10. RealignerTargetCreator

# Therefore, the fasta file needs to be indexed
rule vc_inex_fasta:
    input: "seq_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output: "seq_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
    shell: "samtools faidx {input}"

# java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar
rule vc_realigner_target_creator:
    input: "variants_{assembly}/{sample}.merged.rg.ordered.bam",
           "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.fa",
           "seq_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
    output: "variants_{assembly}/{sample}.merged.rg.ordered.bam.intervals"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/opt/gatk-3.8/GenomeAnalysisTK.jar " + \
           "-T RealignerTargetCreator " + \
           "-R {input[1]} " + \
           "-I {input[0]} " + \
           "-o {output} " + \
           "-nt 8"

### 11. IndelRealigner
rule vc_indel_realigner:
    input: "variants_{assembly}/{sample}.merged.rg.ordered.bam",
           "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.fa",
           "variants_{assembly}/{sample}.merged.rg.ordered.bam.intervals"
    output: "variants_{assembly}/{sample}.indels.bam"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/opt/gatk-3.8/GenomeAnalysisTK.jar " + \
           "-T IndelRealigner " + \
           "-R \"{input[1]}\" " + \
           "-I {input[0]} " + \
           "-targetIntervals {input[2]} " + \
           "-o {output}"

### 12. Base Quality Recalibration
rule vc_base_recalibrator:
    input: "variants_{assembly}/{sample}.indels.bam",
           "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.fa",
           "dbsnp_{assembly}/All_20180418.vcf.gz",
           "dbsnp_{assembly}/All_20180418.vcf.gz.tbi"
    output: "variants_{assembly}/{sample}.indels.recal.csv"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/opt/gatk-3.8/GenomeAnalysisTK.jar " + \
           "-T BaseRecalibrator " + \
           "-R {input[1]} " + \
           "-I {input[0]} " + \
           "-cov ReadGroupCovariate " + \
           "-cov QualityScoreCovariate " + \
           "-cov CycleCovariate " + \
           "-cov ContextCovariate " + \
           "-o {output} " + \
           "-knownSites {input[2]} " + \
           "-nct 8"

### 13. Print Reads
rule vc_print_reads:
    input: "variants_{assembly}/{sample}.indels.bam",
           "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.fa",
           "variants_{assembly}/{sample}.indels.recal.csv"
    output: "variants_{assembly}/{sample}.final.bam"
    conda: "envs/variant_calling.yaml"
    shell: "java -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/opt/gatk-3.8/GenomeAnalysisTK.jar " + \
           "-T PrintReads " + \
           "-R {input[1]} " + \
           "-I {input[0]} " + \
           "-o {output[0]} " + \
           "-BQSR {input[2]} " + \
           "-nct 8"

### 14. variant calling with GATK-HC
# use GATK Haplotypecaller with runtime-optimized settings
# -variant_index_type LINEAR -variant_index_parameter 128000 IW added because
# the GATK program told me so and otherwise would exit with error.
# java -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=4 -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar gatk 
rule vc_snp_calling_with_gatk_hc:
    input: "variants_{assembly}/{sample}.final.bam",
           "seq_{assembly}/Homo_sapiens.{assembly}.dna.primary_assembly.fa",
    output: "variants_{assembly}/{sample}.vcf"
    conda: "envs/variant_calling.yaml"
    shell: "java -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=4 -Xmx35g -Djava.io.tmpdir=/data/lied_egypt_genome/tmp -jar .snakemake/conda/d590255f/opt/gatk-3.8/GenomeAnalysisTK.jar " + \
           "-T HaplotypeCaller " + \
           "-R {input[1]} " + \
           "-I {input[0]} " + \
           "--genotyping_mode DISCOVERY " + \
           "-o {output} " + \
           "-ERC GVCF " + \
           "-variant_index_type LINEAR " + \
           "-variant_index_parameter 128000 " + \
           "-nct 8"

# Doing the variant calling for all 9 samples
rule vc_snp_calling_with_gatk_hc_all:
    input: expand("variants_GRCh38/{sample}.vcf", sample=EGYPT_SAMPLES)

################################################################################
################### Assembly assessment and correction #########################
################################################################################


