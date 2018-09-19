#kate:syntax python;

#######################################
### Analyzing an Egyptian genome
#######################################

from Bio import SeqIO


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
    output: "busco_{assembly}/run_busco_{assembly}_{chr_or_type}/short_summary_busco_{assembly}_{chr_or_type}.txt"
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
# chromosome
# Parameters from the Korean reference genome AK1 (Seo et al. 2016)
# --gapped Perform gapped extension of HSPs after first reducing them to anchor 
#          points
# --gap=600,150 Gap open and gap extension penalty
# --hspthresh=4500 Set the score threshold for the x-drop extension method; HSPs
#                  scoring lower are discarded.
# --seed 12of19 Seeds require a 19bp word with matches in 12 specific positions
# --notransition Don't allow any match positions in seeds to be satisified by 
#                transitions
# --ydrop=15000-chain Set the threshold for terminating gapped extension; this
#                     restricts the endpoints of each local alignment by 
#                     limiting the local region around each anchor in which 
#                     extension is performed
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
rule align_with_lastz_quick:
    input: "repeatmasked_GRCh38/Homo_sapiens.GRCh38.dna.{chr}.fa.masked",
           "repeatmasked_EGYPTREF/Homo_sapiens.EGYPTREF.dna.{scaffold}.fa.masked"
    output: "alignments_lastz/dotplots/{chr}_vs_{scaffold}.maf",
            "alignments_lastz/dotplots/{chr}_vs_{scaffold}.rdotplot"
    conda: "envs/lastz.yaml"
    shell: "lastz {input[0]} {input[1]} " + \
                                  "--notransition " + \
                                  "--step=20 " + \
                                  "--nogapped " + \
                                  "--format=maf " + \
                                  "--hspthresh=4500 " + \
                                  "--seed=12of19 " + \
                                  "--rdotplot={output[1]} " + \
                                  ">{output[0]}"

# Plot the dotplot output of lastz
rule individual_lastz_dotplot:
    input: "alignments_lastz/dotplots/{chr}_vs_{scaffold}.rdotplot"
    output: "alignments_lastz/dotplots/{chr}_vs_{scaffold}.pdf"
    script: "scripts/dotplot.R"

# Plotting for one scaffold the dotplot versus all chromosomes
rule dotplots_scaffold_vs_chromosomes:
    input: expand("alignments_lastz/dotplots/{chr}_vs_{{scaffold}}.rdotplot", \
                  chr=CHR_GRCh38)
    output: "alignments_lastz/dotplots/{scaffold}.pdf"
    script: "scripts/scaffold_vs_grch38.R"            

# Plotting the dotplots for all scaffolds
rule dotplots_scaffold_vs_chromosomes_all:
    input: expand("alignments_lastz/dotplots/{scaffold}.pdf", \
                  scaffold=EGYPT_SCAFFOLDS)

# All versus all comparisons of reference and Egyptian genome
rule align_all_vs_all:
    input: expand("alignments_lastz/{chr}_vs_{scaffold}.maf", \
                  chr=CHR_GRCh38, scaffold=EGYPT_SCAFFOLDS)

