# This script filters Nucmer delta alignment files to contain only alignments
# covering a region specified in a bed file.

# Getting input and output filenames
input_gene = snakemake.input["gene"]
input_align = snakemake.input["align"]
output = snakemake.output[0]

align_write = False
check_coords = False
# Get gene start and end positions from bed file
with open(input_gene,"r") as f_in:
    for line in f_in:
        if line[0] == "#":
            continue
        s = line.split("\t")
        chrom = s[0]
        start = int(s[1])
        end = int(s[2])
        assert(start<end)
# Got through alignment file and write to output those alignments
# covering (part of) the gene
# In this table we remember the chromsome/scaffold combis written so far
written_chrom_header_lines = []
with open(input_align,"r") as f_in, open(output,"w") as f_out:
    header_count = 0
    for line in f_in:
        # Write the two header lines, the first contains the links to 
        # the two input assembly fasta files, the second just contains 
        # 'NUCMER'
        if header_count < 2:
            f_out.write(line)
            header_count += 1
            continue
        s = line.split(" ")
        if line[0] == ">":
            align_write = False
            al_chrom = s[0][1:]
            if al_chrom == chrom:
                chrom_header_line = line
                check_coords = True
            else:
                check_coords = False
                align_write = False
            continue
        if check_coords and not len(s) == 1:
            al_start = int(s[0])
            al_end = int(s[1])
            # There is an alignment covering the desired region
            if start<=al_start<=end or start<=al_end<=end \
               or (al_start <= start and al_end >= end):
                align_write = True
                if not chrom_header_line in written_chrom_header_lines:
                    f_out.write(chrom_header_line)
                    written_chrom_header_lines.append(chrom_header_line)
            else:
                align_write = False
        if align_write:
            f_out.write(line)