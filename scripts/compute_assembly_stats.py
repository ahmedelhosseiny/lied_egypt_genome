# Computing assembly statistics, i.e. N50, N60 etc.
# We compute N50 in accordance with the NCBI N50 values, that is, the size of
# the last scaffold entirely included within 50% of the bases; accordingly,
# the L50 value is the number of contigs entirely included in 50% of bases

# Getting input and output filenames
fname_in = snakemake.input[0]
fname_out = snakemake.output[0]

# Performing computation
lengths = []
with open(fname_in,"r") as f_in, open(fname_out,"w") as f_out:
    for line in f_in:
        if line[:8] == "SCAFFOLD":
            continue
        lengths.append(int(line.split("\t")[1]))
    lengths = sorted(lengths,reverse = True)
    n = sum(lengths)
    stat_type = [50,60,70,80,90]
    for stat in stat_type:
        num_bases_to_sum = stat *n/100
        sum_scaffold_lens = 0
        num_contigs = 0
        for scaffold_len in lengths:
            if sum_scaffold_lens < num_bases_to_sum:
                sum_scaffold_lens += scaffold_len
                last_scaffold_len = scaffold_len
            else:
                break
            num_contigs += 1
        f_out.write("N"+str(stat)+": "+str(last_scaffold_len)+"\n")
        f_out.write("L"+str(stat)+": "+str(num_contigs)+"\n")