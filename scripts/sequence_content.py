# Making a file with sequence information
from Bio import SeqIO

# Getting input and output filenames
fname_in = snakemake.input[0]
fname_out = snakemake.output[0]

# Counting sequence lengths of scaffolds
with open(fname_in, "rU") as handle:
    sum_bases = 0
    for record in SeqIO.parse(handle, "fasta"):
        sum_bases += len(record)
# Computing numbers and writing them to file
with open(fname_in, "rU") as handle, open(fname_out, "w") as f_out:
    header = ["SCAFFOLD","LEN","A","C","G","T","N"]
    header += ["PERCENT_A","PERCENT_C","PERCENT_G","PERCENT_T","PERCENT_BASES"]
    f_out.write("\t".join(header)+"\n")
    for record in SeqIO.parse(handle, "fasta"):
        len_scaf = len(record)
        n_a = record.seq.count("A")
        n_c = record.seq.count("C")
        n_g = record.seq.count("G")
        n_t = record.seq.count("T")
        n_n = record.seq.count("N")
        n_acgtn = sum([n_a,n_c,n_g,n_t,n_n])
        # Make sure that the sequence contains only ACGTN bases
        assert(n_acgtn == len(record))
        percent = [round(100*x/n_acgtn,2) for x in [n_a,n_c,n_g,n_t,n_n]]
        percent_bases = round(100*len(record)/sum_bases,2)
        numbers = [len_scaf,n_a,n_c,n_g,n_t,n_n]+percent+[percent_bases]
        f_out.write(record.id+"\t" \
                             +"\t".join([str(x) for x in numbers]) \
                             +"\n")