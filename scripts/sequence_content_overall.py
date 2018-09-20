from Bio import SeqIO

# Getting input and output filenames
fname_in = snakemake.input[0]
fname_out = snakemake.output[0]

# Counting sequence lengths of scaffolds
with open(fname_in, "rU") as handle, open(fname_out, "w") as f_out:
    num_records = 0
    n = 0
    n_a = 0
    n_c = 0
    n_g = 0
    n_t = 0
    n_n = 0
    for record in SeqIO.parse(handle, "fasta"):
        num_records += 1
        n += len(record)
        n_a += record.seq.count("A")
        n_c += record.seq.count("C")
        n_g += record.seq.count("G")
        n_t += record.seq.count("T")
        n_n += record.seq.count("N")
    f_out.write("NUM_RECORDS: "+str(num_records)+"\n")
    f_out.write("NUM_BASES: "+str(n)+"\n")
    f_out.write("NUM_A: "+str(n_a)+"\n")
    f_out.write("NUM_C: "+str(n_c)+"\n")
    f_out.write("NUM_G: "+str(n_g)+"\n")
    f_out.write("NUM_T: "+str(n_t)+"\n")
    f_out.write("NUM_N: "+str(n_n)+"\n")
    f_out.write("PERCENT_A: "+str(round(100*n_a/n,2))+"\n")
    f_out.write("PERCENT_C: "+str(round(100*n_c/n,2))+"\n")
    f_out.write("PERCENT_G: "+str(round(100*n_g/n,2))+"\n")
    f_out.write("PERCENT_T: "+str(round(100*n_t/n,2))+"\n")
    f_out.write("PERCENT_N: "+str(round(100*n_n/n,2))+"\n")