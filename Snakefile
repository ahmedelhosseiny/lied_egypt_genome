#kate:syntax python;

#######################################
### Analyzing an Egyptian genome
#######################################

from Bio import SeqIO

rule scaffold_names:
    input: "data/pilon.fasta"
    output: "results/scaffold_names.txt"
    shell: "cat {input} | grep '>' > {output}"
    
# Quantifying the sequence content individually for all scaffolds
rule sequence_content:
    input: "data/pilon.fasta"
    output: "results/num_bases.txt",
    run:
        # Counting sequence lengths of scaffolds
        with open(input[0], "rU") as handle:
            sum_bases = 0
            for record in SeqIO.parse(handle, "fasta"):
                sum_bases += len(record)
        # Computing numbers and writing them to file
        with open(input[0], "rU") as handle, open(output[0], "w") as f_out:
            header = ["SCAFFOLD","LEN","A","C","G","T","N"]
            header += ["PERCENT_A","PERCENT_C","PERCENT_G","PERCENT_T","PERCENT_BASES"]
            f_out.write("\t".join(header))
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

# Quantifying the sequence content over all scaffolds
rule sequence_content_overall:
    input: "data/pilon.fasta"
    output: "results/num_all.txt",
    run:
        # Counting sequence lengths of scaffolds
        with open(input[0], "rU") as handle, open(output[0], "w") as f_out:
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
            
# Compute N50 and other related values as statistic for the assembly
rule compute_assembly_stats:
    input: "results/num_bases.txt"
    output: "results/assembly_stats.txt"
    run:
        lengths = []
        with open(input[0],"r") as f_in, open(output[0],"w") as f_out:
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
                stop = False
                for scaffold_len in lengths:
                    if sum_scaffold_lens < num_bases_to_sum:
                        sum_scaffold_lens += scaffold_len
                    else:
                        break
                f_out.write("N"+str(stat)+": "+str(scaffold_len)+"\n")
                    
                

                



