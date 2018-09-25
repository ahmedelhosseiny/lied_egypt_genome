# Making a summary table of the repeatmasker results for an assembly

# Getting input and output filenames
fnames_in = snakemake.input
fname_out = snakemake.output[0]

# scaffold_info is a hash table to which the scaffold info is written for 
# sorting according to scaffold size (largest to smallest) before writing to
# out file
scaffold_info = {}
header = ["filename", "tot_len", "tot_len_excl_xn", "gc_level", \
          "masked_bp", "masked_percent", "tot_interspersed_repeats_len", \
          "tot_interspersed_repeats_percent"]
repeat_types = ["sine", "alu", "mir", "line", "line1", "line2", "l3_cr1", "ltr", \
                "ervl", "ervl_malrs", "erv_class1", "erv_class2", "dna_elements", \
                "hat_chalie", "tcmar_tigger", "unclassified", \
                "small_rna", "sattelites", "simple_repeats", "low_complexity"]
# Append header for all the repeat types if which repeat number, bp length and
# percent is given
for repeat in repeat_types:
    header.append(repeat+"_num")
    header.append(repeat+"_len")
    header.append(repeat+"_percent")
scaffold_info = []
for fname in fnames_in:
    with open(fname,"r") as f_in:
        num_len_perc = []
        for line in f_in:
            # Skip lines wo information
            if line[0] in ["=","\n","-"] \
               or "percentage" in line \
               or " of sequence" in line \
               or "* most" in line \
               or "have been counted" in line:
                continue
            splitted = line.split(":")
            if len(splitted) == 2:
                item,spec = splitted
            else:
                 item = ""
                 spec = splitted[0]
            spec = spec.strip()
            if item == "file name":
                assert(spec == fname.split("/")[1][:-4])
                fname_wo_path = spec
            elif item == "sequences":
                assert(spec == "1")
            elif item == "total length":
                total = [x for x in spec.split(" ") if not x in ["bp","%",""]]
                tot_len = total[0]
                tot_len_no_n = total[1].strip("(")
            elif item == "GC level":
                gc = spec[:-2].strip()
            elif item == "bases masked":
                total = [x for x in spec.split(" ") if not x in ["bp","%","","("]]    
                masked = total[0]
                masked_percent = total[1].strip("(")
            # The next category has only len and perecent
            elif item == "Total interspersed repeats":
                ti_rep = [x for x in spec.split(" ") if not x in ["bp","%",""]]
            # This parses all categories for which three numbers are given,
            # number, length and percent
            else:
                info = [x for x in spec.split(" ") if not x in ["bp","%",""]]
                if not ":" in line:
                    info = info[1:]
                num_len_perc += info
            # Low complexity is the last category to be summarised
            # If this was read, then we can save the entire info for this 
            # scaffold
            if "Low complexity" in line:
                scaffold_info.append([fname_wo_path,tot_len,tot_len_no_n,gc, \
                                        masked,masked_percent, \
                                        ti_rep[0],ti_rep[1]]+num_len_perc)
                continue

# Sort by overall scaffold length and print header and then info to file
sorted_scaffinfo = sorted(scaffold_info,key = lambda x: int(x[1]), reverse=True)
with open(fname_out,"w") as f_out:
    f_out.write("\t".join(header)+"\n")
    for elem in sorted_scaffinfo:
        f_out.write("\t".join(elem)+"\n") 
    
                

                