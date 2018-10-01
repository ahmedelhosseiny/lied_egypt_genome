# Compiling the busco chromosomewise results

# Getting input and output filenames
fnames_in = snakemake.input
fname_out = snakemake.output[0]

allowed_status = ["Complete", "Missing", "Fragmented", "Duplicated"]
all_buscos = {}
buscos_complete = {}
buscos_fragmented = {}

# Read in from Busco result files into hash tables common for all chromosomes
for filename in fnames_in:
    with open(filename, "r") as f_in:
        for line in f_in:
            # Skip comment lines
            if line[0] == '#':
                continue
            s = line.strip("\n").split("\t")
            busco = s[0]
            status = s[1]
            assert(status in allowed_status)
            if not busco in all_buscos:
                all_buscos[busco] = 1
            else: 
                all_buscos[busco] += 1
            if status == "Missing":
                continue
            if status == "Complete" or status == "Duplicated":
                if not busco in buscos_complete:
                    buscos_complete[busco] = 1
                else: 
                    buscos_complete[busco] += 1
            if status == "Fragmented":
                if not busco in buscos_fragmented:
                    buscos_fragmented[busco] = 1
                else: 
                    buscos_fragmented[busco] += 1

# Classifying into complete, single, double, fragmented and missing and writing 
# to result table
double = 0
single = 0
fragmented = 0
missing = 0
num_buscos = 4104
with open(fname_out,"w") as f_out:
    for busco in sorted(all_buscos.keys()):
        if busco in buscos_complete:
            if buscos_complete[busco]>1:
                double += 1
            else:
                single += 1
        elif busco in buscos_fragmented:
                fragmented += 1
        else:
                missing += 1
    assert((double+single+fragmented+missing) == num_buscos)
    c = 100*(double+single)/num_buscos
    s = 100*single/num_buscos
    d = 100*double/num_buscos
    f = 100*fragmented/num_buscos
    m = 100*missing/num_buscos
    f_out.write("C:"+str(round(c,2)) + \
              "%[S:"+str(round(s,2)) + \
              "%,D:"+str(round(d,2)) + \
              "%],F:"+str(round(f,2)) + \
              "%,M:"+str(round(m,2)) + \
              "%,n:"+str(num_buscos))
 