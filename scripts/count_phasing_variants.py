# Because we use the cohort-called SNVs, it is not so straight-forward to count 
# the number of variants Egyptref-individual level variant with respect to 
# certain properties; thus, here we output a binary string for every read 
# denoting whether he fulfills the property (1) or not (0). Main application
# is that we want to use this to correctly reproduce the percentage of variants
# phased which is reported by the 10x genomics longranger phasing summary

# Getting input and output filenames
fname_in = snakemake.input[0]
fname_out = snakemake.output[0]

with open(fname_in,"r") as f_in, open(fname_out,"w") as f_out:
    header_out = ["chrom","ref","alt","al1","al2","phasing"] 
    header_out += ["homozygous","homozygous_ref","homozygous_alt"]
    header_out += ["homozygous_alt1","homozygous_alt_2plus","homozygous_missing","homozygous_indel"]
    header_out += ["homozygous_snv","heterozygous","heterozygous_missing"]
    header_out += ["heterozygous_alt1","heterozygous_ref_alt_multi"]
    header_out += ["heterozygous_ref_alt_multi_indel"]
    header_out += ["heterozygous_alt_alt_multi_snv"]
    header_out += ["multiallelic","multiallelic_incl_indel"]
    header_out += ["multiallelic_no_indel","monoallelic"]
    header_out += ["monoalellic_ref_indel","monoallelic_alt_indel"]
    header_out += ["monoallelic_ref_alt_indel","autosomal"]
    f_out.write(" ".join(header_out)+"\n")
    for line in f_in:
        # Initialize
        autosomal = "0"
        homozygous = "0"
        homozygous_ref = "0"
        homozygous_alt = "0"
        homozygous_alt_2plus = "0"
        homozygous_alt1 = "0"
        homozygous_missing = "0"
        homozygous_indel = "0"
        homozygous_snv = "0"
        heterozygous = "0"
        heterozygous_missing = "0"
        heterozygous_alt1 = "0"
        heterozygous_ref_alt_multi = "0"
        heterozygous_ref_alt_multi_indel = "0"
        heterozygous_alt_alt_multi_snv = "0"
        multiallelic = "0"
        multiallelic_incl_indel = "0"
        multiallelic_no_indel = "0"
        monoallelic = "0"
        monoalellic_ref_indel = "0"
        monoallelic_alt_indel = "0"
        monoallelic_ref_alt_indel = "0"
        # Skip header line
        if line[0] == "#":
            continue
        s = line.split("\t")
        genotype = s[9].split(":")[0]
        ref = s[3]
        assert(len(ref.split(",")) == 1)
        alt = s[4]
        alt_alleles = alt.split(",")
        chrom = s[0]
        # Autosomal or not
        if chrom in [str(x) for x in range(1,23)]:
            autosomal = "1"
        # Phased or not 
        if "/" in genotype:
            phasing = "unphased"
            al1,al2 = genotype.split("/")
        else:
            assert ("|" in genotype)
            phasing = "phased"
            al1,al2 = genotype.split("|")
        # Homozygous
        if al1 == al2:
            homozygous = "1"
            if al1 == "." and al2 == ".":
                homozygous_missing = "1"
            else:
                if al1 == "0" and al2 == "0":
                    homozygous_ref = "1"
                else:
                    homozygous_alt = "1"
                if al1 == "1" and al2 == "1":
                    homozygous_alt1 = "1"
                elif not (al1 == "0" and al2 == "0"): 
                    homozygous_alt_2plus = "1"
                if len(ref)>1 or len(alt_alleles[int(al2)-1])>1:
                    homozygous_indel = "1"
                else:
                    homozygous_snv = "1"
        # Heterozygous
        else:
            heterozygous = "1"
            if al1 == "." or al2 == ".":
                heterozygous_missing = "1"
            elif al1 == "0" and al2 == "1" or al1 == "1" and al2 == "0":
                heterozygous_alt1 = "1"
            elif al1 == "0" or al2 == "0":
                heterozygous_ref_alt_multi = "1"
                if len(ref)>1 or len(alt_alleles[max(int(al1),int(al2))-1])>1:
                    heterozygous_ref_alt_multi_indel = "1"
                else:
                    heterozygous_ref_alt_multi_snv = "1"
            else:
                if len(alt_alleles[int(al1)-1])>1 or len(alt_alleles[int(al2)-1])>1:
                    heterogous_alt_alt_multi_indel = "1"
                else:
                    heterozygous_alt_alt_multi_snv = "1"
        # Cohort-based multi-allelic
        if len(alt_alleles) > 1:
            multiallelic = "1"
            for allele in alt_alleles:
                if len(allele)>1:
                    multiallelic_incl_indel = "1"
            if multiallelic_incl_indel == "0":
                multiallelic_no_indel = "1"
        # Cohort-based indel/SNV
        else:
            monoallelic = "1"
            if len(ref)>1:
                monoalellic_ref_indel = "1"
            if len(alt) > 1:
                monoallelic_alt_indel = "1"
            if len(ref)>1 and len(alt)>1:
                monoallelic_ref_alt_indel = "1"
        out = [ chrom, ref, alt, al1, al2, phasing] 
        out += [ homozygous, homozygous_ref, homozygous_alt]
        out += [ homozygous_alt1,homozygous_alt_2plus,homozygous_missing,homozygous_indel]
        out += [ homozygous_snv,heterozygous,heterozygous_missing]
        out += [ heterozygous_alt1, heterozygous_ref_alt_multi]
        out += [ heterozygous_ref_alt_multi_indel]
        out += [ heterozygous_alt_alt_multi_snv]
        out += [ multiallelic, multiallelic_incl_indel]
        out += [ multiallelic_no_indel, monoallelic]
        out += [ monoalellic_ref_indel, monoallelic_alt_indel]
        out += [ monoallelic_ref_alt_indel,autosomal]
        f_out.write(" ".join(out)+"\n")