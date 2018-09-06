#kate:syntax python;

#######################################
### Analyzing an egyptian genome
#######################################


rule scaffold_names:
    input: "data/pilon.fasta"
    output: "results/num_scaffolds.txt"
    shell: "cat {input} | grep '>' > {output}"
    



