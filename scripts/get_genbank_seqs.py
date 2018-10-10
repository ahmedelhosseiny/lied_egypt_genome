# Getting a Genbank fasta sequence from Entrey using Biopython

from Bio import Entrez, SeqIO

# Getting input and output filenames
genbank_id = wildcards.genbank_id
fname_out = snakemake.output[0]


record = open('als.fasta', 'w')
for seq_id in ids:
    handle = Entrez.efetch(db="nucleotide", id="seq_id", rettype="fasta", retmode="text")
    record = handle.read()
    record.write(record.rstrip('\n'))

