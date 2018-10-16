# Getting a Genbank fasta sequence from Entrez using Biopython

from Bio import Entrez, SeqIO

# It is advised to set the e-mail, so NCBI can contact one before blocking:
# "In case of excessive usage of the E-utilities, NCBI will attempt to contact
# a user at the email address provided before blocking access to the 
# E-utilities."
Entrez.email="Inken.Wohlers@uni-luebeck.de"

# Getting input and output filenames
genbank_id = snakemake.params.genbank_id
fname_out = snakemake.output[0]

# Obtaining the fasta entry and writing it to file
with open(fname_out, 'w') as f_out:
    handle = Entrez.efetch(db="nucleotide", \
                           id=genbank_id, \
                           rettype="fasta", \
                           retmode="text")
    record = handle.read()
    f_out.write(record[:-1]) # There is a newline too much, which is removed

