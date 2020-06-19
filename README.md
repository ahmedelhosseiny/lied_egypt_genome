# An Egyptian personal and population-based Egyptian genome reference (lied\_egypt\_genome)  

**Author:** Inken Wohlers  

**Summary:** These are the workflows used to generate results presented in the Egyptref paper. Workflows are topic-specific.  

**Egyptref Project:**  
A small number of de novo assembled human genomes have been reported to date, and few have been complemented with population-based genetic variation, which is particularly important for North Africa, a region underrepresented in current genome-wide references. Here, we combine long- and short-read whole-genome sequencing data with recent assembly approaches into a de novo assembly of an Egyptian genome. The assembly demonstrates well-balanced quality metrics and is complemented with variant phasing via linked reads into haploblocks, which we associate with gene expression changes in blood. To construct an Egyptian genome reference, we identify genome-wide genetic variation within a cohort of 110 Egyptian individuals. We show that differences in allele frequencies and linkage disequilibrium between Egyptians and Europeans may compromise the transferability of European ancestry-based genetic disease risk and polygenic scores, substantiating the need for multi-ethnic genome references. Thus, the Egyptian genome reference will be a valuable resource for precision medicine.  

**Instructions:**  
Raw base files need to be in place and Bioconda installed (https://bioconda.github.io).  
Please note: Some analyses need a considerable amount of memory (up to 320 GB RAM) and/or CPU (up to 48 cores).  
Run (parts of) the workflow to generate a target file via  
```
snakemake -s Snakefilename --use-conda TARGETFILENAME  
```  
