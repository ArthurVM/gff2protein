# gff_to_protein
A simple python script for taking a gff-fasta pair and extracting protein and cDNA sequences. A. V. Morris 19DEC19.

This tool is exon-aware, stripping out introns to generate cDNA and protein sequence for each gene feature.

Usage:
------
    gff3_to_proteins.py <gff3> <fasta>
    
**Output:**

    gff database        ./<gff_basename>-.db
    protein sequence    ./<gff_basename>-prot.fa
    cDNA sequence       ./<gff_basename>-cdna.fa

Dependencies:
-------------
    Python/Python3
    Biopython
    pyfaidx
    gffutils

To install dependencies:

    pip install -r ./dependencies.txt

N.B. it will query each gene for exons identified by 'CDS' in the feature type column, rather than 'exon'.
