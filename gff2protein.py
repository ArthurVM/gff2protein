#!/usr/bin/env python3
"""AM 19DEC19
Convert GFF3 gene predictions into protein sequences.
This tool is splice-site aware, and will remove introns prior to translation.

Usage:
    gff3_to_proteins.py <gff3> <fasta> <gene_flag>

where <gene_flag> is the name of gene features within the gff (e.g. gene/protein_coding_gene)

Output:
gff database        ./<gff_basename>-.db
protein sequence    ./<gff_basename>-prot.fa
cDNA sequence       ./<gff_basename>-cdna.fa

Dependancies:
Biopython
pyfaidx
gffutils
"""

import sys
import os

from Bio.Seq import Seq
import pyfaidx
import gffutils

def main(gff_file, fa_file, gene_flag):
    fa_rec = pyfaidx.Fasta(fa_file)

    prot_file = "{base}-prot.fa".format(base=os.path.splitext(gff_file)[0])
    cdna_file = "{base}-cdna.fa".format(base=os.path.splitext(gff_file)[0])

    # gff_file = cds_to_exon(gff_file)
    extract_recs(gff_file, fa_rec, prot_file, cdna_file, gene_flag)

def extract_recs(gff_file, fa_rec, prot_file, cdna_file, gene_flag, intron_flag=False):
    prot_f = open(prot_file, "w")
    cdna_f = open(cdna_file, "w")

    db = gff_predictions(gff_file)

    print(f"Extracting {gene_flag} features from {gff_file}")

    # print(db.count_features_of_type(featuretype=None))

    for rec in db.features_of_type(gene_flag, order_by="start"):

        cdna_seq = ""
        c=0

        ## handle introns
        if intron_flag:
            for child in db.children(rec.id, featuretype="exon"):
                c+=1
                cdna_seq += child.sequence(fa_rec, use_strand=True)

            ## check how many children were found, if none check CDSs
            if c==0:
                for child in db.children(rec.id, featuretype="CDS"):
                    c+=1
                    cdna_seq += child.sequence(fa_rec, use_strand=True)
        
        ## handle intronless genomes
        else:
            cdna_seq += rec.sequence(fa_rec, use_strand=True)

        prot_seq = Seq.translate(Seq(cdna_seq), stop_symbol="*", to_stop=False, gap="X")

#         if rec.strand == "+":
#             prot_seq = Seq.translate(Seq(cdna_seq), stop_symbol="*", to_stop=False, gap="X")
#         elif rec.strand == "-":
#             prot_seq = Seq.translate(Seq(cdna_seq).reverse_complement(), stop_symbol="*", to_stop=False, gap="X")

        header = "{id} | {scaff}:{start}-{end} number_exons={nexons} strand={strand}".format(\
        id=rec.id, scaff=rec[0], start=rec.start, end=rec.stop, nexons=c, strand=rec[6])
        cdna_rec_line = ">{header}\n{seq}\n".format(header=header, seq=cdna_seq)
        prot_rec_line = ">{header}\n{seq}\n".format(header=header, seq=prot_seq)

        cdna_f.write(cdna_rec_line)
        prot_f.write(prot_rec_line)

    prot_f.close()
    cdna_f.close()

    print()
    print(f"cDNA FASTA : {cdna_file}")
    print(f"Protein FASTA : {prot_file}")
    print()

def gff_predictions(gff_file):
    """Parse gff3 output and yield SeqFeatures and SeqRecords
    """

    db_handle = "{base}.db".format(base=os.path.splitext(gff_file)[0])

    if os.path.isfile(db_handle) is True:
        print("Found gff db")
        gff_db = gffutils.FeatureDB(db_handle)

    else:
        print("gff db not found. Now creating a new one...")
        gff_db = gffutils.create_db(gff_file, "{base}.db".format(base=os.path.splitext(gff_file)[0]))

    return gff_db

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit()
    main(*sys.argv[1:])
