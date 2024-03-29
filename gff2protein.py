#!/usr/bin/env python3
"""AM 19DEC19
Convert GFF3 gene predictions into protein sequences.
This tool is splice-site aware, and will remove introns prior to translation.

Usage:
    gff3_to_proteins.py <gff3> <fasta> -id <featureID>

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
import argparse

def cds_to_exon(gff_file):
    """Adjusts CDS features to exon for parsing with BCBio GFF module. OBSOLETE! USED DURING TESTING ONLY.
    """
    with open(gff_file, "r") as f:
        newGFF=f.read().replace("CDS", "exon")

    out_handle="{gff}.tmp".format(gff=gff_file)

    with open(out_handle, "w") as f:
        f.write(newGFF)

    return out_handle

def extract_recs(gff_file, fa_rec, featureID, prot_file, cdna_file):
    prot_f = open(prot_file, "w")
    cdna_f = open(cdna_file, "w")

    db = gff_predictions(gff_file)

    # print(db.count_features_of_type(featuretype=None))

    for rec in db.features_of_type(featureID, order_by="start"):

        cdna_seq = ""
        c=0

        for child in db.children(rec.id, featuretype="exon"):
            c+=1
            cdna_seq += child.sequence(fa_rec, use_strand=True)

        ## check how many children were found, if none check CDSs
        if c==0:
            for child in db.children(rec.id, featuretype="CDS"):
                c+=1
                cdna_seq += child.sequence(fa_rec, use_strand=True)

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

def parseArgs(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gff', type=str, help='GFF file containing gene features.')
    parser.add_argument('fasta', type=str, help='Nucleotide sequence file in FASTA format')
    parser.add_argument('-id', '--featureID', type=str, default='gene', help='Name of the feature to be extracted')

    return parser.parse_args()

def main(args):
    fa_rec = pyfaidx.Fasta(args.fasta)

    prot_file = "{base}-prot.fa".format(base=os.path.splitext(args.gff)[0])
    cdna_file = "{base}-cdna.fa".format(base=os.path.splitext(args.gff)[0])

    # gff_file = cds_to_exon(gff_file)
    extract_recs(args.gff, fa_rec, args.featureID, prot_file, cdna_file)

if __name__ == "__main__":
    args = parseArgs(sys.argv)
    main(args)
