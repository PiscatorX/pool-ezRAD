#!/usr/bin/env python 
from Bio import SeqIO
import argparse
import pprint


def contig2clstr(contig_sequences,  clusters):

    cluster_ids_ref = dict([ line.split() for line in clusters])
    #pprint.pprint(cluster_ids)

    #exit(1)
    
    for record in SeqIO.parse(contig_sequences, "fasta"):
        contig_id  = record.id.split('_')[1]
        cluster_id = cluster_ids_ref[contig_id]
        print("{}\t{}\t{}".format(contig_id,cluster_id,record.seq))

if __name__ == '__main__':
    parser = argparse.ArgumentParser("text to fasta")
    parser.add_argument('-s','--sequences',help ="textfile", type=argparse.FileType('r'), required=True)
    parser.add_argument('-c','--cluster_ids', help ="output filename",default = "uniq.fasta", type=argparse.FileType('r'), required=True)
    args = parser.parse_args()
    contig2clstr(args.sequences, args.cluster_ids)

    


