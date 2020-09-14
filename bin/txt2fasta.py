#!/usr/bin/env  python
import argparse

def txt2fasta(text_fp, outfile):

    for contig, line in enumerate(text_fp,1):
        if line[0].isdigit():
            continue
        print(">Contig_{}\n{}".format(contig, line.strip()),file = outfile)
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser("text to fasta")
    parser.add_argument('text',help ="textfile", type=argparse.FileType('r'))
    parser.add_argument('-o','--outfile', help ="output filename",default = "uniq.fasta", type=argparse.FileType('w'))
    args = parser.parse_args()
    txt2fasta(args.text, args.outfile)
