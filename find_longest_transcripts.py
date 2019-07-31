#!/usr/bin/env python3
import sys
from Bio import SeqIO
import re
import collections
import gzip

def get_longest_transcript(fasfile=None):
    if fasfile.endswith('.gz'):
        fh_in = gzip.open(fasfile, 'rt')
    else:
        fh_in = open(fasfile)

    gene_seqid_len = {}
    for rec in SeqIO.parse(fh_in, 'fasta'):
        m = re.search(r'\[gene\=(.+?)\]', rec.description)
        if not m:
            m = re.search(r'\[locus_tag\=(.+?)\]', rec.description)
            if not m:
                print("{} : No '[gene=XXX]' or '[locus_tag=XXX]', will not be output!".format(rec.id), file=sys.stderr)
            continue
        gene = m.group(1)
        transcript_len = len(rec)
        if gene in gene_seqid_len:
            pre_seqid, pre_len = gene_seqid_len[gene]
            if transcript_len > pre_len:
                gene_seqid_len[gene] = (rec.id, transcript_len)
        else:
            gene_seqid_len[gene] = (rec.id, transcript_len)

    fh_in.close()
    seqids = []
    for k, v in gene_seqid_len.items():
        seqid, transcript_len = v
        seqids.append(seqid)

    return seqids


def get_seq(fasfile=None, seqids=None, outfh=sys.stdout):
    if fasfile.endswith('.gz'):
        fh_in = gzip.open(fasfile, 'rt')
    else:
        fh_in = open(fasfile)
    for rec in SeqIO.parse(fh_in, 'fasta'):
        if rec.id in seqids:
            SeqIO.write(rec, outfh, 'fasta')
    fh_in.close()


def main():
    usage = '''
    python3 {0} <in.fasta> <out.fasta>
    '''.format(sys.argv[0])

    if len(sys.argv) != 3:
        sys.exit(usage)

    in_fasta, out_fasta = sys.argv[1:3]

    seqids = get_longest_transcript(fasfile=in_fasta)

    with open(out_fasta, 'w') as fhout:
        get_seq(fasfile=in_fasta, seqids=seqids, outfh=fhout)


if __name__ == '__main__':
    main()









