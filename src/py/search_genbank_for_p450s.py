#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio import Entrez
from Bio import SeqIO


def main():
    # require email from cli
    parser = argparse.ArgumentParser()
    parser.add_argument('-e',
                        help='email address',
                        metavar='email',
                        type=str,
                        required=True)
    args = parser.parse_args()

    # identify myself to Entrez
    Entrez.email = args.e

    # start the search
    search_term = ('(p450 NOT reductase NOT uncharacterized) '
                   'AND txid7460[Organism] '
                   'AND refseq[filter]')
    with Entrez.esearch(db='protein', term=search_term) as handle:
        records = Entrez.read(handle)
        gi_list = records['IdList']

    # download peptide sequences
    genbank_file = 'output/gb/apis_p450.gb'
    with Entrez.efetch(db='protein',
                       id=gi_list,
                       rettype='gb',
                       retmode='text') as handle:
        gb_records = SeqIO.parse(handle, 'gb')
        # write gb to disk
        SeqIO.write(gb_records, genbank_file, 'gb')

    # convert gb to fa
    fa_file = 'output/fa/apis_p450.fa'
    gb_records = SeqIO.parse(genbank_file, 'gb')
    SeqIO.write(gb_records, fa_file, 'fasta')

if __name__ == '__main__':
    main()
