#!/usr/bin/env python

import argparse
import itertools
import random
from pathlib import Path
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio.Seq import Seq


def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--table", type=int, default=11, help="Translation table")
    parser.add_argument("fasta", type=Path)

    return parser.parse_args()


def main():
    args = arguments()

    wobbly_table = create_ambiguous_back_table(args.table)

    with args.fasta.open("r") as f:
        for record in SeqIO.parse(f, "fasta"):

            reverse_translated = wobble_gene(record.seq, wobbly_table, args.table)
            out_record = SeqRecord(
                seq=Seq(reverse_translated),
                id=record.name,
                description="wobbled_third_codon",
            )

            with sys.stdout as sout:
                SeqIO.write(out_record, handle=sout, format="fasta")


def create_ambiguous_back_table(table: int):

    wobbly_table = {}

    t = CodonTable.unambiguous_dna_by_id[table]

    for residue in t.protein_alphabet:
        wobbly_table[residue] = [t.back_table.get(residue)]

    for codon in map("".join, itertools.product("ATCG", repeat=3)):

        aa = t.forward_table.get(codon)
        if aa is None:
            continue
        if codon not in wobbly_table[aa]:
            wobbly_table[aa].append(codon)

    return wobbly_table


def wobble_gene(sequence: Seq, backtable: dict, table_number: int):
    def pick_random_codon(aa):
        if aa == "*":
            return CodonTable.unambiguous_dna_by_id[table_number].stop_codons[0]
        return random.choice(backtable[aa])

    aa_seq = sequence.translate(table=table_number)

    reverse_translated = "".join(map(pick_random_codon, aa_seq))

    return reverse_translated


if __name__ == "__main__":
    main()
