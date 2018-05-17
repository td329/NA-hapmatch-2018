#!/usr/bin/env python
#
# Main segment matching script
# Requires all variant data processed in extract_sites.py
# Requires argweaver output
#
# ===========================================================

from __future__ import print_function
from __future__ import division

# from collections import namedtuple
import argparse
import gzip
import os
import sys
import csv

# ===========================================================

END_OF_CHR = 249250621
OVERLAP = 5e5          # Assumed argweaver was run on 5Mb windows with overlap length OVERLAP

# ===========================================================


def get_positions(filename):
    """
    Opens argweaver output files *smc.gz
    Extracts positions of haplotype segment breakpoints

    :param filename: *smc.gz argweaver output file
    :return: list of (start, end) tuples
    :rtype: list

    """
    positions = []   # few enough rows that this won't cause memory issues
    with gzip.open(filename, 'r') as f:

        f.next()
        regions = f.next().split('\t')
        lower = int(regions[2]) + OVERLAP
        upper = int(regions[3]) - OVERLAP

        for line in f:

            if line.startswith("TREE"):
                fields = line.split("\t")
                start, end = int(fields[1]), int(fields[2])

                if end <= lower or start >= upper:
                    continue
                if start <= 1 or end >= END_OF_CHR:
                    continue

                positions.append((start, end))

    return positions


def test_pos(l):
    """
    Checks for shared private derived alleles at position l[0]

    :param l: line in seqs file
    :return: list of scores for each modern hap allele
    :rtype: list

    """
    scores = [0]*4

    bases = ['A', 'C', 'T', 'G']

    mod = l[1:5]
    ref, aa, anc_pos, anc_neg = l[-4:]
    aa = aa.upper()

    for i, m in enumerate(mod):

        if m == '.':
            continue

        if m == 'R':
            m = ref

        if aa in bases and m != aa:
            if m == anc_pos and m != anc_neg:
                scores[i] = 1
            elif m == anc_neg and m != anc_pos:
                scores[i] = -1

    if scores == [0]*4:
        return None
    else:
        return scores


def get_segment_scores(seqs, positions):
    """
    Counts private derived alleles shared by each modern segment with at most one of the ancient haps

    Input sites file ("seqs") has format:
        > POS    PIMA-1_0    PIMA-1_1    PIMA-2_0    PIMA-2_1    REFERENCE    ANCESTRAL    ANZICK    CK-13
        852875	C	T	C	T	C	C	C	T
        ...
        ...

    :param seqs: name of input site file
    :param positions: tuple of start and end position of segments obtained from argweaver
    :return: position and list of (2) scores
    :rtype: (positions, scores)

    """

    with open(seqs) as s:
        s_it = iter(s)
        s_it.next()
        line = s_it.next().rstrip().split('\t')
        for pos in positions:
            exit_pos = False
            keep = [[] for _ in range(4)]
            seg_score = []

            while int(line[0]) < pos[0]:
                try:
                    line = s_it.next().rstrip().split('\t')
                except StopIteration:
                    exit_pos = True
                    break

            while int(line[0]) <= pos[1]:
                test = test_pos(line)
                if test is not None:
                    for i, t in enumerate(test):
                        keep[i].append(t)
                try:
                    line = s_it.next().rstrip().split('\t')
                except StopIteration:
                    exit_pos = True
                    break

            if exit_pos:
                continue

            for ind_score in keep:
                if not any(ind_score):
                    seg_score.append('.')
                elif 1 in ind_score and -1 in ind_score:    # infinite sites violation
                    seg_score.append('.')
                else:
                    seg_score.append(sum(ind_score))

            yield (pos[0], pos[1]), seg_score


# ===========================================================


p = argparse.ArgumentParser()
p.add_argument('-l', '--lfn', help='list of file names of arg output')
p.add_argument('-f', '--filename', help='file name of arg output, usually directory/out.*.smc.gz')
p.add_argument('-s', '--sequences', help='sequences combined in sites file')
args = p.parse_args()


# ===========================================================


if __name__ == "__main__":

    POP = 'Pima'

    if args.filename and args.lfn:
        sys.exit('Supply only one of filename and lofn')

    if not os.path.exists('{p}_anzick_CK-13.out'.format(p=POP)):
        os.makedirs('{p}_anzick_CK-13.out'.format(p=POP))

    for sample in range(2000, 4010, 10):
        out_name = "{p}_anzick_CK-13.out/{p}_anzick_CK-13.out.{s}.txt".format(p=POP, s=sample)
        print(out_name)

        with open(args.lfn) as lfn:
            for l in lfn:
                arg_file = "/home/td329/projects/NA/S_{p}-1-2_chr1/{l}/out.{s}.smc.gz"\
                    .format(p=POP, l=l.rstrip(), s=sample)
                if not os.path.isfile(arg_file):
                    continue

                with open(out_name, 'a') as out:
                    print('writing:', out_name)
                    out.write("#NAME {}\n".format(arg_file))
                    for sc in get_segment_scores(args.sequences, get_positions(arg_file)):
                        if sc[1] == ['.']*4:
                            continue
                        else:
                            w = csv.writer(out, delimiter='\t')
                            w.writerow([sc[0]] + sc[1])


# =========================================================== 