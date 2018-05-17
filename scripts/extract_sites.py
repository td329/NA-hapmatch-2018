#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse
import gzip


# ===========================================================

END_OF_CHR = 249250621

# ===========================================================


def vcf_to_sites():
    # print(">POS", "SURUI_1-0", "SURUI_1-1", "SURUI_2-0", "SURUI_2-1", sep='\t')

    gen = slice(-2, None)

    for l in sys.stdin:
        alleles = []

        if l.startswith('#'):
            continue

        l = l.rstrip().split('\t')
        POS, REF, ALT, FILTER = l[1], l[3], l[4], l[6]

        if FILTER != 'PASS':
            print(POS, *['.']*4, sep='\t')
            continue

        for i in l[gen]:
            if '.' in i or '/' in i:
                alleles.extend(['.']*2)
            else:
                types = [REF if _ == '0' else ALT for _ in i.split('|')]
                alleles.extend(types)

        print(POS, *alleles, sep='\t')


def fa_to_sites():
    count = 1
    for l in sys.stdin:
        if l.startswith('>1'):
            continue

        if l.startswith('>2') or l.startswith('\n'):
            break

        for a in l.rstrip():
            print(count, a, sep='\t')
            count += 1


def fa_mask_to_bed():
    # extracting chr1 bed mask
    start, pos = 0, 0
    previous = '*'

    for line in sys.stdin:
        if line.startswith(">1"):
            continue
        elif line.startswith("\n") or line.startswith(">2"):
            sys.exit()

        line = line.rstrip()

        for current in line:
            pos += 1
            if previous + current == '01':
                print('chr1', start, pos, sep='\t')
            elif previous + current == '10':
                start = pos

            previous = current


# ===========================================================


def get_source_iterator(file):
    is_vcf = False

    if file.endswith('vcf.gz'):
        is_vcf = True
        f = gzip.open(file, 'r')
    else:
        f = open(file, 'r')

    i = iter(f)
    line = i.next()

    while line.startswith('#') or line.startswith('>'):
            line = i.next()

    while True:
        d = line.rstrip().split('\t')
        if is_vcf:
            als = {'.': '.', '0': d[3], '1': d[4]}
            anz = als[d[-3].split('/')[0]]
            comparison_al = als[d[-1].split('/')[0]]
            yield [d[1], anz, comparison_al]
        else:
            yield d
        line = i.next()

    f.close()


def combine_sites(args):
    mod, ref, anc, ASO, msk = map(get_source_iterator, [args.moderns, args.reference, args.ancestral, args.ASO, args.mask])
    iterators = [mod, ref, anc, ASO, msk]
    flag = {_: False for _ in iterators}
    data = {_: _.next() for _ in iterators}
    order = {mod: slice(0, 4), ref: slice(4, 5), anc: slice(5, 6), ASO: slice(6, None)}

    missing_chars = ['.', '-', 'N']

    i = 17000   # start of chrs missing

    while i <= END_OF_CHR:
        filter_pass = None
        alleles = ['R']*4 + ['.']*4
        i += 1

        if any([flag[key] for key in [ref, anc, ASO]]):
            break

        for it in iterators:
            while int(data[it][0]) <= i and not flag[it]:
                if int(data[it][0]) == i:
                    if it == msk:
                        filter_pass = bool(data[msk][1])
                    else:
                        alleles[order[it]] = data[it][1:]
                try:
                    data[it] = it.next()
                except StopIteration:
                    flag[it] = True

        # Ignoring sites which don't pass Simons x75.fa filter
        if not filter_pass:
            continue

        # Ignoring sites in which no modern has derived allele
        if alleles[:4].count('R') + alleles[:4].count('.') == 4 and alleles[4] == alleles[5].upper():
            continue

        # Ignoring sites where ref, ancestrals or ancients are missing
        # Stronger filtering would ignore segment on which any site is missing
        if any([_ in missing_chars for _ in alleles[4:]]):
            continue

        # Ignoring sites where ancients are the same
        if alleles[-2] == alleles[-1]:
            continue

        print(i, *alleles, sep='\t')

# ===========================================================


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('-m', '--moderns')
    p.add_argument('-r', '--reference')
    p.add_argument('-anc', '--ancestral')
    p.add_argument('-aso', '--ASO')
    p.add_argument('-mask', '--mask')
    args = p.parse_args()

    # print('>POS', 'SURUI-1_0', 'SURUI-1_1', 'SURUI-2_0',\
    # 'SURUI-2_1', 'REFERENCE', 'ANCESTRAL', 'ANZICK', 'RM-83', sep='\t')

    print('>POS', 'PIMA-1_0', 'PIMA-1_1', 'PIMA-2_0', 'PIMA-2_1', 'REFERENCE', 'ANCESTRAL', 'ANZICK', 'CK-13', sep='\t')
    combine_sites(args)

# ===========================================================