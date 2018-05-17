from __future__ import print_function
import sys
import argparse
import gzip

# ['PIMA-1_0', 'PIMA-1_1', 'PIMA-2_0', 'PIMA-2_1', 'REFERENCE', 'ANCESTRAL', 'ANZICK', 'CK-13']

anz = [0]*4
ck = [0]*4

with open('../data/Pima_ref_anc_anzick_CK13_x75-mask_chr1.sites', 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')[1:]
        for i, val in enumerate(line[0:4]):
            if val == '.':
                continue
            elif val == 'R' and line[5]:
                if line[4] == line[6]:
                    anz[i] += 1
                elif line[4] == line[7]:
                    ck[i] += 1
            elif val == line[5]:
                if val == line[6]:
                    anz[i] += 1
                elif val == line[7]:
                    ck[i] += 1

for i, j in zip(anz, ck):
    print(i, j)

print()


anz = [0]*4
ck = [0]*4

with open('../data/Surui_ref_anc_anzick_CK13_x75-mask_chr1.sites', 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')[1:]
        for i, val in enumerate(line[0:4]):
            if val == '.':
                continue
            elif val == 'R' and line[5]:
                if line[4] == line[6]:
                    anz[i] += 1
                elif line[4] == line[7]:
                    ck[i] += 1
            elif val == line[5]:
                if val == line[6]:
                    anz[i] += 1
                elif val == line[7]:
                    ck[i] += 1

for i, j in zip(anz, ck):
    print(i, j)

# output shared private derived alleles:
# anzick ck13
#
# 5648 4226
# 5799 4122
# 5712 4132
# 5856 4060
#
# 5573 3935
# 5438 4034
# 5832 4130
# 5774 4067






# count_difs, count_same = 0, 0
#
# with gzip.open('../data/CK-13.chr1.vcf.gz', 'r') as C, gzip.open('../data/ASO.chr1.vcf.gz', 'r') as A:
#     c_it = iter(C)
#     c = c_it.next()
#     # a_it = iter(A)
#     # a = a_it.next()
#
#     print('>POS\tANZICK\tCK-13')
#
#
#     for a in A:
#         if a.startswith('#'):
#             continue
#
#         pos = int(a.split('\t')[1])
#
#         while c.startswith('#'):
#             c = c_it.next()
#
#         if pos > int(c.split('\t')[1]):
#             c = c_it.next()
#         elif pos == int(c.split('\t')[1]):
#             c_out = c.split('\t')[3]
#             ref_a = a.split('\t')[3]
#             alt_a = a.split('\t')[4]
#             al_anz = a.rstrip().split('\t')[-3]
#
#             if al_anz == '1/1':
#                 anz_out = alt_a
#             else:
#                 anz_out = ref_a
#
#             if c_out != anz_out:
#                 print(pos, anz_out, c_out, sep='\t')
#                 count_difs += 1
#             else:
#                 count_same += 1
#
#             try:
#                 c = c_it.next()
#             except StopIteration:
#                 break

