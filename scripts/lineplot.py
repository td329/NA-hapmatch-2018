from __future__ import division
from __future__ import print_function

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import itertools


for pop, ID in itertools.product(['Pima', 'Surui'], [1,2,3,4]):
#for pop, ID in itertools.product(['Pima'], [1]):
    print(pop, ID)
    scores = {}
    lengths = {}
    for i in range(2000, 4010, 10):
        scores[i] = []
        lengths[i] = []
        with open('../data/{p}_anzick_CK-13.out/{p}_anzick_CK-13.out.{i}.txt'.format(p=pop, i=i)) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.rstrip().split('\t')

                s = line[ID]
                if s == '.':
                    continue
                else:
                    scores[i].append(int(s))

    anzick = []
    CK = []
    for k in scores:
        anzick.append(len(filter(lambda x: x > 0, scores[k])))
        CK.append(len(filter(lambda x: x < 0, scores[k])))

    A = np.mean(anzick)
    C = np.mean(CK)

    if pop == 'Pima':
        plt.scatter(anzick, CK, color='r', alpha=0.1)
        if ID == 1:
            plt.scatter(A, C, color='r', s=100, marker='x', linewidth=3, label='Pima')
        else:
            plt.scatter(A, C, color='r', marker='x', linewidth=3, s=100)
    else:
        plt.scatter(anzick, CK, color='b', alpha=0.1)
        if ID == 1:
            plt.scatter(A, C, color='b', marker='x', s=100, linewidth=3, label='Surui')
        else:
            plt.scatter(A, C, color='b', marker='x', linewidth=3, s=100)

plt.legend()
plt.xlabel('Anzick-1 chr1 segment matches')
plt.ylabel('CK-13 chr1 segment matches')
plt.show()

# Pa = [5648, 5799, 5712, 5856]
# Pc = [4226, 4122, 4132, 4060]
#
# Sa = [5573, 5438, 5832, 5774]
# Sc = [3935, 4034, 4130, 4067]
#
# plt.scatter(Pa, Pc, color='r')
# plt.scatter(Sa, Sc, color='b')
#
# plt.show()
# 5799 4122
# 5712 4132
# 5856 4060
#
# 5573 3935
# 5438 4034
# 5832 4130
# 5774 4067