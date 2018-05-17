from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import math
import random
import sys
import scipy
import itertools


# ===========================================================


exp, log = np.exp, np.log
NUM_ITERATIONS = 20
POST_NORMALISATION = True
RES = 20
MERGE = 800
J1_START = 500      # Anzick
J2_START = 150      # NorthEast
MAX_LENGTH = 50000
IND = None
PI = math.pi
ALPHA = 4
USE_BETA = False

MODEL_ARGS = {'Ne': 0.5e4,
              'join_1': MERGE,
              'join_2': MERGE,
              'merge': MERGE,
              'ancestry_1': 0.5,
              'ancestry_2': 0.5,
              'scale_1': 1.5,
              'scale_2': 1,
              'mu': 1.2e-8,
              'length': 2500}


# ===========================================================


def compare_theoretical_pop_matches(args=MODEL_ARGS, j1=None, j2=None, l=None, pop=None):
    if j1:
        args['join_1'] = j1
    if j2:
        args['join_2'] = j2
    if l:
        args['length'] = l

    lmb1, lmb2 = 1 / args['scale_1'], 1 / args['scale_2']
    t1, t2 = (args['merge'] - args['join_1']) / (2 * args['Ne']), (args['merge'] - args['join_2']) / (2 * args['Ne'])

    if t1 < 0 or t2 < 0:
        raise Exception('Invalid join times: {}, {}'.format(t1, t2))

    if pop:
        if pop == 'Pima':
            p0, p1, p2 = 0, 1, 0
        elif pop == 'Surui':
            p0, p1, p2 = 0, 1, 0
    else:
        p1 = args['ancestry_1']
        if 'ancestry_2' in args:
            p2 = args['ancestry_2']
            if p1 + p2 > 1:
                raise Exception('Ancestry proportions exceed 1')
            p0 = 1 - p1 - p2
        else:
            p2 = 1 - p1
            p0 = 0

    theta = 4 * args['Ne'] * args['mu'] * args['length']
    B = (p0 + p1 * exp(-t1 / lmb1) + p2 * exp(-t2 / lmb2))
    d2 = 1 + theta / 2

    A = {1: p1 * exp(-t1 / lmb1) / (lmb1 + 1),
         2: p2 * exp(-t2 / lmb2) / (lmb2 + 1)}
    d1 = {1: 1 / lmb1 - theta / 2,
          2: 1 / lmb2 - theta / 2}
    t = {1: t1,
         2: t2}
    I0 = {1: p2 * (1 - exp(-t2 / lmb2)) + (2 / 3) * B,
          2: p1 * (1 - exp(-t1 / lmb1)) + (2 / 3) * B}
    lmb = {1: lmb1,
           2: lmb2}

    out = []
    for i in 1, 2:
        I = I0[i]
        I += A[i] * (exp(d1[i] * t[i]) - 1) / d1[i] + A[i] * (exp(-d2 * t[i]) - 1) / d2
        I += A[i] * (exp((1 + 1 / lmb[i]) * t[i]) - 1) * (exp(-(1 + theta / 2) * t[i])) / d2
        I += B / (3 * d2)
        out.append(1 - I)
    return out[0]/(out[0] + out[1])


def extract_random_sample(ind, pop, show_data=False):
    data = []
    sample = random.choice(range(2000, 4010, 10))

    with open('../data/{p}_anzick_CK-13.out/{p}_anzick_CK-13.out.{x}.txt'.format(p=pop, x=sample)) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.rstrip().split('\t')
            pos = line[0][1:-1].split(', ')
            length = int(pos[1]) - int(pos[0])  # Is there supposed to be a +1 here?
            if length > MAX_LENGTH:
                continue
            score = line[ind]
            if score != '.':
                data.append((length, int(score)))

    data = segment_data(data)

    if show_data:
        for k, val in data.items():
            print(k, val)

    return data


def likelihood(p, val):
    k = val[0]
    n = val[0] + val[1]
    if k == n or k == 0:
        log_stirling = 0
    else:
        log_stirling = (1 / 2) * (log(n) - log(2 * PI * k * (n - k)))
        log_stirling += n * log(n) - k * log(k) - (n - k) * log(n - k)

    return exp(log_stirling + k * log(p) + (n - k) * log(1 - p))


def log_stirling(n, k):
    if k == n or k == 0:
        lgs = 0
    else:
        lgs = (1 / 2) * (log(n) - log(2 * PI * k * (n - k)))
        lgs += n * log(n) - k * log(k) - (n - k) * log(n - k)

    return lgs


def log_likelihood(p, val):
    k = val[0]
    n = val[0] + val[1]
    return log_stirling(n, k) + k * log(p) + (n - k) * log(1 - p)


def log_beta(x, y):
    gml = scipy.special.gammaln
    return gml(x) + gml(y) - gml(x + y)


def beta_log_likelihood(p, val):
    k = val[0]
    n = val[0] + val[1]

    a = ALPHA
    b = ALPHA*(1-p)/p

    return log_stirling(n, k) + log_beta(k + a, n - k + b) - log_beta(a, b)


def segment_data(data):
    bins = list(np.linspace(100, MAX_LENGTH, 30))
    data = sorted(data, key=lambda entry: entry[0])
    d_it = iter(data)
    d = d_it.next()
    blocks = []
    out = {}
    # print(bins)
    for i, b in enumerate(bins):
        keep = []
        stay = True
        while d[0] < b and stay:
            keep.append(d[1])
            try:
                d = d_it.next()
            except StopIteration:
                stay = False
        blocks.append(keep)

    for i, length in enumerate(bins):
        dt = blocks[i]
        tmp = [i for i in dt if i > 0]
        obs = (len(tmp), len(dt) - len(tmp))
        if i == 0:
            l = length // 2
        else:
            try:
                l = (length + bins[i-1]) // 2
            except IndexError:
                print(length, bins, i)
                sys.exit()
        out[l] = obs
        # print(l, out[l])

    return out

if __name__ == "__main__":
    proportions = [(1, 0), (0.75, 0.25), (0.50, 0.50), (0.25, 0.75), (0, 1)]
    j1_vls = np.linspace(J1_START, MERGE, RES)
    j2_vls = np.linspace(J2_START, MERGE, RES)
    merge_vls = [500, 550, 600, 650, 700, 750, 800]

    obs =

    keep = []
    for pars in itertools.product(proportions, j1_vls, j2_vls, merge_vls):
        args = MODEL_ARGS
        args['j1'] = pars[1]
        args['j2'] = pars[2]
        args['p1'] = pars[0][0]
        args['p2'] = pars[0][1]
        args['merge'] = pars[3]

        p = compare_theoretical_pop_matches(args=args)

        obs =




