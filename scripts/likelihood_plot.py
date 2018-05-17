#!/usr/bin/env python
#
# Likelihood surface plot using model described in method.pdf
#
# ===========================================================

from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import math
import random
import scipy

# ===========================================================

exp, log, beta_ln = np.exp, np.log, scipy.special.betaln

# ===========================================================

NUM_ITERATIONS = 60         # Num of random samples drawn
RES = 20                    # Num blocks in heat map
GEN_TIME = 25               # Generations in years
MERGE = 15000/GEN_TIME      # Global merge time
J1_START = 12600/GEN_TIME   # Anzick
J2_START = 4630/GEN_TIME    # NorthEast
MAX_LENGTH = 200000         # Exclude the few very long segments
BINS = 30                   # Num bins
IND = None
PI = math.pi
ALPHA = 10                  # Parameter in beta distribution
USE_BETA = True             # Uses binomial if False
CONVERGENCE_PLOT = True     # Plots difference between random iterates if True

MODEL_ARGS = {'Ne': 1e4,
              'merge': MERGE,
              'scale_1': 1.5,
              'scale_2': 1.5,
              'mu': 1.2e-8}

# ===========================================================


def extract_random_sample(ind, pop, show_data=False):
    """
    Method for importing and formatting pymatch output
        - MCMC iterate randomly chosen
        - Data extracted and binned by length using segment_data()

    :param ind: choice of haplotype (1,2,3, or 4)
    :param pop: "Pima" or "Surui"
    :param show_data: if True displays data for testing purposes
    :return: data subdivided into bins by length
    :rtype: dict

    """
    data = []
    sample = random.choice(range(2000, 4010, 10))

    with open('../data/{p}_anzick_CK-13.out/{p}_anzick_CK-13.out.{x}.txt'.format(p=pop, x=sample)) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.rstrip().split('\t')
            pos = line[0][1:-1].split(', ')
            length = int(pos[1]) - int(pos[0]) + 1  # Endpoints included
            score = line[ind]
            if score != '.':
                data.append((length, int(score)))

    data = segment_data(data)

    if show_data:
        for k, val in data.items():
            print(k, val)

    return data


def segment_data(data):
    """
    Bins list of segment-score tuples by length
    Representative length chosen to be midpoint of binning interval

    :param data: list of (length, score) segment tuples
    :return: segment scores binned by length
    :rtype: dict

    """
    bins = list(np.linspace(100, MAX_LENGTH, BINS))
    data = sorted(data, key=lambda entry: entry[0])
    d_it = iter(data)
    d = d_it.next()
    blocks = []
    out = {}
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
        dat = blocks[i]
        tmp = [j for j in dat if j > 0]
        obs = (len(tmp), len(dat) - len(tmp))
        if i == 0:
            l = length // 2
        else:
            l = (length + bins[i - 1]) // 2
        out[l] = obs

    return out


# ===========================================================


def compare_theoretical_pop_matches(args=MODEL_ARGS, j1=None, j2=None, l=None, pop=None):
    """
    Determines segment matching probabilities given model params
    For derivation of probabilities, see method.pdf

    :param args: dict of shared model params
    :param j1: join time population 1
    :param j2: join time population 2
    :param l: length of segment (given by data key)
    :param pop: "Pima" or "Surui"
    :return: probability of segment matching
    :rtype: float

    """
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
            p0, p1, p2 = 0, 0.5, 0.5
        elif pop == 'Surui':
            p0, p1, p2 = 0, 0.5, 0.5

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


# ===========================================================


def log_stirling(n, k):
    """
    :return: log stirling approximation of "n choose k"
    :rtype: float

    """
    if k == n or k == 0:
        lgs = 0
    else:
        lgs = (1 / 2) * (log(n) - log(2 * PI * k * (n - k)))
        lgs += n * log(n) - k * log(k) - (n - k) * log(n - k)

    return lgs


def binomial_log_likelihood(p, val):
    """
    :param val: tuple (k,j) of observed matches with both populations
    :return: binomial log likelihood of observing val given binomial parameter p
    :rtype: float

    """
    k = val[0]
    n = val[0] + val[1]
    return log_stirling(n, k) + k * log(p) + (n - k) * log(1 - p)


def beta_log_likelihood(p, val):
    """
    :param val: tuple (k,j) of observed matches with both populations
    :return: beta log likelihood of val given parameter p and global ALPHA
    :rtype: float

    """
    k = val[0]
    n = val[0] + val[1]

    a = ALPHA
    b = ALPHA*(1-p)/p

    return log_stirling(n, k) + beta_ln(k + a, n - k + b) - beta_ln(a, b)


def int_round(vls):
    """
    :param vls: list of ints to be rounded down
    :return: rounded down list of ints
    """
    vls = list(vls)
    return [int(round(v)) for v in vls]


def update(T, obs, pop):
    """
    Update likelihood surface for each randomly chosen sample

    :param T: pd.DataFrame consisting of current set of log likelihoods
    :param obs:
    :param pop:
    :return: updated likelihoods
    :rtype: pd.DataFrame
    """
    j1_vls = np.linspace(J1_START, MERGE, RES)
    j2_vls = np.linspace(J2_START, MERGE, RES)

    for idx1, row in enumerate(T):
        join_1 = j1_vls[idx1]
        out = []
        for idx2, join_2 in enumerate(j2_vls):
            prob = row[idx2]

            for b, val in obs.items():
                p = compare_theoretical_pop_matches(j1=join_1, j2=join_2, l=b, pop=pop)
                if USE_BETA:
                    prob += beta_log_likelihood(p, val)
                else:
                    prob += binomial_log_likelihood(p, val)
            out.append(prob)

        T[idx1] = out

    return T


# ===========================================================


def plot_heatmap(T):
    """
    :param T: pd.DataFrame of log likelihoods for plotting
    :return: None

    """
    j1_vls = np.linspace(J1_START, MERGE, RES)
    j2_vls = np.linspace(J2_START, MERGE, RES)

    T = pd.DataFrame(T, index=int_round(j1_vls*GEN_TIME), columns=int_round(j2_vls*GEN_TIME))
    axes = sns.heatmap(T, cmap='Blues', cbar_kws={'label': 'Average log likelihood'})
    axes.set_xlabel("CK-13 subpop. join time (years before present)")
    axes.set_ylabel("Anzick subpop. join time (years before present)")
    axes.invert_yaxis()
    axes.set_xticklabels(axes.get_xticklabels(), rotation=45)
    axes.set_yticklabels(axes.get_yticklabels(), rotation=0)
    plt.show()


# ===========================================================


if __name__ == "__main__":

    T = np.zeros((RES, RES))
    dfs = []

    # For each pop and haplotype a random sample chosen
    # That sample is used to update the sample of log likelihoods T
    for j in range(1, NUM_ITERATIONS):
        test = T.copy()
        for pop in ['Surui', 'Pima']:
            for hap in [1, 2, 3, 4]:
                observations = extract_random_sample(hap, pop)
                T = update(T, observations, pop)

        if j > 1:
            # Store absolute point-wise difference between iterates
            # Sum each round for convergence plot
            df = abs(test/(j-1) - T/j)
            dfs.append(np.sum(df)/(RES**2))

        if j % 10 == 0:
            print('Iterate {}/{} completed'.format(j, NUM_ITERATIONS))

    T /= NUM_ITERATIONS

    print(T.max(), T.min())

    if CONVERGENCE_PLOT:
        plt.plot(range(2, NUM_ITERATIONS), dfs)
        plt.xlabel('Number of iterations')
        plt.ylabel('Average absolute difference in likelihood between iterates')
        plt.show()

    plot_heatmap(T)

    # =========================================================