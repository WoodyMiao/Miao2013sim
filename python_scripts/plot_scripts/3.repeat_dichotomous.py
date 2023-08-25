#!/usr/bin/env python

# ----------------------
# Copyright 2023 PMG Lab
# Author: Lin Miao
# Licence: MIT
# Version: 20230825
# ----------------------

import argparse
import pandas as pd
from settings import *
from itertools import product

parser = argparse.ArgumentParser(description='Plot results of 2.repeat_quantitative_gcta.sh')
parser.add_argument('--prefix', type=str, required=True, help='Prefix of the inputs and outputs')
parser.add_argument('--h2g', type=str, required=True, help='Target h2g for the simulation')
parser.add_argument('--neg-alpha', type=str, required=True, help='Power value in the LDAK-Thin Model')
parser.add_argument('--pqtl-vals', type=str, required=True, help='Proprotion values of SNPs to be qtlal')
parser.add_argument('--prevalence-vals', type=str, default=None, help='prevalence values of binary traits')
parser.add_argument('--binary-assoc-test', type=str, required=True, help='chi2|linear|logit, the association test '
                                                                         'performed for a dichotomous phenotype')

args = parser.parse_args()
h2f = float(args.h2g)
pqtl_list = args.pqtl_vals.split(',')
prvl_list = args.prevalence_vals.split(',')
test = args.binary_assoc_test

R = R[1:]
C = C[1:]
legend_label = legend_label[1:]
w = 1 / len(R) - 0.15
s = np.arange(len(R)) * w * 1.2  # shift
s -= s.mean()
bp_style['widths'] = w

results = {q: {prvl: {} for prvl in prvl_list} for q in quantiles}
for q, prvl, pqtl in product(quantiles, prvl_list, pqtl_list):
    results[q][prvl][pqtl] = pd.read_csv(
        f'{args.prefix}.prvl{prvl}_pqtl{pqtl}_alpha-{args.neg_alpha}_h2g{args.h2g}.{test}.one4th_{q}.tsv',
        sep='\t', header=[0, 1, 2], index_col=0)

q_prvl_tup_list = list(product(quantiles, prvl_list))
nrow = len(q_prvl_tup_list)
ncol = len(pqtl_list)
abc = np.array(list('ABCDEFGHIJKLMNOPQRSTUVWX')).reshape(nrow, ncol)

for est, ylabel in est_ylabel.items():
    print(f'Plotting h2={args.h2g}, {est} ...')
    if est == 'se_mrb':
        figsize = (8, 12)
        prog_plot = prog_local_h2
        text_bbox['x'] = 0.03
        text_bbox['y'] = 0.95
    else:
        figsize = (12, 16)
        prog_plot = programs

    fig, ax = plt.subplots(nrow, ncol, figsize=figsize, sharex=True, sharey=True)
    fig.subplots_adjust(left=0.08, bottom=0.03, right=0.99, top=0.99, wspace=0, hspace=0)
    x = np.arange(len(prog_plot))
    bp = dict()
    for i, (q, prvl) in enumerate(q_prvl_tup_list):
        for j, pqtl in enumerate(pqtl_list):
            df = results[q][prvl][pqtl]

            for u in [0, 1]:
                bp[u] = ax[i, j].boxplot(df.loc[:, (prog_plot, R[u], est)], positions=x + s[u],
                                         boxprops=dict(facecolor=C[u], lw=lw), **bp_style)
                if df.loc[:, ('LDER', R[u], est)].isnull().sum() > 0:
                    ax[i, j].boxplot(df.loc[:, ('LDER', R[u], est)].dropna(), positions=x[-1:] + s[u],
                                     boxprops=dict(facecolor=C[u], lw=lw), **bp_style)
                if df.loc[:, ('GBAT', R[u], est)].isnull().sum() > 0:
                    ax[i, j].boxplot(df.loc[:, ('GBAT', R[u], est)].dropna(), positions=x[-1:] + s[u],
                                     boxprops=dict(facecolor=C[u], lw=lw), **bp_style)

            if i == 0 and j == 0:
                ax[i, j].text(s=f'{abc[i, j]}. Pr(QTL)$={pqtl}, α=-0.25$\n'
                                f'$K={prvl}$, ${n_snps[q][0]}≤n$(SNP)$≤{n_snps[q][1]}$',
                              transform=ax[i, j].transAxes, **text_bbox)
            elif i == 0 and j != 0:
                ax[i, j].text(s=f'{abc[i, j]}. Pr(QTL)$={pqtl}, α=-0.25$',
                              transform=ax[i, j].transAxes, **text_bbox)
            elif i != 0 and j == 0:
                ax[i, j].text(s=f'{abc[i, j]}. $K={prvl}$, ${n_snps[q][0]}≤n$(SNP)$≤{n_snps[q][1]}$',
                              transform=ax[i, j].transAxes, **text_bbox)
            else:
                ax[i, j].text(s=f'{abc[i, j]}', transform=ax[i, j].transAxes, **text_bbox)

            if i == nrow - 1:
                ax[i, j].set_xticks(x)
                ax[i, j].set_xticklabels(prog_plot)

    for i in range(nrow):
        for j in range(ncol):
            if est == 'h2_mean':
                ax[i, j].hlines(h2f, xmin=-0.5, xmax=x[-1] + .5, colors='gray', linestyles='dashed', lw=1)
            if est == 'h2_mrb' or est == 'se_mrb':
                ax[i, j].hlines(0, xmin=-0.5, xmax=x[-1] + .5, colors='gray', linestyles='dashed', lw=1)
            if j == 0:
                ax[i, j].set_ylabel(ylabel)
                if est == 'h2_mean' or est == 'h2_sd':
                    if est == 'h2_mean':
                        ax[0, 0].set_ylim(-h2f / 3, h2f * 2.9)
                    elif est == 'h2_sd':
                        ax[0, 0].set_ylim(0, h2f * 1.99)
                    yticks = ax[i, j].get_yticks()
                    yticklabels = [f'{a * 100:.2f}%' for a in yticks]
                    ax[i, j].set_yticks(yticks)
                    ax[i, j].set_yticklabels(yticklabels)

    if est == 'h2_mean':
        ax[0, 0].set_ylim(-h2f / 3, h2f * 2.9)
    elif est == 'h2_sd':
        ax[0, 0].set_ylim(0, h2f * 1.99)
    elif est == 'h2_mrb':
        ax[0, 0].set_ylim(-1.4, 1.9)
    elif est == 'se_mrb':
        ax[0, 0].set_ylim(-1.2, 2.7)

    ax[0, 2].legend([bp[0]['boxes'][0], bp[1]['boxes'][0]], legend_label, title=legend_title,
                    loc='upper right', labelspacing=0.2)
    fig.savefig(f'{args.prefix}.alpha-{args.neg_alpha}_h2g{args.h2g}.{test}.{est}.png')
    plt.close()
