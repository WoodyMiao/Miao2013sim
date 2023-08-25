#!/usr/bin/env python

# ----------------------
# Copyright 2023 PMG Lab
# Author: Lin Miao
# Licence: MIT
# Version: 20230825
# ----------------------

import argparse
import pandas as pd
from itertools import product
from settings import *

parser = argparse.ArgumentParser(description='Plot results of 2.repeat_quantitative_gcta.sh')
parser.add_argument('--prefix', type=str, required=True, help='Prefix of the inputs and outputs')
parser.add_argument('--h2g-vals', type=str, required=True, help='Target h2g values for the simulation')
parser.add_argument('--neg-alpha-vals', type=str, required=True, help='Power values in the LDAK-Thin Model')
parser.add_argument('--pqtl-vals', type=str, required=True, help='Proprotion values of SNPs to be qtlal')

args = parser.parse_args()
h2g_list = args.h2g_vals.split(',')
neg_alpha_list = args.neg_alpha_vals.split(',')
pqtl_list = args.pqtl_vals.split(',')
par_tup_list = list(product(h2g_list, quantiles, neg_alpha_list, pqtl_list))

results = {h2: {q: {neg_alpha: {} for neg_alpha in neg_alpha_list} for q in quantiles} for h2 in h2g_list}
for h2g, q, neg_alpha, pqtl in par_tup_list:
    results[h2g][q][neg_alpha][pqtl] = pd.read_csv(
        f'{args.prefix}.pqtl{pqtl}_alpha-{neg_alpha}_h2g{h2g}.one4th_{q}.tsv', sep='\t', header=[0, 1, 2], index_col=0)

q_alpha_tup_list = list(product(quantiles, neg_alpha_list))
nrow = len(q_alpha_tup_list)
ncol = len(pqtl_list)
abc = np.array(list('ABCDEFGHIJKLMNOPQRSTUVWX')).reshape(nrow, ncol)

for h2g in h2g_list:
    h2f = float(h2g)
    for est, ylabel in est_ylabel.items():
        print(f'Plotting h2={h2g}, {est} ...')
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
        for i, (q, neg_alpha) in enumerate(q_alpha_tup_list):
            for j, pqtl in enumerate(pqtl_list):
                df = results[h2g][q][neg_alpha][pqtl]

                for u in [0, 1, 2]:
                    bp[u] = ax[i, j].boxplot(df.loc[:, (prog_plot, R[u], est)], positions=x + s[u],
                                             boxprops=dict(facecolor=C[u], lw=lw), **bp_style)
                    if df.loc[:, ('LDER', R[u], est)].isnull().sum() > 0:
                        ax[i, j].boxplot(df.loc[:, ('LDER', R[u], est)].dropna(), positions=x[-1:] + s[u],
                                         boxprops=dict(facecolor=C[u], lw=lw), **bp_style)
                    if df.loc[:, ('GBAT', R[u], est)].isnull().sum() > 0:
                        ax[i, j].boxplot(df.loc[:, ('GBAT', R[u], est)].dropna(), positions=x[-1:] + s[u],
                                         boxprops=dict(facecolor=C[u], lw=lw), **bp_style)

                if i == 0 and j == 0:
                    ax[i, j].text(s=f'{abc[i, j]}. Pr(QTL)$={pqtl}$\n'
                                    f'$α=-{neg_alpha}$, ${n_snps[q][0]}≤n$(SNP)$≤{n_snps[q][1]}$',
                                  transform=ax[i, j].transAxes, **text_bbox)
                elif i == 0 and j != 0:
                    ax[i, j].text(s=f'{abc[i, j]}. Pr(QTL)$={pqtl}$',
                                  transform=ax[i, j].transAxes, **text_bbox)
                elif i != 0 and j == 0:
                    ax[i, j].text(s=f'{abc[i, j]}. $α=-{neg_alpha}$, ${n_snps[q][0]}≤n$(SNP)$≤{n_snps[q][1]}$',
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

        ax[0, 2].legend([bp[0]['boxes'][0], bp[1]['boxes'][0], bp[2]['boxes'][0]], legend_label, title=legend_title,
                        loc='upper right', labelspacing=0.2)
        fig.savefig(f'{args.prefix}.h2g{h2g}.{est}.png')
        plt.close()
