#!/usr/bin/env python

# ----------------------
# Copyright 2023 PMG Lab
# Author: Lin Miao
# Licence: MIT
# Version: 20230825
# ----------------------

import logging
import argparse
from multiprocessing.dummy import Pool

import numpy as np
import pandas as pd
from scipy.stats import norm
from pysnptools.snpreader import Bed

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logging.info(f'Start 2.realize_beta_and_gwa.')

parser = argparse.ArgumentParser(description='Realize effect sizes and perform association tests')
parser.add_argument('--nt', type=int, required=True, help='Number of threads running in parallel')
parser.add_argument('--h2g-vals', type=str, required=True, help='Target h2g values for the simulation')
parser.add_argument('--pqtl-vals', type=str, required=True, help='Proprotion values of SNPs to be qtlal')
parser.add_argument('--neg-alpha-vals', type=str, required=True, help='Power values in the LDAK-Thin Model')
parser.add_argument('--n-gwa', type=float, required=True, help='Size of a sample for association tests')
parser.add_argument('--n-rep', type=float, required=True, help='Number of repititions')
parser.add_argument('--gene-list', type=str, default=None, help='File of a list of genes to be included')
parser.add_argument('--out-dir', type=str, required=True, help='Directory for output files, same in all steps')

args = parser.parse_args()
n_gwa = int(args.n_gwa)
n_rep = int(args.n_rep)
h2g_list = args.h2g_vals.split(',')
pqtl_list = args.pqtl_vals.split(',')
neg_alpha_list = args.neg_alpha_vals.split(',')
gene_list = np.loadtxt(args.gene_list, dtype=str)


def process_one_gene(gene):
    dir_gene = f'{args.out_dir}/{gene}'

    # Read plink.bim and plink.fam
    plink_bim = pd.read_csv(f'{dir_gene}/plink.bim', sep='\t', header=None)
    plink_bim.columns = ['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']

    # Read genotypes and calculate MAF
    snp_on_disk = Bed(f'{dir_gene}/plink', count_A1=True)
    snpdata = snp_on_disk.read()
    frq_A1 = np.mean(snpdata.val, axis=0) / 2

    # Standardize genotypes
    snpdata.standardize()
    X = snpdata.val
    n_pop, m = X.shape
    tot_n_gwa = n_rep * n_gwa
    X_gwa = X[:tot_n_gwa].reshape(n_rep, n_gwa, m)
    iid_gwa = snpdata.iid[:tot_n_gwa].reshape(n_rep, n_gwa, 2)

    columns = pd.MultiIndex.from_product((pqtl_list, neg_alpha_list, h2g_list))
    realized_beta = pd.DataFrame(dtype=float, index=np.arange(m), columns=columns)
    realized_beta.columns.names = ['prop_qtl', 'neg_alpha', 'target_h2']
    realized_h2_ = pd.DataFrame(dtype=float, index=[gene], columns=columns)
    realized_h2_.columns.names = ['prop_qtl', 'neg_alpha', 'target_h2']

    columns = pd.MultiIndex.from_product((pqtl_list, neg_alpha_list, h2g_list))
    phenotypes = pd.DataFrame(-9, dtype=np.int8, index=snpdata.iid[:, 0], columns=columns)
    phenotypes.columns.names = ['prop_qtl', 'neg_alpha', 'target_h2']
    for pqtl in pqtl_list:
        m_qtl = int(np.ceil(m * float(pqtl)))
        idx_qtl = np.random.choice(m, size=m_qtl, replace=False)
        for neg_alpha in neg_alpha_list:
            # Realize beta
            beta = np.zeros(m)
            qtl_beta_var = (frq_A1[idx_qtl] * (1 - frq_A1[idx_qtl])) ** (1 - float(neg_alpha))
            beta[idx_qtl] = np.random.multivariate_normal(np.zeros(m_qtl), np.diag(qtl_beta_var))  # m_qtl
            genetic_eff = X @ beta
            genetic_eff_std = genetic_eff.std()
            for h2g in h2g_list:
                # Scale beta to fit target h2g
                scale_factor = float(h2g) ** 0.5 / genetic_eff_std
                genetic_eff_scale_h2g = genetic_eff * scale_factor
                # Calculate phenotypes and realized h2g
                y_norm = genetic_eff_scale_h2g + np.random.normal(0, np.sqrt(1 - float(h2g)), n_pop)
                realized_beta[(pqtl, neg_alpha, h2g)] = beta * scale_factor
                realized_h2_[(pqtl, neg_alpha, h2g)] = genetic_eff_scale_h2g.var() / y_norm.var()

                phenotypes[(pqtl, neg_alpha, h2g)] = y_norm
                # Perform association tests
                y_test = y_norm[:tot_n_gwa].reshape(n_rep, n_gwa, 1)
                z = np.swapaxes(X_gwa, 1, 2) @ y_test / np.sqrt(n_gwa)  # n_rep * m * 1
                p = norm.sf(np.abs(z)) * 2

                for j in range(n_rep):
                    out_pre = f'{dir_gene}/rep{j}/pqtl{pqtl}_alpha-{neg_alpha}_h2g{h2g}'

                    # Write phenotypes of the sample
                    pd.DataFrame(y_test[j], index=pd.MultiIndex.from_arrays(iid_gwa[j].T)).to_csv(
                        f'{out_pre}.pheno', sep=' ', header=False)

                    # Write a summary statistic file for KGGSEE, HESS, and LDSC
                    sumstat = plink_bim.copy()
                    sumstat['Z'] = z[j]
                    sumstat['N'] = n_gwa
                    sumstat['P'] = p[j]
                    sumstat[['CHR', 'BP', 'P', 'SNP', 'A1', 'A2', 'Z', 'N']] \
                        .to_csv(f'{out_pre}.sumstat.gz', sep='\t', index=False)

                    # Write a summary statistic file for LDER
                    sumstat.rename({'SNP': 'snp', 'CHR': 'chr', 'A1': 'a0', 'A2': 'a1', 'Z': 'z'}, axis=1)[
                        ['snp', 'chr', 'a0', 'a1', 'z']].to_csv(
                        f'{out_pre}.lder.sumstat.gz', sep='\t', index=False)

                    # Write a summary statistic file for LDAK
                    sumstat.rename({'SNP': 'Predictor', 'N': 'n'}, axis=1)[['Predictor', 'A1', 'A2', 'n', 'Z']] \
                        .to_csv(f'{out_pre}.ldak.sumstat', sep='\t', index=False)

    realized_beta.to_csv(f'{dir_gene}/realized_effect_sizes.tsv', sep='\t')  # just for debug
    phenotypes.to_csv(f'{dir_gene}/realized_phenotypes.tsv', sep='\t')  # just for debug

    logging.info(f'Done {gene}.')
    return realized_h2_


logging.info(f'Start simulating phenotypes and perform association tests in {args.nt} threads')
realized_h2 = pd.concat(Pool(args.nt).map(process_one_gene, gene_list))
realized_h2.to_csv(f'{args.out_dir}/realized_h2.tsv', sep='\t')
logging.info(f'Done 2.realize_beta_and_gwa.')
