#!/usr/bin/env python

# ----------------------
# Copyright 2023 PMG Lab
# Author: Lin Miao
# Licence: MIT
# Version: 20230825
# ----------------------

import os
import re
import logging
import numpy as np
import pandas as pd

path_of_the_bed_file = 'ref_haplotypes/chr1.codegene_10kbflk.bed'
path_of_the_id_pop_map_file = 'ref_haplotypes/integrated_call_samples_v3.20130502.ALL.panel'
path_of_the_input_vcf_file = 'ref_haplotypes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

logging.info('Read the individual-population map')
pop = pd.read_csv(path_of_the_id_pop_map_file, sep='\t')

logging.info('Read the BED file')
gene_bed = np.loadtxt(path_of_the_bed_file, dtype=int, usecols=[0, 1, 2])
n_genes = gene_bed.shape[0]

logging.info(f'Read the input VCF file.')
meta_info = list()
with os.popen(f'zcat {path_of_the_input_vcf_file}') as I:
    for line in I:
        if re.match(f'##(?:fileformat|FILTER|fileDate|reference|source|contig=<ID=1,|INFO=<ID=EUR_AF)', line):
            meta_info.append(line)
        if re.match('#CHROM', line):
            break
    vcf = pd.read_csv(I, sep='\t', comment='#', header=None)
    vcf.columns = line.strip().split('\t')

logging.info(f'Replace the INFO column with allele frequencies only')
vcf.INFO = vcf.INFO.str.split(';', expand=True)[8]

logging.info('Extract the EUR panel')
vcf = vcf[pd.concat([vcf.columns[:9].to_series(), pop.loc[pop.super_pop == 'EUR', 'sample']])]

logging.info('Drop loci with duplicated positions or rsIDs')
vcf = vcf[(~vcf[['#CHROM', 'POS']].duplicated(keep=False)) &
          (~vcf['ID'].duplicated(keep=False)) & vcf['ID'].str.match(r'^rs\d+$')]

logging.info('Extract only biallelic SNPs')
vcf = vcf[vcf['REF'].isin(['A', 'C', 'G', 'T']) & vcf['ALT'].isin(['A', 'C', 'G', 'T']) &
          np.all(vcf.iloc[:, 9:].isin(['0|0', '0|1', '1|0', '1|1']), axis=1)]

logging.info('Extract SNPs with MAF > 0.01')
allele_frq = vcf.INFO.str.split('=', expand=True)[1].astype(float)
vcf = vcf[(0.01 < allele_frq) & (allele_frq < 0.99)]

logging.info('Extract SNPs in the regions defined by the BED file')
idx_arrays = list()
snp_counts = np.empty(n_genes, dtype=int)
for i in range(n_genes):
    idx = np.where((vcf['POS'] >= gene_bed[i, 1]) & (vcf['POS'] <= gene_bed[i, 2]))[0]
    idx_arrays.append(idx)
    snp_counts[i] = idx.shape[0]
vcf = vcf.iloc[np.sort(np.unique(np.concatenate(idx_arrays)))]

logging.info(f'Write the output VCF file.')
with os.popen(f'gzip -9c >ref_haplotypes/chr1.codegene_10kbflk.vcf.gz', 'w') as O:
    O.writelines(meta_info)
    vcf.to_csv(O, index=False, sep='\t')

logging.info(f'Count SNPs and write SNP counts.')
gene_bed = pd.DataFrame(gene_bed, columns=['CHR', 'START', 'END'])
gene_bed['GENE'] = [f'gene{a + 1:04d}' for a in range(gene_bed.shape[0])]

gene_bed['nSNP'] = snp_counts
q0 = gene_bed['nSNP'].quantile(0)
q1 = gene_bed['nSNP'].quantile(0.25)
q2 = gene_bed['nSNP'].quantile(0.5)
q3 = gene_bed['nSNP'].quantile(0.75)
q4 = gene_bed['nSNP'].quantile(1)

genes_q0_q1 = gene_bed[(q0 < snp_counts) & (snp_counts <= q1)].iloc[::4]
genes_q1_q2 = gene_bed[(q1 < snp_counts) & (snp_counts <= q2)].iloc[::4]
genes_q2_q3 = gene_bed[(q2 < snp_counts) & (snp_counts <= q3)].iloc[::4]
genes_q3_q4 = gene_bed[(q3 < snp_counts) & (snp_counts <= q4)].iloc[::4]
for name_q, df_q in zip(['q0-q1', 'q1-q2', 'q2-q3', 'q3-q4'], [genes_q0_q1, genes_q1_q2, genes_q2_q3, genes_q3_q4]):
    df_q['GENE'].to_csv(f'ref_haplotypes/chr1.one4th_{name_q}.lst', header=False, index=False)

one4th = pd.concat([genes_q0_q1, genes_q1_q2, genes_q2_q3, genes_q3_q4])
one4th['GENE'].to_csv(f'ref_haplotypes/chr1.one4th.lst', header=False, index=False)
one4th.loc[::8, 'GENE'].to_csv(f'ref_haplotypes/chr1.one32th.lst', header=False, index=False)
gene_bed.to_csv('ref_haplotypes/chr1.codegene_10kbflk.snp_counts.tsv', sep='\t', index=False)
