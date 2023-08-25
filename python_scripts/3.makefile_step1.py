#!/usr/bin/env python

# ----------------------
# Copyright 2023 PMG Lab
# Author: Lin Miao
# Licence: MIT
# Version: 20230825
# ----------------------

import os
import logging
import argparse
import numpy as np

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logging.info(f'Start 3.makefile_step1.')

# Set the executables
gcta_exe = 'gcta --threads 1'
ldak_exe = 'ldak'
lder_plinkLD = 'plinkLD.py --thread 1'
ldsc_exe = 'ldsc.py'
plink_exe = 'plink --threads 1'
kggsee_exe = 'java -Xmx4g -jar /app/pmglab/kggsee/kggsee.jar --nt 1'

parser = argparse.ArgumentParser(description='Make shell scripts of h2 estamating programs')
parser.add_argument('--n-ld', type=str, required=True, help='Sizes of samples for LD panels')
parser.add_argument('--n-rep', type=float, required=True, help='Number of repititions')
parser.add_argument('--skip-mkdir', action='store_true', default=False, help='Skip commands of mkdir')
parser.add_argument('--skip-gcta', action='store_true', default=False, help='Skip commands for GCTA')
parser.add_argument('--skip-assoc-ld', action='store_true', default=False, help='Skip commands using LD coefficient '
                                                                                'matrix calculated from samples of '
                                                                                'association tests')
parser.add_argument('--gene-list', type=str, default=None, help='File of a list of genes to be included')
parser.add_argument('--out-dir', type=str, required=True, help='Directory for output files, same in all steps')

args = parser.parse_args()
n_rep = int(args.n_rep)
gene_list = np.loadtxt(args.gene_list, dtype=str)
ld_sfx = [f'ld{n}' for n in args.n_ld.split(',')]
cmd_dict = {method: [] for method in ['kggsee', 'ldak_sumher', 'ldsc', 'lder']}

if not args.skip_assoc_ld:
    ld_sfx = ['gwa'] + ld_sfx
if not args.skip_gcta:
    cmd_dict['gcta'] = []

targets = list()
targets_gcta = list()

for j in range(n_rep):
    logging.info(f'Processing commands for rep{j} ...')

    for gene in gene_list:
        dir_gene = f'{args.out_dir}/{gene}'
        dir_rep = f'{args.out_dir}/{gene}/rep{j}'

        if not args.skip_mkdir:
            os.system(f'mkdir -p {dir_rep}/plinkLD/ldetect-data')
            os.system(f'ln -fs ../../../region.lder {dir_rep}/plinkLD/ldetect-data/fourier_ls-all.bed')

        if not args.skip_gcta:
            targets_gcta.append(f' {gene}rep{j}')
            # Calculate a GRM for GATK
            cmd_dict['gcta'].append(
                f'{gene}rep{j}:\n'
                f'\t{gcta_exe} --bfile {dir_rep}/plink_gwa --make-grm --out {dir_rep}/gcta >/dev/null 2>&1\n'
            )

        for k in ld_sfx:
            targets.append(f' {gene}rep{j}{k}')

            # Make a VCF file for KGGSEE
            cmd_dict['kggsee'].append(
                f'{gene}rep{j}{k}:\n'
                f'\t{plink_exe} --recode vcf-fid bgz --real-ref-alleles '
                f'--out {dir_rep}/plink_{k} --bfile {dir_rep}/plink_{k} >/dev/null 2>&1\n'
                f'\trm {dir_rep}/plink_{k}.log {dir_rep}/plink_{k}.nosex\n'
            )

            # Calculate LD scores for LDAK-SumHer
            cmd_dict['ldak_sumher'].append(
                f'{gene}rep{j}{k}:\n'
                f'\t{ldak_exe} --thin {dir_rep}/ldak_sumher_{k} --bfile {dir_rep}/plink_{k} '
                f'--window-prune .98 --window-kb 100 >/dev/null 2>&1\n'
                "\tawk '{print $$1, 1}' "
                f'{dir_rep}/ldak_sumher_{k}.in >{dir_rep}/ldak_sumher_{k}.weights.thin\n'
                f'\t{ldak_exe} --calc-tagging {dir_rep}/ldak_sumher_{k}.ldak_thin --bfile {dir_rep}/plink_{k} '
                f'--weights {dir_rep}/ldak_sumher_{k}.weights.thin --power -0.25 --window-kb 100 >/dev/null 2>&1\n'
            )

            # Calculate LD scores for LDSC
            cmd_dict['ldsc'].append(
                f'{gene}rep{j}{k}:\n'
                f'\t{ldsc_exe} --yes-really --bfile {dir_rep}/plink_{k} --l2 --ld-wind-kb 100 '
                f'--out {dir_rep}/ldsc_{k} >/dev/null 2>&1\n'
            )

            # Calculate LD scores for LDER
            cmd_dict['lder'].append(
                f'{gene}rep{j}{k}:\n'
                f'\t{lder_plinkLD} --bfile {dir_rep}/plink_{k} --block {dir_gene}/region.lder '
                f'--snplist {dir_gene}/lder.snp --output {dir_rep}/lder_{k}.LD >/dev/null 2>&1\n'
            )

# Write a makefile for each method that calls all shell scripts for every genes and repetitions.
for method, cmd in cmd_dict.items():
    makefile = f'{args.out_dir}/makefile.step1.{method}'

    logging.info(f'Writing {makefile} ...')
    if method == 'gcta':
        first_line = f'{makefile}.done: ' + \
                     ' '.join(targets_gcta) + \
                     f'\n\ttouch {makefile}.done\n'
    else:
        first_line = f'{makefile}.done: ' + \
                     ' '.join(targets) + \
                     f'\n\ttouch {makefile}.done\n'

    with open(f'{makefile}', 'w') as o:
        o.write(first_line)
        o.writelines(cmd)

logging.info(f'Done 3.makefile_step1.')
