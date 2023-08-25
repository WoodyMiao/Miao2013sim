#!/usr/bin/env bash

# ----------------------
# Copyright 2023 PMG Lab
# Author: Lin Miao
# Licence: MIT
# Version: 20230825
# ----------------------

set -e

n_rep=50             # number of repetitions
n_gwa=2e4            # number of the sample size of association tests
n_ld=2e3,5e2         # comma separated LD reference sample sizes
h2g=1e-3             # comma separated target heritability values
pqtl_vals=1,0.1,0.01 # comma separated proportions of QTLs
neg_alpha=1.0        # comma separated negative alpha values
maf_min=0.01         # number of the minor allele frequency threshold
chrom=1              # chromosome ID
jobs=25              # number of jobs (1 thread / job by default) to run in parallel when estimating heritability.

mkdir -p 2.repeat_out

# 1.hapsim_a_population
python_scripts/1.hapsim_a_population.py --nt 5 \
  --n-gwa $n_gwa --n-ld $n_ld --n-rep $n_rep --maf-min $maf_min \
  --region-file ref_haplotypes/chr1.codegene_10kbflk.snp_counts.tsv \
  --vcf-ref ref_haplotypes/chr1.codegene_10kbflk.vcf.gz \
  --gene-list ref_haplotypes/chr1.one32th.lst \
  --out-dir 2.repeat_out

# 2.realize_beta_and_gwa
python_scripts/2.realize_beta_and_gwa.py --nt 5 \
  --h2g-vals $h2g --pqtl-vals $pqtl_vals \
  --neg-alpha-vals $neg_alpha --n-gwa $n_gwa --n-rep $n_rep \
  --gene-list ref_haplotypes/chr1.one32th.lst \
  --out-dir 2.repeat_out

# 3.makefile_step1
python_scripts/3.makefile_step1.py \
  --n-ld $n_ld --n-rep $n_rep \
  --gene-list ref_haplotypes/chr1.one32th.lst \
  --out-dir 2.repeat_out

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

for method in kggsee ldak_sumher ldsc lder; do
  echo "Run $method ..."
  make --jobs $jobs --file 2.repeat_out/makefile.step1.$method
done
echo "Run gcta ..."
make -j 5 -f 2.repeat_out/makefile.step1.gcta

# 4.makefile_step2
python_scripts/4.makefile_step2.py --jobs $jobs --chrom $chrom \
  --n-gwa $n_gwa --n-ld $n_ld --n-rep $n_rep --h2g-vals $h2g \
  --pqtl-vals $pqtl_vals --neg-alpha-vals $neg_alpha --maf-min $maf_min \
  --gcta --kggsee --hess --ldak-gbat --ldak-sumher --ldsc --lder \
  --gene-list ref_haplotypes/chr1.one32th.lst \
  --out-dir $(pwd)/2.repeat_out # Here must use an absolute path, otherwise LDER will fail

for method in gcta kggsee hess ldak_gbat ldak_sumher ldsc lder; do
  echo "Run $method ..."
  sh 2.repeat_out/makefile.step2.$method.sh
done

# 5.harvest_outputs
python_scripts/5.harvest_outputs.py --nt 5 \
  --n-ld $n_ld --n-rep $n_rep --h2g-vals $h2g \
  --pqtl-vals $pqtl_vals --neg-alpha-vals $neg_alpha \
  --gcta --kggsee --ldsc --lder --ldak-gbat --ldak-sumher --hess \
  --gene-list ref_haplotypes/chr1.one32th.lst \
  --out-suffix one32th.tsv --out-dir 2.repeat_out

# 6.plot_results
python_scripts/plot_scripts/2.repeat_quantitative_gcta.py --prefix 2.repeat_out \
  --h2g $h2g --neg-alpha $neg_alpha --pqtl-vals $pqtl_vals
echo "All Done."
