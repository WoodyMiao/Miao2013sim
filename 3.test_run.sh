#!/usr/bin/env bash

# ----------------------
# Copyright 2023 PMG Lab
# Author: Lin Miao
# Licence: MIT
# Version: 20230825
# ----------------------

set -e

# Extract the first five genes in each quantile range of nSNPs in genes for the test
rm -f ref_haplotypes/chr1.one4th.lst.head5concat
for q in q0-q1 q1-q2 q2-q3 q3-q4; do
  head -n 5 ref_haplotypes/chr1.one4th_$q.lst >ref_haplotypes/chr1.one4th_$q.lst.head5
  cat ref_haplotypes/chr1.one4th_$q.lst.head5 >>ref_haplotypes/chr1.one4th.lst.head5concat
done

n_rep=5                  # number of repetitions
n_gwa=5e3                # number of the sample size of association tests
n_ld=2e3,5e2             # comma separated LD reference sample sizes
h2g=1e-2                 # number of the target heritability
pqtl_vals=1,0.1,0.01     # comma separated proportions of QTLs
neg_alpha=1.0            # number of the negative alpha value
maf_min=0.01             # number of the minor allele frequency threshold
prevalence_vals=0.1,0.01 # comma separated prevalence values
n_pop_beta=1e6           # number of the samples size used to rescale the effect sizes
assoc=logit              # the model of which the results will be used to estimate heritability
chrom=1                  # chromosome ID
jobs=25                  # number of jobs (1 thread / job by default) to run in parallel when estimating heritability.

mkdir -p 3.test_out

# 1.hapsim_and_realize
python_scripts/1.hapsim_and_realize_dichotomous.py --nt 5 \
  --n-gwa $n_gwa --n-ld $n_ld --n-rep $n_rep --maf-min $maf_min \
  --h2g $h2g --pqtl-vals $pqtl_vals --neg-alpha $neg_alpha \
  --prevalence-vals $prevalence_vals --n-pop-beta $n_pop_beta \
  --region-file ref_haplotypes/chr1.codegene_10kbflk.snp_counts.tsv \
  --vcf-ref ref_haplotypes/chr1.codegene_10kbflk.vcf.gz \
  --gene-list ref_haplotypes/chr1.one4th.lst.head5concat \
  --out-dir 3.test_out

# 3.makefile_step1
python_scripts/3.makefile_step1.py \
  --n-ld $n_ld --n-rep $n_rep --skip-gcta --skip-assoc-ld \
  --out-dir 3.test_out \
  --gene-list ref_haplotypes/chr1.one4th.lst.head5concat

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

for method in kggsee ldak_sumher ldsc lder; do
  echo "Run $method ..."
  make --jobs $jobs --file 3.test_out/makefile.step1.$method
done

# 4.makefile_step2
python_scripts/4.makefile_step2_dichotomous.py --jobs $jobs --chrom $chrom \
  --n-gwa $n_gwa --n-ld $n_ld --n-rep $n_rep --h2g $h2g --maf-min $maf_min \
  --neg-alpha $neg_alpha --pqtl-vals $pqtl_vals --prevalence-vals $prevalence_vals \
  --assoc-test $assoc --kggsee --ldsc --lder --ldak-gbat --ldak-sumher --hess \
  --gene-list ref_haplotypes/chr1.one4th.lst.head5concat \
  --out-dir $(pwd)/3.test_out # Here must use an absolute path, otherwise LDER will fail.

for method in kggsee hess ldak_gbat ldak_sumher ldsc lder; do
  echo "Run $method ..."
  sh 3.test_out/makefile.step2.$assoc.$method.sh
done

# 5.harvest_outputs
for q in q0-q1 q1-q2 q2-q3 q3-q4; do
  python_scripts/5.harvest_outputs.py --nt 5 \
    --prevalence-vals $prevalence_vals --n-ld $n_ld --n-rep $n_rep \
    --h2g-vals $h2g --pqtl-vals $pqtl_vals --neg-alpha-vals $neg_alpha \
    --kggsee --hess --ldak-gbat --ldak-sumher --ldsc --lder \
    --binary-assoc-test $assoc --skip-assoc-ld \
    --gene-list ref_haplotypes/chr1.one4th_$q.lst.head5 \
    --out-suffix one4th_$q.tsv --out-dir 3.test_out
done

# 6.plot_results
python_scripts/plot_scripts/3.repeat_dichotomous.py --prefix 3.test_out \
  --h2g $h2g --neg-alpha $neg_alpha --pqtl-vals $pqtl_vals \
  --prevalence-vals $prevalence_vals --binary-assoc-test $assoc
echo "All Done."
