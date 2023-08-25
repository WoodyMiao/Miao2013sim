## Materials capable of repeating the simulation analyses in Miao *et al*. 2013

### Repeating the simulation analyses

We created three shell scripts describing how to repeat the simulation analyses (**Figure 1, 2, 3, S1, S2, S3 and S4**)
in our published paper:

Lin Miao, Lin Jiang, Bin Tang, Pak Chung Sham, Miaoxin Li.
Dissecting the high-resolution genetic architecture of complex phenotypes by accurately
estimating gene-based conditional heritability.
The American Journal of Human Genetics (2023). https://doi.org/10.1016/j.ajhg.2023.08.006

The following table is an overview of the simulation in each shell script.

| Shell script in this repository                                                                                                                                                              | Corresponding figure in our paper                                                                                                                                                                                                                                                                                                           |
|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `1.repeat_quantitative.sh`<br/><ul><li>Quantitative Phenotype</li><li>$p_{c} \in \\{1, 0.1, 0.01\\}, \alpha \in \\{-1, -0.25\\}$</li><li>EHE, HESS, GBAT, SumHer, LDSC, LDER</li></ul>       | **Figure1**, **3** and **S4** show $\text{MRB of }h^2$,  $\text{MRB of SE}(h^2)$ and $\text{SD}(h^2)$ respectively calculated from simulating the six possible genetic architectures when $h^2_g=0.1\\%$.<br/>**Figure S2** shows $\text{MRB of}\ h^2$ calculated from simulating the six possible genetic architectures when $h^2_g=1\\%$. |
| `2.repeat_quantitative_gcta.sh`<br/><ul><li>Quantitative Phenotype</li><li>$p_{c} \in \\{1, 0.1, 0.01\\}, \alpha=-1$</li><li>GCTA, EHE, HESS, GBAT, SumHer, LDSC, LDER</li></ul>             | **Figure 2** shows $\text{Mean}(h^2)$, $\text{MRB of }h^2$, $\text{SD}(h^2)$ and $\text{MRB of SE}(h^2)$ calculated from simulating the three possible genetic architectures when $h^2_g=0.1\\%$.                                                                                                                                           |
| `3.repeat_dichotomous.sh`<br/><ul><li>Dichotomous Phenotype</li><li>$K \in \\{0.1,0.01\\}, p_{c} \in \\{1, 0.1, 0.01\\}, \alpha=-0.25$</li><li>EHE, HESS, GBAT, SumHer, LDSC, LDER</li></ul> | **Figure S3** shows $\text{MRB of }h^2$ calculated from simulating the six possible genetic architectures when $h^2_g=0.1\\%$.                                                                                                                                                                                                              |

Note: $p_{c}$ denotes proportion of causal SNPs and $K$ denotes proportion of affected individuals in the
population.

To repeat the analysis using only HapMap3 SNPs (**Figure S1**), filter `ref_haplotypes/chr1.codegene_10kbflk.vcf.gz`
to containing only HapMap3 SNPs and run `1.repeat_quantitative.sh`.
In addition, the script `ref_haplotypes.sh` shows how we made the files in the folder `ref_haplotypes`,
i.e., where we downloaded and how we filtered the chr1 VCF file of the 1000 Genomes Project.

### Setting up a runtime environment

We also created three test scripts (`*.test_run.sh`) for testing a runtime environment. These scripts have the
same structure as the corresponding repeat scripts (`*.repeat_*.sh`), but only run 5 genes with 5 repetitions.
The file `environment.yml` records the conda environment where we run these analyses. The project URLs of the methods
compared in this study are shown in the following table with the version that we run in this study.

| Method[^1] | Project URL                                                                                            | The version we used                                               |
|------------|--------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------|
| KGGSEE     | http://pmglab.top/kggsee/#/                                                                            | v1.1                                                              |
| PLINK      | https://www.cog-genomics.org/plink/                                                                    | v1.90b6.12 (28 Oct 2019)                                          |
| GCTA       | https://yanglab.westlake.edu.cn/software/gcta/#Overview                                                | v1.94.1 (Built at Nov15 2022)                                     |
| HESS[^2]   | https://huwenboshi.github.io/hess/                                                                     | v0.5 (9/October/2017)                                             |
| LDAK       | https://dougspeed.com/                                                                                 | v5.2                                                              |
| LDSC       | https://github.com/bulik/ldsc                                                                          | v1.0.1                                                            |
| LDER       | https://github.com/shuangsong0110/LDER                                                                 | v0.1.0                                                            |
| plinkLD    | https://cloud.tsinghua.edu.cn/f/3f96074ee7ee436895ac/?dl=1 (This URL is provided by the LDER project.) | Created on 2016-6-7. Add shrinkage estimator feature on 2020-8-3. |
| HapSim     | https://cran.r-project.org/web/packages/hapsim/index.html                                              | v0.31                                                             |

[^1]: Paths to the executables of the methods need to be specified in `3.makefile_step1.py`, `4.makefile_step2.py`,
and `4.makefile_step2_dichotomous.py`.

[^2]: The loop of 22 autosomes in step two of HESS code needs to be removed to go through with the result of only one
gene in step one.

### Description of the python scripts

The shell scripts call python scripts in the folder `python_scripts`:

* First, `1.hapsim_a_population.py` and `2.realize_beta_and_gwa.py`
  generate genotypes and quantitative phenotypes, and perform association tests; while
  `1.hapsim_and_realize_dichotomous.py` does the same procedure and is for dichotomous phenotypes.

* Then, `3.makefile_step1.py`, `4.makefile_step2.py`,
  and `4.makefile_step2_dichotomous.py` write GNU makefiles controlling heritability estimation
  by the compared methods. Paths to the executables of the compared methods are specified in these python scripts.

* After the heritability estimation has been done, `5.harvest_outputs.py` collects results and
  calculates $\text{Mean}(h^2)$, $\text{MRB of }h^2$, $\text{SD}(h^2)$ and $\text{MRB of SE}(h^2)$ for each genetic
  architecture. Finally, `plot_scripts/*.py` plot figures.
