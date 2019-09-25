# Evaluating probabilistic programming and fast variational Bayesian inference in phylogenetics

This repository contains the pipeline and data sets supporting the results of the following article:

Fourment Mathieu and Darling E. Aaron. Evaluating probabilistic programming and fast variational Bayesian inference in phylogenetics. [doi:10.1101/702944](https://doi.org/10.1101/702944).

# Requirements

To reproduce the same results the exact version of the following programs should be used.

## [phylostan](https://github.com/4ment/phylostan)
## [physher vb-v1.0](https://github.com/4ment/physher/releases/tag/vb-v1.0)
## [beast v1.10.4](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.4)
## [beast2 v2.5.2](https://github.com/CompEvol/beast2/releases/tag/v2.5.2)

The executables `phylostan`, `physher`, `beast`, `beast2` should be in the `PATH`.

# Running the simulations

``` shell
cd phylostan/examples
scons .
```

or using multithreading (5 threads in this case):

``` shell
cd phylostan/examples
scons -j 5 .
```

The simulations will take several days to complete.