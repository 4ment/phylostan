# Evaluating probabilistic programming and fast variational Bayesian inference in phylogenetics

This repository contains the pipeline and data sets supporting the results of the following article:

Fourment Mathieu and Darling E. Aaron. Evaluating probabilistic programming and fast variational Bayesian inference in phylogenetics. [doi:10.1101/702944](https://doi.org/10.1101/702944).

# Requirements

To reproduce the same results the exact version of the following programs should be used.

## [phylostan](https://github.com/4ment/phylostan)
## [physher vb-v1.0](https://github.com/4ment/physher/releases/tag/vb-v1.0)
## [beast v1.10.4](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.4)
## [beast2 v2.5.2](https://github.com/CompEvol/beast2/releases/tag/v2.5.2)
The `bwest` package must be installed in order to use the Weibull distribution. `phylostan` uses the Weibull distribution to model rate heterogeneity across sites instead of the gamma distribution since automatic differentiation of the gamma distribution is not available in Stan.

The BEAST2 website has a blog [post](https://www.beast2.org/managing-packages/) explaining how to install packages. This package is not part of the official package repository so you need to register this package first.

 - Open BEAUTi, click on `Manage Packages` menu the in the `File` menu. This will open a window showing the package names.
 - Click on the `Package repositories` button, then click on the `Add URL` button.
 - Write the following url `https://4ment.github.io/assets/beast2/packages.xml` in the new window. Click `Done` to close the window.
 - Locate and select the package name `bwest` in the table and click the `Install/Upgrade` button.
 - Close the window et voila the package should be installed.

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

# Generating the script files only

``` shell
cd phylostan/examples
scons dry=1
```

To run `physher` on one the DS1 file:

``` shell
cd examples/DS1/0
physher tree0.json
```
