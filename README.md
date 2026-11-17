# ICSNP: Tools for Multivariate Nonparametrics

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/ICSNP)](https://CRAN.R-project.org/package=ICSNP)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/ICSNP)](https://CRAN.R-project.org/package=ICSNP)
[![R-CMD-check](https://github.com/klauschN/ICSNP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/klauschN/ICSNP/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Classical multivariate methods such as Hotelling's T² test are based on the sample mean vector and the sample covariance matrix, which makes them optimal under multivariate normality but sensitive to outliers and heavy-tailed distributions. The R package **ICSNP** provides nonparametric and robust alternatives based on marginal signs and ranks and on spatial signs and ranks (see, e.g., Nordhausen and Oja, 2018), together with Hotelling's T² test as a parametric reference. The tests can also be applied to invariant coordinates obtained with the package [ICS](https://CRAN.R-project.org/package=ICS), which yields affine invariant tests, for example in the independent component model.

The package has been available on CRAN since 2007, and several other CRAN packages, such as [SpatialNP](https://CRAN.R-project.org/package=SpatialNP), [MNM](https://CRAN.R-project.org/package=MNM) and [ICtest](https://CRAN.R-project.org/package=ICtest), build on it.

## Installation

The released version can be installed from CRAN:

```r
install.packages("ICSNP")
```

The development version can be installed from this repository, which requires a C compiler (on Windows, Rtools):

```r
# install.packages("remotes")
remotes::install_github("klauschN/ICSNP")
```

## Main functionality

| Group | Functions |
|---|---|
| Location tests | `HotellingsT2`, `rank.ctest` (one-, two- and c-sample tests based on marginal signs, ranks and normal scores), `rank.ictest` (tests in the symmetric independent component model), `HP.loc.test` (optimal tests of Hallin and Paindaveine) |
| Tests of independence | `ind.ctest`, `ind.ictest` |
| Location estimates | `spatial.median`, `HR.Mest` (Hettmansperger-Randles estimate), `hl.loc`, `vdw.loc` |
| Shape and scatter estimates | `tyler.shape`, `duembgen.shape`, `symm.huber`, `HP1.shape` and the weighted versions `duembgen.shape.wt` and `symm.huber.wt` |
| Utilities and data | `spatial.sign`, `pair.diff`, `pair.sum`, `pair.prod`, the data set `pulmonary` |

All tests return objects of class `htest`, as for example `t.test` and `wilcox.test` in base R, and the location tests have methods for data matrices, for formulas of the form `Y ~ g` and for `ics` objects. The computationally demanding parts, that is, the loops over all pairs of observations needed for D&uuml;mbgen's shape matrix and the symmetrized Huber scatter matrix, are implemented in C.

## Example

We test whether the location of the changes in three measurements of the pulmonary function of 12 workers exposed to cotton dust is the origin, using the marginal signed-rank test, Hotelling's T² test and the affine invariant signed-rank test of Hallin and Paindaveine:

```r
library("ICSNP")
data("pulmonary")
rank.ctest(pulmonary, scores = "rank")
HotellingsT2(pulmonary)
HP.loc.test(pulmonary, score = "rank")
```

Robust estimates of location and shape are obtained, for example, as

```r
X <- as.matrix(pulmonary)
spatial.median(X)
HR.Mest(X)
duembgen.shape(X)
```

More examples are given in the help pages of the package and in the paper listed below. The R script [paper/ICSNP_paper_code.R](paper/ICSNP_paper_code.R) reproduces all examples of the paper, including its figure.

## Citation

If you use **ICSNP** in your work, please cite the package (see `citation("ICSNP")` in R). A paper describing the package has been submitted:

Nordhausen, K., Oja, H. and Tyler, D. E. ICSNP: An R package for multivariate nonparametric location tests and robust estimates of location and shape. Submitted.

## References

- Nordhausen, K. and Oja, H. (2018). Robust nonparametric inference. *Annual Review of Statistics and Its Application*, 5, 473-500. <https://doi.org/10.1146/annurev-statistics-031017-100247>
- Oja, H. (2010). *Multivariate Nonparametric Methods with R: An Approach Based on Spatial Signs and Ranks*. Springer, New York. <https://doi.org/10.1007/978-1-4419-0468-3>

## Bug reports and contributions

Bug reports, questions and suggestions are welcome as [GitHub issues](https://github.com/klauschN/ICSNP/issues). Pull requests should be accompanied by tests in `tests/testthat`, and the package should pass `R CMD check --as-cran` without errors, warnings or notes.

## License

**ICSNP** is free software released under the GNU General Public License, version 2 or (at your option) any later version, see the file [LICENSE.txt](LICENSE.txt).
