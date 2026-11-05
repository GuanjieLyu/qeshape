# qeshape

Quantile–expectile calibration curves, simultaneous confidence bands, and
bootstrap tests of distributional symmetry in R.

`qeshape` provides the core methods from Guanjie Lyu's manuscript,
*Quantile–expectile calibration curves for distributional shape: estimation,
inference, and symmetry testing* (2026). It uses base R and the standard R
packages, with no additional runtime dependencies.

## Install

Requires R 4.1.0 or later. From the package source directory (the directory
containing `DESCRIPTION`), install with:

```sh
R CMD INSTALL .
```

For GitHub, place the contents of `R_package` at the repository root,
including the hidden `.github` directory. The package name is `qeshape`.

## Quick start

```r
library(qeshape)

set.seed(2026)
x <- rlnorm(150)

curve <- qe_curve(x)
plot(curve)
head(as.data.frame(curve))
predict(curve, probs = c(0.1, 0.5, 0.9))

asymmetry <- qe_asymmetry(x)
plot(asymmetry)

# Small B keeps this example quick. Increase B for an analysis.
band <- qe_band(x, trim = 0.05, conf.level = 0.95, B = 199, seed = 7)
plot(band)
predict(band, probs = c(0.1, 0.5, 0.9))

qe_symmetry_test(x, method = "cvm", B = 199, seed = 7)
qe_symmetry_test(x, method = "ks", B = 199, seed = 7)
```

## Reference manual

The [PDF reference manual](qeshape-manual.pdf) contains the package
overview, statistical assumptions, complete function reference, worked
examples, and index. It is generated from the same help files installed
with the package. To rebuild it from the source directory, with R and
LaTeX available, run:

```sh
Rscript tools/build-manual.R
```

## Functions

| Function | Result |
|:---|:---|
| `qe_curve()` | Empirical calibration curve and sample quantiles |
| `qe_asymmetry()` | Reflected calibration sum minus one |
| `qe_band()` | Simultaneous bootstrap band over a fixed trimmed interval |
| `qe_symmetry_test()` | Centered-bootstrap CvM or KS symmetry test |

Curve, asymmetry, and band objects support `print()`, `plot()`,
`as.data.frame()`, and `predict()`. A test returns an `htest` object.
Predictions evaluate the empirical curve from the retained observations;
they do not interpolate a plotting grid. Band predictions stay within the
interval used to construct the band.

For details, use `help(package = "qeshape")` and the installed methods guide:

```r
file.show(system.file("doc", "methods.md", package = "qeshape"))
```

## What is estimated?

At probability `p`, let `q` be the left empirical quantile, corresponding
to `quantile(x, p, type = 1)`. The calibration estimate is

```text
tau(p) = sum(pmax(q - x, 0)) / sum(abs(x - q)).
```

This describes the expectile level associated with the sample quantile.
It is invariant to translation and positive rescaling. At a sample
endpoint, its value can be zero or one. For a constant sample the package
uses `tau(p) = 0.5`; inference requires a nonconstant sample.

The asymmetry curve is `D(p) = tau(p) + tau(1 - p) - 1`, for
`0 < p <= 0.5`. The default CvM test integrates its squared departure
from zero. The KS option uses its largest absolute departure. The paper
uses CvM as the primary test and reports more conservative finite-sample
behavior for KS.

## Computational choices

- The ordinary empirical bootstrap is implemented through multinomial
  observation counts, avoiding repeated sorting of resampled data.
- The band and KS test evaluate their suprema on exact finite grids.
  The CvM test integrates the empirical step function exactly with unit
  weight. There is no numerical integration grid to tune.
- Bootstrap symmetry statistics are centered at the observed asymmetry
  curve. No symmetry center, density, bandwidth, or variance estimate is
  fitted.
- Band critical values use the left empirical quantile (`type = 1`),
  and band endpoints are clipped to `[0, 1]`.
- Test p-values are `(1 + sum(bootstrap >= observed)) / (B + 1)`.
  This prevents a reported zero caused by a finite bootstrap run. It
  does **not** make the bootstrap test exact in finite samples.
  A rounding tolerance on the asymmetry scale counts numerical ties as
  exceedances; `p.value.threshold` records the cutoff actually used.

The p-value correction is a deliberate difference from the uncorrected
bootstrap proportions in the paper's simulation scripts. Consequently,
this package is an implementation of the methods, not a replacement for
the archived replication programs.

## Using inference responsibly

The theory concerns independent, identically distributed observations,
a finite second moment, and a continuously differentiable distribution
with positive density near the trimmed quantile range. The band and test
results also require the limiting statistic's distribution to be
continuous and to strictly cross the relevant quantile level. The trim
is fixed and chosen independently of the data. Increasing `B` reduces
Monte Carlo error; the asymptotic results concern `B` increasing with
sample size.

The point estimators accept ties, rounded data, and constant samples.
Accepting such input does not establish the regularity assumptions for
inference. A test with fixed trimming detects departures visible on its
trimmed interval; a large p-value does not establish symmetry throughout
the distribution. Separate bands for two samples are descriptive and do
not constitute a calibrated two-sample test.

Missing observations cause an error unless `na.rm = TRUE`; that option
removes `NA` and `NaN`, while infinite values always cause an error. With
`seed = NULL`, bootstrap functions use and advance the caller's random
number stream. A supplied seed is scoped to the call and restores
`.Random.seed` and the RNG kind afterward. R does not expose the cached
normal draw used by its `"Box-Muller"` generator, so that cache cannot be
restored. The default `"Inversion"` normal generator avoids this issue;
use `seed = NULL` when continuation of an existing Box-Muller cache matters.

## Development and citation

From the package source directory:

```sh
R CMD build .
R CMD check --as-cran qeshape_0.1.1.tar.gz
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for contribution guidance, and use
`citation("qeshape")` for the manuscript citation. This package is released
under the MIT license; see [LICENSE.md](LICENSE.md).
