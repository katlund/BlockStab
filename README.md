# BlockStab

The main purpose of this package is to study, verify, and conjecture the stability properties of different versions of Block Gram-Schmidt (BGS) via a skeleton-muscle paradigm.

## Installation

Follow the download options from the Git repository main page.  Then navigate to the repo (`BlockStab`) and run `install_blockstab.m` in MATLAB.  Note that this script only temporarily saves the paths; they will be cleared at the next start-up.  To permanently save `BlockStab` routines to the startup path, run `savepath` after `install_blockstab.m`, which may overwrite paths to other functions with the same names.

## What is new in this version

* [x] [Mixed precision implementations](#mixed-precision)
* [x] Additional low-sync versions of BCGSI+, which help demonstrate finer-grained stability properties
* [x] A Cholesky switch, allowing for users to specify which Cholesky subroutine to use
* [ ] `RunKappaPlot`: a unified, streamlined test engine that avoids redundant runs of skeleton-muscle combinations, simplifies syntax via an options struct, improves display of figure outputs, allows for toggling how and whether figures are saved, and allows for automatic TeX report generation.

To reproduce results from [Carson, et al. 2022](https://doi.org/10.1016/j.laa.2021.12.017), please use release [v1.2022](https://github.com/katlund/BlockStab/releases/tag/v1.2022).

### Mixed precision

Mixed precision routines (i.e., those ending with `_mp`) require one of the additional toolboxes:

* [Advanpix Multiprecision Computing Toolbox](https://www.advanpix.com/), which requires a paid license.  The `mp` subroutine is used.
* [Symbolic Math Toolbox](https://www.mathworks.com/products/symbolic.html), which may also require a paid license.  The subroutine [`vpa`](https://mathworks.com/help/symbolic/vpa.html) is used.

The subroutine `mp_switch` manages which toolbox is called and at what precision via the `param` struct (see below).

## Gram-Schmidt Routines

Each skeleton and muscle can be run individually or via the drivers

```matlab
BGS(XX, s, skel, musc, param)
```

and

```matlab
IntraOrtho(X, musc, param)
```

for skeletons (`BGS`) and muscles (`IntraOrtho`), respectively.

The variable `XX` denotes a block-partitioned matrix with `m` rows, `p` block vectors, and `s` columns per block vector, i.e., $m = ps$.  `X` denotes a single block vector (or a tall-and-skinny matrix) with `m` rows and `s` columns, $s \leq m$. As for the other parameters:

* `skel` - char specifying BGS skeleton
* `musc` - char specifying intra-orthogonalization muscle
* `param`: a struct with the following optional fields:
  * `.chol`: char specifying what type of Cholesky subroutine to call for skeletons that hard-code Cholesky via a block Pythagorean trick (e.g., `bcgs_pip`, `bcgs_pio`, `bcgs_iro_ls`, `bmgs_cwy`, `bmgs_icwy`, and their reorthogonalized and multi-precision versions); default for non-MP versions is `'chol_nan'`, and `'chol_free'` for MP
  * `.mp_package`: char specifying either `'advanpix'` or `'symbolic toolbox'` as the mixed precision package; default: `'advanpix'`
  * `.mp_digits`: int specifiying number of precision digits, e.g., 34 for quadruple precision (in Advanpix) or 32 for quadruple precision in Symbolic Math Toolbox; default: 34
  * `.rpltol`: scalar argument for `cgs_sror` that determines the replacement tolerance; default: 1
  * `.verbose`: boolean for whether to print intermediate loss of orthogonality (LOO) or relative residual (RelRes) per iteration; default: 0

For a list of all currently implemented skeletons and muscles, see `BGS` and`IntraOrtho`. Try, for example:

```matlab
rng(4);
mgs(randn(100,10), true);
```

or equivalently

```matlab
rng(4);
param.verbose = true;
IntraOrtho(randn(100,10), 'MGS', param);
```

Example output:

```matlab
         LOO      |    RelRes
-----------------------------------
  1:  0.0000e+00  |  4.2294e-17
  2:  4.4452e-16  |  4.3673e-17
  3:  4.4469e-16  |  8.0352e-17
  4:  4.4519e-16  |  1.0026e-16
  5:  4.4830e-16  |  1.0289e-16
  6:  5.5750e-16  |  9.8092e-17
  7:  5.5811e-16  |  1.1045e-16
  8:  5.5868e-16  |  1.3219e-16
  9:  5.6484e-16  |  1.3937e-16
 10:  5.6465e-16  |  1.5539e-16
 ```

For block methods, try

```matlab
rng(4);
param.verbose = true;
bmgs(randn(100,20), 2, 'HouseQR', param);
```

or equivalently

```matlab
rng(4);
param.verbose = true;
BGS(randn(100,20), 2, 'BMGS', 'HouseQR', param);
```

Example output:

```matlab
         LOO      |    RelRes
-----------------------------------
  1:  4.4414e-16  |  1.1274e-16
  2:  4.7420e-16  |  1.8449e-16
  3:  7.1500e-16  |  1.8485e-16
  4:  7.3794e-16  |  1.8839e-16
  5:  7.4766e-16  |  1.8543e-16
  6:  7.4959e-16  |  1.9422e-16
  7:  7.6506e-16  |  2.1185e-16
  8:  7.8917e-16  |  2.1527e-16
  9:  7.7993e-16  |  2.2404e-16
 10:  7.9600e-16  |  2.3448e-16
 ```

## Test Routines

There are several test files.  `MakeHeatmap` generates heatmaps comparing loss of orthogonality and residual across many skeleton-muscle combinations for the same test matrix.  `RunKappaPlot` and plots loss of orthogonality and residual trends against matrices with a range of condition numbers. As the Greek letter $\kappa$ is used to denote the 2-norm condition number of a matrix, we refer to these plots as "kappa plots."  See the header for each for full descriptions of their functionalities.  To explore some interesting examples, try the following, and note that `[]` (empty) arguments call default options, which can be quite handy:

* `MakeHeatmap([100 10 2], 'stewart',  {'BCGS', 'BCGS_IRO', 'BCGS_SROR'}, {'CGS', 'HouseQR'}, 1, 1)`
* `KappaPlot([100 10], [], {'MGS', 'MGS_SVL', 'MGS_LTS', 'MGS_CWY', 'MGS_ICWY'})`
* `BlockKappaPlot([100 20 2], [], {'BCGS', 'BCGS_IRO', 'BCGS_IRO_LS'}, {'CGS', 'MGS', 'HouseQR'})`
* `GluedKappaPlot([], [], {'CGS', 'CGS_P', 'MGS', 'CGS_IRO', 'CGS_IRO_LS'})`
* `GluedBlockKappaPlot([], [], {'BCGS', 'BCGS_PIP', 'BCGS_PIO'}, 'HouseQR')`
* `GluedBlockKappaPlotVaryS([], [], [2 5], {'BCGS', 'BMGS'})`
* `MonomialBlockKappaPlot([1000 100],[1 2 5 10],{'BCGS', 'BMGS', 'BCGS_IRO'}, {'CholQR', 'MGS', 'HouseQR'})`

## Documentation

Each file contains a descriptive header.  See especially the following core files:

* `BGS.m` - switches between skeletons
* `IntraOrtho.m` - switches between muscles
* `InnerProd.m` - switches between inner products
* `MakeHeatmap.m` - generates heatmaps for skeleton-muscle combinations; verbose = `true` prints tables to screen
* `RunKappaPlot.m` - generates kappa plots for muscles and skeleton-muscle combinations

## How we cite things

Several papers are foundational for our subroutines.  We provide the full citations here and use abbreviated ones (given as [Author YYYY]) throughout the documentation.

* [Barlow 2019](https://doi.org/10.1137/18M1197400): Barlow, J. Block modified Gram-Schmidt algorithms and their analysis. SIAM Journal on Matrix Analysis and Applications. Vol. 40, 4, pp. 1257--1290, 2019.
* [Barlow & Smoktunowicz 2013](https://doi.org/10.1007/s00211-012-0496-2): Barlow, J. & Smoktunowicz, A. Reorthogonalized block classical Gram-Schmidt. Numerische Mathematik. Vol 123, pp 395--423, 2013.
* [Bielich, et al. 2022](https://doi.org/10.1016/j.parco.2022.102940): Bielich, D, Langou J., Thomas, S., Świrydowicz, K., Yamazaki, I., and Boman, E.G.  Low-synch Gram–Schmidt with delayed reorthogonalization for Krylov solvers.  Parallel Computing. Vol 112, pp 102940, 2022.
* [Carson, et al. 2021](https://doi.org/10.1137/21M1394424): Carson, E., Lund, K., and Rozložník, M.  The stability of block variants of classical Gram-Schmidt.  SIAM Journal on Matrix Analysis and Applications. Vol. 42, 3, pp 1365--1380, 2021.
* [Fukaya, et al. 2020](https://doi.org/10.1137/18M1218212): Fukaya, T., Kannan, R., Nakatsukasa, Y., Yamamoto, Y., & Yanagisawa, Y. Shifted CholeskyQR for computing the QR factorization of ill-conditioned matrices. SIAM Journal on Scientific Computing. Vol. 42, 1, pp A477--A503, 2020.
* [Fukaya, et al. 2014](https://doi.org/10.1109/ScalA.2014.11): Fukaya, T., Nakatsukasa, Y., Yanagisawa, Y., & Yamamoto, Y. CholeskyQR2: A simple and communication-avoiding algorithm for computing a tall-skinny QR factorization on a large-scale parallel system. 2014 5th Workshop on Latest Advances in Scalable Algorithms for Large-Scale Systems, pp 31--38, 2014.
* [Smoktunowicz, et al. 2006](https://doi.org/10.1007/s00211-006-0042-1): Smoktunowicz, A., Barlow, J., and Langou, J. A note on the error analysis of classical Gram-Schmidt. Numerische Mathematik. Vol. 105, 2, pp 299-313, 2006.
* [Stathopoulos & Wu 2002](https://doi.org/10.1137/S1064827500370883): Stathopoulos, A. & Wu, K. A block orthogonalization procedure with constant synchronization requirements. SIAM Journal on Scientific Computing. Vol 23, 6, pp 2165--2182, 2002.
* [Stewart 2008](https://doi.org/10.1137/070682563): Stewart, G. W. Block Gram-Schmidt orthogonalization. SIAM Journal on Scientific Computing. Vol 31, 1, pp 761--775, 2008.
* [Swirydowicz, et al. 2020](https://doi.org/10.1016/j.parco.2022.102940): Świrydowicz, K., Langou, J., Ananthan, S., Yang, U., & Thomas, S. Low synchronization Gram-Schmidt and GMRES algorithms. Technical report, 2020.

## How to cite us

[Carson, et al. 2022](https://doi.org/10.1016/j.laa.2021.12.017) Carson, E., Lund, K, Rozložník, M., and Thomas, S. Block Gram-Schmidt methods and their stability properties. Linear Algebra and its Applications. Vol 638, pp 150--195, 2022.

Please also mention which version of the software you are using by referring, e.g., to a tag or specific commit.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

If you are interested in helping develop an open-source version of this package in either Python or Julia, please contact [Kathryn Lund](mailto:kathryn.d.lund@gmail.com) directly.

## Related projects

See [LowSyncBlockArnoldi](https://gitlab.mpi-magdeburg.mpg.de/lund/low-sync-block-arnoldi) for block Arnoldi versions of these routines and a simple algorithm comparison workflow.
