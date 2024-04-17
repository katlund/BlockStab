# BlockStab

The main purpose of this package is to study, verify, and conjecture the stability properties of different versions of Block Gram-Schmidt (BGS) via a skeleton-muscle paradigm.

## Requirements

MATLAB 2020b or higher is required for saving PDFs of plots with minimal whitespace; see [this link](https://de.mathworks.com/help/matlab/creating_plots/save-figure-with-minimal-white-space.html).  For older versions of MATLAB, be sure to set `run_data.options.save_pdf = 0`, or modify `gen_plots.m` directly to save PDFs via a preferred format.

See [Multiprecision](#multiprecision) for additional package requirements.

We rely on [`linspecer`](https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap) and [`catstruct`](https://www.mathworks.com/matlabcentral/fileexchange/7842-catstruct), which are both included and do not need to be downloaded separately.

Otherwise, the main code base is likely to work with minimal modifications in GNU Octave, but we have not tested this directly.

## Installation

Follow the download options from the Git repository main page.  Then navigate to the repo (`BlockStab`) and run `install_blockstab.m` in MATLAB.  Note that this script only temporarily saves the paths; they will be cleared at the next start-up.  To permanently save `BlockStab` routines to the startup path, run `savepath` after `install_blockstab.m`, which may overwrite paths to other functions with the same names.

## What is new in this version

* [Multiprecision implementations](#multiprecision)
* Additional low-sync versions of BCGSI+, which help demonstrate finer-grained stability properties; in particular, `_a` versions that run an O(eps)-stable `IntraOrtho` on the first block vector for extra stability.
* A Cholesky switch, allowing for users to specify which Cholesky subroutine to use
* [`RunKappaPlot`](#new-test-driver): a unified, streamlined test engine that avoids redundant runs of skeleton-muscle combinations, simplifies syntax via an options struct, improves display of figure outputs, allows for toggling how and whether figures are saved (.eps, .pdf, and .fig formats allowed), and auto-generates a timestamped TeX report, which can be easily compiled and shared with collaborators.
* A new class of test matrices called `piled` matrices: they are similar to `glued`, but can more easily highlight edge-case behavior for some methods.

## Multiprecision

Multiprecision routines (i.e., those ending with `_mp`) require one of the additional toolboxes:

* [Advanpix Multiprecision Computing Toolbox](https://www.advanpix.com/), which requires a paid license.  The `mp` subroutine is used.
* [Symbolic Math Toolbox](https://www.mathworks.com/products/symbolic.html), which may also require a paid license.  The subroutine [`vpa`](https://mathworks.com/help/symbolic/vpa.html) is used.

The subroutine `mp_switch` manages which toolbox is called and at what precision (see below).

We generally aim to simulate mixed-precision with these multiprecision toolboxes, keeping the following in mind:

* We switch to quad (or the user-specified) precision for quantities that feed into Cholesky, and we perform Cholesky in quad.  We aim to study the circumstances under which we can lift condition number bounds on $X$ for routines that use the block Pythagorean theorem.
* To simulate performance-optimized implementations, we only store vectors in double and recast them to quad whenever they are needed for computations feeding into Cholesky.

## Gram-Schmidt Routines

Each skeleton and muscle can be run individually or via the drivers

```matlab
BGS(XX, s, skel, musc, param)
```

and

```matlab
IntraOrtho(X, musc, param)
```

for skeletons ([`BGS`](https://github.com/katlund/BlockStab/blob/master/main/BGS.m)) and muscles ([`IntraOrtho`](https://github.com/katlund/BlockStab/blob/master/main/IntraOrtho.m)), respectively.

The variable `XX` denotes a block-partitioned matrix with `m` rows, `p` block vectors, and `s` columns per block vector, i.e., $m = ps$.  `X` denotes a single block vector (or a tall-and-skinny matrix) with `m` rows and `s` columns, $s \leq m$. As for the other input variables:

* `skel` - char specifying BGS skeleton
* `musc` - char specifying intra-orthogonalization muscle
* `param`: a struct with the following optional fields:
  * `.chol`: char specifying what type of Cholesky subroutine to call for skeletons that hard-code Cholesky via a block Pythagorean trick (e.g., `bcgs_pip`, `bcgs_pio`, `bcgs_iro_ls`, `bmgs_cwy`, `bmgs_icwy`, and their reorthogonalized and multi-precision versions); default for non-MP versions is `'chol_nan'`, and `'chol_aree'` for MP
  * `.mp_package`: char specifying either `'advanpix'` or `'symbolic toolbox'` as the multiprecision package; default: `'advanpix'`
  * `.mp_pair`: a cell specifying the precision pair; the first entry refers to the primary precision, and the second to the (usually) higher secondary precision; when `mp_package` is specified as one of the MP toolboxes, then `{'double', 'quad'}` is the default; otherwise `{'single', 'double'}` is the default (which don't rely on external toolboxes)
  * `.rpltol`: scalar argument for `cgs_sror` that determines the replacement tolerance; default: 1
  * `.verbose`: boolean for whether to print intermediate loss of orthogonality (LOO) or relative residual (RelRes) per iteration; default: 0

For all currently implemented skeletons and muscles, see [`main/skeletons`](https://github.com/katlund/BlockStab/tree/master/main/skeletons) and [`main/muscles`](https://github.com/katlund/BlockStab/tree/master/main/muscles). Try, for example:

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

Expected output (which will vary up to machine-precision errors for different machines and operating systems):

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

Expected output:

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

There are two main test drivers, [`MakeHeatmap`](https://github.com/katlund/BlockStab/blob/master/tests/MakeHeatmap.m) and [`RunKappaPlot`](https://github.com/katlund/BlockStab/blob/master/tests/RunKappaPlot.m).  `MakeHeatmap` generates heatmaps comparing loss of orthogonality and residual across many skeleton-muscle combinations for the same test matrix.  `RunKappaPlot` plots loss of orthogonality and residual trends against matrices with a range of condition numbers. As the Greek letter $\kappa$ is used to denote the 2-norm condition number of a matrix, we refer to these plots as "kappa plots."  See the header for each for full descriptions of their functionalities.  To explore some interesting examples, try the following, and note that `[]` (empty) arguments call default options, which can be quite handy:

* `MakeHeatmap([100 10 2], 'stewart',  {'BCGS', 'BCGS_IRO', 'BCGS_SROR'}, {'CGS', 'HouseQR'}, 1, 1);`
* `RunKappaPlot('laeuchli', [], 'demo.json');`

### New test driver

A major difference compared to the previous version is that all $\kappa$ plots are now managed by a single driver, `RunKappaPlot`, which takes three arguments:

* `mat_type`: `'default'`, `'glued'`, `'laeuchli'`, and `'monomial'`
* `options`: a struct with fields pertaining to the size and scale of trial matrices, as well as flags for saving figures and generating a TeX report
* `config_aile`: a JSON file processed by the subroutine [`alg_config`](https://github.com/katlund/BlockStab/blob/master/tests/auxiliary/alg_config.m)

Setting up the JSON configuration file is a bit tricky, but a number of templates are included.  [`demo.json`](https://github.com/katlund/BlockStab/blob/master/tests/alg_config/demo.json) in particular demonstrates all possible quirks.  There are multiple upsides to the more abstract configuration file:

* Redundancies in skeleton-muscle combinations are avoided, in particular, for skeletons like `bcgs_sror` and `bcgs_iro_ls`, which take only one muscle or none, respectively.
* Direct comparisons between muscles and skeleton-muscle algorithms can be made.  For example, we can plot `cholqr` and `bmgs`$\circ$`houseqr` on the same $\kappa$ plot.
* Direct comparisons between different implementations of the same algorithm are possible.  For example, Cholesky-based routines can be implemented with `chol_nan` or `chol_aree`.  `demo.json` encodes both configurations for several algorithms.
* Direct comparisons between multi-precision and standard double precision implementations of algorithms are also possible.  Again, see `demo.json` for examples.

If you are unfamiliar with JSON, have a look at [JSON formatter](https://jsonformatter.org/json-parser), which provides a nice GUI for parsing and formatting files.

### New plot customization features

`roadmap.m` demonstrates how to use subroutines `gen_plots.m` and `mod_run_data.m` to plot subsets of a previous run, without re-running the tests.  This can be beneficial when running a huge panel of tests that are not easily displayed on the same axis.  (Indeed, choosing colors and symbols for making lines legible is nontrivial.)  We used this feature in [Carson, et al. 2024 B](TBD) to study a progression of method modifications.

## Documentation

Each file contains a descriptive header.  See especially the following core files:

* `BGS.m` - switches between skeletons
* `IntraOrtho.m` - switches between muscles
* `InnerProd.m` - switches between inner products
* `MakeHeatmap.m` - generates heatmaps for skeleton-muscle combinations; verbose = `true` prints tables to screen
* `RunKappaPlot.m` - generates kappa plots for muscles and skeleton-muscle combinations

## How to cite us

Several works are associated with this repository:

* [Carson, et al. 2021](https://doi.org/10.1137/21M1394424): Carson, E., Lund, K., and Rozložník, M.  The stability of block variants of classical Gram-Schmidt.  SIAM Journal on Matrix Analysis and Applications. Vol. 42, 3, pp 1365--1380, 2021. DOI: 10.1137/21M1394424.
* [Carson, et al. 2024 A]: Carson, E., Lund, K., Ma, Y., and Oktay, E.  Reorthogonalized Pythagorean variants of block classical Gram-Schmidt.  In preparation, 2024. DOI: TBD.
* [Carson, et al. 2024 B]: Carson, E., Lund, K., Ma, Y., and Oktay, E.  On the loss of orthogonality of low-synchronization variants of reorthogonalized block Gram-Schmidt.  In preparation, 2024. DOI: TBD.
* [Oktay 2024]: Ph.D. thesis. Faculty of Mathematics and Physics, Charles University, Prague, 2024.
* [Oktay & Carson 2023](https://doi.org/10.1002/pamm.202200060): Okay, E. and Carson, E.  Using mixed precision in low-synchronization reorthogonalized block classical Gram-Schmidt. PAMM. Vol 23, 1, pp e202200060, 2023.  DOI: 10.1002/pamm.202200060

If you are using results from a specific paper, please cite the paper and the version of this software you are using by referring, e.g., to a tag or specific commit.

* [v1.2022](https://github.com/katlund/BlockStab/releases/tag/v1.2022): [Carson, et al. 2021](https://doi.org/10.1137/21M1394424) and [Carson, et al. 2022](https://doi.org/10.1016/j.laa.2021.12.017)
* [v1.2022.mp](https://github.com/katlund/BlockStab/releases/tag/v1.2022.mp): [Oktay & Carson 2023](https://doi.org/10.1002/pamm.202200060).
* [v2.2024-beta](https://github.com/katlund/BlockStab/releases/tag/v2.2024-beta): [Oktay 2024]
* [v2.2024](https://github.com/katlund/BlockStab/releases/tag/v2.2024): [Carson, et al. 2024 A]
* [v3.2024](https://github.com/katlund/BlockStab/releases/tag/v3.2024): [Carson, et al. 2024 B]

To cite this package in general, please use the following format:

```tex
@misc{LunOCetal24,
  title = {{BlockStab}},
  author = {Lund, Kathryn and Oktay, Eda and Carson, Erin C. and Ma, Yuxin},
  year = {2024},
  url = {https://github.com/katlund/BlockStab}
}
```

## How we cite things

Several papers are foundational for our subroutines.  We provide the full citations here and use abbreviated ones (given as [Author YYYY]) throughout the documentation.

* [Barlow 2019](https://doi.org/10.1137/18M1197400): Barlow, J. Block modified Gram-Schmidt algorithms and their analysis. SIAM Journal on Matrix Analysis and Applications. Vol. 40, 4, pp. 1257--1290, 2019. DOI: 10.1137/18M1197400.
* [Barlow & Smoktunowicz 2013](https://doi.org/10.1007/s00211-012-0496-2): Barlow, J. & Smoktunowicz, A. Reorthogonalized block classical Gram-Schmidt. Numerische Mathematik. Vol 123, pp 395--423, 2013. DOI: 10.1007/s00211-012-0496-2.
* [Bielich, et al. 2022](https://doi.org/10.1016/j.parco.2022.102940): Bielich, D, Langou J., Thomas, S., Świrydowicz, K., Yamazaki, I., and Boman, E.G.  Low-synch Gram–Schmidt with delayed reorthogonalization for Krylov solvers.  Parallel Computing. Vol 112, pp 102940, 2022. DOI: 10.1016/j.parco.2022.102940.
* [Fukaya, et al. 2020](https://doi.org/10.1137/18M1218212): Fukaya, T., Kannan, R., Nakatsukasa, Y., Yamamoto, Y., & Yanagisawa, Y. Shifted CholeskyQR for computing the QR factorization of ill-conditioned matrices. SIAM Journal on Scientific Computing. Vol. 42, 1, pp A477--A503, 2020. DOI: 10.1137/18M1218212.
* [Fukaya, et al. 2014](https://doi.org/10.1109/ScalA.2014.11): Fukaya, T., Nakatsukasa, Y., Yanagisawa, Y., & Yamamoto, Y. CholeskyQR2: A simple and communication-avoiding algorithm for computing a tall-skinny QR factorization on a large-scale parallel system. 2014 5th Workshop on Latest Advances in Scalable Algorithms for Large-Scale Systems, pp 31--38, 2014. DOI: 10.1109/ScalA.2014.11.
* [Smoktunowicz, et al. 2006](https://doi.org/10.1007/s00211-006-0042-1): Smoktunowicz, A., Barlow, J., and Langou, J. A note on the error analysis of classical Gram-Schmidt. Numerische Mathematik. Vol. 105, 2, pp 299-313, 2006. DOI: 10.1007/s00211-006-0042-1.
* [Stathopoulos & Wu 2002](https://doi.org/10.1137/S1064827500370883): Stathopoulos, A. & Wu, K. A block orthogonalization procedure with constant synchronization requirements. SIAM Journal on Scientific Computing. Vol 23, 6, pp 2165--2182, 2002. DOI: 10.1137/S1064827500370883.
* [Stewart 2008](https://doi.org/10.1137/070682563): Stewart, G. W. Block Gram-Schmidt orthogonalization. SIAM Journal on Scientific Computing. Vol 31, 1, pp 761--775, 2008. DOI: 10.1137/070682563.
* [Swirydowicz, et al. 2021](https://doi.org/10.1002/nla.2343): Świrydowicz, K., Langou, J., Ananthan, S., Yang, U., & Thomas, S. Low synchronization Gram-Schmidt and generalized minimal residual algorithms. Numerical Linear Algebra with Applications, Vol 28(2), pp e2343, 2021. DOI: 10.1002/nla.2343.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

If you are interested in helping develop an open-source version of this package in either Python or Julia, please contact [Kathryn Lund](mailto:kathryn.d.lund@gmail.com) directly.

## Related projects

See [LowSyncBlockArnoldi](https://gitlab.mpi-magdeburg.mpg.de/lund/low-sync-block-arnoldi) for block Arnoldi and GMRES versions of these routines and a simple algorithm comparison workflow.
