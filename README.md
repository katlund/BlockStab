# BlockStab

## Installation

Follow the download options from the Git repository main page.  Then navigate to the repo and run `install_blockstab.m` in MATLAB.  Note that this script only temporarily saves the paths; they will be cleared at the next start-up.  To permanently save `BlockStab` routines to the startup path, run `savepath`, which may overwrite paths to other functions with the same names.

## Usage

The main purpose of this software is to study, verify, and conjecture the
stability properties of different versions of Block Gram-Schmidt (BGS) via a skeleton-muscle paradigm.  Each skeleton and muscle
can be run individually or via the drivers

```matlab
BGS(XX, s, skel, musc, rpltol, verbose)
```

and

```matlab
IntraOrtho(X, musc, rpltol, verbose)
```

for skeletons and muscles, respectively.

The variable `XX` denotes a matrix with `m` rows, `p` block vectors,
and `s` columns per block vector, i.e., $m = ps$.  `X` denotes a single block vector
(or tall-and-skinny matrix) with `m` rows and `s` columns, $s \leq m$. As for the other parameters:

* `skel` - char specifying BGS skeleton
* `musc` - char specifying intra-orthogonalization muscle
* `rpltol` - replacement tolerance (only required for `cgs_sror` and `bgs_sror`)
* `verbose` - `true` to print the loss of orthogonality (LOO) and relative residual (RelRes) to screen per step of the algorithm; `false` to mute

For a list of all currently implemented skeletons and muscles, see the headers to `BGS` and
`IntraOrtho`. Try, for example:

```matlab
mgs(randn(100,10), true);
```

or

```matlab
IntraOrtho(randn(100,10), 'MGS', [], true);
```

Example output:

```matlab
         LOO      |    RelRes
-----------------------------------
  1:  2.2204e-16  |  3.7634e-17
  2:  2.3515e-16  |  4.1141e-17
  3:  2.3984e-16  |  9.1624e-17
  4:  2.4010e-16  |  8.9259e-17
  5:  2.2531e-16  |  1.1044e-16
  6:  2.3365e-16  |  1.2696e-16
  7:  2.7442e-16  |  1.3645e-16
  8:  2.7454e-16  |  1.4435e-16
  9:  2.6734e-16  |  1.5960e-16
 10:  2.7682e-16  |  1.5283e-16
 ```

For block methods, try

```matlab
bmgs(randn(100,20), 2, 'HouseQR', true);
```

or

```matlab
BGS(randn(100,20), 2, 'BMGS', 'HouseQR', [], true);
```

Example output:

```matlab
         LOO      |    RelRes
-----------------------------------
  1:  6.7663e-16  |  1.2409e-16
  2:  5.8765e-16  |  1.5998e-16
  3:  7.9980e-16  |  1.5196e-16
  4:  8.0337e-16  |  1.9282e-16
  5:  8.0743e-16  |  1.9774e-16
  6:  8.1551e-16  |  2.4859e-16
  7:  8.8906e-16  |  2.5308e-16
  8:  8.7674e-16  |  2.6215e-16
  9:  8.8764e-16  |  2.6683e-16
 10:  9.0044e-16  |  2.5827e-16
 ```

There are several test files.  `MakeHeatmap` generates heatmaps comparing loss of
orthogonality and residual across many skeleton-muscle combinations for the same
test matrix.  `KappaPlot` and similar test files plot loss of orthogonality and
residual trends against matrices with a range of condition numbers. As the Greek
letter $\kappa$ is used to denote the 2-norm condition number of a matrix, we refer to these
plots as "kappa plots."  See the header for each for full descriptions of
their functionalities.  To explore some interesting examples, try the following:

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
* `MatGen.m` - generates matrices used for tests
* `CreateGluedMatrix.m` - generates glued matrices
* `MakeHeatmap.m` - generates heatmaps for skeleton-muscle combinations; verbose = `true` prints tables to screen
* `KappaPlot.m` - generates kappa plots for muscles
* `BlockKappaPlot.m` - generates kappa plots for skeleton-muscle combinations

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

Pull requests are welcome. For major changes, please open an issue first to
discuss what you would like to change.

## Related projects

See [LowSyncBlockArnoldi](https://gitlab.mpi-magdeburg.mpg.de/lund/low-sync-block-arnoldi) for block Arnoldi versions of these routines and a simple algorithm comparison workflow.

## License

[Creative Commons Attribution 4.0 International Public License](https://creativecommons.org/licenses/by/4.0/legalcode)
