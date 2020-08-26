# BlockStab

## Installation
Follow the download options from the Git repository main page.

## Usage
The main purpose of this software is to study, verify, and conjecture the
stability properties of different versions of Block Gram-Schmidt (BGS) and
Block GMRES (BGMRES) via a skeleton-muscle paradigm.

Common parameters:
* `m` - number of rows
* `p` - number of block vectors
* `s` - number of columns per block vector
* `n` - p*s
* `mat` - char specifying matrix type
* `skel` - char specifying BGS skeleton
* `musc` - char specifying intra-orthonormalization muscle
* `rpltol` - replacement tolerance (only for `cgs_sror` and `bgs_sror`)
* `verbose` - true to print information to screen; false to mute

To debug a specific skeleton or muscle, set `verbose = true`.  The loss of
orthogonality (LOO) and relative residual (RelRes) will print to screen per step
of the algorithm.  Try, for example

```
mgs(randn(100,10), true);
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

```
bmgs(randn(100,20), 2, 'HouseQR', true);
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

There are several test files.  See the header for each.  To explore some
interesting examples, try the following:
* `MakeHeatmap([100 10 2], 'stewart',  {'BCGS', 'BCGS_IRO', 'BCGS_SROR'}, {'CGS', 'HouseQR'}, 1, 1)`
* `KappaPlot([100 10], [], {'MGS', 'MGS_SVL', 'MGS_LTS', 'MGS_CWY', 'MGS_ICWY'})`
* `BlockKappaPlot([100 20 2], [], {'BCGS', 'BCGS_IRO', 'BCGS_IRO_LS'}, {'CGS', 'MGS', 'HouseQR'})`
* `GluedKappaPlot([], [], {'CGS', 'CGS_P', 'MGS', 'CGS_IRO', 'CGS_IRO_LS'})`
* `GluedBlockKappaPlot([], [], {'BCGS', 'BCGS_PIP', 'BCGS_PIO'}, 'HouseQR')`
* `GluedBlockKappaPlotVaryS([], [], [2 5], {'BCGS', 'BMGS'})`
* `MonomialBlockKappaPlot([1000 100],[1 2 5 10],{'BCGS', 'BMGS', 'BCGS_IRO'}, {'CholQR', 'MGS', 'HouseQR'})`

## Documentation
Each file contains a descriptive header.  See especially the following core files:
* `MakeHeatmap.m` - generates heatmaps for skeleton-muscle combinations; verbose = `true` prints tables to screen
* `KappaPlot.m` - generates kappa plots for muscles
* `BlockKappaPlot.m` - generates kappa plots for skeleton-muscle combinations
* `MatGen.m` - generates matrices used for tests
* `BGS.m` - switches between skeletons
* `IntraOrtho.m` - switches between muscles
* `InnerProd.m` - switches between inner products

## Contributing
Pull requests are welcome. For major changes, please open an issue first to
discuss what you would like to change.

## How we cite things
Several papers are foundational for our subroutines.  We provide the full citations here and use abbreviated ones (given as [Author YYYY]) throughout the documentation.
* [Barlow 2019]: Barlow, J. Block modified Gram-Schmidt algorithms and their analysis. SIAM Journal on Matrix Analysis and Applications. Vol. 40, 4, pp. 1257--1290, 2019.
* [Barlow & Smoktunowicz 2013]: Barlow, J. & Smoktunowicz, A. Reorthogonalized block classical Gram-Schmidt. Numerische Mathematik. Vol 123, pp 395--423, 2013.
* [Fukaya, et. al. 2018]: Fukaya, T., Kannan, R., Nakatsukasa, Y., Yamamoto, Y., & Yanagisawa, Y. Shifted CholeskyQR for computing the QR factorization of ill-conditioned matrices. arXiv 1809.11085, 2018.
* [Fukaya, et. al. 2014]: Fukaya, T., Nakatsukasa, Y., Yanagisawa, Y., & Yamamoto, Y. CholeskyQR2: A simple and communication-avoiding algorithm for computing a tall-skinny QR factorization on a large-scale parallel system. 2014 5th Workshop on Latest Advances in Scalable Algorithms for Large-Scale Systems, pp 31--38, 2014.
* [Smoktunowicz, et al. 2006]: Smoktunowicz, A., Barlow, J., and Langou, J. A note on the error analysis of classical Gram-Schmidt. Numerische Mathematik. Vol. 105, 2, pp 299-313, 2006.
* [Stathopoulos & Wu 2002]: Stathopoulos, A. & Wu, K. A block orthogonalization procedure with constant synchronization requirements. SIAM Journal on Scientific Computing. Vol 23, 6, pp 2165--2182, 2002.
* [Stewart 2008]: Stewart, G. W. Block Gram-Schmidt orthogonalization. SIAM Journal on Scientific Computing. Vol 31, 1, pp 761--775, 2008.
* [Swirydowicz, et al. 2020]: Świrydowicz, K., Langou, J., Ananthan, S., Yang, U., & Thomas, S. Low synchronization Gram-Schmidt and GMRES algorithms. Technical report, 2020.

## How to cite us
[TBD]

## License
Creative Commons Attribution 4.0 International Public License: https://creativecommons.org/licenses/by/4.0/legalcode.
