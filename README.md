# BlockStab

## Installation
Follow the download options from the Git repository main page.

## Usage
The main purpose of this software is to study stability properties of different versions of Block Gram-Schmidt (BGS) and Block GMRES (BGMRES) via a skeleton-muscle paradigm.

RunTest.m and related variants are functions for running tests, and they allow for specifying different matrices and skeleton-muscle combinations.  For every matrix, RunTest produces a series of heatmaps corresponding to loss of orthonormality, backward error, and runtimes for the specified skeleton-muscle combinations.  Such a format allows one to quickly identify stability patterns for BGS variations.

Basic use-case:
`RunTest([m p s], mat, skel, musc, rpltol, verbose)`
* `m` - number of rows
* `p` - number of block vectors
* `s` - number of columns per block vector
* `mat` - char specifying matrix type
* `skel` - char specifying BGS skeleton
* `musc` - char specifying intra-orthonormalization muscle
* `rpltol` - replacement tolerance (only for `cgs_sror` and `bgs_sror`)
* `verbose` - true to print tables to screen; default is false

Examples:
* `RunTest([1000 20 10], 'laeuchli', {'BCGS_SROR', 'BMGS'}, {'CGS', 'HouseQR'})` will produce and save three 2 x 2 heatmaps: one for loss of orthogonality, one for relative backward error, and one for runtimes, all for computing a QR decomposition of a Läuchli matrix of size 1000 x 200 (partitioned implicitly into 20 block vectors each with 10 columns).
* `RunTest([1000 20 10], 'laeuchli', {'BCGS_SROR', 'BMGS'}, {'CGS', 'HouseQR'}, [], 1)` will do the same, but also print the tables to screen.

See the header of each RunTest script for the specific options that can be set.

## Documentation
Each file contains a descriptive header.  See especially the following core files:
* `RunTest.m` and variants
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
