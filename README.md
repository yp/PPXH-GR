
  PPXH-GR
===========

A heuristic algorithm for the Pure Parsimony Xor Haplotyping problem
(PPXH) based on Graph Realization.

by [Yuri Pirola](http://bimib.disco.unimib.it/index.php/Pirola_Yuri)


------------------------------------------------------------------------



## Introduction ##

This program implements a heuristic for the _Pure Parsimony Xor
Haplotyping_ problem based on a fast algorithm for the _Graph
Realization_ problem.

The heuristic is described in the following papers:

* Paola Bonizzoni, Gianluca Della Vedova, Riccardo Dondi, Yuri Pirola,
  and Romeo Rizzi.
  **Pure Parsimony Xor Haplotyping**.
  _IEEE/ACM Transactions on Computational Biology and Bioinformatics_.
  Vol. 7(4), pp. 598-610 (Oct-Dec, 2010).
  [link](http://dx.doi.org/10.1109/TCBB.2010.52)
  [arXiv version](http://arxiv.org/abs/1001.1210)
* Paola Bonizzoni, Gianluca Della Vedova, Riccardo Dondi, Yuri Pirola,
  and Romeo Rizzi.
  **Pure Parsimony Xor Haplotyping**.
  _Proc. of 5th Int. Symp. Bioinformatics Research and Applications
  (ISBRA)_.
  Springer. Vol. 5542, pp. 186-197 (May, 2009).
  [link](http://dx.doi.org/10.1007/978-3-642-01551-9_19)



## Download ##

PPXH-GR is currently distributed only on source form.
PPXH-GR source code is available at the `yp/PPXH-GR` Git repository
hosted by GitHub.
The repository can be explored using the GitHub web interface at
<https://github.com/yp/PPXH-GR>.

The latest stable version of PPXH-GR can be downloaded in either
[.tar.gz](https://github.com/yp/PPXH-GR/tarball/master) or in
[.zip](https://github.com/yp/PPXH-GR/zipball/master) format.
Previous stable releases can be downloaded from
<https://github.com/yp/PPXH-GR/archives/master>.

It is also possible to clone the entire repository using the following
command:

    $ git clone git://github.com/yp/PPXH-GR.git

Or, if you have a GitHub account, you can fork the project from the
[repository web page](https://github.com/yp/PPXH-GR).



## Dependencies ###

Compilation dependencies:

- [GNU Make](http://www.gnu.org/s/make/)
- [GNU Scientific Library (GSL)](http://www.gnu.org/s/gsl/)


Run-time dependencies:

- [GREAL](http://acgt.cs.tau.ac.il/greal/), a software for the _Graph
  Realization_ problem.
  GREAL is only distributed as executable and needs that the following
  (old) library [`libstdc++-libc6.2-2.so.3`][link] is placed in the working
  directory
- Basic Linux/Unix utilities, such as `bash`, `grep`, and `sed`

[link]: http://rpmfind.net/linux/rpm2html/search.php?query=libstdc%2B%2B-libc6.2-2.so.3



## Compilation ###

The program can be compiled by issuing the command at the command
prompt:

    $ make all

The resulting executable, `ppxh-gr`, can then be copied (if necessary)
in some other folder.



## Usage ##

The program takes as input a genotype matrix and returns a haplotype
matrix that explains the given genotype matrix.
The command line is as follow:

    $ ./ppxh-gr -f FILENAME [-r N]

where `FILENAME` indicates the name of the file containing the genotype
matrix, and `N` indicates the number of solutions that PPXH-GR computes
(default: `10`).

Please note that the following files must be placed in the working
directory:

- `exec-gr.sh` (supplied with PPXH-GR)
- `grealLinux`
- `libstdc++-libc6.2-2.so.3`


The input genotype matrix is described as a sequence of _n_ rows
(representing the genotypes of the _n_ individuals) each one containing
exactly _m_ characters (representing the alleles at the _m_ SNPs/loci
over which the xor-genotypes are defined).
Each character is either `0`/`1` (indicating an homozygous locus) or `2`
(indicating an heterozygous locus).
Please refer to the file `genot-matrix-test.txt` for an example of a
simple 7x4 genotype matrix.



## License ##

PPXH-GR is released under the terms of the GNU General Public License
(GPL) as published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

PPXH-GR is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

Please refer to file `COPYING` or to the
[GNU website](http://www.gnu.org/licenses/) for a copy of the GNU
General Public License.



## Contacts ##

Please contact *Yuri Pirola* for additional information.  
E-mail:   <yuri.pirola@gmail.com>  
Web page: <http://bimib.disco.unimib.it/index.php/Pirola_Yuri>
