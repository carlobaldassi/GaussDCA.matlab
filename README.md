Gaussian Direct Coupling Analysis for protein contacts predicion
================================================================

Overview
--------

This is the code which accompanies the paper ["Fast and accurate multivariate
Gaussian modeling of protein families: Predicting residue contacts and
protein-interaction partners"](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0092721)
by Carlo Baldassi, Marco Zamparo, Christoph Feinauer, Andrea Procaccini,
Riccardo Zecchina, Martin Weigt and Andrea Pagnani, (2014)
PLoS ONE 9(3): e92721. doi:10.1371/journal.pone.0092721

This code, except the `kseq.h` file, is released under the GPL version 3 (or
later) license; see the `LICENSE.md` file for details.

The code is written in [MATLAB](http://www.mathworks.it/products/matlab/); it
provides a function which reads a multiple sequence alignment (in FASTA format)
and returns a ranking of all pairs of residue positions in the aligned
amino-acid sequences.

The code uses some MATLAB MEX files written in C, which need to be compiled; in
case a compiler is not available, it will still work, but it will use (much)
slower fallbacks written in MATLAB code.

Compiling the MEX files
-----------------------

On Unix systems (Linux, Mac OS X, ...) the modules can be compiled following
these steps:

  1. Edit the `Make.config` file in order to set the MATLAB version and/or
     installation directory; in most common cases, setting the MATLAB version
     should be sufficient
  2. Make sure that the zlib library, complete with header files, is installed
     in the system (e.g. on Ubuntu Linux you'll need the `zlib-dev` package).
  3. Run `make` in the top directory

The code can use multiple cores via the OpenMP framework, which however requires
to be supported by the compiler (e.g. GNU gcc supports it, LLVM's clang does not
at the time of writing). You can disable the feature by setting `OMP=0` in the
`Make.config` file.

The code is optimized for 64-bit systems; on 32-bit systems, you may get better
performance my setting `PACKBITS=32` in the `Make.config` file.

Usage
-----

This software provides one main function, `gDCA(filname::String, ...)`. This
function takes the name of a (possibly gzipped) FASTA file, and returns
predicted contact ranking, in the form of a matrix with three columns. The first
two columns containin indices `i` and `j` (with `i < j`), the third one
contaitns a score. The indices start counting from 1, and denote pair of
residue positions in the given alignment; pairs which are separated by less than
a given number of residues (by default 5) are filtered out. The rows of the matrix
are sorted by score in descending order, such that predicted contacts should come
up on top.

The `gDCA` function takes some additional, optional keyword arguments:

 * `pseudocount`: the value of the pseudo-count parameter, between `0` and `1`.
                  the default is `0.8`, which gives good results when the
                  Frobenius norm score is used (see below); a good value for the
                  Direct Information score is `0.2`.
 * `theta`: the value of the similarity threshold. By default it is `'auto'`,
            which means it will be automatically computed (this takes additional
            time); otherwise, a real value between `0` and `1` can be given.
 * `max_gap_fraction`: maximum fraction of gap symbols in a sequence; sequences
                       which exceed this threshold are discarded. The default
                       value is `0.9`.
 * `score`: the scoring function to use. There are two possibilities, `'DI'` for
            the Direct Information, and `'frob'` for the Frobenius norm. The
            default is `'frob'`.
 * `min_separation`: the minimum separation between residues in the output
                     ranking. Must be >= `1`. The default
                     is `5`.

Using gzipped FASTA files is only supported if the MEX module for fasta reading
is compiled.

Examples
--------

Here is a basic usage example, assuming an alignment in FASTA format is found
in the file "alignment.fasta.gz", and the MEX modules were compiled:

  ```
  >> FNR = gDCA('alignment.fasta.gz');
  ```

The above uses the Frobenius norm ranking with default parameters.
This is how to get the Direct Information ranking instead:

  ```
  >> DIR = gDCA('alignment.fasta.gz', 'pseudocount', 0.2, 'score', 'DI');
  ```

