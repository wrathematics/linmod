```
    _ _                           _ 
    | (_)_ __  _ __ ___   ___   __| |
    | | | '_ \| '_ ` _ \ / _ \ / _` |
    | | | | | | | | | | | (_) | (_| |
    |_|_|_| |_|_| |_| |_|\___/ \__,_|

    __     __            _             0.1.0
    \ \   / ___ _ __ ___(_) ___  _ __  
     \ \ / / _ | '__/ __| |/ _ \| '_ \ 
      \ V |  __| |  \__ | | (_) | | | |
       \_/ \___|_|  |___|_|\___/|_| |_|
```


## About 

This package implements linear model and generalized linear model fitters
that behave somewhat in the spirit of R's lm.fit() and glm.fit() routines.
The implementations here are new and not based on the implementations in R,
except that we allow Ross Ihaka's "limited pivoting" strategy as an option
in the fit of linear models in the (possibly) rank degenerate case.  We
also have a much faster (2x+) fitter available if the user wishes to 
assume that the model is of full column rank.


The linmod package is implemented (almost) entirely in Fortran 2003, but a C 
interface is provided with the library.  Additionally, R bindings are provided.  
At the moment, there is some very meager OpenMP usage, which I hope to expand 
over time.


The LM implementations are pretty standard, and the majority of the hard work
is handled by the LAPACK authors (with some modifications).  More specifically, 
the linear model fitter uses a QR factorization of the input data x, and then
generates the solution to the least squares problem by solving the (triangular) 
system of equations:

```
Rb = Q'y
```

See pretty much any book on numerical linear algebra for details; though
I recommend:

    Golub G. H. and Van Loan, C. F. "Matrix Computations", 3rd ed.

The GLM implementation uses iteratively reweighted least squares (IRLS), 
based heavily on the descriptions available in the book: 

    McCallagh, P. and Nelder, J. "Generalized Linear Models", 2nd ed.


Outputs of the available routines are designed to mimic those of R's 
lm.fit() and glm.fit(), and as such should be fully compatible with R.
Additionally, R's outputs were used extensively to test the validity of
the new implementations, which was absolutely invaluable.  And while I
have not studied the implementation details of R's methods, I suspect
that is where the similarities end.  R's implementations use a great deal 
of R code (especially in fitting GLM's), and so I would expect this 
package to manage memory much better.  Also, R uses LINPACK routines for
their linear model fitter; this package uses LAPACK routines (modified
dgels et. al).  As such, I would suspect that its performance gains would
be greater than R's when linking with a high performance LAPACK library 
such as MKL or ACML; but I have not tested this at this time.  Finally, the 
obvious benefit is that the library itself does not even need R, while R's
lm.fit() and glm.fit() do indeed depend quite heavily on R.


Early benchmarking shows a ~15% performance gain over R's lm.fit() at
modestly sized problems when using the rank-revealing method, and a 200% 
performance gain when assuming the model is of full column rank.

Matrices in this library are assumed to be column-major. Use the provided 
transpose functions to convert between the two storage formats as needed.



## License and Copying

Copyright 2013-2014, Drew Schmidt.  All rights reserved.

Currently, the project is licensed and distributed under the Mozilla Public
License 2.0 (MPL).  See the file LICENSE for details.  I anticipate moving
the project to a more permissive license when the library becomes more mature.



## Requirements and Installation

To install, you will need: 

* cmake >= 2.8.1
* A Fortran 2003 compatible compiler with OpenMP support.
* A C99 compatible compiler
* LAPACK and BLAS libraries (not needed if installing the R package)
* R >= 2.14.0 and the RNACI package (if installing the R package)

Both the R package and the standalone library require cmake, because if you
so much as think the word "autotools" around me, I'll punch you in the 
stomach.

To install the R package, simply execute:

```
R CMD INSTALL linmod_0.1.0.tar.gz
```

To build just the shared library, in your terminal, execute:

```
cd linmod/src/linmod/ 
make
```

A static and dynamic library will be placed in the

```
linmod/src/linmod/build
```

tree.



## Usage 

If using the R package, simply use lm_fit() and glm_fit() *exactly* as you
would use lm.fit() and glm.fit().  There are some caveats for the glm fitter;
specifically, the family and link have to each be available to Fortran (R is
slightly more flexible).  Currently, the supported families and 
their associated links are:

* Binomial: cloglog, log, logit, probit, cauchit
* Gamma: identity, inverse, log
* Gaussian: identity, inverse, log
* Poisson: identity, log, sqrt
* Inverse gaussian: inverse, log, identity, "1/mu^2"

Some notable families and links available to R's glm.fit() are missing; those 
will eventually be supported.  But my opinion is that the most important/common
families and links are already supported.

If using the shared library, most of the above applies to you, so you would
do well to read it anyway.  Here, there are both C and Fortran interfaces.

TODO: say more about interfaces



## Contact

Drew Schmidt:

* Project home: https://github.com/wrathematics/linmod
* Bug reports: https://github.com/wrathematics/linmod/issues
* Email: wrathematics .AT. gmail.com
* Twitter: @wrathematics



## Testing

I have rigorously tested this package with the following software/versions:

* gfortran: 4.8.2
* clang: 3.4
* cmake: 2.8.12.2
* R: 3.1.1

Reports from other software/versions (with success or failure) are most welcome.

