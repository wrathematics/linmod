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


Matrices in this library are assumed to be column-major. Use the provided 
transpose functions to convert between the two storage formats if/as needed.



## License and Copying

Copyright 2013-2014, Drew Schmidt.  All rights reserved.

Currently, the project is licensed and distributed under the Mozilla Public
License 2.0 (MPL).  See the file LICENSE for details.  I anticipate moving
the project to a more permissive license when the library becomes more mature.



## System Requirements

To install, you will need: 

* cmake >= 2.8.1
* A Fortran 2003 compatible compiler with OpenMP support.
* A C99 compatible compiler (I really recommend clang).
* LAPACK and BLAS libraries (not needed if installing the R package).
* R >= 2.14.0 and the RNACI R package (if installing the R package).

Both the R package and the standalone library require cmake, because if you
so much as think the word "autotools" around me, I'll punch you in the 
stomach.

Some effort has been made to make the library build with a Fortran
compiler that doesn't have OpenMP support.  I have serious doubts 
that it will actually work, however, and recommend you just
use a remotely recent compiler instead.  

At the moment, the C compiler does not need OpenMP support, 
basically because I don't want to deal with headaches on the Mac.
Over time, this will probably change (primarily affecting the R
bindings); hopefully clang can get its act together by then.



## Installation

To install the R package, assuming you have all of the necessary
system dependencies (stated above), the easiest way to build the
linmod R package is to use the devtools package.

```r
devtools::install_github("wrathematics/RNACI") # dependency
devtools::install_github("wrathematics/linmod")
```

To build just the standalone library without the R package, in your
terminal, execute:

```
cd linmod/src/linmod/ 
make
```

A static and dynamic library will be placed in the 
`linmod/src/linmod/build` tree.



## Usage and Interfaces

For GLM's, the currently supported families and their associated 
links are:

* Binomial: cloglog, log, logit, probit, cauchit
* Gamma: identity, inverse, log
* Gaussian: identity, inverse, log
* Poisson: identity, log, sqrt
* Inverse gaussian: inverse, log, identity, "1/mu^2"

#### R Bindings

Simply use `lm_fit()` and `glm_fit()` as you would R's own
`lm.fit()` and `glm.fit()` respectively.  For a higher level
interface (like `lm()` and `glm()`), see the Rglm2 package.

#### C/Fortran Interface

Ironically, this is fairly unpolished at the moment.  More information
will become available when the interface finalizes.



## Q&A

* Why R bindings?
  - R is what I know and love.  Although I do like to gaze longingly
    at Julia, for the time being, I'm still committed to R.
* Will this ever be on CRAN?
  - I highly doubt it.  The Mac and Windows CRAN machines do
    not support cmake, which is very much needed to build the
    underlying Fortran library.  Additionally, CRAN in general
    frowns on Fortran other than F77.  Honestly, it's not really
    worth the trouble to me.
* Why cmake?
  - Fundamentally, this is a Fortran library, and I make extensive
    use of some of the nice F95 standard features (plus a few from
    F2003).  Using standard make (this ends up including R's 
    Makevars), it becomes an incredibly tedious chore to manage
    Fortran module dependencies.  Because cmake has been
    competently engineered in almost every way, Fortran module
    dependencies are automatically discovered.
* Why Fortran?
  - Because it offers the high-level stuff I want, and it tends to
    perform much better than even other so-called high performance
    languages.  Also, the linear model fitter is a bunch of modifications
    to LAPACK routines, so it was going to include a lot of Fortran
    regardless.
* Why C99?
  - For `stdbool.h`, probably some other stuff.
* Is this production ready?
  - Hell no; it's not even done.  Even the (library) interfaces may change
    somewhat over the next few releases.



## Testing

I have rigorously tested this package with the following software/versions:

* gfortran: 4.8.2
* clang: 3.4
* cmake: 2.8.12.2
* R: 3.1.1

Reports from other software/versions (with success or failure) are most welcome.



## Contact

Drew Schmidt:

* Project home: https://github.com/wrathematics/linmod
* Bug reports: https://github.com/wrathematics/linmod/issues
* Email: wrathematics .AT. gmail.com
* Twitter: @wrathematics

