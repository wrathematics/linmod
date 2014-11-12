# Notes


## lm_fit() internals.

The LAPACK routine `dgels` was deeply modified (into `rdgels`) to use a 
slightly modified modified `dgeqp3` (`rdgeqp3`) for performing a 
rank-revealing QR with column-pivoting (or optionally `dgeqrf` without 
pivoting, but only for full rank problems).  `rdgeqp3` calls `rdlaqp2`,
a highly modified `dlaqp2`, to perform the actual QR+rank detection.

Overwhelmingly, the computation for `rdgels` is in the application
of the elementary reflector via the call to `dlarf` in `rdlaqp2`.

