! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt



! glm_fit
! 
! purpose
! =======
!
! Generalized linear models (glm) using irls.
!
! Implementation is largely based on:
!   McCullagh P. and Nelder, J. A. (1989) "Generalized Linear Models".
!
!
! arguments
! =========
!
! family    (input) character*8
!           the distribution of the response.  choices are:
!           binomial, gamma, gaussian, poisson.
! 
! link      (input) character*8
!           the link function to use.  choice depends on the
!           distribution.  choices are:
!           (binomial) logit, probit, log, cloglog,
!           (gamma) inverse, identity, log,
!           (gaussian) identity, log, inverse
!           (poisson) log, identity, sqrt
! 
! intercept     (input) integer
!           
! 
! stoprule  (input) character*1
!           
! 
! n         (input) integer
!           the number of rows of the data matrix x.  n>=0.
! 
! p         (input) integer
!           the number of columns of the data matrix x.  p>=0.
! 
! x         (input/output) double precision array, dimension (n,p).
!           on entry, the input data matrix.  on exit, the output
!           of a lapack qr factorization is stored for x.
! 
! y         (input) double precision array, dimension (n,1).
!           the response variable.
! 
! offset    (input)
!           
! 
! beta      (output) double precision array, dimension (p,1).
!           the model coefficients.
! 
! wt        (input/output) double precision array, dimension (n,1).
!           on input, a set of starting weights (initialize to 1.0
!           if you have nothing else in mind).  on output, the final
!           set of weights ("working weights") in the irls fit.
! 
! offset    (input) double precision array, dimension (n,1).
!           
! 
! resids    (output)
!           
! 
! maxiter   (input) integer
!           maximum number of iterations.  20 should be sufficient
!           for convergence for most applications if it's going to 
!           happen at all.
! 
! tol       (input) double precision
!           tolerance for the qr decomposition.
! 
! info      (output) integer
!           = 0: successful exit.
!           < 0: if info = -i, the i-th argument had an illegal value.


! eta = x*beta_old
! mu = logit_linkinv(eta)
!
!
! each iteration after the first, we fit the linear model:  z ~ x_tw
! where
!      x_tw = x*wt, and
!      z = sqrt(w) * (eta + 1/w * (y-mu))
!
!
! intercept = 'y', 'n', for whether intercept should be included in null model


! FIXME inputs are:
  ! x, y, weights, start, etastart, mustart, offset, ...

subroutine glm_fit(family, link, intercept, stoprule, n, p, x, y, &
                   beta, wt, offset, resids, maxiter, tol, trace, info) &
  bind(C, name='glm_fit_')
  use, intrinsic :: iso_c_binding
  
  use :: lapack
  use :: glm_check
  use :: glm_loglik_utils
  use :: glm_family_utils
  use :: glm_update_utils
  use :: glm_constants
  implicit none
  
  ! in/out
  integer, intent(in) :: family, link, intercept, stoprule
  integer, intent(in) :: n, p, maxiter, trace
  integer, intent(out) :: info
  double precision, intent(in) :: x(n,p), y(n), offset(n), tol
  double precision, intent(out) :: beta(p), wt(n), resids(n)
  ! local
  integer :: converged, i, j, iter, allocerr, rank, lwork
  double precision :: aic, dev, nulldev, dev_old, tmp
  double precision, allocatable :: beta_old(:)
  double precision, allocatable :: eta(:)
  double precision, allocatable :: mu(:)
  double precision, allocatable :: z(:)
  double precision, allocatable :: rtwt(:)
  double precision, allocatable :: x_tw(:,:)
  double precision, allocatable :: work(:)
  ! intrinsic
  intrinsic :: min, max, dble, dsqrt
  
  
  ! Quick return if possible
  info = glm_check_inputs(n, p, stoprule, maxiter, tol)
  info = glm_check_fam_link(family, link)
  info = glm_check_response(family, n, y)
  if (info < 0) return
  
  
  ! allocate local storage
  allocerr = 0
  
  allocate(beta_old(p), stat=allocerr)
  if (allocerr /= 0) goto 1
  allocate(eta(n), stat=allocerr)
  if (allocerr /= 0) goto 1
  allocate(rtwt(n), stat=allocerr)
  if (allocerr /= 0) goto 1
  allocate(x_tw(n, p), stat=allocerr)
  if (allocerr /= 0) goto 1
  allocate(mu(n), stat=allocerr)
  if (allocerr /= 0) goto 1
  allocate(z(n), stat=allocerr)
  if (allocerr /= 0) goto 1
  
  ! allocate workspace for linear models
  lwork = min(n, p) + max(1, n, p)
  allocate(work(lwork), stat=allocerr)
  if (allocerr /= 0) goto 1
  
  
  ! empty model case
  if (n < 1 .or. p < 1) then
    do i = 1, n
      eta(i) = 0.0d0 + offset(i)
    end do
    
    call glm_linkinv(link, n, eta, mu)
    
    call glm_check_fitted(family, n, mu, info)
    if (info /= 0) return
    
    call glm_variance(family, n, mu, wt)
    
    call glm_residuals(link, n, y, mu, eta, resids)
    
    iter = 0
    goto 1
  end if
  
  
  ! initialize
!!! Parallel block
  !$omp parallel if (n > 5000) private(i) default(shared) 
  !$omp do
  do i = 1, p
    beta_old(i) = 0.0d0
  end do
  !$omp end do
  
  !$omp do
  do i = 1, n
    wt(i) = 1.0d0
  end do
  !$omp end do
  !$omp end parallel
!!! End parallel block
  
  dev = 0.0d0
  
  
  !!! main loop
  IRLS_LOOP: do iter = 0, maxiter
    
    
    ! compute eta = x*beta and mu = linkinv(eta) <==> link(mu) = eta
    if (iter == 0) then
      call glm_initial_mu(family, n, y, wt, mu)
      call glm_link(link, n, mu, eta)
    else 
      call dgemm('n', 'n', n, 1, p, 1.0d0, x, n, beta, p, 0.0d0, eta, n)
      call glm_linkinv(link, n, eta, mu)
    end if
    
    
    ! check for bad fit in the mu's
    call glm_check_fitted(family, n, mu, info)
    if (info /= 0) return
    
    if (stoprule == glm_stoprule_deviance) then
      dev_old = dev
      dev = glm_deviance(family, n, y, mu)
    end if
    
    
    ! Form working response z and x_tw = x*sqrt(wt)
    call glm_update_wt(family, link, n, p, x, x_tw, wt, rtwt, y, mu, eta, z)
    
    call glm_update_z(link, n, rtwt, y, mu, eta, z, offset)
    
    
    ! fit z ~ x_tw
    call glm_update_beta(n, p, beta, beta_old, x_tw, z, work, lwork, info)
    
    
    ! check for convergence
    if (iter > 0) then
      converged = glm_check_convergence(stoprule, p, beta_old, beta, dev, dev_old, tol, iter, maxiter)
      
      if (trace /= 0) then
        if (iter == 1) write (*,*) ""
        write (*, fmt='(a, f0.5, a, i0)') "Deviance = ", dev, " Iterations - ", iter
      end if
    end if
    
    ! converged: break loop and compute likelihood statistics and residuals
    if (converged == glm_convergence_converged) then
      goto 10
    ! infinite parameter values detected: deallocate and return with error
    else if (converged == glm_convergence_infparams) then
      goto 1
    end if
    
  end do IRLS_LOOP
  
  
  10 continue
  
  ! aic, deviance, nulldeviance
  call glm_loglik_stats(family, link, intercept, n, p, x, y, eta, mu, beta, beta_old, dev, aic, nulldev)
  
  ! compute working residuals
  call glm_residuals(link, n, y, mu, eta, resids)
  
  if (converged == glm_convergence_noconvergence) iter = iter - 1
  write (*,*) "iter=",iter
  
  
  1 continue
  if (allocerr /= 0) info = glm_oom
  
  deallocate(beta_old)
  deallocate(eta)
  deallocate(x_tw)
  deallocate(mu)
  deallocate(z)
  deallocate(work)
  
  
  return
end subroutine

