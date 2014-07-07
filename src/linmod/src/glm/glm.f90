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


module glm
  use, intrinsic :: iso_c_binding
  
  use :: lapack
  use :: glm_check
  use :: glm_loglik_utils
  use :: glm_link_utils
  use :: glm_mu_var
  use :: glm_update_utils
  
  implicit none
  
  
  contains
  
  subroutine glm_fit(family, link, intercept, stoprule, n, p, x, y, &
                     beta, wt, offset, resids, maxiter, tol, info) &
  bind(c, name='glm_fit_')
    ! in/out
    integer, intent(in) :: family, link, intercept, stoprule
    integer, intent(in) :: n, p, maxiter
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
    double precision, allocatable :: sqwt(:)
    double precision, allocatable :: x_tw(:,:)
    double precision, allocatable :: work(:)
    ! intrinsic
    intrinsic :: min, max, dble, dsqrt
    
    
    ! quick return if possible
    info = glm_check_fam_link(family, link)
    if (info < 0) return
    
    info = check_response(family, n, y)
    if (info == -8) return
    
    ! allocate local storage
    allocerr = 0
    
    allocate(beta_old(p), stat=allocerr)
    if (allocerr /= 0) stop "out of memory"
    allocate(eta(n), stat=allocerr)
    if (allocerr /= 0) stop "out of memory"
    allocate(sqwt(n), stat=allocerr)
    if (allocerr /= 0) stop "out of memory"
    allocate(x_tw(n, p), stat=allocerr)
    if (allocerr /= 0) stop "out of memory"
    allocate(mu(n), stat=allocerr)
    if (allocerr /= 0) stop "out of memory"
    allocate(z(n), stat=allocerr)
    if (allocerr /= 0) stop "out of memory"
    
    ! allocate workspace for linear models
    lwork = min(n, p) + max(1, n, p)
    allocate(work(lwork), stat=allocerr)
    if (allocerr /= 0) stop "out of memory"
    
    
    ! empty model case
    if (n < 1 .or. p < 1) then
      do i = 1, n
        eta(i) = 0.0d0 + offset(i)
      end do
      
      call glm_linkinv(link, n, eta, mu)
      
      call glm_check_mu(family, n, mu, tol, info)
      if (info /= 0) return
      
      call glm_variance(family, n, mu, wt)
      
      call glm_residuals(link, n, y, mu, eta, resids)
      
      iter = 0
      
      return
    end if
    
    
    ! initialize
!!! Parallel block
    !$omp parallel private(i) default(shared) 
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
      
      
      ! compute eta = x*beta and mu = inverse_link( eta )
      if (iter == 0) then
        call glm_initial_mu(family, n, y, wt, mu)
        call glm_link(link, n, mu, eta)
      else 
        call dgemm('n', 'n', n, 1, p, 1.0d0, x, n, beta, p,  0.0d0, eta, n)
        call glm_linkinv(link, n, eta, mu)
      end if
      
      
      ! check for bad fit in the mu's
      call glm_check_mu(family, n, mu, tol, info)
      if (info /= 0) return
      
      
      ! update wt
      call glm_variance(family, n, mu, wt)
      
      
      if (stoprule == glm_stoprule_deviance) then
        dev_old = dev
        dev = glm_deviance(family, n, y, mu)
      end if
      
      
!!! Parallel block
      !$omp parallel private(i) default(shared) 
      !$omp do private(i)
        do i = 1, n
          sqwt(i) = dsqrt(wt(i))
        end do
      
      
      ! prepare lhs:  x_tw = x*wt
      !$omp do private(i, j)
        do j = 1, p
          do i = 1, n
            x_tw(i,j) = sqwt(i) * x(i,j)
          end do
        end do
      
      
      ! prepare rhs:  z = sqrt(wt) * (x*beta + 1/wt*(y-mu))
      !                 = sqwt * eta + 1/sqwt*(y-mu)
      !$omp do private(i, tmp)
        do i = 1, n
          tmp = sqwt(i)
          z(i) = tmp*eta(i) + 1.0d0/(tmp)*(y(i)-mu(i))
        end do
      !$omp end parallel
!!! End parallel block
      
      
      ! update beta:  fit z ~ x_tw
      call glm_update_beta(n, p, beta, beta_old, x_tw, z, work, lwork, info)
      
      
      ! check for convergence
      if (iter > 0) then
        converged = glm_convergence(stoprule, p, beta_old, beta, dev, dev_old, tol, iter, maxiter)
      end if
      
      if (converged == glm_convergence_converged) then
        ! converged: break loop and compute likelihood statistics and residuals
        goto 10 ! converged
      else if (converged == glm_convergence_infparams) then
        ! infinite parameter values detected: deallocate and return with error
        goto 1
      end if
      
    end do IRLS_LOOP
    
    
    !!! success --- now do all the other stuff
  10   continue
    
    ! aic, deviance, nulldeviance
    call glm_loglik_stats(family, link, intercept, n, p, x, y, eta, mu, beta, beta_old, dev, aic, nulldev)
    
    ! compute working residuals
    call glm_residuals(link, n, y, mu, eta, resids)
    
    write (*,*) "iter=",iter
    
    ! exit subroutine
  1    continue
    
    deallocate(beta_old)
    deallocate(eta)
    deallocate(x_tw)
    deallocate(mu)
    deallocate(z)
    deallocate(work)
    
    
    return
  end
  
  
end module

