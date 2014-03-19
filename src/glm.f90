! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


! glm_fit
! 
! purpose
! =======
!
! generalized linear models (glm) using irls.
!
! implementation is based on "generalized linear models, second 
! edition", by p. mccullagh and j. nelder.
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
! incpt     (input) character*1
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


!
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
! incpt = 'y', 'n', for whether intercept should be included in null model
!
subroutine glm_fit(family, link, incpt, stoprule, n, p, x, y, &
                   beta, wt, offset, resids, maxiter, tol, info)
  use glm_loglik_utils
  use glm_link_utils
  use glm_mu_var
  implicit none
  ! in/out
  character*8         family, link
  character*1         incpt
  integer             stoprule, n, p, maxiter, info
  double precision    x(n,p), y(n), offset(n), beta(p), wt(n), &
                      resids(n), tol
  ! local
  integer             converged, i, j, iter, allocerr, &
                      rank, lwork
  double precision    aic, dev, nulldev, tmp, dev_old
  double precision, allocatable :: beta_old(:)
  double precision, allocatable :: eta(:)
  double precision, allocatable :: mu(:)
  double precision, allocatable :: z(:)
  double precision, allocatable :: sqwt(:)
  double precision, allocatable :: x_tw(:,:)
  double precision, allocatable :: work(:)
  ! parameter
  double precision    zero, one
  parameter ( zero = 0.0d0, one = 1.0d0 )
  ! external
  integer             glm_convergence, check_response
!  external            glm_update_beta
  intrinsic           min, max, dble, dsqrt
  
  
  ! quick return if possible
  info = glm_check_fam_link(family, link)
  if (info.lt.0) return
  
  info = check_response(n, y)
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
  if (n.lt.1 .or. p.lt.1) then
    do i = 1, n
      eta(i) = zero + offset(i)
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
  do i = 1, p
    beta_old(i) = zero
  end do
  
  do i = 1, n
    wt(i) = one
  end do
  
  dev = zero
  
  
  !!! main loop
  main: do iter = 0, maxiter
    
    
    ! compute eta = x*beta and mu = inverse_link( eta )
    if (iter == 0) then
      call glm_initial_mu(family, n, y, wt, mu)
      call glm_link(link, n, mu, eta)
    else 
      call dgemm('n', 'n', n, 1, p, one, x, n, beta, p,  zero, eta, n)
      call glm_linkinv(link, n, eta, mu)
    end if
    
    
    ! check for bad fit in the mu's
    call glm_check_mu(family, n, mu, tol, info)
    if (info /= 0) return
    
    
    ! update wt
    call glm_variance(family, n, mu, wt)
    
    do i = 1, n
      sqwt(i) = dsqrt(wt(i))
    end do
    
    if (stoprule == 3) then
      dev_old = dev
      dev = glm_deviance(family, n, y, mu)
    end if
    
    
    ! prepare lhs:  x_tw = x*wt
    do j = 1, p
      do i = 1, n
        x_tw(i,j) = sqwt(i) * x(i,j)
      end do
    end do
    
    
    ! prepare rhs:  z = sqrt(wt) * (x*beta + 1/wt*(y-mu))
    !                 = sqwt * eta + 1/sqwt*(y-mu)
    do i = 1, n
      tmp = sqwt(i)
      z(i) = tmp*eta(i) + one/(tmp)*(y(i)-mu(i))
    end do
    
    
    ! update beta:  fit z ~ x_tw
    call glm_update_beta(n, p, beta, beta_old, x_tw, z, work, lwork, info)
    
    
    ! check for convergence
    if (iter.gt.0) then
      converged = glm_convergence(stoprule, p, beta_old, beta, dev, dev_old, tol, iter, maxiter, info)
    end if
    
    if (converged == 1) then
      goto 10 ! converged
    else if (converged == 2) then
      goto 1 ! infinite parameter values detected
    end if
    
  end do main
  
  
  !!! success --- now do all the other stuff
10   continue
  
  ! aic, deviance, nulldeviance
  call glm_loglik_stats(family, link, incpt, n, p, x, y, eta, mu, beta, beta_old, dev, aic, nulldev)
  
  ! compute working residuals
  call glm_residuals(link, n, y, mu, eta, resids)
  
  write (*,*) "iter=",iter
  
  goto 1
  
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




! ==================================
! internal functions used in glm_fit
!
! consider these "self-documenting" aka i'm too lazy to explain 
! this mess

subroutine glm_update_beta(n, p, beta, beta_old, x, y, work, lwork, info)
  implicit none
  ! in/out
  integer             n, p, lwork, info
  double precision    beta(p), beta_old(p), x(n,p), y(n), work(lwork)
  ! local
  integer             k, i
  ! external
  external            dgels
  intrinsic           min
  
  
  call dgels('n', n, p, 1, x, n, y, n, work, lwork, info)
  
  k = min(n, p)
  
  do i = 1, k
    beta_old(i) = beta(i)
    beta(i) = y(i)
  end do
  
  return
end










! 1 - converged
! 2 - infinite params
! 3 - no improvement
integer function glm_convergence(stoprule, p, beta_old, beta, dev, dev_old, tol, iter, maxiter, info)
  implicit none
  ! in/out
  integer             stoprule
  integer             p, iter, maxiter, info
  double precision    beta_old(*), beta(*), dev, dev_old, tol
  ! local
  integer             i
  double precision    tmp1, tmp2
  ! intrinsic
  intrinsic           dabs
  logical             disnan
  
  
  glm_convergence = 0
  
  ! check that all parameters are finite and 
  do i = 1, p
    if(disnan(beta(i))) then
      glm_convergence = 2
      info = 1
      goto 1
!        else if (.not.ieee_is_finite(beta(i))) then
!          glm_convergence = 2
!          goto 1
    end if
  end do
  
  if (stoprule == 1) then
    return
  else if (stoprule == 2) then
    do i = 1, p
      tmp1 = dabs(beta(i) - beta_old(i))
      tmp2 = tol*1.0d-6 + tol*dabs(beta_old(i))
      
      if (tmp1.gt.tmp2) then
        glm_convergence = -1
        if (iter == maxiter) info = maxiter
        goto 1
      end if
    end do
    
    glm_convergence = 1
    
1      continue
  else if (stoprule == 3) then
    tmp1 = dabs(dev-dev_old)
    tmp2 = (0.1d0 + dabs(dev))
    
    
    if (tmp1/tmp2.lt.tol) then
      glm_convergence = 1
    else
      glm_convergence = -1
    end if
  end if
  
  return
end



integer function check_response(family, n, y)
  ! in/out
  character*8         family
  integer             n
  double precision    y(*)
  ! local
  integer             i
  ! parameter
  integer             fail
  double precision    zero, one
  parameter ( fail = -8, zero = 0.0d0, one = 1.0d0 )
  
  
  if (family == 'binomial') then
    do i = 1, n
      if (y(i).lt.zero .or. y(i).gt.one) then
        check_response = fail
        return
      end if
    end do
  
  else if (family == 'poisson' .or. family == 'gamma') then
    do i = 1, n
      if (y(i).lt.zero) then
        check_response = fail
        return
      end if
    end do
  
  end if
  
  check_response = 0
  
  return
end





subroutine glm_residuals(link, n, y, mu, eta, resids)
  implicit none
  ! in/out
  character*8         link
  integer             n
  double precision    y(*), mu(*), eta(*), resids(*)
  ! local
  integer             i
  double precision    tmp
  ! parameter
  double precision   one, two, negone
  parameter ( one = 1.0d0, two = 2.0d0, negone = -1.0d0 )
  ! intrinsic
  intrinsic           dexp
  
  
  ! "working" residuals
  
  if (link == 'cloglog') then
    do i = 1, n
      tmp = dexp(eta(i))
      resids(i) = (y(i) - mu(i)) / (tmp * dexp(-tmp))
    end do
  
  else if (link == 'identity') then
    do i = 1, n
      resids(i) = y(i) - mu(i)
    end do
  
  else if (link == 'inverse') then
    do i = 1, n
      tmp = eta(i)
      resids(i) = negone * tmp*tmp * (y(i) - mu(i))
    end do
  
  else if (link == 'log') then
    do i = 1, n
      resids(i) = (y(i) - mu(i)) / mu(i)
    end do
  
  else if (link == 'logit') then
    do i = 1, n
      tmp = mu(i)
      resids(i) = (y(i) - tmp) / tmp / (one - tmp)
    end do
  
  else if (link == 'sqrt') then
    do i = 1, n
      resids(i) = (y(i) - mu(i)) / (two * eta(i))
    end do
  end if
  
  return
end

