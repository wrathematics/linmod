! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module glm_update_utils
  use :: lapack, only : dgels, disnan
  use :: glm_constants
  use :: glm_family_utils
  use :: glm_link_utils
!  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  
  implicit none
  
  
  contains
  
  ! Linear model iteration
  subroutine glm_update_beta(n, p, beta, beta_old, x, y, work, lwork, info)
    ! in/out
    integer, intent(in) :: n, p, lwork
    integer, intent(out) :: info
    double precision, intent(in) :: x(n,p)
    double precision, intent(inout) :: y(n)
    double precision, intent(out) :: beta(p), beta_old(p), work(lwork)
    ! local
    integer :: k, i
    ! external
    intrinsic :: min
    
    
    call dgels('n', n, p, 1, x, n, y, n, work, lwork, info)
    
    k = min(n, p)
    
    !$omp parallel do private(i) default(shared)
      do i = 1, k
        beta_old(i) = beta(i)
        beta(i) = y(i)
      end do
    !$omp end parallel do
    
    return
  end
  
  
  
  
  
  ! rtwt = sqrt( linkinv(eta) / glm_variance(mu) )
  subroutine glm_update_wt(family, link, n, p, x, x_tw, wt, rtwt, y, mu, eta, z)
    ! in/out
    integer, intent(in) :: family, link, n, p
    double precision, intent(in) :: y(n), mu(n), eta(n)
    double precision, intent(out) :: wt(n), rtwt(*), z(n)
    double precision, intent(in) :: x(:,:)
    double precision, intent(out) :: x_tw(:,:)
    ! local
    integer :: i, j
    double precision :: tmp
    
    
    call glm_variance(family, n, mu, wt)
    
    call glm_linkinv_deriv(link, n, eta, rtwt)
    
    !$omp parallel if (n*p > 5000) private(i, j) default(shared) 
    !$omp do
      do i = 1, n
        rtwt(i) = rtwt(i) / dsqrt(wt(i))
      end do
    !$omp end do
    
    
    ! prepare lhs:  x_tw = x*wt
    !$omp do
      do j = 1, p
        do i = 1, n
          x_tw(i, j) = rtwt(i) * x(i, j)
        end do
      end do
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine
  
  
  
  ! Set the working response z
  subroutine glm_update_z(family, link, n, p, x, x_tw, wt, rtwt, y, mu, eta, z)
    ! in/out
    integer, intent(in) :: family, link, n, p
    double precision, intent(in) :: y(n), mu(n), eta(n)
    double precision, intent(out) :: wt(n), rtwt(*), z(n)
    double precision, intent(in) :: x(:,:)
    double precision, intent(out) :: x_tw(:,:)
    ! local
    integer :: i, j
    double precision :: tmp
    
    
    
    
    ! prepare rhs:  z = sqrt(wt) * (x*beta + 1/wt*(y-mu))
    !                 = rtwt * eta + 1/rtwt*(y-mu)
    call glm_linkinv_deriv(link, n, eta, z)
    
    !$omp parallel if (n > 5000) private(i, tmp) default(shared) 
    !$omp do private(i, tmp)
      do i = 1, n
        tmp = rtwt(i)
        z(i) = tmp*eta(i) + tmp*(y(i)-mu(i)) / z(i)
      end do
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine
  
  
  
  ! 1 - converged
  ! 2 - infinite params
  ! 3 - no improvement
  function glm_check_convergence(stoprule, p, beta_old, beta, dev, dev_old, tol, iter, maxiter) &
  result(converged)
    ! in/out
    integer :: converged
    integer, intent(in) :: stoprule
    integer, intent(in) :: p, iter, maxiter
    double precision, intent(in) :: beta_old(*), beta(*), dev, dev_old, tol
    ! local
    integer :: i
    double precision :: tmp1, tmp2
    double precision :: eps = 1.0d-6
    ! intrinsic
    intrinsic :: dabs
    
    
    ! Check parameters
    do i = 1, p
      if(disnan(beta(i))) then
        converged = glm_convergence_infparams
        return
!      else if (.not. ieee_is_finite(beta(i))) then
!        converged = glm_convergence_infparams
        return
      end if
    end do
    
    ! Check stoprule
    if (stoprule == glm_stoprule_maxiter) then
      return
    else if (stoprule == glm_stoprule_coefs) then
      do i = 1, p
        tmp1 = dabs(beta(i) - beta_old(i))
        tmp2 = tol*eps + tol*dabs(beta_old(i))
        
        if (tmp1 > tmp2) then
          converged = glm_convergence_noconvergence
          return
        end if
      end do
      
      converged = glm_convergence_converged
      
    else if (stoprule == glm_stoprule_deviance) then
      tmp1 = dabs(dev-dev_old)
      tmp2 = (0.1d0 + dabs(dev))
      
      if (tmp1/tmp2 < tol) then
        converged = glm_convergence_converged
      else
        converged = glm_convergence_noconvergence
      end if
    end if
    
    return
  end
  
end module
