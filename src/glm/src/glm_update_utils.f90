! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module glm_update_utils
  use lapack
  
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
    
    do i = 1, k
      beta_old(i) = beta(i)
      beta(i) = y(i)
    end do
    
    return
  end
  
  
  
  ! 1 - converged
  ! 2 - infinite params
  ! 3 - no improvement
  function glm_convergence(stoprule, p, beta_old, beta, dev, dev_old, tol, iter, maxiter) &
  result(converged)
    ! in/out
    integer :: converged
    integer, intent(in) :: stoprule
    integer, intent(in) :: p, iter, maxiter
    double precision, intent(in) :: beta_old(*), beta(*), dev, dev_old, tol
    ! local
    integer :: i
    double precision :: tmp1, tmp2
    ! intrinsic
    intrinsic :: dabs
    
    
    converged = 0
    
    ! check that all parameters are finite and 
    do i = 1, p
      if(disnan(beta(i))) then
        converged = 2
!        info = 1
        goto 1
  !        else if (.not.ieee_is_finite(beta(i))) then
  !          converged = 2
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
          converged = -1
!          if (iter == maxiter) info = maxiter
          goto 1
        end if
      end do
      
      converged = 1
      
  1      continue
    else if (stoprule == 3) then
      tmp1 = dabs(dev-dev_old)
      tmp2 = (0.1d0 + dabs(dev))
      
      
      if (tmp1/tmp2.lt.tol) then
        converged = 1
      else
        converged = -1
      end if
    end if
    
    return
  end
  
  
  
  function check_response(family, n, y) &
  result(check)
    ! in/out
    integer :: check
    character*8, intent(in) :: family
    integer, intent(in) :: n
    double precision, intent(in) :: y(*)
    ! local
    integer :: i
    ! parameter
    integer, parameter :: fail = -8
    
    
    if (family == 'binomial') then
      do i = 1, n
        if (y(i) < 0.0d0 .or. y(i) > 1.0d0) then
          check = fail
          return
        end if
      end do
    
    else if (family == 'poisson' .or. family == 'gamma') then
      do i = 1, n
        if (y(i) < 0.0d0) then
          check = fail
          return
        end if
      end do
    
    end if
    
    check = 0
    
    return
  end
  
  
end module
