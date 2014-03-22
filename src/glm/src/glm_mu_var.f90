! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module glm_mu_var
  use lapack
  
  implicit none
  
  
  contains
  
  subroutine glm_check_mu(family, n, mu, tol, info)
    ! in/out
    character*8, intent(in) :: family
    integer, intent(in) :: n
    integer, intent(out) :: info
    double precision, intent(in) :: mu(*), tol
    ! local
    integer :: i
    double precision :: tmp
    
    
    if (family == 'binomial') then
      tmp = 1.0d0 - tol
      do i = 1, n
        if (mu(i) > tmp .or. mu(i) < tol) then
          info = -101
          return
        end if
      end do
    
    else if (family == 'poisson' .or. family == 'gamma') then
      do i = 1, n
        if (mu(i) < tol) then
          info = -101
          return
        end if
      end do
       
    end if
    
    
    return
  end



  ! initialize mu based on the error distribution
  subroutine glm_initial_mu(family, n, y, wt, mu)
    ! in/out
    character*8, intent(in) :: family
    integer, intent(in) :: n
    double precision, intent(in) :: y(*), wt(*)
    double precision, intent(out) :: mu(*)
    ! local
    integer :: i
    
    
    if (family == 'binomial') then
      do i = 1, n
        mu(i) = (wt(i) * y(i) + 0.5d0) / (wt(i) + 1.0d0)
      end do
    
    else if (family == 'gamma' .or. family == 'gaussian') then
      do i = 1, n
        mu(i) = y(i)
      end do
    
    else if (family == 'poisson') then
      !$omp parallel do private(i) default(shared)
        do i = 1, n
          mu(i) = y(i) + 0.1d0
        end do
      !$omp end parallel do
    end if
    
    return
  end



  subroutine glm_variance(family, n, mu, var)
    ! in/out
    character*8 , intent(in) :: family
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n)
    double precision, intent(out) :: var(n)
    ! local
    integer :: i
    double precision :: tmp
    
    
    if (family == 'binomial') then
        do i = 1, n
          tmp = mu(i)
          var(i) = tmp * (1.0d0 - tmp)
        end do
    
    else if (family == 'gamma') then
      do i = 1, n
        tmp = mu(i)
        var(i) = tmp*tmp
      end do
    
    else if (family == 'gaussian') then
      do i = 1, n
        var(i) = 1.0d0
      end do
    
    else if (family == 'poisson') then
      !$omp parallel do private(i) default(shared)
        do i = 1, n
          var(i) = mu(i)
        end do
      !$omp end parallel do
    end if
    
    return
  end
  
  
end module

