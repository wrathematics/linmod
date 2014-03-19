! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module glm_mu_var
  implicit none
  
  
  contains

  subroutine glm_check_mu(family, n, mu, tol, info)
    implicit none
    ! in/out
    character*8         family
    integer             n, info
    double precision    mu(*), tol
    ! local
    integer             i
    double precision    tmp
    ! parameters
    double precision    one, zero
    parameter ( one = 1.0d0, zero = 0.0d0 )
    
    
    if (family == 'binomial') then
      tmp = one-tol
      do i = 1, n
        if (mu(i).gt.tmp .or. mu(i).lt.tol) then
          info = -101
          return
        end if
      end do
    
    else if (family == 'poisson' .or. family == 'gamma') then
      do i = 1, n
        if (mu(i).lt.tol) then
          info = -101
          return
        end if
      end do
    end if
    
    
    return
  end



  ! initialize mu based on the error distribution
  subroutine glm_initial_mu(family, n, y, wt, mu)
    implicit none
    ! in/out
    character*8         family
    integer             n, p, lwork, info
    double precision    y(*), wt(*), mu(*)
    ! local
    integer             i
    ! parameter
    double precision    one, tenth, half
    parameter ( one = 1.0d0, tenth = 0.1d0, half = 0.5d0 )
    ! external
    external            dgels
    
    
    if (family == 'binomial') then
      do i = 1, n
        mu(i) = (wt(i) * y(i) + half) / (wt(i) + one)
      end do
    
    else if (family == 'gamma' .or. family == 'gaussian') then
      do i = 1, n
        mu(i) = y(i)
      end do
    
    else if (family == 'poisson') then
      do i = 1, n
        mu(i) = y(i) + tenth
      end do
    end if
    
    return
  end



  subroutine glm_variance(family, n, mu, var)
    implicit none
    ! in/out
    character*8         family
    integer             n
    double precision    mu(n), var(n)
    ! local
    integer             i
    double precision    tmp
    ! parameter
    double precision    one
    parameter ( one = 1.0d0 )
    
    
    if (family == 'binomial') then
        do i = 1, n
          tmp = mu(i)
          var(i) = tmp * (one - tmp)
        end do
    
    else if (family == 'gamma') then
      do i = 1, n
        tmp = mu(i)
        var(i) = tmp*tmp
      end do
    
    else if (family == 'gaussian') then
      do i = 1, n
        var(i) = one
      end do
    
    else if (family == 'poisson') then
      do i = 1, n
        var(i) = mu(i)
      end do
    end if
    
    return
  end
  
  
end module

