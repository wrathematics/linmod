! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module glm_residuals_utils
  implicit none
  
  
  contains
  
  subroutine glm_residuals(link, n, y, mu, eta, resids)
    ! in/out
    character*8, intent(in) :: link
    integer, intent(in) :: n
    double precision, intent(in) :: y(*), mu(*), eta(*)
    double precision, intent(out) :: resids(*)
    ! local
    integer :: i
    double precision :: tmp
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
        resids(i) = -1.0d0 * tmp*tmp * (y(i) - mu(i))
      end do
    
    else if (link == 'log') then
      do i = 1, n
        resids(i) = (y(i) - mu(i)) / mu(i)
      end do
    
    else if (link == 'logit') then
      do i = 1, n
        tmp = mu(i)
        resids(i) = (y(i) - tmp) / tmp / (1.0d0 - tmp)
      end do
    
    else if (link == 'sqrt') then
      do i = 1, n
        resids(i) = (y(i) - mu(i)) / (2.0d0 * eta(i))
      end do
    end if
    
    return
  end
  
  
end module
