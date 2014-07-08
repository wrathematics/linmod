! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module glm_link_utils
  use :: glm_constants
  implicit none
  
  contains
  
  ! link function
  subroutine glm_link(link, n, x, y)
    ! in/out
    integer, intent(in) :: link
    integer, intent(in) :: n
    double precision, intent(in) :: x(*)
    double precision, intent(out) :: y(*)
    ! local
    integer :: i
    double precision :: tmp
    ! intrinsic
    intrinsic :: dlog, dsqrt
    
    
    if (link == glm_link_cloglog) then
      !$omp parallel do private(i) default(shared)
        do i = 1, n
          y(i) = dlog(-dlog(1.0d0-x(i)))
        end do
      !$omp end parallel do
    
    else if (link == glm_link_identity) then
      !$omp parallel do private(i) default(shared)
        do i = 1, n
          y(i) = x(i)
        end do
      !$omp end parallel do
    
    else if (link == glm_link_inverse) then
      !$omp parallel do private(i) default(shared)
        do i = 1, n
          y(i) = 1.0d0/x(i)
        end do
      !$omp end parallel do
    
    else if (link == glm_link_log) then
      !$omp parallel do private(i) default(shared)
        do i = 1, n
          y(i) = dlog(x(i))
        end do
      !$omp end parallel do
    
    else if (link == glm_link_logit) then
      !$omp parallel do private(i, tmp) default(shared)
        do i = 1, n
          tmp = x(i)
          y(i) = dlog(tmp / (1.0d0 - tmp))
        end do
      !$omp end parallel do
    
    else if (link == glm_link_sqrt) then
      !$omp parallel do private(i) default(shared)
        do i = 1, n
          y(i) = dsqrt(x(i))
        end do
      !$omp end parallel do
    
    end if
    
    return
  end
  
  
  
  ! inverse link function
  subroutine glm_linkinv(link, n, x, y)
    ! in/out
    integer, intent(in) :: link
    integer, intent(in) :: n
    double precision, intent(in) :: x(*)
    double precision, intent(out) :: y(*)
    ! local
    integer :: i
    double precision :: tmp
    ! intrinsic
    intrinsic           dexp, dsqrt
    
    
    if (link == glm_link_cloglog) then
      !$omp parallel do private(i, tmp) default(shared)
        do i = 1, n
          tmp = -dexp(x(i))
          y(i) = -dexp(tmp) - 1.0d0
        end do
      !$omp end parallel do
    
    else if (link == glm_link_identity) then
      !$omp parallel do private(i) default(shared)
        do i = 1, n
          y(i) = x(i)
        end do
      !$omp end parallel do
    
    else if (link == glm_link_inverse) then
      !$omp parallel do private(i) default(shared)
        do i = 1, n
          if (x(i) > 0.0d0) then
            y(i) = 1.0d0/x(i)
          else
            y(i) = 0.0d0
          end if
        end do
      !$omp end parallel do
      
    else if (link == glm_link_log) then
      !$omp parallel do private(i) default(shared)
        do i = 1, n
          y(i) = dexp(x(i))
        end do
      !$omp end parallel do
    
    else if (link == glm_link_logit) then
      !$omp parallel do private(i, tmp) default(shared)
        do i = 1, n
          tmp = dexp(x(i))
          y(i) = tmp / (1.0d0 + tmp)
        end do
      !$omp end parallel do
    
    else if (link == glm_link_sqrt) then
      !$omp parallel do private(i, tmp) default(shared)
        do i = 1, n
          tmp = x(i)
          y(i) = tmp*tmp
        end do
      !$omp end parallel do
    
    end if
    
    return
  end
  
  
  
  subroutine glm_residuals(link, n, y, mu, eta, resids)
    ! in/out
    integer, intent(in) :: link
    integer, intent(in) :: n
    double precision, intent(in) :: y(*), mu(*), eta(*)
    double precision, intent(out) :: resids(*)
    ! local
    integer :: i
    double precision :: tmp
    ! intrinsic
    intrinsic           dexp
    
    
    ! "working" residuals
    
    if (link == glm_link_cloglog) then
      !$omp parallel do private(i, tmp) default(shared)
        do i = 1, n
          tmp = dexp(eta(i))
          resids(i) = (y(i) - mu(i)) / (tmp * dexp(-tmp))
        end do
      !$omp end parallel do
    
    else if (link == glm_link_identity) then
      !$omp parallel do private(i) default(shared)
        do i = 1, n
          resids(i) = y(i) - mu(i)
        end do
      !$omp end parallel do
    
    else if (link == glm_link_inverse) then
      !$omp parallel do private(i, tmp) default(shared)
        do i = 1, n
          tmp = eta(i)
          resids(i) = -1.0d0 * tmp*tmp * (y(i) - mu(i))
        end do
      !$omp end parallel do
    
    else if (link == glm_link_log) then
      !$omp parallel do private(i) default(shared)
        do i = 1, n
          resids(i) = (y(i) - mu(i)) / mu(i)
        end do
      !$omp end parallel do
    
    else if (link == glm_link_logit) then
      !$omp parallel do private(i, tmp) default(shared)
        do i = 1, n
          tmp = mu(i)
          resids(i) = (y(i) - tmp) / tmp / (1.0d0 - tmp)
        end do
      !$omp end parallel do
    
    else if (link == glm_link_sqrt) then
      !$omp parallel do private(i) default(shared)
        do i = 1, n
          resids(i) = (y(i) - mu(i)) / (2.0d0 * eta(i))
        end do
      !$omp end parallel do
    end if
    
    return
  end
  
  
end module
