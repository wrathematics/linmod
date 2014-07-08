! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module glm_loglik_utils
  use :: glm_link_utils
  use :: glm_constants
  
  implicit none
  
  
  contains
  
  ! model deviance calculator
  function glm_deviance(family, n, y, mu) &
  result(dev)
    ! in/out
    double precision :: dev
    integer, intent(in) :: family
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), y(n)
    ! local
    integer :: i
    double precision :: tmp
    
    
    dev = 0.0d0
    
    if (family == glm_family_gaussian) then
      !$omp parallel do private(i, tmp) default(shared) reduction(+:dev)
        do i = 1, n
          tmp = y(i) - mu(i)
          dev = dev + tmp*tmp
        end do
      !$omp end parallel do
      
    
    else if (family == glm_family_poisson) then
      !$omp parallel do private(i) default(shared) reduction(+:dev)
        do i = 1, n
          if (y(i) > 0.0d0) then
            dev = dev + y(i)*dlog(y(i)/mu(i)) 
          end if
          dev = dev + mu(i) - y(i)
        end do
      !$omp end parallel do
      
      dev = 2.0d0 * dev
      
      
    else if (family == glm_family_binomial) then
      !$omp parallel do private(i) default(shared) reduction(+:dev)
        do i = 1, n
          if (y(i) > 0.0d0) then
            dev = dev + dlog(mu(i))
          else
            dev = dev + dlog(1.0d0 - mu(i))
          end if
        end do
      !$omp end parallel do
      
      dev = -2.0d0 * dev
    
    
    else if (family == glm_family_gamma) then
      !$omp parallel do private(i) default(shared) reduction(+:dev)
        do i = 1, n
          if (y(i) > 0.0d0) then
            dev = dev - dlog(y(i)/mu(i))
          end if
          dev = dev + (y(i) - mu(i))/mu(i)
        end do
      !$omp end parallel do
      
      dev = -2.0d0 * dev
    
    end if
    
    return
  end
  
  
  
  ! null model deviance
  subroutine glm_nulldev(family, n, y, mu, dev)
    ! in/out
    integer, intent(in) :: family
    integer, intent(in) :: n
    double precision, intent(in) :: y(*)
    double precision, intent(out) :: mu(*), dev
    ! local
    integer :: i
    double precision :: tmp
    
    
    tmp = 0.0d0
    !$omp parallel private(i) default(shared)
      !$omp do reduction(+:tmp)
        do i = 1, n
          tmp = tmp + y(i)/n
        end do
      
      !$omp do
        do i = 1, n
          mu(i) = tmp
        end do
      !$omp end do
    !$omp end parallel
    
    dev = glm_deviance(family, n, y, mu)
      
    
    return
  end
  
  
  
  ! loglikelihood calculator
  function glm_loglik(family, n, y, mu) &
  result(llik)
    ! in/out
    double precision :: llik
    integer, intent(in) :: family
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n), y(n)
    ! local
    integer :: i
    double precision :: tmp
    
    
    llik = 0.0d0
    
    if (family == glm_family_gaussian) then
      !$omp parallel do private(i, tmp) default(shared) reduction(+:llik)
        do i = 1, n
          tmp = y(i) - mu(i)
          llik = llik - tmp*tmp
        end do
      !$omp end parallel do
      
      llik = llik / 2.0d0
    
    
    else if (family == glm_family_poisson) then
      !$omp parallel do private(i) default(shared) reduction(+:llik)
        do i = 1, n
          if (y(i) > 0.0d0) then
            llik = llik - y(i)*dlog(y(i)/mu(i)) 
          end if
          llik = llik + y(i) - mu(i)
        end do
      !$omp end parallel do
      
      
    else if (family == glm_family_binomial) then
      !$omp parallel do private(i) default(shared) reduction(+:llik)
        do i = 1, n
          if (y(i) > 0.0d0) then
            llik = llik + dlog(mu(i))
          else
            llik = llik + dlog(1.0d0 - mu(i))
          end if
        end do
      !$omp end parallel do
    
    
    else if (family == glm_family_gamma) then
      !$omp parallel do private(i) default(shared) reduction(+:llik)
        do i = 1, n
          if (y(i) > 0.0d0) then
            llik = llik + dlog(y(i)/mu(i))
          end if
          llik = llik - (y(i) - mu(i))/mu(i)
        end do
      !$omp end parallel do
      
    end if
    
    return
  end
  
  
  
  subroutine glm_loglik_stats(family, link, intercept, n, p, x, y, eta, mu, beta, beta_old, dev, aic, nulldev)
    ! in/out
    integer, intent(in) :: intercept, family, link
    integer, intent(in) :: n, p
    double precision, intent(in) :: x(*), y(n), beta(*), beta_old(p)
    double precision, intent(out) :: dev, aic, nulldev
    double precision, intent(out) :: eta(n), mu(*)
    ! local
    integer :: i
    double precision :: tmp
    ! intrinsic
    intrinsic :: dble, sum
    
    
    ! model deviance
  !      dev = -2.0d0 * glm_loglik(family, n, y, mu)
    dev = glm_deviance(family, n, y, mu)
    
    ! model aic
    aic = 2.0d0 * dble(p) + dev
    
    !!! null deviance
    ! null deviance for model with no intercept
    if (intercept == glm_intercept_no) then
      do i = 1, n
        eta(i) = 0.0d0
      end do
      
      call glm_linkinv(link, n, eta, mu)
      nulldev = -2.0d0 * glm_loglik(family, n, y, mu)
    
    ! null deviance (deviance for intercept-only model)
    else
      call glm_nulldev(family, n, y, mu, nulldev)
      
    end if
    
    
    return
  end
  
  
end module
