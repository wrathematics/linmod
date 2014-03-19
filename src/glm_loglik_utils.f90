! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module glm_loglik_utils
  implicit none
  
  contains
  
  ! model deviance calculator
  function glm_deviance(family, n, y, mu) &
  result(dev)
  double precision :: dev
    ! in/out
    character*8         family
    integer             n
    double precision    mu(n), y(n)
    ! local
    integer             i
    double precision    tmp
    ! parameter
    double precision    zero, one, two
    parameter ( zero = 0.0d0, one = 1.0d0, two = 2.0d0 )
    
    
    dev = zero
    
    !!! normal
    if (family == 'gaussian') then
      do i = 1, n
        tmp = y(i) - mu(i)
        dev = dev + tmp*tmp
      end do
      
    !!! poisson
    else if (family == 'poisson') then
      do i = 1, n
        if (y(i).gt.zero) then
          dev = dev + y(i)*dlog(y(i)/mu(i)) 
        end if
        dev = dev + mu(i) - y(i)
      end do
      
      dev = two*dev
      
    !!! binomial
    else if (family == 'binomial') then
      do i = 1, n
        if (y(i).gt.zero) then
          dev = dev + dlog(mu(i))
        else
          dev = dev + dlog(one-mu(i))
        end if
      end do
      
      dev = -two*dev
    
    !!! gamma
    else if (family == 'gamma') then
      do i = 1, n
        if (y(i).gt.zero) then
          dev = dev - dlog(y(i)/mu(i))
        end if
        dev = dev + (y(i) - mu(i))/mu(i)
      end do
      
      dev = -two*dev
    
    end if
    
    return
  end



  ! loglikelihood calculator
  function glm_loglik(family, n, y, mu) &
  result(llik)
  double precision llik
    ! in/out
    character*8         family
    integer             n
    double precision    mu(n), y(n)
    ! local
    integer             i
    double precision    tmp
    ! parameter
    double precision    zero, one, two
    parameter ( zero = 0.0d0, one = 1.0d0, two = 2.0d0 )
    
    
    llik = zero
    
    !!! gaussian
    if (family == 'gaussian') then
      do i = 1, n
        tmp = y(i) - mu(i)
        llik = llik - tmp*tmp
      end do
      
      llik = llik / two
    
    !!! poisson
    else if (family == 'poisson') then
      do i = 1, n
        if (y(i).gt.zero) then
          llik = llik - y(i)*dlog(y(i)/mu(i)) 
        end if
        llik = llik + y(i) - mu(i)
      end do
      
    !!! binomial
    else if (family == 'binomial') then
      do i = 1, n
        if (y(i).gt.zero) then
          llik = llik + dlog(mu(i))
        else
          llik = llik + dlog(one-mu(i))
        end if
      end do
      
    else if (family == 'gamma') then
      do i = 1, n
        if (y(i).gt.zero) then
          llik = llik + dlog(y(i)/mu(i))
        end if
        llik = llik - (y(i) - mu(i))/mu(i)
      end do
      
    end if
    
    return
  end



  function glm_nulldev(family, n, y, mu) &
  result(dev)
  double precision dev
    ! in/out
    character*8         family
    integer             n
    double precision    mu(*), y(*)
    ! local
    integer             i
    double precision    tmp
    ! parameter
    double precision    zero, one, two
    parameter ( zero = 0.0d0, one = 1.0d0, two = 2.0d0 )
    ! functions
    double precision    glm_deviance
    
    
    tmp = zero
    do i = 1, n
      tmp = tmp + y(i)/n
    end do
    
    do i = 1, n
      mu(i) = tmp
    end do
    
    dev = glm_deviance(family, n, y, mu)
      
    
    return
  end



  subroutine glm_loglik_stats(family, link, incpt, n, p, x, y, eta, mu, beta, beta_old, dev, aic, nulldev)
    implicit none
    ! in/out
    character*1         incpt
    character*8         family, link
    integer             n, p
    double precision    x(*), y(n), eta(n), mu(*), beta(*), beta_old(p), dev, aic, nulldev
    ! local
    integer             i
    double precision    tmp
    ! parameter
    double precision    zero, one, two, negtwo
    parameter ( zero = 0.0d0, one = 1.0d0, two = 2.0d0, negtwo = -2.0d0 )
    ! intrinsic
    double precision    glm_loglik, glm_deviance, glm_nulldev
    intrinsic           dble, sum
    
    
    ! model deviance
  !      dev = negtwo * glm_loglik(family, n, y, mu)
    dev = glm_deviance(family, n, y, mu)
    
    ! model aic
    aic = two * dble(p) + dev
    
    !!! null deviance
    ! null deviance for model with no intercept
    if (incpt == 'n') then
      do i = 1, n
        eta(i) = zero
      end do
      
      call glm_linkinv(link, n, eta, mu)
      nulldev = negtwo * glm_loglik(family, n, y, mu)
    
    ! null deviance (deviance for intercept-only model)
    else
      nulldev = glm_nulldev(family, n, y, mu)
      
    end if
    
    
    return
  end
  
  
end module
