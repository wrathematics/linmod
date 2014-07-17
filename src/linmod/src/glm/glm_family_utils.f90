! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt

module glm_family_utils
  use :: lapack
  use :: glm_constants
  
  implicit none
  
  
  contains
  
  subroutine glm_check_fitted(family, n, mu, info)
    ! in/out
    integer, intent(in) :: family
    integer, intent(in) :: n
    integer, intent(out) :: info
    double precision, intent(in) :: mu(*)
    ! local
    integer :: i
    
    
    if (family == glm_family_binomial) then
      do i = 1, n
        if (mu(i) >= 1.0d0 .or. mu(i) <= 0.0d0) then
          info = glm_family_badmu
          return
        end if
      end do
    
    else if (family == glm_family_poisson .or. &
             family == glm_family_gamma) then
      do i = 1, n
        if (mu(i) <= 0.0d0) then
          info = glm_family_badmu
          return
        end if
      end do
       
    end if
    
    
    return
  end



  ! initialize mu based on the error distribution
  subroutine glm_initial_mu(family, n, y, wt, mu)
    ! in/out
    integer, intent(in) :: family
    integer, intent(in) :: n
    double precision, intent(in) :: y(*), wt(*)
    double precision, intent(out) :: mu(*)
    ! local
    integer :: i
    
    
    !$omp parallel if(n > 5000) private(i) default(shared) 
    if (family == glm_family_binomial) then
      !$omp do
        do i = 1, n
          mu(i) = (wt(i) * y(i) + 0.5d0) / (wt(i) + 1.0d0)
        end do
      !$omp end do
    
    else if (family == glm_family_gamma            .or. &
             family == glm_family_gaussian         .or. &
             family == glm_family_inversegaussian) then
      !$omp do
        do i = 1, n
          mu(i) = y(i)
        end do
      !$omp end do
    
    else if (family == glm_family_poisson) then
      !$omp do
        do i = 1, n
          mu(i) = y(i) + 0.1d0
        end do
      !$omp end do
    
    end if
    !$omp end parallel
    
    return
  end



  subroutine glm_variance(family, n, mu, var)
    ! in/out
    integer, intent(in) :: family
    integer, intent(in) :: n
    double precision, intent(in) :: mu(n)
    double precision, intent(out) :: var(n)
    ! local
    integer :: i
    double precision :: tmp
    
    
    !$omp parallel if(n > 5000) private(i, tmp) default(shared) 
    if (family == glm_family_binomial) then
      !$omp do
        do i = 1, n
          tmp = mu(i)
          var(i) = tmp * (1.0d0 - tmp)
        end do
      !$omp end do
    
    else if (family == glm_family_gamma) then
      !$omp do
        do i = 1, n
          tmp = mu(i)
          var(i) = tmp*tmp
        end do
      !$omp end do
    
    else if (family == glm_family_gaussian) then
      !$omp do
        do i = 1, n
          var(i) = 1.0d0
        end do
      !$omp end do
    
    else if (family == glm_family_poisson) then
      !$omp do
        do i = 1, n
          var(i) = mu(i)
        end do
      !$omp end do
    
    else if (family == glm_family_inversegaussian) then
      !$omp do
        do i = 1, n
          tmp = mu(i)
          var(i) = tmp*tmp*tmp
        end do
      !$omp end do
    end if
    !$omp end parallel
    
    return
  end
  
  
end module

