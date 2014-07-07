! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module glm_check
  use :: glm_family_utils
  use :: glm_link_utils
  implicit none
  
  contains
  
  ! check the family and link arguments for valid/supported possibilities
  ! 0: no problem
  ! -1 family is invalid or unsupported
  ! -2 link is invalid or unsupported
  function glm_check_fam_link(family, link) &
  result(check)
    ! in/out
    integer :: check
    integer, intent(in) :: family, link
    
    
    check = 0
    
    if (family == glm_family_binomial) then
        if (link /= glm_link_cloglog  .and. &
            link /= glm_link_log      .and. &
            link /= glm_link_logit)    then
                check = glm_link_unsupported
        end if
    
    else if (family == glm_family_gamma) then
        if (link /= glm_link_identity   .and. &
            link /= glm_link_log        .and. &
            link /= glm_link_inverse)   then
                check = glm_link_unsupported
        end if
    
    else if (family == glm_family_gaussian) then
        if (link /= glm_link_identity      .and. &
            link /= glm_link_log           .and. &
            link /= glm_link_inverse)       then
                check = glm_link_unsupported
        end if
    
    else if (family == glm_family_poisson) then
        if (link /= glm_link_identity     .and. &
            link /= glm_link_log          .and. &
            link /= glm_link_sqrt)         then
                check = glm_link_unsupported
        end if
    
    else
        ! family not supported
        check = glm_family_unsupported
    end if
    
    return
  end
  
  
  
  function check_response(family, n, y) &
  result(check)
    ! in/out
    integer :: check
    integer, intent(in) :: family
    integer, intent(in) :: n
    double precision, intent(in) :: y(*)
    ! local
    integer :: i
    ! parameter
    integer, parameter :: fail = -8
    
    
    check = 0
    
    if (family == glm_family_binomial) then
      do i = 1, n
        if (y(i) < 0.0d0 .or. y(i) > 1.0d0) then
          check = fail
          return
        end if
      end do
    
    else if (family == glm_family_poisson .or. family == glm_family_gamma) then
      do i = 1, n
        if (y(i) < 0.0d0) then
          check = fail
          return
        end if
      end do
    
    end if
    
    return
  end
  
  
end module
