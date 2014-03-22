! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module glm_check
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
    character*8, intent(in) :: family, link
    ! parameters
    integer, parameter :: bad_fam = -1, bad_link = -2
    
    
    check = 0
    
    if (family == 'binomial') then
        if (link /= 'cloglog'  .and. &
            link /= 'log'      .and. &
            link /= 'logit')    then
                check = bad_link
        end if
    
    else if (family == 'gamma') then
        if (link /= 'identity'   .and. &
            link /= 'log'        .and. &
            link /= 'inverse')   then
                check = bad_link
        end if
    
    else if (family == 'gaussian') then
        if (link /= 'identity'      .and. &
            link /= 'log'           .and. &
            link /= 'inverse')       then
                check = bad_link
        end if
    
    else if (family == 'poisson') then
        if (link /= 'identity'     .and. &
            link /= 'log'          .and. &
            link /= 'sqrt')         then
                check = bad_link
        end if
    
    else
        ! family not supported
        check = bad_fam
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
    
    
    check = 0
    
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
    
    return
  end
  
  
end module
