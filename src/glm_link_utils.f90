! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module glm_link_utils
  implicit none
  
  
  contains

  ! link function
  subroutine glm_link(link, n, x, y)
    ! in/out
    character*8         link
    integer             n
    double precision    x(*), y(*)
    ! local
    integer             i
    double precision    tmp
    ! parameter
    double precision    one
    parameter ( one = 1.0d0 )
    ! intrinsic
    intrinsic           dlog, dsqrt
    
    
    if (link == 'cloglog') then
      do i = 1, n
        y(i) = dlog(-dlog(one-x(i)))
      end do
    
    else if (link == 'identity') then
      do i = 1, n
        y(i) = x(i)
      end do
    
    else if (link == 'inverse') then
      do i = 1, n
        y(i) = one/x(i)
      end do
    
    else if (link == 'log') then
      do i = 1, n
        y(i) = dlog(x(i))
      end do
    
    else if (link == 'logit') then
      do i = 1, n
        tmp = x(i)
        y(i) = dlog(tmp / (one-tmp))
      end do
    
    else if (link == 'sqrt') then
      y(i) = dsqrt(x(i))
    
    end if
    
    return
  end



  ! inverse link function
  subroutine glm_linkinv(link, n, x, y)
    ! in/out
    character*8         link
    integer             n
    double precision    x(*), y(*)
    ! local
    integer             i
    double precision    tmp
    ! parameter
    double precision    zero, one
    parameter ( zero = 0.0d0, one = 1.0d0 )
    ! intrinsic
    intrinsic           dexp, dsqrt
    
    
    if (link == 'cloglog') then
      do i = 1, n
        tmp = -dexp(x(i))
        y(i) = -dexp(tmp)-one
      end do
    
    else if (link == 'identity') then
      do i = 1, n
        y(i) = x(i)
      end do
    
    else if (link == 'inverse') then
      do i = 1, n
        if (x(i).gt.zero) then
          y(i) = one/x(i)
        else
          y(i) = zero
        end if
      end do
      
    else if (link == 'log') then
        do i = 1, n
          y(i) = dexp(x(i))
        end do
    
    else if (link == 'logit') then
      do i = 1, n
        tmp = dexp(x(i))
        y(i) = tmp / (one + tmp)
      end do
    
    else if (link == 'sqrt') then
      do i = 1, n
        tmp = x(i)
        y(i) = tmp*tmp
      end do
    
    end if
    
    return
  end



  ! check the family and link arguments for valid/supported possibilities
  ! 0: no problem
  ! -1 family is invalid or unsupported
  ! -2 link is invalid or unsupported

  function glm_check_fam_link(family, link) &
  result(check)
    integer :: check
    ! in/out
    character*8         family, link
    ! parameters
    integer             bad_fam, bad_link
    parameter ( bad_fam = -1, bad_link = -2 )
    
    
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
  
  
end module
