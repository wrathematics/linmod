! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module string_tools
  implicit none
  
  contains
  
  
  ! Makes string all lowercase based on ascii encoding
  function toupper(str) &
  result(ret)
    ! in/out
    character(len=*) :: str
    character(len=len(str)) :: ret
    ! Local
    integer :: i, ascii_i
    integer :: ascii_low = iachar("a")
    integer :: ascii_high = iachar("z")
    
    do i = 1, len(str)
      ascii_i = iachar(str(i:i))
        if (ascii_low <= ascii_i .and. ascii_i <= ascii_high) then
          ret(i:i) = achar(iachar(str(i:i)) - 32)
        else
          ret(i:i) = str(i:i)
        end if
    end do
    
    return
  end function
  
  
  
  ! Makes string all uppercase based on ascii encoding
  function tolower(str) &
  result(ret)
    ! in/out
    character(len=*) :: str
    character(len=len(str)) :: ret
    ! Local
    integer :: i, ascii_i
    integer :: ascii_low = iachar("A")
    integer :: ascii_high = iachar("Z")
    
    do i = 1, len(str)
      ascii_i = iachar(str(i:i))
        if (ascii_low <= ascii_i .and. ascii_i <= ascii_high) then
          ret(i:i) = achar(iachar(str(i:i)) + 32)
        else
          ret(i:i) = str(i:i)
        end if
    end do
    
    return
  end function
  
  
  
  ! Test equality ignoring case
  function equivchar(a, b) &
  result(test)
    ! in/out
    logical :: test
    character(len=1), intent(in) :: a, b
    ! local
    character(len=1) :: a_l, b_l
    
    
    if (a == b) then
      test = .true.
    else
      a_l = tolower(a)
      b_l = tolower(b)
      
      if (a_l == b_l) then
        test = .true.
      else
        test = .false.
      end if
    end if
    
    return
  end function
  
end module


