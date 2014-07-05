! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


if (xlen <= 1) return

do i = 2, xlen
  tmp = x(i)
  j = i
  do
    if (j > 1 .and. x(j-1) > tmp) then
      x(j) = x(j-1)
      j = j - 1
    else
      exit
    end if
  end do
    
    x(j) = tmp
end do

