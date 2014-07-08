! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt

integer :: i, j
integer :: indtmp

if (length <= 1) return

do i = 2, length
  indtmp = ind(i)
  tmp = arr(i)
  j = i
  do
    if (j > 1 .and. ind(j-1) > indtmp) then
      ind(j) = ind(j-1)
      arr(j) = arr(j-1)
      j = j - 1
    else
      exit
    end if
  end do
    
    ind(j) = indtmp
    arr(j) = tmp
end do
