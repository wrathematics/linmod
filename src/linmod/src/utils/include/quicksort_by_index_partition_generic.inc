! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


integer :: i
integer :: pvt

ret = l

pvt = ind(pvtind)

call swap_ind(ind, pvtind, r)
call swap_ind(arr, pvtind, r)

do i = l, r
  if (ind(i) < pvt) then
    call swap_ind(ind, i, ret)
    call swap_ind(arr, i, ret)
    ret = ret + 1
  end if
end do

call swap_ind(ind, ret, r)
call swap_ind(arr, ret, r)

ret = ret - 1 ! for "<=" return
