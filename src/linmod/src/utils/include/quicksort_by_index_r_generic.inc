! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


! choose median of l, r, and (l+r)/2 as pivot
tmp = (/ l, (l+r)/2, r /)
pvt = quicksort_median_of_3(tmp)

! partition ind by pvt
if (tmp(1) == pvt) then
  pvtind = l 
else if (tmp(2) == pvt) then
  pvtind = (l+r)/2
else
  pvtind = r
end if

call swap_ind(ind, pvtind, r)
call swap_ind(arr, pvtind, r)

call quicksort_by_index_partition(arr, ind, l, r, pvtind, itmp)
