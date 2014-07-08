! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module sort_by_index
  implicit none
  
  contains
  
  subroutine insertionsort_by_index(arr, ind, length)
    ! In/out
    integer, intent(in) :: length
    integer, intent(inout) :: ind(*)
    double precision, intent(inout) :: arr(*)
    ! local
    integer :: i, j
    integer :: itmp
    double precision :: dtmp
    
    
    if (length <= 1) return
    
    do i = 2, length
      itmp = ind(i)
      dtmp = arr(i)
      j = i
      do
        if (j > 1 .and. ind(j-1) > itmp) then
          ind(j) = ind(j-1)
          arr(j) = arr(j-1)
          j = j - 1
        else
          exit
        end if
      end do
        
        ind(j) = itmp
        arr(j) = dtmp
    end do
    
    return
  end subroutine
  
  
  
  subroutine quicksort_by_index_partition(arr, ind, l, r, pvtind, ret)
    use swaps
    implicit none
    ! In/out
    integer, intent(in) :: l, r, pvtind
    integer, intent(inout) :: ind(*)
    integer, intent(out) :: ret
    double precision, intent(inout) :: arr(*)
    ! local
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
    
    return
  end subroutine
  
  
  
  recursive subroutine quicksort_by_index_r(arr, ind, l, r)
    use swaps
    use quicksort_utils
    implicit none
    ! inputs
    integer, intent(in) :: l, r
    integer, intent(inout) :: ind(*)
    double precision, intent(inout) :: arr(*)
    ! local
    integer :: itmp, pvtind, pvt
    integer :: tmp(3)
    
    
    if (l < r) then
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
      
      ! finish with insertionsort or quicksort recursively
      if (r-l <= 10) then
        call insertionsort_by_index(arr(l), ind(l), r-l+1)
      else
        call quicksort_by_index_r(arr, ind, l, itmp)
        call quicksort_by_index_r(arr, ind, itmp+2, r)
      end if
    end if
    
    return
  end subroutine
  
  
  
  subroutine quicksort_by_index(arr, ind, length)
    implicit none
    ! inputs
    integer, intent(in) :: length
    integer, intent(inout) :: ind(*)
    double precision, intent(inout) :: arr(*)
    
    
    if (length <= 1) return
    
    call quicksort_by_index_r(arr, ind, 1, length)
    
    return
  end subroutine
  
end module
