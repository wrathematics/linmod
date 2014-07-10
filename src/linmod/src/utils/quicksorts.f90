! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module quicksorts
  use :: swaps
  use :: insertionsorts
  use :: quicksort_utils
  implicit none
  
  
  interface quicksort
    module procedure iquicksort
    module procedure squicksort
    module procedure dquicksort
  end interface
  
  interface quicksort_by_index
    module procedure iquicksort_by_index
    module procedure squicksort_by_index
    module procedure dquicksort_by_index
  end interface
  
  public :: quicksort
  public :: quicksort_by_index
  
  private :: iquicksort, squicksort, dquicksort
  private :: iquicksort_by_index, squicksort_by_index, dquicksort_by_index
  
  contains
  
  
  ! --------------------------------------------------------
  ! Sort
  ! --------------------------------------------------------
  
  recursive subroutine iquicksort_r(x, l, r)
    ! inputs
    integer, intent(in) :: l, r
    integer, intent(inout), dimension(:) :: x
    ! local
    integer :: itmp, ind
    integer :: pvt
    integer :: tmp(3)
    
    
    include 'include/quicksort_r_generic.inc'
    
      if (r-l <= 10) then
        call insertionsort(x(l:), r-l+1)
      else
        ! sort recursively
        call iquicksort_r(x, l, itmp)
        call iquicksort_r(x, itmp+2, r)
      end if
    end if
    
    return
  end subroutine
  
  
  
  subroutine iquicksort(x, xlen)
    ! inputs
    integer, intent(in) :: xlen
    integer, intent(inout) :: x(:)
    
    call iquicksort_r(x, 1, xlen)
    
  end subroutine
  
  
  
  recursive subroutine squicksort_r(x, l, r)
    ! inputs
    integer, intent(in) :: l, r
    real, intent(inout) :: x(:)
    ! local
    integer :: itmp, ind
    real :: pvt
    real :: tmp(3)
    
    
    include 'include/quicksort_r_generic.inc'
    
      if (r-l <= 10) then
        call insertionsort(x(l:), r-l+1)
      else
        ! sort recursively
        call squicksort_r(x, l, itmp)
        call squicksort_r(x, itmp+2, r)
      end if
    end if
    
    return
  end subroutine
  
  
  
  subroutine squicksort(x, xlen)
    ! inputs
    integer, intent(in) :: xlen
    real, intent(inout) :: x(:)
    
    call squicksort_r(x, 1, xlen)
    
  end subroutine
  
  
  
  recursive subroutine dquicksort_r(x, l, r)
    ! inputs
    integer, intent(in) :: l, r
    double precision, intent(inout) :: x(:)
    ! local
    integer :: itmp, ind
    double precision :: pvt
    double precision :: tmp(3)
    
    
    include 'include/quicksort_r_generic.inc'
    
      if (r-l <= 10) then
        call insertionsort(x(l:), r-l+1)
      else
        ! sort recursively
        call dquicksort_r(x, l, itmp)
        call dquicksort_r(x, itmp+2, r)
      end if
    end if
    
    return
  end subroutine
  
  
  
  subroutine dquicksort(x, xlen)
    ! inputs
    integer, intent(in) :: xlen
    double precision, intent(inout) :: x(:)
    
    call dquicksort_r(x, 1, xlen)
    
  end subroutine
  
  
  
  ! --------------------------------------------------------
  ! Sort by index
  ! --------------------------------------------------------
  
  recursive subroutine iquicksort_by_index_r(arr, ind, l, r)
    implicit none
    ! inputs
    integer, intent(in) :: l, r
    integer, intent(inout) :: ind(:)
    integer, intent(inout) :: arr(:)
    ! local
    integer :: itmp, pvtind, pvt
    integer :: tmp(3)
    
    
    if (l < r) then
      include 'include/quicksort_by_index_r_generic.inc'
      
      ! finish with insertionsort or quicksort recursively
      if (r-l <= 10) then
        call insertionsort_by_index(arr(l:), ind(l:), r-l+1)
      else
        call iquicksort_by_index_r(arr, ind, l, itmp)
        call iquicksort_by_index_r(arr, ind, itmp+2, r)
      end if
    end if
    
    return
  end subroutine
  
  
  
  subroutine iquicksort_by_index(arr, ind, length)
    implicit none
    ! inputs
    integer, intent(in) :: length
    integer, intent(inout) :: ind(:)
    integer, intent(inout) :: arr(:)
    
    
    if (length <= 1) return
    
    call iquicksort_by_index_r(arr, ind, 1, length)
    
    return
  end subroutine
  
  
  
  recursive subroutine squicksort_by_index_r(arr, ind, l, r)
    implicit none
    ! inputs
    integer, intent(in) :: l, r
    integer, intent(inout) :: ind(:)
    real, intent(inout) :: arr(:)
    ! local
    integer :: itmp, pvtind, pvt
    integer :: tmp(3)
    
    
    if (l < r) then
      include 'include/quicksort_by_index_r_generic.inc'
      
      ! finish with insertionsort or quicksort recursively
      if (r-l <= 10) then
        call insertionsort_by_index(arr(l:), ind(l:), r-l+1)
      else
        call squicksort_by_index_r(arr, ind, l, itmp)
        call squicksort_by_index_r(arr, ind, itmp+2, r)
      end if
    end if
    
    return
  end subroutine
  
  
  
  subroutine squicksort_by_index(arr, ind, length)
    implicit none
    ! inputs
    integer, intent(in) :: length
    integer, intent(inout) :: ind(:)
    real, intent(inout) :: arr(:)
    
    
    if (length <= 1) return
    
    call squicksort_by_index_r(arr, ind, 1, length)
    
    return
  end subroutine
  
  
  
  recursive subroutine dquicksort_by_index_r(arr, ind, l, r)
    implicit none
    ! inputs
    integer, intent(in) :: l, r
    integer, intent(inout) :: ind(:)
    double precision, intent(inout) :: arr(:)
    ! local
    integer :: itmp, pvtind, pvt
    integer :: tmp(3)
    
    
    if (l < r) then
      include 'include/quicksort_by_index_r_generic.inc'
      
      ! finish with insertionsort or quicksort recursively
      if (r-l <= 10) then
        call insertionsort_by_index(arr(l:), ind(l:), r-l+1)
      else
        call dquicksort_by_index_r(arr, ind, l, itmp)
        call dquicksort_by_index_r(arr, ind, itmp+2, r)
      end if
    end if
    
    return
  end subroutine
  
  
  
  subroutine dquicksort_by_index(arr, ind, length)
    implicit none
    ! inputs
    integer, intent(in) :: length
    integer, intent(inout) :: ind(:)
    double precision, intent(inout) :: arr(:)
    
    
    if (length <= 1) return
    
    call dquicksort_by_index_r(arr, ind, 1, length)
    
    return
  end subroutine
  
end module

