! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module insertionsorts
  use :: swaps
  implicit none
  
  
  interface insertionsort
    module procedure iinsertionsort, sinsertionsort, dinsertionsort
  end interface
  
  interface insertionsort_by_index
    module procedure iinsertionsort_by_index, sinsertionsort_by_index, dinsertionsort_by_index
  end interface
  
  
  contains
  
  
  ! --------------------------------------------------------
  ! Sort
  ! --------------------------------------------------------
  
  subroutine iinsertionsort(x, xlen)
    ! In/out
    integer, intent(in) :: xlen
    integer, intent(inout) :: x(*)
    ! local
    integer :: i, j
    integer :: tmp
    
    include 'include/insertionsort_generic.inc'
    
    return
  end subroutine
  
  
  
  subroutine sinsertionsort(x, xlen)
    ! In/out
    integer, intent(in) :: xlen
    real, intent(inout) :: x(*)
    ! local
    integer :: i, j
    real :: tmp
    
    include 'include/insertionsort_generic.inc'
    
    return
  end subroutine
  
  
  
  subroutine dinsertionsort(x, xlen)
    ! In/out
    integer, intent(in) :: xlen
    double precision, intent(inout) :: x(*)
    ! local
    integer :: i, j
    double precision :: tmp
    
    include 'include/insertionsort_generic.inc'
    
    return
  end subroutine
  
  
  ! --------------------------------------------------------
  ! Sort by index
  ! --------------------------------------------------------
  
  subroutine iinsertionsort_by_index(arr, ind, length)
    ! In/out
    integer, intent(in) :: length
    integer, intent(inout) :: ind(*)
    integer, intent(inout) :: arr(*)
    ! local
    integer :: tmp
    
    include 'include/insertionsort_by_index_generic.inc'
    
    return
  end subroutine
  
  
  
  subroutine sinsertionsort_by_index(arr, ind, length)
    ! In/out
    integer, intent(in) :: length
    integer, intent(inout) :: ind(*)
    real, intent(inout) :: arr(*)
    ! local
    real :: tmp
    
    include 'include/insertionsort_by_index_generic.inc'
    
    return
  end subroutine
  
  
  
  subroutine dinsertionsort_by_index(arr, ind, length)
    ! In/out
    integer, intent(in) :: length
    integer, intent(inout) :: ind(*)
    double precision, intent(inout) :: arr(*)
    ! local
    double precision :: tmp
    
    include 'include/insertionsort_by_index_generic.inc'
    
    return
  end subroutine
  
end module
