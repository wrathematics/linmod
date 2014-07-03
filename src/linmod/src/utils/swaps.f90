! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


! swap two elements of an array
!  inputs:
!   x = array
!     i = index 1
!     j = index 2
module swaps
  implicit none
  
  
  interface swap
    module procedure iswap, sswap, dswap, scswap, dcswap
  end interface
  
  
  interface swap_ind
    module procedure iswap_ind, sswap_ind, dswap_ind, scswap_ind, dcswap_ind
  end interface
  
  
  contains
  
  
  ! Swap values
  subroutine iswap(x, y)
    integer, intent(inout) :: x, y
    integer :: tmp
    
    include 'include/swap_generic.inc'
    
    return
  end subroutine
  
  
  
  subroutine sswap(x, y)
    implicit none
    real, intent(inout) :: x, y
    real :: tmp
    
    include 'include/swap_generic.inc'
    
    return
  end subroutine
  
  
  
  subroutine dswap(x, y)
    implicit none
    double precision, intent(inout) :: x, y
    double precision :: tmp
    
    include 'include/swap_generic.inc'
    
    return
  end subroutine
  
  
  
  subroutine scswap(x, y)
    implicit none
    complex, intent(inout) :: x, y
    complex :: tmp
    
    include 'include/swap_generic.inc'
    
    return
  end subroutine
  
  
  
  subroutine dcswap(x, y)
    implicit none
    double complex, intent(inout) :: x, y
    double complex :: tmp
    
    include 'include/swap_generic.inc'
    
    return
  end subroutine
  
  
  
  ! Swap values at indices
  subroutine iswap_ind(x, i, j)
    integer, intent(inout) :: x(*)
    integer, intent(in) :: i, j
    integer :: tmp
    
    include 'include/swap_generic_indices.inc'
    
    return
  end subroutine
  
  
  
  subroutine sswap_ind(x, i, j)
    implicit none
    real, intent(inout) :: x(*)
    integer, intent(in) :: i, j
    real :: tmp
    
    include 'include/swap_generic_indices.inc'
    
    return
  end subroutine
  
  
  
  subroutine dswap_ind(x, i, j)
    implicit none
    double precision, intent(inout) :: x(*)
    integer, intent(in) :: i, j
    double precision :: tmp
    
    include 'include/swap_generic_indices.inc'
    
    return
  end subroutine
  
  
  
  subroutine scswap_ind(x, i, j)
    implicit none
    complex, intent(inout) :: x(*)
    integer, intent(in) :: i, j
    complex :: tmp
    
    include 'include/swap_generic_indices.inc'
    
    return
  end subroutine
  
  
  
  subroutine dcswap_ind(x, i, j)
    implicit none
    double complex, intent(inout) :: x(*)
    integer, intent(in) :: i, j
    double complex :: tmp
    
    include 'include/swap_generic_indices.inc'
    
    return
  end subroutine
  
end module

