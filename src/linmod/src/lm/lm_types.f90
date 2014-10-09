! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module lm_types
  use, intrinsic :: iso_c_binding
  implicit none
  
  
  type, bind(c) :: qr_t
    integer(c_int) :: n
    integer(c_int) :: p
    real(c_double) :: qr(*)
    real(c_double) :: qraux(*)
    integer(c_int) :: pivot(*)
    real(c_double) :: tol
    integer(c_int) :: rank
  end type myftype
  
  
  contains
  
  
  subroutine recover_qr(qrlist, n, p, qr, qraux, pivot, tol, rank)
    type(qr_t) :: qrlist
    integer(c_int) :: n, p, rank
    real(c_double) :: tol
    real(c_double), pointer, dimension(:) :: qr, qraux
    integer(c_int), pointer, dimension(:) :: pivot
    
    n = qrlist%n
    p = qrlist%p
    tol = qrlist%tol
    rank = qrlist%rank
    
    call c_f_pointer(qrlist%qr, qr, (/n*p/))
    call c_f_pointer(qrlist%qr, qraux, (/p/))
    call c_f_pointer(qrlist%qr, pivot, (/p/))
    
  end subroutine
  
end module
