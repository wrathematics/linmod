! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt



! Purpose
! =======
!
! Linear model fitter using a QR (or LQ if m<n) decomposition.
!
!
! Arguments
! =========
!
! m, n      (input) integer
!           The number of rows and columns of the data matrix a.
! 
! nrhs      (input) integer
!           The number of 'right hand sides', i.e. the number of 
!           columns of b.  nrhs>=0.
! 
! a         (input/output) double precision array, dimension (m,n).
!           On entry, the input data matrix.  On exit, the output
!           of a lapack QR factorization is stored for a.
! 
! b         (input) double precision array, dimension (m,nrhs).
!           The response variable.
! 
! lda, ldb  (input) integer
!           Leading dimension of a and b respectively.
! 
! tol       (input) double precision
!           Numerical tolerance.
! 
! coef      (output) double precision array, dimension (m,nrhs)
!           The coefficients array ("beta").
! 
! eff       (output) double precision array, dimension (m,nrhs)
!           The "effects" array, namely eff := Q^T * b.
! 
! ft        (output) double precision array, dimension (m,nrhs)
!           The "fitted" values, namely ft := Q*(R*b)
! 
! rsd       (output) double precision array, dimension (m,nrhs)
!           The residuals.
! 
! tau       
!           
! 
! jpvt      
!           
! 
! rank      (input/output) integer
!           On input, controls whether numerical rank should be checked
!           (-1 no check, otherwise check). On output, the estimated
!           numerical (columns) rank is returned.
! 
! info      (output) integer
!           = 0: successful exit.
!           < 0: if info = -i, the i-th argument had an illegal value.


module lm
  use, intrinsic :: iso_c_binding
  use :: linmod_omp
  use :: distributions, only : true, false
  implicit none
  
  
  interface
    subroutine rdgels(m, n, nrhs, a, b, tol, coef, eff, ft, rsd, tau, &
                      jpvt, rank, work, lwork, info) &
    bind(c, name='rdgels_')
      integer, intent(in) :: m, n, nrhs, lwork
      integer, intent(out) :: info, jpvt(n)
      integer, intent(inout) :: rank
      double precision, intent(in) :: tol
      double precision, intent(out) :: coef(n, nrhs), tau(*), work(lwork)
      double precision, intent(out), dimension(m, nrhs) :: ft, eff, rsd
      double precision, intent(inout) :: a(m, n), b(m, nrhs)
    end subroutine
  end interface
  
  
  
  contains
  
  subroutine lm_fit(use_offset, n, p, nrhs, x, y, offset, tol, coef, eff, ft, rsd, tau, &
                    jpvt, rank, info) &
  bind(c, name='lm_fit')
    logical(kind=c_bool) :: use_offset
    integer, intent(in) :: n, p, nrhs
    integer, intent(out) :: info, jpvt(p)
    integer, intent(inout) :: rank
    double precision, intent(in) :: offset(n), tol
    double precision, intent(out) :: coef(p, nrhs), tau(*)
    double precision, intent(out), dimension(n, nrhs) :: ft, eff, rsd
    double precision, intent(inout) :: x(n, p), y(n, nrhs)
    ! local
    integer :: np_max, lwork, i, j
    double precision :: tmp(1)
    double precision, allocatable :: work(:)
    
    
    np_max = max(1, n, p) ! FIXME why the HELL is this needed?!
    lwork = -1
    call dgels('n', n, p, nrhs, x, n, y, np_max, tmp, lwork, info)
    lwork = int(tmp(1))
    allocate(work(lwork))
    
    
    if (use_offset .eqv. true) then
      !$omp parallel if(n*nrhs > linmod_omp_minsize) private(i, j) default(shared) 
      !$omp do 
        do j = 1, nrhs
          do i = 1, n
            y(i, j) = y(i, j) - offset(i)
          end do
        end do
      !$omp end do
      !$omp end parallel
    end if
    
    
    call rdgels(n, p, nrhs, x, y, tol, coef, eff, ft, rsd, tau, jpvt, rank, work, lwork, info)
    
    
    if (use_offset .eqv. true) then
      !$omp parallel if(n*nrhs > linmod_omp_minsize) private(i, j) default(shared) 
      !$omp do 
        do j = 1, nrhs
          do i = 1, n
            ft(i, j) = ft(i, j) + offset(i)
          end do
        end do
      !$omp end do
      !$omp end parallel
    end if
    
    if (allocated(work)) deallocate(work)
    
  end subroutine
  
end module

