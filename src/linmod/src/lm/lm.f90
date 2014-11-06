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



!!! TODO: s/m/n/, s/n/p/, s/a/x/ s/b/y/
module lm
  implicit none
  
  
  contains
  
  subroutine lm_fit(m, n, nrhs, a, b, tol, coef, eff, ft, rsd, tau, &
                    jpvt, rank, info) &
  bind(c, name='lm_fit')
    integer, intent(in) :: m, n, nrhs
    integer, intent(out) :: info, jpvt(n)
    integer, intent(inout) :: rank
    double precision, intent(in) :: tol
    double precision, intent(out) :: coef(n, *), tau(*)
    double precision, intent(out), dimension(m, *) :: ft, eff, rsd
    double precision, intent(inout) :: a(m, *), b(m, *)
    ! local
    integer :: lwork
    double precision :: tmp(1)
    double precision, allocatable :: work(:)
    
    
    lwork = -1
    call dgels('n', m, n, nrhs, a, m, b, m, tmp, lwork, info)
    lwork = int(tmp(1))
    allocate(work(lwork))
    
    call rdgels(m, n, nrhs, a, b, tol, coef, eff, ft, rsd, tau, jpvt, rank, work, lwork, info)
    
    if (allocated(work)) deallocate(work)
    
  end subroutine
  
end module

