! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module lmfit
  implicit none
  
  
  interface
    
    subroutine rdgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
      integer, intent(out) :: info
      double precision, intent(in) :: a(*)
      double precision, intent(inout) :: b(*)
      double precision, intent(out) :: work(*)
    end subroutine
    
    
    
    subroutine rdgelqf(m, n, a, lda, tau, work, lwork, info)
      integer, intent(in) :: m, n, lda, lwork
      integer, intent(out) :: info
      double precision, intent(out) :: tau
      double precision, intent(inout) :: a(*), work(*)
    end subroutine
    
    
    
    subroutine rdormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
      character(len=1), intent(in) :: side, trans
      integer, intent(in) :: m, n, k, lda, ldc, lwork
      integer, intent(out) :: info
      double precision, intent(in) :: a(*), tau(*)
      double precision, intent(out) :: c(*), work(*)
    end subroutine
    
    
    
    subroutine rdgeqp3(m, n, a, lda, jpvt, tau, work, lwork, tol, rank, info)
      integer, intent(in) :: m, n, lda, lwork
      integer, intent(out) :: rank, info
      double precision, intent(in) :: tol
      integer, intent(inout) :: jpvt(*)
      double precision, intent(inout) :: a(lda, *), tau(*), work(*)
    end subroutine
    
    
  end interface
  
end module
