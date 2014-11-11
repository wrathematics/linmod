! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module qr_utils
  use :: linmod_omp
  implicit none
  
  
  interface
    
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
    
    
    
    subroutine rdlaqp2(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work, tol, rank)
      integer, intent(in) :: m, n, offset, lda
      integer, intent(out) :: jpvt(*)
      integer, intent(out) :: rank
      double precision, intent(in) :: tol
      double precision, intent(inout) :: a(lda, *), tau(*), vn1(*), vn2(*), work(*)
    end subroutine
    
  end interface
  
  
  
  
  contains
  
  
  subroutine qr_Q(m, n, rank, QR, qraux, Q, info)
    integer, intent(in) :: m, n, rank
    integer, intent(out) :: info
    double precision, intent(in) :: QR(m, n), qraux(*)
    double precision, intent(out) :: Q(m, n)
    ! local
    integer :: i, j, lwork
    double precision :: tmp(1)
    double precision, allocatable :: work(:)
    
    
    call dlacpy('A', m, n, QR, m, Q, m)
    
    lwork = -1
    call dorgqr(m, n, rank, Q, m, qraux, tmp, lwork, info)
    
    lwork = int(tmp(1))
    allocate(work(lwork))
    call dorgqr(m, n, rank, Q, m, qraux, work, lwork, info)
    
    deallocate(work)
    
    return
  end subroutine
  
  
  
  subroutine qr_R(m, n, QR, R)
    integer, intent(in) :: m, n
    double precision, intent(in) :: QR(m, n)
    double precision, intent(out) :: R(n, n)
    ! local
    integer :: i, j
    
    !$omp parallel do if(n*n > linmod_omp_minsize) private(i, j) default(shared)
      do j = 1, n
        do i = 1, n
          R(i, j) = QR(i, j)
        end do
      end do
    !$omp end parallel do
    
    return
  end subroutine
  
end module

