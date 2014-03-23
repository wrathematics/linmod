! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module lapack
  implicit none
  
  
  interface
    
    !!! Legacy lmfit routines
    subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
      integer, intent(out) :: info
      double precision, intent(in) :: a(*)
      double precision, intent(inout) :: b(*)
      double precision, intent(out) :: work(*)
    end subroutine
    
    
    
    subroutine dgelqf(m, n, a, lda, tau, work, lwork, info)
      integer, intent(in) :: m, n, lda, lwork
      integer, intent(out) :: info
      double precision, intent(out) :: tau
      double precision, intent(inout) :: a(*), work(*)
    end subroutine
    
    
    
    subroutine dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
      character(len=1), intent(in) :: side, trans
      integer, intent(in) :: m, n, k, lda, ldc, lwork
      integer, intent(out) :: info
      double precision, intent(in) :: a(*), tau(*)
      double precision, intent(out) :: c(*), work(*)
    end subroutine
    
    
    
    subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
      integer, intent(in) :: m, n, lda, lwork
      integer, intent(out) :: info
      double precision, intent(inout) :: a(*)
      double precision, intent(out) :: tau(*), work(*)
    end subroutine
    
    
    
    !!! Other LAPACK routines
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(len=1), intent(in) :: transa, transb
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      double precision, intent(in) :: alpha, beta
      double precision, intent(in) :: a(*), b(*)
      double precision, intent(out) :: c(*)
    end subroutine
    
    
    
    function disnan(din) &
    result(val)
      logical :: val
      double precision, intent(in) :: din
    end function
    
    
    
    subroutine dlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info)
      character(len=1), intent(in) :: type
      integer, intent(in) :: kl, ku, m, n, lda
      integer, intent(out) :: info
      double precision, intent(in) :: cfrom, cto
      double precision, intent(inout) :: a(*)
    end subroutine
    
    
    
    subroutine dlaset(uplo, m, n, alpha, beta, a, lda)
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: m, n
      double precision, intent(in) :: alpha, beta
      double precision, intent(inout) :: a(*)
    end subroutine
    
    
    
    subroutine dgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, work, lwork, iwork, info)
      character(len=1), intent(in) :: joba, jobu, jobv, jobr, jobt, jobp
      integer, intent(in) :: m, n, lda, ldu, ldv, lwork
      integer, intent(out) :: iwork(*), info
      double precision, intent(in) :: a(*)
      double precision, intent(out) :: sva(*), u(*), v(*), work(*)
    end subroutine
    
    
    
    subroutine dgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
      integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
      integer, intent(in) :: jpvt(*)
      integer, intent(out) :: info
      double precision, intent(in) :: rcond
      double precision, intent(inout) :: a(*), b(*)
    end subroutine
    
    
    
    subroutine dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
      character(len=1), intent(in) :: uplo, trans, diag
      integer, intent(in) :: n, nrhs, lda, ldb
      integer, intent(out) :: info
      double precision, intent(in) :: a(*)
      double precision, intent(inout) :: b(*)
    end subroutine
    
    
    
    subroutine xerbla(srname, info)
      character*(*) :: srname
      integer :: info
    end subroutine
    
    
    
    function lsame(ca, cb) result(val)
      logical :: val
      character(len=1), intent(in) :: ca, cb
    end function
    
    
    
    function ilaenv(ispec, name, opts, n1, n2, n3, n4) result(value)
      integer :: value
      integer :: ispec
      character(len=1) :: name, opts
      integer :: n1, n2, n3, n4
    end function
    
    
    
    subroutine dgeqr2(m, n, a, lda, tau, work, info)
      integer, intent(in) :: m, n, lda
      integer, intent(out) :: info
      double precision, intent(in) :: work(*)
      double precision, intent(inout) :: a(*)
      double precision, intent(out) :: tau(*)
    end subroutine
    
    
    
    subroutine dlarf(side, m, n, v, incv, tau, c, ldc, work)
      character(len=1), intent(in) :: side
      integer, intent(in) :: m, n, incv, ldc
      double precision, intent(in) :: v(*), tau, work(*)
      double precision, intent(inout) :: c(*)
    end subroutine
    
    
    
    subroutine dlarfg(n, alpha, x, incx, tau)
      integer, intent(in) :: n, incx
      double precision, intent(inout) :: alpha, x(*), tau
    end subroutine
    
    
    
    subroutine dorm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
      character(len=1), intent(in) :: side, trans
      integer, intent(in) :: m, n, k, lda, ldc
      integer, intent(out) :: info
      double precision, intent(in) :: a(*), tau(*), work(*)
      double precision, intent(out) :: c(*)
    end subroutine
    
    
    
    subroutine dswap(n, x, incx, y, incy)
      integer, intent(in) :: n, incx, incy
      double precision, intent(inout) :: x(*), y(*)
    end subroutine
    
    
    
    function idamax(n, x, incx) result(index)
      integer :: index
      integer :: n, incx
      double precision :: x(*)
    end function
    
    
    
    function dlamch(cmach) result(val)
      double precision :: val
      character(len=1) :: cmach
    end function
    
    
    
    function dnrm2(n, x, incx) result(res)
      double precision :: res
      integer :: n, incx
      double precision :: x(*)
    end function
    
    
    
    subroutine dlarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork)
      character(len=1), intent(in) :: side, trans, direct, storev
      integer, intent(in) :: m, n, k, ldv, ldt, ldc, ldwork
      double precision, intent(in) :: v(*), t(*), work(*)
      double precision, intent(out) :: c(*)
    end subroutine
    
    
    
    subroutine dlarft(direct, storev, n, k, v, ldv, tau, t, ldt)
      character(len=1), intent(in) :: direct, storev
      integer, intent(in) :: n, k, ldv, ldt
      double precision, intent(in) :: tau(*)
      double precision, intent(out) :: t(*)
      double precision, intent(inout) :: v(*)
    end subroutine
    
    
    
    subroutine dlabad(small, large)
      double precision, intent(inout) :: small, large
    end subroutine
    
    
    
    function dlange(norm, m, n, a, lda, work) result(val)
      double precision :: val
      character(len=1) :: norm
      integer :: m, n, lda
      double precision :: a(*), work(*)
    end function
    
    
    
    subroutine dormlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
      character(len=1), intent(in) :: side, trans
      integer, intent(in) :: m, n, k, lda, ldc, lwork
      integer, intent(out) :: info
      double precision, intent(in) :: a(*), tau(*)
      double precision, intent(inout) :: c(*), work(*)
    end subroutine
    
  end interface
  
end module
