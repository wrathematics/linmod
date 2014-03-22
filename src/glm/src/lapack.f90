module lapack
  implicit none
  
  
  interface
    
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(len=1), intent(in) :: transa, transb
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      double precision, intent(in) :: alpha, beta
      double precision, intent(in) :: a(*), b(*)
      double precision, intent(out) :: c(*)
    end subroutine
    
    
    
    subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
      integer, intent(out) :: info
      double precision, intent(in) :: a(*)
      double precision, intent(inout) :: b(*)
      double precision, intent(out) :: work(*)
    end subroutine
    
    
    
    function disnan(din) &
    result(val)
      logical :: val
      double precision, intent(in) :: din
    end function
    
    
    
!    subroutine dgelqf(m, n, a, lda, tau, work, lwork, info)
!      integer, intent(in) :: m, n, lda, lwork
!      integer, intent(out) :: info
!      double precision, intent(out) :: tau
!      double precision, intent(inout) :: a(*), work(*)
!    end subroutine
!    
!    
!    
!    subroutine dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
!      character(len=1), intent(in) :: side, trans
!      integer, intent(in) :: m, n, k, lda, ldc, lwork
!      integer, intent(out) :: info
!      double precision, intent(in) :: a(*), tau(*)
!      double precision, intent(out) :: c(*), work(*)
!    end subroutine
!    
!    
!    
!    subroutine dlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info)
!      character(len=1), intent(in) :: type
!      integer, intent(in) :: kl, ku, m, n, lda
!      integer, intent(out) :: info
!      double precision, intent(in) :: cfrom(*), cto(*)
!      double precision, intent(inout) :: a(*)
!    end subroutine
!    
!    
!    
!    subroutine dlaset(uplo, m, n, alpha, beta, a, lda)
!      character(len=1), intent(in) :: uplo
!      integer, intent(in) :: m, n
!      double precision, intent(in) :: alpha, beta
!      double precision, intent(inout) :: a(*)
!    end subroutine
!    
!    
!    
!    subroutine dgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, work, lwork, iwork, info)
!      character(len=1), intent(in) :: joba, jobu, jobv, jobr, jobt, jobp
!      integer, intent(in) :: m, n, lda, ldu, ldv, lwork
!      integer, intent(out) :: iwork(*), info
!      double precision, intent(in) :: a(*)
!      double precision, intent(out) :: sva(*), u(*), v(*), work(*)
!    end subroutine
!    
!    
!    
!    subroutine dgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
!      integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
!      integer, intent(in) :: jpvt(*)
!      integer, intent(out) :: info
!      double precision, intent(in) :: rcond
!      double precision, intent(inout) :: a(*), b(*)
!    end subroutine
!    
!    
!    
!    subroutine dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
!      character(len=1), intent(in) :: uplo, trans, diag
!      integer, intent(in) :: n, nrhs, lda, ldb
!      integer, intent(out) :: info
!      double precision, intent(in) :: a(*)
!      double precision, intent(inout) :: b(*)
!    end subroutine
!    
!    
!    
!    subroutine xerbla(srname, info)
!      character*(*) :: srname
!      integer, intent(out) :: info
!    end subroutine
!    
!    
!    
!    function lsame(ca, cb) result(val)
!      logical :: val
!      character(len=1), intent(in) :: ca, cb
!    end function
!    
!    
!    
!    function ilaenv(ispec, name, opts, n1, n2, n3, n4) result(value)
!      
!    end function
    
  end interface
  
end module
