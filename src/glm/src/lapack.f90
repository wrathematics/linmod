module lapack
  implicit none
  
  
  interface
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(len=1), intent(in) :: transa, transb
      integer, intent(in) :: m, n, k
      integer, intent(in) :: lda, ldb, ldc
      double precision, intent(in) :: alpha, beta
      double precision, intent(in) :: a(*), b(*)
      double precision, intent(out) :: c(*)
    end subroutine
  end interface
  
  
  
  interface
    subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, nrhs
      integer, intent(in) :: lda, ldb, lwork
      integer, intent(out) :: info
      double precision, intent(in) :: a(*)
      double precision, intent(inout) :: b(*)
      double precision, intent(out) :: work(*)
    end subroutine
  end interface
  
  
  
  interface
    function disnan(din) &
    result(val)
      logical :: val
      double precision, intent(in) :: din
    end function
  end interface
  
end module
