! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


! OMP statements added to dlacpy, original copyright:
!  -- lapack auxiliary routine (version 3.1) --
!     univ. of tennessee, univ. of california berkeley and nag ltd..
!     november 2006



!  uplo    (input) character*1
!          specifies the part of the matrix a to be copied to b.
!          = 'u':      upper triangular part
!          = 'l':      lower triangular part
!          otherwise:  all of the matrix a
!
!  m       (input) integer
!          the number of rows of the matrix a.  m >= 0.
!
!  n       (input) integer
!          the number of columns of the matrix a.  n >= 0.
!
!  a       (input) double precision array, dimension (lda,n)
!          the m by n matrix a.  if uplo = 'u', only the upper triangle
!          or trapezoid is accessed; if uplo = 'l', only the lower
!          triangle or trapezoid is accessed.
!
!  lda     (input) integer
!          the leading dimension of the array a.  lda >= max(1,m).
!
!  b       (output) double precision array, dimension (ldb,n)
!          on exit, b = a in the locations specified by uplo.
!
!  ldb     (input) integer
!          the leading dimension of the array b.  ldb >= max(1,m).
!
subroutine dlacpy_omp(uplo, m, n, a, lda, b, ldb)
  use string_tools
  implicit none
  
  ! in/out
  character(len=1), intent(in) :: uplo
  integer, intent(in) :: lda, ldb, m, n
  double precision, intent(in) :: a(lda, *)
  double precision, intent(out) :: b(ldb, *)
  ! local
  integer            i, j, itmp
  intrinsic          min
  
  
  if(equivchar(uplo, 'u')) then
    !$omp parallel if(m*n > 5000) private(i, j) default(shared) 
    !$omp do
      do j = 1, n
        do i = 1, min(j, m)
          b(i, j) = a(i, j)
        end do
      end do
    !$omp end do
    !$omp end parallel
  else if(equivchar(uplo, 'l')) then
    !$omp parallel if(m*n > 5000) private(i, j) default(shared) 
    !$omp do
      do j = 1, n
        do i = j, m
          b(i, j) = a(i, j)
        end do
      end do
    !$omp end do
    !$omp end parallel
  else
    !$omp parallel if(m*n > 5000) private(i, j) default(shared) 
    !$omp do
      do j = 1, n
        do i = 1, m
          b(i, j) = a(i, j)
        end do
      end do
    !$omp end do
    !$omp end parallel
  end if
  
  return
end


