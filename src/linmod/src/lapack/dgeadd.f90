! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


!FIXME at the moment, lda == m == ldc

!C = alpha*op(A) + beta*C
!
! m = the number of ROWS of C
! n = number of rows of COLS of C
!
! if trans = 'N', then dim(a) == (m, n)
! if trans = 'T', then dim(a) == (n, m)
subroutine dgeadd(trans, m, n, alpha, a, lda, beta, c, ldc)
  !$ use omp_lib
  
  implicit none
  ! in/out
  character(len=1), intent(in) :: trans
  integer, intent(in) :: m, n, lda, ldc
  double precision, intent(in) :: a(*), alpha, beta
  double precision, intent(inout) :: c(*)
  ! local
  integer :: i, j
  integer :: ind, tind
  
  
  if (trans == 'N' .or. trans == 'n') then
    !$omp parallel if(m*n > 5000) private(i, j, ind) default(shared) 
    !$omp do 
      do j = 1, n
        do i = 1, m
          ind = i + (j-1)*n
          c(ind) = alpha*a(ind) + beta*c(ind)
        end do
      end do
    !$omp end do
    !$omp end parallel
    
  else if (trans == 'T' .or. trans == 't') then
    !$omp parallel if(m*n > 5000) private(i, j, ind, tind) default(shared) 
    !$omp do
      do j = 1, n
        do i = 1, m
          ind = i + (j-1)*n
          tind = j + (i-1)*m
          c(ind) = alpha*a(tind) + beta*c(ind)
        end do
      end do
    !$omp end do
    !$omp end parallel
  end if
  
  return
end subroutine
