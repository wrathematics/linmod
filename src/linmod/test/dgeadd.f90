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
  
  !$ print *, openmp_version
  
  if (trans == 'N' .or. trans == 'n') then
    !$omp parallel if(m*n > 5000) private(i, j, ind) default(shared) 
    !$omp do 
      do j = 1, n
        print *, omp_get_num_threads()
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



!subroutine printer(m, n, a)
!  do i = 1, m
!    do j = 1, n
!      print *, real(a(i, j))
!    end do
!  end do
!end subroutine



program addtwo
  integer, parameter :: m = 3
  integer, parameter :: n = 4
  double precision :: x(m, n)
  double precision :: y(m, n)
  
  do j = 1, n
    do i = 1, m
      y(i, j) = i + (j-1)*m
    end do
  end do
  
  call dgeadd('N', m, n, 0.0d0, x, m, 1.0d0, y, m)
  
  
  print *, int(y)
end program
