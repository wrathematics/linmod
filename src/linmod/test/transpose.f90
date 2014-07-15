subroutine print_matrix(m, n, a)
    integer, intent(in) :: m, n
    double precision, intent(in) :: a(m, n)
    integer :: i, j
    
    
    do j = 1, n
      do i = 1, m
        write(*, '(1000F14.2)', advance="no")( real(a(i, j)) )
      end do
      print *, ""
    end do
    
    return
  end subroutine



subroutine print_tester(m, n, a)
  use :: transposition
  integer :: m, n
  double precision :: a(m, n)
  
  print *, "Input:"
  call print_matrix(m, n, a)
  
  call xpose(m, n, a)
  
  print *, "Transposed:"
  call print_matrix(n, m, a)
  
  call xpose(n, m, a)
  
  print *, "Input again:"
  call print_matrix(m, n, a)
  
  return
end subroutine



subroutine print_blanks(n)
  integer, intent(in) :: n
  integer :: i
  
  do i = 1, n
    print *, ""
  end do
  
  return
end subroutine



program main
  integer, parameter :: m = 5, n = 3
  double precision :: a(m, n), b(n, n)
  a = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], shape(a))
  b = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9], shape(b))
  
  call print_blanks(1)
  print *, "----------------------- Square -----------------------"
  call print_tester(n, n, b)
  
  call print_blanks(3)
  print *, "----------------------- Non-Square -----------------------"
  call print_tester(m, n, a)
  
end program

