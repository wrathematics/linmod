program main
  use :: distributions
  
  double precision :: x
  
  x = 0.78d0
  
  print *, erfinv(x)
  print *, erfinvc(x)
  print *, probit(x)
end program
