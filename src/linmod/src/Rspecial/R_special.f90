module R_special
  use, intrinsic :: iso_c_binding
  implicit none
  
  interface r_set_na
    procedure r_set_na_real, r_set_na_int
  end interface
  
  interface r_set_nan
    procedure r_set_nan_real
  end interface
  
  interface r_set_na_real
    subroutine r_set_na_real(val) &
    bind(C, name="r_set_na_real_")
      import :: C_DOUBLE
      real(kind=C_DOUBLE), intent(inout) :: val
    end subroutine
  end interface
  
  interface r_set_na_int
    subroutine r_set_na_int(val) &
    bind(C, name="r_set_na_int_")
      import :: C_INT
      integer(kind=C_INT), intent(inout) :: val
    end subroutine
  end interface
  
  interface r_set_nan_real
    subroutine r_set_nan_real(val) &
    bind(C, name="r_set_nan_real_")
      import :: C_DOUBLE
      real(kind=C_DOUBLE), intent(inout) :: val
    end subroutine
  end interface
  
end module

