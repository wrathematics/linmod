! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module distributions
  use, intrinsic :: iso_c_binding
  
!  public :: probit
  
  interface
    real(kind=c_double) function erfinv(x) bind(C)
      use, intrinsic :: iso_c_binding
      real(kind=c_double), intent(in), value :: x
    end function
  end interface
  
  interface
    real(kind=c_double) function erfinvc(x) bind(C)
      use, intrinsic :: iso_c_binding
      real(kind=c_double), intent(in), value :: x
    end function
  end interface
  
  interface
    real(kind=c_double) function probit(x) bind(C)
      use, intrinsic :: iso_c_binding
      real(kind=c_double), intent(in), value :: x
    end function
  end interface
end module
