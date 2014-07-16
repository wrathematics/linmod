! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module distributions
  use, intrinsic :: iso_c_binding
  
  logical(kind=c_bool), parameter, public :: true = .true.
  logical(kind=c_bool), parameter, public :: false = .false.
  
  public :: dcauchy, pcauchy, qcauchy
  public :: dnorm, pnorm, qnorm
  
  
  ! ------------------------- Misc -------------------------
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
  
  
  ! ------------------------- Cauchy -------------------------
  interface
    real(kind=c_double) function dcauchy(x, location, scl, log_) bind(C)
      use, intrinsic :: iso_c_binding
      real(kind=c_double), intent(in), value :: x, location, scl
      logical(kind=c_bool), intent(in), value :: log_
    end function
  end interface
  
  interface
    real(kind=c_double) function pcauchy(q, location, scl, lower_tail, log_p) bind(C)
      use, intrinsic :: iso_c_binding
      real(kind=c_double), intent(in), value :: q, location, scl
      logical(kind=c_bool), intent(in), value :: lower_tail, log_p
    end function
  end interface
  
  interface
    real(kind=c_double) function qcauchy(p, location, scl, lower_tail, log_p) bind(C)
      use, intrinsic :: iso_c_binding
      real(kind=c_double), intent(in), value :: p, location, scl
      logical(kind=c_bool), intent(in), value :: lower_tail, log_p
    end function
  end interface
  
  
  ! ------------------------- Normal -------------------------
  interface
    real(kind=c_double) function dnorm(x, mean, sd, log_) bind(C)
      use, intrinsic :: iso_c_binding
      real(kind=c_double), intent(in), value :: x, mean, sd
      logical(kind=c_bool), intent(in), value :: log_
    end function
  end interface
  
  interface
    real(kind=c_double) function pnorm(q, mean, sd, lower_tail, log_p) bind(C)
      use, intrinsic :: iso_c_binding
      real(kind=c_double), intent(in), value :: q, mean, sd
      logical(kind=c_bool), intent(in), value :: lower_tail, log_p
    end function
  end interface
  
  interface
    real(kind=c_double) function qnorm(p, mean, sd, lower_tail, log_p) bind(C)
      use, intrinsic :: iso_c_binding
      real(kind=c_double), intent(in), value :: p, mean, sd
      logical(kind=c_bool), intent(in), value :: lower_tail, log_p
    end function
  end interface
  
end module
