! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


module glm
  implicit none
  
  interface
    subroutine glm_fit(family, link, intercept, stoprule, n, p, x, y, &
                   beta, wt, offset, resids, maxiter, tol, info) &
      bind(C, name='glm_fit_fortran_')
      use, intrinsic :: iso_c_binding
      use :: lapack
      use :: glm_check
      use :: glm_loglik_utils
      use :: glm_link_utils
      use :: glm_family_utils
      use :: glm_update_utils
      implicit none
      
      integer, intent(in) :: family, link, intercept, stoprule
      integer, intent(in) :: n, p, maxiter
      integer, intent(out) :: info
      double precision, intent(in) :: x(n,p), y(n), offset(n), tol
      double precision, intent(out) :: beta(p), wt(n), resids(n)
    end subroutine
  end interface
end module

