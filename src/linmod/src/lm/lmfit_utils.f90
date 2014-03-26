! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module lmfit_utils
  implicit none
  
  
  interface
    
    subroutine dgeadd(trans, m, n, alpha, a, lda, beta, c, ldc)
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, lda, ldc
      double precision, intent(in) :: a(*), alpha, beta
      double precision, intent(inout) :: c(*)
    end subroutine
    
  end interface
  
end module
