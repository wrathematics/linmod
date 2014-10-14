! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module lm
  implicit none
  
  
  contains
  
  subroutine lm_fit(m, n, nrhs, a, b, tol, coef, eff, ft, rsd, tau, &
                    jpvt, rank, info) &
  bind(c, name='lm_fit')
    integer, intent(in) :: m, n, nrhs
    integer, intent(out) :: info, jpvt(n)
    integer, intent(inout) :: rank
    double precision, intent(in) :: tol
    double precision, intent(out) :: coef(n, *), tau(*)
    double precision, intent(out), dimension(m, *) :: ft, eff, rsd
    double precision, intent(inout) :: a(m, *), b(m, *)
    ! local
    integer :: lwork
    double precision :: tmp(1)
    double precision, allocatable :: work(:)
    
    
    lwork = -1
    call dgels('n', m, n, nrhs, a, m, b, m, tmp, lwork, info)
    lwork = int(tmp(1))
    allocate(work(lwork))
    
    call rdgels(m, n, nrhs, a, b, tol, coef, eff, ft, rsd, tau, jpvt, rank, work, lwork, info)
    
    if (allocated(work)) deallocate(work)
    
  end subroutine
  
end module

