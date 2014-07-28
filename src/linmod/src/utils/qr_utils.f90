! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module qr_utils
  implicit none
  
  
  contains
  
  
  subroutine qr_Q(m, n, rank, QR, qraux, Q, info)
    integer, intent(in) :: m, n, rank
    integer, intent(out) :: info
    double precision, intent(in) :: QR(m, n), qraux(*)
    double precision, intent(out) :: Q(m, n)
    ! local
    integer :: i, j, lwork
    double precision :: tmp(1)
    double precision, allocatable :: work(:)
    
    
    call dlacpy('A', m, n, QR, m, Q, m)
    
    lwork = -1
    call dorgqr(m, n, rank, Q, m, qraux, tmp, lwork, info)
    
    lwork = int(tmp(1))
    allocate(work(lwork))
    call dorgqr(m, n, rank, Q, m, qraux, work, lwork, info)
    
    deallocate(work)
    
    return
  end subroutine
  
  
  
  subroutine qr_R(m, n, QR, R)
    integer, intent(in) :: m, n
    double precision, intent(in) :: QR(m, n)
    double precision, intent(out) :: R(n, n)
    ! local
    integer :: i, j
    
    !$omp parallel do if(n*n > 5000) private(i, j) default(shared)
      do j = 1, n
        do i = 1, n
          R(i, j) = QR(i, j)
        end do
      end do
    !$omp end parallel do
    
    return
  end subroutine
end module

