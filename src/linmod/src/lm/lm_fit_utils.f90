! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module lm_fit_utils
  use :: R_special
  use :: quicksorts!, only : quicksort_by_index
  use :: lapack
  use :: qr_utils, only : rdgeqp3
  implicit none
  
  contains
  
  ! m >= n
  subroutine rdgels_qr(m, n, mn, nrhs, a, b, work, lwork, info, &
                    tol, coef, eff, ft, rsd, tau, jpvt, rank, qraux1)
    ! in/out
    integer, intent(in) :: m, n, mn, nrhs, lwork
    integer, intent(out) :: info, jpvt(n)
    integer, intent(inout) :: rank
    double precision, intent(in) :: tol
    double precision, intent(out) :: work(lwork), coef(n, nrhs), tau(*)
    double precision, intent(out), dimension(m, nrhs) :: ft, eff, rsd
    double precision, intent(inout) :: a(m, n), b(m, nrhs)
    double precision, intent(out) :: qraux1
    ! local
    integer i, j
    
    
    ! Assume model matrix is full rank
    if (rank == -1) then 
      rank = n
      call dgeqrf(m, n, a, m, work(1), work(mn+1), lwork-mn, info)
      
      do i = 1, n
        jpvt(i) = i
      end do
      
    ! RRQR
    else
      call rdgeqp3(m, n, a, m, jpvt, work(1), work(mn+1), 3*n+1, tol, rank, info)
    end if
    
    
    ! least-squares problem min || a * x - b ||
    !  b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
    call dormqr('left', 'transpose', m, nrhs, rank, a, m, work(1), b, m, work(mn+1), lwork-mn, info)
    
    ! Store "effects"
!    call dlacpy_omp('All', m, nrhs, b, m, eff, m)
    eff(1:m, 1:nrhs) = b(1:m, 1:nrhs)
    
    
    ! Store qraux(1) == work(1) 
    qraux1 = work(1)
    
    ! workspace at least nrhs, optimally nrhs*nb
    ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
    call dtrtrs('upper', 'no transpose', 'non-unit', rank, nrhs, a, m, b, m, info)
    
    
    !!! Produce fitted.values = Ax = Q*(R*x)
    
    ! Copy over the first RANK elements of numerical solution X
    !$omp parallel if(m*n > 5000) private(i, j) default(shared) 
    !$omp do
      do j = 1, nrhs
        do i = 1, rank
          ft(i, j) = b(i, j)
        end do
        
        do i = rank+1, m
          ft(i, j) = 0.0d0
        end do
      end do
    !$omp end do
    !$omp end parallel
    
    
    ! Compute fitted FT = Q*(R*fitted)
    call dtrmm('L', 'U', 'N', 'N', mn, nrhs, 1.0d0, a, m, ft, m)
    
    tau(1:mn) = work(1:mn)
    call dormqr('L', 'N', m, nrhs, mn, a, m, tau, ft, m, work, lwork, info)
    
    
    !!! FIXME dlacpy_omp doesn't copy beyond first column correctly ?!
    ! Compute residual RSD = B - FT
!    call dgeadd_omp('N', m, nrhs, -1.0d0, ft, m, 1.0d0, rsd, m)
    rsd(1:m, 1:nrhs) = rsd(1:m, 1:nrhs) - ft(1:m, 1:nrhs)
    
    ! Coefficients are stored in the first RANK elements of B
!    call dlacpy_omp('A', n, nrhs, b, 1, coef, 1)
    coef(1:n, 1:nrhs) = b(1:n, 1:nrhs)
    
    call rdgels_fixcoef(m, n, mn, nrhs, rank, jpvt, coef)
    
    return
  end subroutine
  
  
  
  ! m < n
  subroutine rdgels_lq(m, n, mn, nrhs, a, b, work, lwork, info, &
                    tol, coef, eff, ft, rsd, tau, jpvt, rank, qraux1)
    ! in/out
    integer, intent(in) :: m, n, mn, nrhs, lwork
    integer, intent(out) :: info, jpvt(n)
    integer, intent(inout) :: rank
    double precision, intent(in) :: tol
    double precision, intent(out) :: work(lwork), coef(n, nrhs), tau(*)
    double precision, intent(out), dimension(m, nrhs) :: ft, eff, rsd
    double precision, intent(inout) :: a(m, n), b(m, nrhs)
    double precision, intent(out) :: qraux1
    ! local
    integer i, j
    
    
    ! Assume model matrix is full rank
    if (rank == -1) then 
      rank = m
      call dgelqf(m, n, a, m, work(1), work(mn+1), lwork-mn, info)
      
      do i = 1, n
        jpvt(i) = i
      end do
      
    ! RRLQ
    else
!      call rdgeqp3(m, n, a, m, jpvt, work(1), work(mn+1), 3*n+1, tol, rank, info)
    end if
    
    call dormqr('left', 'transpose', m, nrhs, rank, a, m, work(1), b, m, work(mn+1), lwork-mn, info)
    
    ! Store "effects"
!    call dlacpy_omp('All', m, nrhs, b, m, eff, m)
    eff(1:m, 1:nrhs) = b(1:m, 1:nrhs)
    
    
    ! Store qraux(1) == work(1) 
    qraux1 = work(1)
    
    ! workspace at least nrhs, optimally nrhs*nb
    ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
    call dtrtrs('lower', 'no transpose', 'non-unit', rank, nrhs, a, m, b, m, info)
    
    
    !!! Produce fitted.values = Ax = Q*(R*x)
    
    ! Copy over the first RANK elements of numerical solution X
    !$omp parallel if(m*n > 5000) private(i, j) default(shared) 
    !$omp do
      do j = 1, nrhs
        do i = 1, rank
          ft(i, j) = b(i, j)
        end do
        
        do i = rank+1, m
          ft(i, j) = 0.0d0
        end do
      end do
    !$omp end do
    !$omp end parallel
    
    
    ! Compute fitted FT = Q*(R*fitted)
    call dtrmm('L', 'U', 'N', 'N', mn, nrhs, 1.0d0, a, m, ft, m)
    
    tau(1:mn) = work(1:mn)
    call dormlq('L', 'T', m, nrhs, mn, a, m, tau, ft, m, work, lwork, info)
    
    !!! FIXME dlacpy_omp doesn't copy beyond first column correctly ?!
    ! Compute residual RSD = B - FT
!    call dgeadd_omp('N', m, nrhs, -1.0d0, ft, m, 1.0d0, rsd, m)
    rsd(1:m, 1:nrhs) = rsd(1:m, 1:nrhs) - ft(1:m, 1:nrhs)
    
    ! Coefficients are stored in the first RANK elements of B
!    call dlacpy_omp('A', n, nrhs, b, 1, coef, 1)
    coef(1:n, 1:nrhs) = b(1:n, 1:nrhs)
    
    call rdgels_fixcoef(m, n, mn, nrhs, rank, jpvt, coef)
    
    return
  end subroutine
  
  
  
  subroutine rdgels_fixcoef(m, n, mn, nrhs, rank, jpvt, coef)
    ! in/out
    integer, intent(in) :: m, n, mn, nrhs, rank
    integer, intent(in) :: jpvt(n)
    double precision, intent(inout) :: coef(n, nrhs)
    ! local
    integer :: i, j, offset, tail
    integer, allocatable :: pvt(:)
    double precision :: tmpval
    double precision :: na_real
    
    
    call r_set_na(na_real)
    
    ! ------------------  m >= n  ------------------
!    if (m >= n) then
    if (1 > 0) then !!! FIXME TODO
      if (rank == n) return
      
      allocate(pvt(n))
      
      tail = mn
      do j = 1, nrhs
        pvt(1:n) = jpvt(1:n)
        
        ! fill back n-rank values with NA
        do i = 2, mn
          if (jpvt(i) < jpvt(i-1)) then
            tail = i
            exit
          end if
        end do
        
        if (m < n .and. tail == mn) tail = tail + 1
        
        if (tail /= n) then
          do i = tail, n
            coef(i, j) = na_real
          end do
        else if (rank == n-1) then
          coef(n, j) = na_real
        end if
        
        ! reorder
        call quicksort_by_index(coef(1:, j), pvt, n)
      end do
      
      deallocate(pvt)
    ! ------------------  m < n  ------------------
    end if
    
  end subroutine
end module

