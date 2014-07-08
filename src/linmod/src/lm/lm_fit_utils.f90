! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module lm_fit_utils
  use :: R_special
  use :: quicksorts, only : quicksort_by_index
  use :: lapack
  use :: lm, only : rdgeqp3
  implicit none
  
  contains
  
  subroutine rdgels_qr(m, n, mn, nrhs, a, lda, b, ldb, work, lwork, info, &
                    tol, coef, eff, ft, rsd, tau, jpvt, rank, qraux1)
    ! in/out
    integer, intent(in) :: m, n, mn, nrhs, lda, ldb, lwork
    integer, intent(out) :: info, jpvt(n)
    integer, intent(inout) :: rank
    double precision, intent(in) :: tol
    double precision, intent(out) :: work(*), coef(n, *), tau(*)
    double precision, intent(out), dimension(ldb, *) :: ft, eff, rsd
    double precision, intent(inout) :: a(lda, *), b(ldb, *)
    double precision, intent(out) :: qraux1
    ! local
    integer i, j
    
    
    ! Assume model matrix is full rank
    if (rank == -1) then 
      rank = n
      call dgeqrf(m, n, a, lda, work(1), work(mn+1), lwork-mn, info)
    ! RRQR
    else
  !      call rdgeqpf(m, n, a, lda, jpvt, work(1), work(mn+1), tol, rank, info)
      call rdgeqp3(m, n, a, lda, jpvt, work(1), work(mn+1), 3*n+1, tol, rank, info)
    end if
    
    
    ! least-squares problem min || a * x - b ||
    !  b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
    call dormqr('left', 'transpose', m, nrhs, rank, a, lda, work(1), b, ldb, work(mn+1), lwork-mn, info)
    
    ! Store "effects"
    call dlacpy_omp('All', m, nrhs, b, ldb, eff, ldb)
    
    
    ! Store qraux(1) == work(1) 
    qraux1 = work(1)
    
    ! workspace at least nrhs, optimally nrhs*nb
    ! b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
    call dtrtrs('upper', 'no transpose', 'non-unit', rank, nrhs, a, lda, b, ldb, info)
    
    
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
    call dtrmm('L', 'U', 'N', 'N', n, nrhs, 1.0d0, a, lda, ft, ldb)
    
    tau(1:min(m, n)) = work(1:min(m,n))
    call dormqr('L', 'N', m, nrhs, n, a, lda, tau, ft, ldb, work, lwork, info)
    
    ! Compute residual RSD = B - FT
    call dgeadd_omp('N', m, nrhs, -1.0d0, ft, ldb, 1.0d0, rsd, ldb)
    
    ! Coefficients are stored in the first RANK elements of B
!    call dlacpy_omp('A', n, nrhs, b, 1, coef, 1) !!! FIXME doesn't copy beyond first column correctly ?!
    coef(1:n, 1:nrhs) = b(1:n, 1:nrhs)
    
    call rdgels_fixcoef(m, n, mn, nrhs, rank, jpvt, coef)
  end subroutine
  
  
  
  
  subroutine rdgels_lq()
    
  end subroutine
  
  
  
  subroutine rdgels_fixcoef(m, n, mn, nrhs, rank, jpvt, coef)
    ! in/out
    integer, intent(in) :: m, n, mn, nrhs, rank
    integer, intent(in) :: jpvt(n)
    double precision, intent(inout) :: coef(n, *)
    ! local
    integer :: i, j, offset, tail
    integer, allocatable :: pvt(:)
    double precision :: tmpval
    double precision :: na_real
    
    
    call r_set_na(na_real)
    
    if (m >= n) then
      if (rank == n) return
      
      allocate(pvt(n))
      
      
      tail = n
      do j = 1, nrhs
        pvt(1:n) = jpvt(1:n)
        
        ! fill back n-rank values with NA
        do i = 2, n
          if (jpvt(i) < jpvt(i-1)) then
            tail = i
            exit
          end if
        end do
        
        if (tail /= n) then
          do i = tail, n
            coef(i, j) = na_real
          end do
        else if (rank == n-1) then
          coef(n, j) = na_real
        end if
        
        ! reorder
        call dquicksort_by_index(coef(1, j), pvt, n) !!! FIXME
        
        deallocate(pvt)
        
      end do
    end if
    
  end subroutine
end module

