! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


! OMP statements added to dlacpy and dlascl, original copyright:
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006


module lapack_omp
  use :: omp_lib
  use :: linmod_omp
  use :: lapack, only : lsame, dlamch, xerbla
  implicit none
  
  
  contains
  
  ! b = a
  subroutine dlacpy_omp(uplo, m, n, a, lda, b, ldb)
    ! in/out
    character(len=1), intent(in) :: uplo
    integer, intent(in) :: lda, ldb, m, n
    double precision, intent(in) :: a(lda, *)
    double precision, intent(out) :: b(ldb, *)
    ! local
    integer :: i, j, itmp
    intrinsic          min
    
    
    !$omp parallel if (m*n > linmod_omp_minsize) private(i, j) default(shared) 
    if (uplo == 'u' .or. uplo == 'U') then
      !$omp do
        do j = 1, n
          do i = 1, min(j, m)
            b(i, j) = a(i, j)
          end do
        end do
      !$omp end do
    else if (uplo == 'l' .or. uplo == 'L') then
      !$omp do
        do j = 1, n
          do i = j, m
            b(i, j) = a(i, j)
          end do
        end do
      !$omp end do
    else
      !$omp do
        do j = 1, n
          do i = 1, m
            b(i, j) = a(i, j)
          end do
        end do
      !$omp end do
    end if
    !$omp end parallel
    
    return
  end
  
  
  
  ! c = alpha*a + beta*c
  subroutine dgeadd_omp(trans, m, n, alpha, a, lda, beta, c, ldc)
    ! in/out
    character(len=1), intent(in) :: trans
    integer, intent(in) :: m, n, lda, ldc
    double precision, intent(in) :: a(lda, *), alpha, beta
    double precision, intent(inout) :: c(ldc, *)
    ! local
    integer :: i, j
    integer :: ind, tind
    
    
    !$omp parallel if (m*n > linmod_omp_minsize) private(i, j, ind) default(shared) 
    if (trans == 'N' .or. trans == 'n') then
      !$omp do 
        do j = 1, n
          do i = 1, m
            c(i, j) = alpha*a(i, j) + beta*c(i, j)
          end do
        end do
      !$omp end do
      
    else if (trans == 'T' .or. trans == 't') then
      !$omp do
        do j = 1, n
          do i = 1, m
            c(i, j) = alpha*a(j, i) + beta*c(i, j)
          end do
        end do
      !$omp end do
    end if
    !$omp end parallel
    
    return
  end subroutine
  
  
  
  !  Multiplies the M by N real matrix A by the real scalar CTO/CFROM.
  subroutine dlascl_omp(type, kl, ku, cfrom, cto, m, n, a, lda, info)
    ! in/out
    character(len=1), intent(in) :: type
    integer, intent(in) :: kl, ku, lda, m, n
    integer, intent(out) :: info
    double precision, intent(in) :: cfrom, cto
    double precision, intent(inout) :: a(lda, *)
    ! local
    logical :: done
    integer :: i, itype, j, k1, k2, k3, k4
    double precision :: bignum, cfrom1, cfromc, cto1, ctoc, mul, smlnum
    ! intrinsic functions
    intrinsic          abs, max, min
    
    
    ! test the input arguments
    info = 0
    if (lsame(type, 'g')) then
      itype = 0
    else if (lsame(type, 'l')) then
      itype = 1
    else if (lsame(type, 'u')) then
      itype = 2
    else if (lsame(type, 'h')) then
      itype = 3
    else if (lsame(type, 'b')) then
      itype = 4
    else if (lsame(type, 'q')) then
      itype = 5
    else if (lsame(type, 'z')) then
      itype = 6
    else
      itype = -1
    end if
    
    if (itype == -1) then
      info = -1
    else if (cfrom == 0.0d0) then
      info = -4
    else if (m < 0) then
      info = -6
    else if (n < 0 .or. (itype == 4 .and. n /= m) .or. (itype == 5 .and. n /= m)) then
      info = -7
    else if (itype.le.3 .and. lda < max(1, m)) then
      info = -9
    else if (itype.ge.4) then
      if (kl < 0 .or. kl > max(m-1, 0)) then
        info = -2
      else if (ku < 0 .or. ku > max(n-1, 0) .or. &
               ((itype == 4 .or. itype == 5) .and. kl /= ku)) then
        info = -3
      else if ((itype == 4 .and. lda < kl+1) .or. &
               (itype == 5 .and. lda < ku+1) .or. &
               (itype == 6 .and. lda < 2*kl+ku+1)) then
        info = -9
      end if
    end if
    
    if (info /= 0) then
      call xerbla('dlascl', -info)
      return
    end if
    
    ! quick return if possible
    if (n == 0 .or. m == 0) return
    
    ! get machine parameters
    smlnum = dlamch('s')
    bignum = 1.0d0 / smlnum
    
    cfromc = cfrom
    ctoc = cto
    
    done = .false.
    !$omp parallel if (m*n > linmod_omp_minsize) private(i, j) default(shared) 
    do while (.not.done)
      cfrom1 = cfromc*smlnum
      cto1 = ctoc / bignum
      if (abs(cfrom1) > abs(ctoc) .and. ctoc /= 0.0d0) then
        mul = smlnum
        done = .false.
        cfromc = cfrom1
      else if (abs(cto1) > abs(cfromc)) then
        mul = bignum
        done = .false.
        ctoc = cto1
      else
        mul = ctoc / cfromc
        done = .true.
      end if
      
      ! full matrix
      if (itype == 0) then
        !$omp do
          do j = 1, n
            do i = 1, m
              a(i, j) = a(i, j)*mul
            end do
          end do
        !$omp end do
      
      ! lower triangular matrix
      else if (itype == 1) then
        !$omp do
          do j = 1, n
            do i = j, m
              a(i, j) = a(i, j)*mul
            end do
          end do
        !$omp end do
      
      ! upper triangular matrix
      else if (itype == 2) then
        !$omp do
          do j = 1, n
            do i = 1, min(j, m)
              a(i, j) = a(i, j)*mul
            end do
          end do
        !$omp end do
      
      ! upper hessenberg matrix
      else if (itype == 3) then
        !$omp do
          do j = 1, n
            do i = 1, min(j+1, m)
              a(i, j) = a(i, j)*mul
            end do
          end do
        !$omp end do
      
      ! lower half of a symmetric band matrix
      else if (itype == 4) then
        k3 = kl + 1
        k4 = n + 1
        !$omp do
          do j = 1, n
            do i = 1, min(k3, k4-j)
              a(i, j) = a(i, j)*mul
            end do
          end do
        !$omp end do
        
      ! upper half of a symmetric band matrix
      else if (itype == 5) then
        k1 = ku + 2
        k3 = ku + 1
        !$omp do
          do j = 1, n
            do i = max(k1-j, 1), k3
              a(i, j) = a(i, j)*mul
            end do
          end do
        !$omp end do
        
      ! band matrix
      else if (itype == 6) then
        k1 = kl + ku + 2
        k2 = kl + 1
        k3 = 2*kl + ku + 1
        k4 = kl + ku + 1 + m
        !$omp do
          do j = 1, n
            do i = max(k1-j, k2), min(k3, k4-j)
              a(i, j) = a(i, j)*mul
            end do
          end do
        !$omp end do
        
      end if
    
    end do
    !$omp end parallel
    
    
    return
  end subroutine
  
end module
