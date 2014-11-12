! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


!  -- lapack auxiliary routine (version 3.1) --
!     univ. of tennessee, univ. of california berkeley and nag ltd..
!     november 2006


!  dlaqp2 computes a qr factorization with column pivoting of
!  the block a(offset+1:m,1:n).
!  the block a(1:offset,1:n) is accordingly pivoted, but not factorized.
!  arguments


subroutine rdlaqp2(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work, tol, rank)
  use :: lapack
  implicit none
  
  !in/out
  integer, intent(in) :: lda, m, n, offset
  integer, intent(out) :: jpvt(*)
  integer, intent(out) :: rank
  double precision, intent(in) :: tol
  double precision, intent(inout) :: a(lda, *), tau(*), vn1(*), vn2(*), work(*)
  ! local
  integer :: i, itemp, j, mn, offpi, pvt
  double precision :: aii, temp, temp2, tol3z
  ! functions
  intrinsic          abs, max, min, sqrt
  
  
  mn = min(m-offset, n)
!  tol3z = dsqrt(dlamch('epsilon'))
  
  ! compute factorization.
  do i = 1, mn
     offpi = offset + i
    
  ! determine ith pivot column and swap if necessary.
    if (i < mn) then
      do j = i, mn
        if (vn1(j) > tol) then
          pvt = j
          exit
        end if
      end do
    else
      pvt = mn
    end if
    
    ! if (vn1(i) < tol) then
    if (pvt /= i) then
      call dswap(m, a(1, pvt), 1, a(1, i), 1)
      
      itemp = jpvt(pvt)
      jpvt(pvt) = jpvt(i)
      jpvt(i) = itemp
      
      vn1(pvt) = vn1(i)
      vn2(pvt) = vn2(i)
    end if
    
    
    ! generate elementary reflector h(i).
    if (offpi < m) then
      call dlarfg(m-offpi+1, a(offpi, i), a(offpi+1, i), 1, tau(i))
    else
      call dlarfg(1, a(m, i), a(m, i), 1, tau(i))
    end if
    
    
    ! apply h(i)' to a(offset+i:m,i+1:n) from the left.
    if (i < n) then
      aii = a(offpi, i)
      a(offpi, i) = 1.0d0
      call dlarf('left', m-offpi+1, n-i, a(offpi, i), 1, tau(i), a(offpi, i+1), lda, work(1))
      a(offpi, i) = aii
    end if
    
    
    ! update partial column norms.
    do j = i + 1, n
      if (vn1(j) /= 0.0d0) then
        ! note: the following 4 lines follow from the analysis in
        ! lapack working note 176.
        temp = 1.0d0 - (dabs(a(offpi, j)) / vn1(j))**2
        temp = max(temp, 0.0d0)
        temp2 = temp*(vn1(j) / vn2(j))**2
        
        if (temp2 <= tol) then
          if (offpi < m) then
            vn1(j) = dnrm2(m-offpi, a(offpi+1, j), 1)
            vn2(j) = vn1(j)
          else
            vn1(j) = 0.0d0
            vn2(j) = 0.0d0
          end if
        else
          vn1(j) = vn1(j)*dsqrt(temp)
        end if
        
      end if
    end do
  end do
  
  
  ! Estimate numerical rank
  rank = 0
  
  do i = 1, mn-1
    if (jpvt(i) < i) then
      exit
    else
      rank = rank + 1
    end if
  end do
  
  if (jpvt(n) == n .and. vn1(n) > tol) rank = rank + 1
  
  if (rank == 0) rank = 1
  
  return
end subroutine


