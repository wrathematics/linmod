! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


!  -- lapack auxiliary routine (version 3.1) --
!     univ. of tennessee, univ. of california berkeley and nag ltd..
!     november 2006


!  purpose
!  =======
!  dlaqp2 computes a qr factorization with column pivoting of
!  the block a(offset+1:m,1:n).
!  the block a(1:offset,1:n) is accordingly pivoted, but not factorized.
!  arguments
!  =========
!  m       (input) integer
!          the number of rows of the matrix a. m >= 0.
!  n       (input) integer
!          the number of columns of the matrix a. n >= 0.
!  offset  (input) integer
!          the number of rows of the matrix a that must be pivoted
!          but no factorized. offset >= 0.
!  a       (input/output) double precision array, dimension (lda,n)
!          on entry, the m-by-n matrix a.
!          on exit, the upper triangle of block a(offset+1:m,1:n) is 
!          the triangular factor obtained; the elements in block
!          a(offset+1:m,1:n) below the diagonal, together with the
!          array tau, represent the orthogonal matrix q as a product of
!          elementary reflectors. block a(1:offset,1:n) has been
!          accordingly pivoted, but no factorized.
!  lda     (input) integer
!          the leading dimension of the array a. lda >= max(1,m).
!  jpvt    (input/output) integer array, dimension (n)
!          on entry, if jpvt(i)  /=  0, the i-th column of a is permuted
!          to the front of a*p (a leading column); if jpvt(i) = 0,
!          the i-th column of a is a free column.
!          on exit, if jpvt(i) = k, then the i-th column of a*p
!          was the k-th column of a.
!  tau     (output) double precision array, dimension (min(m,n))
!          the scalar factors of the elementary reflectors.
!  vn1     (input/output) double precision array, dimension (n)
!          the vector with the partial column norms.
!  vn2     (input/output) double precision array, dimension (n)
!          the vector with the exact column norms.
!  work    (workspace) double precision array, dimension (n)
!  further details
!  ===============
!  based on contributions by
!    g. quintana-orti, depto. de informatica, universidad jaime i, spain
!    x. sun, computer science dept., duke university, usa
!  partial column norm updating strategy modified by
!    z. drmac and z. bujanovic, dept. of mathematics,
!    university of zagreb, croatia.
!    june 2006.
!  for more details see lapack working note 176.
!  =====================================================================



subroutine rdlaqp2(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work, tol, rank)
  implicit none
  
  !in/out
  integer, intent(in) :: lda, m, n, offset
  integer :: jpvt(*)
  integer, intent(out) :: rank
  double precision, intent(in) :: tol
  double precision :: a(lda, *), tau(*), vn1(*), vn2(*), work(*)
  ! local
  integer :: i, itemp, j, mn, offpi, pvt
  double precision :: aii, temp, temp2, tol3z
  ! subroutines
  external           dlarf, dlarfg, dswap
  ! functions
  intrinsic          abs, max, min, sqrt
  integer            idamax
  double precision   dlamch, dnrm2
  external           idamax, dlamch, dnrm2
  
  
  !     .. executable statements ..
  mn = min(m-offset, n)
  tol3z = sqrt(dlamch('epsilon'))
  
  rank = n
  
!     compute factorization.
  do i = 1, mn
     offpi = offset + i
    
!        determine ith pivot column and swap if necessary.
!     pvt = (i-1) + idamax(n-i+1, vn1(i), 1)
    
!    if(pvt /= i) then
    do j = i, mn
      if (vn1(j) > tol) then
        pvt = j
        exit
      end if
    end do
    
!    if (vn1(i) < tol) then
    if(pvt /= i) then
      call dswap(m, a(1, rank), 1, a(1, i), 1)
      itemp = jpvt(rank)
      jpvt(rank) = jpvt(i)
      jpvt(i) = itemp
      vn1(rank) = vn1(i)
      vn2(rank) = vn2(i)
      
      rank = rank - 1
    end if
    
    
    
    
    
!        generate elementary reflector h(i).
    if(offpi < m) then
      call dlarfg(m-offpi+1, a(offpi, i), a(offpi+1, i), 1, tau(i))
    else
      call dlarfg(1, a(m, i), a(m, i), 1, tau(i))
    end if
    
    if(i < n) then
      ! apply h(i)' to a(offset+i:m,i+1:n) from the left.
      aii = a(offpi, i)
      a(offpi, i) = 1.0d0
      call dlarf('left', m-offpi+1, n-i, a(offpi, i), 1, tau(i), a(offpi, i+1), lda, work(1))
      a(offpi, i) = aii
    end if
    
    ! update partial column norms.
    
    do j = i + 1, n
      if(vn1(j) /= 0.0d0) then
        ! note: the following 4 lines follow from the analysis in
        ! lapack working note 176.
        temp = 1.0d0 - (abs(a(offpi, j)) / vn1(j))**2
        temp = max(temp, 0.0d0)
        temp2 = temp*(vn1(j) / vn2(j))**2
        
        if(temp2 <= tol3z) then
          if(offpi < m) then
            vn1(j) = dnrm2(m-offpi, a(offpi+1, j), 1)
            vn2(j) = vn1(j)
          else
            vn1(j) = 0.0d0
            vn2(j) = 0.0d0
          end if
        else
          vn1(j) = vn1(j)*sqrt(temp)
        end if
        
      end if
    end do
  end do
  
  print *, rank
  
  return
end


