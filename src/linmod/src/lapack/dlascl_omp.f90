! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


! OMP statements added to dlascl, original copyright:
!  -- LAPACK auxiliary routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006

subroutine dlascl( type, kl, ku, cfrom, cto, m, n, a, lda, info )
!     .. scalar arguments ..
      character          type
      integer            info, kl, ku, lda, m, n
      double precision   cfrom, cto
!     ..
!     .. array arguments ..
      double precision   a( lda, * )
!     ..
!  purpose
!  =======
!  dlascl multiplies the m by n real matrix a by the real scalar
!  cto/cfrom.  this is done without over/underflow as long as the final
!  result cto*a(i,j)/cfrom does not over/underflow. type specifies that
!  a may be full, upper triangular, lower triangular, upper hessenberg,
!  or banded.
!  arguments
!  =========
!  type    (input) character*1
!          type indices the storage type of the input matrix.
!          = 'g':  a is a full matrix.
!          = 'l':  a is a lower triangular matrix.
!          = 'u':  a is an upper triangular matrix.
!          = 'h':  a is an upper hessenberg matrix.
!          = 'b':  a is a symmetric band matrix with lower bandwidth kl
!                  and upper bandwidth ku and with the only the lower
!                  half stored.
!          = 'q':  a is a symmetric band matrix with lower bandwidth kl
!                  and upper bandwidth ku and with the only the upper
!                  half stored.
!          = 'z':  a is a band matrix with lower bandwidth kl and upper
!                  bandwidth ku.
!  kl      (input) integer
!          the lower bandwidth of a.  referenced only if type = 'b',
!          'q' or 'z'.
!  ku      (input) integer
!          the upper bandwidth of a.  referenced only if type = 'b',
!          'q' or 'z'.
!  cfrom   (input) double precision
!  cto     (input) double precision
!          the matrix a is multiplied by cto/cfrom. a(i,j) is computed
!          without over/underflow if the final result cto*a(i,j)/cfrom
!          can be represented without over/underflow.  cfrom must be
!          nonzero.
!  m       (input) integer
!          the number of rows of the matrix a.  m >= 0.
!  n       (input) integer
!          the number of columns of the matrix a.  n >= 0.
!  a       (input/output) double precision array, dimension (lda,n)
!          the matrix to be multiplied by cto/cfrom.  see type for the
!          storage type.
!  lda     (input) integer
!          the leading dimension of the array a.  lda >= max(1,m).
!  info    (output) integer
!          0  - successful exit
!          <0 - if info = -i, the i-th argument had an illegal value.
!  =====================================================================
!     .. parameters ..
  double precision   zero, one
  parameter          ( zero = 0.0d0, one = 1.0d0 )
!     ..
!     .. local scalars ..
  logical            done
  integer            i, itype, j, k1, k2, k3, k4
  double precision   bignum, cfrom1, cfromc, cto1, ctoc, mul, smlnum
!     ..
!     .. external functions ..
  logical            lsame
  double precision   dlamch
  external           lsame, dlamch
!     ..
!     .. intrinsic functions ..
  intrinsic          abs, max, min
!     ..
!     .. external subroutines ..
  external           xerbla
!     ..
!     .. executable statements ..
!     test the input arguments
  info = 0
  if( lsame( type, 'g' ) ) then
     itype = 0
  else if( lsame( type, 'l' ) ) then
     itype = 1
  else if( lsame( type, 'u' ) ) then
     itype = 2
  else if( lsame( type, 'h' ) ) then
     itype = 3
  else if( lsame( type, 'b' ) ) then
     itype = 4
  else if( lsame( type, 'q' ) ) then
     itype = 5
  else if( lsame( type, 'z' ) ) then
     itype = 6
  else
     itype = -1
  end if
*
  if( itype == -1 ) then
     info = -1
  else if( cfrom == zero ) then
     info = -4
  else if( m.lt.0 ) then
     info = -6
  else if( n.lt.0 .or. ( itype == 4 .and. n.ne.m ) .or.
 $         ( itype == 5 .and. n.ne.m ) ) then
     info = -7
  else if( itype.le.3 .and. lda.lt.max( 1, m ) ) then
     info = -9
  else if( itype.ge.4 ) then
     if( kl.lt.0 .or. kl.gt.max( m-1, 0 ) ) then
        info = -2
     else if( ku.lt.0 .or. ku.gt.max( n-1, 0 ) .or.
 $            ( ( itype == 4 .or. itype == 5 ) .and. kl.ne.ku ) )
 $             then
        info = -3
     else if( ( itype == 4 .and. lda.lt.kl+1 ) .or.
 $            ( itype == 5 .and. lda.lt.ku+1 ) .or.
 $            ( itype == 6 .and. lda.lt.2*kl+ku+1 ) ) then
        info = -9
     end if
  end if
*
  if( info.ne.0 ) then
     call xerbla( 'dlascl', -info )
     return
  end if
*
!     quick return if possible
*
  if( n == 0 .or. m == 0 )
 $   return
*
!     get machine parameters
*
  smlnum = dlamch( 's' )
  bignum = one / smlnum
*
  cfromc = cfrom
  ctoc = cto
*
1 continue
  cfrom1 = cfromc*smlnum
  cto1 = ctoc / bignum
  if( abs( cfrom1 ).gt.abs( ctoc ) .and. ctoc.ne.zero ) then
     mul = smlnum
     done = .false.
     cfromc = cfrom1
  else if( abs( cto1 ).gt.abs( cfromc ) ) then
     mul = bignum
     done = .false.
     ctoc = cto1
  else
     mul = ctoc / cfromc
     done = .true.
  end if
*
  if( itype == 0 ) then
*
!        full matrix
  *
  do j = 1, n
    do i = 1, m
      a( i, j ) = a( i, j )*mul
    end do
  end do
*
  else if( itype == 1 ) then
*
!        lower triangular matrix
*
  do j = 1, n
    do i = j, m
      a( i, j ) = a( i, j )*mul
    end do
  end do
*
  else if( itype == 2 ) then
*
!        upper triangular matrix
*
  do j = 1, n
    do i = 1, min( j, m )
      a( i, j ) = a( i, j )*mul
    end do
  end do
*
  else if( itype == 3 ) then
*
!        upper hessenberg matrix
*
  do j = 1, n
    do i = 1, min( j+1, m )
      a( i, j ) = a( i, j )*mul
    end do
  end do
*
  else if( itype == 4 ) then
    ! lower half of a symmetric band matrix
    k3 = kl + 1
    k4 = n + 1
    do 110 j = 1, n
      do 100 i = 1, min( k3, k4-j )
        a( i, j ) = a( i, j )*mul
      end do
    end do
    
  else if( itype == 5 ) then
    ! upper half of a symmetric band matrix
    k1 = ku + 2
    k3 = ku + 1
    do j = 1, n
      do i = max( k1-j, 1 ), k3
        a( i, j ) = a( i, j )*mul
      end do
    end do
    
  else if( itype == 6 ) then
    ! band matrix
    k1 = kl + ku + 2
    k2 = kl + 1
    k3 = 2*kl + ku + 1
    k4 = kl + ku + 1 + m
    do j = 1, n
      do i = max( k1-j, k2 ), min( k3, k4-j )
        a( i, j ) = a( i, j )*mul
      end do
    end do
    
  end if
  
  
  if( .not.done ) goto 1
*
  return
*
!     end of dlascl
*
end


