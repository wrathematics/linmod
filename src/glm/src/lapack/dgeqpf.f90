!  -- lapack deprecated driver routine (version 3.1) --
!     univ. of tennessee, univ. of california berkeley and nag ltd..
!     november 2006


subroutine rdgeqpf( m, n, a, lda, jpvt, tau, work, info )
!     .. scalar arguments ..
  integer            info, lda, m, n
!     ..
!     .. array arguments ..
  integer            jpvt( * )
  double precision   a( lda, * ), tau( * ), work( * )
!     ..
!     .. parameters ..
  double precision   zero, one
  parameter          ( zero = 0.0d+0, one = 1.0d+0 )
!     ..
!     .. local scalars ..
  integer            i, itemp, j, ma, mn, pvt
  double precision   aii, temp, temp2, tol3z
!     ..
!     .. external subroutines ..
  external           dgeqr2, dlarf, dlarfg, dorm2r, dswap, xerbla
!     ..
!     .. intrinsic functions ..
  intrinsic          abs, max, min, sqrt
!     ..
!     .. external functions ..
  integer            idamax
  double precision   dlamch, dnrm2
  external           idamax, dlamch, dnrm2
!     ..
!     .. executable statements ..
!
!     test the input arguments
!
  info = 0
  if( m < 0 ) then
     info = -1
  else if( n < 0 ) then
     info = -2
  else if( lda < max( 1, m ) ) then
     info = -4
  end if
  if( info /= 0 ) then
     call xerbla( 'dgeqpf', -info )
     return
  end if
!
  mn = min( m, n )
  tol3z = sqrt(dlamch('epsilon'))
!
!     move initial columns up front
!
  itemp = 1
  do i = 1, n
     if( jpvt( i ) /= 0 ) then
        if( i /= itemp ) then
           call dswap( m, a( 1, i ), 1, a( 1, itemp ), 1 )
           jpvt( i ) = jpvt( itemp )
           jpvt( itemp ) = i
        else
           jpvt( i ) = i
        end if
        itemp = itemp + 1
     else
        jpvt( i ) = i
     end if
  end do
  itemp = itemp - 1
!
!     compute the qr factorization and update remaining columns
!
  if( itemp > 0 ) then
     ma = min( itemp, m )
     call dgeqr2( m, ma, a, lda, tau, work, info )
     if( ma < n ) then
        call dorm2r( 'left', 'transpose', m, n-ma, ma, a, lda, tau, a( 1, ma+1 ), lda, work, info )
     end if
  end if
!
  if( itemp < mn ) then
!
!        initialize partial column norms. the first n elements of
!        work store the exact column norms.
!
     do i = itemp + 1, n
        work( i ) = dnrm2( m-itemp, a( itemp+1, i ), 1 )
        work( n+i ) = work( i )
     end do
!
!        compute factorization
!
     do i = itemp + 1, mn
!
!           determine ith pivot column and swap if necessary
!
        pvt = ( i-1 ) + idamax( n-i+1, work( i ), 1 )
!
        if( pvt /= i ) then
           call dswap( m, a( 1, pvt ), 1, a( 1, i ), 1 )
           itemp = jpvt( pvt )
           jpvt( pvt ) = jpvt( i )
           jpvt( i ) = itemp
           work( pvt ) = work( i )
           work( n+pvt ) = work( n+i )
        end if
!
!           generate elementary reflector h(i)
!
        if( i < m ) then
           call dlarfg( m-i+1, a( i, i ), a( i+1, i ), 1, tau( i ) )
        else
           call dlarfg( 1, a( m, m ), a( m, m ), 1, tau( m ) )
        end if
!
        if( i < n ) then
!
!              apply h(i) to a(i:m,i+1:n) from the left
!
           aii = a( i, i )
           a( i, i ) = one
           call dlarf( 'left', m-i+1, n-i, a( i, i ), 1, tau( i ), a( i, i+1 ), lda, work( 2*n+1 ) )
           a( i, i ) = aii
        end if
!
!           update partial column norms
!
        do j = i + 1, n
           if( work( j ) /= zero ) then
!
!                 note: the following 4 lines follow from the analysis in
!                 lapack working note 176.
!                 
              temp = abs( a( i, j ) ) / work( j )
              temp = max( zero, ( one+temp )*( one-temp ) )
              temp2 = temp*( work( j ) / work( n+j ) )**2
              if( temp2 .le. tol3z ) then 
                 if( m-i > 0 ) then
                    work( j ) = dnrm2( m-i, a( i+1, j ), 1 )
                    work( n+j ) = work( j )
                 else
                    work( j ) = zero
                    work( n+j ) = zero
                 end if
              else
                 work( j ) = work( j )*sqrt( temp )
              end if
           end if
        end do
    end do
  end if
  return
!
!     end of dgeqpf
!
end

