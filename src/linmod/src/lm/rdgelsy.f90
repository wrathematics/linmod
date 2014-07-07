subroutine rdgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
  implicit none
  integer, intent(in) :: m, n, nrhs, lda, ldb
  integer, intent(inout) :: lwork
  integer, intent(out) :: info, rank
  integer, intent(out) :: jpvt(*)
  double precision, intent(out) :: rcond
  double precision, intent(out), dimension(lda, *) :: a, b
  double precision, intent(out) :: work(*)
!     .. parameters ..
  integer            imax, imin
  parameter          (imax = 1, imin = 2)
  double precision   zero, one
  parameter          (zero = 0.0d+0, one = 1.0d+0)
!     ..
!     .. local scalars ..
  logical            lquery
  integer            i, iascl, ibscl, ismax, ismin, j, lwkmin, &
                     lwkopt, mn, nb, nb1, nb2, nb3, nb4
  double precision   anrm, bignum, bnrm, c1, c2, s1, s2, smax, &
                     smaxpr, smin, sminpr, smlnum, wsize
!     ..
!     .. external functions ..
  integer            ilaenv
  double precision   dlamch, dlange
  external           ilaenv, dlamch, dlange
!     ..
!     .. external subroutines ..
  external           dcopy, dgeqp3, dlabad, dlaic1, dlascl, dlaset, &
                     dormqr, dormrz, dtrsm, dtzrzf, xerbla
!     ..
!     .. intrinsic functions ..
  intrinsic          abs, max, min
!     ..
!     .. executable statements ..
!
  mn = min(m, n)
  ismin = mn + 1
  ismax = 2*mn + 1
!
!     test the input arguments.
!
  info = 0
  lquery = (lwork == -1)
  if (m < 0) then
    info = -1
  else if (n < 0) then
    info = -2
  else if (nrhs < 0) then
    info = -3
  else if (lda < max(1, m)) then
    info = -5
  else if (ldb < max(1, m, n)) then
    info = -7
  end if
!
!     figure out optimal block size
!
  if (info == 0) then
    if (mn == 0 .or. nrhs == 0) then
      lwkmin = 1
      lwkopt = 1
    else
      nb1 = ilaenv(1, 'dgeqrf', ' ', m, n, -1, -1)
      nb2 = ilaenv(1, 'dgerqf', ' ', m, n, -1, -1)
      nb3 = ilaenv(1, 'dormqr', ' ', m, n, nrhs, -1)
      nb4 = ilaenv(1, 'dormrq', ' ', m, n, nrhs, -1)
      nb = max(nb1, nb2, nb3, nb4)
      lwkmin = mn + max(2*mn, n + 1, mn + nrhs)
      lwkopt = max(lwkmin, mn + 2*n + nb*(n + 1), 2*mn + nb*nrhs)
    end if
    work(1) = lwkopt
!
    if (lwork < lwkmin .and. .not.lquery) then
      info = -12
    end if
  end if
!
  if (info /= 0) then
    call xerbla('dgelsy', -info)
    return
  else if (lquery) then
    return
  end if
!
!     quick return if possible
!
  if (mn == 0 .or. nrhs == 0) then
    rank = 0
    return
  end if
!
!     get machine parameters
!
  smlnum = dlamch('s') / dlamch('p')
  bignum = one / smlnum
  call dlabad(smlnum, bignum)
!
!     scale a, b if max entries outside range [smlnum,bignum]
!
  anrm = dlange('m', m, n, a, lda, work)
  iascl = 0
  if (anrm > zero .and. anrm < smlnum) then
!
!        scale matrix norm up to smlnum
!
    call dlascl('g', 0, 0, anrm, smlnum, m, n, a, lda, info)
    iascl = 1
  else if (anrm > bignum) then
!
!        scale matrix norm down to bignum
!
    call dlascl('g', 0, 0, anrm, bignum, m, n, a, lda, info)
    iascl = 2
  else if (anrm == zero) then
!
!        matrix all zero. return zero solution.
!
    call dlaset('f', max(m, n), nrhs, zero, zero, b, ldb)
    rank = 0
    go to 70
  end if
!
  bnrm = dlange('m', m, nrhs, b, ldb, work)
  ibscl = 0
  if (bnrm > zero .and. bnrm < smlnum) then
!
!        scale matrix norm up to smlnum
!
    call dlascl('g', 0, 0, bnrm, smlnum, m, nrhs, b, ldb, info)
    ibscl = 1
  else if (bnrm > bignum) then
!
!        scale matrix norm down to bignum
!
    call dlascl('g', 0, 0, bnrm, bignum, m, nrhs, b, ldb, info)
    ibscl = 2
  end if
!
!     compute qr factorization with column pivoting of a:
!        a * p = q * r
!
  call dgeqp3(m, n, a, lda, jpvt, work(1), work(mn+1), lwork-mn, info)
  wsize = mn + work(mn+1)
!
!     workspace: mn+2*n+nb*(n+1).
!     details of householder rotations stored in work(1:mn).
!
!     determine rank using incremental condition estimation
!
  work(ismin) = one
  work(ismax) = one
  smax = abs(a(1, 1))
  smin = smax
  if (abs(a(1, 1)) == zero) then
    rank = 0
    call dlaset('f', max(m, n), nrhs, zero, zero, b, ldb)
    go to 70
  else
    rank = 1
  end if
  
  do
    if (rank < mn) then
      i = rank + 1
      call dlaic1(imin, rank, work(ismin), smin, a(1, i), a(i, i), sminpr, s1, c1)
      call dlaic1(imax, rank, work(ismax), smax, a(1, i), a(i, i), smaxpr, s2, c2)
      
      if (smaxpr*rcond <= sminpr) then
        do i = 1, rank
          work(ismin+i-1) = s1*work(ismin+i-1)
          work(ismax+i-1) = s2*work(ismax+i-1)
        end do
        
        work(ismin+rank) = c1
        work(ismax+rank) = c2
        smin = sminpr
        smax = smaxpr
        rank = rank + 1
      else
        exit
      end if
    end if
  end do
  
!     workspace: 3*mn.
  
!     logically partition r = [ r11 r12 ]
!                             [  0  r22 ]
!     where r11 = r(1:rank,1:rank)
  
!     [r11,r12] = [ t11, 0 ] * y
  
  if (rank < n) then
    call dtzrzf(rank, n, a, lda, work(mn+1), work(2*mn+1), lwork-2*mn, info)
  end if
  
!     workspace: 2*mn.
!     details of householder rotations stored in work(mn+1:2*mn)
  
!     b(1:m,1:nrhs) := q' * b(1:m,1:nrhs)
  
  call dormqr('left', 'transpose', m, nrhs, mn, a, lda, work(1), &
              b, ldb, work(2*mn+1), lwork-2*mn, info)
  wsize = max(wsize, 2*mn+work(2*mn+1))
  
!     workspace: 2*mn+nb*nrhs.
!
!     b(1:rank,1:nrhs) := inv(t11) * b(1:rank,1:nrhs)
  
  call dtrsm('left', 'upper', 'no transpose', 'non-unit', rank, nrhs, one, a, lda, b, ldb)
  
  do j = 1, nrhs
    do i = rank + 1, n
      b(i, j) = zero
    end do
  end do
  
!     b(1:n,1:nrhs) := y' * b(1:n,1:nrhs)
  
  if (rank < n) then
     call dormrz('left', 'transpose', n, nrhs, rank, n-rank, a, &
          lda, work(mn+1), b, ldb, work(2*mn+1), lwork-2*mn, info)
  end if
  
!     workspace: 2*mn+nrhs.
!
!     b(1:n,1:nrhs) := p * b(1:n,1:nrhs)
  
  do j = 1, nrhs
    do i = 1, n
      work(jpvt(i)) = b(i, j)
    end do
    call dcopy(n, work(1), 1, b(1, j), 1)
  end do
  
!     workspace: n.
!
!     undo scaling
  
  if (iascl == 1) then
     call dlascl('g', 0, 0, anrm, smlnum, n, nrhs, b, ldb, info)
     call dlascl('u', 0, 0, smlnum, anrm, rank, rank, a, lda, info)
  else if (iascl == 2) then
     call dlascl('g', 0, 0, anrm, bignum, n, nrhs, b, ldb, info)
     call dlascl('u', 0, 0, bignum, anrm, rank, rank, a, lda, info)
  end if
  if (ibscl == 1) then
     call dlascl('g', 0, 0, smlnum, bnrm, n, nrhs, b, ldb, info)
  else if (ibscl == 2) then
     call dlascl('g', 0, 0, bignum, bnrm, n, nrhs, b, ldb, info)
  end if
  
  
70 continue
  work(1) = lwkopt
  
  return
end subroutine


