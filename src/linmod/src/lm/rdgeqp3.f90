!  purpose
!  =======
  
!  dgeqp3 computes a qr factorization with column pivoting of a
!  matrix a:  a*p = q*r  using level 3 blas.
  
!  arguments
!  =========
  
!  m      (input) integer
  !     the number of rows of the matrix a. m >= 0.
  
!  n      (input) integer
  !     the number of columns of the matrix a.  n >= 0.
  
!  a      (input/output) double precision array, dimension (lda,n)
  !     on entry, the m-by-n matrix a.
  !     on exit, the upper triangle of the array contains the
  !     min(m,n)-by-n upper trapezoidal matrix r; the elements below
  !     the diagonal, together with the array tau, represent the
  !     orthogonal matrix q as a product of min(m,n) elementary
  !     reflectors.
  
!  lda    (input) integer
  !     the leading dimension of the array a. lda >= max(1,m).
  
!  jpvt    (input/output) integer array, dimension (n)
  !     on entry, if jpvt(j) /= 0, the j-th column of a is permuted
  !     to the front of a*p (a leading column); if jpvt(j)=0,
  !     the j-th column of a is a free column.
  !     on exit, if jpvt(j)=k, then the j-th column of a*p was the
  !     the k-th column of a.
  
!  tau    (output) double precision array, dimension (min(m,n))
  !     the scalar factors of the elementary reflectors.
  
!  work    (workspace/output) double precision array, dimension (max(1,lwork))
  !     on exit, if info=0, work(1) returns the optimal lwork.
  
!  lwork   (input) integer
  !     the dimension of the array work. lwork >= 3*n+1.
  !     for optimal performance lwork >= 2*n+(n+1)*nb, where nb
  !     is the optimal blocksize.
  
  !     if lwork = -1, then a workspace query is assumed; the routine
  !     only calculates the optimal size of the work array, returns
  !     this value as the first entry of the work array, and no error
  !     message related to lwork is issued by xerbla.
  
!  info    (output) integer
  !     = 0: successful exit.
  !     < 0: if info = -i, the i-th argument had an illegal value.
  



subroutine rdgeqp3(m, n, a, lda, jpvt, tau, work, lwork, tol, rank, info)
  
!  -- lapack routine (version 3.3.1) --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!  -- april 2011                                            --
  
  integer, intent(out) :: rank
  double precision, intent(in) :: tol
  
  
  ! .. scalar arguments ..
  integer          info, lda, lwork, m, n
  ! ..
  ! .. array arguments ..
  integer          jpvt(*)
  double precision   a(lda, *), tau(*), work(*)
  ! ..
  ! .. parameters ..
  integer          inb, inbmin, ixover
  parameter        (inb = 1, inbmin = 2, ixover = 3)
  ! ..
  ! .. local scalars ..
  logical          lquery
  integer          fjb, iws, j, jb, lwkopt, minmn, minws, na, nb, &
                 nbmin, nfxd, nx, sm, sminmn, sn, topbmn
  ! ..
  ! .. external subroutines ..
  external         dgeqrf, dlaqp2, dlaqps, dormqr, dswap, xerbla
  ! ..
  ! .. external functions ..
  integer          ilaenv
  double precision   dnrm2
  external         ilaenv, dnrm2
  ! ..
  ! .. intrinsic functions ..
  intrinsic        int, max, min
  ! ..
  ! .. executable statements ..
  
  
  ! test input arguments
  info = 0
  lquery = (lwork == -1)
  if(m < 0) then
    info = -1
  else if(n < 0) then
    info = -2
  else if(lda < max(1, m)) then
    info = -4
  end if
  
  if(info == 0) then
    minmn = min(m, n)
    if(minmn == 0) then
      iws = 1
      lwkopt = 1
    else
      iws = 3*n + 1
      nb = ilaenv(inb, 'dgeqrf', ' ', m, n, -1, -1)
      lwkopt = 2*n + (n + 1)*nb
    end if
    work(1) = lwkopt
    
    if((lwork < iws) .and. .not.lquery) then
      info = -8
    end if
  end if
  
  if(info /= 0) then
    call xerbla('dgeqp3', -info)
    return
  else if(lquery) then
    return
  end if
  
  ! quick return if possible.
  
  if(minmn == 0) then
    return
  end if
  
  ! move initial columns up front.
  
  
  nfxd = 1
  do j = 1, n
    if(jpvt(j) /= 0) then
      if(j /= nfxd) then
        call dswap(m, a(1, j), 1, a(1, nfxd), 1)
        jpvt(j) = jpvt(nfxd)
        jpvt(nfxd) = j
      else
        jpvt(j) = j
      end if
      nfxd = nfxd + 1
    else
      jpvt(j) = j
    end if
  end do
  nfxd = nfxd - 1
  
  
  ! factorize fixed columns
  ! =======================
  
  ! compute the qr factorization of fixed columns and update
  ! remaining columns.
  
  if(nfxd > 0) then
    na = min(m, nfxd)
!*cc     call dgeqr2(m, na, a, lda, tau, work, info)
    call dgeqrf(m, na, a, lda, tau, work, lwork, info)
    iws = max(iws, int(work(1)))
    if(na < n) then
!*cc        call dorm2r('left', 'transpose', m, n-na, na, a, lda,
!*cc  $                tau, a(1, na+1), lda, work, info)
       call dormqr('left', 'transpose', m, n-na, na, a, lda, tau, a(1, na+1), lda, work, lwork, info)
       iws = max(iws, int(work(1)))
    end if
  end if
  
  ! factorize free columns
  ! ======================
  
  if(nfxd < minmn) then
  
    sm = m - nfxd
    sn = n - nfxd
    sminmn = minmn - nfxd
  
    ! determine the block size.
    nb = ilaenv(inb, 'dgeqrf', ' ', sm, sn, -1, -1)
    nbmin = 2
    nx = 0
  
    if((nb > 1) .and. (nb < sminmn)) then
      ! determine when to cross over from blocked to unblocked code.
      nx = max(0, ilaenv(ixover, 'dgeqrf', ' ', sm, sn, -1, -1))
      
      if(nx < sminmn) then
      ! determine if workspace is large enough for blocked code.
        
        minws = 2*sn + (sn+1)*nb
        iws = max(iws, minws)
        if(lwork < minws) then
          ! not enough workspace to use optimal nb: reduce nb and
          ! determine the minimum value of nb.
          nb = (lwork-2*sn) / (sn+1)
          nbmin = max(2, ilaenv(inbmin, 'dgeqrf', ' ', sm, sn, -1, -1))
        end if
      end if
    end if
  
    ! initialize partial column norms. the first n elements of work
    ! store the exact column norms.
    do j = nfxd + 1, n
      work(j) = dnrm2(sm, a(nfxd+1, j), 1)
      work(n+j) = work(j)
    end do
  
    if((nb >= nbmin) .and. (nb < sminmn) .and. (nx < sminmn)) then
      ! use blocked code initially.
      j = nfxd + 1
      
      ! compute factorization: while loop.
      topbmn = minmn - nx
      do while(j <= topbmn)
         jb = min(nb, topbmn-j+1)
  
        ! factorize jb columns among columns j:n.
  
         call dlaqps(m, n-j+1, j-1, jb, fjb, a(1, j), lda, &
                    jpvt(j), tau(j), work(j), work(n+j), &
                    work(2*n+1), work(2*n+jb+1), n-j+1)
  
         j = j + fjb
      end do
    else
       j = nfxd + 1
    end if
  
    ! use unblocked code to factor the last or only block.
    if(j <= minmn) then
      call rdlaqp2(m, n-j+1, j-1, a(1, j), lda, jpvt(j), tau(j), work(j), work(n+j), work(2*n+1), tol, rank)
    end if
  
  end if
  
  work(1) = iws
  
  
!!!  !!! Estimate numerical rank
!!!  rank = 0
!!!  
!!!  do i = 1, n
!!!    if (i /= jpvt(i)) then
!!!     rank = rank + 1
!!!    end if
!!!  end do
  
  return
end
