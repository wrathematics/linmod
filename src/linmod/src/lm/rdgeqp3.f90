! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt

! Modified from LAPACK routine dgeqp3, made to perform column
! pivoting as R does.  Original copyright:
!  -- lapack routine (version 3.3.1) --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!  -- april 2011                                            --


! computes a qr factorization with column pivoting of a using level 3 blas.
subroutine rdgeqp3(m, n, a, lda, jpvt, tau, work, lwork, tol, rank, info)
  use :: qr_utils, only : rdlaqp2
  use :: lapack
  
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
  ! .. intrinsic functions ..
  intrinsic        int, max, min
  
  
  ! test input arguments
  info = 0
  lquery = (lwork == -1)
  if (m < 0) then
    info = -1
  else if (n < 0) then
    info = -2
  else if (lda < max(1, m)) then
    info = -4
  end if
  
  if (info == 0) then
    minmn = min(m, n)
    if (minmn == 0) then
      iws = 1
      lwkopt = 1
    else
      iws = 3*n + 1
      nb = ilaenv(inb, 'dgeqrf', ' ', m, n, -1, -1)
      lwkopt = 2*n + (n + 1)*nb
    end if
    work(1) = lwkopt
    
    if ((lwork < iws) .and. .not.lquery) then
      info = -8
    end if
  end if
  
  if (info /= 0) then
    call xerbla('dgeqp3', -info)
    return
  else if (lquery) then
    return
  end if
  
  ! quick return if possible.
  if (minmn == 0) then
    return
  end if
  
  ! move initial columns up front.
  nfxd = 1
  do j = 1, n
    if (jpvt(j) /= 0) then
      if (j /= nfxd) then
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
  if (nfxd > 0) then
    na = min(m, nfxd)
    call dgeqrf(m, na, a, lda, tau, work, lwork, info)
    iws = max(iws, int(work(1)))
    if (na < n) then
       call dormqr('left', 'transpose', m, n-na, na, a, lda, tau, a(1, na+1), lda, work, lwork, info)
       iws = max(iws, int(work(1)))
    end if
  end if
  
  ! factorize free columns
  ! ======================
  if (nfxd < minmn) then
  
    sm = m - nfxd
    sn = n - nfxd
    sminmn = minmn - nfxd
  
    ! determine the block size.
    nb = ilaenv(inb, 'dgeqrf', ' ', sm, sn, -1, -1)
    nbmin = 2
    nx = 0
  
    if ((nb > 1) .and. (nb < sminmn)) then
      ! determine when to cross over from blocked to unblocked code.
      nx = max(0, ilaenv(ixover, 'dgeqrf', ' ', sm, sn, -1, -1))
      
      if (nx < sminmn) then
      ! determine if workspace is large enough for blocked code.
        minws = 2*sn + (sn+1)*nb
        iws = max(iws, minws)
        if (lwork < minws) then
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
  
    if ((nb >= nbmin) .and. (nb < sminmn) .and. (nx < sminmn)) then
      ! use blocked code initially.
      j = nfxd + 1
      
      ! compute factorization: while loop.
      topbmn = minmn - nx
      do while (j <= topbmn)
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
    if (j <= minmn) then
      call rdlaqp2(m, n-j+1, j-1, a(1, j), lda, jpvt(j), tau(j), work(j), work(n+j), work(2*n+1), tol, rank)
    end if
  
  end if
  
  work(1) = iws
  
  return
end subroutine
