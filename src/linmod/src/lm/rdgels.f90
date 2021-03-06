! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt

! Extremely modified from LAPACK routine dgels, original copyright:
!  -- lapack driver routine (version 3.3.1) --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!  -- april 2011


subroutine rdgels(m, n, nrhs, a, b, tol, coef, eff, ft, rsd, tau, &
                  jpvt, rank, work, lwork, info) &
bind(c, name='rdgels_')
  use :: lapack
  use :: lapack_omp
  use :: lm_fit_utils
  use :: qr_utils, only : rdgelqf, rdormqr
  
  implicit none
  
  ! in/out
!  type(qr_t) :: qrlist
  integer, intent(in) :: m, n, nrhs, lwork
  integer, intent(out) :: info, jpvt(n)
  integer, intent(inout) :: rank
  double precision, intent(in) :: tol
  double precision, intent(out) :: coef(n, nrhs), tau(*), work(lwork)
  double precision, intent(out), dimension(m, nrhs) :: ft, eff, rsd
  double precision, intent(inout) :: a(m, n), b(m, nrhs)
  ! FIXME
  double precision :: qraux1
  ! local
  integer :: lda, ldb
  logical :: lquery
  integer :: brow, i, iascl, ibscl, j, mn, nb, scllen, wsize
  double precision :: anrm, bignum, bnrm, smlnum
  double precision :: rwork(1)
  ! functions
  intrinsic :: dble, max, min, int
  
  
  lda = m
  ldb = max(m, n)
  
  ! test the input arguments.
  info = 0
  mn = min(m, n)
  lquery = (lwork == -1)
  if (m < 0) then
     info = -2
  else if (n < 0) then
     info = -3
  else if (nrhs < 0) then
     info = -4
  else if (lda < max(1, m)) then
     info = -6
  else if (ldb < max(1, m, n)) then
     info = -8
  else if (lwork < max(1, mn+max(mn, nrhs)) .and. .not.lquery) then
     info = -10
  end if
  
  
  
  ! figure out optimal block size
  if (info == 0 .or. info == -10) then
    
    if(m >= n) then
      nb = ilaenv(1, 'dgeqrf', ' ', m, n, -1, -1)
      nb = max(nb, ilaenv(1, 'dormqr', 'lt', m, nrhs, n, -1))
    else
      nb = ilaenv(1, 'dgelqf', ' ', m, n, -1, -1)
      nb = max(nb, ilaenv(1, 'dormlq', 'ln', n, nrhs, m, -1))
    end if
    
    wsize = max(1, mn+max(mn, nrhs)*nb)
    work(1) = dble(wsize)
  end if
  
  if (info /= 0) then
     call xerbla('dgels ', -info)
     goto 50 ! deallocate temp storage and return
  else if(lquery) then
    goto 50 ! deallocate temp storage and return
  end if
  
  
  ! quick return if possible
  if (min(m, n, nrhs) == 0) then
     call dlaset('full', max(m, n), nrhs, 0.0d0, 0.0d0, b, ldb)
     goto 50 ! deallocate temp storage and return
  end if
  
  
  ! get machine parameters
  smlnum = dlamch('s') / dlamch('p')
  bignum = 1.0d0 / smlnum
  call dlabad(smlnum, bignum)
  
  
  ! scale a, b if max element outside range [smlnum,bignum]
  anrm = dlange('m', m, n, a, lda, rwork)
  iascl = 0
  
  ! scale matrix norm up to smlnum
  if (anrm > 0.0d0 .and. anrm < smlnum) then
    call dlascl_omp('g', 0, 0, anrm, smlnum, m, n, a, lda, info)
    iascl = 1
  else if(anrm > bignum) then
    ! scale matrix norm down to bignum
    call dlascl_omp('g', 0, 0, anrm, bignum, m, n, a, lda, info)
    iascl = 2
  else if(anrm == 0.0d0) then
  ! matrix all zero. return zero solution.
    call dlaset('f', max(m, n), nrhs, 0.0d0, 0.0d0, b, ldb)
    go to 50
  end if
  
  
  brow = m
  
  bnrm = dlange('m', brow, nrhs, b, ldb, rwork)
  ibscl = 0
  
  ! scale matrix norm up to smlnum
  if (bnrm > 0.0d0 .and. bnrm < smlnum) then
    call dlascl_omp('g', 0, 0, bnrm, smlnum, brow, nrhs, b, ldb, info)
    ibscl = 1
  ! scale matrix norm down to bignum
  else if(bnrm > bignum) then
    call dlascl_omp('g', 0, 0, bnrm, bignum, brow, nrhs, b, ldb, info)
    ibscl = 2
  end if
  
  
  ! Copy B over to RSD for later residual calculation
  call dlacpy_omp('All', m, nrhs, b, m, rsd, m)
  
  !$omp parallel do if(n > linmod_omp_minsize) private(i) default(shared) 
    do i = 1, n
      jpvt(i) = 0
    end do
  !$omp end parallel do
  
  
  ! Fit the linear model, do all the extra value wrangling
  if (m >= n) then
    call rdgels_qr(m, n, mn, nrhs, a, b, work, lwork, info, &
                   tol, coef, eff, ft, rsd, tau, jpvt, rank, qraux1)
  else
    call rdgels_lq(m, n, mn, nrhs, a, b, work, lwork, info, &
                   tol, coef, eff, ft, rsd, tau, jpvt, rank, qraux1)
  end if
  
  if (info > 0) return
  
  scllen = n
  
  
  ! undo scaling
  if (iascl == 1) then
    call dlascl_omp('g', 0, 0, anrm, smlnum, scllen, nrhs, b, ldb, info)
  else if(iascl == 2) then
    call dlascl_omp('g', 0, 0, anrm, bignum, scllen, nrhs, b, ldb, info)
  end if
  if (ibscl == 1) then
    call dlascl_omp('g', 0, 0, smlnum, bnrm, scllen, nrhs, b, ldb, info)
  else if(ibscl == 2) then
    call dlascl_omp('g', 0, 0, bignum, bnrm, scllen, nrhs, b, ldb, info)
  end if
  
  
  50 continue
!!!!!!  work(1) = dble(wsize)
  work(1) = qraux1
  
  return
end subroutine

