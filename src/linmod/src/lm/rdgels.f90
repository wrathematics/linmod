! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


!FIXME change to dgelsy ???



! Modified LAPACK routine dgels, original copyright:
!  -- lapack driver routine (version 3.3.1) --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!  -- april 2011



! Purpose
! =======
!
! Linear model fitter using a QR (or LQ if m<n) decomposition.
!
!
! Arguments
! =========
!
! m         (input) integer
!           The number of rows of the data matrix a.  m>=0.
! 
! n         (input) integer
!           The number of columns of the data matrix a.  n>=0.
! 
! nrhs      (input) integer
!           The number of 'right hand sides', i.e. the number of 
!           columns of b.  nrhs>=0.
! 
! a         (input/output) double precision array, dimension (m,n).
!           On entry, the input data matrix.  On exit, the output
!           of a lapack QR factorization is stored for a.
! 
! b         (input) double precision array, dimension (m,nrhs).
!           The response variable.
! 
! resids    (output)
!           
! 
! tol       (input) double precision
!           Tolerance for the qr decomposition.
! 
! info      (output) integer
!           = 0: successful exit.
!           < 0: if info = -i, the i-th argument had an illegal value.

subroutine rdgels(m, n, nrhs, a, lda, b, ldb, work, lwork, info, &
                  tol, coef, eff, ft, rsd, tau, jpvt, rank)
  use lapack
  use lapack_omp
  use rdgels_utils
  
  implicit none
  
  ! in/out
  integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
  integer, intent(out) :: info, jpvt(n)
  integer, intent(inout) :: rank
  double precision, intent(in) :: tol
  double precision, intent(out) :: work(*), coef(n, *), tau(*)
  double precision, intent(out), dimension(ldb, *) :: ft, eff, rsd
  double precision, intent(inout) :: a(lda, *), b(ldb, *)
  ! FIXME
  double precision :: qraux1
  ! local
  logical :: lquery
  integer :: brow, i, iascl, ibscl, j, mn, nb, scllen, wsize
  double precision :: anrm, bignum, bnrm, smlnum
  double precision :: rwork(1)
  ! functions
  intrinsic :: dble, max, min
  
  
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
  if(info == 0 .or. info == -10) then
    
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
  
  if(info /= 0) then
     call xerbla('dgels ', -info)
     return
  else if(lquery) then
     return
  end if
  
  
  ! quick return if possible
  if(min(m, n, nrhs) == 0) then
     call dlaset('full', max(m, n), nrhs, 0.0d0, 0.0d0, b, ldb)
     return
  end if
  
  
  ! get machine parameters
  smlnum = dlamch('s') / dlamch('p')
  bignum = 1.0d0 / smlnum
  call dlabad(smlnum, bignum)
  
  
  ! scale a, b if max element outside range [smlnum,bignum]
  anrm = dlange('m', m, n, a, lda, rwork)
  iascl = 0
  
  ! scale matrix norm up to smlnum
  if(anrm > 0.0d0 .and. anrm < smlnum) then
    call dlascl('g', 0, 0, anrm, smlnum, m, n, a, lda, info)
    iascl = 1
  else if(anrm > bignum) then
    ! scale matrix norm down to bignum
    call dlascl('g', 0, 0, anrm, bignum, m, n, a, lda, info)
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
  if(bnrm > 0.0d0 .and. bnrm < smlnum) then
    call dlascl('g', 0, 0, bnrm, smlnum, brow, nrhs, b, ldb, info)
    ibscl = 1
  ! scale matrix norm down to bignum
  else if(bnrm > bignum) then
    call dlascl('g', 0, 0, bnrm, bignum, brow, nrhs, b, ldb, info)
    ibscl = 2
  end if
  
  
!--------------------------------------------------------------------
  ! Copy B over to RSD for later residual calculation
  call dlacpy_omp('All', m, nrhs, b, ldb, rsd, ldb)
  
  !$omp parallel do if(n > 5000) private(i) default(shared) 
    do i = 1, n
      jpvt(i) = 0
    end do
  !$omp end parallel do
  
  
  
  if (m >= n) then
    call rdgels_qr(m, n, mn, nrhs, a, lda, b, ldb, work, lwork, info, &
                  tol, coef, eff, ft, rsd, tau, jpvt, rank, qraux1)
    
    if (info > 0) then
       return
    end if
    
    scllen = n
  
  
  ! ----------- Compute LQ factorization of A
  else
  
     call dgelqf(m, n, a, lda, work(1), work(mn+1), lwork-mn, info)
!
!        workspace at least m, optimally m*nb.
!
!
!           underdetermined system of equations a * x = b
!
!           b(1:m,1:nrhs) := inv(l) * b(1:m,1:nrhs)
!
      call dtrtrs('lower', 'no transpose', 'non-unit', m, nrhs, a, lda, b, ldb, info)
!
      if(info > 0) then
         return
      end if
!
!           b(m+1:n,1:nrhs) = 0
!
      do j = 1, nrhs
         do i = m + 1, n
            b(i, j) = 0.0d0
         end do
      end do
!
!           b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
!
      call dormlq('l', 't', n, nrhs, m, a, lda, work(1), b, ldb, work(mn+1), lwork-mn, info)
!
!           workspace at least nrhs, optimally nrhs*nb
!
      scllen = n
!
  end if
  
  
  ! undo scaling
  if(iascl == 1) then
     call dlascl('g', 0, 0, anrm, smlnum, scllen, nrhs, b, ldb, info)
  else if(iascl == 2) then
     call dlascl('g', 0, 0, anrm, bignum, scllen, nrhs, b, ldb, info)
  end if
  if(ibscl == 1) then
     call dlascl('g', 0, 0, smlnum, bnrm, scllen, nrhs, b, ldb, info)
  else if(ibscl == 2) then
     call dlascl('g', 0, 0, bignum, bnrm, scllen, nrhs, b, ldb, info)
  end if
  
50 continue
!!!!!!  work(1) = dble(wsize)
  work(1) = qraux1
  
  return
  
  
  
  
  
  
  
  contains
  subroutine rdgelqf(m, n, a, lda, tau, work, lwork, info)
    integer, intent(in) :: m, n, lda, lwork
    integer, intent(out) :: info
    double precision, intent(out) :: tau
    double precision, intent(inout) :: a(*), work(*)
  end subroutine
  
  
  
  subroutine rdormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    character(len=1), intent(in) :: side, trans
    integer, intent(in) :: m, n, k, lda, ldc, lwork
    integer, intent(out) :: info
    double precision, intent(in) :: a(*), tau(*)
    double precision, intent(out) :: c(*), work(*)
  end subroutine
  
end subroutine



