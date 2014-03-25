! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt



! Modified LAPACK routine dgels, original copyright:
!  -- lapack driver routine (version 3.3.1) --
!  -- lapack is a software package provided by univ. of tennessee,    --
!  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
!  -- april 2011


subroutine rdgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info, &
                  tol, coef, eff, ft, rsd, tau, rank)
  use lapack
  use lmfit_utils
  
  implicit none
  
  ! new arguments
  integer, intent(out) :: rank
  double precision, intent(in) :: tol
  double precision, intent(out) :: coef(n, *)
  double precision, intent(out), dimension(ldb, *) :: ft, eff, rsd, tau
  
  double precision :: qraux1
  
!     .. scalar arguments ..
  character          trans
  integer            info, lda, ldb, lwork, m, n, nrhs
!     ..
!     .. array arguments ..
  double precision   a(lda, *), b(ldb, *), work(*)
!     .. parameters ..
  double precision   zero, one
  parameter          (zero = 0.0d0, one = 1.0d0)
!     ..
!     .. local scalars ..
  logical            lquery, tpsd
  integer            brow, i, iascl, ibscl, j, mn, nb, scllen, wsize
  double precision   anrm, bignum, bnrm, smlnum
!     ..
!     .. local arrays ..
  double precision   rwork(1)
!     ..
!!     .. external functions ..
!  logical            lsame
!  integer            ilaenv
!  double precision   dlamch, dlange
!  external           lsame, ilaenv, dlabad, dlamch, dlange
!!     ..
!!     .. external subroutines ..
!  external           dgelqf, dgeqrf, dlascl, dlaset, dormlq, dormqr, dtrtrs, xerbla
!     ..
!     .. intrinsic functions ..
  intrinsic          dble, max, min
!     ..
!     .. executable statements ..
!
!     test the input arguments.
!
  info = 0
  mn = min(m, n)
  lquery = (lwork == -1)
  if(.not.(lsame(trans, 'n') .or. lsame(trans, 't'))) then
     info = -1
  else if(m < 0) then
     info = -2
  else if(n < 0) then
     info = -3
  else if(nrhs < 0) then
     info = -4
  else if(lda < max(1, m)) then
     info = -6
  else if(ldb < max(1, m, n)) then
     info = -8
  else if(lwork < max(1, mn+max(mn, nrhs)) .and. .not.lquery) then
     info = -10
  end if
!
!     figure out optimal block size
!
  if(info == 0 .or. info == -10) then
!
     tpsd = .true.
     if(lsame(trans, 'n')) tpsd = .false.
!
     if(m >= n) then
        nb = ilaenv(1, 'dgeqrf', ' ', m, n, -1, -1)
        if(tpsd) then
           nb = max(nb, ilaenv(1, 'dormqr', 'ln', m, nrhs, n, -1))
        else
           nb = max(nb, ilaenv(1, 'dormqr', 'lt', m, nrhs, n, -1))
        end if
     else
        nb = ilaenv(1, 'dgelqf', ' ', m, n, -1, -1)
        if(tpsd) then
           nb = max(nb, ilaenv(1, 'dormlq', 'lt', n, nrhs, m, -1))
        else
           nb = max(nb, ilaenv(1, 'dormlq', 'ln', n, nrhs, m, -1))
        end if
     end if
!
     wsize = max(1, mn+max(mn, nrhs)*nb)
     work(1) = dble(wsize)
!
  end if
!
  if(info /= 0) then
     call xerbla('dgels ', -info)
     return
  else if(lquery) then
     return
  end if
!
!     quick return if possible
!
  if(min(m, n, nrhs) == 0) then
     call dlaset('full', max(m, n), nrhs, zero, zero, b, ldb)
     return
  end if
!
!     get machine parameters
!
  smlnum = dlamch('s') / dlamch('p')
  bignum = one / smlnum
  call dlabad(smlnum, bignum)
!
!     scale a, b if max element outside range [smlnum,bignum]
!
  anrm = dlange('m', m, n, a, lda, rwork)
  iascl = 0
  if(anrm > zero .and. anrm < smlnum) then
!
!        scale matrix norm up to smlnum
!
     call dlascl('g', 0, 0, anrm, smlnum, m, n, a, lda, info)
     iascl = 1
  else if(anrm > bignum) then
!
!        scale matrix norm down to bignum
!
     call dlascl('g', 0, 0, anrm, bignum, m, n, a, lda, info)
     iascl = 2
  else if(anrm == zero) then
!
!        matrix all zero. return zero solution.
!
     call dlaset('f', max(m, n), nrhs, zero, zero, b, ldb)
     go to 50
  end if
!
  brow = m
  if(tpsd) brow = n
  bnrm = dlange('m', brow, nrhs, b, ldb, rwork)
  ibscl = 0
  if(bnrm > zero .and. bnrm < smlnum) then
!
!        scale matrix norm up to smlnum
!
     call dlascl('g', 0, 0, bnrm, smlnum, brow, nrhs, b, ldb, info)
     ibscl = 1
  else if(bnrm > bignum) then
!
!        scale matrix norm down to bignum
!
     call dlascl('g', 0, 0, bnrm, bignum, brow, nrhs, b, ldb, info)
     ibscl = 2
  end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Copy B over to RSD for later residual calculation
  call dlacpy('All', m, nrhs, b, ldb, rsd, ldb)
  
  
  if(m >= n) then
!
!        compute qr factorization of a
!
!!!!!!!!!!!!!!!!!!!!!!!!!
    call dgeqrf(m, n, a, lda, work(1), work(mn+1), lwork-mn, info)
    
     if(.not.tpsd) then
!
!           least-squares problem min || a * x - b ||
!
!           b(1:m,1:nrhs) := q**t * b(1:m,1:nrhs)
!
        call rdormqr('left', 'transpose', m, nrhs, n, a, lda, work(1), b, ldb, work(mn+1), lwork-mn, info)
        
        call dlacpy('All', n, nrhs, b, ldb, coef, ldb)
!!!!!!!        ! Effects array
!!!!!!!        call dlacpy('All', m, nrhs, b, ldb, eff, ldb)
!
!           workspace at least nrhs, optimally nrhs*nb
!
!           b(1:n,1:nrhs) := inv(r) * b(1:n,1:nrhs)
!
        ! Store qraux(1) == work(1) 
        qraux1 = work(1)
        
        call dtrtrs('upper', 'no transpose', 'non-unit', n, nrhs, a, lda, b, ldb, info)
        
! Produce fitted.values = Ax = Q*(R*x)
        ! Copy over the first RANK elements of numerical soln X
        call dlacpy('All', m, nrhs, b, ldb, ft, ldb)
        
        ! Pretend A="QR" is the upper triangular R and compute R*x
        call dtrmm('L', 'U', 'N', 'N', n, nrhs, 1.0d0, a, lda, ft, ldb)
        
        ! Compute fitted FT = Q*(R*x)
        call dormqr('L', 'N', m, nrhs, n, a, lda, tau, ft, ldb, work(1), lwork, info)
        
        ! Compute residual RSD = FT-b
        call dgeadd('N', m, nrhs, -1.0d0, ft, ldb, 1.0d0, rsd, ldb)
        
!
        if(info > 0) then
           return
        end if
!
        scllen = n
!
     else
!
!           overdetermined system of equations a**t * x = b
!
!           b(1:n,1:nrhs) := inv(r**t) * b(1:n,1:nrhs)
!
        call dtrtrs('upper', 'transpose', 'non-unit', n, nrhs, a, lda, b, ldb, info)
!
        if(info > 0) then
           return
        end if
!
!           b(n+1:m,1:nrhs) = zero
!
        do j = 1, nrhs
           do i = n + 1, m
              b(i, j) = zero
           end do
        end do
!
!           b(1:m,1:nrhs) := q(1:n,:) * b(1:n,1:nrhs)
!
        call rdormqr('left', 'no transpose', m, nrhs, n, a, lda, work(1), b, ldb, work(mn+1), lwork-mn, info)
!
!           workspace at least nrhs, optimally nrhs*nb
!
        scllen = m
!
     end if
!
  else
!
!        compute lq factorization of a
!
     call rdgelqf(m, n, a, lda, work(1), work(mn+1), lwork-mn, info)
!
!        workspace at least m, optimally m*nb.
!
     if(.not.tpsd) then
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
              b(i, j) = zero
           end do
        end do
!
!           b(1:n,1:nrhs) := q(1:n,:)**t * b(1:m,1:nrhs)
!
        call dormlq('left', 'transpose', n, nrhs, m, a, lda, work(1), b, ldb, work(mn+1), lwork-mn, info)
!
!           workspace at least nrhs, optimally nrhs*nb
!
        scllen = n
!
     else
!
!           overdetermined system min || a**t * x - b ||
!
!           b(1:n,1:nrhs) := q * b(1:n,1:nrhs)
!
        call dormlq('left', 'no transpose', n, nrhs, m, a, lda, work(1), b, ldb, work(mn+1), lwork-mn, info)
        
        ! workspace at least nrhs, optimally nrhs*nb
        
        ! b(1:m,1:nrhs) := inv(l**t) * b(1:m,1:nrhs)
        
        call dtrtrs('lower', 'transpose', 'non-unit', m, nrhs, a, lda, b, ldb, info)
        
        if(info > 0) then
           return
        end if
        scllen = m
     end if
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
end



