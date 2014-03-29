! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


!  -- lapack routine (version 3.1) --
  ! univ. of tennessee, univ. of california berkeley and nag ltd..
  ! november 2006


subroutine rdormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
  use lapack
  implicit none
  
  ! in/out
  character(len=1), intent(in) :: side, trans
  integer, intent(in) :: m, n, k, lda, ldc, lwork
  integer, intent(out) :: info
  double precision, intent(out) :: tau(*), work(*)
  double precision, intent(inout) :: a(lda, *), c(ldc, *)
  ! local
  integer, parameter :: nbmax = 64, ldt = nbmax+1
  logical :: left, lquery, notran
  integer :: i, i1, i2, i3, ib, ic, iinfo, iws, jc, ldwork, &
             lwkopt, mi, nb, nbmin, ni, nq, nw
  double precision :: t(ldt, nbmax)
  ! intrinsics
  intrinsic          max, min
  
  
  ! test the input arguments
  
      info = 0
      left = lsame(side, 'l')
      notran = lsame(trans, 'n')
      lquery = (lwork == -1)
  
  ! nq is the order of q and nw is the minimum dimension of work
  if(left) then
     nq = m
     nw = n
  else
     nq = n
     nw = m
  end if
  
  if(.not.left .and. .not.lsame(side, 'r')) then
     info = -1
  else if(.not.notran .and. .not.lsame(trans, 't')) then
     info = -2
  else if(m < 0) then
     info = -3
  else if(n < 0) then
     info = -4
  else if(k < 0 .or. k > nq) then
     info = -5
  else if(lda < max(1, nq)) then
     info = -7
  else if(ldc < max(1, m)) then
     info = -10
  else if(lwork < max(1, nw) .and. .not.lquery) then
     info = -12
  end if
  
  if(info == 0) then
    ! determine the block size.  nb may be at most nbmax, where nbmax
    ! is used to define the local array t.
    nb = min(nbmax, ilaenv(1, 'dormqr', side // trans, m, n, k, -1))
    lwkopt = max(1, nw)*nb
    work(1) = lwkopt
  end if
  
  if(info /= 0) then
    call xerbla('dormqr', -info)
    return
  else if(lquery) then
    return
  end if
  
  ! quick return if possible
  if(m == 0 .or. n == 0 .or. k == 0) then
    work(1) = 1
    return
  end if
  
  nbmin = 2
  ldwork = nw
  if(nb > 1 .and. nb < k) then
    iws = nw*nb
    if(lwork < iws) then
      nb = lwork / ldwork
      nbmin = max(2, ilaenv(2, 'dormqr', side // trans, m, n, k, -1))
    end if
  else
    iws = nw
  end if
  
  if(nb < nbmin .or. nb.ge.k) then
    ! use unblocked code
    call dorm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, iinfo)
  else
    ! use blocked code
    if((left .and. .not.notran) .or. (.not.left .and. notran)) then
      i1 = 1
      i2 = k
      i3 = nb
    else
      i1 = ((k-1) / nb)*nb + 1
      i2 = 1
      i3 = -nb
    end if
    
    if(left) then
      ni = n
      jc = 1
    else
      mi = m
      ic = 1
    end if
    
    do i = i1, i2, i3
      ib = min(nb, k-i+1)
      
      ! form the triangular factor of the block reflector
      ! h = h(i) h(i+1) . . . h(i+ib-1)
      call dlarft('forward', 'columnwise', nq-i+1, ib, a(i, i), lda, tau(i), t, ldt)
      
      if(left) then
        ! h or h' is applied to c(i:m,1:n)
        mi = m - i + 1
        ic = i
      else
        ! h or h' is applied to c(1:m,i:n)
        
        ni = n - i + 1
        jc = i
      end if
      
      ! apply h or h'
      call dlarfb(side, trans, 'forward', 'columnwise', mi, ni, ib, a(i, i), lda, t, ldt, c(ic, jc), ldc, work, ldwork)
    end do
  end if
  work(1) = lwkopt
  return
  
  return
end subroutine

