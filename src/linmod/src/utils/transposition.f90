! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2014, Schmidt


module transposition
  use, intrinsic :: iso_c_binding
  implicit none
  
  private :: factor, xpose_square, xpose_nonsquare
  public :: xpose
  
  contains
  
  
  subroutine factor(n, ifact, ipower, nexp, npower)
        ! in/out
        integer, intent(in) :: n
        integer, intent(out) :: ifact(8), ipower(8), nexp(8), npower
        ! local
        integer :: ip = 0, ifcur = 0, idiv = 2, npart, iquot
        
        
        npart = n
        
  10    iquot = npart / idiv
        if (npart  - idiv*iquot) 60, 20, 60
  20    if (idiv - ifcur) 40, 40, 30
  30    ip = ip+1
        ifact(ip) = idiv
        ipower(ip) = idiv
        ifcur = idiv
        nexp(ip) = 1
        goto 50
  40    ipower(ip) = idiv*ipower(ip)
        nexp(ip) = nexp(ip) + 1
  50    npart = iquot
        goto 10
  60    if (iquot - idiv) 100, 100, 70
  70    if (idiv - 2) 80, 80, 90
  80    idiv = 3
        goto 10
  90    idiv = idiv + 2
        goto 10
  100   if (npart - 1) 140, 140, 110
  110   if (npart - ifcur) 130, 130, 120
  120   ip = ip + 1
        ifact(ip) = npart
        ipower(ip) = npart
        nexp(ip) = 1
  130   ipower(ip) = npart*ipower(ip)
        nexp(ip) = nexp(ip) + 1
  140   npower = ip
        return
  end subroutine
  
  
  
  subroutine xpose_square(n, a)
    ! in/out
    integer, intent(in) :: n
    double precision, intent(inout) :: a(n, n)
    ! local
    integer :: i, j
    double precision :: tmp
    
    
    !$omp parallel if(n*n > 5000) private(i, j, tmp) default(shared) 
    !$omp do 
    do j = 1, n
      do i = 1, j-1
        tmp = a(j, i)
        a(j, i) = a(i, j)
        a(i, j) = tmp
      end do
    end do
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine
  
  
  ! ACM Algorithm 467 -- Transposition of a rectangular matrix in situ.
  ! by Norman Brenner, mit, 1/72.  cf. alg. 380, cacm, 5/70.
  ! Modified to allocate its own workspace, moved square transposition
  subroutine xpose_nonsquare(n1, n2, a) &
    bind(c, name='xpose_')
        ! in/out
        integer(kind=c_int), intent(in) :: n1, n2
        real(kind=c_double), intent(inout) :: a(*)
        ! local
        logical, allocatable :: moved(:)
        integer :: n12, nwork
        double precision atemp, btemp
        integer i, ia1, ia2, idiv, ip, isoid, istart, &
                itest, m, mmia1, mmia2, mmist, n, ncount, npower
        integer :: iexp(8) = 0, ifact(8) = 0, ipower(8) = 0, nexp(8) = 0
        intrinsic mod
        
        
        n12 = n1*n2
        nwork = (n1+n2) / 2
        allocate(moved(nwork))
        
        n = n1
        m = n1 * n2 - 1
        
        !  modulus m is factored into prime powers.  eight factors
        !  suffice up to m = 2*3*5*7*11*13*17*19 = 9,767,520.
        call factor (m, ifact, ipower, nexp, npower)
        do 40 ip = 1, npower
          iexp(ip) = 0
  40    continue
        !  generate every divisor of m less than m/2.
        idiv = 1
  50    if (idiv >= m / 2) go to 190
        !  the number of elements whose index is divisible by idiv
        !  and by no other divisor of m is the euler totient
        !  function, phi(m/idiv).
        ncount = m / idiv
        do 60 ip = 1, npower
          if (iexp(ip) == nexp(ip)) go to 60
          ncount = (ncount / ifact(ip)) * (ifact(ip) - 1)
  60    continue
        do 70 i = 1, nwork
          moved(i) = .false.
  70    continue
        !  the starting point of a subcycle is divisible only by idiv
        !  and must not appear in any other subcycle.
        istart = idiv
  80    mmist = m - istart
        if (istart == idiv) go to 120
        if (istart > nwork) go to 90
        if (moved(istart)) go to 160
  90    isoid = istart / idiv
        do 100 ip = 1, npower
          if (iexp(ip) == nexp(ip)) go to 100
          if (mod (isoid, ifact(ip)) == 0) go to 160
  100   continue
        if (istart < nwork) go to 120
        itest = istart
  110   itest = mod (n * itest, m)
        if (itest < istart .or. itest > mmist) go to 160
        if (itest > istart .and. itest < mmist) go to 110
  120   atemp = a(istart+1)
        btemp = a(mmist+1)
        ia1 = istart
  130   ia2 = mod (n * ia1, m)
        mmia1 = m - ia1
        mmia2 = m - ia2
        if (ia1 < nwork) moved(ia1) = .true.
        if (mmia1 < nwork) moved(mmia1) = .true.
        ncount = ncount - 2
        !  move two elements, the second from the negative
        !  subcycle.  check first for subcycle closure.
        if (ia2 == istart) go to 140
        if (mmia2 == istart) go to 150
        a(ia1+1) = a(ia2+1)
        a(mmia1+1) = a(mmia2+1)
        ia1 = ia2
        go to 130
  140   a(ia1+1) = atemp
        a(mmia1+1) = btemp
        go to 160
  150   a(ia1+1) = btemp
        a(mmia1+1) = atemp
  160   istart = istart + idiv
        if (ncount > 0) go to 80
        do 180 ip = 1, npower
          if (iexp(ip) == nexp(ip)) go to 170
          iexp(ip) = iexp(ip) + 1
          idiv = idiv * ifact(ip)
          go to 50
  170     iexp(ip) = 0
          idiv = idiv / ipower(ip)
  180   continue
  190   continue
        
        
        deallocate(moved)
        return
        
  end subroutine
  
  
  
  subroutine xpose(m, n, a)
    ! in/out
    integer, intent(in) :: m, n
    double precision, intent(inout) :: a(m, n)
    
    if (m < 2 .or. n < 2) then
      return
    else if (m == n) then
      call xpose_square(n, a)
    else
      call xpose_nonsquare(m, n, a)
    end if
    
    return
  end subroutine
end module

