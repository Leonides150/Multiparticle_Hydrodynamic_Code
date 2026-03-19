! Main program for project multipolar_solutions
! Sagnik Singha
! Created: Mon Apr  4 11:02:57 CDT 2016
!
!
!
!
! module parameters
!   implicit none
!   real*8 :: pi = 4.d0 * atan(1.d0), 
!   real*8 :: EPS=3.d0*(10**(-7)), FPMIN=(10**(0))
!   integer :: ITMAX=100
! end module parameters

! program main
!   implicit none
!   integer :: m,sig
!   real*8 :: a,ro(2),rn(2),x,rel_ro,tau,pi
!   real*8, external :: gammq

!   pi = 4.d0 * atan(1.d0)
!   m=2
!   !sig = 1
!   a = dble(m)
!   !rel_ro = (ro(1)-rn(1))**2 + (ro(2)-rn(2))**2
!   !x = pi*(1.d0/dble(sig)**2)*rel_ro
!   x = 1.d0
!   tau = gammq(a,x)

!   write(*,*)'tau=',tau

! end program main

!!!!!
real*8 function gammq(a,x)
  implicit none 

  real*8 :: a,x
  real*8 :: gammcf,gamser,gln

  if(x.lt.0.d0 .OR. a.le.0.d0)write(*,*) 'bad arguments'
  if(x.lt.a+1.d0)then
     call gser(gamser,a,x,gln)
     !write(*,*)'gamser=',gamser
     gammq = 1.d0 - gamser
  else
     call gcf(gammcf,a,x,gln)
     !write(*,*)gammcf,a,x,gln
     gammq = gammcf
  endif

  !write(*,*)'gammq=',gammq
  
  return 
end function gammq

!!!!!!
subroutine gser(gamser,a,x,gln)
  implicit none
  integer :: ITMAX
  real *8 :: EPS
  real*8 :: a,gamser,gln,x
  integer :: n
  real*8 :: ap,del,sum
  real*8, external :: gammln
  parameter(ITMAX=100, EPS=3.e-7)

  gln = gammln(a)
  !write(*,*)'gln=',gln

  if(x.le.0)then
     if(x.lt.0)write(*,*) 'x<0 in gser'
     gamser = 0
     return
  endif
  ap = a
  sum = 1.d0/a
  del = sum
  do n=1,ITMAX
     ap = ap + 1
     del = (del * x) / ap
     sum = sum + del
     if(abs(del).lt.abs(sum)*EPS)goto 1
  enddo
  write(*,*) 'a#$ too large, ITMAX too small in gser'
1 gamser = sum*exp(-x+a*log(x)-gln)

  return
end subroutine gser

!!!!!!
subroutine gcf(gammcf,a,x,gln)
  implicit none
  integer :: ITMAX
  real*8 :: EPS, FPMIN
  real*8 :: a,gammcf,gln,x
  integer :: i
  real*8 :: an,b,d,del,h
  real*8 :: c
  real*8, external :: gammln
  parameter(ITMAX=100,EPS=3E-7, FPMIN=1.E-30)
!!!
  !write(*,*)FPMIN
!!!
  gln = gammln(a)
  b = x + 1.d0 - a
  c = 1.d0/FPMIN
  d = 1.d0/b
  h = d
  do i = 1,ITMAX
     an = -i*(i-a)
     b = b + 2.d0
     d = an*d + b
     if(abs(d).lt.FPMIN)d=FPMIN
     c = b + an/c
     if(abs(c).lt.FPMIN)c=FPMIN
     d = 1.d0/d
     del = d*c
     h = h*del
     if(abs(del-1.d0).lt.EPS)goto 1
  enddo
  write(*,*) ' a too large, ITMAX too small in gcf'
  stop!!!!_!!!!!!!
1 gammcf = exp(-x+a*log(x)-gln)*h

  !write(*,*)'gammcf=',gammcf 
  return
end subroutine gcf

!!!!!

real*8 function gammln(xx)
  implicit none

  real *8 :: xx
  integer :: j
  real*8 :: ser, stp, tmp, x, y, cof(6)
  save cof,stp
  data cof,stp/76.18009172947146d0,&
       -86.50532032941677d0,24.01409824083091d0,&
       -1.231739572450155d0,0.1208650973866179d-2,&
       -0.5395239384953d-5,2.5066282746310005d0/

  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  do j=1,6
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gammln=tmp+log(stp*ser/x)

  !write(*,*)'gammln=',gammln

  return  
end function gammln






