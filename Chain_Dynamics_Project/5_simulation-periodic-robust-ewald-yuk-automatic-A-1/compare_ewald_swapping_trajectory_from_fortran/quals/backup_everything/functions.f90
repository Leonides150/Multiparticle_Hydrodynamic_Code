

!_________Angle finder_______________!
real*8 function theta(a)
  use parameters
  implicit none

  real*8 :: a(2),ratio,mag,eps,angle,abs_ratio
  real*8,external::theta_cos

  mag = sqrt((a(1)**2)+(a(2)**2))
  eps = 0.0001d0
  ratio = a(2)/a(1)
  abs_ratio = abs(a(1)/mag)

  !write(*,*)a
  
  if(mag.lt.eps)then
     write(*,*)'STOP-Angle Undefined, a,mag=',a,mag
     stop
  elseif(abs_ratio.lt.eps)then
     !write(*,*)'*'
     angle = theta_cos(a,mag)                    
  elseif((a(2).ge.0.d0) .and. a(1).gt.0.d0)then
     !!!write(*,*)'@'
     angle = atan(ratio)         !say theta/1st
  elseif((a(2).gt.0.d0) .and. a(1).lt.0.d0)then
     !!!write(*,*)'@@'
     angle = pi + atan(ratio)    !without adjustment atan gives -theta/2nd 
  elseif((a(2).lt.0.d0) .and. a(1).lt.0.d0)then
     !!!write(*,*)'@@@'
     angle = pi + atan(ratio)    !without adjustment atan gives theta/3rd
  elseif((a(2).lt.0.d0) .and. a(1).gt.0.d0)then
     !!!write(*,*)'@@@@'
     angle = (2*pi) + atan(ratio)  !4th
  elseif(a(2) .eq. 0.d0 .and. a(1).lt.0.d0)then
     angle = pi                  !when angle is pi
  end if

  theta = angle
end function theta
!______for angles closer to 90degree____!
real*8 function theta_cos(a,mag)
  use parameters
  implicit none

  real*8 :: a(2),mag,angle

  if(a(2) .gt. 0.d0)then          !first & second quadrant

     angle = acos(a(1)/mag)
 
  elseif(a(2) .lt. 0.d0)then      !third & fourth quadrant

     angle = (2.d0*pi) - acos(a(1)/mag)
 
  elseif(a(2) .eq. 0.d0 .and. a(1) .lt. 0.d0)then

     angle = pi                   !when angle is pi & not '0'

  endif

  theta_cos = angle
end function theta_cos


!_______Factorial______________!

integer function factorial(nnmax)
  use parameters
  implicit none
  integer :: nn,res,nnmax
  res = 1

  do nn=1, nnmax
     res = res*nn
  enddo

  factorial = res
  !write(*,*)'factorial=',factorial
  
end function factorial







































!#______________Non-Multipolar Phi_Minus___________
! complex*16 function phi_minus(xx,yy)
!   use parameters
!   implicit none

!   real*8 :: i_re,i_im,ro(2),ro_n(2),rel_ro(2),p_m1,p_m2,angle
!   complex*16 :: i,p_m3
!   real*8 :: xx,yy,rel_ro_mag,nx,ny
!   real*8,external :: theta

!   i_re = 0.d0
!   i_im = 1.d0
!   i = cmplx(i_re,i_im)

!   ro(1) = xx
!   ro(2) = yy
 
!   nx = 1.d0
!   ny = 1.d0
!   ro_n(1) = nx*Lx
!   ro_n(2) = ny*Ly
!   rel_ro(1) = ro(1)-ro_n(1)
!   rel_ro(2) = ro(2)-ro_n(2)

!   p_m1 = 0.5d0/dble(m)
!   rel_ro_mag = sqrt((rel_ro(1)**2)+(rel_ro(2)**2))
!   p_m2 = 1.d0/(rel_ro_mag**m)
!   angle = theta(rel_ro)

!   p_m3 = exp(i*dble(m)*angle)


!   phi_minus = p_m1*p_m2*p_m3

! end function phi_minus

!___________Alternate Phi_Minus Non-Polar_____________
! complex*16 function phi_alternate(xx,yy)
!   use parameters
!   implicit none

!   real*8 :: i_re,i_im,xx,yy,ro(2),nx,ny,ro_n(2),rel_ro(2),rel_ro_mag
!   real*8 :: p_a_1,p_a_2
!   complex*16 :: i,p_a_3

!   i_re = 0.d0
!   i_im = 1.d0
!   i = cmplx(i_re,i_im)

!   ro(1) = xx
!   ro(2) = yy

!   nx = 1.d0
!   ny = 1.d0
!   ro_n(1) = nx*Lx
!   ro_n(2) = ny*Ly
!   rel_ro(1) = ro(1)-ro_n(1)
!   rel_ro(2) = ro(2)-ro_n(2)
!   rel_ro_mag = sqrt((rel_ro(1)**2)+(rel_ro(2)**2))

!   p_a_1 = 1.d0/(2.d0*m)
!   p_a_2 = 1.d0/(rel_ro_mag**2)
!   p_a_3 = (rel_ro(1) + (i*rel_ro(2)))**m

!   phi_alternate = p_a_1*p_a_2*p_a_3

! end function phi_alternate
