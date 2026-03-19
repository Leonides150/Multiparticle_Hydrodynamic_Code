! Main program for project wigner
! Sagnik Singha
! Created: Wed Feb 24 13:55:35 CST 2016
!
!
!
!
module parameters
  implicit none
  real*8 :: pi = 4.d0 * atan(1.d0),Cw = 1.3105329259d0,sig=16.d0
  integer :: depth=30
end module parameters
!!!
program main
  use parameters
  implicit none
  real*8 :: ro(2)

  ro(1) = 0.5d0
  ro(2)= 0.5d0
  call wigner(ro)

end program main

!!!wigner function
subroutine wigner(ro)
  use parameters
  implicit none
  real*8 :: ro(2),ro_n(2), k_n(2),Lx,Ly,sumint
  integer :: i,j,nx,ny
  real*8, external :: sum_int
  complex*16, external :: sum_exp
  complex*16 :: w,sumexp

  Lx = 1.d0
  Ly = 1.d0
  sumint = 0.d0
  sumexp = 0.d0

  do i=1,100
     do j=1,100
        nx = i-50
        ny = j-50

        ro_n(1) = nx*Lx
        ro_n(2)= ny*Ly
        k_n(1) = (nx/Lx)
        k_n(2) = (ny/Ly)

        sumint = sumint + sum_int(ro,ro_n)

        if(nx.ne.0)then
           if(ny.ne.0)then
              sumexp = sumexp + sum_exp(ro,k_n)
           end if
        end if

     enddo
  enddo
  
  write(*,*)'sigma=16'
  write(*,*)'sumint=',sumint
  write(*,*)'sumexp=',sumexp
  
  w = (0.5d0*sumint) + (0.5d0*(1.d0/(pi*Lx*Ly))*sumexp) - ((sig**2)/&
       (2*Lx*Ly))+Cw
  
  write(*,*)'w=',w

end subroutine wigner

!!!exponential integral sum
real*8 function sum_int(ro,ro_n)
  use parameters
  integer :: nx,ny
  real*8 :: x,Lx,y,Ly,int
  real*8 :: ro(2), ro_n(2)

  ! nx=1
  ! ny=1
  ! Lx=1.d0
  ! Ly=1.d0
  ! sig=1.d0
  !sum1 = 0.5*rep
  
  xx = (pi*(((ro(1)-ro_n(1))**2 +(ro(2)-ro_n(2))**2)))/sig**2
  call cf(depth,xx)
  sum_int = int
    
end function sum_int

!!!complex sum
complex*16 function sum_exp(ro,k_n)
  use parameters
 
  real*8 :: kn_sq,dot_pr
  real*8 :: ro(2), k_n(2),func1,func2
  complex*16 :: func

  kn_sq = ((k_n(1))**2) + ((k_n(2))**2)
  dot_pr = ((x*nx)/Lx)+((y*ny)/Ly)
  ! func = exp((-pi*(sig**2)*((kn_sq))+(2*pi*(dot_pr))))
  ! sum_exp = func
  func1 = (1/kn_sq)*exp(-pi*(sig**2)*kn_sq)*cos(2*pi*dot_pr)
  func2 = (1/kn_sq)*exp(-pi*(sig**2)*kn_sq)*sin(2*pi*dot_pr)
  func = CMPLX(SNGL(func1), SNGL(func2))
  sum_exp = func
     
end function sum_exp



  
