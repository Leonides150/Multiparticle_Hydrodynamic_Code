! Main program for project gradient_correct_Ewald
! Sagnik Singha
! Created: Mon Sep  5 16:20:23 CDT 2016
!
!
!
!

! program main
!   use parameters
!   use parameters2
!   implicit none
  
!   real*8 :: x,y,relpos(2)
!   !real*8 :: grad(2)
!   complex*16 :: grad(2),grad_ewald(2)
!   integer :: ig

!   Open(1,FILE='grad_ewald_xy3.dat',STATUS='UNKNOWN')
!   !!Open(1,FILE='grad_gamma_comp12.dat',STATUS='UNKNOWN')
!   !!Open(2,FILE='grad_gamma_comp22.dat',STATUS='UNKNOWN')
!   write(1,*)'#','sigma=',sigma,'imax=',imax,'n1max=',n1max
    
!   do ig = 1,imax
!      write(*,*)ig
!      relpos(1) = ini_valx + incre_x*dble(ig)
!      relpos(2) = ini_valx + incre_x*dble(ig)

!      call ewald(relpos,grad_ewald)
!      !call grad_xy(relpos,grad)
!      write(1,*)relpos(1),REALPART(grad_ewald(1)),REALPART(grad_ewald(2))
!   end do
  
!   close(1)
!   close(2)
! end program main

!____________________________________!
subroutine grad_xy(relpos,grad,flag)
  use parameters
  use parameters2
  use parameters_kn
  implicit none

  real*8 :: relpos(2),coordx_p(2),coordx_m(2),coordy_p(2),coordy_m(2)
  complex*16 :: x_comp,ewald_sum_1,ewald_sum_2,grad(2)
  !real*8 :: grad(2)
  complex*16 :: func_plusx,func_minx,func_plusy,func_miny
  !real*8 :: func_plusx,func_minx,func_plusy,func_miny
  character(5) :: flag
  real*8,external :: gammq
  real*8 :: angle1,angle2,angle3,angle4

  coordx_p(1) = relpos(1) + h
  coordx_p(2) = relpos(2)

  coordx_m(1) = relpos(1) - h
  coordx_m(2) = relpos(2)

  coordy_p(1) = relpos(1)
  coordy_p(2) = relpos(2) + h

  coordy_m(1) = relpos(1)
  coordy_m(2) = relpos(2) - h

  !write(*,*)'coordx_p',coordx_p
  !write(*,*)'coordx_m',coordx_m
  !write(*,*)'coordy_p',coordy_p
  !write(*,*)'coordy_m',coordy_m
  !write(*,*)'-~-~-~'

  if(flag.eq."ewald")then

     call ewald(coordx_p,func_plusx)
     call ewald(coordx_m,func_minx)
     call ewald(coordy_p,func_plusy)
     call ewald(coordy_m,func_miny)

  elseif(flag.eq."phinp")then

     call phi_minus_np(coordx_p,func_plusx,angle1)
     call phi_minus_np(coordx_m,func_minx,angle2)
     call phi_minus_np(coordy_p,func_plusy,angle3)
     call phi_minus_np(coordy_m,func_miny,angle4)

  elseif(flag.eq."gamma")then

     !write(*,*)'coordx_p',coordx_p
     !write(*,*)'coordx_m',coordx_m
     !write(*,*)'coordy_p',coordy_p
     !write(*,*)'coordy_m',coordy_m
     !write(*,*)'+++++++'

     call cplx_gamma(coordx_p,func_plusx)
     call cplx_gamma(coordx_m,func_minx)
     call cplx_gamma(coordy_p,func_plusy)
     call cplx_gamma(coordy_m,func_miny)

  elseif(flag.eq."exsum")then
     !think and write expo part in second sum, ro and k_n creating problem
     call get_fourier_term(coordx_p,func_plusx)
     call get_fourier_term(coordx_m,func_minx)
     call get_fourier_term(coordy_p,func_plusy)
     call get_fourier_term(coordy_m,func_miny)

  end if

  grad(1) = (func_plusx-func_minx)/(2.d0*h)
  grad(2) = (func_plusy-func_miny)/(2.d0*h)

end subroutine grad_xy


