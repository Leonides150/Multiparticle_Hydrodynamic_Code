module swap_traj_para
  implicit none
  real*8 :: kappa=1.d0!2.5d0!4.2!3.2d0
  real*8 :: b=0.01d0!0.65d0
  real*8 :: Aswap
end module swap_traj_para
!__________________________!
module file_para
  integer :: num=307
end module file_para
!__________________________!
program main
  use file_para
  implicit none
  complex *16 :: ewald_sum(2)
  real*8 :: v_s,ro(2),incr
  integer :: i
  !real*8 :: t(11852),y(11852),dy(11852)
  real*8 :: t(614),y(614),dy(614)

  open(2,file='swap_traj_fit_2-7_1.dat',status='unknown')

  call read_para(t,y,dy)

  call get_sigma
  call swap_traj_amp_auto_calc
    
  do i=1,num
     incr = 0.025d0
     ro(1) = 0.5d0*y(2*i-1)!incr*dble(i)
     ro(2) = 0.d0!0.5d0*y(2*i)  !0.d0
     
     call get_sigma
     call swap_trajec_vel(ro,v_s)
     call ewald(ro,ewald_sum)
     write(2,*)y(2*i-1),y(2*i),dy(2*i-1),dy(2*i),REALPART(ewald_sum(1)),REALPART(ewald_sum(2)),v_s
     !write(2,*)ro(1),ro(2),REALPART(ewald_sum(1)),REALPART(ewald_sum(2))
     
  end do
  
  close(2)
  
end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine swap_trajec_vel(ro,v_s)
  use swap_traj_para
  implicit none
  real*8 :: g_ro,v_s
  real*8 :: ro(2)
  real*8,external :: swap_repul_strength

  !Aswap = 87.d0!35.d0!8.d0!3.325d0
  
  write(*,*)'from_swap_trajec_vel=',Aswap
  g_ro = swap_repul_strength(ro)
  v_s = Aswap*g_ro * ro(1)

  !write(*,*)'subroutine_swap_trajec_vel_g_ro',g_ro
  !write(*,*)'subroutine_swap_trajec_vel_v_s=',v_s  
    
end subroutine swap_trajec_vel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function swap_repul_strength(ro)
  use swap_traj_para
  implicit none
  real*8 :: numerator,denominator,ro(2),ro_mag
!  real*8 :: b,kappa

 ! b=1.d0
 ! kappa=1.d0
 ! kappa=3.2d0!2.8d0!2.1d0!1.8d0!0.8d0
 ! b=0.65d0  

  ro_mag = sqrt( (ro(1)**2) + (ro(2)**2) )
  numerator = exp(-kappa*ro_mag)
  denominator = ((b**2)+(ro_mag**2))**2

  swap_repul_strength = numerator/denominator
  
  !write(*,*)'func_swap_repul_strength=',swap_repul_strength
  
end function swap_repul_strength
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!To read initial parameters!!!Not related to swapping velocity!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_sigma
  use parameters
  use main_module
  implicit none

  sigma = f_sigma*Lx
  !write(*,*)'sigma=',sigma

end subroutine get_sigma
!____________________________!
subroutine read_para(t,y,dy)
  use file_para
  implicit none
  integer :: k
  !real*8 :: t(11852),y(11852),dy(11852)
  real*8 :: t(614),y(614),dy(614)

  open(1,file='long_box_dx2_7.dat',status='unknown')
  !open(2,file='test_output_file.dat',status='unknown')
  do k=1,num
     
     read(1,*)t(k),y(2*k-1),y(2*k),dy(2*k-1),dy(2*k)
     !write(2,*)t,y(2*k-1),y(2*k),dy(2*k-1),dy(2*k)
     
  end do
 
  close(1)
  !close(2)
  
end subroutine read_para
!______________________________!
subroutine swap_traj_amp_auto_calc
  use swap_traj_para
  implicit none
  real*8 :: ro_st(2),ro_st_mag
  real*8 :: p,g,q
  complex *16 :: ewald_sum_st(2)

  !___preferred stopping distance___!
  ro_st(1) = 0.5d0*4.3398227692000004d0 
  ro_st(2) = 0.d0
  
  ro_st_mag = sqrt( (ro_st(1)**2) + (ro_st(2)**2) )
  g = ( ro_st(1)*exp(-kappa*ro_st_mag) )/(b**2 + ro_st_mag**2)**2
  call ewald(ro_st,ewald_sum_st)
  q = REALPART(ewald_sum_st(1))
  write(*,*)'g=',g
  write(*,*)'q=',q
  Aswap = -q/g
  write(*,*)'from_swap_traj_amp_auto=',Aswap
    
end subroutine swap_traj_amp_auto_calc



