module swap_traj_para
  implicit none
  !real*8 :: kappa=8.50d0
  real*8 :: kappa=6.d0!!!5.8d0!5.5d0!4.2!3.2d0
!Sagnik  real*8 :: b=1.5d0,x0=1.75d0!b=2.8d0,x0=1.65d0!!b=2.d0,x0=1.782d0!0.65d0
!  real*8 :: b=0.d0,x0=1.95d0!b=2.8d0,x0=1.65d0!!b=2.d0,x0=1.782d0!0.65d0
  !real*8 :: b=0.d0,x0=1.95d0,rev_h=1.d0,z_shift=0.d0
  !!real*8 :: b=1.5d0,x0=1.95d0,rev_h=0.6d0, z_shift = 1.d0
  real*8 :: b=3.2d0,x0=1.95d0,rev_h=0.2d0, z_shift = 1.d0
!!!Note: Decreasing rev_h to make the plot rise in x[1.5-2.15] range
!!!      Increasing b to mke the plot dip back down in x[1-1.5] range - also increases A_swap value
!!!      Decreasing kappa lowers the tip and - also decreases A_swap value
  
  real*8 :: Aswap
end module swap_traj_para
!__________________________!
module file_para
  integer :: num1=307,num2=5926
end module file_para
!__________________________!
program main
  use file_para
  implicit none
  complex *16 :: ewald_sum(2)
  real*8 :: v_s,ro(2),incr,vq(2),ro_z
  integer :: i
  !real*8 :: t(11852),y(11852),dy(11852)
  real*8, allocatable :: t(:),y(:),dy(:),z(:),dz(:)

  !open(1,file='swap_traj_fit_1.dat',status='unknown')
  !Near-field!
  open(2,file='swap_traj_fit_2-7_2.dat',status='unknown')

  allocate(t(2*num1),y(2*num1),dy(2*num1),z(num1),dz(num1))
  call read_para(t,y,dy,num1,z,dz)

  call get_sigma
  call swap_traj_amp_auto_calc
    
  do i=1,num1!5926
     incr = 0.025d0
     ro(1) = 0.5d0*y(2*i-1)!incr*dble(i)
     ro(2) = 0.d0!0.5d0*y(2*i)  !0.d0
     ro_z = 0.5d0*z(i) !1.d0
     
     call get_sigma
     call swap_trajec_vel(ro,v_s,ro_z)
     call ewald(ro,ewald_sum)
     write(2,*)y(2*i-1),y(2*i),dy(2*i-1),dy(2*i),REALPART(ewald_sum(1)),REALPART(ewald_sum(2)),v_s
     
  end do

  deallocate(t,y,dy,z,dz)
  !Far-field!
  open(3,file='swap_traj_fit_8_1.dat',status='unknown')
  allocate(t(2*num2),y(2*num2),dy(2*num2),z(num2),dz(num2))
  call read_para(t,y,dy,num2,z,dz)

  call get_sigma
  call swap_traj_amp_auto_calc
    
  do i=1,num2!5926
     incr = 0.025d0
     ro(1) = 0.5d0*y(2*i-1)!incr*dble(i)
     ro(2) = 0.d0!0.5d0*y(2*i)  !0.d0
     ro_z = 0.5d0*z(i) !1.d0
     
     call get_sigma
     call swap_trajec_vel(ro,v_s,ro_z)
     call ewald(ro,ewald_sum)
     write(3,*)y(2*i-1),y(2*i),dy(2*i-1),dy(2*i),REALPART(ewald_sum(1)),REALPART(ewald_sum(2)),v_s
     !call quadrupole(ro,vq)
     !write(3,*)y(2*i-1),y(2*i),dy(2*i-1),dy(2*i),vq(1),vq(2),v_s
     !vq(1),vq(2)
  end do

  !Extends to 1!
  open(4,file='extend-to-1.dat',status='unknown')

  call get_sigma
  call swap_traj_amp_auto_calc

  do i=1,100
     incr = -0.0035d0
     ro(1) =  0.5d0*2.6999988556000001 + dble(i)*incr 
     ro(2) = 0.d0
     ro_z =  1.d0
     
     call get_sigma
     call swap_trajec_vel(ro,v_s,ro_z)
     call ewald(ro,ewald_sum)
     write(4,*)ro(1),ro(2),0,0,REALPART(ewald_sum(1)),REALPART(ewald_sum(2)),v_s
  end do
  
  !close(1)
  close(2)
  close(3)
  close(4)
  
end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine swap_trajec_vel(ro,v_s,ro_z)
  use swap_traj_para
  implicit none
  real*8 :: g_ro,v_s
  real*8 :: ro(2),ro_z
  real*8,external :: swap_repul_strength

  !Aswap = 87.d0!35.d0!8.d0!3.325d0
  
  write(*,*)'from_swap_trajec_vel=',Aswap
  g_ro = swap_repul_strength(ro)
  v_s = Aswap*g_ro * ro(1) * (ro_z + z_shift)

  !write(*,*)'subroutine_swap_trajec_vel_g_ro',g_ro
  !write(*,*)'subroutine_swap_trajec_vel_v_s=',v_s  
    
end subroutine swap_trajec_vel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function swap_repul_strength(ro)
  use swap_traj_para
  implicit none
  real*8 :: numerator,denominator,ro(2),ro_mag,ro_mag_shift

  ro_mag = sqrt( (ro(1)**2) + (ro(2)**2) )
  !ro_mag_shift = sqrt( (ro(1)-x0)**2 + (ro(2)**2) )

  ro_mag_shift=ro_mag-x0
  !numerator = exp(-kappa*ro_mag)
  numerator = 1.d0/(rev_h + exp(kappa*ro_mag_shift))
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
subroutine read_para(t,y,dy,n,z,dz)
  use file_para
  implicit none
  integer :: k,n
  real*8 :: t(2*n),y(2*n),dy(2*n),z(n),dz(n)
  !real*8 :: t(614),y(614),dy(614)
  y=0.d0
  dy=0.d0

  if(n.eq.307)then
     open(1,file='long_box_dx2_7.dat',status='unknown')
  elseif(n.eq.5926)then
     open(1,file='long_box_dx8_0.dat',status='unknown')
  end if
  !open(2,file='test_output_file.dat',status='unknown')
  do k=1,n!5926

     read(1,*)t(k),y(2*k-1),z(k),dy(2*k-1),dz(k)
     !!!!old!!!!read(1,*)t(k),y(2*k-1),y(2*k),dy(2*k-1),dy(2*k)
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
  complex *16 :: ewald_sum_st(2),ro_st_mag_shift

  !___preferred stopping distance___!
  ro_st(1) = 0.5d0*4.3398227692000004d0 
  ro_st(2) = 0.d0
  
  ro_st_mag = sqrt( (ro_st(1)**2) + (ro_st(2)**2) )
  !ro_mag_shift = sqrt( (ro_st(1)-x0)**2 + (ro_st(2)**2) )
  ro_st_mag_shift=ro_st_mag-x0
  
  !g = ( ro_st(1)*exp(-kappa*ro_st_mag) )/(b**2 + ro_st_mag**2)**2
  g = (ro_st(1) * (ro_z + z_shift))/((rev_h + exp(kappa*ro_st_mag_shift))*(b**2 + ro_st_mag**2)**2)
  call ewald(ro_st,ewald_sum_st)
  q = REALPART(ewald_sum_st(1))
  write(*,*)'g=',g
  write(*,*)'q=',q
  Aswap = -q/g
  write(*,*)'from_swap_traj_amp_auto=',Aswap
    
end subroutine swap_traj_amp_auto_calc


subroutine quadrupole(ro,vq)
  implicit none
  real*8 :: ro(2),vq(2),ro_mag

  
  ro_mag = sqrt( (ro(1)**2) + (ro(2)**2) )
  vq(1) = (-2*(ro(1)**3) + 6*(ro(2)**2)*ro(1) )/(ro_mag**6)
  vq(2) = ( 2*(ro(2)**3) - 6*(ro(1)**2)*ro(2) )/(ro_mag**6)

end subroutine quadrupole
