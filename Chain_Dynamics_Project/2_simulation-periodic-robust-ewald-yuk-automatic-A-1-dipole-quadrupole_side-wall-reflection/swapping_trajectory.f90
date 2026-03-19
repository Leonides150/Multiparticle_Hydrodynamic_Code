
! program main
!   implicit none
!   complex *16 :: ewald_sum(2)
!   real*8 :: v_s,ro(2),incr
!   integer :: i

!   open(1,file='swapping_vel.dat',status='unknown')

!   do i=1,500
!      incr = 0.01d0
!      ro(1) = incr*dble(i)
!      ro(2) = incr*dble(i)
!      ro(2) = 1.5d0

!      call get_sigma
!      call swap_trajec_vel(ro,v_s)
!      call ewald(ro,ewald_sum)
!      write(1,*)ro(1),ro(2),v_s,REALPART(ewald_sum(1)),REALPART(ewald_sum(2))
     
!   end do
  
!   close(1)
  
! end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine swap_traj_activation(i,j,v_s)
  use new_cood_yukawa
  use main_module
  implicit none
  
  integer :: i,j
  real*8 :: v_s(2),ro(2)

  ro(1) = new_relpos(1)/dia
  ro(2) = new_relpos(2)/dia

  if(i.eq.j)then
     v_s = 0.d0
  else
     call swap_trajec_vel(ro,v_s)
  end if

  
end subroutine swap_traj_activation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine swap_trajec_vel(ro,v_s)
  use swap_traj_para
  implicit none
  real*8 :: g_ro,v_s(2)
  real*8 :: ro(2)
  real*8,external :: swap_repul_strength

  !!!call swap_traj_amp_auto_calc
  g_ro = swap_repul_strength(ro)
  v_s(1) = Aswap*g_ro * ro(1)
  v_s(2) = Aswap*g_ro * ro(2)
  !write(*,*)"swap_traj_subroutine_from_sub",v_s
  !write(*,*)'Aswap=',Aswap,'g_ro=',g_ro,'from swap_traj'

!!!write(*,*)'subroutine_swap_trajec_vel_g_ro',g_ro
!!!write(*,*)'subroutine_swap_trajec_vel_v_s=',v_s  
  
end subroutine swap_trajec_vel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function swap_repul_strength(ro)
  use swap_traj_para
  implicit none
  real*8 :: numerator,denominator,ro(2),ro_mag,ro_mag_shift
  !real*8 :: b,kappa

  ro_mag = sqrt( (ro(1)**2) + (ro(2)**2) )
  ro_mag_shift=ro_mag-x0
  numerator = 1.d0/(1.d0 + exp(kappa*ro_mag_shift))
  denominator = ((b**2)+(ro_mag**2))**2

  swap_repul_strength = numerator/denominator
  
!!!write(*,*)'func_swap_repul_strength=',swap_repul_strength
  
end function swap_repul_strength
!_________________________________!
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
  ro_st_mag_shift=ro_st_mag-x0
  
  g = (ro_st(1))/((1.d0 + exp(kappa*ro_st_mag_shift))*(b**2 + ro_st_mag**2)**2)
  call ewald(ro_st,ewald_sum_st)
  q = REALPART(ewald_sum_st(1))
  !write(*,*)'g=',g
  !write(*,*)'q=',q
  Aswap = -q/g
  !write(*,*)'from_swap_traj_amp_auto=',Aswap
    
end subroutine swap_traj_amp_auto_calc
