
subroutine trailing_end_potential(y,dy)
  use main_module
  use para_bk
  implicit none
  real*8 :: y(2*numpar),dy(2*numpar)
  real*8 :: x,grad_pot_bk
  integer :: i

!!!  call pre_cal_g
  
  do i=1,numpar

     x = y(2*i-1) - x_st_bk                            !!!Particle x_dist from rear edge!!!
     call potential_vel_bk(x,grad_pot_bk)
     dy(2*i-1) = dy(2*i-1) + grad_pot_bk - v_front        !!!Add potential and front velocity!!!
     !write(*,*)'i=',i,'x=',x,'pot_vel=',grad_pot_bk
     
  end do

end subroutine trailing_end_potential

subroutine leading_end_potential(y,dy)
  use main_module
  use para_bk
  implicit none
  real*8 :: y(2*numpar),dy(2*numpar)
  real*8 :: x,grad_pot_bk
  integer :: i

  do i=1,numpar

     x = x_st_ft - y(2*i-1)                           !!!Particle x_dist from rear edge!!!
     call potential_vel_bk(x,grad_pot_bk)
     dy(2*i-1) = dy(2*i-1) - frac_ft*grad_pot_bk !- v_front        !!!Add potential and front velocity!!!
     
  end do
    
end subroutine leading_end_potential
  
subroutine potential_vel_bk(x,grad_pot_bk)
  use para_bk
  implicit none
  real*8 :: func1,func2,x,grad_pot_bk

  func1 = ( -1.d0 - (alpha_bk*m_bk*x) )/(x**2)
  func2 = exp(-alpha_bk*m_bk*x)
   grad_pot_bk = -g_sq*func1*func2

   ! write(*,*)'func1=',func1,'func2=',func2,'g_sq=',g_sq

end subroutine potential_vel_bk

subroutine pre_cal_g
  use para_bk
  implicit none
  real*8 :: num,den,frac_bk=0.88d0

  g_sq = frac_bk*num/den
  num = exp(alpha_bk*m_bk*x_st_bk)
  den = (-1-alpha_bk*m_bk*x_st_bk)/(x_st_bk)**2

  end subroutine pre_cal_g
