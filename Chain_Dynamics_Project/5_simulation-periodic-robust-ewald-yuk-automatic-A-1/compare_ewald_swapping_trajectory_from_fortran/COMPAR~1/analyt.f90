! Main program for project analyt_grad_correct
! Sagnik Singha
! Created: Mon Sep 12 10:42:51 CDT 2016
!
!
!
!
!____________________________________!
subroutine deriv_gamma(ii,rel_ro,d_gamma)
  use parameters
  use parameters2
  use main_module
  implicit none

  integer :: ii
  real*8 :: rel_ro(2),rel_ro_sq,nx,gam_arg1,gam_arg2
  real*8 ::  dg1,dg2,dg3,dg4,dg5,d_gamma,sim_gam
  real*8, external :: gammln

  nx=dble(ii)
  rel_ro_sq = ((rel_ro(1))**2)+((rel_ro(2))**2)
  gam_arg1 = dble(m)
  gam_arg2 = (pi*rel_ro_sq)/(sigma**2)
  sim_gam = exp(gammln(dble(m)))

  dg1 = exp(-gam_arg2)
  dg2 = -(gam_arg2)**(m-1)
  dg3 = pi/(sigma**2)
  dg4 = 2.d0*(rel_ro(1)-(nx*Lx))
  dg5 = 1.d0/sim_gam

  d_gamma = dg1*dg2*dg3*dg4*dg5

end subroutine deriv_gamma

!___________________________!
subroutine deriv_phi_np(ii,jj,rel_ro,d_phi_np)
  use parameters
  use parameters2
  implicit none

  integer :: ii,jj
  real*8 :: rel_ro(2),rel_ro_mag,angle_phi,nx,ny
  real*8 ::  dp1,re_part,cm_part
  complex*16 :: dp2,dp3,d_phi_np
  real*8,external :: theta

  nx=dble(ii)
  ny=dble(jj)

  rel_ro_mag = sqrt((rel_ro(1)**2)+(rel_ro(2)**2))
  angle_phi = theta(rel_ro)
  
  dp1 = rel_ro_mag**(-m-2)
  dp2 = exp(i*dble(m)*angle_phi)
  re_part = rel_ro(1)             !(rel_ro(1)-(nx*Lx))
  cm_part = rel_ro(2)             !(rel_ro(2)-(ny*Ly))
  dp3 = dcmplx(re_part,cm_part)
  d_phi_np = (-0.5d0)*dp1*dp2*dp3

end subroutine deriv_phi_np
!___________________________!
subroutine deriv_exp(ii,jj,rel_ro,f_t)
  use parameters
  use parameters_kn
  use main_module
  implicit none

  integer :: ii,jj
  real*8 :: rel_ro(2),nx,ny,f_t1
  complex*16 :: d_exp,de1,de2,f_t2,f_t3,f_t
  complex*16,external :: phi_plus

  nx=dble(ii)
  ny=dble(jj)
  de1 = (2*pi*i*nx)/Lx
  de2 = exp((2*pi*i)*((rel_ro(1)*(nx/Lx))+(rel_ro(2)*(ny/Ly))))
  d_exp = de1*de2
  
  f_t1 = 1.d0/((k_nmod)**2)
  f_t2 = phi_plus(k_n,k_nmod)
  f_t3 = exp(-pi*(sigma**2)*(k_nmod**2))
  f_t = f_t1*f_t2*f_t3*d_exp
  
end subroutine deriv_exp
