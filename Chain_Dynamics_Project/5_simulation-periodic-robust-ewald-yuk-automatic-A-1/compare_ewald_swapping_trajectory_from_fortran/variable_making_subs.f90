
subroutine get_kn(ii,jj)
  use parameters
  use parameters_kn
  use main_module
  implicit none
  
  integer :: ii,jj
  real*8 :: nx,ny
  
  nx=dble(ii)                            
  ny=dble(jj)                           
  k_n(1) = (nx/Lx)
  k_n(2) = (ny/Ly)
  k_nmod = sqrt((k_n(1))**2 + (k_n(2))**2)

end subroutine get_kn

!________________________________________!
subroutine get_fourier_term(relpos,prod)
  use parameters
  use parameters_kn

  real*8 :: relpos(2),s_s_1
  complex*16 :: prod,s_s_3,s_s_2,s_s_4
  complex*16,external :: phi_plus,expo_sum
  character(5) :: flag

  s_s_1 = 1.d0/((k_nmod)**2)
  s_s_2 = phi_plus(k_n,k_nmod)
  s_s_3 = expo_sum(relpos,k_n)
  s_s_4 = exp(-pi*(sigma**2)*(k_nmod**2))
  
  prod = s_s_1*s_s_2*s_s_3*s_s_4
  
end subroutine get_fourier_term

!_________________________________!
subroutine get_rel_ro(ii,jj,ro,rel_ro,rel_ro_mag)
  use parameters
  use main_module
  implicit none
  
  integer :: ii,jj
  real*8 :: nx,ny,ro(2),ro_n(2),rel_ro(2),rel_ro_mag
  
  nx = dble(ii)
  ny = dble(jj)
  
  ro_n(1) = nx*Lx
  ro_n(2)= ny*Ly
  rel_ro(1) = ro(1)-ro_n(1)
  rel_ro(2) = ro(2)-ro_n(2)
  rel_ro_mag = sqrt((rel_ro(1)**2)+(rel_ro(2)**2))
 
end subroutine get_rel_ro
!_________________________________!
subroutine get_in_gamma(relpos,incom_gamma)
  use parameters
  implicit none

  real*8 :: gam_arg1,gam_arg2,incom_gamma,relpos(2),relpos_mag_sq
  real*8,external :: gammq

  relpos_mag_sq = ((relpos(1))**2)+((relpos(2))**2)
  gam_arg1 = dble(m)
  gam_arg2 = pi*(1.d0/sigma**2)*(relpos_mag_sq)
  incom_gamma = gammq(gam_arg1,gam_arg2)  !Incomplete Gamma!Q(x)

end subroutine get_in_gamma
!_________________________________!
subroutine cplx_gamma(relpos,c_incom_gamma)
  implicit none

  real*8 :: relpos(2),incom_gamma,relpos_mag
  complex*16 :: c_incom_gamma

  call get_in_gamma(relpos,incom_gamma)
  !! write(*,*)relpos,incom_gamma,'___'
  c_incom_gamma = dcmplx(incom_gamma,0.d0)
  !! write(*,*)relpos,c_incom_gamma,'-~-~'

end subroutine cplx_gamma
