
!____________________________________!
subroutine ewald(relpos,ewald_sum)
  use parameters
  use parameters2
  use main_module
  implicit none

  real*8 :: relpos(2),angle_np
  complex*16 :: ewald_sum(2),sum1(2),fourier_sum(2),mult,sec_term(2)
  complex*16 :: phi_min_np
  integer,external :: factorial

  call first_sum(relpos,sum1)
  call second_sum(relpos,fourier_sum)
  mult = ((pi**(m-1))*(i**m))/(2*(factorial(m))*Lx*Ly)
  sec_term(1) = fourier_sum(1)*mult
  sec_term(2) = fourier_sum(2)*mult

  ewald_sum(1) = sum1(1) + sec_term(1)
  ewald_sum(2) = sum1(2) + sec_term(2)

end subroutine ewald

!____________________________________!
subroutine first_sum(ro,sum1)
  use parameters
  use flag_self_int
  implicit none
  real*8 :: ro(2),xx,yy,incom_gamma,rel_ro(2),rel_ro_mag
  real*8 :: angle1,p_m1,p_m2,gam_arg1,gam_arg2
  real*8 :: incom_gam_compl,grad_gam_a
  complex*16 :: p_m3,phi_min_np,sum1(2),prod(2),grad_phinp(2),grad_phinp_a
  complex*16 :: grad_gam(2),product_a,sum1_a
  integer :: ii,jj
  character(5) :: flag
  real*8,external :: theta

  sum1 = dcmplx(0.d0,0.d0)

  do ii= -n1max,n1max          
     do jj=-m1max,m1max

        if(flag_self .AND. (ii.eq.0 .AND. jj.eq.0) )then
           continue

        else

           call get_rel_ro(ii,jj,ro,rel_ro,rel_ro_mag)

           if(rel_ro_mag .le. (3.d0*sigma))then!!!!

              call phi_minus_np(rel_ro,phi_min_np,angle1)
              call get_in_gamma(rel_ro,incom_gamma)

              flag = "gamma"
              call grad_xy(rel_ro,grad_gam,flag)
              !call deriv_gamma(ii,rel_ro,grad_gam_a)

              flag = "phinp"
              call grad_xy(rel_ro,grad_phinp,flag)
              !call deriv_phi_np(ii,jj,rel_ro,grad_phinp_a)

              prod(1) = (grad_gam(1)*phi_min_np) + (grad_phinp(1)&
                   *incom_gamma)
              prod(2) = (grad_gam(2)*phi_min_np) + (grad_phinp(2)&
                   *incom_gamma)
              product_a = (grad_gam_a*phi_min_np)+(grad_phinp_a&
                   *incom_gamma)
              !product = phi_min_np*incom_gamma
              sum1(1) = sum1(1) + prod(1)
              sum1(2) = sum1(2) + prod(2) 
              sum1_a = sum1_a + product_a

           else!!!!
              continue!!!!
           end if!!!!

        endif

     end do
  end do
  
end subroutine first_sum

!___________________________________________________________!
subroutine second_sum(relpos,fourier_sum)
  use parameters
  use parameters_kn
  implicit none

  complex*16 :: sum2(2),fourier_sum(2),product2,grad_fp(2),grad_fp_a,sum2_a
  complex*16 :: fourier_sum_a
  real*8 :: relpos(2)
  integer :: ii,jj
  character(5) :: flag
  
  sum2 = dcmplx(0.d0,0.d0)
  
  do ii=-n1max,n1max                                                       
     do jj=-m1max,m1max
        
        call get_kn(ii,jj)
        
        if( (ii.ne.0) .OR.(jj.ne.0) )then

           flag = "exsum"
           call grad_xy(relpos,grad_fp,flag)
           !!!call deriv_exp(ii,jj,relpos,grad_fp_a)
           
           sum2(1) = sum2(1) + grad_fp(1)
           sum2(2) = sum2(2) + grad_fp(2)
           sum2_a = sum2_a + grad_fp_a
           !call get_fourier_term(relpos,product2)
           !sum2 = sum2 + product2
           
        endif
     enddo
  enddo
  
  fourier_sum(1) = sum2(1)
  fourier_sum(2) = sum2(2)
  fourier_sum_a = sum2_a
  
end subroutine second_sum

!__________________________________!
complex*16 function phi_plus(vect,vect_mag)
  use parameters
  use parameters_kn
  implicit none

  real*8 :: vect(2), phi,vect_mag
  real*8 :: parta
  complex*16 :: partb
  real*8 :: theta

  !!i_re = 0.d0
  !!i_im = 1.d0
  !!i = cmplx(i_re,i_im) !i = dcmplx(i_re,i_im) *error source*
  
  phi  = theta(vect)
  parta = vect_mag**m  
  partb = exp(i*m*phi)
  phi_plus = parta*partb

end function phi_plus

!________________________!
complex*16 function expo_sum(ro,k_n)
  use parameters
  implicit none

  real*8 :: dot_pro,ro(2),k_n(2)
  real*8 :: nx,ny
  integer :: ii,jj
  
  dot_pro=dot_product(ro,k_n)
  expo_sum = exp(2.d0*i*(pi)*(dot_pro))

end function expo_sum

!_____________________!
subroutine phi_minus_np(ro,phi_min_np,angle1)
  use parameters
  implicit none
  real*8 :: ro(2),nx,ny,ro_n(2),gamma_handler
  real*8 :: rel_ro(2),ro_mag,angle1,p_m1,p_m2,gam_arg1,gam_arg2
  complex*16 :: p_m3,phi_min_np,phi_sum
  integer :: ii,jj
  real*8,external :: theta

  phi_sum = dcmplx(0.d0,0.d0)

  ro_mag = sqrt((ro(1)**2)+(ro(2)**2))
  angle1 = theta(ro)   
  p_m1 = 0.5d0/dble(m)
  p_m2 = 1.d0/(ro_mag**m)
  p_m3 = exp(i*dble(m)*angle1)
  phi_min_np = p_m1*p_m2*p_m3   !Phi_Minus Product!

end subroutine phi_minus_np
