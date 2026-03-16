!Comment: Quadrupolar velocity spit out at particle contact i.e. rel. positio!n = 1, is 1. As shown below
!velocity array= -0.49999996108164402        2.0568846035183744E-015  0.49999996108164335        6.8113034229438106E-016

module parameters_main
  implicit none
  !Input values of these constants separately
  integer :: nosteps=80000!99900!1024!2048!4096!8192!16384!32768
  real*8, parameter :: Ca=1.d0, eta=1.d0, vin=1.d0, mu=1.d0, W=1.d0,Q=1.d0
  real*8, parameter :: pp=1.d0, qq=1.d0
  
  !integer ::  m_yuk=2
  !real*8 ::  A=-1.d0, ka=3.d0
  
end module parameters_main

!For animation
subroutine exportbeadcoords(y,dy,n)!added V1,V2
  use velocities
  use main_module
  implicit none
  integer :: i, n
  integer, save :: j
  real*8 :: rad, y(2*n),dy(2*n)
  character(len=10), save :: numeral, prefix
  character(len=20) :: filename
  logical, save :: flag =.TRUE.
  integer :: res1,res2
  
  rad=0.5d0*dia
  if(flag) then
     !Starting Point
     j=0
     prefix='A'
     flag=.FALSE.
  end if
  !!!Commenting out the line below so that file system starts from A00000.dat
  !j=j+1
  write(numeral,'(i5.5)') j
  filename=trim(prefix)//trim(numeral)//'.dat'

  !To only print every 15th file till 1500 & every 50th afterwards!
  res1 = mod(j,5)
  res2 = mod(j,5)
  if(j.ne.1)then
     if(j.le.1500)then
        if(res1.eq.0)then
           
           open(1,file=filename,status='unknown')!!
           do i=1,n,1
              write(1,*) i, y(2*i-1), y(2*i),dy(2*i-1),dy(2*i),rad,Lx,Ly
           end do
           close(1)
           
        end if
     else
        if(res2.eq.0)then
           
           open(1,file=filename,status='unknown')!!
           do i=1,n,1
              write(1,*) i, y(2*i-1), y(2*i),dy(2*i-1),dy(2*i),rad,Lx,Ly
           end do
           close(1)
           
        end if
     end if
  else
     
     open(1,file=filename,status='unknown')!!
     do i=1,n,1
        write(1,*) i, y(2*i-1), y(2*i),dy(2*i-1),dy(2*i),rad,Lx,Ly
     end do
     close(1)
     
  end if
  
  j=j+1
  
              
end subroutine exportbeadcoords
!___________________!
subroutine importbeadcoords(y,dy)  
  use main_module
  implicit none
  integer :: k
  real*8 :: y(2*numpar),dy(2*numpar)
  real*8 :: a1,a2,a3     !dummy variables

  open(2,file='A26520.dat',status='unknown')
  
  do k=1,numpar
     !To restart from already created files
     read(2,*)a1,y(2*k-1), y(2*k),dy(2*k-1),dy(2*k),a2,a3
     !To read initial config from a file
     !read(2,*)y(2*k-1),y(2*k)
  end do
 
  close(2)
 
end subroutine importbeadcoords
!___________________!
! subroutine write_one_file
!   implicit none
!   logical, save :: ifopen =.TRUE.

!   if(ifopen)then
!      open(15,file='particle_pos_evol.dat',status='unknown')
!      ifopen=.FALSE.
!   end if
!   write(15,*)relpos,y
  
! end subroutine write_one_file
!___________________!
program main
  use yukawa_para
  use velocities
  use parameters_main
  use main_module
  use rectangle
  use random_parameters
  use gauss_para
  !use configuration
  implicit none
  integer :: n, l, config,i
  real*8, allocatable :: y(:), yo(:), dy(:)
  real*8 :: t=0.d0 !!! dt= 0.2d0 !!! 0.1d0 !!! 0.05d0!!! 0.025d0 !0.1d0!0.00625d0!0.0125d0!0.00625d0!0.0125d0!0.3992!0.05d0!0.2!0.1!0.05!0.025!0.0125d0!0.00625d0

  
  interface
     subroutine derives(t,y,dy,n)
       implicit none
       integer :: n
       real*8 :: t, Y(n), dY(n)
     end subroutine derives
  end interface

!!!For two particle simulations switch on generate af!!!
  call read_parameters
  call random_seed(put=seed)
  !call generate_LxLy
  write(*,*)'Lx,Ly=',Lx,Ly
  call generate_af
  write(*,*)'area_fraction=',af
  call get_sigma

  write(*,*)'dt=',dt
      
  !Creating File
  !open(3,file='V1-V2-comparison.dat',status='unknown')
  open(2,file='frames.dat',status='unknown')
  write(2,*) nosteps
  allocate(y(2*numpar),yo(2*numpar),dy(2*numpar))

  !!!GAUSS VELOCITY ALLOCATION!!!
  allocate(gauss_vel(2*numpar))

   
!!!  call setIC(y)
!!!  call setIC_sinusoid(y) !!!Use this together with setIC
  !call setSQAR(y,2*numpar)
  ! call rand_seq_adsorp(y)
  !!! call importbeadcoords(y,dy)
  !!! call adhoc_perturb(y) !!! Used to perturb selected particles a little bit more
  !!!call two_particles(y)
  call linear_varying_array(y)
  !!! call test_pos(y)

  !!!GAUSS VELOCITY DISTR!!!
  !!! call gauss_assign
  
  !!! do i=1,numpar
  !!!   write(*,*)gauss_vel(2*i-1),gauss_vel(2*i)
  !!! end do
  
  call exportbeadcoords(y,dy,numpar)

  open(20,file='time_and_positions.dat',status='unknown')
  write(20,*)t,y(1),y(2),y(3),y(4),dy(1),dy(2),dy(3),dy(4)

 
  do l=1,nosteps
!!!write(*,*)'nostep',l!!!!
!!!!call rk4(y,yo,2*numpar,t,dt,derives)
     call euler(y,yo,2*numpar,t,dt,derives)!!!!
     y=yo; t=t+dt
!!!!For outputting velocity!!!!
     
     
     
     !!!y(1)= 3.d0!!!0.d0
     !!!y(2)=5.d0!!!0.d0
     !!!y(3)=7.d0!!!0.1d0*l
     !!!y(4)=5.d0!!!0.d0
     !!!y(5)=3.d0+10.d0
     !!!y(6)=5.d0
     !!!y(7)=7.d0+10.d0
     !!!y(8)=5.d0
     
     call derives(t,y,dy,2*numpar)     
     call exportbeadcoords(y,dy,numpar)
  
     write(20,*)t,y(1),y(2),y(3),y(4),dy(1),dy(2),dy(3),dy(4)
     write(20,*)l,y(1),y(2),y(3),y(4),dy(1),dy(2),dy(3),dy(4)
     !write(*,*)'first_particle_velocity=',dy(1),dy(2)
     !write(*,*)'second_particle_velocity=',dy(3),dy(4)
     !write(*,*)'main_program_vel=',dy(1),dy(3)
     !if(y(3)-y(1) .le. 1.d0)then
     !   STOP
     !end if
     
  enddo
  

  close(2)
  close(20)
  
end program main

!subroutine for RK4
subroutine derives(t,y,dy,nn)
  use velocities
  use parameters_main
  use flag_self_int
  use main_module
  use gauss_para
  implicit none 
  integer :: i,j,nn
  real*8 :: t,y(nn),dy(nn),r1(2),r2(2),V(2),relpos(2)
  
  dy=0.d0
  !!!GAUSS VELOCITY TRANSFER!!!
  !!! call dy_initialize(dy) 
  
  !Calculating velocity
  do i=1,numpar
     !!! dy((2*i)-1)=0.d0!0.75d0!(mu*vin)
     do j=1,numpar

        !write(*,*)'####','particle ij',i,j!!!!
        
        if(i .EQ. j)then
           flag_self = .true.
        else
           flag_self = .false.
        end if

        r1(1)=y(2*i-1)
        r1(2)=y(2*i)
        r2(1)=y(2*j-1)
        r2(2)=y(2*j)
        !relpos(1)=r2(1)-r1(1)
        !relpos(2)=r2(2)-r1(2)
        relpos(1)=r1(1)-r2(1)           
        relpos(2)=r1(2)-r2(2)
        
        pos_diff(1) = relpos(1)
        pos_diff(2) = relpos(2)
        
        call combined_VEL(i,j,V,relpos)
        
        dy(2*i-1)=dy(2*i-1)+V(1)
        dy(2*i)=dy(2*i)+V(2)
                
     enddo
  enddo

!!!  write(*,*)'derives'
!!!  do i=1,numpar
!!!     write(*,*)dy(2*i-1),dy(2*i)
!!!  end do

end subroutine derives

!!!GAUSS VELOCITY TRANSFER
subroutine dy_initialize(dy)
  use gauss_para
  use main_module
  implicit none
  real*8 :: dy(2*numpar)
  integer :: i
  
  do i=1,numpar
     dy(2*i-1)=gauss_vel(2*i-1) * g_stren
     dy(2*i)=gauss_vel(2*i) * g_stren
  end do
!  write(*,*)'inside derives'
!  do i=1,numpar
!     write(*,*)dy(2*i-1),dy(2*i)
!  end do
  
end subroutine dy_initialize

!Gives combined velocity due to Quadrupolar & Yukava potential
subroutine combined_VEL(i,j,V,relpos)
  use velocities
  use parameters_main
  use new_cood_yukawa
  implicit none
  integer :: i,j
!!!real*8 :: V1(2), V2(2),
  real*8 :: V(2),relpos(2),v_s(2),alpha,beta(2)
  alpha = 0.5d0!1.d0
  beta(1) = 1.d0!!! 0.75d0
  beta(2) = 0.d0
  
  !new_relpos(1) = relpos(1)
  !new_relpos(2) = relpos(2)
  !write(*,*)'with_in_combined_VEL'
  call getVsingle(V_d,V_q,relpos)

  !the underlying calls are dirty, make them robust
  call short_dist(relpos(1),new_relpos(1),1)
  call short_dist(relpos(2),new_relpos(2),2)
  !new_relpos(1)=relpos(1)
  !new_relpos(1)=relpos(1)
  !write(*,*)'new_relpos',new_relpos
  
  call yukava_activation(i,j,V2)
  !!!V2=0.d0!!!! 
  !write(*,*)'Repul_Vel=',V2
  !write(*,*)'Hydro_Vel=',V1
!!!SWAPPING TRAJECTORY!!!
  call swap_traj_activation(i,j,v_s)
  !call swap_trajec_vel(relpos,v_s)
  
  !V2=0.d0
  v_s = 0.d0
  !V_d=0.d0
  !V_q=0.d0

!!!MAKING SELF-TERM ZERO!!!
  if(i .EQ. j)then
     V_d = 0.d0
  end if
    
  V(1)= alpha*V_d(1) + (1-alpha)*V_q(1) + V2(1) + (1-alpha)*v_s(1)
!!!v_s added to only x component since along x dir!!!
  V(2)= alpha*V_d(2) + (1-alpha)*V_q(2) + V2(2) + (1-alpha)*v_s(2)
  !write(*,*)'Ewald_vel=',V_q(1),V_q(2)
  !write(*,*)'Swap_traj_vel=',v_s
  
end subroutine combined_VEL

!Neighbour Finding!
subroutine short_dist(cood,new_cood,cn)
  use main_module
  implicit none
  real*8 :: cood,cood_rem,new_cood,LL
  integer :: cn

  if(cn.eq.1) then
     LL=Lx
  else if(cn.eq.2) then
     LL=Ly
  else
     write(*,*) 'cn=',cn,'  something wrong in short_dist'
  end if

  cood_rem = modulo(cood,LL)

  if( cood_rem .gt. (0.5d0*LL) )then
     new_cood = cood_rem - LL
  else
     new_cood = cood_rem
  end if

end subroutine short_dist

!_Activates Yukawa if different particles__!
subroutine yukava_activation(i,j,V2)
  use parameters_main
  use new_cood_yukawa
  implicit none
  
  integer :: i,j
  real*8 ::V2(2),V_yuk(2)

  if(i.eq.j)then
     V2(1) = 0.d0
     V2(2) = 0.d0
  else
     call potential_vel(V_yuk)
     V2(1) = V_yuk(1)
     V2(2) = V_yuk(2)
     !write(*,*)'yuk_inside',V2!!!!
     !STOP!!!!
  end if

  !write(*,*)'Yuk-activation-sub',relpos,V2!!!!
  !write(*,*)
  
end subroutine yukava_activation

!Quadrupole-Hydrodynamics Velocity
subroutine getVsingle(V_m1,V_m2,relpos)
  use parameters_main
  use parameters
  use main_module
  implicit none
  real*8 :: dVx,dVy,V_m1(2),V_m2(2),const,relpos(2),x,y
  complex*16 :: grad_ewald(2)
  

  x=relpos(1)
  y=relpos(2)
  
  !Dipole call
  m=1
  call ewald(relpos,grad_ewald)
  V_m1(1) = (dia**(m+1))*REALPART(grad_ewald(1))
  V_m1(2) = (dia**(m+1))*REALPART(grad_ewald(2))
  
  !!! write(*,*)'m=',m
  !!! write(*,*)'dipolar_getVsingle=',V_m1
  !Quadrupole call
  m=2
  call ewald(relpos,grad_ewald)
  V_m2(1) = (dia**(m+1))*REALPART(grad_ewald(1))
  V_m2(2) = (dia**(m+1))*REALPART(grad_ewald(2))

  !write(*,*)'with_in_getVsingle=',m
  !Incorporating mobility matrix 
  !V(2) = 1.3d0 * (dia**(m+1))*REALPART(grad_ewald(2))

  !write(*,*)'Ewald_vel=',REALPART(grad_ewald(1)),REALPART(grad_ewald(2))
  !write(*,*)'hydro-stuff-relpos,V',relpos,V!!!!
  

end subroutine getVsingle

! !subroutine for allocating intial values for linear arrangement
! subroutine setIC(y)
!   use parameters_main
!   use parameters
!   use main_module
!   implicit none
!   integer mm,z,ii,jj
!   real*8 :: y(numpar)
  
!   do z=1,numpar
!      y((2*z)-1)= 0.d0!!!!*p_sp  Make a similar setup for linear arrangement as you have made for rectangular
!      y(2*z)= (z-1)*dia 
!   enddo

!   !write(*,*)'rad-setIC',dia!!!!
  
! end subroutine setIC

! !To initialize as a square array
! subroutine setSQAR(y,mm)
!   use parameters_main
!   use parameters
!   use main_module
!   use rectangle
!   implicit none
!   integer mm,z,ii,jj,numpar_x,numpar_y
!   real*8 :: y(mm),rad

!   z=1  
!   y=0.d0
!   numpar_x = floor(sqrt(dfloat(numpar)))
!   numpar_y = numpar_x
!   write(*,*)'numpar_x&y',numpar_x,numpar_y

!   do ii=1,numpar_x!!!!floor(sqrt(dfloat(numpar)))
!      do jj=1,numpar_y!!!!floor(sqrt(dfloat(numpar)))

!         y((2*z)-1)=(ii-1)*dia*p_sf_x
!         y(2*z)=(jj-1)*dia*p_sf_y

!         z=z+1
!      enddo
!   enddo

!   call randomize(y)

! end subroutine setSQAR

!__________________________!
subroutine potential_vel(grad_fxy)
  use parameters_main
  use parameters
  use main_module
  use new_cood_yukawa
  implicit none
  real*8 :: x,y,ro,grad_f,grad_fxy(2)
  
  x=new_relpos(1)/(dia) !!!!****
  y=new_relpos(2)/(dia) 
  ro=sqrt((x**2)+(y**2))
  
  call grad_potential(grad_f,x,y,ro)
  call gradxy_potential(grad_fxy,grad_f,x,y,ro)
  !write(*,*)'from pot_vel',grad_fxy!!!!
  
end subroutine potential_vel
  
!__________________________!
subroutine grad_potential(grad_f,x,y,ro)
  use yukawa_para
  use parameters_main
  use parameters
  use main_module
  implicit none
  real*8 :: relpos(2),x,y,fexp,ro,grad_f,fac

  fexp=exp(-ka*ro)
  fac = ((-ka/(ro**m_yuk))-(m_yuk/(ro**(m_yuk+1))))
  call back_calc_A
  !A = -10.d0
  !write(*,*)'A=',A
  grad_f = A*fexp*fac*(1.d0/(dia)) !!!!**** !1/dia due to differentiation chain rule
  

end subroutine grad_potential
!_____________________________!
subroutine back_calc_A
  use main_module
  use yukawa_para
  implicit none
  real*8 :: r_yuk,fac,f0

  r_yuk = 1.d0

!!!NOTE!!
!!!Value of quadrupole is not 1 at x=1, so to compensate for that fact, f0=2 makes summation of quadrupole and yukawa=0 at x=1 in non-periodic systems. In periodic systems, there would be a slight deviation in the numbers as a result the sum=0 at x=1.002 for f0=2.
  f0=2.d0!1.d0 !correction for vaue of quadrupole at 1
  
  fac = f0*((-ka/(r_yuk**m_yuk))-(m_yuk/(r_yuk**(m_yuk+1))))

  A = (dia)/( exp(-ka*r_yuk)*fac ) !!!!****
!!!!!---
  !A = (6.0383273022452132E-002)*A

end subroutine back_calc_A

!_____________________________!
subroutine gradxy_potential(grad_fxy,grad_f,x,y,ro)
  use parameters_main
  implicit none
  real*8 :: grad_f,grad_fx,grad_fy,grad_fxy(2),x,y,ro,unit_vect(2)

  unit_vect(1) = x/ro
  unit_vect(2) = y/ro

  grad_fx = grad_f*unit_vect(1)
  grad_fy = grad_f*unit_vect(2)

  grad_fxy(1) = grad_fx
  grad_fxy(2) = grad_fy

  !write(*,*)'gradxy_pot',x,y,grad_f,grad_fxy!!!!
  
end subroutine gradxy_potential
!_____________________________!

subroutine read_rectangle
  use rectangle
  implicit none
  
  write(*,*)
  write(*,*) 'READING RECTANGLE DATA'
  read(*,*)Nx
  read(*,*)Ny
  write(*,*) 'Nx,Ny=',Nx,Ny

  !p_sf_x = 2.d0!!!!!
  !p_sf_y = 2.d0!!!!!

end subroutine read_rectangle

subroutine get_rectangle
  use main_module
  use rectangle
  implicit none

  write(*,*)
  write(*,*) 'CREATING  RECTANGLE CONFIGURATION'
  numpar = Nx*Ny
  Lx = p_sf_x*dia*Nx + del_x
  Ly = p_sf_y*dia*Ny + del_y
  write(*,*)'Nx-Ny-numpar',Nx,Ny,numpar
  write(*,*) 'lx,Ly=',Lx,Ly

  !call setSQAR(y,2*numpar)
  
end subroutine get_rectangle
  
subroutine randomize(y)
  use main_module
  use random_parameters
  implicit none
  !integer, dimension(12) :: seed
  real*8 :: x
  real*8 :: y(2*numpar)
  integer :: i

  !write(*,*)'randomize_subroutine_seed=',seed
  call random_seed(put=seed) 

  do i=1,2*numpar
     call random_number(x)
!!!!write(*,*)'x=',x
     y(i)=y(i)+(x*factor)
  end do
end subroutine randomize

!Reads remaining parameters from the bash script
subroutine read_parameters
  use yukawa_para
  use random_parameters
  use main_module
  use gauss_para
  !use swap_traj_para
  !use rectangle
  implicit none

!!!for two particle simulations switch off af and switch on Lx,Ly!!!
  
  !read(*,*)af
  read(*,*)Lx
  read(*,*)Ly
  read(*,*)dt
  read(*,*)numpar
  read(*,*)seed_val
  seed = seed_val
  read(*,*)ka
  read(*,*)m_yuk
  !read(*,*)g_stren
  !read(*,*)p_sf_x
  !read(*,*)p_sf_y
  !read(*,*)factor
  !write(*,*)'p_sf_x,p_sf_y,factor=',p_sf_x,p_sf_y,factor
  write(*,*)'numpar=',numpar
  !write(*,*)'area_fraction=',af
  write(*,*)'kappa=',ka
  write(*,*)'m_yuk=',m_yuk
 
end subroutine read_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine generate_LxLy
  use parameters
  use main_module
  implicit none
  real*8 :: num,den

  num = numpar*pi*(dia**2)*0.25d0
  den = af
  Lx = sqrt(num/den)
  Ly = Lx
  
end subroutine generate_LxLy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine generate_af
  use parameters
  use main_module
  real*8 :: num,den


  num = numpar*pi*0.25d0*(dia**2)
  !Ly = Lx
  den = Lx * Ly
  af = num/den
  write(*,*)'From generate_af:Lx=',Lx,'Ly=',Ly
end subroutine generate_af
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_sigma
  use parameters
  use main_module
  implicit none

  sigma = f_sigma*Lx
  write(*,*)'sigma=',sigma

end subroutine get_sigma

!Two particle placement!
subroutine two_particles(y)
  use main_module
  use two_particles_para
  implicit none
  real*8 :: y(2*numpar)
  !real*8 :: drho(2)
  if(numpar.ne.2) then
     write(*,*) 'numpar.ne.2'
     stop
  else
     read(*,*)drho(1)
     read(*,*)drho(2)
     write(*,*)'drho=',drho
     y(1)=(Lx-drho(1))/2.d0
     y(2)=(Ly-drho(2))/2.d0
     y(3)=(Lx+drho(1))/2.d0
     y(4)=(Ly+drho(2))/2.d0
  end if
  write(*,*) 'configuration'
  write(*,*) y

end subroutine two_particles

!subroutine for allocating intial values for linear arrangement
subroutine setIC(y)
  use two_particles_para
  use main_module
  implicit none
  integer i,k,part_pair,imax
  real*8 :: y(numpar)
  k=1
  part_pair = numpar/2

  !drho(1) = 3.d0
  !drho(2) = 0.d0
  read(*,*)drho(1)
  read(*,*)drho(2)
  write(*,*)'particle separations=',drho
  write(*,*)'numpar=',numpar
  
  !do i=1,part_pair+1,2
  !do i=1,numpar-1,2
  !   write(*,*)'i=',i,'k=',k
  !   y(2*i-1) = (Lx - (2*k-1)*drho(1))/2.d0
  !   y(2*i) = (Ly - (2*k-1)*drho(2))/2.d0
  !   y(2*i+1) = (Lx + (2*k-1)*drho(1))/2.d0
  !   y(2*i+2) = (Ly + (2*k-1)*drho(2))/2.d0
  !   k=k+1
  !   
  !enddo

  imax = Lx/drho(1) ! + 2
  do i=1,numpar

     !###First particle at half standard inter-particle separation from cell wall
     y(2*i-1) = dble( (i-1)+0.5)*drho(1)
     y(2*i) = 0.5d0*Ly
     
  end do
  
  !write(*,*)'rad-setIC',dia!!!!
  do i=1,numpar-1
     write(*,*)y(2*i-1),y(2*i) 
  end do
end subroutine setIC


!Linear Chain with 2 varying interparticle gaps!
subroutine linear_varying_array(y_new)
  use main_module
  use two_particles_para
  implicit none
  integer :: i,j,k,l,particle_count,num_extra_spaces
  real*8 :: y(2*numpar),y_new(2*numpar),free_space_x
  integer :: left_lim, middle_lim, right_lim
  particle_count=0
  k=1
  l=1  

  write(*,*)'Creating a linear array of particles around center of cell'
  
  read(*,*)drho(1)
  read(*,*)drho(2)
  write(*,*)"Initial Interparticle distance",drho(1),drho(2)
  
  read(*,*)pair_num
  write(*,*)"Pair_num",pair_num
  read(*,*)numpar_extra
  write(*,*)"Numpar_extra",numpar_extra

!!!Creating an array of particles starting from centre
  !  do i = 1,numpar-numpar_extra-1,2
  !     
  !     particle_count = particle_count + 2
  !     !write(*,*)particle_count
  !     y(2*i-1) = (Lx - (2*k-1)*drho(1))/2.d0
  !     y(2*i) = (Ly - (2*k-1)*drho(2))/2.d0
  !     y(2*i+1) = (Lx + (2*k-1)*drho(1))/2.d0
  !     y(2*i+2) = (Ly + (2*k-1)*drho(2))/2.d0
  !     k=k+1
  !     
  !  end do
  

  do i=1,numpar-numpar_extra
     
     !###First particle at half standard inter-particle separation from cell wall
     y(2*i-1) = dble( (i-1)+0.5)*drho(1)
     y(2*i) = 0.5d0*Ly
     
  end do

  do i=1,numpar-1
     write(*,*)y(2*i-1),y(2*i)
  end do
  
!!!No. of extra particles allowed with given drho & drho_new & dia
  free_space_x = ( (2*pair_num) ) * drho(1)  
  drho_new(1) = free_space_x/(numpar_extra + 2*pair_num )
  drho_new(2) = drho(2)
  
  write(*,*)"Free space created=",free_space_x
  write(*,*)"New particle separation=",drho_new(1)
  
  if(drho_new(1) .lt. 1.d0)then
     write(*,*)"Warning: Particle spacing less than diameter"
     STOP
  end if

!!!Shifting particle position array to accomodate newer particles
  ! numpar = numpar + numpar_extra
!  write(*,*)"No. of extra particles=",numpar_extra
!  do i = 1,(2*pair_num) + numpar_extra - 1,2
!     y_new(2*i-1) = (Lx - (2*l-1)*drho_new(1))/2.d0
!     y_new(2*i) = (Ly - (2*l-1)*drho_new(2))/2.d0
!     y_new(2*i+1) = (Lx + (2*l-1)*drho_new(1))/2.d0
!     y_new(2*i+2) = (Ly + (2*l-1)*drho_new(2))/2.d0
!     l=l+1
     !write(*,*)"new placement",i
!  end do

!!!Populating other particle coordinates into y_new
!  do i = 1,2*(numpar-numpar_extra-2*pair_num)
!     y_new(2*(2*pair_num + numpar_extra) + i) = y(2*(2*pair_num) + i)
!  end do

!  do i = 1, (2*pair_num) + numpar_extra

!     y_new(2*i-1) = y(2*i-1)
     
!  end do

  left_lim = (numpar - 2*pair_num - numpar_extra)/2
  middle_lim = 2*pair_num + numpar_extra
  right_lim = (numpar - 2*pair_num - numpar_extra)/2
  
  do i=1,left_lim
     y_new(2*i-1) = y(2*i-1)
     y_new(2*i) = 0.5d0*Ly
  end do

  do i=left_lim + 1,left_lim + middle_lim
     y_new(2*i-1) = left_lim*drho(1) + dble( (l-1)+0.5)*drho_new(1)
     y_new(2*i) = 0.5d0*Ly
     l = l + 1
  end do

  do i=left_lim + middle_lim + 1,left_lim + middle_lim + right_lim
     y_new(2*i-1) = left_lim*drho(1) + middle_lim*drho_new(1) + dble( (k-1)+0.5)*drho(1)
     y_new(2*i) = 0.5d0*Ly
     k = k + 1
  end do

  do i=1,numpar-1
     write(*,*)y_new(2*i-1),y_new(2*i)
  end do
  
  
end subroutine linear_varying_array

!!!_____________________________________!!!
!subroutine for allocating intial values for sinusoidal arrangement
subroutine setIC_sinusoid(y)
  use two_particles_para
  use main_module
  implicit none
  integer:: i,part_pair
  real*8 :: y(2*numpar),pi = 4.d0 * atan(1.d0)
  real*8 :: k,incr,fac1,fac2,waves
  
  fac1 = 1.d0!0.5d0!1.d0
  fac2 = 1.d0!0.05d0! 0.0125d0! 0.3d0! 0.5d0
  waves = 1.d0!6.d0! 1.d0
  write(*,*)'particle separations=',drho
  write(*,*)'numpar=',numpar

  k=0
  incr=(2*pi*waves)/(numpar)!0.2d0!0.4d0!0.2d0!0.1d0
  do i=1,2*numpar,2

     ! y(i) = y(i) + 1.d0/( 2.d0 + fac1 * ( sin(dble(k)) ) )!##sin density##!
     y(i) = y(i) + fac2 * sin(dble(k) )!+ (pi/2.d0))!##sin position perturbation##!
     
     ! y(2*i-1) = y(2*i-1) + 0.5d0 * sin(2*pi*k)
     ! y(2*i)= y(2*i) !+ sin(2*pi*k*(drho(1)/Lx))
     ! y(2*i+1) = y(2*i+1) + 0.5d0 * sin(2*pi*k)
     ! y(2*i+2) = y(2*i+2)
      
     k=k+incr
     write(*,*)'k=',k
  enddo
  
end subroutine setIC_sinusoid

subroutine adhoc_perturb(y)
  use two_particles_para
  use main_module
  implicit none
  real*8 :: y(2*numpar)

  y(19) = y(19) - 0.176d0
  y(23) = y(23) + 0.176d0
  !y(31) = y(31) - 0.176d0 !0.125d0 !0.25d0 !0.5d0
  !y(35) = y(35) + 0.176d0 !0.125d0 !0.25d0 !0.5d0
  
end subroutine adhoc_perturb

subroutine test_pos(y)
  use main_module
  implicit none
  real*8 :: y(2*numpar)
  
  y(1)=5.d0
  y(2)=5.d0
  
  ! y(3)=5.d0+10.d0
  ! y(4)=5.d0
  
  ! y(5)=5.d0+10.d0
  ! y(6)=5.d0+10.d0
  
  ! y(7)=5.d0
  ! y(8)=5.d0+10.d0
  
  ! !!!!!!!!!right-up-cell
  ! y(9)=3.d0+10.d0
  ! y(10)=5.d0+10.d0
  ! y(11)=7.d0+10.d0
  ! y(12)=5.d0+10.d0
  ! !!!!!!!!!up-cell
  ! y(13)=3.d0
  ! y(14)=5.d0+10.d0
  ! y(15)=7.d0
  ! y(16)=5.d0+10.d0
  
end subroutine test_pos
