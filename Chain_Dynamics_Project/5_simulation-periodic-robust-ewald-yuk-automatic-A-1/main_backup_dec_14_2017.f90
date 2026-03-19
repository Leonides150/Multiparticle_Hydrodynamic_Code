
module parameters_main
  implicit none
  !Input values of these constants separately
  integer :: nosteps=10!99000
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
  
  rad=0.5d0*dia
  if(flag) then
     !Starting Point
     j=0
     prefix='A'
     flag=.FALSE.
  end if
  j=j+1
  write(numeral,'(i5.5)') j
  filename=trim(prefix)//trim(numeral)//'.dat'
  open(1,file=filename,status='unknown')
  do i=1,n,1
     write(1,*) i, y(2*i-1), y(2*i),dy(2*i-1),dy(2*i),rad,Lx     
  end do
  close(1)
end subroutine exportbeadcoords
!___________________!
subroutine importbeadcoords(y,dy)  
  use main_module
  implicit none
  integer :: k
  real*8 :: y(2*numpar),dy(2*numpar)
  real*8 :: a1,a2,a3     !dummy variables

  open(2,file='A04531.dat',status='unknown')
  
  do k=1,numpar
     read(2,*)a1,y(2*k-1), y(2*k),dy(2*k-1),dy(2*k),a2,a3
  end do
 
  close(2)
 
end subroutine importbeadcoords
!___________________!
program main
  use yukawa_para
  use velocities
  use parameters_main
  use main_module
  use rectangle
  use random_parameters
  !use configuration
  implicit none
  integer :: n, l, config
  real*8, allocatable :: y(:), yo(:), dy(:)
  real*8 :: dt=0.00625d0, t=0.d0
  
  interface
     subroutine derives(t,y,dy,n)
       implicit none
       integer :: n
       real*8 :: t, Y(n), dY(n)
     end subroutine derives
  end interface

  call read_parameters
  call random_seed(put=seed)
  call generate_LxLy
  write(*,*)'Lx,Ly=',Lx,Ly
  call get_sigma
      
  !Creating File
  !open(3,file='V1-V2-comparison.dat',status='unknown')
  open(2,file='frames.dat',status='unknown')
  write(2,*) nosteps
  allocate(y(2*numpar),yo(2*numpar),dy(2*numpar))

  !call setIC(y)
  !call setSQAR(y,2*numpar)
  call rand_seq_adsorp(y)
  !call importbeadcoords(y,dy)

 
  do l=1,nosteps
!!!!write(*,*)'nostep',l!!!!
!!!!call rk4(y,yo,2*numpar,t,dt,derives)
     call euler(y,yo,2*numpar,t,dt,derives)!!!!
     
!!!!For outputting velocity!!!!
     call derives(t,yo,dy,2*numpar)
     call exportbeadcoords(y,dy,numpar)
     !write(3,*)l,pos_diff,V2,V1
     y=yo; t=t+dt
     
  enddo
  
  !close(3)
  close(2)
  !close(4)
  
end program

!subroutine for RK4
subroutine derives(t,y,dy,nn)
  use velocities
  use parameters_main
  use flag_self_int
  use main_module
  implicit none 
  integer :: i,j,nn
  real*8 :: t,y(nn),dy(nn),r1(2),r2(2),V(2),relpos(2)
  dy=0.d0

  !Calculating velocity
  do i=1,numpar
     !dy((2*i)-1)=(mu*vin)
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

end subroutine derives

!Gives combined velocity due to Quadrupolar & Yukava potential
subroutine combined_VEL(i,j,V,relpos)
  use velocities
  use parameters_main
  use new_cood_yukawa
  implicit none
  integer :: i,j
!!!real*8 :: V1(2), V2(2),
  real*8 :: V(2), relpos(2)
  
  !new_relpos(1) = relpos(1)
  !new_relpos(2) = relpos(2)
  call getVsingle(V1,relpos)

  call short_dist(relpos(1),new_relpos(1))
  call short_dist(relpos(2),new_relpos(2))
  !new_relpos(1)=relpos(1)
  !new_relpos(1)=relpos(1)
  !write(*,*)'new_relpos',new_relpos

call yukava_activation(i,j,V2)
  !!!V2=0.d0!!!!
  write(*,*)'Repul_Vel=',V2
  write(*,*)'Hydro_Vel=',V1
  
  V(1)=V1(1)+V2(1)
  V(2)=V1(2)+V2(2)

end subroutine combined_VEL

!Neighbour Finding!
subroutine short_dist(cood,new_cood)
  use main_module
  implicit none
  real*8 :: cood,cood_rem,new_cood

  cood_rem = modulo(cood,Lx)

  if( cood_rem .gt. (0.5d0*Lx) )then
     new_cood = cood_rem - Lx
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
subroutine getVsingle(V,relpos)
  use parameters_main
  use parameters
  use main_module
  implicit none
  real*8 :: dVx,dVy,V(2),const,relpos(2),x,y
  complex*16 :: grad_ewald(2)
  

  x=relpos(1)
  y=relpos(2)
  
  call ewald(relpos,grad_ewald)
  V(1) = (dia**(m+1))*REALPART(grad_ewald(1))
  V(2) = (dia**(m+1))*REALPART(grad_ewald(2))

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
  
  x=new_relpos(1)/dia 
  y=new_relpos(2)/dia 
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
  grad_f = A*fexp*fac*(1.d0/dia)  !1/dia due to differentiation chain rule
  

end subroutine grad_potential
!_____________________________!
subroutine back_calc_A
  use main_module
  use yukawa_para
  implicit none
  real*8 :: r_yuk,fac

  r_yuk = 1.d0
  
  fac = ((-ka/(r_yuk**m_yuk))-(m_yuk/(r_yuk**(m_yuk+1))))

  A = dia/( exp(-ka*r_yuk)*fac )
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
  !use rectangle
  implicit none

  read(*,*)af
  read(*,*)numpar
  read(*,*)seed_val
  seed = seed_val
  read(*,*)ka
  read(*,*)m_yuk
  !read(*,*)p_sf_x
  !read(*,*)p_sf_y
  !read(*,*)factor
  !write(*,*)'p_sf_x,p_sf_y,factor=',p_sf_x,p_sf_y,factor
  write(*,*)'numpar=',numpar
  write(*,*)'area_fraction=',af
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
  implicit none
  real*8 :: y(2*numpar)

end subroutine two_particles
