! module constants
!   implicit none

!   real*8 :: diameter = 1.d0,Lx=64.d0
!   integer ::  n=1024
!   integer, dimension(12) :: seed
!   integer :: seed_val
  
! end module constants
!!!!!____!!!!!!!!!!!!!!!____!!!!!!!!!
!program main
!  use constants
!  implicit none
!  real*8, allocatable :: y(:)
!
!  seed_val=0
!  seed=seed_val
!  allocate(y(2*n))
!  call random_seed(put=seed)  
!
!  Open(101,file='rand_adsorp.dat',status='unknown')
!
!  call rand_seq_adsorp(y)
!  
!
!  close(101)
!  
!end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rand_seq_adsorp(y)
  use main_module
  use random_parameters
  implicit none
  real*8 :: y(2*numpar)
  
  call gen_rand_array(y)
  
end subroutine rand_seq_adsorp
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gen_rand_array(y)
  use main_module
  implicit none
  integer :: i,z,overlap_count0,overlap_count,ii
  real*8 :: y(2*numpar),rand_cood(2)
  logical :: overlap

  overlap_count=0 
  do i = 1,numpar
     write(*,*) 'i=',i

     overlap =.TRUE.
     overlap_count0=-1

     do while(overlap)

        !if(i .gt. 1)then
        !   do ii = 1,i-1
        !      write(101,*)i-1,y((2*ii)-1),y(2*ii),diameter*0.5d0,Lx
        !   end do
        !   write(101,*)
        !end if

        overlap_count0=overlap_count0+1
        call random_cood(rand_cood)
        y((2*i)-1) = rand_cood(1)
        y(2*i) = rand_cood(2)

        !if(i .gt. 1)then
        !   write(101,*)i-1,rand_cood(1),rand_cood(2),diameter*0.5d0,Lx
        !   write(101,*)
        !end if

!!!Avoids check for i=1, forces to run the loop second time
        if(i.gt.1)then 
           call check_overlap(y,i,overlap)
        else
           overlap=.FALSE.
        end if
     end do

     overlap_count=max(overlap_count,overlap_count0)

     !write(*,*)i,'th particle',y((2*i)-1),y(2*i)
     !write(*,*)'overlap_count0=',overlap_count0
     !Write(*,*)

  end do
  write(*,*) 'overlap_count=',overlap_count

end subroutine gen_rand_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_overlap(y,i,overlap)
  use main_module
  implicit none
  integer :: j,i
  real*8 :: rel_cood(2),rel_dist
  real*8 :: y(2*numpar),prev_par_cood(2),new_prev_par_cood(2),new_rel_cood(2)
  character(len=11) :: flag
  logical :: overlap
  
  do j = 1,i-1
        
     rel_cood(1) = y((2*j)-1) - y((2*i)-1)
     rel_cood(2) = y(2*j) - y(2*i)

     call short_dist(rel_cood(1),new_rel_cood(1),1)
     call short_dist(rel_cood(2),new_rel_cood(2),2)
     
!!!rel_dist = sqrt(rel_cood(1)**2 + rel_cood(2)**2)
     rel_dist = sqrt(new_rel_cood(1)**2 + new_rel_cood(2)**2)
     if(rel_dist .lt. dia)then
        overlap=.TRUE.
        return
     end if
     overlap=.FALSE.

  end do

end subroutine check_overlap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Neighbour Finding!
! subroutine short_dist(cood,new_cood)
!   use main_module
!   implicit none
!   real*8 :: cood,cood_rem,new_cood

!   cood_rem = modulo(cood,Lx)

!   if( cood_rem .gt. (0.5d0*Lx) )then
!      new_cood = cood_rem - Lx
!   else
!      new_cood = cood_rem
!   end if

! end subroutine short_dist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine random_cood(rand_cood)
  use main_module
  use random_parameters
  implicit none
  real*8 :: x,rand_cood(2)
  integer :: i
   
  do i = 1,2
     
     call random_number(x)
     rand_cood(i) = x*Lx
     
  end do

end subroutine random_cood
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

