program main
  implicit none
  
  real*8 :: vel,grad_pot
  integer :: i

  open(1,file='test_trailing.dat',status='unknown')
  
  do i=1,10
     call potential_vel(dble(i),grad_pot)
     write(1,*)i,grad_pot
  end do
  
  close(1)
end program main

subroutine potential_vel(x,grad_pot)
  implicit none

  real*8 :: g=1.d0, alpha=1.d0, m=1.d0
  real*8 :: func1,func2,x,grad_pot

  grad_pot = g**2 * func1 * func2
  func1 = ( -1.d0 - (alpha*m*x) )/(x**2)
  func2 = exp(-alpha*m*x)


end subroutine potential_vel
