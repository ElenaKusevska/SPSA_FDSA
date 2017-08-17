
!-----------------------------------------
! Implementation of the fdSA algorithm:
!-----------------------------------------

module fdsa_routine
implicit none
contains

subroutine fdsa (f, theta_k, filename, a, c, convergence, niter)
   
   ! interface for the abstract function:
   interface
      real(kind=8) function f(var)
         real(kind=8), allocatable, dimension (:), intent(in) :: var
         real(kind=8) :: x, y
      end function
   end interface

   ! variables for the fdsa algorithm:
   real(kind=8), intent(inout), allocatable, dimension(:) :: theta_k
   real(kind=8), allocatable, dimension(:) :: gradient
   real(kind=8), allocatable, dimension(:) :: coords_plus, coords_minus
   real(kind=8), intent(in) :: a, c, convergence
   real(kind=8) :: L_plus, L_minus, L_new, L_old, ck, ak, test
   integer, intent(in) :: niter
   integer :: p ! number of varibles to be optimized
   integer :: k ! counter for the steps in the iteration
   integer :: i, j ! counters

   ! for naming the output files:
   character(len=5), intent(in) :: filename
   character(len=10) :: filename1
   
   open (unit=3, file=filename, status='new', action='write')
   
   p = size(theta_k)
   allocate (gradient(p), coords_plus(p), coords_minus(p))
   
   L_old = f(theta_k)
   write(3,*) theta_k, L_old
   k = 1 ! counter
   test = 1

   !compute:
   do while (test .gt. convergence)

      ! compute the k-th gain coefficients (step size):
      ck = c / ((dble(k))**(1.0d0/6.0d0)) ! c = x2 - x1
      ak = a / (dble(k)) ! SA coefficient

      ! perturb the coordinates one by one:
      do j = 1, p
         do i = 1, p
            coords_plus(i) = theta_k(i)
            coords_minus(i) = theta_k(i)
         end do
         coords_plus(j) = theta_k(j) + ck
         coords_minus(j) = theta_k(j) - ck

         ! evaluate the function with the j-th coordinate perturbed:
         L_plus = f(coords_plus)
         L_minus = f(coords_minus)

         ! compute the j-th component of the gradient:
         gradient(j) = (L_plus - L_minus) / (2*ck)
      end do
      write(9,*) k, gradient

      ! update the theta estimate:
      do i = 1, p
         theta_k(i) = theta_k(i) - ak*gradient(i)
      end do
      L_new = f(theta_k)
      write(3,*) theta_k, L_new
         
      ! test for convergence, or termination
      if (k .ge. niter) exit
      test = dabs(L_new - L_old)
      L_old = L_new
      k = k + 1
   
   end do
   write(9,*)
   write(9,*) "  k  ", " test "
   write(9,*) k, test

   close(3)

   ! write last point in optimization to a separate file
   filename1 = filename // '_last'
   open (unit=2, file=filename1, status='new', action='write')
   write(2,*) theta_k, L_new
   close(2)

   deallocate(gradient, coords_plus, coords_minus)

end subroutine fdsa

end module fdsa_routine
