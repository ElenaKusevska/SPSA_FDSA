
!-----------------------------------------
! Implementation of the SPSA algorithm:
!-----------------------------------------

module spsa_routine
use other_routines
implicit none
contains

subroutine spsa (f, theta_k, random_number_distribution, filename, &
      a, c, alpha, gama, convergence, niter, define_a_using_B, B_multiplier)
   
   ! interface for the abstract function:
   interface
      real(kind=8) function f(var)
         real(kind=8), allocatable, dimension (:), intent(in) :: var
         real(kind=8) :: x, y
      end function
   end interface

   ! variables for the spsa algorithm:
   real(kind=8), allocatable, dimension(:), intent(inout) :: theta_k
   real(kind=8), intent(inout) :: a
   real(kind=8), intent(in) :: c, convergence
   real(kind=8), intent(in) :: B_multiplier ! B = B_multiplier * niter
   real(kind=8), intent(in) :: alpha, gama
   character(len=9), intent(in) :: random_number_distribution
   character(len=1), intent(in) :: define_a_using_B
   integer, intent(in) :: niter
   real(kind=8), allocatable, dimension(:) :: del_k
   real(kind=8), allocatable, dimension(:) :: coords_plus, coords_minus
   real(kind=8), allocatable, dimension(:) :: gradient
   real(kind=8) :: ak, ck, B, mfact, test
   real(kind=8) :: L_plus, L_minus, L_new, L_old ! objective function
   integer :: p ! number of varibles to be optimized
   integer :: k ! counter for the steps in the spsa procedure

   ! for seeding the random number generator:
   integer, allocatable, dimension(:) :: seed
   integer :: ns

   ! for naming the output files:
   character(len=5), intent(in) :: filename
   character(len=10) :: filename1

   !counter:
   integer :: i

   !---------------------------------
   ! Initial operations:
   !---------------------------------

   B = B_multiplier*dble(niter)
   if ( define_a_using_B .eq. 'Y') then
      a = a * (( B + 1.0d0)**alpha)
   end if

   ! allocate the random perturbation vector, the vectors to hold the two 
   ! sets of perturbed coordinates, and the gradient vector:
   p = size(theta_k)
   allocate (del_k(p), coords_plus(p), coords_minus(p), gradient(p))

   ! initialize the random number generator:
   call random_seed(size=ns)
   allocate (seed(ns))
   do i = 1, ns
      seed(i) = i
   end do
   write(*,*) "seed:", seed
   call random_seed(put=seed)

   !------------------------
   ! The SPSA peocedure:
   !------------------------

   open (unit=2, file=filename, status='new', action='write')

   k = 1 ! counter
   test = 1.0d0
   L_old = f(theta_k)
   write (2,*) theta_k, L_old

   do while (test .gt. convergence)
   
      ! compute the k-th gain coefficients (step size)
      ak = a / ((B + k + 1)**alpha)
      ck = c / ((k + 1)**gama)

      ! compute the simultaneous perturbation vector:
      if (random_number_distribution == 'bern_dist') then
         call bernoulli_distribution (del_k)
      else if (random_number_distribution == 'split_uni') then
         call split_uniform (del_k)
      else
         write(*,*) 'no distribution function selected'
         stop
      end if

      ! compute the perturbed coordinates:
      do i = 1, p
         coords_plus(i) = theta_k(i) + ck*del_k(i)
         coords_minus(i) = theta_k(i) - ck*del_k(i)
      end do
   
      ! evaluate the function with the perturbed coordinates:
      L_plus = f(coords_plus)
      L_minus = f(coords_minus)

      ! determine the gradient:
      mfact = (L_plus - L_minus) / (2*ck) ! multiplicative factor
      do i = 1, p
         gradient(i) = mfact * (1/del_k(i))
      end do
      write(9,*) k, gradient

      ! update the theta estimate:
      do i = 1, p
         theta_k(i) = theta_k(i) - ak*gradient(i)
      end do
      L_new = f(theta_k)
      write (2,*) theta_k, L_new

      ! test for convergence, or termination
      if (k .ge. (niter + niter/3)) exit
      test = dabs(L_new - L_old)
      L_old = L_new
      k = k + 1
   
   end do
   write(9,*)
   write(9,*) "  k  ", " test "
   write(9,*) k, test

   ! write last point in optimization to a separate file
   filename1 = filename // '_last'
   open (unit=2, file=filename1, status='new', action='write')
   write(2,*) theta_k, L_new
   close(2)

   deallocate (del_k, coords_plus, coords_minus, gradient, seed)

end subroutine spsa

end module spsa_routine
