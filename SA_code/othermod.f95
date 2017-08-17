module other_routines
implicit none
contains

!----------------------------
! Functions to optimize:
!----------------------------

! function 1:
real(kind=8) function f1 (var)
   real(kind=8), allocatable, dimension (:), intent(in) :: var
   real(kind=8) :: x, y
   x = var(1)
   y = var(2)
   f1 = ( sin(x)**2 + cos(x)**2 ) / ( 5.0d0 + x**2 + y**2 )
end function f1

! number of variables in function 1:
integer function f1_nvar ()
   f1_nvar = 2
end function f1_nvar

! function 2:
real(kind=8) function f2 (var)
   real(kind=8), allocatable, dimension (:), intent(in) :: var
   real(kind=8) :: x, y
   x = var(1)
   y = var(2)
   f2 =-(y + 47.0d0)*sin(sqrt(abs(y + x/2.0d0 + 47.0d0))) - x*sin(sqrt(abs(x - (y + 47.0d0)))) 
end function f2

! number of variables in function 2:
integer function f2_nvar ()
   f2_nvar = 2
end function f2_nvar

!-----------------------------------
! Random number distributions:
!-----------------------------------

subroutine bernoulli_distribution (del_k)
   real(kind=8), allocatable, dimension(:), intent(inout) :: del_k
   integer :: i, n

   n = size(del_k)
   do i = 1, n
      call random_number(del_k(i))
      if (del_k(i) .gt. 0.5) then
         del_k(i) = 1.0d0
      else if (del_k(i) .le. 0.5) then
         del_k(i) = -1.0d0
      end if
   end do
end subroutine bernoulli_distribution

subroutine split_uniform (del_k)
   real(kind=8), allocatable, dimension(:), intent(inout) :: del_k
   integer :: i, n

   ! mean magnitude = 0.95 
   ! variance = 1
   ! del_k(i) e [-1.4908,-0.4092] ^ [0.4092,1.4908]

   n = size(del_k)
   i = 1
   do while (i .le. n)
      call random_number(del_k(i))
      del_k(i) = (del_k(i)*2.0d0*1.4908d0) - 1.4908d0
      if (dabs(del_k(i)) .ge. 0.4092d0) then
         i = i + 1
      end if
   end do
end subroutine split_uniform

end module other_routines
