program main
   use spsa_routine
   use fdsa_routine
   use other_routines
   implicit none

   ! variables passed on to the spsa and fdsa routines:
   real(kind=8), allocatable, dimension(:) :: guess
   real(kind=8) a, c, alpha, gama ! the gain coefficients
   real(kind=8) :: convergence ! |E_old - E_new|
   real(kind=8) :: B_multiplier ! B = B_multiplier * niter
   integer :: niter ! maximum number of iterations
   integer :: select_function, run_number
   integer :: p ! number of independent variables to be optimized

   ! for determining the number of variables to be optimized by 
   ! reading from the inputfile - read(1,*,iostat=io):
   real(kind=8) :: dummy
   integer :: io

   ! control variables:
   character(len=9) :: random_number_dist
   character(len=4) :: select_method
   character(len=1) :: define_a_using_B ! use stability constant B
!                                                 or not? [Y/N]
   
   ! input/output file names:
   character(len=13) :: run_outputfile
   character(len=7) :: inputname
   character(len=5) :: filename

   ! for a better structured input file:
   ! (these variables will not be used anywhere 
   ! in the program)
   character(len=2) :: name1
   character(len=5) :: name9
   character(len=6) :: name8
   character(len=7) :: name2
   character(len=9) :: name3
   character(len=12) :: name4
   character(len=13) :: name5
   character(len=27) :: name6
   character(len=25) :: name7

   ! counter:
   integer :: i

   ! the program can produce filenames ending in *00-*99. This is good
   ! for systematically testing the effect of different choices of 
   ! coefficients for the algorithms on the convergence of the 
   ! optimization procedure, while saving everything in the same folder 
   ! for easy comparison:
   write(*,*) 'what run is this? [0-99]'
   read(*,*) run_number

   if ( (run_number .gt. 99) .or. (run_number .lt. 0) ) then
      write(*,*) 'invalid run number'
      stop
   end if

   !--------------------------------------
   ! read the initial guess from
   ! the input file
   !--------------------------------------

   ! define the name of the input file that contains the initial guess:
   ! (guess00-guess99)
   if ( (run_number .ge. 0) .and. (run_number .le. 9) ) then
      write (inputname, '(A5,I1,I1)') 'guess', 0, run_number
   else if ( (run_number .ge. 10) .and. (run_number .le. 99) ) then
      write (*, '(A5,I2)') 'guess', run_number
      write (inputname, '(A5,I2)') 'guess', run_number
   end if

   open (unit=1, file=inputname, status='old', action='read')

   ! determine the number of independent variables to be optimized:
   p = 0
   do
      read(1,*,iostat=io) dummy
      if (io/=0) exit
      p = p + 1
      write(*,*) io
   end do
   allocate (guess(p))
   
   ! read the initial guess from the input file:
   rewind(1)
   do i = 1, p
      read(1,*) guess(i)
   end do
   close(1)

   !--------------------------------------
   ! read the coefficients and the 
   ! control variables for the 
   ! procedures from the input file:
   !--------------------------------------

   ! define the name of the input file that contains the coefficients
   ! and the control variables: 
   ! (coeff00-coeff99)
   if ( (run_number .ge. 0) .and. (run_number .le. 9) ) then
      write (inputname, "(A5,I1,I1)") 'coeff', 0, run_number
   else if ( (run_number .ge. 10) .and. (run_number .le. 99) ) then
      write (inputname, "(A5,I2)") 'coeff', run_number
   end if

   ! read from that file:
   open (unit=1, file=inputname, status='old', action='read')

   read(1,*) name1, a
   read(1,*) name1, c 
   read(1,*) name2, niter 
   read(1,*) name2, select_method
   if ( (select_method .ne. 'spsa') .and. (select_method .ne. 'fdsa') ) then
      write(*,*) ' no valid optimization procedure selected'
      stop
   end if
   read(1,*) name3, select_function
   read(1,*) name4, convergence 
   if ( select_method .eq. 'spsa') then
      read(1,*) name5, B_multiplier
      read(1,*) name6, random_number_dist 
      read(1,*) name7, define_a_using_B
      read(1,*) name8, alpha
      read(1,*) name9, gama
   end if
   close (1)

   !-------------------------------------
   ! call optimization procedure:
   !-------------------------------------

   ! define the name of the output file where the gradient at each
   ! iteration, as well as the total number of iterations
   ! and the final value L_new - L_old are written:
   if ( (run_number .ge. 0) .and. (run_number .le. 9) ) then
      write (run_outputfile, "(A3,I1,I1,A8)") 'run', 0, run_number, &
         '_details'
   else if ( (run_number .ge. 10) .and. (run_number .le. 99) ) then
      write (run_outputfile, "(A3,I2,A8)") 'run', run_number, '_details'
   end if
   
   open (unit=9, file=run_outputfile, status='new', action='write')

   write(9,*) 'k   gradient vector, at each step of the optimization (k):'
   write(9,*)

   ! define the name of the outputfile that can be used to plot
   ! the result of the optimization procedure:
   if ( (run_number .ge. 0) .and. (run_number .le. 9) ) then
      write (filename, "(A1,I1,A1,I1,I1)") 'f', select_function, '_', 0, &
         run_number
   else if ( (run_number .ge. 10) .and. (run_number .le. 99) ) then
      write (filename, "(A1,I1,A1,I2)") 'f', select_function, '_', &
         run_number
   end if

   ! call the appropriate optimization procedure
   if ( select_function == 1 ) then
      if ( p .ne. f1_nvar()) then
         write(*,*) 'initial guess does not match number of variables ', &
            'for f1'
         stop
      end if
      if ( select_method == 'spsa' ) then
         call spsa (f1, guess, random_number_dist, filename, a, c, alpha, &
            gama, convergence, niter, define_a_using_B, B_multiplier)
      else if ( select_method == 'fdsa' ) then
         call fdsa (f1, guess, filename, a, c, convergence, niter)
      end if
   else if ( select_function == 2 ) then
      if ( p .ne. f2_nvar()) then
         write(*,*) 'initial guess does not match number of variables ', &
            'for f1'
         stop
      end if
      if ( select_method == 'spsa' ) then
         call spsa (f2, guess, random_number_dist, filename, a, c, alpha, &
            gama, convergence, niter, define_a_using_B, B_multiplier)
      else if ( select_method == 'fdsa' ) then
         call fdsa (f2, guess, filename, a, c, convergence, niter)
      end if
   else
      write(*,*) 'no valid function selected'
      stop
   end if

   close(9)
   deallocate( guess )
end program main
