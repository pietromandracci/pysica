!	COPYRIGHT 2008 Pietro Mandracci
!
!	This program is free software: you can redistribute it and/or modify
!	it under the terms of the GNU General Public License as published by
!	the Free Software Foundation, either version 3 of the License, or
!	(at your option) any later version.
!
!	This program is distributed in the hope that it will be useful,
!	but WITHOUT ANY WARRANTY; without even the implied warranty of
!	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!	GNU General Public License for more details.
!
!	You should have received a copy of the GNU General Public License
!	along with this program.  If not, see <http://www.gnu.org/licenses/>.


module f_random

   use f_precision

   implicit none

contains

   subroutine init_random_seed()

      ! Initialises the random number generator, using a seed based on the system time
      ! This subroutine was copied from the GNU95 compiler documentation for the
      ! RANDOM_SEED function

      integer              :: i, n, clock
      integer, allocatable :: seed(:)
         
      call random_seed(size = n)
      allocate(seed(n))
          
      call system_clock(COUNT=clock)
          
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(PUT = seed)
          
      deallocate(seed)

   end subroutine init_random_seed


   integer function random_integer(lower, upper)

      ! Returns an integer random number bewtween lower and upper, both included  

      !f2py intent(in) :: lower, upper

      ! Parameters
      integer, intent(in) ::  lower, upper

      ! Local variables
      real(sp) :: rand01

      call random_number(rand01)

      random_integer = lower + nint( rand01 * (upper - lower) )
   
   end function random_integer


   subroutine random_integers_array(array, n_elements, lower, upper, check_repetitions)

      ! Returns an array of integer random numbers, all of them between lower and upper, both included  
      ! If requested, repetitions are avoided. 
      ! However, if the range of integers (upper-lower+1) is lower than the number of array elements,
      ! it is impossible to avoid repetitions, so the check is automatically turned off

      !f2py intent(hide)  :: n_elements
      !f2py intent(in)    :: lower, upper, check_repetitions
      !f2py intent(inout) :: array

      ! Parameters
      integer,   intent(inout) :: array(1:n_elements)
      integer,   intent(in)    :: lower, upper, n_elements
      logical,   intent(in)    :: check_repetitions

      ! Local variables
      integer :: i, j
      integer :: r 
      logical :: check, duplicate, regenerate


      ! If the range of integers (upper-lower) is lower than the number of array elements,
      ! it is impossible to avoid repetitions, so the check is turned off
      check = .true.
      if (check_repetitions) then
         if (abs(upper-lower+1) < n_elements) check = .false.
      else
         check = .false.
      endif

      ! generate the first index as a random integer in range lower, upper
      r = random_integer(lower,upper)
      array(1) = r
      ! starting from the second index, if required, check that duplicates are avoided
      do i = 2, n_elements
         ! generate a random integer in range lower, upper
         r = random_integer(lower,upper)
         ! if check is required, regenerate the random integer until it is different from all the previous ones
         if (check) then
            regenerate = .true.
            do while (regenerate)
               ! compare random number with all previously generated to check duplicates
               duplicate = .false.
               do j = 1, i-1                          
                  if (r == array(j)) then 
                     duplicate = .true.
                     exit
                  endif
               enddo
               ! in case of duplication, regenerate and recheck
               if (duplicate) then
                  call no_operation !necessary for unknown reason !!!!!!!!!!!
                  r = random_integer(lower,upper)
               else
                  regenerate = .false.
               endif
            enddo ! while (regenerate)
         endif ! (check)
         array(i) = r
      enddo 

   end subroutine random_integers_array


   subroutine no_operation

   end subroutine no_operation

end module f_random
