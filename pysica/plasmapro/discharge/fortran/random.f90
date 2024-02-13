! COPYRIGHT (c) 2020-2024 Pietro Mandracci

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module f_random
     
   use f_precision

   use f_constants   

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

      ! Returns an integer random number between lower and upper, both included  

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


   subroutine random_maxwell_velocity(vx, vy, vz, v, v_mean)

      ! Generate a random velocity vector following maxwell distribution

      !f2py intent(out) :: vx, vy, vz, v
      !f2py intent(in)  :: v_mean

      ! Parameters
      real(dp), intent(out) :: vx, vy, vz, v
      real(dp), intent(in)  :: v_mean           ! Mean speed of the gas molecules

      ! Local variables
      real(dp) :: costheta, sintheta, phi, vxy
      real(sp) :: rand1, rand2

      ! Generate a random velocity module according to maxwell distribution
      v = random_maxwell_speed(v_mean)
      
      ! Generate a random cosinus of theta in [0,1[
      call random_number(rand1)
      costheta = 2 * rand1 - 1
      sintheta = sqrt(1 - costheta * costheta)      

      ! Generate a random phi in [0, 2*PI[
      call random_number(rand2)      
      phi = 2 * PI * rand2

      ! Calculate the velocity components
      vxy = v      * sintheta
      vx  = vxy    * cos(phi)
      vy  = vxy    * sin(phi)
      vz  = v      * costheta
      
   end subroutine random_maxwell_velocity


   real(dp) function random_maxwell_speed(v_mean)

     ! Return a random speed (modulus of the velocity vector)
     ! extracted from an ensamble of particles following the maxwell distribution
     ! for which the mean speed is given 

     ! Parameters
     real(dp), intent(in) :: v_mean ! Mean value of the speed distribution

     ! Local variables
     real(dp) :: y, x, f, xmax, ymax
     real(sp) :: rand_x, rand_y

     xmax = 4.0
     ymax = 4.0 / (SQRT_PI * E_EULER)
     
     do
        call random_number(rand_x)
        call random_number(rand_y)
        x = rand_x * xmax
        y = rand_y * ymax
        f = maxwell_speed_reduced(x)
        if (y <= f) exit
     enddo
     
     random_maxwell_speed = SQRT_PI / 2 * v_mean * x

   end function    

   
   real(dp) pure function maxwell_speed_reduced(x)

      ! Returns the value of the proability density of the Maxwell ditribution
      ! for the adimensional variable x = 2/sqrt(PI) * v / v_mean = sqrt(m/2kT)*v
      ! where v_mean is the mean speed (m, T, k, are mass, temperature and Boltzmann constant)

      ! Parameters
      real(dp), intent(in) :: x

      maxwell_speed_reduced = 4 / SQRT_PI * x*x * exp(-x*x)
      
   end function maxwell_speed_reduced

   
   real(dp) pure function maxwell_vmean(temperature, mass, squared)

      ! Return the mean speed of an ensamble of particles following the maxwell distribution,
      ! optionally, it is given squared to avoid calculation of the square root and then squaring

      ! Parameters
      real(dp), intent(in) :: temperature
      real(dp), intent(in) :: mass
      logical,  intent(in) :: squared

      ! Local variables
      real(dp) :: v2

      v2 = 8 * K_BOLTZMANN * temperature / (PI * mass)
      if (squared) then
         maxwell_vmean = v2
      else
         maxwell_vmean = sqrt(v2)
      endif

   end function maxwell_vmean
        

end module f_random
