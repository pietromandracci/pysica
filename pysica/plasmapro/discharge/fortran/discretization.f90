! COPYRIGHT (c) 2020-2022 Pietro Mandracci

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

module f_discretization

   use f_precision

   use f_constants

   implicit none

contains

   subroutine generate_trig_table(table, n, debug)

      ! Generate a table of pre-calculated values of the functions sin(x), cos(x) and tan(x)
      ! A pre-defined 2-dimensional array must be passed to the function
      ! The array columns will be in the order: theta/rad; theta/deg; sin(theta), cos(theta), tan(theta)
      ! WARNING: if the array is passed by a python function as a numpy array, it MUST be defined using the flag
      !          order='fortran' (or order='F').  Eg. x = numpy.zeros([n,5], order='F')

      ! Parameters

      integer,  intent(in)    :: n                  ! Number of values in range [0,pi[ for which functions should be
                                                            !   calculated
      real(dp), intent(inout) :: table(1:n, 1:5)    ! Array of calculated values
      logical,  intent(in)    :: debug

      ! Local variables

      integer  :: i
      real(dp) :: delta                             ! Step for the calculation of the trig fuctions
      real(dp) :: theta

      delta = PI / n
      if (debug) then
         print *
         print *, "N        = ", n
         print *, "delta    = ", delta
         print *
      endif
      theta = 0.0
      do i=1, n
         table(i,1) = theta                  ! Angle in radiants
         table(i,2) = theta*180.0/PI         ! Angle in degrees
         table(i,3) = sin(theta)             ! 
         table(i,4) = cos(theta)
         table(i,5) = tan(theta)
         theta = theta + delta
      enddo

      if (debug) then
         do i=1, n   
            print *, table(i,1), int(table(i,2)), table(i,3), table(i,4), table(i,5)
         enddo
      endif
   
   end subroutine generate_trig_table


   real(dp) function trig_function(table, n, angle, input_index, output_index)

      ! Calculate trigonometric functions using a pre-calculated table of their values
      ! The table rows contain: angle in rad, angle in deg, sin(angle), cos(angle), tan(angle)

      ! Parameters

      integer,  intent(in) :: n                ! Number of values in range [0,pi[ for which functions have been
                                               !    calculated
      real(dp), intent(in) :: table(1:n, 1:5)  ! Array of calculated values
      real(dp), intent(in) :: angle            ! Value of the angle given as input, in rad or in deg, depending on
                                               !    input_index
      integer,  intent(in) :: input_index      ! Index that specifies if the angle is given in rad or deg
                                               !    1-> angle in rad; 2-> angle in deg
      integer,  intent(in) :: output_index     ! Index of the function to give as output
                                               !    1-> angle in rad; 2-> angle in deg; 3-> sin; 4-> cos; 5-> tan
      ! Local variables

      integer :: i

      do i=1,n
         if (angle < table(i,input_index)) exit
      enddo
      trig_function = table(i-1, output_index)

   end function trig_function

end module f_discretization
