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

module f_math

   use f_precision

!   use f_constants

!   use f_io

   implicit none

contains

   subroutine moving_average_weighted(a, a_smooth, weights, na, nw, fill, skip, normalize)

      ! Returns the moving average of an array, using a window of 2w+1 points

      !f2py intent(in)   :: a, weights, skip, normalize
      !f2py intent(out)  :: a_smooth     
      !f2py intent(hide) :: na, nw

      ! Parameters
      real(dp), dimension(1:na),  intent(in)             :: a            ! Array of which the moving average is requested
      real(dp), dimension(1:na),  intent(out)            :: a_smooth     ! Smoothed array
      real(dp), dimension(1:nw),  intent(in)             :: weights      ! Weights
      integer,                    intent(in)             :: na
      integer,                    intent(in)             :: nw
      real(dp),                   intent(in),  optional  :: fill         ! Number used to fill the unsmoothed
      integer,                    intent(in),  optional  :: skip         ! Number of elements of the array to be skipped    
      logical,                    intent(in),  optional  :: normalize
      

      ! Local variables
      logical                                 :: norm  
      real(dp), dimension(1:nw)               :: w_norm           ! Normalized weights      
      real(dp)                                :: c_norm           ! Normalization coefficient
      real(dp)                                :: smooth           ! Value of smoothed array element
      real(dp)                                :: f                ! Fill value
      integer                                 :: window, i, j, sk

      if (present(fill)) then
         f = fill
      else
         f = 0
      endif
      
      if (present(normalize)) then
         norm = normalize
      else
         norm = .false.
      endif

      if (present(skip)) then
         sk = skip
      else
         sk = 0
      endif

      window = nw - 1 
      a_smooth = f      
      
      if (norm) then
         c_norm = 0
         do i = 1, nw
            c_norm = c_norm + weights(i)
         enddo
         c_norm = 2 * c_norm - weights(1)
         w_norm = weights / c_norm
      else
         w_norm = weights         
      endif

      ! If the array length is not enough to calculate the average for at least one point, return
      if (na < 2 * window + 1) return

      do i = window + 1 + sk, na - window - sk
         smooth = 0.0
         do j = i-window, i+window
            smooth = smooth + a(j) * w_norm(abs(j-i)+1)
         enddo
         a_smooth(i) = smooth
      enddo
      
   end subroutine moving_average_weighted


end module f_math
