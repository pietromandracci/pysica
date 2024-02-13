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

module f_statistics
     
   use f_precision

   use f_constants   

   implicit none

contains
     

   subroutine generate_histogram_mask(mask, bin_limits, n_bins, n_min_data, x_array, n_data, n_max_bins)

      ! Generate an histogram from an array of univariate data, the bins are determined in such a way that
      ! their number does not exceed n_max_bins and each bin, except the last one, contains et least n_min_data
      ! A logical mask is created as a 2-dimensional array, each column of which defines an histogram bin
      ! in each column the true values correspond to the elements of x_array that belong to that bin
      
      !f2py intent(inout)  :: mask, bin_limits, n_min_data
      !f2py intent(out)    :: n_bins
      !f2py intent(in)     :: x_array, n_max_bins, debug
      !f2py intent(hide)   :: n_data

      ! Parameters
      logical,  dimension(n_data, n_max_bins), intent(inout) :: mask        ! Each column of the array refers to a bin and contains
                                                                            !    a logical mask with .true. values
                                                                            !    for the indexes of x_array components
                                                                            !    that belong to that bin
      real(dp), dimension(0:n_max_bins),       intent(inout) :: bin_limits  ! Upper limits of the bins
                                                                            !    except the zero element which is the lower limit
                                                                            !    of the first bin
      integer,                                 intent(out)   :: n_bins      ! Actual number of generated bins 
      integer,                                 intent(inout) :: n_min_data  ! Minimum number of elements in each bin
      integer,                                 intent(in)    :: n_max_bins  ! Maximum number of bins to be generated      
      real(dp), dimension(n_data),             intent(in)    :: x_array     ! Array of data values to be binned
      integer,                                 intent(in)    :: n_data      ! Number of data values

      ! Local variables
      integer, parameter         :: factor = 1000
      logical, dimension(n_data) :: unbinned                             ! True values identify x_array elements
                                                                          ! which do not belong to a bin yet
      real(dp)                   :: x_min, x_max, x_limit, x_delta, bonus
      integer                    :: i_delta, i_bin, i_data, count_data

      ! Check that the minimun number requested for each bin does not exceed the number of data
      if (n_min_data > n_data) n_min_data = n_data

      ! This factor is required to avoid that, due to the cumulation of rounding errors, 
      ! the upper limit of the last bin is lower then the maximum data value
      bonus         = 1.0 + 1.0 / real(factor * n_max_bins)
      x_min         = minval(x_array) 
      x_max         = maxval(x_array) 
      x_delta       = bonus * (x_max - x_min) / real(n_max_bins)
      unbinned      = .true.            ! No data elements have been binned yet
      bin_limits(0) = x_min             ! Lower limit of the first bin
      i_bin         = 1                 ! Bin counter
      count_data    = 0                 ! Counter of data belonging to the bin
      x_limit       = x_min
      do i_delta = 1, n_max_bins
         x_limit = x_limit + x_delta    ! Try this Upper limit for this bin         
         do i_data = 1, n_data
            ! If this datum is in the bin, set the corresponding element of the mask
            if ( (unbinned(i_data)) .and. (x_array(i_data) <= x_limit) ) then
               mask(i_data, i_bin) = .true.
               unbinned(i_data)    = .false.
               count_data          = count_data + 1  ! Number of data in this bin
            endif
         enddo
         if ( (count_data >= n_min_data) .or. (x_limit >= x_max) ) then 
            bin_limits(i_bin) = x_limit
            i_bin             = i_bin + 1
            count_data        = 0
         endif
      enddo
      n_bins = i_bin - 1
      
   end subroutine generate_histogram_mask


   subroutine generate_histogram_mask_reversed(mask, bin_limits, n_bins, n_min_data, x_array, n_data, n_max_bins)

      ! Generate an histogram from an array of univariate data, the bins are determined in such a way that
      ! their number does not exceed n_max_bins and each bin, except the first one, contains et least n_min_data
      ! A logical mask is created as a 2-dimensional array, each column of which defines an histogram bin
      ! in each column the true values correspond to the elements of x_array that belong to that bin
      
      !f2py intent(inout)  :: mask, bin_limits, n_min_data
      !f2py intent(out)    :: n_bins
      !f2py intent(in)     :: x_array, n_max_bins, debug
      !f2py intent(hide)   :: n_data

      ! Parameters
      logical,  dimension(n_data, n_max_bins), intent(inout) :: mask        ! Each column of the array refers to a bin and contains
                                                                            !    a logical mask with .true. values
                                                                            !    for the indexes of x_array components
                                                                            !    that belong to that bin
      real(dp), dimension(0:n_max_bins),       intent(inout) :: bin_limits  ! Upper limits of the bins
                                                                            !    except the zero element which is the lower limit
                                                                            !    of the first bin
      integer,                                 intent(out)   :: n_bins      ! Actual number of generated bins 
      integer,                                 intent(inout) :: n_min_data  ! Minimum number of elements in each bin
      integer,                                 intent(in)    :: n_max_bins  ! Maximum number of bins to be generated      
      real(dp), dimension(n_data),             intent(in)    :: x_array     ! Array of data values to be binned
      integer,                                 intent(in)    :: n_data      ! Number of data values

      ! Local variables
      integer, parameter         :: factor = 1000
      logical, dimension(n_data) :: unbinned                             ! True values identify x_array elements
                                                                          ! which do not belong to a bin yet
      real(dp)                   :: x_min, x_max, x_limit, x_delta, bonus
      integer                    :: i_delta, i_bin, i_data, count_data

      ! Check that the minimun number requested for each bin does not exceed the number of data
      if (n_min_data > n_data) n_min_data = n_data

      ! This factor is required to avoid that, due to the cumulation of rounding errors, 
      ! the upper limit of the last bin is lower then the maximum data value
      bonus         = 1.0 + 1.0 / real(factor * n_max_bins)
      x_min         = minval(x_array) 
      x_max         = maxval(x_array) 
      x_delta       = bonus * (x_max - x_min) / real(n_max_bins)
      unbinned      = .true.            ! No data elements have been binned yet
      bin_limits(0) = x_max             ! Upper limit of the first bin
      i_bin         = 1                 ! Bin counter
      count_data    = 0                 ! Counter of data belonging to the bin
      x_limit       = x_max
      do i_delta = 1, n_max_bins
         x_limit = x_limit - x_delta    ! Try this lower limit for this bin         
         do i_data = 1, n_data
            ! If this datum is in the bin, set the corresponding element of the mask
            if ( (unbinned(i_data)) .and. (x_array(i_data) >= x_limit) ) then
               mask(i_data, i_bin) = .true.
               unbinned(i_data)    = .false.
               count_data          = count_data + 1  ! Number of data in this bin
            endif
         enddo
         if ( (count_data >= n_min_data) .or. (x_limit <= x_min) ) then 
            bin_limits(i_bin) = x_limit
            i_bin             = i_bin + 1
            count_data        = 0
         endif
      enddo
      n_bins = i_bin - 1
      
   end subroutine generate_histogram_mask_reversed   
   

end module f_statistics
