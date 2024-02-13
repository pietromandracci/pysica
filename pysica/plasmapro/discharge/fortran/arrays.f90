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

module f_arrays

   use f_precision

   implicit none

contains

!   subroutine tofalse_minima(data, mask, factor, n_data)

      ! Given the data array and an associated logical mask, set to false the number of mask elements
      ! required to have a number of true elements that is a multiple of factor
      ! the removed elements are the ones with the lower values of corresponding data elements

      ! Parameters
!      real(dp), dimension(n_data), intent(in)    :: data
!      logical,  dimension(n_data), intent(inout) :: mask
!      integer,                     intent(in)    :: factor
!      integer,                     intent(in)    :: n_data

      ! Local variables
!      integer :: n_remove, i, j

!      n_remove = mod(n_data, factor)
!      if (n_remove > 0) then
!         do i = 1, n_remove
!            j = minloc(data, mask)
!            mask(j) = .false.
!         enddo
!      endif
   
!      end subroutine tofalse_minima

   real(dp) function array_average_masked(array, mask, n_rows, n_elements, row)

      ! Returns the average of a row of a two dimensional array using a mask
      ! only the array elements corresponding to true elements of the mask are counted

      !f2py intent(in)   :: array, mask, row
      !f2py intent(hide) :: n_rows, n_elements

      ! Parameters
      real(dp), intent(in) :: array(1:n_rows, 1:n_elements)
      logical,  intent(in) :: mask(1:n_rows,  1:n_elements)
      integer,  intent(in) :: n_rows, n_elements, row

      ! Local variables
      real(dp) :: array_sum, array_size
      integer  :: i

      if ( (row >= 1) .and. (row <= n_rows) ) then
         array_size = 0
         array_sum  = 0
         do i = 1, n_elements
            if (mask(row,i)) then
               array_size = array_size + 1
               array_sum  = array_sum  + array(row,i)
            endif
         enddo
         if (array_size > 0) then
            array_average_masked = array_sum / array_size
         else
            array_average_masked = 0.0_dp
         endif
      else 
         array_average_masked = 0.0_dp
      endif

   end function array_average_masked


   integer function search_index(array, n_elements, value, nearest)

      ! Returns the index of the first array element that is greater than the value provided
      ! if the value provided is greater than the maximum of the array, the index of the latter is returned
      ! if nearest is set to .true. then instead of the index of the first greater value, 
      ! the index of the nearest value is provided

      !f2py intent(hide) :: n_elements
      !f2py intent(in)   :: value, array, nearest
      
      ! Parameters
      real(dp), intent(in) :: array(1:n_elements), value
      integer,  intent(in) :: n_elements
      logical,  intent(in) :: nearest

      ! Local variables
      integer :: array_index

      do array_index = 1, n_elements
         if (array(array_index) > value) exit
      enddo

      if (array_index > n_elements) then 
         search_index = n_elements
      else
         if ((nearest) .and. (array_index > 1)) then
            if ( (array(array_index) - value) <= (value - array(array_index-1)) ) then
               search_index = array_index
            else
               search_index = array_index -1 
            endif
         else
            search_index = array_index
         endif !((nearest) .and. (array_index > 1))
      endif !(array_index > n_elements)

   end function search_index


   integer function n_state_elements(logical_array, n_rows, n_elements, row, logical_state)

      ! Returns the number of true (or false) elements of a row of a logical bidimensional array
      
      !f2py intent(hide) :: n_rows, n_elements
      !f2py intent(in)   :: row, logical_array, logical_state
 
      ! Parameters
      logical, intent(in) :: logical_array(1:n_rows, 1:n_elements)
      logical, intent(in) :: logical_state
      integer, intent(in) :: row, n_rows, n_elements

      ! Local variables 
      integer  :: i, n_count

      if ((row > 0) .and. (row <= n_rows)) then 
         n_count = 0
         do i = 1, n_elements
            if (logical_array(row, i).eqv.logical_state) n_count = n_count + 1
         enddo
         n_state_elements = n_count
      else
         n_state_elements = 0
      endif

   end function n_state_elements


   subroutine find_absolute_indexes(array, n_rows, n_elements, row, relative_indexes, absolute_indexes, n_indexes, state_tofind)

      ! Find the absolute indexes of .true. (.false.) items of a row of a logical bidimensional array
      ! that are identified by relative indexes, i.e. counting only the .true. (.false.) items
      ! Example
      ! given the array [T,F,T,T,F,F,F,T,T,F,T,F,T,F,T]
      ! the absolute indexes corresponding to the relative indexes 1,3,7 of .true. items will be 
      ! 1,4,13
      ! the absolute indexes corresponding to the relative indexes 1,3,7 of .false. items will be 
      ! 2,6,14

      !f2py intent(hide) :: n_rows, n_elements, n_indexes
      !f2py intent(in)   :: array, row, relative_indexes, state_tofind
      !f2py intent(out)  :: absolute_indexes

      ! Parameters
      logical, intent(in)  :: array(1:n_rows, 1:n_elements)  
      integer, intent(in)  :: relative_indexes(1:n_indexes)  
      integer, intent(out) :: absolute_indexes(1:n_indexes)  
      integer, intent(in)  :: n_rows, n_elements, n_indexes
      integer, intent(in)  :: row                      
      logical, intent(in)  :: state_tofind                   ! logical state of the elements to be found

      ! Local variables
      integer :: i, j, k

      if ((row > 0).and.(row <= n_rows)) then            ! if the selected row is out of range do nothing
         j = 1                                           ! j counts the .true. (.false.) items of the array only (relative index)
         do i = 1, n_elements                            ! i counts all the array items (absolute index)
            if (array(row, i).eqv.state_tofind) then     ! if the item is .true. (.false.) increase only the absolute index
               do k = 1, n_indexes                       ! if the relative index equals one of the set given, store the absolute index
                  if (j == relative_indexes(k)) then 
                     absolute_indexes(k) = i
                     exit
                  endif 
               enddo ! k = 1, n_indexes
               j = j + 1
            endif
         enddo ! i = 1, n_items
      endif
 
   end subroutine find_absolute_indexes


   subroutine set_elements(array, n_rows, n_elements, row, indexes_array, n_indexes, newstate)

      ! Set the logical state of some items of a row of a logical bidimensional array
      ! the indexes of the elements to change are passed as elements of an integer array

      !f2py intent(hide)  :: n_rows, n_elements, n_indexes
      !f2py intent(in)    :: row, indexes_array, newstate
      !f2py intent(inout) :: array

      ! Parameters
      logical, intent(inout) :: array(1:n_rows, 1:n_elements)   ! some of the True elements of a row of this array must be changed
      integer, intent(in)    :: indexes_array(1:n_indexes)      ! each element identifies one of the True items to set False
      integer, intent(in)    :: n_rows, n_elements, n_indexes
      integer, intent(in)    :: row                             ! row on which to operate
      logical, intent(in)    :: newstate                        ! state to set

      ! Local variables
      integer :: i

      if ((row > 0).and.(row <= n_rows)) then ! if the selected row is out of range do nothing
         do i = 1, n_indexes
            array(row, indexes_array(i)) = newstate
         enddo
      endif
 
   end subroutine set_elements


   subroutine invert_elements(array, n_rows, n_elements, row, indexes_array, n_indexes)

      ! Invert (from .true. to .false. or viceversa) the logical state of some items of a row of a logical bidimensional array
      ! the indexes of the elements to change are passed as elements of an integer array

      !f2py intent(hide)  :: n_rows, n_elements, n_indexes
      !f2py intent(in)    :: row, indexes_array
      !f2py intent(inout) :: array

      ! Parameters
      logical, intent(inout) :: array(1:n_rows, 1:n_elements)   ! some of the True elements of a row of this array must be changed
      integer, intent(in)    :: indexes_array(1:n_indexes)      ! each element identifies one of the True items to set False
      integer, intent(in)    :: n_rows, n_elements, n_indexes
      integer, intent(in)    :: row                             ! row on which to operate

      ! Local variables
      integer :: i
      logical :: newstate

      if ((row > 0).and.(row <= n_rows)) then ! if the selected row is out of range do nothing
         do i = 1, n_indexes                 
            newstate = .not.(array(row, indexes_array(i)))
            array(row, indexes_array(i)) = newstate 
         enddo
      endif
 
   end subroutine invert_elements


   subroutine invert_elements_relative(array, n_rows, n_elements, row, relative_indexes, n_indexes, state_tochange)

      ! Invert (from .true. to .false. or viceversa) the logical state of a given number of .true. (.false.) items 
      ! of a row of a logical bidimensional array
      ! the .true. (.false.) elements to change are identified by the elements of an integer array
      ! each of the elements of the "indexes_array" array identifies a .true. (.false.) item to be changed
      ! Example
      ! if any of the indexes_array elements is 5, then the 5th .true. (.false.) item 
      ! of the selected row of the array (which is NOT usally the 5th element of the array) 
      ! will be set to .false. (.true.)

      !f2py intent(hide)  :: n_rows, n_elements, n_indexes
      !f2py intent(in)    :: row, indexes_array, state_tochange
      !f2py intent(inout) :: array

      ! Parameters
      logical, intent(inout) :: array(1:n_rows, 1:n_elements)   ! some of the True elements of a row of this array must be changed
      integer, intent(in)    :: relative_indexes(1:n_indexes)   ! each element identifies one of the True items to set False
      integer, intent(in)    :: n_rows, n_elements, n_indexes
      integer, intent(in)    :: row                             ! row on which to operate
      logical, intent(in)    :: state_tochange                  ! logical state of the elements to be changed 

      ! Local variables
      integer :: absolute_indexes(1:n_indexes)

      call find_absolute_indexes(array, n_rows, n_elements, row, relative_indexes, absolute_indexes, n_indexes, state_tochange)

      call invert_elements(array, n_rows, n_elements, row, absolute_indexes, n_indexes)

   end subroutine invert_elements_relative


   integer function true_size(logical_array, n_elements)

      ! Returns the number of true elements of a logical array

      !f2py intent(in)   :: logical_array
      !f2py intent(hide) :: n_elements

      ! Parameters
      logical, intent(in) :: logical_array(1:n_elements)
      integer, intent(in) :: n_elements

      ! Local variables
      integer:: i, array_size 
      
      array_size = 0 
      do i=1, n_elements
         if (logical_array(i)) array_size = array_size + 1
      enddo
      true_size = array_size

   end function true_size


end module f_arrays
