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


module f_io

   implicit none

contains

   subroutine pause()

      ! Local variables

      character :: c

      ! Stops the program execution and waits until "RETURN" is pressed

      print *, ""
      print *, "*** Press RETURN to continue ***"
      !   read (*,'()')
      call fget(c)   ! WARNING: this is a GNU95 extension and is not portable

   end subroutine pause

end module f_io
