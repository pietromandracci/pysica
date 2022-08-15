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

module f_constants

   use f_precision

   implicit none
   
   real(dp), parameter :: E_EULER              = 2.718281828459045      ! exp(1.0_dp)
   real(dp), parameter :: PI                   = 3.1415926535897931_dp
   real(dp), parameter :: SQRT_PI              = 1.772453850905515      ! sqrt(PI)
   real(dp), parameter :: K_BOLTZMANN          = 1.3806503E-23_dp       ! Boltzmann's constant / J K**-1
   real(dp), parameter :: ATOMIC_ENERGY_EV     = 27.21_dp               ! Atomic unit of energy / eV
   real(dp), parameter :: BOHR_RADIUS          = 0.29E-10_dp            ! Bohr radius / m
   real(dp), parameter :: E_CHARGE             = -1.602176E-19_dp       ! Electron charge / C
   real(dp), parameter :: E_MASS               = 9.109382E-31_dp        ! Electron mass / kg
   real(dp), parameter :: E_CHARGETOMASSRATIO  = -1.7588196E11_dp       ! Electron charge to mass ratio / C*kg**-1
   real(dp), parameter :: E_V2TOENERGY_AU      = 1.0447687E-13_dp       ! Multiply this by the squared velocity to obtain
                                                                        !   electron kinetic energy in atomic units
   real(dp), parameter :: EPSILONZERO          = 8.854187817E-12_dp     ! Electric permittivity of vacuum


end module f_constants
