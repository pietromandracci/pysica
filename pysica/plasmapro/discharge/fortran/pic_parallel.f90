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

module f_pic

   use f_precision

   use f_constants

   use f_io

   implicit none

contains

   subroutine calculate_velocity_increments(dv_z, rho, psi, z, isactive, weight, cm_ratio, &
                                           &dt, psi_0, distance, area, N_cells, N_types, N_particles)

      ! Calculate the velocity increments (z-component) of electrons and ions, based on the values of potential at the electrodes
      ! and on the charge distribution, using the Particle In Cell method.
      ! a) charge density is calculated at grid points;
      ! b) electric potential is calculated at grid points;
      ! c) electric field z-components are calculated in the center of each cell;
      ! d) acceleration and velocity increment in the z direction are calculated for each charge,
      !    considering the contribution of the electric field at each cell
      ! In steps (a) and (d) the "cell size" is considered: a function that weights the charge of each particle inside the cell
    
      !f2py intent(out)  :: dv_z, rho, psi
      !f2py intent(in)   :: z, isactive, weight, cm_ratio, dt, psi_0, distance, area
      !f2py intent(hide) :: N_cells, N_types, N_particles

      ! Parameters
      real(dp), dimension(0:N_types,1:N_particles),  intent(out) :: dv_z         ! Velocity increment during timestep, z-component
      real(dp), dimension(0:N_cells),                intent(out) :: rho          ! Charge density
      real(dp), dimension(0:N_cells),                intent(out) :: psi          ! Electric potential
      real(dp), dimension(0:N_types,1:N_particles),  intent(in)  :: z            ! Particles positions, z-component
      logical,  dimension(0:N_types,1:N_particles),  intent(in)  :: isactive     ! Select which particles are active
      real(dp), dimension(0:N_types),                intent(in)  :: weight       ! Simulation weight of each particle type
      real(dp), dimension(0:N_types),                intent(in)  :: cm_ratio     ! Charge to mass ratio of each particle type (signed)
      real(dp),                                      intent(in)  :: dt           ! Simulation timestep
      real(dp),                                      intent(in)  :: psi_0        ! Actual value of the electric potential
                                                                                 !    at the biased electrode
      real(dp),                                      intent(in)  :: distance     ! Distance between the electrodes
      real(dp),                                      intent(in)  :: area         ! Area of the electrodes
      integer,                                       intent(in)  :: N_cells
      integer,                                       intent(in)  :: N_types
      integer,                                       intent(in)  :: N_particles

      ! Local variables
      real(dp), dimension(1:N_cells) :: E            ! Electric field, z component
      integer                        :: i, j
      real(dp)                       :: delta, psi_N     

      delta = distance / real(N_cells)
      psi_N = 0.0

      ! Calculate the charge density at the grid points
      call calculate_charge_density(rho, z, isactive, weight, delta, area, N_types, N_particles, N_cells)
      
      ! Calculate the electric potential at the grid points
      call calculate_potential(E, psi, rho, N_cells, psi_0, psi_N, delta) 
      
      ! Calculate the velocity increment for each active particle
      do i = 0, N_types
         !$omp parallel do default(shared), private(j)         
         do j = 1, N_particles
            if (isactive(i,j)) then
               dv_z(i,j) = electric_acceleration(z(i,j), cm_ratio(i), delta, E, N_cells) * dt 
            endif
         enddo
         !$omp end parallel do
      enddo

   end subroutine calculate_velocity_increments


   subroutine calculate_charge_density(rho, z, isactive, weight, delta, area, N_types, N_particles, N_cells)

      !f2py intent(out)  :: rho
      !f2py intent(in)   :: z, isactive, weight, delta, area, N_cells
      !f2py intent(hide) :: N_types, N_particles

      ! Parameters
      real(dp), dimension(0:N_cells),                intent(out) :: rho          ! Electric charge density
      real(dp), dimension(0:N_types,1:N_particles),  intent(in)  :: z            ! Particles positions, z-component
      logical,  dimension(0:N_types,1:N_particles),  intent(in)  :: isactive     ! Select which particles are active
      real(dp), dimension(0:N_types),                intent(in)  :: weight       ! Simulation weight of each particle type
      real(dp),                                      intent(in)  :: delta        ! Spatial dimension of each interval
      real(dp),                                      intent(in)  :: area         ! Area of the electrodes
      integer,                                       intent(in)  :: N_types
      integer,                                       intent(in)  :: N_particles
      integer,                                       intent(in)  :: N_cells

      ! Local variables
      integer  :: grid_index, type_index, particle_index
      real(dp) :: z_grid, distance, multiplier

      do grid_index = 0, N_cells
         z_grid = real(grid_index) * delta
         rho(grid_index) = 0.0
         do type_index = 0, N_types              
            if (type_index == 0) then
               multiplier =   weight(type_index)
            else
               multiplier = - weight(type_index)
            endif
!            !$omp parallel do default(shared) private(particle_index, distance) reduction(+:rho) ! PEGGIORATIVO
            do particle_index = 1, N_particles
               if (isactive(type_index, particle_index)) then
                  distance = z_grid - z(type_index, particle_index) 
                  rho(grid_index) = rho(grid_index) +  multiplier * particle_size(distance, delta)
               endif
            enddo ! particle_index = 1, N_particles
!            !$omp end parallel do
         enddo ! type_index = 0, N_types
      enddo ! grid_index = 0, N_cells
      rho = rho * E_CHARGE / (area*delta)

   end subroutine calculate_charge_density


   pure subroutine calculate_potential(E_field, psi, rho, N_cells, psi_0, psi_N, delta)

      !f2py intent(out)  :: psi, E_field
      !f2py intent(in)   :: rho, psi_0, psi_N, delta
      !f2py intent(hide) :: N_cells

      ! Parameters
      real(dp), dimension(1:N_cells), intent(out) :: E_field             ! Electric field
      real(dp), dimension(0:N_cells), intent(out) :: psi                 ! Electric potential
      real(dp), dimension(0:N_cells), intent(in)  :: rho                 ! Electric charge density
      real(dp),                       intent(in)  :: psi_0, psi_N
      real(dp),                       intent(in)  :: delta               ! Spatial dimension of each interval
      integer,                        intent(in)  :: N_cells  

      ! Local variables
      integer                            :: i
      real(dp)                           :: psi_1, C_sum
      real(dp), dimension(0:N_cells) :: constant            ! Electric charge density

      ! Calculate constants C = rho * delta**2 / epsilon0
      do i = 0, N_cells-1
         constant(i) = - rho(i) * delta*delta / EPSILONZERO
      enddo

      ! Calculate summation 
      C_sum = 0
      do i = 1, N_cells-1
         C_sum = C_sum + real(i) * constant(N_cells-i)
      enddo

      ! Calculate potential at i = 1
      psi_1 = psi_N + real(N_cells-1) * psi_0 - C_sum
      psi_1 = psi_1 / N_cells

      ! Calculate all potentials
      psi(0) = psi_0
      psi(1) = psi_1
      do i = 2, N_cells
         psi(i) = 2 * psi(i-1) - psi(i-2) + constant(i-1)
      enddo 

      ! Calculate electric field
      do i = 1, N_cells
         E_field(i) = ( psi(i-1) - psi(i) ) / delta
      enddo

   end subroutine calculate_potential


   real(dp) pure function electric_acceleration(z_charge, cm_ratio, delta, E_field, N_cells)

      ! Calculate the z-component of the acceleration of a charge, accordiing to the PIC scheme
     
      !f2py intent(in)   :: z_charge, cm_ratio, electric_field, delta
      !f2py intent(hide) :: N_cells

      ! Parameters
      real(dp),                       intent(in) :: z_charge      ! Position of the charge on which the force acts
      real(dp),                       intent(in) :: cm_ratio      ! Charge to mass ratio (whith sign)
      real(dp),                       intent(in) :: delta         ! Separation between the points of the grid
                                                                  !    on which E_field is calculated
      real(dp), dimension(1:N_cells), intent(in) :: E_field       ! Electric field intensity at each cell 
      integer,                        intent(in) :: N_cells

      ! Local variables
      integer  :: i
      real(dp) :: z_field, distance, sum

      sum = 0.0
!      !$omp parallel do default(shared) private(i, z_field, distance)  
      do i = 1, N_cells
         z_field = ( real(i) - 0.5 ) * delta ! Position of cell center
         distance = z_field - z_charge
         sum = sum + E_field(i) * particle_size(distance, delta)
      enddo
!      !$omp end parallel do 
      electric_acceleration = cm_ratio * sum

   end function electric_acceleration 


   real(dp) pure function particle_size(x, delta)

      ! Size function for the PIC scheme

      !f2py intent(in)  :: x, delta
         
      ! Parameters
      real(dp), intent(in)  :: x, delta

      ! Local variables
      if (delta == 0) then 
         particle_size = 0.0
      else
         if ( abs(x) > abs(delta) ) then
            particle_size = 0.0
         else
            particle_size = 1.0 - abs(x/delta)
         endif
      endif

   end function particle_size

end module f_pic
