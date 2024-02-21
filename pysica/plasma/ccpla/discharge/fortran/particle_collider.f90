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

module f_particle_collider

   use f_precision

   use f_io

   use f_random

   use f_arrays

   use f_constants

   use f_scattering

   implicit none

contains

   subroutine check_collisions(x, y, z, vx, vy, vz, v, v2, &
                              &isactive, restart, &
                              &cm_ratio, weight,  &        
                              &e_loss, e_loss_ions, &                             
                              &p_collision, p_collision_ions, &
                              &p_limits, n_limits, &
                              &limit1_exc, limit1_diss, &
                              &coll_f_tot_ions, p_limits_ions, &
                              &kdt_recomb, ion_density, isactive_recomb, &
                              &min_scattered, &
                              &n_totpro, &
                              &en_ion, en_exc, en_diss, n_excpro, n_disspro, mean_speed, &    
                              &coll_null, coll_ela, coll_ion, coll_exc, coll_dis, coll_rec, &              
                              &v_values, v_values_ions, n_v_values, n_v_values_ions, &
                              & x_part_added,  y_part_added,  z_part_added, &
                              &vx_part_added, vy_part_added, vz_part_added, &
                              &w_part_added, n_part_added, &
                              &n_types, n_particles, nmax_excpro, nmax_disspro, &
                              &debug_level)

      ! Parameters
      real(dp), dimension(0:n_types, 1:n_particles), intent(inout) :: x, y, z,&            ! Positions and velocities of electrons and ions
                                                                     &vx, vy, vz, v, v2
      logical,  dimension(0:n_types, 1:n_particles), intent(inout) :: isactive             ! Select which elctrons/ions are active
      logical,  dimension(0:n_types, 1:n_particles), intent(inout) :: restart              ! Select which particles must restart
                                                                                           !    leap-frog scheme at the next iteration 
      real(dp), dimension(0:n_types),                intent(in)    :: cm_ratio             ! Charge/mass ratio (signed) of each
                                                                                           !    particle type
      real(dp), dimension(0:n_types),                intent(inout) :: weight               ! Computational weight of each particle type 
                                                                                           !    (number of real particles
                                                                                           !    represented by each simulated one)
      real(dp), dimension(1:n_types),                intent(in)    :: e_loss               ! Mean loss of e- energy in elastic scattering
      real(dp), dimension(1:n_types,1:n_types),      intent(in)    :: e_loss_ions          ! Mean loss of ion energy in elastic scattering      
      real(dp), dimension(1:n_v_values),             intent(in)    :: v_values             ! Velocities for which electrons
                                                                                           !    collision frequency values are given
      real(dp), dimension(1:n_v_values, &                                                  ! Velocities for which ions collision frequency 
                         &1:n_v_values_ions),        intent(in)    :: v_values_ions        !    values are given
      integer,                                       intent(in)    :: n_v_values           ! Number of tabulated values of el collision freq.
      integer,                                       intent(in)    :: n_v_values_ions      ! Number of tabulated values of ion collision freq.
      real(dp),                                      intent(in)    :: p_collision          ! Collision probability (including null collisions)
                                                                                           !    for electrons
      real(dp), dimension(1:n_types),                intent(in)    :: p_collision_ions     ! Collision probability (including null collisions)
                                                                                           !    for ions
      real(dp), dimension(1:n_totpro,1:n_v_values),  intent(in)    :: p_limits             ! Values used to select collision type for electrons
      integer,                                       intent(in)    :: n_limits             ! Number of probability limits for el collisions
      integer,                                       intent(in)    :: limit1_exc           ! Index of prob limit for e- impact excitation
      integer,                                       intent(in)    :: limit1_diss          ! Index of prob limit for e- impact dissociation
      real(dp), dimension(1:n_types, &
                         &1:n_types, &
                         &1:n_v_values_ions),        intent(in)    :: coll_f_tot_ions      ! Tolal collision frequency (all processes)
                                                                                           !    of an ion on a specific neutral   
      real(dp), dimension(1:n_types, &
                         &1:n_types, &
                         &1:2, &
                         &1:n_v_values_ions),        intent(in)    :: p_limits_ions        ! Probability limits for elastic, charge exchange
                                                                                           !   and null  scattering for a
                                                                                           !   specific ion/neutral combination
      real(dp), dimension(1:n_types,1:nmax_disspro, &
                         &1:n_v_values),             intent(in)    :: kdt_recomb           ! Collision rate coefficients for e-/ion recombination
      integer,                                       intent(in)    :: min_scattered        ! Minimum number of scattered particles to apply
                                                                                           !    the many particles method 
                                                                                           !    multiplied by timestep / m**3
      real(dp), dimension(1:n_types),                intent(in)    :: ion_density          ! Actual number density of each ion type
      logical,                                       intent(in)    :: isactive_recomb      ! Set if ion-electron recombination is active
      integer,                                       intent(in)    :: n_totpro             ! Total maximum possible number of scattering
                                                                                           !   processes for electrons:
                                                                                           !    elastic + ionization + dissociation
      real(dp),                                      intent(inout) :: coll_null            ! Number of observed null collisions
      real(dp), dimension(1:n_types),                intent(inout) :: coll_ela             ! Number of observed elastic collisions
                                                                                           !    for each neutral gas type
      real(dp), dimension(1:n_types),                intent(inout) :: coll_ion             ! Number of observed ionization collisions for each 
                                                                                           !    neutral gas type
      real(dp), dimension(1:n_types,1:nmax_excpro),  intent(inout) :: coll_exc             ! Number of observed excitation collisions for each
                                                                                           !    neutral gas type and excitation process     
      real(dp), dimension(1:n_types,1:nmax_disspro), intent(inout) :: coll_dis             ! Number of observed dissociation collisions for each
                                                                                           !    neutral gas type and dissociation process 
      real(dp), dimension(1:n_types,1:nmax_disspro), intent(inout) :: coll_rec             ! Number of observed electron-ion recombinations
                                                                                           !    for each neutral gas type
                                                                                           !    and dissociative recombination process 
      real(dp), dimension(1:n_types),                intent(in)    :: en_ion               ! Ionization energy of each neutral type
      real(dp), dimension(1:n_types,1:nmax_excpro),  intent(in)    :: en_exc               ! Excitation energy of each neutral type and 
                                                                                           !    each excitation process      
      real(dp), dimension(1:n_types,1:nmax_disspro), intent(in)    :: en_diss              ! Dissociation energy of each neutral type and 
                                                                                           !    each dissociation process
      integer,  dimension(1:n_types),                intent(in)    :: n_excpro             ! Number of excitation processes available 
                                                                                           !    for each neutral type      
      integer,  dimension(1:n_types),                intent(in)    :: n_disspro            ! Number of dissociation processes available 
                                                                                           !    for each neutral type
      integer,                                       intent(in)    :: nmax_excpro          ! Maximum number of excitation types      
      integer,                                       intent(in)    :: nmax_disspro         ! Maximum number of dissociation types
      real(dp), dimension(1:n_types),                intent(in)    :: mean_speed           ! Mean speed of neutrals molecules/ m*s**-1

      real(dp), dimension(0:n_types, 1:n_particles), intent(inout) :: x_part_added, &
                                                                     &y_part_added, &
                                                                     &z_part_added         ! Positions of particles that need to be added
      real(dp), dimension(0:n_types, 1:n_particles), intent(inout) :: vx_part_added, &
                                                                     &vy_part_added, &
                                                                     &vz_part_added        ! Velocities of particles that need to be added
      integer,  dimension(0:n_types),                intent(inout) :: n_part_added         ! Number of particles that need to be added
      integer,  dimension(0:n_types, 1:n_particles), intent(inout) :: w_part_added         ! Number of computational particles to add
                                                                                           !   for each real one
      integer,                                       intent(in)    :: n_types, n_particles ! Number of ion types and max number of
                                                                                           !   active electron/ions (for each type)
      integer,                                       intent(in)    :: debug_level


      ! Local variables
      integer  :: type_index, n_active, n_scattered 
      real(dp) :: collision_probability

      do type_index = 0, n_types
         ! Calculate the number of active particles of this type
         n_active =  count( isactive(type_index,:) ) 

         ! Select the correct collision probability (including null collisions), depending on the particle type (electron or ion)
         if (type_index == 0) then
            collision_probability  = p_collision                     ! electrons
         else
            collision_probability =  p_collision_ions(type_index)    ! ions
         endif

         ! Calculate the mean (expected) number of scattering particles (including null collisions)
         !    -> collision probability by number of active particles
         n_scattered = nint(collision_probability * n_active)

!         if (type_index == 0) print *,'***** type ', type_index, ' --> ', &
!                                    & nint(collision_probability*100), n_active, n_scattered !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! If the number of scattered particles is high enough, apply the faster method, otherwise the slower one
         if ( (n_scattered > min_scattered) ) then  
            call manage_collisions_fast_method(type_index, n_active, n_scattered)
         else 
            call manage_collisions_slow_method(type_index, collision_probability)          
         endif
      enddo

   contains

     
      ! +-----------------------------+
      ! | General collision managment |
      ! +-----------------------------+

      subroutine manage_collisions_slow_method(type_index, collision_probability)

         ! For each particle, compare the collision probability (including null collisions) with a random number
         ! if there is a collision, call the proper subroutine that will decide which type of collision is (maybe null)
         ! this method requires calculating a random number for each particle (scattered or not), so it is slower the the other one

         ! Parameters
         integer,  intent(in) :: type_index             ! type of scattering particles
         real(dp), intent(in) :: collision_probability  ! probability that a single particle has a collision (including null ones)
                                                        ! during the timestep dt

         ! Local variables
         integer  :: particle_index
         logical  :: collision
         real(sp) :: rand

         do particle_index = 1, n_particles                                   
            ! Operate only on existing particles
            if (isactive(type_index,particle_index)) then                      
               call random_number(rand)
               ! If there is a collision (maybe also a null one)
               if (rand < collision_probability) then
                  if ( (debug_level > 1) ) print *, 'COLLISION DETECTED (SLOW METHOD)' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! If the particle is an electron (type_index=0)
                  if (type_index == 0) then 
                     collision = .false.
                     call electron_collisions(particle_index, collision)
                     ! If it was a null collision, check if a recombination happened (if recombination processes have been set to on)
                     if (collision .eqv. .false.) then
                        coll_null = coll_null + weight(0)
                        if (isactive_recomb) call electron_ion_recombination(particle_index)
                     endif
                  ! If the particle is an ion (type_index>0)  
                  else ! (type_index == 0)
                     call ion_collisions(type_index,particle_index)
                  endif ! (type_index == 0)
               ! If there was no collision    
               else ! (rand < collision_probability)
                  ! If the particle is an electron, and there was no collision, check if a recombination happened
                  ! (if recombination processes have been set to on)
                  if ( (type_index == 0) .and. (isactive_recomb) ) call electron_ion_recombination(particle_index)
               endif ! (rand < collision_probability)
            endif ! (isactive(particle_index))
         enddo ! particle_index

      end subroutine manage_collisions_slow_method


      subroutine manage_collisions_fast_method(type_index, n_active, n_scattered)

         ! Generate a set of random indexes, each representing a scattered particle
         ! the number of indexes is equal to the mean number of scattered particles
         ! this method requires calculating a random number for each *scattered* particle only (including null scattering events)

         ! Parameters
         integer, intent(in) :: type_index   ! type of scattering particles
         integer, intent(in) :: n_active     ! number of active particles of the type considered
         integer, intent(in) :: n_scattered  ! mean number of scattered particles, equal to the number of random indexes to generate

         ! Local variables
         integer                            :: i, particle_index
         logical                            :: collision
         logical, dimension(1:n_particles)  :: check_electron_recombination
         integer, allocatable               :: scattered_indexes_relative(:), scattered_indexes_absolute(:)

         ! Allocate arrays to deal with the indexes of scattered particles
         allocate(scattered_indexes_relative(1:n_scattered))
         allocate(scattered_indexes_absolute(1:n_scattered))
      
         ! Determine randomly the indexes of particles that were subject to scattering
         ! will return a set of relative indexes, identifying the scattered particles
         ! each relative index identifies the position of the particle among the active ones
         !   e.g. 1 is the first active particle, i.e. the first position for which the isactive array has a .true. value
         call random_integers_array(scattered_indexes_relative, n_scattered, 1, n_active, .true.)
         
         ! Transform relative indexes to absolute indexes
         ! each absolute index identify the position of the particle in the position and velocity arrays (x,y,z,vx,vy,vz,v, etc...)
         ! e.g. if the first three positions are occupied by inactive particles, the relative index 1 corresponds to the absolute index 4
         call find_absolute_indexes(isactive, n_types+1, n_particles, type_index+1, &
                                   &scattered_indexes_relative, scattered_indexes_absolute, n_scattered, .true.)   

         ! Manage electron collisions
         if (type_index == 0) then
            ! Check electron collisions, including electron/ion recombination
            if (isactive_recomb) then
               ! Initialize the array used to keep track of scattered electrons, since only unscattered ones may recombine
               check_electron_recombination = isactive(0, :) ! Only active electrons can recombine
               do i = 1, n_scattered
                  if ( (debug_level > 1) ) print *, 'EL COLLISION DETECTED (FAST METHOD) WITH RECOMBINATION' !!!!!!!!!!!!!!!!!!!!!
                  particle_index = scattered_indexes_absolute(i)
                  collision = .false.
                  call electron_collisions(particle_index, collision)
                  ! If the electron had a non-null collision, will not be checked for recombination
                  if (collision) then 
                     check_electron_recombination(particle_index) = .false.  ! Since the electron has scattered, it can not recombine
                  else
                     coll_null = coll_null + weight(0)
                  endif
               enddo ! i
               ! Manage electron-ion recombinations (of unscattered electrons)
               do particle_index = 1, n_particles
                  if (check_electron_recombination(particle_index)) call electron_ion_recombination(particle_index)
               enddo          
            ! Check electron collisions, excluding electron/ion recombination
            else
               do i = 1, n_scattered
                  if ( (debug_level > 1) ) print *, 'EL COLLISION DETECTED (FAST METHOD) WITHOUT RECOMBINATION' !!!!!!!!!!!!!!!!!!
                  particle_index = scattered_indexes_absolute(i)
                  collision = .false.
                  call electron_collisions(particle_index, collision)
                  if (.not.(collision)) coll_null = coll_null + weight(0)
               enddo ! i 
            endif
         ! Manage ions collisions
         else
            do i = 1, n_scattered
               if ( (debug_level > 1) ) print *, 'ION COLLISION DETECTED (FAST METHOD)' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               particle_index = scattered_indexes_absolute(i)
               call ion_collisions(type_index,particle_index)
            enddo ! particle_index 
         endif !(type_index == 0)

         deallocate(scattered_indexes_relative, scattered_indexes_absolute)

      end subroutine manage_collisions_fast_method


      ! +----------------------------------------+
      ! | Electron-neutral collisions management |
      ! +----------------------------------------+

      subroutine electron_collisions(electron_index, collision_happened)

         ! Once a collision between an electron and a neutral has been detected,
         ! select which is the collision type: elastic, ionization, excitation, dissociation
         ! or a null collision (i.e. no collision at all) and call the proper subroutine to manage it
 
         ! Parameters
         integer, intent(in)    :: electron_index         ! Index identifying the scattered eletron
         logical, intent(inout) :: collision_happened     ! Will be set to true if the collision is not null

         ! Local variables
         character*12 :: text
         real(dp)     :: v2_0
         integer      :: l, v_index, process_index, neutral_type, exc_type, diss_type
         real(sp)     :: rand2 

         ! Initialize some variables for debug purposes
         text = '------------'
         neutral_type = -1
         exc_type     = -1
         diss_type    = -1
         v2_0         = v2(0, electron_index)

         ! Select the nearest velocity value for which collision prob. are tabulated
         v_index = search_index(v_values,  n_v_values, v(0,electron_index), .true.)

         ! Decide which type of scattering happened
         call random_number(rand2)
         do process_index = 1, n_totpro
            if ( rand2 <= p_limits(process_index, v_index) ) exit
         enddo         
         ! plimits from 1 to n_types are for elastic scattering
         if (process_index <= n_types) then 
            neutral_type = process_index
            call elastic_scattering_electron(electron_index, neutral_type)        
            collision_happened = .true.
            text = 'ELASTIC     '            
         ! plimits from n_types+1 to limit1_exc are for ionization scattering
         ! note that limit1_exc is the first excitation limit in the python array, which starts from 0
         ! in the Fortran array, starting from 1, it is the last ionization index
         elseif (process_index <= limit1_exc) then
            neutral_type = process_index - n_types
            call ionization_scattering_electron(electron_index, neutral_type)
            collision_happened = .true.
            text = 'IONIZATION  '            
         ! plimits from limit1_exc+1 to limit_diss are for excitation scattering
         ! note that limit1_diss is the first dissociation limit in the python array, which starts from 0
         ! in the Fortran array, starting from 1, it is the last excitation index            
         !elseif (process_index <= n_limits) then
         elseif (process_index <= limit1_diss) then
            ! Identify which type of neutral and which excitation process are involved
            l = process_index - limit1_exc
            do neutral_type = 1, n_types
               if ( l > n_excpro(neutral_type) ) then
                  l = l - n_excpro(neutral_type)
               else 
                  exc_type = l
                  exit
               endif 
            enddo ! neutral_type
            call excitation_scattering_electron(electron_index, neutral_type, exc_type)
            collision_happened = .true.
            text = 'EXCITATION  '            
         ! plimits from limit_diss+1 to n_limits are for dissociation scattering    
         elseif (process_index <= n_limits) then
            ! Identify which type of neutral and which dissociation process are involved
            l = process_index - limit1_diss
            do neutral_type = 1, n_types
               if ( l > n_disspro(neutral_type) ) then
                  l = l - n_disspro(neutral_type)
               else 
                  diss_type = l
                  exit
               endif 
            enddo ! neutral_type
            call dissociation_scattering_electron(electron_index, neutral_type, diss_type)
            collision_happened = .true.
            text = 'DISSOCIATION'
         ! plimits greater than n_limits are for null collisions   
         else            ! Null collision: do nothing  
            collision_happened = .false.
            text = 'NULL       '
         endif  ! (process_index <= n_types) 
   
         if ( (debug_level > 1).and.collision_happened ) then
            call print_collision_report_electron(text, rand2, electron_index, v2_0, v_index, &
                                               &neutral_type, process_index, diss_type)
         endif 

      end subroutine electron_collisions


      ! +-------------------------------------+
      ! | Electron-neutral elastic collisions |
      ! +-------------------------------------+

      subroutine elastic_scattering_electron(electron_index, neutral_type)

         ! Parameters
         integer, intent(in) :: electron_index, neutral_type

         ! Local variables
         real(dp) :: verx_scatter      ! Projection of scattered velocity versor on x-axis in the frame of reference of the target
         real(dp) :: very_scatter      ! Projection of scattered velocity versor on y-axis in the frame of reference of the target
         real(dp) :: verz_scatter      ! Projection of scattered velocity versor on z-axis in the frame of reference of the target
         real(dp) :: verx1,very1,verz1 ! Projections on x,y,z axis of the velocity versor after scattering
         logical  :: debug

         if (debug_level > 2) then
            debug = .true.
         else
            debug = .false.
         endif
 
         ! Calculate scattering angles (respect to particle velocity before scattering)
         call anisotropic_scattering_electron(verx_scatter, very_scatter, verz_scatter, v2(0, electron_index), .true., debug)
         ! call isotropic_scattering(verx_scatter, very_scatter, verz_scatter, debug)

         ! Calculate the direction of particle velocity after scattering
         call scattering_tarfor2disfor(verx1,very1,verz1, &
                                      &vx(0,electron_index), &
                                      &vy(0,electron_index), &
                                      &vz(0,electron_index), &
                                      &v( 0,electron_index),  &
                                      &verx_scatter, very_scatter, verz_scatter, &
                                      &debug)

         ! Calculate new electron velocity module (due to kinetic energy transferred to the neutral particle)
         ! NOTE: verz_scatter is the cosinus of the scattering angle in the frame of reference of the target
         v2(0,electron_index) = v2(0,electron_index) * (1 - e_loss(neutral_type) * (1 - verz_scatter))
         v(0,electron_index)  = sqrt(v2(0,electron_index))
         
         ! Calculate the new components of velocity vector
         vx(0,electron_index) = v(0,electron_index) * verx1
         vy(0,electron_index) = v(0,electron_index) * very1
         vz(0,electron_index) = v(0,electron_index) * verz1

         ! Increase collisions counter
         coll_ela(neutral_type) = coll_ela(neutral_type) + weight(0)

      end subroutine elastic_scattering_electron


      ! +----------------------------------------+
      ! | Electron-neutral ionization collisions |
      ! +----------------------------------------+      
      
      subroutine ionization_scattering_electron(electron_index, neutral_type)

         ! Parameters
         integer, intent(in) :: electron_index, neutral_type

         ! Local variables
         real(dp) :: res_energy_j      ! residual energy: energy of incoming electron less ionization energy (in joules)
         real(dp) :: e_sc_energy_j     ! scattered electron energy (in joules)
         real(dp) :: e_ej_energy_j     ! ejected electron energy (in joules)
         real(dp) :: v2_ej             ! ejected electron squared speed
         real(dp) :: v_ej              ! ejected electron speed
         real(dp) :: verx_scatter      ! Projection of scattered velocity versor on x-axis in the electron velocity frame (*)
         real(dp) :: very_scatter      ! Projection of scattered velocity versor on y-axis in the electron velocity frame (*)
         real(dp) :: verz_scatter      ! Projection of scattered velocity versor on z-axis in the electron velocity frame (*)
         real(dp) :: verx1,very1,verz1 ! Projections on x,y,z axis of the velocity versor after scattering
         real(dp) :: vx_neu            ! Velocity x-component of the (randomly generated) ionized atom 
         real(dp) :: vy_neu            ! Velocity y-component of the (randomly generated) ionized atom 
         real(dp) :: vz_neu            ! Velocity z-component of the (randomly generated) ionized atom 
         real(dp) :: v_neu             ! Velocity x-component of the (randomly generated) ionized atom 
         logical  :: debug

         ! (*) i.e. the frame in which the z-axis has the same direction as the impinging electron velocity vector before scattering

         if (debug_level > 2) then
            debug = .true.
         else
            debug = .false.
         endif

         ! Calculate the residual energy after ionization (incoming electron energy less ionization energy) in joules
         ! note that E_CHARGE is negative
         res_energy_j = 0.5 * E_MASS * v2(0, electron_index) + en_ion(neutral_type) * E_CHARGE
         
         ! If the residual energy is negative, the ionization process is not possible
         if (res_energy_j < 0.0) then
            if (debug) print *, 'ERROR: negative residual energy found in ionization process! Residual energy =', res_energy_j
            return
         endif
         
         ! Calculate new direction of electron velocity after scattering

         ! Generate a new random direction for scattered electron (respect to impinging electron direction)
         call anisotropic_scattering_electron(verx_scatter, very_scatter, verz_scatter, v2(0, electron_index), .true., debug)

         ! Calculate the direction of scattered electron in the laboratory frame
         call scattering_tarfor2disfor(verx1, very1, verz1, &
                                      &vx(0,electron_index), &
                                      &vy(0,electron_index), &
                                      &vz(0,electron_index), &
                                      &v(0,electron_index),  &
                                      &verx_scatter, very_scatter, verz_scatter, &
                                      &debug)

         ! Calculate new modulus of electron velocity after scattering

         ! Residual energy is equally shared between the scattered and ejected electron
         e_sc_energy_j = res_energy_j / 2.0         ! energy of the scattered electron            
         e_ej_energy_j = e_sc_energy_j              ! energy of the electron ejected from the atom
         
         ! Calculate scattered electron speed
         v2(0,electron_index) = 2.0 * e_sc_energy_j / E_MASS
         v( 0,electron_index) = sqrt(v2(0,electron_index))

         ! Calculate the new components of the scattered electron velocity vector after scattering
         ! from modulus and direction
         vx(0,electron_index) = v(0,electron_index) * verx1
         vy(0,electron_index) = v(0,electron_index) * very1
         vz(0,electron_index) = v(0,electron_index) * verz1

         ! Calculate the ejected electron speed
         v2_ej = 2.0 * e_ej_energy_j / E_MASS
         v_ej  = sqrt(v2_ej)

         ! The ejected electron velocity versor is kept simmetric respect to the scattered electron one
         ! the z component of the versor is the same (z is the direction of the impinging electron velocity)
         ! the x and y component are opposite
         verx_scatter = - verx_scatter
         very_scatter = - very_scatter

         ! Calculate the direction of ejected electron in the laboratory frame         
         call scattering_tarfor2disfor(verx1, very1, verz1, &
                                      &vx(0,electron_index), &
                                      &vy(0,electron_index), &
                                      &vz(0,electron_index), &
                                      &v(0,electron_index),  &
                                      &verx_scatter, very_scatter, verz_scatter, &
                                      &debug)

         ! Generate a random velocity vector for the neutral, which will become the ion velocity
         call random_maxwell_velocity(vx_neu, vy_neu, vz_neu, v_neu, mean_speed(neutral_type))
         
         ! Add the ejected electron  and the ion  considering their weights
         call ionization_add_particles(electron_index, neutral_type, &
                                      &v_ej*verx1, v_ej*very1, v_ej*verz1, &
                                      &vx_neu, vy_neu, vz_neu)

         ! Increase collisions counter
         coll_ion(neutral_type) = coll_ion(neutral_type) + weight(0)
         
      end subroutine ionization_scattering_electron

      
      subroutine ionization_add_particles(electron_index, neutral_type, &
                                         &vx_new_el, vy_new_el, vz_new_el,&
                                         &vx_new_io, vy_new_io, vz_new_io)

         ! Parameters
         integer,  intent(in) :: electron_index                   ! Index of the ionizing electron
         integer,  intent(in) :: neutral_type                     ! Type of ionized neutral
         real(dp), intent(in) :: vx_new_el, vy_new_el, vz_new_el  ! Velocity components of the ejected electron
         real(dp), intent(in) :: vx_new_io, vy_new_io, vz_new_io  ! Velocity components of the new ion         

         ! Local variables
         integer  :: n_ions_added
         real(dp) :: weight_ratio
         real(sp) :: rand

         ! Store a computational electron to be added
         n_part_added(0)                  = n_part_added(0) + 1   ! increase the number of events that added electrons (ionization in this case)
         w_part_added(0,n_part_added(0))  = 1                     ! will add a single computational electron for each ionization event
         x_part_added(0,n_part_added(0))  = x(0,electron_index)   ! store the position of the ejected e- (x-component)
         y_part_added(0,n_part_added(0))  = y(0,electron_index)   ! store the position of the ejected e- (y-component)
         z_part_added(0,n_part_added(0))  = z(0,electron_index)   ! store the position of the ejected e- (z-component)
         vx_part_added(0,n_part_added(0)) = vx_new_el             ! store the velocity of the ejected e- (x-component)
         vy_part_added(0,n_part_added(0)) = vy_new_el             ! store the velocity of the ejected e- (y-component)
         vz_part_added(0,n_part_added(0)) = vz_new_el             ! store the velocity of the ejected e- (z-component)         

         ! Calculate the number of computational ions to add as a result of the ionization
         ! It will be equal to the ratio between the weights of electron and ion type
         ! if it is lower then one, then it becomes the probability that an ion is added
         weight_ratio = weight(0) / weight(neutral_type)
         if (weight_ratio >= 1) then
            n_ions_added = nint(weight_ratio)  
         else
            call random_number(rand)
            if (rand < weight_ratio) then
               n_ions_added = 1
            else
               n_ions_added = 0
            endif
         endif
         if (n_ions_added > 0) then
            ! Store a computational ion to be added
            n_part_added(neutral_type)                             = n_part_added(neutral_type) + 1  ! increase the number of events 
                                                                                                     !   that added ions
            w_part_added( neutral_type,n_part_added(neutral_type)) = n_ions_added                    ! number of computational ions
                                                                                                     !   added by this event
            x_part_added( neutral_type,n_part_added(neutral_type)) = x(0,electron_index)             ! ion position is the one of the
            y_part_added( neutral_type,n_part_added(neutral_type)) = y(0,electron_index)             !    impinging electron 
            z_part_added( neutral_type,n_part_added(neutral_type)) = z(0,electron_index)             ! 
            vx_part_added(neutral_type,n_part_added(neutral_type)) = vx_new_io                       ! ion velocity is the one of the 
            vy_part_added(neutral_type,n_part_added(neutral_type)) = vy_new_io                       !    ionized neutral
            vz_part_added(neutral_type,n_part_added(neutral_type)) = vz_new_io                       ! 
         endif

      end subroutine ionization_add_particles


      ! +----------------------------------------+
      ! | Electron-neutral excitation collisions |
      ! +----------------------------------------+      

      subroutine excitation_scattering_electron(electron_index, neutral_type, exc_type)

         ! Parameters
         integer, intent(in) :: electron_index, neutral_type, exc_type

         ! Local variables
         real(dp) :: e_energy_j        ! electron energy in joules
         real(dp) :: verx_scatter      ! Projection of scattered velocity versor on x-axis in the frame of reference of the target
         real(dp) :: very_scatter      ! Projection of scattered velocity versor on y-axis in the frame of reference of the target
         real(dp) :: verz_scatter      ! Projection of scattered velocity versor on z-axis in the frame of reference of the target
         real(dp) :: verx1,very1,verz1 ! Projections on x,y,z axis of the velocity versor after scattering
         logical  :: debug

         if (debug_level > 2) then
            debug = .true.
         else
            debug = .false.
         endif

         ! Calculate particle energy after scattering in joules; E_CHARGE is negative
         e_energy_j = 0.5 * E_MASS * v2(0, electron_index) + en_exc(neutral_type, exc_type) * E_CHARGE

         ! If the electron energy after collision is negative, the excitation process was not possible
         if (e_energy_j < 0.0) then
            if (debug) print *, 'ERROR: negative residual energy found in excitation process! Residual energy =', e_energy_j
            return
         endif         
         
         ! Generate a new random direction for particle velocity after scattering
         call anisotropic_scattering_electron(verx_scatter, very_scatter, verz_scatter, v2(0, electron_index), .true., debug)

         ! Calculate the direction of particle velocity after scattering
         call scattering_tarfor2disfor(verx1,very1,verz1, &
                                      &vx(0,electron_index), &
                                      &vy(0,electron_index), &
                                      &vz(0,electron_index), &
                                      &v(0,electron_index), &
                                      &verx_scatter, very_scatter, verz_scatter, &
                                      &debug)

         ! Calculate particle velocity after scattering
         v2(0,electron_index) = 2.0 * e_energy_j / E_MASS
         v( 0,electron_index) = sqrt(v2(0,electron_index))

         ! Calculate the new components of velocity vector
         vx(0,electron_index) = v(0,electron_index) * verx1
         vy(0,electron_index) = v(0,electron_index) * very1
         vz(0,electron_index) = v(0,electron_index) * verz1

         ! Increase collisions counter
         coll_exc(neutral_type, exc_type) = coll_exc(neutral_type, exc_type) + weight(0)

      end subroutine excitation_scattering_electron

      
      ! +------------------------------------------+
      ! | Electron-neutral dissociation collisions |
      ! +------------------------------------------+      
      
      subroutine dissociation_scattering_electron(electron_index, neutral_type, diss_type)

         ! Parameters
         integer, intent(in) :: electron_index, neutral_type, diss_type

         ! Local variables
         real(dp) :: e_energy_j        ! electron energy in joules
         real(dp) :: verx_scatter      ! Projection of scattered velocity versor on x-axis in the frame of reference of the target
         real(dp) :: very_scatter      ! Projection of scattered velocity versor on y-axis in the frame of reference of the target
         real(dp) :: verz_scatter      ! Projection of scattered velocity versor on z-axis in the frame of reference of the target
         real(dp) :: verx1,very1,verz1 ! Projections on x,y,z axis of the velocity versor after scattering
         logical  :: debug
         
         if (debug_level > 2) then
            debug = .true.
         else
            debug = .false.
         endif

         ! Calculate particle energy after scattering in joules; E_CHARGE is negative
         e_energy_j = 0.5 * E_MASS * v2(0, electron_index) + en_diss(neutral_type, diss_type) * E_CHARGE
         
         ! If the electron energy after collision is negative, the dissocitaion process was not possible
         if (e_energy_j < 0) then
            if(debug) print *, 'ERROR: negative residual energy found in dissociation process! Residual energy =', e_energy_j
            return
         endif         
         
         ! Generate a new random direction for particle velocity after scattering
         call anisotropic_scattering_electron(verx_scatter, very_scatter, verz_scatter, v2(0, electron_index), .true., debug)

         ! Calculate the direction of particle velocity after scattering
         call scattering_tarfor2disfor(verx1,very1,verz1, &
                                      &vx(0,electron_index), &
                                      &vy(0,electron_index), &
                                      &vz(0,electron_index), &
                                      &v(0,electron_index), &
                                      &verx_scatter, very_scatter, verz_scatter, &
                                      &debug)

         ! Calculate particle velocity after scattering
         v2(0,electron_index) = 2.0 * e_energy_j / E_MASS
         v( 0,electron_index) = sqrt(v2(0,electron_index))

         ! Calculate the new components of velocity vector
         vx(0,electron_index) = v(0,electron_index) * verx1
         vy(0,electron_index) = v(0,electron_index) * very1
         vz(0,electron_index) = v(0,electron_index) * verz1

         ! Increase collisions counter
         coll_dis(neutral_type, diss_type) = coll_dis(neutral_type, diss_type) + weight(0)

      end subroutine dissociation_scattering_electron


      ! +-----------------------------------+
      ! | Ion-neutral collisions management |
      ! +-----------------------------------+
      
      subroutine ion_collisions(ion_type, ion_index)

         ! Once a collision between an ion and a neutral has been detected,
         ! select which is the collision type: elastic, charge exchange or a null collision
         ! and call the proper subroutine to manage it        

         ! Parameters
         integer, intent(in) :: ion_type, ion_index

         ! Local variables
         character*12                     :: text
         real(dp), dimension(1:n_types)   :: v_neu, vx_neu, vy_neu, vz_neu     ! Modulus and components of the neutral velocity
         real(dp), dimension(1:n_types)   :: v_rel, vx_rel, vy_rel, vz_rel     ! Modulus and components of the relative velocity ion-neutral
         real(dp)                         :: v_0,   vx_0,   vy_0,   vz_0, v2_0 ! Velocity before collision (stored for debbuging)
         real(dp)                         :: coll_f_sum                        ! Sum of the total collision frequencies for each neutral
         real(dp)                         :: coll_f_ratio                      ! Collision freq for a neutral divided by the total
         real(dp), dimension(1:n_types)   :: p_limits_neutral
         integer,  dimension(1:n_types)   :: v_index
         integer                          :: process_index, neutral_type   
         real(sp)                         :: rand2, rand3
         real(dp)                         :: m_ratio                           ! Ratio of masses: m_ion / m_neutral

         ! Initialize some variables for debug purposes
         text = '------------'
         neutral_type = -1
         vx_0         = vx(ion_type, ion_index)
         vy_0         = vy(ion_type, ion_index)
         vz_0         = vz(ion_type, ion_index)
         v_0          = v(ion_type, ion_index)
         v2_0         = v2(ion_type, ion_index)         

         ! Calculate the relative speed of the ion respect to each neutral type
         do neutral_type = 1, n_types
            ! Generate a random velocity vector for the neutral using a Maxwell distribution
            ! with the specific mean speed of this neutral (mean energy is the same for all neutrals,
            ! since it depends on temperatute only, but mean speed depends on mass also)
            call random_maxwell_velocity(vx_neu(neutral_type), vy_neu(neutral_type), vz_neu(neutral_type), v_neu(neutral_type), &
                                        &mean_speed(neutral_type))            
            ! Calculate the relative speed between the ion and the (radomly generated) neutral
            vx_rel(neutral_type) = vx(ion_type,ion_index) - vx_neu(neutral_type)
            vy_rel(neutral_type) = vy(ion_type,ion_index) - vy_neu(neutral_type)
            vz_rel(neutral_type) = vz(ion_type,ion_index) - vz_neu(neutral_type)
            v_rel(neutral_type)  = sqrt(   vx_rel(neutral_type) * vx_rel(neutral_type) + &
                                         & vy_rel(neutral_type) * vy_rel(neutral_type) + &
                                         & vz_rel(neutral_type) * vz_rel(neutral_type)     ) 
            ! For each neutral type, select the nearest velocity value for which collision prob. are tabulated
            v_index(neutral_type) = search_index(v_values_ions, n_v_values_ions, v_rel(neutral_type), .true.)
         enddo

         ! Calculate the sum of all collision frequencies for the different neutrals
         ! for each neutral, the frequency is taken considering the proper relative ion-neutral speed         
         coll_f_sum = 0
         do neutral_type = 1, n_types
            coll_f_sum = coll_f_sum + coll_f_tot_ions(ion_type, neutral_type, v_index(neutral_type))
         enddo
         ! Create an array of the probability limits for the scattering with each neutral type
         ! for each neutral, the frequency is taken considering the proper relative ion-neutral speed         
         do neutral_type = 1, n_types
            coll_f_ratio = coll_f_tot_ions(ion_type, neutral_type, v_index(neutral_type)) / coll_f_sum
            if (neutral_type == 1) then
               p_limits_neutral(neutral_type) = coll_f_ratio
            else
               p_limits_neutral(neutral_type) = p_limits_neutral(neutral_type-1) + coll_f_ratio
            endif
         enddo
         
         ! Decide the type of neutral target (it could still be a null collision, but this will be tested later)
         call random_number(rand2)
         do neutral_type = 1, n_types
            if ( rand2 <= p_limits_neutral(neutral_type) ) exit
         enddo ! neutral_type

         ! Decide if the collision was elastic, charge exchange or null
         call random_number(rand3)         
         if (rand3 < p_limits_ions(ion_type, neutral_type, 1, v_index(neutral_type))) then 
            process_index = 1
            text = 'ELASTIC ION '
            ! Calculate mass ratio: m_ion / m_neutral = (q / m_neutral) / (q / m_ion)
            m_ratio = cm_ratio(neutral_type) / cm_ratio(ion_type)
            call elastic_scattering_ion(ion_type, ion_index, neutral_type, m_ratio, &
                                       &vx_neu(neutral_type), vy_neu(neutral_type), vz_neu(neutral_type), &
                                       &vx_rel(neutral_type), vy_rel(neutral_type), vz_rel(neutral_type), &                 
                                       &v_rel(neutral_type))
         elseif (rand3 < p_limits_ions(ion_type, neutral_type, 2, v_index(neutral_type))) then 
            process_index = 2
            text = 'CHARGE EX   '
            call charge_exchange_scattering_ion(ion_type, ion_index, neutral_type, &
                                               &vx_neu(neutral_type), vy_neu(neutral_type), vz_neu(neutral_type), &
                                               &v_neu(neutral_type))            
         else            ! Null collision
            process_index = -1
            text = 'NULL ION   '
         endif 

         if ( (debug_level > 1) ) then ! .and. (text .ne. 'NULL ION   ') ) then
            call print_collision_report_ion(text, rand2, rand3, ion_type, ion_index, &
                                           &v_0,   vx_0,   vy_0,   vz_0, v2_0, &                 
                                           &v_neu, vx_neu, vy_neu, vz_neu, &
                                           &v_rel, vx_rel, vy_rel, vz_rel, &
                                           &v_index, &
                                           &neutral_type, p_limits_neutral, process_index)
         endif

      end subroutine ion_collisions

      
      ! +--------------------------------+
      ! | Ion-neutral elastic collisions |
      ! +--------------------------------+
      
      subroutine elastic_scattering_ion(ion_type, ion_index, neutral_type, m_ratio, &
                                       &vx_neu, vy_neu, vz_neu, &
                                       &vx_rel_0, vy_rel_0, vz_rel_0, vrel_0)

         ! To manage ion elastic collision we need to express velocity vectors in 3 different frames of reference (FOR)
         ! DISFOR) discharge frame of reference, in which both ion and neutral are moving,
         !         and the z axis is parallel to the electric field direction,
         !         the initial and final ion velocity vectors will be expresses in this FOR
         ! TARFOR) target frame of reference, in which the target (neutral particle) is at rest,
         !         and the z axis is parallel to the direction of the relative (ion-neutral) velocity vector before scattering.
         ! COMFOR) center of mass (COM) frame of reference, in which the COM of the ion-neutral system is at rest,
         !         and the z axis is parallel to the direction of the relative velocity between ion and neutral before scattering
         !         the elastic ion scattering is considered isotropic in this FOR
         ! I have avoided the use of the term "Laboratory frame of reference" because it may refer to both the DISFOR or the TARFOR,
         ! depending on the context

         ! Parameters
         integer,  intent(in) :: ion_type, ion_index, neutral_type    ! Indexes identifying the ion and the neutral type
         real(dp), intent(in) :: vx_neu, vy_neu, vz_neu               ! Velocity vector of the neutral particle in the DISFOR
         real(dp), intent(in) :: vx_rel_0, vy_rel_0, vz_rel_0, vrel_0 ! Velocity vector and speed of the ion in the TARFOR
         real(dp), intent(in) :: m_ratio                              ! Ratio of masses ion/target

         ! Local variables
         real(dp) :: costheta_comfor              ! Cosinus of the scattering angle in the COMFOR
         real(dp) :: theta_scatter_tarfor         ! Scattering angle (ion direction) in the TARFOR
         real(dp) :: theta_recoil_tarfor          ! Recoil angles (neutral direction) in the TARFOR
         real(dp) :: phi_scatter_tarfor           ! Scattering anomaly in the TARFOR
         real(dp) :: verx_tarfor                  ! Projection of scattered relative velocity versor on x-axis in the TARFOR
         real(dp) :: very_tarfor                  ! Projection of scattered relative velocity versor on y-axis in the TARFOR
         real(dp) :: verz_tarfor                  ! Projection of scattered relative velocity versor on z-axis in the TARFOR
         real(dp) :: verx_disfor                  ! Projection of scattered relative velocity versor on x-axis in the DISFOR
         real(dp) :: very_disfor                  ! Projection of scattered relative velocity versor on y-axis in the DISFOR
         real(dp) :: verz_disfor                  ! Projection of scattered relative velocity versor on z-axis in the DISFOR
         real(dp) :: vx_rel_1, vy_rel_1, vz_rel_1 ! Relative velocity vector ion-neutral after scattering in the DISFOR
         real(dp) :: vrel_1                       ! Relative speed ion-neutral after scattering in the DISFOR
         real(sp) :: rand
        
         ! Generate an isotropically random scattering angle in the COMFOR
         call random_number(rand)
         costheta_comfor = 2 * rand - 1

         ! Transform the scattering and recoil angles from the COMFOR to the TARFOR
         ! recoil angle in the TARFOR will not be used because we do not follow the neutral fate
         call scattering_elastic_comfor2tarfor(theta_scatter_tarfor, theta_recoil_tarfor, costheta_comfor, m_ratio)

         ! Generate a random scattering anomaly in the TARFOR
         call random_number(rand)
         phi_scatter_tarfor = 2 * PI * rand

         ! Calculate the components of the scattered relative velocity versor in the TARFOR
         verx_tarfor = sin(theta_scatter_tarfor) * cos(phi_scatter_tarfor)
         very_tarfor = sin(theta_scatter_tarfor) * sin(phi_scatter_tarfor)
         verz_tarfor = cos(theta_scatter_tarfor)
         
         ! Trasform the scattered relative velocity versor from TARFOR to DISFOR coordianates
         call scattering_tarfor2disfor(verx_disfor, very_disfor, verz_disfor, &
                                      &vx_rel_0, vy_rel_0, vz_rel_0, vrel_0, &
                                      &verx_tarfor, very_tarfor, verz_tarfor, &
                                      &.false.)

         ! Calculate the modulus of the scattered relative velocity vector
         vrel_1   = vrel_0 * sqrt((1 - e_loss_ions(ion_type,neutral_type) * (1 - costheta_comfor)))         

         ! Calculate the scattered relative velocity vector in the DISFOR
         vx_rel_1 = vrel_1 * verx_disfor
         vy_rel_1 = vrel_1 * very_disfor         
         vz_rel_1 = vrel_1 * verz_disfor
         
         ! Calculate the new components of the ion velocity vector in the DISFOR     
         vx(ion_type,ion_index) = vx_rel_1 + vx_neu
         vy(ion_type,ion_index) = vy_rel_1 + vy_neu
         vz(ion_type,ion_index) = vz_rel_1 + vz_neu
         v(ion_type,ion_index)  = sqrt( vx(ion_type,ion_index) * vx(ion_type,ion_index) + &
                                       &vy(ion_type,ion_index) * vy(ion_type,ion_index) + &
                                       &vz(ion_type,ion_index) * vz(ion_type,ion_index)     )

      end subroutine elastic_scattering_ion

      
      ! +----------------------------------------+
      ! | Ion-neutral charge-exchange collisions |
      ! +----------------------------------------+
      
      subroutine charge_exchange_scattering_ion(ion_type, ion_index, neutral_type, vx_neu, vy_neu, vz_neu, v_neu)

         ! Parameters
         integer,  intent(in) :: ion_type, ion_index, neutral_type
         real(dp), intent(in) :: vx_neu, vy_neu, vz_neu, v_neu

         ! Reduce ion velocity to the one of the neutral (atom / molecule)
         v( ion_type,ion_index) = v_neu
         v2(ion_type,ion_index) = v_neu * v_neu
         vx(ion_type,ion_index) = vx_neu
         vy(ion_type,ion_index) = vy_neu
         vz(ion_type,ion_index) = vz_neu

      end subroutine charge_exchange_scattering_ion


      ! +---------------------------------------+
      ! | Electron-ion recombination management |
      ! +---------------------------------------+     

      subroutine electron_ion_recombination(electron_index)

         ! Parameters
         integer, intent(in) :: electron_index

         ! Local variables
         character*13  :: text
         real(dp)      :: v2_0
         real(dp)      :: fdt_recomb_total  ! Total recombination frequency multiplied by timestep
         real(dp)      :: p_recomb_abs      ! Absolute recombination probability 
         real(dp)      :: p_recomb_rel      ! Relative recombination probability, conditioned by the absence of electron collision
         real(dp)      :: p_realcoll_abs    ! Absolute probability of a real collision (not null)
         real(dp)      :: p_realcoll_rel    ! Relative probability of a real collision (not null),
                                            !     conditioned by the presence of electron collision (including null one)
         real(dp)      :: p_limit 
         real(sp)      :: rand1, rand2
         integer       :: ion_type, neutral_type, process_index, diss_type, v_index 

         ! Initialize some variables for debug purposes
         text = '------------'
         neutral_type  = -1
         diss_type     = -1
         process_index = -1
         v2_0          = v2(0, electron_index)      

         ! Find the index for the nearest velocity for which cross sections are tabulated
         v_index = search_index(v_values, n_v_values, v(0,electron_index), .true.)

         ! Calculate total recombination frequency multiplied by timestep
         ! f(v) * dt = n * k(v) * dt = n * s(v) * v * dt
         fdt_recomb_total = 0
         do ion_type = 1, n_types
            do diss_type = 1, n_disspro(ion_type)
               fdt_recomb_total = fdt_recomb_total + ion_density(ion_type) * &!1.0D16 * &
                                                    &kdt_recomb(ion_type, diss_type, v_index)
            enddo
         enddo

         ! Get the relative probability that a real scattering (not null) happened,
         ! when a scattering event (including null ones) was detected
         p_realcoll_rel = p_limits(n_limits, v_index)
         ! Calculate the absolute probability of a real collision (not null)
         p_realcoll_abs = p_realcoll_rel *  p_collision
         ! Calculate the absolute probability of recombination
         ! p = 1 - exp(- n * s(v) * v * dt) = 1 - exp(- f(v) * dt)
         p_recomb_abs = 1.D0 - exp(- fdt_recomb_total)
         ! Calculate the probability of a recombination process, conditioned by the absence of a real collision
         p_recomb_rel = p_recomb_abs / (1.D0 - p_realcoll_abs)

         ! Decide if the electron recombines during this timestep
         call random_number(rand1)
         if (rand1 < p_recomb_rel) then
            if (debug_level > 1) then
               print *, 'n_totpro                       :', n_totpro
               print *, 'p_collision                    :', p_collision
               print *, 'p_realcoll_abs, p_realcoll_rel :', p_realcoll_abs, p_realcoll_rel
               print *, 'fdt_recomb_total               :', fdt_recomb_total
               print *, 'p_recomb_abs,   p_recomb_rel   :', p_recomb_abs,   p_recomb_rel
               print *, 'rand1                          :', rand1
            endif
            ! Check the recombination channel (type of recombined ion and of recombination process)
            call random_number(rand2)
            p_limit = 0
            do ion_type = 1, n_types
               do diss_type = 1, n_disspro(ion_type)
                  p_limit = p_limit + ion_density(ion_type) * kdt_recomb(ion_type, diss_type, v_index) / fdt_recomb_total
                  if (rand2 <= p_limit) exit                     
               enddo
            enddo
            ! Increase the proper collision counter
            coll_rec(ion_type, diss_type) = coll_rec(ion_type, diss_type) + weight(0)

            if (debug_level > 1) then
               print *, 'rand2                          :', rand2
               print *, 'ion_type                       :', ion_type
               print *, 'diss_type                      :', diss_type
            endif

            call recombination_remove_particles(electron_index, ion_type)

            if (debug_level > 1) then
               text = 'RECOMBINATION'
               call print_collision_report_electron(text, rand2, electron_index, v2_0, v_index, &
                                                   &neutral_type, process_index, diss_type)
            endif ! (debug_level > 1)
         endif ! (rand1 < p_recomb_rel)

      end subroutine electron_ion_recombination


      subroutine recombination_remove_particles(electron_index, ion_type)

         ! Parameters
         integer, intent(in) :: electron_index, ion_type

         ! Local variables
         integer  :: n, n_ions_removed, ion_index
         real(dp) :: weight_ratio
         real(sp) :: rand

         ! Remove the recombined electron
         call remove_particle(0, electron_index)

         ! Calculate the number of simulation ions to remove
         ! It will be equal to the ratio between the weights of electron and ion type
         ! if it is lower then one, then it becomes the probability that an ion is removed
         weight_ratio = weight(0)/weight(ion_type)
         if (weight_ratio >= 1) then
            n_ions_removed = nint(weight_ratio)  
         else
            call random_number(rand)
            if (rand < weight_ratio) then
               n_ions_removed = 1
            else
               n_ions_removed = 0
            endif
         endif

         ! Remove the proper number of ions
         if (n_ions_removed > 0) then
            n = 0
            ! Search for active ions to remove
            do ion_index = 1, n_particles
               if (isactive(ion_type, ion_index).eqv..true.) then 
                  call remove_particle(ion_type, ion_index) 
                  n = n + 1
                  if (n >= n_ions_removed) exit
               endif
            enddo ! ion_index = 1, n_particles
            ! enddo ! i = 1, n_ions_removed
            if (debug_level > 1) then
               print *, 'ion_index                                  :', ion_index
               print *, 'n_ions                                     :', count( isactive(ion_type,:))
            endif
         endif ! (n_ions_removed > 0)

       end subroutine  recombination_remove_particles

      
      ! +-------------------------------+
      ! | Remove single electron or ion |
      ! +-------------------------------+

      subroutine remove_particle(t, m)

         ! Parameters
         integer, intent(in) :: t, m

         if ( (t<=n_types) .and. (m<=n_particles) ) then
            isactive(t, m) = .false.
            restart(t, m)  = .true.
            x(t, m)        = 0
            y(t, m)        = 0
            z(t, m)        = 0
            vx(t, m)       = 0
            vy(t, m)       = 0
            vz(t, m)       = 0
            v(t, m)        = 0
            v2(t, m)       = 0
         endif

      end subroutine remove_particle

      
      ! +----------------------------------------------------+
      ! | Print a report of the collisions (debug mode only) |
      ! +----------------------------------------------------+

      subroutine print_collision_report_electron(s, rand2, particle_index, v2_0, v_index, &
                                                &neutral_type, process_index, diss_type)      

         ! Parameters
         character*12, intent(in) :: s
         real(dp),     intent(in) :: v2_0
         real(sp),     intent(in) :: rand2
         integer,      intent(in) :: particle_index, v_index, neutral_type, process_index, diss_type

         ! Local variables
         real(dp) :: e0, e1
         integer  :: l

         ! Calculate electron energy (in eV) before and after scattering
         e0 = 0.5 * v2_0 / abs(cm_ratio(0))
         e1 = 0.5 * v2(0,particle_index) / abs(cm_ratio(0))

         print *, " "
         print *, "ELECTRON COLLISION DETAILS"
         print *, "Collision type    = ", s         
         print *, "Electron index    = ", particle_index
         if (neutral_type >= 0) then 
            print *, "Target type       = ", neutral_type
            print *, "Energy loss       = ", e_loss(neutral_type)
         endif
         print *, "rand2             = ", rand2
         print *, "p_limits          = "
         do l = 1, n_totpro
            print *, p_limits(l, v_index)
         enddo
         print *, "process_index     = ", process_index, "( of", n_totpro,")"
         print *, "plimit1_exc       = ", limit1_exc
         print *, "plimit1_diss      = ", limit1_diss            

         if (diss_type >=0) print *, "Dissoc  type      = ", diss_type
         print *, "energy0 / eV      = ", e0
         print *, "energy1 / eV      = ", e1

         print *, "vx, vy, vz, v, v2 = ", vx(type_index,particle_index),  &
                                      &vy(type_index,particle_index), &
                                      &vz(type_index,particle_index), &
                                      &v(type_index,particle_index), &
                                      &v2(type_index,particle_index) 
         if (debug_level >1) call pause

      end subroutine print_collision_report_electron


      subroutine print_collision_report_ion(s, rand2, rand3, type_index, particle_index, &
                                           &v_0  , vx_0,   vy_0,   vz_0, v2_0, &
                                           &v_neu, vx_neu, vy_neu, vz_neu, &
                                           &v_rel, vx_rel, vy_rel, vz_rel, &
                                           &v_index, &
                                           &neutral_type, p_limits_neutral, process_index)

         ! Parameters
         character*12,                   intent(in) :: s
         real(dp),                       intent(in) :: v_0  , vx_0,   vy_0,   vz_0, v2_0
         real(dp), dimension(1:n_types), intent(in) :: v_neu, vx_neu, vy_neu, vz_neu
         real(dp), dimension(1:n_types), intent(in) :: v_rel, vx_rel, vy_rel, vz_rel
         integer,  dimension(1:n_types), intent(in) :: v_index
         real(dp), dimension(1:n_types), intent(in) :: p_limits_neutral
         real(sp),                       intent(in) :: rand2, rand3
         integer,                        intent(in) :: type_index, particle_index, neutral_type, process_index

         ! Local variables
         real(dp) :: e0, e1
         integer  :: l

         ! Calculate ion energy (in eV) before and after scattering
         e0 = 0.5 * v2_0 / abs(cm_ratio(type_index))
         e1 = 0.5 * v2(type_index,particle_index) / abs(cm_ratio(type_index))

         print *, " "
         print *, "ION COLLISION DETAILS"
         print *, "Collision type    = ", s         
         print *, "Ion type, index   = ", type_index, particle_index
         if (neutral_type >= 0) then 
            print *, "Target type       = ", neutral_type
            print *, "Energy loss       = ", e_loss_ions(type_index, neutral_type)
         endif
         print *, "rand2             = ", rand2
         print *, "p_limits neutral  = ", p_limits_neutral
         print *, "neutral           = ", neutral_type
         print *, "v_rel             = ", v_rel(neutral_type)
         print *, "v_index           = ", v_index(neutral_type)
         print *, "rand3             = ", rand3
         print *, "p_limits proc     = ", p_limits_ions(type_index, neutral_type, 1, v_index(neutral_type)), &
                                         &p_limits_ions(type_index, neutral_type, 2, v_index(neutral_type))
         print *, "proc index        = ", process_index
         print *, ""
         print *, "Neutral and relative (ion-neutral) velocity vectors"
         do l = 1, n_types
            print *, "type: ", l,        "v_mean:", mean_speed(l), "v_index:", v_index(l)
            print *, "v_neu:", v_neu(l), "vn_x:",   vx_neu(l),     "vn_y:",    vy_neu(l), "vn_z:", vz_neu(l)
            print *, "v_rel:", v_rel(l), "vr_x:",   vx_rel(l),     "vr_y:",    vy_rel(l), "vr_z:", vz_rel(l)
         enddo
         print *, ""
         print *, "Initial and final ion velocity vector"
         print *, "v0:   ", v_0,      "v0_x:",   vx_0,          "v0_y:",    vy_0,     "v0_z:", vz_0,     "v0^2:", v2_0
         print *, "v1:   ", v(type_index,particle_index), &
                 &"v1_x:", vx(type_index,particle_index), &
                 &"v1_y:", vy(type_index,particle_index), &
                 &"v1_z:", vz(type_index,particle_index), &
                 &"v1^2:", v2(type_index,particle_index)         
         if (debug_level >1) call pause

      end subroutine print_collision_report_ion
      

   end subroutine check_collisions


end module f_particle_collider

