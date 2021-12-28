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
                              &isactive, restart,&
                              &cm_ratio, weight,  &         
                              &v_ratio, v_ratio_ions, &    
                              &p_collision, p_collision_ions, &
                              &p_limits, p_limits_ions, & 
                              &n_limits, n_limits_ions, &
                              &kdt_recomb, ion_density, isactive_recomb, &
                              &min_scattered, &
                              &n_totpro, n_totpro_ions, &
                              &en_ion, en_diss, n_disspro, mean_speed, &    
                              &coll_null, coll_ela, coll_ion, coll_dis, coll_rec, &              
                              &v_values, v_values_ions, n_v_values, n_v_values_ions, &
                              &x_part_added, y_part_added, z_part_added, w_part_added, n_part_added, &
                              &n_types, n_particles, nmax_disspro, &
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
      real(dp), dimension(1:n_types),                intent(in)    :: v_ratio              ! Ratio between electron velocity values
                                                                                           !    measured before and after an elastic scattering
      real(dp), dimension(1:n_types,1:n_types),      intent(in)    :: v_ratio_ions         ! Ratio between ion velocity values
                                                                                           !    measured before and after an elastic scattering 
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
      real(dp), dimension(1:n_types,1:n_totpro_ions, &
                         &1:n_v_values_ions),        intent(in)    :: p_limits_ions        ! Values used to select collision type for ions
      integer,                                       intent(in)    :: n_limits             ! Number of probability limits for el collisions
      integer,                                       intent(in)    :: n_limits_ions        ! Number of probability limits for ions collisions

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
      integer,                                       intent(in)    :: n_totpro_ions        ! Total maximum possible number of scattering
                                                                                           !  processes for ions: elastic + charge exchange
      real(dp),                                      intent(inout) :: coll_null            ! Number of observed null collisions
      real(dp), dimension(1:n_types),                intent(inout) :: coll_ela             ! Number of observed elastic collisions
                                                                                           !    for each neutral gas type
      real(dp), dimension(1:n_types),                intent(inout) :: coll_ion             ! Number of observed ionization collisions for each 
                                                                                           !    neutral gas type
      real(dp), dimension(1:n_types,1:nmax_disspro), intent(inout) :: coll_dis             ! Number of observed dissociation collisions for each
                                                                                           !    neutral gas type and dissociation process 
      real(dp), dimension(1:n_types,1:nmax_disspro), intent(inout) :: coll_rec             ! Number of observed electron-ion recombinations
                                                                                           !    for each neutral gas type
                                                                                           !    and dissociative recombination process 
      real(dp), dimension(1:n_types),                intent(in)    :: en_ion               ! Ionization energy of each neutral type
      real(dp), dimension(1:n_types,1:nmax_disspro), intent(in)    :: en_diss              ! Dissociation energy of each neutral type and 
                                                                                           !    each dissociation process
      integer,  dimension(1:n_types),                intent(in)    :: n_disspro            ! Number of dissociation processes available 
                                                                                           !    for each neutral type
      integer,                                       intent(in)    :: nmax_disspro         ! Maximum number of dissociation types
      real(dp), dimension(1:n_types),                intent(in)    :: mean_speed           ! Mean speed of neutrals molecules/ m*s**-1

      real(dp), dimension(0:n_types, 1:n_particles), intent(inout) :: x_part_added, &
                                                                     &y_part_added, &
                                                                     &z_part_added         ! Positions of particles that need to be added
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
         !    -> collision probablility by number of active particles
         n_scattered = nint(collision_probability * n_active)

!         if (type_index == 0) print *,'***** type ', type_index, ' --> ', &
!                                    & nint(collision_probability*100), n_active, n_scattered !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! If the number of scattered particles is high enough, apply the faster method, otherwise the slower one
         if ( (n_scattered > min_scattered) ) then  
            call manage_collisions_fast_method(type_index, n_active, n_scattered)
         else 
            call manage_collisions_slow_method(type_index, collision_probability)          
         endif
      enddo

   contains

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
                  if ( (debug_level > 1) ) print *, 'COLLISION DETECTED (FAST METHOD) WITH RECOMBINATION' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
               do particle_index = 1, n_particles
                  if (check_electron_recombination(particle_index)) call electron_ion_recombination(particle_index)
               enddo          
            ! Check electron collisions, excluding electron/ion recombination
            else
               do i = 1, n_scattered
                  if ( (debug_level > 1) ) print *, 'COLLISION DETECTED (FAST METHOD) WITHOUT RECOMBINATION' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  particle_index = scattered_indexes_absolute(i)
                  collision = .false.
                  call electron_collisions(particle_index, collision)
                  if (.not.(collision)) coll_null = coll_null + weight(0)
               enddo ! i 
            endif

            
!            ! Initialize the array used to keep track of scattered electrons, since only unscattered ones may recombine
!            check_electron_recombination = isactive(0, :) ! Only active electrons can recombine
!            do i = 1, n_scattered
!               if ( (debug_level > 1) ) print *, 'COLLISION DETECTED (FAST METHOD)' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               particle_index = scattered_indexes_absolute(i)
!               collision = .false.
!               call electron_collisions(particle_index, collision)
!               ! If the electron had a non-null collision, will not be checked for recombination
!               if (collision) then 
!                  check_electron_recombination(particle_index) = .false.  ! Since the electron has scattered, it can not recombine
!               else
!                  coll_null = coll_null + weight(0)
!               endif
!            enddo ! i 
!            ! Check unscattered electrons for recombination with ions
!            if (isactive_recomb) then
!               do particle_index = 1, n_particles
!                  if (isactive_recomb .and. check_electron_recombination(particle_index)) then
!                  if (check_electron_recombination(particle_index)) call electron_ion_recombination(particle_index)
!               enddo
!            endif   

         ! Manage ions collisions
         else
            do i = 1, n_scattered
               if ( (debug_level > 1) ) print *, 'COLLISION DETECTED (FAST METHOD)' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               particle_index = scattered_indexes_absolute(i)
               call ion_collisions(type_index,particle_index)
            enddo ! particle_index 
         endif !(type_index == 0)

         deallocate(scattered_indexes_relative, scattered_indexes_absolute)

      end subroutine manage_collisions_fast_method


      subroutine electron_collisions(electron_index, collision_happened)
 
         ! Parameters
         integer, intent(in)    :: electron_index
         logical, intent(inout) :: collision_happened

         ! Local variables
         character*12 :: string
         real(dp)     :: v2_0
         integer      :: l, v_index, process_index, neutral_type, diss_type
         real(sp)     :: rand2 

         ! Initialize some variables for debug purposes
         string = '------------'
         neutral_type = -1
         diss_type    = -1
         v2_0         = v2(0, electron_index)      

         ! Select the nearest velocity value for which collision prob. are tabulated
         v_index = search_index(v_values,  n_v_values, v(0,electron_index), .true.)

         ! Decide which type of scattering happened
         call random_number(rand2)
         do process_index = 1, n_totpro
            if ( rand2 <= p_limits(process_index, v_index) ) exit
         enddo

         if (process_index <= n_types) then 
            neutral_type = process_index
            call elastic_scattering_electron(electron_index, neutral_type)        
            collision_happened = .true.
            string = 'ELASTIC     '
         elseif (process_index <= 2 * n_types) then
            neutral_type = process_index - n_types
            call ionization_scattering_electron(electron_index, neutral_type)
            collision_happened = .true.
            string = 'IONIZATION  '
         elseif (process_index <= n_limits) then
            ! Identify which type of neutral and which dissociation process are involved
            l = process_index - n_types - n_types
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
            string = 'DISSOCIATION'
         else            ! Null collision: do nothing  
            collision_happened = .false.
            string = 'NULL       '
         endif  ! (process_index <= n_types) 
   
         if ( (debug_level > 1) ) then
            call print_collision_report(string, rand2, 0, electron_index, v2_0, v_index, &
                                        &neutral_type, process_index, diss_type)
         endif 

      end subroutine electron_collisions
 

      subroutine electron_ion_recombination(electron_index)

         ! Parameters
         integer, intent(in) :: electron_index

         ! Local variables
         character*13  :: string
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
         string = '------------'
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
               string = 'RECOMBINATION'
               call print_collision_report(string, rand2, 0, electron_index, v2_0, v_index, &
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


      subroutine ion_collisions(ion_type, ion_index)

         ! Parameters
         integer, intent(in) :: ion_type, ion_index

         ! Local variables
         character*12 :: string
         real(dp)     :: v2_0
         integer      :: v_index, process_index, neutral_type
         real(sp)     :: rand2

         ! Initialize some variables for debug purposes
         string = '------------'
         neutral_type = -1
         v2_0         = v2(ion_type, ion_index)

         ! Select the nearest velocity value for which collision prob. are tabulated
         v_index = search_index(v_values_ions, n_v_values_ions, v(ion_type,ion_index), .true.)

         ! Decide which type of scattering happened
         call random_number(rand2)
         do process_index = 1, n_totpro_ions
            if ( rand2 <= p_limits_ions(ion_type, process_index, v_index) ) exit
         enddo ! process_index

         if (process_index <= n_types ) then 
            neutral_type = process_index  
            call elastic_scattering_ion(ion_type, ion_index, neutral_type)
            string = 'ELASTIC ION '   
         elseif (process_index < n_limits_ions) then
            neutral_type = process_index - n_types
            call charge_exchange_scattering_ion(ion_type, ion_index, neutral_type)
            string = 'CHARGE EX   '  
         else            ! Null collision  
            neutral_type = -1
            string = 'NULL ION   '
         endif  ! (process_index < 1 )

         if ( (debug_level > 1) ) then ! .and. (string .ne. 'NULL ION   ') ) then
            call print_collision_report(string, rand2, ion_type, ion_index, v2_0, &
                                       &v_index, neutral_type, process_index, -1)
         endif

      end subroutine ion_collisions


      subroutine elastic_scattering_electron(electron_index, neutral_type)

         ! Parameters
         integer, intent(in) :: electron_index, neutral_type

         ! Local variables
         real(dp) :: e_energy_au       ! electron energy in atomic units
         real(dp) :: verx_scatter      ! Projection of scattered velocity versor on x-axis in the frame of reference of the target
         real(dp) :: very_scatter      ! Projection of scattered velocity versor on y-axis in the frame of reference of the target
         real(dp) :: verz_scatter      ! Projection of scattered velocity versor on z-axis in the frame of reference of the target
         real(dp) :: verx1,very1,verz1 ! Projections on x,y,z axis of the velocity versor after scattering
         logical  :: debug

         if (debug_level > 1) then
            debug = .true.
         else
            debug = .false.
         endif

         ! Calculate particle energy before scattering (in atomic units)
         e_energy_au = E_V2TOENERGY_AU * v2(0, electron_index)             
 
         ! Calculate scattering angles (respect to particle velocity before scattering)
         call elastic_scattering_electron_neutral(verx_scatter, very_scatter, verz_scatter, e_energy_au, .true., debug)
         ! call isotropic_scattering(verx_scatter, very_scatter, verz_scatter, debug)

         ! Calculate the direction of particle velocity after scattering
         call scattering_3d(verx1,very1,verz1, &
                           &vx(0,electron_index), &
                           &vy(0,electron_index), &
                           &vz(0,electron_index), &
                           &v( 0,electron_index),  &
                           &verx_scatter, very_scatter, verz_scatter, &
                           &debug)

         ! Calculate new electron velocity (due to kinetic energy transferred to the neutral particle)
         v (0,electron_index) = v_ratio(neutral_type) * v(0,electron_index)
         v2(0,electron_index) = v(0,electron_index)   * v(0,electron_index)  

         ! Calculate the new components of velocity vector
         vx(0,electron_index) = v(0,electron_index) * verx1
         vy(0,electron_index) = v(0,electron_index) * very1
         vz(0,electron_index) = v(0,electron_index) * verz1

         ! Increase collisions counter
         coll_ela(neutral_type) = coll_ela(neutral_type) + weight(0)

      end subroutine elastic_scattering_electron


      subroutine ionization_scattering_electron(electron_index, neutral_type)

         ! Parameters
         integer, intent(in) :: electron_index, neutral_type

         ! Local variables
         real(dp) :: e_energy_j        ! electron energy in joules
         real(dp) :: verx_scatter      ! Projection of scattered velocity versor on x-axis in the frame of reference of the target
         real(dp) :: very_scatter      ! Projection of scattered velocity versor on y-axis in the frame of reference of the target
         real(dp) :: verz_scatter      ! Projection of scattered velocity versor on z-axis in the frame of reference of the target
         real(dp) :: verx1,very1,verz1 ! Projections on x,y,z axis of the velocity versor after scattering
!         integer  :: i, n_ions_added
!         real(dp) :: weight_ratio
!         real(sp) :: rand
!         real(dp) :: p
         logical  :: debug

         if (debug_level > 1) then
            debug = .true.
         else
            debug = .false.
         endif
 
         ! Generate a new random direction for particle velocity after scattering
         call isotropic_scattering(verx_scatter, very_scatter, verz_scatter, debug)

         ! Calculate the direction of particle velocity after scattering
         call scattering_3d(verx1,very1,verz1, &
                           &vx(0,electron_index), &
                           &vy(0,electron_index), &
                           &vz(0,electron_index), &
                           &v(0,electron_index),  &
                           &verx_scatter, very_scatter, verz_scatter, &
                           &debug)

         ! Calculate particle energy after scattering in joules; E_CHARGE is negative
         e_energy_j = 0.5 * E_MASS * v2(0, electron_index) + en_ion(neutral_type) * E_CHARGE  
         if (e_energy_j < 0) e_energy_j = 0

         ! Calculate particle velocity after scattering
         v2(0,electron_index) = 2.0 * e_energy_j / E_MASS
         v( 0,electron_index) = sqrt(v2(0,electron_index))

         ! Calculate the new components of velocity vector
         vx(0,electron_index) = v(0,electron_index) * verx1
         vy(0,electron_index) = v(0,electron_index) * very1
         vz(0,electron_index) = v(0,electron_index) * verz1

         ! Increase collisions counter
         coll_ion(neutral_type) = coll_ion(neutral_type) + weight(0)

         call ionization_add_particles(electron_index, neutral_type)

      end subroutine ionization_scattering_electron


      subroutine ionization_add_particles(electron_index, neutral_type)

         ! Parameters
         integer, intent(in) :: electron_index, neutral_type

         ! Local variables
         integer  :: n_ions_added
         real(dp) :: weight_ratio
         real(sp) :: rand
!         logical  :: debug

         ! Add an electron
         n_part_added(0)                 = n_part_added(0) + 1    ! will add an electron
         w_part_added(0,n_part_added(0)) = 1                      ! will add a single computational electron for each ionization event
         x_part_added(0,n_part_added(0)) = x(0,electron_index)    ! save the position (x-component)
         y_part_added(0,n_part_added(0)) = y(0,electron_index)    ! save the position (y-component)
         z_part_added(0,n_part_added(0)) = z(0,electron_index)    ! save the position (z-component)

         ! Calculate the number of ions to add
         ! It will be equal to the ratio between the weights of electron and ion type
         ! if it is lower then one, then it becomes the probability that an ion is added
         weight_ratio = weight(0)/weight(neutral_type)
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
            n_part_added(neutral_type)                            = n_part_added(neutral_type) + 1  ! will add a group of comput. ions 
                                                                                                    !   of type neutral_type
            w_part_added(neutral_type,n_part_added(neutral_type)) = n_ions_added                    ! number of computational ions of this type
                                                                                                    !   to add
            x_part_added(neutral_type,n_part_added(neutral_type)) = x(0,electron_index)             ! save the position (x-component) 
            y_part_added(neutral_type,n_part_added(neutral_type)) = y(0,electron_index)             ! save the position (y-component) 
            z_part_added(neutral_type,n_part_added(neutral_type)) = z(0,electron_index)             ! save the position (z-component) 
         endif

      end subroutine ionization_add_particles


      subroutine dissociation_scattering_electron(electron_index, neutral_type, diss_type)

         ! Parameters
         integer, intent(in) :: electron_index, neutral_type, diss_type

         ! Local variables
         real(dp) :: e_energy_j        ! electron energy in joules
         real(dp) :: verx_scatter      ! Projection of scattered velocity versor on x-axis in the frame of reference of the target
         real(dp) :: very_scatter      ! Projection of scattered velocity versor on y-axis in the frame of reference of the target
         real(dp) :: verz_scatter      ! Projection of scattered velocity versor on z-axis in the frame of reference of the target
         real(dp) :: verx1,very1,verz1 ! Projections on x,y,z axis of the velocity versor after scattering

         ! Generate a new random direction for particle velocity after scattering
         call isotropic_scattering(verx_scatter, very_scatter, verz_scatter, .false.)

         ! Calculate the direction of particle velocity after scattering
         call scattering_3d(verx1,very1,verz1, &
                           &vx(0,electron_index), &
                           &vy(0,electron_index), &
                           &vz(0,electron_index), &
                            &v(0,electron_index), &
                           &verx_scatter, very_scatter, verz_scatter, &
                           &.false.)

         ! Calculate particle energy after scattering in joules; E_CHARGE is negative
         e_energy_j = 0.5 * E_MASS * v2(0, electron_index) + en_diss(neutral_type, diss_type) * E_CHARGE
         if (e_energy_j < 0) e_energy_j = 0 

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


      subroutine elastic_scattering_ion(ion_type, ion_index, neutral_type)

         ! Parameters
         integer, intent(in) :: ion_type, ion_index, neutral_type

         ! Local variables
         real(dp) :: verx_scatter      ! Projection of scattered velocity versor on x-axis in the frame of reference of the target
         real(dp) :: very_scatter      ! Projection of scattered velocity versor on y-axis in the frame of reference of the target
         real(dp) :: verz_scatter      ! Projection of scattered velocity versor on z-axis in the frame of reference of the target
         real(dp) :: verx1,very1,verz1 ! Projections on x,y,z axis of the velocity versor after scattering

         ! Calculate elastic scattering angles
         call isotropic_scattering(verx_scatter, very_scatter, verz_scatter, .false.)

         ! Calculate the direction of particle velocity after scattering
         call scattering_3d(verx1,very1,verz1, &
                           &vx(ion_type,ion_index), &
                           &vy(ion_type,ion_index), &
                           &vz(ion_type,ion_index), &
                           &v(ion_type,ion_index),  &
                           &verx_scatter, very_scatter, verz_scatter, &
                           &.false.)

         ! Reduce particle velocity (due to kinetic energy transferred to the neutral particle)
         v( ion_type,ion_index) = v_ratio_ions(ion_type,neutral_type) * v(ion_type,ion_index)
         v2(ion_type,ion_index) = v(ion_type,ion_index) * v(ion_type,ion_index)  

         ! Calculate the new components of velocity vector
         vx(ion_type,ion_index) = v(ion_type,ion_index) * verx1
         vy(ion_type,ion_index) = v(ion_type,ion_index) * very1
         vz(ion_type,ion_index) = v(ion_type,ion_index) * verz1

      end subroutine elastic_scattering_ion


      subroutine charge_exchange_scattering_ion(ion_type, ion_index, neutral_type)

         ! Parameters
         integer, intent(in) :: ion_type, ion_index, neutral_type

         ! Local variables
         real(dp) :: vratio  ! Ratio of ion velocities before and after charge exchange

         ! Reduce ion velocity to the thermal speed of the atom / molecule
         vratio = mean_speed(neutral_type) / v(ion_type,ion_index)

         ! Calculate velocity modulus and components 
         v( ion_type,ion_index) = vratio * v(ion_type,ion_index)
         v2(ion_type,ion_index) = v(ion_type,ion_index) * v(ion_type,ion_index)

         vx(ion_type,ion_index) = vratio * vx(ion_type,ion_index)
         vy(ion_type,ion_index) = vratio * vy(ion_type,ion_index) 
         vz(ion_type,ion_index) = vratio * vz(ion_type,ion_index)

      end subroutine charge_exchange_scattering_ion


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


      subroutine print_collision_report(s, rand2, type_index, particle_index, v2_0, v_index, &
                                       &neutral_type, process_index, diss_type)      

         ! Parameters
         character*12, intent(in) :: s
         real(dp),     intent(in) :: v2_0
         real(sp),     intent(in) :: rand2
         integer,      intent(in) :: type_index, particle_index, v_index, neutral_type, process_index, diss_type

         ! Loval variables
         real(dp) :: e0, e1
         integer  :: l

         ! Calculate particle energy (in eV) before and after scattering
         e0 = 0.5 * v2_0 / abs(cm_ratio(type_index))
         e1 = 0.5 * v2(type_index,particle_index) / abs(cm_ratio(type_index))

         print *, " "
         print *, "COLLISION DETAILS"
         print *, "Particle type  = ", type_index
         print *, "Particle index = ", particle_index
         print *, "Collision type:  ", s
         if (neutral_type >= 0) then 
            print *, "Target type   = ", neutral_type
            if (type_index == 0) then 
               print *, "v ratio        = ", v_ratio(neutral_type)
            else
               print *, "v ratio        = ", v_ratio_ions(type_index, neutral_type)
            endif
         endif
         print *, "rand2          = ", rand2
         print *, "p_limits       = "
         if (type_index == 0) then
            do l = 1, n_totpro
               print *, p_limits(l, v_index)
            enddo
            print *, "process_index  = ", process_index, "( of", n_totpro,")"
         else
            do l = 1, n_totpro_ions
               print *, p_limits_ions(type_index, l, v_index)
            enddo
            print *, "process_index  = ", process_index, "( of", n_totpro_ions,")"
         endif
         if (diss_type >=0) print *, "Dissoc  type  = ", diss_type
         print *, "energy0 / eV   = ", e0
         print *, "energy1 / eV   = ", e1

         print *, "vx, vy, vz, v, v2 = ", vx(type_index,particle_index),  &
                                      &vy(type_index,particle_index), &
                                      &vz(type_index,particle_index), &
                                      &v(type_index,particle_index), &
                                      &v2(type_index,particle_index) 
         if (debug_level >1) call pause

      end subroutine print_collision_report


   end subroutine check_collisions


end module f_particle_collider

