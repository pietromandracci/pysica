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

module f_main

   use f_precision

   use f_io

   use f_pic

   use f_particle_mover

   use f_particle_collider

   use f_particle_manager

!   use omp_lib

   implicit none

contains

   subroutine simccp(x, y, z, &                                         ! particles positions (inout)
                   &vx, vy, vz, v,&                                     ! particles velocity components and modulus (inout)
                   &isactive, &                                         ! array used to identify active particles (inout) 
                   &restart, &                                          ! array used to restart LP scheme (inout) 
                   &cm_ratio, &                                         ! charge/mass ratio (signed) of each particle type (in)
                   &weight, &                                           ! computational weights (inout)
                   &rescale_factor, &                                   ! factor used for weight rescaling (in)                   
                   &v_values, &                                         ! speeds for which e- coll xsec are known (in)
                   &coll_freq, &                                        ! max value of e- collision frequencies over all speeds (in)
                   &p_limits, &                                         ! probability limits for electron collisions (in)
                   &n_limits, &                                         ! number of probability limits for electron collisions (in)
                   &limit1_exc, &                                       ! index of first probability limit for excitation (in)
                   &limit1_diss, &                                      ! index of first probability limit for dissociation (in)                 
                   &v_values_ions, &                                    ! speeds for which ion coll xsec are known (in) 
                   &coll_freq_ions, &                                   ! max value of ion coll frequencies (in)
                   &coll_f_tot_ions, &                                  ! toal collision frequency for each neutral type (all processes) (in)
                   &p_limits_ions, &                                    ! probability limits for ion scattering (in)      
                   &k_recomb, &                                         ! coll rate coeff for e-/ion dissociative recombination (in)
                   &isactive_recomb, &                                  ! set if electron/ion recombination processes are activated
                   &min_scattered, &                                    ! min mumber of scattering needed to apply many particles method (in)
                   &e_loss, e_loss_ions, &                              ! mean energy loss of e-, ions during e- scattering (in)
                   &en_ion, en_exc, en_diss, &                          ! ionization, excitation, dissociation energies (in)
                   &n_excpro, n_disspro, &                              ! number of excitation, dissociation processes (in)
                   &mean_speed, &                                       ! mean speed (in)                   
                   &se_coefficient, &                                   ! secondary emission coefficients for each neutral type (in)
                   &distance, length, Vmax, omega, phi, &               ! reactor and electric potential characteristics (in)
                   &lateral_loss, &                                     ! if true, electons and ions are lost atthe borded (in)
                   &rho, psi, &                                         ! electric charge density and electric potential (out) 
                   &average_current, &                                  ! average electric current measured at the electrodes (out)
                   &dt, duration, &                                     ! timestep (in), required duration of simulation (in)
                   &coll_null,coll_ela, coll_ion, &
                   &coll_exc, coll_dis, coll_rec, &                     ! collisions (inout)
                   &n_particles, n_v_values, n_totpro, &                ! arrays dimensions (hide)
                   &n_v_values_ions, &                                  ! arrays dimensions (hide)
                   &n_types, nmax_excpro, nmax_disspro, &               ! arrays dimensions (hide)
                   &n_cells, &                                          ! arrays dimensions (hide)
                   &debug_level)                                        ! used for debug purposes (in)

      !  Simulates evolution of the velocities of charged particles 
      !  (e.g. particles) moving in a variable electric field E(t) = Emax * sin( omega*t + phi )
      !  when they are subject to scattering with some collision
      !  centers uniformly distributed in space.
      !  Several types of scattering center may be considered and the dependence of scattering cross section
      !  from particle velocity must be given to the function.
      !  The densities of particles and collision centers are considered fixed in time.


      !f2py intent(in)        :: cm_ratio, rescale-factor, v_values, coll_freq, p_limits, n_limits, limit1_exc, limit1_diss
      !f2py intent(in)        :: v_values_ions, coll_freq_ions, coll_f_tot_ions, p_limists_ions
      !f2py intent(in)        :: k_recomb, isactive_recomb, min_scattered, e_loss, e_loss_ions, en_ion, en_exc, en_diss, n_excpro, n_disspro
      !f2py intent(in)        :: mean_speed, se_coefficient, distance, length, vmax, omega, phi, lateral_loss, dt, duration, debug_level
      !f2py intent(inout)     :: x, y, z, vx, vy, vz, v, weight
      !f2py intent(inout)     :: isactive, restart, coll_null, coll_ela, coll_ion, coll_exc, coll_dis, coll_rec
      !f2py intent(inout)     :: rho, psi
      !f2py intent(overwrite) :: average_current
      !f2py intent(hide)      :: n_particles, n_v_values, n_v_values_ions, n_totpro, n_types, nmax_excpro, nmax_disspro, n_cells      

      !  Parameters
      real(dp), dimension(0:n_types,1:n_particles),     intent(inout) :: x, y, z         ! Particles positions components
      real(dp), dimension(0:n_types,1:n_particles),     intent(inout) :: vx, vy, vz, v   ! Particles velocities: components and modulus
      logical,  dimension(0:n_types,1:n_particles),     intent(inout) :: isactive        ! Selects which particlesa are active
      logical,  dimension(0:n_types,1:n_particles),     intent(inout) :: restart         ! Selects which particles must restart 
                                                                                         !    leap-frog scheme at the next iteration
      real(dp), dimension(0:n_types),                   intent(in)    :: cm_ratio        ! Charge/mass ratios (signed) of each particle type
      real(dp), dimension(0:n_types),                   intent(inout) :: weight          ! Weight of each particle type (number of real 
                                                                                         !    particle for each simulated one)
      real(dp), dimension(0:n_types),                   intent(in)    :: rescale_factor  ! Factor by which the number of computational particles
                                                                                         !    will be divided when they exceed
                                                                                         !    the max allowed number
      real(dp), dimension(1:n_v_values),                intent(in)    :: v_values        ! Velocities for which electron collision frequencies
                                                                                         !    are given
      real(dp),                                         intent(in)    :: coll_freq       ! Total global collision frequency for electrons
                                                                                         !    (maximum value over all speeds)
      real(dp), dimension(1:n_totpro,1:n_v_values),     intent(in)    :: p_limits        ! Values used to select collision type for electrons
      integer,                                          intent(in)    :: n_limits        ! Number of probability limits for electron collisions
      integer,                                          intent(in)    :: limit1_exc      ! Index of prob limit for e- impact excitation
      integer,                                          intent(in)    :: limit1_diss     ! Index of prob limit for e- impact dissociation
      real(dp), dimension(1:n_types,1:n_v_values_ions), intent(in)    :: v_values_ions   ! Velocities for which ion collision frequencies
                                                                                         !    are given
      real(dp), dimension(1:n_types),                   intent(in)    :: coll_freq_ions  ! Total global collision frequency
                                                                                         !    for ions (maximum value over all speeds)
      real(dp), dimension(1:n_types, &
                         &1:n_types, &
                         &1:n_v_values_ions),           intent(in)    :: coll_f_tot_ions ! Total collision frequency (all processes) for 
                                                                                         !    collision of an ion on a specific neutral   
      real(dp), dimension(1:n_types, &
                         &1:n_types, &
                         &1:2, &
                         &1:n_v_values_ions),           intent(in)    :: p_limits_ions   ! Probability limits for ion scattering
                                                                                         !    (elastic/charge exchange/neutral)
                                                                                         !    for ech ion/neutral combination
      real(dp), dimension(1:n_types, &
                         &1:nmax_disspro, &
                         &1:n_v_values),                intent(in)    :: k_recomb        ! Collision rate. coeff for e-/ion recomb. / m**3*s**-1
      logical,                                          intent(in)    :: isactive_recomb ! Set if electron/ion recombiantion is active
      integer,                                          intent(in)    :: min_scattered   ! Minimum number of scatteredparticles needed
                                                                                         !    to apply the many particles method
      real(dp), dimension(1:n_types),                   intent(in)    :: e_loss          ! Mean loss of e- energy in elastic scattering
      real(dp), dimension(1:n_types,1:n_types),         intent(in)    :: e_loss_ions     ! Mean loss of ion energy in elastic scattering      
      real(dp), dimension(1:n_types),                   intent(in)    :: en_ion          ! Ionization energy of each neutral type
      real(dp), dimension(1:n_types,1:nmax_excpro),     intent(in)    :: en_exc          ! Excitation energy for each neutral type and 
                                                                                         !    each excitation process      
      real(dp), dimension(1:n_types,1:nmax_disspro),    intent(in)    :: en_diss         ! Dissociation energy for each neutral type and 
                                                                                         !    each dissociation process
      integer,  dimension(1:n_types),                   intent(in)    :: n_excpro        ! Number of available excitation processes 
                                                                                         !    for each neutral type      
      integer,  dimension(1:n_types),                   intent(in)    :: n_disspro       ! Number of available dissociation processes 
                                                                                         !    for each neutral type
      real(dp), dimension(1:n_types),                   intent(in)    :: mean_speed      ! Mean speed of neutral molecules / m*s**-1
      real(dp), dimension(1:n_types),                   intent(in)    :: se_coefficient  ! Secondary emission coefficient of each gas molecule
      real(dp),                                         intent(in)    :: distance        ! Distance between electrodes / m
      real(dp),                                         intent(in)    :: length          ! Lateral length of the electrodes / m
      real(dp),                                         intent(in)    :: Vmax            ! Electric bias (peak value for RF field) / V
      real(dp),                                         intent(in)    :: omega           ! Electric bias pulsation / rad*s**-1 
      real(dp),                                         intent(inout) :: phi             ! Electric bias phase / rad
      logical,                                          intent(in)    :: lateral_loss    ! Decide if e-/i+ are lost at the border, or recovered 
      real(dp), dimension(0:n_cells),                   intent(inout) :: rho             ! Charge density
      real(dp), dimension(0:n_cells),                   intent(inout) :: psi             ! Electric potential
      real(dp),                                         intent(in)    :: dt              ! Simulation timestep / s 
      real(dp),                                         intent(in)    :: duration        ! Required duration of the simulation time / s
      real(dp),                                         intent(inout) :: average_current ! Average electric current passed between the electrodes
                                                                                         !    (positive direction is the one following z axis)
      real(dp),                                         intent(inout) :: coll_null       ! Number of observed null collisions
      real(dp), dimension(1:n_types),                   intent(inout) :: coll_ela        ! Number of observed elastic collisions 
                                                                                         !    for each neutral type
      real(dp), dimension(1:n_types),                   intent(inout) :: coll_ion        ! Number of observed ionization collisions 
                                                                                         !    for each neutral gas type
      real(dp), dimension(1:n_types,1:nmax_excpro),     intent(inout) :: coll_exc        ! Number of observed excitation collisions
                                                                                         !    for each neutral type and excitation proc. 
      real(dp), dimension(1:n_types,1:nmax_disspro),    intent(inout) :: coll_dis        ! Number of observed dissociation collisions
                                                                                         !    for each neutral type and dissociation proc. 
      real(dp), dimension(1:n_types,1:nmax_disspro),    intent(inout) :: coll_rec        ! Number of observed electron-ion recombinations 
                                                                                         !    for each neutral gas type  
                                                                                         !    and dissociative recomb. process 
      integer,                                          intent(in)    :: n_particles     ! Maximum number of particles (same for each type)
      integer,                                          intent(in)    :: n_v_values      ! Number of tabulated values of electrons collision 
                                                                                         !    frequency
      integer,                                          intent(in)    :: n_v_values_ions ! Number of tabulated values of ion collision frequency
      integer,                                          intent(in)    :: n_totpro        ! Total maximum possible number of scattering processes 
                                                                                         !    for electrons: elastic + ionization + dissociation
      integer,                                          intent(in)    :: n_types         ! Number of neutral types (and number of ion types)
      integer,                                          intent(in)    :: nmax_excpro     ! Maximum number of excitation types      
      integer,                                          intent(in)    :: nmax_disspro    ! Maximum number of dissociation types
      integer,                                          intent(in)    :: n_cells         ! Number of cells used in the PIC scheme
      integer,                                          intent(in)    :: debug_level     ! Amount of output given for debugging purposes

      !  Local variables
      real(dp)                                     :: time               ! Simulation time
      real(dp), dimension(0:n_types,1:n_particles) :: dv_z               ! Increment of vz during time dt
      real(dp)                                     :: Vactual            ! Actual value of electric potential at the electrode
      real(dp)                                     :: area               ! Area of the electrodes
      real(dp)                                     :: volume             ! Volume embedded between the electrodes
      real(dp)                                     :: sinphi             ! sin(omega * t)
      real(dp)                                     :: d_phi              ! Phase increment in time dt
      real(dp), dimension(0:n_types,1:n_particles) :: v2                 ! Square of the velocity module      
      real(dp), dimension(1:n_types)               :: ion_density        ! Actual number density of each ion type
      real(dp)                                     :: p_collision        ! Collision probability for a single electron
      real(dp), dimension(1:n_types)               :: p_collision_ions   ! Collision probabilities for single ions
      real(dp), dimension(1:n_types, &
                         &1:nmax_disspro, &
                         &1:n_v_values)            :: kdt_recomb         ! Collision rate coeff. for e-/ion recombination,
                                                                         !    multiplied by timestep / m**3
      integer                                      :: n_ele_abs_0        ! Number of electrons absorbed by electrode at z=0
      integer                                      :: n_ele_abs_d        ! Number of electrons absorbed by electrode at z=d
      integer                                      :: n_ion_abs_0        ! Number of ions absorbed by electrode at z=0
      integer                                      :: n_ion_abs_d        ! Number of ions absorbed by electrode at z=d
      integer,  dimension(0:n_types,1:n_particles) :: indexes_removed    ! Indexes of removed particles of each type
      integer,  dimension(0:n_types)               :: n_part_removed     ! Number of removed particles of each type
      integer,  dimension(0:n_types)               :: n_part_added       ! Number of particles of each type (electrons and each ion type)
                                                                         !    added during an iteration due to ionization or SE processes
      integer,  dimension(0:n_types,1:n_particles) :: w_part_added       ! Number of computational particles to add for each real one
      real(dp), dimension(0:n_types,1:n_particles) :: x_part_added       ! Positions of the particles of each type to be added (x-component)
      real(dp), dimension(0:n_types,1:n_particles) :: y_part_added       ! Positions of the particles of each type to be added (y-component)
      real(dp), dimension(0:n_types,1:n_particles) :: z_part_added       ! Positions of the particles of each type to be added (z-component)
      real(dp), dimension(0:n_types,1:n_particles) :: vx_part_added      ! Velocities of the particles of each type to be added (x-component)
      real(dp), dimension(0:n_types,1:n_particles) :: vy_part_added      ! Velocities of the particles of each type to be added (y-component)
      real(dp), dimension(0:n_types,1:n_particles) :: vz_part_added      ! Velocities of the particles of each type to be added (z-component)      
      integer                                      :: i

      ! Variables used for performance checking
      integer                                      :: n_iterations
      real                                         :: tcpu_start, tcpu_end
      real,            dimension(0:6)              :: tcpu, tcpu_sum, tcpu_min, tcpu_max
      integer(kind=8)                              :: clock_rate, clock_start, clock_end
      integer(kind=8), dimension(0:6)              :: clock, clock_sum, clock_min, clock_max
      
      tcpu       = 0
      tcpu_sum   = 0
      tcpu_min   = 1
      tcpu_max   = 0
      clock      = 0
      clock_sum  = 0
      clock_min  = 1000000000
      clock_max  = 0      

      ! Initialize internal clock, since cloc_rate is kind=8, rate should be 1 us
      call system_clock(count_rate=clock_rate)    
      ! if (debug_level >0) print *, 'CLOCK RATE = ', clock_rate

      ! Initialize random number generator
      call init_random_seed()

      ! Calculate electrodes area and plasma volume
      area = length * length
      volume = area * distance

      ! Multiply the recombination rate coefficient by timestep
      kdt_recomb = k_recomb * dt

      ! Erase the counters of particles absorbed by the two electrodes
      ! the average current will be given by :
      ! i_average = (n_ele_abs_d - n_ion_abs_d - n_ele_abs_0 + n_ion_abs_0) * electron_charge / time
      n_ele_abs_0 = 0
      n_ele_abs_d = 0
      n_ion_abs_0 = 0
      n_ion_abs_d = 0
 
      ! Calculate the phase increment for the electric potential
      d_phi = omega * dt 

      ! Calculate the probability that a single electron has at least 1 scattering event
      ! of any type (including null ones) during time dt
      ! dt  -> time interval; coll_freq -> collision frequency (maximum value over all electron speeds)
      p_collision = 1.D0 - exp(- coll_freq * dt)
      
      ! Calculate, for each  ion type,  the probability that a single ion has
      ! at least 1 scattering event of any type (including null ones)
      ! during time dt
      ! dt  -> time interval; coll_freq -> collision frequencies for each ion type (maximum value over all ion speeds)
      do i = 1, n_types
         p_collision_ions(i) = 1.D0 - exp(- coll_freq_ions(i) * dt )
      enddo 

      ! Erase time counter (counts only the time elapsed during the call to simccp fuction, so it is erases at each new call)
      time = 0.0
      ! Reset the iteration counter, used for debug purposes only
      n_iterations = 0 

      if (debug_level > 0) call print_header_function

      ! +-----------+
      ! | MAIN LOOP |
      ! +-----------+
      
      do while (time < duration)  ! Repeat the main loop for the requested duration

         n_iterations = n_iterations + 1         

         if (debug_level > 1) call print_header_iteration

         ! +----------------------------------------------------+
         ! | Initialize counters of added and removed particles |
         ! +----------------------------------------------------+                
         
         !do i = 0, n_types      
         !   n_part_added(i)   = 0  ! Inizialize the counter of particles to add due to ionization events or secondary emission
         !   n_part_removed(i) = 0  ! Inizialize the counter of particles to remove since they left the plasma volume
         !enddo

         n_part_added   = 0
         n_part_removed = 0

         ! +--------------------------------------------------------+
         ! | Calculate electric potential and particle acceleration |
         ! +--------------------------------------------------------+
         
         if (debug_level > 0) then
            if (debug_level > 1) print *, "-> PIC"
            call cpu_time(tcpu_start)
            call system_clock(clock_start)
         endif

         ! Calculate the variable part of the electric potential at the electrode -> V(t) = V0 * sin[phi(t)];
         ! where  phi(t) = phi0 + omega*t 
         ! In case of DC bias the electric potential will be constant ->  V(t) = V0 -> sin(phi) = 1
         if (omega == 0) then  
            sinphi = 1.0_dp
         else
            sinphi = sin(phi)
         endif
         Vactual = Vmax * sinphi         
                  
         ! Calculate the actual velocity increment of e- and ions, based on the potential at the biased electrode
         ! and on the charge distribution         
         call calculate_velocity_increments(dv_z, rho, psi, z, isactive, weight, cm_ratio, &
                                            &dt, Vactual, distance, area, n_cells, n_types, n_particles)
         
         if (debug_level > 0) then
            call cpu_time(tcpu_end)
            tcpu(1) = tcpu_end - tcpu_start
            call system_clock(clock_end)
            clock(1) = clock_end - clock_start
         endif

         ! +-----------------------------------------------+
         ! | Move the particles using the leap-frog method |
         ! +-----------------------------------------------+
         
         if (debug_level > 0) then
            if (debug_level > 1) print *, "-> LEAP FROG"
            call cpu_time(tcpu_start)
            call system_clock(clock_start)            
         endif
         
         ! Move the charged particles         
         call leap_frog(x, y, z, vx, vy, vz, v, v2, &
                       &isactive, restart, &
                       &weight, & 
                       &n_types, n_particles, &
                       &dv_z, dt, &
                       &debug_level)
         
         if (debug_level > 0) then
            call cpu_time(tcpu_end)
            tcpu(2) = tcpu_end - tcpu_start
            call system_clock(clock_end)
            clock(2) = clock_end - clock_start
         endif

         ! +------------------------+
         ! |Check plasma boundaries |
         ! +------------------------+

         if (debug_level > 0) then
            if (debug_level > 1) print *, "-> CHECK BOUNDARIES"
            call cpu_time(tcpu_start)
            call system_clock(clock_start)            
         endif
         
         ! Check if some particles left the plasma volume, check also secondary emission processes         
         call check_boundaries(x, y, z, & 
                             &isactive, &
                             &weight, &
                             &x_part_added, y_part_added, z_part_added, w_part_added, n_part_added, &
                             &indexes_removed, n_part_removed, &
                             &n_types, n_particles, &
                             &distance, length, lateral_loss, &
                             &se_coefficient, &
                             &n_ele_abs_0, n_ele_abs_d, n_ion_abs_0, n_ion_abs_d, &
                             &debug_level)
         
         if (debug_level > 0) then
            call cpu_time(tcpu_end)
            tcpu(3) = tcpu_end - tcpu_start
            call system_clock(clock_end)
            clock(3) = clock_end - clock_start
         endif
        
         ! +------------------+
         ! | Remove particles |
         ! +------------------+

         if (debug_level > 0) then
            if (debug_level > 1) print *, "-> REMOVE PARTICLES"
            call cpu_time(tcpu_start)
            call system_clock(clock_start)            
         endif
         
         ! Remove particles that left the plasma volume         
         call remove_particles(x, y, z, vx, vy, vz, v, v2, &
                              &isactive, restart, &
                              &indexes_removed, n_part_removed, n_types, n_particles, &
                              &debug_level)
         
         if (debug_level > 0) then
            call cpu_time(tcpu_end)
            tcpu(4) = tcpu_end - tcpu_start
            call system_clock(clock_end)
            clock(4) = clock_end - clock_start
         endif

         ! +----------------------+
         ! | Collision Management |
         ! +----------------------+        

         if (debug_level > 0) then
            if (debug_level > 1) print *, "-> COLLISIONS CHECK"
            call cpu_time(tcpu_start)
            call system_clock(clock_start)            
         endif
         
         ! Calculate ion number density for each ion type, which is required to calculate recombination probability
         do i = 1, n_types
            ion_density(i) = real( count( isactive(i,:) ) ) * weight(i) / volume
         enddo
         ! Check which particles experienced collisions and act properly depending on collision type
         call check_collisions(x, y, z, vx, vy, vz, v, v2, &
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
         
         if (debug_level > 0) then
            call cpu_time(tcpu_end)
            tcpu(5) = tcpu_end - tcpu_start
            call system_clock(clock_end)
            clock(5) = clock_end - clock_start
         endif

         ! +---------------+
         ! | Add particles |
         ! +---------------+        
         
         if (debug_level > 0) then
            if (debug_level > 1) print *, "-> ADD PARTICLES"
            call cpu_time(tcpu_start)
            call system_clock(clock_start)            
         endif
         
         ! Add electrons and ions generated by ionization processes
         call add_particles(x, y, z, vx, vy, vz, v, v2, &
                           &isactive, restart, &
                           &weight, &
                           &rescale_factor, &
                           & x_part_added,  y_part_added,  z_part_added, &
                           &vx_part_added, vy_part_added, vz_part_added, &
                           &w_part_added, n_part_added, &
                           &n_types, n_particles, &
                           &debug_level)
         if (debug_level > 0) then
            call cpu_time(tcpu_end)
            tcpu(6) = tcpu_end - tcpu_start
            call system_clock(clock_end)
            clock(6) = clock_end - clock_start
         endif

         ! Calculate the new time and new phase -> phi(t+dt) = phi0 + omega * (t+dt) = phi(t) + omega * dt = phi(t) + dphi
         time = time + dt
         phi  = phi  + d_phi
         
         if (debug_level > 0) then
            tcpu(0)  = sum(tcpu(1:6))  ! Total CPU time in this iteration
            clock(0) = sum(clock(1:6)) ! Total clock time in this iteration
            ! Update the counter of CPU time, min and max, among all iterations
            do i = 0, 6
               tcpu_sum(i)  = tcpu_sum(i)  + tcpu(i)
               if (tcpu(i) < tcpu_min(i)) tcpu_min(i) = tcpu(i)
               if (tcpu(i) > tcpu_max(i)) tcpu_max(i) = tcpu(i)
               clock_sum(i) = clock_sum(i)  + clock(i)
               if (clock(i) < clock_min(i)) clock_min(i) = clock(i)
               if (clock(i) > clock_max(i)) clock_max(i) = clock(i)
            enddo
         endif   
         
      enddo ! (time < duration)

      ! Calculate the average electric current that was flowing between the electrodes during the simulated time
      average_current = (n_ele_abs_d - n_ion_abs_d - n_ele_abs_0 + n_ion_abs_0) * E_CHARGE / time

      if (debug_level > 1) call collision_report
      if (debug_level > 0) then
         call cpu_time_report         
         call print_end_signature      
      endif

    contains

      subroutine print_header_function

         print *, ""
         print *, "*** BEGIN OF FORTRAN FUNCTION SIMCCP ***"
         print *, ""         
         print *, "Electric field pulsation [rad/s]        = ", omega
         if (omega > 0) then
            print *, "Electric field phase [rad]              = ", phi
            print *, "Max electric potential [V]              = ", Vmax
            print *, "Start electric potential [V/m]          = ", Vmax*sin(phi)
         else
            print *, "Electric field intensity (static) [V/m] = ", Vmax
         endif
         print *, "Time step dt [s]                        = ", dt
         print *, "Requested duration [s]                  = ", duration
!         print *, "Maximum velocity increase in dt [m/s] "
         do i = 0, n_types
            print *, "Particle type, cm ratio          ", i, cm_ratio(i)
         enddo
!        print *, "isactive", isactive
!        print *, "restart ", restart
         if (debug_level > 1) call pause

      end subroutine print_header_function


      subroutine print_header_iteration
         print *, " "
         print *, "ITERATION NUMBER                    = ", n_iterations
         print *, "Time, requested duration [s]        = ", time, duration
         print *, "Phase [rad]                         = ", phi
         !print *, " "         
         if (debug_level > 1) call pause
      end subroutine print_header_iteration

      subroutine collision_report

         integer :: i,j 
         
         print *, " "
         print *, "ELECTRON COLLISION REPORT"
         print *, "Elastic:      ", coll_ela
         print *, "Ionization:   ", coll_ion
         print *, "Excitation:   ", sum(coll_exc)
         do i = 1, n_types
            do j = 1, nmax_excpro
               print *, i, j, coll_exc(i,j)
            enddo
         enddo         
         print *, "Dissociation: ", sum(coll_dis)
         do i = 1, n_types
            do j = 1, nmax_disspro
               print *, i, j, coll_dis(i,j)
            enddo
         enddo
         print *, "Null:         ", coll_null
         if (debug_level > 1) call pause         
         
      end subroutine collision_report
      
      subroutine print_end_signature

         if (debug_level > 1) then
            print *, " "
            print *, "CURRENT BALANCE"
            print *, "Electrons absobed at z=0: ", n_ele_abs_0
            print *, "Electrons absobed at z=d: ", n_ele_abs_d
            print *, "Difference:               ", n_ele_abs_d - n_ele_abs_0
            print *, "Ions absobed at z=0:      ", n_ion_abs_0
            print *, "Ions absobed at z=d:      ", n_ion_abs_d
            print *, "Difference:               ", n_ion_abs_d - n_ion_abs_0
            print *, "Net electric current:     ", average_current, " A"
         endif
         
         print *, " "
         print *, "Time                                    = ", time
         print *, "Phase                                   = ", phi
!        print *, "isactive", isactive
!        print *, "restart ", restart
         print *, ""
         print *, "*** END OF FORTRAN FUNCTION SIMCCP ***"
         print *, ""         
         if (debug_level > 1) call pause

      end subroutine print_end_signature

      
      subroutine cpu_time_report
       
         ! Local variables
         real, dimension(0:6) ::  tcpu_mean, clock_mean

         tcpu_mean =  tcpu_sum / n_iterations
         clock_mean = real(clock_sum) / real(n_iterations)

         print *, ''         
         print *, 'CPU TIME REPORT'
         print *, 'Iterations       ->', n_iterations
!         print *, ''
         print *, 'CPU TIME (s)       ', char(9), char(9),'Mean', char(9), char(9), 'delta', char(9), char(9),'%'         
         print *, 'All steps        ->', tcpu_mean(0), tcpu_max(0)-tcpu_min(0)  
         print *, 'Particle in cell ->', tcpu_mean(1), tcpu_max(1)-tcpu_min(1), tcpu_sum(1)/tcpu_sum(0)*100
         print *, 'Leap frog        ->', tcpu_mean(2), tcpu_max(2)-tcpu_min(2), tcpu_sum(2)/tcpu_sum(0)*100
         print *, 'Check boundaries ->', tcpu_mean(3), tcpu_max(3)-tcpu_min(3), tcpu_sum(3)/tcpu_sum(0)*100
         print *, 'Remove particles ->', tcpu_mean(4), tcpu_max(4)-tcpu_min(4), tcpu_sum(4)/tcpu_sum(0)*100
         print *, 'Check collisions ->', tcpu_mean(5), tcpu_max(5)-tcpu_min(5), tcpu_sum(5)/tcpu_sum(0)*100
         print *, 'Add particles    ->', tcpu_mean(6), tcpu_max(6)-tcpu_min(6), tcpu_sum(6)/tcpu_sum(0)*100
!         print *, ''
         print *, 'CLOCK TIME (s)     ', char(9), char(9),'Mean', char(9), char(9), 'delta', char(9), char(9),'%' 
         print *, 'All steps        ->', clock_mean(0)*1.0E-9, (clock_max(0)-clock_min(0))*1.0E-9  
         print *, 'Particle in cell ->', clock_mean(1)*1.0E-9, (clock_max(1)-clock_min(1))*1.0E-9, &
                                         real(clock_sum(1))/real(clock_sum(0))*100.0
         print *, 'Leap frog        ->', clock_mean(2)*1.0E-9, (clock_max(2)-clock_min(2))*1.0E-9, &
                                         real(clock_sum(2))/real(clock_sum(0))*100.0
         print *, 'Check boundaries ->', clock_mean(3)*1.0E-9, (clock_max(3)-clock_min(3))*1.0E-9, &
                                         real(clock_sum(3))/real(clock_sum(0))*100.0
         print *, 'Remove particles ->', clock_mean(4)*1.0E-9, (clock_max(4)-clock_min(4))*1.0E-9, &
                                         real(clock_sum(4))/real(clock_sum(0))*100.0
         print *, 'Check collisions ->', clock_mean(5)*1.0E-9, (clock_max(5)-clock_min(5))*1.0E-9, &
                                         real(clock_sum(5))/real(clock_sum(0))*100.0
         print *, 'Add particles    ->', clock_mean(6)*1.0E-9, (clock_max(6)-clock_min(6))*1.0E-9, &
                                         real(clock_sum(6))/real(clock_sum(0))*100.0
         if (debug_level > 1) call pause
         
      end subroutine cpu_time_report
      
   end subroutine simccp

end module f_main

