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


module f_particle_mover

   use f_precision

   use f_io

   use f_random

   use f_arrays

!   use f_discretization

   use f_constants

   implicit none

contains

   subroutine leap_frog(x, y, z, vx, vy, vz, v, v2, &
                       &isactive, restart, &
                       &weight, &
                       &n_types, n_particles, &
                       &dv_z, dt,&
                       &debug_level)

      !f2py intent(hide)      :: n_types, n_particles
      !f2py intent(in)        :: sinphi, dt, dv_z, debug_level
      !f2py intent(inout)     :: x, y, z, v, vx, vy, vz, v2
      !f2py intent(inout)     :: isactive, restart

      ! Parameters
      real(dp), dimension(0:n_types, 1:n_particles), intent(inout) :: x, y, z, vx, vy, vz, v, v2    ! Positions and velocities of the particles
      logical,  dimension(0:n_types, 1:n_particles), intent(inout) :: isactive, restart             ! Select which particles are active,
                                                                                                    !    and which must restart leap-frog scheme
                                                                                                    !    at the next iteration
      real(dp), dimension(0:n_types),                intent(in)    :: weight                        ! Weight of each particle type 
                                                                                                    !    (number of real particles represented
                                                                                                    !    by each simulated one)
      real(dp),                                      intent(in)    :: dt
      real(dp), dimension(0:n_types, 1:n_particles), intent(in)    :: dv_z
      integer,                                       intent(in)    :: n_types, n_particles
      integer,                                       intent(in)    :: debug_level

      ! Local variables
      integer  :: particle_type, particle
 
      do particle_type = 0, n_types
!         dv_z = dvmax_z(particle_type) * sinphi
         do particle = 1, n_particles
            if (isactive(particle_type, particle)) then 
               call leap_frog_particle
            endif
         enddo
      enddo

   contains
      
      subroutine leap_frog_particle

         ! Local variables
         logical  :: restarted_particle

         ! Calculates new particle velocity v(t-dt/2)
         ! If restart is true, then stored values of velocity and position are referred to the same time
         ! so it is needed to shift velocity of dt/2 in order to apply leap-frog scheme
         if (restart(particle_type, particle)) then                                         ! restart = .true. -->  vz is vz(0)
            vz(particle_type, particle) = vz(particle_type, particle) + &
                                         &dv_z(particle_type, particle) / 2     ! vz becomes vz(dt/2) = vz(0) + az(0) * dt/2
            restart(particle_type, particle) = .false.
            restarted_particle = .true.
         else                                                                      ! restart = .false. --> vz is vz(t-3*dt/2)
            vz(particle_type, particle) = vz(particle_type, particle) + &
                                          &dv_z(particle_type, particle)         ! vz becomes vz(t-dt/2) = vz(t-3*dt/2) + az(t-dt) * dt
            restarted_particle = .false.
         endif ! (restart(particle_type, particle))

         ! Calculate the modulus of velocity (v) and store its squared value (to be used later)
         v2(particle_type, particle) = vx(particle_type,particle) * vx(particle_type,particle) + &
                             &vy(particle_type,particle) * vy(particle_type,particle) + &
                             &vz(particle_type,particle) * vz(particle_type,particle)
         v(particle_type, particle)  = sqrt(v2(particle_type, particle))

         ! Calculates new position x(t), y(t), z(t)
         x(particle_type,particle) = x(particle_type,particle) + vx(particle_type,particle) * dt  ! x(i*dt) = x((i-1)*dt) + vx((i-1/2)*dt) * dt
         y(particle_type,particle) = y(particle_type,particle) + vy(particle_type,particle) * dt  ! y(i*dt) = y((i-1)*dt) + vy((i-1/2)*dt) * dt
         z(particle_type,particle) = z(particle_type,particle) + vz(particle_type,particle) * dt  ! z(i*dt) = z((i-1)*dt) + vz((i-1/2)*dt) * dt

         if (debug_level > 2) then 
            print *, ""
            if (restarted_particle) print *, " -> Restarted !"
            print *, "Particle type, max type, weight     = ", particle_type, n_types, weight(particle_type)
            print *, "Particle number, max particles      = ", particle, n_particles
            print *, "dv_z [m/s]                          = ", dv_z(particle_type, particle)
            print *, "vx, vy, vz, v [m/s]                 = ", vx(particle_type,particle), &
                                                              &vy(particle_type,particle), &
                                                              &vz(particle_type,particle), &
                                                              &v(particle_type,particle)
            print *, "x,y,z [m]                           = ", x(particle_type,particle), &
                                                              &y(particle_type,particle), &
                                                              &z(particle_type,particle)
         endif

      end subroutine leap_frog_particle

   end subroutine leap_frog
 

   subroutine check_boundaries(x, y, z, & 
                             &isactive, &
                             &weight, &
                             &x_part_added, y_part_added, z_part_added, w_part_added, n_part_added, &
                             &indexes_removed, n_removed, &
                             &n_types, n_particles, &
                             &distance, length, lateral_loss, &
                             &se_coefficient, &
                             &n_ele_abs_0, n_ele_abs_d, n_ion_abs_0, n_ion_abs_d, &
                             &debug_level)

      ! Parameters
      real(dp), dimension(0:n_types, 1:n_particles), intent(inout) :: x, y, z                ! Positions and velocities
                                                                                             !    of the particles
      logical,  dimension(0:n_types, 1:n_particles), intent(inout) :: isactive               ! Select which particles
                                                                                             !    are active, and which ones
                                                                                             !    must restart leap-frog
                                                                                             !    scheme at the next iteration
      real(dp), dimension(0:n_types),                intent(inout) :: weight                 ! Weight of each particle type 
                                                                                             !    (number of real particles
                                                                                             !     represented by each !
                                                                                             !    simulated one)
      real(dp), dimension(0:n_types, 1:n_particles), intent(inout) :: x_part_added, &        ! Positions of particles that need to be added
                                                                     &y_part_added, &
                                                                     &z_part_added         
      integer,  dimension(0:n_types, 1:n_particles), intent(inout) :: w_part_added           ! Number of comp. particles to add for each real one
      integer,  dimension(0:n_types),                intent(inout) :: n_part_added           ! Number of particles that need to be added
      integer,  dimension(0:n_types, 1:n_particles), intent(inout) :: indexes_removed        ! Indexes of removed particles of each type
      integer,  dimension(0:n_types),                intent(inout) :: n_removed              ! Number of removed particles of each type
      integer,                                       intent(inout) :: n_ele_abs_0            ! Number of electrons absorbed by electrode at z=0
      integer,                                       intent(inout) :: n_ele_abs_d            ! Number of electrons absorbed by electrode at z=d
      integer,                                       intent(inout) :: n_ion_abs_0            ! Number of ions absorbed by electrode at z=0
      integer,                                       intent(inout) :: n_ion_abs_d            ! Number of ions absorbed by electrode at z=d 
      real(dp),                                      intent(in)    :: distance               ! Distance between electrodes / m
      real(dp),                                      intent(in)    :: length                 ! Lateral length of the electrodes / m
      logical,                                       intent(in)    :: lateral_loss           ! Decide if e-/i+ are lost at the border,
                                                                                             !    or recovered       
      real(dp), dimension(1:n_types),                intent(in)    :: se_coefficient
      integer,                                       intent(in)    :: n_types, n_particles   
      integer,                                       intent(in)    :: debug_level


      ! Local variables
      integer :: particle_type, particle

      do particle_type = 0, n_types
         do particle = 1, n_particles
            if (isactive(particle_type,particle)) then 
               call check_boundaries_particle
            endif
         enddo
      enddo

   contains

      subroutine check_boundaries_particle

         ! If the particle has reached one of the electrodes, remove it and add to the electric current
         ! if it is an ion reaching an electrode, activate secondary emission process
         if (z(particle_type, particle) > distance) then
            if (particle_type == 0) then
                n_ele_abs_d = n_ele_abs_d + int(weight(particle_type))
            else
                n_ion_abs_d = n_ion_abs_d + int(weight(particle_type))
            endif
            if (particle_type > 0) call secondary_emission(particle_type, x(particle_type,particle), y(particle_type,particle), &
                                                          &distance)
            n_removed(particle_type) = n_removed(particle_type) + 1
            indexes_removed(particle_type, n_removed(particle_type)) = particle 
         elseif (z(particle_type, particle) < 0) then
            if (particle_type == 0) then
               n_ele_abs_0 = n_ele_abs_0 + int(weight(particle_type))
            else
               n_ion_abs_0 = n_ion_abs_0 + int(weight(particle_type))
            endif
            if (particle_type > 0)  call secondary_emission(particle_type, x(particle_type,particle), y(particle_type,particle), &
                                                           &0.0_dp)
            n_removed(particle_type) = n_removed(particle_type) + 1
            indexes_removed(particle_type, n_removed(particle_type)) = particle
         endif ! (z(particle_type, particle) > distance) 

        ! What to do if the particle has reached the edge of the electrode
         if (lateral_loss) then
            ! If the particle has reached the edge of the electrode, remove it
            if ( (x(particle_type,particle) > length) .or. &
                &(y(particle_type,particle) > length) .or. &
                &(x(particle_type,particle) < 0     ) .or. &
                &(y(particle_type,particle) < 0     )       ) then 
               n_removed(particle_type) = n_removed(particle_type) + 1
               indexes_removed(particle_type, n_removed(particle_type)) = particle
            endif
         else   
            ! If the particle has reached the edge of the electrode, re-enters at the opposite border (infinite length electrodes)
            if (x(particle_type,particle) > length)     then
               x(particle_type,particle) = 0.0
            elseif (x(particle_type,particle) < 0     ) then
               x(particle_type,particle) = length
            endif
            if     (y(particle_type,particle) > length) then
               y(particle_type,particle) = 0.0
            elseif (y(particle_type,particle) < 0     ) then
               y(particle_type,particle) = length
            endif   
         endif

      end subroutine check_boundaries_particle


      subroutine secondary_emission(ion_type, x, y, z)

         ! Test if an ion approaching an electrode will emit a secondary electron

         ! Parameters
         integer,  intent(in) :: ion_type ! Type of the ion approaching the electrode 
         real(dp), intent(in) :: x,y,z    ! Position at which electrons will be emitted

         ! Local variables
         real(sp) :: rand1!, rand2
         integer  :: i, n_secondaries
         real(dp) :: f, nw_secondaries


         ! Calculate the number of computational electrons to be emitted by each computational ion
         nw_secondaries = se_coefficient(ion_type) * weight(ion_type) / weight(0)

         i = int(nw_secondaries)  ! -> number of comp. electrons to be surely emitted for each comp. ion (may be zero)
         f = nw_secondaries - i   ! -> probability to emit a further comp. electron for each comp. ion

         if (f > 0) then
            call random_number(rand1)
            if (rand1 < f) then
               i = i + 1
            endif
         endif 

         n_secondaries = i
          
         if (debug_level > 0) then
            print *, ""
            print *, "SECONDARY EMISSION"
            print *, "ion                = ", ion_type 
            print *, "SE                 = ", se_coefficient(ion_type)
            if (f > 0 ) print *, "rand1              = ", rand1
            print *, "weight ratio       = ", weight(ion_type) / weight(0)
            print *, "nw_secondaries     = ", nw_secondaries
            print *, "n_secondaries      = ", n_secondaries
            if (debug_level > 1) call pause
            call pause
         endif

         ! Check that the number of SE does not exceed the capacity of arrays used to store the positions at which they must be generated
         !if (n_secondaries > n_particles) n_secondaries = n_particles  

         ! For the electric current, emission of secondaries electrons is equivalent to ions absorpion
         if (z >= distance) then
            n_ion_abs_d = n_ion_abs_d + n_secondaries * int(weight(0))
         else
            n_ion_abs_0 = n_ion_abs_0 + n_secondaries * int(weight(0))
         endif

         if (n_secondaries > 0) then
            n_part_added(0)                 = n_part_added(0) + 1  ! will add some electrons from this ion
            w_part_added(0,n_part_added(0)) = n_secondaries        ! save the number of computational electrons to be added
            x_part_added(0,n_part_added(0)) = x                    ! save the position of emission secondary electron emitted (x-component) 
            y_part_added(0,n_part_added(0)) = y                    ! save the position of emission secondary electron emitted (y-component) 
            z_part_added(0,n_part_added(0)) = z                    ! save the position of emission secondary electron emitted (z-component) 
         endif

      end subroutine secondary_emission

   end subroutine check_boundaries

end module f_particle_mover

