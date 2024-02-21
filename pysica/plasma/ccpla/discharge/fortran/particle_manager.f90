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

module f_particle_manager

   use f_precision

   use f_io

   use f_random

   use f_arrays

   use f_constants

   implicit none

contains

   subroutine add_particles(x, y, z, vx, vy, vz, v, v2, &
                           &isactive, restart, &
                           &weight, &
                           &rescale_factor, &
                           & x_part_added,  y_part_added,  z_part_added, &
                           &vx_part_added, vy_part_added, vz_part_added, &                           
                           &w_part_added, n_part_added, &
                           &n_types, n_particles, &
                           &debug_level)

      ! Parameters
      real(dp), dimension(0:n_types, 1:n_particles), intent(inout) :: x, y, z              ! Positions of particles
      real(dp), dimension(0:n_types, 1:n_particles), intent(inout) :: vx, vy, vz, v, v2    ! Velocities of particles
      logical,  dimension(0:n_types, 1:n_particles), intent(inout) :: isactive, restart    ! Select which particles are active, and which ones
                                                                                           !    must restart leap-frog scheme
                                                                                           !    at the next iteration
      real(dp), dimension(0:n_types),                intent(inout) :: weight               ! Weight of each particle type 
      real(dp), dimension(0:n_types),                intent(in)    :: rescale_factor       ! Factor by which the number of particles
                                                                                           !    will be divided
      real(dp), dimension(0:n_types, 1:n_particles), intent(inout) :: x_part_added, &
                                                                     &y_part_added, &
                                                                     &z_part_added         ! Positions of particles that need to be added
      real(dp), dimension(0:n_types, 1:n_particles), intent(inout) :: vx_part_added, &
                                                                     &vy_part_added, &
                                                                     &vz_part_added        ! Velocities of particles that need to be added
      integer,  dimension(0:n_types),                intent(inout) :: n_part_added         ! Number of events (e.g. ionizations) leading to
                                                                                           !    particle addition      
      integer,  dimension(0:n_types, 1:n_particles), intent(inout) :: w_part_added         ! Number of comp. particles to add for each event

      integer,                                       intent(in)    :: n_types, n_particles
      integer,                                       intent(in)    :: debug_level

      ! Local variables
      integer  :: n, m, t, n_active, n_min, n_add, ntot_part_added
      real(dp) :: x0, y0, z0, vx0, vy0, vz0, wratio, rand
      logical  :: rescaled
   
      do t = 0, n_types
         ! Calculate the active number of particles of this type
         ! n_active = n_state_elements(isactive, n_types+1, n_particles, t+1, .true.)
         n_active = count(isactive(t,:))
!         print *, 'Type, n_active:',t, n_active,n_state_elements(isactive, n_types+1, n_particles, t+1, .true.)
         ! If there have been events that added particles of this type, manage particle addition
         if (n_part_added(t) > 0) then
!            print *, 'Type, events:  ',t, n_part_added(t)
            ! Calculate the total number of computational particles to be added for this type
!            ntot_part_added = 0
!            do n = 1, n_part_added(t)
!               ntot_part_added = ntot_part_added + w_part_added(t, n)
!            enddo
            ntot_part_added = sum(w_part_added(t,:))
!            print *, 'Type, added particles:',t, ntot_part_added
            ! Check that, adding the required number of particles, the max allowed number is not exceeded
            ! otherwise rescale the number of particles and their weight
            rescaled = .false.
            if ( ((n_active + ntot_part_added) > n_particles) .and. (rescale_factor(t) > 1.0) ) then
               call rescale_particles(t, rescale_factor(t), n_active)
               rescaled = .true.                 
            endif
            ! Add the required number of computational particles of type t
            do n = 1, n_part_added(t)
               ! If the ensable has been rescaled, the number of comp. particles to add for each event must be rescaled as well
               if (rescaled) then
                  ! The number of particles to add must be reduced of the same factor used to increase their weight
                  wratio = w_part_added(t, n) / rescale_factor(t)
                  ! If the result is less than one particle to add, it becomes the probabilty that a particle is added
                  if (wratio >= 1) then
                     w_part_added(t, n) = nint(wratio)
                  else
                     call random_number(rand)
                     if (rand < wratio) then
                        w_part_added(t, n) = 1
                     else
                        w_part_added(t, n) = 0
                     endif ! (rand < wratio)   
                  endif ! (wratio >= 1)
               endif ! (rescaled)
               ! Add the proper number of computational particles
               do m = 1, w_part_added(t, n)
                  call add_particle(t, &
                                  &  x_part_added(t, n) , y_part_added(t, n),  z_part_added(t, n), &
                                  & vx_part_added(t, n), vy_part_added(t, n), vz_part_added(t, n))
               enddo !m
            enddo !n
            if (debug_level > 1) print *,'ADDED PARTICLES: TYPE', t, 'n. ', n_part_added(t)
            
!            print *, 'ADDING PARTICLES'
!            print *, 'type            =', t
!            print *, 'n_active        =', n_active
!            print *, 'n_part_added    =', n_part_added(t)
!            print *, 'ntot_part_added =', ntot_part_added
!            call pause
  
            n_part_added(t) = 0 ! Erase the number of particles of type t to add
         endif 

         ! If the number of computational particles is too low, then
         ! repopulate the particles ensamble and reduce their weight
         if ( (n_particles >= 100 * rescale_factor(t)) .and. (weight(t) >= rescale_factor(t)) ) then
            n_min = int(n_particles / 100)
            if (n_active < n_min) then
               ! Location and velocity of the new particles will be equal to the average of the existing ones
               x0  = array_average_masked(x,  isactive, n_types+1, n_particles, t+1)
               y0  = array_average_masked(y,  isactive, n_types+1, n_particles, t+1)
               z0  = array_average_masked(z,  isactive, n_types+1, n_particles, t+1)
               vx0 = array_average_masked(vx, isactive, n_types+1, n_particles, t+1)
               vy0 = array_average_masked(vy, isactive, n_types+1, n_particles, t+1)
               vz0 = array_average_masked(vz, isactive, n_types+1, n_particles, t+1)
               n_add = nint(n_active * rescale_factor(t)) - n_active
               do n = 1, n_add ! Add the required number of particles of type t
                  call add_particle(t, x0, y0, z0, vx0, vy0, vz0)
               enddo
               weight(t) = weight(t) / rescale_factor(t)
            endif ! n_active < n_min
         endif ! n_particles >= 1000        
      enddo ! t

   contains

      subroutine add_particle(t, x0, y0, z0, vx0, vy0, vz0)
      
         ! Add a particle of type t at the position (x0,y0,z0), with velocity (vx0, vy0, vz0)

         ! Parameters
         integer,  intent(in) :: t
         real(dp), intent(in) :: x0, y0, z0, vx0, vy0, vz0

         ! Local variables
         integer :: m

         do m = 1, n_particles
            if (isactive(t, m).eqv..false.) exit
         enddo
         if (m <= n_particles) then
            isactive(t, m) = .true.
            restart(t, m)  = .true.
            x(t, m)        = x0
            y(t, m)        = y0
            z(t, m)        = z0
            vx(t, m)       = vx0
            vy(t, m)       = vy0
            vz(t, m)       = vz0
            v2(t,m)        = vx(t, m) * vx(t, m) + vy(t, m) * vy(t, m) + vz(t, m) * vz(t, m)
            v(t, m)        = sqrt(v2(t,m))     
          endif

      end subroutine add_particle


      subroutine rescale_particles(t, f, n_active)
 
          ! Rescale the number of particles of type t by a factor f

          ! Parameters
          integer,  intent(in) :: t, n_active
          real(dp), intent(in) :: f       

          ! Local variables
          integer              :: n_remove
          integer, allocatable :: tofalse_indexes_relative(:), tofalse_indexes_absolute(:)
      
          ! Calculate the number of particles to remove
          n_remove = n_active - nint(abs(n_active / f))

          if (debug_level > 1) then
             print *
             print *, "RESCALE PARTICLES of type", t
             print *, "Scaling factor           ", f
             print *, "Number of particles      ", n_active
             print *, "Number to be removed     ", n_remove
          endif

          allocate(tofalse_indexes_relative(1:n_remove))
          allocate(tofalse_indexes_absolute(1:n_remove))
 
          ! Determine randomly the particles to remove
          ! will return a set of indexes, identifying the particles that must be removed
          ! each particle is represented by a .true. item in a row of the isactive array
          call random_integers_array(tofalse_indexes_relative, n_remove, 1, n_active, .true.)

          ! Remove the particles
          call find_absolute_indexes(isactive, n_types+1, n_particles, t+1, tofalse_indexes_relative, tofalse_indexes_absolute, &
                                    &n_remove, .true.)
          call set_elements(isactive, n_types+1, n_particles, t+1, tofalse_indexes_absolute, n_remove, .false.)
          call set_elements(restart,  n_types+1, n_particles, t+1, tofalse_indexes_absolute, n_remove, .true.)

          deallocate(tofalse_indexes_relative, tofalse_indexes_absolute)

          ! Increase proportionally the weight of this type of particle
          weight(t) = abs(weight(t) * f)

          if (debug_level > 1) then
             print *, "New number of particles", n_state_elements(isactive, n_types+1, n_particles, t+1, .true.)
             call pause       
          endif

      end subroutine rescale_particles

   end subroutine add_particles


   subroutine remove_particles(x, y, z, vx, vy, vz, v, v2, &
                              &isactive, restart, &
                              &indexes_removed, n_part_removed, &
                              &n_types, n_particles, &
                              &debug_level)

      ! Parameters
      real(dp), dimension(0:n_types,1:n_particles), intent(inout) :: x, y, z            ! Positions of particles
      real(dp), dimension(0:n_types,1:n_particles), intent(inout) :: vx, vy, vz, v, v2  ! Velocities of particles
      logical,  dimension(0:n_types,1:n_particles), intent(inout) :: isactive, restart  ! Select which particles are active, and which ones
                                                                                        !    must restart leap-frog scheme at the next iteration
      integer,  dimension(0:n_types,1:n_particles), intent(in)    :: indexes_removed    ! Indexes of removed particles of each type
      integer,  dimension(0:n_types),               intent(in)    :: n_part_removed     ! Number of removed particles of each type
      integer,                                      intent(in)    :: n_types, n_particles
      integer,                                      intent(in)    :: debug_level

      ! Local variables
      integer   :: t, i

      do t = 0, n_types
         if (n_part_removed(t) > 0) then
            do i = 1, n_part_removed(t)
               call remove_particle(t, indexes_removed(t,i))
            enddo
         endif
      enddo

   contains

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

   end subroutine remove_particles

end module f_particle_manager

