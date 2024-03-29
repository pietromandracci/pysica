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

module f_scattering

   use f_precision

   use f_constants

   use f_io

   implicit none

contains


   subroutine isotropic_scattering(vers_x, vers_y, vers_z, debug)

      ! Calculates the scattering angles for a completely isotropic scattering process

      !f2py intent(out) :: vers_x, vers_y, vers_z
      !f2py intent(in)  :: debug

      ! Parameters
      real(dp), intent(out) :: vers_x   ! Projection on x axis of velocity vector after scattering (sinthetacosphi)
      real(dp), intent(out) :: vers_y   ! Projection on y axis of velocity vector after scattering (sinthetasinphi)
      real(dp), intent(out) :: vers_z   ! Projection on z axis of velocity vector after scattering (costheta)
      logical,  intent(in)  :: debug    ! Set debug mode 

      ! Local variables
      real(dp) :: sintheta, phi !theta
      real(sp) :: r1,r2             ! Random numbers in range [0,1]

      if (debug) then
         print *,""
         print *,"**SUBROUTINE ISOTROPIC SCATTERING**"
      endif

      ! Calculate scattering angles and projections of scattered velocity versor
      call random_number(r1)
 
      vers_z = 2 * r1 - 1                    ! cos(theta) uniformly distribuited in range [-1,1]
      sintheta = sqrt(1 - vers_z * vers_z)

      call random_number(r2)                 ! phi uniformly distribuited in range [0,2*PI]
      phi = 2 * PI * r2
      
      vers_x = sintheta * cos(phi)
      vers_y = sintheta * sin(phi)

      if (debug) then
         print *,"IN"
         print *,"r1             = ", r1
         print *,"r2             = ", r2
         print *,"sinthetasinphi = ", vers_x
         print *,"sinthetacosphi = ", vers_y
         print *,"costheta       = ", vers_z
         call pause
      endif
  
   end subroutine isotropic_scattering
   

   subroutine anisotropic_scattering_electron(vers_x, vers_y, vers_z, xi, isatom, debug)

      ! Calculates the scattering angles for elastic collision of an electron with an atom or a nonpolar molecule
      ! the used formula is correct for elastic scattering with atoms (screened Coulomb potential) and is taken from the paper
      ! PHYSICAL REVIEW E 65 (2002) 037402

      !f2py intent(out) :: vers_x, vers_y, vers_z
      !f2py intent(in)  :: xi, atom, debug

      ! Parameters
      real(dp), intent(out) :: vers_x   ! Projection on x axis of velocity vector after scattering (sinthetacosphi)
      real(dp), intent(out) :: vers_y   ! Projection on y axis of velocity vector after scattering (sinthetasinphi)
      real(dp), intent(out) :: vers_z   ! Projection on z axis of velocity vector after scattering (costheta)
      real(dp), intent(in)  :: xi       ! Parameter the meaning of which depends on the type of target particle
                                        !   atom     -> squared electron speed
                                        !   molecule -> value of the xi function for the current electron energy 
      logical,  intent(in)  :: isatom   ! if .true. the scattering is calculated for an atom, if .false. for a molecule
      logical,  intent(in)  :: debug    ! Set debug mode 

      ! Local variables
      real(dp) :: e                     ! electron energy espressed in atomic energy units (au = 27.21 eV)
      real(dp) :: sintheta, phi
      real(sp) :: r1,r2 

      if (debug) then
         print *,""
         print *,"**SUBROUTINE ELASTIC SCATTERING NEUTRAL**"
      endif

      ! Calculate scattering angles
      call random_number(r1)

      if (isatom) then
         e = E_V2TOENERGY_AU * xi                     ! Calculate energy in atomic units
         vers_z   = 1 - 2 * r1 / (1 + 8 * e * (1 - r1))   
      else 
         vers_z   = 1 - 2 * r1 * ( 1 - xi ) / ( 1 + xi * (1 - 2 * r1) )
      endif
      sintheta = sqrt(1 - vers_z * vers_z)     

      call random_number(r2)
      phi = 2 * PI * r2
      vers_x = sintheta * cos(phi)
      vers_y = sintheta * sin(phi)
   
      if (debug) then
         print *,"IN"
         print *,"isatom          = ", isatom
         if (isatom) then
            print *,"speed / m s**-1 = ", xi
            print *,"energy / eu     = ", e
         else 
            print *,"xi              = ", xi
         endif
         print *,"r1              = ", r1
         print *,"r2              = ", r2
         print *,"OUT"
         print *,"sinthetasinphi  = ", vers_x
         print *,"sinthetacosphi  = ", vers_y
         print *,"costheta        = ", vers_z
         !call pause
      endif
  
   end subroutine anisotropic_scattering_electron


   subroutine scattering_elastic_comfor2tarfor(theta_tarfor, phi_tarfor, costheta_comfor, m_ratio)

      ! Transform the scattering and recoil angles of an elastic scattering from the COMFOR to the TARFOR
      ! TARFOR) target frame of reference, in which the target is at rest,
      !         and the z axis is parallel to the direction of the relative (projectile-target) velocity vector before scattering.
      ! COMFOR) center of mass frame of reference, in which the COM of projectile and target is at rest,
      !         and the z axis is parallel to the direction of the relative (projectile-target) velocity vector before scattering.     
      ! See e.g.: Lieberman M.A., Lichtenberg A.J. "Principles of Plasma Discharges and Materials Processing"
      ! (Wiley 2005)(ISBN 0471720011) eq. 3.2.10-3.2.11 , pag. 51-52

      !f2py intent(out) :: theta_tarfor, phi_tarfor
      !f2py intent(in)  :: costheta_comfor, m_ratio

      ! Parameters
      real(dp), intent(out) :: theta_tarfor    ! Scattering angle in the TARFOR
      real(dp), intent(out) :: phi_tarfor      ! Recoil angle     in the TARFOR
      real(dp), intent(in)  :: costheta_comfor ! Cosinus of the scattering angle in the COMFOR      
      real(dp), intent(in)  :: m_ratio         ! Mass ratio: projectile mass / target mass

      ! Local variables
      real(dp) :: sintheta_comfor ! Sinus of the scattering angle in the COMFOR
      real(dp) :: theta_comfor    ! Scattering angle in the COMFOR

      sintheta_comfor = sqrt(1 - costheta_comfor * costheta_comfor)
      theta_comfor    = acos(costheta_comfor)
      theta_tarfor    = atan( sintheta_comfor / ( m_ratio + costheta_comfor )  )
      phi_tarfor      = (PI - theta_comfor) / 2.0

   end subroutine scattering_elastic_comfor2tarfor

   
   subroutine scattering_tarfor2disfor(versorx1, versory1, versorz1, vx, vy, vz, v, &
                                      &sinthetacosphi_scatter, sinthetasinphi_scatter, costheta_scatter, &
                                      &debug)

      ! Calculates the new components of velocity after collision in the discharge frame of reference starting from:
      ! - the direction of projectile velocity before scattering in the discharge frame of reference
      ! - the scattering angles in the target frame of reference
      ! The discharge frame of reference (DISFOR) has
      ! - the z axis parallel to the electric field direction
      ! - the x and y axes parallel to the borders of the electrodes
      ! The target frame or reference (TARFOR) has
      ! - the z axis parallel to the direction of the relative projectile-target velocity vector before scattering
      ! - the x and y' axes directions are not important since scattering anomaly angle is always random

      !f2py intent(out) :: versorx1, versory1, versorz1
      !f2py intent(in)  :: vx, vy, vz, v
      !f2py intent(in)  :: sinthetacosphi_scatter, sinthetasinphi_scatter, costheta_scatter, debug


      ! Parameters
      real(dp), intent(out) :: versorx1               ! Projection on the DISFOR x axis of velocity versor after scattering
      real(dp), intent(out) :: versory1               ! Projection on the DISFOR y axis of velocity versor after scattering
      real(dp), intent(out) :: versorz1               ! Projection on the DISFOR z axis of velocity versor after scattering
      real(dp), intent(in)  :: vx                     ! Projection on the DISFOR x axis of velocity before scattering
      real(dp), intent(in)  :: vy                     ! Projection on the DISFOR y axis of velocity before scattering
      real(dp), intent(in)  :: vz                     ! Projection on the DISFOR z axis of velocity before scattering
      real(dp), intent(in)  :: v                      ! Modulus of velocity before scattering
      real(dp), intent(in)  :: sinthetacosphi_scatter ! Projection on the TARFOR x' axis of velocity versor after scattering
      real(dp), intent(in)  :: sinthetasinphi_scatter ! Projection on the TARFOR y' axis of velocity versor after scattering
      real(dp), intent(in)  :: costheta_scatter       ! Projection on the TARFOR z' axis of velocity versor after scattering
      logical,  intent(in)  :: debug                  ! Set debug mode 

      ! Local variables
      real(dp) :: costheta0   ! Projection on z axis   of velocity versor before scattering
      real(dp) :: sintheta0   ! Projection on x-y plane of velocity versor before scattering
      real(dp) :: cosphi0     ! Cosinus of the anomaly of velocity versor before scattering
      real(dp) :: sinphi0     ! Sinus of the anomaly of velocity versor before scattering
      real(dp) :: vr

      if (debug) then
         print *,""
         print *,"**SUBROUTINE SCATTERING 3D**"
         print *,""
         print *,"IN"
         print *,"vx                      = ", vx
         print *,"vy                      = ", vy
         print *,"vz                      = ", vz
         print *,"v                       = ", v
         print *,"sinthetacosphi_scatter  = ", sinthetacosphi_scatter
         print *,"sinthetasinphi_scatter  = ", sinthetasinphi_scatter
         print *,"costheta_scatter        = ", costheta_scatter
         print *,""
      endif

      ! Calculate angles before scattering
      ! vr is the projection of the velocity vector in the x-y plane
      if (v > 0) then
         costheta0 = vz / v
         vr = sqrt(vx * vx + vy * vy) 
         if (vr > 0) then 
            sintheta0 = vr / v
            cosphi0   = vx / vr
            sinphi0   = vy / vr
         else
            sintheta0 = 0
            cosphi0   = 1
            sinphi0   = 0
         endif 
      else
         costheta0 = 1
         sintheta0 = 0
         cosphi0   = 1
         sinphi0   = 0
      endif      

      ! Calculate the projections of velocity versor after scattering
      ! see for example Plasma Sources Sci. Technol.21 (2012) 055028
      versorx1 = - sinthetacosphi_scatter * costheta0 * cosphi0 - &
                  &sinthetasinphi_scatter             * sinphi0 + &
                  &costheta_scatter       * sintheta0 * cosphi0

      versory1 =   sinthetacosphi_scatter * costheta0 * sinphi0 + &
                  &sinthetasinphi_scatter             * cosphi0 + &
                  &costheta_scatter       * sintheta0 * sinphi0

      versorz1 =  - sinthetacosphi_scatter  * sintheta0 + costheta_scatter * costheta0 

      if (debug) then
         print *,"INTERMEDIATE"
         print *,"costheta0               = ", costheta0
         print *,"sintheta0               = ", sintheta0
         print *,"cosphi0                 = ", cosphi0 
         print *,"sinphi0                 = ", sinphi0
         print *,"vr                      = ", vr
         print *,""
         print *,"OUT"
         print *,"versorx1                = ", versorx1
         print *,"versory1                = ", versory1
         print *,"versorz1                = ", versorz1
         call pause
      endif
  
   end subroutine scattering_tarfor2disfor

end module f_scattering
