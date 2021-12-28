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


module f_scattering

   use f_precision

   use f_constants

   use f_io

   implicit none

contains

   subroutine scattering_3d(versorx1, versory1, versorz1, vx, vy, vz, v, &
                           &sinthetacosphi_scatter, sinthetasinphi_scatter, costheta_scatter, &
                           &debug)

      ! Calculates the new components of velocity after collision in the laboratory frame of reference 
      ! starting from the angles describing the new velocity in the frame of reference of old velocity

      !f2py intent(out) :: versorx1, versory1, versorz1
      !f2py intent(in)  :: vx, vy, vz, v
      !f2py intent(in)  :: sinthetacosphi_scatter, sinthetasinphi_scatter, costheta_scatter, debug


      ! Parameters
      real(dp), intent(out) :: versorx1               ! Projection on x axis of velocity versor after scattering
      real(dp), intent(out) :: versory1               ! Projection on y axis of velocity versor after scattering
      real(dp), intent(out) :: versorz1               ! Projection on z axis of velocity versor after scattering
      real(dp), intent(in)  :: vx                     ! Projection on x axis of velocity before scattering
      real(dp), intent(in)  :: vy                     ! Projection on y axis of velocity before scattering
      real(dp), intent(in)  :: vz                     ! Projection on z axis of velocity before scattering
      real(dp), intent(in)  :: v                      ! Modulus of velocity before scattering
      real(dp), intent(in)  :: sinthetacosphi_scatter ! Projection on x axis of velocity versor after scattering
                                                      !    in the reference of the particle velocity
      real(dp), intent(in)  :: sinthetasinphi_scatter ! Projection on y axis of velocity versor after scattering
                                                      !    in the reference of the particle velocity
      real(dp), intent(in)  :: costheta_scatter       ! Projection on z axis of velocity versor after scattering
                                                      !    in the reference of the particle velocity
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
  
   end subroutine scattering_3d


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


   subroutine elastic_scattering_electron_neutral(vers_x, vers_y, vers_z, e, isatom, debug)

      ! Calculates the scattering angles for elastic collision of an electron with an atom or a nonpolar molecule
      ! the used formula is correct for elastic scattering with atoms (screened Coulomb potential) and is taken from the paper
      ! PHYSICAL REVIEW E 65 (2002) 037402

      !f2py intent(out) :: vers_x, vers_y, vers_z
      !f2py intent(in)  :: e, atom, debug

      ! Parameters
      real(dp), intent(out) :: vers_x   ! Projection on x axis of velocity vector after scattering (sinthetacosphi)
      real(dp), intent(out) :: vers_y   ! Projection on y axis of velocity vector after scattering (sinthetasinphi)
      real(dp), intent(out) :: vers_z   ! Projection on z axis of velocity vector after scattering (costheta)
      real(dp), intent(in)  :: e        ! Parameter the meaning of which depends on the type of target particle
                                        !   atom     -> electron energy espressed in atomic energy units (au = 27.21 eV)
                                        !   molecule -> value of the xi function for the current electron energy 
      logical,  intent(in)  :: isatom   ! if .true. the scattering is calculated for an atom, if .false. for a molecule
      logical,  intent(in)  :: debug    ! Set debug mode 

      ! Local variables
      real(dp) :: sintheta, phi
      real(sp) :: r1,r2 

      if (debug) then
         print *,""
         print *,"**SUBROUTINE ELASTIC SCATTERING NEUTRAL**"
      endif

      ! Calculate scattering angles
      call random_number(r1)

      if (isatom) then 
         vers_z   = 1 - 2 * r1 / (1 + 8 * e * (1 - r1))   
      else 
         vers_z   = 1 - 2 * r1 * ( 1 - e ) / ( 1 + e * (1 - 2 * r1) )
      endif
      sintheta = sqrt(1 - vers_z * vers_z)     

      call random_number(r2)
      phi = 2 * PI * r2
      vers_x = sintheta * cos(phi)
      vers_y = sintheta * sin(phi)
   
      if (debug) then
         print *,"IN"
         print *,"isatom         = ", isatom
         print *,"energy / xi    = ", e
         print *,"r1             = ", r1
         print *,"r2             = ", r2
         print *,"OUT"
         print *,"sinthetasinphi = ", vers_x
         print *,"sinthetacosphi = ", vers_y
         print *,"costheta       = ", vers_z
         call pause
      endif
  
   end subroutine elastic_scattering_electron_neutral


end module f_scattering
