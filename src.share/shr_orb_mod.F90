!===============================================================================
! SVN $Id: shr_orb_mod.F90 25434 2010-11-04 22:46:24Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/release_tags/cesm1_2_x_n00_share3_130528/shr/shr_orb_mod.F90 $
!===============================================================================
!NOTES:
!05-17-24: XJ added "do_exo_higheccen" and "do_exo_MFA" adapted from Adams et al 2019
!          All changes commented by !ada:  
!06-11-24: XJ deleted the longterm evolution for Earth's orbit.

MODULE shr_orb_mod

   use shr_kind_mod
   use shr_sys_mod
   use shr_cal_mod
   use shr_const_mod
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit
   use exoplanet_mod,  only: exo_eccen, exo_obliq, exo_mvelp, &       		     
                             do_exo_synchronous,do_exo_higheccen, &
                             exo_porb, exo_ndays, do_exo_MFA

   IMPLICIT none

   !----------------------------------------------------------------------------
   ! PUBLIC: Interfaces and global data
   !----------------------------------------------------------------------------
   public :: shr_orb_cosz
   public :: shr_orb_params
   public :: shr_orb_decl
   public :: shr_orb_print

   real   (SHR_KIND_R8),public,parameter :: SHR_ORB_UNDEF_REAL = 1.e36_SHR_KIND_R8 ! undefined real 
   integer(SHR_KIND_IN),public,parameter :: SHR_ORB_UNDEF_INT  = 2000000000        ! undefined int

   !----------------------------------------------------------------------------
   ! PRIVATE: by default everything else is private to this module
   !----------------------------------------------------------------------------
   private

   real   (SHR_KIND_R8),parameter :: pi                 = SHR_CONST_PI
   real   (SHR_KIND_R8),parameter :: SHR_ORB_ECCEN_MIN  =   0.0_SHR_KIND_R8 ! min value for eccen
   real   (SHR_KIND_R8),parameter :: SHR_ORB_ECCEN_MAX  =  0.99_SHR_KIND_R8 ! max value for eccen
   real   (SHR_KIND_R8),parameter :: SHR_ORB_OBLIQ_MIN  = -90.0_SHR_KIND_R8 ! min value for obliq
   real   (SHR_KIND_R8),parameter :: SHR_ORB_OBLIQ_MAX  = +90.0_SHR_KIND_R8 ! max value for obliq
   real   (SHR_KIND_R8),parameter :: SHR_ORB_MVELP_MIN  =   0.0_SHR_KIND_R8 ! min value for mvelp
   real   (SHR_KIND_R8),parameter :: SHR_ORB_MVELP_MAX  = 360.0_SHR_KIND_R8 ! max value for mvelp

!===============================================================================
CONTAINS
!===============================================================================

real(SHR_KIND_R8) FUNCTION shr_orb_cosz(jday,lat,lon,declin)

   !----------------------------------------------------------------------------
   !
   ! FUNCTION to return the cosine of the solar zenith angle.
   ! Assumes 365.0 days/year.
   !
   !--------------- Code History -----------------------------------------------
   !
   ! Original Author: Brian Kauffman
   ! Date:            Jan/98
   ! History:         adapted from statement FUNCTION in share/orb_cosz.h
   !
   !----------------------------------------------------------------------------

   real   (SHR_KIND_R8),intent(in) :: jday   ! Julian cal day (1.xx to 365.xx)
   real   (SHR_KIND_R8),intent(in) :: lat    ! Centered latitude (radians)
   real   (SHR_KIND_R8),intent(in) :: lon    ! Centered longitude (radians)
   real   (SHR_KIND_R8),intent(in) :: declin ! Solar declination (radians)

   !ada:  what follows is a borrowing of the Kepler's equation code, with all
   !the appropriate variable declarations.
   real   (SHR_KIND_R8) :: ndays = exo_ndays  ! Length of rotation period in Earth days   
   real   (SHR_KIND_R8) :: dayspy = exo_porb  ! days per year
   real   (SHR_KIND_R8),parameter :: ve     = 0._SHR_KIND_R8   ! Calday of vernal equinox
   real   (SHR_KIND_R8),parameter :: eccen     = exo_eccen
                                                     ! assumes Jan 1 = calday 1
 
   real   (SHR_KIND_R8) ::   lambm  ! Lambda m, mean long of perihelion (rad)
   real   (SHR_KIND_R8) ::   lmm    ! Intermediate argument involving lambm
   real   (SHR_KIND_R8) ::   lamb   ! Lambda, the earths long of perihelion
   real   (SHR_KIND_R8) ::   invrho ! Inverse normalized sun/earth distance
   real   (SHR_KIND_R8) ::   sinl   ! Sine of lmm

   real   (SHR_KIND_R8) ::   E0     ! Initialization of eccentric anomaly
   real   (SHR_KIND_R8) ::   M0     ! Initialization of mean anomaly
   real   (SHR_KIND_R8) ::   E      ! Eccentric anomaly
   real   (SHR_KIND_R8) ::   toler  ! Tolerance for convergence
   integer(SHR_KIND_IN) ::   n      ! index   

   !----------------------------------------------------------------------------
   lambm = (jday - ve)*2._SHR_KIND_R8*pi/dayspy
   lmm   = lambm - pi

   ! ada: Use a numerical scheme to relate the mean anomaly (time fraction of
   ! the orbit) to eccentric anomaly (needed for the orbital longitude), as a
   ! function of eccentricity.
   E0 = 0.0_SHR_KIND_R8
   E = lmm
   toler = 1.0E-4_SHR_KIND_R8
   n = 0_SHR_KIND_IN

   DO WHILE ((abs(E - E0) .gt. toler) .and. (n .lt. 100_SHR_KIND_IN))
      E0 = E
      M0 = E0 - eccen*sin(E0)
      E = E0 + (lmm - M0)
      n = n + 1_SHR_KIND_IN
   END DO

!----------------------------------------------------------------------------


   if (do_exo_synchronous) then    !synchronous rotation
     shr_orb_cosz = sin(lat)*sin(declin) - &
     &          cos(lat)*cos(declin)*cos(lon)
   else   ! default, Earth present day orbit
!ada: Here "ndays" is the pseudosynchronous day length in Earth days, from
!exoplanet_mod. This sets the rotation period for the subsolar motion.
     shr_orb_cosz = sin(lat)*sin(declin) - &
     &      cos(lat)*cos(declin)*cos(jday/ndays*2.0_SHR_KIND_R8*pi - E +lon)
!     &              cos(lat)*cos(declin)*cos(jday*2.0_SHR_KIND_R8*pi + lon)
!      &              cos(lat)*cos(declin)*cos((jday-1.0)*2.0_SHR_KIND_R8*pi + lon)
   endif

END FUNCTION shr_orb_cosz

!===============================================================================

SUBROUTINE shr_orb_params( iyear_AD , eccen  , obliq , mvelp     ,     &
           &               obliqr   , lambm0 , mvelpp, log_print )

!-------------------------------------------------------------------------------
!
! Calculate earths orbital parameters using Dave Threshers formula which 
! came from Berger, Andre.  1978  "A Simple Algorithm to Compute Long-Term 
! Variations of Daily Insolation".  Contribution 18, Institute of Astronomy 
! and Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium
!
!------------------------------Code history-------------------------------------
!
! Original Author: Erik Kluzek
! Date:            Oct/97
!
!-------------------------------------------------------------------------------

   !----------------------------- Arguments ------------------------------------
   integer(SHR_KIND_IN),intent(in)    :: iyear_AD  ! Year to calculate orbit for
   real   (SHR_KIND_R8),intent(inout) :: eccen     ! orbital eccentricity
   real   (SHR_KIND_R8),intent(inout) :: obliq     ! obliquity in degrees
   real   (SHR_KIND_R8),intent(inout) :: mvelp     ! moving vernal equinox long
   real   (SHR_KIND_R8),intent(out)   :: obliqr    ! Earths obliquity in rad
   real   (SHR_KIND_R8),intent(out)   :: lambm0    ! Mean long of perihelion at
                                                   ! vernal equinox (radians)
   real   (SHR_KIND_R8),intent(out)   :: mvelpp    ! moving vernal equinox long
                                                   ! of perihelion plus pi (rad)
   logical             ,intent(in)    :: log_print ! Flags print of status/error
   real   (SHR_KIND_R8) :: degrad = pi/180._SHR_KIND_R8   ! degree to radian conversion factor
   character(len=*),parameter :: subname = '(shr_orb_params)'

   
   !-------------------------- Formats -----------------------------------------
   character(*),parameter :: svnID  = "SVN " // &
   "$Id: shr_orb_mod.F90 25434 2010-11-04 22:46:24Z tcraig $"
   character(*),parameter :: svnURL = "SVN <unknown URL>" 
!  character(*),parameter :: svnURL = "SVN " // &
!  "$URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/release_tags/cesm1_2_x_n02_share3_130715/shr/shr_orb_mod.F90 $"
   character(len=*),parameter :: F00 = "('(shr_orb_params) ',4a)"
   character(len=*),parameter :: F01 = "('(shr_orb_params) ',a,i9)"
   character(len=*),parameter :: F02 = "('(shr_orb_params) ',a,f6.3)"
   character(len=*),parameter :: F03 = "('(shr_orb_params) ',a,es14.6)"


   !! exoplanet_mod  linking for orbital property specification
   obliq = exo_obliq
   eccen = exo_eccen
   mvelp = exo_mvelp

 
   do while (mvelp .lt. 0.0_SHR_KIND_R8)
     mvelp = mvelp + 360.0_SHR_KIND_R8
   end do
   do while (mvelp .ge. 360.0_SHR_KIND_R8)
     mvelp = mvelp - 360.0_SHR_KIND_R8
   end do


   ! Orbit needs the obliquity in radians
 
   obliqr = obliq*degrad
 
   ! 180 degrees must be added to mvelp since observations are made from the
   ! earth and the sun is considered (wrongly for the algorithm) to go around
   ! the earth. For a more graphic explanation see Appendix B in:
   !
   ! A. Berger, M. Loutre and C. Tricot. 1993.  Insolation and Earth Orbital
   ! Periods.  J. of Geophysical Research 98:10,341-10,362.
   !
   ! Additionally, orbit will need this value in radians. So mvelp becomes
   ! mvelpp (mvelp plus pi)
 
   mvelpp = (mvelp + 180._SHR_KIND_R8)*degrad
   lambm0 = 0

   if ( log_print ) then
     write(s_logunit,F03) '------ Computed Orbital Parameters ------'
     write(s_logunit,F03) 'Eccentricity      = ',eccen
     write(s_logunit,F03) 'Obliquity (deg)   = ',obliq
     write(s_logunit,F03) 'Obliquity (rad)   = ',obliqr
     write(s_logunit,F03) 'Long of perh(deg) = ',mvelp
     write(s_logunit,F03) 'Long of perh(rad) = ',mvelpp
     write(s_logunit,F03) 'Long at v.e.(rad) = ',lambm0
     write(s_logunit,F03) '-----------------------------------------'
   end if
 
END SUBROUTINE shr_orb_params

!===============================================================================

SUBROUTINE shr_orb_decl(calday ,eccen ,mvelpp ,lambm0 ,obliqr ,delta ,eccf)

!-------------------------------------------------------------------------------
!
! Compute earth/orbit parameters using formula suggested by
! Duane Thresher.
!
!---------------------------Code history----------------------------------------
!
! Original version:  Erik Kluzek
! Date:              Oct/1997
!
!-------------------------------------------------------------------------------

   !------------------------------Arguments--------------------------------
   real   (SHR_KIND_R8),intent(in)  :: calday ! Calendar day, including fraction
   real   (SHR_KIND_R8),intent(in)  :: eccen  ! Eccentricity
   real   (SHR_KIND_R8),intent(in)  :: obliqr ! Earths obliquity in radians
   real   (SHR_KIND_R8),intent(in)  :: lambm0 ! Mean long of perihelion at the 
                                              ! vernal equinox (radians)
   real   (SHR_KIND_R8),intent(in)  :: mvelpp ! moving vernal equinox longitude
                                              ! of perihelion plus pi (radians)
   real   (SHR_KIND_R8),intent(out) :: delta  ! Solar declination angle in rad
   real   (SHR_KIND_R8),intent(out) :: eccf   ! Earth-sun distance factor (ie. (1/r)**2)
 
   !---------------------------Local variables-----------------------------

   real   (SHR_KIND_R8),parameter :: dayspy = exo_porb                   ! days per year
   real   (SHR_KIND_R8),parameter :: ve     = 1.0_SHR_KIND_R8   ! Calday of vernal equinox

   !!real   (SHR_KIND_R8),parameter :: dayspy = 365.0_SHR_KIND_R8  ! days per year
   !!real   (SHR_KIND_R8),parameter :: ve     = 80.5_SHR_KIND_R8   ! Calday of vernal equinox
                                                     ! assumes Jan 1 = calday 1
 
   real   (SHR_KIND_R8) ::   lambm  ! Lambda m, mean long of perihelion (rad)
   real   (SHR_KIND_R8) ::   lmm    ! Intermediate argument involving lambm
   real   (SHR_KIND_R8) ::   lamb   ! Lambda, the earths long of perihelion
   real   (SHR_KIND_R8) ::   invrho ! Inverse normalized sun/earth distance
   real   (SHR_KIND_R8) ::   sinl   ! Sine of lmm
   
   ! ada: Variables added for Kepler equation  
   real   (SHR_KIND_R8) ::   E0     ! Initialization of eccentric anomaly
   real   (SHR_KIND_R8) ::   M0     ! Initialization of mean anomaly
   real   (SHR_KIND_R8) ::   E      ! Eccentric anomaly
   real   (SHR_KIND_R8) ::   toler  ! Tolerance for convergence
   integer(SHR_KIND_IN) ::   n      ! index
 
 
   ! Compute eccentricity factor and solar declination using
   ! day value where a round day (such as 213.0) refers to 0z at
   ! Greenwich longitude.
   !
   ! Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
   ! Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
   ! 35:2362-2367.
   !
   ! To get the earths true longitude (position in orbit; lambda in Berger 
   ! 1978) which is necessary to find the eccentricity factor and declination,
   ! must first calculate the mean longitude (lambda m in Berger 1978) at
   ! the present day.  This is done by adding to lambm0 (the mean longitude
   ! at the vernal equinox, set as March 21 at noon, when lambda=0; in radians)
   ! an increment (delta lambda m in Berger 1978) that is the number of
   ! days past or before (a negative increment) the vernal equinox divided by
   ! the days in a model year times the 2*pi radians in a complete orbit.
 
   lambm = lambm0 + (calday - ve)*2._SHR_KIND_R8*pi/dayspy
   lmm   = lambm  - mvelpp


   ! Adams: If do_exo_higheccen, use a numerical scheme to relate the mean
   ! anomaly (time fraction of the orbit) to eccentric anomaly (needed for the
   ! orbital longitude), as a function of eccentricity.

   if (do_exo_higheccen) then
      E0 = 0.0_SHR_KIND_R8
      E = lmm
      toler = 1.0E-4_SHR_KIND_R8
      n = 0_SHR_KIND_IN

      DO WHILE ((abs(E - E0) .gt. toler) .and. (n .lt. 100_SHR_KIND_IN))
         E0 = E
         M0 = E0 - eccen*sin(E0)
         E = E0 + (lmm - M0)
         n = n + 1_SHR_KIND_IN
      END DO

      lamb = pi - 2.0_SHR_KIND_R8 * atan2(sqrt(1.0_SHR_KIND_R8 -eccen)*cos(E/2.0_SHR_KIND_R8), sqrt(1.0_SHR_KIND_R8 + eccen)*sin(E/2.0_SHR_KIND_R8))

   else        
       ! The earths true longitude, in radians, is then found from
       ! the formula in Berger 1978:
 
       sinl  = sin(lmm)
       lamb  = lambm  + eccen*(2._SHR_KIND_R8*sinl + eccen*(1.25_SHR_KIND_R8*sin(2._SHR_KIND_R8*lmm)  &
       &     + eccen*((13.0_SHR_KIND_R8/12.0_SHR_KIND_R8)*sin(3._SHR_KIND_R8*lmm) - 0.25_SHR_KIND_R8*sinl)))
   endif

   ! Using the obliquity, eccentricity, moving vernal equinox longitude of
   ! perihelion (plus), and earths true longitude, the declination (delta)
   ! and the normalized earth/sun distance (rho in Berger 1978; actually inverse
   ! rho will be used), and thus the eccentricity factor (eccf), can be 
   ! calculated from formulas given in Berger 1978.
 
   invrho = (1._SHR_KIND_R8 + eccen*cos(lamb)) / (1._SHR_KIND_R8 - eccen*eccen)
 
   ! Set solar declination and eccentricity factor
 
   delta  = asin(sin(obliqr)*sin(lamb))
   
   if (do_exo_MFA) then
      eccf   = invrho*invrho*sqrt(1._SHR_KIND_R8 - eccen*eccen)
   else
      eccf   = invrho*invrho
   endif
   
   return
 
END SUBROUTINE shr_orb_decl

!===============================================================================

SUBROUTINE shr_orb_print( iyear_AD, eccen, obliq, mvelp )

!-------------------------------------------------------------------------------
!
! Print out the information on the Earths input orbital characteristics
!
!---------------------------Code history----------------------------------------
!
! Original version:  Erik Kluzek
! Date:              Oct/1997
!
!-------------------------------------------------------------------------------

   !---------------------------Arguments----------------------------------------
   integer(SHR_KIND_IN),intent(in) :: iyear_AD ! requested Year (AD)
   real   (SHR_KIND_R8),intent(in) :: eccen    ! eccentricity (unitless) 
                                               ! (typically 0 to 0.1)
   real   (SHR_KIND_R8),intent(in) :: obliq    ! obliquity (-90 to +90 degrees) 
                                               ! typically 22-26
   real   (SHR_KIND_R8),intent(in) :: mvelp    ! moving vernal equinox at perhel
                                               ! (0 to 360 degrees)
   !-------------------------- Formats -----------------------------------------
   character(len=*),parameter :: F00 = "('(shr_orb_print) ',4a)"
   character(len=*),parameter :: F01 = "('(shr_orb_print) ',a,i9.4)"
   character(len=*),parameter :: F02 = "('(shr_orb_print) ',a,f6.3)"
   character(len=*),parameter :: F03 = "('(shr_orb_print) ',a,es14.6)"
   !----------------------------------------------------------------------------
 
   if (s_loglev > 0) then
   if ( iyear_AD .ne. SHR_ORB_UNDEF_INT ) then
     if ( iyear_AD > 0 ) then
       write(s_logunit,F01) 'Orbital parameters calculated for year: AD ',iyear_AD
     else
       write(s_logunit,F01) 'Orbital parameters calculated for year: BC ',iyear_AD
     end if
   else if ( obliq /= SHR_ORB_UNDEF_REAL ) then
     write(s_logunit,F03) 'Orbital parameters: '
     write(s_logunit,F03) 'Obliquity (degree):              ', obliq
     write(s_logunit,F03) 'Eccentricity (unitless):         ', eccen
     write(s_logunit,F03) 'Long. of moving Perhelion (deg): ', mvelp
   else
     write(s_logunit,F03) 'Orbit parameters not set!'
   end if
   endif
 
END SUBROUTINE shr_orb_print
!===============================================================================


END MODULE shr_orb_mod
