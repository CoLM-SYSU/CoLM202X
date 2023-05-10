MODULE MOD_Turbulence_LEddy

!-----------------------------------------------------------------------
 use precision
 IMPLICIT NONE
 SAVE

! PUBLIC MEMBER FUNCTIONS:
  public :: moninobuk_leddy
  public :: moninobukm_leddy


! PRIVATE MEMBER FUNCTIONS:
  private :: psi


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

 subroutine moninobuk_leddy(hu,ht,hq,displa,z0m,z0h,z0q,obu,um, hpbl, &
                            ustar,fh2m,fq2m,fm10m,fm,fh,fq)

! ======================================================================
!
! Implement the LZD2022 scheme (Liu et al., 2022), which accounts for large
! eddy effects by inlcuding the boundary layer height in the phim function,
! to compute friction velocity, relation for potential temperature and
! humidity profiles of surface boundary layer.
!
! References:
! [1] Zeng et al., 1998: Intercomparison of bulk aerodynamic algorithms
!     for the computation of sea surface fluxes using TOGA CORE and TAO data. 
!     J. Climate, 11: 2628-2644.
! [2] Liu et al., 2022: A surface flux estimation scheme accounting for
!     large-eddy effects for land surface modeling. GRL, 49, e2022GL101754.
!
! Created by Shaofeng Liu, May 5, 2023
!
! ======================================================================

  use precision
  use PhysicalConstants, only : vonkar
  implicit none

! ---------------------- dummy argument --------------------------------

  real(r8), INTENT(in) :: hu       ! observational height of wind [m]
  real(r8), INTENT(in) :: ht       ! observational height of temperature [m]
  real(r8), INTENT(in) :: hq       ! observational height of humidity [m]
  real(r8), INTENT(in) :: displa   ! displacement height [m]
  real(r8), INTENT(in) :: z0m      ! roughness length, momentum [m]
  real(r8), INTENT(in) :: z0h      ! roughness length, sensible heat [m]
  real(r8), INTENT(in) :: z0q      ! roughness length, latent heat [m]
  real(r8), INTENT(in) :: obu      ! monin-obukhov length (m)
  real(r8), INTENT(in) :: um       ! wind speed including the stablity effect [m/s]
  real(r8), INTENT(in) :: hpbl     ! atmospheric boundary layer height [m]

  real(r8), INTENT(out) :: ustar   ! friction velocity [m/s]
  real(r8), INTENT(out) :: fh2m    ! relation for temperature at 2m
  real(r8), INTENT(out) :: fq2m    ! relation for specific humidity at 2m
  real(r8), INTENT(out) :: fm10m   ! integral of profile function for momentum at 10m
  real(r8), INTENT(out) :: fm      ! integral of profile function for momentum
  real(r8), INTENT(out) :: fh      ! integral of profile function for heat
  real(r8), INTENT(out) :: fq      ! integral of profile function for moisture

!------------------------ local variables ------------------------------

  real(r8) zldis  ! reference height "minus" zero displacement heght [m]
  real(r8) zetam, &
           zetam2 ! transition point of flux-gradient relation (wind profile)
  real(r8) zetat  ! transition point of flux-gradient relation (temp. profile)
  real(r8) zeta   ! dimensionless height used in Monin-Obukhov theory
  real(r8) Bm     ! Coefficient of the LZD2022 scheme: Bm = 0.0047*(-hpbl/L) + 0.1854
  real(r8) Bm2    ! max(Bm, 0.2722)

! real(r8), external :: psi    ! stability function for unstable case
!-----------------------------------------------------------------------
! adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.

! wind profile
        zldis=hu-displa
        zeta=zldis/obu
!
! Begin: Shaofeng Liu, 2023.05.05
!
		Bm     = 0.0047 * (-hpbl/obu) + 0.1854
		zetam  = 0.5*Bm**4 * ( -16. - (256. + 4.*Bm**(-4)**0.5) )
		Bm2    = max(Bm, 0.2722)
		zetam2 = min(zetam, -0.13)

        if(zeta < zetam2)then           ! zeta < zetam2
          fm    = log(zetam2*obu/z0m) - psi(1,zetam2) &
                + psi(1,z0m/obu) - 2.*Bm2 * ( (-zeta)**(-0.5)-(-zetam2)**(-0.5) )
          ustar = vonkar*um / fm
!
! End: Shaofeng Liu, 2023.05.05
!
        else if(zeta < 0.)then          ! zetam2 <= zeta < 0
          fm    = log(zldis/z0m) - psi(1,zeta) + psi(1,z0m/obu)
          ustar = vonkar*um / fm
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fm    = log(zldis/z0m) + 5.*zeta - 5.*z0m/obu
          ustar = vonkar*um / fm
        else                            !  1 < zeta, phi=5+zeta
          fm    = log(obu/z0m) + 5. - 5.*z0m/obu + (5.*log(zeta)+zeta-1.)
          ustar = vonkar*um / fm
        endif

        ! for 10 meter wind-velocity
        zldis=10.+z0m
        zeta=zldis/obu
        zetam=1.574
        if(zeta < -zetam)then           ! zeta < -1
          fm10m  = log(-zetam*obu/z0m) - psi(1,-zetam) &
                + psi(1,z0m/obu) + 1.14*((-zeta)**0.333-(zetam)**0.333)
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          fm10m  = log(zldis/z0m) - psi(1,zeta) + psi(1,z0m/obu)
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fm10m  = log(zldis/z0m) + 5.*zeta - 5.*z0m/obu
        else                            !  1 < zeta, phi=5+zeta
          fm10m  = log(obu/z0m) + 5. - 5.*z0m/obu + (5.*log(zeta)+zeta-1.)
        endif

! temperature profile
        zldis=ht-displa
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
          fh    = log(-zetat*obu/z0h)-psi(2,-zetat) &
                + psi(2,z0h/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333))
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          fh    = log(zldis/z0h) - psi(2,zeta) + psi(2,z0h/obu)
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fh    = log(zldis/z0h) + 5.*zeta - 5.*z0h/obu
        else                            !  1 < zeta, phi=5+zeta
          fh    = log(obu/z0h) + 5. - 5.*z0h/obu + (5.*log(zeta)+zeta-1.)
        endif

        ! for 2 meter screen temperature
        zldis=2.+z0h  ! ht-displa
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
          fh2m = log(-zetat*obu/z0h)-psi(2,-zetat) &
                + psi(2,z0h/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333))
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          fh2m = log(zldis/z0h) - psi(2,zeta) + psi(2,z0h/obu)
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fh2m = log(zldis/z0h) + 5.*zeta - 5.*z0h/obu
        else                            !  1 < zeta, phi=5+zeta
          fh2m = log(obu/z0h) + 5. - 5.*z0h/obu + (5.*log(zeta)+zeta-1.)
        endif

! humidity profile
        zldis=hq-displa
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
          fq    = log(-zetat*obu/z0q) - psi(2,-zetat) &
                + psi(2,z0q/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333))
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          fq    = log(zldis/z0q) - psi(2,zeta) + psi(2,z0q/obu)
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fq    = log(zldis/z0q) + 5.*zeta - 5.*z0q/obu
        else                            !  1 < zeta, phi=5+zeta
          fq    = log(obu/z0q) + 5. - 5.*z0q/obu + (5.*log(zeta)+zeta-1.)
        endif

      ! for 2 meter screen humidity
        zldis=2.+z0h
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
           fq2m = log(-zetat*obu/z0q)-psi(2,-zetat) &
                 + psi(2,z0q/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333))
        elseif (zeta < 0.) then         ! -1 <= zeta < 0
           fq2m = log(zldis/z0q)-psi(2,zeta)+psi(2,z0q/obu)
        else if (zeta <= 1.) then       !  0 <= zeta <= 1
           fq2m = log(zldis/z0q)+5.*zeta-5.*z0q/obu
        else                            ! 1 < zeta, phi=5+zeta
           fq2m = log(obu/z0q)+5.-5.*z0q/obu+(5.*log(zeta)+zeta-1.)
        endif

 end subroutine moninobuk_leddy


 subroutine moninobukm_leddy(hu,ht,hq,displa,z0m,z0h,z0q,obu,um,displat,z0mt, hpbl, &
                             ustar,fh2m,fq2m,htop,fmtop,fm,fh,fq,fht,fqt,phih)

! ======================================================================
!
! !DESCRIPTION:
!
!
! Original author : Yongjiu Dai, September 15, 1999
!
! calculation of friction velocity, relation for potential temperatur
! and humidity profiles of surface boundary layer.
! the scheme is based on the work of Zeng et al. (1998):
! Intercomparison of bulk aerodynamic algorithms for the computation
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate, Vol. 11: 2628-2644
!
! REVISIONS:
! Hua Yuan, 09/2017: adapted from moninobuk FUNCTION to calculate canopy top
!                    fm, fq and phih for roughness sublayer u/k profile calculation
! Shaofeng Liu, 05/2023: implement the LZD2022 scheme (Liu et al., 2022), which 
!						 accounts for large eddy effects by including the 
!                        boundary leyer height in the phim function. 
! ======================================================================

  use precision
  use PhysicalConstants, only : vonkar
  implicit none

! ---------------------- dummy argument --------------------------------

  real(r8), INTENT(in) :: hu       ! observational height of wind [m]
  real(r8), INTENT(in) :: ht       ! observational height of temperature [m]
  real(r8), INTENT(in) :: hq       ! observational height of humidity [m]
  real(r8), INTENT(in) :: displa   ! displacement height [m]
  real(r8), INTENT(in) :: displat  ! displacement height of the top layer [m]
  real(r8), INTENT(in) :: z0m      ! roughness length, momentum [m]
  real(r8), INTENT(in) :: z0h      ! roughness length, sensible heat [m]
  real(r8), INTENT(in) :: z0q      ! roughness length, latent heat [m]
  real(r8), INTENT(in) :: z0mt     ! roughness length of the top layer, latent heat [m]
  real(r8), INTENT(in) :: htop     ! canopy top height of the top layer [m]
  real(r8), INTENT(in) :: obu      ! monin-obukhov length (m)
  real(r8), INTENT(in) :: um       ! wind speed including the stablity effect [m/s]
  real(r8), INTENT(in) :: hpbl     ! atmospheric boundary layer height [m]

  real(r8), INTENT(out) :: ustar   ! friction velocity [m/s]
  real(r8), INTENT(out) :: fh2m    ! relation for temperature at 2m
  real(r8), INTENT(out) :: fq2m    ! relation for specific humidity at 2m
  real(r8), INTENT(out) :: fmtop   ! integral of profile function for momentum at 10m
  real(r8), INTENT(out) :: fm      ! integral of profile function for momentum
  real(r8), INTENT(out) :: fh      ! integral of profile function for heat
  real(r8), INTENT(out) :: fq      ! integral of profile function for moisture
  real(r8), INTENT(out) :: fht     ! integral of profile function for heat at the top layer
  real(r8), INTENT(out) :: fqt     ! integral of profile function for moisture at the top layer
  real(r8), INTENT(out) :: phih    ! phi(h), similarity function for sensible heat

!------------------------ local variables ------------------------------

  real(r8) zldis  ! reference height "minus" zero displacement heght [m]
  real(r8) zetam, &
           zetam2 ! transition point of flux-gradient relation (wind profile)
  real(r8) zetat  ! transition point of flux-gradient relation (temp. profile)
  real(r8) zeta   ! dimensionless height used in Monin-Obukhov theory
  real(r8) Bm     ! Coefficient of the LZD2022 scheme: Bm = 0.0047*(-hpbl/L) + 0.1854
  real(r8) Bm2    ! max(Bm, 0.2722)

! real(r8), external :: psi    ! stability function for unstable case
!-----------------------------------------------------------------------
! adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.

! wind profile
        zldis=hu-displa
        zeta=zldis/obu
!
! Begin: Shaofeng Liu, 2023.05.05
!
		Bm     = 0.0047 * (-hpbl/obu) + 0.1854
		zetam  = 0.5*Bm**4 * ( -16. - (256. + 4.*Bm**(-4)**0.5) )
		Bm2    = max(Bm, 0.2722)
		zetam2 = min(zetam, -0.13)

        if(zeta < zetam2)then           ! zeta < zetam2
          fm    = log(zetam2*obu/z0m) - psi(1,zetam2) &
                + psi(1,z0m/obu) - 2.*Bm2 * ( (-zeta)**(-0.5)-(-zetam2)**(-0.5) )
          ustar = vonkar*um / fm
!
! End: Shaofeng Liu, 2023.05.05
!
        else if(zeta < 0.)then          ! zetam2 <= zeta < 0
          fm    = log(zldis/z0m) - psi(1,zeta) + psi(1,z0m/obu)
          ustar = vonkar*um / fm
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fm    = log(zldis/z0m) + 5.*zeta - 5.*z0m/obu
          ustar = vonkar*um / fm
        else                            !  1 < zeta, phi=5+zeta
          fm    = log(obu/z0m) + 5. - 5.*z0m/obu + (5.*log(zeta)+zeta-1.)
          ustar = vonkar*um / fm
        endif

        ! for canopy top wind-velocity
        !NOTE: changed for canopy top wind-velocity (no wake assumed)
        zldis=htop-displa
        zeta=zldis/obu
        zetam=1.574
        if(zeta < -zetam)then           ! zeta < -1
          fmtop  = log(-zetam*obu/z0m) - psi(1,-zetam) &
                + psi(1,z0m/obu) + 1.14*((-zeta)**0.333-(zetam)**0.333)
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          fmtop  = log(zldis/z0m) - psi(1,zeta) + psi(1,z0m/obu)
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fmtop  = log(zldis/z0m) + 5.*zeta - 5.*z0m/obu
        else                            !  1 < zeta, phi=5+zeta
          fmtop  = log(obu/z0m) + 5. - 5.*z0m/obu + (5.*log(zeta)+zeta-1.)
        endif

! temperature profile
        zldis=ht-displa
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
          fh    = log(-zetat*obu/z0h)-psi(2,-zetat) &
                + psi(2,z0h/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333))
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          fh    = log(zldis/z0h) - psi(2,zeta) + psi(2,z0h/obu)
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fh    = log(zldis/z0h) + 5.*zeta - 5.*z0h/obu
        else                            !  1 < zeta, phi=5+zeta
          fh    = log(obu/z0h) + 5. - 5.*z0h/obu + (5.*log(zeta)+zeta-1.)
        endif

        ! for 2 meter screen temperature
        zldis=2.+z0h  ! ht-displa
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
          fh2m = log(-zetat*obu/z0h)-psi(2,-zetat) &
                + psi(2,z0h/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333))
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          fh2m = log(zldis/z0h) - psi(2,zeta) + psi(2,z0h/obu)
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fh2m = log(zldis/z0h) + 5.*zeta - 5.*z0h/obu
        else                            !  1 < zeta, phi=5+zeta
          fh2m = log(obu/z0h) + 5. - 5.*z0h/obu + (5.*log(zeta)+zeta-1.)
        endif

        ! for top layer temperature
        zldis=displat+z0mt-displa  ! ht-displa
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
          fht = log(-zetat*obu/z0h)-psi(2,-zetat) &
                + psi(2,z0h/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333))
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          fht = log(zldis/z0h) - psi(2,zeta) + psi(2,z0h/obu)
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fht = log(zldis/z0h) + 5.*zeta - 5.*z0h/obu
        else                            !  1 < zeta, phi=5+zeta
          fht = log(obu/z0h) + 5. - 5.*z0h/obu + (5.*log(zeta)+zeta-1.)
        endif

        ! for canopy top phi(h)
        ! CESM TECH NOTE EQ. (5.31)
        zldis=htop-displa  ! ht-displa
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
           phih = 0.9*vonkar**(1.333)*(-zeta)**(-0.333)
        else if(zeta < 0.)then          ! -1 <= zeta < 0
           phih = (1. - 16.*zeta)**(-0.5)
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
           phih = 1. + 5.*zeta
        else                            !  1 < zeta, phi=5+zeta
           phih = 5. + zeta
        endif

! humidity profile
        zldis=hq-displa
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
          fq    = log(-zetat*obu/z0q) - psi(2,-zetat) &
                + psi(2,z0q/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333))
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          fq    = log(zldis/z0q) - psi(2,zeta) + psi(2,z0q/obu)
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fq    = log(zldis/z0q) + 5.*zeta - 5.*z0q/obu
        else                            !  1 < zeta, phi=5+zeta
          fq    = log(obu/z0q) + 5. - 5.*z0q/obu + (5.*log(zeta)+zeta-1.)
        endif

        ! for 2 meter screen humidity
        zldis=2.+z0h
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
           fq2m = log(-zetat*obu/z0q)-psi(2,-zetat) &
                 + psi(2,z0q/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333))
        elseif (zeta < 0.) then         ! -1 <= zeta < 0
           fq2m = log(zldis/z0q)-psi(2,zeta)+psi(2,z0q/obu)
        else if (zeta <= 1.) then       !  0 <= zeta <= 1
           fq2m = log(zldis/z0q)+5.*zeta-5.*z0q/obu
        else                            ! 1 < zeta, phi=5+zeta
           fq2m = log(obu/z0q)+5.-5.*z0q/obu+(5.*log(zeta)+zeta-1.)
        endif

        ! for top layer humidity
        zldis=displat+z0mt-displa  ! ht-displa
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
           fqt = log(-zetat*obu/z0q)-psi(2,-zetat) &
                 + psi(2,z0q/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333))
        elseif (zeta < 0.) then         ! -1 <= zeta < 0
           fqt = log(zldis/z0q)-psi(2,zeta)+psi(2,z0q/obu)
        else if (zeta <= 1.) then       !  0 <= zeta <= 1
           fqt = log(zldis/z0q)+5.*zeta-5.*z0q/obu
        else                            ! 1 < zeta, phi=5+zeta
           fqt = log(obu/z0q)+5.-5.*z0q/obu+(5.*log(zeta)+zeta-1.)
        endif

 end subroutine moninobukm_leddy


 
 real(r8) function psi(k,zeta)

!=======================================================================
! stability function for unstable case (rib < 0)

  use precision
  implicit none

  integer k
  real(r8) zeta  ! dimensionless height used in Monin-Obukhov theory
  real(r8) chik  !

      chik = (1.-16.*zeta)**0.25
      if(k == 1)then
        psi = 2.*log((1.+chik)*0.5)+log((1.+chik*chik)*0.5)-2.*atan(chik)+2.*atan(1.)
      else
        psi = 2.*log((1.+chik*chik)*0.5)
      endif

 end function psi


END MODULE MOD_Turbulence_LEddy
! --------- EOP ------------

