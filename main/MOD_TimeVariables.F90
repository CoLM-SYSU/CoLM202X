#include <define.h>

MODULE MOD_TimeVariables
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

use precision
use timemanager
#ifdef PFT_CLASSIFICATION
USE MOD_PFTimeVars
#endif
#ifdef PC_CLASSIFICATION
USE MOD_PCTimeVars
#endif
#ifdef BGC
USE MOD_BGCTimeVars
#endif
#ifdef LATERAL_FLOW
USE MOD_HydroTimeVars
#endif
#ifdef URBAN_MODEL
USE MOD_UrbanTimeVars
#endif

IMPLICIT NONE
SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run
      real(r8), allocatable :: z_sno    (:,:)   ! node depth [m]
      real(r8), allocatable :: dz_sno   (:,:)   ! interface depth [m]
      real(r8), allocatable :: t_soisno (:,:)   ! soil temperature [K]
      real(r8), allocatable :: wliq_soisno(:,:) ! liquid water in layers [kg/m2]
      real(r8), allocatable :: wice_soisno(:,:) ! ice lens in layers [kg/m2]
      real(r8), allocatable :: h2osoi (:,:)     ! volumetric soil water in layers [m3/m3]
      real(r8), allocatable :: smp(:,:)         ! soil matrix potential [mm]
      real(r8), allocatable :: hk (:,:)         ! hydraulic conductivity [mm h2o/s]
      real(r8), allocatable :: rootr(:,:)       ! water exchange between soil and root. Positive: soil->root [?]
#ifdef PLANT_HYDRAULIC_STRESS
      real(r8), allocatable :: vegwp(:,:)       ! vegetation water potential [mm]
      real(r8), allocatable :: gs0sun   (:)     ! working copy of sunlit stomata conductance
      real(r8), allocatable :: gs0sha   (:)     ! working copy of shalit stomata conductance
#endif
#ifdef OzoneStress
      real(r8), allocatable :: o3coefv_sun(:) ! Ozone stress factor for photosynthesis on sunlit leaf
      real(r8), allocatable :: o3coefv_sha(:) ! Ozone stress factor for photosynthesis on shaded leaf
      real(r8), allocatable :: o3coefg_sun(:) ! Ozone stress factor for stomata on sunlit leaf
      real(r8), allocatable :: o3coefg_sha(:) ! Ozone stress factor for stomata on shaded leaf
      real(r8), allocatable :: lai_old    (:) ! lai in last time step
      real(r8), allocatable :: o3uptakesun(:) ! Ozone does, sunlit leaf (mmol O3/m^2)
      real(r8), allocatable :: o3uptakesha(:) ! Ozone does, shaded leaf (mmol O3/m^2)
#endif
      real(r8), allocatable :: rstfacsun_out(:) ! factor of soil water stress on sunlit leaf
      real(r8), allocatable :: rstfacsha_out(:) ! factor of soil water stress on shaded leaf
      real(r8), allocatable :: gssun_out(:)     ! stomata conductance on sunlit leaf
      real(r8), allocatable :: gssha_out(:)     ! stomata conductance on shaded leaf
      real(r8), allocatable :: t_grnd   (:)     ! ground surface temperature [K]

#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
      real(r8), allocatable :: assim_RuBP_sun_out   (:) !1
      real(r8), allocatable :: assim_RuBP_sha_out   (:) !1
      real(r8), allocatable :: assim_Rubisco_sun_out(:) !1
      real(r8), allocatable :: assim_Rubisco_sha_out(:) !1
      real(r8), allocatable :: assimsun_out         (:) !1
      real(r8), allocatable :: assimsha_out         (:) !1
      real(r8), allocatable :: etrsun_out           (:) !1
      real(r8), allocatable :: etrsha_out           (:) !1
      real(r8), allocatable :: cisun_out            (:) !1
      real(r8), allocatable :: cisha_out            (:) !1
      real(r8), allocatable :: Dsun_out             (:) !1
      real(r8), allocatable :: Dsha_out             (:) !1
      real(r8), allocatable :: gammasun_out         (:) !1
      real(r8), allocatable :: gammasha_out         (:) !1
      real(r8), allocatable :: lambdasun_out        (:)
      real(r8), allocatable :: lambdasha_out        (:)
      real(r8), allocatable :: lambda_out           (:)
#endif
#endif

      real(r8), allocatable :: tleaf    (:)     ! leaf temperature [K]
      real(r8), allocatable :: ldew     (:)     ! depth of water on foliage [mm]
      real(r8), allocatable :: ldew_rain(:)     ! depth of rain on foliage [mm]
      real(r8), allocatable :: ldew_snow(:)     ! depth of rain on foliage [mm]
      real(r8), allocatable :: sag      (:)     ! non dimensional snow age [-]
      real(r8), allocatable :: scv      (:)     ! snow cover, water equivalent [mm]
      real(r8), allocatable :: snowdp   (:)     ! snow depth [meter]
      real(r8), allocatable :: fveg     (:)     ! fraction of vegetation cover
      real(r8), allocatable :: fsno     (:)     ! fraction of snow cover on ground
      real(r8), allocatable :: sigf     (:)     ! fraction of veg cover, excluding snow-covered veg [-]
      real(r8), allocatable :: green    (:)     ! leaf greenness
      real(r8), allocatable :: tlai     (:)     ! leaf area index
      real(r8), allocatable :: lai      (:)     ! leaf area index
      real(r8), allocatable :: laisun   (:)     ! leaf area index for sunlit leaf
      real(r8), allocatable :: laisha   (:)     ! leaf area index for shaded leaf
      real(r8), allocatable :: tsai     (:)     ! stem area index
      real(r8), allocatable :: sai      (:)     ! stem area index
      real(r8), allocatable :: coszen   (:)     ! cosine of solar zenith angle
      real(r8), allocatable :: alb  (:,:,:)     ! averaged albedo [-]
      real(r8), allocatable :: ssun (:,:,:)     ! sunlit canopy absorption for solar radiation (0-1)
      real(r8), allocatable :: ssha (:,:,:)     ! shaded canopy absorption for solar radiation (0-1)
      real(r8), allocatable :: thermk   (:)     ! canopy gap fraction for tir radiation
      real(r8), allocatable :: extkb    (:)     ! (k, g(mu)/mu) direct solar extinction coefficient
      real(r8), allocatable :: extkd    (:)     ! diffuse and scattered diffuse PAR extinction coefficient
      real(r8), allocatable :: zwt      (:)     ! the depth to water table [m]
      real(r8), allocatable :: wa       (:)     ! water storage in aquifer [mm]
      real(r8), allocatable :: wat      (:)     ! total water storage [mm]
      real(r8), allocatable :: dpond    (:)     ! depth of ponding water [mm]

      real(r8), allocatable :: t_lake(:,:)      ! lake layer teperature [K]
      real(r8), allocatable :: lake_icefrac(:,:)! lake mass fraction of lake layer that is frozen
      real(r8), allocatable :: savedtke1(:)     ! top level eddy conductivity (W/m K)

      REAL(r8), allocatable :: snw_rds     (:,:) !effective grain radius (col,lyr) [microns, m-6]
      REAL(r8), allocatable :: mss_bcpho   (:,:) !mass of hydrophobic BC in snow  (col,lyr) [kg]
      REAL(r8), allocatable :: mss_bcphi   (:,:) !mass of hydrophillic BC in snow (col,lyr) [kg]
      REAL(r8), allocatable :: mss_ocpho   (:,:) !mass of hydrophobic OC in snow  (col,lyr) [kg]
      REAL(r8), allocatable :: mss_ocphi   (:,:) !mass of hydrophillic OC in snow (col,lyr) [kg]
      REAL(r8), allocatable :: mss_dst1    (:,:) !mass of dust species 1 in snow  (col,lyr) [kg]
      REAL(r8), allocatable :: mss_dst2    (:,:) !mass of dust species 2 in snow  (col,lyr) [kg]
      REAL(r8), allocatable :: mss_dst3    (:,:) !mass of dust species 3 in snow  (col,lyr) [kg]
      REAL(r8), allocatable :: mss_dst4    (:,:) !mass of dust species 4 in snow  (col,lyr) [kg]
      REAL(r8), allocatable :: ssno    (:,:,:,:) !snow layer absorption [-]

      real(r8), allocatable :: trad     (:) ! radiative temperature of surface [K]
      real(r8), allocatable :: tref     (:) ! 2 m height air temperature [kelvin]
      real(r8), allocatable :: qref     (:) ! 2 m height air specific humidity
      real(r8), allocatable :: rst      (:) ! canopy stomatal resistance (s/m)
      real(r8), allocatable :: emis     (:) ! averaged bulk surface emissivity
      real(r8), allocatable :: z0m      (:) ! effective roughness [m]
      real(r8), allocatable :: displa   (:) ! zero displacement height [m]
      real(r8), allocatable :: zol      (:) ! dimensionless height (z/L) used in Monin-Obukhov theory
      real(r8), allocatable :: rib      (:) ! bulk Richardson number in surface layer
      real(r8), allocatable :: ustar    (:) ! u* in similarity theory [m/s]
      real(r8), allocatable :: qstar    (:) ! q* in similarity theory [kg/kg]
      real(r8), allocatable :: tstar    (:) ! t* in similarity theory [K]
      real(r8), allocatable :: fm       (:) ! integral of profile function for momentum
      real(r8), allocatable :: fh       (:) ! integral of profile function for heat
      real(r8), allocatable :: fq       (:) ! integral of profile function for moisture

      ! PUBLIC MEMBER FUNCTIONS:
      public :: allocate_TimeVariables
      public :: deallocate_TimeVariables
      public :: READ_TimeVariables
      public :: WRITE_TimeVariables
#ifdef CLMDEBUG
      public :: check_TimeVariables
#endif


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_TimeVariables
! --------------------------------------------------------------------
! Allocates memory for CLM 1d [numpatch] variables
! ------------------------------------------------------

  use precision
  USE GlobalVars
  use spmd_task
  use mod_landpatch, only : numpatch
  IMPLICIT NONE


  if (p_is_worker) then

     if (numpatch > 0) then

        allocate (z_sno      (maxsnl+1:0,      numpatch))
        allocate (dz_sno     (maxsnl+1:0,      numpatch))
        allocate (t_soisno   (maxsnl+1:nl_soil,numpatch))
        allocate (wliq_soisno(maxsnl+1:nl_soil,numpatch))
        allocate (wice_soisno(maxsnl+1:nl_soil,numpatch))
        allocate (smp        (1:nl_soil,numpatch))
        allocate (hk         (1:nl_soil,numpatch))
        allocate (h2osoi     (1:nl_soil,numpatch))
        allocate (rootr      (1:nl_soil,numpatch))
#ifdef PLANT_HYDRAULIC_STRESS
        allocate (vegwp      (1:nvegwcs,numpatch))
        allocate (gs0sun               (numpatch))
        allocate (gs0sha               (numpatch))
#endif
#ifdef OzoneStress
        allocate (o3coefv_sun          (numpatch)) ! Ozone stress factor for photosynthesis on sunlit leaf
        allocate (o3coefv_sha          (numpatch)) ! Ozone stress factor for photosynthesis on shaded leaf
        allocate (o3coefg_sun          (numpatch)) ! Ozone stress factor for stomata on sunlit leaf
        allocate (o3coefg_sha          (numpatch)) ! Ozone stress factor for stomata on shaded leaf
        allocate (lai_old              (numpatch)) ! lai in last time step
        allocate (o3uptakesun          (numpatch)) ! Ozone does, sunlit leaf (mmol O3/m^2)
        allocate (o3uptakesha          (numpatch)) ! Ozone does, shaded leaf (mmol O3/m^2)
#endif
        allocate (rstfacsun_out        (numpatch))
        allocate (rstfacsha_out        (numpatch))
        allocate (gssun_out            (numpatch))
        allocate (gssha_out            (numpatch))
        allocate (t_grnd               (numpatch))
#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
            allocate ( assim_RuBP_sun_out        (numpatch) )
            allocate ( assim_RuBP_sha_out        (numpatch) )
            allocate ( assim_Rubisco_sun_out        (numpatch) )
            allocate ( assim_Rubisco_sha_out        (numpatch) )
            allocate ( assimsun_out        (numpatch) )
            allocate ( assimsha_out        (numpatch) )
            allocate ( etrsun_out        (numpatch) )
            allocate ( etrsha_out        (numpatch) )
            allocate ( cisun_out        (numpatch) )
            allocate ( cisha_out        (numpatch) )
            allocate ( Dsun_out        (numpatch) )
            allocate ( Dsha_out        (numpatch) )
            allocate ( gammasun_out        (numpatch) )
            allocate ( gammasha_out        (numpatch) )
            allocate ( lambdasun_out        (numpatch) )
            allocate ( lambdasha_out        (numpatch) )
            allocate ( lambda_out                   (numpatch) )
#endif
#endif
        allocate (tleaf                (numpatch))
        allocate (ldew                 (numpatch))
        allocate (ldew_rain            (numpatch))
        allocate (ldew_snow            (numpatch))
        allocate (sag                  (numpatch))
        allocate (scv                  (numpatch))
        allocate (snowdp               (numpatch))
        allocate (fveg                 (numpatch))
        allocate (fsno                 (numpatch))
        allocate (sigf                 (numpatch))
        allocate (green                (numpatch))
        allocate (tlai                 (numpatch))
        allocate (lai                  (numpatch))
        allocate (laisun               (numpatch))
        allocate (laisha               (numpatch))
        allocate (tsai                 (numpatch))
        allocate (sai                  (numpatch))
        allocate (coszen               (numpatch))
        allocate (alb              (2,2,numpatch))
        allocate (ssun             (2,2,numpatch))
        allocate (ssha             (2,2,numpatch))
        allocate (thermk               (numpatch))
        allocate (extkb                (numpatch))
        allocate (extkd                (numpatch))
        allocate (zwt                  (numpatch))
        allocate (wa                   (numpatch))
        allocate (wat                  (numpatch))
        allocate (dpond                (numpatch))

        allocate (t_lake       (nl_lake,numpatch))    !new lake scheme
        allocate (lake_icefrac (nl_lake,numpatch))    !new lake scheme
        allocate (savedtke1            (numpatch))    !new lake scheme

        allocate (snw_rds           (maxsnl+1:0,numpatch))
        allocate (mss_bcpho         (maxsnl+1:0,numpatch))
        allocate (mss_bcphi         (maxsnl+1:0,numpatch))
        allocate (mss_ocpho         (maxsnl+1:0,numpatch))
        allocate (mss_ocphi         (maxsnl+1:0,numpatch))
        allocate (mss_dst1          (maxsnl+1:0,numpatch))
        allocate (mss_dst2          (maxsnl+1:0,numpatch))
        allocate (mss_dst3          (maxsnl+1:0,numpatch))
        allocate (mss_dst4          (maxsnl+1:0,numpatch))
        allocate (ssno          (2,2,maxsnl+1:1,numpatch))

        allocate (trad                 (numpatch))
        allocate (tref                 (numpatch))
        allocate (qref                 (numpatch))
        allocate (rst                  (numpatch))
        allocate (emis                 (numpatch))
        allocate (z0m                  (numpatch))
        allocate (displa               (numpatch))
        allocate (zol                  (numpatch))
        allocate (rib                  (numpatch))
        allocate (ustar                (numpatch))
        allocate (qstar                (numpatch))
        allocate (tstar                (numpatch))
        allocate (fm                   (numpatch))
        allocate (fh                   (numpatch))
        allocate (fq                   (numpatch))

     end if
  end if

#ifdef PFT_CLASSIFICATION
     CALL allocate_PFTimeVars
#endif

#ifdef PC_CLASSIFICATION
     CALL allocate_PCTimeVars
#endif

#ifdef BGC
     CALL allocate_BGCTimeVars
#endif

#ifdef LATERAL_FLOW
     CALL allocate_HydroTimeVars
#endif

#ifdef URBAN_MODEL
     CALL allocate_UrbanTimeVars
#endif

  END SUBROUTINE allocate_TimeVariables



  SUBROUTINE deallocate_TimeVariables ()

     use spmd_task
     use mod_landpatch, only : numpatch
     implicit none

     ! --------------------------------------------------
     ! Deallocates memory for CLM 1d [numpatch] variables
     ! --------------------------------------------------

     if (p_is_worker) then

        if (numpatch > 0) then

           deallocate (z_sno    )
           deallocate (dz_sno   )
           deallocate (t_soisno    )
           deallocate (wliq_soisno )
           deallocate (wice_soisno )
           deallocate (smp )
           deallocate (hk  )
           deallocate (h2osoi )
           deallocate (rootr  )
           deallocate (rstfacsun_out )
           deallocate (rstfacsha_out )
           deallocate (gssun_out )
           deallocate (gssha_out )
#ifdef WUEdiag
#ifdef PFT_CLASSIFICATION
           deallocate ( assim_RuBP_sun_out        )
           deallocate ( assim_RuBP_sha_out        )
           deallocate ( assim_Rubisco_sun_out        )
           deallocate ( assim_Rubisco_sha_out        )
           deallocate ( assimsun_out        )
           deallocate ( assimsha_out        )
           deallocate ( etrsun_out        )
           deallocate ( etrsha_out        )
           deallocate ( cisun_out        )
           deallocate ( cisha_out        )
           deallocate ( Dsun_out        )
           deallocate ( Dsha_out        )
           deallocate ( gammasun_out        )
           deallocate ( gammasha_out        )
           deallocate ( lambdasun_out        )
           deallocate ( lambdasha_out        )
           deallocate ( lambda_out                   )
#endif
#endif
#ifdef PLANT_HYDRAULIC_STRESS
           deallocate (vegwp  )
           deallocate (gs0sun )
           deallocate (gs0sha )
#endif
#ifdef OzoneStress
           deallocate (o3coefv_sun) ! Ozone stress factor for photosynthesis on sunlit leaf
           deallocate (o3coefv_sha) ! Ozone stress factor for photosynthesis on shaded leaf
           deallocate (o3coefg_sun) ! Ozone stress factor for stomata on sunlit leaf
           deallocate (o3coefg_sha) ! Ozone stress factor for stomata on shaded leaf
           deallocate (lai_old    ) ! lai in last time step
           deallocate (o3uptakesun) ! Ozone does, sunlit leaf (mmol O3/m^2)
           deallocate (o3uptakesha) ! Ozone does, shaded leaf (mmol O3/m^2)
#endif
           deallocate (t_grnd )
           deallocate (tleaf  )
           deallocate (ldew   )
           deallocate (ldew_rain)
           deallocate (ldew_snow)
           deallocate (sag    )
           deallocate (scv    )
           deallocate (snowdp )
           deallocate (fveg   )
           deallocate (fsno   )
           deallocate (sigf   )
           deallocate (green  )
           deallocate (tlai   )
           deallocate (lai    )
           deallocate (laisun )
           deallocate (laisha )
           deallocate (tsai   )
           deallocate (sai    )
           deallocate (coszen )
           deallocate (alb    )
           deallocate (ssun   )
           deallocate (ssha   )
           deallocate (thermk )
           deallocate (extkb  )
           deallocate (extkd  )
           deallocate (zwt    )
           deallocate (wa     )
           deallocate (wat    )
           deallocate (dpond  )

           deallocate (t_lake )      ! new lake scheme
           deallocate (lake_icefrac) ! new lake scheme
           deallocate (savedtke1)    ! new lake scheme

           deallocate (snw_rds  )
           deallocate (mss_bcpho)
           deallocate (mss_bcphi)
           deallocate (mss_ocpho)
           deallocate (mss_ocphi)
           deallocate (mss_dst1 )
           deallocate (mss_dst2 )
           deallocate (mss_dst3 )
           deallocate (mss_dst4 )
           deallocate (ssno     )

           deallocate (trad   )
           deallocate (tref   )
           deallocate (qref   )
           deallocate (rst    )
           deallocate (emis   )
           deallocate (z0m    )
           deallocate (displa )
           deallocate (zol    )
           deallocate (rib    )
           deallocate (ustar  )
           deallocate (qstar  )
           deallocate (tstar  )
           deallocate (fm     )
           deallocate (fh     )
           deallocate (fq     )

        end if
     end if

#if (defined PFT_CLASSIFICATION)
     CALL deallocate_PFTimeVars
#endif

#if (defined PC_CLASSIFICATION)
     CALL deallocate_PCTimeVars
#endif

#if (defined BGC)
     CALL deallocate_BGCTimeVars
#endif

#ifdef LATERAL_FLOW
     CALL deallocate_HydroTimeVars
#endif

#if (defined URBAN_MODEL)
     CALL deallocate_UrbanTimeVars
#endif

  END SUBROUTINE deallocate_TimeVariables


  !---------------------------------------
  function save_to_restart (idate, deltim, itstamp, ptstamp) result(rwrite)

     use mod_namelist
     implicit none

     logical :: rwrite

     integer,  intent(in) :: idate(3)
     real(r8), intent(in) :: deltim
     type(timestamp), intent(in) :: itstamp, ptstamp


     ! added by yuan, 08/31/2014
     select case (trim(DEF_WRST_FREQ))
     case ('TIMESTEP')
        rwrite = .true.
     case ('HOURLY')
        rwrite = isendofhour (idate, deltim)
     case ('DAILY')
        rwrite = isendofday(idate, deltim)
     case ('MONTHLY')
        rwrite = isendofmonth(idate, deltim)
     case ('YEARLY')
        rwrite = isendofyear(idate, deltim)
     end select

     if (rwrite) then
        rwrite = (ptstamp < itstamp)
     end if

  end function save_to_restart

  !---------------------------------------
  SUBROUTINE WRITE_TimeVariables (idate, site, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use mod_namelist, only : DEF_REST_COMPRESS_LEVEL
     USE mod_landpatch
     use ncio_vector
     USE GlobalVars
     IMPLICIT NONE

     integer, INTENT(in) :: idate(3)
     character(LEN=*), intent(in) :: site
     character(LEN=*), intent(in) :: dir_restart

     ! Local variables
     character(LEN=256) :: file_restart
     character(len=14)  :: cdate
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL

     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_'//trim(cdate)//'.nc'

     call ncio_create_file_vector (file_restart, landpatch)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')

     CALL ncio_define_dimension_vector (file_restart, landpatch, 'snow',     -maxsnl       )
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'snowp1',   -maxsnl+1     )
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'soilsnow', nl_soil-maxsnl)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'soil',     nl_soil)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'lake',     nl_lake)

#ifdef PLANT_HYDRAULIC_STRESS
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'vegnodes', nvegwcs)
#endif

     CALL ncio_define_dimension_vector (file_restart, landpatch, 'band', 2)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'rtyp', 2)

     ! Time-varying state variables which reaquired by restart run
     call ncio_write_vector (file_restart, 'z_sno   '   , 'snow', -maxsnl, 'patch', landpatch, z_sno , compress) !  node depth [m]
     call ncio_write_vector (file_restart, 'dz_sno  '   , 'snow', -maxsnl, 'patch', landpatch, dz_sno, compress) !  interface depth [m]
     call ncio_write_vector (file_restart, 't_soisno'   , 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, t_soisno   , compress) !  soil temperature [K]
     call ncio_write_vector (file_restart, 'wliq_soisno', 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, wliq_soisno, compress) !  liquid water in layers [kg/m2]
     call ncio_write_vector (file_restart, 'wice_soisno', 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, wice_soisno, compress) !  ice lens in layers [kg/m2]
     call ncio_write_vector (file_restart, 'smp',         'soil', nl_soil, 'patch', landpatch, smp, compress) !  soil matrix potential [mm]
     call ncio_write_vector (file_restart, 'hk',          'soil', nl_soil, 'patch', landpatch, hk, compress) !  hydraulic conductivity [mm h2o/s]
#ifdef PLANT_HYDRAULIC_STRESS
     call ncio_write_vector (file_restart, 'vegwp',   'vegnodes', nvegwcs, 'patch', landpatch, vegwp, compress) !  vegetation water potential [mm]
     call ncio_write_vector (file_restart, 'gs0sun',    'patch', landpatch, gs0sun, compress) !  working copy of sunlit stomata conductance
     call ncio_write_vector (file_restart, 'gs0sha',    'patch', landpatch, gs0sha, compress) !  working copy of shalit stomata conductance
#endif
#ifdef OzoneStress
     call ncio_write_vector (file_restart, 'lai_old    ', 'patch', landpatch, lai_old    , compress)
     call ncio_write_vector (file_restart, 'o3uptakesun', 'patch', landpatch, o3uptakesun, compress)
     call ncio_write_vector (file_restart, 'o3uptakesha', 'patch', landpatch, o3uptakesha, compress)
#endif
     call ncio_write_vector (file_restart, 't_grnd  '   , 'patch', landpatch, t_grnd    , compress) !  ground surface temperature [K]
     call ncio_write_vector (file_restart, 'tleaf   '   , 'patch', landpatch, tleaf     , compress) !  leaf temperature [K]
     call ncio_write_vector (file_restart, 'ldew    '   , 'patch', landpatch, ldew      , compress) !  depth of water on foliage [mm]
     call ncio_write_vector (file_restart, 'ldew_rain'  , 'patch', landpatch, ldew_rain , compress) !  depth of water on foliage [mm]
     call ncio_write_vector (file_restart, 'ldew_snow'  , 'patch', landpatch, ldew_snow , compress) !  depth of water on foliage [mm]
     call ncio_write_vector (file_restart, 'sag     '   , 'patch', landpatch, sag       , compress) !  non dimensional snow age [-]
     call ncio_write_vector (file_restart, 'scv     '   , 'patch', landpatch, scv       , compress) !  snow cover, water equivalent [mm]
     call ncio_write_vector (file_restart, 'snowdp  '   , 'patch', landpatch, snowdp    , compress) !  snow depth [meter]
     call ncio_write_vector (file_restart, 'fveg    '   , 'patch', landpatch, fveg      , compress) !  fraction of vegetation cover
     call ncio_write_vector (file_restart, 'fsno    '   , 'patch', landpatch, fsno      , compress) !  fraction of snow cover on ground
     call ncio_write_vector (file_restart, 'sigf    '   , 'patch', landpatch, sigf      , compress) !  fraction of veg cover, excluding snow-covered veg [-]
     call ncio_write_vector (file_restart, 'green   '   , 'patch', landpatch, green     , compress) !  leaf greenness
     call ncio_write_vector (file_restart, 'lai     '   , 'patch', landpatch, lai       , compress) !  leaf area index
     call ncio_write_vector (file_restart, 'tlai    '   , 'patch', landpatch, tlai      , compress) !  leaf area index
     call ncio_write_vector (file_restart, 'sai     '   , 'patch', landpatch, sai       , compress) !  stem area index
     call ncio_write_vector (file_restart, 'tsai    '   , 'patch', landpatch, tsai      , compress) !  stem area index
     call ncio_write_vector (file_restart, 'coszen  '   , 'patch', landpatch, coszen    , compress) !  cosine of solar zenith angle
     call ncio_write_vector (file_restart, 'alb     '   , 'band', 2, 'rtyp', 2, 'patch', landpatch, alb , compress) !  averaged albedo [-]
     call ncio_write_vector (file_restart, 'ssun    '   , 'band', 2, 'rtyp', 2, 'patch', landpatch, ssun, compress) !  sunlit canopy absorption for solar radiation (0-1)
     call ncio_write_vector (file_restart, 'ssha    '   , 'band', 2, 'rtyp', 2, 'patch', landpatch, ssha, compress) !  shaded canopy absorption for solar radiation (0-1)
     call ncio_write_vector (file_restart, 'thermk  '   , 'patch', landpatch, thermk    , compress) !  canopy gap fraction for tir radiation
     call ncio_write_vector (file_restart, 'extkb   '   , 'patch', landpatch, extkb     , compress) !  (k, g(mu)/mu) direct solar extinction coefficient
     call ncio_write_vector (file_restart, 'extkd   '   , 'patch', landpatch, extkd     , compress) !  diffuse and scattered diffuse PAR extinction coefficient
     call ncio_write_vector (file_restart, 'zwt     '   , 'patch', landpatch, zwt       , compress) !  the depth to water table [m]
     call ncio_write_vector (file_restart, 'wa      '   , 'patch', landpatch, wa        , compress) !  water storage in aquifer [mm]
     call ncio_write_vector (file_restart, 'dpond   '   , 'patch', landpatch, dpond     , compress) ! depth of ponding water

     call ncio_write_vector (file_restart, 't_lake  '   , 'lake', nl_lake, 'patch', landpatch, t_lake      , compress) !
     call ncio_write_vector (file_restart, 'lake_icefrc', 'lake', nl_lake, 'patch', landpatch, lake_icefrac, compress) !
     call ncio_write_vector (file_restart, 'savedtke1  ', 'patch', landpatch, savedtke1   , compress) !
     call ncio_write_vector (file_restart, 'snw_rds  ', 'snow', -maxsnl, 'patch', landpatch, snw_rds  , compress)
     call ncio_write_vector (file_restart, 'mss_bcpho', 'snow', -maxsnl, 'patch', landpatch, mss_bcpho, compress)
     call ncio_write_vector (file_restart, 'mss_bcphi', 'snow', -maxsnl, 'patch', landpatch, mss_bcphi, compress)
     call ncio_write_vector (file_restart, 'mss_ocpho', 'snow', -maxsnl, 'patch', landpatch, mss_ocpho, compress)
     call ncio_write_vector (file_restart, 'mss_ocphi', 'snow', -maxsnl, 'patch', landpatch, mss_ocphi, compress)
     call ncio_write_vector (file_restart, 'mss_dst1 ', 'snow', -maxsnl, 'patch', landpatch, mss_dst1 , compress)
     call ncio_write_vector (file_restart, 'mss_dst2 ', 'snow', -maxsnl, 'patch', landpatch, mss_dst2 , compress)
     call ncio_write_vector (file_restart, 'mss_dst3 ', 'snow', -maxsnl, 'patch', landpatch, mss_dst3 , compress)
     call ncio_write_vector (file_restart, 'mss_dst4 ', 'snow', -maxsnl, 'patch', landpatch, mss_dst4 , compress)
     call ncio_write_vector (file_restart, 'ssno', 'band', 2, 'rtyp', 2, 'snowp1', -maxsnl+1, 'patch', landpatch, ssno, compress)

     ! Additional va_vectorriables required by reginal model (such as WRF ) RSM)
     call ncio_write_vector (file_restart, 'trad ', 'patch', landpatch, trad , compress) !     radiative temperature of surface [K]
     call ncio_write_vector (file_restart, 'tref ', 'patch', landpatch, tref , compress) !     2 m height air temperature [kelvin]
     call ncio_write_vector (file_restart, 'qref ', 'patch', landpatch, qref , compress) !     2 m height air specific humidity
     call ncio_write_vector (file_restart, 'rst  ', 'patch', landpatch, rst  , compress) !     canopy stomatal resistance (s/m)
     call ncio_write_vector (file_restart, 'emis ', 'patch', landpatch, emis , compress) !     averaged bulk surface emissivity
     call ncio_write_vector (file_restart, 'z0m  ', 'patch', landpatch, z0m  , compress) !     effective roughness [m]
     call ncio_write_vector (file_restart, 'zol  ', 'patch', landpatch, zol  , compress) !     dimensionless height (z/L) used in Monin-Obukhov theory
     call ncio_write_vector (file_restart, 'rib  ', 'patch', landpatch, rib  , compress) !     bulk Richardson number in surface layer
     call ncio_write_vector (file_restart, 'ustar', 'patch', landpatch, ustar, compress) !     u* in similarity theory [m/s]
     call ncio_write_vector (file_restart, 'qstar', 'patch', landpatch, qstar, compress) !     q* in similarity theory [kg/kg]
     call ncio_write_vector (file_restart, 'tstar', 'patch', landpatch, tstar, compress) !     t* in similarity theory [K]
     call ncio_write_vector (file_restart, 'fm   ', 'patch', landpatch, fm   , compress) !     integral of profile function for momentum
     call ncio_write_vector (file_restart, 'fh   ', 'patch', landpatch, fh   , compress) !     integral of profile function for heat
     call ncio_write_vector (file_restart, 'fq   ', 'patch', landpatch, fq   , compress) !     integral of profile function for moisture

#if (defined PFT_CLASSIFICATION)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pft_'//trim(cdate)//'.nc'
     CALL WRITE_PFTimeVars (file_restart)
#endif

#if (defined PC_CLASSIFICATION)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pc_'//trim(cdate)//'.nc'
     CALL WRITE_PCTimeVars (file_restart)
#endif

#if (defined BGC)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_bgc_'//trim(cdate)//'.nc'
     CALL WRITE_BGCTimeVars (file_restart)
#endif

#if (defined LATERAL_FLOW)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_basin_'//trim(cdate)//'.nc'
     CALL WRITE_HydroTimeVars (file_restart)
#endif

#if (defined URBAN_MODEL)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_urban_'//trim(cdate)//'.nc'
     CALL WRITE_UrbanTimeVars (file_restart)
#endif
  end subroutine WRITE_TimeVariables

  !---------------------------------------
  SUBROUTINE READ_TimeVariables (idate, site, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use mod_namelist
     use spmd_task
     use ncio_vector
#ifdef CLMDEBUG
     USE mod_colm_debug
#endif
     USE mod_landpatch
     USE GlobalVars

     IMPLICIT NONE

     integer, INTENT(in) :: idate(3)
     character(LEN=*), intent(in) :: site
     character(LEN=*), intent(in) :: dir_restart

     ! Local variables
     character(LEN=256) :: file_restart
     character(len=14)  :: cdate

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif

     if (p_is_master) then
        write(*,'(/,A26)') 'Loading Time Variables ...'
     end if

     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
     file_restart = trim(dir_restart) // '/' // trim(site) //'_restart_'//trim(cdate)//'.nc'

     ! Time-varying state variables which reaquired by restart run
     call ncio_read_vector (file_restart, 'z_sno   '   , -maxsnl, landpatch, z_sno ) !  node depth [m]
     call ncio_read_vector (file_restart, 'dz_sno  '   , -maxsnl, landpatch, dz_sno) !  interface depth [m]
     call ncio_read_vector (file_restart, 't_soisno'   , nl_soil-maxsnl, landpatch, t_soisno   ) !  soil temperature [K]
     call ncio_read_vector (file_restart, 'wliq_soisno', nl_soil-maxsnl, landpatch, wliq_soisno) !  liquid water in layers [kg/m2]
     call ncio_read_vector (file_restart, 'wice_soisno', nl_soil-maxsnl, landpatch, wice_soisno) !  ice lens in layers [kg/m2]
     call ncio_read_vector (file_restart, 'smp',         nl_soil,        landpatch, smp        ) !  soil matrix potential [mm]
     call ncio_read_vector (file_restart, 'hk',          nl_soil,        landpatch, hk         ) !  hydraulic conductivity [mm h2o/s]
#ifdef PLANT_HYDRAULIC_STRESS
     call ncio_read_vector (file_restart, 'vegwp',       nvegwcs,        landpatch, vegwp      ) !  vegetation water potential [mm]
     call ncio_read_vector (file_restart, 'gs0sun  ',    landpatch, gs0sun     ) !  working copy of sunlit stomata conductance
     call ncio_read_vector (file_restart, 'gs0sha  ',    landpatch, gs0sha     ) !  working copy of shalit stomata conductance
#endif
#ifdef OzoneStress
     call ncio_read_vector (file_restart, 'lai_old    ', landpatch, lai_old    )
     call ncio_read_vector (file_restart, 'o3uptakesun', landpatch, o3uptakesun)
     call ncio_read_vector (file_restart, 'o3uptakesha', landpatch, o3uptakesha)
#endif
     call ncio_read_vector (file_restart, 't_grnd  '   , landpatch, t_grnd     ) !  ground surface temperature [K]
     call ncio_read_vector (file_restart, 'tleaf   '   , landpatch, tleaf      ) !  leaf temperature [K]
     call ncio_read_vector (file_restart, 'ldew    '   , landpatch, ldew       ) !  depth of water on foliage [mm]
     call ncio_read_vector (file_restart, 'ldew_rain    '   , landpatch, ldew_rain       ) !  depth of rain on foliage [mm]
     call ncio_read_vector (file_restart, 'ldew_snow    '   , landpatch, ldew_snow       ) !  depth of snow on foliage [mm]
     call ncio_read_vector (file_restart, 'sag     '   , landpatch, sag        ) !  non dimensional snow age [-]
     call ncio_read_vector (file_restart, 'scv     '   , landpatch, scv        ) !  snow cover, water equivalent [mm]
     call ncio_read_vector (file_restart, 'snowdp  '   , landpatch, snowdp     ) !  snow depth [meter]
     call ncio_read_vector (file_restart, 'fveg    '   , landpatch, fveg       ) !  fraction of vegetation cover
     call ncio_read_vector (file_restart, 'fsno    '   , landpatch, fsno       ) !  fraction of snow cover on ground
     call ncio_read_vector (file_restart, 'sigf    '   , landpatch, sigf       ) !  fraction of veg cover, excluding snow-covered veg [-]
     call ncio_read_vector (file_restart, 'green   '   , landpatch, green      ) !  leaf greenness
     call ncio_read_vector (file_restart, 'lai     '   , landpatch, lai        ) !  leaf area index
     call ncio_read_vector (file_restart, 'tlai    '   , landpatch, tlai       ) !  leaf area index
     call ncio_read_vector (file_restart, 'sai     '   , landpatch, sai        ) !  stem area index
     call ncio_read_vector (file_restart, 'tsai    '   , landpatch, tsai       ) !  stem area index
     call ncio_read_vector (file_restart, 'coszen  '   , landpatch, coszen     ) !  cosine of solar zenith angle
     call ncio_read_vector (file_restart, 'alb     '   , 2, 2, landpatch, alb  ) !  averaged albedo [-]
     call ncio_read_vector (file_restart, 'ssun    '   , 2, 2, landpatch, ssun ) !  sunlit canopy absorption for solar radiation (0-1)
     call ncio_read_vector (file_restart, 'ssha    '   , 2, 2, landpatch, ssha ) !  shaded canopy absorption for solar radiation (0-1)
     call ncio_read_vector (file_restart, 'thermk  '   , landpatch, thermk     ) !  canopy gap fraction for tir radiation
     call ncio_read_vector (file_restart, 'extkb   '   , landpatch, extkb      ) !  (k, g(mu)/mu) direct solar extinction coefficient
     call ncio_read_vector (file_restart, 'extkd   '   , landpatch, extkd      ) !  diffuse and scattered diffuse PAR extinction coefficient
     call ncio_read_vector (file_restart, 'zwt     '   , landpatch, zwt        ) !  the depth to water table [m]
     call ncio_read_vector (file_restart, 'wa      '   , landpatch, wa         ) !  water storage in aquifer [mm]
     call ncio_read_vector (file_restart, 'dpond   '   , landpatch, dpond      ) ! depth of ponding water

     call ncio_read_vector (file_restart, 't_lake  '   , nl_lake, landpatch, t_lake      ) !
     call ncio_read_vector (file_restart, 'lake_icefrc', nl_lake, landpatch, lake_icefrac) !
     call ncio_read_vector (file_restart, 'savedtke1', landpatch, savedtke1) !

     call ncio_read_vector (file_restart, 'snw_rds  ', -maxsnl, landpatch, snw_rds  ) !
     call ncio_read_vector (file_restart, 'mss_bcpho', -maxsnl, landpatch, mss_bcpho) !
     call ncio_read_vector (file_restart, 'mss_bcphi', -maxsnl, landpatch, mss_bcphi) !
     call ncio_read_vector (file_restart, 'mss_ocpho', -maxsnl, landpatch, mss_ocpho) !
     call ncio_read_vector (file_restart, 'mss_ocphi', -maxsnl, landpatch, mss_ocphi) !
     call ncio_read_vector (file_restart, 'mss_dst1 ', -maxsnl, landpatch, mss_dst1 ) !
     call ncio_read_vector (file_restart, 'mss_dst2 ', -maxsnl, landpatch, mss_dst2 ) !
     call ncio_read_vector (file_restart, 'mss_dst3 ', -maxsnl, landpatch, mss_dst3 ) !
     call ncio_read_vector (file_restart, 'mss_dst4 ', -maxsnl, landpatch, mss_dst4 ) !
     call ncio_read_vector (file_restart, 'ssno', 2,2, -maxsnl+1, landpatch, ssno) !

     ! Additional variables required by reginal model (such as WRF ) RSM)
     call ncio_read_vector (file_restart, 'trad ', landpatch, trad ) ! radiative temperature of surface [K]
     call ncio_read_vector (file_restart, 'tref ', landpatch, tref ) ! 2 m height air temperature [kelvin]
     call ncio_read_vector (file_restart, 'qref ', landpatch, qref ) ! 2 m height air specific humidity
     call ncio_read_vector (file_restart, 'rst  ', landpatch, rst  ) ! canopy stomatal resistance (s/m)
     call ncio_read_vector (file_restart, 'emis ', landpatch, emis ) ! averaged bulk surface emissivity
     call ncio_read_vector (file_restart, 'z0m  ', landpatch, z0m  ) ! effective roughness [m]
     call ncio_read_vector (file_restart, 'zol  ', landpatch, zol  ) ! dimensionless height (z/L) used in Monin-Obukhov theory
     call ncio_read_vector (file_restart, 'rib  ', landpatch, rib  ) ! bulk Richardson number in surface layer
     call ncio_read_vector (file_restart, 'ustar', landpatch, ustar) ! u* in similarity theory [m/s]
     call ncio_read_vector (file_restart, 'qstar', landpatch, qstar) ! q* in similarity theory [kg/kg]
     call ncio_read_vector (file_restart, 'tstar', landpatch, tstar) ! t* in similarity theory [K]
     call ncio_read_vector (file_restart, 'fm   ', landpatch, fm   ) ! integral of profile function for momentum
     call ncio_read_vector (file_restart, 'fh   ', landpatch, fh   ) ! integral of profile function for heat
     call ncio_read_vector (file_restart, 'fq   ', landpatch, fq   ) ! integral of profile function for moisture

#if (defined PFT_CLASSIFICATION)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pft_'//trim(cdate)//'.nc'
     CALL READ_PFTimeVars (file_restart)
#endif

#if (defined PC_CLASSIFICATION)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pc_'//trim(cdate)//'.nc'
     CALL READ_PCTimeVars (file_restart)
#endif

#if (defined BGC)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_bgc_'//trim(cdate)//'.nc'
     CALL READ_BGCTimeVars (file_restart)
#endif

#if (defined LATERAL_FLOW)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_basin_'//trim(cdate)//'.nc'
     CALL READ_HydroTimeVars (file_restart)
#endif

#if (defined URBAN_MODEL)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_urban_'//trim(cdate)//'.nc'
     CALL READ_UrbanTimeVars (file_restart)
#endif

#ifdef CLMDEBUG
     call check_TimeVariables
#endif

     if (p_is_master) then
        write(*,*) 'Loading Time Variables done.'
     end if

  end subroutine READ_TimeVariables

  !---------------------------------------
#ifdef CLMDEBUG
  SUBROUTINE check_TimeVariables ()

     use spmd_task
     use mod_colm_debug

     IMPLICIT NONE

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif
     if (p_is_master) then
        write(*,'(/,A27)') 'Checking Time Variables ...'
     end if

     call check_vector_data ('z_sno       ', z_sno )      !  node depth [m]
     call check_vector_data ('dz_sno      ', dz_sno)      !  interface depth [m]
     call check_vector_data ('t_soisno    ', t_soisno   ) !  soil temperature [K]
     call check_vector_data ('wliq_soisno ', wliq_soisno) !  liquid water in layers [kg/m2]
     call check_vector_data ('wice_soisno ', wice_soisno) !  ice lens in layers [kg/m2]
     call check_vector_data ('smp         ', smp        ) !  soil matrix potential [mm]
     call check_vector_data ('hk          ', hk         ) !  hydraulic conductivity [mm h2o/s]
#ifdef PLANT_HYDRAULIC_STRESS
     call check_vector_data ('vegwp       ', vegwp      ) !  vegetation water potential [mm]
     call check_vector_data ('gs0sun      ', gs0sun     ) !  working copy of sunlit stomata conductance
     call check_vector_data ('gs0sha      ', gs0sha     ) !  working copy of shalit stomata conductance
#endif
#ifdef OzoneStress
     call check_vector_data ('o3coefv_sun', o3coefv_sun)
     call check_vector_data ('o3coefv_sha', o3coefv_sha)
     call check_vector_data ('o3coefg_sun', o3coefg_sun)
     call check_vector_data ('o3coefg_sha', o3coefg_sha)
     call check_vector_data ('lai_old    ', lai_old    )
     call check_vector_data ('o3uptakesun', o3uptakesun)
     call check_vector_data ('o3uptakesha', o3uptakesha)
#endif
     call check_vector_data ('t_grnd      ', t_grnd     ) !  ground surface temperature [K]
     call check_vector_data ('tleaf       ', tleaf      ) !  leaf temperature [K]
     call check_vector_data ('ldew        ', ldew       ) !  depth of water on foliage [mm]
     call check_vector_data ('ldew_rain   ', ldew_rain       ) !  depth of rain on foliage [mm]
     call check_vector_data ('ldew_snow   ', ldew_snow       ) !  depth of snow on foliage [mm]
     call check_vector_data ('sag         ', sag        ) !  non dimensional snow age [-]
     call check_vector_data ('scv         ', scv        ) !  snow cover, water equivalent [mm]
     call check_vector_data ('snowdp      ', snowdp     ) !  snow depth [meter]
     call check_vector_data ('fveg        ', fveg       ) !  fraction of vegetation cover
     call check_vector_data ('fsno        ', fsno       ) !  fraction of snow cover on ground
     call check_vector_data ('sigf        ', sigf       ) !  fraction of veg cover, excluding snow-covered veg [-]
     call check_vector_data ('green       ', green      ) !  leaf greenness
     call check_vector_data ('lai         ', lai        ) !  leaf area index
     call check_vector_data ('tlai        ', tlai       ) !  leaf area index
     call check_vector_data ('sai         ', sai        ) !  stem area index
     call check_vector_data ('tsai        ', tsai       ) !  stem area index
     call check_vector_data ('coszen      ', coszen     ) !  cosine of solar zenith angle
     call check_vector_data ('alb         ', alb  ) !  averaged albedo [-]
     call check_vector_data ('ssun        ', ssun ) !  sunlit canopy absorption for solar radiation (0-1)
     call check_vector_data ('ssha        ', ssha ) !  shaded canopy absorption for solar radiation (0-1)
     call check_vector_data ('thermk      ', thermk     ) !  canopy gap fraction for tir radiation
     call check_vector_data ('extkb       ', extkb      ) !  (k, g(mu)/mu) direct solar extinction coefficient
     call check_vector_data ('extkd       ', extkd      ) !  diffuse and scattered diffuse PAR extinction coefficient
     call check_vector_data ('zwt         ', zwt        ) !  the depth to water table [m]
     call check_vector_data ('wa          ', wa         ) !  water storage in aquifer [mm]
     call check_vector_data ('dpond       ', dpond      ) !  depth of ponding water

     call check_vector_data ('t_lake      ', t_lake      ) !
     call check_vector_data ('lake_icefrc ', lake_icefrac) !
     call check_vector_data ('savedtke1   ', savedtke1   ) !

#if (defined PFT_CLASSIFICATION)
     CALL check_PFTimeVars
#endif

#if (defined PC_CLASSIFICATION)
     CALL check_PCTimeVars
#endif

#if (defined BGC)
     CALL check_BGCTimeVars
#endif

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif

  end subroutine check_TimeVariables
#endif


END MODULE MOD_TimeVariables
! ------ EOP --------------
