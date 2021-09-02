#include <define.h> 

MODULE MOD_TimeVariables
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

  USE precision
  USE timemanager
  IMPLICIT NONE
  SAVE

! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run
  REAL(r8), allocatable :: z_sno       (:,:) !node depth [m]
  REAL(r8), allocatable :: dz_sno      (:,:) !interface depth [m]
  REAL(r8), allocatable :: t_soisno    (:,:) !soil temperature [K]
  REAL(r8), allocatable :: wliq_soisno (:,:) !liquid water in layers [kg/m2]
  REAL(r8), allocatable :: wice_soisno (:,:) !ice lens in layers [kg/m2]
  REAL(r8), allocatable :: h2osoi      (:,:) !volumetric soil water in layers [m3/m3]
  REAL(r8), allocatable :: rstfac        (:) !factor of soil water stress 
  REAL(r8), allocatable :: t_grnd        (:) !ground surface temperature [K]

  REAL(r8), allocatable :: tleaf         (:) !leaf temperature [K]
  REAL(r8), allocatable :: ldew          (:) !depth of water on foliage [mm]
  REAL(r8), allocatable :: sag           (:) !non dimensional snow age [-]
  REAL(r8), allocatable :: scv           (:) !snow cover, water equivalent [mm]
  REAL(r8), allocatable :: snowdp        (:) !snow depth [meter]
  REAL(r8), allocatable :: fveg          (:) !fraction of vegetation cover
  REAL(r8), allocatable :: fsno          (:) !fraction of snow cover on ground
  REAL(r8), allocatable :: sigf          (:) !fraction of veg cover, excluding snow-covered veg [-]
  REAL(r8), allocatable :: green         (:) !leaf greenness
  REAL(r8), allocatable :: tlai          (:) !leaf area index
  REAL(r8), allocatable :: lai           (:) !leaf area index
  REAL(r8), allocatable :: laisun        (:) !leaf area index
  REAL(r8), allocatable :: laisha        (:) !leaf area index
  REAL(r8), allocatable :: tsai          (:) !stem area index
  REAL(r8), allocatable :: sai           (:) !stem area index
  REAL(r8), allocatable :: coszen        (:) !cosine of solar zenith angle
  REAL(r8), allocatable :: alb       (:,:,:) !averaged albedo [-]
  REAL(r8), allocatable :: ssun      (:,:,:) !sunlit canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: ssha      (:,:,:) !shaded canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: thermk        (:) !canopy gap fraction for tir radiation
  REAL(r8), allocatable :: extkb         (:) !(k, g(mu)/mu) direct solar extinction coefficient
  REAL(r8), allocatable :: extkd         (:) !diffuse and scattered diffuse PAR extinction coefficient
  REAL(r8), allocatable :: zwt           (:) !the depth to water table [m]
  REAL(r8), allocatable :: wa            (:) !water storage in aquifer [mm]
  REAL(r8), allocatable :: wat           (:) !total water storage [mm]

  REAL(r8), allocatable :: t_lake      (:,:) !lake layer teperature [K]
  REAL(r8), allocatable :: lake_icefrac(:,:) !lake mass fraction of lake layer that is frozen

  REAL(r8), allocatable :: trad          (:) !radiative temperature of surface [K]
  REAL(r8), allocatable :: tref          (:) !2 m height air temperature [kelvin]
  REAL(r8), allocatable :: qref          (:) !2 m height air specific humidity
  REAL(r8), allocatable :: rst           (:) !canopy stomatal resistance (s/m)
  REAL(r8), allocatable :: emis          (:) !averaged bulk surface emissivity
  REAL(r8), allocatable :: z0m           (:) !effective roughness [m]
  REAL(r8), allocatable :: displa        (:) !zero displacement height [m]
  REAL(r8), allocatable :: zol           (:) !dimensionless height (z/L) used in Monin-Obukhov theory
  REAL(r8), allocatable :: rib           (:) !bulk Richardson number in surface layer
  REAL(r8), allocatable :: ustar         (:) !u* in similarity theory [m/s]
  REAL(r8), allocatable :: qstar         (:) !q* in similarity theory [kg/kg]
  REAL(r8), allocatable :: tstar         (:) !t* in similarity theory [K]
  REAL(r8), allocatable :: fm            (:) !integral of profile function for momentum
  REAL(r8), allocatable :: fh            (:) !integral of profile function for heat
  REAL(r8), allocatable :: fq            (:) !integral of profile function for moisture

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_TimeVariables
  PUBLIC :: READ_TimeVariables
  PUBLIC :: WRITE_TimeVariables
  PUBLIC :: deallocate_TimeVariables

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_TimeVariables 
! --------------------------------------------------------------------
! Allocates memory for CLM 1d [numpatch] variables
! ------------------------------------------------------

     USE precision
     USE GlobalVars
     USE MOD_PFTimeVars
     USE MOD_PCTimeVars
     IMPLICIT NONE

     allocate (z_sno             (maxsnl+1:0,numpatch))
     allocate (dz_sno            (maxsnl+1:0,numpatch))
     allocate (t_soisno    (maxsnl+1:nl_soil,numpatch))
     allocate (wliq_soisno (maxsnl+1:nl_soil,numpatch))
     allocate (wice_soisno (maxsnl+1:nl_soil,numpatch))
     allocate (h2osoi             (1:nl_soil,numpatch))
     allocate (rstfac                       (numpatch))
     allocate (t_grnd                       (numpatch))
     allocate (tleaf                        (numpatch))
     allocate (ldew                         (numpatch))
     allocate (sag                          (numpatch))
     allocate (scv                          (numpatch))
     allocate (snowdp                       (numpatch))
     allocate (fveg                         (numpatch))
     allocate (fsno                         (numpatch))
     allocate (sigf                         (numpatch))
     allocate (green                        (numpatch))
     allocate (tlai                         (numpatch))
     allocate (lai                          (numpatch))
     allocate (laisun                       (numpatch))
     allocate (laisha                       (numpatch))
     allocate (tsai                         (numpatch))
     allocate (sai                          (numpatch))
     allocate (coszen                       (numpatch))
     allocate (alb                      (2,2,numpatch))
     allocate (ssun                     (2,2,numpatch))
     allocate (ssha                     (2,2,numpatch))
     allocate (thermk                       (numpatch))
     allocate (extkb                        (numpatch))
     allocate (extkd                        (numpatch))
     allocate (zwt                          (numpatch))
     allocate (wa                           (numpatch))
     allocate (wat                          (numpatch))

     allocate (t_lake               (nl_lake,numpatch)) 
     allocate (lake_icefrac         (nl_lake,numpatch))

     allocate (trad                         (numpatch))
     allocate (tref                         (numpatch))
     allocate (qref                         (numpatch))
     allocate (rst                          (numpatch))
     allocate (emis                         (numpatch))
     allocate (z0m                          (numpatch))
     allocate (displa                       (numpatch))
     allocate (zol                          (numpatch))
     allocate (rib                          (numpatch))
     allocate (ustar                        (numpatch))
     allocate (qstar                        (numpatch))
     allocate (tstar                        (numpatch))
     allocate (fm                           (numpatch))
     allocate (fh                           (numpatch))
     allocate (fq                           (numpatch))

#ifdef PFT_CLASSIFICATION
     CALL allocate_PFTimeVars
#endif

#ifdef PC_CLASSIFICATION
     CALL allocate_PCTimeVars
#endif
 
  END SUBROUTINE allocate_TimeVariables


  SUBROUTINE READ_TimeVariables (idate,dir_restart_hist,casename)
! --------------------------------------------------------------------
! Read the model variables for restart run [histTimeVar]
! ...............................................................

     USE precision
     USE MOD_PFTimeVars
     USE MOD_PCTimeVars
     IMPLICIT NONE

     CHARACTER(LEN=255), intent(in) :: dir_restart_hist
     CHARACTER(LEN=256), intent(in) :: casename
     INTEGER, intent(in) :: idate(3)     !calendar (year, julian day, seconds)

     INTEGER :: lhistTimeVar             !logical unit number of restart time-varying file
     INTEGER :: id(3)                    !calendar (year, julian day, seconds)
     CHARACTER(LEN=255) :: cdate         !character for date
     CHARACTER(LEN=256) :: fhistTimeVar  !file name of time-varying file

     ! the model variables for restart run
     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1),idate(2),idate(3)

     lhistTimeVar = 100
     fhistTimeVar = trim(dir_restart_hist)//trim(casename)//'-'//'rstTimeVar'//'-'//trim(cdate)
     print*,trim(fhistTimeVar)
     open(unit=lhistTimeVar,file=trim(fhistTimeVar),status='unknown',&
                            form='unformatted',action='read')

     read (lhistTimeVar) id, & !
         ! Time-varying state variables which reaquired by restart run
           z_sno,           &! node depth [m]
           dz_sno,          &! interface depth [m]
           t_soisno,        &! soil temperature [K]
           wliq_soisno,     &! liquid water in layers [kg/m2]
           wice_soisno,     &! ice lens in layers [kg/m2]
           t_grnd,          &! ground surface temperature [K]
           tleaf,           &! leaf temperature [K]
           ldew,            &! depth of water on foliage [mm]
           sag,             &! non dimensional snow age [-]
           scv,             &! snow cover, water equivalent [mm]
           snowdp,          &! snow depth [meter]
           fveg,            &! fraction of vegetation cover
           fsno,            &! fraction of snow cover on ground
           sigf,            &! fraction of veg cover, excluding snow-covered veg [-]
           green,           &! leaf greenness
           lai,             &! leaf area index
           tlai,            &! leaf area index
           sai,             &! stem area index
           tsai,            &! stem area index
           coszen,          &! cosine of solar zenith angle
           alb,             &! averaged albedo [-]
           ssun,            &! sunlit canopy absorption for solar radiation (0-1)
           ssha,            &! shaded canopy absorption for solar radiation (0-1)
           thermk,          &! canopy gap fraction for tir radiation
           extkb,           &! (k, g(mu)/mu) direct solar extinction coefficient
           extkd,           &! diffuse and scattered diffuse PAR extinction coefficient
           zwt,             &! the depth to water table [m]
           wa,              &! water storage in aquifer [mm]
 
           t_lake,          &! lake layer teperature [K]
           lake_icefrac,    &! lake mass fraction of lake layer that is frozen

         ! Additional variables required by reginal model (such as WRF & RSM) 
           trad,            &! radiative temperature of surface [K]
           tref,            &! 2 m height air temperature [kelvin]
           qref,            &! 2 m height air specific humidity
           rst,             &! canopy stomatal resistance (s/m)
           emis,            &! averaged bulk surface emissivity
           z0m,             &! effective roughness [m]
           zol,             &! dimensionless height (z/L) used in Monin-Obukhov theory
           rib,             &! bulk Richardson number in surface layer
           ustar,           &! u* in similarity theory [m/s]
           qstar,           &! q* in similarity theory [kg/kg]
           tstar,           &! t* in similarity theory [K]
           fm,              &! integral of profile function for momentum
           fh,              &! integral of profile function for heat
           fq                ! integral of profile function for moisture

     ! PFT/PC time variabls
#ifdef PFT_CLASSIFICATION
     read (lhistTimeVar)    &!
           tleaf_p,         &! shaded leaf temperature [K]
           ldew_p,          &! depth of water on foliage [mm]
           sigf_p,          &! fraction of veg cover, excluding snow-covered veg [-]
           lai_p,           &! leaf area index
           tlai_p,          &! true leaf area index
           sai_p,           &! stem area index
           tsai_p,          &! true stem area index
           ssun_p,          &! sunlit canopy absorption for solar radiation (0-1)      
           ssha_p,          &! shaded canopy absorption for solar radiation (0-1)
           thermk_p,        &! canopy gap fraction for tir radiation
           extkb_p,         &! (k, g(mu)/mu) direct solar extinction coefficient
           extkd_p,         &! diffuse and scattered diffuse PAR extinction coefficient
           tref_p,          &! 2 m height air temperature [kelvin]
           qref_p,          &! 2 m height air specific humidity
           rst_p,           &! canopy stomatal resistance (s/m)
           z0m_p             ! effective roughness [m]                                 
#endif

#ifdef PC_CLASSIFICATION
     read (lhistTimeVar)    &!
           tleaf_c,         &! leaf temperature [K]
           ldew_c,          &! depth of water on foliage [mm]
           sigf_c,          &! fraction of veg cover, excluding snow-covered veg [-]
           lai_c,           &! leaf area index
           tlai_c,          &! true leaf area index
           sai_c,           &! stem area index
           tsai_c,          &! true stem area index
           ssun_c,          &! sunlit canopy absorption for solar radiation (0-1)
           ssha_c,          &! shaded canopy absorption for solar radiation (0-1)
           thermk_c,        &! canopy gap fraction for tir radiation
           fshade_c,        &! canopy gap fraction for tir radiation
           extkb_c,         &! (k, g(mu)/mu) direct solar extinction coefficient
           extkd_c,         &! diffuse and scattered diffuse PAR extinction coefficient
           rst_c,           &! canopy stomatal resistance (s/m)
           z0m_c             ! effective roughness [m]                                 
#endif
      
           IF (id(1) /= idate(1) .or. id(2) /= idate(2) .or. id(3) /= idate(3)) THEN
              print*, 'id = ', id, 'idate = ', idate
              print*, 'The date of initial data is NOT IDENTICAL TO initial set-up'
              CALL abort
           ENDIF

     close(lhistTimeVar)

  END SUBROUTINE READ_TimeVariables


  SUBROUTINE WRITE_TimeVariables (idate,dir_restart_hist,casename)
! --------------------------------------------------------------------
! Write out the model variables for restart run [histTimeVar]
! --------------------------------------------------------------------

     USE precision
     USE MOD_PFTimeVars
     USE MOD_PCTimeVars
     IMPLICIT NONE
     
     INTEGER, intent(in) :: idate(3)     !calendar (year, julian day, seconds)
     CHARACTER(LEN=255), intent(in) :: dir_restart_hist
     CHARACTER(LEN=256), intent(in) :: casename

     INTEGER :: lhistTimeVar             !logical unit number of restart time-varying file
     INTEGER :: id(3)                    !calendar (year, julian day, seconds), temporal
     CHARACTER(LEN=255) :: cdate         !character for date
     CHARACTER(LEN=256) :: fhistTimeVar  !file name of time-varying file

! ...............................................................

     id(:) = idate(:)
     CALL adj2begin(id)

     ! the model variables for restart run 
     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') id(1), id(2), id(3)

     lhistTimeVar = 100
     fhistTimeVar = trim(dir_restart_hist)//trim(casename)//'-'//'rstTimeVar'//'-'//trim(cdate)
     print*,trim(fhistTimeVar)
     open(unit=lhistTimeVar,file=trim(fhistTimeVar),status='unknown',&
                            form='unformatted',action='write')

     write(lhistTimeVar) id, & !
         ! Time-varying state variables which reaquired by restart run
           z_sno,           &! node depth [m]
           dz_sno,          &! interface depth [m]
           t_soisno,        &! soil temperature [K]
           wliq_soisno,     &! liquid water in layers [kg/m2]
           wice_soisno,     &! ice lens in layers [kg/m2]
           t_grnd,          &! ground surface temperature [K]
           tleaf,           &! leaf temperature [K]
           ldew,            &! depth of water on foliage [mm]
           sag,             &! non dimensional snow age [-]
           scv,             &! snow cover, water equivalent [mm]
           snowdp,          &! snow depth [meter]
           fveg,            &! fraction of vegetation cover
           fsno,            &! fraction of snow cover on ground
           sigf,            &! fraction of veg cover, excluding snow-covered veg [-]
           green,           &! leaf greenness
           lai,             &! leaf area index
           tlai,            &! leaf area index
           sai,             &! stem area index
           tsai,            &! stem area index
           coszen,          &! cosine of solar zenith angle
           alb,             &! averaged albedo [-]
           ssun,            &! sunlit canopy absorption for solar radiation (0-1)
           ssha,            &! shaded canopy absorption for solar radiation (0-1)
           thermk,          &! canopy gap fraction for tir radiation
           extkb,           &! (k, g(mu)/mu) direct solar extinction coefficient
           extkd,           &! diffuse and scattered diffuse PAR extinction coefficient
           zwt,             &! the depth to water table [m]
           wa,              &! water storage in aquifer [mm]
 
           t_lake,          &! lake layer teperature [K]
           lake_icefrac,    &! lake mass fraction of lake layer that is frozen

         ! Additional variables required by reginal model (such as WRF & RSM) 
           trad,            &! radiative temperature of surface [K]
           tref,            &! 2 m height air temperature [kelvin]
           qref,            &! 2 m height air specific humidity
           rst,             &! canopy stomatal resistance (s/m)
           emis,            &! averaged bulk surface emissivity
           z0m,             &! effective roughness [m]
           zol,             &! dimensionless height (z/L) used in Monin-Obukhov theory
           rib,             &! bulk Richardson number in surface layer
           ustar,           &! u* in similarity theory [m/s]
           qstar,           &! q* in similarity theory [kg/kg]
           tstar,           &! t* in similarity theory [K]
           fm,              &! integral of profile function for momentum
           fh,              &! integral of profile function for heat
           fq                ! integral of profile function for moisture

      ! PFT/PC time variabls
#ifdef PFT_CLASSIFICATION
     write(lhistTimeVar)    &!
           tleaf_p,         &! shaded leaf temperature [K]
           ldew_p,          &! depth of water on foliage [mm]
           sigf_p,          &! fraction of veg cover, excluding snow-covered veg [-]
           lai_p,           &! leaf area index
           tlai_p,          &! true leaf area index
           sai_p,           &! stem area index
           tsai_p,          &! true stem area index
           ssun_p,          &! sunlit canopy absorption for solar radiation (0-1)       
           ssha_p,          &! shaded canopy absorption for solar radiation (0-1)
           thermk_p,        &! canopy gap fraction for tir radiation
           extkb_p,         &! (k, g(mu)/mu) direct solar extinction coefficient
           extkd_p,         &! diffuse and scattered diffuse PAR extinction coefficient
           tref_p,          &! 2 m height air temperature [kelvin]
           qref_p,          &! 2 m height air specific humidity
           rst_p,           &! canopy stomatal resistance (s/m)
           z0m_p             ! effective roughness [m]                                 
#endif

#ifdef PC_CLASSIFICATION
     WRITE(lhistTimeVar)    &!
           tleaf_c,         &! leaf temperature [K]
           ldew_c,          &! depth of water on foliage [mm]
           sigf_c,          &! fraction of veg cover, excluding snow-covered veg [-]
           lai_c,           &! leaf area index
           tlai_c,          &! true leaf area index
           sai_c,           &! stem area index
           tsai_c,          &! true stem area index
           ssun_c,          &! sunlit canopy absorption for solar radiation (0-1)
           ssha_c,          &! shaded canopy absorption for solar radiation (0-1)
           thermk_c,        &! canopy gap fraction for tir radiation
           fshade_c,        &! canopy gap fraction for tir radiation
           extkb_c,         &! (k, g(mu)/mu) direct solar extinction coefficient
           extkd_c,         &! diffuse and scattered diffuse PAR extinction coefficient
           rst_c,           &! canopy stomatal resistance (s/m)
           z0m_c             ! effective roughness [m]                                 
#endif
 
     close(lhistTimeVar)

  END SUBROUTINE WRITE_TimeVariables

  SUBROUTINE deallocate_TimeVariables
! --------------------------------------------------
! Deallocates memory for CLM 1d [numpatch] variables
! --------------------------------------------------
     USE MOD_PFTimeVars
     USE MOD_PCTimeVars

     deallocate (z_sno        )
     deallocate (dz_sno       )
     deallocate (t_soisno     )
     deallocate (wliq_soisno  )
     deallocate (wice_soisno  )
     deallocate (h2osoi       )
     deallocate (rstfac       )
     deallocate (t_grnd       )
     deallocate (tleaf        )
     deallocate (ldew         )
     deallocate (sag          )
     deallocate (scv          )
     deallocate (snowdp       )
     deallocate (fveg         )
     deallocate (fsno         )
     deallocate (sigf         )
     deallocate (green        )
     deallocate (tlai         )
     deallocate (lai          )
     deallocate (laisun       )
     deallocate (laisha       )
     deallocate (tsai         )
     deallocate (sai          )
     deallocate (coszen       )
     deallocate (alb          )
     deallocate (ssun         )
     deallocate (ssha         )
     deallocate (thermk       )
     deallocate (extkb        )
     deallocate (extkd        )
     deallocate (zwt          )
     deallocate (wa           )
     deallocate (wat          )

     deallocate (t_lake       ) 
     deallocate (lake_icefrac )

     deallocate (trad         )
     deallocate (tref         )
     deallocate (qref         )
     deallocate (rst          )
     deallocate (emis         )
     deallocate (z0m          )
     deallocate (displa       )
     deallocate (zol          )
     deallocate (rib          )
     deallocate (ustar        )
     deallocate (qstar        )
     deallocate (tstar        )
     deallocate (fm           )
     deallocate (fh           )
     deallocate (fq           )
     
#ifdef PFT_CLASSIFICATION
     CALL deallocate_PFTimeVars
#endif

#ifdef PC_CLASSIFICATION
     CALL deallocate_PCTimeVars
#endif
  
  END SUBROUTINE deallocate_TimeVariables

END MODULE MOD_TimeVariables
! ---------- EOP ------------
