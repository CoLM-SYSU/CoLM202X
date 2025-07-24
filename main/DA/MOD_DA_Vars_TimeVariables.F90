#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_TimeVariables
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Process time-varying state variables for data assimilation
!
! AUTHOR:
!   Lu Li, 12/2024: Initial version
!   Lu Li, 07/2025: Remove unused variables and clean codes
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_TimeManager
   USE MOD_Namelist, only: DEF_DA_ENS
   IMPLICIT NONE
   SAVE

   ! public functions
   PUBLIC :: allocate_TimeVariables_ens
   PUBLIC :: deallocate_TimeVariables_ens
   PUBLIC :: READ_TimeVariables_ens
   PUBLIC :: WRITE_TimeVariables_ens
#ifdef RangeCheck
   PUBLIC :: check_TimeVariables_ens
#endif

   ! define variables
   ! Time-varying state variables which required by restart run
   real(r8), allocatable :: z_sno_ens       (:,:,:) ! node depth [m]
   real(r8), allocatable :: dz_sno_ens      (:,:,:) ! interface depth [m]
   real(r8), allocatable :: t_soisno_ens    (:,:,:) ! soil temperature [K]
   real(r8), allocatable :: wliq_soisno_ens (:,:,:) ! liquid water in layers [kg/m2]
   real(r8), allocatable :: wice_soisno_ens (:,:,:) ! ice lens in layers [kg/m2]
   real(r8), allocatable :: smp_ens         (:,:,:) ! soil matrix potential [mm]
   real(r8), allocatable :: hk_ens          (:,:,:) ! hydraulic conductivity [mm h2o/s]
   real(r8), allocatable :: t_grnd_ens        (:,:) ! ground surface temperature [K]
   real(r8), allocatable :: tleaf_ens         (:,:) ! leaf temperature [K]
   real(r8), allocatable :: ldew_ens          (:,:) ! depth of water on foliage [mm]
   real(r8), allocatable :: ldew_rain_ens     (:,:) ! depth of rain on foliage [mm]
   real(r8), allocatable :: ldew_snow_ens     (:,:) ! depth of rain on foliage [mm]
   real(r8), allocatable :: fwet_snow_ens     (:,:) ! vegetation snow fractional cover [-]
   real(r8), allocatable :: sag_ens           (:,:) ! non dimensional snow age [-]
   real(r8), allocatable :: scv_ens           (:,:) ! snow cover, water equivalent [mm]
   real(r8), allocatable :: snowdp_ens        (:,:) ! snow depth [meter]
   real(r8), allocatable :: fveg_ens          (:,:) ! fraction of vegetation cover
   real(r8), allocatable :: fsno_ens          (:,:) ! fraction of snow cover on ground
   real(r8), allocatable :: sigf_ens          (:,:) ! fraction of veg cover, excluding snow-covered veg [-]
   real(r8), allocatable :: green_ens         (:,:) ! leaf greenness
   real(r8), allocatable :: tlai_ens          (:,:) ! leaf area index
   real(r8), allocatable :: lai_ens           (:,:) ! leaf area index
   real(r8), allocatable :: tsai_ens          (:,:) ! stem area index 
   real(r8), allocatable :: sai_ens           (:,:) ! stem area index
   real(r8), allocatable :: alb_ens       (:,:,:,:) ! averaged albedo [-]
   real(r8), allocatable :: ssun_ens      (:,:,:,:) ! sunlit canopy absorption for solar radiation (0-1)
   real(r8), allocatable :: ssha_ens      (:,:,:,:) ! shaded canopy absorption for solar radiation (0-1)
   real(r8), allocatable :: ssoi_ens      (:,:,:,:) ! soil absorption for solar radiation (0-1)
   real(r8), allocatable :: ssno_ens      (:,:,:,:) ! snow absorption for solar radiation (0-1)
   real(r8), allocatable :: thermk_ens        (:,:) ! canopy gap fraction for tir radiation
   real(r8), allocatable :: extkb_ens         (:,:) ! (k, g(mu)/mu) direct solar extinction coefficient
   real(r8), allocatable :: extkd_ens         (:,:) ! diffuse and scattered diffuse PAR extinction coefficient
   real(r8), allocatable :: zwt_ens           (:,:) ! the depth to water table [m]
   real(r8), allocatable :: wdsrf_ens         (:,:) ! depth of surface water [mm]
   real(r8), allocatable :: wa_ens            (:,:) ! water storage in aquifer [mm]
   real(r8), allocatable :: wetwat_ens        (:,:) ! water storage in wetland [mm]
   real(r8), allocatable :: t_lake_ens      (:,:,:) ! lake layer teperature [K]
   real(r8), allocatable :: lake_icefrac_ens(:,:,:) ! lake mass fraction of lake layer that is frozen
   real(r8), allocatable :: savedtke1_ens     (:,:) ! top level eddy conductivity (W/m K)

   ! diagnostic variables for DA
   real(r8), allocatable :: h2osoi_ens      (:,:,:) ! volumetric soil water in layers [m3/m3]
   real(r8), allocatable :: t_brt_ens       (:,:,:) ! brightness temperature for radiance calculation [K]
   real(r8), allocatable :: t_brt             (:,:) ! brightness temperature for radiance calculation [K]
   real(r8), allocatable :: trad_ens          (:,:) ! radiative temperature of surface [K]
   real(r8), allocatable :: tref_ens          (:,:) ! 2 m height air temperature [kelvin]
   real(r8), allocatable :: qref_ens          (:,:) ! 2 m height air specific humidity
   real(r8), allocatable :: ustar_ens         (:,:) ! u* in similarity theory [m/s]
   real(r8), allocatable :: qstar_ens         (:,:) ! q* in similarity theory [kg/kg]
   real(r8), allocatable :: tstar_ens         (:,:) ! t* in similarity theory [K]
   real(r8), allocatable :: fm_ens            (:,:) ! integral of profile FUNCTION for momentum
   real(r8), allocatable :: fh_ens            (:,:) ! integral of profile FUNCTION for heat
   real(r8), allocatable :: fq_ens            (:,:) ! integral of profile FUNCTION for moisture



!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE allocate_TimeVariables_ens()

!-----------------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Vars_Global
      USE MOD_SPMD_Task
      USE MOD_LandPatch, only: numpatch
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            ! allocate all time-varying state variables
            allocate (z_sno_ens      (maxsnl+1:0,      DEF_DA_ENS,numpatch)); z_sno_ens       (:,:,:) = spval
            allocate (dz_sno_ens     (maxsnl+1:0,      DEF_DA_ENS,numpatch)); dz_sno_ens      (:,:,:) = spval
            allocate (t_soisno_ens   (maxsnl+1:nl_soil,DEF_DA_ENS,numpatch)); t_soisno_ens    (:,:,:) = spval
            allocate (wliq_soisno_ens(maxsnl+1:nl_soil,DEF_DA_ENS,numpatch)); wliq_soisno_ens (:,:,:) = spval
            allocate (wice_soisno_ens(maxsnl+1:nl_soil,DEF_DA_ENS,numpatch)); wice_soisno_ens (:,:,:) = spval
            allocate (smp_ens               (1:nl_soil,DEF_DA_ENS,numpatch)); smp_ens         (:,:,:) = spval
            allocate (hk_ens                (1:nl_soil,DEF_DA_ENS,numpatch)); hk_ens          (:,:,:) = spval
            allocate (t_grnd_ens                      (DEF_DA_ENS,numpatch)); t_grnd_ens        (:,:) = spval
            allocate (tleaf_ens                       (DEF_DA_ENS,numpatch)); tleaf_ens         (:,:) = spval
            allocate (ldew_ens                        (DEF_DA_ENS,numpatch)); ldew_ens          (:,:) = spval
            allocate (ldew_rain_ens                   (DEF_DA_ENS,numpatch)); ldew_rain_ens     (:,:) = spval
            allocate (ldew_snow_ens                   (DEF_DA_ENS,numpatch)); ldew_snow_ens     (:,:) = spval
            allocate (fwet_snow_ens                   (DEF_DA_ENS,numpatch)); fwet_snow_ens     (:,:) = spval
            allocate (sag_ens                         (DEF_DA_ENS,numpatch)); sag_ens           (:,:) = spval
            allocate (scv_ens                         (DEF_DA_ENS,numpatch)); scv_ens           (:,:) = spval
            allocate (snowdp_ens                      (DEF_DA_ENS,numpatch)); snowdp_ens        (:,:) = spval
            allocate (fveg_ens                        (DEF_DA_ENS,numpatch)); fveg_ens          (:,:) = spval
            allocate (fsno_ens                        (DEF_DA_ENS,numpatch)); fsno_ens          (:,:) = spval
            allocate (sigf_ens                        (DEF_DA_ENS,numpatch)); sigf_ens          (:,:) = spval
            allocate (green_ens                       (DEF_DA_ENS,numpatch)); green_ens         (:,:) = spval
            allocate (tlai_ens                        (DEF_DA_ENS,numpatch)); tlai_ens          (:,:) = spval
            allocate (lai_ens                         (DEF_DA_ENS,numpatch)); lai_ens           (:,:) = spval
            allocate (tsai_ens                        (DEF_DA_ENS,numpatch)); tsai_ens          (:,:) = spval
            allocate (sai_ens                         (DEF_DA_ENS,numpatch)); sai_ens           (:,:) = spval
            allocate (alb_ens                     (2,2,DEF_DA_ENS,numpatch)); alb_ens       (:,:,:,:) = spval
            allocate (ssun_ens                    (2,2,DEF_DA_ENS,numpatch)); ssun_ens      (:,:,:,:) = spval
            allocate (ssha_ens                    (2,2,DEF_DA_ENS,numpatch)); ssha_ens      (:,:,:,:) = spval
            allocate (ssoi_ens                    (2,2,DEF_DA_ENS,numpatch)); ssoi_ens      (:,:,:,:) = spval
            allocate (ssno_ens                    (2,2,DEF_DA_ENS,numpatch)); ssno_ens      (:,:,:,:) = spval
            allocate (thermk_ens                      (DEF_DA_ENS,numpatch)); thermk_ens        (:,:) = spval
            allocate (extkb_ens                       (DEF_DA_ENS,numpatch)); extkb_ens         (:,:) = spval
            allocate (extkd_ens                       (DEF_DA_ENS,numpatch)); extkd_ens         (:,:) = spval
            allocate (zwt_ens                         (DEF_DA_ENS,numpatch)); zwt_ens           (:,:) = spval
            allocate (wdsrf_ens                       (DEF_DA_ENS,numpatch)); wdsrf_ens         (:,:) = spval
            allocate (wa_ens                          (DEF_DA_ENS,numpatch)); wa_ens            (:,:) = spval
            allocate (wetwat_ens                      (DEF_DA_ENS,numpatch)); wetwat_ens        (:,:) = spval
            allocate (t_lake_ens              (nl_lake,DEF_DA_ENS,numpatch)); t_lake_ens      (:,:,:) = spval
            allocate (lake_icefrac_ens        (nl_lake,DEF_DA_ENS,numpatch)); lake_icefrac_ens(:,:,:) = spval
            allocate (savedtke1_ens                   (DEF_DA_ENS,numpatch)); savedtke1_ens     (:,:) = spval
            
            ! diagnostic variables for DA
            allocate (h2osoi_ens            (1:nl_soil,DEF_DA_ENS,numpatch)); h2osoi_ens      (:,:,:) = spval
            allocate (t_brt_ens                     (2,DEF_DA_ENS,numpatch)); t_brt_ens       (:,:,:) = spval
            allocate (t_brt                                    (2,numpatch)); t_brt             (:,:) = spval
            allocate (trad_ens                        (DEF_DA_ENS,numpatch)); trad_ens          (:,:) = spval
            allocate (tref_ens                        (DEF_DA_ENS,numpatch)); tref_ens          (:,:) = spval
            allocate (qref_ens                        (DEF_DA_ENS,numpatch)); qref_ens          (:,:) = spval
            allocate (ustar_ens                       (DEF_DA_ENS,numpatch)); ustar_ens         (:,:) = spval
            allocate (qstar_ens                       (DEF_DA_ENS,numpatch)); qstar_ens         (:,:) = spval
            allocate (tstar_ens                       (DEF_DA_ENS,numpatch)); tstar_ens         (:,:) = spval
            allocate (fm_ens                          (DEF_DA_ENS,numpatch)); fm_ens            (:,:) = spval
            allocate (fh_ens                          (DEF_DA_ENS,numpatch)); fh_ens            (:,:) = spval
            allocate (fq_ens                          (DEF_DA_ENS,numpatch)); fq_ens            (:,:) = spval
         ENDIF
      ENDIF

   END SUBROUTINE allocate_TimeVariables_ens

!-----------------------------------------------------------------------------

   SUBROUTINE deallocate_TimeVariables_ens()

!-----------------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_LandPatch, only: numpatch
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            deallocate (z_sno_ens                  )
            deallocate (dz_sno_ens                 )
            deallocate (t_soisno_ens               )
            deallocate (wliq_soisno_ens            )
            deallocate (wice_soisno_ens            )
            deallocate (smp_ens                    )
            deallocate (hk_ens                     )
            deallocate (t_grnd_ens                 )
            deallocate (tleaf_ens                  )
            deallocate (ldew_ens                   )
            deallocate (ldew_rain_ens              )
            deallocate (ldew_snow_ens              )
            deallocate (fwet_snow_ens              )
            deallocate (sag_ens                    )
            deallocate (scv_ens                    )
            deallocate (snowdp_ens                 )
            deallocate (fveg_ens                   )
            deallocate (fsno_ens                   )
            deallocate (sigf_ens                   )
            deallocate (green_ens                  )
            deallocate (tlai_ens                   )
            deallocate (lai_ens                    )
            deallocate (tsai_ens                   )
            deallocate (sai_ens                    )
            deallocate (alb_ens                    )
            deallocate (ssun_ens                   )
            deallocate (ssha_ens                   )
            deallocate (ssoi_ens                   )
            deallocate (ssno_ens                   )
            deallocate (thermk_ens                 )
            deallocate (extkb_ens                  )
            deallocate (extkd_ens                  )
            deallocate (zwt_ens                    )
            deallocate (wdsrf_ens                  )
            deallocate (wa_ens                     )
            deallocate (wetwat_ens                 )
            deallocate (t_lake_ens                 )
            deallocate (lake_icefrac_ens           )
            deallocate (savedtke1_ens              )

            deallocate (h2osoi_ens                 )
            deallocate (t_brt_ens                  )
            deallocate (t_brt                      )
            deallocate (trad_ens                   )
            deallocate (tref_ens                   )
            deallocate (qref_ens                   )
            deallocate (ustar_ens                  )
            deallocate (qstar_ens                  )
            deallocate (tstar_ens                  )
            deallocate (fm_ens                     )
            deallocate (fh_ens                     )
            deallocate (fq_ens                     )
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_TimeVariables_ens

!-----------------------------------------------------------------------------

   SUBROUTINE WRITE_TimeVariables_ens(idate, lc_year, site, dir_restart)

!-----------------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_Namelist, only : DEF_REST_CompressLevel, DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS, &
         DEF_USE_IRRIGATION, DEF_USE_Dynamic_Lake
      USE MOD_LandPatch
      USE MOD_NetCDFVector
      USE MOD_Vars_Global
      USE MOD_Vars_TimeInvariants, only : dz_lake
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer, intent(in) :: idate(3)
      integer, intent(in) :: lc_year      !year of land cover type data
      character(len=*), intent(in) :: site
      character(len=*), intent(in) :: dir_restart

!------------------------ Local Variables ------------------------------------
      character(len=256) :: file_restart
      character(len=14)  :: cdate
      character(len=256) :: cyear         !character for lc_year
      integer :: compress
      integer :: i

!-----------------------------------------------------------------------------
      compress = DEF_REST_CompressLevel

      ! land cover type year
      write(cyear,'(i4.4)') lc_year
      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)

      IF (p_is_master) THEN
         CALL system('mkdir -p ' // trim(dir_restart)//'/'//trim(cdate))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      file_restart = trim(dir_restart)// '/'//trim(cdate)//'/' // trim(site) //'_restart_'//trim(cdate)//'_lc'//trim(cyear)//'_ens'//'.nc'

      CALL ncio_create_file_vector      (file_restart, landpatch)

      CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'snow',     -maxsnl       )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'snowp1',   -maxsnl+1     )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'soilsnow', nl_soil-maxsnl)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'soil',     nl_soil)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'lake',     nl_lake)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'band', 2)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'rtyp', 2)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'ens', DEF_DA_ENS)

      ! Time-varying state variables which reaquired by restart run
      CALL ncio_write_vector (file_restart, 'z_sno   '   , 'snow',     -maxsnl,        'ens', DEF_DA_ENS, 'patch', landpatch, z_sno_ens,       compress) ! node depth [m]
      CALL ncio_write_vector (file_restart, 'dz_sno  '   , 'snow',     -maxsnl,        'ens', DEF_DA_ENS, 'patch', landpatch, dz_sno_ens,      compress) ! interface depth [m]
      CALL ncio_write_vector (file_restart, 't_soisno'   , 'soilsnow', nl_soil-maxsnl, 'ens', DEF_DA_ENS, 'patch', landpatch, t_soisno_ens,    compress) ! soil temperature [K]
      CALL ncio_write_vector (file_restart, 'wliq_soisno', 'soilsnow', nl_soil-maxsnl, 'ens', DEF_DA_ENS, 'patch', landpatch, wliq_soisno_ens, compress) ! liquid water in layers [kg/m2]
      CALL ncio_write_vector (file_restart, 'wice_soisno', 'soilsnow', nl_soil-maxsnl, 'ens', DEF_DA_ENS, 'patch', landpatch, wice_soisno_ens, compress) ! ice lens in layers [kg/m2]
      CALL ncio_write_vector (file_restart, 'smp',         'soil',     nl_soil,        'ens', DEF_DA_ENS, 'patch', landpatch, smp_ens,         compress) ! soil matrix potential [mm]
      CALL ncio_write_vector (file_restart, 'hk',          'soil',     nl_soil,        'ens', DEF_DA_ENS, 'patch', landpatch, hk_ens,          compress) ! hydraulic conductivity [mm h2o/s]
      CALL ncio_write_vector (file_restart, 't_grnd',    'ens', DEF_DA_ENS, 'patch', landpatch, t_grnd_ens,    compress) ! ground surface temperature [K]
      CALL ncio_write_vector (file_restart, 'tleaf',     'ens', DEF_DA_ENS, 'patch', landpatch, tleaf_ens,     compress) ! leaf temperature [K]
      CALL ncio_write_vector (file_restart, 'ldew',      'ens', DEF_DA_ENS, 'patch', landpatch, ldew_ens,      compress) ! depth of water on foliage [mm]
      CALL ncio_write_vector (file_restart, 'ldew_rain', 'ens', DEF_DA_ENS, 'patch', landpatch, ldew_rain_ens, compress) ! depth of water on foliage [mm]
      CALL ncio_write_vector (file_restart, 'ldew_snow', 'ens', DEF_DA_ENS, 'patch', landpatch, ldew_snow_ens, compress) ! depth of water on foliage [mm]
      CALL ncio_write_vector (file_restart, 'fwet_snow', 'ens', DEF_DA_ENS, 'patch', landpatch, fwet_snow_ens, compress) ! vegetation snow fractional cover [-]
      CALL ncio_write_vector (file_restart, 'sag',       'ens', DEF_DA_ENS, 'patch', landpatch, sag_ens,       compress) ! non dimensional snow age [-]
      CALL ncio_write_vector (file_restart, 'scv',       'ens', DEF_DA_ENS, 'patch', landpatch, scv_ens,       compress) ! snow cover, water equivalent [mm]
      CALL ncio_write_vector (file_restart, 'snowdp',    'ens', DEF_DA_ENS, 'patch', landpatch, snowdp_ens,    compress) ! snow depth [meter]
      CALL ncio_write_vector (file_restart, 'fveg',      'ens', DEF_DA_ENS, 'patch', landpatch, fveg_ens,      compress) ! fraction of vegetation cover
      CALL ncio_write_vector (file_restart, 'fsno',      'ens', DEF_DA_ENS, 'patch', landpatch, fsno_ens,      compress) ! fraction of snow cover on ground
      CALL ncio_write_vector (file_restart, 'sigf',      'ens', DEF_DA_ENS, 'patch', landpatch, sigf_ens,      compress) ! fraction of veg cover, excluding snow-covered veg [-]
      CALL ncio_write_vector (file_restart, 'green',     'ens', DEF_DA_ENS, 'patch', landpatch, green_ens,     compress) ! leaf greenness
      CALL ncio_write_vector (file_restart, 'tlai',      'ens', DEF_DA_ENS, 'patch', landpatch, tlai_ens,      compress) ! leaf area index
      CALL ncio_write_vector (file_restart, 'lai',       'ens', DEF_DA_ENS, 'patch', landpatch, lai_ens,       compress) ! leaf area index
      CALL ncio_write_vector (file_restart, 'tsai',      'ens', DEF_DA_ENS, 'patch', landpatch, tsai_ens,      compress) ! stem area index
      CALL ncio_write_vector (file_restart, 'sai',       'ens', DEF_DA_ENS, 'patch', landpatch, sai_ens,       compress) ! stem area index
      CALL ncio_write_vector (file_restart, 'alb     '   , 'band', 2, 'rtyp', 2, 'ens', DEF_DA_ENS, 'patch', landpatch, alb_ens,  compress) ! averaged albedo [-]
      CALL ncio_write_vector (file_restart, 'ssun    '   , 'band', 2, 'rtyp', 2, 'ens', DEF_DA_ENS, 'patch', landpatch, ssun_ens, compress) ! sunlit canopy absorption for solar radiation (0-1)
      CALL ncio_write_vector (file_restart, 'ssha    '   , 'band', 2, 'rtyp', 2, 'ens', DEF_DA_ENS, 'patch', landpatch, ssha_ens, compress) ! shaded canopy absorption for solar radiation (0-1)
      CALL ncio_write_vector (file_restart, 'ssoi    '   , 'band', 2, 'rtyp', 2, 'ens', DEF_DA_ENS, 'patch', landpatch, ssoi_ens, compress) ! shaded canopy absorption for solar radiation (0-1)
      CALL ncio_write_vector (file_restart, 'ssno    '   , 'band', 2, 'rtyp', 2, 'ens', DEF_DA_ENS, 'patch', landpatch, ssno_ens, compress) ! shaded canopy absorption for solar radiation (0-1)
      CALL ncio_write_vector (file_restart, 'thermk  ',  'ens', DEF_DA_ENS, 'patch', landpatch, thermk_ens,    compress) ! canopy gap fraction for tir radiation
      CALL ncio_write_vector (file_restart, 'extkb',     'ens', DEF_DA_ENS, 'patch', landpatch, extkb_ens,     compress) ! (k, g(mu)/mu) direct solar extinction coefficient
      CALL ncio_write_vector (file_restart, 'extkd',     'ens', DEF_DA_ENS, 'patch', landpatch, extkd_ens,     compress) ! diffuse and scattered diffuse PAR extinction coefficient
      CALL ncio_write_vector (file_restart, 'zwt',       'ens', DEF_DA_ENS, 'patch', landpatch, zwt_ens,       compress) ! the depth to water table [m]
      CALL ncio_write_vector (file_restart, 'wdsrf',     'ens', DEF_DA_ENS, 'patch', landpatch, wdsrf_ens,     compress) ! depth of surface water [mm]
      CALL ncio_write_vector (file_restart, 'wa',        'ens', DEF_DA_ENS, 'patch', landpatch, wa_ens,        compress) ! water storage in aquifer [mm]
      CALL ncio_write_vector (file_restart, 'wetwat',    'ens', DEF_DA_ENS, 'patch', landpatch, wetwat_ens,    compress) ! water storage in wetland [mm]
      CALL ncio_write_vector (file_restart, 't_lake  '   , 'lake', nl_lake, 'ens', DEF_DA_ENS, 'patch', landpatch, t_lake_ens,       compress)
      CALL ncio_write_vector (file_restart, 'lake_icefrc', 'lake', nl_lake, 'ens', DEF_DA_ENS, 'patch', landpatch, lake_icefrac_ens, compress)
      CALL ncio_write_vector (file_restart, 'savedtke1  ', 'ens', DEF_DA_ENS, 'patch', landpatch, savedtke1_ens, compress)

   END SUBROUTINE WRITE_TimeVariables_ens

!-----------------------------------------------------------------------------

   SUBROUTINE READ_TimeVariables_ens(idate, lc_year, site, dir_restart)

!-----------------------------------------------------------------------------
      USE MOD_Namelist
      USE MOD_SPMD_Task
      USE MOD_NetCDFVector
#ifdef RangeCheck
      USE MOD_RangeCheck
#endif
      USE MOD_LandPatch
      USE MOD_Vars_Global
      USE MOD_Vars_TimeInvariants, only: dz_lake
      IMPLICIT NONE

!------------------------ Dummy Arguments ------------------------------------
      integer, intent(in) :: idate(3)
      integer, intent(in) :: lc_year      !year of land cover type data
      character(len=*), intent(in) :: site
      character(len=*), intent(in) :: dir_restart

!------------------------ Local Variables ------------------------------------
      character(len=256) :: file_restart
      character(len=14)  :: cdate, cyear

!-----------------------------------------------------------------------------
#ifdef USEMPI
      CALL mpi_barrier(p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         write (*, *) 'Loading Ensemble Time Variables ...'
      END IF

      ! land cover type year
      write (cyear, '(i4.4)') lc_year

      write (cdate, '(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
      file_restart = trim(dir_restart)//'/'//trim(cdate)//'/'//trim(site)//'_restart_'//trim(cdate)//'_lc'//trim(cyear)//'_ens'//'.nc'

      ! Time-varying state variables which reaquired by restart run
      CALL ncio_read_vector(file_restart, 'z_sno   ', -maxsnl, DEF_DA_ENS, landpatch, z_sno_ens)             ! node depth [m]
      CALL ncio_read_vector(file_restart, 'dz_sno  ', -maxsnl, DEF_DA_ENS, landpatch, dz_sno_ens)            ! interface depth [m]
      CALL ncio_read_vector(file_restart, 't_soisno', nl_soil - maxsnl, DEF_DA_ENS, landpatch, t_soisno_ens)   ! soil temperature [K]
      CALL ncio_read_vector(file_restart, 'wliq_soisno', nl_soil - maxsnl, DEF_DA_ENS, landpatch, wliq_soisno_ens)! liquid water in layers [kg/m2]
      CALL ncio_read_vector(file_restart, 'wice_soisno', nl_soil - maxsnl, DEF_DA_ENS, landpatch, wice_soisno_ens)! ice lens in layers [kg/m2]
      CALL ncio_read_vector(file_restart, 'smp', nl_soil, DEF_DA_ENS, landpatch, smp_ens)        ! soil matrix potential [mm]
      CALL ncio_read_vector(file_restart, 'hk', nl_soil, DEF_DA_ENS, landpatch, hk_ens)         ! hydraulic conductivity [mm h2o/s]
      CALL ncio_read_vector(file_restart, 't_grnd  ', DEF_DA_ENS, landpatch, t_grnd_ens) ! ground surface temperature [K]
      CALL ncio_read_vector(file_restart, 'tleaf   ', DEF_DA_ENS, landpatch, tleaf_ens) ! leaf temperature [K]
      CALL ncio_read_vector(file_restart, 'ldew    ', DEF_DA_ENS, landpatch, ldew_ens) ! depth of water on foliage [mm]
      CALL ncio_read_vector(file_restart, 'ldew_rain', DEF_DA_ENS, landpatch, ldew_rain_ens) ! depth of water on foliage [mm]
      CALL ncio_read_vector(file_restart, 'ldew_snow', DEF_DA_ENS, landpatch, ldew_snow_ens) ! depth of water on foliage [mm]
      CALL ncio_read_vector(file_restart, 'fwet_snow', DEF_DA_ENS, landpatch, fwet_snow_ens) ! vegetation snow fractional cover [-]
      CALL ncio_read_vector(file_restart, 'sag     ', DEF_DA_ENS, landpatch, sag_ens) ! non dimensional snow age [-]
      CALL ncio_read_vector(file_restart, 'scv     ', DEF_DA_ENS, landpatch, scv_ens) ! snow cover, water equivalent [mm]
      CALL ncio_read_vector(file_restart, 'snowdp  ', DEF_DA_ENS, landpatch, snowdp_ens) ! snow depth [meter]
      CALL ncio_read_vector(file_restart, 'fveg    ', DEF_DA_ENS, landpatch, fveg_ens) ! fraction of vegetation cover
      CALL ncio_read_vector(file_restart, 'fsno    ', DEF_DA_ENS, landpatch, fsno_ens) ! fraction of snow cover on ground
      CALL ncio_read_vector(file_restart, 'sigf    ', DEF_DA_ENS, landpatch, sigf_ens) ! fraction of veg cover, excluding snow-covered veg [-]
      CALL ncio_read_vector(file_restart, 'green   ', DEF_DA_ENS, landpatch, green_ens) ! leaf greenness
      CALL ncio_read_vector(file_restart, 'lai     ', DEF_DA_ENS, landpatch, lai_ens) ! leaf area index
      CALL ncio_read_vector(file_restart, 'tlai    ', DEF_DA_ENS, landpatch, tlai_ens) ! leaf area index
      CALL ncio_read_vector(file_restart, 'sai     ', DEF_DA_ENS, landpatch, sai_ens) ! stem area index
      CALL ncio_read_vector(file_restart, 'tsai    ', DEF_DA_ENS, landpatch, tsai_ens) ! stem area index
      CALL ncio_read_vector(file_restart, 'alb     ', 2, 2, DEF_DA_ENS, landpatch, alb_ens) ! averaged albedo [-]
      CALL ncio_read_vector(file_restart, 'ssun    ', 2, 2, DEF_DA_ENS, landpatch, ssun_ens) ! sunlit canopy absorption for solar radiation (0-1)
      CALL ncio_read_vector(file_restart, 'ssha    ', 2, 2, DEF_DA_ENS, landpatch, ssha_ens) ! shaded canopy absorption for solar radiation (0-1)
      CALL ncio_read_vector(file_restart, 'ssoi    ', 2, 2, DEF_DA_ENS, landpatch, ssoi_ens) ! soil absorption for solar radiation (0-1)
      CALL ncio_read_vector(file_restart, 'ssno    ', 2, 2, DEF_DA_ENS, landpatch, ssno_ens) ! snow absorption for solar radiation (0-1)
      CALL ncio_read_vector(file_restart, 'thermk  ', DEF_DA_ENS, landpatch, thermk_ens) ! canopy gap fraction for tir radiation
      CALL ncio_read_vector(file_restart, 'extkb   ', DEF_DA_ENS, landpatch, extkb_ens) ! (k, g(mu)/mu) direct solar extinction coefficient
      CALL ncio_read_vector(file_restart, 'extkd   ', DEF_DA_ENS, landpatch, extkd_ens) ! diffuse and scattered diffuse PAR extinction coefficient
      CALL ncio_read_vector(file_restart, 'zwt     ', DEF_DA_ENS, landpatch, zwt_ens) ! the depth to water table [m]
      CALL ncio_read_vector(file_restart, 'wdsrf   ', DEF_DA_ENS, landpatch, wdsrf_ens) ! depth of surface water [mm]
      CALL ncio_read_vector(file_restart, 'wa      ', DEF_DA_ENS, landpatch, wa_ens) ! water storage in aquifer [mm]
      CALL ncio_read_vector(file_restart, 'wetwat  ', DEF_DA_ENS, landpatch, wetwat_ens) ! water storage in wetland [mm]
      CALL ncio_read_vector(file_restart, 't_lake  ', nl_lake, DEF_DA_ENS, landpatch, t_lake_ens) ! lake temperature [K]
      CALL ncio_read_vector(file_restart, 'lake_icefrc', nl_lake, DEF_DA_ENS, landpatch, lake_icefrac_ens) ! lake ice fraction [-]
      CALL ncio_read_vector(file_restart, 'savedtke1  ', DEF_DA_ENS, landpatch, savedtke1_ens) ! saved tke1 [m2/s2]

#ifdef RangeCheck
      CALL check_TimeVariables_ens
#endif

      IF (p_is_master) THEN
         write (*, *) 'Loading Ensemble Time Variables done.'
      END IF

   END SUBROUTINE READ_TimeVariables_ens

#ifdef RangeCheck
!-----------------------------------------------------------------------

   SUBROUTINE check_TimeVariables_ens()

!-----------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_RangeCheck
      USE MOD_Namelist, only: DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS, DEF_USE_IRRIGATION, &
         DEF_USE_SNICAR, DEF_USE_Dynamic_Lake
      USE MOD_Vars_TimeInvariants, only: dz_lake
      IMPLICIT NONE

#ifdef USEMPI
      CALL mpi_barrier(p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write (*, *) 'Checking Ensemble Time Variables ...'
      END IF

      CALL check_vector_data ('z_sno       [m]    ', z_sno_ens       ) ! node depth [m]
      CALL check_vector_data ('dz_sno      [m]    ', dz_sno_ens      ) ! interface depth [m]
      CALL check_vector_data ('t_soisno    [K]    ', t_soisno_ens    ) ! soil temperature [K]
      CALL check_vector_data ('wliq_soisno [kg/m2]', wliq_soisno_ens ) ! liquid water in layers [kg/m2]
      CALL check_vector_data ('wice_soisno [kg/m2]', wice_soisno_ens ) ! ice lens in layers [kg/m2]
      CALL check_vector_data ('smp         [mm]   ', smp_ens         ) ! soil matrix potential [mm]
      CALL check_vector_data ('hk          [mm/s] ', hk_ens          ) ! hydraulic conductivity [mm h2o/s]
      CALL check_vector_data ('t_grnd      [K]    ', t_grnd_ens      ) ! ground surface temperature [K]
      CALL check_vector_data ('tleaf       [K]    ', tleaf_ens       ) ! leaf temperature [K]
      CALL check_vector_data ('ldew        [mm]   ', ldew_ens        ) ! depth of water on foliage [mm]
      CALL check_vector_data ('ldew_rain   [mm]   ', ldew_rain_ens   ) ! depth of rain on foliage [mm]
      CALL check_vector_data ('ldew_snow   [mm]   ', ldew_snow_ens   ) ! depth of snow on foliage [mm]
      CALL check_vector_data ('fwet_snow   [mm]   ', fwet_snow_ens   ) ! vegetation snow fractional cover [-]
      CALL check_vector_data ('sag         [-]    ', sag_ens         ) ! non dimensional snow age [-]
      CALL check_vector_data ('scv         [mm]   ', scv_ens         ) ! snow cover, water equivalent [mm]
      CALL check_vector_data ('snowdp      [m]    ', snowdp_ens      ) ! snow depth [meter]
      CALL check_vector_data ('fveg        [-]    ', fveg_ens        ) ! fraction of vegetation cover
      CALL check_vector_data ('fsno        [-]    ', fsno_ens        ) ! fraction of snow cover on ground
      CALL check_vector_data ('sigf        [-]    ', sigf_ens        ) ! fraction of veg cover, excluding snow-covered veg [-]
      CALL check_vector_data ('green       [-]    ', green_ens       ) ! leaf greenness
      CALL check_vector_data ('lai         [-]    ', lai_ens         ) ! leaf area index
      CALL check_vector_data ('tlai        [-]    ', tlai_ens        ) ! leaf area index
      CALL check_vector_data ('sai         [-]    ', sai_ens         ) ! stem area index
      CALL check_vector_data ('tsai        [-]    ', tsai_ens        ) ! stem area index
      CALL check_vector_data ('alb         [-]    ', alb_ens         ) ! averaged albedo [-]
      CALL check_vector_data ('ssun        [-]    ', ssun_ens        ) ! sunlit canopy absorption for solar radiation (0-1)
      CALL check_vector_data ('ssha        [-]    ', ssha_ens        ) ! shaded canopy absorption for solar radiation (0-1)
      CALL check_vector_data ('ssoi        [-]    ', ssoi_ens        ) ! soil absorption for solar radiation (0-1)
      CALL check_vector_data ('ssno        [-]    ', ssno_ens        ) ! snow absorption for solar radiation (0-1)
      CALL check_vector_data ('thermk      [-]    ', thermk_ens      ) ! canopy gap fraction for tir radiation
      CALL check_vector_data ('extkb       [-]    ', extkb_ens       ) ! (k, g(mu)/mu) direct solar extinction coefficient
      CALL check_vector_data ('extkd       [-]    ', extkd_ens       ) ! diffuse and scattered diffuse PAR extinction coefficient
      CALL check_vector_data ('zwt         [m]    ', zwt_ens         ) ! the depth to water table [m]
      CALL check_vector_data ('wdsrf       [mm]   ', wdsrf_ens       ) ! depth of surface water [mm]
      CALL check_vector_data ('wa          [mm]   ', wa_ens          ) ! water storage in aquifer [mm]
      CALL check_vector_data ('wetwat      [mm]   ', wetwat_ens      ) ! water storage in wetland [mm]
      CALL check_vector_data ('t_lake      [K]    ', t_lake_ens      ) ! lake temperature [K]
      CALL check_vector_data ('lake_icefrc [-]    ', lake_icefrac_ens) ! lake ice fraction [-]
      CALL check_vector_data ('savedtke1   [W/m K]', savedtke1_ens   ) ! saved tke1 [m2/s2]

#ifdef USEMPI
      CALL mpi_barrier(p_comm_glb, p_err)
#endif

   END SUBROUTINE check_TimeVariables_ens
!-----------------------------------------------------------------------------
#endif

END MODULE MOD_DA_Vars_TimeVariables
#endif