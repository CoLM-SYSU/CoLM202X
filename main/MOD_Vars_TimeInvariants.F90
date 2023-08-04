#include <define.h>

! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
MODULE MOD_Vars_PFTimeInvariants
! -----------------------------------------------------------------
! !DESCRIPTION:
! Define PFT time invariables
!
! Added by Hua Yuan, 08/2019
! -----------------------------------------------------------------

  USE MOD_Precision
  USE MOD_Vars_Global
  IMPLICIT NONE
  SAVE

  ! for LULC_IGBP_PFT and LULC_IGBP_PC
  INTEGER , allocatable :: pftclass    (:)    !PFT type
  REAL(r8), allocatable :: pftfrac     (:)    !PFT fractional cover
  REAL(r8), allocatable :: htop_p      (:)    !canopy top height [m]
  REAL(r8), allocatable :: hbot_p      (:)    !canopy bottom height [m]

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_PFTimeInvariants
  PUBLIC :: READ_PFTimeInvariants
  PUBLIC :: WRITE_PFTimeInvariants
  PUBLIC :: deallocate_PFTimeInvariants
#ifdef RangeCheck
  PUBLIC :: check_PFTimeInvariants
#endif

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_PFTimeInvariants
  ! --------------------------------------------------------------------
  ! Allocates memory for CoLM PFT 1d [numpft] variables
  ! --------------------------------------------------------------------

     USE MOD_SPMD_Task
     USE MOD_LandPFT,   only : numpft
     USE MOD_Precision
     IMPLICIT NONE

     IF (p_is_worker) THEN
        IF (numpft > 0) THEN
           allocate (pftclass      (numpft))
           allocate (pftfrac       (numpft))
           allocate (htop_p        (numpft))
           allocate (hbot_p        (numpft))
        ENDIF
     ENDIF

  END SUBROUTINE allocate_PFTimeInvariants

  SUBROUTINE READ_PFTimeInvariants (file_restart)

     use MOD_NetCDFVector
     USE MOD_LandPFT
     IMPLICIT NONE

     character(LEN=*), intent(in) :: file_restart

     call ncio_read_vector (file_restart, 'pftclass', landpft, pftclass) !
     call ncio_read_vector (file_restart, 'pftfrac ', landpft, pftfrac ) !
     call ncio_read_vector (file_restart, 'htop_p  ', landpft, htop_p  ) !
     call ncio_read_vector (file_restart, 'hbot_p  ', landpft, hbot_p  ) !

  end subroutine READ_PFTimeInvariants

  SUBROUTINE WRITE_PFTimeInvariants (file_restart)

     use MOD_NetCDFVector
     use MOD_LandPFT
     USE MOD_Namelist
     USE MOD_Vars_Global
     IMPLICIT NONE

     ! Local variables
     character(len=*), intent(in) :: file_restart
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL

     call ncio_create_file_vector (file_restart, landpft)
     CALL ncio_define_dimension_vector (file_restart, landpft, 'pft')

     call ncio_write_vector (file_restart, 'pftclass', 'pft', landpft, pftclass, compress) !
     call ncio_write_vector (file_restart, 'pftfrac ', 'pft', landpft, pftfrac , compress) !
     call ncio_write_vector (file_restart, 'htop_p  ', 'pft', landpft, htop_p  , compress) !
     call ncio_write_vector (file_restart, 'hbot_p  ', 'pft', landpft, hbot_p  , compress) !

  end subroutine WRITE_PFTimeInvariants

  SUBROUTINE deallocate_PFTimeInvariants
! --------------------------------------------------
! Deallocates memory for CoLM PFT 1d [numpft] variables
! --------------------------------------------------
     USE MOD_SPMD_Task
     USE MOD_LandPFT

     IF (p_is_worker) THEN
        IF (numpft > 0) THEN
           deallocate (pftclass)
           deallocate (pftfrac )
           deallocate (htop_p  )
           deallocate (hbot_p  )
        ENDIF
     ENDIF

  END SUBROUTINE deallocate_PFTimeInvariants

#ifdef RangeCheck
  SUBROUTINE check_PFTimeInvariants ()

     use MOD_RangeCheck
     IMPLICIT NONE

     call check_vector_data ('pftfrac', pftfrac) !
     call check_vector_data ('htop_p ', htop_p ) !
     call check_vector_data ('hbot_p ', hbot_p ) !

  end subroutine check_PFTimeInvariants
#endif

END MODULE MOD_Vars_PFTimeInvariants
#endif

MODULE MOD_Vars_TimeInvariants
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

  USE MOD_Precision
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
  USE MOD_Vars_PFTimeInvariants
#endif
#ifdef BGC
  USE MOD_BGC_Vars_TimeInvariants
#endif
#ifdef URBAN_MODEL
  USE MOD_Urban_Vars_TimeInvariants
#endif
  IMPLICIT NONE
  SAVE

! -----------------------------------------------------------------
! surface classification and soil information
  INTEGER,  allocatable :: patchclass     (:)  !index of land cover type of the patches at the fraction > 0
  INTEGER,  allocatable :: patchtype      (:)  !land water type
  LOGICAL,  allocatable :: patchmask      (:)  !patch mask

  REAL(r8), allocatable :: patchlatr      (:)  !latitude in radians
  REAL(r8), allocatable :: patchlonr      (:)  !longitude in radians

  REAL(r8), allocatable :: lakedepth      (:)  !lake depth
  REAL(r8), allocatable :: dz_lake      (:,:)  !new lake scheme

  REAL(r8), allocatable :: soil_s_v_alb   (:)  !albedo of visible of the saturated soil
  REAL(r8), allocatable :: soil_d_v_alb   (:)  !albedo of visible of the dry soil
  REAL(r8), allocatable :: soil_s_n_alb   (:)  !albedo of near infrared of the saturated soil
  REAL(r8), allocatable :: soil_d_n_alb   (:)  !albedo of near infrared of the dry soil

  REAL(r8), allocatable :: vf_quartz    (:,:)  !volumetric fraction of quartz within mineral soil
  REAL(r8), allocatable :: vf_gravels   (:,:)  !volumetric fraction of gravels
  REAL(r8), allocatable :: vf_om        (:,:)  !volumetric fraction of organic matter
  REAL(r8), allocatable :: vf_sand      (:,:)  !volumetric fraction of sand
  REAL(r8), allocatable :: wf_gravels   (:,:)  !gravimetric fraction of gravels
  REAL(r8), allocatable :: wf_sand      (:,:)  !gravimetric fraction of sand
  REAL(r8), allocatable :: OM_density   (:,:)  !OM density (kg/m3)
  REAL(r8), allocatable :: BD_all       (:,:)  !bulk density of soil (GRAVELS + ORGANIC MATTER + Mineral Soils,kg/m3)

  REAL(r8), allocatable :: wfc          (:,:)  !field capacity
  REAL(r8), allocatable :: porsl        (:,:)  !fraction of soil that is voids [-]
  REAL(r8), allocatable :: psi0         (:,:)  !minimum soil suction [mm] (NOTE: "-" valued)
  REAL(r8), allocatable :: bsw          (:,:)  !clapp and hornbereger "b" parameter [-]
#ifdef vanGenuchten_Mualem_SOIL_MODEL
  REAL(r8), allocatable :: theta_r      (:,:)  ! residual moisture content [-]
  REAL(r8), allocatable :: alpha_vgm    (:,:)  ! a parameter corresponding approximately to the inverse of the air-entry value
  REAL(r8), allocatable :: L_vgm        (:,:)  ! pore-connectivity parameter [dimensionless]
  REAL(r8), allocatable :: n_vgm        (:,:)  ! a shape parameter [dimensionless]
  REAL(r8), allocatable :: sc_vgm       (:,:)  ! saturation at the air entry value in the classical vanGenuchten model [-]
  REAL(r8), allocatable :: fc_vgm       (:,:)  ! a scaling factor by using air entry value in the Mualem model [-]
#endif
  REAL(r8), allocatable :: hksati       (:,:)  !hydraulic conductivity at saturation [mm h2o/s]
  REAL(r8), allocatable :: csol         (:,:)  !heat capacity of soil solids [J/(m3 K)]
  REAL(r8), allocatable :: k_solids     (:,:)  !thermal conductivity of soil solids [W/m-K]
  REAL(r8), allocatable :: dksatu       (:,:)  !thermal conductivity of saturated soil [W/m-K]
  real(r8), allocatable :: dksatf       (:,:)  !thermal conductivity of saturated frozen soil [W/m-K]
  REAL(r8), allocatable :: dkdry        (:,:)  !thermal conductivity for dry soil  [W/(m-K)]
  REAL(r8), allocatable :: BA_alpha     (:,:)  !alpha in Balland and Arp(2005) thermal conductivity scheme
  REAL(r8), allocatable :: BA_beta      (:,:)  !beta in Balland and Arp(2005) thermal conductivity scheme
  REAL(r8), allocatable :: htop           (:)  !canopy top height [m]
  REAL(r8), allocatable :: hbot           (:)  !canopy bottom height [m]

  real(r8), allocatable :: dbedrock       (:)  !depth to bedrock
  integer , allocatable :: ibedrock       (:)  !bedrock level


  REAL(r8) :: zlnd                             !roughness length for soil [m]
  REAL(r8) :: zsno                             !roughness length for snow [m]
  REAL(r8) :: csoilc                           !drag coefficient for soil under canopy [-]
  REAL(r8) :: dewmx                            !maximum dew
  REAL(r8) :: wtfact                           !fraction of model area with high water table
  REAL(r8) :: capr                             !tuning factor to turn first layer T into surface T
  REAL(r8) :: cnfac                            !Crank Nicholson factor between 0 and 1
  REAL(r8) :: ssi                              !irreducible water saturation of snow
  REAL(r8) :: wimp                             !water impremeable if porosity less than wimp
  REAL(r8) :: pondmx                           !ponding depth (mm)
  REAL(r8) :: smpmax                           !wilting point potential in mm
  REAL(r8) :: smpmin                           !restriction for min of soil poten. (mm)
  REAL(r8) :: trsmx0                           !max transpiration for moist soil+100% veg.  [mm/s]
  REAL(r8) :: tcrit                            !critical temp. to determine rain or snow

! PUBLIC MEMBER FUNCTIONS:
  public :: allocate_TimeInvariants
  public :: deallocate_TimeInvariants
  public :: READ_TimeInvariants
  public :: WRITE_TimeInvariants

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_TimeInvariants ()
  ! --------------------------------------------------------------------
  ! Allocates memory for CoLM 1d [numpatch] variables
  ! --------------------------------------------------------------------

     use MOD_Precision
     USE MOD_Vars_Global
     use MOD_SPMD_Task
     use MOD_LandPatch, only: numpatch
     IMPLICIT NONE

     if (p_is_worker) then

        if (numpatch > 0) then

           allocate (patchclass           (numpatch))
           allocate (patchtype            (numpatch))
           allocate (patchmask            (numpatch))

           allocate (patchlonr            (numpatch))
           allocate (patchlatr            (numpatch))

           allocate (lakedepth            (numpatch))
           allocate (dz_lake      (nl_lake,numpatch))

           allocate (soil_s_v_alb         (numpatch))
           allocate (soil_d_v_alb         (numpatch))
           allocate (soil_s_n_alb         (numpatch))
           allocate (soil_d_n_alb         (numpatch))

           allocate (vf_quartz    (nl_soil,numpatch))
           allocate (vf_gravels   (nl_soil,numpatch))
           allocate (vf_om        (nl_soil,numpatch))
           allocate (vf_sand      (nl_soil,numpatch))
           allocate (wf_gravels   (nl_soil,numpatch))
           allocate (wf_sand      (nl_soil,numpatch))
           allocate (OM_density   (nl_soil,numpatch))
           allocate (BD_all       (nl_soil,numpatch))
           allocate (wfc          (nl_soil,numpatch))
           allocate (porsl        (nl_soil,numpatch))
           allocate (psi0         (nl_soil,numpatch))
           allocate (bsw          (nl_soil,numpatch))
#ifdef vanGenuchten_Mualem_SOIL_MODEL
           allocate (theta_r      (nl_soil,numpatch))
           allocate (alpha_vgm    (nl_soil,numpatch))
           allocate (L_vgm        (nl_soil,numpatch))
           allocate (n_vgm        (nl_soil,numpatch))
           allocate (sc_vgm       (nl_soil,numpatch))
           allocate (fc_vgm       (nl_soil,numpatch))
#endif
           allocate (hksati       (nl_soil,numpatch))
           allocate (csol         (nl_soil,numpatch))
           allocate (k_solids     (nl_soil,numpatch))
           allocate (dksatu       (nl_soil,numpatch))
           allocate (dksatf       (nl_soil,numpatch))
           allocate (dkdry        (nl_soil,numpatch))
           allocate (BA_alpha     (nl_soil,numpatch))
           allocate (BA_beta      (nl_soil,numpatch))
           allocate (htop                 (numpatch))
           allocate (hbot                 (numpatch))
           allocate (dbedrock             (numpatch))
           allocate (ibedrock             (numpatch))
     end if

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
     CALL allocate_PFTimeInvariants
#endif

#ifdef BGC
     CALL allocate_BGCTimeInvariants
#endif

#ifdef URBAN_MODEL
     CALL allocate_UrbanTimeInvariants
#endif

  end if

  END SUBROUTINE allocate_TimeInvariants

  !---------------------------------------
  SUBROUTINE READ_TimeInvariants (lc_year, casename, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use MOD_Namelist
     use MOD_SPMD_Task
     use MOD_NetCDFVector
     use MOD_NetCDFSerial
#ifdef RangeCheck
     USE MOD_RangeCheck
#endif
     USE MOD_LandPatch
     USE MOD_Vars_Global

     IMPLICIT NONE

     INTEGER         , intent(in) :: lc_year
     character(LEN=*), intent(in) :: casename
     character(LEN=*), intent(in) :: dir_restart

     ! Local variables
     character(LEN=256) :: file_restart, cyear

     write(cyear,'(i4.4)') lc_year
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_const' // '_lc' // trim(cyear) // '.nc'

     call ncio_read_vector (file_restart, 'patchclass',   landpatch, patchclass)          !
     call ncio_read_vector (file_restart, 'patchtype' ,   landpatch, patchtype )          !
     call ncio_read_vector (file_restart, 'patchmask' ,   landpatch, patchmask )          !

     call ncio_read_vector (file_restart, 'patchlonr' ,   landpatch, patchlonr )          !
     call ncio_read_vector (file_restart, 'patchlatr' ,   landpatch, patchlatr )          !

     call ncio_read_vector (file_restart, 'lakedepth',    landpatch, lakedepth)           !
     call ncio_read_vector (file_restart, 'dz_lake' ,     nl_lake, landpatch, dz_lake)    !

     call ncio_read_vector (file_restart, 'soil_s_v_alb', landpatch, soil_s_v_alb)        ! albedo of visible of the saturated soil
     call ncio_read_vector (file_restart, 'soil_d_v_alb', landpatch, soil_d_v_alb)        ! albedo of visible of the dry soil
     call ncio_read_vector (file_restart, 'soil_s_n_alb', landpatch, soil_s_n_alb)        ! albedo of near infrared of the saturated soil
     call ncio_read_vector (file_restart, 'soil_d_n_alb', landpatch, soil_d_n_alb)        ! albedo of near infrared of the dry soil

     call ncio_read_vector (file_restart, 'vf_quartz ',   nl_soil, landpatch, vf_quartz ) ! volumetric fraction of quartz within mineral soil
     call ncio_read_vector (file_restart, 'vf_gravels',   nl_soil, landpatch, vf_gravels) ! volumetric fraction of gravels
     call ncio_read_vector (file_restart, 'vf_om     ',   nl_soil, landpatch, vf_om     ) ! volumetric fraction of organic matter
     call ncio_read_vector (file_restart, 'vf_sand   ',   nl_soil, landpatch, vf_sand   ) ! volumetric fraction of sand
     call ncio_read_vector (file_restart, 'wf_gravels',   nl_soil, landpatch, wf_gravels) ! gravimetric fraction of gravels
     call ncio_read_vector (file_restart, 'wf_sand   ',   nl_soil, landpatch, wf_sand   ) ! gravimetric fraction of sand
     call ncio_read_vector (file_restart, 'OM_density',   nl_soil, landpatch, OM_density) ! OM density
     call ncio_read_vector (file_restart, 'BD_all    ',   nl_soil, landpatch, BD_all    ) ! bulk density of soil
     call ncio_read_vector (file_restart, 'wfc       ',   nl_soil, landpatch, wfc       ) ! field capacity
     call ncio_read_vector (file_restart, 'porsl  ' ,     nl_soil, landpatch, porsl     ) ! fraction of soil that is voids [-]
     call ncio_read_vector (file_restart, 'psi0   ' ,     nl_soil, landpatch, psi0      ) ! minimum soil suction [mm] (NOTE: "-" valued)
     call ncio_read_vector (file_restart, 'bsw    ' ,     nl_soil, landpatch, bsw       ) ! clapp and hornbereger "b" parameter [-]
#ifdef vanGenuchten_Mualem_SOIL_MODEL
     call ncio_read_vector (file_restart, 'theta_r  ' ,   nl_soil, landpatch, theta_r   ) ! residual moisture content [-]
     call ncio_read_vector (file_restart, 'alpha_vgm' ,   nl_soil, landpatch, alpha_vgm ) ! a parameter corresponding approximately to the inverse of the air-entry value
     call ncio_read_vector (file_restart, 'L_vgm    ' ,   nl_soil, landpatch, L_vgm     ) ! pore-connectivity parameter [dimensionless]
     call ncio_read_vector (file_restart, 'n_vgm    ' ,   nl_soil, landpatch, n_vgm     ) ! a shape parameter [dimensionless]
     call ncio_read_vector (file_restart, 'sc_vgm   ' ,   nl_soil, landpatch, sc_vgm    ) ! saturation at the air entry value in the classical vanGenuchten model [-]
     call ncio_read_vector (file_restart, 'fc_vgm   ' ,   nl_soil, landpatch, fc_vgm    ) ! a scaling factor by using air entry value in the Mualem model [-]
#endif
     call ncio_read_vector (file_restart, 'hksati ' ,     nl_soil, landpatch, hksati )    ! hydraulic conductivity at saturation [mm h2o/s]
     call ncio_read_vector (file_restart, 'csol   ' ,     nl_soil, landpatch, csol   )    ! heat capacity of soil solids [J/(m3 K)]
     call ncio_read_vector (file_restart, 'k_solids',     nl_soil, landpatch, k_solids)   ! thermal conductivity of soil solids [W/m-K]
     call ncio_read_vector (file_restart, 'dksatu ' ,     nl_soil, landpatch, dksatu )    ! thermal conductivity of unfrozen saturated soil [W/m-K]
     call ncio_read_vector (file_restart, 'dksatf ' ,     nl_soil, landpatch, dksatf )    ! thermal conductivity of frozen saturated soil [W/m-K]
     call ncio_read_vector (file_restart, 'dkdry  ' ,     nl_soil, landpatch, dkdry  )    ! thermal conductivity for dry soil  [W/(m-K)]
     call ncio_read_vector (file_restart, 'BA_alpha',     nl_soil, landpatch, BA_alpha)   ! alpha in Balland and Arp(2005) thermal conductivity scheme
     call ncio_read_vector (file_restart, 'BA_beta' ,     nl_soil, landpatch, BA_beta )   ! beta in Balland and Arp(2005) thermal conductivity scheme
     call ncio_read_vector (file_restart, 'htop' ,    landpatch, htop)                    !
     call ncio_read_vector (file_restart, 'hbot' ,    landpatch, hbot)                    !

     IF(DEF_USE_BEDROCK)THEN
        call ncio_read_vector (file_restart, 'debdrock' ,    landpatch, dbedrock)         !
        call ncio_read_vector (file_restart, 'ibedrock' ,    landpatch, ibedrock)         !
     ENDIF

     call ncio_read_bcast_serial (file_restart, 'zlnd  ', zlnd  ) ! roughness length for soil [m]
     call ncio_read_bcast_serial (file_restart, 'zsno  ', zsno  ) ! roughness length for snow [m]
     call ncio_read_bcast_serial (file_restart, 'csoilc', csoilc) ! drag coefficient for soil under canopy [-]
     call ncio_read_bcast_serial (file_restart, 'dewmx ', dewmx ) ! maximum dew
     call ncio_read_bcast_serial (file_restart, 'wtfact', wtfact) ! fraction of model area with high water table
     call ncio_read_bcast_serial (file_restart, 'capr  ', capr  ) ! tuning factor to turn first layer T into surface T
     call ncio_read_bcast_serial (file_restart, 'cnfac ', cnfac ) ! Crank Nicholson factor between 0 and 1
     call ncio_read_bcast_serial (file_restart, 'ssi   ', ssi   ) ! irreducible water saturation of snow
     call ncio_read_bcast_serial (file_restart, 'wimp  ', wimp  ) ! water impremeable if porosity less than wimp
     call ncio_read_bcast_serial (file_restart, 'pondmx', pondmx) ! ponding depth (mm)
     call ncio_read_bcast_serial (file_restart, 'smpmax', smpmax) ! wilting point potential in mm
     call ncio_read_bcast_serial (file_restart, 'smpmin', smpmin) ! restriction for min of soil poten. (mm)
     call ncio_read_bcast_serial (file_restart, 'trsmx0', trsmx0) ! max transpiration for moist soil+100% veg.  [mm/s]
     call ncio_read_bcast_serial (file_restart, 'tcrit ', tcrit ) ! critical temp. to determine rain or snow

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_pft_const' // '_lc' // trim(cyear) // '.nc'
     CALL READ_PFTimeInvariants (file_restart)
#endif

#if (defined BGC)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_bgc_const' // '_lc' // trim(cyear) // '.nc'
     CALL READ_BGCTimeInvariants (file_restart)
#endif

#if (defined URBAN_MODEL)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_urb_const' // '_lc' // trim(cyear) // '.nc'
     CALL READ_UrbanTimeInvariants (file_restart)
#endif

#ifdef RangeCheck
     call check_TimeInvariants ()
#endif

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif

     if (p_is_master) then
        write(*,'(A29)') 'Loading Time Invariants done.'
     end if

  end subroutine READ_TimeInvariants

  !---------------------------------------
  SUBROUTINE WRITE_TimeInvariants (lc_year, casename, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use MOD_Namelist, only : DEF_REST_COMPRESS_LEVEL, DEF_USE_BEDROCK
     use MOD_SPMD_Task
     use MOD_NetCDFSerial
     use MOD_NetCDFVector
     use MOD_LandPatch
     USE MOD_Vars_Global

     IMPLICIT NONE

     INTEGER         , intent(in) :: lc_year
     character(len=*), intent(in) :: casename
     character(len=*), intent(in) :: dir_restart

     ! Local Variables
     character(len=256) :: file_restart, cyear
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL

     write(cyear,'(i4.4)') lc_year
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_const' //'_lc'// trim(cyear) // '.nc'

     call ncio_create_file_vector (file_restart, landpatch)

     CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'soil', nl_soil)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'lake', nl_lake)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'band', 2)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'rtyp', 2)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'snow',     -maxsnl       )
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'snowp1',   -maxsnl+1     )
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'soilsnow', nl_soil-maxsnl)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'soil',     nl_soil)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'lake',     nl_lake)

     call ncio_write_vector (file_restart, 'patchclass', 'patch', landpatch, patchclass)                            !
     call ncio_write_vector (file_restart, 'patchtype' , 'patch', landpatch, patchtype )                            !
     call ncio_write_vector (file_restart, 'patchmask' , 'patch', landpatch, patchmask )                            !

     call ncio_write_vector (file_restart, 'patchlonr' , 'patch', landpatch, patchlonr )                            !
     call ncio_write_vector (file_restart, 'patchlatr' , 'patch', landpatch, patchlatr )                            !

     call ncio_write_vector (file_restart, 'lakedepth' , 'patch', landpatch, lakedepth , compress)                  !
     call ncio_write_vector (file_restart, 'dz_lake'   ,  'lake', nl_lake, 'patch', landpatch, dz_lake, compress)   !

     call ncio_write_vector (file_restart, 'soil_s_v_alb', 'patch', landpatch, soil_s_v_alb, compress)              ! albedo of visible of the saturated soil
     call ncio_write_vector (file_restart, 'soil_d_v_alb', 'patch', landpatch, soil_d_v_alb, compress)              ! albedo of visible of the dry soil
     call ncio_write_vector (file_restart, 'soil_s_n_alb', 'patch', landpatch, soil_s_n_alb, compress)              ! albedo of near infrared of the saturated soil
     call ncio_write_vector (file_restart, 'soil_d_n_alb', 'patch', landpatch, soil_d_n_alb, compress)              ! albedo of near infrared of the dry soil

     call ncio_write_vector (file_restart, 'vf_quartz ', 'soil', nl_soil, 'patch', landpatch, vf_quartz , compress) ! volumetric fraction of quartz within mineral soil
     call ncio_write_vector (file_restart, 'vf_gravels', 'soil', nl_soil, 'patch', landpatch, vf_gravels, compress) ! volumetric fraction of gravels
     call ncio_write_vector (file_restart, 'vf_om     ', 'soil', nl_soil, 'patch', landpatch, vf_om     , compress) ! volumetric fraction of organic matter
     call ncio_write_vector (file_restart, 'vf_sand   ', 'soil', nl_soil, 'patch', landpatch, vf_sand   , compress) ! volumetric fraction of sand
     call ncio_write_vector (file_restart, 'wf_gravels', 'soil', nl_soil, 'patch', landpatch, wf_gravels, compress) ! gravimetric fraction of gravels
     call ncio_write_vector (file_restart, 'wf_sand   ', 'soil', nl_soil, 'patch', landpatch, wf_sand   , compress) ! gravimetric fraction of sand
     call ncio_write_vector (file_restart, 'OM_density', 'soil', nl_soil, 'patch', landpatch, OM_density, compress) ! OM_density
     call ncio_write_vector (file_restart, 'BD_all    ', 'soil', nl_soil, 'patch', landpatch, BD_all    , compress) ! bulk density of soil
     call ncio_write_vector (file_restart, 'wfc       ', 'soil', nl_soil, 'patch', landpatch, wfc       , compress) ! field capacity
     call ncio_write_vector (file_restart, 'porsl     ', 'soil', nl_soil, 'patch', landpatch, porsl     , compress) ! fraction of soil that is voids [-]
     call ncio_write_vector (file_restart, 'psi0      ', 'soil', nl_soil, 'patch', landpatch, psi0      , compress) ! minimum soil suction [mm] (NOTE: "-" valued)
     call ncio_write_vector (file_restart, 'bsw       ', 'soil', nl_soil, 'patch', landpatch, bsw       , compress) ! clapp and hornbereger "b" parameter [-]

#ifdef vanGenuchten_Mualem_SOIL_MODEL
     call ncio_write_vector (file_restart, 'theta_r  ' , 'soil', nl_soil, 'patch', landpatch, theta_r   , compress) ! residual moisture content [-]
     call ncio_write_vector (file_restart, 'alpha_vgm' , 'soil', nl_soil, 'patch', landpatch, alpha_vgm , compress) ! a parameter corresponding approximately to the inverse of the air-entry value
     call ncio_write_vector (file_restart, 'L_vgm    ' , 'soil', nl_soil, 'patch', landpatch, L_vgm     , compress) ! pore-connectivity parameter [dimensionless]
     call ncio_write_vector (file_restart, 'n_vgm    ' , 'soil', nl_soil, 'patch', landpatch, n_vgm     , compress) ! a shape parameter [dimensionless]
     call ncio_write_vector (file_restart, 'sc_vgm   ' , 'soil', nl_soil, 'patch', landpatch, sc_vgm    , compress) ! saturation at the air entry value in the classical vanGenuchten model [-]
     call ncio_write_vector (file_restart, 'fc_vgm   ' , 'soil', nl_soil, 'patch', landpatch, fc_vgm    , compress) ! a scaling factor by using air entry value in the Mualem model [-]
#endif
     call ncio_write_vector (file_restart, 'hksati   ' , 'soil', nl_soil, 'patch', landpatch, hksati    , compress) ! hydraulic conductivity at saturation [mm h2o/s]
     call ncio_write_vector (file_restart, 'csol     ' , 'soil', nl_soil, 'patch', landpatch, csol      , compress) ! heat capacity of soil solids [J/(m3 K)]
     call ncio_write_vector (file_restart, 'k_solids ' , 'soil', nl_soil, 'patch', landpatch, k_solids  , compress) ! thermal conductivity of soil solids [W/m-K]
     call ncio_write_vector (file_restart, 'dksatu   ' , 'soil', nl_soil, 'patch', landpatch, dksatu    , compress) ! thermal conductivity of saturated soil [W/m-K]
     call ncio_write_vector (file_restart, 'dksatf   ' , 'soil', nl_soil, 'patch', landpatch, dksatf    , compress) ! thermal conductivity of saturated soil [W/m-K]
     call ncio_write_vector (file_restart, 'dkdry    ' , 'soil', nl_soil, 'patch', landpatch, dkdry     , compress) ! thermal conductivity for dry soil  [W/(m-K)]
     call ncio_write_vector (file_restart, 'BA_alpha ' , 'soil', nl_soil, 'patch', landpatch, BA_alpha  , compress) ! alpha in Balland and Arp(2005) thermal conductivity scheme
     call ncio_write_vector (file_restart, 'BA_beta  ' , 'soil', nl_soil, 'patch', landpatch, BA_beta   , compress) ! beta in Balland and Arp(2005) thermal conductivity scheme

     call ncio_write_vector (file_restart, 'htop' , 'patch', landpatch, htop)                                       !
     call ncio_write_vector (file_restart, 'hbot' , 'patch', landpatch, hbot)                                       !

     IF(DEF_USE_BEDROCK)THEN
        call ncio_write_vector (file_restart, 'debdrock' , 'patch', landpatch, dbedrock)                            !
        call ncio_write_vector (file_restart, 'ibedrock' , 'patch', landpatch, ibedrock)                            !
     ENDIF

     if (p_is_master) then

        call ncio_create_file (file_restart)

        call ncio_write_serial (file_restart, 'zlnd  ', zlnd  ) ! roughness length for soil [m]
        call ncio_write_serial (file_restart, 'zsno  ', zsno  ) ! roughness length for snow [m]
        call ncio_write_serial (file_restart, 'csoilc', csoilc) ! drag coefficient for soil under canopy [-]
        call ncio_write_serial (file_restart, 'dewmx ', dewmx ) ! maximum dew
        call ncio_write_serial (file_restart, 'wtfact', wtfact) ! fraction of model area with high water table
        call ncio_write_serial (file_restart, 'capr  ', capr  ) ! tuning factor to turn first layer T into surface T
        call ncio_write_serial (file_restart, 'cnfac ', cnfac ) ! Crank Nicholson factor between 0 and 1
        call ncio_write_serial (file_restart, 'ssi   ', ssi   ) ! irreducible water saturation of snow
        call ncio_write_serial (file_restart, 'wimp  ', wimp  ) ! water impremeable if porosity less than wimp
        call ncio_write_serial (file_restart, 'pondmx', pondmx) ! ponding depth (mm)
        call ncio_write_serial (file_restart, 'smpmax', smpmax) ! wilting point potential in mm
        call ncio_write_serial (file_restart, 'smpmin', smpmin) ! restriction for min of soil poten. (mm)
        call ncio_write_serial (file_restart, 'trsmx0', trsmx0) ! max transpiration for moist soil+100% veg.  [mm/s]
        call ncio_write_serial (file_restart, 'tcrit ', tcrit ) ! critical temp. to determine rain or snow

     end if

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_pft_const' //'_lc'// trim(cyear) // '.nc'
     CALL WRITE_PFTimeInvariants (file_restart)
#endif

#if (defined BGC)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_bgc_const' //'_lc'// trim(cyear) // '.nc'
     CALL WRITE_BGCTimeInvariants (file_restart)
#endif

#if (defined URBAN_MODEL)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_urb_const' //'_lc'// trim(cyear) // '.nc'
     CALL WRITE_UrbanTimeInvariants (file_restart)
#endif

   end subroutine WRITE_TimeInvariants

  SUBROUTINE deallocate_TimeInvariants ()

     use MOD_SPMD_Task
     use MOD_LandPatch, only: numpatch

     implicit none

     ! --------------------------------------------------
     ! Deallocates memory for CoLM 1d [numpatch] variables
     ! --------------------------------------------------

     if (p_is_worker) then

        if (numpatch > 0) then

           deallocate (patchclass     )
           deallocate (patchtype      )
           deallocate (patchmask      )

           deallocate (patchlonr      )
           deallocate (patchlatr      )

           deallocate (lakedepth      )
           deallocate (dz_lake        )

           deallocate (soil_s_v_alb   )
           deallocate (soil_d_v_alb   )
           deallocate (soil_s_n_alb   )
           deallocate (soil_d_n_alb   )

           deallocate (vf_quartz      )
           deallocate (vf_gravels     )
           deallocate (vf_om          )
           deallocate (vf_sand        )
           deallocate (wf_gravels     )
           deallocate (wf_sand        )
           deallocate (OM_density     )
           deallocate (BD_all         )
           deallocate (wfc            )
           deallocate (porsl          )
           deallocate (psi0           )
           deallocate (bsw            )
#ifdef vanGenuchten_Mualem_SOIL_MODEL
           deallocate (theta_r        )
           deallocate (alpha_vgm      )
           deallocate (L_vgm          )
           deallocate (n_vgm          )
           deallocate (sc_vgm         )
           deallocate (fc_vgm         )
#endif
           deallocate (hksati         )
           deallocate (csol           )
           deallocate (k_solids       )
           deallocate (dksatu         )
           deallocate (dksatf         )
           deallocate (dkdry          )
           deallocate (BA_alpha       )
           deallocate (BA_beta        )

           deallocate (htop           )
           deallocate (hbot           )

           deallocate (dbedrock       )
           deallocate (ibedrock       )

        end if
     end if

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
     CALL deallocate_PFTimeInvariants
#endif

#ifdef BGC
     CALL deallocate_BGCTimeInvariants
#endif

#ifdef URBAN_MODEL
     CALL deallocate_UrbanTimeInvariants
#endif
  END SUBROUTINE deallocate_TimeInvariants

#ifdef RangeCheck
   !---------------------------------------
   SUBROUTINE check_TimeInvariants ()

      use MOD_SPMD_Task
      use MOD_RangeCheck
      use MOD_Namelist, only : DEF_USE_BEDROCK

      IMPLICIT NONE

      if (p_is_master) then
         write(*,'(/,A29)') 'Checking Time Invariants ...'
      end if

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif

     call check_vector_data ('lakedepth    [m]     ', lakedepth   ) !
     call check_vector_data ('dz_lake      [m]     ', dz_lake     ) ! new lake scheme

     call check_vector_data ('soil_s_v_alb [-]     ', soil_s_v_alb) ! albedo of visible of the saturated soil
     call check_vector_data ('soil_d_v_alb [-]     ', soil_d_v_alb) ! albedo of visible of the dry soil
     call check_vector_data ('soil_s_n_alb [-]     ', soil_s_n_alb) ! albedo of near infrared of the saturated soil
     call check_vector_data ('soil_d_n_alb [-]     ', soil_d_n_alb) ! albedo of near infrared of the dry soil
     call check_vector_data ('vf_quartz    [m3/m3] ', vf_quartz   ) ! volumetric fraction of quartz within mineral soil
     call check_vector_data ('vf_gravels   [m3/m3] ', vf_gravels  ) ! volumetric fraction of gravels
     call check_vector_data ('vf_om        [m3/m3] ', vf_om       ) ! volumetric fraction of organic matter
     call check_vector_data ('vf_sand      [m3/m3] ', vf_sand     ) ! volumetric fraction of sand
     call check_vector_data ('wf_gravels   [kg/kg] ', wf_gravels  ) ! gravimetric fraction of gravels
     call check_vector_data ('wf_sand      [kg/kg] ', wf_sand     ) ! gravimetric fraction of sand
     call check_vector_data ('OM_density   [kg/m3] ', OM_density  ) ! OM density
     call check_vector_data ('BD_all       [kg/m3] ', BD_all      ) ! bulk density of soils
     call check_vector_data ('wfc          [m3/m3] ', wfc         ) ! field capacity
     call check_vector_data ('porsl        [m3/m3] ', porsl       ) ! fraction of soil that is voids [-]
     call check_vector_data ('psi0         [mm]    ', psi0        ) ! minimum soil suction [mm] (NOTE: "-" valued)
     call check_vector_data ('bsw          [-]     ', bsw         ) ! clapp and hornbereger "b" parameter [-]
#ifdef vanGenuchten_Mualem_SOIL_MODEL
     call check_vector_data ('theta_r      [m3/m3] ', theta_r     ) ! residual moisture content [-]
     call check_vector_data ('alpha_vgm    [-]     ', alpha_vgm   ) ! a parameter corresponding approximately to the inverse of the air-entry value
     call check_vector_data ('L_vgm        [-]     ', L_vgm       ) ! pore-connectivity parameter [dimensionless]
     call check_vector_data ('n_vgm        [-]     ', n_vgm       ) ! a shape parameter [dimensionless]
     call check_vector_data ('sc_vgm       [-]     ', sc_vgm      ) ! saturation at the air entry value in the classical vanGenuchten model [-]
     call check_vector_data ('fc_vgm       [-]     ', fc_vgm      ) ! a scaling factor by using air entry value in the Mualem model [-]
#endif
     call check_vector_data ('hksati       [mm/s]  ', hksati      ) ! hydraulic conductivity at saturation [mm h2o/s]
     call check_vector_data ('csol         [J/m3/K]', csol        ) ! heat capacity of soil solids [J/(m3 K)]
     call check_vector_data ('k_solids     [W/m/K] ', k_solids    ) ! thermal conductivity of soil solids [W/m-K]
     call check_vector_data ('dksatu       [W/m/K] ', dksatu      ) ! thermal conductivity of unfrozen saturated soil [W/m-K]
     call check_vector_data ('dksatf       [W/m/K] ', dksatf      ) ! thermal conductivity of frozen saturated soil [W/m-K]
     call check_vector_data ('dkdry        [W/m/K] ', dkdry       ) ! thermal conductivity for dry soil  [W/(m-K)]
     call check_vector_data ('BA_alpha     [-]     ', BA_alpha    ) ! alpha in Balland and Arp(2005) thermal conductivity scheme
     call check_vector_data ('BA_beta      [-]     ', BA_beta     ) ! beta in Balland and Arp(2005) thermal conductivity scheme

     call check_vector_data ('htop         [m]     ', htop        )
     call check_vector_data ('hbot         [m]     ', hbot        )

     IF(DEF_USE_BEDROCK)THEN
        call check_vector_data ('dbedrock     [m]     ', dbedrock    ) !
     ENDIF

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif

     if (p_is_master) then
        write(*,'(/,A)') 'Checking Constants ...'
        write(*,'(A,E20.10)') 'zlnd   [m]    ', zlnd   ! roughness length for soil [m]
        write(*,'(A,E20.10)') 'zsno   [m]    ', zsno   ! roughness length for snow [m]
        write(*,'(A,E20.10)') 'csoilc [-]    ', csoilc ! drag coefficient for soil under canopy [-]
        write(*,'(A,E20.10)') 'dewmx  [mm]   ', dewmx  ! maximum dew
        write(*,'(A,E20.10)') 'wtfact [-]    ', wtfact ! fraction of model area with high water table
        write(*,'(A,E20.10)') 'capr   [-]    ', capr   ! tuning factor to turn first layer T into surface T
        write(*,'(A,E20.10)') 'cnfac  [-]    ', cnfac  ! Crank Nicholson factor between 0 and 1
        write(*,'(A,E20.10)') 'ssi    [-]    ', ssi    ! irreducible water saturation of snow
        write(*,'(A,E20.10)') 'wimp   [m3/m3]', wimp   ! water impremeable if porosity less than wimp
        write(*,'(A,E20.10)') 'pondmx [mm]   ', pondmx ! ponding depth (mm)
        write(*,'(A,E20.10)') 'smpmax [mm]   ', smpmax ! wilting point potential in mm
        write(*,'(A,E20.10)') 'smpmin [mm]   ', smpmin ! restriction for min of soil poten. (mm)
        write(*,'(A,E20.10)') 'trsmx0 [mm/s] ', trsmx0 ! max transpiration for moist soil+100% veg.  [mm/s]
        write(*,'(A,E20.10)') 'tcrit  [K]    ', tcrit  ! critical temp. to determine rain or snow
     end if

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
     CALL check_PFTimeInvariants
#endif

#ifdef BGC
     CALL check_BGCTimeInvariants
#endif

   end subroutine check_TimeInvariants
#endif

END MODULE MOD_Vars_TimeInvariants
! ---------- EOP ------------
