#include <define.h>

MODULE MOD_TimeInvariants 
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

use precision
IMPLICIT NONE
SAVE
! -----------------------------------------------------------------
! surface classification and soil information
  INTEGER,  allocatable :: patchclass     (:)  !index of land cover type of the patches at the fraction > 0
  INTEGER,  allocatable :: patchtype      (:)  !land water type

  REAL(r8), allocatable :: patchlatr      (:)  !latitude in radians
  REAL(r8), allocatable :: patchlonr      (:)  !longitude in radians

  REAL(r8), allocatable :: lakedepth      (:)  !lake depth
  REAL(r8), allocatable :: dz_lake      (:,:)  !new lake scheme

  REAL(r8), allocatable :: soil_s_v_alb   (:)  !albedo of visible of the saturated soil
  REAL(r8), allocatable :: soil_d_v_alb   (:)  !albedo of visible of the dry soil
  REAL(r8), allocatable :: soil_s_n_alb   (:)  !albedo of near infrared of the saturated soil
  REAL(r8), allocatable :: soil_d_n_alb   (:)  !albedo of near infrared of the dry soil
  REAL(r8), allocatable :: porsl        (:,:)  !fraction of soil that is voids [-]
  REAL(r8), allocatable :: psi0         (:,:)  !minimum soil suction [mm] (NOTE: "-" valued)
#ifdef Campbell_SOIL_MODEL
  REAL(r8), allocatable :: bsw          (:,:)  !clapp and hornbereger "b" parameter [-]
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
  REAL(r8), allocatable :: theta_r      (:,:)
  REAL(r8), allocatable :: alpha_vgm    (:,:)
  REAL(r8), allocatable :: L_vgm        (:,:)
  REAL(r8), allocatable :: n_vgm        (:,:)
  REAL(r8), allocatable :: sc_vgm       (:,:)
  REAL(r8), allocatable :: fc_vgm       (:,:)
#endif
  REAL(r8), allocatable :: hksati       (:,:)  !hydraulic conductivity at saturation [mm h2o/s]
  REAL(r8), allocatable :: csol         (:,:)  !heat capacity of soil solids [J/(m3 K)]
  REAL(r8), allocatable :: dksatu       (:,:)  !thermal conductivity of saturated soil [W/m-K]
  REAL(r8), allocatable :: dkdry        (:,:)  !thermal conductivity for dry soil  [W/(m-K)]

  REAL(r8), allocatable :: htop           (:)  !canopy top height [m]
  REAL(r8), allocatable :: hbot           (:)  !canopy bottom height [m]

#ifdef USE_DEPTH_TO_BEDROCK
  real(r8), allocatable :: dbedrock       (:)  ! depth to bedrock
  integer , allocatable :: ibedrock       (:)  ! bedrock level
#endif

  REAL(r8) :: zlnd         ! roughness length for soil [m]
  REAL(r8) :: zsno         ! roughness length for snow [m]
  REAL(r8) :: csoilc       ! drag coefficient for soil under canopy [-]
  REAL(r8) :: dewmx        ! maximum dew
  REAL(r8) :: wtfact       ! fraction of model area with high water table
  REAL(r8) :: capr         ! tuning factor to turn first layer T into surface T
  REAL(r8) :: cnfac        ! Crank Nicholson factor between 0 and 1
  REAL(r8) :: ssi          ! irreducible water saturation of snow
  REAL(r8) :: wimp         ! water impremeable if porosity less than wimp
  REAL(r8) :: pondmx       ! ponding depth (mm)
  REAL(r8) :: smpmax       ! wilting point potential in mm
  REAL(r8) :: smpmin       ! restriction for min of soil poten. (mm)
  REAL(r8) :: trsmx0       ! max transpiration for moist soil+100% veg.  [mm/s]
  REAL(r8) :: tcrit        ! critical temp. to determine rain or snow

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
  ! Allocates memory for CLM 1d [numpatch] variables
  ! --------------------------------------------------------------------

     use precision
     USE GlobalVars
     use spmd_task
     use mod_landpatch, only : numpatch
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeInvars
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeInvars
#endif
     IMPLICIT NONE

  if (p_is_worker) then

     if (numpatch > 0) then

        allocate (patchclass           (numpatch))
        allocate (patchtype            (numpatch))
        
        allocate (patchlonr            (numpatch))
        allocate (patchlatr            (numpatch))

        allocate (lakedepth            (numpatch))
        allocate (dz_lake      (nl_lake,numpatch))

        allocate (soil_s_v_alb         (numpatch))
        allocate (soil_d_v_alb         (numpatch))
        allocate (soil_s_n_alb         (numpatch))
        allocate (soil_d_n_alb         (numpatch))

        allocate (porsl        (nl_soil,numpatch))
        allocate (psi0         (nl_soil,numpatch))
#ifdef Campbell_SOIL_MODEL
        allocate (bsw          (nl_soil,numpatch))
#endif
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
        allocate (dksatu       (nl_soil,numpatch))
        allocate (dkdry        (nl_soil,numpatch))
     
        allocate (htop                 (numpatch))
        allocate (hbot                 (numpatch))

#ifdef USE_DEPTH_TO_BEDROCK
        allocate (dbedrock             (numpatch))
        allocate (ibedrock             (numpatch))
#endif

     end if

#ifdef PFT_CLASSIFICATION
     CALL allocate_PFTimeInvars
#endif

#ifdef PC_CLASSIFICATION
     CALL allocate_PCTimeInvars
#endif

  end if

  END SUBROUTINE allocate_TimeInvariants

  !---------------------------------------
  SUBROUTINE READ_TimeInvariants (casename, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use mod_namelist
     use spmd_task
     use ncio_vector
     use ncio_serial
     use mod_colm_debug
     USE mod_landpatch
     USE GlobalVars
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeInvars
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeInvars
#endif

     IMPLICIT NONE

     character(LEN=*), intent(in) :: casename
     character(LEN=*), intent(in) :: dir_restart
  
     ! Local variables
     character(LEN=256) :: file_restart

     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_const.nc'

     call ncio_read_vector (file_restart, 'patchclass',   landpatch, patchclass) !
     call ncio_read_vector (file_restart, 'patchtype' ,   landpatch, patchtype ) !
     
     call ncio_read_vector (file_restart, 'patchlonr' ,   landpatch, patchlonr ) !
     call ncio_read_vector (file_restart, 'patchlatr' ,   landpatch, patchlatr ) !
     
     call ncio_read_vector (file_restart, 'lakedepth',    landpatch, lakedepth) !
     call ncio_read_vector (file_restart, 'dz_lake' ,     nl_lake, landpatch, dz_lake) !
     
     call ncio_read_vector (file_restart, 'soil_s_v_alb', landpatch, soil_s_v_alb) ! albedo of visible of the saturated soil
     call ncio_read_vector (file_restart, 'soil_d_v_alb', landpatch, soil_d_v_alb) ! albedo of visible of the dry soil
     call ncio_read_vector (file_restart, 'soil_s_n_alb', landpatch, soil_s_n_alb) ! albedo of near infrared of the saturated soil
     call ncio_read_vector (file_restart, 'soil_d_n_alb', landpatch, soil_d_n_alb) ! albedo of near infrared of the dry soil

     call ncio_read_vector (file_restart, 'porsl  ' ,     nl_soil, landpatch, porsl  ) ! fraction of soil that is voids [-]
     call ncio_read_vector (file_restart, 'psi0   ' ,     nl_soil, landpatch, psi0   ) ! minimum soil suction [mm] (NOTE: "-" valued)
#ifdef Campbell_SOIL_MODEL                                         
     call ncio_read_vector (file_restart, 'bsw    ' ,     nl_soil, landpatch, bsw    ) ! clapp and hornbereger "b" parameter [-]
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
     call ncio_read_vector (file_restart, 'theta_r  ' ,   nl_soil, landpatch, theta_r   ) 
     call ncio_read_vector (file_restart, 'alpha_vgm' ,   nl_soil, landpatch, alpha_vgm ) 
     call ncio_read_vector (file_restart, 'L_vgm    ' ,   nl_soil, landpatch, L_vgm     ) 
     call ncio_read_vector (file_restart, 'n_vgm    ' ,   nl_soil, landpatch, n_vgm     ) 
     call ncio_read_vector (file_restart, 'sc_vgm   ' ,   nl_soil, landpatch, sc_vgm    ) 
     call ncio_read_vector (file_restart, 'fc_vgm   ' ,   nl_soil, landpatch, fc_vgm    ) 
#endif                                                             
     call ncio_read_vector (file_restart, 'hksati ' ,     nl_soil, landpatch, hksati ) ! hydraulic conductivity at saturation [mm h2o/s]
     call ncio_read_vector (file_restart, 'csol   ' ,     nl_soil, landpatch, csol   ) ! heat capacity of soil solids [J/(m3 K)]
     call ncio_read_vector (file_restart, 'dksatu ' ,     nl_soil, landpatch, dksatu ) ! thermal conductivity of saturated soil [W/m-K]
     call ncio_read_vector (file_restart, 'dkdry  ' ,     nl_soil, landpatch, dkdry  ) ! thermal conductivity for dry soil  [W/(m-K)]
     
     call ncio_read_vector (file_restart, 'htop' ,    landpatch, htop) !
     call ncio_read_vector (file_restart, 'hbot' ,    landpatch, hbot) !

#ifdef USE_DEPTH_TO_BEDROCK
     call ncio_read_vector (file_restart, 'debdrock' ,    landpatch, dbedrock) !
     call ncio_read_vector (file_restart, 'ibedrock' ,    landpatch, ibedrock) !
#endif

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

#if (defined PFT_CLASSIFICATION)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_pft_const.nc'
     CALL READ_PFTimeInvars (file_restart)
#endif 

#if (defined PC_CLASSIFICATION)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_pc_const.nc'
     CALL READ_PCTimeInvars (file_restart)
#endif 

#ifdef CLMDEBUG 
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
  SUBROUTINE WRITE_TimeInvariants (casename, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use mod_namelist, only : DEF_REST_COMPRESS_LEVEL 
     use spmd_task
     use ncio_serial
     use ncio_vector
     use mod_landpatch
     USE GlobalVars
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeInvars
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeInvars
#endif

     IMPLICIT NONE

     character(len=*), intent(in) :: casename
     character(len=*), intent(in) :: dir_restart
     
     ! Local Variables
     character(len=256) :: file_restart
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL 

     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_const.nc'

     call ncio_create_file_vector (file_restart, landpatch)

     CALL ncio_define_pixelset_dimension (file_restart, landpatch)
     CALL ncio_define_dimension_vector (file_restart, 'soil', nl_soil)
     CALL ncio_define_dimension_vector (file_restart, 'lake', nl_lake)
     CALL ncio_define_dimension_vector (file_restart, 'band',   2)
     CALL ncio_define_dimension_vector (file_restart, 'wetdry', 2)
     
     call ncio_write_vector (file_restart, 'patchclass', 'vector', landpatch, patchclass) !
     call ncio_write_vector (file_restart, 'patchtype' , 'vector', landpatch, patchtype ) !
     
     call ncio_write_vector (file_restart, 'patchlonr' , 'vector', landpatch, patchlonr ) !
     call ncio_write_vector (file_restart, 'patchlatr' , 'vector', landpatch, patchlatr ) !
     
     call ncio_write_vector (file_restart, 'lakedepth', 'vector', landpatch, lakedepth ,      compress) !
     call ncio_write_vector (file_restart, 'dz_lake' ,  'lake', nl_lake, 'vector', landpatch, dz_lake, compress) !
     
     call ncio_write_vector (file_restart, 'soil_s_v_alb', 'vector', landpatch, soil_s_v_alb, compress) ! albedo of visible of the saturated soil
     call ncio_write_vector (file_restart, 'soil_d_v_alb', 'vector', landpatch, soil_d_v_alb, compress) ! albedo of visible of the dry soil
     call ncio_write_vector (file_restart, 'soil_s_n_alb', 'vector', landpatch, soil_s_n_alb, compress) ! albedo of near infrared of the saturated soil
     call ncio_write_vector (file_restart, 'soil_d_n_alb', 'vector', landpatch, soil_d_n_alb, compress) ! albedo of near infrared of the dry soil

     call ncio_write_vector (file_restart, 'porsl  ' , 'soil', nl_soil, 'vector', landpatch, porsl , compress) ! fraction of soil that is voids [-]
     call ncio_write_vector (file_restart, 'psi0   ' , 'soil', nl_soil, 'vector', landpatch, psi0  , compress) ! minimum soil suction [mm] (NOTE: "-" valued)
#ifdef Campbell_SOIL_MODEL
     call ncio_write_vector (file_restart, 'bsw    ' , 'soil', nl_soil, 'vector', landpatch, bsw   , compress) ! clapp and hornbereger "b" parameter [-]
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
     call ncio_write_vector (file_restart, 'theta_r  ' , 'soil', nl_soil, 'vector', landpatch, theta_r  , compress) 
     call ncio_write_vector (file_restart, 'alpha_vgm' , 'soil', nl_soil, 'vector', landpatch, alpha_vgm, compress) 
     call ncio_write_vector (file_restart, 'L_vgm    ' , 'soil', nl_soil, 'vector', landpatch, L_vgm    , compress) 
     call ncio_write_vector (file_restart, 'n_vgm    ' , 'soil', nl_soil, 'vector', landpatch, n_vgm    , compress) 
     call ncio_write_vector (file_restart, 'sc_vgm   ' , 'soil', nl_soil, 'vector', landpatch, sc_vgm   , compress) 
     call ncio_write_vector (file_restart, 'fc_vgm   ' , 'soil', nl_soil, 'vector', landpatch, fc_vgm   , compress) 
#endif
     call ncio_write_vector (file_restart, 'hksati ' , 'soil', nl_soil, 'vector', landpatch, hksati, compress) ! hydraulic conductivity at saturation [mm h2o/s]
     call ncio_write_vector (file_restart, 'csol   ' , 'soil', nl_soil, 'vector', landpatch, csol  , compress) ! heat capacity of soil solids [J/(m3 K)]
     call ncio_write_vector (file_restart, 'dksatu ' , 'soil', nl_soil, 'vector', landpatch, dksatu, compress) ! thermal conductivity of saturated soil [W/m-K]
     call ncio_write_vector (file_restart, 'dkdry  ' , 'soil', nl_soil, 'vector', landpatch, dkdry , compress) ! thermal conductivity for dry soil  [W/(m-K)]

     call ncio_write_vector (file_restart, 'htop' , 'vector', landpatch, htop) !
     call ncio_write_vector (file_restart, 'hbot' , 'vector', landpatch, hbot) !

#ifdef USE_DEPTH_TO_BEDROCK
     call ncio_write_vector (file_restart, 'debdrock' , 'vector', landpatch, dbedrock) !
     call ncio_write_vector (file_restart, 'ibedrock' , 'vector', landpatch, ibedrock) !
#endif

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

#if (defined PFT_CLASSIFICATION)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_pft_const.nc'
     CALL WRITE_PFTimeInvars (file_restart)
#endif 

#if (defined PC_CLASSIFICATION)
     file_restart = trim(dir_restart) // '/' // trim(casename) //'_restart_pc_const.nc'
     CALL WRITE_PCTimeInvars (file_restart)
#endif 

   end subroutine WRITE_TimeInvariants 

  SUBROUTINE deallocate_TimeInvariants ()

     use spmd_task
     use mod_landpatch, only : numpatch
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeInvars
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeInvars
#endif
     implicit none

     ! --------------------------------------------------
     ! Deallocates memory for CLM 1d [numpatch] variables
     ! --------------------------------------------------

     if (p_is_worker) then

        if (numpatch > 0) then

           deallocate (patchclass)
           deallocate (patchtype )
           
           deallocate (patchlonr )
           deallocate (patchlatr )

           deallocate (lakedepth)
           deallocate (dz_lake  )

           deallocate (soil_s_v_alb)
           deallocate (soil_d_v_alb)
           deallocate (soil_s_n_alb)
           deallocate (soil_d_n_alb)

           deallocate (porsl  )
           deallocate (psi0   )
#ifdef Campbell_SOIL_MODEL
           deallocate (bsw    )
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
           deallocate (theta_r  )
           deallocate (alpha_vgm)
           deallocate (L_vgm    )
           deallocate (n_vgm    )
           deallocate (sc_vgm   )
           deallocate (fc_vgm   )
#endif
           deallocate (hksati )
           deallocate (csol   )
           deallocate (dksatu )
           deallocate (dkdry  )

           deallocate (htop)
           deallocate (hbot)
          
#ifdef USE_DEPTH_TO_BEDROCK
           deallocate (dbedrock)
           deallocate (ibedrock)
#endif

        end if
     end if

#ifdef PFT_CLASSIFICATION
     CALL deallocate_PFTimeInvars
#endif

#ifdef PC_CLASSIFICATION
     CALL deallocate_PCTimeInvars
#endif

  END SUBROUTINE deallocate_TimeInvariants

#ifdef CLMDEBUG
   !---------------------------------------
   SUBROUTINE check_TimeInvariants ()

      use spmd_task
      use mod_colm_debug
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeInvars
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeInvars
#endif

      IMPLICIT NONE
      
      if (p_is_master) then
         write(*,'(/,A29)') 'Checking Time Invariantes ...'
      end if

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif

      call check_vector_data ('lakedepth   ', lakedepth   ) !
      call check_vector_data ('dz_lake     ', dz_lake     ) ! new lake scheme
      
      call check_vector_data ('soil_s_v_alb', soil_s_v_alb) ! albedo of visible of the saturated soil
      call check_vector_data ('soil_d_v_alb', soil_d_v_alb) ! albedo of visible of the dry soil
      call check_vector_data ('soil_s_n_alb', soil_s_n_alb) ! albedo of near infrared of the saturated soil
      call check_vector_data ('soil_d_n_alb', soil_d_n_alb) ! albedo of near infrared of the dry soil
      call check_vector_data ('porsl       ', porsl       ) ! fraction of soil that is voids [-]
      call check_vector_data ('psi0        ', psi0        ) ! minimum soil suction [mm] (NOTE: "-" valued)
#ifdef Campbell_SOIL_MODEL
      call check_vector_data ('bsw         ', bsw         ) ! clapp and hornbereger "b" parameter [-]
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      call check_vector_data ('theta_r     ', theta_r     ) 
      call check_vector_data ('alpha_vgm   ', alpha_vgm   ) 
      call check_vector_data ('L_vgm       ', L_vgm       ) 
      call check_vector_data ('n_vgm       ', n_vgm       ) 
      call check_vector_data ('sc_vgm      ', sc_vgm      ) 
      call check_vector_data ('fc_vgm      ', fc_vgm      ) 
#endif
      call check_vector_data ('hksati      ', hksati      ) ! hydraulic conductivity at saturation [mm h2o/s]
      call check_vector_data ('csol        ', csol        ) ! heat capacity of soil solids [J/(m3 K)]
      call check_vector_data ('dksatu      ', dksatu      ) ! thermal conductivity of saturated soil [W/m-K]
      call check_vector_data ('dkdry       ', dkdry       ) ! thermal conductivity for dry soil  [W/(m-K)]
      
      call check_vector_data ('htop        ', htop        ) 
      call check_vector_data ('hbot        ', hbot        ) 

#ifdef USE_DEPTH_TO_BEDROCK
      call check_vector_data ('dbedrock    ', dbedrock    ) !
#endif

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif
     
      if (p_is_master) then
         write(*,'(A7,E20.10)') 'zlnd  ', zlnd   ! roughness length for soil [m]
         write(*,'(A7,E20.10)') 'zsno  ', zsno   ! roughness length for snow [m]
         write(*,'(A7,E20.10)') 'csoilc', csoilc ! drag coefficient for soil under canopy [-]
         write(*,'(A7,E20.10)') 'dewmx ', dewmx  ! maximum dew
         write(*,'(A7,E20.10)') 'wtfact', wtfact ! fraction of model area with high water table
         write(*,'(A7,E20.10)') 'capr  ', capr   ! tuning factor to turn first layer T into surface T
         write(*,'(A7,E20.10)') 'cnfac ', cnfac  ! Crank Nicholson factor between 0 and 1
         write(*,'(A7,E20.10)') 'ssi   ', ssi    ! irreducible water saturation of snow
         write(*,'(A7,E20.10)') 'wimp  ', wimp   ! water impremeable if porosity less than wimp
         write(*,'(A7,E20.10)') 'pondmx', pondmx ! ponding depth (mm)
         write(*,'(A7,E20.10)') 'smpmax', smpmax ! wilting point potential in mm
         write(*,'(A7,E20.10)') 'smpmin', smpmin ! restriction for min of soil poten. (mm)
         write(*,'(A7,E20.10)') 'trsmx0', trsmx0 ! max transpiration for moist soil+100% veg.  [mm/s]
         write(*,'(A7,E20.10)') 'tcrit ', tcrit  ! critical temp. to determine rain or snow
      end if

#ifdef PFT_CLASSIFICATION
     CALL check_PFTimeInvars
#endif

#ifdef PC_CLASSIFICATION
     CALL check_PCTimeInvars
#endif

   end subroutine check_TimeInvariants 
#endif

END MODULE MOD_TimeInvariants
