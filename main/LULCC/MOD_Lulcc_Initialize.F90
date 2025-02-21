#include <define.h>
#ifdef LULCC
MODULE MOD_Lulcc_Initialize

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: LulccInitialize

CONTAINS

   SUBROUTINE LulccInitialize (casename,dir_landdata,dir_restart,&
                               idate,greenwich)

!-----------------------------------------------------------------------
!
! !DESCRIPTION:
!  Initialization routine for Land-use-Land-cover-change (Lulcc) case
!
!  Created by Hua Yuan, 04/08/2022
!
! !REVISIONS:
!  08/2023, Wenzong Dong: Porting to MPI version and share the same code
!           with MOD_Initialize:initialize()
!
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_Mesh
   USE MOD_LandElm
   USE MOD_LandPatch
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
#endif
   USE MOD_LandUrban
   USE MOD_Const_LC
   USE MOD_Const_PFT
   USE MOD_TimeManager
   USE MOD_Lulcc_Vars_TimeInvariants
   USE MOD_Lulcc_Vars_TimeVariables
   USE MOD_SrfdataRestart
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables
   USE MOD_Initialize
#ifdef SrfdataDiag
   USE MOD_SrfdataDiag, only: gdiag, srfdata_diag_init
#endif

   IMPLICIT NONE

!-------------------------- Dummy Arguments ----------------------------
   character(len=*), intent(in) :: casename      ! case name
   character(len=*), intent(in) :: dir_landdata
   character(len=*), intent(in) :: dir_restart

   integer, intent(inout) :: idate(3)   ! year, julian day, seconds of the starting time
   logical, intent(in)    :: greenwich  ! true: greenwich time, false: local time

!-------------------------- Local Variables ----------------------------
   integer :: year, jday

!-----------------------------------------------------------------------

      ! initial time of model run
      ! ............................
      CALL adj2begin(idate)

      year = idate(1)
      jday = idate(2)

      CALL Init_GlobalVars
      CAll Init_LC_Const
      CAll Init_PFT_Const

      ! deallocate pixelset and mesh data of previous year
      CALL mesh_free_mem
      CALL landelm%forc_free_mem
      CALL landpatch%forc_free_mem
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL landpft%forc_free_mem
#endif
#ifdef URBAN_MODEL
      CALL landurban%forc_free_mem
#endif

      ! load pixelset and mesh data of next year
      ! CALL pixel%load_from_file  (dir_landdata)
      ! CALL gblock%load_from_file (dir_landdata)
      CALL mesh_load_from_file     (dir_landdata, year)
      CALL pixelset_load_from_file (dir_landdata, 'landelm'  , landelm  , numelm  , year)

      ! load CATCHMENT of next year
#ifdef CATCHMENT
      CALL pixelset_load_from_file (dir_landdata, 'landhru'  , landhru  , numhru  , year)
#endif

      ! load landpatch data of next year
      CALL pixelset_load_from_file (dir_landdata, 'landpatch', landpatch, numpatch, year)

      ! load pft data of PFT/PC of next year
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL pixelset_load_from_file (dir_landdata, 'landpft'  , landpft  , numpft  , year)
      CALL map_patch_to_pft
#endif

      ! load urban data of next year
#ifdef URBAN_MODEL
      CALL pixelset_load_from_file (dir_landdata, 'landurban', landurban, numurban, year)
      CALL map_patch_to_urban
#endif

      ! initialize for data associated with land element
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CALL elm_vector_init ()
#ifdef CATCHMENT
      CALL hru_vector_init ()
#endif
#endif

      ! build element subfraction of next year which it's needed in the MOD_Lulcc_TransferTrace
      IF (p_is_worker) THEN
         CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
      ENDIF

      ! initialize for SrfdataDiag, it is needed in the MOD_Lulcc_TransferTrace for outputing transfer_matrix
#ifdef SrfdataDiag
#ifdef GRIDBASED
      CALL init_gridbased_mesh_grid ()
      CALL gdiag%define_by_copy (gridmesh)
#else
      CALL gdiag%define_by_ndims(3600,1800)
#endif
      CALL srfdata_diag_init (dir_landdata)
#endif

      ! --------------------------------------------------------------------
      ! Deallocates memory for CoLM 1d [numpatch] variables
      ! --------------------------------------------------------------------
      CALL deallocate_TimeInvariants
      CALL deallocate_TimeVariables

      ! initialize all state variables of next year
      CALL initialize (casename,dir_landdata,dir_restart,&
                       idate,year,greenwich,lulcc_call=.true.)

   END SUBROUTINE LulccInitialize

END MODULE MOD_Lulcc_Initialize
#endif
! ---------- EOP ------------
