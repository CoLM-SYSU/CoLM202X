#include <define.h>
#ifdef LULCC
MODULE MOD_Lulcc_Initialize

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: LulccInitialize

!-----------------------------------------------------------------------
   CONTAINS
!-----------------------------------------------------------------------

 SUBROUTINE LulccInitialize (casename,dir_landdata,dir_restart,&
                             idate,greenwich)

! ======================================================================
!
! Initialization routine for Land-use-Land-cover-change (Lulcc) case
!
! ======================================================================

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_Mesh
   USE MOD_LandElm
   USE MOD_LandPatch
   USE MOD_Const_LC
   USE MOD_Const_PFT
   use MOD_TimeManager
   USE MOD_Lulcc_Vars_TimeInvariants
   USE MOD_Lulcc_Vars_TimeVariables
   USE MOD_SrfdataRestart
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables
   USE MOD_Initialize

   IMPLICIT NONE

   ! ----------------------------------------------------------------------
   character(len=*), intent(in) :: casename      ! case name
   character(len=*), intent(in) :: dir_landdata
   character(len=*), intent(in) :: dir_restart

   integer, intent(inout) :: idate(3)   ! year, julian day, seconds of the starting time
   logical, intent(in)    :: greenwich  ! true: greenwich time, false: local time

   ! local vars
   INTEGER :: year, jday
   ! ----------------------------------------------------------------------

   ! initial time of model run
   ! ............................
   CALL adj2begin(idate)

   year = idate(1)
   jday = idate(2)

   CALL Init_GlobalVars
   CAll Init_LC_Const
   CAll Init_PFT_Const

   ! deallocate pixelset and mesh data of previous
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
   ! call pixel%load_from_file  (dir_landdata)
   ! call gblock%load_from_file (dir_landdata)
   call mesh_load_from_file     (dir_landdata, year)
   CALL pixelset_load_from_file (dir_landdata, 'landelm'  , landelm  , numelm  , year)

#ifdef CATCHMENT
   CALL pixelset_load_from_file (dir_landdata, 'landhru'  , landhru  , numhru  , year)
#endif

   call pixelset_load_from_file (dir_landdata, 'landpatch', landpatch, numpatch, year)

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   call pixelset_load_from_file (dir_landdata, 'landpft'  , landpft  , numpft  , year)
   CALL map_patch_to_pft
#endif

#ifdef URBAN_MODEL
   CALL pixelset_load_from_file (dir_landdata, 'landurban', landurban, numurban, year)
   CALL map_patch_to_urban
#endif

#if (defined UNSTRUCTURED || defined CATCHMENT)
   CALL elm_vector_init ()
#ifdef CATCHMENT
   CALL hru_vector_init ()
#endif
#endif

   ! --------------------------------------------------------------------
   ! Deallocates memory for CoLM 1d [numpatch] variables
   ! --------------------------------------------------------------------
   CALL deallocate_TimeInvariants
   CALL deallocate_TimeVariables

   CALL initialize (casename,dir_landdata,dir_restart,&
                    idate,year,greenwich,lulcc_call=.true.)

 END SUBROUTINE LulccInitialize

END MODULE MOD_Lulcc_Initialize
#endif
