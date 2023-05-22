#include <define.h>

MODULE MOD_OzoneData

   USE mod_grid
   USE mod_data_type
   USE mod_mapping_grid2pset
   use MOD_Vars_1DForcing, only: forc_ozone
   IMPLICIT NONE

   CHARACTER(len=256) :: file_ozone
   TYPE(grid_type) :: grid_ozone

   TYPE(block_data_real8_2d) :: f_ozone
   
   type (mapping_grid2pset_type) :: mg2p_ozone

CONTAINS

   ! ----------
   SUBROUTINE init_ozone_data (time, idate)
      
      USE spmd_task
      USE mod_namelist
      USE timemanager
      USE mod_grid
      USE ncio_serial
      USE ncio_block
      USE mod_landpatch
      USE mod_colm_debug
      IMPLICIT NONE
      
      type(timestamp), intent(in) :: time
      integer,         intent(in) :: idate(3)

      ! Local Variables
      REAL(r8), allocatable :: lat(:), lon(:)
      INTEGER :: itime
      INTEGER :: iyear, month, mday
      CHARACTER(LEN=8) :: syear, smonth

      call julian2monthday(idate(1),idate(2),month,mday)
      iyear = idate(1)
      if(idate(1) .lt. 2013)iyear = 2013
      if(idate(1) .gt. 2021)iyear = 2021
      write(syear,"(I4.4)")  iyear
      write(smonth,"(I2.2)") month
      file_ozone = trim(DEF_dir_rawdata) // '/Ozone/China/'//trim(syear)//trim(smonth)//'_O3_v2.nc'

      CALL ncio_read_bcast_serial (file_ozone, 'latitude', lat)
      CALL ncio_read_bcast_serial (file_ozone, 'longitude', lon)

      CALL grid_ozone%define_by_center (lat, lon)
      
      CALL allocate_block_data (grid_ozone, f_ozone)  
      
      call mg2p_ozone%build (grid_ozone, landpatch)

      itime = mday

      CALL ncio_read_block_time (file_ozone, 'O3', grid_ozone, itime, f_ozone)
#ifdef CoLMDEBUG
      CALL check_block_data ('Ozone', f_ozone)
#endif

!      IF (p_is_worker) THEN
!         IF (numpatch > 0) THEN
!            allocate (lnfm (numpatch))
!         ENDIF
!      ENDIF

   END SUBROUTINE init_ozone_data 

   ! ----------
   SUBROUTINE update_ozone_data (time, deltim)
      
      USE timemanager
      USE mod_namelist
      USE ncio_block
      USE mod_colm_debug
      IMPLICIT NONE
      
      type(timestamp), intent(in) :: time
      REAL(r8), intent(in) :: deltim

      ! Local Variables
      type(timestamp) :: time_next
      INTEGER :: month, mday
      INTEGER :: iyear, imonth, imonth_next, iday, iday_next
      CHARACTER(LEN=8) :: syear, smonth

      call julian2monthday(time%year,time%day,month,mday)
      imonth = month
      iday   = mday

      time_next = time + int(deltim)
      call julian2monthday(time_next%year,time_next%day,month,mday)
      imonth_next = month
      iday_next   = mday

      iyear = time_next%year
      if(time_next%year .lt. 2013)iyear=2013
      if(time_next%year .gt. 2021)iyear=2021
      if(imonth_next /= imonth)then
         write(syear,"(I4.4)")  iyear
         write(smonth,"(I2.2)") month
         file_ozone = trim(DEF_dir_rawdata) // '/Ozone/China/'//trim(syear)//trim(smonth)//'_O3_v2.nc'
      end if

      IF (iday_next /= iday .and. .not.(month .eq. 2 .and. iday_next .eq. 29 .and. .not.(isleapyear(iyear)))) THEN
         CALL ncio_read_block_time (file_ozone, 'O3', grid_ozone, iday_next, f_ozone)
#ifdef CoLMDEBUG
         CALL check_block_data ('Ozone', f_ozone)
#endif         
      
         call mg2p_ozone%map_aweighted (f_ozone, forc_ozone) 
         forc_ozone = forc_ozone * 1.e-9 
#ifdef CoLMDEBUG
         call check_vector_data ('Ozone', forc_ozone)
#endif         
      ENDIF

   END SUBROUTINE update_ozone_data 

END MODULE MOD_OzoneData
