#include <define.h>

#ifdef SinglePoint
module MOD_HistSingle

   !----------------------------------------------------------------------------
   ! DESCRIPTION:
   ! 
   !     Write out model results at sites to history files.
   !
   ! Created by Shupeng Zhang, July 2023
   !
   ! TODO...(need complement)
   !----------------------------------------------------------------------------

contains
! ----- subroutines ------

   ! -- write history time --
   subroutine hist_single_write_time (filename, dataname, time, itime)

      use MOD_Namelist
      USE MOD_SingleSrfData
      USE MOD_NetCDFSerial
      USE MOD_Landpatch, only : numpatch
#ifdef URBAN_MODEL
      USE MOD_Landurban, only : numurban
#endif
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname

      integer, intent(in)  :: time(3)
      integer, intent(out) :: itime

      ! Local variables
      logical :: fexists

      inquire (file=filename, exist=fexists)
      if (.not. fexists) then
         call ncio_create_file (trim(filename))
         CALL ncio_define_dimension(filename, 'time', 0)
         call ncio_define_dimension(filename, 'patch', numpatch)
#ifdef URBAN_MODEL
         call ncio_define_dimension(filename, 'urban', numurban)
#endif

         call ncio_write_serial (filename, 'lat', SITE_lat_location)
         CALL ncio_put_attr (filename, 'lat', 'long_name', 'latitude')
         CALL ncio_put_attr (filename, 'lat', 'units', 'degrees_north')

         call ncio_write_serial (filename, 'lon', SITE_lon_location)
         CALL ncio_put_attr (filename, 'lon', 'long_name', 'longitude')
         CALL ncio_put_attr (filename, 'lon', 'units', 'degrees_east')
      endif

      call ncio_write_time (filename, dataname, time, itime, DEF_HIST_FREQ)

   END SUBROUTINE hist_single_write_time 

end module MOD_HistSingle
#endif
