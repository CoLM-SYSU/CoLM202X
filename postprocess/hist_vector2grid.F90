#include <define.h>

program hist_vector2grid

   use mod_namelist
   use mod_vector2grid
   USE netcdf
   USE ncio_serial

   implicit none

   ! Local variables 
   LOGICAL :: fexists
   character(len=256) :: meshfile, input, output
   INTEGER :: ifile, ncid, nvars, ivar, timelen
   INTEGER, allocatable :: varids(:)
   CHARACTER(len=256), allocatable :: varnames(:)

   call getarg (1, meshfile)

   ifile = 0
   DO WHILE (.true.)
      read(*,*) input
      inquire (file=trim(input), exist=fexists)
      IF (fexists) THEN

         IF (ifile == 0) THEN
            CALL hist_vector2grid_init (meshfile, input)
         ENDIF
         
         CALL get_output_filename (input, output)

         CALL ncio_create_file (output)

         CALL copy_dimension_nc (input, output, 'time')
         CALL copy_variable_nc  (input, output, 'time')
         CALL ncio_define_dimension (output, 'latitude',  size(latitude))
         CALL ncio_define_dimension (output, 'longitude', size(longitude))
         CALL ncio_write_serial (output, 'latitude',  latitude,  'latitude')
         CALL ncio_write_serial (output, 'longitude', longitude, 'longitude')

         CALL ncio_inquire_length (input, 'time', timelen)

         CALL nccheck( nf90_open (trim(input), NF90_NOWRITE, ncid) )
         CALL nccheck( nf90_inquire (ncid, nVariables = nvars) ) 

         allocate (varids   (nvars))
         allocate (varnames (nvars))
         
         CALL nccheck( nf90_inq_varids (ncid, nvars, varids) ) 
         DO ivar = 1, nvars
            CALL nccheck( nf90_inquire_variable (ncid, varids(ivar), varnames(ivar)) )
         ENDDO
         CALL nccheck( nf90_close (ncid) )
      
         DO ivar = 1, nvars
            CALL hist_vector2grid_one_var (input, output, varnames(ivar), timelen)
         ENDDO

         deallocate (varids)
         deallocate (varnames)

         ifile = ifile + 1
      ELSE
         exit
      ENDIF
   ENDDO

   deallocate (addr2d)

end program hist_vector2grid
