#include <define.h>

PROGRAM post_vector2grid

   USE MOD_Namelist
   USE MOD_Vector2Grid
   USE netcdf
   USE MOD_NetCDFSerial

   IMPLICIT NONE

   ! Local variables 
   logical :: fexists
   character(len=256) :: meshfile, inputlist, input, output
   integer :: ifile, ncid, nvars, ivar, timelen
   integer, allocatable :: varids(:)
   character(len=256), allocatable :: varnames(:)
   character(len=256) :: svarname
   logical :: onevar

   integer :: iost
   
   IF (COMMAND_ARGUMENT_COUNT() == 0) THEN
      write(*,*)  'Usage    : PATH/post_vector2grid.x' 
      write(*,*)  ' <Arg 1> : meshfile'
      write(*,*)  ' <Arg 2> : vector_file_list'
      
      STOP
   ENDIF

   CALL get_command_argument (1, meshfile )
   CALL get_command_argument (2, inputlist)

   onevar = .false.
   IF (COMMAND_ARGUMENT_COUNT() > 2) THEN
      CALL get_command_argument (3, svarname)
      onevar = .true.
   ENDIF

   open(unit=18, file=trim(inputlist))
   ifile = 0
   DO WHILE (.true.)
      read(18,*,iostat=iost) input
      IF (iost /= 0) EXIT
   
      inquire (file=trim(input), exist=fexists)
      IF (fexists) THEN

         IF (ifile == 0) THEN
            CALL hist_vector2grid_init (meshfile, input)
         ENDIF
         
         CALL get_output_filename (input, output)
      
         write(*,'(A," --> ",A)') trim(input), trim(output)

         CALL ncio_create_file (output)

         CALL copy_dimension_nc (input, output, 'time')
         CALL copy_variable_nc  (input, output, 'time')
         CALL ncio_define_dimension (output, 'latitude',  size(latitude))
         CALL ncio_define_dimension (output, 'longitude', size(longitude))
         CALL ncio_write_serial (output, 'latitude',  latitude,  'latitude')
         CALL ncio_write_serial (output, 'longitude', longitude, 'longitude')

         IF (.not. onevar) THEN

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
               CALL hist_vector2grid_one_var_time (input, output, varnames(ivar), timelen)
            ENDDO

            deallocate (varids)
            deallocate (varnames)

         ELSE

            CALL hist_vector2grid_one_var (input, output, svarname)

         ENDIF

         ifile = ifile + 1
      ENDIF
   ENDDO
   close(18)

   IF (allocated(addr2d)) deallocate (addr2d)

END PROGRAM post_vector2grid
