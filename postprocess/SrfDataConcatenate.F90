PROGRAM srfdata_concatenate

   USE MOD_NetCDFSerial
   USE MOD_Utils
   IMPLICIT NONE

   ! Local variables 
   character(len=256) :: dirlanddata, dirvar, prefix, varname    
   character(len=256) :: level, typefilter, output, rshp
   character(len=256) :: tmpfile, file_list_cmd
   integer :: timevals (8)

   logical   :: dim1to2
   integer   :: filter, nfile, ifile, nthis, dsp, ntotal
   integer*8 :: bsnmax
   character(len=256) :: line, blockinfo, varfile, levfile, landfile

   real(r8), allocatable :: longitude(:), latitude(:)
   integer  :: nlon, nlat

   type :: varreal
      real(r8), allocatable :: val(:)
   END type
   type :: varint
      integer*8, allocatable :: val(:)
   END type

   type(varreal), allocatable :: varvec(:)
   type(varint ), allocatable :: bsnvec(:)

   real(r8),  allocatable :: varcache1(:), varcache2(:,:)
   integer*8, allocatable :: bsncache1(:), bsncache2(:,:)
   real(r8),  allocatable :: vardata(:)
   integer*8, allocatable :: eindex (:)
   integer  , allocatable :: settyp (:)
   integer  , allocatable :: order  (:)

   integer :: stat, i, j, ibasin
   real(r8), parameter :: spval = -1.e36_r8  !missing value

   IF (COMMAND_ARGUMENT_COUNT() == 0) THEN
      write(*,*)  'Usage    : PATH/srfdata_contenate' 
      write(*,*)  ' <Arg 1> : path_to_landdata'
      write(*,*)  ' <Arg 2> : dir_var'
      write(*,*)  ' <Arg 3> : prefix'
      write(*,*)  ' <Arg 4> : varname'
      write(*,*)  ' <Arg 5> : patch | pft | urban'
      write(*,*)  ' <Arg 6> : typefilter'
      write(*,*)  ' <Arg 7> : output'
      write(*,*)  ' <Arg 8> : reshape or not (optional, value T/F, reshape from 1D to 2D)'

      STOP
   ENDIF

   CALL get_command_argument (1, dirlanddata)
   CALL get_command_argument (2, dirvar     )
   CALL get_command_argument (3, prefix     )
   CALL get_command_argument (4, varname    )
   CALL get_command_argument (5, level      )
   CALL get_command_argument (6, typefilter )
   CALL get_command_argument (7, output     )

   dim1to2 = .false.
   IF (COMMAND_ARGUMENT_COUNT() > 7) THEN
      CALL get_command_argument (8, rshp)
      IF (trim(rshp) == 'T') THEN
         dim1to2 = .true.
      ENDIF
   ENDIF

   read(typefilter,*) filter

   CALL date_and_time (values = timevals)
   write(tmpfile,'("temp_", I4.4,5I2.2,I3.3, ".txt")') timevals(1:3), timevals(5:8)

   file_list_cmd = 'ls ' // trim(dirlanddata) // '/' // trim(dirvar) &
      // '/' // trim(prefix) // '*.nc > ' // trim(tmpfile)
   CALL system(file_list_cmd)
   
   nfile = 0
   open(unit=10, file=trim(tmpfile))
   DO WHILE (.true.)
      read(10, '(A)', IOSTAT=stat) line
      IF(IS_IOSTAT_END(stat)) EXIT

      nfile = nfile + 1
   ENDDO
   close(10)

   allocate (varvec (nfile))
   allocate (bsnvec (nfile))

   open(unit=10, file=trim(tmpfile))
   DO ifile = 1, nfile

      read(10, '(A)', IOSTAT=stat) line
      IF(IS_IOSTAT_END(stat)) EXIT

      blockinfo = line(len_trim(line)-10:len_trim(line)-3)

      varfile = trim(dirlanddata) // '/' // trim(dirvar) &
         // '/' // trim(prefix) // '_' // trim(blockinfo) // '.nc'

      CALL ncio_read_serial (varfile, varname, vardata)

      IF (trim(level) == 'patch') THEN
         levfile = trim(dirlanddata) // '/landpatch/landpatch_' //  trim(blockinfo) // '.nc'
      ELSEIF (trim(level) == 'pft') THEN
         levfile = trim(dirlanddata) // '/landpft/landpft_' //  trim(blockinfo) // '.nc'
      ELSEIF (trim(level) == 'urban') THEN
         levfile = trim(dirlanddata) // '/landurban/landurban_' //  trim(blockinfo) // '.nc'
      ENDIF
         
      CALL ncio_read_serial (levfile, 'eindex', eindex)
      CALL ncio_read_serial (levfile, 'settyp', settyp)

      nthis = count(settyp == filter)
      IF (nthis > 0) THEN
         allocate (varvec(ifile)%val(nthis))
         allocate (bsnvec(ifile)%val(nthis))

         varvec(ifile)%val = pack(vardata, settyp == filter)
         bsnvec(ifile)%val = pack(eindex,  settyp == filter)
      ENDIF

   ENDDO
   close(10)

   ntotal = 0
   bsnmax = -1
   DO ifile = 1, nfile
      IF (allocated(varvec(ifile)%val)) THEN
         ntotal = ntotal + size(varvec(ifile)%val)
      
         IF (dim1to2) THEN
            bsnmax = max(bsnmax, maxval(bsnvec(ifile)%val))
         ENDIF
      ENDIF
   ENDDO

   IF (ntotal > 0) THEN

      IF (dim1to2) THEN
         landfile = trim(dirlanddata) // '/mesh/mesh.nc'
         CALL ncio_read_serial (landfile, 'latitude' , latitude )
         CALL ncio_read_serial (landfile, 'longitude', longitude)

         nlat = size(latitude)
         nlon = size(longitude)

         allocate (varcache2 (nlon, nlat))
         allocate (bsncache2 (nlon, nlat))

         varcache2(:,:) = spval
         bsncache2(:,:) = -1
         
         DO ifile = 1, nfile
            IF (allocated(varvec(ifile)%val)) THEN
               nthis = size(varvec(ifile)%val)

               DO ibasin = 1, nthis
                  i = mod(bsnvec(ifile)%val(ibasin), nlon)
                  IF (i == 0) i = nlon 
                  j = (bsnvec(ifile)%val(ibasin)-1) / nlon + 1

                  varcache2(i,j) = varvec(ifile)%val(ibasin)
                  bsncache2(i,j) = bsnvec(ifile)%val(ibasin)

               ENDDO
            ENDIF
         ENDDO
      ELSE
         allocate (varcache1 (ntotal))
         allocate (bsncache1 (ntotal))
         
         dsp = 0
         DO ifile = 1, nfile
            IF (allocated(varvec(ifile)%val)) THEN
               nthis = size(varvec(ifile)%val)
                  
               varcache1(dsp+1:dsp+nthis) = varvec(ifile)%val
               bsncache1(dsp+1:dsp+nthis) = bsnvec(ifile)%val

               dsp = dsp + nthis
            ENDIF
         ENDDO
         
         allocate (order (ntotal))
         order = (/(i, i=1,ntotal)/)

         CALL quicksort (ntotal, bsncache1, order)

         varcache1 = varcache1(order)
      ENDIF

      CALL ncio_create_file (output)
   
      IF (dim1to2) THEN
         CALL ncio_define_dimension (output, 'longitude', nlon)
         CALL ncio_define_dimension (output, 'latitude' , nlat)
         CALL ncio_write_serial (output, varname,    varcache2, 'longitude', 'latitude')
         CALL ncio_write_serial (output, 'elmindex', bsncache2, 'longitude', 'latitude')
               
         CALL ncio_write_serial (output, 'latitude', latitude, 'latitude')
         CALL ncio_put_attr (output, 'latitude' , 'long_name', 'latitude')
         CALL ncio_put_attr (output, 'latitude' , 'units', 'degrees_north')
         
         CALL ncio_write_serial (output, 'longitude', longitude, 'longitude')
         CALL ncio_put_attr (output, 'longitude', 'long_name', 'longitude')
         CALL ncio_put_attr (output, 'longitude', 'units', 'degrees_east')
      ELSE
         CALL ncio_define_dimension (output, 'vec', ntotal)
         CALL ncio_write_serial (output, 'elmindex', bsncache1, 'vec')
         CALL ncio_write_serial (output, varname,    varcache1, 'vec')
      ENDIF
         
      CALL ncio_put_attr (output, varname, 'missing_value', spval)
   ENDIF

   IF (allocated(varvec)) THEN
      DO ifile = 1, size(varvec)
         IF (allocated(varvec(ifile)%val)) deallocate(varvec(ifile)%val)
      ENDDO
      deallocate(varvec)
   ENDIF
   IF (allocated(bsnvec)) THEN
      DO ifile = 1, size(bsnvec)
         IF (allocated(bsnvec(ifile)%val)) deallocate(bsnvec(ifile)%val)
      ENDDO
      deallocate(bsnvec)
   ENDIF
   
   IF (allocated(longitude)) deallocate(longitude)
   IF (allocated(latitude )) deallocate(latitude )

   IF (allocated(varcache1)) deallocate(varcache1)
   IF (allocated(varcache2)) deallocate(varcache2)
   IF (allocated(bsncache1)) deallocate(bsncache1)
   IF (allocated(bsncache2)) deallocate(bsncache2)

   IF (allocated(vardata)) deallocate(vardata)
   IF (allocated(eindex )) deallocate(eindex )
   IF (allocated(settyp )) deallocate(settyp )
   IF (allocated(order  )) deallocate(order  )

   file_list_cmd = 'rm ' // tmpfile
   CALL system(file_list_cmd)

END PROGRAM srfdata_concatenate
