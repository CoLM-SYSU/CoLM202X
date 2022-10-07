#include <define.h>
#define OPENMP 24

module mod_vector2grid

   use precision
   USE ncio_serial
   use netcdf
   USE GlobalVars, only : spval
   implicit none

   INTEGER :: ndim1out, ndim2out
   INTEGER, allocatable :: addr2d(:,:)
   REAL(r8), allocatable :: latitude(:), longitude(:)
   CHARACTER(len=256) :: coord1name, coord2name

contains

   ! -------
   SUBROUTINE hist_vector2grid_init (meshfile, input)
      
      USE mod_utils
      USE mod_namelist
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: meshfile, input

      ! Local Variables
      INTEGER :: ndim1, ndim2, ndim0, i1, i2, i1n, i2n, iloc
      LOGICAL :: found
      INTEGER, allocatable :: basin(:)
#ifdef CATCHMENT
      INTEGER, allocatable :: cat(:,:), hru(:,:)
      INTEGER, allocatable :: htype(:)
#endif
#ifdef UNSTRUCTURED
      INTEGER, allocatable :: polygon(:,:)
#endif

#ifdef CATCHMENT
      CALL ncio_read_serial (meshfile, 'icat',  cat)
      CALL ncio_read_serial (meshfile, 'hunit', hru)
      ndim1 = size(cat,1); ndim2 = size(cat,2)

      CALL ncio_read_serial (meshfile, 'latitude',  latitude)
      CALL ncio_read_serial (meshfile, 'longitude', longitude)
      coord1name = 'latitude'
      coord2name = 'longitude'
#endif
#ifdef UNSTRUCTURED
      CALL ncio_read_serial (DEF_file_landgrid, 'patchtypes', polygon)
      ndim1 = size(polygon,1); ndim2 = size(polygon,2)
      
      longitude = (/(-180 + 360.0/43200/2 + (i1-1)*360.0/43200, i1=1,43200)/)
      latitude  = (/(  90 - 180.0/21600/2 - (i2-1)*180.0/21600, i2=1,21600)/)
      coord1name = 'longitude'
      coord2name = 'latitude'
#endif

      allocate (addr2d (ndim1, ndim2))
      addr2d(:,:) = -1

      CALL ncio_read_serial (input, 'basin', basin)
#ifdef CATCHMENT
      CALL ncio_read_serial (input, 'htype', htype)
#endif
      
      ndim0 = size(basin)

      DO i1 = 1, ndim1
         DO i2 = 1, ndim2

#ifdef CATCHMENT
            IF (cat(i1,i2) <= 0) cycle
#ENDIF
#ifdef UNSTRUCTURED
            IF (polygon(i1,i2) <= 0) cycle
#ENDIF

            found = .false.
            DO i1n = max(i1-1,1), i1
               DO i2n = max(i2-1,1), i2

#ifdef CATCHMENT
                  IF ((cat(i1n,i2n) == cat(i1,i2)) .and. (hru(i1n,i2n) == hru(i1,i2)) &
#ENDIF
#ifdef UNSTRUCTURED
                  IF ((polygon(i1n,i2n) == polygon(i1,i2)) &
#ENDIF
                     .and. (addr2d(i1n,i2n) /= -1)) THEN
                     addr2d(i1,i2) = addr2d(i1n,i2n)
                     found = .true.
                     exit
                  ENDIF
               ENDDO
            ENDDO

            IF (.not. found) THEN
#ifdef CATCHMENT
               iloc = find_in_sorted_list2 (hru(i1,i2),cat(i1,i2), ndim0, htype,basin)
#ENDIF
#ifdef UNSTRUCTURED
               iloc = find_in_sorted_list1 (polygon(i1,i2), ndim0, basin)
#ENDIF
               IF (iloc > 0) addr2d(i1,i2) = iloc
            ENDIF

         ENDDO
      ENDDO

      ndim1out = ndim1
      ndim2out = ndim2

#ifdef CATCHMENT
      deallocate (cat)
      deallocate (hru)
      deallocate (htype)
#ENDIF
#ifdef UNSTRUCTURED
      deallocate (polygon)
#ENDIF
      deallocate (basin)

      write(*,*) 'Init done.'

   END SUBROUTINE hist_vector2grid_init

   ! -------
   SUBROUTINE hist_vector2grid_one_var (input, output, varname, timelen)
      
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: input, output
      CHARACTER(len=*), intent(in) :: varname
      INTEGER, intent(in) :: timelen

      ! Local Variables
      INTEGER :: ncid, varid, ndims, id
      INTEGER :: i1, i2
      INTEGER, allocatable :: dimids(:), dimlens(:)
      CHARACTER(len=256), allocatable :: dimnames (:)
      REAL(r8), allocatable :: data2in(:,:),     data2out(:,:,:)
      REAL(r8), allocatable :: data3in(:,:,:),   data3out(:,:,:,:)
      REAL(r8), allocatable :: data4in(:,:,:,:), data4out(:,:,:,:,:)

      CALL nccheck( nf90_open (trim(input), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid (ncid, trim(varname), varid) )
      CALL nccheck( nf90_inquire_variable (ncid, varid, ndims = ndims) )
      allocate (dimids  (ndims))
      CALL nccheck( nf90_inquire_variable (ncid, varid, dimids = dimids) )
      allocate (dimnames (ndims))
      allocate (dimlens  (ndims))
      DO id = 1, ndims
         CALL nccheck( nf90_inquire_dimension (ncid, dimids(id), dimnames(id), dimlens(id)) )
      ENDDO
      CALL nccheck( nf90_close (ncid) )
               
      IF (ndims >=2) THEN
         IF ((trim(dimnames(ndims-1)) == 'hydrounit') .and. (trim(dimnames(ndims)) == 'time')) THEN
            IF (ndims == 2) THEN

               CALL ncio_read_serial (input, varname, data2in)

               allocate (data2out (ndim1out, ndim2out, timelen))
               data2out(:,:,:) = spval

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i1,i2)
#endif
               DO i1 = 1, ndim1out
                  DO i2 = 1, ndim2out
                     IF (addr2d(i1,i2) > 0) THEN
                        data2out(i1,i2,:) = data2in(addr2d(i1,i2),:)
                     ENDIF
                  ENDDO
               ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

               CALL ncio_write_serial (output, varname, data2out, &
                  coord1name, coord2name, 'time', compress = 1)

               deallocate(data2in )
               deallocate(data2out)

            ELSEIF (ndims == 3) THEN

               CALL ncio_read_serial (input, varname, data3in)

               CALL copy_dimension_nc (input, output, dimnames(1))

               allocate (data3out (dimlens(1), ndim1out, ndim2out, timelen))
               data3out(:,:,:,:) = spval
               DO i1 = 1, ndim1out
                  DO i2 = 1, ndim2out
                     IF (addr2d(i1,i2) > 0) THEN
                        data3out(:,i1,i2,:) = data3in(:,addr2d(i1,i2),:)
                     ENDIF
                  ENDDO
               ENDDO

               CALL ncio_write_serial (output, varname, data3out, &
                  dimnames(1), coord1name, coord2name, 'time', compress = 1)

               deallocate(data3in )
               deallocate(data3out)

            ELSEIF (ndims == 4) THEN

               CALL ncio_read_serial (input, varname, data4in)

               CALL copy_dimension_nc (input, output, dimnames(1))
               CALL copy_dimension_nc (input, output, dimnames(2))

               allocate (data4out (dimlens(1), dimlens(2), ndim1out, ndim2out, timelen))
               data4out(:,:,:,:,:) = spval
               DO i1 = 1, ndim1out
                  DO i2 = 1, ndim2out
                     IF (addr2d(i1,i2) > 0) THEN
                        data4out(:,:,i1,i2,:) = data4in(:,:,addr2d(i1,i2),:)
                     ENDIF
                  ENDDO
               ENDDO

               CALL ncio_write_serial (output, varname, data4out, &
                  dimnames(1), dimnames(2), coord1name, coord2name, 'time', compress = 1)

               deallocate(data4in )
               deallocate(data4out)

            ENDIF
      
            CALL copy_variable_attr_nc (input, output, varname)
            
            write(*,*) 'From <', trim(input), '> to <', trim(output), '> : ', trim(varname)
         ENDIF
      ENDIF

      deallocate (dimids)
      deallocate (dimlens)
      deallocate (dimnames)

   END SUBROUTINE hist_vector2grid_one_var

   ! ------
   subroutine copy_dimension_nc (input, output, dimname)
      
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: input, output
      CHARACTER(len=*), intent(in) :: dimname

      ! Local Variables
      INTEGER :: ncid, varid, dimid, dimlen
      LOGICAL :: dim_exist

      CALL nccheck( nf90_open (trim(output), NF90_NOWRITE, ncid) )
      dim_exist = (nf90_inq_dimid (ncid, trim(dimname), dimid) == NF90_NOERR)

      IF (.not. dim_exist) THEN
         CALL nccheck( nf90_close (ncid) )

         CALL nccheck( nf90_open (trim(input), NF90_NOWRITE, ncid) )
         CALL nccheck( nf90_inq_dimid (ncid, trim(dimname), dimid) )
         CALL nccheck( nf90_inquire_dimension (ncid, dimid, len = dimlen) )
         CALL nccheck( nf90_close (ncid) )

         CALL nccheck( nf90_open (trim(output), NF90_WRITE, ncid) )
         CALL nccheck( nf90_redef (ncid) )
         CALL nccheck( nf90_def_dim (ncid, trim(dimname), dimlen, dimid) )
         CALL nccheck( nf90_enddef (ncid) )
         CALL nccheck( nf90_close (ncid) )
      ELSE
         CALL nccheck( nf90_close (ncid) )
      ENDIF

   END subroutine copy_dimension_nc 

   ! ------
   subroutine copy_variable_nc (input, output, varname)
      
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: input, output
      CHARACTER(len=*), intent(in) :: varname

      ! Local Variables
      INTEGER :: ncid, ncidin, ncidout
      INTEGER :: varid, varidin, varidout
      INTEGER :: xtype, ndims, dimlen, id, natts, iattr
      LOGICAL :: var_exist
      INTEGER, allocatable :: dimids (:), dimlens (:)
      CHARACTER(len=256), allocatable :: dimnames (:)
      CHARACTER(len=256) :: attrname

      REAL(r8), allocatable :: data1 (:)
      REAL(r8), allocatable :: data2 (:,:)
      REAL(r8), allocatable :: data3 (:,:,:)
      REAL(r8), allocatable :: data4 (:,:,:,:)
      REAL(r8), allocatable :: data5 (:,:,:,:,:)
      REAL(r8), allocatable :: data6 (:,:,:,:,:,:)

      CALL nccheck( nf90_open (trim(output), NF90_NOWRITE, ncid) )
      var_exist = (nf90_inq_varid(ncid, trim(varname), varid) == NF90_NOERR)
      CALL nccheck( nf90_close (ncid) )

      IF (var_exist) RETURN

      CALL nccheck( nf90_open (trim(input), NF90_NOWRITE, ncid) )
      CALL nccheck( nf90_inq_varid (ncid, trim(varname), varid) )
      CALL nccheck( nf90_inquire_variable (ncid, varid, xtype = xtype, ndims = ndims) )
      allocate (dimids   (ndims))
      CALL nccheck( nf90_inquire_variable (ncid, varid, dimids = dimids) )
      allocate (dimnames (ndims))
      allocate (dimlens  (ndims))
      DO id = 1, ndims
         CALL nccheck( nf90_inquire_dimension (ncid, dimids(id), dimnames(id), dimlens(id)) )
      ENDDO
      CALL nccheck( nf90_close (ncid) )

      DO id = 1, ndims
         CALL copy_dimension_nc (input, output, dimnames(id))
      ENDDO

      CALL nccheck( nf90_open (trim(output), NF90_WRITE, ncidout) )
      DO id = 1, ndims
         CALL nccheck( nf90_inq_dimid (ncidout, dimnames(id), dimids(id)) )
      ENDDO
      CALL nccheck( nf90_redef (ncidout) )
      CALL nccheck( nf90_def_var (ncidout, trim(varname), xtype, dimids, varidout) )
      CALL nccheck( nf90_enddef (ncidout) )

      CALL nccheck( nf90_open (trim(input),  NF90_NOWRITE, ncidin)  )
      CALL nccheck( nf90_inq_varid (ncidin,  trim(varname), varidin ) )
      select case (ndims)
      case (1)
         allocate (data1 (dimlens(1)))
         CALL nccheck( nf90_get_var (ncidin,  varidin , data1) )
         CALL nccheck( nf90_put_var (ncidout, varidout, data1) )
         deallocate (data1)
      case (2)
         allocate (data2 (dimlens(1),dimlens(2)))
         CALL nccheck( nf90_get_var (ncidin,  varidin , data2) )
         CALL nccheck( nf90_put_var (ncidout, varidout, data2) )
         deallocate (data2)
      case (3)
         allocate (data3 (dimlens(1),dimlens(2),dimlens(3)))
         CALL nccheck( nf90_get_var (ncidin,  varidin , data3) )
         CALL nccheck( nf90_put_var (ncidout, varidout, data3) )
         deallocate (data3)
      case (4)
         allocate (data4 (dimlens(1),dimlens(2),dimlens(3),dimlens(4)))
         CALL nccheck( nf90_get_var (ncidin,  varidin , data4) )
         CALL nccheck( nf90_put_var (ncidout, varidout, data4) )
         deallocate (data4)
      case (5)
         allocate (data5 (dimlens(1),dimlens(2),dimlens(3),dimlens(4),dimlens(5)))
         CALL nccheck( nf90_get_var (ncidin,  varidin , data5) )
         CALL nccheck( nf90_put_var (ncidout, varidout, data5) )
         deallocate (data5)
      case (6)
         allocate (data6 (dimlens(1),dimlens(2),dimlens(3),dimlens(4),dimlens(5),dimlens(6)))
         CALL nccheck( nf90_get_var (ncidin,  varidin , data6) )
         CALL nccheck( nf90_put_var (ncidout, varidout, data6) )
         deallocate (data6)
      END select 
       
      CALL nccheck( nf90_close (ncidin)  )
      CALL nccheck( nf90_close (ncidout) )

      CALL copy_variable_attr_nc (input, output, varname)

   END subroutine copy_variable_nc 

   ! ------
   subroutine copy_variable_attr_nc (input, output, varname)
      
      IMPLICIT NONE
      CHARACTER(len=*), intent(in) :: input, output
      CHARACTER(len=*), intent(in) :: varname

      ! Local Variables
      INTEGER :: ncidin, ncidout, varidin, varidout, natts, iattr
      CHARACTER(len=256) :: attrname

      CALL nccheck( nf90_open (trim(output), NF90_WRITE, ncidout) )
      CALL nccheck( nf90_inq_varid (ncidout,  trim(varname), varidout) )
      CALL nccheck( nf90_open (trim(input),  NF90_NOWRITE, ncidin)  )
      CALL nccheck( nf90_inq_varid (ncidin,  trim(varname), varidin ) )
       
      CALL nccheck( nf90_inquire_variable (ncidin, varidin, natts = natts) )
      IF (natts > 0) THEN
         CALL nccheck( nf90_redef (ncidout) )
         DO iattr = 1, natts
            CALL nccheck( nf90_inq_attname (ncidin, varidin, iattr, attrname) )
            CALL nccheck( nf90_copy_att (ncidin, varidin, attrname, ncidout, varidout) )
         ENDDO
         CALL nccheck( nf90_enddef (ncidout) )
      ENDIF

      CALL nccheck( nf90_close (ncidin)  )
      CALL nccheck( nf90_close (ncidout) )

   END subroutine copy_variable_attr_nc 

   ! ------
   SUBROUTINE get_output_filename (input, output)
      
      IMPLICIT NONE

      CHARACTER(len=*), intent(in)  :: input
      CHARACTER(len=*), intent(out) :: output

      ! Local variables
      INTEGER :: i

      i = len_trim (input) 
      DO while (i > 0)
         IF (input(i:i) == '.') exit
         i = i - 1
      ENDDO

      IF (i > 0) THEN
         output = input(1:i-1) // '_grid' // '.nc'
      ELSE
         output = input // '_grid' // '.nc'
      ENDIF

   END SUBROUTINE get_output_filename

END module mod_vector2grid
