#include <define.h>
#define OPENMP 24

MODULE mod_vector2grid

   USE MOD_Precision
   USE MOD_NetCDFSerial
   USE netcdf
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   integer :: ndim1out, ndim2out
   integer, allocatable :: addr2d(:,:)
   real(r8), allocatable :: latitude(:), longitude(:)
   character(len=256) :: coord1name, coord2name

CONTAINS

   ! -------
   SUBROUTINE hist_vector2grid_init (meshfile, input)
      
   USE MOD_Utils
   IMPLICIT NONE
   character(len=*), intent(in) :: meshfile, input

   ! Local Variables
   integer :: ndim1, ndim2, ndim0, i1, i2, i1n, i2n, iloc
   integer :: i1s, i1e, i2s, i2e
   logical :: found
   integer, allocatable :: elm(:)
#ifdef CATCHMENT
   integer, allocatable :: cat(:,:), hru(:,:)
   integer, allocatable :: htype(:)
#endif
#ifdef UNSTRUCTURED
   integer, allocatable :: polygon(:,:)
#endif
   integer, allocatable  :: addr2d_g(:,:)
   real(r8), allocatable :: lon(:), lat(:)

#ifdef CATCHMENT
      CALL ncio_read_serial (meshfile, 'icatchment2d', cat)
      CALL ncio_read_serial (meshfile, 'ihydrounit2d', hru)
      ndim1 = size(cat,1); ndim2 = size(cat,2)
#endif
#ifdef UNSTRUCTURED
      CALL ncio_read_serial (meshfile, 'elmindex', polygon)
      ndim1 = size(polygon,1); ndim2 = size(polygon,2)
#endif

      allocate (addr2d_g (ndim1, ndim2))
      addr2d_g(:,:) = -1

#ifdef CATCHMENT
      CALL ncio_read_serial (input, 'bsn_hru', elm  )
      CALL ncio_read_serial (input, 'typ_hru', htype)
#else
      CALL ncio_read_serial (input, 'elmindex', elm)
#endif
      
      ndim0 = size(elm)

      DO i1 = 1, ndim1
         DO i2 = 1, ndim2

#ifdef CATCHMENT
            IF (cat(i1,i2) <= 0) CYCLE
#endif
#ifdef UNSTRUCTURED
            IF (polygon(i1,i2) <= 0) CYCLE
#endif

            found = .false.
            DO i1n = max(i1-1,1), i1
               DO i2n = max(i2-1,1), i2

#ifdef CATCHMENT
                  IF ((cat(i1n,i2n) == cat(i1,i2)) .and. (hru(i1n,i2n) == hru(i1,i2)) &
#endif
#ifdef UNSTRUCTURED
                  IF ((polygon(i1n,i2n) == polygon(i1,i2)) &
#endif
                     .and. (addr2d_g(i1n,i2n) /= -1)) THEN
                     addr2d_g(i1,i2) = addr2d_g(i1n,i2n)
                     found = .true.
                     EXIT
                  ENDIF
               ENDDO
            ENDDO

            IF (.not. found) THEN
#ifdef CATCHMENT
               iloc = find_in_sorted_list2 (hru(i1,i2),cat(i1,i2), ndim0, htype,elm)
#endif
#ifdef UNSTRUCTURED
               iloc = find_in_sorted_list1 (polygon(i1,i2), ndim0, elm)
#endif
               IF (iloc > 0) addr2d_g(i1,i2) = iloc
            ENDIF

         ENDDO
      ENDDO

      DO i1 = 1, ndim1
         IF (any(addr2d_g(i1,:) > 0)) THEN
            i1s = i1
            EXIT
         ENDIF
      ENDDO
      
      DO i1 = ndim1, 1, -1
         IF (any(addr2d_g(i1,:) > 0)) THEN
            i1e = i1
            EXIT
         ENDIF
      ENDDO

      DO i2 = 1, ndim2
         IF (any(addr2d_g(:,i2) > 0)) THEN
            i2s = i2
            EXIT
         ENDIF
      ENDDO
      
      DO i2 = ndim2, 1, -1
         IF (any(addr2d_g(:,i2) > 0)) THEN
            i2e = i2
            EXIT
         ENDIF
      ENDDO

      ndim1out = i1e - i1s + 1
      ndim2out = i2e - i2s + 1

      allocate(addr2d (ndim1out,ndim2out))
      addr2d = addr2d_g(i1s:i1e,i2s:i2e)

      CALL ncio_read_serial (meshfile, 'latitude',  lat)
      CALL ncio_read_serial (meshfile, 'longitude', lon)

#ifdef UNSTRUCTURED
      allocate(longitude(ndim1out))
      allocate(latitude (ndim2out))

      longitude = lon(i1s:i1e)
      latitude  = lat(i2s:i2e)

      coord1name = 'longitude'
      coord2name = 'latitude'
#endif
#ifdef CATCHMENT
      allocate(latitude (ndim1out))
      allocate(longitude(ndim2out))

      latitude  = lat(i1s:i1e)
      longitude = lon(i2s:i2e)

      coord1name = 'latitude'
      coord2name = 'longitude'
#endif

      deallocate (addr2d_g)
#ifdef CATCHMENT
      deallocate (cat)
      deallocate (hru)
      deallocate (htype)
#endif
#ifdef UNSTRUCTURED
      deallocate (polygon)
#endif
      deallocate (elm)

      deallocate(lat, lon)

      write(*,*) 'Init done.'

   END SUBROUTINE hist_vector2grid_init

   ! -------
   SUBROUTINE hist_vector2grid_one_var_time (input, output, varname, timelen)
      
   IMPLICIT NONE
   character(len=*), intent(in) :: input, output
   character(len=*), intent(in) :: varname
   integer, intent(in) :: timelen

   ! Local Variables
   integer :: ncid, varid, ndims, id
   integer :: i1, i2
   integer, allocatable :: dimids(:), dimlens(:)
   character(len=256), allocatable :: dimnames (:)
   real(r8), allocatable :: data2in(:,:),     data2out(:,:,:)
   real(r8), allocatable :: data3in(:,:,:),   data3out(:,:,:,:)
   real(r8), allocatable :: data4in(:,:,:,:), data4out(:,:,:,:,:)

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
         IF (((trim(dimnames(ndims-1)) == 'hydrounit') .or. (trim(dimnames(ndims-1)) == 'element')) &
            .and. (trim(dimnames(ndims)) == 'time')) THEN

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
            
            write(*,*) '     ', trim(varname), ' done.'
         ENDIF
      ENDIF

      deallocate (dimids)
      deallocate (dimlens)
      deallocate (dimnames)

   END SUBROUTINE hist_vector2grid_one_var_time

   ! -------
   SUBROUTINE hist_vector2grid_one_var (input, output, varname)
      
   IMPLICIT NONE
   character(len=*), intent(in) :: input, output
   character(len=*), intent(in) :: varname

   ! Local Variables
   integer :: ncid, varid, ndims, id
   integer :: i1, i2
   integer, allocatable :: dimids(:), dimlens(:)
   character(len=256), allocatable :: dimnames (:)
   real(r8), allocatable :: data1in(:),   data1out(:,:)
   real(r8), allocatable :: data2in(:,:), data2out(:,:,:)

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
               
      IF (ndims == 1) THEN

         CALL ncio_read_serial (input, varname, data1in)

         allocate (data1out (ndim1out, ndim2out))
         data1out(:,:) = spval

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i1,i2)
#endif
         DO i1 = 1, ndim1out
            DO i2 = 1, ndim2out
               IF (addr2d(i1,i2) > 0) THEN
                  data1out(i1,i2) = data1in(addr2d(i1,i2))
               ENDIF
            ENDDO
         ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

         CALL ncio_write_serial (output, varname, data1out, &
            coord1name, coord2name, compress = 1)

         deallocate(data1in )
         deallocate(data1out)

      ELSEIF (ndims == 2) THEN

         CALL ncio_read_serial (input, varname, data2in)

         CALL copy_dimension_nc (input, output, dimnames(1))

         allocate (data2out (dimlens(1), ndim1out, ndim2out))
         data2out(:,:,:) = spval
         DO i1 = 1, ndim1out
            DO i2 = 1, ndim2out
               IF (addr2d(i1,i2) > 0) THEN
                  data2out(:,i1,i2) = data2in(:,addr2d(i1,i2))
               ENDIF
            ENDDO
         ENDDO

         CALL ncio_write_serial (output, varname, data2out, &
            dimnames(1), coord1name, coord2name, compress = 1)

         deallocate(data2in )
         deallocate(data2out)

      ENDIF

      deallocate (dimids)
      deallocate (dimlens)
      deallocate (dimnames)

   END SUBROUTINE hist_vector2grid_one_var

   ! ------
   SUBROUTINE copy_dimension_nc (input, output, dimname)
      
   IMPLICIT NONE
   character(len=*), intent(in) :: input, output
   character(len=*), intent(in) :: dimname

   ! Local Variables
   integer :: ncid, varid, dimid, dimlen
   logical :: dim_exist

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

   END SUBROUTINE copy_dimension_nc 

   ! ------
   SUBROUTINE copy_variable_nc (input, output, varname)
      
   IMPLICIT NONE
   character(len=*), intent(in) :: input, output
   character(len=*), intent(in) :: varname

   ! Local Variables
   integer :: ncid, ncidin, ncidout
   integer :: varid, varidin, varidout
   integer :: xtype, ndims, dimlen, id, natts, iattr
   logical :: var_exist
   integer, allocatable :: dimids (:), dimlens (:)
   character(len=256), allocatable :: dimnames (:)
   character(len=256) :: attrname

   real(r8), allocatable :: data1 (:)
   real(r8), allocatable :: data2 (:,:)
   real(r8), allocatable :: data3 (:,:,:)
   real(r8), allocatable :: data4 (:,:,:,:)
   real(r8), allocatable :: data5 (:,:,:,:,:)
   real(r8), allocatable :: data6 (:,:,:,:,:,:)

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
      select CASE (ndims)
      CASE (1)
         allocate (data1 (dimlens(1)))
         CALL nccheck( nf90_get_var (ncidin,  varidin , data1) )
         CALL nccheck( nf90_put_var (ncidout, varidout, data1) )
         deallocate (data1)
      CASE (2)
         allocate (data2 (dimlens(1),dimlens(2)))
         CALL nccheck( nf90_get_var (ncidin,  varidin , data2) )
         CALL nccheck( nf90_put_var (ncidout, varidout, data2) )
         deallocate (data2)
      CASE (3)
         allocate (data3 (dimlens(1),dimlens(2),dimlens(3)))
         CALL nccheck( nf90_get_var (ncidin,  varidin , data3) )
         CALL nccheck( nf90_put_var (ncidout, varidout, data3) )
         deallocate (data3)
      CASE (4)
         allocate (data4 (dimlens(1),dimlens(2),dimlens(3),dimlens(4)))
         CALL nccheck( nf90_get_var (ncidin,  varidin , data4) )
         CALL nccheck( nf90_put_var (ncidout, varidout, data4) )
         deallocate (data4)
      CASE (5)
         allocate (data5 (dimlens(1),dimlens(2),dimlens(3),dimlens(4),dimlens(5)))
         CALL nccheck( nf90_get_var (ncidin,  varidin , data5) )
         CALL nccheck( nf90_put_var (ncidout, varidout, data5) )
         deallocate (data5)
      CASE (6)
         allocate (data6 (dimlens(1),dimlens(2),dimlens(3),dimlens(4),dimlens(5),dimlens(6)))
         CALL nccheck( nf90_get_var (ncidin,  varidin , data6) )
         CALL nccheck( nf90_put_var (ncidout, varidout, data6) )
         deallocate (data6)
      END select 
       
      CALL nccheck( nf90_close (ncidin)  )
      CALL nccheck( nf90_close (ncidout) )

      CALL copy_variable_attr_nc (input, output, varname)

   END SUBROUTINE copy_variable_nc 

   ! ------
   SUBROUTINE copy_variable_attr_nc (input, output, varname)
      
   IMPLICIT NONE
   character(len=*), intent(in) :: input, output
   character(len=*), intent(in) :: varname

   ! Local Variables
   integer :: ncidin, ncidout, varidin, varidout, natts, iattr
   character(len=256) :: attrname

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

   END SUBROUTINE copy_variable_attr_nc 

   ! ------
   SUBROUTINE get_output_filename (input, output)
      
   IMPLICIT NONE

      character(len=*), intent(in)  :: input
      character(len=*), intent(out) :: output

      ! Local variables
      integer :: i

      i = len_trim (input) 
      DO WHILE (i > 0)
         IF (input(i:i) == '.') EXIT
         i = i - 1
      ENDDO

      IF (i > 0) THEN
         output = input(1:i-1) // '_grid' // '.nc'
      ELSE
         output = input // '_grid' // '.nc'
      ENDIF

   END SUBROUTINE get_output_filename

END MODULE mod_vector2grid
