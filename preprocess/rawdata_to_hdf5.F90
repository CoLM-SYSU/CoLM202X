PROGRAM bin_to_hdf5

   USE MOD_Precision
   USE colm_io_serial
   IMPLICIT NONE

   integer, parameter :: nlat = 21600
   integer, parameter :: nlon = 43200

   integer, parameter :: nxblk = 5
   integer, parameter :: nyblk = 5

   character (len=256) :: bindir
   character (len=256) :: h5dir
   character (len=256) :: lndname

   character, allocatable :: a_chr1 (:,:)
   integer*1, allocatable :: a_int8 (:,:)
   integer*2, allocatable :: a_int16 (:,:)
   real(r8),  allocatable :: a_real8 (:,:)

   integer :: length
   integer :: irow
   integer :: iunit

   integer :: n8, Julian_day
   character(len=256) :: c

   integer, parameter :: compress = 1

   CALL getarg (1, bindir)
   CALL getarg (2, h5dir)

   !-------------------------------
   lndname = trim(bindir) // 'forest_height/Forest_Height.bin'

   allocate (a_chr1 (nlon,nlat))
   iunit = 100
   inquire (iolength=length) a_chr1 (:,1)
   open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
   DO irow = 1, nlat
      read (iunit, rec=irow) a_chr1 (:,irow)
   ENDDO
   close (iunit)

   allocate (a_int8 (nlon, nlat))
   a_int8 = ichar(a_chr1)

   lndname = trim(h5dir) // 'Forest_Height.h5'
   CALL colm_create_file (lndname)
   CALL colm_write_serial (lndname, '/forest_height', a_int8, &
      compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

   write(*,*) 'Forest height done'

   deallocate (a_chr1)
   deallocate (a_int8)

   !-------------------------------
   lndname = trim(bindir) // 'glacier/glacier.bin'

   allocate (a_int16 (nlon,nlat))
   iunit = 100
   inquire (iolength=length) a_int16 (:,1)
   open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
   DO irow = 1, nlat
      read (iunit, rec=irow) a_int16 (:,irow)
   ENDDO
   close (iunit)

   lndname = trim(h5dir) // 'glacier.h5'
   CALL colm_create_file (lndname)
   CALL colm_write_serial (lndname, '/glacier', a_int16, &
      compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

   write(*,*) 'Glacier done'

   deallocate (a_int16)

   !-------------------------------

   allocate (a_chr1 (nlon,nlat))
   allocate (a_int8 (nlon, nlat))

   CALL execute_command_line ('mkdir -p ' // trim(h5dir) // '/lai/global_30s_10_year_avg')

   DO n8 = 1, 46
      Julian_day = 1 + (N8-1)*8
      write(c,'(i3.3)') Julian_day

      lndname = trim(bindir) // 'lai/global_30s_10_year_avg/LAI_BNU_' // trim(c)

      iunit = 100
      inquire (iolength=length) a_chr1 (:,1)
      open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
      DO irow = 1, nlat
         read (iunit, rec=irow) a_chr1 (:,irow)
      ENDDO
      close (iunit)

      a_int8 = ichar(a_chr1)

      lndname = trim(h5dir) // 'lai/global_30s_10_year_avg/LAI_BNU_' // trim(c) // '.h5'
      CALL colm_create_file (lndname)
      CALL colm_write_serial (lndname, '/lai', a_int8, &
         compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

      write(*,*) 'lai ' // trim(c) // ' done'
   ENDDO

   deallocate (a_chr1)
   deallocate (a_int8)

   !-------------------------------
   lndname = trim(bindir) // 'lake_depth/GlobalLakeDepth.bin'

   allocate (a_int16 (nlon,nlat))
   iunit = 100
      inquire (iolength=length) a_int16 (:,1)
   open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
   DO irow = 1, nlat
      read (iunit, rec=irow) a_int16 (:,irow)
   ENDDO
   close (iunit)

   lndname = trim(h5dir) // 'lake_depth.h5'
   CALL colm_create_file (lndname)
   CALL colm_write_serial (lndname, '/lake_depth', a_int16, &
      compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

   write(*,*) 'Lake depth done'

   deallocate (a_int16)

   !-------------------------------
   lndname = trim(bindir) // 'lake_wetland/glwd.bin'

   allocate (a_chr1 (nlon,nlat))
   iunit = 100
      inquire (iolength=length) a_chr1 (:,1)
   open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
   DO irow = 1, nlat
      read (iunit, rec=irow) a_chr1 (:,irow)
   ENDDO
   close (iunit)

   allocate (a_int8 (nlon, nlat))
   a_int8 = ichar(a_chr1)

   lndname = trim(h5dir) // 'lake_wetland.h5'
   CALL colm_create_file (lndname)
   CALL colm_write_serial (lndname, '/lake_wetland', a_int8, &
      compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

   write(*,*) 'Lake wetland done'

   deallocate (a_chr1)
   deallocate (a_int8)

   !-------------------------------
   lndname = trim(bindir) // 'RAW_DATA_updated/landtypes_usgs_update.bin'

   allocate (a_chr1 (nlon,nlat))
   iunit = 100
      inquire (iolength=length) a_chr1 (:,1)
   open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
   DO irow = 1, nlat
      read (iunit, rec=irow) a_chr1 (:,irow)
   ENDDO
   close (iunit)

   allocate (a_int8 (nlon, nlat))
   a_int8 = ichar(a_chr1)

   lndname = trim(h5dir) // 'landtype_usgs_update.h5'
   CALL colm_create_file (lndname)
   CALL colm_write_serial (lndname, '/landtype', a_int8, &
      compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

   write(*,*) 'Landtype done'

   deallocate (a_chr1)
   deallocate (a_int8)

   !-------------------------------
   lndname = trim(bindir) // 'soil_brightness/soilcol_clm_30s.bin'

   allocate (a_chr1 (nlon,nlat))
   iunit = 100
      inquire (iolength=length) a_chr1 (:,1)
   open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
   DO irow = 1, nlat
      read (iunit, rec=irow) a_chr1 (:,irow)
   ENDDO
   close (iunit)

   allocate (a_int8 (nlon, nlat))
   a_int8 = ichar(a_chr1)

   lndname = trim(h5dir) // 'soil_brightness.h5'
   CALL colm_create_file (lndname)
   CALL colm_write_serial (lndname, '/soil_brightness', a_int8, &
      compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

   write(*,*) 'Soil brightness done'

   deallocate (a_chr1)
   deallocate (a_int8)

   !-------------------------------
   allocate (a_real8 (nlon,nlat))

   CALL execute_command_line ('mkdir -p ' // trim(h5dir) // '/soil')

   lndname = trim(h5dir) // 'soil/theta_s.h5'
   CALL colm_create_file (lndname)

   lndname = trim(h5dir) // 'soil/psi_s.h5'
   CALL colm_create_file (lndname)

   lndname = trim(h5dir) // 'soil/lambda.h5'
   CALL colm_create_file (lndname)

   lndname = trim(h5dir) // 'soil/k_s.h5'
   CALL colm_create_file (lndname)

   lndname = trim(h5dir) // 'soil/csol.h5'
   CALL colm_create_file (lndname)

   lndname = trim(h5dir) // 'soil/tksatu.h5'
   CALL colm_create_file (lndname)

   lndname = trim(h5dir) // 'soil/tkdry.h5'
   CALL colm_create_file (lndname)

   DO n8 = 1, 8
      write(c,'(i1)') n8

      ! (1) Read in the saturated water content [cm3/cm3]
      lndname = trim(bindir) // 'RAW_DATA_updated/theta_s_l'//trim(c)

      iunit = 100
      inquire (iolength=length) a_real8 (:,1)
      open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
      DO irow = 1, nlat
         read (iunit, rec=irow) a_real8 (:,irow)
      ENDDO
      close (iunit)

      lndname = trim(h5dir) // 'soil/theta_s.h5'
      CALL colm_write_serial (lndname, '/theta_s_l'//trim(c), a_real8, &
         compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

      write(*,*) 'Theta_s_l' // trim(c) // ' done'

      ! (2) Read in the matric potential at saturation [cm]
      lndname = trim(bindir) // 'RAW_DATA_updated/psi_s_l'//trim(c)

      iunit = 100
      inquire (iolength=length) a_real8 (:,1)
      open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
      DO irow = 1, nlat
         read (iunit, rec=irow) a_real8 (:,irow)
      ENDDO
      close (iunit)

      lndname = trim(h5dir) // 'soil/psi_s.h5'
      CALL colm_write_serial (lndname, 'psi_s_l'//trim(c), a_real8, &
         compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

      write(*,*) 'psi_s_l' // trim(c) // ' done'

      ! (3) Read in the pore size distribution index [dimensionless]
      lndname = trim(bindir) // 'RAW_DATA_updated/lambda_l'//trim(c)

      iunit = 100
      inquire (iolength=length) a_real8 (:,1)
      open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
      DO irow = 1, nlat
         read (iunit, rec=irow) a_real8 (:,irow)
      ENDDO
      close (iunit)

      lndname = trim(h5dir) // 'soil/lambda.h5'
      CALL colm_write_serial (lndname, 'lambda_l'//trim(c), a_real8, &
         compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

      write(*,*) 'lambda_l' // trim(c) // ' done'

      ! (4) Read in the saturated hydraulic conductivity [cm/day]
      lndname = trim(bindir) // 'RAW_DATA_updated/k_s_l'//trim(c)

      iunit = 100
      inquire (iolength=length) a_real8 (:,1)
      open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
      DO irow = 1, nlat
         read (iunit, rec=irow) a_real8 (:,irow)
      ENDDO
      close (iunit)

      lndname = trim(h5dir) // 'soil/k_s.h5'
      CALL colm_write_serial (lndname, 'k_s_l'//trim(c), a_real8, &
         compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

      write(*,*) 'k_s_l' // trim(c) // ' done'

      ! (5) Read in the heat capacity of soil solids [J/(m3 K)]
      lndname = trim(bindir) // 'RAW_DATA_updated/csol_l'//trim(c)

      iunit = 100
      inquire (iolength=length) a_real8 (:,1)
      open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
      DO irow = 1, nlat
         read (iunit, rec=irow) a_real8 (:,irow)
      ENDDO
      close (iunit)

      lndname = trim(h5dir) // 'soil/csol.h5'
      CALL colm_write_serial (lndname, 'csol_l'//trim(c), a_real8, &
         compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

      write(*,*) 'csol_l' // trim(c) // ' done'

      ! (6) Read in the thermal conductivity of saturated soil [W/m-K]
      lndname = trim(bindir) // 'RAW_DATA_updated/tksatu_l'//trim(c)

      iunit = 100
      inquire (iolength=length) a_real8 (:,1)
      open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
      DO irow = 1, nlat
         read (iunit, rec=irow) a_real8 (:,irow)
      ENDDO
      close (iunit)

      lndname = trim(h5dir) // 'soil/tksatu.h5'
      CALL colm_write_serial (lndname, 'tksatu_l'//trim(c), a_real8, &
         compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

      write(*,*) 'tksatu_l' // trim(c) // ' done'

      ! (7) Read in the thermal conductivity for dry soil [W/(m-K)]
      lndname = trim(bindir) // 'RAW_DATA_updated/tkdry_l'//trim(c)

      iunit = 100
      inquire (iolength=length) a_real8 (:,1)
      open (iunit, file=trim(lndname), access='direct', recl=length, form='unformatted', status='old')
      DO irow = 1, nlat
         read (iunit, rec=irow) a_real8 (:,irow)
      ENDDO
      close (iunit)

      lndname = trim(h5dir) // 'soil/tkdry.h5'
      CALL colm_write_serial (lndname, 'tkdry_l'//trim(c), a_real8, &
         compress, chunk = (/ nlon/nxblk, nlat/nyblk /))

      write(*,*) 'tkdry_l' // trim(c) // ' done'

   ENDDO

   deallocate (a_real8)

END PROGRAM bin_to_hdf5
