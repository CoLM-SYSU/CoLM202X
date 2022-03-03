#include <define.h>
SUBROUTINE rd_soil_properties(dir_rawdata)
! ----------------------------------------------------------------------
! => Read in soil characteristic dataset from original "raw" data files -
!     data with 30 arc seconds resolution
! => Fill the missing data 
! => Estimate the soil hydraulic and thermal parameters at the resolution of 30 arc seconds
!
! 6. Global Soil Characteristics 
!    (http://globalchange.bnu.edu.cn)
! 6.1 percentage of gravel (% volume)
! 6.2 percentage of sand   (% weight)
! 6.3 percentage of clay   (% weight)
! 6.4 organic Carbon (SOC) (% weight)
! 6.5 bulk density (BD)    (g/cm3)
! 6.6 ...
!
! Reference: 
! (1) http://globalchange.bnu.edu.cn
! (2) Shangguan et al., 2014: 
!     A global soil data set for earth system modeling. 
!     J. of Advances in Modeling Earth Systems, DOI: 10.1002/2013MS000293
! (3) Dai et al.,2014: Implementation of a New Global Soil Dataset in the Common Land Model.
!
! Created by Yongjiu Dai, 12/2013
! ----------------------------------------------------------------------
use precision
use spmd_TM
#if (defined usempi)
use spmd_io
#endif

IMPLICIT NONE

! arguments:
      character(len=256), intent(in) :: dir_rawdata 

! local variables:
      integer, parameter :: nlat=21600    ! 180*(60*2)
      integer, parameter :: nlon=43200    ! 360*(60*2)

      character(len=256) lndname

#if (defined usempi)
      character, allocatable :: land_chr1(:,:)
      integer(kind=1), allocatable ::  land_int1(:,:)
      integer(kind=2), allocatable ::  land_int2(:,:)
#else
      character(len=1) land_chr1(nlon)
      character(len=2) land_chr2(nlon)
      integer(kind=1)  land_int1(nlon)
      integer(kind=2)  land_int2(nlon)
#endif

      ! (1) global land cover characteristics
      ! ---------------------------------
      integer, allocatable :: landtypes(:,:)  ! GLCC USGS/MODIS IGBP land cover types 

      ! (6) global soil characteristcs
      ! --------------------------
!     integer, allocatable :: nonsoil (:,:)   !

      integer(kind=1), allocatable :: int_soil_grav_l (:,:) ! gravel content   (% of volume)
      integer(kind=1), allocatable :: int_soil_sand_l (:,:) ! sand percentage  (% of weight)
      integer(kind=1), allocatable :: int_soil_clay_l (:,:) ! clay percentage  (% of weight)
      integer(kind=2), allocatable :: int_soil_oc_l   (:,:) ! organic carbon percentage (% of weight)
      integer(kind=2), allocatable :: int_soil_bd_l   (:,:) ! bulk density     (g/cm3)

!   ---------------------------------------------------------------
      integer i, j
      integer nrow, ncol
      integer iunit    
      integer length
      integer nsl, MODEL_SOIL_LAYER

! soil hydraulic parameters
      real(r8), allocatable :: theta_s_l   (:,:) ! saturated water content (cm3/cm3)
      real(r8), allocatable :: psi_s_l     (:,:) ! matric potential at saturation (cm)
      real(r8), allocatable :: lambda_l    (:,:) ! pore size distribution index (dimensionless)
      real(r8), allocatable :: k_s_l       (:,:) ! saturated hydraulic conductivity (cm/day)

! soil thermal parameters
      real(r8), allocatable :: csol_l      (:,:) ! heat capacity of soil solids [J/(m3 K)]
      real(r8), allocatable :: tksatu_l    (:,:) ! thermal conductivity of saturated soil [W/m-K]
      real(r8), allocatable :: tkdry_l     (:,:) ! thermal conductivity for dry soil  [W/(m-K)]

! CLM soil layer thickiness and depths
      integer nl_soil 
      real(r8), allocatable ::  zsoi(:)  ! soil layer depth [m]
      real(r8), allocatable ::  dzsoi(:) ! soil node thickness [m]
      real(r8), allocatable ::  zsoih(:) ! interface level below a zsoi level [m]

! soil hydraulic and thermal parameters
      real(r8) soil_grav_l  ! gravel content   (% of volume)
      real(r8) soil_sand_l  ! sand percentage  (% of weight)
      real(r8) soil_clay_l  ! clay percentage  (% of weight)
      real(r8) soil_oc_l    ! organic carbon percentage (% of weight)
      real(r8) soil_bd_l    ! bulk density     (g/cm3)

      real(r8) theta_s      ! saturated water content (cm3/cm3)
      real(r8) psi_s        ! matric potential at saturation (cm)
      real(r8) lambda       ! pore size distribution index (dimensionless)
      real(r8) k_s          ! saturated hydraulic conductivity (cm/day)
      real(r8) csol         ! heat capacity of soil solids [J/(m3 K)]
      real(r8) tksatu       ! thermal conductivity of saturated soil [W/m-K]
      real(r8) tkdry        ! thermal conductivity for dry soil  [W/(m-K)]

      character c
      real(r8) a
      real(r8) soildepth
      integer ii, iii, iiii, jj, jjj, jjjj

#if(defined FAO_STATSGO_SOILMAP)
! ----------------------------------------------------------------------
! Creates land model surface dataset from original "raw" data files -
!     USGS data with 30 arc seconds resolution:
!  -  soil texture (FAO+STATSGO)
! ----------------------------------------------------------------------
      character(LEN=1), allocatable :: landsola_chr(:,:) ! soil texture type in upper 30cm
      character(LEN=1), allocatable :: landsolb_chr(:,:) ! soil texture type in 30-100cm

    ! relative amounts of sand (s), and clay (c) in the < 2 mm fraction of
    ! the component layer was then estimated using table:
      real(r8), dimension(17) ::  s_ = (/92.,82.,58.,17.,10.,43.,58.,10.,32.,52.,&
                                          6.,22., 0., 0., 0., 0., 0./)
      real(r8), dimension(17) ::  c_ = (/ 3., 6.,10.,13., 5.,18.,27.,34.,34.,42.,&
                                         47.,58., 0., 0., 0., 0., 0./)
      !integer ia
      real(r8) :: a_
      integer  :: L1, L2
#endif

      real(r8) r8_min,   r8_max
      real(r8) r8_min_g, r8_max_g

      integer  :: nrow_start
      integer  :: nrow_end

#if (defined usempi)
      integer(kind=MPI_OFFSET_KIND) :: fdisp
#endif

! Initialize MPI tasks
#if (defined usempi)
      nrow_start = fine_lat_map%bdisp(1) + 1
      nrow_end   = fine_lat_map%bdisp(1) + fine_lat_map%bstrd(1) 
#else
      nrow_start = 1
      nrow_end   = nlat
#endif

! ........................................
! ... (1) gloabl land cover characteristics  
! ........................................
      allocate (landtypes(nlon,nrow_start:nrow_end))
#if(defined USGS_CLASSIFICATION)
     ! GLCC USGS classification
      lndname = trim(dir_rawdata)//'RAW_DATA_updated/landtypes_usgs_update.bin'
#endif
#if(defined IGBP_CLASSIFICATION)
     ! MODIS IGBP classification
      lndname = trim(dir_rawdata)//'RAW_DATA_updated/landtypes_igbp_update.bin'
#endif

      if (p_master) print*,trim(lndname)

#if (defined usempi)
      allocate (land_chr1(nlon,nrow_start:nrow_end))
      fdisp = 0
      call mpi_rdwr_data (lndname, 'read', fdisp, fine_lat_map, nlon, land_chr1)
      landtypes(:,:) = ichar(land_chr1(:,:))
      deallocate (land_chr1)
#else
      inquire(iolength=length) land_chr1
      iunit = 100
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      do nrow = 1, nlat
         read(iunit,rec=nrow,err=100) land_chr1
! modifiedy by yuan, 06/02/2016
         !do ncol = 1, nlon
         !   landtypes(ncol,nrow) = ichar(land_chr1(ncol))
         !enddo
         landtypes(:,nrow) = ichar(land_chr1(:))
      enddo
      close (iunit)
#endif 

! .................................
! ... (6) global soil charateristics
! .................................
      nl_soil = 10
      allocate ( zsoi(1:nl_soil), dzsoi(1:nl_soil), zsoih(0:nl_soil) )

      ! ----------------------------------
      ! soil layer thickness, depths (m)
      ! ----------------------------------
      do nsl = 1, nl_soil
        zsoi(nsl) = 0.025*(exp(0.5*(nsl-0.5))-1.)  ! node depths
      end do

      dzsoi(1) = 0.5*(zsoi(1)+zsoi(2))         ! =zsoih(1)
      dzsoi(nl_soil) = zsoi(nl_soil)-zsoi(nl_soil-1)
      do nsl = 2, nl_soil-1
         dzsoi(nsl) = 0.5*(zsoi(nsl+1)-zsoi(nsl-1))  ! thickness b/n two interfaces
      end do

      zsoih(0) = 0.
      zsoih(nl_soil) = zsoi(nl_soil) + 0.5*dzsoi(nl_soil)
      do nsl = 1, nl_soil-1
         zsoih(nsl) = 0.5*(zsoi(nsl)+zsoi(nsl+1))    ! interface depths
      enddo

! -----------------------
! non-soil classification
!     NONSOIL
!     ----------------------
!     CODE   VALUE
!     -19    Inland water
!     -18    Urban
!     -17    Salt flats
!     -16    Rock debris
!     -15    No data
!     -14    Island
!     -13    Humanly disturbed
!     -12    Glaciers & permanent snow
!     -11    Fishponds
!     -10    Dunes & shifting sands
!      1     Soil
!     ----------------------
!
!     iunit = 100
!     inquire(iolength=length) land_int1
!     lndname = trim(dir_rawdata)//'soil/NONSOIL'
!     print*,lndname
!
!     allocate ( nonsoil(nlon,nlat) )
!
!     ii = 0
!     iii = 0
!     iiii = 0
!
!     open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
!     do nrow = 1, nlat
!        read(iunit,rec=nrow,err=100) land_int1
!        do ncol = 1, nlon
!           nonsoil(ncol,nrow) = land_int1(ncol)
!           if(nonsoil(ncol,nrow) == -16)then
!              ii = ii + 1
!           endif
!           if(nonsoil(ncol,nrow) == -14)then
!              iii = iii + 1
!           endif
!           if(nonsoil(ncol,nrow) == -18)then
!              iiii = iiii + 1
!           endif
!        enddo
!     enddo
!     print*, minval(nonsoil), maxval(nonsoil)
!     print*,'Rock debris =', ii, 'island = ', iii, 'Urban =', iiii
!     close (iunit)
! -----------------------

#if(defined FAO_STATSGO_SOILMAP)
! -----------------FAO/+STATSGO soil map-----------------------------------------
! .. (1) soil 0 - 30 cm;
! .. (2) soil 30 - 100 cm;

      allocate ( landsola_chr (nlon,nrow_start:nrow_end) )
      allocate ( landsolb_chr (nlon,nrow_start:nrow_end) )

      lndname = trim(dir_rawdata)//'USGS_soil/soilcat.30s'
      if (p_master) print*,trim(lndname)
#if (defined usempi)
      fdisp = 0
      call mpi_rdwr_data (lndname, 'read', fdisp, fine_lat_map, nlon, landsola_chr)
#else
      inquire(iolength=length) landsola_chr1
      iunit = 100
      open(iunit, file=lndname,access='direct',recl=length,&
                    form='unformatted',status='old')
      do nrow = 1, nlat
         read(iunit,rec=nrow,err=100) landsola_chr(1:nlon,nrow)
      enddo
      close(iunit)
#endif

      lndname = trim(dir_rawdata)//'USGS_soil/soilcatb.30s'
      if (p_master) print*,trim(lndname)
#if (defined usempi)
      fdisp = 0
      call mpi_rdwr_data (lndname, 'read', fdisp, fine_lat_map, nlon, landsolb_chr)
#else
      inquire(iolength=length) landsolb_chr1
      iunit = 100
      open(iunit, file=lndname,access='direct',recl=length,&
                    form='unformatted',status='old')
      do nrow = 1, nlat
         read(iunit,rec=nrow,err=100) landsolb_chr(1:nlon,nrow)
      enddo
      close(iunit)

#endif
#endif

! modifiedy by yuan, 06/02/2016
#if (!defined FAO_STATSGO_SOILMAP)
      allocate ( int_soil_grav_l (nlon,nrow_start:nrow_end) ,&
                 int_soil_sand_l (nlon,nrow_start:nrow_end) ,&
                 int_soil_clay_l (nlon,nrow_start:nrow_end) ,&
                 int_soil_oc_l   (nlon,nrow_start:nrow_end) ,&
                 int_soil_bd_l   (nlon,nrow_start:nrow_end)  )
#endif

      allocate ( theta_s_l       (nlon,nrow_start:nrow_end) ,&
                 psi_s_l         (nlon,nrow_start:nrow_end) ,&
                 lambda_l        (nlon,nrow_start:nrow_end) ,&
                 k_s_l           (nlon,nrow_start:nrow_end)  )

      allocate ( csol_l          (nlon,nrow_start:nrow_end) ,&
                 tksatu_l        (nlon,nrow_start:nrow_end) ,&
                 tkdry_l         (nlon,nrow_start:nrow_end)  )

      DO nsl = 1, 8
         MODEL_SOIL_LAYER = nsl
         write(c,'(i1)') MODEL_SOIL_LAYER 

! modifiedy by yuan, 06/02/2016
         !allocate ( int_soil_grav_l (nlon,nlat) ,&
         !           int_soil_sand_l (nlon,nlat) ,&
         !           int_soil_clay_l (nlon,nlat) ,&
         !           int_soil_oc_l   (nlon,nlat) ,&
         !           int_soil_bd_l   (nlon,nlat)  )

         !allocate ( theta_s_l       (nlon,nlat) ,&
         !           psi_s_l         (nlon,nlat) ,&
         !           lambda_l        (nlon,nlat) ,&
         !           k_s_l           (nlon,nlat)  )

         !allocate ( csol_l          (nlon,nlat) ,&
         !           tksatu_l        (nlon,nlat) ,&
         !           tkdry_l         (nlon,nlat)  )

#if (!defined FAO_STATSGO_SOILMAP)
         ! ------------------------------------
         ! (6.1) precentage of gravel (% volume)
         ! ------------------------------------
#if (defined SoilGrid_Rock_Fragments)
         lndname = trim(dir_rawdata)//'SoilGrids_Gravel/GRAV_L'//trim(c)
#else
         lndname = trim(dir_rawdata)//'soil/GRAV_L'//trim(c)
#endif
         if (p_master) print*,trim(lndname)
         
#if (defined usempi)
         fdisp = 0
         call mpi_rdwr_data (lndname, 'read', fdisp, fine_lat_map, nlon, int_soil_grav_l)
#else
         inquire(iolength=length) land_int1
         iunit = 100
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = 1, nlat
            read(iunit,rec=nrow,err=100) land_int1
! modifiedy by yuan, 06/02/2016
            !do ncol = 1, nlon
            !   int_soil_grav_l(ncol,nrow) = land_int1(ncol)
            !enddo
            int_soil_grav_l(:,nrow) = land_int1(:)
         enddo
         close (iunit)
#endif

         ! ----------------------------------
         ! (6.2) percentage of sand (% weight)
         ! ----------------------------------
         lndname = trim(dir_rawdata)//'soil/SAND_L'//trim(c)
         if (p_master) print*,trim(lndname)

#if (defined usempi)
         fdisp = 0
         call mpi_rdwr_data (lndname, 'read', fdisp, fine_lat_map, nlon, int_soil_sand_l)
#else
         inquire(iolength=length) land_int1
         iunit = 100
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = 1, nlat
            read(iunit,rec=nrow,err=100) land_int1
! modifiedy by yuan, 06/02/2016
            !do ncol = 1, nlon
            !   int_soil_sand_l(ncol,nrow) = land_int1(ncol)
            !enddo
            int_soil_sand_l(:,nrow) = land_int1(:)
         enddo
         close (iunit)
#endif

         ! ----------------------------------
         ! (6.3) percentage of clay (% weight)
         ! ----------------------------------
         lndname = trim(dir_rawdata)//'soil/CLAY_L'//trim(c)
         if (p_master) print*,trim(lndname)

#if (defined usempi)
         fdisp = 0
         call mpi_rdwr_data (lndname, 'read', fdisp, fine_lat_map, nlon, int_soil_clay_l)
#else
         inquire(iolength=length) land_int1
         iunit = 100
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = 1, nlat 
            read(iunit,rec=nrow,err=100) land_int1
! modifiedy by yuan, 06/02/2016
            !do ncol = 1, nlon
            !   int_soil_clay_l(ncol,nrow) = land_int1(ncol)
            !enddo
            int_soil_clay_l(:,nrow) = land_int1(:)
         enddo
         close (iunit)
#endif

         ! -------------------------------------
         ! (6.4) percentage of organic carbon (%)
         ! -------------------------------------
         lndname = trim(dir_rawdata)//'soil/OC_L'//trim(c)
         if (p_master) print*,trim(lndname)

#if (defined usempi)
         fdisp = 0
         call mpi_rdwr_data (lndname, 'read', fdisp, fine_lat_map, nlon, int_soil_oc_l)
#else
         inquire(iolength=length) land_int2
         iunit = 100
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = 1, nlat 
            read(iunit,rec=nrow,err=100) land_int2
! modifiedy by yuan, 06/02/2016
            !do ncol = 1, nlon
            !   int_soil_oc_l(ncol,nrow) = land_int2(ncol)
            !enddo
            int_soil_oc_l(:,nrow) = land_int2(:)
         enddo
         close (iunit)
#endif

         ! -------------------------
         ! (6.5) bulk density (g/cm3)
         ! -------------------------
         lndname = trim(dir_rawdata)//'soil/BD_L'//trim(c)
         if (p_master) print*,trim(lndname)

#if (defined usempi)
         fdisp = 0
         call mpi_rdwr_data (lndname, 'read', fdisp, fine_lat_map, nlon, int_soil_bd_l)
#else
         inquire(iolength=length) land_int2
         iunit = 100
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
         do nrow = 1, nlat 
            read(iunit,rec=nrow,err=100) land_int2 
! modifiedy by yuan, 06/02/2016
            !do ncol = 1, nlon 
            !   int_soil_bd_l(ncol,nrow) = land_int2(ncol)
            !enddo 
            int_soil_bd_l(:,nrow) = land_int2(:)
         enddo
         close (iunit)
#endif

         ! ---------
         ! (6.6) ...
         ! ---------
#endif



         ! ---------------------------------------
         ! calculate the soil hydraulic parameters
         ! ---------------------------------------
         do j = nrow_start, nrow_end
            do i = 1, nlon

#if(defined FAO_STATSGO_SOILMAP)
! -----------------FAO/+STATSGO soil map-----------------------------------------
            L1 = ichar(landsola_chr(i,j))
            L2 = ichar(landsolb_chr(i,j))

            soil_grav_l = 0.0
            soil_sand_l = 0.0
            soil_clay_l = 0.0
            soil_oc_l = 0.0
            soil_bd_l = 0.0

            if(zsoih(nsl+1) .lt. 0.3)then
               if(L1.ge.1 .and. L1.le.17) then
                  soil_grav_l = 0.0
                  soil_sand_l = s_(L1)
                  soil_clay_l = c_(L1)
                  soil_oc_l = 0.0
                                   a_ = 0.489 - 0.00126*s_(L1)
                  soil_bd_l = (1.- a_)*2.7
               endif
            else
               if(L2.ge.1 .and. L2.le.17) then
                  soil_grav_l = 0.0
                  soil_sand_l = s_(L2)
                  soil_clay_l = c_(L2)
                  soil_oc_l = 0.0
                                   a_ = 0.489 - 0.00126*s_(L2)
                  soil_bd_l = (1.- a_)*2.7
               endif
            endif
#if(defined USGS_CLASSIFICATION)
             if(landtypes(i,j)/=0 .and. landtypes(i,j)/=16 .and. & ! NOT OCEAN(0) / LAND WATER BODIES(16)
                 landtypes(i,j)/=24) then ! NOT SNOW or ICE (24)
#endif
#if(defined IGBP_CLASSIFICATION)
             if(landtypes(i,j)/=0 .and. & !NOT OCEAN(0)
                landtypes(i,j)/=17 .and. landtypes(i,j)/=15)then !NOT LAND WATER BODIES(17)/ SNOW OR ICE(15)
#endif
                if (soil_sand_l*soil_clay_l < 0.01) then
                   soil_sand_l = 43.
                   soil_clay_l = 18.
                end if

             end if
! -----------------FAO/+STATSGO soil map-----------------------------------------
#else
               soil_grav_l = int_soil_grav_l(i,j)
               soil_sand_l = int_soil_sand_l(i,j)
               soil_clay_l = int_soil_clay_l(i,j)
               soil_oc_l   = int_soil_oc_l  (i,j) * 0.01
               soil_bd_l   = int_soil_bd_l  (i,j) * 0.01

               if(soil_grav_l < 0.0) soil_grav_l = 0.0  ! missing value = -100
#endif

!#if(defined USGS_CLASSIFICATION)
!              if(landtypes(i,j)/=0 .and. & !NOT OCEAN(0)/WATER BODIES(16)/GLACIER and ICESHEET(24)
!                 landtypes(i,j)/=16 .and. landtypes(i,j)/=24)then
!#endif
!#if(defined IGBP_CLASSIFICATION)
!              if(landtypes(i,j)/=0 .and. & !NOT OCEAN(0)/WATER BODIES(17)/GLACIER and ICE SHEET(15)
!                 landtypes(i,j)/=17 .and. landtypes(i,j)/=15)then
!#endif

#if(defined USGS_CLASSIFICATION)
               if(landtypes(i,j)==16)then   !WATER BODIES(16)
#endif
#if(defined IGBP_CLASSIFICATION)
               if(landtypes(i,j)==17)then   !WATER BODIES(17)
#endif
                  soil_grav_l = 0.
                  soil_sand_l = 10.
                  soil_clay_l = 45.
                  soil_oc_l   = 3.0
                  soil_bd_l   = 1.2
               endif

#if(defined USGS_CLASSIFICATION)
               if(landtypes(i,j)==24)then   !GLACIER and ICESHEET(24)
#endif
#if(defined IGBP_CLASSIFICATION)
               if(landtypes(i,j)==15)then   !GLACIER and ICE SHEET(15)
#endif
                  soil_grav_l = 90.
                  soil_sand_l = 89.
                  soil_clay_l = 1.
                  soil_oc_l   = 0.
                  soil_bd_l   = 2.0
               endif

               if(landtypes(i,j)/=0)then    !NOT OCEAN(0)
                  ! checking the soil physical properties
                  ! ------------------------------------
                  if( soil_sand_l < 0.0 ) soil_sand_l = 43.   ! missing value = -100
                  if( soil_clay_l < 0.0 ) soil_clay_l = 18.   ! missing value = -100
                  if( soil_oc_l   < 0.0 ) soil_oc_l = 1.0     ! missing value = -999
                  if( soil_bd_l   < 0.0 ) soil_bd_l = 1.2     ! missing value = -999

                  if( soil_sand_l < 1.0 ) soil_sand_l = 1.
                  if( soil_clay_l < 1.0 ) soil_clay_l = 1.

                  a = soil_sand_l + soil_clay_l
                  if( a >= 96. ) then
                      soil_sand_l = soil_sand_l * 96. / a
                      soil_clay_l = soil_clay_l * 96. / a
                  endif
                  if( soil_oc_l < 0.01 ) soil_oc_l = 0.01
                  if( soil_oc_l > 58.0 ) soil_oc_l = 58.0
                  if( soil_bd_l < 0.1  ) soil_bd_l = 0.1 

                  soildepth = zsoih(MODEL_SOIL_LAYER)*100.0 

                 ! estimating soil hydraulic properties
                 ! ------------------------------------
                  CALL soil_hydraulic_parameters( soil_sand_l,soil_clay_l, &
                       soil_oc_l,soil_bd_l,soildepth, &
                       theta_s,psi_s,lambda,k_s )

                 ! estimating soil thermal properties
                 ! ----------------------------------
                  CALL soil_thermal_parameters(soil_grav_l,soil_sand_l,soil_clay_l, &
                       soil_oc_l,soil_bd_l,theta_s,soildepth,&
                       csol,tksatu,tkdry)

                 ! updating the hydraulic properties of soil-gravel mixtures
                 ! -------------------------------------------------------
                  theta_s = theta_s * (1.-soil_grav_l/100.)
                  k_s = k_s * (1.-soil_grav_l/100.)
                 ! k_s = k_s * 2.0*(1.-soil_grav_l/100.)/(2.0+soil_grav_l/100.)  ! (Peck and Waston (1979)
                 ! k_s = k_s * (1.-soil_grav_mass_fraction/100.)  ! (Brakensiek et al., 1986; Bagarello and Iovino, 2007)

               else                         !OCEAN
                  theta_s = -1.0e36
                  psi_s   = -1.0e36
                  lambda  = -1.0e36
                  k_s     = -1.0e36

                  csol    = -1.0e36
                  tksatu  = -1.0e36
                  tkdry   = -1.0e36
               endif

               theta_s_l(i,j) = theta_s
               psi_s_l  (i,j) = psi_s
               lambda_l (i,j) = lambda
               k_s_l    (i,j) = k_s

               csol_l   (i,j) = csol
               tksatu_l (i,j) = tksatu
               tkdry_l  (i,j) = tkdry

            enddo 
         enddo

         call print_minmax ('theta  =', theta_s_l, theta_s_l .gt. -1.0e30)
         call print_minmax ('psi    =', psi_s_l  , psi_s_l   .gt. -1.0e30)
         call print_minmax ('lambda =', lambda_l , lambda_l  .gt. -1.0e30)
         call print_minmax ('Ks     =', k_s_l    , k_s_l     .gt. -1.0e30)
         call print_minmax ('csol   =', csol_l   , csol_l    .gt. -1.0e30)
         call print_minmax ('tksatu =', tksatu_l , tksatu_l  .gt. -1.0e30)
         call print_minmax ('tkdry  =', tkdry_l  , tkdry_l   .gt. -1.0e30)

! (1) Write out the saturated water content [cm3/cm3]
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/theta_s_l'//trim(c)
         if (p_master) print*,trim(lndname)
#if (defined usempi)
         fdisp = 0
         call mpi_rdwr_data (lndname, 'write', fdisp, fine_lat_map, nlon, theta_s_l)
#else
         inquire(iolength=length) theta_s_l(:,1)
         iunit = 100
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) theta_s_l(:,j)
         enddo
         close(iunit)
#endif

! (2) Write out the matric potential at saturation [cm]
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/psi_s_l'//trim(c)
         if (p_master) print*,trim(lndname)
#if (defined usempi)
         fdisp = 0
         call mpi_rdwr_data (lndname, 'write', fdisp, fine_lat_map, nlon, psi_s_l)
#else
         inquire(iolength=length) psi_s_l(:,1)
         iunit = 100
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) psi_s_l(:,j)
         enddo
         close(iunit)
#endif

! (3) Write out the pore size distribution index [dimensionless]
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/lambda_l'//trim(c)
         if (p_master) print*,trim(lndname)
#if (defined usempi)
         fdisp = 0
         call mpi_rdwr_data (lndname, 'write', fdisp, fine_lat_map, nlon, lambda_l)
#else
         inquire(iolength=length) lambda_l(:,1)
         iunit = 100
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) lambda_l(:,j)
         enddo
         close(iunit)
#endif

! (4) Write out the saturated hydraulic conductivity [cm/day]
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/k_s_l'//trim(c)
         if (p_master) print*,trim(lndname)
#if (defined usempi)
         fdisp = 0
         call mpi_rdwr_data (lndname, 'write', fdisp, fine_lat_map, nlon, k_s_l)
#else
         inquire(iolength=length) k_s_l(:,1)
         iunit = 100
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) k_s_l(:,j)
         enddo
         close(iunit)
#endif

! (5) Write out the heat capacity of soil solids [J/(m3 K)]
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/csol_l'//trim(c)
         if (p_master) print*,trim(lndname)
#if (defined usempi)
         fdisp = 0
         call mpi_rdwr_data (lndname, 'write', fdisp, fine_lat_map, nlon, csol_l)
#else
         inquire(iolength=length) csol_l(:,1)
         iunit = 100
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) csol_l(:,j)
         enddo
         close(iunit)
#endif

! (6) Write out the thermal conductivity of saturated soil [W/m-K]
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/tksatu_l'//trim(c)
         if (p_master) print*,trim(lndname)
#if (defined usempi)
         fdisp = 0
         call mpi_rdwr_data (lndname, 'write', fdisp, fine_lat_map, nlon, tksatu_l)
#else
         inquire(iolength=length) tksatu_l(:,1)
         iunit = 100
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) tksatu_l(:,j)
         enddo
         close(iunit)
#endif

! (7) Write out the thermal conductivity for dry soil [W/(m-K)]
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/tkdry_l'//trim(c)
         if (p_master) print*,trim(lndname)
#if (defined usempi)
         fdisp = 0
         call mpi_rdwr_data (lndname, 'write', fdisp, fine_lat_map, nlon, tkdry_l)
#else
         inquire(iolength=length) tkdry_l(:,1)
         iunit = 100
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) tkdry_l(:,j)
         enddo
         close(iunit)
#endif

! modifiedy by yuan, 06/02/2016
         !deallocate ( int_soil_grav_l ,&
         !             int_soil_sand_l ,&
         !             int_soil_clay_l ,&
         !             int_soil_oc_l   ,&
         !             int_soil_bd_l    )

         !deallocate ( theta_s_l       ,&
         !             psi_s_l         ,&
         !             lambda_l        ,&
         !             k_s_l            )

         !deallocate ( csol_l          ,&
         !             tksatu_l        ,&
         !             tkdry_l          )

      ENDDO
         
! modifiedy by yuan, 06/02/2016
#if(defined FAO_STATSGO_SOILMAP)
      deallocate ( landsola_chr, &
                   landsolb_chr     )
#else
      deallocate ( int_soil_grav_l ,&
                   int_soil_sand_l ,&
                   int_soil_clay_l ,&
                   int_soil_oc_l   ,&
                   int_soil_bd_l    )
#endif

      deallocate ( theta_s_l       ,&
                   psi_s_l         ,&
                   lambda_l        ,&
                   k_s_l            )

      deallocate ( csol_l          ,&
                   tksatu_l        ,&
                   tkdry_l          )


      deallocate ( landtypes )
      deallocate ( zsoi, dzsoi, zsoih )

      go to 1000
100   print 102,nrow,lndname
101   print 102,j,lndname
102   format(' record =',i8,',  error occured on file: ',a50)
1000  continue

if (p_master) print*,'------ END rd_soil_properties ------'

END SUBROUTINE rd_soil_properties



!#if(defined FAO_STATSGO_SOILMAP)
!  integer function ia(chr,n,ispval)
!
!!  purpose: to convert a n-bytes character (chr) to integer ia.
!!        ** the integer data file is saved as a n-byte character
!!           data file. this function is used to recover the
!!           character data to the integer data.
!!
!!  n      --- the number of bytes in chr
!!  ispval --- default value for the negative integer.
!
!      character*(*) chr
!      integer bit_1, bit_2
!
!      bit_1 = '200'O     ! BINARY '10000000'
!      bit_2 = '377'O     ! BINARY '11111111'
!      ia    = 0
!
!      ii1 = ichar(chr(1:1))
!! .. get the sign -- isn=0 positive, isn=1 negative:
!      jj  = iand(ii1,bit_1)
!      isn = ishft(jj,-7)
!
!! .. for negative number:
!!    because the negative integers are represented by the supplementary
!!    binary code inside machine.
!
!        if (isn.eq.1) then
!          do m = n+1,4
!             nbit = (m-1)*8
!             jj = ishft(bit_2,nbit)
!             ia = ieor(jj,ia)
!          end do
!        endif
!
!!   .. get the byte from chr:
!         do m = 1,n
!           ii2 = ichar(chr(m:m))
!! new IBM xlf 8.1 compiler fix: thanks to Jim Edwards
!           if (ii2.lt.0) ii2 = ii2 + 256
!           mshft = (n-m)*8
!           ia2   = ishft(ii2,mshft)
!!   .. the abs(integer):
!           ia = ieor(ia,ia2)
!         end do
!
!      if (ia.lt.0) ia = ispval
!
!      return
!  end function ia
!
!!-----------------------------------------------------------------------
!!EOP
!#endif
