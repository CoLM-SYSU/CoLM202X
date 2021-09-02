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
USE GlobalVars

IMPLICIT NONE

! arguments:
      character(len=256), intent(in) :: dir_rawdata 
      character(len=256) lndname
      character(len=1) land_chr1(nlon)
      integer(kind=1)  land_int1(nlon30s)
      integer(kind=2)  land_int2(nlon30s)

      ! (1) global land cover characteristics
      ! ---------------------------------
      integer, allocatable :: landtypes(:,:)  ! GLCC USGS/MODIS IGBP land cover types 

      ! (6) global soil characteristcs
      ! --------------------------
      integer(kind=1), allocatable :: int_soil_grav_l (:,:) ! gravel content   (% of volume)
      integer(kind=1), allocatable :: int_soil_sand_l (:,:) ! sand percentage  (% of weight)
      integer(kind=1), allocatable :: int_soil_clay_l (:,:) ! clay percentage  (% of weight)
      integer(kind=2), allocatable :: int_soil_oc_l   (:,:) ! organic carbon percentage (% of weight)
      integer(kind=2), allocatable :: int_soil_bd_l   (:,:) ! bulk density     (g/cm3)

      ! ---------------------------------------------------------------
      integer i, j, i1, j1
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
      CHARACTER(len=256) suffix
      real(r8) a
      real(r8) soildepth
      integer ii, iii, iiii, jj, jjj, jjjj

! ........................................
! ... (1) gloabl land cover characteristics  
! ........................................
      iunit = 100 
      inquire(iolength=length) land_chr1 
      allocate ( landtypes(nlon,nlat) )

#if(defined USGS_CLASSIFICATION)
      ! GLCC USGS classification
      ! -------------------
      lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/landtypes_usgs_update.bin'
#else
      lndname = trim(dir_rawdata)//'landtypes/landtypes-modis-igbp-2005.bin'
#endif
      
      print*,lndname
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
      do nrow = 1, nlat 
         read(iunit,rec=nrow,err=100) land_chr1 
         landtypes(:,nrow) = ichar(land_chr1(:)) 
      enddo 
      close (iunit)

#ifdef USGS_CLASSIFICATION
      suffix = ''
#else
      suffix = '.igbp'
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

      allocate ( int_soil_grav_l (nlon30s,nlat30s) ,&
                 int_soil_sand_l (nlon30s,nlat30s) ,&
                 int_soil_clay_l (nlon30s,nlat30s) ,&
                 int_soil_oc_l   (nlon30s,nlat30s) ,&
                 int_soil_bd_l   (nlon30s,nlat30s)  )

      allocate ( theta_s_l       (nlon,nlat) ,&
                 psi_s_l         (nlon,nlat) ,&
                 lambda_l        (nlon,nlat) ,&
                 k_s_l           (nlon,nlat)  )

      allocate ( csol_l          (nlon,nlat) ,&
                 tksatu_l        (nlon,nlat) ,&
                 tkdry_l         (nlon,nlat)  )

      iunit = 100
      DO nsl = 1, 8
         MODEL_SOIL_LAYER = nsl
         write(c,'(i1)') MODEL_SOIL_LAYER 

         ! ------------------------------------
         ! (6.1) precentage of gravel (% volume)
         ! ------------------------------------
         inquire(iolength=length) land_int1
         lndname = trim(dir_rawdata)//'soil/GRAV_L'//trim(c)
         print*,lndname
         
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = 1, nlat30s
            read(iunit,rec=nrow,err=100) land_int1
            int_soil_grav_l(:,nrow) = land_int1(:)
         enddo
         close (iunit)

         ! ----------------------------------
         ! (6.2) percentage of sand (% weight)
         ! ----------------------------------
         inquire(iolength=length) land_int1
         lndname = trim(dir_rawdata)//'soil/SAND_L'//trim(c)
         print*,lndname

         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = 1, nlat30s
            read(iunit,rec=nrow,err=100) land_int1
            int_soil_sand_l(:,nrow) = land_int1(:)
         enddo
         close (iunit)

         ! ----------------------------------
         ! (6.3) percentage of clay (% weight)
         ! ----------------------------------
         inquire(iolength=length) land_int1
         lndname = trim(dir_rawdata)//'soil/CLAY_L'//trim(c)
         print*,lndname

         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = 1, nlat30s
            read(iunit,rec=nrow,err=100) land_int1
            int_soil_clay_l(:,nrow) = land_int1(:)
         enddo
         close (iunit)

         ! -------------------------------------
         ! (6.4) percentage of organic carbon (%)
         ! -------------------------------------
         inquire(iolength=length) land_int2
         lndname = trim(dir_rawdata)//'soil/OC_L'//trim(c)
         print*,lndname

         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = 1, nlat30s 
            read(iunit,rec=nrow,err=100) land_int2
            int_soil_oc_l(:,nrow) = land_int2(:)
         enddo
         close (iunit)

         ! -------------------------
         ! (6.5) bulk density (g/cm3)
         ! -------------------------
         inquire(iolength=length) land_int2
         lndname = trim(dir_rawdata)//'soil/BD_L'//trim(c)
         print*,lndname

         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
         do nrow = 1, nlat30s 
            read(iunit,rec=nrow,err=100) land_int2 
            int_soil_bd_l(:,nrow) = land_int2(:)
         enddo
         close (iunit)

         ! ---------
         ! (6.6) ...
         ! ---------



         ! ---------------------------------------
         ! calculate the soil hydraulic parameters
         ! ---------------------------------------

! added by yuan, 06/02/2016
#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP, "soil hydraulic parameters..."
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,i1,j1) &
!$OMP PRIVATE(soil_grav_l,soil_sand_l,soil_clay_l,soil_oc_l,soil_bd_l) &
!$OMP PRIVATE(theta_s,psi_s,lambda,k_s,csol,tksatu,tkdry) &
!$OMP PRIVATE(soildepth,a)
#endif
         do j = 1, nlat
            do i = 1, nlon

               i1 = i
               j1 = j

#ifndef USGS_CLASSIFICATION
               i1 = int((i+1)/2)
               j1 = int((j+1)/2)
#endif
               soil_grav_l = int_soil_grav_l(i1,j1)
               soil_sand_l = int_soil_sand_l(i1,j1)
               soil_clay_l = int_soil_clay_l(i1,j1)
               soil_oc_l   = int_soil_oc_l  (i1,j1) * 0.01
               soil_bd_l   = int_soil_bd_l  (i1,j1) * 0.01

               if(soil_grav_l < 0.0) soil_grav_l = 0.0  ! missing value = -100

#if(defined USGS_CLASSIFICATION)
               if(landtypes(i,j)==16)then   !WATER BODIES(16)
#else
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
#else
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
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

         print*,'theta  =', minval(theta_s_l, mask = theta_s_l .gt. -1.0e30), maxval(theta_s_l, mask = theta_s_l .gt. -1.0e30)
         print*,'psi    =', minval(psi_s_l,   mask = psi_s_l   .gt. -1.0e30), maxval(psi_s_l,   mask = psi_s_l   .gt. -1.0e30)
         print*,'lambda =', minval(lambda_l,  mask = lambda_l  .gt. -1.0e30), maxval(lambda_l,  mask = lambda_l  .gt. -1.0e30)
         print*,'Ks     =', minval(k_s_l,     mask = k_s_l     .gt. -1.0e30), maxval(k_s_l,     mask = k_s_l     .gt. -1.0e30)
         print*,'csol   =', minval(csol_l,    mask = csol_l    .gt. -1.0e30), maxval(csol_l,    mask = csol_l    .gt. -1.0e30)
         print*,'tksatu =', minval(tksatu_l,  mask = tksatu_l  .gt. -1.0e30), maxval(tksatu_l,  mask = tksatu_l  .gt. -1.0e30)
         print*,'tkdry  =', minval(tkdry_l,   mask = tkdry_l   .gt. -1.0e30), maxval(tkdry_l,   mask = tkdry_l   .gt. -1.0e30)

! (1) Write out the saturated water content [cm3/cm3]
         inquire(iolength=length) theta_s_l(:,1)
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/theta_s_l'//trim(c)//trim(suffix)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) theta_s_l(:,j)
         enddo
         close(iunit)

! (2) Write out the matric potential at saturation [cm]
         inquire(iolength=length) psi_s_l(:,1)
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/psi_s_l'//trim(c)//trim(suffix)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) psi_s_l(:,j)
         enddo
         close(iunit)

! (3) Write out the pore size distribution index [dimensionless]
         inquire(iolength=length) lambda_l(:,1)
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/lambda_l'//trim(c)//trim(suffix)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) lambda_l(:,j)
         enddo
         close(iunit)

! (4) Write out the saturated hydraulic conductivity [cm/day]
         inquire(iolength=length) k_s_l(:,1)
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/k_s_l'//trim(c)//trim(suffix)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) k_s_l(:,j)
         enddo
         close(iunit)

! (5) Write out the heat capacity of soil solids [J/(m3 K)]
         inquire(iolength=length) csol_l(:,1)
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/csol_l'//trim(c)//trim(suffix)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) csol_l(:,j)
         enddo
         close(iunit)

! (6) Write out the thermal conductivity of saturated soil [W/m-K]
         inquire(iolength=length) tksatu_l(:,1)
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/tksatu_l'//trim(c)//trim(suffix)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) tksatu_l(:,j)
         enddo
         close(iunit)

! (7) Write out the thermal conductivity for dry soil [W/(m-K)]
         inquire(iolength=length) tkdry_l(:,1)
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/tkdry_l'//trim(c)//trim(suffix) 
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
         do j = 1, nlat
            write(iunit,rec=j,err=101) tkdry_l(:,j)
         enddo
         close(iunit)

      ENDDO
         
      deallocate ( int_soil_grav_l ,&
                   int_soil_sand_l ,&
                   int_soil_clay_l ,&
                   int_soil_oc_l   ,&
                   int_soil_bd_l    )

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

print*,'------ END rd_soil_properties ------'

END SUBROUTINE rd_soil_properties
