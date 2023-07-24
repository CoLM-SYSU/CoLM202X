#include <define.h>

MODULE MOD_UserSpecifiedForcing

!DESCRIPTION
!===========
   !---This MODULE is used for read atmospheric forcing dataset from various sources.
   ! ------------------------------------------------------------
   !     Read forcing data from :
   !     1)  PRINCETON     2)  GSWP2         3)  GSWP3
   !     4)  QIAN          5)  CRUNCEPV4     6)  CRUNCEPV7
   !     7)  ERA5LAND      8)  ERA5          9)  MSWX
   !     10) WFDE5         11) CRUJRA        12) WFDEI
   !     13) JRA55         14) GDAS          15) CLDAS
   !     16) CMFD          17) TPMFD         18) CMIP6
   !     19) POINT
   !
   !     PLEASE modify the following codes when specified forcing used
   ! ------------------------------------------------------------
!Original Author:
!-------------------
   !---Shupeng Zhang and Zhongwang Wei

!References:
!-------------------
   !---In preparation


!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------
   !* :SUBROUTINE:"init_user_specified_forcing" : initialization of the selected forcing dataset
   !* :SUBROUTINE:"metfilename"  :  identify the forcing file name
   !* :SUBROUTINE:"metpreprocess" :  preprocess the forcing data

!REVISION HISTORY
   !----------------
   ! 2023.05.01  Shupeng Zhang and Zhongwang Wei @ SYSU
   ! 2021.12.02  Shupeng Zhang and Zhongwang Wei @ SYSU
   ! Siguang Zhu and Nan Wei, 10/2014: metpreprocess for forc_q calibration
   ! Hua Yuan, 04/2014: initial code of forcing structure for CoLM2014

   use MOD_Precision

   implicit none

   character(len=256) :: dataset

   logical  :: solarin_all_band   ! whether solar radiation in all bands is available
   real(r8) :: HEIGHT_V           ! observation height of wind speed
   real(r8) :: HEIGHT_T           ! observation height of air temperature
   real(r8) :: HEIGHT_Q           ! observation height of specific humidity

   integer  :: NVAR      ! variable number of forcing data
   integer  :: startyr   ! start year of forcing data        <MARK #1>
   integer  :: startmo   ! start month of forcing data
   integer  :: endyr     ! end year of forcing data
   integer  :: endmo     ! end month of forcing data

   integer, allocatable :: dtime(:)          ! time interval of forcing data
   integer, allocatable :: offset(:)         ! offset of forcing data

   logical :: leapyear   ! leapyear calendar
   logical :: data2d     ! data in 2 dimension (lon, lat)
   logical :: hightdim   ! have "z" dimension
   logical :: dim2d      ! lat/lon value in 2 dimension (lon, lat)

   character(len=256) :: latname                   ! dimension name of latitude
   character(len=256) :: lonname                   ! dimension name of longitude

   character(len=256) :: groupby                   ! file grouped by year/month

   character(len=256), allocatable :: fprefix(:)   ! file prefix
   character(len=256), allocatable :: vname(:)     ! variable name
   character(len=256), allocatable :: tintalgo(:)  ! interpolation algorithm

   ! ----- public subroutines -----
   public :: init_user_specified_forcing ! initialization of the selected forcing dataset
   public :: metfilename                 ! identify the forcing file name
   public :: metpreprocess               ! preprocess the forcing data

CONTAINS

   ! ----------------
   subroutine init_user_specified_forcing

      use MOD_Namelist
      implicit none

      ! Local variables
      integer :: ivar,NVAR_default

      NVAR = DEF_forcing%NVAR
      NVAR_default=NVAR
      if (DEF_USE_CBL_HEIGHT) then
         NVAR=NVAR+1
      endif

      IF (allocated(dtime )) deallocate(dtime)
      IF (allocated(offset)) deallocate(offset)
      allocate (dtime  (NVAR))
      allocate (offset (NVAR))

      IF (allocated(fprefix )) deallocate(fprefix )
      IF (allocated(vname   )) deallocate(vname   )
      IF (allocated(tintalgo)) deallocate(tintalgo)
      allocate (fprefix  (NVAR))
      allocate (vname    (NVAR))
      allocate (tintalgo (NVAR))

      solarin_all_band = DEF_forcing%solarin_all_band ! whether solar radiation in all bands is available
      HEIGHT_V         = DEF_forcing%HEIGHT_V         ! observation height of wind speed
      HEIGHT_T         = DEF_forcing%HEIGHT_T         ! observation height of air temperature
      HEIGHT_Q         = DEF_forcing%HEIGHT_Q         ! observation height of specific humidity

      startyr          = DEF_forcing%startyr          ! start year of forcing data
      startmo          = DEF_forcing%startmo          ! start month of forcing data
      endyr            = DEF_forcing%endyr            ! end year of forcing data
      endmo            = DEF_forcing%endmo            ! end month of forcing data
      dtime(:)         = DEF_forcing%dtime(:)         ! time interval of forcing data
      offset(:)        = DEF_forcing%offset(:)        ! offset of forcing data

      leapyear         = DEF_forcing%leapyear         ! whether leapyear calendar
      data2d           = DEF_forcing%data2d           ! whether data in 2 dimension (lon, lat)
      hightdim         = DEF_forcing%hightdim         ! whether have "z" dimension (height)
      dim2d            = DEF_forcing%dim2d            ! whether lat/lon value in 2 dimension (lon, lat)

      latname          = DEF_forcing%latname          ! dimension name of latitude
      lonname          = DEF_forcing%lonname          ! dimension name of longitude

      groupby          = DEF_forcing%groupby          ! file grouped by year/month

      do ivar = 1, NVAR_default
         fprefix (ivar) = DEF_forcing%fprefix(ivar)  ! file prefix
         vname   (ivar) = DEF_forcing%vname(ivar)    ! variable name
         tintalgo(ivar) = DEF_forcing%tintalgo(ivar) ! interpolation algorithm
      end do
      if (DEF_USE_CBL_HEIGHT) then
         fprefix (NVAR) = DEF_forcing%CBL_fprefix
         vname   (NVAR) = DEF_forcing%CBL_vname
         tintalgo(NVAR) = DEF_forcing%CBL_tintalgo
         dtime(NVAR)    = DEF_forcing%CBL_dtime
         offset(NVAR)   = DEF_forcing%CBL_offset
      endif
   end subroutine init_user_specified_forcing

   ! ----------------
   FUNCTION metfilename(year, month, day, var_i)

      use MOD_Namelist
      implicit none

      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: day
      integer, intent(in) :: var_i
      character(len=256)  :: metfilename
      character(len=256)  :: yearstr
      character(len=256)  :: monthstr

      write(yearstr, '(I4.4)') year
      write(monthstr, '(I2.2)') month

      select case (trim(DEF_forcing%dataset))
      case ('PRINCETON') ! Princeton forcing data
      !DESCRIPTION
      !===========
         !---Princeton Global Meteorological Forcing Dataset for Land Surface Modeling

      !data source:
      !-------------------
         !---https://rda.ucar.edu/datasets/ds314.0/

      !References:
      !-------------------
         !---Sheffield, J., G. Goteti, and E. F. Wood, 2006: Development of a 50-year high-resolution
         !   global dataset of meteorological forcings for land surface modeling J. Climate, 19(13),
         !   3088-3111.

      !REVISION HISTORY
      !----------------
         !---2022.05.01   Zhongwang Wei @ SYSU: remove the "z" dimension

         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//'-'//trim(yearstr)//'.nc'
      case ('GSWP3')     ! GSWP3 forcing data
      !DESCRIPTION
      !===========
         !---Global Meteorological Forcing Dataset for Global Soil Wetness Project Phase 3

      !data source:
      !-------------------
         !---http://hydro.iis.u-tokyo.ac.jp/GSWP3/
         !---https://www.isimip.org/gettingstarted/input-data-bias-adjustment/details/4/

      !References:
      !-------------------
         !---Dirmeyer, P. A., Gao, X., Zhao, M., Guo, Z., Oki, T. and Hanasaki, N. (2006) GSWP-2:
         !   Multimodel Analysis and Implications for Our Perception of the Land Surface. Bulletin
         !   of the American Meteorological Society, 87(10), 1381–98.

      !REVISION HISTORY
      !----------------
         !---
         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//'-'//trim(monthstr)//'.nc'
      case ('QIAN')      ! Qian forcing data
      !DESCRIPTION
      !===========
         !---Qian Global Meteorological Forcing Dataset from 1948 to 2004

      !data source:
      !-------------------
         !---Not available now!

      !References:
      !-------------------
         !---Qian T., and co-authors, 2006: Simulation of Global Land Surface Conditions from 1948 to 2004.
         !   Part I: Forcing Data and Evaluations. J. Hydrometeorol., 7, 953-975.

      !REVISION HISTORY
      !----------------
         !---

         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//'-'//trim(monthstr)//'.nc'
      case ('CRUNCEPV4') ! CRUNCEP V4 forcing data
      !DESCRIPTION
      !===========
         !---CRUNCEP Version 4 - Atmospheric Forcing Data for the Community Land Model

      !data source:
      !-------------------
         !---http://dods.extra.cea.fr/data/p529viov/cruncep/V5_1901_2013/

      !References:
      !-------------------
         !---Viovy, N. (2010), CRU‐NCEP dataset.
         !   [Available at: http://dods.extra.cea.fr/data/p529viov/cruncep/readme.htm.]

      !REVISION HISTORY
      !----------------
         !---

         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//'-'//trim(monthstr)//'.nc'
      case ('CRUNCEPV7') ! CRUNCEP V7 forcing data
      !DESCRIPTION
      !===========
         !---CRUNCEP Version 7 - Atmospheric Forcing Data for the Community Land Model

      !data source:
      !-------------------
         !---https://rda.ucar.edu/datasets/ds314.3/

      !References:
      !-------------------
         !---Viovy, Nicolas. (2018). CRUNCEP Version 7 -
         !   Atmospheric Forcing Data for the Community Land Model.
         !   Research Data Archive at the National Center for Atmospheric Research,
         !   Computational and Information Systems Laboratory.
         !   https://doi.org/10.5065/PZ8F-F017. Accessed 05 May 2023.

      !REVISION HISTORY
      !----------------
         !---

         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//'-'//trim(monthstr)//'.nc'
      case ('ERA5LAND') ! ERA5-Land forcing data
      !DESCRIPTION
      !===========
         !---enhanced global dataset for the land component of the fifth
         !   generation of European ReAnalysis (ERA5)

      !data source:
      !-------------------
         !---https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form

      !References:
      !-------------------
         !---Muñoz-Sabater, J., Dutra, E., Agustí-Panareda, A., Albergel, C., Arduini, G., Balsamo, G., Boussetta, S.,
         !   Choulga, M., Harrigan, S., Hersbach, H. and Martens, B., 2021. ERA5-Land:
         !   A state-of-the-art global reanalysis dataset for land applications. Earth System
         !   Science Data, 13(9), pp.4349-4383.

      !REVISION HISTORY
      !----------------
         !---2021.11.01   Zhongwang Wei @ SYSU: zip file to reduce the size of the data; remove offset and scale_factor

         metfilename = '/'//trim(fprefix(var_i))//'_'//trim(yearstr)//'_'//trim(monthstr)
         select case (var_i)
         case (1)
            metfilename = trim(metfilename) // '_2m_temperature.nc'
         case (2)
            metfilename = trim(metfilename) //'_specific_humidity.nc'
         case (3)
            metfilename = trim(metfilename) //'_surface_pressure.nc'
         case (4)
            metfilename = trim(metfilename) //'_total_precipitation_m_hr.nc'
         case (5)
            metfilename = trim(metfilename) //'_10m_u_component_of_wind.nc'
         case (6)
            metfilename = trim(metfilename) //'_10m_v_component_of_wind.nc'
         case (7)
            metfilename = trim(metfilename) //'_surface_solar_radiation_downwards_w_m2.nc'
         case (8)
            metfilename = trim(metfilename) //'_surface_thermal_radiation_downwards_w_m2.nc'
         END select
      case ('ERA5') ! ERA5 forcing data
      !DESCRIPTION
      !===========
         !---The fifth generation of European ReAnalysis (ERA5)

      !data source:
      !-------------------
         !---https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

      !References:
      !-------------------
         !---Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz‐Sabater, J.,
         !   Nicolas, J., Peubey, C., Radu, R., Schepers, D. and Simmons, A., 2020.
         !   The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730), pp.1999-2049.

      !REVISION HISTORY
      !----------------
         !---2021.11.01   Zhongwang Wei @ SYSU: zip file to reduce the size of the data; remove offset and scale_factor

         metfilename = '/'//trim(fprefix(var_i))//'_'//trim(yearstr)//'_'//trim(monthstr)
         select case (var_i)
         case (1)
            metfilename = trim(metfilename) // '_2m_temperature.nc4'
         case (2)
            metfilename = trim(metfilename) //'_q.nc4'
         case (3)
            metfilename = trim(metfilename) //'_surface_pressure.nc4'
         case (4)
            metfilename = trim(metfilename) //'_mean_total_precipitation_rate.nc4'
         case (5)
            metfilename = trim(metfilename) //'_10m_u_component_of_wind.nc4'
         case (6)
            metfilename = trim(metfilename) //'_10m_v_component_of_wind.nc4'
         case (7)
            metfilename = trim(metfilename) //'_mean_surface_downward_short_wave_radiation_flux.nc4'
         case (8)
            metfilename = trim(metfilename) //'_mean_surface_downward_long_wave_radiation_flux.nc4'
         END select
      case ('MSWX') ! MSWX forcing data
      !DESCRIPTION
      !===========
         !---Multi-Source Weather forcing data

      !data source:
      !-------------------
         !---https://www.gloh2o.org/mswx/

      !References:
      !-------------------
         !---Beck, H.E., van Dijk, A.I., Larraondo, P.R., McVicar, T.R., Pan, M., Dutra, E. and Miralles, D.G., 2022.
         !   MSWX: Global 3-hourly 0.1 bias-corrected meteorological data including near-real-time
         !   updates and forecast ensembles. Bulletin of the American Meteorological Society, 103(3), pp.E710-E732.

      !REVISION HISTORY
      !----------------
         !---2021.11.01   Zhongwang Wei @ SYSU: Regroup data into monthly

         metfilename = '/'//trim(fprefix(var_i))//'_'//trim(yearstr)//'_'//trim(monthstr)//'.nc'
      case ('WFDE5')
      !DESCRIPTION
      !===========
         !---WATCH Forcing Data methodology applied to ERA5 reanalysis data

      !data source:
      !-------------------
         !---https://doi.org/10.24381/cds.20d54e34

      !References:
      !-------------------
         !---Cucchi, M., Weedon, G.P., Amici, A., Bellouin, N., Lange, S., Müller Schmied, H.,
         !   Hersbach, H. and Buontempo, C., 2020. WFDE5: bias-adjusted ERA5 reanalysis data
         !   for impact studies. Earth System Science Data, 12(3), pp.2097-2120.

      !REVISION HISTORY
      !----------------
         !---2021.11.01   Zhongwang Wei @ SYSU: zip file to reduce the size of the data; remove offset and scale_factor

         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//trim(monthstr)//'_v2.0.nc'
      case ('CRUJRA')
      !DESCRIPTION
      !===========
         !---Collection of CRU JRA forcing datasets of gridded land surface blend
         !   of Climatic Research Unit (CRU) and Japanese reanalysis (JRA) data

      !data source:
      !-------------------
         !---https://catalogue.ceda.ac.uk/uuid/863a47a6d8414b6982e1396c69a9efe8

      !References:
      !-------------------
         !---University of East Anglia Climatic Research Unit; Harris, I.C. (2019):
         !   CRU JRA: Collection of CRU JRA forcing datasets of gridded land surface blend
         !   of Climatic Research Unit (CRU) and Japanese reanalysis (JRA) data..
         ! Centre for Environmental Data Analysis, date of citation.
         !http://catalogue.ceda.ac.uk/uuid/863a47a6d8414b6982e1396c69a9efe8

      !REVISION HISTORY
      !----------------
         !---2021.11.01   Zhongwang Wei @ SYSU: zip file to reduce the size of the data; remove offset and scale_factor

         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//'.365d.noc.nc'


      case ('WFDEI')
      !DESCRIPTION
      !===========
         !---WATCH Forcing Data methodology applied to ERA-Interim reanalysis data

      !data source:
      !-------------------
         !---https://doi.org/10.24381/cds.20d54e34

      !References:
      !-------------------
         !---Weedon, G.P., Balsamo, G., Bellouin, N., Gomes, S., Best, M.J. and Viterbo, P., 2014.
         !   The WFDEI meteorological forcing data set: WATCH Forcing Data methodology applied
         !   to ERA‐Interim reanalysis data. Water Resources Research, 50(9), pp.7505-7514.

      !REVISION HISTORY
      !----------------
         !---2021.11.01   Zhongwang Wei @ SYSU: zip file to reduce the size of the data; remove offset and scale_factor

         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//'-'//trim(monthstr)//'.nc'
      case ('JRA55')
      !DESCRIPTION
      !===========
         !---the Japanese 55-year Reanalysis

      !data source:
      !-------------------
         !---https://jra.kishou.go.jp/JRA-55/index_en.html

      !References:
      !-------------------
         !---Kobayashi, S., Y. Ota, Y. Harada, A. Ebita, M. Moriya, H. Onoda, K. Onogi,
         !   H. Kamahori, C. Kobayashi, H. Endo, K. Miyaoka, and K. Takahashi , 2015:
         !   The JRA-55 Reanalysis: General specifications and basic characteristics.
         !   J. Meteor. Soc. Japan, 93, 5-48, doi:10.2151/jmsj.2015-001.
         !---Harada, Y., H. Kamahori, C. Kobayashi, H. Endo, S. Kobayashi, Y. Ota, H. Onoda,
         !   K. Onogi, K. Miyaoka, and K. Takahashi, 2016: The JRA-55 Reanalysis:
         !   Representation of atmospheric circulation and climate variability, J. Meteor. Soc. Japan,
         !   94, 269-302, doi:10.2151/jmsj.2016-015.

      !REVISION HISTORY
      !----------------
         !---2021.11.01   Zhongwang Wei @ SYSU: zip file to reduce the size of the data; remove offset and scale_factor


         metfilename = '/'//trim(fprefix(var_i))//'_'//trim(yearstr)//'.nc'
      case ('GDAS')
      !DESCRIPTION
      !===========
         !--- Forcing Data From Global Data Assimilation System

      !data source:
      !-------------------
         !--https://disc.sci.gsfc.nasa.gov/datasets/GLDAS_NOAH025_3H_V2.1/summary

      !References:
      !-------------------
         !---Beaudoing, H. and M. Rodell, NASA/GSFC/HSL (2020), GLDAS Noah Land Surface Model L4 3 hourly 0.25 x 0.25
         !   degree V2.1, Greenbelt, Maryland, USA, Goddard Earth Sciences Data and Information Services
         !   Center (GES DISC), Accessed: [Data Access Date], 10.5067/E7TYRXPJKWOQ
         !---Rodell, M., P.R. Houser, U. Jambor, J. Gottschalck, K. Mitchell, C. Meng, K. Arsenault, B. Cosgrove,
         !   J. Radakovich, M. Bosilovich, J.K. Entin, J.P. Walker, D. Lohmann, and D. Toll, 2004:
         !   The Global Land Data Assimilation System, Bull. Amer. Meteor. Soc., 85, 381-394,
         !   doi:10.1175/BAMS-85-3-381

      !REVISION HISTORY
      !----------------
         !---2021.11.01   Zhongwang Wei @ SYSU: merge the data into monthly file; zip file to reduce the size of the data; remove offset and scale_factor

         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//trim(monthstr)//'.nc4'
      case ('CLDAS')

      !DESCRIPTION
      !===========
         !---The Real-Time Product Dataset Of The China Meteorological Administration
         !   Land Data Assimilation System

      !data source:
      !-------------------
         !--CMA, not pulicly available

      !References:
      !-------------------
         !---Xia, Y.L.; Hao, Z.C.; Shi, C.X.; Li, Y.H.; Meng, J.; Xu, T.R.; Wu, X.Y.; Zhang, B.Q.
         !    Regional and global land data assimilation systems: Innovations, challenges, and
         !    prospects. J. Meteorol. Res. 2019, 33, 159–189.

      !REVISION HISTORY
      !----------------
         !---2021.11.01   Zhongwang Wei @ SYSU: gap filling for the missing data; zip file to reduce the size of the data; remove offset and scale_factor

         metfilename = '/'//trim(fprefix(var_i))//'-'//trim(yearstr)//trim(monthstr)//'.nc'
      case ('CMFD')
      !DESCRIPTION
      !===========
         !--- The China Meteorological Forcing Dataset

      !data source:
      !-------------------
         !--https://data.tpdc.ac.cn/en/data/8028b944-daaa-4511-8769-965612652c49/

      !References:
      !-------------------
         !---He, J., Yang, K., Tang, W., Lu, H., Qin, J., Chen, Y. and Li, X., 2020.
         !   The first high-resolution meteorological forcing dataset for land process
         !   studies over China. Scientific data, 7(1), p.25.

         !REVISION HISTORY
      !----------------
         !---

         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//trim(monthstr)//'.nc4'
      case ('CMIP6')
      !DESCRIPTION
      !===========
         !---the Climate Model Intercomparison Project Phase 6 (CMIP6) forcing data sets

      !data source:
      !-------------------
         !---https://esgf-node.llnl.gov/projects/cmip6/

      !References:
      !-------------------
         !---

         !REVISION HISTORY
      !----------------
         !---2021.11.01   Zhongwang Wei @ SYSU: regroup the data into monthly file;
         !   zip file to reduce the size of the data; remove offset and scale_factor


         metfilename = '/'//trim(fprefix(var_i))//'_'//trim(yearstr)//'.nc'
      case ('TPMFD')
      !DESCRIPTION
      !===========
         !---A high-resolution near-surface meteorological forcing dataset for the Third Pole region

      !data source:
      !-------------------
         !---https://data.tpdc.ac.cn/zh-hans/data/44a449ce-e660-44c3-bbf2-31ef7d716ec7

      !References:
      !-------------------
         !---Yang, K., Jiang, Y., Tang, W., He, J., Shao, C., Zhou, X., Lu, H.,
         !   Chen, Y., Li, X., Shi, J. (2023). A high-resolution near-surface
         !   meteorological forcing dataset for the Third Pole region （TPMFD, 1979-2020）.
         !   National Tibetan Plateau/Third Pole Environment Data Center,
         !   https://doi.org/10.11888/Atmos.tpdc.300398. https://cstr.cn/18406.11.Atmos.tpdc.300398.

      !REVISION HISTORY
      !----------------
         !---2021.11.01   Zhongwang Wei @ SYSU: regroup the data into monthly file;
         !   zip file to reduce the size of the data; remove offset and scale_factor

         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//trim(monthstr)//'.nc'
      case ('POINT')
         metfilename = '/'//trim(fprefix(1))
      end select
   if (DEF_USE_CBL_HEIGHT) then
      select case (var_i)
      case (9)
         metfilename = '/'//trim(fprefix(9))//'_'//trim(yearstr)//'_'//trim(monthstr)//'_boundary_layer_height.nc4'
      END select
   endif
   END FUNCTION metfilename

 ! preprocess for forcing data [not applicable yet for PRINCETON]
 ! ------------------------------------------------------------
   SUBROUTINE metpreprocess(grid, forcn)

      use MOD_Const_Physical
      use MOD_Namelist
      use MOD_SPMD_Task
      use MOD_Block
      use MOD_Grid
      use MOD_DataType
      USE MOD_Qsadv
      implicit none
      type(grid_type), intent(in) :: grid
      type(block_data_real8_2d), intent(inout) :: forcn(:)

      integer  :: iblkme, ib, jb, i, j
      real(r8) :: es, esdT, qsat_tmp, dqsat_tmpdT, e, ea

      !----------------------------------------------------------------------------
      ! use polynomials to calculate saturation vapor pressure and derivative with
      ! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
      ! required to convert relative humidity to specific humidity
      !----------------------------------------------------------------------------
      if (trim(DEF_forcing%dataset) == 'POINT') then
#ifdef SinglePoint
         call qsadv(forcn(1)%blk(gblock%xblkme(1),gblock%yblkme(1))%val(1,1), &
                    forcn(3)%blk(gblock%xblkme(1),gblock%yblkme(1))%val(1,1), &
                    es,esdT,qsat_tmp,dqsat_tmpdT)
         if (qsat_tmp < forcn(2)%blk(gblock%xblkme(1),gblock%yblkme(1))%val(1,1)) THEN
            forcn(2)%blk(gblock%xblkme(1),gblock%yblkme(1))%val(1,1) = qsat_tmp
         ENDIF
#endif
      else
         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)

            do j = 1, grid%ycnt(jb)
               do i = 1, grid%xcnt(ib)

                  select case (trim(DEF_forcing%dataset))

                  case ('PRINCETON')

                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('GSWP2')

                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('GSWP3')
                     if (forcn(1)%blk(ib,jb)%val(i,j)<212.0) forcn(1)%blk(ib,jb)%val(i,j) = 212.0
                     if (forcn(4)%blk(ib,jb)%val(i,j)<0.0) forcn(4)%blk(ib,jb)%val(i,j) = 0.0
                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('QIAN')

                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                     e  = forcn(3)%blk(ib,jb)%val(i,j) * forcn(2)%blk(ib,jb)%val(i,j) &
                        / (0.622_R8 + 0.378_R8 * forcn(2)%blk(ib,jb)%val(i,j))
                     ea = 0.70_R8 + 5.95e-05_R8 * 0.01_R8 * e * exp(1500.0_R8/forcn(1)%blk(ib,jb)%val(i,j))
                     forcn(8)%blk(ib,jb)%val(i,j) = ea * stefnc * forcn(1)%blk(ib,jb)%val(i,j)**4

                  case ('CRUNCEPV4')

                     if (forcn(1)%blk(ib,jb)%val(i,j) < 212.0) forcn(1)%blk(ib,jb)%val(i,j) = 212.0
                     if (forcn(4)%blk(ib,jb)%val(i,j) < 0.0)   forcn(4)%blk(ib,jb)%val(i,j) = 0.0
                     if (forcn(7)%blk(ib,jb)%val(i,j) < 0.0)   forcn(7)%blk(ib,jb)%val(i,j) = 0.0
                     ! 12th grade of Typhoon 32.7-36.9 m/s
                     if (abs(forcn(5)%blk(ib,jb)%val(i,j)) > 40.0) forcn(5)%blk(ib,jb)%val(i,j) = &
                        40.0*forcn(5)%blk(ib,jb)%val(i,j)/abs(forcn(5)%blk(ib,jb)%val(i,j))
                     if (abs(forcn(6)%blk(ib,jb)%val(i,j)) > 40.0) forcn(6)%blk(ib,jb)%val(i,j) = &
                        40.0*forcn(6)%blk(ib,jb)%val(i,j)/abs(forcn(6)%blk(ib,jb)%val(i,j))
                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('CRUNCEPV7')

                     if (forcn(1)%blk(ib,jb)%val(i,j) < 212.0) forcn(1)%blk(ib,jb)%val(i,j) = 212.0
                     if (forcn(4)%blk(ib,jb)%val(i,j) < 0.0)   forcn(4)%blk(ib,jb)%val(i,j) = 0.0
                     if (forcn(7)%blk(ib,jb)%val(i,j) < 0.0)   forcn(7)%blk(ib,jb)%val(i,j) = 0.0
                     ! 12th grade of Typhoon 32.7-36.9 m/s
                     ! NOTE by Wenzong: This is a problem when running a GNU-compiled program, because there is
                     ! no data of forcn(5), temporarily comment the code below
                     ! if (abs(forcn(5)%blk(ib,jb)%val(i,j)) > 40.0) forcn(5)%blk(ib,jb)%val(i,j) = &
                     !    40.0*forcn(5)%blk(ib,jb)%val(i,j)/abs(forcn(5)%blk(ib,jb)%val(i,j))
                     if (abs(forcn(6)%blk(ib,jb)%val(i,j)) > 40.0) forcn(6)%blk(ib,jb)%val(i,j) = &
                        40.0*forcn(6)%blk(ib,jb)%val(i,j)/abs(forcn(6)%blk(ib,jb)%val(i,j))
                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('ERA5LAND')

                     forcn(4)%blk(ib,jb)%val(i,j)=forcn(4)%blk(ib,jb)%val(i,j) * 1000./3600.
                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('ERA5')

                     if (forcn(4)%blk(ib,jb)%val(i,j) < 0.0)   forcn(4)%blk(ib,jb)%val(i,j) = 0.0
                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif
                     if (forcn(4)%blk(ib,jb)%val(i,j) < 0.0)   forcn(4)%blk(ib,jb)%val(i,j) = 0.0

                  case ('MSWX')

                     forcn(1)%blk(ib,jb)%val(i,j)=forcn(1)%blk(ib,jb)%val(i,j)+273.15
                     forcn(4)%blk(ib,jb)%val(i,j)=forcn(4)%blk(ib,jb)%val(i,j)/10800.
                     if (forcn(4)%blk(ib,jb)%val(i,j)>1000.0) forcn(4)%blk(ib,jb)%val(i,j) = 0.0
                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('WFDE5')

                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('WFDEI')

                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('CLDAS') ! CLDAS forcing

                     forcn(4)%blk(ib,jb)%val(i,j)=forcn(4)%blk(ib,jb)%val(i,j)/3600.
                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('CMFD') ! CMFD forcing

                     forcn(4)%blk(ib,jb)%val(i,j)=forcn(4)%blk(ib,jb)%val(i,j)/3600.
                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('CRUJRA') ! CRUJRA forcing
                     forcn(4)%blk(ib,jb)%val(i,j)=forcn(4)%blk(ib,jb)%val(i,j)/21600.
                     forcn(7)%blk(ib,jb)%val(i,j)=forcn(7)%blk(ib,jb)%val(i,j)/21600.
                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('GDAS') ! GDAS forcing

                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('JRA55') ! JRA55 forcing

                     forcn(4)%blk(ib,jb)%val(i,j)=forcn(4)%blk(ib,jb)%val(i,j)/86400.
                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('TPMFD') ! TPMFD forcing

                     forcn(4)%blk(ib,jb)%val(i,j)=forcn(4)%blk(ib,jb)%val(i,j)/3600.! convert to mm/s
                     forcn(3)%blk(ib,jb)%val(i,j)=forcn(3)%blk(ib,jb)%val(i,j)*100. ! convert to pa

                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif

                  case ('CMIP6') ! CMIP6 forcing

                     if (forcn(4)%blk(ib,jb)%val(i,j) < 0.0)   forcn(4)%blk(ib,jb)%val(i,j) = 0.0
                     call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     endif
                     if (forcn(4)%blk(ib,jb)%val(i,j) < 0.0)   forcn(4)%blk(ib,jb)%val(i,j) = 0.0
                  end select

               end do
            end do
         end do
      end if

   END SUBROUTINE metpreprocess

END MODULE MOD_UserSpecifiedForcing
! ----------- EOP ---------------
