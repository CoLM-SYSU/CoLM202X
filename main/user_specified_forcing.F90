
#include <define.h>

MODULE user_specified_forcing

! ------------------------------------------------------------
! MODULE NANE: 
!     PRIMCETON GSWP2 QIAN POINT
!
! PURPOSE :
!     Read PRINCETON/GSWP2/QIAN/CRUNCEP/POINT forcing data     
!
!     PLEASE modify the following codes when specified forcing used
!     metpreprocess modified by siguang & weinan for forc_q calibration
! ------------------------------------------------------------

#if (defined USE_PRINCETON_DATA)

   use precision
   implicit none

 ! parameter setting
 ! ------------------------------------------------------------
   integer, parameter :: NVAR    = 8              ! variable number of forcing data
   integer, parameter :: nlats   = 180            ! number of latitudes
   integer, parameter :: nlons   = 360            ! number of longitudes
   integer, parameter :: startyr = 1948           ! start year of forcing data
   integer, parameter :: startmo = 1              ! start month of forcing data
   integer, parameter :: endyr   = 2006           ! end year of forcing data
   integer, parameter :: endmo   = 12             ! end month of forcing data
   integer, parameter :: dtime(NVAR)  = 10800     ! temporal resolution
   integer, parameter :: offset(NVAR) = 0         ! time offset (seconds)
   integer, parameter :: nlands  = 1              ! land grid number in 1d

   logical, parameter :: leapyear = .true.        ! leapyear calendar
   logical, parameter :: data2d   = .true.        ! data in 2 dimension (lon, lat)
   logical, parameter :: hightdim = .true.        ! have "z" dimension
   logical, parameter :: dim2d    = .false.       ! lat/lon value in 2 dimension (lon, lat)
   logical, parameter :: latrev   = .true.        ! need to reverse latitudes
   logical, parameter :: lonadj   = .true.        ! need to adjust longitude, 0~360 -> -180~180

   character(len=256), parameter :: latname = 'latitude'  ! dimension name of latitude
   character(len=256), parameter :: lonname = 'longitude' ! dimension name of longitude

 ! file grouped by year/month
   character(len=256), parameter :: groupby = 'year'
   
 ! prefix of forcing data file
   character(len=256), parameter :: fprefix(NVAR) = [character(len=256) :: &
      'tas/tas_3hourly_', 'shum/shum_3hourly_', 'pres/pres_3hourly_', &
      'prcp/prcp_3hourly_', 'NULL', 'wind/wind_3hourly_', &
      'dswrf/dswrf_3hourly_', 'dlwrf/dlwrf_3hourly_']

 ! variable name of forcing data file
   character(len=256), parameter :: vname(NVAR) = [character(len=256) :: &
      'tas', 'shum', 'pres', 'prcp', 'NULL', 'wind', 'dswrf', 'dlwrf']
 
 ! interpolation method
   character(len=256), parameter :: tintalgo(NVAR) = [character(len=256) :: &
      'linear', 'linear', 'linear', 'nearest', 'NULL', 'linear', 'coszen', 'linear']

   INTERFACE getfilename
      MODULE procedure getfilename
   END INTERFACE

   INTERFACE metpreprocess
      MODULE procedure metpreprocess
   END INTERFACE

   public metpreprocess
 
CONTAINS

   FUNCTION getfilename(year, month, var_i)

      implicit none
      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: var_i
      character(len=256)  :: getfilename
      character(len=256)  :: yearstr

      write(yearstr, '(I4.4)') year
      getfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//'-'//trim(yearstr)//'.nc'
      return
   END FUNCTION getfilename

 ! preprocess for forcing data [not applicable yet for PRINCETON]
 ! ------------------------------------------------------------
   SUBROUTINE metpreprocess(forcn)

      implicit none
      real(r8), intent(inout) :: forcn(:,:,:)

      integer  :: i, j
      real(r8) :: es, esdT, qsat_tmp, dqsat_tmpdT

!----------------------------------------------------------------------------
! use polynomials to calculate saturation vapor pressure and derivative with
! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
! required to convert relative humidity to specific humidity
!----------------------------------------------------------------------------
      do i = 1, nlats 
         do j = 1, nlons
            call qsadv(forcn(j,i,1),forcn(j,i,3),es,esdT,qsat_tmp,dqsat_tmpdT)
            if (qsat_tmp < forcn(j,i,2)) then
               forcn(j,i,2) = qsat_tmp
            endif
         end do
      end do

   END SUBROUTINE metpreprocess

#endif
 


#if(defined USE_GSWP2_DATA)

   use precision
   implicit none

 ! parameter setting
 ! ------------------------------------------------------------
   integer, parameter :: NVAR    = 8              ! variable number of forcing data
   integer, parameter :: nlats   = 150            ! number of latitudes
   integer, parameter :: nlons   = 360            ! number of longitudes
   integer, parameter :: startyr = 1982           ! start year of forcing data
   integer, parameter :: startmo = 7              ! start month of forcing data
   integer, parameter :: endyr   = 1995           ! end year of forcing data
   integer, parameter :: endmo   = 12             ! end month of forcing data
   integer, parameter :: dtime(NVAR)  = 10800     ! temporal resolution
   integer, parameter :: offset(NVAR) = 10800     ! time offset (seconds)
   integer, parameter :: nlands  = 15238          ! land grid number in 1d

   logical, parameter :: leapyear = .true.        ! leapyear calendar
   logical, parameter :: data2d   = .false.       ! data in 2 dimension (lon, lat)
   logical, parameter :: hightdim = .false.       ! have "z" dimension
   logical, parameter :: dim2d    = .true.        ! lat/lon value in 2 dimension (lon, lat)
   logical, parameter :: latrev   = .false.       ! need to reverse latitudes
   logical, parameter :: lonadj   = .false.       ! need to adjust longitude, 0~360 -> -180~180

   character(len=256), parameter :: latname = 'nav_lat'  ! dimension name of latitude
   character(len=256), parameter :: lonname = 'nav_lon'  ! dimension name of longitude

 ! file grouped by year/month
   character(len=256), parameter :: groupby = 'month'
   
 ! prefix of forcing data file
   character(len=256), parameter :: fprefix(NVAR) = [character(len=256) :: &
      'Tair_cru/Tair_cru', 'Qair_cru/Qair_cru', 'PSurf_ecor/PSurf_ecor', &
      'Rainf_gswp/Rainf_gswp', 'Rainf_C_gswp/Rainf_C_gswp', 'Wind_ncep/Wind_ncep', &
      'SWdown_srb/SWdown_srb', 'LWdown_srb/LWdown_srb']

 ! variable name of forcing data file
   character(len=256), parameter :: vname(NVAR) = [character(len=256) :: &
      'Tair', 'Qair', 'PSurf', 'Rainf', 'Rainf_C', 'Wind', 'SWdown', 'LWdown']
 
 ! interpolation method
   character(len=256), parameter :: tintalgo(NVAR) = [character(len=256) :: &
      'linear', 'linear', 'linear', 'nearest', 'nearest', 'linear', 'coszen', 'linear']
   
   INTERFACE getfilename
      MODULE procedure getfilename
   END INTERFACE

   INTERFACE metpreprocess
      MODULE procedure metpreprocess
   END INTERFACE
 
   public metpreprocess

CONTAINS

   FUNCTION getfilename(year, month, var_i)

      implicit none
      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: var_i
      character(len=256)  :: getfilename
      character(len=256)  :: yearstr
      character(len=256)  :: monthstr
      
      write(yearstr, '(I4.4)') year
      write(monthstr, '(I2.2)') month
      getfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//trim(monthstr)//'.nc'
      return
   END FUNCTION getfilename

 ! preprocess for forcing data [not applicable yet for GSWP2]
 ! ------------------------------------------------------------
   SUBROUTINE metpreprocess(forcn)

      implicit none
      real(r8), intent(inout) :: forcn(:,:,:)

      integer  :: i, j
      real(r8) :: es, esdT, qsat_tmp, dqsat_tmpdT

!----------------------------------------------------------------------------
! use polynomials to calculate saturation vapor pressure and derivative with
! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
! required to convert relative humidity to specific humidity
!----------------------------------------------------------------------------
      do i = 1, nlats 
         do j = 1, nlons
            call qsadv(forcn(j,i,1),forcn(j,i,3),es,esdT,qsat_tmp,dqsat_tmpdT)
            if (qsat_tmp < forcn(j,i,2)) then
               forcn(j,i,2) = qsat_tmp
            endif
         end do
      end do

   END SUBROUTINE metpreprocess

#endif


 
#if(defined USE_QIAN_DATA)

   use precision
   use PhysicalConstants
   implicit none

 ! parameter setting
 ! ------------------------------------------------------------
   integer, parameter :: NVAR    = 8              ! variable number of forcing data
   integer, parameter :: nlats   = 94             ! number of latitudes
   integer, parameter :: nlons   = 192            ! number of longitudes
   integer, parameter :: startyr = 1972           ! start year of forcing data       <MARK #2>
   integer, parameter :: startmo = 1              ! start month of forcing data
   integer, parameter :: endyr   = 2004           ! end year of forcing data
   integer, parameter :: endmo   = 12             ! end month of forcing data
   integer, parameter :: dtime(NVAR)  = (/ &      ! temporal resolution
      10800, 10800, 10800, 21600, 0, 10800, 21600, 0/)
   integer, parameter :: offset(NVAR) = (/ &      ! time offset (seconds)
      5400, 5400, 5400, 10800, 0, 5400, 0, 0/)    ! ..
   integer, parameter :: nlands  = 1              ! land grid number in 1d

   logical, parameter :: leapyear = .false.       ! leapyear calendar
   logical, parameter :: data2d   = .true.        ! data in 2 dimension (lon, lat)
   logical, parameter :: hightdim = .false.       ! have "z" dimension
   logical, parameter :: dim2d    = .true.        ! lat/lon value in 2 dimension (lon, lat)
   logical, parameter :: latrev   = .true.        ! need to reverse latitudes
   logical, parameter :: lonadj   = .true.        ! need to adjust longitude, 0~360 -> -180~180

   character(len=256), parameter :: latname = 'LATIXY'  ! dimension name of latitude
   character(len=256), parameter :: lonname = 'LONGXY'  ! dimension name of longitude

 ! file grouped by year/month
   character(len=256), parameter :: groupby = 'month'
   
 ! prefix of forcing data file
   character(len=256), parameter :: fprefix(NVAR) = [character(len=256) :: &
      'TmpPrsHumWnd3Hrly/clmforc.Qian.c2006.T62.TPQW.', &
      'TmpPrsHumWnd3Hrly/clmforc.Qian.c2006.T62.TPQW.', &
      'TmpPrsHumWnd3Hrly/clmforc.Qian.c2006.T62.TPQW.', &
      'Precip6Hrly/clmforc.Qian.c2006.T62.Prec.', &
      'NULL', &
      'TmpPrsHumWnd3Hrly/clmforc.Qian.c2006.T62.TPQW.', &
      'Solar6Hrly/clmforc.Qian.c2006.T62.Solr.', &
      'NULL']

 ! variable name of forcing data file
   character(len=256), parameter :: vname(NVAR) = [character(len=256) :: &
      'TBOT', 'QBOT', 'PSRF', 'PRECTmms', 'NULL', 'WIND', 'FSDS', 'NULL']
 
 ! interpolation method
   character(len=256), parameter :: tintalgo(NVAR) = [character(len=256) :: &
      'linear', 'linear', 'linear', 'nearest', 'NULL', 'linear', 'coszen', 'NULL']
   
   INTERFACE getfilename
      MODULE procedure getfilename
   END INTERFACE

   INTERFACE metpreprocess
      MODULE procedure metpreprocess
   END INTERFACE
 
   public metpreprocess

CONTAINS

   FUNCTION getfilename(year, month, var_i)

      implicit none
      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: var_i
      character(len=256)  :: getfilename
      character(len=256)  :: yearstr
      character(len=256)  :: monthstr
      
      write(yearstr, '(I4.4)') year
      write(monthstr, '(I2.2)') month
      getfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//'-'//trim(monthstr)//'.nc'
      return
   END FUNCTION getfilename

 ! preprocess for forcing data: calculate LW
 ! ------------------------------------------------------------
   SUBROUTINE metpreprocess(forcn)

      implicit none
      real(r8), intent(inout) :: forcn(:,:,:)

      real(r8) :: e, ea
      integer  :: i, j
      real(r8) :: es, esdT, qsat_tmp, dqsat_tmpdT

!----------------------------------------------------------------------------
! use polynomials to calculate saturation vapor pressure and derivative with
! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
! required to convert relative humidity to specific humidity
!----------------------------------------------------------------------------
      do i = 1, nlats 
         do j = 1, nlons
            call qsadv(forcn(j,i,1),forcn(j,i,3),es,esdT,qsat_tmp,dqsat_tmpdT)
            if (qsat_tmp < forcn(j,i,2)) then
               forcn(j,i,2) = qsat_tmp
            endif
         end do
      end do
      
      do i = 1, nlats
         do j = 1, nlons
            e  = forcn(j,i,3) * forcn(j,i,2) / (0.622_R8 + 0.378_R8 * forcn(j,i,2))
            ea = 0.70_R8 + 5.95e-05_R8 * 0.01_R8 * e * exp(1500.0_R8/forcn(j,i,1))
            forcn(j,i,8) = ea * stefnc * forcn(j,i,1)**4
         end do
      end do

   END SUBROUTINE metpreprocess

#endif



#if(defined USE_CRUNCEP_DATA)

   use precision
   use PhysicalConstants
   implicit none

 ! parameter setting
 ! ------------------------------------------------------------
   integer, parameter :: NVAR    = 8              ! variable number of forcing data
   integer, parameter :: nlats   = 360            ! number of latitudes
   integer, parameter :: nlons   = 720            ! number of longitudes
   integer, parameter :: startyr = 2000           ! start year of forcing data        <MARK #1>
   integer, parameter :: startmo = 1              ! start month of forcing data
   integer, parameter :: endyr   = 2003           ! end year of forcing data
   integer, parameter :: endmo   = 12             ! end month of forcing data
   integer, parameter :: dtime(NVAR)  = (/ &      ! temporal resolution
      21600, 21600, 21600, 21600, 0, 21600, 21600, 21600/)
   integer, parameter :: offset(NVAR) = (/ &      ! time offset (seconds)
      10800, 10800, 10800, 10800, 0, 10800, 0, 10800/)    ! ..
   integer, parameter :: nlands  = 1              ! land grid number in 1d

   logical, parameter :: leapyear = .false.       ! leapyear calendar
   logical, parameter :: data2d   = .true.        ! data in 2 dimension (lon, lat)
   logical, parameter :: hightdim = .false.       ! have "z" dimension
   logical, parameter :: dim2d    = .true.        ! lat/lon value in 2 dimension (lon, lat)
   logical, parameter :: latrev   = .true.        ! need to reverse latitudes
   logical, parameter :: lonadj   = .true.        ! need to adjust longitude, 0~360 -> -180~180

   character(len=256), parameter :: latname = 'LATIXY'  ! dimension name of latitude
   character(len=256), parameter :: lonname = 'LONGXY'  ! dimension name of longitude

 ! file grouped by year/month
   character(len=256), parameter :: groupby = 'month'
   
 ! prefix of forcing data file
   character(len=256), parameter :: fprefix(NVAR) = [character(len=256) :: &
      'TPHWL6Hrly/clmforc.cruncep.V4.c2011.0.5d.TPQWL.', &
      'TPHWL6Hrly/clmforc.cruncep.V4.c2011.0.5d.TPQWL.', &
      'TPHWL6Hrly/clmforc.cruncep.V4.c2011.0.5d.TPQWL.', &
      'Precip6Hrly/clmforc.cruncep.V4.c2011.0.5d.Prec.', &
      'NULL', &
      'TPHWL6Hrly/clmforc.cruncep.V4.c2011.0.5d.TPQWL.', &
      'Solar6Hrly/clmforc.cruncep.V4.c2011.0.5d.Solr.', &
      'TPHWL6Hrly/clmforc.cruncep.V4.c2011.0.5d.TPQWL.']

 ! variable name of forcing data file
   character(len=256), parameter :: vname(NVAR) = [character(len=256) :: &
      'TBOT', 'QBOT', 'PSRF', 'PRECTmms', 'NULL', 'WIND', 'FSDS', 'FLDS']
 
 ! interpolation method
   character(len=256), parameter :: tintalgo(NVAR) = [character(len=256) :: &
      'linear', 'linear', 'linear', 'nearest', 'NULL', 'linear', 'coszen', 'linear']
   
   INTERFACE getfilename
      MODULE procedure getfilename
   END INTERFACE

   INTERFACE metpreprocess
      MODULE procedure metpreprocess
   END INTERFACE
 
   public metpreprocess

CONTAINS

   FUNCTION getfilename(year, month, var_i)

      implicit none
      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: var_i
      character(len=256)  :: getfilename
      character(len=256)  :: yearstr
      character(len=256)  :: monthstr
      
      write(yearstr, '(I4.4)') year
      write(monthstr, '(I2.2)') month
      getfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//'-'//trim(monthstr)//'.nc'
      return
   END FUNCTION getfilename

 ! preprocess for forcing data: calculate LW
 ! ------------------------------------------------------------
   SUBROUTINE metpreprocess(forcn)

      implicit none
      real(r8), intent(inout) :: forcn(:,:,:)

      real(r8) :: e, ea
      integer  :: i, j
      real(r8) :: es, esdT, qsat_tmp, dqsat_tmpdT

!----------------------------------------------------------------------------
! use polynomials to calculate saturation vapor pressure and derivative with
! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
! required to convert relative humidity to specific humidity
!----------------------------------------------------------------------------
      do i = 1, nlats 
         do j = 1, nlons
            call qsadv(forcn(j,i,1),forcn(j,i,3),es,esdT,qsat_tmp,dqsat_tmpdT)
            if (qsat_tmp < forcn(j,i,2)) then
               forcn(j,i,2) = qsat_tmp
            endif
         end do
      end do
      
   END SUBROUTINE metpreprocess

#endif



#if(defined USE_POINT_DATA)

   use precision
   implicit none
 
 ! parameter setting
 ! ------------------------------------------------------------
   integer, parameter :: NVAR    = 8              ! variable number of forcing data
 ! not applied for POINT
   integer, parameter :: nlats   = 2              ! number of latitudes
   integer, parameter :: nlons   = 2              ! number of longitudes
   integer, parameter :: startyr = 1948           ! start year of forcing data
   integer, parameter :: startmo = 7              ! start month of forcing data
   integer, parameter :: endyr   = 2006           ! end year of forcing data
   integer, parameter :: endmo   = 12             ! end month of forcing data
   integer, parameter :: dtime(NVAR)  = 1800      ! temporal resolution
   integer, parameter :: offset(NVAR) = 0         ! time offset (seconds)
   integer, parameter :: nlands  = 1              ! land grid number in 1d

 ! not applied for POINT
   logical, parameter :: leapyear = .true.        ! leapyear calendar
   logical, parameter :: data2d   = .true.        ! data in 2 dimension (lon, lat)
   logical, parameter :: hightdim = .false.       ! have "z" dimension
   logical, parameter :: dim2d    = .false.       ! lat/lon value in 2 dimension (lon, lat)
   logical, parameter :: latrev   = .true.        ! need to reverse latitudes
   logical, parameter :: lonadj   = .true.        ! need to adjust longitude, 0~360 -> -180~180

 ! not applied for POINT
   character(len=256), parameter :: latname = 'NULL' ! dimension name of latitude
   character(len=256), parameter :: lonname = 'NULL' ! dimension name of longitude

 ! not applied for POINT
   character(len=256), parameter :: groupby = 'NULL'
   
 ! prefix of forcing data file
   character(len=256), parameter :: fprefix(NVAR) = 'VAL.DAT.CTRL.INT'

 ! not applied for POINT
   character(len=256), parameter :: vname(NVAR) = 'NULL'
 
 ! not applied for POINT
   character(len=256), parameter :: tintalgo(NVAR) = 'NULL'

   INTERFACE getfilename
      MODULE procedure getfilename
   END INTERFACE

   INTERFACE metpreprocess
      MODULE procedure metpreprocess
   END INTERFACE
 
   public metpreprocess

CONTAINS

   FUNCTION getfilename(year, month, var_i)

      implicit none
      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: var_i
      character(len=256)  :: getfilename

      getfilename = '/'//trim(fprefix(1))
      return
   END FUNCTION getfilename

 ! preprocess for forcing data [not applicable yet for POINT]
 ! ------------------------------------------------------------
   SUBROUTINE metpreprocess(forcn)

      implicit none
      real(r8), intent(inout) :: forcn(:,:,:)

      real(r8) :: es, esdT, qsat_tmp, dqsat_tmpdT

!----------------------------------------------------------------------------
! use polynomials to calculate saturation vapor pressure and derivative with
! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
! required to convert relative humidity to specific humidity
!----------------------------------------------------------------------------

      call qsadv(forcn(1,1,1),forcn(1,1,3),es,esdT,qsat_tmp,dqsat_tmpdT)
      if (qsat_tmp < forcn(1,1,2)) forcn(1,1,2) = qsat_tmp


   END SUBROUTINE metpreprocess

#endif


END MODULE user_specified_forcing
! ----------- EOP ---------------
