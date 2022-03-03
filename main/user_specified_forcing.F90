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

   use precision
   implicit none

   character(len=256) :: dataset

   logical  :: solarin_all_band   
   real(r8) :: HEIGHT_V            
   real(r8) :: HEIGHT_T           
   real(r8) :: HEIGHT_Q           

   integer  :: NVAR      ! variable number of forcing data
   integer  :: startyr   ! start year of forcing data        <MARK #1>
   integer  :: startmo   ! start month of forcing data
   integer  :: endyr     ! end year of forcing data
   integer  :: endmo     ! end month of forcing data
   integer, allocatable :: dtime(:)          
   integer, allocatable :: offset(:)        
   
   integer :: nlands     ! land grid number in 1d

   logical :: leapyear   ! leapyear calendar
   logical :: data2d     ! data in 2 dimension (lon, lat)
   logical :: hightdim   ! have "z" dimension
   logical :: dim2d      ! lat/lon value in 2 dimension (lon, lat)

   character(len=256) :: latname                   ! dimension name of latitude
   character(len=256) :: lonname                   ! dimension name of longitude

   character(len=256) :: groupby                   ! file grouped by year/month

   character(len=256), allocatable :: fprefix(:) 
   character(len=256), allocatable :: vname(:) 
   character(len=256), allocatable :: tintalgo(:) 

   ! ----- public subroutines -----
   public :: init_user_specified_forcing
   public :: metfilename
   public :: metpreprocess

CONTAINS

   ! ----------------
   subroutine init_user_specified_forcing 

      use mod_namelist
      implicit none

      ! Local variables
      integer :: ivar

      NVAR = DEF_forcing%NVAR

      allocate (dtime  (NVAR))
      allocate (offset (NVAR))

      allocate (fprefix  (NVAR))
      allocate (vname    (NVAR))
      allocate (tintalgo (NVAR))
   
     solarin_all_band = DEF_forcing%solarin_all_band   
     HEIGHT_V         = DEF_forcing%HEIGHT_V            
     HEIGHT_T         = DEF_forcing%HEIGHT_T           
     HEIGHT_Q         = DEF_forcing%HEIGHT_Q           
                     
     startyr          = DEF_forcing%startyr  
     startmo          = DEF_forcing%startmo   
     endyr            = DEF_forcing%endyr     
     endmo            = DEF_forcing%endmo     
     dtime(:)         = DEF_forcing%dtime(:)          
     offset(:)        = DEF_forcing%offset(:)        
                     
     nlands           = DEF_forcing%nlands     
                     
     leapyear         = DEF_forcing%leapyear   
     data2d           = DEF_forcing%data2d     
     hightdim         = DEF_forcing%hightdim   
     dim2d            = DEF_forcing%dim2d      
                     
     latname          = DEF_forcing%latname                   
     lonname          = DEF_forcing%lonname                  
                     
     groupby          = DEF_forcing%groupby                   

     do ivar = 1, NVAR
        fprefix (ivar) = DEF_forcing%fprefix(ivar) 
        vname   (ivar) = DEF_forcing%vname(ivar) 
        tintalgo(ivar) = DEF_forcing%tintalgo(ivar) 
     end do

   end subroutine init_user_specified_forcing 

   ! ----------------
   FUNCTION metfilename(year, month, day, var_i)

      use mod_namelist
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
      case ('PRINCETON')  
         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//'-'//trim(yearstr)//'.nc'
      case ('GSWP2')
         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//trim(monthstr)//'.nc'
      case ('QIAN')
         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//trim(monthstr)//'.nc'
      case ('CRUNCEP')
         metfilename = '/'//trim(fprefix(var_i))//trim(yearstr)//'-'//trim(monthstr)//'.nc'
      case ('ERA5LAND')
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
      case ('POINT')
         metfilename = '/'//trim(fprefix(1))
      end select

   END FUNCTION metfilename

 ! preprocess for forcing data [not applicable yet for PRINCETON]
 ! ------------------------------------------------------------
   SUBROUTINE metpreprocess(grid, forcn)

      use PhysicalConstants
      use mod_namelist
      use spmd_task
      use mod_block
      use mod_grid
      use mod_data_type
      implicit none
      type(grid_type), intent(in) :: grid
      type(block_data_real8_2d), intent(inout) :: forcn(:)

      integer  :: ib, jb, i, j
      real(r8) :: es, esdT, qsat_tmp, dqsat_tmpdT, e, ea

!----------------------------------------------------------------------------
! use polynomials to calculate saturation vapor pressure and derivative with
! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
! required to convert relative humidity to specific humidity
!----------------------------------------------------------------------------
      if (trim(DEF_forcing%dataset) == 'POINT') then

         call qsadv(forcn(1)%blk(1,1)%val(1,1),forcn(3)%blk(1,1)%val(1,1), &
            es,esdT,qsat_tmp,dqsat_tmpdT)
         if (qsat_tmp < forcn(2)%blk(1,1)%val(1,1)) forcn(2)%blk(1,1)%val(1,1) = qsat_tmp

      else
         do jb = 1, gblock%nyblk
            do ib = 1, gblock%nxblk
               if (gblock%pio(ib,jb) == p_iam_glb) then

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

                        case ('CRUNCEP')  

                           call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                              es,esdT,qsat_tmp,dqsat_tmpdT)
                           if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                              forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                           endif

                        case ('ERA5LAND')

                           call qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                              es,esdT,qsat_tmp,dqsat_tmpdT)
                           if (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) then
                              forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                           endif

                        end select

                     end do
                  end do
               end if
            end do
         end do
      end if

   END SUBROUTINE metpreprocess

END MODULE user_specified_forcing
! ----------- EOP ---------------
