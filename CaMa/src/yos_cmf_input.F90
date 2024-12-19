MODULE YOS_CMF_INPUT
!==========================================================
!* PURPOSE: Shared variables for CaMa-Flood model configulation
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
!

! Modified by Zhongwang Wei @ SYSU 2022.11.20: add water re-infiltration calculation 

! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
   USE PARKIND1,                ONLY: JPIM, JPRB, JPRM
   IMPLICIT NONE
   SAVE 
   !================================================
   !*** CMF default files
   logical                         :: LLOGOUT                 !! true: log output to file
   integer(KIND=JPIM)              :: LOGNAM                  !! default log    file FID
   integer(KIND=JPIM)              :: NSETFILE                !! input namelist file FID
   integer(KIND=JPIM)              :: TMPNAM                  !! temporal I/O   file FIG
   character(LEN=256)              :: CLOGOUT                 !! default log    file name
   character(LEN=256)              :: CSETFILE                !! input namelist file name

   DATA LLOGOUT       /.TRUE./
   DATA CLOGOUT       /'./log_CaMa.txt'/
   DATA CSETFILE      /'../run/cama_flood.nml'/

   !================================================
   !*** NAMELIST/NRUNVER/
   logical                         :: LADPSTP                 !! true: use adaptive time step

   logical                         :: LFPLAIN                 !! true: consider floodplain (false: only river channel)
   logical                         :: LKINE                   !! true: use kinematic wave
   logical                         :: LFLDOUT                 !! true: floodplain flow (high-water channel flow) active
   logical                         :: LPTHOUT                 !! true: activate bifurcation scheme
   logical                         :: LDAMOUT                 !! true: activate dam operation (under development)
   logical                         :: LLEVEE                  !! true: activate levee scheme  (under development)

   !~~ used in ECMWF
   logical                         :: LROSPLIT                !! true: input if surface (Qs) and sub-surface (Qsb) runoff
   logical                         :: LWEVAP                  !! true: input water evaporation to extract from floodplain
   logical                         :: LWEVAPFIX               !! true: water balance closure extracting water from evap when available
   logical                         :: LWINFILT                !! true: input water infiltration to extract from floodplain
   logical                         :: LWINFILTFIX             !! true: water balance closure extracting water from Infiltration when available
   logical                         :: LWEXTRACTRIV            !! true: also extract water from rivers 
   logical                         :: LSLOPEMOUTH             !! true: prescribe water level slope == elevation slope on river month
   logical                         :: LGDWDLY                 !! true: Activate ground water reservoir and delay
   logical                         :: LSLPMIX                 !! true: activate mixed kinematic and local inertia based on slope

   logical                         :: LMEANSL                 !! true : boundary condition for mean sea level
   logical                         :: LSEALEV                 !! true : boundary condition for variable sea level

   logical                         :: LOUTINS                 !! true: diagnose instantaneous discharge 

   logical                         :: LRESTART                !! true: initial condition from restart file
   logical                         :: LSTOONLY                !! true: storage only restart (mainly for data assimilation)

   logical                         :: LOUTPUT                 !! true: use standard output (to file)
   logical                         :: LOUTINI                 !! true: output initial storage (netCDF only)

   logical                         :: LGRIDMAP                !! true: for standard XY gridded 2D map
   logical                         :: LLEAPYR                 !! true: neglect leap year (Feb29 skipped)
   logical                         :: LMAPEND                 !! true: for map data endian conversion
   logical                         :: LBITSAFE                !! true: for Bit Identical (removed from v410, set in Mkinclude)
   logical                         :: LSTG_ES                 !! true: for Vector Processor optimization (CMF_OPT_FLDSTG_ES) 

   logical                         :: LSEDOUT                 !! true: sediment scheme                 

   !================================================
   !*** NAMELIST/NCONF/
   character(LEN=256)              :: CDIMINFO                !! Dimention Information

   real(KIND=JPRB)                 :: DT                      !! Time Step Length [SEC] (should be multiple of 60)

   integer(KIND=JPIM)              :: IFRQ_OUT                !! [hour]: frequency to write output     e.g. (1,2,3,6,12,24) hour
   integer(KIND=JPIM)              :: IFRQ_INP                !! [hour]: frequency to update runoff    e.g. (1,2,3,6,12,24) hour
   integer(KIND=JPIM)              :: IFRQ_SL                 !! [min]:  frequency to update sea level e.g. (1,2,5,10,15,20,30,60) min

   !*** set by CDIMINFO
   integer(KIND=JPIM)              :: NX                      !! NUMBER OF GRIDS IN HORIZONTAL
   integer(KIND=JPIM)              :: NY                      !! NUMBER OF GRIDS IN VERTICAL
   integer(KIND=JPIM)              :: NLFP                    !! NUMBER OF VERTICAL LEVELS DEFINING FLOODPLAIN

   integer(KIND=JPIM)              :: NXIN                    !! NUMBER OF GRIDS IN HORIZONTAL
   integer(KIND=JPIM)              :: NYIN                    !! NUMBER OF GRIDS IN VERTICAL
   integer(KIND=JPIM)              :: INPN                    !! MAX INPUT NUMBER

   real(KIND=JPRB)                 :: WEST                    !! west, east, north, south edge of the domain [deg]
   real(KIND=JPRB)                 :: EAST
   real(KIND=JPRB)                 :: NORTH
   real(KIND=JPRB)                 :: SOUTH

   !*** calculated from IFRQ & DT
   real(KIND=JPRB)                 :: DTIN                    !! SECOND IN INPUT TIME STEP [SEC]
   real(KIND=JPRB)                 :: DTSL                    !! SECOND IN TIME STEP [SEC]

   !================================================
   !*** NAMELIST/PARAM/
   !* parameters
   real(KIND=JPRB)                 :: PMANRIV              !! manning roughness (river)
   real(KIND=JPRB)                 :: PMANFLD              !! manning roughness (floodplain)
   real(KIND=JPRB)                 :: PGRV                 !! gravity acceleration [m/s2]
   real(KIND=JPRB)                 :: PDSTMTH              !! downstream distance at river mouth [m]
   real(KIND=JPRB)                 :: PCADP                !! CFL coefficient
   real(KIND=JPRB)                 :: PMINSLP              !! minimum topographic slope (kinematic wave) [m/m]
   !* missing values
   integer(KIND=JPIM)              :: IMIS                !! integer undefined
   real(KIND=JPRM)                 :: RMIS                !! real    undefined
   real(KIND=JPRB)                 :: DMIS                !! double  undefined
   !* file suffix
   character(LEN=256)              :: CSUFBIN            ! .bin suffix for binary (2D map)
   character(LEN=256)              :: CSUFVEC            ! .vec suffix for binary (1D vector)
   character(LEN=256)              :: CSUFPTH            ! .pth suffix for binary (1D bifurcation channel)
   character(LEN=256)              :: CSUFCDF            ! .nc  suffix for netCDF
!================================================
#ifdef IFS_CMF
   ! Fluxes buffers for IFS coupling
   real(KIND=JPRB), ALLOCATABLE :: ZBUFFO(:,:,:)
   real(KIND=JPRB), ALLOCATABLE :: ZBUFFI(:,:,:)
   real(KIND=JPRB), ALLOCATABLE :: ZACC0(:,:)
   real(KIND=JPRB), ALLOCATABLE :: ZACC1(:,:)
   !Time step to be advanced within DRV_ADVANCE used for IFS coupling
   integer(KIND=JPIM)           :: ISTEPADV
   !===============================================
#endif

END MODULE YOS_CMF_INPUT
