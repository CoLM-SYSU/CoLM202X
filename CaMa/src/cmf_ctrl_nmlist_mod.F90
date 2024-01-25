MODULE CMF_CTRL_NMLIST_MOD
!==========================================================
!* PURPOSE: read CaMa-Flood Model configulations from namelist ("input_flood.nam" as default)
!
!* CONTAINS:
! -- CMF_CONFIG_NAMELIST  : read namelist for CaMa-Flood 
! -- CMF_CONFIG_CHECK     : check config conflict
!
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  29Jul 2019
!             Adapted mostly from CMF v362 CONTROL0.F90

  ! Modified by Zhongwang Wei @ SYSU 2022.11.20: add water re-infiltration calculation 

! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
! shared variables in module
   USE PARKIND1,                only: JPIM, JPRB, JPRM
   USE YOS_CMF_INPUT,           only: LOGNAM
   USE YOS_CMF_MAP,             only: REGIONTHIS, REGIONALL

   IMPLICIT NONE
CONTAINS
   !####################################################################
   ! -- CMF_CONFIG_NAMELIST  : read namelist for CaMa-Flood 
   ! -- CMF_CONFIG_CHECK     : check config conflict
   !
   !####################################################################
   SUBROUTINE CMF_CONFIG_NMLIST
   USE YOS_CMF_INPUT,      only: TMPNAM,   NSETFILE,   CSETFILE
   ! run version
   USE YOS_CMF_INPUT,      only: LADPSTP,  LFPLAIN,  LKINE,    LFLDOUT,  LPTHOUT,  LDAMOUT,  &
                              & LROSPLIT, LGDWDLY,  LSLPMIX,  LMEANSL,  LSEALEV,  LOUTPUT,  &
                              & LRESTART, LSTOONLY, LGRIDMAP, LLEAPYR,  LMAPEND,  LBITSAFE, &
                              & LSTG_ES,  LLEVEE,   LOUTINS,  LOUTINI,  LSEDOUT, &
                              & LSLOPEMOUTH,LWEVAP,LWINFILT, LWEVAPFIX,LWINFILTFIX,LWEXTRACTRIV
   ! dimention & time
   USE YOS_CMF_INPUT,      only: CDIMINFO, DT,       NX,NY,    NLFP,     NXIN,NYIN,    INPN, &
                              & IFRQ_INP, DTIN,     WEST,EAST,NORTH,SOUTH
   ! parameters
   USE YOS_CMF_INPUT,      only: PMANRIV,  PMANFLD,  PDSTMTH,  PMINSLP,  PGRV, PCADP, &
                              & IMIS, RMIS, DMIS,   CSUFBIN,  CSUFVEC,  CSUFPTH,  CSUFCDF
   USE CMF_UTILS_MOD,      only: INQUIRE_FID
   IMPLICIT NONE
   !* local
   character(LEN=8)              :: CREG                 !! 
   !
   NAMELIST/NRUNVER/  LADPSTP,  LFPLAIN,  LKINE,    LFLDOUT,  LPTHOUT,  LDAMOUT,  &
                     LROSPLIT, LGDWDLY,  LSLPMIX,  LMEANSL,  LSEALEV,  LOUTPUT,  &
                     LRESTART, LSTOONLY, LGRIDMAP, LLEAPYR,  LMAPEND,  LBITSAFE, &
                     LSTG_ES,  LLEVEE,   LSEDOUT,  LOUTINS,  LSLOPEMOUTH,        &
                     LWEVAP, LWINFILT,  LWEVAPFIX,LWINFILTFIX,LWEXTRACTRIV,       LOUTINI

   NAMELIST/NDIMTIME/ CDIMINFO, DT, IFRQ_INP

   NAMELIST/NPARAM/   PMANRIV, PMANFLD, PGRV,    PDSTMTH, PCADP,   PMINSLP, &
                     IMIS, RMIS, DMIS, CSUFBIN, CSUFVEC, CSUFPTH, CSUFCDF
      !================================================
      write(LOGNAM,*) ""
      write(LOGNAM,*) "!--------------------"

      ! *** 0. SET INPUT UNIT AND open FILE 
      NSETFILE=INQUIRE_FID()               !!  for namelist
      open(NSETFILE,FILE=CSETFILE,STATUS="OLD")
      write(LOGNAM,*) "CMF::CONFIG_NMLIST: namelist opened: ", TRIM(CSETFILE), NSETFILE 

      !============================
      !*** 1. basic simulation run version

      ! * defaults
      LADPSTP  = .TRUE.            !! true: use adaptive time step
      LFPLAIN  = .TRUE.            !! true: consider floodplain (false: only river channel)
      LKINE    = .FALSE.           !! true: use kinematic wave
      LFLDOUT  = .TRUE.            !! true: floodplain flow (high-water channel flow) active
      LPTHOUT  = .FALSE.           !! true: activate bifurcation scheme
      LDAMOUT  = .FALSE.           !! true: activate dam operation (under development)
      LLEVEE   = .FALSE.           !! true: activate levee scheme  (under development)
      LSEDOUT  = .FALSE.           !! true: activate sediment transport (under development)
      LOUTINS  = .FALSE.           !! true: diagnose instantaneous discharge
      !!=== this part is used by ECMWF
      LROSPLIT = .FALSE.           !! true: input if surface (Qs) and sub-surface (Qsb) runoff
      LWEVAP   = .FALSE.           !! true: input evaporation to extract from river 
      LWINFILT   = .FALSE.         !! true: input infiltration to extract from river 
      LWEVAPFIX= .FALSE.           !! true: water balance closure extracting water from evap when available
      LWINFILTFIX= .FALSE.         !! true: water balance closure extracting water from infiltration when available
      LGDWDLY  = .FALSE.           !! true: Activate ground water reservoir and delay
      LSLPMIX  = .FALSE.           !! true: activate mixed kinematic and local inertia based on slope
      LWEXTRACTRIV=.FALSE.         !! true: also extract water from rivers 
      LSLOPEMOUTH =.FALSE.         !! true: prescribe water level slope == elevation slope on river month
      !!===

      !! dinamic sea level
      LMEANSL  = .FALSE.           !! true : boundary condition for mean sea level
      LSEALEV  = .FALSE.           !! true : boundary condition for variable sea level

      !! restaer & output
      LRESTART = .FALSE.           !! true: initial condition from restart file
      LSTOONLY = .FALSE.           !! true: storage only restart (mainly for data assimilation)
      LOUTPUT  = .TRUE.            !! true: use standard output (to file)
      LOUTINI  = .FALSE.           !! true: output initial storage (netCDF only)

      LGRIDMAP = .TRUE.            !! true: for standard XY gridded 2D map
      LLEAPYR  = .FALSE.            !! true: neglect leap year (Feb29 skipped)
      LMAPEND  = .FALSE.           !! true: for map data endian conversion
      LBITSAFE = .FALSE.           !! true: for Bit Identical (not used from v410, set in Mkinclude)
      LSTG_ES  = .FALSE.           !! true: for Vector Processor optimization (CMF_OPT_FLDSTG_ES) 

      !* change
      rewind(NSETFILE)
      read(NSETFILE,NML=NRUNVER)

      write(LOGNAM,*) ""
      write(LOGNAM,*) "=== NAMELIST, NRUNVER ==="
      write(LOGNAM,*) "LADPSTP ",  LADPSTP
      write(LOGNAM,*) "LFPLAIN ",  LFPLAIN
      write(LOGNAM,*) "LKINE   ",  LKINE
      write(LOGNAM,*) "LFLDOUT ",  LFLDOUT
      write(LOGNAM,*) "LPTHOUT ",  LPTHOUT
      write(LOGNAM,*) "LDAMOUT ",  LDAMOUT
      write(LOGNAM,*) "LLEVEE  ",  LLEVEE
      write(LOGNAM,*) "LSEDOUT ",  LSEDOUT
      write(LOGNAM,*) "LOUTINS ",  LOUTINS
      write(LOGNAM,*) ""
      write(LOGNAM,*) "LROSPLIT ", LROSPLIT
      write(LOGNAM,*) "LWEVAP   ", LWEVAP
      write(LOGNAM,*) "LWEVAPFIX", LWEVAPFIX
      write(LOGNAM,*) "LWINFILT   ", LWINFILT
      write(LOGNAM,*) "LWINFILTFIX", LWINFILTFIX
      write(LOGNAM,*) "LWEXTRACTRIV", LWEXTRACTRIV
      write(LOGNAM,*) "LGDWDLY  ", LGDWDLY
      write(LOGNAM,*) "LSLPMIX  ", LSLPMIX
      write(LOGNAM,*) "LSLOPEMOUTH ", LSLOPEMOUTH
      write(LOGNAM,*) ""
      write(LOGNAM,*) "LMEANSL: ", LSEALEV
      write(LOGNAM,*) "LSEALEV: ", LSEALEV
      write(LOGNAM,*) ""
      write(LOGNAM,*) "LRESTART ", LRESTART
      write(LOGNAM,*) "LSTOONLY ", LSTOONLY
      write(LOGNAM,*) "LOUTPUT  ", LOUTPUT
      write(LOGNAM,*) "LOUTINI  ", LOUTINI
      write(LOGNAM,*) ""
      write(LOGNAM,*) "LGRIDMAP ", LGRIDMAP
      write(LOGNAM,*) "LLEAPYR  ", LLEAPYR
      write(LOGNAM,*) "LMAPEND  ", LMAPEND
      write(LOGNAM,*) "LBITSAFE ", LBITSAFE
      write(LOGNAM,*) "LSTG_ES " , LSTG_ES

      !============================
      !*** 2. set model dimention & time

      !* defaults (from namelist)
      CDIMINFO ="NONE"
      DT       = 24*60*60          !! dt = 1day (automatically set by adaptive time step)
      IFRQ_INP = 24                !! daily (24h) input

      !* change
      rewind(NSETFILE)
      read(NSETFILE,NML=NDIMTIME)

      DTIN  = IFRQ_INP*60*60       !! hour -> second

      write(LOGNAM,*) ""
      write(LOGNAM,*) "=== NAMELIST, NCONF ==="
      write(LOGNAM,*) "CDIMINFO  ", TRIM(CDIMINFO)
      write(LOGNAM,*) "DT        ", DT
      write(LOGNAM,*) "DTIN      ", DTIN
      write(LOGNAM,*) "IFRQ_INP  ", IFRQ_INP

      !==========
      !* default (from diminfo)
      NX    = 1440                 !! 15 minute resolution
      NY    = 720
      NLFP  = 10                   !! 10 floodplain layer
      NXIN  = 360                  !! 1 degree input
      NYIN  = 180
      INPN  = 1                    !! maximum number of input grids corresponding to one CaMa-Flood grid
      WEST  = -180._JPRB              !! west, east, north, south edges of the domain
      EAST  =  180._JPRB
      NORTH =  90._JPRB
      SOUTH = -90._JPRB

      !* value from CDIMINFO
      IF( CDIMINFO/="NONE" )THEN
         write(LOGNAM,*) "CMF::CONFIG_NMLIST: read DIMINFO ", TRIM(CDIMINFO)

         TMPNAM=INQUIRE_FID()
         open(TMPNAM,FILE=CDIMINFO,FORM='FORMATTED')
         read(TMPNAM,*) NX
         read(TMPNAM,*) NY
         read(TMPNAM,*) NLFP
         read(TMPNAM,*) NXIN
         read(TMPNAM,*) NYIN
         read(TMPNAM,*) INPN
         read(TMPNAM,*) 
         IF( LGRIDMAP )THEN
            read(TMPNAM,*) WEST
            read(TMPNAM,*) EAST
            read(TMPNAM,*) NORTH
            read(TMPNAM,*) SOUTH
         ENDIF
         close(TMPNAM)
      ENDIF

      !* check
      write(LOGNAM,*) ""
      write(LOGNAM,*) "=== DIMINFO ==="
      write(LOGNAM,*) "NX,NY,NLFP     ", NX,  NY,  NLFP
      write(LOGNAM,*) "NXIN,NYIN,INPN ", NXIN,NYIN,INPN
      IF( LGRIDMAP ) THEN
         write(LOGNAM,*) "WEST,EAST,NORTH,SOUTH ", WEST,EAST,NORTH,SOUTH
      ENDIF

      !============================
      !*** 3. set PARAM: parameters
      ! * defaults
      PMANRIV=0.03_JPRB                              !! manning coefficient river
      PMANFLD=0.10_JPRB                              !! manning coefficient floodplain
      PGRV   =9.8_JPRB                               !! gravity accerelation
      PDSTMTH=10000._JPRB                            !! downstream distance at river mouth [m]
      PCADP  =0.7_JPRB                               !! CFL coefficient
      PMINSLP=1.E-5                                  !! minimum slope (kinematic wave)

      IMIS=-9999_JPIM
      RMIS=1.E20_JPRM
      DMIS=1.E20_JPRB

      CSUFBIN='.bin'
      CSUFVEC='.vec'
      CSUFPTH='.pth'
      CSUFCDF='.nc'

      ! * change
      rewind(NSETFILE)
      read(NSETFILE,NML=NPARAM)

      write(LOGNAM,*) ""
      write(LOGNAM,*) "=== NAMELIST, NPARAM ==="
      write(LOGNAM,*) "PMANRIV  ", PMANRIV
      write(LOGNAM,*) "PMANRIV  ", PMANFLD
      write(LOGNAM,*) "PGRV     ", PGRV
      write(LOGNAM,*) "PDSTMTH  ", PDSTMTH
      write(LOGNAM,*) "PCADP    ", PCADP
      write(LOGNAM,*) "PMINSLP  ", PMINSLP
      write(LOGNAM,*) ""
      write(LOGNAM,*) "IMIS     ", IMIS
      write(LOGNAM,*) "RMIS     ", RMIS
      write(LOGNAM,*) "DMIS     ", DMIS
      write(LOGNAM,*) ""
      write(LOGNAM,*) "CSUFBIN  ", TRIM(CSUFBIN)
      write(LOGNAM,*) "CSUFVEC  ", TRIM(CSUFVEC)
      write(LOGNAM,*) "CSUFPTH  ", TRIM(CSUFPTH)
      write(LOGNAM,*) "CSUFCDF  ", TRIM(CSUFCDF)

      !===============================
      !*** close FILE 
      close(NSETFILE)

      write(LOGNAM,*) "CMF::CONFIG_NMLIST: end "

      write(LOGNAM,*) "--------------------!"
      write(LOGNAM,*) ""

      IF (REGIONALL>=2 )THEN 
         write(CREG,'(I0)') REGIONTHIS                !! Regional Output for MPI run
         CSUFVEC=TRIM(CSUFVEC)//'-'//TRIM(CREG)       !! Change suffix of output file for each MPI node (only vector output)
      ENDIF

   END SUBROUTINE CMF_CONFIG_NMLIST
   !####################################################################


   !####################################################################
   SUBROUTINE CMF_CONFIG_CHECK
   USE YOS_CMF_INPUT,      only: LADPSTP,  LFPLAIN,  LKINE,    LPTHOUT,     &
                              & LROSPLIT, LGDWDLY,  LSEALEV, &
                              & LWEVAP,   LWEVAPFIX, LWINFILT,   LWINFILTFIX,  LWEXTRACTRIV
   USE YOS_CMF_INPUT,      only: DT, DTIN, DTSL
   IMPLICIT NONE
   !================================================

      write(LOGNAM,*) "CMF::CONFIG_CHECK: check setting conflicts"

      !*** 1. check for time step
      IF ( DT<60 .or. MOD( INT(DT),60 )/=0 ) THEN
         write(LOGNAM,*) "DT= ", DT
         write(LOGNAM,*) "DT should be multiple of 60. CaMa-Flood controls time by MINUTE"
         write(LOGNAM,*) "stop"
         STOP 9
      ENDIF

      IF ( MOD( INT(DTIN), INT(DT) )/=0 ) THEN
         write(LOGNAM,*) "DTIN, DT= ", DTIN, DT
         write(LOGNAM,*) "DTIN should be multiple of DT"
         write(LOGNAM,*) "stop"
         STOP 9
      ENDIF

      IF ( LSEALEV .and. MOD( INT(DTSL), INT(DT) )/=0 ) THEN
         write(LOGNAM,*) "DTSL, DT= ", DTIN, DT
         write(LOGNAM,*) "DTSL should be multiple of DT"
         write(LOGNAM,*) "stop"
         STOP 9
      ENDIF

      !*** 2. check for physics options

      IF ( .not.LFPLAIN .and. .not.LKINE ) THEN
         write(LOGNAM,*) "LFPLAIN=.false. & LKINE=.false."
         write(LOGNAM,*) "CAUTION: NO FLOODPLAIN OPTION reccomended to be used with kinematic wave (LKINE=.true.)"
      ENDIF

      IF ( LKINE .and. LADPSTP ) THEN
         write(LOGNAM,*) "LKINE=.true. & LADPSTP=.true."
         write(LOGNAM,*) "adaptive time step reccoomended only with local inertial equation (LKINE=.false.)"
         write(LOGNAM,*) "Set appropriate fixed time step for Kinematic Wave"
      ENDIF

      IF ( LKINE .and. LPTHOUT ) THEN
         write(LOGNAM,*) "LKINE=.true. & LPATHOUT=.true."
         write(LOGNAM,*) "bifurcation channel flow only available with local inertial equation (LKINE=.false.)"
         write(LOGNAM,*) "STOP"
         STOP 9
      ENDIF

      IF ( LGDWDLY .and. .not. LROSPLIT ) THEN
         write(LOGNAM,*) "LGDWDLY=true and LROSPLIT=false"
         write(LOGNAM,*) "Ground water reservoir can only be active when runoff splitting is on"
      ENDIF

      IF ( LWEVAPFIX .and. .not. LWEVAP ) THEN
         write(LOGNAM,*) "LWEVAPFIX=true and LWEVAP=false"
         write(LOGNAM,*) "LWEVAPFIX can only be active if LWEVAP is active"
      ENDIF 

      !  add water re-infiltration calculation 
      IF ( LWINFILTFIX .and. .not. LWINFILT ) THEN
         write(LOGNAM,*) "LWINFILTFIX=true and LWINFILT=false"
         write(LOGNAM,*) "LWINFILTFIX can only be active if LWINFILT is active"
      ENDIF

      IF ( LWEXTRACTRIV .and. .not. LWEVAP ) THEN
         write(LOGNAM,*) "LWEXTRACTRIV=true and LWEVAP=false"
         write(LOGNAM,*) "LWEXTRACTRIV can only be active if LWEVAP is active"
      ENDIF 

      write(LOGNAM,*) "CMF::CONFIG_CHECK: end"

   END SUBROUTINE CMF_CONFIG_CHECK
!####################################################################

END MODULE CMF_CTRL_NMLIST_MOD
