MODULE CMF_CTRL_OUTPUT_MOD
!==========================================================
!* PURPOSE: Control CaMa-Flood standard output file (binary / netCDF)  
!
!* CONTAINS:
! -- CMF_OUTPUT_NMLIST : Read output file info from namelist
! -- CMF_OUTPUT_INIT   : Create & Open standard output files
! -- CMF_OUTPUT_WRITE  : Write output to files
! -- CMF_OUTPUT_END    : Close standard output files
!
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
! shared variables in MODULE
   USE PARKIND1,                only: JPIM, JPRB, JPRM
   USE YOS_CMF_INPUT,           only: LOGNAM,  IFRQ_OUT
   USE YOS_CMF_INPUT,           only: CSUFBIN, CSUFVEC, CSUFPTH, CSUFCDF
   USE YOS_CMF_INPUT,           only: LPTHOUT, LDAMOUT, LLEVEE,  LWEVAP, LWINFILT,LGDWDLY, LOUTINS,LROSPLIT
   IMPLICIT NONE
   !============================
   SAVE
   !*** NAMELIST/NOUTPUT/ from inputnam
   character(LEN=256)              ::  COUTDIR           ! OUTPUT DIRECTORY
   character(LEN=256)              ::  CVARSOUT          ! Comma-separated list of output variables to save 
   character(LEN=256)              ::  COUTTAG           ! Output Tag Name for each experiment
   !
   logical                         ::  LOUTVEC           ! TRUE FOR VECTORIAL OUTPUT, FALSE FOR NX,NY OUTPUT
   logical                         ::  LOUTCDF           ! true for netcdf outptu false for binary
   integer(KIND=JPIM)              ::  NDLEVEL           ! NETCDF DEFLATION LEVEL 
   !
   logical                         ::  LOUTTXT           ! TRUE FOR Text output for some gauges
   character(LEN=256)              ::  CGAUTXT           ! List of Gauges (ID, IX, IY)
   !
   NAMELIST/NOUTPUT/ COUTDIR,CVARSOUT,COUTTAG,LOUTCDF,NDLEVEL,LOUTVEC,IFRQ_OUT,LOUTTXT,CGAUTXT
   !
   !*** local variables
   integer(KIND=JPIM)              :: NVARS              ! temporal output var number
   parameter                         (NVARS=30)          ! actual   output var number
   integer(KIND=JPIM)              :: NVARSOUT
   integer(KIND=JPIM)              :: IRECOUT            ! Output file irec

   !*** type for output file    
   type TVAROUT
      character(LEN=256)              :: CVNAME             ! output variable name
      character(LEN=256)              :: CVLNAME            ! output variable long name
      character(LEN=256)              :: CVUNITS            ! output units
      character(LEN=256)              :: CFILE              ! output full path file name 
      integer(KIND=JPIM)              :: BINID              ! output binary output file ID
      integer(KIND=JPIM)              :: NCID               ! output netCDF output file ID
      integer(KIND=JPIM)              :: VARID              ! output netCDF output variable ID
      integer(KIND=JPIM)              :: TIMID              ! output netCDF time   variable ID 
      integer(KIND=JPIM)              :: IRECNC               ! Current time record for writting 
   END type TVAROUT 
   type(TVAROUT),ALLOCATABLE       :: VAROUT(:)          ! output variable type set

CONTAINS
   !####################################################################
   ! -- CMF_OUTPUT_NMLIST : Read output file info from namelist
   ! -- CMF_OUTPUT_INIT   : Create & Open standard output files
   ! -- CMF_OUTPUT_WRITE  : Write output to files
   ! -- CMF_OUTPUT_END    : Close standard output files
   ! --
   !####################################################################
   SUBROUTINE CMF_OUTPUT_NMLIST
   ! reed setting from namelist
   ! -- Called from CMF_DRV_NMLIST
   USE YOS_CMF_INPUT,      only: CSETFILE,NSETFILE
   USE CMF_UTILS_MOD,      only: INQUIRE_FID
   IMPLICIT NONE
      !================================================
      write(LOGNAM,*) ""
      write(LOGNAM,*) "!---------------------!"

      !*** 1. open namelist
      NSETFILE=INQUIRE_FID()
      open(NSETFILE,FILE=CSETFILE,STATUS="OLD")
      write(LOGNAM,*) "CMF::OUTPUT_NMLIST: namelist open in unit:", TRIM(CSETFILE), NSETFILE 

      !*** 2. default value
      COUTDIR="./"
      CVARSOUT="outflw,storge,rivdph"
      COUTTAG="_cmf"
      LOUTCDF=.FALSE.
      NDLEVEL=0
      LOUTVEC=.FALSE.
      IFRQ_OUT = 24                !! daily (24h) output
      !
      LOUTTXT=.FALSE.
      CGAUTXT="None"

      !*** 3. read namelist
      rewind(NSETFILE)
      read(NSETFILE,NML=NOUTPUT)

      write(LOGNAM,*)   "=== NAMELIST, NOUTPUT ==="
      write(LOGNAM,*)   "COUTDIR:  ", TRIM(COUTDIR)
      write(LOGNAM,*)   "CVARSOUT: ", TRIM(CVARSOUT)
      write(LOGNAM,*)   "COUTTAG:  ", TRIM(COUTTAG)

      write(LOGNAM,*)   "LOUTCDF:  ", LOUTCDF
      IF( LOUTCDF )THEN
         write(LOGNAM,*) "NDLEVEL:  ", NDLEVEL
      ENDIF
      IF( LOUTVEC )THEN
         write(LOGNAM,*) "LOUTVEC:  ", LOUTVEC
      ENDIF
      write(LOGNAM,*)   "IFRQ_OUT  ", IFRQ_OUT

      write(LOGNAM,*)   "IFRQ_OUT  ", LOUTTXT
      write(LOGNAM,*)   "CGAUTXRT  ", CGAUTXT

      close(NSETFILE)

      write(LOGNAM,*) "CMF::OUTPUT_NMLIST: end"

   END SUBROUTINE CMF_OUTPUT_NMLIST
   !####################################################################





   !####################################################################
   SUBROUTINE CMF_OUTPUT_INIT
   ! Initialize output module (create/open files)
   ! -- Called from CMF_DRV_INIT
   USE YOS_CMF_INPUT,           only: NX,NY
   USE YOS_CMF_TIME,            only: ISYYYY, ISMM,   ISDD,   ISHOUR, ISMIN
   USE YOS_CMF_MAP,             only: NSEQMAX,NPTHOUT,NPTHLEV,REGIONTHIS
   USE CMF_UTILS_MOD,           only: INQUIRE_FID
   IMPLICIT NONE
   !* Local variables 
   character(LEN=256)              :: CTIME, CTMP
   integer(KIND=JPIM)              :: JF,J,J0
   character(LEN=256)              :: CVNAMES(NVARS)
   !================================================
      write(LOGNAM,*) ""
      write(LOGNAM,*) "!---------------------!"

      write(LOGNAM,*) "CMF::OUTPUT_INIT: check output variables"
      !! Start by finding out # of output variables 
      NVARSOUT=0
      J0=1
      DO J=1,LEN(TRIM(CVARSOUT))
         IF( (J>J0) .and. (CVARSOUT(J:J) .eq. ',') ) THEN
            CTMP=TRIM(ADJUSTL(CVARSOUT(J0:J-1)))
            IF (LEN(CTMP) > 0 ) THEN
               NVARSOUT=NVARSOUT+1
               CVNAMES(NVARSOUT)=CTMP
            ENDIF
            J0=J+1
         ENDIF
      ENDDO
      ! Last one 
      IF ( J0 < LEN(TRIM(CVARSOUT)) ) THEN
         J=LEN(TRIM(CVARSOUT))
         CTMP=TRIM(ADJUSTL(CVARSOUT(J0:J)))
            IF (LEN(CTMP) > 0 ) THEN
               NVARSOUT=NVARSOUT+1
               CVNAMES(NVARSOUT)=CTMP
            ENDIF
      ENDIF 

      IF ( NVARSOUT == 0 ) THEN
         write(LOGNAM,*) "CMF::OUTPUT_INIT: No output files will be produced!"
         RETURN
      ENDIF 

      allocate(VAROUT(NVARSOUT))
      write(CTIME,'(A14,I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') 'seconds since ',ISYYYY,'-',ISMM,'-',ISDD,' ',ISHOUR,":",ISMIN

      !* Loop on variables and create files 
      !  add water re-infiltration calculation
      ! currently was not used in colm-cama coupling model

      DO JF=1,NVARSOUT
         write(LOGNAM,*) "Creating output for variable:", TRIM( CVNAMES(JF) )
         SELECT CASE (CVNAMES(JF))
            CASE ('rivout')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='river discharge'
               VAROUT(JF)%CVUNITS='m3/s'
            CASE ('rivsto')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='river storage'
               VAROUT(JF)%CVUNITS='m3'
            CASE ('rivdph')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='river depth'
               VAROUT(JF)%CVUNITS='m'
            CASE ('rivvel')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='river velocity'
               VAROUT(JF)%CVUNITS='m/s'
            CASE ('fldout')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='floodplain discharge'
               VAROUT(JF)%CVUNITS='m3/s'
            CASE ('fldsto')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='floodplain storage'
               VAROUT(JF)%CVUNITS='m3'
            CASE ('flddph')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='floodplain depth'
               VAROUT(JF)%CVUNITS='m'  
            CASE ('fldfrc')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='flooded fraction'
               VAROUT(JF)%CVUNITS='0-1'  
            CASE ('fldare')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='flooded area'
               VAROUT(JF)%CVUNITS='m2'
            CASE ('sfcelv')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='water surface elevation'
               VAROUT(JF)%CVUNITS='m'
            CASE ('totout')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='discharge (river+floodplain)'
               VAROUT(JF)%CVUNITS='m3/s'
            CASE ('outflw')                   !! comparability for previous output name
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='discharge (river+floodplain)'
               VAROUT(JF)%CVUNITS='m3/s'
            CASE ('totsto')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='total storage (river+floodplain)'
               VAROUT(JF)%CVUNITS='m3'
            CASE ('storge')                   !! comparability for previous output name
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='total storage (river+floodplain)'
               VAROUT(JF)%CVUNITS='m3'
            CASE ('pthflw')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='bifurcation channel discharge'
               VAROUT(JF)%CVUNITS='m3/s'
            CASE ('pthout')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='net bifurcation discharge'
               VAROUT(JF)%CVUNITS='m3/s'
            CASE ('maxsto')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='daily maximum storage'
               VAROUT(JF)%CVUNITS='m3'  
            CASE ('maxflw')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='daily maximum discharge'
               VAROUT(JF)%CVUNITS='m3/s' 
            CASE ('maxdph')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='daily maximum river depth'
               VAROUT(JF)%CVUNITS='m' 
            CASE ('runoff')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='Surface runoff'
               VAROUT(JF)%CVUNITS='m3/s' 
            CASE ('runoffsub')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='sub-surface runoff'
               VAROUT(JF)%CVUNITS='m3/s' 
            CASE ('damsto')   !!! added
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='reservoir storage'
               VAROUT(JF)%CVUNITS='m3' 
            CASE ('daminf')   !!! added
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='reservoir inflow'
               VAROUT(JF)%CVUNITS='m3/s' 
            CASE ('levsto')   !!! added
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='protected area storage'
               VAROUT(JF)%CVUNITS='m3' 
            CASE ('levdph')   !!! added
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='protected area depth'
               VAROUT(JF)%CVUNITS='m' 
            CASE ('gdwsto')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='ground water storage'
               VAROUT(JF)%CVUNITS='m3'
            CASE ('gwsto')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='ground water storage'
               VAROUT(JF)%CVUNITS='m3'
            CASE ('gwout')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='ground water discharge'
               VAROUT(JF)%CVUNITS='m3/s'  
            CASE ('wevap')
               IF (LWEVAP) THEN
                  VAROUT(JF)%CVNAME=CVNAMES(JF)
                  VAROUT(JF)%CVLNAME='water evaporation'
                  VAROUT(JF)%CVUNITS='m3/s'
               ENDIF
            CASE ('winfilt')
               IF (LWINFILT) THEN
                  VAROUT(JF)%CVNAME=CVNAMES(JF)
                  VAROUT(JF)%CVLNAME='water infiltration'
                  VAROUT(JF)%CVUNITS='m3/s'
               ENDIF
            CASE ('outins')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='instantaneous discharge'
               VAROUT(JF)%CVUNITS='m3/s' 
            CASE ('outflw_ocean')
               VAROUT(JF)%CVNAME=CVNAMES(JF)
               VAROUT(JF)%CVLNAME='discharge to ocean'
               VAROUT(JF)%CVUNITS='m3/s'
            CASE DEFAULT
         write(LOGNAM,*) trim(CVNAMES(JF)), ' Not defined in CMF_CREATE_OUTCDF_MOD'

         END SELECT
         VAROUT(JF)%BINID=INQUIRE_FID()
      ENDDO

      IRECOUT=0  ! Initialize Output record to 1 (shared in netcdf & binary)

   CONTAINS
      !==========================================================
      !+ CREATE_OUTBIN
      !+ CREATE_OUTCDF
      !==========================================================
      SUBROUTINE CREATE_OUTBIN
      IMPLICIT NONE
         !================================================
         IF( TRIM(VAROUT(JF)%CVNAME)=='pthflw' ) THEN   !! bifurcation channel
            IF( REGIONTHIS==1 )THEN
               VAROUT(JF)%CFILE=TRIM(COUTDIR)//TRIM(VAROUT(JF)%CVNAME)//TRIM(COUTTAG)//TRIM(CSUFPTH)
               open(VAROUT(JF)%BINID,FILE=VAROUT(JF)%CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NPTHOUT*NPTHLEV)
               write(LOGNAM,*) "output file opened in unit: ", TRIM(VAROUT(JF)%CFILE), VAROUT(JF)%BINID
            ENDIF
         ELSEIF( LOUTVEC )THEN   !!  1D land only output
            VAROUT(JF)%CFILE=TRIM(COUTDIR)//TRIM(VAROUT(JF)%CVNAME)//TRIM(COUTTAG)//TRIM(CSUFVEC)
            open(VAROUT(JF)%BINID,FILE=VAROUT(JF)%CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NSEQMAX)
            write(LOGNAM,*) "output file opened in unit: ", TRIM(VAROUT(JF)%CFILE), VAROUT(JF)%BINID
         ELSE                   !!  2D default map output
            IF( REGIONTHIS==1 )THEN
               VAROUT(JF)%CFILE=TRIM(COUTDIR)//TRIM(VAROUT(JF)%CVNAME)//TRIM(COUTTAG)//TRIM(CSUFBIN)
               write(LOGNAM,*) "  -- ", TRIM(VAROUT(JF)%CFILE)
               open(VAROUT(JF)%BINID,FILE=VAROUT(JF)%CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
               write(LOGNAM,*) "output file opened in unit: ", TRIM(VAROUT(JF)%CFILE), VAROUT(JF)%BINID
            ENDIF
         ENDIF
      END SUBROUTINE CREATE_OUTBIN
      !==========================================================
      !+
      !+
      !+
      !==========================================================
      SUBROUTINE CREATE_OUTCDF
#ifdef UseCDF_CMF
      USE YOS_CMF_INPUT,           only: RMIS
      USE YOS_CMF_MAP,             only: D1LON, D1LAT
      USE CMF_UTILS_MOD,           only: NCERROR
      USE NETCDF
      IMPLICIT NONE
      integer(KIND=JPIM)  :: TIMEID,VARID,LATID,LONID
         !============
         VAROUT(JF)%IRECNC=1 ! initialize record current writting record to 1 

         !============
         VAROUT(JF)%CFILE=TRIM(COUTDIR)//'o_'//TRIM(VAROUT(JF)%CVNAME)//TRIM(COUTTAG)//TRIM(CSUFCDF)
         ! Create file 
         CALL NCERROR( NF90_CREATE(VAROUT(JF)%CFILE,NF90_NETCDF4,VAROUT(JF)%NCID),&
                     'CREATING FILE:'//TRIM(VAROUT(JF)%CFILE) )
         !=== set dimension ===
         CALL NCERROR( NF90_DEF_DIM(VAROUT(JF)%NCID, 'time', NF90_UNLIMITED, TIMEID) )
         CALL NCERROR( NF90_DEF_DIM(VAROUT(JF)%NCID, 'lat', NY, LATID) )
         CALL NCERROR( NF90_DEF_DIM(VAROUT(JF)%NCID, 'lon', NX, LONID) )

         !=== define variables ===
         CALL NCERROR( NF90_DEF_VAR(VAROUT(JF)%NCID, 'lat', NF90_FLOAT, (/LATID/), VARID) )
         CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VARID, 'long_name','latitude') )
         CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VARID, 'units','degrees_north') )

         CALL NCERROR( NF90_DEF_VAR(VAROUT(JF)%NCID, 'lon', NF90_FLOAT, (/LONID/), VARID) )
         CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VARID, 'long_name','longitude') )
         CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VARID, 'units','degrees_east') )

         CALL NCERROR( NF90_DEF_VAR(VAROUT(JF)%NCID, 'time', NF90_DOUBLE, (/TIMEID/), VAROUT(JF)%TIMID) ) 
         CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VAROUT(JF)%TIMID, 'long_name','time') )
         CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VAROUT(JF)%TIMID, 'units',CTIME) )

         !===
         CALL NCERROR( NF90_DEF_VAR(VAROUT(JF)%NCID, VAROUT(JF)%CVNAME, NF90_FLOAT, &
                     (/LONID,LATID,TIMEID/), VAROUT(JF)%VARID,DEFLATE_LEVEL=NDLEVEL),     &
                     'Creating Variable')

         CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VAROUT(JF)%VARID, 'long_name', TRIM(VAROUT(JF)%CVLNAME)) )
         CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VAROUT(JF)%VARID, 'units',     TRIM(VAROUT(JF)%CVUNITS)) )
         CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VAROUT(JF)%VARID, '_FillValue',RMIS) )

         CALL NCERROR( NF90_ENDDEF(VAROUT(JF)%NCID) )

         !=== put lon lat info ===
         CALL NCERROR ( NF90_INQ_VARID(VAROUT(JF)%NCID,'lon',VARID),'getting id' )
         CALL NCERROR( NF90_PUT_VAR(VAROUT(JF)%NCID,VARID,D1LON))

         CALL NCERROR ( NF90_INQ_VARID(VAROUT(JF)%NCID,'lat',VARID),'getting id' )
         CALL NCERROR( NF90_PUT_VAR(VAROUT(JF)%NCID,VARID,D1LAT))

         write(LOGNAM,*) 'CFILE: ',TRIM(VAROUT(JF)%CFILE),' CVAR:',TRIM(VAROUT(JF)%CVNAME),&
                        ' CLNAME: ',TRIM(VAROUT(JF)%CVLNAME),' CUNITS: ',TRIM(VAROUT(JF)%CVUNITS)
         write(LOGNAM,*) 'open in UNIT: ',VAROUT(JF)%NCID
#endif
      END SUBROUTINE CREATE_OUTCDF
   !==========================================================

   END SUBROUTINE CMF_OUTPUT_INIT
   !####################################################################


   !####################################################################
   SUBROUTINE CMF_OUTPUT_WRITE
   !======
   USE CMF_UTILS_MOD,           only: vecD2mapR
   ! save results to output files
   ! -- Called either from "MAIN/Coupler" or CMF_DRV_ADVANCE
   USE YOS_CMF_INPUT,      only: NX, NY, LOUTINI
   USE YOS_CMF_MAP,        only: NSEQMAX, NPTHOUT, NPTHLEV, REGIONTHIS
   USE YOS_CMF_TIME,       only: JYYYYMMDD, JHHMM, JHOUR, JMIN, KSTEP
   USE YOS_CMF_PROG,       only: P2RIVSTO,     P2FLDSTO,     P2GDWSTO, &
                              & P2DAMSTO,     P2LEVSTO,     D2COPY       !!! added
   USE YOS_CMF_DIAG,       only: D2RIVDPH,     D2FLDDPH,     D2FLDFRC,     D2FLDARE,     D2SFCELV,     D2STORGE, &
                              & D2OUTFLW_AVG, D2RIVOUT_AVG, D2FLDOUT_AVG, D2PTHOUT_AVG, D1PTHFLW_AVG,  &
                              & D2RIVVEL_AVG, D2GDWRTN_AVG, D2RUNOFF_AVG, D2ROFSUB_AVG, D2WEVAPEX_AVG,D2WINFILTEX_AVG, &
                              & D2OUTFLW_MAX, D2STORGE_MAX, D2RIVDPH_MAX, &
                              & D2DAMINF_AVG, D2OUTINS, D2LEVDPH   !!! added
#ifdef UseMPI_CMF
   USE CMF_CTRL_MPI_MOD,   only: CMF_MPI_AllReduce_R2MAP, CMF_MPI_AllReduce_R1PTH
#endif
   IMPLICIT NONE
   integer(KIND=JPIM)          :: JF
   real(KIND=JPRB),POINTER     :: D2VEC(:,:) ! point data location to output
   !*** LOCAL
   real(KIND=JPRM)             :: R2OUT(NX,NY)
   real(KIND=JPRM)             :: R1POUT(NPTHOUT,NPTHLEV)
      !================================================
      write(LOGNAM,*) ""
      write(LOGNAM,*) "!---------------------!"

      !*** 0. check date:hour with output frequency
      IF ( MOD(JHOUR,IFRQ_OUT)==0 .and. JMIN==0 ) THEN             ! JHOUR: end of time step , NFPPH: output frequency (hour)

         !*** 1. update IREC & calc average variable
         IRECOUT=IRECOUT+1 
         write(LOGNAM,*) 'CMF::OUTPUT_WRITE: write at time: ', JYYYYMMDD, JHHMM, IRECOUT

         !*** 2. check variable name & allocate data to pointer DVEC
         DO JF=1,NVARSOUT
            SELECT CASE (VAROUT(JF)%CVNAME)
               CASE ('rivsto')
                  D2COPY=P2RIVSTO  !! convert Double to Single precision when using SinglePrecisionMode 
                  D2VEC => D2COPY  !!   (Storage variables are kept as Float64 in SinglePrecisionMode)
               CASE ('fldsto')
                  D2COPY=P2FLDSTO
                  D2VEC => D2COPY
               CASE ('rivout')
                  D2VEC => D2RIVOUT_AVG
               CASE ('rivdph')
                  D2VEC => D2RIVDPH
               CASE ('rivvel')
                  D2VEC => D2RIVVEL_AVG
               CASE ('fldout')
                  D2VEC => D2FLDOUT_AVG
               CASE ('flddph')
                  D2VEC => D2FLDDPH
               CASE ('fldfrc')
                  D2VEC => D2FLDFRC
               CASE ('fldare')
                  D2VEC => D2FLDARE
               CASE ('sfcelv')
                  D2VEC => D2SFCELV
               CASE ('totout')
                  D2VEC => D2OUTFLW_AVG
               CASE ('outflw')            !!  compatibility for previous file name
                  D2VEC => D2OUTFLW_AVG
               CASE ('totsto')
                  D2VEC => D2STORGE
               CASE ('storge')            !!  compatibility for previous file name
                  D2VEC => D2STORGE
               CASE ('pthout')
                  IF( .not. LPTHOUT ) CYCLE
                  D2VEC => D2PTHOUT_AVG
               CASE ('pthflw')
                  IF( .not. LPTHOUT ) CYCLE
               CASE ('maxflw')
                  D2VEC =>  D2OUTFLW_MAX
               CASE ('maxdph')
                  D2VEC =>  D2RIVDPH_MAX
               CASE ('maxsto')
                  D2VEC =>  D2STORGE_MAX
               CASE ('outins')
                  IF( .not. LOUTINS ) CYCLE
                  D2VEC =>  D2OUTINS
               CASE ('gwsto')
                  IF( .not. LGDWDLY ) CYCLE
                  D2COPY=P2GDWSTO
                  D2VEC =>  D2COPY
               CASE ('gdwsto')
                  IF( .not. LGDWDLY ) CYCLE
                  D2COPY=P2GDWSTO
                  D2VEC =>  D2COPY
               CASE ('gwout')
                  IF( .not. LGDWDLY ) CYCLE
                  D2VEC =>  D2GDWRTN_AVG
               CASE ('gdwrtn')
               IF( .not. LGDWDLY ) CYCLE
                  D2VEC =>  D2GDWRTN_AVG
                  CASE ('runoff')             !!  compatibility for previous file name
                  D2VEC =>  D2RUNOFF_AVG  
               CASE ('runoffsub')           !!  compatibility for previous file name
                  IF( .not. LROSPLIT ) CYCLE
                  D2VEC =>  D2ROFSUB_AVG  
               CASE ('rofsfc')
                  D2VEC =>  D2RUNOFF_AVG
               CASE ('rofsub')
                  D2VEC =>  D2ROFSUB_AVG
               CASE ('wevap')
                  IF( .not. LWEVAP ) CYCLE
                  D2VEC => D2WEVAPEX_AVG
               CASE ('winfilt')
                  IF( .not. LWINFILT ) CYCLE
                  D2VEC => D2WINFILTEX_AVG
               CASE ('damsto')   !!! added
                  IF( .not. LDAMOUT ) CYCLE
                  D2COPY=P2DAMSTO
                  D2VEC => D2COPY
               CASE ('daminf')   !!! added
                  IF( .not. LDAMOUT ) CYCLE
                  D2VEC =>  d2daminf_avg
               CASE ('levsto')   !!! added
                  IF( .not. LLEVEE ) CYCLE
                  D2COPY=P2LEVSTO
                  D2VEC => D2COPY
               CASE ('levdph')   !!! added
                  IF( .not. LLEVEE ) CYCLE
                  D2VEC =>  D2LEVDPH
               CASE DEFAULT
            END SELECT   !! variable name select

            IF( KSTEP==0 .and. LOUTINI )THEN  !! write storage only when LOUTINI specified
               IF ( .not. LOUTCDF ) CYCLE
               IF ( VAROUT(JF)%CVNAME/='rivsto' .and. VAROUT(JF)%CVNAME/='fldsto' .and. VAROUT(JF)%CVNAME/='gwsto' ) CYCLE
            ENDIF

            !! convert 1Dvector to 2Dmap
            IF( VAROUT(JF)%CVNAME/='pthflw' ) THEN  !! usual 2D map variable
               CALL vecD2mapR(D2VEC,R2OUT)             !! MPI node data is gathered by vecP2mapR
#ifdef UseMPI_CMF
               CALL CMF_MPI_AllReduce_R2MAP(R2OUT)
#endif
            ELSE
               IF( .not. LPTHOUT ) CYCLE
               R1POUT(:,:)=real(D1PTHFLW_AVG(:,:))
#ifdef UseMPI_CMF
               CALL CMF_MPI_AllReduce_R1PTH(R1POUT)
#endif
            ENDIF

            !*** 3. write D2VEC to output file
            IF ( LOUTCDF ) THEN
               IF ( REGIONTHIS==1 ) CALL WRTE_OUTCDF  !! netCDFG
            ELSE
               IF( VAROUT(JF)%CVNAME=='pthflw' ) THEN
                  IF ( REGIONTHIS==1 ) CALL WRTE_OUTPTH(VAROUT(JF)%BINID,IRECOUT,R1POUT)        !! 1D bifu channel
               ELSE
                  IF( LOUTVEC )THEN
                     CALL WRTE_OUTVEC(VAROUT(JF)%BINID,IRECOUT,D2VEC)         !! 1D vector (optional)
                  ELSE
                     IF ( REGIONTHIS==1 ) CALL WRTE_OUTBIN(VAROUT(JF)%BINID,IRECOUT,R2OUT)         !! 2D map
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

         write(LOGNAM,*) 'CMF::OUTPUT_WRITE: end'


      ENDIF



   !==========================================================
   CONTAINS
      !+ WRTE_OUTBIN
      !+ WRTE_OUTPTH
      !+ WRTE_OUTVEC
      !+ WRTE_OUTCDF
      !==========================================================
      SUBROUTINE WRTE_OUTBIN(IFN,IREC,R2OUTDAT)
      IMPLICIT NONE
      !*** INPUT
      integer(KIND=JPIM),intent(in)   :: IFN                 !! FILE NUMBER
      integer(KIND=JPIM),intent(in)   :: IREC                !! RECORD
      real(KIND=JPRM)                 :: R2OUTDAT(NX,NY)
         !================================================
         write(IFN,REC=IREC) R2OUTDAT

      END SUBROUTINE WRTE_OUTBIN
      !==========================================================
      !+
      !+
      !+
      !==========================================================
      SUBROUTINE WRTE_OUTPTH(IFN,IREC,R2OUTDAT)
      IMPLICIT NONE
      !*** INPUT
      integer(KIND=JPIM),intent(in)   :: IFN                 !! FILE NUMBER
      integer(KIND=JPIM),intent(in)   :: IREC                !! RECORD
      real(KIND=JPRM)                 :: R2OUTDAT(NPTHOUT,NPTHLEV)
      !================================================
         write(IFN,REC=IREC) R2OUTDAT

      END SUBROUTINE WRTE_OUTPTH
      !==========================================================
      !+
      !+
      !+
      !==========================================================
      SUBROUTINE WRTE_OUTVEC(IFN,IREC,D2OUTDAT)
      IMPLICIT NONE
      !*** INPUT
      integer(KIND=JPIM),intent(in)   :: IFN                 !! FILE NUMBER
      integer(KIND=JPIM),intent(in)   :: IREC                !! RECORD
      real(KIND=JPRB),intent(in)      :: D2OUTDAT(NSEQMAX,1) !! OUTPUT DATA
      !*** LOCAL
      real(KIND=JPRM)                 :: R2OUTDAT(NSEQMAX,1)
         !================================================
         R2OUTDAT(:,:)=real(D2OUTDAT(:,:))
         write(IFN,REC=IREC) R2OUTDAT

      END SUBROUTINE WRTE_OUTVEC
      !==========================================================
      !+
      !+
      !+
      !==========================================================
      SUBROUTINE WRTE_OUTCDF
#ifdef UseCDF_CMF
      USE NETCDF 
      USE YOS_CMF_TIME,            only: KMINSTART,KMINNEXT
      USE CMF_UTILS_MOD,           only: NCERROR
      IMPLICIT NONE
      real(KIND=JPRB)                 :: XTIME ! seconds since start of the run !

         !================================================
         XTIME=real( (KMINNEXT-KMINSTART),JPRB) *60._JPRB      !! for netCDF
         CALL NCERROR( NF90_PUT_VAR(VAROUT(JF)%NCID,VAROUT(JF)%TIMID,XTIME,(/VAROUT(JF)%IRECNC/)) )

         CALL NCERROR( NF90_PUT_VAR(VAROUT(JF)%NCID,VAROUT(JF)%VARID,R2OUT(1:NX,1:NY),(/1,1,VAROUT(JF)%IRECNC/),(/NX,NY,1/)) )

         ! update IREC
         VAROUT(JF)%IRECNC=VAROUT(JF)%IRECNC+1

         ! Comment out this as it slows down significantly the writting in the cray  
         !CALL NCERROR( NF90_SYNC(VAROUT(JF)%NCID) )  
#endif
      END SUBROUTINE WRTE_OUTCDF 
   !==========================================================

   END SUBROUTINE CMF_OUTPUT_WRITE
   !####################################################################





   !####################################################################
   SUBROUTINE CMF_OUTPUT_END
   ! Finalize output module (close files)
   ! -- Called from CMF_DRV_END
#ifdef UseCDF_CMF
   USE NETCDF
   USE CMF_UTILS_MOD,           only: NCERROR
#endif
   USE YOS_CMF_MAP,             only: REGIONTHIS
   IMPLICIT NONE
   ! Local variables
   integer(KIND=JPIM)              :: JF
      !================================================
      write(LOGNAM,*) ""
      write(LOGNAM,*) "!---------------------!"
      write(LOGNAM,*) "CMF::OUTPUT_END: finalize output module"

      IF( REGIONTHIS==1 )THEN
         IF (LOUTCDF) THEN
#ifdef UseCDF_CMF
            DO JF=1,NVARSOUT
               CALL NCERROR( NF90_CLOSE(VAROUT(JF)%NCID))
               write(LOGNAM,*) "Output netcdf output unit closed:",VAROUT(JF)%NCID
            ENDDO
#endif
         ELSE !! binary output
            DO JF=1,NVARSOUT
               close(VAROUT(JF)%BINID)
               write(LOGNAM,*) "Output binary output unit closed:",VAROUT(JF)%BINID
            ENDDO
            IF( LOUTVEC )THEN
               CALL WRTE_mapR2vecD  !! write map-vector conversion file
            ENDIF
         ENDIF
      ENDIF

      write(LOGNAM,*) "CMF::OUTPUT_END: end"


   CONTAINS
      !==========================================================
      !+ WRTE_mapR2vecD
      !+
      !+
      !==========================================================
      SUBROUTINE WRTE_mapR2vecD       !! 1D sequence vector informtion required to convert MPI distributed vector output to 2D map
      USE YOS_CMF_INPUT,      only: TMPNAM
      USE YOS_CMF_MAP,        only: I1SEQX, I1SEQY, NSEQMAX
      USE CMF_UTILS_MOD,      only: INQUIRE_FID
      IMPLICIT NONE
      !* local variable
      character(LEN=256)         :: CFILE1
      !================================================
         IF( LOUTVEC )THEN
            CFILE1='./ind_xy'//TRIM(CSUFVEC)

            write(LOGNAM,*) "LOUTVEC: write mapR2vecD conversion table", TRIM(CFILE1)

            TMPNAM=INQUIRE_FID()
            open(TMPNAM,FILE=CFILE1,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NSEQMAX)
            write(TMPNAM,REC=1) I1SEQX
            write(TMPNAM,REC=2) I1SEQY
            close(TMPNAM)
         ENDIF

      END SUBROUTINE WRTE_mapR2vecD
      !================================================


   END SUBROUTINE CMF_OUTPUT_END
      !####################################################################



   !####################################################################
   SUBROUTINE CMF_OUTTXT_WRTE
   USE YOS_CMF_DIAG,       only: D2OUTFLW
   USE YOS_CMF_TIME,       only: IYYYYMMDD,ISYYYY
   USE YOS_CMF_MAP,        only: I2VECTOR
   USE CMF_UTILS_MOD,      only: INQUIRE_FID

   ! local
   integer(KIND=JPIM)                  :: GID, GIX, GIY, GISEQ
   character(len=256),SAVE             :: GNAME

   integer(KIND=JPIM),SAVE             :: IGAUGE, NGAUGE, NGAUGEX
   integer(KIND=JPIM),ALLOCATABLE,SAVE :: WriteID(:), WriteISEQ(:)
   character(len=9),ALLOCATABLE,SAVE   :: WriteName(:)
   real(KIND=JPRB),ALLOCATABLE,SAVE    :: WriteOut(:)

   ! File IO
   integer(KIND=JPIM),SAVE             :: LOGOUTTXT
   character(len=4),SAVE               :: cYYYY
   character(len=256),SAVE             :: CLEN, CFMT
   character(len=256),SAVE             :: COUTTXT
   logical,SAVE                        :: IsOpen
   DATA IsOpen       /.FALSE./
      ! ======
      IF( LOUTTXT )THEN

         IF( .not. IsOpen)THEN
            IsOpen=.TRUE.

            NGAUGEX=0
            LOGOUTTXT=INQUIRE_FID()
            open(LOGOUTTXT,FILE=CGAUTXT,FORM='formatted',STATUS='old')
            read(LOGOUTTXT,*) NGAUGE
            DO IGAUGE=1, NGAUGE
               read(LOGOUTTXT,*) GID, GNAME, GIX, GIY
               IF( I2VECTOR(GIX,GIY)>0 )THEN
                  NGAUGEX=NGAUGEX+1
               ENDIF
            ENDDO
            close(LOGOUTTXT)

            allocate( WriteID(NGAUGEX),WriteISEQ(NGAUGEX),WriteOut(NGAUGEX),WriteName(NGAUGEX))

            NGAUGEX=0
            open(LOGOUTTXT,FILE=CGAUTXT,FORM='formatted',STATUS='old')
            read(LOGOUTTXT,*) NGAUGE
            DO IGAUGE=1, NGAUGE
               read(LOGOUTTXT,*) GID, GNAME, GIX, GIY
               IF( I2VECTOR(GIX,GIY)>0 )THEN
                  NGAUGEX=NGAUGEX+1
                  WriteID(NGAUGEX)  =GID
                  WriteName(NGAUGEX)=TRIM(GNAME)
                  WriteISEQ(NGAUGEX)=I2VECTOR(GIX,GIY)
               ENDIF
            ENDDO
            close(LOGOUTTXT)

            ! ============
            write(CYYYY,'(i4.4)') ISYYYY
            COUTTXT='./outtxt-'//TRIM(cYYYY)//'.txt'

            LOGOUTTXT=INQUIRE_FID()
            open(LOGOUTTXT,FILE=COUTTXT,FORM='formatted')

            write(CLEN,'(i0)') NGAUGE
            CFMT="(i10,"//TRIM(CLEN)//"(i10))"
            write(LOGOUTTXT,CFMT) NGAUGEX, ( WriteID(IGAUGE),IGAUGE=1,NGAUGEX )

            CFMT="(i10,"//TRIM(CLEN)//"(x,a9))"
            write(LOGOUTTXT,CFMT) NGAUGEX, ( WriteName(IGAUGE),IGAUGE=1,NGAUGEX )


            CFMT="(i10,"//TRIM(CLEN)//"(f10.2))"
         ENDIF

         DO IGAUGE=1, NGAUGEX
            GISEQ=WriteISEQ(IGAUGE)
            WriteOut(IGAUGE) = D2OUTFLW(GISEQ,1)
         ENDDO
            
         write(LOGOUTTXT,CFMT) IYYYYMMDD, ( WriteOUT(IGAUGE),IGAUGE=1,NGAUGEX )

      ENDIF

   END SUBROUTINE CMF_OUTTXT_WRTE
      !####################################################################

END MODULE CMF_CTRL_OUTPUT_MOD
