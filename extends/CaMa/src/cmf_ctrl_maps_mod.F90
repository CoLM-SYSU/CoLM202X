MODULE CMF_CTRL_MAPS_MOD
!==========================================================
!* PURPOSE: Control CaMa-Flood map/topography data
!
!* CONTAINS:
! -- CMF_MAPS_NMLIST   : configuration from namelist
! -- CMF_RIVMAP_INIT  : read & set river network map 
! -- CMF_TOPO_INIT    : read & set topography
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
   ! shared variables in module
   USE PARKIND1,                ONLY: JPIM, JPRB, JPRD, JPRM
   USE YOS_CMF_INPUT,           only: LOGNAM
   IMPLICIT NONE
   SAVE
   !*** NAMELIST/NMAP/ from inputnam
   character(LEN=256)              :: CNEXTXY         !! river network nextxy
   character(LEN=256)              :: CGRAREA         !! catchment area
   character(LEN=256)              :: CELEVTN         !! bank top elevation
   character(LEN=256)              :: CNXTDST         !! distance to next outlet
   character(LEN=256)              :: CRIVLEN         !! river channel length
   character(LEN=256)              :: CFLDHGT         !! floodplain elevation profile
   !* river channel parameters
   character(LEN=256)              :: CRIVWTH         !! channel width
   character(LEN=256)              :: CRIVHGT         !! channel depth
   character(LEN=256)              :: CRIVMAN         !! river manning coefficient
   !* optional maps
   character(LEN=256)              :: CPTHOUT         !! bifurcation channel table
   character(LEN=256)              :: CGDWDLY         !! Groundwater Delay Parameter
   character(LEN=256)              :: CMEANSL         !! mean sea level
   !* MPI parallelization
   character(LEN=256)              :: CMPIREG         !! MPI region map
   !* netCDF map
   logical                         :: LMAPCDF         !! true for netCDF map input
   character(LEN=256)              :: CRIVCLINC       !! river map netcdf
   character(LEN=256)              :: CRIVPARNC       !! river parameter netcdf (WIDTH,HEIGHT, Manning, ground wateer delay)
   character(LEN=256)              :: CMEANSLNC       !! mean sea level netCDF
   character(LEN=256)              :: CMPIREGNC       !! MPI region map in netcdf

   NAMELIST/NMAP/     CNEXTXY,  CGRAREA,  CELEVTN,  CNXTDST, CRIVLEN, CFLDHGT, &
                     CRIVWTH,  CRIVHGT,  CRIVMAN,  CPTHOUT, CGDWDLY, CMEANSL, &
                     CMPIREG,  LMAPCDF,  CRIVCLINC,CRIVPARNC,CMEANSLNC,CMPIREGNC


CONTAINS
!####################################################################
! -- CMF_MAP_NMLIST   : configuration from namelist
! -- CMF_RIVMAP_INIT  : read & set river network map 
! -- CMF_TOPO_INIT    : read & set topography
!
!
!####################################################################
   SUBROUTINE CMF_MAPS_NMLIST
   ! reed setting from namelist
   ! -- Called from CMF_DRV_NMLIST
   USE YOS_CMF_INPUT,      only: CSETFILE,NSETFILE,LMEANSL,LGDWDLY
   USE CMF_UTILS_MOD,      only: INQUIRE_FID
   IMPLICIT NONE
      !================================================
      !*** 1. open namelist
      write(LOGNAM,*) ""
      write(LOGNAM,*) "!---------------------!"

      NSETFILE=INQUIRE_FID()
      open(NSETFILE,FILE=CSETFILE,STATUS="OLD")
      write(LOGNAM,*) "CMF::MAP_NMLIST: namelist open in unit: ", TRIM(CSETFILE), NSETFILE 

      !*** 2. default value
      CNEXTXY="./nextxy.bin"
      CGRAREA="./ctmare.bin"
      CELEVTN="./elevtn.bin"
      CNXTDST="./nxtdst.bin"
      CRIVLEN="./rivlen.bin"
      CFLDHGT="./fldhgt.bin"

      CRIVWTH="./rivwth.bin"
      CRIVHGT="./rivhgt.bin"
      CRIVMAN="./rivman.bin"

      CPTHOUT="./bifprm.txt"
      CGDWDLY="NONE"
      CMEANSL="NONE"

      CMPIREG="NONE"

      LMAPCDF=.FALSE.
      CRIVCLINC="NONE"
      CRIVPARNC="NONE"
      CMEANSLNC="NONE"
      CMPIREGNC="NONE"

      !*** 3. read namelist
      rewind(NSETFILE)
      read(NSETFILE,NML=NMAP)

      write(LOGNAM,*)     "=== NAMELIST, NMAP ==="
      write(LOGNAM,*)     "LMAPCDF:   ", LMAPCDF
      IF( LMAPCDF )THEN
         write(LOGNAM,*)   "CRIVCLINC: ", TRIM(CRIVCLINC)
         write(LOGNAM,*)   "CRIVPARNC: ", TRIM(CRIVPARNC)
         IF( LMEANSL ) THEN
            write(LOGNAM,*) "CMEANSLNC: ", TRIM(CMEANSLNC)
         ENDIF
#ifdef UseMPI_CMF
         write(LOGNAM,*) "CMPIREGNC:   ", TRIM(CMPIREGNC)
#endif
      ELSE
         write(LOGNAM,*)   "CNEXTXY:   ", TRIM(CNEXTXY)
         write(LOGNAM,*)   "CGRAREA:   ", TRIM(CGRAREA)
         write(LOGNAM,*)   "CELEVTN:   ", TRIM(CELEVTN)
         write(LOGNAM,*)   "CNXTDST:   ", TRIM(CNXTDST)
         write(LOGNAM,*)   "CRIVLEN:   ", TRIM(CRIVLEN)
         write(LOGNAM,*)   "CFLDHGT:   ", TRIM(CFLDHGT)

         write(LOGNAM,*)   "CRIVWTH:   ", TRIM(CRIVWTH)
         write(LOGNAM,*)   "CRIVHGT:   ", TRIM(CRIVHGT)
         write(LOGNAM,*)   "CRIVMAN:   ", TRIM(CRIVMAN)

         write(LOGNAM,*)   "CPTHOUT:   ", TRIM(CPTHOUT)
         IF( LGDWDLY )THEN
            write(LOGNAM,*) "CGDWDLY:    ",TRIM(CGDWDLY)
         ENDIF
         IF( LMEANSL )THEN
            write(LOGNAM,*) "CMEANSL:   ", TRIM(CMEANSL)
         ENDIF
#ifdef UseMPI_CMF
         write(LOGNAM,*) "CMPIREG:   ", TRIM(CMPIREG)
#endif
      ENDIF

      close(NSETFILE)

      write(LOGNAM,*) "CMF::MAP_NMLIST: end"

   END SUBROUTINE CMF_MAPS_NMLIST
!####################################################################





!####################################################################
   SUBROUTINE CMF_RIVMAP_INIT
   ! read & set river network map 
   ! -- CALL from CMF_DRV_INIT
   USE YOS_CMF_INPUT,      only: TMPNAM, NX,NY,NLFP, LPTHOUT
   USE YOS_CMF_MAP,        only: I2NEXTX,I2NEXTY, I2REGION, REGIONALL,REGIONTHIS, &
                              & I1SEQX, I1SEQY,  I1NEXT,  I2VECTOR, D1LON,    D1LAT,      &
                              & NSEQRIV,  NSEQALL,  NSEQMAX
   USE CMF_UTILS_MOD,      only: INQUIRE_FID
   IMPLICIT NONE
      !================================================
      write(LOGNAM,*) ""
      write(LOGNAM,*) "!---------------------!"
      write(LOGNAM,*) 'CMF::RIVMAP_INIT: river network initialization'

      ! *** 1. allocate ARRAYS
      allocate( I2NEXTX(NX,NY) )
      allocate( I2NEXTY(NX,NY) )
      allocate( I2REGION(NX,NY) )
      allocate( D1LON(NX) )
      allocate( D1LAT(NY) )

      !============================
      !*** 2a. read river network map
      write(LOGNAM,*) 'CMF::RIVMAP_INIT: read nextXY & set lat lon'
      IF( LMAPCDF )THEN
#ifdef UseCDF_CMF
         CALL READ_MAP_CDF
#endif
      ELSE
         CALL READ_MAP_BIN
      ENDIF

      !*** 2b. calculate river sequence & regions
      write(LOGNAM,*) 'CMF::RIVMAP_INIT: calc region'
      CALL CALC_REGION

      !============================
      !*** 3. conversion 2D map -> 1D vector
      write(LOGNAM,*) 'CMF::RIVMAP_INIT: calculate 1d river sequence'

      CALL CALC_1D_SEQ                                  !! 2D map to 1D vector conversion. for faster calculation

      write(LOGNAM,*) '  NSEQRIV=',NSEQRIV
      write(LOGNAM,*) '  NSEQALL=',NSEQALL

      !*** 3c. Write Map Data                                       !! used for combining mpi distributed output into one map
      IF( REGIONTHIS==1 )THEN
         TMPNAM=INQUIRE_FID()
         open(TMPNAM,FILE='./mapdata.txt',FORM='FORMATTED')
         write(TMPNAM,*) 'NX',        NX
         write(TMPNAM,*) 'NY',        NY
         write(TMPNAM,*) 'NLFP',      NLFP
         write(TMPNAM,*) 'REGIONALL', REGIONALL
         write(TMPNAM,*) 'NSEQMAX',   NSEQMAX
         close(TMPNAM)
      ENDIF

      !============================
      !*** 4.  bifurcation channel parameters
      IF( LPTHOUT )THEN
         write(LOGNAM,*) 'CMF::RIVMAP_INIT: read bifurcation channel setting'
         CALL READ_BIFPARAM
      ENDIF

      deallocate(I2NEXTY,I2REGION )

      write(LOGNAM,*) 'CMF::RIVMAP_INIT: end'

   CONTAINS
      !==========================================================
      !+ READ_MAP_BIN
      !+ READ_MAP_CDF
      !+ CALC_REGION
      !+ CALC_1D_SEQ
      !+ READ_BIFPRM
      !==========================================================
      SUBROUTINE READ_MAP_BIN
      USE YOS_CMF_INPUT,      only: TMPNAM, LMAPEND
      USE YOS_CMF_INPUT,      only: WEST,EAST,NORTH,SOUTH
      USE CMF_UTILS_MOD,      only: INQUIRE_FID, CONV_ENDI
      IMPLICIT NONE
      !* local variables
      integer(KIND=JPIM),SAVE    :: IX,IY
         !==========================================================
         !*** read river map
         write(LOGNAM,*)'RIVMAP_INIT: nextxy binary: ',TRIM(CNEXTXY)
         TMPNAM=INQUIRE_FID()
         open(TMPNAM,FILE=CNEXTXY,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
         read(TMPNAM,REC=1) I2NEXTX
         read(TMPNAM,REC=2) I2NEXTY
         close(TMPNAM)

         IF ( LMAPEND )THEN
            CALL CONV_ENDI(I2NEXTX,NX,NY)
            CALL CONV_ENDI(I2NEXTY,NX,NY)
         ENDIF

         !*** calculate lat, lon
         IF( WEST>=-180._JPRB .and. EAST<=360._JPRB .and. SOUTH>=-180._JPRB .and. NORTH<=180._JPRB )THEN  !! bugfix_v396a
   !$OMP PARALLEL DO
            DO IX=1,NX
               D1LON(IX)=WEST +(DBLE(IX)-0.5D0)*(EAST-WEST)  /DBLE(NX)
            ENDDO
   !$OMP END PARALLEL DO
   !$OMP PARALLEL DO
            DO IY=1,NY
               D1LAT(IY)=NORTH-(DBLE(IY)-0.5D0)*(NORTH-SOUTH)/DBLE(NY)
            ENDDO
   !$OMP END PARALLEL DO
         ENDIF

      END SUBROUTINE READ_MAP_BIN
      !==========================================================
      !+
      !+
      !+
      !==========================================================
#ifdef UseCDF_CMF
      SUBROUTINE READ_MAP_CDF
      USE CMF_UTILS_MOD  ,only: NCERROR
      USE NETCDF
      IMPLICIT NONE
      !* local variables
      integer(KIND=JPIM)              :: NCID,VARID
         !================================================
         write(LOGNAM,*)'RIVMAP_INIT: nextxy netCDF: ', TRIM(CRIVCLINC)

         CALL NCERROR (NF90_OPEN(CRIVCLINC,NF90_NOWRITE,NCID),'opening '//TRIM(CRIVCLINC) )

         !*** next xy
         CALL NCERROR ( NF90_INQ_VARID(NCID,'nextx',VARID),'getting id' )
         CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2NEXTX),'reading data' ) 

         CALL NCERROR ( NF90_INQ_VARID(NCID,'nexty',VARID),'getting id' )
         CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2NEXTY),'reading data' )

         !*** lat, lon
         CALL NCERROR ( NF90_INQ_VARID(NCID,'lat',VARID),'getting id' )
         CALL NCERROR ( NF90_GET_VAR(NCID,VARID,D1LAT),'reading data' )

         CALL NCERROR ( NF90_INQ_VARID(NCID,'lon',VARID),'getting id' )
         CALL NCERROR ( NF90_GET_VAR(NCID,VARID,D1LON),'reading data' )

         CALL NCERROR( NF90_CLOSE(NCID))

      END SUBROUTINE READ_MAP_CDF
#endif
   !==========================================================
   !+
   !+
   !+
   !==========================================================
      SUBROUTINE CALC_REGION    !! evenly allocate pixels to mpi nodes (updated in v4.03. MPI region given from file)
      USE YOS_CMF_INPUT,           only: IMIS
#ifdef UseCDF_CMF
      USE CMF_UTILS_MOD,           only: NCERROR
      USE NETCDF
#endif
      IMPLICIT NONE
      !* local variables
      integer(KIND=JPIM),ALLOCATABLE  :: REGIONGRID(:)
      !
      integer(KIND=JPIM),SAVE         :: IX,IY
      integer(KIND=JPIM),SAVE         :: IREGION
#ifdef UseMPI_CMF
#ifdef UseCDF_CMF
INTEGER(KIND=JPIM)              :: NCID
INTEGER(KIND=JPIM)              :: VARID
#endif
#endif
!$OMP THREADPRIVATE               (IX)
      !================================================
         write(LOGNAM,*) 'RIVMAP_INIT: region code'

         !*** read MPI region map
         REGIONALL=1
         I2REGION(:,:)=IMIS
!$OMP PARALLEL DO
         DO IY=1, NY
            DO IX=1, NX
               IF( I2NEXTX(IX,IY)/=IMIS ) THEN
                  I2REGION(IX,IY)=1
               ENDIF
            ENDDO
         ENDDO
!$OMP END PARALLEL DO

   !! Use MPI: read MPI region map, allocate regions to MPI nodes
#ifdef UseMPI_CMF
         IF ( LMAPCDF ) THEN
#ifdef UseCDF_CMF
            CALL NCERROR (NF90_OPEN(CMPIREGNC,NF90_NOWRITE,NCID),'opening '//TRIM(CMPIREGNC) )
            CALL NCERROR (NF90_INQ_VARID(NCID, 'mpireg',VARID),'getting id' )
            CALL NCERROR (NF90_GET_VAR(NCID,VARID,I2REGION),'reading data' )
            CALL NCERROR (NF90_CLOSE(NCID))
#endif
         ELSE
            write(LOGNAM,*)'RIVMAP_INIT: read MPI region: ',TRIM(CNEXTXY)
            TMPNAM=INQUIRE_FID()
            open(TMPNAM,FILE=CMPIREG,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
            read(TMPNAM,REC=1) I2REGION
            close(TMPNAM)
         ENDIF

         REGIONALL=1
!$OMP PARALLEL DO REDUCTION(max:REGIONALL)
         DO IY=1, NY
            DO IX=1, NX
               REGIONALL=MAX( REGIONALL, I2REGION(IX,IY) )
            ENDDO
         ENDDO
!$OMP END PARALLEL DO
#endif


         write(LOGNAM,*)'RIVMAP_INIT: count number of grid in each region: '
         allocate(REGIONGRID(REGIONALL))
         REGIONGRID(:)=0
         !! OMP reduction operation for array might not be available in some environment
         DO IY=1, NY
            DO IX=1, NX
               IF( I2REGION(IX,IY)>0 ) THEN
                  IREGION=I2REGION(IX,IY)
                  REGIONGRID(IREGION)=REGIONGRID(IREGION)+1
               ENDIF
            ENDDO
         ENDDO

         NSEQMAX=0
         DO IREGION=1, REGIONALL
            NSEQMAX=MAX(NSEQMAX,REGIONGRID(IREGION))  !! maximum nseqall among all MPI region
         ENDDO

         write(LOGNAM,*) 'CALC_REGION: REGIONALL= ', REGIONALL
         write(LOGNAM,*) 'CALC_REGION: NSEQMAX='   , NSEQMAX
         WRITE(LOGNAM,*) 'CALC_REGION: NSEQALL='   , NSEQALL
      END SUBROUTINE CALC_REGION
      !==========================================================
      !+
      !+
      !+
      !==========================================================
      SUBROUTINE CALC_1D_SEQ
      ! OpenMP is not used, because results of this subroutine highly depents on calculation order
      USE YOS_CMF_INPUT,           only: IMIS
      IMPLICIT NONE
      !* local variables
      integer(KIND=JPIM)              :: IX,IY,JX,JY,ISEQ,JSEQ,ISEQ1,ISEQ2,AGAIN

      integer(KIND=JPIM),ALLOCATABLE  :: NUPST(:,:), UPNOW(:,:)
         !================================================
         write(LOGNAM,*) 'RIVMAP_INIT: convert 2D map to 1D sequence'

         allocate( NUPST(NX,NY) )
         allocate( UPNOW(NX,NY) )

         allocate( I1SEQX(NSEQMAX) )
         allocate( I1SEQY(NSEQMAX) )
         allocate( I1NEXT(NSEQMAX) )
         allocate( I2VECTOR(NX,NY) )
         I1SEQX(:)=0
         I1SEQY(:)=0
         I1NEXT(:)=0
         I2VECTOR(:,:)=0

         ! count number of upstream 
         NUPST(:,:)=0
         UPNOW(:,:)=0
         DO IY=1, NY
            DO IX=1, NX
               IF( I2NEXTX(IX,IY).gt.0 .and. I2REGION(IX,IY)==REGIONTHIS )THEN
                  JX=I2NEXTX(IX,IY)
                  JY=I2NEXTY(IX,IY)
                  NUPST(JX,JY)=NUPST(JX,JY)+1
               ENDIF
            ENDDO
         ENDDO

         ! register upmost grid in 1d sequence
         ISEQ=0
         DO IY=1, NY
            DO IX=1, NX
               IF( I2NEXTX(IX,IY).gt.0 .and. I2REGION(IX,IY)==REGIONTHIS )THEN
                  IF( NUPST(IX,IY)==UPNOW(IX,IY) )THEN
                  ISEQ=ISEQ+1
                  I1SEQX(ISEQ)=IX
                  I1SEQY(ISEQ)=IY
                  I2VECTOR(IX,IY)=ISEQ
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         ISEQ1=1
         ISEQ2=ISEQ

         AGAIN=1
         DO WHILE( AGAIN==1 )
            AGAIN=0
            JSEQ=ISEQ2
            DO ISEQ=ISEQ1, ISEQ2
               IX=I1SEQX(ISEQ)
               IY=I1SEQY(ISEQ)
               JX=I2NEXTX(IX,IY)
               JY=I2NEXTY(IX,IY)
               UPNOW(JX,JY)=UPNOW(JX,JY)+1
               IF( UPNOW(JX,JY)==NUPST(JX,JY) .and. I2NEXTX(JX,JY)>0 )THEN !! if all upstream calculated, register to 1D sequence
                  JSEQ=JSEQ+1
                  I1SEQX(JSEQ)=JX
                  I1SEQY(JSEQ)=JY
                  I2VECTOR(JX,JY)=JSEQ
                  AGAIN=1
               ENDIF
            ENDDO
            ISEQ1=ISEQ2+1
            ISEQ2=JSEQ
         ENDDO
         NSEQRIV=JSEQ

         ISEQ=NSEQRIV
         DO IY=1, NY
            DO IX=1, NX
               IF( I2NEXTX(IX,IY).lt.0 .and. I2NEXTX(IX,IY).NE.IMIS .and. I2REGION(IX,IY)==REGIONTHIS )THEN
                  ISEQ=ISEQ+1
                  I1SEQX(ISEQ)=IX
                  I1SEQY(ISEQ)=IY
                  I2VECTOR(IX,IY)=ISEQ
               ENDIF
            ENDDO
         ENDDO
         NSEQALL=ISEQ

         DO ISEQ=1, NSEQALL
            IX=I1SEQX(ISEQ)
            IY=I1SEQY(ISEQ)
            IF( I2NEXTX(IX,IY)>0 )THEN
               JX=I2NEXTX(IX,IY)
               JY=I2NEXTY(IX,IY)
               I1NEXT(ISEQ)=I2VECTOR(JX,JY)
            ELSE
               I1NEXT(ISEQ)=I2NEXTX(IX,IY)
            ENDIF
         ENDDO

         deallocate(NUPST,UPNOW)
            
      END SUBROUTINE CALC_1D_SEQ
      !==========================================================
      !+
      !+
      !+
      !==========================================================
      SUBROUTINE READ_BIFPARAM    !! evenly allocate pixels to mpi nodes (not used in vcurrent version)
      USE YOS_CMF_INPUT,      only: PMANRIV, PMANFLD
      USE YOS_CMF_MAP,        only: NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN,&
                                 & PTH_DST, PTH_ELV, PTH_WTH,  PTH_MAN
      USE CMF_UTILS_MOD,      only: INQUIRE_FID
      IMPLICIT NONE
      !* local variables
      integer(KIND=JPIM)         :: IX,IY, JX,JY
      integer(KIND=JPIM)         :: IPTH,  ILEV,  NPTHOUT1
      real(KIND=JPRB)            :: PELV,  PWTH,  PDPH
         !================================================
         write(LOGNAM,*)"RIVMAP_INIT: Bifuraction channel:", TRIM(CPTHOUT)

         TMPNAM=INQUIRE_FID()
         open(TMPNAM,FILE=CPTHOUT,FORM='FORMATTED')
         read(TMPNAM,*) NPTHOUT,NPTHLEV

         write(LOGNAM,*) "Bifurcation channel dimantion", NPTHOUT, NPTHLEV

         allocate( PTH_UPST(NPTHOUT) )
         allocate( PTH_DOWN(NPTHOUT) )
         allocate( PTH_DST(NPTHOUT)  )
         allocate( PTH_ELV(NPTHOUT,NPTHLEV) )
         allocate( PTH_WTH(NPTHOUT,NPTHLEV) )
         allocate( PTH_MAN(NPTHLEV)  )

         NPTHOUT1=0
         DO IPTH=1, NPTHOUT
            READ(TMPNAM,*) IX, IY, JX, JY, PTH_DST(IPTH), PELV, PDPH, (PTH_WTH(IPTH,ILEV),ILEV=1,NPTHLEV)
            PTH_UPST(IPTH)=I2VECTOR(IX,IY)
            PTH_DOWN(IPTH)=I2VECTOR(JX,JY)
            IF (PTH_UPST(IPTH) > 0 .and. PTH_DOWN(IPTH) > 0) THEN
               NPTHOUT1=NPTHOUT1+1
            ENDIF
            DO ILEV=1, NPTHLEV
               IF( ILEV==1 )THEN            !!ILEV=1: water channel bifurcation. consider bifurcation channel depth
                  PWTH=PTH_WTH(IPTH,ILEV)
                  IF( PWTH>0 )THEN
                     PTH_ELV(IPTH,ILEV)=PELV - PDPH
                  ELSE
                     PTH_ELV(IPTH,ILEV)=1.E20
                  ENDIF
               ELSE
                  PWTH=PTH_WTH(IPTH,ILEV)
                  IF( PWTH>0 )THEN
                     PTH_ELV(IPTH,ILEV)=PELV + ILEV - 2.0    !! ILEV=2: bank top level 
                  ELSE
                     PTH_ELV(IPTH,ILEV)=1.E20
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         close(TMPNAM)

         DO ILEV=1, NPTHLEV
            IF( ILEV==1 )THEN
               PTH_MAN(ILEV)=PMANRIV
            ELSE
               PTH_MAN(ILEV)=PMANFLD
            ENDIF
         ENDDO

         IF (NPTHOUT /= NPTHOUT1) THEN
            write(LOGNAM,*)"Bifuraction channel outside of domain. Only valid:", NPTHOUT1
         ENDIF

      END SUBROUTINE READ_BIFPARAM
      !==========================================================

   END SUBROUTINE CMF_RIVMAP_INIT
!####################################################################





   !####################################################################
   SUBROUTINE CMF_TOPO_INIT
   ! read & set topography map 
   ! -- call from CMF_DRV_INIT
   USE YOS_CMF_INPUT,  only: TMPNAM,   NX,NY,NLFP, LMAPEND,  &
                           & LFPLAIN,  LMEANSL,  LGDWDLY,  LSLPMIX, LSLOPEMOUTH
   USE YOS_CMF_MAP,    only: D2NXTDST, D2GRAREA, D2ELEVTN, D2RIVLEN, &
                           & D2RIVWTH, D2RIVHGT, D2FLDHGT, D2RIVELV, &
                           & D2FLDGRD, D2RIVMAN, P2RIVSTOMAX, P2FLDSTOMAX,  &
                           & DFRCINC,  NSEQALL,  NSEQMAX, D2MEANSL, D2DWNELV, &
                           & D2GDWDLY, I2MASK
   IMPLICIT NONE
   !================================================
      write(LOGNAM,*) ""
      write(LOGNAM,*) "!---------------------!"

      write(LOGNAM,*) 'CMF::TOPO_INIT: topography map initialization'

      ! *** 1. allocate ARRAYS
      allocate( D2GRAREA(NSEQMAX,1) )
      allocate( D2ELEVTN(NSEQMAX,1) )
      allocate( D2NXTDST(NSEQMAX,1) )
      allocate( D2RIVLEN(NSEQMAX,1) )
      allocate( D2RIVWTH(NSEQMAX,1) )
      allocate( D2RIVHGT(NSEQMAX,1) )
      allocate( D2FLDHGT(NSEQMAX,1,NLFP) )
      allocate( D2RIVMAN(NSEQMAX,1) )
      allocate( D2MEANSL(NSEQMAX,1) )
      allocate( D2DWNELV(NSEQMAX,1) )
      allocate( D2GDWDLY(NSEQMAX,1) )
      allocate( I2MASK(NSEQMAX,1) )

      D2GRAREA(:,:)  =0._JPRB
      D2ELEVTN(:,:)  =0._JPRB
      D2NXTDST(:,:)  =0._JPRB
      D2RIVLEN(:,:)  =0._JPRB
      D2RIVWTH(:,:)  =0._JPRB
      D2RIVHGT(:,:)  =0._JPRB
      D2FLDHGT(:,:,:)=0._JPRB
      D2RIVMAN(:,:)  =0._JPRB
      D2MEANSL(:,:)  =0._JPRB
      D2DWNELV(:,:)  =0._JPRB
      D2GDWDLY(:,:)  =0._JPRB
      I2MASK(:,:)    =0._JPIM     !! mask for calculation (IFS slopemix: Kinemacti Wave for Mask=1; Reservoir: dam=2, dam upstream=1)

      !============================
      ! *** 2. Read topo map
      write(LOGNAM,*) 'CMF::TOPO_INIT: read topography maps'
      IF ( .not. LMAPCDF ) THEN
         CALL READ_TOPO_BIN
      ELSE
         CALL READ_TOPO_CDF
      ENDIF

      !============================
      ! *** 3a. Calc Channel Parameters
      write(LOGNAM,*) 'TOPO_INIT: calc river channel parameters'

      ALLOCATE(P2RIVSTOMAX(NSEQMAX,1))
      allocate(D2RIVELV(NSEQMAX,1))

      IF ( LFPLAIN ) THEN
         P2RIVSTOMAX(:,:) = D2RIVLEN(:,:) * D2RIVWTH(:,:) * D2RIVHGT(:,:)
      ELSE
         write(LOGNAM,*) 'TOPO_INIT: no floodplain (rivstomax=1.D18)'
         P2RIVSTOMAX(:,:) = 1.E18
      ENDIF
      D2RIVELV(:,:) = D2ELEVTN(:,:) - D2RIVHGT(:,:)

      !*** 3b. Calc Channel Parameters
      write(LOGNAM,*) 'TOPO_INIT: calc floodplain parameters'

      ALLOCATE(P2FLDSTOMAX(NSEQMAX,1,NLFP))
      allocate(D2FLDGRD(NSEQMAX,1,NLFP))
      CALL SET_FLDSTG

      !*** 3c. Calc downstream boundary
      write(LOGNAM,*) 'TOPO_INIT: calc downstream boundary elevation'
      D2DWNELV(:,:)=D2ELEVTN(:,:)
      IF( LMEANSL ) THEN
         D2DWNELV(:,:)=D2ELEVTN(:,:)+D2MEANSL(:,:)
      ENDIF

   CONTAINS
      !==========================================================
      !+ READ_TOPO_BIN
      !+ READ_TOPO_CDF
      !+ SET_FLDSTG
      !+ SET_SLOPEMIX
      !==========================================================
   SUBROUTINE READ_TOPO_BIN
   USE CMF_UTILS_MOD,       only: mapR2vecD, CONV_END,  INQUIRE_FID
   IMPLICIT NONE
   !* local variables
   integer(KIND=JPIM)          :: ILFP
   real(KIND=JPRM),ALLOCATABLE :: R2TEMP(:,:)
   real(KIND=JPRB),ALLOCATABLE :: D2TEMP(:,:)
      !================================================
      allocate(R2TEMP(NX,NY))
      allocate(D2TEMP(NSEQMAX,1))

      TMPNAM=INQUIRE_FID()

      write(LOGNAM,*)'TOPO_INIT: unit-catchment area : ',TRIM(CGRAREA) 
      open(TMPNAM,FILE=CGRAREA,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
      read(TMPNAM,REC=1) R2TEMP(:,:)
      IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
      CALL mapR2vecD(R2TEMP,D2GRAREA)
      close(TMPNAM)

      write(LOGNAM,*)'TOPO_INIT: ground elevation : ',TRIM(CELEVTN)
      open(TMPNAM,FILE=CELEVTN,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
      read(TMPNAM,REC=1) R2TEMP(:,:)
      IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
      CALL mapR2vecD(R2TEMP,D2ELEVTN)
      close(TMPNAM)

      write(LOGNAM,*)'TOPO_INIT: downstream distance : ',TRIM(CNXTDST)
      open(TMPNAM,FILE=CNXTDST,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
      read(TMPNAM,REC=1) R2TEMP(:,:)
      IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
      CALL mapR2vecD(R2TEMP,D2NXTDST)
      close(TMPNAM)

      write(LOGNAM,*)'TOPO_INIT: river channel length : ',TRIM(CRIVLEN)
      open(TMPNAM,FILE=CRIVLEN,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
      read(TMPNAM,REC=1) R2TEMP(:,:)
      IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
      CALL mapR2vecD(R2TEMP,D2RIVLEN)
      close(TMPNAM)

      write(LOGNAM,*)'TOPO_INIT: floodplain elevation profile : ',TRIM(CFLDHGT)
      open(TMPNAM,FILE=TRIM(CFLDHGT),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
      DO ILFP=1,NLFP
         READ(TMPNAM,REC=ILFP) R2TEMP
         IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
         CALL mapR2vecD(R2TEMP,D2TEMP)
         D2FLDHGT(:,:,ILFP)= D2TEMP(:,:)
      ENDDO
      close(TMPNAM)

      !*** river channel / groundwater parameters)

      write(LOGNAM,*)'TOPO_INIT: river channel depth : ',TRIM(CRIVHGT)
      open(TMPNAM,FILE=CRIVHGT,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
      read(TMPNAM,REC=1) R2TEMP(:,:)
      IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
      CALL mapR2vecD(R2TEMP,D2RIVHGT)
      close(TMPNAM)

      write(LOGNAM,*)'TOPO_INIT: river channel width : ',TRIM(CRIVWTH)
      open(TMPNAM,FILE=CRIVWTH,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
      read(TMPNAM,REC=1) R2TEMP(:,:)
      IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
      CALL mapR2vecD(R2TEMP,D2RIVWTH)
      close(TMPNAM)

      write(LOGNAM,*)'TOPO_INIT: manning coefficient river: ',TRIM(CRIVMAN)
      open(TMPNAM,FILE=TRIM(CRIVMAN),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
      read(TMPNAM,REC=1) R2TEMP(:,:)
      IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
      CALL mapR2vecD(R2TEMP,D2RIVMAN)
      close(TMPNAM)

      IF( LGDWDLY )THEN
         write(LOGNAM,*)'TOPO_INIT: groundwater delay parameter: ',TRIM(CGDWDLY)
         open(TMPNAM,FILE=TRIM(CGDWDLY),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
         read(TMPNAM,REC=1) R2TEMP(:,:)
         IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
         CALL mapR2vecD(R2TEMP,D2GDWDLY)
         close(TMPNAM)
      ENDIF

      IF( LSLPMIX )THEN
         write(LOGNAM,*)'TOPO_INIT: LSLPMIX only used in IFS, not availabke with binary map'
      ENDIF
      IF( LSLOPEMOUTH )THEN
         write(LOGNAM,*)'TOPO_INIT: LSLOPEMOUTH only used in IFS, not availabke with binary map'
      ENDIF

      ! ==========

      IF( LMEANSL ) THEN
         write(LOGNAM, *)'TOPO_INIT: mean sea level: ', TRIM(CMEANSL)
         open(TMPNAM, FILE=CMEANSL, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=4*NX*NY)
         read(TMPNAM, REC=1) R2TEMP(:,:)
         IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
         CALL mapR2vecD(R2TEMP, D2MEANSL)
         close(TMPNAM)
      ENDIF

      deallocate(R2TEMP)
      deallocate(D2TEMP)

   END SUBROUTINE READ_TOPO_BIN
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE READ_TOPO_CDF
#ifdef UseCDF_CMF
   USE NETCDF 
   USE CMF_UTILS_MOD,            only: NCERROR,mapR2vecD
   USE YOS_CMF_MAP,              only: D2ELEVSLOPE     !! only used in ECMWF
   IMPLICIT NONE
   !* local variables
   integer(KIND=JPIM)               :: NCID,VARID,STATUS
   integer(KIND=JPIM)               :: ILEV
   real(KIND=JPRM),ALLOCATABLE      :: R2TEMP(:,:)
   real(KIND=JPRB),ALLOCATABLE      :: D2TEMP(:,:)
      !================================================
      allocate(R2TEMP(NX,NY))
      allocate(D2TEMP(NSEQMAX,1))

      !! CLIM FILE
      CALL NCERROR (NF90_OPEN(CRIVCLINC,NF90_NOWRITE,NCID),'opening '//TRIM(CRIVCLINC) )

      write(LOGNAM,*)'TOPO_INIT: ctmare:',TRIM(CRIVCLINC)
      CALL NCERROR ( NF90_INQ_VARID(NCID,'ctmare',VARID),'getting id' )
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' )
      CALL mapR2vecD(R2TEMP,D2GRAREA)

      write(LOGNAM,*)'TOPO_INIT: elevtn:',TRIM(CRIVCLINC)
      CALL NCERROR ( NF90_INQ_VARID(NCID,'elevtn',VARID),'getting id' )
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
      CALL mapR2vecD(R2TEMP,D2ELEVTN)

      write(LOGNAM,*)'TOPO_INIT: nxtdst:',TRIM(CRIVCLINC)
      CALL NCERROR ( NF90_INQ_VARID(NCID,'nxtdst',VARID),'getting id' )
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
      CALL mapR2vecD(R2TEMP,D2NXTDST)

      write(LOGNAM,*)'TOPO_INIT: rivlen:',TRIM(CRIVCLINC)
      CALL NCERROR ( NF90_INQ_VARID(NCID,'rivlen',VARID),'getting id' )
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
      CALL mapR2vecD(R2TEMP,D2RIVLEN)

      write(LOGNAM,*)'TOPO_INIT: fldhgt:',TRIM(CRIVCLINC)
      CALL NCERROR ( NF90_INQ_VARID(NCID,'fldhgt',VARID),'getting id' )
      DO ILEV=1,NLFP
         CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP,(/1,1,ILEV/),(/NX,NY,1/)),'reading data' ) 
         CALL mapR2vecD(R2TEMP,D2TEMP)
         D2FLDHGT(:,:,ILEV)=D2TEMP(:,:)
      ENDDO

      CALL NCERROR( NF90_CLOSE(NCID))

      IF ( LSLOPEMOUTH ) THEN
         allocate( D2ELEVSLOPE(NSEQMAX,1) )
         write(LOGNAM,*)'TOPO_INIT: elevslope:',TRIM(CRIVPARNC)
         STATUS = NF90_INQ_VARID(NCID,'elevslope',VARID)
         IF (STATUS /= 0 ) THEN
            write(LOGNAM,*)'TOPO_INIT: elevslope: not present, aborting'
            STOP 9 
         ELSE
            CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
         ENDIF 
         CALL mapR2vecD(R2TEMP,D2ELEVSLOPE)
      ENDIF

      !!========== 
      !! PAR FILE (river channel / groundwater parameters)
      CALL NCERROR (NF90_OPEN(CRIVPARNC,NF90_NOWRITE,NCID),'opening '//TRIM(CRIVPARNC) )

      write(LOGNAM,*)'TOPO_INIT: rivwth:',TRIM(CRIVPARNC)
      CALL NCERROR ( NF90_INQ_VARID(NCID,'rivwth',VARID),'getting id' )
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
      CALL mapR2vecD(R2TEMP,D2RIVWTH)

      write(LOGNAM,*)'TOPO_INIT: rivhgt:',TRIM(CRIVPARNC)
      CALL NCERROR ( NF90_INQ_VARID(NCID,'rivhgt',VARID),'getting id' )
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
      CALL mapR2vecD(R2TEMP,D2RIVHGT)

      write(LOGNAM,*)'TOPO_INIT: rivman:',TRIM(CRIVPARNC)
      CALL NCERROR ( NF90_INQ_VARID(NCID,'rivman',VARID),'getting id' )
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
      CALL mapR2vecD(R2TEMP,D2RIVMAN)

      IF ( LGDWDLY ) THEN
         write(LOGNAM,*)'TOPO_INIT: GDWDLY:',TRIM(CRIVPARNC)
         STATUS = NF90_INQ_VARID(NCID,'gdwdly',VARID)
         IF (STATUS /= 0 ) THEN
            write(LOGNAM,*)'TOPO_INIT: GDWDLY: not present, setting to zero'
            R2TEMP(:,:) = 0._JPRB
         ELSE
            CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
         ENDIF 
         CALL mapR2vecD(R2TEMP,D2GDWDLY)
      ENDIF

      I2MASK(:,:)=0_JPIM
      IF ( LSLPMIX ) THEN
         CALL SET_SLOPEMIX
      ENDIF

      CALL NCERROR( NF90_CLOSE(NCID))

      !!========== 
      !! MEAN SEA LEVEL FILE
      IF( LMEANSL ) THEN
         CALL NCERROR (NF90_OPEN(CMEANSLNC,NF90_NOWRITE,NCID),'opening '//TRIM(CMEANSLNC) )
         write(LOGNAM,*)'TOPO_INIT: rivhgt:',TRIM(CMEANSLNC)
         CALL NCERROR ( NF90_INQ_VARID(NCID,'meansl',VARID),'getting id' )
         CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
         CALL mapR2vecD ( R2TEMP,D2MEANSL  )
         CALL NCERROR ( NF90_CLOSE(NCID) )
      ENDIF 

      deallocate(R2TEMP)
      deallocate(D2TEMP)
#endif
   END SUBROUTINE READ_TOPO_CDF
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE SET_FLDSTG
   IMPLICIT NONE
   !* local variables
   integer(KIND=JPIM),SAVE  ::  ISEQ, I
   real(KIND=JPRB),SAVE     ::  DSTONOW
   real(KIND=JPRB),SAVE     ::  DSTOPRE
   real(KIND=JPRB),SAVE     ::  DHGTPRE
   real(KIND=JPRB),SAVE     ::  DWTHINC
!$OMP THREADPRIVATE               (I,DSTONOW,DSTOPRE,DHGTPRE,DWTHINC)
      !================================================
      P2FLDSTOMAX(:,:,:) = 0._JPRD
      D2FLDGRD(:,:,:)    = 0._JPRB
      DFRCINC=dble(NLFP)**(-1.)
      !
!$OMP PARALLEL DO
      DO ISEQ=1, NSEQALL
         DSTOPRE = P2RIVSTOMAX(ISEQ,1)
         DHGTPRE = 0._JPRB
         DWTHINC = D2GRAREA(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * DFRCINC
         DO I=1, NLFP
            DSTONOW = D2RIVLEN(ISEQ,1) * ( D2RIVWTH(ISEQ,1) + DWTHINC*(DBLE(I)-0.5) ) * (D2FLDHGT(ISEQ,1,I)-DHGTPRE)
            P2FLDSTOMAX(ISEQ,1,I) = DSTOPRE + DSTONOW
            D2FLDGRD(ISEQ,1,I) = (D2FLDHGT(ISEQ,1,I)-DHGTPRE) * DWTHINC**(-1.)
            DSTOPRE = P2FLDSTOMAX(ISEQ,1,I)
            DHGTPRE = D2FLDHGT(ISEQ,1,I)
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

   !
   END SUBROUTINE SET_FLDSTG
   !==========================================================
   !+
   !+
   !+
   !==========================================================
   SUBROUTINE SET_SLOPEMIX    !! only used in IFS0
#ifdef UseCDF_CMF
   USE NETCDF 
   USE CMF_UTILS_MOD,           only: NCERROR,mapI2vecI

   IMPLICIT NONE
   integer(KIND=JPIM),ALLOCATABLE  :: I2TEMP(:,:)
   integer(KIND=JPIM)              :: ISEQ, I0, I1
   integer(KIND=JPIM)              :: NCID,VARID,STATUS

      allocate(I2TEMP(NX,NY))
      write(LOGNAM,*)'TOPO_INIT: mask_slope:',TRIM(CRIVPARNC)
      STATUS =  NF90_INQ_VARID(NCID,'mask_slope',VARID)
      IF (STATUS /= 0 ) THEN
         write(LOGNAM,*)'TOPO_INIT: mask_slope: LSLPMIX should be set to FALSE: ABORTING!'
         STOP 9
      ENDIF 
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2TEMP),'reading data' ) 
      CALL mapI2vecI(I2TEMP,I2MASK)
      I0=0
      I1=0
      DO ISEQ=1,NSEQALL
         IF (I2MASK(ISEQ,1) == 1 ) THEN  !! kinematic wave applied
            I1=I1+1
         ENDIF
         IF (I2MASK(ISEQ,1) == 0 ) THEN 
            I0=I0+1
         ENDIF
      ENDDO
      write(LOGNAM,*)'TOPO_INIT: sum(mask==0), sum(mask==1)',I0,I1
      IF ( I0+I1 .NE. NSEQALL ) THEN 
         write(LOGNAM,*)'TOPO_INIT: mask==0 + mask == 1 does not match NSEQALL.. something wrong, aborting'
         STOP 9
      ENDIF 

      deallocate(I2TEMP)
#endif
   END SUBROUTINE SET_SLOPEMIX

   END SUBROUTINE CMF_TOPO_INIT
!####################################################################


END MODULE CMF_CTRL_MAPS_MOD
