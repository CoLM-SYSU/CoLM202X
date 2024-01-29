MODULE CMF_CTRL_DAMOUT_MOD
!==========================================================
!* PURPOSE: CaMa-Flood reservoir operation scheme (under development)
!
! (C) R. Hanazaki & D.Yamazaki (U-Tokyo)  Feb 2020
!
!* CONTAINS:
! -- CMF_DEM_NMLIST  : Read setting from namelist
! -- CMF_DAM_INIT    : Initialize dam data
! -- CMF_CALC_DAMOUT : Calculate inflow and outflow at dam
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
! 
!* Updated: May 2023, by Shulei Zhang(zhangshlei@mail.sysu.edu.cn)
!  -- Introducing three additional reservoir operation schemes (H06, V13, and LIS) along with their respective input reading modules.
!  -- To access the reservoir parameter preprocessing code, please refer to the "preprocessing" folder.
!==========================================================
   USE PARKIND1,                only: JPIM, JPRB, JPRM, JPRD
   USE CMF_UTILS_MOD,           only: INQUIRE_FID
   USE YOS_CMF_INPUT,           only: LOGNAM, IMIS, LDAMOUT, LPTHOUT, LRESTART
   USE YOS_CMF_INPUT,           only: NX, NY, DT
   USE YOS_CMF_MAP,             only: I2VECTOR, I1NEXT, NSEQALL, NSEQRIV, NSEQMAX
   USE YOS_CMF_MAP,             only: NPTHOUT,  NPTHLEV, PTH_UPST, PTH_DOWN, PTH_ELV, I2MASK!! bifurcation pass
   USE YOS_CMF_PROG,            only: D2RIVOUT, D2FLDOUT, P2RIVSTO, P2FLDSTO, P2DAMSTO, P2DAMINF, D2RUNOFF
   USE YOS_CMF_DIAG,            only: D2RIVINF, D2FLDINF
   USE YOS_CMF_TIME,            only: ISYYYY        

   !============================
   IMPLICIT NONE
   SAVE
   !*** NAMELIST/NDAMOUT/
   character(LEN=256)              :: CDAMFILE         !! dam parameter files
   character(LEN=3)                :: LDAMOPT          !! dam scheme
   logical                         :: LDAMTXT          !! true: dam inflow-outflw txt output
   NAMELIST/NDAMOUT/   CDAMFILE, LDAMOPT, LDAMTXT
   !*** CDAMFILE contains:
   ! damloc.csv: basic information, provided by GRAND
   ! damflow.csv: flow characteristics, estimated by simulated natural flow
   ! damsto.csv: storage characteristics, estimated using GRSAD & ReGeom datasets
   ! damfcperiod.csv: flood control period, estimated by simulated natural flow (for V13)
   ! /water_use_data: water use grid and daily water use for each dam (for H06 and V13)
   !*** LDAMOPT contains:
   ! H06: use Hanasaki 2006 scheme
   ! V13: use Voisin 2013 scheme
   ! LIS: use LISFLOOD scheme
   ! H22: use Hanazaki 2022 scheme

   !*** dam map
   integer(KIND=JPIM),ALLOCATABLE  :: DamSeq(:)            !! coresponding ISEQ of each dam
   integer(KIND=JPIM),ALLOCATABLE  :: I1DAM(:)             !! dam map: 1=dam, 10=upstream of dam, 11: dam grid & downstream is also dam, 0=other

   !*** dam basic information
   integer(KIND=JPIM)              :: IDAM, NDAM           !! number of dams
   integer(KIND=JPIM)              :: NDAMX                !! exclude dams
   integer(KIND=JPIM),ALLOCATABLE  :: GRanD_ID(:)          !! GRanD ID
   character(LEN=256)              :: DamName              !! dam name
   integer(KIND=JPIM)              :: IX, IY               !! IX,IY of dam grid
   real(KIND=JPRB)                 :: DamLon, DamLat       !! longitude, latitude of dam body
   real(KIND=JPRB)                 :: totalsto             !! total storage capacity of reservoir (mcm)
   real(KIND=JPRB),ALLOCATABLE     :: upreal(:)            !! observed drainage area of reservoir (km2)
   character(LEN=256),ALLOCATABLE  :: MainUse(:)           !! main use of dam
   integer(KIND=JPIM)              :: CYear                !! construction year

   !*** dam parameters
   real(KIND=JPRB),ALLOCATABLE     :: Qn(:), Qf(:)         !! Qn: normal discharge; Qf: flood discharge (m3/s)
   real(KIND=JPRB),ALLOCATABLE     :: TotVol(:)            !! total storage capacity of reservoir (mcm)
   real(KIND=JPRB),ALLOCATABLE     :: FldVol(:)            !! flood control storage (mcm)
   real(KIND=JPRB),ALLOCATABLE     :: NorVol(:)            !! normal storage (mcm)
   real(KIND=JPRB),ALLOCATABLE     :: ConVol(:)            !! conservative storage (mcm)

   !*** water use data 
   real(KIND=JPRB),ALLOCATABLE     :: WUSE_DD(:,:)         !! daily water demand (m3/s)
   real(KIND=JPRB),ALLOCATABLE     :: WUSE_AD(:)           !! average daily demand (m3/s)

   !*** dam parameters for different schemes
   !*** H22
   real(KIND=JPRB),ALLOCATABLE     :: H22_EmeVol(:)        !! storage volume to start emergency operation (m3)
   real(KIND=JPRB),ALLOCATABLE     :: H22_k(:)             !! release coefficient (-)
   !*** H06
   real(KIND=JPRB),ALLOCATABLE     :: H06_c(:)             !! the ratio between capacity and mean annual inflow (-)
   real(KIND=JPRB),ALLOCATABLE     :: H06_DPI(:)           !! the ratio between annual mean demand and annual mean inflow (-) 
   !*** V13
   integer(KIND=JPIM),ALLOCATABLE  :: StFC_Mth(:)          !! the start of the flood control period
   integer(KIND=JPIM),ALLOCATABLE  :: NdFC_Mth(:)          !! the end of the flood control period
   integer(KIND=JPIM),ALLOCATABLE  :: StOP_Mth(:)          !! the start of the operational year

CONTAINS
   !####################################################################
   !* CONTAINS:
   ! -- CMF_DEMOUT_NMLIST  : Read setting from namelist
   ! -- CMF_DAMOUT_INIT    : Initialize dam data
   ! -- CMF_DAMOUT_CALC    : Calculate inflow and outflow at dam
   !####################################################################
   SUBROUTINE CMF_DAMOUT_NMLIST
   ! reed setting from namelist
   ! -- Called from CMF_DRV_NMLIST
   USE YOS_CMF_INPUT,      only: CSETFILE,NSETFILE,LDAMOUT
   USE CMF_UTILS_MOD,      only: INQUIRE_FID
   IMPLICIT NONE
   !================================================
      write(LOGNAM,*) ""
      write(LOGNAM,*) "!---------------------!"

      !*** 1. open namelist
      NSETFILE=INQUIRE_FID()
      open(NSETFILE,FILE=CSETFILE,STATUS="OLD")
      write(LOGNAM,*) "CMF::DAMOUT_NMLIST: namelist open in unit: ", TRIM(CSETFILE), NSETFILE 

      !*** 2. read namelist
      rewind(NSETFILE)
      read(NSETFILE,NML=NDAMOUT)

      IF( LDAMOUT )THEN
         write(LOGNAM,*)   "=== NAMELIST, NDAMOUT ==="
         write(LOGNAM,*)   "CDAMFILE: ", CDAMFILE
         write(LOGNAM,*)   "LDAMOPT: " , LDAMOPT
         write(LOGNAM,*)   "LDAMTXT: " , LDAMTXT
      ENDIF

      close(NSETFILE)

      write(LOGNAM,*) "CMF::DAMOUT_NMLIST: end" 

   END SUBROUTINE CMF_DAMOUT_NMLIST
   !####################################################################


   !####################################################################
   SUBROUTINE CMF_DAMOUT_INIT
   USE CMF_UTILS_MOD,      only: INQUIRE_FID
   USE YOS_CMF_INPUT,      only: NX, NY, LRESTART, LPTHOUT
   USE YOS_CMF_MAP,        only: I2VECTOR, I1NEXT, NSEQALL, NSEQMAX
   USE YOS_CMF_PROG,       only: P2RIVSTO, P2DAMSTO, P2DAMINF
   USE YOS_CMF_MAP,        only: NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN, PTH_ELV , I2MASK !! bifurcation pass
   USE YOS_CMF_TIME,       only: ISYYYY        

   ! read setting from CDAMFILE
   IMPLICIT NONE
   integer(KIND=JPIM)         :: NDAMFILE
   integer(KIND=JPIM)         :: ISEQ, JSEQ
   integer(KIND=JPIM)         :: IPTH, ILEV, ISEQP, JSEQP
   character(LEN=256)         :: CDAMFILE_tmp

      !####################################################################
      !! ================ read dam basic information ================
      write(LOGNAM,*) "!---------------------!"
      CDAMFILE_tmp = trim(CDAMFILE)//'damloc.csv'
      write(LOGNAM,*) "CMF::DAMOUT_INIT: initialize dam", trim(CDAMFILE_tmp) 

      NDAMFILE=INQUIRE_FID()
      open(NDAMFILE,FILE=CDAMFILE_tmp,STATUS="OLD")
      read(NDAMFILE,*) NDAM
      read(NDAMFILE,*)        

      write(LOGNAM,*) "CMF::DAMOUT_INIT: number of dams", NDAM

      !! --- allocate ---
      allocate(DamSeq(NDAM))
      allocate(I1DAM(NSEQMAX))
      allocate(GRanD_ID(NDAM))
      allocate(MainUse(NDAM))
      allocate(upreal(NDAM))

      !! --------------
      DamSeq(:)= IMIS
      I1DAM(:) = 0
      NDAMX    = 0

      DO IDAM = 1, NDAM
         read(NDAMFILE,*) GRanD_ID(IDAM), DamName, DamLon, DamLat, totalsto, upreal(IDAM), MainUse(IDAM), CYear, IX, IY

         !! --------------
         IF (CYear > ISYYYY) CYCLE        !! check construction year
         IF (IX<=0 .or. IX > NX .or. IY<=0 .or. IY > NY ) CYCLE
         ISEQ=I2VECTOR(IX,IY)  
         IF( I1NEXT(ISEQ)==-9999 .or. ISEQ<=0 ) CYCLE
         NDAMX=NDAMX+1

         !! --------------
         DamSeq(IDAM)=ISEQ
         I1DAM(ISEQ)=1
         I2MASK(ISEQ,1)=2   !! reservoir grid. skipped for adaptive time step

      ENDDO
      close(NDAMFILE)

      write(LOGNAM,*) "CMF::DAMOUT_INIT: allocated dams:", NDAMX 

      !! ================ read dam flow parameters ================
      write(LOGNAM,*) "!---------------------!"
      CDAMFILE_tmp = trim(CDAMFILE)//'damflow.csv'
      write(LOGNAM,*) "CMF::DAMOUT_INIT: dam flow parameters:", trim(CDAMFILE_tmp) 

      NDAMFILE=INQUIRE_FID()
      open(NDAMFILE,FILE=CDAMFILE_tmp,STATUS="OLD")
      read(NDAMFILE,*) 

      !! --- allocate ---
      allocate(Qn(NDAM), Qf(NDAM))

      DO IDAM = 1, NDAM
         read(NDAMFILE,*) GRanD_ID(IDAM), Qn(IDAM), Qf(IDAM)
      ENDDO
      close(NDAMFILE)

      !! ================ read dam storage parameters ================
      write(LOGNAM,*) "!---------------------!"
      CDAMFILE_tmp = trim(CDAMFILE)//'damsto.csv'
      write(LOGNAM,*) "CMF::DAMOUT_INIT: dam storage parameters:", trim(CDAMFILE_tmp) 

      NDAMFILE=INQUIRE_FID()
      open(NDAMFILE,FILE=CDAMFILE_tmp,STATUS="OLD")
      read(NDAMFILE,*) 

      !! --- allocate ---
      allocate(TotVol(NDAM), FldVol(NDAM), NorVol(NDAM), ConVol(NDAM))

      DO IDAM = 1, NDAM
         read(NDAMFILE,*) GRanD_ID(IDAM), TotVol(IDAM), FldVol(IDAM), NorVol(IDAM), ConVol(IDAM)
         TotVol(IDAM) = TotVol(IDAM) * 1.E6    !! from Million Cubic Meter to m3
         FldVol(IDAM) = FldVol(IDAM) * 1.E6
         NorVol(IDAM) = NorVol(IDAM) * 1.E6
         ConVol(IDAM) = ConVol(IDAM) * 1.E6
      ENDDO
      close(NDAMFILE)

      !! ================ Read parameters for different schemes ================
      IF(LDAMOPT == "H06")THEN  
         !! --- allocate ---
         allocate(H06_DPI(NDAM), H06_c(NDAM))

         H06_DPI = WUSE_AD / Qn    
         H06_c = TotVol / (Qn * 86400. * 365.) !! from m3/s to m3

         !! --- Read water use data ---
         write(LOGNAM,*) "CMF:: READ_WATER_USE "
         CALL READ_WATER_USE
      ENDIF

      IF(LDAMOPT == "H22")THEN 
         !! --- allocate ---
         allocate(H22_EmeVol(NDAM), H22_k(NDAM))

         H22_EmeVol = FldVol + (TotVol - FldVol)*0.2  
         H22_k =  MAX((1. - (TotVol - FldVol) / upreal /0.2),0._JPRB)
      ENDIF

      IF(LDAMOPT == "V13")THEN  
         !! --- allocate ---
         allocate(H06_DPI(NDAM), H06_c(NDAM))

         H06_DPI = WUSE_AD / Qn    
         H06_c = TotVol / (Qn * 86400. * 365.) !! from m3/s to m3

         !! --- read fcperiod.csv data ---
         write(LOGNAM,*) "!---------------------!"
         CDAMFILE_tmp = trim(CDAMFILE)//'damfcperiod.csv'
         write(LOGNAM,*) "CMF::DAMOUT_INIT: dam flow parameters:", trim(CDAMFILE_tmp) 

         NDAMFILE=INQUIRE_FID()
         open(NDAMFILE,FILE=CDAMFILE_tmp,STATUS="OLD")
         read(NDAMFILE,*) 

         !! --- allocate ---
         allocate(StFC_Mth(NDAM), NdFC_Mth(NDAM),StOP_Mth(NDAM))

         DO IDAM = 1, NDAM
            read(NDAMFILE,*) GRanD_ID(IDAM), StFC_Mth(IDAM), NdFC_Mth(IDAM),StOP_Mth(IDAM)
         ENDDO
         close(NDAMFILE)

         !! --- Read water use data ---
         write(LOGNAM,*) "CMF:: READ_WATER_USE "
         CALL READ_WATER_USE
      ENDIF
      !####################################################################

      !! mark upstream of dam grid, for applying kinematic wave routine to suppress storage buffer effect.
      DO ISEQ=1, NSEQALL
         IF( I1DAM(ISEQ)==0 .and. I1NEXT(ISEQ)>0 )THEN !! if target is non-dam grid
            JSEQ=I1NEXT(ISEQ)
            IF( I1DAM(JSEQ)==1 .or. I1DAM(JSEQ)==11 )THEN !! if downstream is dam
               I1DAM(ISEQ)=10            !! mark upstream of dam grid by "10"
               I2MASK(ISEQ,1)=1   !! reservoir upstream grid. skipped for adaptive time step
            ENDIF
         ENDIF

         IF( I1DAM(ISEQ)==1 .and. I1NEXT(ISEQ)>0 )THEN !! if target is dam grid
            JSEQ=I1NEXT(ISEQ)
               IF( I1DAM(JSEQ)==1 .or. I1DAM(JSEQ)==11 )THEN !! if downstream is dam
                  I1DAM(ISEQ)=11            !! mark upstream of dam grid by "11"
                  I2MASK(ISEQ,1)=2   !! reservoir grid (cascading). skipped for adaptive time step
               ENDIF
         ENDIF
      ENDDO

      !! Initialize dam storage
      IF( .not. LRESTART )THEN
         P2DAMSTO(:,1)=0._JPRD
         DO IDAM=1, NDAM
               IF( DamSeq(IDAM)>0 )THEN
                  ISEQ=DamSeq(IDAM)
                  P2DAMSTO(ISEQ,1)=NorVol(IDAM)*0.5  !! set initial storage to Normal Storage Volume
                  P2RIVSTO(ISEQ,1)=MAX(P2RIVSTO(ISEQ,1),NorVol(IDAM)*0.5_JPRD) !! also set initial river storage, in order to keep consistency
               ENDIF
         ENDDO
      ENDIF

      !! Initialize dam inflow
      DO ISEQ=1, NSEQALL
         P2DAMINF(ISEQ,1)=0._JPRD
      ENDDO

      !! Stop bifurcation at dam & dam-upstream grids
      IF( LPTHOUT )THEN
         DO IPTH=1, NPTHOUT
            ISEQP=PTH_UPST(IPTH)
            JSEQP=PTH_DOWN(IPTH)
            IF( ISEQP<=0 .or. JSEQP<=0) CYCLE
            IF( I1DAM(ISEQP)>0 .or. I1DAM(JSEQP)>0 )THEN
               DO ILEV=1, NPTHLEV
                  PTH_ELV(IPTH,ILEV)=1.E20  !! no bifurcation
               ENDDO
            ENDIF
         ENDDO
      ENDIF

   END SUBROUTINE CMF_DAMOUT_INIT



   SUBROUTINE READ_WATER_USE
   USE YOS_CMF_TIME,       only: NSTEPS   

   IMPLICIT NONE
   character(LEN=256)         :: CDAMFILE_WUSE_YEAR
   character(LEN=4)           :: CYYYY
   character(LEN=16)          :: tmp_i, tmp_name  
   integer(KIND=JPIM)         :: NDAMFILE

      !=======================================  
      write(CYYYY,'(I4)') ISYYYY
      CDAMFILE_WUSE_YEAR = trim(CDAMFILE)//'water_use_data/'//trim(adjustl(CYYYY))//'.txt'
      write(LOGNAM,*) "CMF::DAMOUT_INIT: dam water use file:", trim(CDAMFILE_WUSE_YEAR)

      !! --- allocate ---
      !! from dam water use data
      allocate(WUSE_DD(NDAM, NSTEPS))    
      allocate(WUSE_AD(NDAM))          

      !! read dam water use
      NDAMFILE=INQUIRE_FID()
      open(NDAMFILE,FILE=trim(CDAMFILE_WUSE_YEAR),STATUS="OLD")

      DO IDAM = 1, NDAM         
         read(NDAMFILE,*) tmp_i, tmp_name, WUSE_DD(IDAM,:)  
         WUSE_AD(IDAM) = sum(WUSE_DD(IDAM,:))/NSTEPS
      ENDDO
      close(NDAMFILE)

   END SUBROUTINE READ_WATER_USE

   !####################################################################
   SUBROUTINE CMF_DAMOUT_WATBAL
   IMPLICIT NONE

   ! SAVE for OMP
   integer(KIND=JPIM),SAVE    :: ISEQD
   !*** water balance
   real(KIND=JPRB),SAVE       :: DamInflow
   real(KIND=JPRB),SAVE       :: DamOutflw           !! Total outflw 
   real(KIND=JPRD),SAVE       :: GlbDAMSTO, GlbDAMSTONXT, GlbDAMINF, GlbDAMOUT, DamMiss
   
!$OMP THREADPRIVATE    (ISEQD,DamInflow,DamOutflw)
      ! ==========================================
      !* 4) update reservoir storage and check water DamMiss --------------------------
      GlbDAMSTO    = 0._JPRB
      GlbDAMSTONXT = 0._JPRB
      GlbDAMINF    = 0._JPRB
      GlbDAMOUT    = 0._JPRB
      
!$OMP PARALLEL DO REDUCTION(+:GlbDAMSTO, GlbDAMSTONXT, GlbDAMINF, GlbDAMOUT)
      DO IDAM=1, NDAM
         IF( DamSeq(IDAM)<=0 ) CYCLE
         ISEQD = DamSeq(IDAM)
         
         DamInflow = D2RIVINF(ISEQD,1) + D2FLDINF(ISEQD,1) + D2RUNOFF(ISEQD,1)
         DamOutflw = D2RIVOUT(ISEQD,1) + D2FLDOUT(ISEQD,1)
         !!P2DAMINF(ISEQD,1)=DamInflow   !! if water balance needs to be checked in the output file, P2DAMINF should be updated.
         
         GlbDAMSTO = GlbDAMSTO + P2DAMSTO(ISEQD,1)
         GlbDAMINF = GlbDAMINF + DamInflow*DT
         GlbDAMOUT = GlbDAMOUT + DamOutflw*DT
         
         P2DAMSTO(ISEQD,1) = P2DAMSTO(ISEQD,1) + DamInflow * DT - DamOutflw * DT
         
         GlbDAMSTONXT = GlbDAMSTONXT + P2DAMSTO(ISEQD,1)
      ENDDO
!$OMP END PARALLEL DO
  
      DamMiss = GlbDAMSTO-GlbDAMSTONXT+GlbDAMINF-GlbDAMOUT
      write(LOGNAM,*) "CMF::DAM_CALC: DamMiss at all dams:", DamMiss*1.D-9
  
   END SUBROUTINE CMF_DAMOUT_WATBAL
!####################################################################


!####################################################################
   SUBROUTINE CMF_DAMOUT_CALC
   USE YOS_CMF_INPUT,      only: DT
   USE YOS_CMF_MAP,        only: I1NEXT,   NSEQALL,  NSEQRIV
   USE YOS_CMF_PROG,       only: D2RIVOUT, D2FLDOUT, P2RIVSTO, P2FLDSTO
   USE YOS_CMF_PROG,       only: P2DAMSTO, P2DAMINF, D2RUNOFF   
   USE YOS_CMF_DIAG,       only: D2RIVINF, D2FLDINF
   USE YOS_CMF_TIME,       only: KSTEP, ISMM       

   ! local
   IMPLICIT NONE
   ! SAVE for OMP
   integer(KIND=JPIM),SAVE    :: ISEQD
   !** dam variables
   real(KIND=JPRB),SAVE       :: DamVol
   real(KIND=JPRB),SAVE       :: DamInflow
   real(KIND=JPRB),SAVE       :: DamOutflw           !! Total outflw 
   ! real(KIND=JPRB),SAVE       :: DamOutTmp           !! Total outflw 

   !*** water balance
   real(KIND=JPRD),SAVE       :: GlbDAMSTO, GlbDAMSTONXT, GlbDAMINF, GlbDAMOUT, DamMiss

!$OMP THREADPRIVATE    (ISEQD,DamVol,DamInflow,DamOutflw)
!====================
!CONTAINS
!+ UPDATE_INFLOW: replace dam upstream with kinamatic wave, calculate inflow to dam
!+ MODIFY_OUTFLW: modify outflw to avoid negative storage
!+ DAM_OPERATION_H06
!+ DAM_OPERATION_V13
!+ DAM_OPERATION_LIS
!+ DAM_OPERATION_H22
!==========================================================

!* (1) Replace discharge in upstream grids with kinematic outflow
!     to avoid storage buffer effect (Shin et al., 2019, WRR)
! ------  rivout at upstream grids of dam, rivinf to dam grids are updated.
      CALL UPDATE_INFLOW        

!* (2) Reservoir Operation
!====================================
!$OMP PARALLEL DO
      DO IDAM=1, NDAM
         IF( DamSeq(IDAM)<=0 ) CYCLE
         ISEQD=DamSeq(IDAM)

         !! *** 2a update dam volume and inflow -----------------------------------
         DamVol    = P2DAMSTO(ISEQD,1)    
         DamInflow = P2DAMINF(ISEQD,1)

         !================================ 2b Reservoir Operation ================================ 
         !! option: Hanasaki 2006 scheme
         IF( LDAMOPT == "H06" )THEN
            CALL DAM_OPERATION_H06(DamVol, DamInflow, DamOutflw, &
                                 MainUse(IDAM), TotVol(IDAM), Qn(IDAM),H06_DPI(IDAM), H06_c(IDAM), &
                                 WUSE_DD(IDAM, KSTEP), WUSE_AD(IDAM))
         ENDIF

         !! option: Voisin 2013 scheme
         IF( LDAMOPT == "V13" )THEN
            CALL DAM_OPERATION_V13(DamVol, DamInflow, DamOutflw, &
                                 MainUse(IDAM), TotVol(IDAM), Qn(IDAM),H06_DPI(IDAM), H06_c(IDAM), &
                                 WUSE_DD(IDAM, KSTEP), WUSE_AD(IDAM),&
                                 StFC_Mth(IDAM), NdFC_Mth(IDAM), StOP_Mth(IDAM))
         ENDIF

         !! option: LISFLOOD scheme
         IF( LDAMOPT == "LIS" )THEN
            CALL DAM_OPERATION_LISFLOOD(DamVol, DamInflow, DamOutflw, &
                                       ConVol(IDAM), NorVol(IDAM), FldVol(IDAM), TotVol(IDAM), Qn(IDAM), Qf(IDAM))
         ENDIF

         !! option: Hanazaki 2022 scheme
         IF( LDAMOPT == "H22" )THEN 
            CALL DAM_OPERATION_H22(DamVol, DamInflow, DamOutflw, &
                                    FldVol(IDAM)*0.5, FldVol(IDAM), H22_EmeVol(IDAM), Qn(IDAM), Qf(IDAM),H22_k(IDAM))
         ENDIF

         !! *** 2c flow limitter
         DamOutflw = MIN( DamOutflw, DamVol/DT, real(P2RIVSTO(ISEQD,1)+P2FLDSTO(ISEQD,1),JPRB)/DT )
         DamOutflw = MAX( DamOutflw, 0._JPRB ) 

         !! update CaMa variables  (treat all outflow as RIVOUT in dam grid, no fldout)
         D2RIVOUT(ISEQD,1) = DamOutflw
         D2FLDOUT(ISEQD,1) = 0._JPRB
      ENDDO
!$OMP END PARALLEL DO
!====================================


!* 3) modify outflow to suppless negative discharge, update RIVOUT,FLDOUT,RIVINF,FLDINF
      CALL MODIFY_OUTFLW


!* 4) update reservoir storage and check water DamMiss --------------------------
      GlbDAMSTO    = 0._JPRB
      GlbDAMSTONXT = 0._JPRB
      GlbDAMINF    = 0._JPRB
      GlbDAMOUT    = 0._JPRB

!$OMP PARALLEL DO REDUCTION(+:GlbDAMSTO, GlbDAMSTONXT, GlbDAMINF, GlbDAMOUT)
      DO IDAM=1, NDAM
         IF( DamSeq(IDAM)<=0 ) CYCLE
         ISEQD = DamSeq(IDAM)

         DamInflow = D2RIVINF(ISEQD,1) + D2FLDINF(ISEQD,1) + D2RUNOFF(ISEQD,1)
         DamOutflw = D2RIVOUT(ISEQD,1) + D2FLDOUT(ISEQD,1)
         !!P2DAMINF(ISEQD,1)=DamInflow   !! if water balance needs to be checked in the output file, P2DAMINF should be updated.

         GlbDAMSTO = GlbDAMSTO + P2DAMSTO(ISEQD,1)
         GlbDAMINF = GlbDAMINF + DamInflow*DT
         GlbDAMOUT = GlbDAMOUT + DamOutflw*DT

         P2DAMSTO(ISEQD,1) = P2DAMSTO(ISEQD,1) + DamInflow * DT - DamOutflw * DT

         GlbDAMSTONXT = GlbDAMSTONXT + P2DAMSTO(ISEQD,1)
      ENDDO
!$OMP END PARALLEL DO

      DamMiss = GlbDAMSTO-GlbDAMSTONXT+GlbDAMINF-GlbDAMOUT
      write(LOGNAM,*) "CMF::DAM_CALC: DamMiss at all dams:", DamMiss*1.D-9


CONTAINS
!==========================================================
!+ UPDATE_INFLOW: replace dam upstream with kinamatic wave, calculate inflow to dam
!+ MODIFY_OUTFLW: modify outflw to avoid negative storage
!+ DAM_OPERATION_H06
!+ DAM_OPERATION_V13
!+ DAM_OPERATION_LIS
!+ DAM_OPERATION_H22
!==========================================================

   SUBROUTINE DAM_OPERATION_H06(Vol_dam, inflw, outflw, &
                             dam_purpose, Vol_tot, Q_n, DPI, c, &
                             dam_wuse, dam_wuse_avg)
   !! reference: Hanasaki, N., Kanae, S., & Oki, T. (2006).
   !! A reservoir operation scheme for global river routing models. 
   !! Journal of Hydrology, 327(1-2), 22-41.

   IMPLICIT NONE
   real(KIND=JPRB), intent(in)             :: Vol_dam       ! Current volume of the dam (m3)
   real(KIND=JPRB), intent(in)             :: inflw         ! Inflow rate to the dam (m3/s)
   real(KIND=JPRB), intent(out)            :: outflw        ! Outflow rate from the dam (m3/s)
   !*** parameter
   character(LEN=256), intent(in)          :: dam_purpose   ! main use of dam
   real(KIND=JPRB), intent(in)             :: Vol_tot       ! total storage capacity (m3)
   real(KIND=JPRB), intent(in)             :: Q_n           ! normal discharge (m3/s)
   real(KIND=JPRB), intent(in)             :: DPI           ! the ratio between annual mean demand and annual mean inflo (-)
   real(KIND=JPRB), intent(in)             :: c             ! the ratio between capacity and mean annual inflow (-)
   real(KIND=JPRB), intent(in)             :: dam_wuse      ! daily water demand (m3/s)
   real(KIND=JPRB), intent(in)             :: dam_wuse_avg  ! average daily demand (m3/s)
   !*** local
   real(KIND=JPRB)                         :: a=0.85        ! adjustment factor for Krls
   real(KIND=JPRB)                         :: M=0.5         ! the ratio between minimum release and long‐term annual mean inflow
   real(KIND=JPRB)                         :: R             ! demand‐controlled release ratio = MIN(1, a*c)
   real(KIND=JPRB)                         :: Krls          ! the ratio between initial storage and the long‐term target storage = S/(a*C)
   real(KIND=JPRB)                         :: OutTmp        ! temporary outflow
   !=====================================================
      ! !! parameter
      ! M = 0.5
      ! a = 0.85 
      !! parameter calculation
      Krls =  Vol_dam / (Vol_tot * a)
      R    =  MIN(1., a * c)

      !! calculate DamOutTmp
      IF( dam_purpose == "Irrigation" )THEN
         !! irrigation reservoir
         IF ( DPI < (1. - M) ) THEN
            OutTmp = Q_n + dam_wuse - dam_wuse_avg
         ELSE
            OutTmp = Q_n * (M + (1. - M) * dam_wuse / dam_wuse_avg)
         ENDIF
      ELSE
         !! not irrigation reservoir
         OutTmp = Q_n 
      ENDIF

      !! calculate DamOutflw
      IF ( c >= 0.5 ) THEN
         outflw = Krls * OutTmp 
      ELSE
         outflw = R * Krls * OutTmp + (1. - R) * inflw
      ENDIF

   END SUBROUTINE DAM_OPERATION_H06


   SUBROUTINE DAM_OPERATION_V13(Vol_dam, inflw, outflw, &
                             dam_purpose, Vol_tot, Q_n, DPI, c, &
                             dam_wuse, dam_wuse_avg, &
                             MthStFC,MthNdFC,MthStOP)
   !! reference: Voisin, N., Li, H., Ward, D., Huang, M., Wigmosta, M., & Leung, L. R. (2013). 
   !! On an improved sub-regional water resources management representation for integration into earth system models. 
   !! Hydrology and Earth System Sciences, 17(9), 3605-3622.
   !! reference: https://github.com/IMMM-SFA/mosartwmpy

   IMPLICIT NONE
   real(KIND=JPRB), intent(in)             :: Vol_dam       ! Current volume of the dam (m3)
   real(KIND=JPRB), intent(in)             :: inflw         ! Inflow rate to the dam (m3/s)
   real(KIND=JPRB), intent(out)            :: outflw        ! Outflow rate from the dam (m3/s)
   !*** parameter
   character(LEN=256), intent(in)          :: dam_purpose   ! main use of dam
   real(KIND=JPRB), intent(in)             :: Vol_tot       ! total storage capacity (m3)
   real(KIND=JPRB), intent(in)             :: Q_n           ! normal discharge (m3/s)
   real(KIND=JPRB), intent(in)             :: DPI           ! the ratio between annual mean demand and annual mean inflo (-)
   real(KIND=JPRB), intent(in)             :: c             ! the ratio between capacity and mean annual inflow (-)
   real(KIND=JPRB), intent(in)             :: dam_wuse      ! daily water demand (m3/s)
   real(KIND=JPRB), intent(in)             :: dam_wuse_avg  ! average daily demand (m3/s)
   integer(KIND=JPIM), intent(in)          :: MthStFC       ! the start of the flood control period
   integer(KIND=JPIM), intent(in)          :: MthNdFC       ! the end of the flood control period
   integer(KIND=JPIM), intent(in)          :: MthStOP       ! the start of the operational year
   !*** local
   real(KIND=JPRB)                         :: a=0.85        ! adjustment factor for Krls
   real(KIND=JPRB)                         :: M=0.5         ! the ratio between minimum release and long‐term annual mean inflow
   real(KIND=JPRB)                         :: R             ! demand‐controlled release ratio = MIN(1, a*c)
   real(KIND=JPRB)                         :: Krls          ! the ratio between initial storage and the long‐term target storage = S/(a*C)
   real(KIND=JPRB)                         :: OutTmp        ! temporary outflow
   real(KIND=JPRB)                         :: drop          ! temporary drop flow
   !=====================================================
      !! parameter calculation
      Krls =  Vol_dam / (Vol_tot * a)
      R    =  MIN(1., a * c)

      IF( dam_purpose == "Irrigation" .or. dam_purpose == "Flood-control") THEN
         !! irrigation reservoir
         !! calculate OutTmp
         IF ( DPI < (1. - M) ) THEN
            OutTmp = Q_n + dam_wuse - dam_wuse_avg
         ELSE
            OutTmp = Q_n * (M + (1. - M) * dam_wuse / dam_wuse_avg)
         ENDIF

         !! conmbined irrigation & flood control
         !*** stage-1: from MthStFC to MthNdFC
         IF ( MthStFC <= MthNdFC ) THEN
            IF ( ISMM >= MthStFC .and. ISMM < MthNdFC) THEN
               IF ( inflw < Q_n ) THEN
                  drop = abs(inflw - Q_n)
                  OutTmp = OutTmp + drop
               ENDIF
            ENDIF
         ELSEIF ( MthStFC > MthNdFC ) THEN
            IF ( ISMM >= MthStFC .or. ISMM < MthNdFC) THEN
               IF ( inflw < Q_n ) THEN
                  drop = abs(inflw - Q_n)
                  OutTmp = OutTmp + drop
               ENDIF
            ENDIF
         ENDIF

         !*** stage-2: from MthNdFC to MthStOP
         IF ( MthNdFC <= MthStOP ) THEN
            IF ( ISMM >= MthNdFC .and. ISMM < MthStOP ) THEN
               ! IF ( inflw > Q_n ) THEN
               !   fill = abs(inflw - Q_n)
               ! end
               IF (OutTmp > Q_n) THEN
                  OutTmp = Q_n
               ENDIF
            ENDIF
         ELSEIF ( MthNdFC > MthStOP ) THEN
            IF ( ISMM >= MthNdFC .or. ISMM < MthStOP ) THEN
               ! IF ( inflw > Q_n ) THEN
               !   fill = abs(inflw - Q_n)
               ! end
               IF (OutTmp > Q_n) THEN
                  OutTmp = Q_n
               ENDIF
            ENDIF
         ENDIF

      !! not irrigation reservoir
      ELSE    
         OutTmp = Q_n   
      ENDIF

      !! calculate DamOutflw
      IF ( c >= 0.5 ) THEN
         outflw = Krls * OutTmp 
      ELSE
         outflw = R * Krls * OutTmp + (1. - R) * inflw
      ENDIF

  ! the end of the flood control period (NDFC) is defined as the first month of the wet period preceding the start of the operational year
  ! The start of the flood control period (STFC) is defined as the month with the lowest flow within the dry period preceding the start of the operational period.
  ! the first month at which the long-term mean monthly flow falls below the long- term mean annual flow

   END SUBROUTINE DAM_OPERATION_V13


   SUBROUTINE DAM_OPERATION_LISFLOOD(Vol_dam, inflw, outflw, &
                                    Vol_con, Vol_nor, Vol_fld, Vol_tot, Q_n, Q_f)
   !! reference: https://ec-jrc.github.io/lisflood-model/3_03_optLISFLOOD_reservoirs/
   !! reference: https://github.com/ec-jrc/lisflood-code

   IMPLICIT NONE
   ! Input variables
   real(KIND=JPRB), intent(in)             :: Vol_dam     ! Current volume of the dam (m3)
   real(KIND=JPRB), intent(in)             :: inflw       ! Inflow rate to the dam (m3/s)
   real(KIND=JPRB), intent(out)            :: outflw      ! Outflow rate from the dam (m3/s)
   ! Parameters
   real(KIND=JPRB), intent(in)             :: Vol_con     ! volume of water at the conservation level of the dam (m3)
   real(KIND=JPRB), intent(in)             :: Vol_nor     ! volume of water at the normal level of the dam (m3)
   real(KIND=JPRB), intent(in)             :: Vol_fld     ! volume of water at the flood controld level of the dam (m3)
   real(KIND=JPRB), intent(in)             :: Vol_tot     ! total storage capacity (m3)
   real(KIND=JPRB), intent(in)             :: Q_n         ! normal discharge (m3/s)
   real(KIND=JPRB), intent(in)             :: Q_f         ! flood control discharge  (m3/s)

   ! Local variables
   real(KIND=JPRB)                         :: F            ! fractional level of dam
   real(KIND=JPRB)                         :: L_c          ! fractional level at conservation level
   real(KIND=JPRB)                         :: L_n          ! fractional level at normal operating level
   real(KIND=JPRB)                         :: L_f          ! fractional level at flood control level
   real(KIND=JPRB)                         :: Q_min        ! minimum outflow rate
   real(KIND=JPRB)                         :: Q_max        ! maximum outflow rate
   real(KIND=JPRB)                         :: AdjLn        ! adjustment factor for L_n,in range 0.01 to 0.99
   real(KIND=JPRB)                         :: AdjQnor      ! adjustment factor for Q_nor, in range 0.25 to 2
   real(KIND=JPRB)                         :: Ladj_f       ! adjusted fractional level for maximum outflow
   real(KIND=JPRB)                         :: Qadj_nor     ! adjusted normal outflow rate

      !=====================================================
      AdjLn = 0.5 ! adjustment factor for L_n,in range 0.01 to 0.99
      AdjQnor = 1 ! adjustment factor for Q_n, in range 0.25 to 2
      
      !*** Parameters calculation  
      F = Vol_dam / Vol_tot
      L_c = Vol_con / Vol_tot
      L_n = Vol_nor / Vol_tot
      L_f = Vol_fld / Vol_tot
      Ladj_f  = L_n + AdjLn * (L_f - L_n)

      Q_min = Q_n * 0.1
      Qadj_nor = MAX(Q_min, MIN(AdjQnor * Q_n, Q_f)) 
      
      !*** operation scheme
      IF (F <= 2* L_c) THEN
         outflw = MIN(Q_min, Vol_dam/86400.)
      
      ELSEIF (F > 2* L_c .and. F <= L_n) THEN
         outflw = Q_min + (Qadj_nor - Q_min) * ((F - 2*L_c)/(L_n - 2*L_c))

      ELSEIF (F > L_n .and. F <= Ladj_f) THEN
         outflw = Qadj_nor 

      ELSEIF (F > Ladj_f .and. F <= L_f) THEN
         outflw = Qadj_nor + (Q_f - Qadj_nor) * ((F - Ladj_f)/(L_f - Ladj_f))
      
      ELSEIF (F > L_f) THEN
         Q_max = MIN(Q_f, MAX(1.2 * inflw, Qadj_nor))   
         outflw = MAX(Q_max, (F - L_f - 0.01) * (Vol_tot/86400.))
      ENDIF
      
      ! the condition described below is applied in order to prevent outflow values that are too large compared to the inflow value.
      ! reference: https://github.com/ec-jrc/lisflood-code
      IF (F < L_f .and. outflw > MIN(1.2 * inflw, Qadj_nor)) THEN
         outflw = MIN(outflw, MAX(inflw, Qadj_nor))
      ENDIF
      
   END SUBROUTINE DAM_OPERATION_LISFLOOD


   SUBROUTINE DAM_OPERATION_H22(Vol_dam, inflw, outflw, &
                              Vol_con, Vol_fld, Vol_eme, Q_n, Q_f, k)
   !! reference: Hanazaki, R., Yamazaki, D., & Yoshimura, K. (2022). 
   !! Development of a reservoir flood control scheme for global flood models. 
   !! Journal of Advances in Modeling Earth Systems, 14(3), e2021MS002944.
   
   IMPLICIT NONE
   real(KIND=JPRB), intent(in)                :: Vol_dam     ! Current volume of the dam (m3)
   real(KIND=JPRB), intent(in)                :: inflw       ! Inflow rate to the dam (m3/s)
   real(KIND=JPRB), intent(out)               :: outflw      ! Outflow rate from the dam (m3/s)
   !*** parameter
   real(KIND=JPRB), intent(in)                :: Vol_con     ! conservative storage (m3) 
   real(KIND=JPRB), intent(in)                :: Vol_fld     ! flood control storage (m3) 
   real(KIND=JPRB), intent(in)                :: Vol_eme     ! emergency storage (m3)   
   real(KIND=JPRB), intent(in)                :: Q_n         ! normal discharge (m3/s)
   real(KIND=JPRB), intent(in)                :: Q_f         ! flood discharge (m3/s)
   real(KIND=JPRB), intent(in)                :: k           ! release coefficient
   !=====================================================

      !! case1: impoundment
      IF( Vol_dam <=  Vol_con)THEN
         outflw = Q_n * (Vol_dam / Vol_fld )

      !! case2: water supply
      ELSEIF( Vol_con < Vol_dam .and. Vol_dam <= Vol_fld )THEN
         IF( Q_f <= inflw )THEN
            outflw = Q_n * 0.5 +   (Vol_dam-Vol_con)/( Vol_fld-Vol_con)      * (Q_f - Q_n)
         ELSE
            outflw = Q_n * 0.5 + (((Vol_dam-Vol_con)/( Vol_eme-Vol_con))**2) * (Q_f - Q_n)
         ENDIF  

      !! case3: flood control
      ELSEIF( Vol_fld < Vol_dam .and. Vol_dam <= Vol_eme ) THEN
         IF( Q_f <= inflw ) THEN
            outflw = Q_f +  k * (Vol_dam-Vol_fld)/(Vol_eme-Vol_fld) * (inflw-Q_f)
         ELSE
            outflw = Q_n*0.5 + (((Vol_dam-Vol_con)/(Vol_eme-Vol_con))**2)* (Q_f - Q_n)
         ENDIF

      !! case4: emergency operation
      ELSE
         outflw = MAX(inflw, Q_f)
      ENDIF

   END SUBROUTINE DAM_OPERATION_H22


   SUBROUTINE UPDATE_INFLOW
   USE YOS_CMF_INPUT,      only: PMINSLP, PMANFLD
   USE YOS_CMF_MAP,        only: D2RIVLEN, D2RIVMAN, D2ELEVTN, D2NXTDST, D2RIVWTH
   USE YOS_CMF_PROG,       only: D2RIVOUT_PRE, D2FLDOUT_PRE
   USE YOS_CMF_DIAG,       only: D2RIVDPH, D2RIVVEL, D2FLDDPH
   IMPLICIT NONE
   ! SAVE for OpenMP
   integer(KIND=JPIM),SAVE    :: ISEQ, JSEQ
   real(KIND=JPRB),SAVE       :: DSLOPE,DAREA,DVEL,DSLOPE_F,DARE_F,DVEL_F
!$OMP THREADPRIVATE     (JSEQ,DSLOPE,DAREA,DVEL,DSLOPE_F,DARE_F,DVEL_F)
   !============================

   !*** 1a. reset outflw & dam inflow
!$OMP PARALLEL DO
      DO ISEQ=1, NSEQALL
         IF( I1DAM(ISEQ)>0 )THEN  !! if dam grid or upstream of dam, reset variables
            D2RIVOUT(ISEQ,1) = 0._JPRB
            D2FLDOUT(ISEQ,1) = 0._JPRB
            P2DAMINF(ISEQ,1) = 0._JPRD
         ENDIF
      ENDDO
!$OMP END PARALLEL DO

   !*** 1b. calculate dam inflow, using previous tstep discharge
#ifndef NoAtom_CMF
!$OMP PARALLEL DO  !! No OMP Atomic for bit-identical simulation (set in Mkinclude)
#endif
      DO ISEQ=1, NSEQALL
         IF( I1DAM(ISEQ)==10 .or. I1DAM(ISEQ)==11 )THEN  !! if dam grid or upstream of dam
            JSEQ=I1NEXT(ISEQ)
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
            P2DAMINF(JSEQ,1) = P2DAMINF(JSEQ,1) + D2RIVOUT_PRE(ISEQ,1) + D2FLDOUT_PRE(ISEQ,1) 
         ENDIF
      ENDDO
#ifndef NoAtom_CMF
!$OMP END PARALLEL DO
#endif

   !*** 1c. discharge for upstream grids of dams
!$OMP PARALLEL DO  !! No OMP Atomic for bit-identical simulation (set in Mkinclude)
      DO ISEQ=1, NSEQRIV
         IF( I1DAM(ISEQ)==10 )THEN  !! if downstream is DAM
            JSEQ   = I1NEXT(ISEQ)
            ! === river flow
            DSLOPE = (D2ELEVTN(ISEQ,1)-D2ELEVTN(JSEQ,1)) * D2NXTDST(ISEQ,1)**(-1.)
            DSLOPE = MAX(DSLOPE,PMINSLP)

            DVEL   = D2RIVMAN(ISEQ,1)**(-1.) * DSLOPE**0.5 * D2RIVDPH(ISEQ,1)**(2./3.)
            DAREA  = D2RIVWTH(ISEQ,1) * D2RIVDPH(ISEQ,1)

            D2RIVVEL(ISEQ,1) = DVEL
            D2RIVOUT(ISEQ,1) = DAREA * DVEL
            D2RIVOUT(ISEQ,1) = MIN( D2RIVOUT(ISEQ,1), real(P2RIVSTO(ISEQ,1),JPRB)/DT )
            !=== floodplain flow
            DSLOPE_F = MIN( 0.005_JPRB,DSLOPE )    !! set MIN [instead of using weirequation for efficiency]
            DVEL_F   = PMANFLD**(-1.) * DSLOPE_F**0.5 * D2FLDDPH(ISEQ,1)**(2./3.)
            DARE_F   = P2FLDSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.)
            DARE_F   = MAX( DARE_F - D2FLDDPH(ISEQ,1)*D2RIVWTH(ISEQ,1), 0._JPRB )   !!remove above river channel     area

            D2FLDOUT(ISEQ,1) = DARE_F * DVEL_F
            D2FLDOUT(ISEQ,1) = MIN(  D2FLDOUT(ISEQ,1)*1._JPRD, P2FLDSTO(ISEQ,1)/DT )
         ENDIF
      ENDDO
!$OMP END PARALLEL DO

   END SUBROUTINE UPDATE_INFLOW
!==========================================================
!+
!+
!+
!==========================================================
   SUBROUTINE MODIFY_OUTFLW
   ! modify outflow in order to avoid negative storage
   USE YOS_CMF_MAP,        only: NSEQMAX
   IMPLICIT NONE

   real(KIND=JPRD)            :: P2STOOUT(NSEQMAX,1)                      !! total outflow from a grid     [m3]
   real(KIND=JPRD)            :: P2RIVINF(NSEQMAX,1)                      !! 
   real(KIND=JPRD)            :: P2FLDINF(NSEQMAX,1)                      !! 

   real(KIND=JPRB)            :: D2RATE(NSEQMAX,1)                        !! outflow correction
   ! SAVE for OpenMP
   integer(KIND=JPIM),SAVE    :: ISEQ, JSEQ
   real(KIND=JPRB),SAVE       :: OUT_R1, OUT_R2, OUT_F1, OUT_F2, DIUP, DIDW
   !$OMP THREADPRIVATE     (JSEQ,OUT_R1, OUT_R2, OUT_F1, OUT_F2, DIUP, DIDW)
   !================================================
  
!*** 1. initialize & calculate P2STOOUT for normal cells

!$OMP PARALLEL DO
      DO ISEQ=1, NSEQALL
         P2RIVINF(ISEQ,1) = 0._JPRD
         P2FLDINF(ISEQ,1) = 0._JPRD
         P2STOOUT(ISEQ,1) = 0._JPRD
         D2RATE(ISEQ,1) = 1._JPRB
      ENDDO
!$OMP END PARALLEL DO

!! for normal cells ---------
#ifndef NoAtom_CMF
!$OMP PARALLEL DO
#endif
      DO ISEQ=1, NSEQRIV                                                    !! for normalcells
         JSEQ=I1NEXT(ISEQ) ! next cell's pixel
         OUT_R1 = MAX(  D2RIVOUT(ISEQ,1),0._JPRB )
         OUT_R2 = MAX( -D2RIVOUT(ISEQ,1),0._JPRB )
         OUT_F1 = MAX(  D2FLDOUT(ISEQ,1),0._JPRB )
         OUT_F2 = MAX( -D2FLDOUT(ISEQ,1),0._JPRB )
         DIUP=(OUT_R1+OUT_F1)*DT
         DIDW=(OUT_R2+OUT_F2)*DT
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
         P2STOOUT(ISEQ,1) = P2STOOUT(ISEQ,1) + DIUP 
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
         P2STOOUT(JSEQ,1) = P2STOOUT(JSEQ,1) + DIDW 
      ENDDO
#ifndef NoAtom_CMF
!$OMP END PARALLEL DO
#endif

!! for river mouth grids ------------
!$OMP PARALLEL DO
      DO ISEQ=NSEQRIV+1, NSEQALL
         OUT_R1 = MAX( D2RIVOUT(ISEQ,1), 0._JPRB )
         OUT_F1 = MAX( D2FLDOUT(ISEQ,1), 0._JPRB )
         P2STOOUT(ISEQ,1) = P2STOOUT(ISEQ,1) + OUT_R1*DT + OUT_F1*DT
      ENDDO
!$OMP END PARALLEL DO

!============================
!*** 2. modify outflow

!$OMP PARALLEL DO
      DO ISEQ=1, NSEQALL
         IF ( P2STOOUT(ISEQ,1) > 1.E-8 ) THEN
            D2RATE(ISEQ,1) = MIN( (P2RIVSTO(ISEQ,1)+P2FLDSTO(ISEQ,1)) * P2STOOUT(ISEQ,1)**(-1.), 1._JPRD )
         ENDIF
      ENDDO
!$OMP END PARALLEL DO

!! normal pixels------
#ifndef NoAtom_CMF
!$OMP PARALLEL DO  !! No OMP Atomic for bit-identical simulation (set in Mkinclude)
#endif
      DO ISEQ=1, NSEQRIV ! for normal pixels
         JSEQ=I1NEXT(ISEQ)
         IF( D2RIVOUT(ISEQ,1) >= 0._JPRB )THEN
            D2RIVOUT(ISEQ,1) = D2RIVOUT(ISEQ,1)*D2RATE(ISEQ,1)
            D2FLDOUT(ISEQ,1) = D2FLDOUT(ISEQ,1)*D2RATE(ISEQ,1)
         ELSE
            D2RIVOUT(ISEQ,1) = D2RIVOUT(ISEQ,1)*D2RATE(JSEQ,1)
            D2FLDOUT(ISEQ,1) = D2FLDOUT(ISEQ,1)*D2RATE(JSEQ,1)
         ENDIF
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
         P2RIVINF(JSEQ,1) = P2RIVINF(JSEQ,1) + D2RIVOUT(ISEQ,1)             !! total inflow to a grid (from upstream)
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
         P2FLDINF(JSEQ,1) = P2FLDINF(JSEQ,1) + D2FLDOUT(ISEQ,1)
      ENDDO
#ifndef NoAtom_CMF
!$OMP END PARALLEL DO
#endif

      D2RIVINF(:,:)=P2RIVINF(:,:)  !! needed for SinglePrecisionMode
      D2FLDINF(:,:)=P2FLDINF(:,:)

!! river mouth-----------------
!$OMP PARALLEL DO
      DO ISEQ=NSEQRIV+1, NSEQALL
         D2RIVOUT(ISEQ,1) = D2RIVOUT(ISEQ,1)*D2RATE(ISEQ,1)
         D2FLDOUT(ISEQ,1) = D2FLDOUT(ISEQ,1)*D2RATE(ISEQ,1)
      ENDDO
!$OMP END PARALLEL DO

   END SUBROUTINE MODIFY_OUTFLW
!==========================================================

   END SUBROUTINE CMF_DAMOUT_CALC
!####################################################################



!####################################################################
   SUBROUTINE CMF_DAMOUT_WRTE
   USE YOS_CMF_PROG,       only: P2DAMSTO, P2DAMINF, D2RIVOUT
   USE YOS_CMF_TIME,       only: IYYYYMMDD,ISYYYY
   USE CMF_UTILS_MOD,      only: INQUIRE_FID

   ! local
   character(LEN=36)          :: WriteTXT(NDAMX) !, WriteTXT2(NDAMX)

   ! File IO
   integer(KIND=JPIM),SAVE    :: ISEQD, JDAM
   integer(KIND=JPIM),SAVE    :: LOGDAM
   character(LEN=4),SAVE      :: CYYYY
   character(LEN=256),SAVE    :: CLEN, CFMT
   character(LEN=256),SAVE    :: DAMTXT
   logical,SAVE               :: IsOpen
   DATA IsOpen       /.FALSE./

      ! ======

      IF( LDAMTXT .and. LDAMTXT)THEN

         IF( .not. IsOpen)THEN
            IsOpen=.TRUE.
            write(CYYYY,'(i4.4)') ISYYYY
            DAMTXT='./damtxt-'//trim(CYYYY)//'.txt'

            LOGDAM=INQUIRE_FID()
            open(LOGDAM,FILE=DAMTXT,FORM='formatted')

            write(CLEN,'(i0)') NDAMX
            CFMT="(i10,"//TRIM(CLEN)//"(a36))"

            JDAM=0
            DO IDAM=1, NDAM
               IF( DamSeq(IDAM)<=0 ) CYCLE
               JDAM=JDAM+1
               
               write(WriteTxt(JDAM), '(i12)') GRanD_ID(IDAM)
            ENDDO

            write(LOGDAM,CFMT) NDAMX, (WriteTXT(JDAM) ,JDAM=1, NDAMX)

         ENDIF

         JDAM=0
         DO IDAM=1, NDAM
            IF( DamSeq(IDAM)<=0 ) CYCLE
            JDAM=JDAM+1
            ISEQD=DamSeq(IDAM)
            
            write(WriteTxt(JDAM), '(3f12.3)') P2DAMSTO(ISEQD,1)*1.E-6, P2DAMINF(ISEQD,1), D2RIVOUT(ISEQD,1) !！ P2DAMSTO m3 to Million Cubic Meter
         ENDDO

         CFMT="(i10,"//TRIM(CLEN)//"(a36))"  
         write(LOGDAM,CFMT) IYYYYMMDD, (WriteTXT(JDAM),JDAM=1, NDAMX)

      ENDIF

   END SUBROUTINE CMF_DAMOUT_WRTE
!####################################################################


END MODULE CMF_CTRL_DAMOUT_MOD

