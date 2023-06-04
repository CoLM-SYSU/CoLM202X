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
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM, JPRD
USE YOS_CMF_INPUT,           ONLY: LOGNAM, IMIS, LDAMOUT
!============================
IMPLICIT NONE
SAVE
!*** NAMELIST/NDAM/
CHARACTER(LEN=256)              :: CDAMFILE         !! dam parameter files
LOGICAL                         :: LDAMTXT          !! true: dam inflow-outflw txt output
CHARACTER(LEN=3)                :: LDAMOPT          !! dam scheme
NAMELIST/NDAMOUT/   CDAMFILE, LDAMTXT, LDAMOPT 
!*** CDAMFILE contains:
! damloc.csv: basic information, provided by GRAND
! damflow.csv: flow characteristics, estimated by simulated natural flow
! damsto.csv: storage characteristics, estimated using GRSAD & ReGeom datasets
! damfcperiod.csv: flood control period, estimated by simulated natural flow (for V13)
! /water_use: water use grid and daily water use for each dam (for H06 and V13)
!*** LDAMOPT contains:
! H06: use Hanasaki 2006 scheme
! V13: use Voisin 2013 scheme
! LIS: use LISFLOOD scheme
! H22: use Hanazaki 2022 scheme

!*** dam map
INTEGER(KIND=JPIM),ALLOCATABLE  :: DamSeq(:)   !! coresponding ISEQ of each dam
INTEGER(KIND=JPIM),ALLOCATABLE  :: I1DAM(:)    !! dam map: 1=dam, 10=upstream of dam, 11: dam grid & downstream is also dam, 0=other

!*** dam basic information
INTEGER(KIND=JPIM)              :: IDAM, NDAM           !! number of dams
INTEGER(KIND=JPIM)              :: NDAMX                !! exclude dams
INTEGER(KIND=JPIM),ALLOCATABLE  :: GRanD_ID(:)          !! GRanD ID
CHARACTER(LEN=256)              :: DamName              !! dam name
INTEGER(KIND=JPIM)              :: IX, IY               !! IX,IY of dam grid
REAL(KIND=JPRB)                 :: DamLon, DamLat       !! longitude, latitude of dam body
REAL(KIND=JPRB)                 :: totalsto             !! total storage capacity of reservoir (mcm)
REAL(KIND=JPRB),ALLOCATABLE     :: upreal(:)            !! observed drainage area of reservoir (km2)
CHARACTER(LEN=256),ALLOCATABLE  :: MainUse(:)           !! main use of dam
INTEGER(KIND=JPIM)              :: CYear                !! construction year

!*** dam parameters
REAL(KIND=JPRB),ALLOCATABLE     :: Qn(:), Qf(:)         !! Qn: normal discharge; Qf: flood discharge (m3/s)
REAL(KIND=JPRB),ALLOCATABLE     :: TotVol(:)            !! total storage capacity of reservoir (mcm)
REAL(KIND=JPRB),ALLOCATABLE     :: FldVol(:)            !! flood control storage (mcm)
REAL(KIND=JPRB),ALLOCATABLE     :: NorVol(:)            !! normal storage (mcm)
REAL(KIND=JPRB),ALLOCATABLE     :: ConVol(:)            !! conservative storage (mcm)

!*** water use data 
REAL(KIND=JPRB),ALLOCATABLE     :: WUSE_DD(:,:)         !! daily water demand (m3/s)
REAL(KIND=JPRB),ALLOCATABLE     :: WUSE_AD(:)           !! average daily demand (m3/s)

!*** dam parameters for different schemes
!*** H22
REAL(KIND=JPRB),ALLOCATABLE     :: H22_EmeVol(:)        !! storage volume to start emergency operation (m3)
REAL(KIND=JPRB),ALLOCATABLE     :: H22_k(:)             !! release coefficient (-)
!*** H06
REAL(KIND=JPRB),ALLOCATABLE     :: H06_c(:)             !! the ratio between capacity and mean annual inflow (-)
REAL(KIND=JPRB),ALLOCATABLE     :: H06_DPI(:)           !! the ratio between annual mean demand and annual mean inflow (-) 
!*** V13
INTEGER(KIND=JPIM),ALLOCATABLE  :: StFC_Mth(:)          !! the start of the flood control period
INTEGER(KIND=JPIM),ALLOCATABLE  :: NdFC_Mth(:)          !! the end of the flood control period
INTEGER(KIND=JPIM),ALLOCATABLE  :: StOP_Mth(:)          !! the start of the operational year

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
USE YOS_CMF_INPUT,      ONLY: CSETFILE,NSETFILE,LDAMOUT
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"

!*** 1. open namelist
NSETFILE=INQUIRE_FID()
OPEN(NSETFILE,FILE=CSETFILE,STATUS="OLD")
WRITE(LOGNAM,*) "CMF::DAMOUT_NMLIST: namelist OPEN in unit: ", TRIM(CSETFILE), NSETFILE 

!*** 2. default value
! CDAMFILE="./dam_params.csv"
! LDAMTXT=.TRUE.
! SDAM_H22=.TRUE.

!*** 3. read namelist
REWIND(NSETFILE)
READ(NSETFILE,NML=NDAMOUT)

IF( LDAMOUT )THEN
  WRITE(LOGNAM,*)   "=== NAMELIST, NDAMOUT ==="
  WRITE(LOGNAM,*)   "CDAMFILE: ", CDAMFILE
  WRITE(LOGNAM,*)   "LDAMTXT: " , LDAMTXT
  WRITE(LOGNAM,*)   "LDAMOPT: " , LDAMOPT
ENDIF

CLOSE(NSETFILE)

WRITE(LOGNAM,*) "CMF::DAMOUT_NMLIST: end" 

END SUBROUTINE CMF_DAMOUT_NMLIST
!####################################################################


!####################################################################
SUBROUTINE CMF_DAMOUT_INIT
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
USE YOS_CMF_INPUT,      ONLY: NX, NY, LRESTART, LPTHOUT
USE YOS_CMF_MAP,        ONLY: I2VECTOR, I1NEXT, NSEQALL, NSEQMAX
USE YOS_CMF_PROG,       ONLY: P2RIVSTO, P2DAMSTO, P2DAMINF
USE YOS_CMF_MAP,        ONLY: NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN, PTH_ELV , I2MASK !! bifurcation pass
USE YOS_CMF_TIME,       ONLY: ISYYYY        

! read setting from CDAMFILE
IMPLICIT NONE
INTEGER(KIND=JPIM)         :: NDAMFILE
INTEGER(KIND=JPIM)         :: ISEQ, JSEQ
INTEGER(KIND=JPIM)         :: IPTH, ILEV, ISEQP, JSEQP
CHARACTER(LEN=256)         :: CDAMFILE_tmp

!####################################################################
!! ================ READ dam basic information ================
WRITE(LOGNAM,*) "!---------------------!"
CDAMFILE_tmp = trim(CDAMFILE)//'damloc.csv'
WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: initialize dam", trim(CDAMFILE_tmp) 

NDAMFILE=INQUIRE_FID()
OPEN(NDAMFILE,FILE=CDAMFILE_tmp,STATUS="OLD")
READ(NDAMFILE,*) NDAM
READ(NDAMFILE,*)        

WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: number of dams", NDAM

!! --- ALLOCATE ---
!! calculate from CDAMFILE
ALLOCATE(DamSeq(NDAM))
!! dam map, dam variable
ALLOCATE(I1DAM(NSEQMAX))
!! from CDAMFILE
ALLOCATE(GRanD_ID(NDAM))
ALLOCATE(MainUse(NDAM))
ALLOCATE(upreal(NDAM))

!! --------------
DamSeq(:)= IMIS
I1DAM(:) = 0
NDAMX    = 0

DO IDAM = 1, NDAM
  READ(NDAMFILE,*) GRanD_ID(IDAM), DamName, DamLon, DamLat, totalsto, upreal(IDAM), MainUse(IDAM), CYear, IX, IY

  !! --------------
  IF (CYear > ISYYYY) cycle        !! check construction year
  IF (IX<=0 .or. IX > NX .or. IY<=0 .or. IY > NY ) cycle
  ISEQ=I2VECTOR(IX,IY)  
  IF( I1NEXT(ISEQ)==-9999 .or. ISEQ<=0 ) cycle
  NDAMX=NDAMX+1

  !! --------------
  DamSeq(IDAM)=ISEQ
  I1DAM(ISEQ)=1
  I2MASK(ISEQ,1)=2   !! reservoir grid. skipped for adaptive time step

END DO
CLOSE(NDAMFILE)

WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: allocated dams:", NDAMX 

!! ================ READ dam flow parameters ================
WRITE(LOGNAM,*) "!---------------------!"
CDAMFILE_tmp = trim(CDAMFILE)//'damflow.csv'
WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: dam flow parameters:", trim(CDAMFILE_tmp) 

NDAMFILE=INQUIRE_FID()
OPEN(NDAMFILE,FILE=CDAMFILE_tmp,STATUS="OLD")
READ(NDAMFILE,*) 

!! --- ALLOCATE ---
ALLOCATE(Qn(NDAM), Qf(NDAM))

DO IDAM = 1, NDAM
  READ(NDAMFILE,*) GRanD_ID(IDAM), Qn(IDAM), Qf(IDAM)
END DO
CLOSE(NDAMFILE)

!! ================ READ dam storage parameters ================
WRITE(LOGNAM,*) "!---------------------!"
CDAMFILE_tmp = trim(CDAMFILE)//'damsto.csv'
WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: dam storage parameters:", trim(CDAMFILE_tmp) 

NDAMFILE=INQUIRE_FID()
OPEN(NDAMFILE,FILE=CDAMFILE_tmp,STATUS="OLD")
READ(NDAMFILE,*) 

!! --- ALLOCATE ---
ALLOCATE(TotVol(NDAM), FldVol(NDAM), NorVol(NDAM), ConVol(NDAM))

DO IDAM = 1, NDAM
  READ(NDAMFILE,*) GRanD_ID(IDAM), TotVol(IDAM), FldVol(IDAM), NorVol(IDAM), ConVol(IDAM)
  TotVol(IDAM) = TotVol(IDAM) * 1.E6    !! from Million Cubic Meter to m3
  FldVol(IDAM) = FldVol(IDAM) * 1.E6
  NorVol(IDAM) = NorVol(IDAM) * 1.E6
  ConVol(IDAM) = ConVol(IDAM) * 1.E6
END DO
CLOSE(NDAMFILE)

!! ================ READ water use data ================
IF(LDAMOPT == "H06" .OR. LDAMOPT == "V13")THEN  
  WRITE(LOGNAM,*) "CMF:: READ_WATER_USE "
  CALL READ_WATER_USE 
ENDIF

!! ================ parameters for different schemes ================
IF(LDAMOPT == "H06" .OR. LDAMOPT == "V13")THEN  
  !! --- ALLOCATE ---
  ALLOCATE(H06_DPI(NDAM), H06_c(NDAM))

  H06_DPI = WUSE_AD / Qn    
  H06_c = TotVol / (Qn * 86400. * 365.) !! Qn m3/s to m3
ENDIF

IF(LDAMOPT == "H22")THEN 
  !! --- ALLOCATE ---
  ALLOCATE(H22_EmeVol(NDAM), H22_k(NDAM))

  H22_EmeVol = FldVol + (TotVol - FldVol)*0.2  
  H22_k =  max((1. - (TotVol - FldVol) / upreal /0.2),0._JPRB)
ENDIF

IF(LDAMOPT == "V13")THEN  
 WRITE(LOGNAM,*) "!---------------------!"
  CDAMFILE_tmp = trim(CDAMFILE)//'damfcperiod.csv'
  WRITE(LOGNAM,*) "CMF::DAMOUT_INIT: dam flow parameters:", trim(CDAMFILE_tmp) 

  NDAMFILE=INQUIRE_FID()
  OPEN(NDAMFILE,FILE=CDAMFILE_tmp,STATUS="OLD")
  READ(NDAMFILE,*) 

  !! --- ALLOCATE ---
  ALLOCATE(StFC_Mth(NDAM), NdFC_Mth(NDAM),StOP_Mth(NDAM))

  DO IDAM = 1, NDAM
    READ(NDAMFILE,*) GRanD_ID(IDAM), StFC_Mth(IDAM), NdFC_Mth(IDAM),StOP_Mth(IDAM)
  END DO
  CLOSE(NDAMFILE)
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
END DO

!! Initialize dam storage
IF( .not. LRESTART )THEN
  P2DAMSTO(:,1)=0._JPRD
  DO IDAM=1, NDAM
    IF( DamSeq(IDAM)>0 )THEN
      ISEQ=DamSeq(IDAM)
      P2DAMSTO(ISEQ,1)=NorVol(IDAM)*0.5  !! set initial storage to Normal Storage Volume
      P2RIVSTO(ISEQ,1)=max(P2RIVSTO(ISEQ,1),NorVol(IDAM)*0.5_JPRD) !! also set initial river storage, in order to keep consistency
    ENDIF
  END DO
ENDIF

!! Initialize dam inflow
DO ISEQ=1, NSEQALL
  P2DAMINF(ISEQ,1)=0._JPRD
END DO

!! Stop bifurcation at dam & dam-upstream grids
IF( LPTHOUT )THEN
  DO IPTH=1, NPTHOUT
    ISEQP=PTH_UPST(IPTH)
    JSEQP=PTH_DOWN(IPTH)
    IF( ISEQP<=0 .or. JSEQP<=0) CYCLE
    IF( I1DAM(ISEQP)>0 .or. I1DAM(JSEQP)>0 )THEN
      DO ILEV=1, NPTHLEV
        PTH_ELV(IPTH,ILEV)=1.E20  !! no bifurcation
      END DO
    ENDIF
  END DO
ENDIF
!####################################################################

CONTAINS
SUBROUTINE READ_WATER_USE
  USE YOS_CMF_TIME,       ONLY: NSTEPS   

  IMPLICIT NONE
  CHARACTER(LEN=256)         :: CDAMFILE_WUSE_YEAR
  CHARACTER(len=4)           :: CYYYY
  CHARACTER(LEN=16)          :: tmp_i, tmp_name  
  !=======================================  
  write(CYYYY,'(I4)') ISYYYY
  CDAMFILE_WUSE_YEAR = trim(CDAMFILE)//'water_use/'//trim(adjustl(CYYYY))//'.txt'
  write(LOGNAM,*) "CMF::DAMOUT_INIT: dam water use file:", trim(CDAMFILE_WUSE_YEAR)

  !! --- ALLOCATE ---
  !! from dam water use data
  ALLOCATE(WUSE_DD(NDAM, NSTEPS))    
  ALLOCATE(WUSE_AD(NDAM))          

  !! read dam water use
  NDAMFILE=INQUIRE_FID()
  OPEN(NDAMFILE,FILE=trim(CDAMFILE_WUSE_YEAR),STATUS="OLD")

  DO IDAM = 1, NDAM
    READ(NDAMFILE,*) tmp_i, tmp_name, WUSE_DD(IDAM,:)  
    WUSE_AD(IDAM) = sum(WUSE_DD(IDAM,:))/NSTEPS     
  END DO
  CLOSE(NDAMFILE)

END SUBROUTINE READ_WATER_USE

END SUBROUTINE CMF_DAMOUT_INIT
!####################################################################



!####################################################################
SUBROUTINE CMF_DAMOUT_CALC
USE YOS_CMF_INPUT,      ONLY: DT
USE YOS_CMF_MAP,        ONLY: I1NEXT,   NSEQALL,  NSEQRIV
USE YOS_CMF_PROG,       ONLY: D2RIVOUT, D2FLDOUT, P2RIVSTO, P2FLDSTO
USE YOS_CMF_PROG,       ONLY: P2DAMSTO, P2DAMINF, D2RUNOFF   
USE YOS_CMF_DIAG,       ONLY: D2RIVINF, D2FLDINF
USE YOS_CMF_TIME,       ONLY: KSTEP, ISMM       

! local
IMPLICIT NONE
! SAVE for OMP
INTEGER(KIND=JPIM),SAVE    :: ISEQD
!** dam variables
REAL(KIND=JPRB),SAVE       :: DamVol
REAL(KIND=JPRB),SAVE       :: DamInflow
REAL(KIND=JPRB),SAVE       :: DamOutflw           !! Total outflw 
! REAL(KIND=JPRB),SAVE       :: DamOutTmp           !! Total outflw 

!*** water balance
REAL(KIND=JPRD),SAVE       :: GlbDAMSTO, GlbDAMSTONXT, GlbDAMINF, GlbDAMOUT, DamMiss

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
  DamOutflw = min( DamOutflw, DamVol/DT, real(P2RIVSTO(ISEQD,1)+P2FLDSTO(ISEQD,1),JPRB)/DT )
  DamOutflw = max( DamOutflw, 0._JPRB ) 

  !! update CaMa variables  (treat all outflow as RIVOUT in dam grid, no fldout)
  D2RIVOUT(ISEQD,1) = DamOutflw
  D2FLDOUT(ISEQD,1) = 0._JPRB
END DO
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
END DO
!$OMP END PARALLEL DO

DamMiss = GlbDAMSTO-GlbDAMSTONXT+GlbDAMINF-GlbDAMOUT
WRITE(LOGNAM,*) "CMF::DAM_CALC: DamMiss at all dams:", DamMiss*1.D-9


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
  REAL(KIND=JPRB), INTENT(IN)             :: Vol_dam       ! Current volume of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)             :: inflw         ! Inflow rate to the dam (m3/s)
  REAL(KIND=JPRB), INTENT(OUT)            :: outflw        ! Outflow rate from the dam (m3/s)
  !*** parameter
  CHARACTER(LEN=256), INTENT(IN)          :: dam_purpose   ! main use of dam
  REAL(KIND=JPRB), INTENT(IN)             :: Vol_tot       ! total storage capacity (m3)
  REAL(KIND=JPRB), INTENT(IN)             :: Q_n           ! normal discharge (m3/s)
  REAL(KIND=JPRB), INTENT(IN)             :: DPI           ! the ratio between annual mean demand and annual mean inflo (-)
  REAL(KIND=JPRB), INTENT(IN)             :: c             ! the ratio between capacity and mean annual inflow (-)
  REAL(KIND=JPRB), INTENT(IN)             :: dam_wuse      ! daily water demand (m3/s)
  REAL(KIND=JPRB), INTENT(IN)             :: dam_wuse_avg  ! average daily demand (m3/s)
  !*** local
  REAL(KIND=JPRB)                         :: a=0.85            !! adjustment factor for Krls
  REAL(KIND=JPRB)                         :: M=0.5             !! the ratio between minimum release and long‐term annual mean inflow
  REAL(KIND=JPRB)                         :: R             !! demand‐controlled release ratio = min(1, a*c)
  REAL(KIND=JPRB)                         :: Krls          !! the ratio between initial storage and the long‐term target storage = S/(a*C)
  REAL(KIND=JPRB)                         :: OutTmp        !! temporary outflow
  !=====================================================
  ! !! parameter
  ! M = 0.5
  ! a = 0.85 
  !! parameter calculation
  Krls =  Vol_dam / (Vol_tot * a)
  R    =  min(1., a * c)

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
  REAL(KIND=JPRB), INTENT(IN)             :: Vol_dam       ! Current volume of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)             :: inflw         ! Inflow rate to the dam (m3/s)
  REAL(KIND=JPRB), INTENT(OUT)            :: outflw        ! Outflow rate from the dam (m3/s)
  !*** parameter
  CHARACTER(LEN=256), INTENT(IN)          :: dam_purpose   ! main use of dam
  REAL(KIND=JPRB), INTENT(IN)             :: Vol_tot       ! total storage capacity (m3)
  REAL(KIND=JPRB), INTENT(IN)             :: Q_n           ! normal discharge (m3/s)
  REAL(KIND=JPRB), INTENT(IN)             :: DPI           ! the ratio between annual mean demand and annual mean inflo (-)
  REAL(KIND=JPRB), INTENT(IN)             :: c             ! the ratio between capacity and mean annual inflow (-)
  REAL(KIND=JPRB), INTENT(IN)             :: dam_wuse      ! daily water demand (m3/s)
  REAL(KIND=JPRB), INTENT(IN)             :: dam_wuse_avg  ! average daily demand (m3/s)
  INTEGER(KIND=JPIM), INTENT(IN)          :: MthStFC       ! the start of the flood control period
  INTEGER(KIND=JPIM), INTENT(IN)          :: MthNdFC       ! the end of the flood control period
  INTEGER(KIND=JPIM), INTENT(IN)          :: MthStOP       ! the start of the operational year
  !*** local
  REAL(KIND=JPRB)                         :: a=0.85            !! adjustment factor for Krls
  REAL(KIND=JPRB)                         :: M=0.5             !! the ratio between minimum release and long‐term annual mean inflow
  REAL(KIND=JPRB)                         :: R             !! demand‐controlled release ratio = min(1, a*c)
  REAL(KIND=JPRB)                         :: Krls          !! the ratio between initial storage and the long‐term target storage = S/(a*C)
  REAL(KIND=JPRB)                         :: OutTmp        !! temporary outflow
  REAL(KIND=JPRB)                         :: drop          !! temporary drop flow
  !=====================================================
  !! parameter calculation
  Krls =  Vol_dam / (Vol_tot * a)
  R    =  min(1., a * c)

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
    if ( MthStFC <= MthNdFC ) then
        if ( ISMM >= MthStFC .and. ISMM < MthNdFC) then
          if ( inflw < Q_n ) then
            drop = abs(inflw - Q_n)
            OutTmp = OutTmp + drop
          endif
        endif
    elseif ( MthStFC > MthNdFC ) then
        if ( ISMM >= MthStFC .or. ISMM < MthNdFC) then
          if ( inflw < Q_n ) then
            drop = abs(inflw - Q_n)
            OutTmp = OutTmp + drop
          endif
        endif
    endif

    !*** stage-2: from MthNdFC to MthStOP
    if ( MthNdFC <= MthStOP ) then
      if ( ISMM >= MthNdFC .and. ISMM < MthStOP ) then
          ! if ( inflw > Q_n ) then
          !   fill = abs(inflw - Q_n)
          ! end
          if (OutTmp > Q_n) then
            OutTmp = Q_n
          endif
      endif
    elseif ( MthNdFC > MthStOP ) then
      if ( ISMM >= MthNdFC .or. ISMM < MthStOP ) then
          ! if ( inflw > Q_n ) then
          !   fill = abs(inflw - Q_n)
          ! end
          if (OutTmp > Q_n) then
            OutTmp = Q_n
          endif
      endif
    endif

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
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_dam     ! Current volume of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)                :: inflw       ! Inflow rate to the dam (m3/s)
  REAL(KIND=JPRB), INTENT(OUT)               :: outflw      ! Outflow rate from the dam (m3/s)
  ! Parameters
  REAL(KIND=JPRB), INTENT(IN)  :: Vol_con          ! volume of water at the conservation level of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)  :: Vol_nor          ! volume of water at the normal level of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)  :: Vol_fld          ! volume of water at the flood controld level of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)  :: Vol_tot          ! total storage capacity (m3)
  REAL(KIND=JPRB), INTENT(IN)  :: Q_n              ! normal discharge (m3/s)
  REAL(KIND=JPRB), INTENT(IN)  :: Q_f              ! flood control discharge  (m3/s)

  ! Local variables
  REAL(KIND=JPRB)              :: F                ! fractional level of dam
  REAL(KIND=JPRB)              :: L_c              ! fractional level at conservation level
  REAL(KIND=JPRB)              :: L_n              ! fractional level at normal operating level
  REAL(KIND=JPRB)              :: L_f              ! fractional level at flood control level
  REAL(KIND=JPRB)              :: Q_min            ! minimum outflow rate
  REAL(KIND=JPRB)              :: Q_max            ! maximum outflow rate
  REAL(KIND=JPRB)              :: AdjLn            ! adjustment factor for L_n,in range 0.01 to 0.99
  REAL(KIND=JPRB)              :: AdjQnor          ! adjustment factor for Q_nor, in range 0.25 to 2
  REAL(KIND=JPRB)              :: Ladj_f           ! adjusted fractional level for maximum outflow
  REAL(KIND=JPRB)              :: Qadj_nor         ! adjusted normal outflow rate

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
  Qadj_nor = max(Q_min, min(AdjQnor * Q_n, Q_f)) 
  
  !*** operation scheme
  IF (F <= 2* L_c) THEN
    outflw = min(Q_min, Vol_dam/86400.)
  
  ELSEIF (F > 2* L_c .and. F <= L_n) THEN
    outflw = Q_min + (Qadj_nor - Q_min) * ((F - 2*L_c)/(L_n - 2*L_c))

  ELSEIF (F > L_n .and. F <= Ladj_f) THEN
    outflw = Qadj_nor 

  ELSEIF (F > Ladj_f .and. F <= L_f) THEN
    outflw = Qadj_nor + (Q_f - Qadj_nor) * ((F - Ladj_f)/(L_f - Ladj_f))
  
  ELSEIF (F > L_f) THEN
    Q_max = min(Q_f, max(1.2 * inflw, Qadj_nor))   
    outflw = max(Q_max, (F - L_f - 0.01) * (Vol_tot/86400.))
  ENDIF
  
  ! the condition described below is applied in order to prevent outflow values that are too large compared to the inflow value.
  ! reference: https://github.com/ec-jrc/lisflood-code
  IF (F < L_f .and. outflw > min(1.2 * inflw, Qadj_nor)) THEN
    outflw = min(outflw, max(inflw, Qadj_nor))
  ENDIF
  
END SUBROUTINE DAM_OPERATION_LISFLOOD


SUBROUTINE DAM_OPERATION_H22(Vol_dam, inflw, outflw, &
                             Vol_con, Vol_fld, Vol_eme, Q_n, Q_f, k)
  !! reference: Hanazaki, R., Yamazaki, D., & Yoshimura, K. (2022). 
  !! Development of a reservoir flood control scheme for global flood models. 
  !! Journal of Advances in Modeling Earth Systems, 14(3), e2021MS002944.
  
  IMPLICIT NONE
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_dam     ! Current volume of the dam (m3)
  REAL(KIND=JPRB), INTENT(IN)                :: inflw       ! Inflow rate to the dam (m3/s)
  REAL(KIND=JPRB), INTENT(OUT)               :: outflw      ! Outflow rate from the dam (m3/s)
  !*** parameter
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_con     ! conservative storage (m3) 
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_fld     ! flood control storage (m3) 
  REAL(KIND=JPRB), INTENT(IN)                :: Vol_eme     ! emergency storage (m3)   
  REAL(KIND=JPRB), INTENT(IN)                :: Q_n         ! normal discharge (m3/s)
  REAL(KIND=JPRB), INTENT(IN)                :: Q_f         ! flood discharge (m3/s)
  REAL(KIND=JPRB), INTENT(IN)                :: k           ! release coefficient
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
    outflw = max(inflw, Q_f)
  ENDIF

END SUBROUTINE DAM_OPERATION_H22


SUBROUTINE UPDATE_INFLOW
USE YOS_CMF_INPUT,      ONLY: PMINSLP, PMANFLD
USE YOS_CMF_MAP,        ONLY: D2RIVLEN, D2RIVMAN, D2ELEVTN, D2NXTDST, D2RIVWTH
USE YOS_CMF_PROG,       ONLY: D2RIVOUT_PRE, D2FLDOUT_PRE
USE YOS_CMF_DIAG,       ONLY: D2RIVDPH, D2RIVVEL, D2FLDDPH
IMPLICIT NONE
! SAVE for OpenMP
INTEGER(KIND=JPIM),SAVE    :: ISEQ, JSEQ
REAL(KIND=JPRB),SAVE       :: DSLOPE,DAREA,DVEL,DSLOPE_F,DARE_F,DVEL_F
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
END DO
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
END DO
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
    DSLOPE = max(DSLOPE,PMINSLP)

    DVEL   = D2RIVMAN(ISEQ,1)**(-1.) * DSLOPE**0.5 * D2RIVDPH(ISEQ,1)**(2./3.)
    DAREA  = D2RIVWTH(ISEQ,1) * D2RIVDPH(ISEQ,1)

    D2RIVVEL(ISEQ,1) = DVEL
    D2RIVOUT(ISEQ,1) = DAREA * DVEL
    D2RIVOUT(ISEQ,1) = MIN( D2RIVOUT(ISEQ,1), real(P2RIVSTO(ISEQ,1),JPRB)/DT )
    !=== floodplain flow
    DSLOPE_F = min( 0.005_JPRB,DSLOPE )    !! set min [instead of using weirequation for efficiency]
    DVEL_F   = PMANFLD**(-1.) * DSLOPE_F**0.5 * D2FLDDPH(ISEQ,1)**(2./3.)
    DARE_F   = P2FLDSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.)
    DARE_F   = MAX( DARE_F - D2FLDDPH(ISEQ,1)*D2RIVWTH(ISEQ,1), 0._JPRB )   !!remove above river channel     area

    D2FLDOUT(ISEQ,1) = DARE_F * DVEL_F
    D2FLDOUT(ISEQ,1) = MIN(  D2FLDOUT(ISEQ,1)*1._JPRD, P2FLDSTO(ISEQ,1)/DT )
  ENDIF
END DO
!$OMP END PARALLEL DO

END SUBROUTINE UPDATE_INFLOW
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE MODIFY_OUTFLW
! modify outflow in order to avoid negative storage
USE YOS_CMF_MAP,        ONLY: NSEQMAX
IMPLICIT NONE

REAL(KIND=JPRD)            :: P2STOOUT(NSEQMAX,1)                      !! total outflow from a grid     [m3]
REAL(KIND=JPRD)            :: P2RIVINF(NSEQMAX,1)                      !! 
REAL(KIND=JPRD)            :: P2FLDINF(NSEQMAX,1)                      !! 

REAL(KIND=JPRB)            :: D2RATE(NSEQMAX,1)                        !! outflow correction
! SAVE for OpenMP
INTEGER(KIND=JPIM),SAVE    :: ISEQ, JSEQ
REAL(KIND=JPRB),SAVE       :: OUT_R1, OUT_R2, OUT_F1, OUT_F2, DIUP, DIDW
!$OMP THREADPRIVATE     (JSEQ,OUT_R1, OUT_R2, OUT_F1, OUT_F2, DIUP, DIDW)
!================================================
  
!*** 1. initialize & calculate P2STOOUT for normal cells

!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
  P2RIVINF(ISEQ,1) = 0._JPRD
  P2FLDINF(ISEQ,1) = 0._JPRD
  P2STOOUT(ISEQ,1) = 0._JPRD
  D2RATE(ISEQ,1) = 1._JPRB
END DO
!$OMP END PARALLEL DO

!! for normal cells ---------
#ifndef NoAtom_CMF
!$OMP PARALLEL DO
#endif
DO ISEQ=1, NSEQRIV                                                    !! for normalcells
  JSEQ=I1NEXT(ISEQ) ! next cell's pixel
  OUT_R1 = max(  D2RIVOUT(ISEQ,1),0._JPRB )
  OUT_R2 = max( -D2RIVOUT(ISEQ,1),0._JPRB )
  OUT_F1 = max(  D2FLDOUT(ISEQ,1),0._JPRB )
  OUT_F2 = max( -D2FLDOUT(ISEQ,1),0._JPRB )
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
END DO
#ifndef NoAtom_CMF
!$OMP END PARALLEL DO
#endif

!! for river mouth grids ------------
!$OMP PARALLEL DO
DO ISEQ=NSEQRIV+1, NSEQALL
  OUT_R1 = max( D2RIVOUT(ISEQ,1), 0._JPRB )
  OUT_F1 = max( D2FLDOUT(ISEQ,1), 0._JPRB )
  P2STOOUT(ISEQ,1) = P2STOOUT(ISEQ,1) + OUT_R1*DT + OUT_F1*DT
END DO
!$OMP END PARALLEL DO

!============================
!*** 2. modify outflow

!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
  IF ( P2STOOUT(ISEQ,1) > 1.E-8 ) THEN
    D2RATE(ISEQ,1) = min( (P2RIVSTO(ISEQ,1)+P2FLDSTO(ISEQ,1)) * P2STOOUT(ISEQ,1)**(-1.), 1._JPRD )
  ENDIF
END DO
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
END DO
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
END DO
!$OMP END PARALLEL DO

END SUBROUTINE MODIFY_OUTFLW
!==========================================================

END SUBROUTINE CMF_DAMOUT_CALC
!####################################################################



!####################################################################
SUBROUTINE CMF_DAMOUT_WRTE
USE YOS_CMF_PROG,       ONLY: P2DAMSTO, P2DAMINF, D2RIVOUT
USE YOS_CMF_TIME,       ONLY: IYYYYMMDD,ISYYYY
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID

! local
CHARACTER(len=36)          :: WriteTXT(NDAMX) !, WriteTXT2(NDAMX)

! File IO
INTEGER(KIND=JPIM),SAVE    :: ISEQD, JDAM
INTEGER(KIND=JPIM),SAVE    :: LOGDAM
CHARACTER(len=4),SAVE      :: CYYYY
CHARACTER(len=256),SAVE    :: CLEN, CFMT
CHARACTER(len=256),SAVE    :: DAMTXT
LOGICAL,SAVE               :: IsOpen
DATA IsOpen       /.FALSE./

! ======

IF( LDAMTXT .and. LDAMTXT)THEN

  IF( .not. IsOpen)THEN
    IsOpen=.TRUE.
    WRITE(CYYYY,'(i4.4)') ISYYYY
    DAMTXT='./damtxt-'//trim(CYYYY)//'.txt'

    LOGDAM=INQUIRE_FID()
    OPEN(LOGDAM,FILE=DAMTXT,FORM='formatted')

    WRITE(CLEN,'(i0)') NDAMX
    CFMT="(i10,"//TRIM(CLEN)//"(a36))"

    JDAM=0
    DO IDAM=1, NDAM
      IF( DamSeq(IDAM)<=0 ) CYCLE
      JDAM=JDAM+1
      ISEQD=DamSeq(IDAM)

      ! WRITE(WriteTxt(JDAM), '(i12,2f12.2)') GRanD_ID(IDAM), (FldVol(IDAM)+ConVol(IDAM))*1.E-9, ConVol(IDAM)*1.E-9
      ! WRITE(WriteTxt2(JDAM),'(3f12.2)') upreal(IDAM),   Qf(IDAM), Qn(IDAM)
      WRITE(WriteTxt(JDAM), '(i12)') GRanD_ID(IDAM)
    END DO

    WRITE(LOGDAM,CFMT) NDAMX, (WriteTXT(JDAM) ,JDAM=1, NDAMX)

    ! CFMT="(a10,"//TRIM(CLEN)//"(a36))"
    ! WRITE(LOGDAM,CFMT)  "Date", (WriteTXT2(JDAM),JDAM=1, NDAMX)
  ENDIF

  JDAM=0
  DO IDAM=1, NDAM
    IF( DamSeq(IDAM)<=0 ) CYCLE
    JDAM=JDAM+1
    ISEQD=DamSeq(IDAM)
    ! P2DAMSTO m3 to Million Cubic Meter
    WRITE(WriteTxt(JDAM), '(3f12.3)') P2DAMSTO(ISEQD,1)*1.E-6, P2DAMINF(ISEQD,1), D2RIVOUT(ISEQD,1)
  END DO

  CFMT="(i10,"//TRIM(CLEN)//"(a36))"  
  WRITE(LOGDAM,CFMT) IYYYYMMDD, (WriteTXT(JDAM),JDAM=1, NDAMX)

ENDIF

END SUBROUTINE CMF_DAMOUT_WRTE
!####################################################################


END MODULE CMF_CTRL_DAMOUT_MOD
