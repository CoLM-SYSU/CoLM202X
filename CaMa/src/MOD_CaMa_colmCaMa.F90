#include <define.h>
MODULE MOD_CaMa_colmCaMa
#if(defined CaMa_Flood)
!DESCRIPTION
!===========
! This MODULE is the coupler for the colm and CaMa-Flood model.

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------
!* :SUBROUTINE:"colm_CaMa_init" :  Initialization of the coupler
!* :SUBROUTINE:"colm_CaMa_drv"  :  Coupling between colm and CaMa-Flood
!* :SUBROUTINE:"colm_CaMa_exit" :  Finalization of the coupler
!* :SUBROUTINE:"get_fldinfo"    :  Get floodplain information from CaMa-Flood model
!* :SUBROUTINE:"get_fldevp"     :  calculate floodplain evaporation

!REVISION HISTORY
!----------------
! 2023.02.21  Zhongwang Wei @ SYSU
! 2021.12.02  Zhongwang Wei @ SYSU
! 2020.10.01  Zhongwang Wei @ SYSU

use MOD_Namelist
USE MOD_CaMa_Vars
USE PARKIND1,                  ONLY: JPRB, JPRM, JPIM
USE CMF_DRV_CONTROL_MOD,       ONLY: CMF_DRV_INPUT,   CMF_DRV_INIT,    CMF_DRV_END
USE CMF_DRV_ADVANCE_MOD,       ONLY: CMF_DRV_ADVANCE
USE CMF_CTRL_FORCING_MOD,      ONLY: CMF_FORCING_GET, CMF_FORCING_PUT
USE CMF_CTRL_OUTPUT_MOD,       ONLY: CMF_OUTPUT_INIT,CMF_OUTPUT_END,NVARSOUT,VAROUT
USE YOS_CMF_INPUT,             ONLY: NXIN, NYIN, DT,DTIN,IFRQ_INP,LLEAPYR,NX,NY,RMIS,DMIS
USE MOD_Precision,                 ONLY: r8,r4
USE YOS_CMF_INPUT ,            ONLY: LROSPLIT,LWEVAP,LWINFILT
USE MOD_SPMD_Task
USE CMF_CTRL_TIME_MOD
USE MOD_Vars_Global,                ONLY : spval
USE MOD_Vars_1DFluxes
USE MOD_Qsadv

IMPLICIT NONE
!----------------------- Dummy argument --------------------------------
INTEGER I,J
INTEGER(KIND=JPIM)              :: ISTEPX              ! total time step
INTEGER(KIND=JPIM)              :: ISTEPADV            ! time step to be advanced within DRV_ADVANCE
REAL(KIND=JPRB),ALLOCATABLE     :: ZBUFF(:,:,:)        ! Buffer to store forcing runoff

INTERFACE colm_CaMa_init
   MODULE PROCEDURE colm_CaMa_init
END INTERFACE

INTERFACE colm_CaMa_drv
   MODULE PROCEDURE colm_CaMa_drv
END INTERFACE

   INTERFACE colm_CaMa_exit
      MODULE PROCEDURE colm_CaMa_exit
   end INTERFACE
CONTAINS

SUBROUTINE colm_CaMa_init
   USE MOD_LandPatch
   USE YOS_CMF_TIME,          ONLY: YYYY0
   IMPLICIT NONE
   !** local variables

   INTEGER i,j
   INTEGER(KIND=JPIM)          :: JF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF(p_is_master)THEN
      !Namelist handling
      CALL CMF_DRV_INPUT
      !get the time information from colm namelist
      DT       = IFRQ_INP*3600                                              ! time step of model simulation [sec]
      DTIN     = IFRQ_INP*3600                                              ! time step of input data [sec]
      SYEAR    = DEF_simulation_time%start_year                             ! start year
      SMON     = DEF_simulation_time%start_month                            ! start month
      SDAY     = DEF_simulation_time%start_day                              ! start day
      SHOUR    = DEF_simulation_time%start_sec/3600                         ! start hour
      EYEAR    = DEF_simulation_time%end_year                               ! end year
      EMON     = DEF_simulation_time%end_month                              ! end month
      EDAY     = DEF_simulation_time%end_day                                ! end day
      EHOUR    = DEF_simulation_time%end_sec/3600                           ! end hour
      LLEAPYR  = DEF_forcing%leapyear                                       ! leap year flag

      !----------------------- Dummy argument --------------------------------
      YYYY0    = SYEAR
      RMIS     = spval
      DMIS     = spval

      CALL CMF_DRV_INIT       !INITIALIZATION

      !Initialize varialbes to be outputed from variable list
      DO JF=1,NVARSOUT
         SELECT CASE (VAROUT(JF)%CVNAME)
         CASE ('rivout') ! river discharge [m3/s]
            DEF_hist_cama_vars%rivout=.true.
         CASE ('rivsto') ! river storage  [m3]
            DEF_hist_cama_vars%rivsto=.true.
         CASE ('rivdph') ! river depth   [m]
            DEF_hist_cama_vars%rivdph=.true.
         CASE ('rivvel') ! river velocity [m/s]
            DEF_hist_cama_vars%rivvel=.true.
         CASE ('fldout') ! floodplain discharge [m3/s]
            DEF_hist_cama_vars%fldout=.true.
         CASE ('fldsto') ! floodplain storage [m3]
            DEF_hist_cama_vars%fldsto=.true.
         CASE ('flddph') ! floodplain depth [m]
            DEF_hist_cama_vars%flddph=.true.
         CASE ('fldfrc') ! floodplain fraction (0-1)
            DEF_hist_cama_vars%fldfrc=.true.
         CASE ('fldare') ! floodplain area [m2]
            DEF_hist_cama_vars%fldare=.true.
         CASE ('sfcelv') ! water surface elevation [m]
            DEF_hist_cama_vars%sfcelv=.true.
         CASE ('totout') ! total discharge(river+floodplain) [m3/s]
            DEF_hist_cama_vars%totout=.true.
         CASE ('outflw') ! compatibility for previous file name =totout
            DEF_hist_cama_vars%outflw=.true.
         CASE ('totsto') ! total storage(river+floodplain) [m3]
            DEF_hist_cama_vars%totsto=.true.
         CASE ('storge') ! compatibility for previous file name =totsto
            DEF_hist_cama_vars%storge=.true.
         CASE ('pthflw') ! bifurcation channel discharge [m3/s]
            DEF_hist_cama_vars%pthflw=.true.
         CASE ('pthout') ! net bifurcation discharge [m3/s]
            DEF_hist_cama_vars%pthout=.true.
         CASE ('maxflw') ! daily maximum discharge [m3/s]
            DEF_hist_cama_vars%maxflw=.true.
         CASE ('maxdph') ! daily maximum depth [m]
            DEF_hist_cama_vars%maxdph=.true.
         CASE ('maxsto') ! daily maximum storage [m3]
            DEF_hist_cama_vars%maxsto=.true.
         !TODO: check the difference between gwsto and gdwsto
         CASE ('gwsto')  ! ground water storage [m3]
            DEF_hist_cama_vars%gwsto=.true.
         CASE ('gdwsto') ! ground water storage [m3]
            DEF_hist_cama_vars%gdwsto=.true.
         CASE ('gwout')  ! ground water discharge [m3/s]
            DEF_hist_cama_vars%gwout=.true.
         !TODO: check the difinition of gdwrtn, runoff, rofsfc, rofsub
         CASE ('gdwrtn') ! Ground water return flow [m3/s]
            DEF_hist_cama_vars%gdwrtn=.true.
         CASE ('runoff') ! total runoff [m3/s]                !!  compatibility for previous file name
            DEF_hist_cama_vars%runoff=.true.
         CASE ('runoffsub') ! subsurface runoff [m3/s]        !!  compatibility for previous file name
            DEF_hist_cama_vars%runoffsub=.true.
         CASE ('rofsfc') ! surface runoff [m3/s]              !!  compatibility for previous file name
            DEF_hist_cama_vars%rofsfc=.true.
         CASE ('rofsub')  ! input sub-surface runoff [m3/s]
            DEF_hist_cama_vars%rofsub=.true.
         CASE ('damsto')   ! reservoir storage [m3]
            DEF_hist_cama_vars%damsto=.true.
         CASE ('daminf')   ! reservoir inflow [m3/s]
               DEF_hist_cama_vars%daminf=.true.
         CASE ('levsto')   !flood storage in protected side (storage betwen river & levee) [m3]
               DEF_hist_cama_vars%levsto=.true.
         CASE ('levdph')   !flood depth in protected side [m]
               DEF_hist_cama_vars%levdph=.true.
         CASE ('wevap')    ! input inundation Evaporation [m]
            IF (LWEVAP) then
               DEF_hist_cama_vars%wevap=.true.
            ELSE
               DEF_hist_cama_vars%wevap=.false.
            ENDIF
         CASE ('winfilt')  ! input inundation re-infiltrition [m]
            IF (LWINFILT) then
               DEF_hist_cama_vars%winfilt=.true.
            else
               DEF_hist_cama_vars%winfilt=.false.
            ENDIF
         CASE DEFAULT
            stop
         END SELECT
      end do
   ENDIF

      !Broadcast the variables to all the processors
      CALL mpi_bcast (NX      ,   1, MPI_INTEGER,   p_root, p_comm_glb, p_err) ! number of grid points in x-direction of CaMa-Flood
      CALL mpi_bcast (NY      ,   1, MPI_INTEGER,   p_root, p_comm_glb, p_err) ! number of grid points in y-direction of CaMa-Flood
      CALL mpi_bcast (IFRQ_INP ,   1, MPI_INTEGER,  p_root, p_comm_glb, p_err) ! input frequency of CaMa-Flood (hour)
      CALL mpi_bcast (LWEVAP ,   1, MPI_LOGICAL,  p_root, p_comm_glb, p_err)   ! switch for inundation evaporation
      CALL mpi_bcast (LWINFILT ,   1, MPI_LOGICAL,  p_root, p_comm_glb, p_err) ! switch for inundation re-infiltration

      !Allocate the data structure for cama
      CALL gcama%define_by_ndims (NX, NY)  !define the data structure for cama
      CALL mp2g_cama%build (landpatch, gcama) !build the mapping between cama and mpi
      CALL mg2p_cama%build (gcama, landpatch)

      CALL cama_gather%set (gcama)

      !Allocate the cama-flood related variable for accumulation
      CALL allocate_2D_cama_Fluxes  (gcama) !allocate the 2D variables
      CALL allocate_acc_cama_Fluxes () !allocate the accumulation variables
      CALL FLUSH_acc_cama_fluxes    () !initialize the accumulation variables

      !Only master processor allocate the 2D variables
      IF (p_is_master) THEN
         ALLOCATE (runoff_2d (NX,NY))
         ALLOCATE (fevpg_2d  (NX,NY))
         ALLOCATE (finfg_2d  (NX,NY))
         !Allocate data buffer for input forcing, flood fraction and flood depth
         ALLOCATE (ZBUFF(NX,NY,4))
         ALLOCATE (fldfrc_tmp(NX,NY))
         ALLOCATE (flddepth_tmp(NX,NY))
         !Initialize the data buffer for input forcing, flood fraction and flood depth
         runoff_2d(:,:)    = 0.0D0 !runoff in master processor
         fevpg_2d(:,:)     = 0.0D0 !evaporation in master processor
         finfg_2d(:,:)     = 0.0D0 !re-infiltration in master processor
         ZBUFF(:,:,:)      = 0.0D0 !input forcing in master processor
         fldfrc_tmp(:,:)   = 0.0D0 !flood fraction in master processor
         flddepth_tmp(:,:) = 0.0D0 !flood depth in master processor
      ENDIF
      !Allocate the cama-flood related variable in worker processors
      IF (p_is_worker) THEN
         ALLOCATE (flddepth_cama(numpatch)) !flood depth in worker processors
         ALLOCATE (fldfrc_cama(numpatch))   !flood fraction in worker processors
         ALLOCATE (fevpg_fld(numpatch))     !evaporation in worker processors
         ALLOCATE (finfg_fld(numpatch))     !re-infiltration in worker processors
         flddepth_cama(:)     =  0.0D0
         fldfrc_cama(:)       =  0.0D0
         fevpg_fld(:)         =  0.0D0
         finfg_fld(:)         =  0.0D0
end IF
end SUBROUTINE colm_CaMa_init

!####################################################################
SUBROUTINE colm_cama_drv(idate_sec)
   IMPLICIT NONE
   INTEGER, intent(in)  :: idate_sec     ! calendar (year, julian day, seconds)
   !Accumulate cama-flood related flux variables
   CALL accumulate_cama_fluxes
   ! If the time is the same as the input time step of cama-flood
   IF  (MOD(idate_sec,3600*int(IFRQ_INP))==0) THEN
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      ! Prepare sending the accumulated runoff flux varilble to cama model (master processor to worker processors)
      CALL colm2cama_real8 (a_rnof_cama, f_rnof_cama, runoff_2d)

      ! Prepare sending the accumulated inundation evaporation flux to cama model (master processor to worker processors)
      ! only if the inundation evaporation is turned on
      IF (LWEVAP) THEN
#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
         CALL colm2cama_real8 (a_fevpg_fld, f_fevpg_fld, fevpg_2d)
      ENDIF
      ! Prepare sending the accumulated inundation re-infiltrition flux to cama model (master processor to worker processors)
      ! only if the inundation re-infiltrition is turned on
      IF (LWINFILT) THEN
#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
         CALL colm2cama_real8 (a_finfg_fld, f_finfg_fld, finfg_2d)
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      ! Reset the accumulation variables
      CALL flush_acc_cama_fluxes
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      ! Initialize the variables unit for cama-flood input
      IF(p_is_master)THEN
         DO i = 1,NX           ! cama_gather%ginfo%nlon
            DO j = 1, NY       ! cama_gather%ginfo%nlat
               ZBUFF(i,j,1)=runoff_2d(i,j)/1000.0D0   ! mm/s -->m/s
               ZBUFF(i,j,2)=0.0D0
               IF (LWEVAP) THEN
                  ZBUFF(i,j,3)=fevpg_2d(i,j)/1000.0D0 ! mm/s -->m/s
               ELSE
                  ZBUFF(i,j,3)=0.0D0
               ENDIF
               IF (LWINFILT) THEN
                  ZBUFF(i,j,4)=finfg_2d(i,j)/1000.0D0  !mm/s -->m/s
               ELSE
                  ZBUFF(i,j,4)=0.0D0
               ENDIF
            ENDDO
         ENDDO

         ! Simulating the hydrodynamics in continental-scale rivers
         ! ----------------------------------------------------------------------

         ! Get the time step of cama-flood simulation
         ISTEPADV=INT(DTIN/DT,JPIM)
         ! Interporlate variables & send to CaMa-Flood
         CALL CMF_FORCING_PUT(ZBUFF)
         ! Advance CaMa-Flood model for ISTEPADV
         CALL CMF_DRV_ADVANCE(ISTEPADV)
         ! Get the flood depth and flood fraction from cama-flood model
         IF (LWINFILT .or. LWEVAP) THEN
            CALL get_fldinfo()
         ENDIF
      ENDIF
      ! Send the flood depth and flood fraction from master processors to worker processors
      IF (LWINFILT .or. LWEVAP) THEN
#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
         CALL cama2colm_real8 (flddepth_tmp, f_flddepth_cama, flddepth_cama)! unit [m]
#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
         CALL cama2colm_real8 (fldfrc_tmp,   f_fldfrc_cama,   fldfrc_cama  ) ! unit [%]
#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
         flddepth_cama=flddepth_cama*1000.D0 !m --> mm
         fldfrc_cama=fldfrc_cama/100.D0     !% --> [0-1]
      ENDIF
   ENDIF
END SUBROUTINE colm_cama_drv

SUBROUTINE colm_cama_exit
   IMPLICIT NONE
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   ! finalize CaMa-Flood
   CALL deallocate_acc_cama_Fluxes ()
   IF(p_is_master)THEN
      ! finalize CaMa-Flood
      DEALLOCATE(ZBUFF)
      DEALLOCATE (runoff_2d)
      DEALLOCATE (fevpg_2d)
      DEALLOCATE (finfg_2d)
   ENDIF
   IF (p_is_worker) THEN
      DEALLOCATE (flddepth_cama)
      DEALLOCATE (fldfrc_cama)
      DEALLOCATE (fevpg_fld)
      DEALLOCATE (finfg_fld)
   end IF
END SUBROUTINE colm_cama_exit

!####################################################################
SUBROUTINE get_fldinfo()
   !DESCRIPTION
   !===========
   ! This subrountine prepare cama output variables for inundation evaporation and re-infiltration

   !ANCILLARY FUNCTIONS AND SUBROUTINES
   !-------------------
   !* :SUBROUTINE:"vecD2mapD" : convert 1D vector data -> 2D map data (REAL*8), see CAMA/CMF_UTILS_MOD.F90


   !REVISION HISTORY
   !----------------
   ! 2020.10.01  Zhongwang Wei @ SYSU

   USE YOS_CMF_INPUT,      ONLY:  NX, NY !grid number
   USE YOS_CMF_DIAG,       ONLY:  D2FLDDPH,D2FLDFRC !1D vector data of flood depth and flood fraction
   USE CMF_UTILS_MOD,      ONLY:  vecD2mapD          !convert 1D vector data -> 2D map data (REAL*8)

   IMPLICIT NONE

   !----------------------- Dummy argument --------------------------------
   INTEGER i,j

   !================================================
   !! convert 1Dvector to 2Dmap
   CALL vecD2mapD(D2FLDFRC,flddepth_tmp)             !! MPI node data is gathered by VEC2MAP
   CALL vecD2mapD(D2FLDDPH,fldfrc_tmp)               !! MPI node data is gathered by VEC2MAP

   do i    = 1, NX
      do j = 1, NY
         IF (flddepth_tmp(i,j) .LT.    0.0)        flddepth_tmp(i,j) = 0.0
         IF (fldfrc_tmp(i,j)   .LT.    0.0)        fldfrc_tmp(i,j)   = 0.0
         IF (fldfrc_tmp(i,j)   .GT.    100.0)      fldfrc_tmp(i,j)   = 100.0    !!If fraction is larger than 100%, it is set to 100%.
      ENDDO
   ENDDO
END SUBROUTINE get_fldinfo

SUBROUTINE get_fldevp (hu,ht,hq,us,vs,tm,qm,rhoair,psrf,tssea,&
	hpbl, &
    taux,tauy,fseng,fevpg,tref,qref,z0m,zol,rib,ustar,qstar,tstar,fm,fh,fq)
   !DESCRIPTION
   !===========
   ! This subrountine compute surface fluxes, derviatives, and exchange coefficiants
   ! This is the main SUBROUTINE to execute the calculation of thermal processes
   ! and surface fluxes

   !ANCILLARY FUNCTIONS AND SUBROUTINES
   !-------------------
   !* :SUBROUTINE:"qsadv"         :     !  computes saturation mixing ratio and change in saturation
   !* :SUBROUTINE:" moninobukini" :     !  initialzation of Monin-Obukhov length, see MOD_FrictionVelocity.F90
   !* :SUBROUTINE:" moninobuk"    :     !  calculation of friction velocity, relation for potential temperature
                                        !  and humidity profiles of surface boundary layer,see MOD_FrictionVelocity.F90
   !REVISION HISTORY
   !----------------
   ! 2023.05.05  Shaofeng Liu @ SYSU:
   !             add option to call moninobuk_leddy, the LargeEddy
   !             surface turbulence scheme (LZD2022); make a proper update of um.
   ! 2020.10.01  Zhongwang Wei @ SYSU
   ! 2002.08.30  Yongjiu Dai   @ BNU
   ! 1999.09.15  Yongjiu Dai   @ BNU
   USE MOD_Precision
   USE MOD_Const_Physical, ONLY : cpair,rgas,vonkar,grav
   USE MOD_FrictionVelocity
   USE MOD_TurbulenceLEddy
   IMPLICIT NONE
   REAL(r8), INTENT(in)  :: hu      ! agcm reference height of wind [m]
   REAL(r8), INTENT(in)  :: ht      ! agcm reference height of temperature [m]
   REAL(r8), INTENT(in)  :: hq      ! agcm reference height of humidity [m]
   REAL(r8), INTENT(in)  :: us      ! wind component in eastward direction [m/s]
   REAL(r8), INTENT(in)  :: vs      ! wind component in northward direction [m/s]
   REAL(r8), INTENT(in)  :: tm      ! temperature at agcm reference height [kelvin]
   REAL(r8), INTENT(in)  :: qm      ! specific humidity at agcm reference height [kg/kg]
   REAL(r8), INTENT(in)  :: rhoair  ! density air [kg/m3]
   REAL(r8), INTENT(in)  :: psrf    ! atmosphere pressure at the surface [pa] [not used]
   REAL(r8), INTENT(in)  :: tssea   ! inundation surface temperature [K]-->set to tgrnd
   REAL(r8), INTENT(in)  :: hpbl    ! atmospheric boundary layer height [m]
   REAL(r8), INTENT(out) :: taux    ! wind stress: E-W [kg/m/s**2]
   REAL(r8), INTENT(out) :: tauy    ! wind stress: N-S [kg/m/s**2]
   REAL(r8), INTENT(out) :: fseng   ! sensible heat flux from ground [mm/s]
   REAL(r8), INTENT(out) :: fevpg   ! evaporation heat flux from ground [mm/s]
   REAL(r8), INTENT(out) :: tref    ! 2 m height air temperature [kelvin]
   REAL(r8), INTENT(out) :: qref    ! 2 m height air humidity [?]
   REAL(r8), INTENT(out) :: z0m     ! effective roughness [m]
   REAL(r8), INTENT(out) :: zol     ! dimensionless height (z/L) used in Monin-Obukhov theory
   REAL(r8), INTENT(out) :: rib     ! bulk Richardson number in surface layer
   REAL(r8), INTENT(out) :: ustar   ! friction velocity [m/s]
   REAL(r8), INTENT(out) :: tstar   ! temperature scaling parameter
   REAL(r8), INTENT(out) :: qstar   ! moisture scaling parameter
   REAL(r8), INTENT(out) :: fm      ! integral of profile function for momentum
   REAL(r8), INTENT(out) :: fh      ! integral of profile function for heat
   REAL(r8), INTENT(out) :: fq      ! integral of profile function for moisture

   !----------------------- Dummy argument --------------------------------
   INTEGER i
   INTEGER niters            ! maximum number of iterations for surface temperature
   INTEGER iter              ! iteration index
   INTEGER nmozsgn           ! number of times moz changes sign

   REAL(r8) :: beta          ! coefficient of conective velocity [-]
   REAL(r8) :: displax       ! zero-displacement height [m]
   REAL(r8) :: dth           ! diff of virtual temp. between ref. height and surface
   REAL(r8) :: dqh           ! diff of humidity between ref. height and surface
   REAL(r8) :: dthv          ! diff of vir. poten. temp. between ref. height and surface
   REAL(r8) :: eg            ! water vapor pressure at temperature T [Pa]
   REAL(r8) :: degdT         ! d(eg)/dT
   REAL(r8) :: obu           ! monin-obukhov length [m]
   REAL(r8) :: obuold        ! monin-obukhov length from previous iteration
   REAL(r8) :: qsatg         ! ground saturated specific humidity [kg/kg]
   REAL(r8) :: qsatgdT       ! d(qsatg)/dT
   REAL(r8) :: ram           ! aerodynamical resistance [s/m]
   REAL(r8) :: rah           ! thermal resistance [s/m]
   REAL(r8) :: raw           ! moisture resistance [s/m]
   REAL(r8) :: raih          ! temporary variable [kg/m2/s]
   REAL(r8) :: raiw          ! temporary variable [kg/m2/s]
   REAL(r8) :: fh2m          ! relation for temperature at 2m
   REAL(r8) :: fq2m          ! relation for specific humidity at 2m
   REAL(r8) :: fm10m         ! integral of profile function for momentum at 10m
   REAL(r8) :: thm           ! intermediate variable (tm+0.0098*ht)
   REAL(r8) :: th            ! potential temperature (kelvin)
   REAL(r8) :: thv           ! virtual potential temperature (kelvin)
   REAL(r8) :: thvstar       ! virtual potential temperature scaling parameter
   REAL(r8) :: um            ! wind speed including the stablity effect [m/s]
   REAL(r8) :: ur            ! wind speed at reference height [m/s]
   REAL(r8) :: visa          ! kinematic viscosity of dry air [m2/s]
   REAL(r8) :: wc            ! convective velocity [m/s]
   REAL(r8) :: wc2           ! wc**2
   REAL(r8) :: xt            !  ---->  temporary variables
   REAL(r8) :: xq            !  ---->  temporary variables
   REAL(r8) :: zii           ! convective boundary height [m]
   REAL(r8) :: zldis         ! reference height "minus" zero displacement heght [m]
   REAL(r8) :: z0mg          ! roughness length over ground, momentum [m]
   REAL(r8) :: z0hg          ! roughness length over ground, sensible heat [m]
   REAL(r8) :: z0qg          ! roughness length over ground, latent heat [m]

   REAL, parameter :: zsice = 0.04  ! sea ice aerodynamic roughness length [m]

   !-----------------------------------------------------------------------
   ! Potential temperatur at the reference height
   beta = 1.      ! -  (in computing W_*)
   zii  = 1000.    ! m  (pbl height)

   !-----------------------------------------------------------------------
   ! Compute sensible and latent fluxes and their derivatives with respect
   ! to surface temperature using surface temperatures from previous time step.
   !-----------------------------------------------------------------------
   ! Initialization variables
   nmozsgn = 0
   obuold  = 0.
   ! Calculate saturation mixing ratio and change in saturation
   CALL qsadv(tssea,psrf,eg,degdT,qsatg,qsatgdT)

   ! Potential temperatur at the reference height
   thm = tm + 0.0098*ht                   ! intermediate variable equivalent to
                                          ! tm*(pgcm/psrf)**(rgas/cpair)
   th = tm*(100000./psrf)**(rgas/cpair)   ! potential T
   thv = th*(1.+0.61*qm)                  ! virtual potential T
   ur = max(0.1,sqrt(us*us+vs*vs))        ! limit set to 0.1

   dth   = thm-tssea                      ! diff of potential temp. between ref. height and surface
   dqh   = qm-qsatg                       ! diff of humidity between ref. height and surface
   dthv  = dth*(1.+0.61*qm)+0.61*th*dqh   ! diff of vir. poten. temp. between ref. height and surface
   !TODO: check if this is correct, inundation may occur over vegetated surface
   zldis = hu-0.                          ! reference height "minus" zero displacement heght

   ! Kinematic viscosity of dry air (m2/s)- Andreas (1989) CRREL Rep. 89-11
   visa=1.326e-5*(1.+6.542e-3*tm + 8.301e-6*tm**2 - 4.84e-9*tm**3)

   ! Loop to obtain initial and good ustar and zo
   ustar=0.06    ! initial value of ustar
   wc=0.5        ! initial value of wc
   !initial value of um
   IF(dthv.ge.0.) THEN
   um=max(ur,0.1)
   else
   um=sqrt(ur*ur+wc*wc)
   ENDIF
   ! initial value of z0mg and ustar
   do i=1,5
   z0mg=0.013*ustar*ustar/grav+0.11*visa/ustar
   ustar=vonkar*um/log(zldis/z0mg)
   ENDDO
   !
   CALL moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu)

   ! Evaluated stability-dependent variables using moz from prior iteration
   niters  = 10
   displax = 0.

   !----------------------------------------------------------------
   ITERATION : do iter = 1, niters         ! begin stability iteration
   !----------------------------------------------------------------
   ! Compute stability-dependent variables
   z0mg = 0.013*ustar*ustar/grav + 0.11*visa/ustar
   xq   = 2.67*(ustar*z0mg/visa)**0.25 - 2.57
   xt= xq
   z0qg=z0mg/exp(xq)
   z0hg=z0mg/exp(xt)
   !calculation of friction velocity, relation for potential temperature
   if (DEF_USE_CBL_HEIGHT) then
      CALL moninobuk_leddy(hu,ht,hq,displax,z0mg,z0hg,z0qg,obu,um, hpbl, &
                           ustar,fh2m,fq2m,fm10m,fm,fh,fq)
   else
      CALL moninobuk(hu,ht,hq,displax,z0mg,z0hg,z0qg,obu,um,&
                     ustar,fh2m,fq2m,fm10m,fm,fh,fq)
   endif
   !get qstar and tstar
   tstar = vonkar/fh*dth
   qstar = vonkar/fq*dqh

   thvstar=tstar*(1.+0.61*qm)+0.61*th*qstar
   zol=zldis*vonkar*grav*thvstar/(ustar**2*thv) !z/L
   IF(zol >= 0.) THEN       ! stable
   zol = min(2.,max(zol,1.e-6))
   else                     ! unstable
   zol = max(-100.,min(zol,-1.e-6))
   ENDIF
   obu = zldis/zol

   IF(zol >= 0.)THEN
   um = max(ur,0.1) !wind speed at reference height
   else
   if (DEF_USE_CBL_HEIGHT) then !//TODO: Shaofeng, 2023.05.18
     zii = max(5.*hu,hpbl)
   endif !//TODO: Shaofeng, 2023.05.18
   wc = (-grav*ustar*thvstar*zii/thv)**(1./3.) !convective velocity scale
   wc2 = beta*beta*(wc*wc)                     !convective velocity scale squared
   um = sqrt(ur*ur+wc2)                        !wind speed with convective velocity scale
   ENDIF

   IF (obuold*obu < 0.) nmozsgn = nmozsgn+1
   IF(nmozsgn >= 4) EXIT

   obuold = obu

   !----------------------------------------------------------------
   ENDDO ITERATION                         ! end stability iteration
   !----------------------------------------------------------------

   ! Get derivative of fluxes with repect to ground temperature
   ram    = 1./(ustar*ustar/um)
   rah    = 1./(vonkar/fh*ustar)
   raw    = 1./(vonkar/fq*ustar)

   raih   = rhoair*cpair/rah
   raiw   = rhoair/raw
   !cgrnds = raih
   !cgrndl = raiw*qsatgdT

   rib = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

   ! Surface fluxes of momentum, sensible and latent
   ! using ground temperatures from previous time step
   taux   = -rhoair*us/ram
   tauy   = -rhoair*vs/ram

   fseng  = -raih*dth
   fevpg  = -raiw*dqh
   !fsena  = fseng
   !fevpa  = fevpg

   ! 2 m height air temperature
   tref   = thm + vonkar/fh*dth * (fh2m/vonkar - fh/vonkar)
   qref   = qm + vonkar/fq*dqh * (fq2m/vonkar - fq/vonkar)
   z0m   = z0mg
   end SUBROUTINE get_fldevp

#endif
end MODULE MOD_CaMa_colmCaMa
