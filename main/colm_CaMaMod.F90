#include <define.h>
module colm_CaMaMod
#if(defined CaMa_Flood)
   use mod_namelist
   USE MOD_CaMa_Variables
   USE PARKIND1,                ONLY: JPRB, JPRM, JPIM
   USE CMF_DRV_CONTROL_MOD,     ONLY: CMF_DRV_INPUT,   CMF_DRV_INIT,    CMF_DRV_END
   USE CMF_DRV_ADVANCE_MOD,     ONLY: CMF_DRV_ADVANCE
   USE CMF_CTRL_FORCING_MOD,    ONLY: CMF_FORCING_GET, CMF_FORCING_PUT
   USE CMF_CTRL_OUTPUT_MOD,     ONLY: CMF_OUTPUT_INIT,CMF_OUTPUT_END
   USE YOS_CMF_INPUT,           ONLY: NXIN, NYIN, DT,DTIN,IFRQ_INP,LLEAPYR,NX,NY,RMIS,DMIS
   use precision,               only: r8,r4
   USE YOS_CMF_INPUT,           ONLY: LROSPLIT,LWEVAP,LWINFILT
   use spmd_task
   use CMF_CTRL_TIME_MOD    
   use GlobalVars, only : spval
   use MOD_1D_Fluxes

   IMPLICIT NONE
   !** local variables
   INTEGER i,j
   INTEGER(KIND=JPIM)              :: ISTEPX              ! total time step
   INTEGER(KIND=JPIM)              :: ISTEPADV            ! time step to be advanced within DRV_ADVANCE
   REAL(KIND=JPRB),ALLOCATABLE     :: ZBUFF(:,:,:)        ! Buffer to store forcing runoff

   real(r8), allocatable :: Effarea     (:,:)  
   real(r8), allocatable :: Effdepth     (:,:)  
   interface colm_CaMa_init
      module procedure colm_CaMa_init
   end interface

   interface colm_CaMa_drv
      module procedure colm_CaMa_drv
   end interface

   interface colm_CaMa_exit
      module procedure colm_CaMa_exit
   end interface
CONTAINS

   subroutine colm_CaMa_init  

      use CMF_CTRL_OUTPUT_MOD
      USE mod_landpatch
      USE MOD_CaMa_Variables
      USE YOS_CMF_TIME,       ONLY: YYYY0
      implicit none

      INTEGER nlon_cama, nlat_cama
      integer dtime
      integer nnn,i,j,IX,IY,mmm
      !*** local variables
      INTEGER(KIND=JPIM)          :: JF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      if(p_is_master)then 
         !*** 1a. Namelist handling
         CALL CMF_DRV_INPUT
         DT       = IFRQ_INP*3600!DEF_simulation_time%timestep
         DTIN     = IFRQ_INP*3600
         SYEAR    = DEF_simulation_time%start_year                              !  start year  
         SMON     = DEF_simulation_time%start_month                            !  month       
         SDAY     = DEF_simulation_time%start_day                             !  day        
         SHOUR    = DEF_simulation_time%start_sec/3600                        !  hour       
         EYEAR    = DEF_simulation_time%end_year                     		!  end year    
         EMON     = DEF_simulation_time%end_month                        	!  month       
         EDAY     = DEF_simulation_time%end_day                        	!  day         
         EHOUR    = DEF_simulation_time%end_sec/3600                        	!  hour 
         LLEAPYR  = DEF_forcing%leapyear
         YYYY0    = SYEAR
         RMIS     = spval
         DMIS     = spval

         !*** 1b. INITIALIZATION
         !CALL CMF_OUTPUT_INIT
         CALL CMF_DRV_INIT
         !*** 2. check variable name & allocate data to pointer DVEC
         DO JF=1,NVARSOUT
            SELECT CASE (VAROUT(JF)%CVNAME)
            CASE ('rivout')
               DEF_hist_cama_vars%rivout=.true.
            CASE ('rivsto')
               DEF_hist_cama_vars%rivsto=.true.
            CASE ('rivdph')
               DEF_hist_cama_vars%rivdph=.true.
            CASE ('rivvel')
               DEF_hist_cama_vars%rivvel=.true.
            CASE ('fldout')
               DEF_hist_cama_vars%fldout=.true.
            CASE ('fldsto')
               DEF_hist_cama_vars%fldsto=.true.
            CASE ('flddph')
               DEF_hist_cama_vars%flddph=.true.
            CASE ('fldfrc')
               DEF_hist_cama_vars%fldfrc=.true.
            CASE ('fldare')
               DEF_hist_cama_vars%fldare=.true.
            CASE ('sfcelv')
               DEF_hist_cama_vars%sfcelv=.true.
            CASE ('totout')
               DEF_hist_cama_vars%totout=.true.
            CASE ('outflw')            !!  compatibility for previous file name
               DEF_hist_cama_vars%outflw=.true.
            CASE ('totsto')
               DEF_hist_cama_vars%totsto=.true.
            CASE ('storge')            !!  compatibility for previous file name
               DEF_hist_cama_vars%storge=.true.
            CASE ('pthout')
               DEF_hist_cama_vars%pthout=.true.
            CASE ('maxflw')
               DEF_hist_cama_vars%maxflw=.true.
            CASE ('maxdph')
               DEF_hist_cama_vars%maxdph=.true.
            CASE ('maxsto')
               DEF_hist_cama_vars%maxsto=.true.
            CASE ('gwsto')
               DEF_hist_cama_vars%gwsto=.true.
            CASE ('gdwsto')
               DEF_hist_cama_vars%gdwsto=.true.
            CASE ('gwout')
               DEF_hist_cama_vars%gwout=.true.
            CASE ('gdwrtn')
               DEF_hist_cama_vars%gdwrtn=.true.
            CASE ('runoff')             !!  compatibility for previous file name
               DEF_hist_cama_vars%runoff=.true. 
            CASE ('runoffsub')           !!  compatibility for previous file name
               DEF_hist_cama_vars%runoffsub=.true. 
            CASE ('rofsfc')
               DEF_hist_cama_vars%rofsfc=.true. 
            CASE ('rofsub')
               DEF_hist_cama_vars%rofsub=.true. 
            CASE ('damsto')   !!! added
               DEF_hist_cama_vars%damsto=.true. 
            CASE ('daminf')   !!! added
               DEF_hist_cama_vars%daminf=.true. 
            CASE ('levsto')   !!! added
               DEF_hist_cama_vars%levsto=.true. 
            CASE ('levdph')   !!! added
               DEF_hist_cama_vars%levdph=.true. 
            CASE ('wevap')   !!! added
               DEF_hist_cama_vars%wevap=.true. 
            CASE ('winfilt')   !!! added
               DEF_hist_cama_vars%winfilt=.true.             
            CASE DEFAULT
               stop
            END SELECT
         end do
      endif
      CALL mpi_bcast (NX      ,   1, MPI_LOGICAL,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (NY      ,   1, MPI_LOGICAL,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (IFRQ_INP ,   1, MPI_LOGICAL,  p_root, p_comm_glb, p_err)
      
      CALL mpi_bcast (LWEVAP ,   1, MPI_LOGICAL,  p_root, p_comm_glb, p_err)
      CALL mpi_bcast (LWINFILT ,   1, MPI_LOGICAL,  p_root, p_comm_glb, p_err)

      CALL gcama%define_by_ndims (NX, NY) 
      call mp2g_cama%build (landpatch, gcama)
      call mg2p_cama%build (gcama, landpatch)

      CALL cama_gather%set (gcama)

      call allocate_2D_cama_Fluxes  (gcama)
      call allocate_acc_cama_Fluxes ()
      call FLUSH_acc_cama_fluxes    ()

      if(p_is_master)then 
         ALLOCATE (runoff_2d (NX,NY))
         ALLOCATE (fevpg_2d  (NX,NY))
         ALLOCATE (finfg_2d  (NX,NY))
         !*** 1c. allocate data buffer for input forcing
         ALLOCATE (ZBUFF(NX,NY,4))
         ALLOCATE (fldfrc_tmp(NX,NY))
         ALLOCATE (flddepth_tmp(NX,NY))
         runoff_2d(:,:)    = 0.0D0
         fevpg_2d(:,:)     = 0.0D0
         finfg_2d(:,:)     = 0.0D0
         ZBUFF(:,:,:)      = 0.0D0
         fldfrc_tmp(:,:)   = 0.0D0
         flddepth_tmp(:,:) = 0.0D0
      endif
      if (p_is_worker) then
         ALLOCATE (flddepth_cama(numpatch))
         ALLOCATE (fldfrc_cama(numpatch))
         ALLOCATE (fevpg_fld(numpatch))
         ALLOCATE (finfg_fld(numpatch))
         flddepth_cama(:)     =  0.0D0
         fldfrc_cama(:)       =  0.0D0
         fevpg_fld(:)         =  0.0D0
         finfg_fld(:)         =  0.0D0
end if  
   end subroutine colm_CaMa_init

   ! -----
   subroutine colm_cama_drv(idate_sec)
      implicit none
      integer, intent(in)  :: idate_sec     ! calendar (year, julian day, seconds)
      call accumulate_cama_fluxes
      if  (MOD(idate_sec,3600*int(IFRQ_INP))==0) then
#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
         call colm2cama_real8 (a_rnof_cama, f_rnof_cama, runoff_2d)
IF (LWEVAP) THEN
#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
         call colm2cama_real8 (a_fevpg_fld, f_fevpg_fld, fevpg_2d)
ENDIF

IF (LWINFILT) THEN
#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
         call colm2cama_real8 (a_finfg_fld, f_finfg_fld, finfg_2d)
ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
         CALL flush_acc_cama_fluxes
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
         if(p_is_master)then
               do i = 1,NX ! cama_gather%ginfo%nlon
                  do j = 1, NY !cama_gather%ginfo%nlat
                     ZBUFF(i,j,1)=runoff_2d(i,j)/1000.0D0 !mm/s -->m/s
                     ZBUFF(i,j,2)=0.0D0
                     IF (LWEVAP) THEN
                        ZBUFF(i,j,3)=fevpg_2d(i,j)/1000.0D0 !mm/s -->m/s
                     ELSE
                        ZBUFF(i,j,3)=0.0D0
                     ENDIF
                     IF (LWINFILT) THEN
                        !if(finfg_2d(i,j)<0.0d0) finfg_2d(i,j)=0.0d0
                        ZBUFF(i,j,4)=finfg_2d(i,j)/1000.0D0  !mm/s -->m/s
                     ELSE                  
                        ZBUFF(i,j,4)=0.0D0
                     ENDIF
                  enddo
               enddo
   
         ! Simulating the hydrodynamics in continental-scale rivers
         ! ----------------------------------------------------------------------
         ISTEPADV=INT(DTIN/DT,JPIM)
         !*  2a Read forcing from file, This is only relevant in Stand-alone mode 
         !CALL CMF_FORCING_GET(ZBUFF(:,:,:))
         !*  2b Interporlate runoff & send to CaMa-Flood 
         
               CALL CMF_FORCING_PUT(ZBUFF)
         !*  2c  Advance CaMa-Flood model for ISTEPADV
         CALL CMF_DRV_ADVANCE(ISTEPADV)
               IF (LWINFILT .or. LWEVAP) THEN
                  call get_fldinfo()
               endif
         endif
         IF (LWINFILT .or. LWEVAP) THEN

#ifdef USEMPI
            CALL mpi_barrier (p_comm_glb, p_err)
#endif   
            call cama2colm_real8 (flddepth_tmp, f_flddepth_cama, flddepth_cama)
#ifdef USEMPI
            CALL mpi_barrier (p_comm_glb, p_err)
#endif         
            call cama2colm_real8 (fldfrc_tmp,   f_fldfrc_cama,   fldfrc_cama  )
#ifdef USEMPI
            CALL mpi_barrier (p_comm_glb, p_err)
#endif               
            flddepth_cama=flddepth_cama*1000.D0 !m --> mm
         endif
   endif
   end subroutine colm_cama_drv

   subroutine colm_cama_exit
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      !*** 3a. finalize CaMa-Flood 
      call deallocate_acc_cama_Fluxes ()
      if(p_is_master)then
         !*** 3a. finalize CaMa-Flood 
         DEALLOCATE(ZBUFF)
         DEALLOCATE (runoff_2d)
         DEALLOCATE (fevpg_2d)
         DEALLOCATE (finfg_2d)

         !  CALL CMF_DRV_END
      endif
      if (p_is_worker) then
         DEALLOCATE (flddepth_cama)
         DEALLOCATE (fldfrc_cama)
         DEALLOCATE (fevpg_fld)
         DEALLOCATE (finfg_fld)
      end if
   end subroutine colm_cama_exit



   !####################################################################
   SUBROUTINE get_fldinfo()
      ! save results to master process
      USE YOS_CMF_INPUT,      ONLY:  NX, NY
      USE YOS_CMF_PROG,       ONLY:  D2FLDSTO 
      USE YOS_CMF_DIAG,       ONLY:  D2FLDDPH,D2FLDFRC
      USE CMF_UTILS_MOD,      ONLY:  VEC2MAPD

      IMPLICIT NONE

      integer i,j

      !================================================
      !! convert 1Dvector to 2Dmap
      CALL VEC2MAPD(D2FLDFRC,flddepth_tmp)             !! MPI node data is gathered by VEC2MAP m
      CALL VEC2MAPD(D2FLDDPH,fldfrc_tmp)             !! MPI node data is gathered by VEC2MAP m
      do i    = 1, NX
         do j = 1, NY
            if (flddepth_tmp(i,j) .lt.    0.0)      flddepth_tmp(i,j)=0.0
            if (fldfrc_tmp(i,j)   .lt.    0.0)      fldfrc_tmp(i,j)=0.0
         enddo
      enddo
   END SUBROUTINE get_fldinfo



   subroutine fldfluxes (hu,ht,hq,&
      us,vs,tm,qm,rhoair,psrf,tssea,&
      taux,tauy,fseng,fevpg,tref,qref,&
      z0m,zol,rib,ustar,qstar,tstar,fm,fh,fq)
      ! compute surface fluxes, derviatives, and exchange coefficiants

      !=======================================================================
      ! this is the main subroutine to execute the calculation of thermal processes
      ! and surface fluxes
      !
      ! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
      !=======================================================================
   
      use precision
      use PhysicalConstants, only : cpair,rgas,vonkar,grav
      use FRICTION_VELOCITY
      implicit none
   
      !----------------------- Dummy argument --------------------------------
      real(r8), INTENT(in) :: hu      ! agcm reference height of wind [m]
      real(r8), INTENT(in) :: ht      ! agcm reference height of temperature [m]
      real(r8), INTENT(in) :: hq      ! agcm reference height of humidity [m]
      real(r8), INTENT(in) :: us      ! wind component in eastward direction [m/s]
      real(r8), INTENT(in) :: vs      ! wind component in northward direction [m/s]
      real(r8), INTENT(in) :: tm      ! temperature at agcm reference height [kelvin]
      real(r8), INTENT(in) :: qm      ! specific humidity at agcm reference height [kg/kg]
      real(r8), INTENT(in) :: rhoair  ! density air [kg/m3]
      real(r8), INTENT(in) :: psrf    ! atmosphere pressure at the surface [pa] [not used]
      real(r8), INTENT(in) :: tssea   ! inundation surface temperature [K]-->set to tgrnd
   
      real(r8), INTENT(out) :: &
      taux,     &! wind stress: E-W [kg/m/s**2]
      tauy,     &! wind stress: N-S [kg/m/s**2]
      fseng,    &! sensible heat flux from ground [mm/s]
      fevpg,    &! evaporation heat flux from ground [mm/s]
   
      tref,     &! 2 m height air temperature [kelvin]
      qref,     &! 2 m height air humidity
      z0m,      &! effective roughness [m]
      zol,      &! dimensionless height (z/L) used in Monin-Obukhov theory
      rib,      &! bulk Richardson number in surface layer
      ustar,    &! friction velocity [m/s]
      tstar,    &! temperature scaling parameter
      qstar,    &! moisture scaling parameter
      fm,       &! integral of profile function for momentum
      fh,       &! integral of profile function for heat
      fq!,       &! integral of profile function for moisture
!      cgrndl,   &! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
!      cgrnds     ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
   
   !------------------------ LOCAL VARIABLES ------------------------------
   integer i
   integer niters, &! maximum number of iterations for surface temperature
   iter,      &! iteration index
   nmozsgn     ! number of times moz changes sign
   
   real(r8) :: &
   beta,      &! coefficient of conective velocity [-]
   displax,   &! zero-displacement height [m]
   dth,       &! diff of virtual temp. between ref. height and surface
   dqh,       &! diff of humidity between ref. height and surface
   dthv,      &! diff of vir. poten. temp. between ref. height and surface
   eg,        &! water vapor pressure at temperature T [Pa]
   degdT,     &! d(eg)/dT
   obu,       &! monin-obukhov length (m)
   obuold,    &! monin-obukhov length from previous iteration
   qsatg,     &! ground saturated specific humidity [kg/kg]
   qsatgdT,   &! d(qsatg)/dT
   ram,       &! aerodynamical resistance [s/m]
   rah,       &! thermal resistance [s/m]
   raw,       &! moisture resistance [s/m]
   raih,      &! temporary variable [kg/m2/s]
   raiw,      &! temporary variable [kg/m2/s]
   fh2m,      &! relation for temperature at 2m
   fq2m,      &! relation for specific humidity at 2m
   fm10m,     &! integral of profile function for momentum at 10m
   thm,       &! intermediate variable (tm+0.0098*ht)
   th,        &! potential temperature (kelvin)
   thv,       &! virtual potential temperature (kelvin)
   thvstar,   &! virtual potential temperature scaling parameter
   um,        &! wind speed including the stablity effect [m/s]
   ur,        &! wind speed at reference height [m/s]
   visa,      &! kinematic viscosity of dry air [m2/s]
   wc,        &! convective velocity [m/s]
   wc2,       &! wc**2
   xt,        &!
   xq,        &!
   zii,       &! convective boundary height [m]
   zldis,     &! reference height "minus" zero displacement heght [m]
   z0mg,      &! roughness length over ground, momentum [m]
   z0hg,      &! roughness length over ground, sensible heat [m]
   z0qg        ! roughness length over ground, latent heat [m]
   
   real, parameter :: zsice = 0.04  ! sea ice aerodynamic roughness length [m]
   
   !-----------------------------------------------------------------------
   ! potential temperatur at the reference height
   beta = 1.      ! -  (in computing W_*)
   zii = 1000.    ! m  (pbl height)
   
   !-----------------------------------------------------------------------
   !     Compute sensible and latent fluxes and their derivatives with respect 
   !     to ground temperature using ground temperatures from previous time step.
   !-----------------------------------------------------------------------
   ! Initialization variables
   nmozsgn = 0
   obuold = 0.
   
   call qsadv(tssea,psrf,eg,degdT,qsatg,qsatgdT)
   
   ! potential temperatur at the reference height
   thm = tm + 0.0098*ht              ! intermediate variable equivalent to
                       ! tm*(pgcm/psrf)**(rgas/cpair)
   th = tm*(100000./psrf)**(rgas/cpair) ! potential T
   thv = th*(1.+0.61*qm)             ! virtual potential T
   ur = max(0.1,sqrt(us*us+vs*vs))   ! limit set to 0.1
   
   dth   = thm-tssea
   dqh   = qm-qsatg
   dthv  = dth*(1.+0.61*qm)+0.61*th*dqh
   zldis = hu-0.
   
   ! Kinematic viscosity of dry air (m2/s)- Andreas (1989) CRREL Rep. 89-11
   visa=1.326e-5*(1.+6.542e-3*tm + 8.301e-6*tm**2 - 4.84e-9*tm**3)
   
   ! loop to obtain initial and good ustar and zo
   ustar=0.06
   wc=0.5
   if(dthv.ge.0.) then
   um=max(ur,0.1)
   else
   um=sqrt(ur*ur+wc*wc)
   endif
   
   do i=1,5
   z0mg=0.013*ustar*ustar/grav+0.11*visa/ustar
   ustar=vonkar*um/log(zldis/z0mg)
   enddo

   
   call moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu)
   
   ! Evaluated stability-dependent variables using moz from prior iteration
   niters=10
   displax = 0.
   
   !----------------------------------------------------------------
   ITERATION : do iter = 1, niters         ! begin stability iteration
   !----------------------------------------------------------------
   
   !if(nint(oro).eq.0)then   ! ocean
   z0mg=0.013*ustar*ustar/grav + 0.11*visa/ustar
   xq=2.67*(ustar*z0mg/visa)**0.25 - 2.57
   xt= xq
   z0qg=z0mg/exp(xq)
   z0hg=z0mg/exp(xt)
  ! endif
   
   call moninobuk(hu,ht,hq,displax,z0mg,z0hg,z0qg,obu,um,&
       ustar,fh2m,fq2m,fm10m,fm,fh,fq)
   
   tstar = vonkar/fh*dth
   qstar = vonkar/fq*dqh
   
   thvstar=tstar+0.61*th*qstar
   zol=zldis*vonkar*grav*thvstar/(ustar**2*thv)
   if(zol >= 0.) then       ! stable
   zol = min(2.,max(zol,1.e-6))
   else                     ! unstable
   zol = max(-100.,min(zol,-1.e-6))
   endif
   obu = zldis/zol
   
   if(zol >= 0.)then
   um = max(ur,0.1)
   else
   wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
   wc2 = beta*beta*(wc*wc)
   um = sqrt(ur*ur+wc2)
   endif
   
   if (obuold*obu < 0.) nmozsgn = nmozsgn+1
   if(nmozsgn >= 4) EXIT
   
   obuold = obu
   
   !----------------------------------------------------------------
   enddo ITERATION                         ! end stability iteration
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
   
   ! surface fluxes of momentum, sensible and latent 
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
   end subroutine fldfluxes

#endif
end module colm_camaMod
