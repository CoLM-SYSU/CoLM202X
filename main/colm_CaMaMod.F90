#include <define.h>
module colm_CaMaMod
#if(defined CaMa_Flood)
   use mod_namelist
   USE MOD_CaMa_Variables
   USE PARKIND1,                ONLY: JPRB, JPRM, JPIM
   !USE YOS_CMF_TIME,            ONLY: NSTEPS
   USE CMF_DRV_CONTROL_MOD,     ONLY: CMF_DRV_INPUT,   CMF_DRV_INIT,    CMF_DRV_END
   USE CMF_DRV_ADVANCE_MOD,     ONLY: CMF_DRV_ADVANCE
   USE CMF_CTRL_FORCING_MOD,    ONLY: CMF_FORCING_GET, CMF_FORCING_PUT
   USE CMF_CTRL_OUTPUT_MOD,     ONLY: CMF_OUTPUT_INIT,CMF_OUTPUT_END
   !USE CMF_CTRL_RESTART_MOD,    ONLY: CMF_RESTART_WRITE
   USE YOS_CMF_INPUT,           ONLY: NXIN, NYIN, DT,DTIN,IFRQ_INP,LLEAPYR,NX,NY,RMIS,DMIS

   !use YOS_CMF_MAP,             only: I2NEXTX
   use precision,               only: r8,r4
   !use YOS_CMF_MAP,             ONLY: D1LON, D1LAT ,D2GRAREA
   !USE CMF_UTILS_MOD,           ONLY: VEC2MAPD
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

   subroutine colm_CaMa_init !(nlon_cama, nlat_cama)

      use CMF_CTRL_OUTPUT_MOD
      USE mod_landpatch
      USE MOD_CaMa_Variables
      USE YOS_CMF_TIME,       ONLY: YYYY0
      implicit none

      !INTEGER, intent(in) :: nlon_cama, nlat_cama
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
         DT       = DEF_simulation_time%timestep
         DTIN     = DEF_simulation_time%timestep
         IFRQ_INP = DTIN/3600.0
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
            CASE DEFAULT
               stop
            END SELECT
         end do
      endif
      CALL mpi_bcast (NX      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (NY      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)

      nlon_cama=NX
      nlat_cama=NY
      CALL gcama%define_by_ndims (nlon_cama, nlat_cama) 
      call mp2g_cama%build (landpatch, gcama)
      call mg2p_cama%build (gcama, landpatch)

      CALL cama_gather%set (gcama)

      call allocate_2D_cama_Fluxes  (gcama)
      call allocate_acc_cama_Fluxes ()
      call FLUSH_acc_cama_fluxes    ()

      if(p_is_master)then 
         allocate (runoff_2d (cama_gather%ginfo%nlon,cama_gather%ginfo%nlat))
         !*** 1c. allocate data buffer for input forcing
         ALLOCATE (ZBUFF(NXIN,NYIN,2))
         ALLOCATE (Effarea(NX,NY))
         ALLOCATE (Effdepth(NX,NY))
      endif

#ifdef USEMPI  
      CALL mpi_bcast (DEF_hist_cama_vars%rivout      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%rivsto      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%rivdph      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%rivvel      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%fldout      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%fldsto      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%flddph      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%fldfrc      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%fldare      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%sfcelv      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%totout      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%outflw      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%totsto      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%storge      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%pthout      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%maxflw      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%maxdph      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%maxsto      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%gwsto       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%gdwsto      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%gwout       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%gdwrtn      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%runoff      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%runoffsub   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%rofsfc      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%rofsub      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%damsto      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_cama_vars%daminf      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
#endif  
   end subroutine colm_CaMa_init

   ! -----
   subroutine colm_cama_drv
      implicit none
    !  call colm2cama_real8 (rnof, f_rnof_cama, runoff_2d)
      CALL flush_acc_cama_fluxes 
      call accumulate_cama_fluxes
      call colm2cama_real8 (a_rnof_cama, f_rnof_cama, runoff_2d)
      !CALL flush_acc_cama_fluxes !may not need
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      if(p_is_master)then
         do j = 1, cama_gather%ginfo%nlat
            do i = 1, cama_gather%ginfo%nlon
               if(runoff_2d(i,j) < 1.e-10) runoff_2d(i,j) = 0.
               ZBUFF(i,j,1)=runoff_2d(i,j)/1000.0
               ZBUFF(i,j,2)=0.D0
            enddo
         enddo
   
         ! Simulating the hydrodynamics in continental-scale rivers
         ! ----------------------------------------------------------------------
         ISTEPADV=INT(DTIN/DT,JPIM)
         !*  2a Read forcing from file, This is only relevant in Stand-alone mode 
         !CALL CMF_FORCING_GET(ZBUFF(:,:,:))
         !*  2b Interporlate runoff & send to CaMa-Flood 
         
         CALL CMF_FORCING_PUT(ZBUFF(:,:,:))
         !*  2c  Advance CaMa-Flood model for ISTEPADV
         CALL CMF_DRV_ADVANCE(ISTEPADV)
         Effarea(:,:)=0.0
         Effdepth(:,:)=0.0
         call get_flddepth()
      endif
      !     call master2IO_2d_real8(Effdepth,IO_Effdepth)
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
         deallocate (runoff_2d)
         deallocate (Effarea)
         deallocate (Effdepth)
         !  CALL CMF_DRV_END
      endif
   end subroutine colm_cama_exit



   !####################################################################
   SUBROUTINE get_flddepth
      USE CMF_UTILS_MOD,           ONLY: VEC2MAP
      ! save results to master process
      USE YOS_CMF_INPUT,      ONLY: NX, NY
      USE YOS_CMF_PROG,       ONLY:   D2FLDSTO    
      USE YOS_CMF_MAP,        ONLY:  D2GRAREA,D2RIVWTH,D2RIVLEN


      !!!!-----------will be used when master has mutiple nodes
      !#ifdef UseMPI
      !   USE CMF_CTRL_MPI_MOD,   ONLY: CMF_MPI_REDUCE_R2MAP, CMF_MPI_REDUCE_R1PTH
      !#endif
      !!!!-----------will be used when master has mutiple nodes

      IMPLICIT NONE
      !*** LOCAL
      REAL(KIND=JPRM)             :: D2RIVWTHVEC(NX,NY)
      REAL(KIND=JPRM)             :: D2RIVLENVEC(NX,NY)
      REAL(KIND=JPRM)             :: D2FLDSTOVEC(NX,NY)
      REAL(KIND=JPRM)             :: D2GRAREAVEC(NX,NY)
      integer i,j

      !================================================
      !! convert 1Dvector to 2Dmap
      CALL VEC2MAP(D2RIVWTH,D2RIVWTHVEC)             !! MPI node data is gathered by VEC2MAP
      CALL VEC2MAP(D2RIVLEN,D2RIVLENVEC)             !! MPI node data is gathered by VEC2MAP
      CALL VEC2MAP(D2FLDSTO,D2FLDSTOVEC)             !! MPI node data is gathered by VEC2MAP
      CALL VEC2MAP(D2GRAREA,D2GRAREAVEC)             !! MPI node data is gathered by VEC2MAP
      do j = 1,NY
         do i = 1,NX
            Effarea(i,j)=max(D2GRAREAVEC(i,j)-D2RIVWTHVEC(i,j)*D2RIVLENVEC(i,j),0.0d0)
            if (Effarea(i,j)>0.0 .and. D2FLDSTOVEC(i,j)>1.e-5) then
               Effdepth(i,j)=max(D2FLDSTOVEC(i,j)/Effarea(i,j),0.0d0)
            endif
         enddo
      enddo


      !!!!-----------will be used when master has mutiple nodes
      !#ifdef UseMPI
      !      CALL CMF_MPI_REDUCE_R2MAP(R2OUT)
      !#endif
   END SUBROUTINE get_flddepth

#endif
end module colm_camaMod
