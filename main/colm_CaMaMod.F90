#include <define.h>
module colm_CaMaMod
#if(defined CaMa_Flood)
   use mod_namelist
   use MOD_1D_cama_Fluxes
   use MOD_2D_cama_Fluxes
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
   use mod_hist,                only:ghist,var_out   !,mp2g_hist     
   use CMF_CTRL_TIME_MOD    
   use GlobalVars, only : spval
   IMPLICIT NONE
   !** local variables
   INTEGER i,j
   INTEGER(KIND=JPIM)              :: ISTEPX              ! total time step
   INTEGER(KIND=JPIM)              :: ISTEPADV            ! time step to be advanced within DRV_ADVANCE
   REAL(KIND=JPRB),ALLOCATABLE     :: ZBUFF(:,:,:)        ! Buffer to store forcing runoff
   real(r8), allocatable           :: runoff_2d(:,:)         ! total runoff [mm/s]

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
      implicit none
      integer dtime
      integer nnn,i,j,IX,IY,mmm
      !*** local variables
      INTEGER(KIND=JPIM)          :: JF
      call allocate_2D_cama_Fluxes(ghist)
      call allocate_1D_cama_Fluxes ()
      if(p_is_master)then 
         allocate (runoff_2d(DEF_nlon_hist,DEF_nlat_hist))
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
         print *,EHOUR       
         LLEAPYR  = DEF_forcing%leapyear
         RMIS     = spval
         DMIS     = spval
         print *,DTIN,DT,IFRQ_INP

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

   subroutine colm_cama_drv
      implicit none
      call var_out (runoff_2d)
      if(p_is_master)then
         do j = 1, DEF_nlat_hist
            do i = 1, DEF_nlon_hist
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
      !*** 3a. finalize CaMa-Flood 
      call deallocate_1D_cama_Fluxes ()
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

SUBROUTINE master2IO_2d_real8 (InVar,OutVar)
   !=======================================================================
   ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
   !=======================================================================     
   use precision
   use mod_namelist
   use timemanager
   use spmd_task
   use mod_2d_cama_fluxes
   use MOD_1D_cama_Fluxes
   use mod_block
   use mod_data_type
   use mod_landpatch
   use mod_mapping_pset2grid
   use mod_colm_debug
   USE MOD_TimeInvariants, only : patchtype
   USE MOD_1D_Acc_cama_Fluxes
   USE mod_hist
   use mod_grid
   !use GlobalVars, only : spval
   IMPLICIT NONE
   real(r8), INTENT(in) ::  Invar (DEF_nlon_hist, DEF_nlat_hist)
   type(block_data_real8_2d), INTENT(inout) :: OutVar
   real(r8), allocatable     ::  vectmp(:)  
   logical,  allocatable     ::  filter(:)
   integer :: xblk, yblk, xloc, yloc
   integer :: iblk, jblk, idata, ixseg, iyseg
   integer :: rmesg(2), smesg(2), isrc, iproc
   real(r8), allocatable :: rbuf(:,:), sbuf(:,:), vdata(:,:)
   integer :: xdsp, ydsp, xcnt, ycnt
   character(len=256) :: fileblock
   !call ghist%define_by_ndims (NXIN, NYIN)
   !call mp2g_hist%build (landpatch, ghist)
   !call set_segment_info (ghist)
   if (p_is_master) then
      do iyseg = 1, hist_block_info%nyseg
         do ixseg = 1, hist_block_info%nxseg
             iblk = hist_block_info%xsegs(ixseg)%blk
             jblk = hist_block_info%ysegs(iyseg)%blk
             xdsp = hist_block_info%xsegs(ixseg)%gdsp
             ydsp = hist_block_info%ysegs(iyseg)%gdsp
             xcnt = hist_block_info%xsegs(ixseg)%cnt
             ycnt = hist_block_info%ysegs(iyseg)%cnt
             
             allocate (sbuf (xcnt,ycnt))
             sbuf = InVar (xdsp+1:xdsp+xcnt, ydsp+1:ydsp+ycnt)
             smesg = (/ixseg, iyseg/)
             call mpi_send (smesg, 2, MPI_INTEGER, &
                gblock%pio(iblk,jblk), 10000, p_comm_glb, p_err) 
             call mpi_send (sbuf, xcnt*ycnt, MPI_DOUBLE, &
                gblock%pio(iblk,jblk), 10000, p_comm_glb, p_err)
             deallocate (sbuf)
          end do
       end do

       DO iproc = 0, p_np_io-1
          smesg = (/0, 0/)
          CALL mpi_send(smesg, 2, MPI_INTEGER, p_address_io(iproc), 10000, p_comm_glb, p_err)
       ENDDO
   elseif  (p_is_io) then
      DO WHILE (.true.)
         call mpi_recv (rmesg, 2, MPI_INTEGER, p_root, 10000, p_comm_glb, p_stat, p_err)
         ixseg = rmesg(1)
         iyseg = rmesg(2)
             
         IF ((ixseg > 0) .and. (iyseg > 0)) THEN
            iblk = hist_block_info%xsegs(ixseg)%blk
            jblk = hist_block_info%ysegs(iyseg)%blk
            xdsp = hist_block_info%xsegs(ixseg)%bdsp
            ydsp = hist_block_info%ysegs(iyseg)%bdsp
            xcnt = hist_block_info%xsegs(ixseg)%cnt
            ycnt = hist_block_info%ysegs(iyseg)%cnt

            allocate (rbuf(xcnt,ycnt))
            call mpi_recv (rbuf, xcnt*ycnt, MPI_DOUBLE, &
               p_root, 10000, p_comm_glb, p_stat, p_err)
            OutVar%blk(iblk,jblk)%val(xdsp+1:xdsp+xcnt,ydsp+1:ydsp+ycnt)= rbuf
            deallocate (rbuf)
         ELSE
            exit
         ENDIF
     end do
   endif
   END SUBROUTINE master2IO_2d_real8

   

   

         
#endif
end module colm_camaMod
