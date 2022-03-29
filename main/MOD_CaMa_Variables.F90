#include <define.h>

module MOD_CaMa_Variables
#if(defined CaMa_Flood)

   use precision
   USE mod_grid
   use mod_data_type
   USE mod_mapping_pset2grid
   USE mod_mapping_grid2pset

   real(r8) :: nacc              ! number of accumulation
   real(r8), allocatable     :: a_rnof_cama (:) ! on worker : total runoff [mm/s]
   type(block_data_real8_2d) :: f_rnof_cama     ! on IO     : total runoff [mm/s]
   real(r8), allocatable     :: runoff_2d (:,:) ! on Master : total runoff [mm/s]

   TYPE(grid_type) :: gcama
   
   TYPE (mapping_pset2grid_type) :: mp2g_cama
   TYPE (mapping_grid2pset_type) :: mg2p_cama

   TYPE (grid_concat_type) :: cama_gather
 
   type(block_data_real8_2d) :: IO_Effdepth    ! inundation to water depth [m]
   type(block_data_real8_2d) :: IO_Effarea

   TYPE history_var_cama_type
      LOGICAL :: rivout       = .false.
      LOGICAL :: rivsto       = .false.
      LOGICAL :: rivdph       = .false.
      LOGICAL :: rivvel       = .false.
      LOGICAL :: fldout       = .false.
      LOGICAL :: fldsto       = .false.
      LOGICAL :: flddph       = .false.
      LOGICAL :: fldfrc       = .false.
      LOGICAL :: fldare       = .false.
      LOGICAL :: sfcelv       = .false.
      LOGICAL :: totout       = .false.
      LOGICAL :: outflw       = .false.
      LOGICAL :: totsto       = .false.
      LOGICAL :: storge       = .false.
      LOGICAL :: pthout       = .false.
      LOGICAL :: maxflw       = .false.
      LOGICAL :: maxdph       = .false.
      LOGICAL :: maxsto       = .false.
      LOGICAL :: gwsto        = .false.
      LOGICAL :: gdwsto       = .false.
      LOGICAL :: gwout        = .false.
      LOGICAL :: gdwrtn       = .false.
      LOGICAL :: runoff       = .false.
      LOGICAL :: runoffsub    = .false.
      LOGICAL :: rofsfc       = .false.
      LOGICAL :: rofsub       = .false.
      LOGICAL :: damsto       = .false.
      LOGICAL :: daminf       = .false.
   END TYPE history_var_cama_type

   TYPE (history_var_cama_type) :: DEF_hist_cama_vars

   ! --- subroutines ---
   public :: allocate_acc_cama_fluxes
   public :: deallocate_acc_cama_fluxes
   public :: flush_acc_cama_fluxes
   public :: accumulate_cama_fluxes
   
   public :: allocate_2D_cama_Fluxes
   
   PUBLIC :: colm2cama_real8
   PUBLIC :: cama2colm_real8

   PUBLIC :: hist_out_cama

contains

   ! -----
   subroutine allocate_acc_cama_fluxes 

      use spmd_task
      USE GlobalVars
      use mod_landpatch, only : numpatch
      implicit none

      if (p_is_worker) then
         if (numpatch > 0) then
            allocate (a_rnof_cama(numpatch))
         end if
      end if

   end subroutine allocate_acc_cama_fluxes

   ! -----
   subroutine deallocate_acc_cama_fluxes()

      use spmd_task
      use mod_landpatch, only : numpatch
      implicit none

      if (p_is_worker) then
         if (numpatch > 0) then
            deallocate (a_rnof_cama)
         end if
      end if

   end subroutine deallocate_acc_cama_fluxes

   ! -----
   SUBROUTINE FLUSH_acc_cama_fluxes ()

      use spmd_task
      use mod_landpatch, only : numpatch
      use GlobalVars,    only : spval 
      implicit none

      if (p_is_worker) then

         nacc = 0

         if (numpatch > 0) then
            ! flush the Fluxes for accumulation
            a_rnof_cama (:) = spval
         end if
      end if

   END SUBROUTINE FLUSH_acc_cama_fluxes

   ! -----
   SUBROUTINE accumulate_cama_fluxes 

      use precision
      use spmd_task
      USE MOD_1D_Fluxes, only : rnof
      use mod_landpatch, only : numpatch

      IMPLICIT NONE

      if (p_is_worker) then
         if (numpatch > 0) then
            nacc = nacc + 1
            call acc1d_cama (rnof, a_rnof_cama)  
         end if
      end if

   END SUBROUTINE accumulate_cama_fluxes

   ! -----
   SUBROUTINE acc1d_cama (var, s)

      use precision
      use GlobalVars, only: spval

      IMPLICIT NONE

      real(r8), intent(in)    :: var(:)
      real(r8), intent(inout) :: s  (:)
      ! Local variables
      integer :: i

      do i = lbound(var,1), ubound(var,1)
         if (var(i) /= spval) then
            if (s(i) /= spval) then
               s(i) = s(i) + var(i)
            else
               s(i) = var(i)
            end if
         end if
      end do
      
   END SUBROUTINE acc1d_cama

   ! -----
   SUBROUTINE allocate_2D_cama_Fluxes (grid)
      ! --------------------------------------------------------------------
      ! Allocates memory for CLM 2d [lon_points,lat_points] variables
      ! --------------------------------------------------------------------

      use spmd_task
      use mod_grid
      use mod_data_type
      implicit none

      type(grid_type), intent(in) :: grid

      if (p_is_io) then
         call allocate_block_data (grid, f_rnof_cama)  ! total runoff [mm/s]
         call allocate_block_data (grid, IO_Effdepth)  ! inundation to water depth [m]
         call allocate_block_data (grid, IO_Effarea)  ! inundation to water depth [m]
      end if

   END SUBROUTINE allocate_2D_cama_Fluxes

   ! -----
   SUBROUTINE hist_out_cama (file_hist, itime_in_file)

      USE spmd_task
      USE CMF_CALC_DIAG_MOD,  ONLY: CMF_DIAG_AVERAGE, CMF_DIAG_RESET
      USE YOS_CMF_PROG,       ONLY: D2RIVSTO,     D2FLDSTO,     D2GDWSTO, &
         d2damsto !!! added
      USE YOS_CMF_DIAG,       ONLY: D2RIVDPH,     D2FLDDPH,     D2FLDFRC,     D2FLDARE,     &
         D2SFCELV,     D2STORGE, &
         D2OUTFLW_AVG, D2RIVOUT_AVG, D2FLDOUT_AVG, D2PTHOUT_AVG, D1PTHFLW_AVG, &
         D2RIVVEL_AVG, D2GDWRTN_AVG, D2RUNOFF_AVG, D2ROFSUB_AVG,               &
         D2OUTFLW_MAX, D2STORGE_MAX, D2RIVDPH_MAX, &
         d2daminf_avg   !!! added

      IMPLICIT NONE
      
      character(LEN=*), intent(in) :: file_hist
      integer, intent(in) :: itime_in_file

      if (p_is_master) then
         !*** average variable
         CALL CMF_DIAG_AVERAGE

         !*** write output data
         call flux_map_and_write_2d_cama (DEF_hist_cama_vars%rivout, &
            D2RIVOUT_AVG, file_hist, 'rivout', itime_in_file,'river discharge','m3/s')

         call flux_map_and_write_2d_cama (DEF_hist_cama_vars%rivsto, &
            D2RIVSTO, file_hist, 'rivsto', itime_in_file,'river storage','m3')

         call flux_map_and_write_2d_cama (DEF_hist_cama_vars%rivdph, &
            D2RIVDPH, file_hist, 'rivdph', itime_in_file,'river depth','m')

         call flux_map_and_write_2d_cama (DEF_hist_cama_vars%rivvel, &
            D2RIVVEL_AVG, file_hist, 'rivvel', itime_in_file,'river velocity','m/s')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%fldout, &
            D2FLDOUT_AVG, file_hist, 'fldout', itime_in_file,'floodplain discharge','m3/s')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%fldsto, &
            D2FLDSTO, file_hist, 'fldsto', itime_in_file,'floodplain storage','m3')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%flddph, &
            D2FLDDPH, file_hist, 'flddph', itime_in_file,'floodplain depth','m')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%fldfrc, & 
            D2FLDFRC, file_hist, 'fldfrc', itime_in_file,'flooded fraction','0-1')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%fldare, &
            D2FLDARE, file_hist, 'fldare', itime_in_file,'flooded area','m2')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%sfcelv, &
            D2SFCELV, file_hist, 'sfcelv', itime_in_file,'water surface elevation','m')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%totout, &
            D2OUTFLW_AVG, file_hist, 'totout', itime_in_file,'discharge (river+floodplain)','m3/s')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%outflw, &
            D2OUTFLW_AVG, file_hist, 'outflw', itime_in_file,'discharge (river+floodplain)','m3/s')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%totsto, &
            D2STORGE, file_hist, 'totsto', itime_in_file,'total storage (river+floodplain)','m3')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%storge, &
            D2STORGE, file_hist, 'storge', itime_in_file,'total storage (river+floodplain)','m3')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%pthout, &
            D2PTHOUT_AVG, file_hist, 'pthout', itime_in_file,'net bifurcation discharge','m3/s')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%gdwsto, &
            D2GDWSTO, file_hist, 'gdwsto', itime_in_file,'ground water storage','m3')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%gwsto, &
            D2GDWSTO, file_hist, 'gwsto', itime_in_file,'ground water storage','m3')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%gwout, &
            D2GDWRTN_AVG, file_hist, 'gwout', itime_in_file,'ground water discharge','m3/s')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%runoff, &
            D2RUNOFF_AVG, file_hist, 'runoff', itime_in_file,'Surface runoff','m3/s')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%runoffsub, &
            D2ROFSUB_AVG  , file_hist, 'runoffsub', itime_in_file,'sub-surface runoff','m3/s')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%maxsto, &
            D2STORGE_MAX, file_hist, 'maxsto', itime_in_file,'daily maximum storage','m3')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%maxflw, &
            D2OUTFLW_MAX, file_hist, 'maxflw', itime_in_file,'daily maximum discharge','m3/s')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%maxdph, &
            D2RIVDPH_MAX, file_hist, 'maxdph', itime_in_file,'daily maximum river depth','m')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%damsto, &
            d2damsto, file_hist, 'damsto', itime_in_file,'reservoir storage','m3')

         call flux_map_and_write_2d_cama(DEF_hist_cama_vars%daminf, &
            d2daminf_avg, file_hist, 'daminf', itime_in_file,'reservoir inflow','m3/s')

         !*** reset variable
         CALL CMF_DIAG_RESET
      ENDIF

   END SUBROUTINE hist_out_cama 

   ! -----
   SUBROUTINE hist_write_cama_time (filename, dataname, time, itime)

      use spmd_task
      use ncio_serial
      USE YOS_CMF_INPUT, ONLY: NX, NY
      USE YOS_CMF_MAP,   ONLY: D1LON, D1LAT
      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      integer, intent(in)  :: time(3)
      integer, intent(out) :: itime

      ! Local variables
      logical :: fexists

      if (p_is_master) then
         inquire (file=filename, exist=fexists)
         if (.not. fexists) then
            call ncio_create_file (trim(filename))
            CALL ncio_define_dimension(filename, 'time', 0)
            call ncio_define_dimension(filename,'lat_cama', NY)
            call ncio_define_dimension(filename,'lon_cama', NX)
            call ncio_write_serial (filename, 'lat_cama', D1LAT,'lat_cama')
            call ncio_write_serial (filename, 'lon_cama', D1LON,'lon_cama')
         endif

         call ncio_write_time (filename, dataname, time, itime)
      ENDIF

   END SUBROUTINE hist_write_cama_time

   ! -----
   subroutine flux_map_and_write_2d_cama (is_hist, &
         var_in, file_hist, varname, itime_in_file,longname,units)

      USE mod_namelist
      USE YOS_CMF_INPUT,  ONLY: NX, NY
      USE YOS_CMF_MAP,    ONLY: NSEQALL
      USE PARKIND1,       ONLY: JPRM
      USE CMF_UTILS_MOD,  ONLY: VEC2MAPD
      use ncio_serial

      IMPLICIT NONE
      logical, intent(in) :: is_hist
      real(r8), INTENT(in) ::  var_in (NSEQALL, 1) 
      character(len=*), intent(in) :: file_hist
      character(len=*), intent(in) :: varname
      integer, intent(in) :: itime_in_file
      character (len=*), intent(in),optional :: longname
      character (len=*), intent(in),optional :: units

      REAL(r8) :: R2OUT(NX,NY)
      integer  :: compress

      if (.not. is_hist) return

      CALL VEC2MAPD(var_in,R2OUT)

      compress = DEF_HIST_COMPRESS_LEVEL 
      call ncio_write_serial_time (file_hist, varname,  &
         itime_in_file, R2OUT, 'lon_cama', 'lat_cama', 'time',compress,longname,units)

   end subroutine flux_map_and_write_2d_cama

   ! -----
   SUBROUTINE colm2cama_real8 (WorkerVar, IOVar, MasterVar)

      !=======================================================================
      ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
      !=======================================================================

      use precision
      use mod_namelist
      use timemanager
      use spmd_task
      use mod_block
      use mod_data_type
      use mod_landpatch
      use mod_mapping_pset2grid
      use mod_colm_debug
      USE MOD_TimeInvariants, only : patchtype

      !use GlobalVars, only : spval
      IMPLICIT NONE

      real(r8),                  intent(inout) :: WorkerVar(:)
      TYPE(block_data_real8_2d), intent(inout) :: IOVar
      real(r8),                  INTENT(inout) :: MasterVar(:,:)

      type(block_data_real8_2d) :: sumwt
      real(r8), allocatable     :: vectmp(:)  
      logical,  allocatable     :: filter(:)
      integer :: xblk, yblk, xloc, yloc
      integer :: iblk, jblk, idata, ixseg, iyseg
      integer :: rmesg(3), smesg(3), isrc
      real(r8), allocatable :: rbuf(:,:), sbuf(:,:), vdata(:,:)
      integer :: xdsp, ydsp, xcnt, ycnt

      if(p_is_master)then 
         MasterVar(:,:) = spval
      endif 

      IF (p_is_worker) THEN
         where (WorkerVar /= spval) 
            WorkerVar = WorkerVar / nacc
         endwhere

         if (numpatch > 0) then
            allocate (filter (numpatch))
            allocate (vectmp (numpatch))

            filter(:) = patchtype < 99
            vectmp (:) = 1.
         end if
      ENDIF

      CALL mp2g_cama%map (WorkerVar, IOVar, spv = spval, msk = filter) 

      if (p_is_io) then
         call allocate_block_data (gcama, sumwt)
      end if

      call mp2g_cama%map (vectmp, sumwt, spv = spval, msk = filter)

      if (p_is_io) then
         do yblk = 1, gblock%nyblk
            do xblk = 1, gblock%nxblk
               if (gblock%pio(xblk,yblk) == p_iam_glb) then
                  do yloc = 1, gcama%ycnt(yblk) 
                     do xloc = 1, gcama%xcnt(xblk) 

                        if (sumwt%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) then
                           IF (IOVar%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                              IOVar%blk(xblk,yblk)%val(xloc,yloc) &
                                 = IOVar%blk(xblk,yblk)%val(xloc,yloc) &
                                 / sumwt%blk(xblk,yblk)%val(xloc,yloc)
                           ENDIF
                        else
                           IOVar%blk(xblk,yblk)%val(xloc,yloc) = spval
                        end if

                     end do
                  end do

               end if
            end do
         end do
      end if

      if (p_is_master) then
         do idata = 1, cama_gather%ndatablk
            call mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, 10011, p_comm_glb, p_stat, p_err)
            isrc  = rmesg(1)
            ixseg = rmesg(2)
            iyseg = rmesg(3)

            xdsp = cama_gather%xsegs(ixseg)%gdsp
            ydsp = cama_gather%ysegs(iyseg)%gdsp
            xcnt = cama_gather%xsegs(ixseg)%cnt
            ycnt = cama_gather%ysegs(iyseg)%cnt

            allocate (rbuf(xcnt,ycnt))
            call mpi_recv (rbuf, xcnt * ycnt, MPI_DOUBLE, &
               isrc, 10011, p_comm_glb, p_stat, p_err)
            MasterVar (xdsp+1:xdsp+xcnt,ydsp+1:ydsp+ycnt) = rbuf
            deallocate (rbuf)
         end do

      elseif (p_is_io) then
         do iyseg = 1, cama_gather%nyseg
            do ixseg = 1, cama_gather%nxseg

               iblk = cama_gather%xsegs(ixseg)%blk
               jblk = cama_gather%ysegs(iyseg)%blk

               if (gblock%pio(iblk,jblk) == p_iam_glb) then
                  xdsp = cama_gather%xsegs(ixseg)%bdsp
                  ydsp = cama_gather%ysegs(iyseg)%bdsp
                  xcnt = cama_gather%xsegs(ixseg)%cnt
                  ycnt = cama_gather%ysegs(iyseg)%cnt

                  allocate (sbuf (xcnt,ycnt))
                  sbuf = IOVar%blk(iblk,jblk)%val(xdsp+1:xdsp+xcnt,ydsp+1:ydsp+ycnt)

                  smesg = (/p_iam_glb, ixseg, iyseg/)
                  call mpi_send (smesg, 3, MPI_INTEGER, &
                     p_root, 10011, p_comm_glb, p_err) 
                  call mpi_send (sbuf, xcnt*ycnt, MPI_DOUBLE, &
                     p_root, 10011, p_comm_glb, p_err)

                  deallocate (sbuf)

               end if
            end do
         end do
      end if

      if (allocated(filter)) deallocate(filter)
      if (allocated(vectmp)) deallocate(vectmp)

   END SUBROUTINE colm2cama_real8

   ! -----
   SUBROUTINE cama2colm_real8 (MasterVar, IOVar, WorkerVar)
      !=======================================================================
      ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
      !=======================================================================     
      use precision
      use mod_namelist
      use timemanager
      use spmd_task
      use mod_block
      use mod_data_type
      use mod_landpatch
      use mod_mapping_pset2grid
      use mod_colm_debug
      USE MOD_TimeInvariants, only : patchtype
      use mod_grid
      !use GlobalVars, only : spval
      IMPLICIT NONE
      real(r8),                  INTENT(in)    :: MasterVar (:,:)
      type(block_data_real8_2d), INTENT(inout) :: IOVar
      REAL(r8),                  intent(inout) :: WorkerVar (:)

      integer :: xblk, yblk, xloc, yloc
      integer :: iblk, jblk, idata, ixseg, iyseg
      integer :: rmesg(2), smesg(2), isrc, iproc
      real(r8), allocatable :: rbuf(:,:), sbuf(:,:), vdata(:,:)
      integer :: xdsp, ydsp, xcnt, ycnt

      !call ghist%define_by_ndims (NXIN, NYIN)
      !call mp2g_hist%build (landpatch, ghist)
      !call set_segment_info (ghist)
      if (p_is_master) then
         do iyseg = 1, cama_gather%nyseg
            do ixseg = 1, cama_gather%nxseg
               iblk = cama_gather%xsegs(ixseg)%blk
               jblk = cama_gather%ysegs(iyseg)%blk

               IF (gblock%pio(iblk,jblk) >= 0) THEN
                  xdsp = cama_gather%xsegs(ixseg)%gdsp
                  ydsp = cama_gather%ysegs(iyseg)%gdsp
                  xcnt = cama_gather%xsegs(ixseg)%cnt
                  ycnt = cama_gather%ysegs(iyseg)%cnt

                  allocate (sbuf (xcnt,ycnt))
                  sbuf = MasterVar (xdsp+1:xdsp+xcnt, ydsp+1:ydsp+ycnt)
                  smesg = (/ixseg, iyseg/)
                  call mpi_send (smesg, 2, MPI_INTEGER, &
                     gblock%pio(iblk,jblk), 10000, p_comm_glb, p_err) 
                  call mpi_send (sbuf, xcnt*ycnt, MPI_DOUBLE, &
                     gblock%pio(iblk,jblk), 10000, p_comm_glb, p_err)
                  deallocate (sbuf)
               ENDIF
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
               iblk = cama_gather%xsegs(ixseg)%blk
               jblk = cama_gather%ysegs(iyseg)%blk
               xdsp = cama_gather%xsegs(ixseg)%bdsp
               ydsp = cama_gather%ysegs(iyseg)%bdsp
               xcnt = cama_gather%xsegs(ixseg)%cnt
               ycnt = cama_gather%ysegs(iyseg)%cnt

               allocate (rbuf(xcnt,ycnt))
               call mpi_recv (rbuf, xcnt*ycnt, MPI_DOUBLE, &
                  p_root, 10000, p_comm_glb, p_stat, p_err)
               IOVar%blk(iblk,jblk)%val(xdsp+1:xdsp+xcnt,ydsp+1:ydsp+ycnt)= rbuf
               deallocate (rbuf)
            ELSE
               exit
            ENDIF
         end do
      endif

      CALL mg2p_cama%map_aweighted (IOVar, WorkerVar)

   END SUBROUTINE cama2colm_real8

#endif

END module MOD_CaMa_Variables
! ----- EOP ---------
