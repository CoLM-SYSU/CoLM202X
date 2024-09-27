#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_HillslopeNetwork
!--------------------------------------------------------------------------------
! DESCRIPTION:
! 
!    Surface networks (hillslope bands): data and communication subroutines.
!
! Created by Shupeng Zhang, May 2023
!--------------------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE
   
   ! -- data type --
   type :: hillslope_network_info_type
      integer :: nhru
      integer , pointer :: ihru (:) ! location of HRU in global vector "landhru"
      integer , pointer :: indx (:) ! index of HRU
      real(r8), pointer :: area (:) ! area of HRU [m^2]
      real(r8), pointer :: agwt (:) ! water area only including (patchtype <= 2) [m^2]
      real(r8), pointer :: hand (:) ! height above nearest drainage [m]
      real(r8), pointer :: elva (:) ! elevation [m]
      real(r8), pointer :: plen (:) ! average drainage path length to downstream HRU [m]
      real(r8), pointer :: flen (:) ! interface length between this and downstream HRU [m]
      integer , pointer :: inext(:) ! location of next HRU in this basin
   END type hillslope_network_info_type

   ! -- Instance --
   type(hillslope_network_info_type), pointer :: hillslope_network (:)
      
CONTAINS
   
   ! ----------
   SUBROUTINE hillslope_network_init ()

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Mesh
   USE MOD_Pixel
   USE MOD_LandHRU
   USE MOD_LandPatch
   USE MOD_Vars_TimeInvariants, only : patchtype
   USE MOD_Utils
   USE MOD_UserDefFun
   IMPLICIT NONE

   ! Local Variables
   character(len=256) :: hillslope_network_file

   integer :: numbasin, maxnumhru, ibasin, nhru, hs, he, ihru, ipatch, ps, pe, i, j, ipxl
   integer :: iworker, mesg(2), nrecv, irecv, isrc, idest
   
   integer , allocatable :: indxhru (:,:)
   real(r8), allocatable :: areahru (:,:)
   real(r8), allocatable :: handhru (:,:)
   real(r8), allocatable :: elvahru (:,:)
   real(r8), allocatable :: plenhru (:,:)
   real(r8), allocatable :: lfachru (:,:)
   integer , allocatable :: nexthru (:,:)
   
   integer , allocatable :: basinindex (:)
   integer , allocatable :: icache (:,:)
   real(r8), allocatable :: rcache (:,:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      numbasin = numelm

      hillslope_network => null()

      hillslope_network_file = DEF_CatchmentMesh_data 

      IF (p_is_master) THEN
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_index',      indxhru)
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_area',       areahru)
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_hand',       handhru)
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_elva',       elvahru)
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_pathlen',    plenhru)
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_facelen',    lfachru)
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_downstream', nexthru)
      ENDIF

      IF (p_is_master) maxnumhru = size(indxhru,1) 
#ifdef USEMPI
      CALL mpi_bcast (maxnumhru, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
#ifdef USEMPI
         DO iworker = 1, p_np_worker
            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            nrecv = mesg(2)
            isrc  = mesg(1)

            IF (nrecv > 0) THEN
               
               allocate (basinindex (nrecv))
               allocate (icache (maxnumhru,nrecv))
               allocate (rcache (maxnumhru,nrecv))

               CALL mpi_recv (basinindex, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               idest = isrc

               DO irecv = 1, nrecv
                  icache(:,irecv) = indxhru(:,basinindex(irecv))
               ENDDO
               CALL mpi_send (icache, maxnumhru*nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache(:,irecv) = areahru(:,basinindex(irecv))
               ENDDO
               CALL mpi_send (rcache, maxnumhru*nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache(:,irecv) = handhru(:,basinindex(irecv))
               ENDDO
               CALL mpi_send (rcache, maxnumhru*nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache(:,irecv) = elvahru(:,basinindex(irecv))
               ENDDO
               CALL mpi_send (rcache, maxnumhru*nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache(:,irecv) = plenhru(:,basinindex(irecv))
               ENDDO
               CALL mpi_send (rcache, maxnumhru*nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache(:,irecv) = lfachru(:,basinindex(irecv))
               ENDDO
               CALL mpi_send (rcache, maxnumhru*nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  icache(:,irecv) = nexthru(:,basinindex(irecv))
               ENDDO
               CALL mpi_send (icache, maxnumhru*nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               deallocate (basinindex)
               deallocate (icache)
               deallocate (rcache)

            ENDIF
         ENDDO
#else
         IF (numbasin > 0) THEN

            allocate (basinindex (numbasin))
            allocate (icache (maxnumhru,numbasin))
            allocate (rcache (maxnumhru,numbasin))
            
            DO ibasin = 1, numbasin
               basinindex(ibasin) = mesh(ibasin)%indx
            ENDDO

            DO ibasin = 1, numbasin
               icache(:,ibasin) = indxhru(:,basinindex(ibasin))
            ENDDO
            indxhru = icache

            DO ibasin = 1, numbasin
               rcache(:,ibasin) = areahru(:,basinindex(ibasin))
            ENDDO
            areahru = rcache

            DO ibasin = 1, numbasin
               rcache(:,ibasin) = handhru(:,basinindex(ibasin))
            ENDDO
            handhru = rcache

            DO ibasin = 1, numbasin
               rcache(:,ibasin) = elvahru(:,basinindex(ibasin))
            ENDDO
            elvahru = rcache

            DO ibasin = 1, numbasin
               rcache(:,ibasin) = plenhru(:,basinindex(ibasin))
            ENDDO
            plenhru = rcache

            DO ibasin = 1, numbasin
               rcache(:,ibasin) = lfachru(:,basinindex(ibasin))
            ENDDO
            lfachru = rcache

            DO ibasin = 1, numbasin
               icache(:,ibasin) = nexthru(:,basinindex(ibasin))
            ENDDO
            nexthru = icache
               
            deallocate (basinindex)
            deallocate (icache)
            deallocate (rcache)

         ENDIF
#endif
      ENDIF

      IF (p_is_worker) THEN
               
#ifdef USEMPI
         mesg(1:2) = (/p_iam_glb,numbasin/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err) 

         IF (numbasin > 0) THEN
            allocate (basinindex (numbasin))
            DO ibasin = 1, numbasin
               basinindex(ibasin) = mesh(ibasin)%indx
            ENDDO

            CALL mpi_send (basinindex, numbasin, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err) 

            allocate (indxhru (maxnumhru,numbasin))
            CALL mpi_recv (indxhru, maxnumhru*numbasin, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (areahru (maxnumhru,numbasin))
            CALL mpi_recv (areahru, maxnumhru*numbasin, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (handhru (maxnumhru,numbasin))
            CALL mpi_recv (handhru, maxnumhru*numbasin, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (elvahru (maxnumhru,numbasin))
            CALL mpi_recv (elvahru, maxnumhru*numbasin, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (plenhru (maxnumhru,numbasin))
            CALL mpi_recv (plenhru, maxnumhru*numbasin, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (lfachru (maxnumhru,numbasin))
            CALL mpi_recv (lfachru, maxnumhru*numbasin, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (nexthru (maxnumhru,numbasin))
            CALL mpi_recv (nexthru, maxnumhru*numbasin, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
#endif

         IF (numbasin > 0) THEN
            allocate( hillslope_network (numbasin))
         ENDIF

         DO ibasin = 1, numbasin
               
            nhru = count(indxhru(:,ibasin) >= 0)
            hillslope_network(ibasin)%nhru = nhru

            IF (nhru > 0) THEN

               allocate (hillslope_network(ibasin)%ihru  (nhru))
               allocate (hillslope_network(ibasin)%indx  (nhru))
               allocate (hillslope_network(ibasin)%area  (nhru))
               allocate (hillslope_network(ibasin)%agwt  (nhru))
               allocate (hillslope_network(ibasin)%hand  (nhru))
               allocate (hillslope_network(ibasin)%elva  (nhru))
               allocate (hillslope_network(ibasin)%plen  (nhru))
               allocate (hillslope_network(ibasin)%flen  (nhru))
               allocate (hillslope_network(ibasin)%inext (nhru))

               hillslope_network(ibasin)%indx = indxhru(1:nhru,ibasin) 
               hillslope_network(ibasin)%area = areahru(1:nhru,ibasin) * 1.0e6 ! km^2 to m^2
               hillslope_network(ibasin)%hand = handhru(1:nhru,ibasin)         ! m
               hillslope_network(ibasin)%elva = elvahru(1:nhru,ibasin)         ! m
               hillslope_network(ibasin)%plen = plenhru(1:nhru,ibasin) * 1.0e3 ! km to m      
               hillslope_network(ibasin)%flen = lfachru(1:nhru,ibasin) * 1.0e3 ! km to m

               hs = basin_hru%substt(ibasin)
               he = basin_hru%subend(ibasin)
               hillslope_network(ibasin)%ihru = (/ (i, i = hs, he) /)

               DO i = 1, nhru
                  IF (nexthru(i,ibasin) >= 0) THEN
                     j = findloc_ud(indxhru(1:nhru,ibasin) == nexthru(i,ibasin))
                     hillslope_network(ibasin)%inext(i) = j 
                  ELSE
                     hillslope_network(ibasin)%inext(i) = -1
                  ENDIF
               ENDDO

               DO i = 1, nhru
                  hillslope_network(ibasin)%agwt(i) = 0
                  ps = hru_patch%substt(i+hs-1)
                  pe = hru_patch%subend(i+hs-1)
                  DO ipatch = ps, pe
                     IF (patchtype(ipatch) <= 2) THEN
                        DO ipxl = landpatch%ipxstt(ipatch), landpatch%ipxend(ipatch)
                           hillslope_network(ibasin)%agwt(i) = hillslope_network(ibasin)%agwt(i) &
                              + 1.0e6 * areaquad ( &
                              pixel%lat_s(mesh(ibasin)%ilat(ipxl)), pixel%lat_n(mesh(ibasin)%ilat(ipxl)), &
                              pixel%lon_w(mesh(ibasin)%ilon(ipxl)), pixel%lon_e(mesh(ibasin)%ilon(ipxl)) )
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO

            ELSE
               hillslope_network(ibasin)%ihru  => null() 
               hillslope_network(ibasin)%indx  => null()
               hillslope_network(ibasin)%area  => null()
               hillslope_network(ibasin)%agwt  => null()
               hillslope_network(ibasin)%hand  => null()
               hillslope_network(ibasin)%elva  => null()
               hillslope_network(ibasin)%plen  => null()
               hillslope_network(ibasin)%flen  => null()
               hillslope_network(ibasin)%inext => null()
            ENDIF
         ENDDO

      ENDIF 
         
      IF (allocated(indxhru)) deallocate(indxhru)  
      IF (allocated(areahru)) deallocate(areahru)
      IF (allocated(handhru)) deallocate(handhru)
      IF (allocated(elvahru)) deallocate(elvahru)
      IF (allocated(plenhru)) deallocate(plenhru)
      IF (allocated(lfachru)) deallocate(lfachru)
      IF (allocated(nexthru)) deallocate(nexthru)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
      IF (p_is_master) write(*,'(A)') 'Read surface network information done.'
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      
   END SUBROUTINE hillslope_network_init
   
   ! ----------
   SUBROUTINE hillslope_network_final ()

   IMPLICIT NONE

   ! Local Variables
   integer :: ibasin

      IF (associated(hillslope_network)) THEN
         DO ibasin = 1, size(hillslope_network)
            IF (associated(hillslope_network(ibasin)%ihru )) deallocate(hillslope_network(ibasin)%ihru )
            IF (associated(hillslope_network(ibasin)%indx )) deallocate(hillslope_network(ibasin)%indx )
            IF (associated(hillslope_network(ibasin)%area )) deallocate(hillslope_network(ibasin)%area )
            IF (associated(hillslope_network(ibasin)%agwt )) deallocate(hillslope_network(ibasin)%agwt )
            IF (associated(hillslope_network(ibasin)%hand )) deallocate(hillslope_network(ibasin)%hand )
            IF (associated(hillslope_network(ibasin)%elva )) deallocate(hillslope_network(ibasin)%elva )
            IF (associated(hillslope_network(ibasin)%plen )) deallocate(hillslope_network(ibasin)%plen )
            IF (associated(hillslope_network(ibasin)%flen )) deallocate(hillslope_network(ibasin)%flen )
            IF (associated(hillslope_network(ibasin)%inext)) deallocate(hillslope_network(ibasin)%inext)
         ENDDO

         deallocate(hillslope_network)
      ENDIF

   END SUBROUTINE hillslope_network_final

END MODULE MOD_Catch_HillslopeNetwork
#endif
