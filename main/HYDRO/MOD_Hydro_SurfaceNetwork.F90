#include <define.h>

#ifdef LATERAL_FLOW
MODULE MOD_Hydro_SurfaceNetwork
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
   TYPE :: surface_network_info_type
      INTEGER :: nhru
      INTEGER , pointer :: ihru (:)
      INTEGER , pointer :: indx (:)
      REAL(r8), pointer :: area (:)
      REAL(r8), pointer :: hand (:)
      REAL(r8), pointer :: elva (:)
      REAL(r8), pointer :: plen (:)
      REAL(r8), pointer :: flen (:)
      INTEGER , pointer :: inext(:)
   END TYPE surface_network_info_type

   ! -- Instance --
   TYPE(surface_network_info_type), pointer :: surface_network (:)
      
CONTAINS
   
   ! ----------
   SUBROUTINE surface_network_init ()

      USE MOD_SPMD_Task
      USE MOD_Namelist
      USE MOD_NetCDFSerial
      USE MOD_Mesh
      USE MOD_Pixel
      USE MOD_LandHRU
      USE MOD_Utils
      IMPLICIT NONE

      ! Local Variables
      CHARACTER(len=256) :: surface_network_file

      INTEGER :: numbasin, maxnumhru, ibasin, nhru, istt, iend, ihru, i, j, ipxl
      INTEGER :: iworker, mesg(2), nrecv, irecv, isrc, idest
   
      INTEGER , allocatable :: indxhru (:,:)
      REAL(r8), allocatable :: areahru (:,:)
      REAL(r8), allocatable :: handhru (:,:)
      REAL(r8), allocatable :: elvahru (:,:)
      REAL(r8), allocatable :: plenhru (:,:)
      REAL(r8), allocatable :: lfachru (:,:)
      INTEGER , allocatable :: nexthru (:,:)
      
      INTEGER , allocatable :: basinindex (:)
      INTEGER , allocatable :: icache (:,:)
      REAL(r8), allocatable :: rcache (:,:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      numbasin = numelm

      surface_network_file = DEF_CatchmentMesh_data 

      IF (p_is_master) THEN
         CALL ncio_read_serial (surface_network_file, 'hydrounit_index',      indxhru)
         CALL ncio_read_serial (surface_network_file, 'hydrounit_area',       areahru)
         CALL ncio_read_serial (surface_network_file, 'hydrounit_hand',       handhru)
         CALL ncio_read_serial (surface_network_file, 'hydrounit_elva',       elvahru)
         CALL ncio_read_serial (surface_network_file, 'hydrounit_pathlen',    plenhru)
         CALL ncio_read_serial (surface_network_file, 'hydrounit_facelen',    lfachru)
         CALL ncio_read_serial (surface_network_file, 'hydrounit_downstream', nexthru)
      ENDIF

      IF (p_is_master) maxnumhru = size(indxhru,1) 
#ifdef USEMPI
      CALL mpi_bcast (maxnumhru, 1, MPI_INTEGER, p_root, p_comm_glb, p_err)
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
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_err) 

         IF (numbasin > 0) THEN
            allocate (basinindex (numbasin))
            DO ibasin = 1, numbasin
               basinindex(ibasin) = mesh(ibasin)%indx
            ENDDO

            CALL mpi_send (basinindex, numbasin, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_err) 

            allocate (indxhru (maxnumhru,numbasin))
            CALL mpi_recv (indxhru, maxnumhru*numbasin, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (areahru (maxnumhru,numbasin))
            CALL mpi_recv (areahru, maxnumhru*numbasin, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (handhru (maxnumhru,numbasin))
            CALL mpi_recv (handhru, maxnumhru*numbasin, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (elvahru (maxnumhru,numbasin))
            CALL mpi_recv (elvahru, maxnumhru*numbasin, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (plenhru (maxnumhru,numbasin))
            CALL mpi_recv (plenhru, maxnumhru*numbasin, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (lfachru (maxnumhru,numbasin))
            CALL mpi_recv (lfachru, maxnumhru*numbasin, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (nexthru (maxnumhru,numbasin))
            CALL mpi_recv (nexthru, maxnumhru*numbasin, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
#endif

         IF (numbasin > 0) THEN
            allocate( surface_network (numbasin))
         ENDIF

         DO ibasin = 1, numbasin
               
            nhru = count(indxhru(:,ibasin) >= 0)
            surface_network(ibasin)%nhru = nhru

            IF (nhru > 0) THEN

               allocate (surface_network(ibasin)%indx  (nhru))
               allocate (surface_network(ibasin)%area  (nhru))
               allocate (surface_network(ibasin)%hand  (nhru))
               allocate (surface_network(ibasin)%elva  (nhru))
               allocate (surface_network(ibasin)%plen  (nhru))
               allocate (surface_network(ibasin)%flen  (nhru))
               allocate (surface_network(ibasin)%inext (nhru))

               surface_network(ibasin)%indx = indxhru(1:nhru,ibasin) 
               surface_network(ibasin)%area = areahru(1:nhru,ibasin) * 1.0e6 ! km^2 to m^2
               surface_network(ibasin)%hand = handhru(1:nhru,ibasin)         ! m
               surface_network(ibasin)%elva = elvahru(1:nhru,ibasin)         ! m
               surface_network(ibasin)%plen = plenhru(1:nhru,ibasin) * 1.0e3 ! km to m      
               surface_network(ibasin)%flen = lfachru(1:nhru,ibasin) * 1.0e3 ! km to m

               allocate (surface_network(ibasin)%ihru (nhru))
               istt = basin_hru%substt(ibasin)
               iend = basin_hru%subend(ibasin)
               surface_network(ibasin)%ihru = (/ (i, i = istt, iend) /)

               DO i = 1, nhru
                  IF (nexthru(i,ibasin) >= 0) THEN
                     j = findloc(indxhru(1:nhru,ibasin), nexthru(i,ibasin), dim=1)
                     surface_network(ibasin)%inext(i) = j 
                  ELSE
                     surface_network(ibasin)%inext(i) = -1
                  ENDIF
               ENDDO

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
      
   END SUBROUTINE surface_network_init
   
   ! ----------
   SUBROUTINE surface_network_final ()

      IMPLICIT NONE

      ! Local Variables
      INTEGER :: ibasin

      IF (associated(surface_network)) THEN
         DO ibasin = 1, size(surface_network)
            IF (associated(surface_network(ibasin)%ihru )) deallocate(surface_network(ibasin)%ihru )
            IF (associated(surface_network(ibasin)%indx )) deallocate(surface_network(ibasin)%indx )
            IF (associated(surface_network(ibasin)%area )) deallocate(surface_network(ibasin)%area )
            IF (associated(surface_network(ibasin)%hand )) deallocate(surface_network(ibasin)%hand )
            IF (associated(surface_network(ibasin)%elva )) deallocate(surface_network(ibasin)%elva )
            IF (associated(surface_network(ibasin)%plen )) deallocate(surface_network(ibasin)%plen )
            IF (associated(surface_network(ibasin)%flen )) deallocate(surface_network(ibasin)%flen )
            IF (associated(surface_network(ibasin)%inext)) deallocate(surface_network(ibasin)%inext)
         ENDDO

         deallocate(surface_network)
      ENDIF

   END SUBROUTINE surface_network_final

END MODULE MOD_Hydro_SurfaceNetwork
#endif
