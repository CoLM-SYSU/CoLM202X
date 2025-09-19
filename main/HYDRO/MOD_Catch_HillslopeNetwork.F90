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
   type :: hillslope_network_type
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
      real(r8), pointer :: fldprof (:,:) ! flood area profile
   CONTAINS
      final :: hillslope_network_free_mem
   END type hillslope_network_type

   integer :: nfldstep

CONTAINS

   ! ----------
   SUBROUTINE hillslope_network_init (ne, elmindex, hillslope_network)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_UserDefFun
   IMPLICIT NONE

   integer, intent(in) :: ne
   integer, intent(in) :: elmindex (:)
   type(hillslope_network_type), pointer :: hillslope_network(:)

   ! Local Variables
   character(len=256) :: hillslope_network_file

   integer :: maxnumhru, ie, nhru, hs, i, j
   integer :: iworker, mesg(2), nrecv, irecv, isrc, idest

   integer , allocatable :: eid (:)

   integer , allocatable :: nhru_all(:), nhru_in_bsn(:)

   integer , allocatable :: indxhru (:,:)
   real(r8), allocatable :: handhru (:,:)
   real(r8), allocatable :: elvahru (:,:)
   real(r8), allocatable :: plenhru (:,:)
   real(r8), allocatable :: lfachru (:,:)
   integer , allocatable :: nexthru (:,:)
   real(r8), allocatable :: fldstep (:,:,:)

   integer , allocatable :: icache (:,:)
   real(r8), allocatable :: rcache (:,:), rcache2(:,:,:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      hillslope_network_file = DEF_CatchmentMesh_data

      IF (p_is_master) THEN

         CALL ncio_read_serial (hillslope_network_file, 'basin_numhru',        nhru_all)
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_index',      indxhru)
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_hand',       handhru)
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_elva',       elvahru)
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_pathlen',    plenhru)
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_facelen',    lfachru)
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_downstream', nexthru)
         CALL ncio_read_serial (hillslope_network_file, 'hydrounit_flood_step', fldstep)

      ENDIF

      IF (p_is_master) maxnumhru = size(indxhru,1)
      IF (p_is_master) nfldstep  = size(fldstep,1)
#ifdef USEMPI
      CALL mpi_bcast (maxnumhru, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (nfldstep,  1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      hillslope_network => null()

      IF (p_is_master) THEN
#ifdef USEMPI
         DO iworker = 1, p_np_worker
            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            nrecv = mesg(2)
            isrc  = mesg(1)

            IF (nrecv > 0) THEN

               allocate (eid (nrecv))
               allocate (nhru_in_bsn (nrecv))
               allocate (icache (maxnumhru,nrecv))
               allocate (rcache (maxnumhru,nrecv))

               CALL mpi_recv (eid, nrecv, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               idest = isrc

               nhru_in_bsn = nhru_all(eid)
               CALL mpi_send (nhru_in_bsn, nrecv, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

               DO irecv = 1, nrecv
                  icache(:,irecv) = indxhru(:,eid(irecv))
               ENDDO
               CALL mpi_send (icache, maxnumhru*nrecv, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

               DO irecv = 1, nrecv
                  rcache(:,irecv) = handhru(:,eid(irecv))
               ENDDO
               CALL mpi_send (rcache, maxnumhru*nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

               DO irecv = 1, nrecv
                  rcache(:,irecv) = elvahru(:,eid(irecv))
               ENDDO
               CALL mpi_send (rcache, maxnumhru*nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

               DO irecv = 1, nrecv
                  rcache(:,irecv) = plenhru(:,eid(irecv))
               ENDDO
               CALL mpi_send (rcache, maxnumhru*nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

               DO irecv = 1, nrecv
                  rcache(:,irecv) = lfachru(:,eid(irecv))
               ENDDO
               CALL mpi_send (rcache, maxnumhru*nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

               DO irecv = 1, nrecv
                  icache(:,irecv) = nexthru(:,eid(irecv))
               ENDDO
               CALL mpi_send (icache, maxnumhru*nrecv, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

               allocate (rcache2 (nfldstep,maxnumhru,nrecv))

               DO irecv = 1, nrecv
                  rcache2(:,:,irecv) = fldstep(:,:,eid(irecv))
               ENDDO
               CALL mpi_send (rcache2, nfldstep*maxnumhru*nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (eid)
               deallocate (nhru_in_bsn)
               deallocate (icache)
               deallocate (rcache)
               deallocate (rcache2)
            ENDIF
         ENDDO
#else
         IF (ne > 0) THEN

            allocate (nhru_in_bsn      (ne))
            allocate (icache (maxnumhru,ne))
            allocate (rcache (maxnumhru,ne))

            nhru_in_bsn = nhru_all(elmindex)

            DO ie = 1, ne
               icache(:,ie) = indxhru(:,elmindex(ie))
            ENDDO
            indxhru = icache

            DO ie = 1, ne
               rcache(:,ie) = handhru(:,elmindex(ie))
            ENDDO
            handhru = rcache

            DO ie = 1, ne
               rcache(:,ie) = elvahru(:,elmindex(ie))
            ENDDO
            elvahru = rcache

            DO ie = 1, ne
               rcache(:,ie) = plenhru(:,elmindex(ie))
            ENDDO
            plenhru = rcache

            DO ie = 1, ne
               rcache(:,ie) = lfachru(:,elmindex(ie))
            ENDDO
            lfachru = rcache

            DO ie = 1, ne
               icache(:,ie) = nexthru(:,elmindex(ie))
            ENDDO
            nexthru = icache

            deallocate (icache)
            deallocate (rcache)

            allocate (rcache2 (nfldstep,maxnumhru,ne))

            DO ie = 1, ne
               rcache2(:,:,ie) = fldstep(:,:,elmindex(ie))
            ENDDO
            fldstep = rcache2

            deallocate (rcache2)

         ENDIF
#endif
      ENDIF

      IF (p_is_worker) THEN

#ifdef USEMPI
         mesg(1:2) = (/p_iam_glb, ne/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (ne > 0) THEN

            CALL mpi_send (elmindex, ne, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)

            allocate (nhru_in_bsn (ne))
            CALL mpi_recv (nhru_in_bsn, ne, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (indxhru (maxnumhru,ne))
            CALL mpi_recv (indxhru, maxnumhru*ne, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (handhru (maxnumhru,ne))
            CALL mpi_recv (handhru, maxnumhru*ne, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (elvahru (maxnumhru,ne))
            CALL mpi_recv (elvahru, maxnumhru*ne, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (plenhru (maxnumhru,ne))
            CALL mpi_recv (plenhru, maxnumhru*ne, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (lfachru (maxnumhru,ne))
            CALL mpi_recv (lfachru, maxnumhru*ne, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (nexthru (maxnumhru,ne))
            CALL mpi_recv (nexthru, maxnumhru*ne, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (fldstep (nfldstep,maxnumhru,ne))
            CALL mpi_recv (fldstep, nfldstep*maxnumhru*ne, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
#endif

         IF (ne > 0) THEN
            allocate( hillslope_network (ne))
         ENDIF

         hs = 0

         DO ie = 1, ne

            nhru = count(indxhru(:,ie) >= 0)

            IF (nhru > 0) THEN

               IF (nhru /= nhru_in_bsn(ie)) THEN
                  write(*,*) 'Warning : numbers of hydro units from file mismatch!'
               ENDIF

               allocate (hillslope_network(ie)%ihru  (nhru))
               allocate (hillslope_network(ie)%indx  (nhru))
               allocate (hillslope_network(ie)%area  (nhru))
               allocate (hillslope_network(ie)%agwt  (nhru))
               allocate (hillslope_network(ie)%hand  (nhru))
               allocate (hillslope_network(ie)%elva  (nhru))
               allocate (hillslope_network(ie)%plen  (nhru))
               allocate (hillslope_network(ie)%flen  (nhru))
               allocate (hillslope_network(ie)%inext (nhru))

               hillslope_network(ie)%indx = indxhru(1:nhru,ie)
               hillslope_network(ie)%hand = handhru(1:nhru,ie)         ! m
               hillslope_network(ie)%elva = elvahru(1:nhru,ie)         ! m
               hillslope_network(ie)%plen = plenhru(1:nhru,ie) * 1.0e3 ! km to m
               hillslope_network(ie)%flen = lfachru(1:nhru,ie) * 1.0e3 ! km to m

               hillslope_network(ie)%ihru = (/ (i, i = hs+1, hs+nhru) /)

               allocate (hillslope_network(ie)%fldprof (nfldstep,nhru))

               DO i = 1, nhru
                  DO j = 1, nfldstep
                     IF (j == 1) THEN
                        hillslope_network(ie)%fldprof(j,i) = max(fldstep(j,i,ie) * (j-0.5)/nfldstep, 1.e-3)
                     ELSE
                        hillslope_network(ie)%fldprof(j,i) = hillslope_network(ie)%fldprof(j-1,i) &
                           + (fldstep(j,i,ie) - fldstep(j-1,i,ie)) * (j-0.5)/nfldstep
                     ENDIF
                  ENDDO
               ENDDO

               DO i = 1, nhru
                  IF (nexthru(i,ie) >= 0) THEN
                     j = findloc_ud(indxhru(1:nhru,ie) == nexthru(i,ie))
                     hillslope_network(ie)%inext(i) = j
                  ELSE
                     hillslope_network(ie)%inext(i) = -1
                  ENDIF
               ENDDO

            ELSE
               ! for lake
               hillslope_network(ie)%ihru  => null()
               hillslope_network(ie)%indx  => null()
               hillslope_network(ie)%area  => null()
               hillslope_network(ie)%agwt  => null()
               hillslope_network(ie)%hand  => null()
               hillslope_network(ie)%elva  => null()
               hillslope_network(ie)%plen  => null()
               hillslope_network(ie)%flen  => null()
               hillslope_network(ie)%inext => null()
            ENDIF

            hillslope_network(ie)%nhru = nhru_in_bsn(ie)
            hs = hs + nhru_in_bsn(ie)

         ENDDO

      ENDIF

      IF (allocated(nhru_all   )) deallocate(nhru_all   )
      IF (allocated(nhru_in_bsn)) deallocate(nhru_in_bsn)

      IF (allocated(indxhru)) deallocate(indxhru)
      IF (allocated(handhru)) deallocate(handhru)
      IF (allocated(elvahru)) deallocate(elvahru)
      IF (allocated(plenhru)) deallocate(plenhru)
      IF (allocated(lfachru)) deallocate(lfachru)
      IF (allocated(nexthru)) deallocate(nexthru)
      IF (allocated(fldstep)) deallocate(fldstep)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
      IF (p_is_master) write(*,'(A)') 'Read hillslope network information done.'
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE hillslope_network_init

   ! ---------
   SUBROUTINE hillslope_network_free_mem (this)

   IMPLICIT NONE
   type(hillslope_network_type) :: this

      IF (associated(this%ihru   )) deallocate(this%ihru   )
      IF (associated(this%indx   )) deallocate(this%indx   )
      IF (associated(this%area   )) deallocate(this%area   )
      IF (associated(this%agwt   )) deallocate(this%agwt   )
      IF (associated(this%hand   )) deallocate(this%hand   )
      IF (associated(this%elva   )) deallocate(this%elva   )
      IF (associated(this%plen   )) deallocate(this%plen   )
      IF (associated(this%flen   )) deallocate(this%flen   )
      IF (associated(this%inext  )) deallocate(this%inext  )
      IF (associated(this%fldprof)) deallocate(this%fldprof)

   END SUBROUTINE hillslope_network_free_mem

END MODULE MOD_Catch_HillslopeNetwork
#endif
