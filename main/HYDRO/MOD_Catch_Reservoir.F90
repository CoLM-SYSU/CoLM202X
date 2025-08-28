#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_Reservoir
!-----------------------------------------------------------------------
! DESCRIPTION:
!
!    Reservoir module in catchment mesh.
!
! Created by Shupeng Zhang, July 2025
!-----------------------------------------------------------------------

   integer, allocatable :: reservoir_id (:)

CONTAINS

   ! -------
   SUBROUTINE readin_reservoir_data ( )

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Utils
   USE MOD_Catch_BasinNetwork,     only : numbasin
   USE MOD_Catch_RiverLakeNetwork, only : lake_id

   IMPLICIT NONE

   ! Local variables
   character(len=256) :: filename
   integer :: funit, nline, iline, maxlakeid, thisid, thislake, laktyp, thisresv, numlake, ierr
   integer :: mesg(2), isrc, iwork

   integer, allocatable :: allresv (:), senddata(:), recvdata(:)
   logical, allocatable :: lakemask(:)

      IF (p_is_master) THEN

         filename = trim(DEF_dir_runtime)//'HydroLAKES_Reservoir.txt'

         funit = 101
         open (funit, file = trim(filename), action = 'READ')

         nline = 0
         maxlakeid = 0
         rewind (funit)
         read(funit,*,iostat=ierr) ! skip first header line
         DO WHILE (.true.)
            read(funit,*,iostat=ierr) thisid
            IF (ierr /= 0) exit
            nline = nline + 1
            maxlakeid = max(maxlakeid, thisid)
         ENDDO

         allocate (allresv (maxlakeid))
         allresv(:) = 0

         rewind (funit)
         read(funit,*,iostat=ierr)
         DO iline = 1, nline
            read(funit,*) thislake, laktyp, thisresv
            allresv(thislake) = thisresv
         ENDDO

         close (funit)

      ENDIF

      IF (p_is_worker) THEN
         allocate (reservoir_id (numbasin))
         reservoir_id(:) = 0
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN

         IF (numbasin > 0) THEN
            allocate (lakemask (numbasin))
            lakemask = lake_id > 0
            numlake = count(lakemask)
         ELSE
            numlake = 0
         ENDIF

         mesg = (/p_iam_glb, numlake/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (numlake > 0) THEN

            allocate (senddata (numlake))
            senddata = pack(lake_id, lakemask)
            CALL mpi_send (senddata, numlake, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)

            allocate (recvdata (numlake))
            CALL mpi_recv (recvdata, numlake, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            CALL unpack_inplace (recvdata, lakemask, reservoir_id)

            deallocate (senddata)
            deallocate (recvdata)
         ENDIF

      ENDIF

      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1

            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc    = mesg(1)
            numlake = mesg(2)
            IF (numlake > 0) THEN

               allocate(recvdata (numlake))
               CALL mpi_recv (recvdata, numlake, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               allocate(senddata (numlake))
               senddata = allresv(recvdata)
               CALL mpi_send (senddata, numlake, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_err)

               deallocate (recvdata)
               deallocate (senddata)
            ENDIF

         ENDDO
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else

      DO ibasin = 1, numbasin
         IF (lake_id(ibasin) > 0) THEN
            reservoir_id(ibasin) = allresv(lake_id(ibasin))
         ENDIF
      ENDDO

#endif

      IF (allocated(lakemask)) deallocate (lakemask)
      IF (allocated(allresv )) deallocate (allresv )

   END SUBROUTINE readin_reservoir_data

END MODULE MOD_Catch_Reservoir
#endif
