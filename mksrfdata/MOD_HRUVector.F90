#include <define.h>

#if (defined CATCHMENT)
MODULE MOD_HRUVector

!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    Address of Data associated with HRU.
!
!    To output a vector, Data is gathered from worker processes directly to
!    master.  "hru_data_address" stores information on how to reorganize data
!    gathered.  The output data in vector is sorted by global element index
!    (i.e. catchment index)
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_DataType
   IMPLICIT NONE

   integer :: totalnumhru
   type(pointer_int32_1d), allocatable :: hru_data_address (:)

   integer*8, allocatable :: eindx_hru (:)
   integer,   allocatable :: htype_hru (:)

CONTAINS

   ! --------
   SUBROUTINE hru_vector_init

   USE MOD_SPMD_Task
   USE MOD_Utils
   USE MOD_Mesh
   USE MOD_LandElm
   USE MOD_LandHRU
   USE MOD_LandPatch
   USE MOD_ElmVector
#ifdef CROP
   USE MOD_LandCrop
#endif
   IMPLICIT NONE

   ! Local Variables
   integer   :: mesg(2), iwork, isrc, ndata

   integer, allocatable :: nhru_bsn (:)
   integer, allocatable :: nhru_bsn_glb (:)
   integer, allocatable :: rbuff (:)

   integer, allocatable :: hru_dsp_glb (:)
   integer :: ielm, i, ielm_glb

   integer :: nhru, nelm, hru_dsp_loc

      IF (p_is_worker) THEN

         CALL elm_hru%build (landelm, landhru,   use_frac = .true.)

         CALL hru_patch%build (landhru, landpatch, use_frac = .true.)

         IF (numelm > 0) THEN
            allocate (nhru_bsn (numelm))
            nhru_bsn = elm_hru%subend - elm_hru%substt + 1
         ENDIF

#ifdef USEMPI
         mesg = (/p_iam_glb, numelm/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numelm > 0) THEN
            CALL mpi_send (nhru_bsn, numelm, MPI_INTEGER, p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
#endif
      ENDIF

      IF (p_is_master) THEN

         allocate (hru_data_address (0:p_np_worker-1))

         allocate (nhru_bsn_glb (totalnumelm))

#ifdef USEMPI
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate (rbuff (ndata))

               CALL mpi_recv (rbuff, ndata, MPI_INTEGER, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

               nhru_bsn_glb(elm_data_address(p_itis_worker(isrc))%val) = rbuff

               IF (sum(rbuff) > 0) THEN
                  allocate(hru_data_address(p_itis_worker(isrc))%val (sum(rbuff)))
               ENDIF

               deallocate(rbuff)
            ENDIF
         ENDDO
#else
         nhru_bsn_glb(elm_data_address(0)%val) = nhru_bsn
         IF (sum(nhru_bsn) > 0) THEN
            allocate(hru_data_address(0)%val (sum(nhru_bsn)))
         ENDIF
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN

         totalnumhru = sum(nhru_bsn_glb)

         allocate (hru_dsp_glb (totalnumelm))
         hru_dsp_glb(1) = 0
         DO ielm = 2, totalnumelm
            hru_dsp_glb(ielm) = hru_dsp_glb(ielm-1) + nhru_bsn_glb(ielm-1)
         ENDDO

         DO iwork = 0, p_np_worker-1
            IF (allocated(elm_data_address(iwork)%val)) THEN
               nelm = size(elm_data_address(iwork)%val)
               hru_dsp_loc = 0
               DO ielm = 1, nelm
                  ielm_glb = elm_data_address(iwork)%val(ielm)
                  nhru = nhru_bsn_glb(ielm_glb)
                  IF (nhru > 0) THEN
                     hru_data_address(iwork)%val (hru_dsp_loc+1:hru_dsp_loc+nhru) = &
                        (/ (i, i = hru_dsp_glb(ielm_glb)+1, hru_dsp_glb(ielm_glb)+nhru) /)
                     hru_dsp_loc = hru_dsp_loc + nhru
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (totalnumhru, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

#ifdef USEMPI
      IF (p_is_worker) THEN
         mesg = (/p_iam_glb, numhru/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numhru > 0) THEN
            CALL mpi_send (landhru%settyp, numhru, MPI_INTEGER, p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
      ENDIF
#endif

      IF (p_is_master) THEN

         allocate (eindx_hru (totalnumhru))

         DO ielm = 1, totalnumelm
            eindx_hru(hru_dsp_glb(ielm)+1:hru_dsp_glb(ielm)+nhru_bsn_glb(ielm)) = &
               eindex_glb(ielm)
         ENDDO

         allocate (htype_hru (totalnumhru))

#ifdef USEMPI
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate (rbuff (ndata))

               CALL mpi_recv (rbuff, ndata, MPI_INTEGER, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
               htype_hru(hru_data_address(p_itis_worker(isrc))%val) = rbuff

               deallocate(rbuff)
            ENDIF
         ENDDO
#else
         htype_hru(hru_data_address(0)%val) = landhru%settyp
#endif

         ! To distinguish between lake HRUs and hillslopes, the program sets the
         ! type of lake HRUs as a negative number.
         ! Set it as a positive number for output.
         htype_hru = abs(htype_hru)

      ENDIF

#ifdef USEMPI
      CALL mpi_bcast (totalnumhru, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

   END SUBROUTINE hru_vector_init

   ! ----------
   SUBROUTINE hru_vector_final ()

   IMPLICIT NONE

      IF (allocated(hru_data_address))   deallocate (hru_data_address)
      IF (allocated(eindx_hru)) deallocate (eindx_hru)
      IF (allocated(htype_hru)) deallocate (htype_hru)

   END SUBROUTINE hru_vector_final

END MODULE MOD_HRUVector
#endif
