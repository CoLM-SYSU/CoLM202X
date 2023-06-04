#include <define.h>

#ifdef LULC_IGBP_PC

MODULE MOD_LandPC

   !------------------------------------------------------------------------------------
   ! DESCRIPTION:
   !
   !    Build pixelset "landpc" (Plant Community).
   !
   !    In CoLM, the global/regional area is divided into a hierarchical structure:
   !    1. If GRIDBASED or UNSTRUCTURED is defined, it is
   !       ELEMENT >>> PATCH
   !    2. If CATCHMENT is defined, it is
   !       ELEMENT >>> HRU >>> PATCH
   !    If Plant Function Type classification is used, PATCH is further divided into PFT.
   !    If Plant Community classification is used,     PATCH is further divided into PC.
   ! 
   !    "landpc" refers to pixelset PC.
   !
   ! Created by Shupeng Zhang, May 2023
   !    porting codes from Hua Yuan's OpenMP version to MPI parallel version.
   !------------------------------------------------------------------------------------

   USE MOD_Pixelset
   IMPLICIT NONE

   ! ---- Instance ----
   INTEGER :: numpc
   TYPE(pixelset_type) :: landpc

   INTEGER, allocatable :: patch2pc(:)    !projection from patch to PC
   INTEGER, allocatable :: pc2patch(:)    !projection from PC to patch

CONTAINS

   ! -------------------------------
   SUBROUTINE landpc_build (lc_year)

      USE MOD_Precision
      USE MOD_SPMD_Task
      USE MOD_LandPatch
      USE MOD_Const_LC
      IMPLICIT NONE

      INTEGER, intent(in) :: lc_year

      ! Local Variables
      INTEGER  :: ipatch, ipc, npc, m, npc_glb

      IF (p_is_master) THEN
         write(*,'(A)') 'Making land plant community tiles :'
      ENDIF

#ifdef SinglePoint
      IF (patchtypes(SITE_landtype) == 0) THEN
         numpc = 1
         allocate (patch2pc (1))
         allocate (pc2patch (1))
         patch2pc(1) = 1
         pc2patch(1) = 1

         allocate (landpc%eindex (numpc))
         allocate (landpc%settyp (numpc))
         allocate (landpc%ipxstt (numpc))
         allocate (landpc%ipxend (numpc))
         allocate (landpc%ielm   (numpc))

         landpc%eindex(1) = 1
         landpc%ipxstt(1) = 1
         landpc%ipxend(1) = 1
         landpc%settyp(1) = SITE_landtype
         landpc%ielm  (1) = 1
      ELSE
         numpc = 0
      ENDIF

      landpc%nset = numpc
      CALL landpc%set_vecgs

      RETURN
#endif

      if (p_is_worker) then

         numpc = 0

         DO ipatch = 1, numpatch
            m = landpatch%settyp(ipatch)
            IF (patchtypes(m) == 0) THEN
               numpc = numpc + 1
            ENDIF
         ENDDO

         IF (numpc > 0) THEN

            allocate (pc2patch (numpc   ))
            allocate (patch2pc (numpatch))

            patch2pc(:) = -1

            allocate (landpc%eindex (numpc))
            allocate (landpc%settyp (numpc))
            allocate (landpc%ipxstt (numpc))
            allocate (landpc%ipxend (numpc))
            allocate (landpc%ielm   (numpc))

            npc = 0
            DO ipatch = 1, numpatch
               m = landpatch%settyp(ipatch)
               IF (patchtypes(m) == 0) THEN

                  npc = npc + 1

                  landpc%ielm  (npc) = landpatch%ielm  (ipatch)
                  landpc%eindex(npc) = landpatch%eindex(ipatch)
                  landpc%ipxstt(npc) = landpatch%ipxstt(ipatch)
                  landpc%ipxend(npc) = landpatch%ipxend(ipatch)
                  landpc%settyp(npc) = m

                  pc2patch(npc) = ipatch
                  patch2pc(ipatch) = npc

               ENDIF
            ENDDO
         ENDIF
      ENDIF

      landpc%nset = numpc
      CALL landpc%set_vecgs

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numpc, npc_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', npc_glb, ' plant community tiles.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numpc, ' plant community tiles.'
#endif

   END SUBROUTINE landpc_build

   ! ----------------------
   SUBROUTINE map_patch_to_pc

      USE MOD_SPMD_Task
      USE MOD_LandPatch
      USE MOD_Const_LC
      IMPLICIT NONE

      INTEGER :: ipatch, ipc

      IF (p_is_worker) THEN

         IF ((numpatch <= 0) .or. (numpc <= 0)) return

         IF (allocated(pc2patch)) deallocate(pc2patch)
         IF (allocated(patch2pc)) deallocate(patch2pc)

         allocate (pc2patch (numpc   ))
         allocate (patch2pc (numpatch))

         ipc = 0
         DO ipatch = 1, numpatch
            IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
               ipc = ipc + 1
               patch2pc(ipatch) = ipc
               pc2patch(ipc) = ipatch
            ELSE
               patch2pc(ipatch) = -1
            ENDIF
         ENDDO

      ENDIF

   END SUBROUTINE map_patch_to_pc

END MODULE MOD_LandPC

#endif
