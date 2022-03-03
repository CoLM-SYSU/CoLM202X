#include <define.h>

#ifdef PC_CLASSIFICATION

MODULE mod_landpc

   USE mod_pixelset
   IMPLICIT NONE

   ! ---- Instance ----
   INTEGER :: numpc
   TYPE(pixelset_type) :: landpc

   INTEGER, allocatable :: patch2pc(:)    !projection from patch to PC
   INTEGER, allocatable :: pc2patch(:)    !projection from PC to patch

CONTAINS
   
   ! -------------------------------
   SUBROUTINE landpc_build ()

      USE precision
      USE spmd_task
      USE mod_landpatch
      USE LC_Const
      IMPLICIT NONE

      ! Local Variables
      INTEGER  :: ipatch, ipc, npc, m
      
      IF (p_is_master) THEN
         write(*,*) 'Making land plant communities :'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      if (p_is_worker) then

         numpc = 0

         DO ipatch = 1, numpatch
            m = landpatch%ltyp(ipatch)
            IF (patchtypes(m) == 0) THEN
               numpc = numpc + 1
            ENDIF
         ENDDO

         IF (numpc > 0) THEN

            allocate (pc2patch (numpc   ))
            allocate (patch2pc (numpatch))

            allocate (landpc%unum (numpc))
            allocate (landpc%iunt (numpc))
            allocate (landpc%ltyp (numpc))
            allocate (landpc%istt (numpc))
            allocate (landpc%iend (numpc))

            npc = 0
            DO ipatch = 1, numpatch
               m = landpatch%ltyp(ipatch)
               IF (patchtypes(m) == 0) THEN

                  npc = npc + 1
                        
                  landpc%iunt(npc) = landpatch%iunt(ipatch)
                  landpc%unum(npc) = landpatch%unum(ipatch)
                  landpc%istt(npc) = landpatch%istt(ipatch)
                  landpc%iend(npc) = landpatch%iend(ipatch)
                  landpc%ltyp(npc) = m

                  pc2patch(npc) = ipatch
                  patch2pc(ipatch) = npc

               ENDIF
            ENDDO
         ENDIF
      ENDIF

      landpc%nset = numcell
      CALL landpc%set_vecgs

   END SUBROUTINE landpc_build

   ! ----------------------
   SUBROUTINE map_patch_to_pc

      USE spmd_task
      USE mod_landpatch
      USE LC_Const
      IMPLICIT NONE

      INTEGER :: ipatch, ipc

      IF (p_is_worker) THEN

         IF ((numpatch <= 0) .or. (numpc <= 0)) return
         
         allocate (pc2patch (numpc   ))
         allocate (patch2pc (numpatch))
      
         ipc = 0
         DO ipatch = 1, numpatch
            IF (patchtypes(landpatch%ltyp(ipatch)) == 0) THEN
               ipc = ipc + 1
               patch2pc(ipatch) = ipc
               pc2patch(ipc) = ipatch
            ELSE
               patch2pc(ipatch) = -1
            ENDIF
         ENDDO
         
         write(*,'(I10,A14,I5)') numpft, ' pcs on worker', p_iam_glb

      ENDIF

   END SUBROUTINE map_patch_to_pc

END MODULE mod_landpc

#endif
