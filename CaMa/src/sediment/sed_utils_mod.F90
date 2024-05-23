MODULE sed_utils_mod
!==========================================================
!* PURPOSE: Shared ulitity functions/subroutines for CaMa-Flood sediment scheme
! (C) M. Hatono (Hiroshima Univ)  Jan 2023
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
CONTAINS
   !####################################################################
   !-- splitchar           : same function as splitting characters in CaMa
   !-- sed_diag_average    : 
   !-- sed_diag_reset      :
   !####################################################################
   SUBROUTINE splitchar(allvars,vnames)
   ! same function as splitting characters in CaMa
   USE PARKIND1,                only: JPIM
   IMPLICIT NONE
   SAVE
   character(len=256), intent(in)  :: allvars
   character(len=256), intent(out) :: vnames(:)
   integer(kind=JPIM)              :: nvarsout, j0, j
   character(len=256)              :: ctmp

      nvarsout = 0
      j0 = 1
      DO j = 1, len(trim(allvars))
         IF ( (j>j0) .and. (allvars(j:j).eq.',') ) THEN
            ctmp = trim(adjustl(allvars(j0:j-1)))
            IF ( len(ctmp) > 0 ) THEN
               nvarsout = nvarsout + 1
               vnames(nvarsout) = ctmp
            ENDIF
            j0 = j + 1
         ENDIF
      ENDDO

      ! last one
      IF ( j0 < len(trim(allvars)) ) THEN
         j = len(trim(allvars))
         ctmp = trim(adjustl(allvars(j0:j)))
         IF ( len(ctmp) > 0 ) THEN
            nvarsout = nvarsout + 1
            vnames(nvarsout) = ctmp
         ENDIF
      ENDIF
   END SUBROUTINE splitchar
   !==========================================================
   !+
   !==========================================================
   SUBROUTINE sed_diag_average
   USE yos_cmf_sed,             only: d2rivout_sed, d2rivvel_sed, sadd_riv
   IMPLICIT NONE

      d2rivout_sed(:) = d2rivout_sed(:) /dble(sadd_riv)
      d2rivvel_sed(:) = d2rivvel_sed(:) /dble(sadd_riv)
   END SUBROUTINE sed_diag_average
   !==========================================================
   !+
   !==========================================================
   SUBROUTINE sed_diag_reset
   USE PARKIND1,                only: JPRB
   USE YOS_CMF_PROG,            only: P2RIVSTO
   USE yos_cmf_sed,             only: d2rivsto_pre, d2rivout_sed, d2rivvel_sed, &
                                       sadd_riv, sadd_out, sedDT
   IMPLICIT NONE

      sadd_riv = 0
      d2rivout_sed(:) = 0._JPRB
      d2rivvel_sed(:) = 0._JPRB
      d2rivsto_pre(:) = P2RIVSTO(:,1)

      sadd_out = sadd_out + sedDT
   END SUBROUTINE sed_diag_reset
!####################################################################
END MODULE sed_utils_mod
