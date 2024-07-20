MODULE CMF_CALC_LAKEIN_MOD
!==========================================================
!* PURPOSE: Lake-river coupling in ILS
!
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
   USE PARKIND1,                only: JPIM, JPRM, JPRB
CONTAINS
!####################################################################
! -- CMF_CALC_LAKEIN
! -- CMF_LAKEIN_AVE
! -- CMF_LAKEIN_AVERAGE
! -- CMF_RESET_LAKEIN
!####################################################################
   SUBROUTINE CMF_CALC_LAKEIN
   USE YOS_CMF_INPUT,      only: DT
   USE YOS_CMF_MAP,        only: NSEQALL, D2GRAREA
   USE YOS_CMF_PROG,       only: P2RIVSTO, P2FLDSTO
   USE YOS_CMF_DIAG,       only: D2STORGE
   USE YOS_CMF_ICI,        only: D2LAKEFRC, D2RUNIN
   IMPLICIT NONE
   !*** Local
   integer(KIND=JPIM)         :: ISEQ
   real(KIND=JPRB)            :: DRIVRIN, DFLDRIN
   !*** Lake Parameter
   real(KIND=JPRB)            :: RINDMP
      !================================================
      RINDMP = 1.D0

      DO ISEQ=1, NSEQALL
         DRIVRIN = P2RIVSTO(ISEQ,1) * D2LAKEFRC(ISEQ,1) / (RINDMP * 8.64D4) * DT !! m3
         DFLDRIN = P2FLDSTO(ISEQ,1) * D2LAKEFRC(ISEQ,1) / (RINDMP * 8.64D4) * DT !! m3
         P2RIVSTO(ISEQ,1) = P2RIVSTO(ISEQ,1) - DRIVRIN
         P2FLDSTO(ISEQ,1) = P2FLDSTO(ISEQ,1) - DFLDRIN
         D2RUNIN(ISEQ,1) = (DRIVRIN + DFLDRIN) / DT / D2GRAREA(ISEQ,1) * 1.D3 !! kg/m2/s

         D2STORGE(ISEQ,1)=P2RIVSTO(ISEQ,1)+P2FLDSTO(ISEQ,1)
      ENDDO

   END SUBROUTINE CMF_CALC_LAKEIN
   !####################################################################


   !####################################################################
   SUBROUTINE CMF_LAKEIN_AVE
   USE YOS_CMF_INPUT,      only: DT
   USE YOS_CMF_MAP,        only: NSEQALL
   USE YOS_CMF_ICI,        only: D2RUNIN, D2RUNIN_AVG
   IMPLICIT NONE
   !*** Local
   integer(KIND=JPIM)         :: ISEQ
   !================================================
      DO ISEQ=1, NSEQALL
         D2RUNIN_AVG(ISEQ,1)=D2RUNIN_AVG(ISEQ,1)+D2RUNIN(ISEQ,1)*DT
      ENDDO

   END SUBROUTINE CMF_LAKEIN_AVE
   !####################################################################


   !####################################################################
   SUBROUTINE CMF_LAKEIN_AVERAGE
   USE YOS_CMF_DIAG,       only: NADD
   USE YOS_CMF_MAP,        only: NSEQALL
   USE YOS_CMF_ICI,        only: D2RUNIN_AVG
   IMPLICIT NONE
   !*** Local
   integer(KIND=JPIM)         :: ISEQ
   !================================================
   DO ISEQ=1, NSEQALL
      D2RUNIN_AVG(ISEQ,1)=D2RUNIN_AVG(ISEQ,1)/DBLE(NADD)
   ENDDO

   END SUBROUTINE CMF_LAKEIN_AVERAGE
   !####################################################################

   !####################################################################
   SUBROUTINE CMF_RESET_LAKEIN
   USE YOS_CMF_ICI,        only: D2RUNIN_AVG
   IMPLICIT NONE
      !================================================
      D2RUNIN_AVG(:,:)=0._JPRB

   END SUBROUTINE CMF_RESET_LAKEIN
   !####################################################################

END MODULE CMF_CALC_LAKEIN_MOD
