MODULE cmf_calc_sedpar_mod
!==========================================================
!* PURPOSE: parameters for sediment transport
! (C) M.Hatono  (Hiroshima-U)  May 2021
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not USE this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
   USE PARKIND1,                only: JPIM, JPRB

CONTAINS
   !####################################################################
   ! -- csVel
   ! -- sVel
   ! -- susVel
   !####################################################################
   FUNCTION calc_settingVelocity() result(setVel)
   USE YOS_CMF_INPUT,           only:  PGRV
   USE yos_cmf_sed,             only:  nsed, psedD, pset, pwatD, sDiam, visKin
   IMPLICIT NONE
   SAVE
   real(kind=JPRB)                 ::  setVel(nsed) ! setting velocity [m/s]
   real(kind=JPRB)                 ::  sTmp(nsed)

      sTmp(:) = 6.d0 * visKin / sDiam(:)
      setVel(:) = pset * ( sqrt( 2.d0/3.d0*(psedD-pwatD)/pwatD*PGRV*sDiam(:) &
                                 + sTmp(:)*sTmp(:) ) - sTmp(:) )
   END FUNCTION calc_settingVelocity
   !=====================================================

   FUNCTION calc_criticalShearVelocity(diam) result(csVel)
   IMPLICIT NONE
   SAVE
   real(kind=JPRB)                 ::  csVel ! critical shear velocity[(cm/s)^2]
   real(kind=JPRB), intent(in)     ::  diam ![m]
   real(kind=JPRB)                 ::  cA, cB
      !========
      cB = 1.d0
      IF ( diam >= 0.00303d0 ) THEN
         cA = 80.9d0
      ELSE IF ( diam >= 0.00118d0 ) THEN
         cA = 134.6d0
         cB = 31.d0 / 32.d0
      ELSE IF ( diam >= 0.000565d0 ) THEN
         cA = 55.d0
      ELSE IF ( diam >= 0.000065d0 ) THEN
         cA = 8.41d0
         cB = 11.d0 / 32.d0
      ELSE
         cA = 226.d0
      ENDIF
      
      csVel = cA * ( diam*100.d0 ) ** cB
      RETURN
   END FUNCTION calc_criticalShearVelocity
   !=====================================================

   FUNCTION calc_shearVelocity(rivvel,rivdph) result(sVel)
   USE YOS_CMF_INPUT,           only:  PGRV, PMANRIV
   IMPLICIT NONE
   SAVE
   real(kind=JPRB)                 ::  sVel ! shear velocity[m/s]
   real(kind=JPRB), intent(in)     ::  rivvel, rivdph 
   !========

      sVel = sqrt ( PGRV * PMANRIV**2.d0 * rivvel**2.d0 * rivdph**(-1.d0/3.d0) )  !bug fix 2022/11/22
      RETURN
   END FUNCTION calc_shearVelocity
   !=====================================================

   FUNCTION calc_suspendVelocity(csVel,sVel,setVel) result(susVel) ! Uchida and Fukuoka (2019) Eq.44
   USE yos_cmf_sed,             only: lambda, nsed, vonKar
   IMPLICIT NONE
   SAVE
   real(kind=JPRB)                 ::  susVel(nsed)   ! suspend velocity [m/s]
   !!  real(kind=JPRB), intent(in)     ::  csVel(nsed), sVel, setVel(nsed)    ! crit shear, shear velocity, setting velocity [m/s]
   real(kind=JPRB), intent(in)     ::  csVel(:), sVel, setVel(:)    ! crit shear, shear velocity, setting velocity [m/s]
   integer(kind=JPIM)              ::  ised
   real(kind=JPRB)                 ::  alpha, a, cB, sTmp
   !========
      !--------------!
      ! set constant !
      !--------------!
      alpha = vonKar / 6.d0
      a = 0.08d0
      cB = 1.d0 - lambda
  
      !-----------------------!
      ! calc suspend velocity !
      !-----------------------!
      susVel(:) = 0.d0
      DO ised = 1, nsed
         IF ( csVel(ised) > sVel ) CYCLE
         sTmp = setVel(ised) / alpha / sVel
         susVel(ised) = max( setVel(ised) * cB / (1.d0+sTmp) * (1.d0-a*sTmp) / (1.d0+(1.d0-a)*sTmp), 0.d0 )
      ENDDO
   END FUNCTION calc_suspendVelocity

!####################################################################

END MODULE cmf_calc_sedpar_mod
