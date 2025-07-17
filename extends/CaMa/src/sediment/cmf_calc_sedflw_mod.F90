MODULE cmf_calc_sedflw_mod
!==========================================================
!* PURPOSE: physics for sediment transport
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
CONTAINS
!####################################################################
! -- CMF_CALC_SEDFLW
! --
! --
!####################################################################
   SUBROUTINE cmf_calc_sedflw
   USE PARKIND1,                only: JPIM, JPRB
   USE YOS_CMF_INPUT,           only: PGRV
   USE YOS_CMF_MAP,             only: D2RIVLEN, D2RIVWTH, NSEQALL
   USE YOS_CMF_PROG,            only: P2RIVSTO
   USE YOS_CMF_DIAG,            only: D2RIVDPH
   USE yos_cmf_sed,             only: lambda, nsed, sedDT, setVel, &
                                       d2layer, d2sedcon, d2rivsto_pre
   USE sed_utils_mod,           only: sed_diag_average, sed_diag_reset

   IMPLICIT NONE
   !$ SAVE
   SAVE
   integer(kind=JPIM)              :: ISEQ
   real(kind=JPRB)                 :: sedsto(NSEQALL,nsed)
   real(kind=JPRB)                 :: shearVel(NSEQALL)
   real(kind=JPRB)                 :: critShearVel(NSEQALL,nsed), dMean(NSEQALL), susVel(NSEQALL,nsed)
   real(kind=JPRB), parameter      :: IGNORE_DPH = 0.05d0
   !================================================

      CALL sed_diag_average 
  !$omp parallel do
      DO iseq = 1, NSEQALL
         sedsto(iseq,:) = d2sedcon(iseq,:) * max(d2rivsto_pre(iseq), 0.d0)
      ENDDO
  !$omp END parallel DO

      CALL calc_params
      CALL calc_advection
      CALL calc_entrainment
      CALL calc_exchange

  !$omp parallel DO
      DO iseq = 1, NSEQALL
         IF ( P2RIVSTO(iseq,1) < D2RIVWTH(iseq,1)*D2RIVLEN(iseq,1)*IGNORE_DPH ) CYCLE
         d2sedcon(iseq,:) = sedsto(iseq,:) / P2RIVSTO(iseq,1)
      ENDDO
  !$omp END parallel DO

      CALL sed_diag_reset

   CONTAINS
      !==========================================================
      !+ calc_params
      !+ calc_advection
      !+ calc_entrainment
      !+ calc_exchange
      !==========================================================
      SUBROUTINE calc_params
      USE yos_cmf_sed,           only: pset, revEgia, sDiam, visKin, d2rivvel_sed
      USE cmf_calc_sedpar_mod,   only: calc_criticalShearVelocity, calc_shearVelocity, calc_suspendVelocity
      IMPLICIT NONE
      SAVE
      integer(kind=JPIM)            ::  ised, iseq
      real(kind=JPRB)               ::  csVel0, sTmp, sTmp1(nsed)
      !=====================================================
      
         DO iseq = 1, NSEQALL
            !-------------------------!
            ! critical shear velocity !
            !-------------------------!
            
            IF ( sum(d2layer(iseq,:)) <= 0.d0 ) THEN
               critShearVel(iseq,:) = 1e20
            ELSE IF ( revEgia ) THEN
               dMean(iseq) = 0.d0
               DO ised = 1, nsed
                  dMean(iseq) = dMean(iseq) + sDiam(ised)*d2layer(iseq,ised)/sum(d2layer(iseq,:))
               ENDDO
               csVel0 = calc_criticalShearVelocity(dMean(iseq))
               DO ised = 1, nsed
                  IF ( sDiam(ised) / dMean(iseq) >= 0.4d0 ) THEN
                     critShearVel(iseq,ised) = sqrt( csVel0*sDiam(ised)/dMean(iseq) ) * &
                     & ( log10(19.d0)/log10(19.d0*sDiam(ised)/dMean(iseq)) ) * 0.01d0
                  ELSE
                     critShearVel(iseq,ised) = sqrt( 0.85*csVel0 ) * 0.01d0
                  ENDIF
               ENDDO      
            ELSE
               DO ised = 1, nsed
                  critShearVel(iseq,ised) = sqrt( calc_criticalShearVelocity(sDiam(ised)) ) * 0.01d0
               ENDDO
            ENDIF
         
            !------------------------------------------------------!
            ! shear velocity, suspend velocity, Karman coefficient !
            !------------------------------------------------------!
            IF ( d2rivvel_sed(iseq) == 0.d0 .or. D2RIVDPH(iseq,1) < IGNORE_DPH ) THEN
               shearVel(iseq) = 0.d0
               susVel(iseq,:) = 0.d0
            ELSE
               shearVel(iseq) = calc_shearVelocity(d2rivvel_sed(iseq), D2RIVDPH(iseq,1))
               susVel(iseq,:) = calc_suspendVelocity(critShearVel(iseq,:), shearVel(iseq), setVel(:))
            ENDIF
         ENDDO
      END SUBROUTINE calc_params
   !=====================================================

      SUBROUTINE calc_advection
      USE YOS_CMF_MAP,           only:  I1NEXT
      USE yos_cmf_sed,           only:  d2rivout_sed, d2bedout, d2sedout, &
                                       d2bedout_avg, d2sedout_avg, psedD, pwatD
      IMPLICIT NONE
      real(kind=JPRB)               ::  bOut(NSEQALL,nsed), brate(NSEQALL,nsed)
      real(kind=JPRB)               ::  sOut(NSEQALL,nsed), srate(NSEQALL,nsed)
      integer(kind=JPIM)            ::  ised, iseq
      ! SAVE for omop
      real(kind=JPRB), SAVE         ::  plusVel, minusVel
      integer(kind=JPIM), SAVE      ::  iseq0, iseq1
!$omp threadprivate ( plusVel, minusVel, iseq0, iseq1 )
         !========
         bOut(:,:) = 0.d0
         sOut(:,:) = 0.d0
!$omp parallel DO
         DO iseq = 1, NSEQALL
            
            IF ( d2rivout_sed(iseq) >= 0.d0 ) THEN
               iseq0 = iseq
               iseq1 = I1NEXT(iseq)
            ELSE
               iseq0 = I1NEXT(iseq)
               iseq1 = iseq
            ENDIF

            IF ( d2rivout_sed(iseq) == 0.d0 ) THEN
               d2sedout(iseq,:) = 0.d0
               d2bedout(iseq,:) = 0.d0
               CYCLE
            ENDIF

            !-------------------!
            ! calc suspend flow !
            !-------------------!
            IF ( iseq0 < 0 ) THEN
               d2sedout(iseq,:) = d2sedcon(iseq1,:) * d2rivout_sed(iseq)
            ELSE
               d2sedout(iseq,:) = d2sedcon(iseq0,:) * d2rivout_sed(iseq)
               sOut(iseq0,:) = sOut(iseq0,:) + abs(d2sedout(iseq,:))*sedDT
            ENDIF

            !--------------!
            ! calc bedflow !
            !--------------!
            IF ( minval(critShearVel(iseq,:)) >= shearVel(iseq) .or. sum(d2layer(iseq,:)) == 0.d0 .or. iseq0 < 0  ) THEN
               d2bedout(iseq,:) = 0.d0
            ELSE 
               DO ised = 1, nsed
                  IF ( critShearVel(iseq,ised) >= shearVel(iseq) .or. d2layer(iseq,ised) == 0.d0 ) THEN
                     d2bedout(iseq,ised) = 0.d0
                     CYCLE
                  ENDIF
                  plusVel = shearVel(iseq) + critShearVel(iseq,ised)
                  minusVel = shearVel(iseq) - critShearVel(iseq,ised)
                  d2bedout(iseq,ised) = 17.d0 * D2RIVWTH(iseq,1) * plusVel * minusVel * minusVel & 
                  & / ((psedD-pwatD)/pwatD) / PGRV * d2layer(iseq,ised) / sum(d2layer(iseq,:)) 
                  bOut(iseq0,ised) = bOut(iseq0,ised) + d2bedout(iseq,ised)*sedDT
               ENDDO
            ENDIF
         ENDDO
!$omp END parallel DO

         !--------------------------------------------!
         ! adjust outflow IF larget than sedsto/layer !
         !--------------------------------------------!
         brate(:,:) = 1.d0
         srate(:,:) = 1.d0
!$omp parallel DO
         DO iseq = 1, NSEQALL
            IF ( minval(sOut(iseq,:)) <= 1e-8 ) THEN
               DO ised = 1, nsed
                  IF ( sOut(iseq,ised) > 1e-8 ) THEN
                     srate(iseq,ised) = min ( sedsto(iseq,ised) / sOut(iseq,ised), 1.d0 )
                  ENDIF
               ENDDO
            ELSE
               srate(iseq,:) = min ( sedsto(iseq,:) / sOut(iseq,:), 1.d0 )
            ENDIF
            IF ( minval(bOut(iseq,:)) <= 1e-8 ) THEN
               DO ised = 1, nsed
                  IF ( bOut(iseq,ised) > 1e-8 ) THEN
                     brate(iseq,ised) = min( d2layer(iseq,ised) / bOut(iseq,ised), 1.d0 )
                  ENDIF
               ENDDO
            ELSE
               brate(iseq,:) = min( d2layer(iseq,:) / bOut(iseq,:), 1.d0 )
            ENDIF
      ENDDO
!$omp END parallel DO
    
         DO iseq = 1, NSEQALL
            IF ( d2rivout_sed(iseq) >= 0.d0 ) THEN
               iseq0 = iseq
               iseq1 = I1NEXT(iseq)
            ELSE
               iseq0 = I1NEXT(iseq)
               iseq1 = iseq
            ENDIF

            IF ( iseq0 > 0 ) THEN
               d2sedout(iseq,:) = d2sedout(iseq,:) * srate(iseq0,:)
               sedsto(iseq0,:) = max( sedsto(iseq0,:)-abs(d2sedout(iseq,:))*sedDT, 0.d0 )
               d2bedout(iseq,:) = d2bedout(iseq,:) * brate(iseq0,:)
               d2layer(iseq0,:) = max( d2layer(iseq0,:)-abs(d2bedout(iseq,:))*sedDT, 0.d0 )
            ENDIF
            IF ( iseq1 > 0 ) THEN
               sedsto(iseq1,:) = max( sedsto(iseq1,:)+abs(d2sedout(iseq,:))*sedDT, 0.d0 )
               d2layer(iseq1,:) = max( d2layer(iseq1,:)+abs(d2bedout(iseq,:))*sedDT, 0.d0 ) 
            ENDIF

            d2bedout_avg(iseq,:) = d2bedout_avg(iseq,:) + d2bedout(iseq,:)*sedDT
            d2sedout_avg(iseq,:) = d2sedout_avg(iseq,:) + d2sedout(iseq,:)*sedDT
         ENDDO

      END SUBROUTINE calc_advection
  !=====================================================

      SUBROUTINE calc_entrainment
      USE yos_cmf_sed,           only:  vonKar, d2netflw, d2netflw_avg, d2sedinp, d2seddep, totlyrnum
      
      IMPLICIT NONE
      real(kind=JPRB)               ::  dTmp(NSEQALL,nsed), D(NSEQALL,nsed), Es(NSEQALL,nsed), Zd(NSEQALL,nsed)
      integer(kind=JPIM)            ::  ilyr, ised, iseq
      real(kind=JPRB),SAVE          ::  dTmp1, layerP
!$omp threadprivate  ( dTmp1, layerP )
!========

!$omp parallel DO
         DO iseq = 1, NSEQALL
            IF ( D2RIVDPH(iseq,1) < IGNORE_DPH ) THEN
               d2netflw(iseq,:) = 0.d0
               CYCLE
            ENDIF

            !----------------------!
            ! calculate suspension !
            !----------------------!
            IF ( sum(d2layer(iseq,:)) == 0.d0 .or. all(susVel(iseq,:)==0.d0) ) THEN
               Es(iseq,:) = 0.d0
            ELSE
               Es(iseq,:) = susVel(iseq,:) * (1.d0-lambda) * D2RIVWTH(iseq,1) * D2RIVLEN(iseq,1) * d2layer(iseq,:) / sum(d2layer(iseq,:))
               Es(iseq,:) = max( Es(iseq,:), 0.d0 )
            ENDIF

            !----------------------!
            ! calculate deposition !
            !----------------------!
            IF ( shearVel(iseq) == 0.d0 .or. all(setVel(:)==0.d0) ) THEN
               D(iseq,:) = 0.d0
            ELSE
               Zd(iseq,:) = 6.d0 * setVel(:) / vonKar / shearVel(iseq)
               D(iseq,:) = setVel(:) * D2RIVWTH(iseq,1) * D2RIVLEN(iseq,1) * d2sedcon(iseq,:) * Zd(iseq,:) / (1.d0-exp(-Zd(iseq,:))) 
               D(iseq,:) = max( D(iseq,:), 0.d0 )
            ENDIF
            d2netflw(iseq,:) = Es(iseq,:) - D(iseq,:)
      
            !-------------------------------------------!
            ! IF >0, suspension ; IF <0, deposition     !
            ! adjust netflw IF larger than sedsto/layer !
            !-------------------------------------------!
            DO ised = 1, nsed
               IF ( d2netflw(iseq,ised) == 0.d0 ) THEN
                  CYCLE
               ELSE IF ( d2netflw(iseq,ised) > 0.d0 ) THEN
                  dTmp1 = d2netflw(iseq,ised)*sedDT/(1.d0-lambda)
                  IF ( dTmp1 < d2layer(iseq,ised) ) THEN
                     d2layer(iseq,ised) = d2layer(iseq,ised) - dTmp1
                  ELSE
                     d2netflw(iseq,ised) = d2layer(iseq,ised) * (1.d0-lambda) / sedDT
                     d2layer(iseq,ised) = 0.d0
                  ENDIF
                  sedsto(iseq,ised) = sedsto(iseq,ised) + d2netflw(iseq,ised) * sedDT
               ELSE
                  IF ( abs(d2netflw(iseq,ised))*sedDT < sedsto(iseq,ised) ) THEN
                     sedsto(iseq,ised) = max (sedsto(iseq,ised) - abs(d2netflw(iseq,ised))*sedDT, 0.d0 )
                  ELSE
                     d2netflw(iseq,ised) = - sedsto(iseq,ised) / sedDT
                     sedsto(iseq,ised) = 0.d0
                  ENDIF
                  d2layer(iseq,ised) = d2layer(iseq,ised) + abs(d2netflw(iseq,ised))*sedDT/(1.d0-lambda)
               ENDIF
            ENDDO

            sedsto(iseq,:) = sedsto(iseq,:) + d2sedinp(iseq,:)*sedDT    
            IF ( sum(sedsto(iseq,:)) > P2RIVSTO(iseq,1) * 0.01d0 ) THEN
               dTmp(iseq,:) = ( sum(sedsto(iseq,:)) - P2RIVSTO(iseq,1)*0.01d0 ) * sedsto(iseq,:)/sum(sedsto(iseq,:))
               d2netflw(iseq,:) = d2netflw(iseq,:) - dTmp(iseq,:)/sedDT
               sedsto(iseq,:) = sedsto(iseq,:) - dTmp(iseq,:)
               d2layer(iseq,:) = d2layer(iseq,:) + dTmp(iseq,:)/(1.d0-lambda)
            ENDIF
            
            d2netflw_avg(iseq,:) = d2netflw_avg(iseq,:) + d2netflw(iseq,:)*sedDT
         ENDDO
!$omp END parallel DO
      END SUBROUTINE calc_entrainment
  !=====================================================

      SUBROUTINE calc_exchange
      ! redistribute into vertical bed layers
      USE yos_cmf_sed,           only:  d2seddep, lyrdph, totlyrnum
      
      IMPLICIT NONE
      integer(kind=JPIM)            ::  ilyr, ised, iseq, jlyr, slyr
      real(kind=JPRB)               ::  diff, lyrvol, layerP(nsed), seddepP(totlyrnum+1,nsed), tmp(nsed)

         DO iseq = 1, NSEQALL
            lyrvol = lyrdph * D2RIVWTH(iseq,1) * D2RIVLEN(iseq,1)

            IF ( minval(d2layer(iseq,:)) < 0.d0 ) d2layer(iseq,:) = max( d2layer(iseq,:), 0.d0 )
            IF ( minval(d2seddep(iseq,:,:)) < 0.d0 ) d2seddep(iseq,:,:) = max( d2seddep(iseq,:,:), 0.d0 )

            !---------------------------------------!
            ! IF bed storage less than layer volume !
            !---------------------------------------!
            IF ( sum(d2layer(iseq,:)) + sum(d2seddep(iseq,:,:)) <= lyrvol ) THEN
               d2layer(iseq,:) = d2layer(iseq,:) + sum(d2seddep(iseq,:,:),dim=1)
               d2seddep(iseq,:,:) = 0.d0
               CYCLE
            ENDIF

            !------------------------------------!
            ! distribute into top exchange layer !
            !------------------------------------!
            layerP(:) = d2layer(iseq,:)
            IF ( sum(layerP(:)) >= lyrvol ) THEN
               d2layer(iseq,:) = layerP(:) * min( lyrvol/sum(layerP(:)), 1.d0 )
               layerP(:) = max( layerP(:) - d2layer(iseq,:), 0.d0 )
               slyr = 0
            ELSE IF ( sum(d2seddep(iseq,:,:)) > 0.d0 ) THEN
               layerP(:) = 0.d0
               DO ilyr = 1, totlyrnum
                  diff = lyrvol - sum(d2layer(iseq,:))
                  IF ( diff <= 0.d0 ) EXIT
                  IF ( sum(d2seddep(iseq,ilyr,:)) <= diff ) THEN
                     d2layer(iseq,:) = d2layer(iseq,:) + d2seddep(iseq,ilyr,:)
                     d2seddep(iseq,ilyr,:) = 0.d0
                     slyr = ilyr + 1
                  ELSE
                     tmp(:) = diff * d2seddep(iseq,ilyr,:) / sum(d2seddep(iseq,ilyr,:))
                     d2layer(iseq,:) = d2layer(iseq,:) + tmp(:)
                     d2seddep(iseq,ilyr,:) = max( d2seddep(iseq,ilyr,:) - tmp(:), 0.d0 )
                     slyr = ilyr
                     EXIT
                  ENDIF
               ENDDO
            ELSE
               d2seddep(iseq,:,:) = 0.d0
               CYCLE
            ENDIF
            IF ( sum(d2seddep(iseq,:,:)) == 0.d0 ) CYCLE

            !-----------------------------------!
            ! distribute remaining bedload into !
            ! vertical deposition layers        !
            !-----------------------------------!
            seddepP(1,:) = layerP(:)
            seddepP(2:,:) = d2seddep(iseq,:,:)
            d2seddep(iseq,:,:) = 0.d0
            DO ilyr = 1, totlyrnum - 1
               IF ( sum(d2seddep(iseq,ilyr,:)) == lyrvol ) CYCLE
               DO jlyr = slyr+1, totlyrnum + 1
                  diff = lyrvol - sum(d2seddep(iseq,ilyr,:))
                  IF ( diff <= 0.d0 ) EXIT
                  IF ( sum(seddepP(jlyr,:)) <= diff ) THEN
                     d2seddep(iseq,ilyr,:) = d2seddep(iseq,ilyr,:) + seddepP(jlyr,:)
                     seddepP(jlyr,:) = 0.d0
                  ELSE
                     tmp(:) = diff * seddepP(jlyr,:) / sum(seddepP(jlyr,:))
                     d2seddep(iseq,ilyr,:) = d2seddep(iseq,ilyr,:) + tmp(:)
                     seddepP(jlyr,:) = max(seddepP(jlyr,:) - tmp(:), 0.d0)
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
               
               IF ( sum(seddepP) > 0.d0 ) THEN
                  d2seddep(iseq,totlyrnum,:) = sum(seddepP, dim=1)
               ENDIF
         ENDDO

      END SUBROUTINE calc_exchange
   END SUBROUTINE cmf_calc_sedflw
   !####################################################################

END MODULE cmf_calc_sedflw_mod
