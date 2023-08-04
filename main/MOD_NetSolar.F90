#include <define.h>

MODULE MOD_NetSolar

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: netsolar


!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------


   SUBROUTINE netsolar (ipatch,idate,deltim,dlon,patchtype,&
                        forc_sols,forc_soll,forc_solsd,forc_solld,&
                        alb,ssun,ssha,lai,sai,rho,tau,ssno,&
                        parsun,parsha,sabvsun,sabvsha,sabg,sabg_lyr,sr,&
                        solvd,solvi,solnd,solni,srvd,srvi,srnd,srni,&
                        solvdln,solviln,solndln,solniln,srvdln,srviln,srndln,srniln)
!
! !DESCRIPTION:
! Net solar absorbed by surface
!
! Original author : Yongjiu Dai, 09/15/1999; 09/11/2001
!
! REVISIONS:
! Hua Yuan, 05/2014: added for solar radiation output [vars: so*, sr*]
!
! Hua Yuan, 08/2014: added for local noon calculation
!
! Hua Yuan, 08/2020: added for PFT and PC calculation
!
! Hua Yuan, 12/2022: calculated snow layer absorption by SNICAR model
!
! !USES:
   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_Namelist, only: DEF_USE_SNICAR
   USE MOD_TimeManager, only: isgreenwich
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT, only: patch_pft_s, patch_pft_e
   USE MOD_Vars_PFTimeInvariants
   USE MOD_Vars_PFTimeVariables
   USE MOD_Vars_1DPFTFluxes
#endif

   IMPLICIT NONE

! Dummy argument
   INTEGER,  intent(in) :: ipatch     !patch index
   INTEGER,  intent(in) :: idate(3)   !model time
   INTEGER,  intent(in) :: patchtype  !land water TYPE (99-sea)

   REAL(r8), intent(in) :: dlon       !logitude in radians
   REAL(r8), intent(in) :: deltim     !seconds in a time step [second]

   REAL(r8), intent(in) :: &
         forc_sols,  &! atm vis direct beam solar rad onto srf [W/m2]
         forc_soll,  &! atm nir direct beam solar rad onto srf [W/m2]
         forc_solsd, &! atm vis diffuse solar rad onto srf [W/m2]
         forc_solld   ! atm nir diffuse solar rad onto srf [W/m2]

   REAL(r8), dimension(1:2,1:2), intent(in) :: &
         alb,      &! averaged albedo [-]
         ssun,     &! sunlit canopy absorption for solar radiation
         ssha       ! shaded canopy absorption for solar radiation

   REAL(r8), dimension(1:2,1:2,maxsnl+1:1), intent(inout) :: &
         ssno       ! snow layer absorption

   REAL(r8), intent(in) :: &
         lai,      &! leaf area index
         sai,      &! stem area index
         rho(2,2), &! leaf reflectance (iw=iband, il=life and dead)
         tau(2,2)   ! leaf transmittance (iw=iband, il=life and dead)

   REAL(r8), intent(out) :: &
         parsun,   &! PAR absorbed by sunlit vegetation [W/m2]
         parsha,   &! PAR absorbed by shaded vegetation [W/m2]
         sabvsun,  &! solar absorbed by sunlit vegetation [W/m2]
         sabvsha,  &! solar absorbed by shaded vegetation [W/m2]
         sabg,     &! solar absorbed by ground  [W/m2]
         sr,       &! total reflected solar radiation (W/m2)
         solvd,    &! incident direct beam vis solar radiation (W/m2)
         solvi,    &! incident diffuse beam vis solar radiation (W/m2)
         solnd,    &! incident direct beam nir solar radiation (W/m2)
         solni,    &! incident diffuse beam nir solar radiation (W/m2)
         srvd,     &! reflected direct beam vis solar radiation (W/m2)
         srvi,     &! reflected diffuse beam vis solar radiation (W/m2)
         srnd,     &! reflected direct beam nir solar radiation (W/m2)
         srni,     &! reflected diffuse beam nir solar radiation (W/m2)
         solvdln,  &! incident direct beam vis solar radiation at local noon(W/m2)
         solviln,  &! incident diffuse beam vis solar radiation at local noon(W/m2)
         solndln,  &! incident direct beam nir solar radiation at local noon(W/m2)
         solniln,  &! incident diffuse beam nir solar radiation at local noon(W/m2)
         srvdln,   &! reflected direct beam vis solar radiation at local noon(W/m2)
         srviln,   &! reflected diffuse beam vis solar radiation at local noon(W/m2)
         srndln,   &! reflected direct beam nir solar radiation at local noon(W/m2)
         srniln     ! reflected diffuse beam nir solar radiation at local noon(W/m2)

   REAL(r8), intent(out) :: &
         sabg_lyr(maxsnl+1:1)   ! solar absorbed by snow layers [W/m2]

! ----------------local variables ---------------------------------
   INTEGER  :: local_secs
   REAL(r8) :: radpsec, sabvg

   INTEGER ps, pe, pc
!=======================================================================

      sabvsun = 0.
      sabvsha = 0.
      parsun  = 0.
      parsha  = 0.

      sabg  = 0.
      sabg_lyr(:) = 0.

      IF (patchtype == 0) THEN
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         ps = patch_pft_s(ipatch)
         pe = patch_pft_e(ipatch)
         sabvsun_p(ps:pe) = 0.
         sabvsha_p(ps:pe) = 0.
         parsun_p(ps:pe)  = 0.
         parsha_p(ps:pe)  = 0.
#endif
      ENDIF


      IF (forc_sols+forc_soll+forc_solsd+forc_solld > 0.) THEN
         IF (patchtype < 4) THEN        !non lake and ocean
          ! Radiative fluxes onto surface
            parsun  = forc_sols*ssun(1,1) + forc_solsd*ssun(1,2)
            parsha  = forc_sols*ssha(1,1) + forc_solsd*ssha(1,2)
            sabvsun = forc_sols*ssun(1,1) + forc_solsd*ssun(1,2) &
                    + forc_soll*ssun(2,1) + forc_solld*ssun(2,2)
            sabvsha = forc_sols*ssha(1,1) + forc_solsd*ssha(1,2) &
                    + forc_soll*ssha(2,1) + forc_solld*ssha(2,2)
            sabvg   = forc_sols *(1.-alb(1,1)) + forc_solsd*(1.-alb(1,2)) &
                    + forc_soll *(1.-alb(2,1)) + forc_solld*(1.-alb(2,2))
            sabg    = sabvg - sabvsun - sabvsha

            IF (patchtype == 0) THEN
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
               parsun_p(ps:pe)  = forc_sols*ssun_p(1,1,ps:pe) + forc_solsd*ssun_p(1,2,ps:pe)
               parsha_p(ps:pe)  = forc_sols*ssha_p(1,1,ps:pe) + forc_solsd*ssha_p(1,2,ps:pe)
               sabvsun_p(ps:pe) = forc_sols*ssun_p(1,1,ps:pe) + forc_solsd*ssun_p(1,2,ps:pe) &
                                + forc_soll*ssun_p(2,1,ps:pe) + forc_solld*ssun_p(2,2,ps:pe)
               sabvsha_p(ps:pe) = forc_sols*ssha_p(1,1,ps:pe) + forc_solsd*ssha_p(1,2,ps:pe) &
                                + forc_soll*ssha_p(2,1,ps:pe) + forc_solld*ssha_p(2,2,ps:pe)
#endif
            ENDIF

         ELSE               !lake or ocean
            sabvg = forc_sols *(1.-alb(1,1)) + forc_soll *(1.-alb(2,1)) &
                  + forc_solsd*(1.-alb(1,2)) + forc_solld*(1.-alb(2,2))
            sabg = sabvg
         ENDIF

         IF (DEF_USE_SNICAR) THEN

            IF (patchtype < 4) THEN        !non lake and ocean
               ! normalization
               IF(sum(ssno(1,1,:))>0.) ssno(1,1,:) = (1-alb(1,1)-ssun(1,1)-ssha(1,1)) * ssno(1,1,:)/sum(ssno(1,1,:))
               IF(sum(ssno(1,2,:))>0.) ssno(1,2,:) = (1-alb(1,2)-ssun(1,2)-ssha(1,2)) * ssno(1,2,:)/sum(ssno(1,2,:))
               IF(sum(ssno(2,1,:))>0.) ssno(2,1,:) = (1-alb(2,1)-ssun(2,1)-ssha(2,1)) * ssno(2,1,:)/sum(ssno(2,1,:))
               IF(sum(ssno(2,2,:))>0.) ssno(2,2,:) = (1-alb(2,2)-ssun(2,2)-ssha(2,2)) * ssno(2,2,:)/sum(ssno(2,2,:))
            ELSE                           !lake case
               ! normalization
               IF(sum(ssno(1,1,:))>0.) ssno(1,1,:) = (1-alb(1,1)) * ssno(1,1,:)/sum(ssno(1,1,:))
               IF(sum(ssno(1,2,:))>0.) ssno(1,2,:) = (1-alb(1,2)) * ssno(1,2,:)/sum(ssno(1,2,:))
               IF(sum(ssno(2,1,:))>0.) ssno(2,1,:) = (1-alb(2,1)) * ssno(2,1,:)/sum(ssno(2,1,:))
               IF(sum(ssno(2,2,:))>0.) ssno(2,2,:) = (1-alb(2,2)) * ssno(2,2,:)/sum(ssno(2,2,:))
            ENDIF

            ! snow layer absorption
            sabg_lyr(:) = forc_sols*ssno(1,1,:) + forc_solsd*ssno(1,2,:) &
                        + forc_soll*ssno(2,1,:) + forc_solld*ssno(2,2,:)
         ENDIF
      ENDIF

      solvd = forc_sols
      solvi = forc_solsd
      solnd = forc_soll
      solni = forc_solld
      srvd  = solvd*alb(1,1)
      srvi  = solvi*alb(1,2)
      srnd  = solnd*alb(2,1)
      srni  = solni*alb(2,2)
      sr    = srvd + srvi + srnd + srni

      ! calculate the local secs
      radpsec = pi/12./3600.
      IF ( isgreenwich ) THEN
         local_secs = idate(3) + nint((dlon/radpsec)/deltim)*deltim
         local_secs = mod(local_secs,86400)
      ELSE
         local_secs = idate(3)
      ENDIF

      IF (local_secs == 86400/2) THEN
         solvdln = forc_sols
         solviln = forc_solsd
         solndln = forc_soll
         solniln = forc_solld
         srvdln  = solvdln*alb(1,1)
         srviln  = solviln*alb(1,2)
         srndln  = solndln*alb(2,1)
         srniln  = solniln*alb(2,2)
      ELSE
         solvdln = spval
         solviln = spval
         solndln = spval
         solniln = spval
         srvdln  = spval
         srviln  = spval
         srndln  = spval
         srniln  = spval
      ENDIF

   END SUBROUTINE netsolar

END MODULE MOD_NetSolar
