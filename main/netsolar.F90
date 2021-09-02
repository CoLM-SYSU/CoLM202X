#include <define.h>

 SUBROUTINE netsolar (ipatch,idate,deltim,dlon,patchtype,&
                      forc_sols,forc_soll,forc_solsd,forc_solld,&
                      alb,ssun,ssha,lai,sai,rho,tau,&
                      parsun,parsha,sabvsun,sabvsha,sabg,sabvg,sr,&
                      solvd,solvi,solnd,solni,srvd,srvi,srnd,srni,&
                      solvdln,solviln,solndln,solniln,srvdln,srviln,srndln,srniln)
!=======================================================================
! Net solar absorbed by surface
! Original author : Yongjiu Dai, 09/15/1999; 09/11/2001
!=======================================================================

  USE precision
  USE GlobalVars
  USE timemanager, only: isgreenwich
  USE MOD_PFTimeInvars
  USE MOD_PFTimeVars
  USE MOD_1D_PFTFluxes
  USE MOD_PCTimeInvars
  USE MOD_PCTimeVars
  USE MOD_1D_PCFluxes

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
        sabvg,    &! solar absorbed by ground + vegetation [W/m2]
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

! ----------------local variables ---------------------------------
   INTEGER  :: local_secs
   REAL(r8) :: radpsec

   INTEGER ps, pe, pc
!=======================================================================
       
        sabvsun = 0.
        sabvsha = 0.
        parsun  = 0.
        parsha  = 0.

        sabg  = 0.
        sabvg = 0.

IF (patchtype == 0) THEN

#ifdef PFT_CLASSIFICATION
        ps = patch_pft_s(ipatch)      
        pe = patch_pft_e(ipatch)
        sabvsun_p(ps:pe) = 0.
        sabvsha_p(ps:pe) = 0.
        parsun_p(ps:pe)  = 0.
        parsha_p(ps:pe)  = 0.
#endif

#ifdef PC_CLASSIFICATION
        pc = patch2pc(ipatch)
        sabvsun_c(:,pc) = 0.
        sabvsha_c(:,pc) = 0.
        parsun_c(:,pc)  = 0.
        parsha_c(:,pc)  = 0.
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
              sabvg   = forc_sols *(1.-alb(1,1)) + forc_soll *(1.-alb(2,1)) &
                      + forc_solsd*(1.-alb(1,2)) + forc_solld*(1.-alb(2,2))
              sabg    = sabvg - sabvsun - sabvsha

              !TODO: bug exist
              ! 08/17/2021, yuan: LAI PAR, 区别lai和sai的吸收
              !parsun  = parsun * lai*(1.-rho(1,1)-tau(1,1)) / &
              !   ( lai*(1.-rho(1,1)-tau(1,1)) + sai*(1.-rho(1,2)-tau(1,2)) )

              !parsha  = parsha * lai*(1.-rho(1,1)-tau(1,1)) / &
              !   ( lai*(1.-rho(1,1)-tau(1,1)) + sai*(1.-rho(1,2)-tau(1,2)) )

IF (patchtype == 0) THEN

#ifdef PFT_CLASSIFICATION
              parsun_p(ps:pe)  = forc_sols*ssun_p(1,1,ps:pe) + forc_solsd*ssun_p(1,2,ps:pe)
              parsha_p(ps:pe)  = forc_sols*ssha_p(1,1,ps:pe) + forc_solsd*ssha_p(1,2,ps:pe)
              sabvsun_p(ps:pe) = forc_sols*ssun_p(1,1,ps:pe) + forc_solsd*ssun_p(1,2,ps:pe) &
                               + forc_soll*ssun_p(2,1,ps:pe) + forc_solld*ssun_p(2,2,ps:pe)
              sabvsha_p(ps:pe) = forc_sols*ssha_p(1,1,ps:pe) + forc_solsd*ssha_p(1,2,ps:pe) &
                               + forc_soll*ssha_p(2,1,ps:pe) + forc_solld*ssha_p(2,2,ps:pe)
#endif

#ifdef PC_CLASSIFICATION
              parsun_c(:,pc)  = forc_sols*ssun_c(1,1,:,pc) + forc_solsd*ssun_c(1,2,:,pc)
              parsha_c(:,pc)  = forc_sols*ssha_c(1,1,:,pc) + forc_solsd*ssha_c(1,2,:,pc)
              sabvsun_c(:,pc) = forc_sols*ssun_c(1,1,:,pc) + forc_solsd*ssun_c(1,2,:,pc) &
                              + forc_soll*ssun_c(2,1,:,pc) + forc_solld*ssun_c(2,2,:,pc) 
              sabvsha_c(:,pc) = forc_sols*ssha_c(1,1,:,pc) + forc_solsd*ssha_c(1,2,:,pc) &
                              + forc_soll*ssha_c(2,1,:,pc) + forc_solld*ssha_c(2,2,:,pc) 
#endif

ENDIF
           ELSE               !lake or ocean
              sabvg = forc_sols *(1.-alb(1,1)) + forc_soll *(1.-alb(2,1)) &
                    + forc_solsd*(1.-alb(1,2)) + forc_solld*(1.-alb(2,2))
              sabg = sabvg
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

        !print *, "solar radiation balance check:", forc_sols+forc_soll+forc_solsd+forc_solld-&
        !   sabg-sabvsun-sabvsha-sr

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
