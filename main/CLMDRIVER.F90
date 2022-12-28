#include <define.h>

SUBROUTINE CLMDRIVER (idate,deltim,dolai,doalb,dosst,oro)


!=======================================================================
!
! CLM MODEL DRIVER
!
! Original author : Yongjiu Dai, 09/30/1999; 08/30/2002, 03/2014
!
!=======================================================================

 use precision
 use PhysicalConstants, only: tfrz, rgas, vonkar
 USE GlobalVars
 USE LC_Const
 use MOD_TimeInvariants
 use MOD_TimeVariables
 use MOD_1D_Forcing
 use MOD_1D_Fluxes
 USE mod_landpatch, only : numpatch
 USE mod_namelist, only : DEF_forcing
 USE mod_forcing, only : forcmask
 use omp_lib
#ifdef CaMa_Flood
 use MOD_CaMa_Variables, only : flddepth_cama,fldfrc_cama,fevpg_fld,finfg_fld
#endif
 IMPLICIT NONE

  integer,  INTENT(in) :: idate(3) ! model calendar for next time step (year, julian day, seconds)
  real(r8), INTENT(in) :: deltim   ! seconds in a time-step

  logical,  INTENT(in) :: dolai    ! true if time for time-varying vegetation paramter
  logical,  INTENT(in) :: doalb    ! true if time for surface albedo calculation
  logical,  INTENT(in) :: dosst    ! true if time for update sst/ice/snow

  real(r8), INTENT(inout) :: oro(numpatch)  ! ocean(0)/seaice(2)/ flag

! -------------- Local varaibles -------------------
  real(r8), allocatable :: z_soisno (:,:)
  real(r8), allocatable :: dz_soisno(:,:)

  integer :: i, m

! ======================================================================

  !TODO: can be removed below
  allocate ( z_soisno  (maxsnl+1:nl_soil,numpatch) )
  allocate ( dz_soisno (maxsnl+1:nl_soil,numpatch) )

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,m) &
!$OMP SCHEDULE(STATIC, 1)
#endif
      DO i = 1, numpatch

         IF (DEF_forcing%has_missing_value) THEN
            IF (.not. forcmask(i)) cycle
         ENDIF
         
         m = patchclass(i)
         !TODO: can be removed
         z_soisno (maxsnl+1:0,i) = z_sno (maxsnl+1:0,i)
         z_soisno (1:nl_soil ,i) = z_soi (1:nl_soil)
         dz_soisno(maxsnl+1:0,i) = dz_sno(maxsnl+1:0,i)
         dz_soisno(1:nl_soil ,i) = dz_soi(1:nl_soil)

         !TODO: 整理变量次序
         CALL CLMMAIN (i, idate,           coszen(i),       deltim,          &
         patchlonr(i),    patchlatr(i),    patchclass(i),   patchtype(i),    &
         doalb,           dolai,           dosst,           oro(i),          &

       ! SOIL INFORMATION AND LAKE DEPTH
         soil_s_v_alb(i), soil_d_v_alb(i), soil_s_n_alb(i), soil_d_n_alb(i), &
         vf_quartz(1:,i), vf_gravels(1:,i),vf_om(1:,i),     vf_sand(1:,i),   &
         wf_gravels(1:,i),wf_sand(1:,i),   porsl(1:,i),     psi0(1:,i),      &
         bsw(1:,i),                                                          &
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         theta_r(1:,i),   alpha_vgm(1:,i), n_vgm(1:,i),     L_vgm(1:,i),     &
         sc_vgm (1:,i),   fc_vgm   (1:,i),                                   &
#endif
         hksati(1:,i),    csol(1:,i),      k_solids(1:,i),  dksatu(1:,i),    &
         dksatf(1:,i),    dkdry(1:,i),                                       &
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
         BA_alpha(1:,i),  BA_beta(1:,i),                                     &
#endif
         rootfr(1:,m),    lakedepth(i),    dz_lake(1:,i),                    &  
#if(defined CaMa_Flood)
                flddepth_cama(i),fldfrc_cama(i),fevpg_fld(i),  finfg_fld(i),        &!
#endif

       ! VEGETATION INFORMATION
         htop(i),         hbot(i),         sqrtdi(m),                        &
         effcon(m),       vmax25(m),                                         &
#ifdef PLANT_HYDRAULIC_STRESS
         kmax_sun(m),     kmax_sha(m),     kmax_xyl(m),     kmax_root(m),    &
         psi50_sun(m),    psi50_sha(m),    psi50_xyl(m),    psi50_root(m),   &
         ck(m),                                                              &
#endif
         slti(m),         hlti(m),                                           &
         shti(m),         hhti(m),         trda(m),         trdm(m),         &
         trop(m),         gradm(m),        binter(m),       extkn(m),        &
         chil(m),         rho(1:,1:,m),    tau(1:,1:,m),                     &

       ! ATMOSPHERIC FORCING
         forc_pco2m(i),   forc_po2m(i),    forc_us(i),      forc_vs(i),      &
         forc_t(i),       forc_q(i),       forc_prc(i),     forc_prl(i),     &
         forc_rain(i),    forc_snow(i),    forc_psrf(i),    forc_pbot(i),    &
         forc_sols(i),    forc_soll(i),    forc_solsd(i),   forc_solld(i),   &
         forc_frl(i),     forc_hgt_u(i),   forc_hgt_t(i),   forc_hgt_q(i),   &
         forc_rhoair(i),                                                     &

       ! LAND SURFACE VARIABLES REQUIRED FOR RESTART
         z_soisno(maxsnl+1:,i),            dz_soisno(maxsnl+1:,i),           &
         t_soisno(maxsnl+1:,i),            wliq_soisno(maxsnl+1:,i),         &
         wice_soisno(maxsnl+1:,i),         smp(1:,i),          hk(1:,i),     &
         t_grnd(i),       tleaf(i),        ldew(i),     ldew_rain(i),      ldew_snow(i),             &
         sag(i),          scv(i),          snowdp(i),       fveg(i),         &
         fsno(i),         sigf(i),         green(i),        lai(i),          &
         sai(i),          alb(1:,1:,i),    ssun(1:,1:,i),   ssha(1:,1:,i),   &
         thermk(i),       extkb(i),        extkd(i),                         &
#ifdef PLANT_HYDRAULIC_STRESS
         vegwp(1:,i),     gs0sun(i),       gs0sha(i),                        &
#endif
#ifdef OzoneStress
         lai_old(i),      o3uptakesun(i),  o3uptakesha(i)  ,forc_ozone(i),   &
#endif
         zwt(i),          dpond(i),        wa(i),                            &
         t_lake(1:,i),    lake_icefrac(1:,i),               savedtke1(i),    & 

       ! additional diagnostic variables for output
         laisun(i),       laisha(i),       rootr(1:,i),                      &
         rstfacsun(i),    rstfacsha(i),    h2osoi(1:,i),    wat(i),          &

       ! FLUXES
         taux(i),         tauy(i),         fsena(i),        fevpa(i),        &
         lfevpa(i),       fsenl(i),        fevpl(i),        etr(i),          &
         fseng(i),        fevpg(i),        olrg(i),         fgrnd(i),        &
         trad(i),         tref(i),         qref(i),         rsur(i),         &
         rnof(i),         qintr(i),        qinfl(i),        qdrip(i),        &
         rst(i),          assim(i),        respc(i),        sabvsun(i),      &
         sabvsha(i),      sabg(i),         sr(i),           solvd(i),        &
         solvi(i),        solnd(i),        solni(i),        srvd(i),         &
         srvi(i),         srnd(i),         srni(i),         solvdln(i),      &
         solviln(i),      solndln(i),      solniln(i),      srvdln(i),       &
         srviln(i),       srndln(i),       srniln(i),       qcharge(i),      &
         xerr(i),         zerr(i),                                           &

       ! TUNABLE modle constants
         zlnd,            zsno,            csoilc,          dewmx,           &
         wtfact,          capr,            cnfac,           ssi,             &
         wimp,            pondmx,          smpmax,          smpmin,          &
         trsmx0,          tcrit,                                             &

       ! additional variables required by coupling with WRF model
         emis(i),         z0m(i),          zol(i),          rib(i),          &
         ustar(i),        qstar(i),        tstar(i),                         &
         fm(i),           fh(i),           fq(i) )


         z_sno (maxsnl+1:0,i) = z_soisno (maxsnl+1:0,i)
         dz_sno(maxsnl+1:0,i) = dz_soisno(maxsnl+1:0,i)

#if(defined BGC)
         if(patchtype(i) .eq. 0)then
            CALL bgc_driver (i,idate(1:3),deltim, patchlatr(i)*180/PI,patchlonr(i)*180/PI)
         end if
#endif
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

  deallocate ( z_soisno  )
  deallocate ( dz_soisno )

END SUBROUTINE CLMDRIVER
! ---------- EOP ------------
