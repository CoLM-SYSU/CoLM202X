#include <define.h>

module MOD_1D_Acc_Fluxes

   use precision

   real(r8) :: nac              ! number of accumulation
   real(r8), allocatable :: nac_ln   (:)

   real(r8), allocatable :: a_us     (:)  
   real(r8), allocatable :: a_vs     (:)
   real(r8), allocatable :: a_t      (:)
   real(r8), allocatable :: a_q      (:)
   real(r8), allocatable :: a_prc    (:)
   real(r8), allocatable :: a_prl    (:)
   real(r8), allocatable :: a_pbot   (:)
   real(r8), allocatable :: a_frl    (:)
   real(r8), allocatable :: a_solarin(:)

   real(r8), allocatable :: a_taux   (:)
   real(r8), allocatable :: a_tauy   (:)
   real(r8), allocatable :: a_fsena  (:)
   real(r8), allocatable :: a_lfevpa (:)
   real(r8), allocatable :: a_fevpa  (:)
   real(r8), allocatable :: a_fsenl  (:)
   real(r8), allocatable :: a_fevpl  (:)
   real(r8), allocatable :: a_etr    (:)
   real(r8), allocatable :: a_fseng  (:)
   real(r8), allocatable :: a_fevpg  (:)
   real(r8), allocatable :: a_fgrnd  (:)
   real(r8), allocatable :: a_sabvsun(:)  
   real(r8), allocatable :: a_sabvsha(:) 
   real(r8), allocatable :: a_sabg   (:)
   real(r8), allocatable :: a_olrg   (:)
   real(r8), allocatable :: a_rnet   (:)
   real(r8), allocatable :: a_xerr   (:)
   real(r8), allocatable :: a_zerr   (:)
   real(r8), allocatable :: a_rsur   (:)
   real(r8), allocatable :: a_rnof   (:)
   real(r8), allocatable :: a_qintr  (:)
   real(r8), allocatable :: a_qinfl  (:)
   real(r8), allocatable :: a_qdrip  (:)
   real(r8), allocatable :: a_rstfacsun (:)
   real(r8), allocatable :: a_rstfacsha (:)
#ifdef VARIABLY_SATURATED_FLOW
   real(r8), allocatable :: a_dpond  (:)
#endif
   real(r8), allocatable :: a_zwt    (:)
   real(r8), allocatable :: a_wa     (:)
   real(r8), allocatable :: a_wat    (:)
   real(r8), allocatable :: a_assim  (:)
   real(r8), allocatable :: a_respc  (:)

   real(r8), allocatable :: a_qcharge(:)

   real(r8), allocatable :: a_t_grnd(:) 
   real(r8), allocatable :: a_tleaf (:) 
   real(r8), allocatable :: a_ldew  (:) 
   real(r8), allocatable :: a_scv   (:) 
   real(r8), allocatable :: a_snowdp(:) 
   real(r8), allocatable :: a_fsno  (:) 
   real(r8), allocatable :: a_sigf  (:) 
   real(r8), allocatable :: a_green (:) 
   real(r8), allocatable :: a_lai   (:) 
   real(r8), allocatable :: a_laisun(:) 
   real(r8), allocatable :: a_laisha(:) 
   real(r8), allocatable :: a_sai   (:) 

   real(r8), allocatable :: a_alb(:,:,:)    

   real(r8), allocatable :: a_emis (:)
   real(r8), allocatable :: a_z0m  (:)
   real(r8), allocatable :: a_trad (:)
   real(r8), allocatable :: a_tref (:)
   real(r8), allocatable :: a_qref (:)
   real(r8), allocatable :: a_rain (:)
   real(r8), allocatable :: a_snow (:)  

   real(r8), allocatable :: a_t_soisno    (:,:)    
   real(r8), allocatable :: a_wliq_soisno (:,:)
   real(r8), allocatable :: a_wice_soisno (:,:)
   real(r8), allocatable :: a_h2osoi      (:,:)
   real(r8), allocatable :: a_rootr       (:,:)
#ifdef PLANT_HYDRAULIC_STRESS
   real(r8), allocatable :: a_vegwp       (:,:)
#endif
   real(r8), allocatable :: a_t_lake      (:,:) 
   real(r8), allocatable :: a_lake_icefrac(:,:) 

   real(r8), allocatable :: a_ustar(:) 
   real(r8), allocatable :: a_tstar(:)
   real(r8), allocatable :: a_qstar(:)
   real(r8), allocatable :: a_zol  (:)
   real(r8), allocatable :: a_rib  (:)
   real(r8), allocatable :: a_fm   (:)
   real(r8), allocatable :: a_fh   (:)
   real(r8), allocatable :: a_fq   (:)

   real(r8), allocatable :: a_us10m(:) 
   real(r8), allocatable :: a_vs10m(:) 
   real(r8), allocatable :: a_fm10m(:) 

   real(r8), allocatable :: a_sr     (:)
   real(r8), allocatable :: a_solvd  (:)
   real(r8), allocatable :: a_solvi  (:)
   real(r8), allocatable :: a_solnd  (:)
   real(r8), allocatable :: a_solni  (:)
   real(r8), allocatable :: a_srvd   (:)
   real(r8), allocatable :: a_srvi   (:)
   real(r8), allocatable :: a_srnd   (:)
   real(r8), allocatable :: a_srni   (:)
   real(r8), allocatable :: a_solvdln(:)
   real(r8), allocatable :: a_solviln(:)
   real(r8), allocatable :: a_solndln(:)
   real(r8), allocatable :: a_solniln(:)
   real(r8), allocatable :: a_srvdln (:)
   real(r8), allocatable :: a_srviln (:)
   real(r8), allocatable :: a_srndln (:)
   real(r8), allocatable :: a_srniln (:)

   public :: allocate_acc_fluxes
   public :: deallocate_acc_fluxes
   public :: flush_acc_fluxes
   public :: accumulate_fluxes

contains

   subroutine allocate_acc_fluxes 

      use spmd_task
      USE GlobalVars
      use mod_landpatch, only : numpatch
      implicit none

      if (p_is_worker) then
         if (numpatch > 0) then

            allocate (a_us     (numpatch))  
            allocate (a_vs     (numpatch))
            allocate (a_t      (numpatch))
            allocate (a_q      (numpatch))
            allocate (a_prc    (numpatch))
            allocate (a_prl    (numpatch))
            allocate (a_pbot   (numpatch))
            allocate (a_frl    (numpatch))
            allocate (a_solarin(numpatch))

            allocate (a_taux      (numpatch))
            allocate (a_tauy      (numpatch))
            allocate (a_fsena     (numpatch))
            allocate (a_lfevpa    (numpatch))
            allocate (a_fevpa     (numpatch))
            allocate (a_fsenl     (numpatch))
            allocate (a_fevpl     (numpatch))
            allocate (a_etr       (numpatch))
            allocate (a_fseng     (numpatch))
            allocate (a_fevpg     (numpatch))
            allocate (a_fgrnd     (numpatch))
            allocate (a_sabvsun   (numpatch))  
            allocate (a_sabvsha   (numpatch)) 
            allocate (a_sabg      (numpatch))
            allocate (a_olrg      (numpatch))
            allocate (a_rnet      (numpatch))
            allocate (a_xerr      (numpatch))
            allocate (a_zerr      (numpatch))
            allocate (a_rsur      (numpatch))
            allocate (a_rnof      (numpatch))
            allocate (a_qintr     (numpatch))
            allocate (a_qinfl     (numpatch))
            allocate (a_qdrip     (numpatch))
            allocate (a_rstfacsun (numpatch))
            allocate (a_rstfacsha (numpatch))
#ifdef VARIABLY_SATURATED_FLOW
            allocate (a_dpond     (numpatch))
#endif
            allocate (a_zwt       (numpatch))
            allocate (a_wa        (numpatch))
            allocate (a_wat       (numpatch))
            allocate (a_assim     (numpatch))
            allocate (a_respc     (numpatch))

            allocate (a_qcharge   (numpatch))

            allocate (a_t_grnd    (numpatch)) 
            allocate (a_tleaf     (numpatch)) 
            allocate (a_ldew      (numpatch)) 
            allocate (a_scv       (numpatch)) 
            allocate (a_snowdp    (numpatch)) 
            allocate (a_fsno      (numpatch)) 
            allocate (a_sigf      (numpatch)) 
            allocate (a_green     (numpatch)) 
            allocate (a_lai       (numpatch)) 
            allocate (a_laisun    (numpatch)) 
            allocate (a_laisha    (numpatch)) 
            allocate (a_sai       (numpatch)) 

            allocate (a_alb   (2,2,numpatch))    

            allocate (a_emis      (numpatch))
            allocate (a_z0m       (numpatch))
            allocate (a_trad      (numpatch))
            allocate (a_tref      (numpatch))
            allocate (a_qref      (numpatch))
            allocate (a_rain      (numpatch))
            allocate (a_snow      (numpatch))  

            allocate (a_t_soisno    (maxsnl+1:nl_soil,numpatch))    
            allocate (a_wliq_soisno (maxsnl+1:nl_soil,numpatch))
            allocate (a_wice_soisno (maxsnl+1:nl_soil,numpatch))
            allocate (a_h2osoi      (1:nl_soil,       numpatch))
            allocate (a_rootr       (1:nl_soil,       numpatch))
#ifdef PLANT_HYDRAULIC_STRESS
            allocate (a_vegwp       (1:nvegwcs,       numpatch))
#endif
            allocate (a_t_lake      (nl_lake,numpatch)) 
            allocate (a_lake_icefrac(nl_lake,numpatch)) 

            allocate (a_ustar     (numpatch)) 
            allocate (a_tstar     (numpatch))
            allocate (a_qstar     (numpatch))
            allocate (a_zol       (numpatch))
            allocate (a_rib       (numpatch))
            allocate (a_fm        (numpatch))
            allocate (a_fh        (numpatch))
            allocate (a_fq        (numpatch))

            allocate (a_us10m     (numpatch)) 
            allocate (a_vs10m     (numpatch)) 
            allocate (a_fm10m     (numpatch)) 

            allocate (a_sr        (numpatch))
            allocate (a_solvd     (numpatch))
            allocate (a_solvi     (numpatch))
            allocate (a_solnd     (numpatch))
            allocate (a_solni     (numpatch))
            allocate (a_srvd      (numpatch))
            allocate (a_srvi      (numpatch))
            allocate (a_srnd      (numpatch))
            allocate (a_srni      (numpatch))
            allocate (a_solvdln   (numpatch))
            allocate (a_solviln   (numpatch))
            allocate (a_solndln   (numpatch))
            allocate (a_solniln   (numpatch))
            allocate (a_srvdln    (numpatch))
            allocate (a_srviln    (numpatch))
            allocate (a_srndln    (numpatch))
            allocate (a_srniln    (numpatch))

            allocate (nac_ln      (numpatch))

         end if
      end if

   end subroutine allocate_acc_fluxes

   subroutine deallocate_acc_fluxes ()

      use spmd_task
      use mod_landpatch, only : numpatch
      implicit none

      if (p_is_worker) then
         if (numpatch > 0) then

            deallocate (a_us     )  
            deallocate (a_vs     )
            deallocate (a_t      )
            deallocate (a_q      )
            deallocate (a_prc    )
            deallocate (a_prl    )
            deallocate (a_pbot   )
            deallocate (a_frl    )
            deallocate (a_solarin)

            deallocate (a_taux      )
            deallocate (a_tauy      )
            deallocate (a_fsena     )
            deallocate (a_lfevpa    )
            deallocate (a_fevpa     )
            deallocate (a_fsenl     )
            deallocate (a_fevpl     )
            deallocate (a_etr       )
            deallocate (a_fseng     )
            deallocate (a_fevpg     )
            deallocate (a_fgrnd     )
            deallocate (a_sabvsun   )  
            deallocate (a_sabvsha   ) 
            deallocate (a_sabg      )
            deallocate (a_olrg      )
            deallocate (a_rnet      )
            deallocate (a_xerr      )
            deallocate (a_zerr      )
            deallocate (a_rsur      )
            deallocate (a_rnof      )
            deallocate (a_qintr     )
            deallocate (a_qinfl     )
            deallocate (a_qdrip     )
            deallocate (a_rstfacsun )
            deallocate (a_rstfacsha )
#ifdef VARIABLY_SATURATED_FLOW
            deallocate (a_dpond     )
#endif
            deallocate (a_zwt       )
            deallocate (a_wa        )
            deallocate (a_wat       )
            deallocate (a_assim     )
            deallocate (a_respc     )

            deallocate (a_qcharge   )

            deallocate (a_t_grnd    ) 
            deallocate (a_tleaf     ) 
            deallocate (a_ldew      ) 
            deallocate (a_scv       ) 
            deallocate (a_snowdp    ) 
            deallocate (a_fsno      ) 
            deallocate (a_sigf      ) 
            deallocate (a_green     ) 
            deallocate (a_lai       ) 
            deallocate (a_laisun    ) 
            deallocate (a_laisha    ) 
            deallocate (a_sai       ) 

            deallocate (a_alb  ) 

            deallocate (a_emis      )
            deallocate (a_z0m       )
            deallocate (a_trad      )
            deallocate (a_tref      )
            deallocate (a_qref      )
            deallocate (a_rain      )
            deallocate (a_snow      )  

            deallocate (a_t_soisno    )    
            deallocate (a_wliq_soisno )
            deallocate (a_wice_soisno )
            deallocate (a_h2osoi      ) 
            deallocate (a_rootr       )
#ifdef PLANT_HYDRAULIC_STRESS
            deallocate (a_vegwp       )
#endif
            deallocate (a_t_lake      ) 
            deallocate (a_lake_icefrac) 

            deallocate (a_ustar     ) 
            deallocate (a_tstar     )
            deallocate (a_qstar     )
            deallocate (a_zol       )
            deallocate (a_rib       )
            deallocate (a_fm        )
            deallocate (a_fh        )
            deallocate (a_fq        )

            deallocate (a_us10m     ) 
            deallocate (a_vs10m     ) 
            deallocate (a_fm10m     ) 

            deallocate (a_sr        )
            deallocate (a_solvd     )
            deallocate (a_solvi     )
            deallocate (a_solnd     )
            deallocate (a_solni     )
            deallocate (a_srvd      )
            deallocate (a_srvi      )
            deallocate (a_srnd      )
            deallocate (a_srni      )
            deallocate (a_solvdln   )
            deallocate (a_solviln   )
            deallocate (a_solndln   )
            deallocate (a_solniln   )
            deallocate (a_srvdln    )
            deallocate (a_srviln    )
            deallocate (a_srndln    )
            deallocate (a_srniln    )

            deallocate (nac_ln      )

         end if
      end if

   end subroutine deallocate_acc_fluxes

   !-----------------------
   SUBROUTINE FLUSH_acc_fluxes ()

      use spmd_task
      use mod_landpatch, only : numpatch
      use GlobalVars,    only : spval 
      implicit none

      if (p_is_worker) then

         nac = 0

         if (numpatch > 0) then

            ! flush the Fluxes for accumulation
            a_us     (:) = spval
            a_vs     (:) = spval
            a_t      (:) = spval
            a_q      (:) = spval
            a_prc    (:) = spval
            a_prl    (:) = spval
            a_pbot   (:) = spval
            a_frl    (:) = spval
            a_solarin(:) = spval

            a_taux    (:) = spval
            a_tauy    (:) = spval
            a_fsena   (:) = spval
            a_lfevpa  (:) = spval
            a_fevpa   (:) = spval
            a_fsenl   (:) = spval
            a_fevpl   (:) = spval
            a_etr     (:) = spval
            a_fseng   (:) = spval
            a_fevpg   (:) = spval
            a_fgrnd   (:) = spval
            a_sabvsun (:) = spval
            a_sabvsha (:) = spval
            a_sabg    (:) = spval
            a_olrg    (:) = spval
            a_rnet    (:) = spval
            a_xerr    (:) = spval
            a_zerr    (:) = spval
            a_rsur    (:) = spval
            a_rnof    (:) = spval
            a_qintr   (:) = spval
            a_qinfl   (:) = spval
            a_qdrip   (:) = spval
            a_rstfacsun(:) = spval
            a_rstfacsha(:) = spval
#ifdef VARIABLY_SATURATED_FLOW
            a_dpond   (:) = spval
#endif
            a_zwt     (:) = spval
            a_wa      (:) = spval
            a_wat     (:) = spval
            a_assim   (:) = spval
            a_respc   (:) = spval

            a_qcharge (:) = spval

            a_t_grnd  (:) = spval
            a_tleaf   (:) = spval
            a_ldew    (:) = spval
            a_scv     (:) = spval
            a_snowdp  (:) = spval
            a_fsno    (:) = spval
            a_sigf    (:) = spval
            a_green   (:) = spval
            a_lai     (:) = spval
            a_laisun  (:) = spval
            a_laisha  (:) = spval
            a_sai     (:) = spval
            
            a_alb   (:,:,:) = spval

            a_emis      (:) = spval
            a_z0m       (:) = spval
            a_trad      (:) = spval
            a_tref      (:) = spval
            a_qref      (:) = spval 
            a_rain      (:) = spval
            a_snow      (:) = spval

            a_t_soisno     (:,:) = spval 
            a_wliq_soisno  (:,:) = spval
            a_wice_soisno  (:,:) = spval
            a_h2osoi       (:,:) = spval
            a_rootr        (:,:) = spval
#ifdef PLANT_HYDRAULIC_STRESS
            a_vegwp        (:,:) = spval
#endif
            a_t_lake       (:,:) = spval
            a_lake_icefrac (:,:) = spval

            a_ustar (:) = spval 
            a_tstar (:) = spval
            a_qstar (:) = spval
            a_zol   (:) = spval
            a_rib   (:) = spval
            a_fm    (:) = spval
            a_fh    (:) = spval
            a_fq    (:) = spval

            a_us10m (:) = spval
            a_vs10m (:) = spval
            a_fm10m (:) = spval

            a_sr       (:) = spval
            a_solvd    (:) = spval
            a_solvi    (:) = spval
            a_solnd    (:) = spval
            a_solni    (:) = spval
            a_srvd     (:) = spval
            a_srvi     (:) = spval
            a_srnd     (:) = spval
            a_srni     (:) = spval
            a_solvdln  (:) = spval
            a_solviln  (:) = spval
            a_solndln  (:) = spval
            a_solniln  (:) = spval
            a_srvdln   (:) = spval
            a_srviln   (:) = spval
            a_srndln   (:) = spval
            a_srniln   (:) = spval

            nac_ln  (:) = 0

         end if
      end if

   END SUBROUTINE FLUSH_acc_fluxes

   SUBROUTINE accumulate_fluxes 
      ! ----------------------------------------------------------------------
      ! perfrom the grid average mapping: average a subgrid input 1d vector 
      ! of length numpatch to a output 2d array of length [ghist%xcnt,ghist%ycnt]
      !
      ! Created by Yongjiu Dai, 03/2014
      !---------------------------------------------------------------------

      use precision
      use spmd_task
      use mod_landpatch,     only : numpatch
      use PhysicalConstants, only : vonkar, stefnc, cpair, rgas, grav
      use MOD_TimeInvariants
      use MOD_TimeVariables
      use MOD_1D_Forcing
      use MOD_1D_Fluxes
      use FRICTION_VELOCITY
      use mod_colm_debug
      use GlobalVars

      IMPLICIT NONE

      ! Local Variables

      real(r8), allocatable :: r_trad  (:)

      real(r8), allocatable :: r_ustar (:)
      real(r8), allocatable :: r_tstar (:)
      real(r8), allocatable :: r_qstar (:)
      real(r8), allocatable :: r_zol   (:)
      real(r8), allocatable :: r_rib   (:)
      real(r8), allocatable :: r_fm    (:)
      real(r8), allocatable :: r_fh    (:)
      real(r8), allocatable :: r_fq    (:)

      real(r8), allocatable :: r_us10m (:)
      real(r8), allocatable :: r_vs10m (:)
      real(r8), allocatable :: r_fm10m (:)

      !---------------------------------------------------------------------
      integer  ib, jb, i, j
      real(r8) rhoair,thm,th,thv,ur,displa_av,zldis,hgt_u,hgt_t,hgt_q
      real(r8) z0m_av,z0h_av,z0q_av,us,vs,tm,qm,psrf
      real(r8) obu,fh2m,fq2m
      real(r8) um,thvstar,beta,zii,wc,wc2

      if (p_is_worker) then
         if (numpatch > 0) then

            nac = nac + 1

            call acc1d (forc_us  , a_us  )
            call acc1d (forc_vs  , a_vs  )
            call acc1d (forc_t   , a_t   )
            call acc1d (forc_q   , a_q   )
            call acc1d (forc_prc , a_prc )
            call acc1d (forc_prl , a_prl )
            call acc1d (forc_pbot, a_pbot)
            call acc1d (forc_frl , a_frl )

            call acc1d (forc_sols,  a_solarin)
            call acc1d (forc_soll,  a_solarin)
            call acc1d (forc_solsd, a_solarin)
            call acc1d (forc_solld, a_solarin)

            call acc1d (taux    , a_taux   )
            call acc1d (tauy    , a_tauy   )
            call acc1d (fsena   , a_fsena  )
            call acc1d (lfevpa  , a_lfevpa )
            call acc1d (fevpa   , a_fevpa  )
            call acc1d (fsenl   , a_fsenl  )
            call acc1d (fevpl   , a_fevpl  )
            call acc1d (etr     , a_etr    )
            call acc1d (fseng   , a_fseng  )
            call acc1d (fevpg   , a_fevpg  )
            call acc1d (fgrnd   , a_fgrnd  )
            call acc1d (sabvsun , a_sabvsun)
            call acc1d (sabvsha , a_sabvsha)
            call acc1d (sabg    , a_sabg   )
            call acc1d (olrg    , a_olrg   )

            rnet = sabg + sabvsun + sabvsha - olrg + forc_frl
            call acc1d (rnet    , a_rnet   )

            call acc1d (xerr   , a_xerr   )
            call acc1d (zerr   , a_zerr   )
            call acc1d (rsur   , a_rsur   )
            call acc1d (rnof   , a_rnof   )
            call acc1d (qintr  , a_qintr  )
            call acc1d (qinfl  , a_qinfl  )
            call acc1d (qdrip  , a_qdrip  )
            call acc1d (rstfacsun , a_rstfacsun )
            call acc1d (rstfacsha , a_rstfacsha )
#ifdef VARIABLY_SATURATED_FLOW
            call acc1d (dpond  , a_dpond  )
#endif
            call acc1d (zwt    , a_zwt    )
            call acc1d (wa     , a_wa     )
            call acc1d (wat    , a_wat    )
            call acc1d (assim  , a_assim  )
            call acc1d (respc  , a_respc  )

            call acc1d (qcharge, a_qcharge)

            call acc1d (t_grnd , a_t_grnd )
            call acc1d (tleaf  , a_tleaf  )
            call acc1d (ldew   , a_ldew   )
            call acc1d (scv    , a_scv    )
            call acc1d (snowdp , a_snowdp )
            call acc1d (fsno   , a_fsno   )
            call acc1d (sigf   , a_sigf   )
            call acc1d (green  , a_green  )
            call acc1d (lai    , a_lai    )
            call acc1d (laisun , a_laisun )
            call acc1d (laisha , a_laisha )
            call acc1d (sai    , a_sai    )

            call acc3d (alb    , a_alb    )

            call acc1d (emis   , a_emis   )
            call acc1d (z0m    , a_z0m    )

            allocate (r_trad (numpatch))
            do i = 1, numpatch
               r_trad(i) = (olrg(i)/stefnc)**0.25
            end do
            call acc1d (r_trad , a_trad   )
            deallocate (r_trad )

            call acc1d (tref   , a_tref   )
            call acc1d (qref   , a_qref   )

            call acc1d (forc_rain, a_rain )
            call acc1d (forc_snow, a_snow )

            call acc2d (t_soisno   , a_t_soisno   )
            call acc2d (wliq_soisno, a_wliq_soisno)
            call acc2d (wice_soisno, a_wice_soisno)

            call acc2d (h2osoi     , a_h2osoi     )
            call acc2d (rootr      , a_rootr      )
#ifdef PLANT_HYDRAULIC_STRESS
            call acc2d (vegwp      , a_vegwp      )
#endif
            call acc2d (t_lake      , a_t_lake      )
            call acc2d (lake_icefrac, a_lake_icefrac)

            allocate (r_ustar (numpatch))
            allocate (r_tstar (numpatch))
            allocate (r_qstar (numpatch))
            allocate (r_zol   (numpatch))
            allocate (r_rib   (numpatch))
            allocate (r_fm    (numpatch))
            allocate (r_fh    (numpatch))
            allocate (r_fq    (numpatch))

            allocate (r_us10m (numpatch))
            allocate (r_vs10m (numpatch))
            allocate (r_fm10m (numpatch))

            do i = 1, numpatch

               z0m_av = z0m(i) 
               z0h_av = z0m(i)
               z0q_av = z0m(i)

               displa_av = 2./3.*z0m_av/0.07

               hgt_u = max(forc_hgt_u(i), 5.+displa_av)
               hgt_t = max(forc_hgt_t(i), 5.+displa_av)
               hgt_q = max(forc_hgt_q(i), 5.+displa_av)
               zldis = hgt_u-displa_av

               us = forc_us(i)
               vs = forc_vs(i)
               tm = forc_t (i)
               qm = forc_q (i)
               psrf = forc_psrf(i)
               rhoair = (psrf - 0.378*qm*psrf/(0.622+0.378*qm)) / (rgas*tm)

               r_ustar(i) = sqrt(max(1.e-6,sqrt(taux(i)**2+tauy(i)**2))/rhoair)
               r_tstar(i) = -fsena(i)/(rhoair*r_ustar(i))/cpair
               r_qstar(i) = -fevpa(i)/(rhoair*r_ustar(i))

               thm = tm + 0.0098*hgt_t
               th  = tm*(100000./psrf)**(rgas/cpair)
               thv = th*(1.+0.61*qm)

               r_zol(i) = zldis*vonkar*grav * (r_tstar(i)+0.61*th*r_qstar(i)) &
                  / (r_ustar(i)**2*thv)

               if(r_zol(i) >= 0.)then   !stable
                  r_zol(i) = min(2.,max(r_zol(i),1.e-6))
               else                       !unstable
                  r_zol(i) = max(-100.,min(r_zol(i),-1.e-6))
               endif

               beta = 1.
               zii = 1000.
               thvstar=r_tstar(i)+0.61*th*r_qstar(i)
               ur = sqrt(us*us+vs*vs)
               if(r_zol(i) >= 0.)then
                  um = max(ur,0.1)
               else
                  wc = (-grav*r_ustar(i)*thvstar*zii/thv)**(1./3.)
                  wc2 = beta*beta*(wc*wc)
                  um = max(0.1,sqrt(ur*ur+wc2))
               endif

               obu = zldis/r_zol(i)
               call moninobuk(hgt_u,hgt_t,hgt_q,displa_av,z0m_av,z0h_av,z0q_av,&
                  obu,um,r_ustar(i),fh2m,fq2m,r_fm10m(i),r_fm(i),r_fh(i),r_fq(i))

               ! bug found by chen qiying 2013/07/01 
               r_rib(i) = r_zol(i) /vonkar * r_ustar(i)**2 / (vonkar/r_fh(i)*um**2)
               r_rib(i) = min(5.,r_rib(i))

               r_us10m(i) = us/um * r_ustar(i) /vonkar * r_fm10m(i)
               r_vs10m(i) = vs/um * r_ustar(i) /vonkar * r_fm10m(i)

            end do

            call acc1d (r_ustar, a_ustar)
            call acc1d (r_tstar, a_tstar)
            call acc1d (r_qstar, a_qstar)
            call acc1d (r_zol  , a_zol  )
            call acc1d (r_rib  , a_rib  )
            call acc1d (r_fm   , a_fm   )
            call acc1d (r_fh   , a_fh   )
            call acc1d (r_fq   , a_fq   )

            call acc1d (r_us10m, a_us10m)
            call acc1d (r_vs10m, a_vs10m)
            call acc1d (r_fm10m, a_fm10m)

            deallocate (r_ustar )
            deallocate (r_tstar )
            deallocate (r_qstar )
            deallocate (r_zol   )
            deallocate (r_rib   )
            deallocate (r_fm    )
            deallocate (r_fh    )
            deallocate (r_fq    )

            deallocate (r_us10m )
            deallocate (r_vs10m )
            deallocate (r_fm10m )

            call acc1d (sr     , a_sr     )
            call acc1d (solvd  , a_solvd  )
            call acc1d (solvi  , a_solvi  )
            call acc1d (solnd  , a_solnd  )
            call acc1d (solni  , a_solni  )
            call acc1d (srvd   , a_srvd   )
            call acc1d (srvi   , a_srvi   )
            call acc1d (srnd   , a_srnd   )
            call acc1d (srni   , a_srni   )
            call acc1d (solvdln, a_solvdln)
            call acc1d (solviln, a_solviln)
            call acc1d (solndln, a_solndln)
            call acc1d (solniln, a_solniln)
            call acc1d (srvdln , a_srvdln )
            call acc1d (srviln , a_srviln )
            call acc1d (srndln , a_srndln )
            call acc1d (srniln , a_srniln )

            do i = 1, numpatch
               if (solvdln(i) /= spval) then
                  nac_ln(i) = nac_ln(i) + 1
               end if
            end do

         end if
      end if

   END SUBROUTINE accumulate_fluxes


   !------
   SUBROUTINE acc1d (var, s)

      use precision
      use GlobalVars, only: spval

      IMPLICIT NONE

      real(r8), intent(in)    :: var(:)
      real(r8), intent(inout) :: s  (:)
      ! Local variables
      integer :: i

      do i = lbound(var,1), ubound(var,1)
         if (var(i) /= spval) then
            if (s(i) /= spval) then
               s(i) = s(i) + var(i)
            else
               s(i) = var(i)
            end if
         end if
      end do
      
   END SUBROUTINE acc1d

   !------
   SUBROUTINE acc2d (var, s)

      use precision
      use GlobalVars, only: spval

      IMPLICIT NONE

      real(r8), intent(in)    :: var(:,:)
      real(r8), intent(inout) :: s  (:,:)
      ! Local variables
      integer :: i1, i2

      do i1 = lbound(var,1), ubound(var,1)
         do i2 = lbound(var,2), ubound(var,2)
            if (var(i1,i2) /= spval) then
               if (s(i1,i2) /= spval) then
                  s(i1,i2) = s(i1,i2) + var(i1,i2)
               else
                  s(i1,i2) = var(i1,i2)
               end if
            end if
         end do
      end do

   END SUBROUTINE acc2d

   !------
   SUBROUTINE acc3d (var, s)

      use precision
      use GlobalVars, only: spval

      IMPLICIT NONE

      real(r8), intent(in)    :: var(:,:,:)
      real(r8), intent(inout) :: s  (:,:,:)
      ! Local variables
      integer :: i1, i2, i3

      do i1 = lbound(var,1), ubound(var,1)
         do i2 = lbound(var,2), ubound(var,2)
            do i3 = lbound(var,3), ubound(var,3)
               if (var(i1,i2,i3) /= spval) then
                  if (s(i1,i2,i3) /= spval) then
                     s(i1,i2,i3) = s(i1,i2,i3) + var(i1,i2,i3)
                  else
                     s(i1,i2,i3) = var(i1,i2,i3)
                  end if
               end if
            end do
         end do
      end do

   END SUBROUTINE acc3d

end module MOD_1D_Acc_Fluxes
! ----- EOP ---------
