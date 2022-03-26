#include <define.h>
module bgc_veg_CNPhenologyMod

  use PFT_Const, only: &
      isevg  , issed  , isstd  , leaf_long, woody  , leafcn , frootcn, livewdcn, deadwdcn, &
      lflitcn, lf_flab, lf_fcel, lf_flig  , fr_flab, fr_fcel, fr_flig, &

! crop variables
      minplanttemp, minplantjday, maxplantjday, hybgdd, manunitro, planttemp, lfemerg, mxmat, grnfill, mxtmp, &
      baset  , gddmin 

  use MOD_TimeInvariants, only: &
      ndays_on        , ndays_off      , fstor2tran, crit_dayl  , crit_onset_fdd, crit_onset_swi, &
      crit_offset_fdd , crit_offset_swi, soilpsi_on, soilpsi_off, lwtop

  use GlobalVars, only: &
 !crop variables
      nswheat         , nirrig_swheat     , nsugarcane  , nirrig_sugarcane  , &
      nwwheat         , nirrig_wwheat     , ntmp_corn   , nirrig_tmp_corn   , &
      ntrp_corn       , nirrig_trp_corn   , nmiscanthus , nirrig_miscanthus , &
      nswitchgrass    , nirrig_switchgrass, ncotton     , nirrig_cotton     , &
      nrice           , nirrig_rice       , ntmp_soybean, nirrig_tmp_soybean, &
      ntrp_soybean    , nirrig_trp_soybean, &
      spval

  use MOD_TimeVariables, only: &
      t_soisno, smp, dayl, prev_dayl, prec10, prec60, prec365, prec_today, prec_daily, accumnstep

  use MOD_PFTimeVars, only: &
      tref_p       , tempavg_tref_p , annavg_tref_p  , gdd0_p        , gdd8_p            , &
      gdd10_p      , gdd020_p       , gdd820_p       , gdd1020_p     , nyrs_crop_active_p, &
      bglfr_p      , bgtr_p         , lgsf_p         , offset_flag_p , offset_counter_p  , &
      onset_flag_p , onset_counter_p, onset_gddflag_p, onset_gdd_p   , onset_fdd_p       , &
      onset_swi_p  , offset_fdd_p   , offset_swi_p   , dormant_flag_p, tlai_p            , &

      prev_leafc_to_litter_p        , prev_frootc_to_litter_p        , days_active_p     , &

      leafc_p            , frootc_p            , livestemc_p         , &
      livestemn_p        , livecrootc_p        , grainc_p, grainn_p  , &

      leafc_storage_p    , frootc_storage_p    , livestemc_storage_p , & 
      deadstemc_storage_p, livecrootc_storage_p, deadcrootc_storage_p, &
      leafn_storage_p    , frootn_storage_p    , livestemn_storage_p , &
      deadstemn_storage_p, livecrootn_storage_p, deadcrootn_storage_p, &

      leafc_xfer_p       , frootc_xfer_p       , livestemc_xfer_p    , &
      deadstemc_xfer_p   , livecrootc_xfer_p   , deadcrootc_xfer_p   , &
      leafn_xfer_p       , frootn_xfer_p       , livestemn_xfer_p    , &
      deadstemn_xfer_p   , livecrootn_xfer_p   , deadcrootn_xfer_p   , &
      gresp_storage_p    , &

      leaf_prof_p        , froot_prof_p        , &
      cropseedc_deficit_p, cropseedn_deficit_p , &

! crop variables
      cropplant_p       , idop_p              , a5tmin_p            , a10tmin_p        , t10_p          , &
      cumvd_p           , hdidx_p             , vf_p                , cphase_p         , fert_counter_p , &
      croplive_p        , gddplant_p          , gddtsoi_p           , harvdate_p       , gddmaturity_p  , &
      huigrain_p        , huileaf_p        , peaklai_p      , &
      fert_p            , tref_min_p          , tref_max_p          , tref_min_inst_p  , tref_max_inst_p, &
      fertnitro_p       , latbaset_p    ! input from files
      
  use MOD_1D_PFTFluxes, only: &
      livestemc_to_deadstemc_p       , livecrootc_to_deadcrootc_p   , &

      leafc_storage_to_xfer_p        , frootc_storage_to_xfer_p     , &
      livestemc_storage_to_xfer_p    , deadstemc_storage_to_xfer_p  , &
      livecrootc_storage_to_xfer_p   , deadcrootc_storage_to_xfer_p , &
      gresp_storage_to_xfer_p        , &
      
      leafc_xfer_to_leafc_p          , frootc_xfer_to_frootc_p        , &
      livestemc_xfer_to_livestemc_p  , deadstemc_xfer_to_deadstemc_p  , &
      livecrootc_xfer_to_livecrootc_p, deadcrootc_xfer_to_deadcrootc_p, &

      livestemn_to_deadstemn_p       , livecrootn_to_deadcrootn_p     , &
      livestemn_to_retransn_p        , livecrootn_to_retransn_p       , &

      leafn_storage_to_xfer_p        , frootn_storage_to_xfer_p       , &
      livestemn_storage_to_xfer_p    , deadstemn_storage_to_xfer_p    , &
      livecrootn_storage_to_xfer_p   , deadcrootn_storage_to_xfer_p   , &
      
      leafn_xfer_to_leafn_p          , frootn_xfer_to_frootn_p        , &
      livestemn_xfer_to_livestemn_p  , deadstemn_xfer_to_deadstemn_p  , &
      livecrootn_xfer_to_livecrootn_p, deadcrootn_xfer_to_deadcrootn_p, &
      cpool_to_leafc_p               , cpool_to_frootc_p              , &
      leafc_to_litter_p              , frootc_to_litter_p             , &
      leafn_to_litter_p              , frootn_to_litter_p             , &
      leafn_to_retransn_p            , &

      crop_seedc_to_leaf_p           , crop_seedn_to_leaf_p           , &
      grainc_to_seed_p               , grainn_to_seed_p               , &
      grainc_to_food_p               , grainn_to_food_p               , &
      cpool_to_grainc_p              , npool_to_grainn_p              , &
      livestemc_to_litter_p          , livestemn_to_litter_p          , &
      cpool_to_livestemc_p                

  use MOD_PFTimeInvars, only: pftclass, pftfrac

  use MOD_1D_Fluxes, only: &
      phenology_to_met_c , phenology_to_cel_c , phenology_to_lig_c, &
      phenology_to_met_n , phenology_to_cel_n , phenology_to_lig_n, &
      grainc_to_cropprodc, grainn_to_cropprodn

  use MOD_1D_forcing, only: forc_prc, forc_prl

  use timemanager
  use precision
  use bgc_DaylengthMod, only: daylength

implicit none

public CNPhenology

integer, parameter :: NOT_Planted   = 999 ! If not planted   yet in year
integer, parameter :: NOT_Harvested = 999 ! If not harvested yet in year

contains

subroutine CNPhenology(i,ps,pe,nl_soil,idate,dz_soi,deltim,dlat,npcropmin,phase)

implicit none

integer ,intent(in) :: i
integer ,intent(in) :: ps
integer ,intent(in) :: pe
integer ,intent(in) :: nl_soil
integer ,intent(in) :: idate(3)
real(r8),intent(in) :: deltim
real(r8),intent(in) :: dlat
integer ,intent(in) :: npcropmin
real(r8),intent(in) :: dz_soi(nl_soil)
integer ,intent(in) :: phase

real(r8) dayspyr
integer  h,m    ! 1 for north hemsiphere; 2 for south hemisphere

   if(isleapyear(idate(1)))then
      dayspyr = 366
   else
      dayspyr = 365
   end if

   if(dlat > 0)then
      h = 1
   else 
      h = 2
   end if
 if ( phase == 1 ) then
    call CNPhenologyClimate    (i,ps,pe,idate(1:3),deltim,dayspyr,npcropmin,nl_soil,dz_soi,dlat)

    call CNEvergreenPhenology  (i,ps,pe,deltim,dayspyr)

    call CNSeasonDecidPhenology(i,ps,pe,idate(1:3),deltim,dayspyr,dlat)

    call CNStressDecidPhenology(i,ps,pe,deltim,dayspyr)

#ifdef CROP
!   if (doalb .and. num_pcropp > 0 ) then
    if(dlat >= 0)then
       h = 1
    else
       h = 2
    end if
    call CropPhenology(i,ps,pe,idate(1:3),h,deltim,dayspyr,npcropmin)
!   end if
#endif
 else if ( phase == 2 ) then
    ! the same onset and offset routines are called regardless of
    ! phenology type - they depend only on onset_flag, offset_flag, bglfr, and bgtr

    call CNOnsetGrowth(i,ps,pe,deltim)

    call CNOffsetLitterfall(i,ps,pe,deltim,npcropmin)

    call CNBackgroundLitterfall(i,ps,pe)

    call CNLivewoodTurnover(i,ps,pe)

    call CNGrainToProductPools(i,ps,pe)

    call CNLitterToColumn(i,ps,pe,nl_soil,npcropmin)
 else
    write(*,*) 'bad phenology phase'
 end if

end subroutine CNPhenology

subroutine CNPhenologyClimate (i,ps,pe,idate,deltim,dayspyr,npcropmin,nl_soil,dz_soi,dlat)
   !
   ! !DESCRIPTION:
   integer ,intent(in) :: i
   integer ,intent(in) :: ps
   integer ,intent(in) :: pe
   integer ,intent(in) :: idate(3)
   real(r8),intent(in) :: deltim
   real(r8),intent(in) :: dayspyr ! days per year (days)
   integer ,intent(in) :: npcropmin
   integer ,intent(in) :: nl_soil
   real(r8),intent(in) :: dz_soi(nl_soil)
   real(r8),intent(in) :: dlat
   ! For coupled carbon-nitrogen code (CN).
   !
   ! !LOCAL VARIABLES:
   real(r8), parameter :: yravg   = 20.0_r8      ! length of years to average for gdd
   real(r8), parameter :: yravgm1 = yravg-1.0_r8 ! minus 1 of above
   integer :: m,ivt
   logical , parameter :: isconst_baset = .true. ! .true. for constant base temperature
                                                 ! .false. for latidinal varied base temperature
   real(r8) stepperday,  nsteps
   integer month, mday
   !-----------------------------------------------------------------------

     ! set time steps

   stepperday = 86400._r8 / deltim
   do m = ps , pe
      tempavg_tref_p(m) = tempavg_tref_p(m) + tref_p(m) * (deltim/86400._r8/dayspyr)
#ifdef CROP
      if(idate(3) .eq. 1800 .or. tref_max_inst_p(m) .eq. spval)then
         tref_max_inst_p(m) = tref_p(m)
      else
         tref_max_inst_p(m) = max(tref_max_inst_p(m) , tref_p(m))
      end if

      if(idate(3) .eq. 1800 .or. tref_min_inst_p(m) .eq. spval)then
         tref_min_inst_p(m) = tref_p(m)
      else
         tref_min_inst_p(m) = min(tref_min_inst_p(m) , tref_p(m))
      end if
      if(idate(3) .eq. 84600)then
         tref_max_p(m) = tref_max_inst_p(m)
         tref_min_p(m) = tref_min_inst_p(m)
      end if
#endif
   end do

   accumnstep(i) = accumnstep(i) + 1
   prec_today(i) = forc_prc(i) + forc_prl(i)
#ifdef CROP

   nsteps = amin1(5._r8 * stepperday, accumnstep(i))
   do m = ps , pe
      a5tmin_p (m) = (a5tmin_p(m)  * (nsteps - 1) + tref_min_p(m) ) / nsteps
   end do
#endif
   nsteps = amin1(10._r8 * stepperday, accumnstep(i))
   prec10     (i)  = ( prec10  (i) * (nsteps - 1) + prec_today(i) ) / nsteps
#ifdef CROP
   do m = ps , pe
      t10_p    (m) = ( t10_p   (m) * (nsteps - 1) + tref_p    (m) ) / nsteps
      a10tmin_p(m) = (a10tmin_p(m) * (nsteps - 1) + tref_min_p(m) ) / nsteps
   end do
#endif

   nsteps = amin1(60._r8 * stepperday, accumnstep(i))
   prec60     (i) = ( prec60 (i) * (nsteps - 1) + prec_today(i) ) / nsteps

   nsteps = amin1(365._r8 * stepperday, accumnstep(i))
   prec365    (i) = ( prec365(i) * (nsteps - 1) + prec_today(i) ) / nsteps

   call julian2monthday(idate(1),idate(2),month,mday)
   do m = ps , pe
      ivt = pftclass(m)
      if(((month .ge. 4 .and. month .le. 9) .and. dlat .ge. 0) .or. &
         ((month .gt. 9 .or.  month .lt. 4) .and. dlat .lt. 0)) then
         gdd0_p (m) = gdd0_p (m) + max(0._r8, min(26._r8, &
                            tref_p(m) - 273.15)) * deltim / 86400._r8
         gdd8_p (m) = gdd8_p (m) + max(0._r8, min(30._r8, &
                            tref_p(m) - 273.15 - 8)) * deltim / 86400._r8
         gdd10_p(m) = gdd10_p(m) + max(0._r8, min(30._r8, &
                            tref_p(m) - 273.15 - 10)) * deltim / 86400._r8
      end if
#ifdef CROP
      if(croplive_p(m))then
         if(ivt == nwwheat .or. ivt == nirrig_wwheat)then
            gddplant_p(m) = gddplant_p(m) + vf_p(m)
            gddtsoi_p(m)  = gddtsoi_p (m) + vf_p(m)
         else
            if((.not. isconst_baset) .and. ((ivt == nswheat) .or. (ivt == nirrig_swheat) .or. &
                                            (ivt == nsugarcane) .or. (ivt == nirrig_sugarcane)))then
               gddplant_p(m) = gddplant_p(m) + max(0._r8, min(mxtmp(ivt), &
                               tref_p(m) - (273.15 + latbaset_p(m)))) * deltim / 86400._r8
            else
               gddplant_p(m) = gddplant_p(m) + max(0._r8, min(mxtmp(ivt), &
                               tref_p(m) - (273.15 + baset(ivt)))) * deltim / 86400._r8
            end if
            gddtsoi_p(m) = gddtsoi_p(m) +  max(0._r8, min(mxtmp(ivt), &
                           (t_soisno(1,i) * dz_soi(1) + t_soisno(2,i) * dz_soi(2)) &
                            / (dz_soi(1) + dz_soi(2)))) * deltim / 86400._r8 
         end if
      else
         gddplant_p(m) = 0._r8
         gddtsoi_p(m) = 0._r8
      end if
#endif
   end do

!   if(idate(2) .eq. dayspyr .and. idate(3) .eq. 86400)then
!      annavg_tref(i) = tempavg_tref(i)
!      tempavg_tref(i) = 0
!   end if

     !
     ! The following crop related steps are done here rather than CropPhenology
     ! so that they will be completed each time-step rather than with doalb.
     !
     ! The following lines come from ibis's climate.f + stats.f
     ! gdd SUMMATIONS ARE RELATIVE TO THE PLANTING DATE (see subr. updateAccFlds)
   do m = ps , pe
      ivt = pftclass(m)
      if (idate(2) == 1 .and. idate(3) ==1800 .and. nyrs_crop_active_p(m) == 0) then ! YR 1:
         gdd020_p(m)  = 0._r8                             ! set gdd..20 variables to 0
         gdd820_p(m)  = 0._r8                             ! and crops will not be planted
         gdd1020_p(m) = 0._r8
      end if
      if (isendofyear(idate,deltim)) then        ! <-- END of EVERY YR:
         nyrs_crop_active_p(m) = nyrs_crop_active_p(m) + 1
         if (nyrs_crop_active_p(m) == 1) then                     ! <-- END of YR 1
            gdd020_p(m)  = gdd0_p(m)                                ! <-- END of YR 1
            gdd820_p(m)  = gdd8_p(m)                                ! <-- END of YR 1
            gdd1020_p(m) = gdd10_p(m)                               ! <-- END of YR 1
         end if                                                 ! <-- END of YR 1
         gdd020_p(m)  = (yravgm1* gdd020_p(m)  + gdd0_p(m))  / yravg  ! gdd..20 must be long term avgs
         gdd820_p(m)  = (yravgm1* gdd820_p(m)  + gdd8_p(m))  / yravg  ! so ignore results for yrs 1 & 2
         gdd1020_p(m) = (yravgm1* gdd1020_p(m) + gdd10_p(m)) / yravg
         gdd0_p (m) = 0._r8
         gdd8_p (m) = 0._r8
         gdd10_p(m) = 0._r8
      end if
      
   end do

end subroutine CNPhenologyClimate


subroutine CNEvergreenPhenology (i,ps,pe,deltim,dayspyr)
     !
     ! !DESCRIPTION:
     ! For coupled carbon-nitrogen code (CN).
     !
     !
     ! !ARGUMENTS:
   integer ,intent(in) :: i
   integer ,intent(in) :: ps
   integer ,intent(in) :: pe
   real(r8),intent(in) :: deltim
   real(r8),intent(in) :: dayspyr                    ! Days per year
     !
     ! !LOCAL VARIABLES:

   real(r8):: tranr
   real(r8):: t1                         ! temporary variable
   integer :: ivt, m
     !-----------------------------------------------------------------------

      do m = ps , pe
         ivt = pftclass(m)
         if (isevg(ivt)) then
            bglfr_p(m) = 1._r8/(leaf_long(ivt) * dayspyr * 86400._r8)
            bgtr_p(m)  = 0._r8
            lgsf_p(m)  = 0._r8
         end if
      end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    if (CN_evergreen_phenology_opt == 1) then
!    do fp = 1,num_soilp
!       p = filter_soilp(fp)
!       print*,'isevg',isevg(ivt),ivt
     do m = ps , pe
        ivt = pftclass(m)
        if (isevg(ivt)) then

           tranr=0.0002_r8
          ! set carbon fluxes for shifting storage pools to transfer pools
!          if (use_matrixcn) then
!             matrix_phtransfer(p,ileafst_to_ileafxf_phc)   =  matrix_phtransfer(p,ileafst_to_ileafxf_phc) + tranr/dt
!             matrix_phtransfer(p,ifrootst_to_ifrootxf_phc) =  matrix_phtransfer(p,ifrootst_to_ifrootxf_phc) + tranr/dt
!             if (woody(ivt(p)) == 1.0_r8) then
!                matrix_phtransfer(p,ilivestemst_to_ilivestemxf_phc)   = matrix_phtransfer(p,ilivestemst_to_ilivestemxf_phc) + t     ranr/dt
!                matrix_phtransfer(p,ideadstemst_to_ideadstemxf_phc)   = matrix_phtransfer(p,ideadstemst_to_ideadstemxf_phc) + t     ranr/dt
!                matrix_phtransfer(p,ilivecrootst_to_ilivecrootxf_phc) = matrix_phtransfer(p,ilivecrootst_to_ilivecrootxf_phc) +      tranr/dt
!                matrix_phtransfer(p,ideadcrootst_to_ideadcrootxf_phc) = matrix_phtransfer(p,ideadcrootst_to_ideadcrootxf_phc) +      tranr/dt
!             end if
!          end if !use_matrixcn
           leafc_storage_to_xfer_p(m)  = tranr * leafc_storage_p(m)/deltim
           frootc_storage_to_xfer_p(m) = tranr * frootc_storage_p(m)/deltim
           if (woody(ivt) == 1) then
              livestemc_storage_to_xfer_p(m)  = tranr * livestemc_storage_p(m)/deltim
              deadstemc_storage_to_xfer_p(m)  = tranr * deadstemc_storage_p(m)/deltim
              livecrootc_storage_to_xfer_p(m) = tranr * livecrootc_storage_p(m)/deltim
              deadcrootc_storage_to_xfer_p(m) = tranr * deadcrootc_storage_p(m)/deltim
              gresp_storage_to_xfer_p(m)      = tranr * gresp_storage_p(m)/deltim
           end if

         ! set nitrogen fluxes for shifting storage pools to transfer pools
!         if (use_matrixcn) then
!            matrix_nphtransfer(p,ileafst_to_ileafxf_phn)   =  matrix_nphtransfer(p,ileafst_to_ileafxf_phn) + tranr/dt
!            matrix_nphtransfer(p,ifrootst_to_ifrootxf_phn) =  matrix_nphtransfer(p,ifrootst_to_ifrootxf_phn) + tranr/dt
!            if (woody(ivt(p)) == 1.0_r8) then
!               matrix_nphtransfer(p,ilivestemst_to_ilivestemxf_phn) = matrix_nphtransfer(p,ilivestemst_to_ilivestemxf_phn) + tr     anr/dt
!               matrix_nphtransfer(p,ideadstemst_to_ideadstemxf_phn) = matrix_nphtransfer(p,ideadstemst_to_ideadstemxf_phn) + tr     anr/dt
!               matrix_nphtransfer(p,ilivecrootst_to_ilivecrootxf_phn) = matrix_nphtransfer(p,ilivecrootst_to_ilivecrootxf_phn)      + tranr/dt
!               matrix_nphtransfer(p,ideadcrootst_to_ideadcrootxf_phn) = matrix_nphtransfer(p,ideadcrootst_to_ideadcrootxf_phn)      + tranr/dt
!            end if
!         end if !use_matrixcn
          leafn_storage_to_xfer_p(m)  = tranr * leafn_storage_p(m)/deltim
          frootn_storage_to_xfer_p(m) = tranr * frootn_storage_p(m)/deltim
          if (woody(ivt) == 1) then
             livestemn_storage_to_xfer_p(m)  = tranr * livestemn_storage_p(m)/deltim
             deadstemn_storage_to_xfer_p(m)  = tranr * deadstemn_storage_p(m)/deltim
             livecrootn_storage_to_xfer_p(m) = tranr * livecrootn_storage_p(m)/deltim
             deadcrootn_storage_to_xfer_p(m) = tranr * deadcrootn_storage_p(m)/deltim
          end if

          t1 = 1.0_r8 / deltim

!         if (use_matrixcn) then
!            matrix_phtransfer(p,ileafxf_to_ileaf_phc)   = matrix_phtransfer(p,ileafxf_to_ileaf_phc) + t1
!            matrix_phtransfer(p,ifrootxf_to_ifroot_phc) = matrix_phtransfer(p,ifrootxf_to_ifroot_phc) + t1
!
!            matrix_nphtransfer(p,ileafxf_to_ileaf_phn)   = matrix_nphtransfer(p,ileafxf_to_ileaf_phn) + t1
!            matrix_nphtransfer(p,ifrootxf_to_ifroot_phn) = matrix_nphtransfer(p,ifrootxf_to_ifroot_phn) + t1
!            if (woody(ivt(p)) == 1.0_r8) then
!               matrix_phtransfer(p,ilivestemxf_to_ilivestem_phc)   = matrix_phtransfer(p,ilivestemxf_to_ilivestem_phc) + t1
!               matrix_phtransfer(p,ideadstemxf_to_ideadstem_phc)   = matrix_phtransfer(p,ideadstemxf_to_ideadstem_phc) + t1
!               matrix_phtransfer(p,ilivecrootxf_to_ilivecroot_phc) = matrix_phtransfer(p,ilivecrootxf_to_ilivecroot_phc) + t1
!               matrix_phtransfer(p,ideadcrootxf_to_ideadcroot_phc) = matrix_phtransfer(p,ideadcrootxf_to_ideadcroot_phc) + t1
!
!               matrix_nphtransfer(p,ilivestemxf_to_ilivestem_phn)   = matrix_nphtransfer(p,ilivestemxf_to_ilivestem_phn) + t1
!               matrix_nphtransfer(p,ideadstemxf_to_ideadstem_phn)   = matrix_nphtransfer(p,ideadstemxf_to_ideadstem_phn) + t1
!               matrix_nphtransfer(p,ilivecrootxf_to_ilivecroot_phn) = matrix_nphtransfer(p,ilivecrootxf_to_ilivecroot_phn) + t1
!               matrix_nphtransfer(p,ideadcrootxf_to_ideadcroot_phn) = matrix_nphtransfer(p,ideadcrootxf_to_ideadcroot_phn) + t1
!            end if
!         end if !use_matrixcn

          leafc_xfer_to_leafc_p(m)   = t1 * leafc_xfer_p(m)
          frootc_xfer_to_frootc_p(m) = t1 * frootc_xfer_p(m)

          leafn_xfer_to_leafn_p(m)   = t1 * leafn_xfer_p(m)
          frootn_xfer_to_frootn_p(m) = t1 * frootn_xfer_p(m)
          if (woody(ivt) == 1) then
             livestemc_xfer_to_livestemc_p(m)   = t1 * livestemc_xfer_p(m)
             deadstemc_xfer_to_deadstemc_p(m)   = t1 * deadstemc_xfer_p(m)
             livecrootc_xfer_to_livecrootc_p(m) = t1 * livecrootc_xfer_p(m)
             deadcrootc_xfer_to_deadcrootc_p(m) = t1 * deadcrootc_xfer_p(m)

             livestemn_xfer_to_livestemn_p(m)   = t1 * livestemn_xfer_p(m)
             deadstemn_xfer_to_deadstemn_p(m)   = t1 * deadstemn_xfer_p(m)
             livecrootn_xfer_to_livecrootn_p(m) = t1 * livecrootn_xfer_p(m)
             deadcrootn_xfer_to_deadcrootn_p(m) = t1 * deadcrootn_xfer_p(m)
          end if

       end if ! end of if (isevg(ivt(p)) == 1._r8) then

    end do ! end of pft loop

!    end if ! end of if (CN_evergreen_phenology_opt == 1) then
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end subroutine CNEvergreenPhenology

   subroutine CNSeasonDecidPhenology(i,ps,pe,idate,deltim,dayspyr,dlat)

    !
    ! !DESCRIPTION:
    ! For coupled carbon-nitrogen code (CN).
    ! This routine handles the seasonal deciduous phenology code (temperate
    ! deciduous vegetation that has only one growing season per year).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer ,intent(in) :: i
    integer ,intent(in) :: ps
    integer ,intent(in) :: pe
    integer ,intent(in) :: idate(3)
    real(r8),intent(in) :: deltim
    real(r8),intent(in) :: dayspyr
    real(r8),intent(in) :: dlat

    !
    ! !LOCAL VARIABLES:
    real(r8):: ws_flag        !winter-summer solstice flag (0 or 1)
    real(r8):: crit_onset_gdd !critical onset growing degree-day sum
    real(r8):: soilt
    integer :: idate2_last
    integer :: ivt, m
    !-----------------------------------------------------------------------

      ! start patch loop

      idate2_last = idate(2) - 1
      if(idate2_last .le. 0)idate2_last=idate2_last+365
      prev_dayl(i)=daylength(dlat,idate2_last)
      dayl(i)     =daylength(dlat,idate(2))
!      print*,'issed',issed(ivt),ivt

      do m = ps , pe
         ivt = pftclass(m)
         if (issed(ivt)) then

         ! set background litterfall rate, background transfer rate, and
         ! long growing season factor to 0 for seasonal deciduous types
            bglfr_p(m) = 0._r8
            bgtr_p(m) = 0._r8
            lgsf_p(m) = 0._r8

         ! onset gdd sum from Biome-BGC, v4.1.2
            crit_onset_gdd = exp(4.8_r8 + 0.13_r8*(annavg_tref_p(m) - 273.15_r8))

         ! set flag for solstice period (winter->summer = 1, summer->winter = 0)
            if (dayl(i) >= prev_dayl(i)) then
               ws_flag = 1._r8
            else
               ws_flag = 0._r8
            end if

         ! update offset_counter and test for the end of the offset period
            if (offset_flag_p(m) == 1.0_r8) then
            ! decrement counter for offset period
               offset_counter_p(m) = offset_counter_p(m) - deltim

            ! if this is the end of the offset_period, reset phenology
            ! flags and indices
               if (offset_counter_p(m) == 0.0_r8) then
               ! this code block was originally handled by call cn_offset_cleanup(i)
               ! inlined during vectorization

                  offset_flag_p(m) = 0._r8
                  offset_counter_p(m) = 0._r8
                  dormant_flag_p(m) = 1._r8
                  days_active_p(m) = 0._r8
!               if (use_cndv) then
!                  pftmayexist(i) = .true.
!               end if

               ! reset the previous timestep litterfall flux memory
                  prev_leafc_to_litter_p(m) = 0._r8
                  prev_frootc_to_litter_p(m) = 0._r8
               end if
            end if

            ! update onset_counter and test for the end of the onset period
            if (onset_flag_p(m) == 1.0_r8) then
            ! decrement counter for onset period
               onset_counter_p(m) = onset_counter_p(m) - deltim

            ! if this is the end of the onset period, reset phenology
            ! flags and indices
               if (onset_counter_p(m) == 0.0_r8) then
               ! this code block was originally handled by call cn_onset_cleanup(i)
               ! inlined during vectorization

                  onset_flag_p(m) = 0.0_r8
                  onset_counter_p(m) = 0.0_r8
               ! set all transfer growth rates to 0.0
                  leafc_xfer_to_leafc_p(m)   = 0.0_r8
                  frootc_xfer_to_frootc_p(m) = 0.0_r8
                  leafn_xfer_to_leafn_p(m)   = 0.0_r8
                  frootn_xfer_to_frootn_p(m) = 0.0_r8
                  if (woody(ivt) == 1) then
                     livestemc_xfer_to_livestemc_p(m)   = 0.0_r8
                     deadstemc_xfer_to_deadstemc_p(m)   = 0.0_r8
                     livecrootc_xfer_to_livecrootc_p(m) = 0.0_r8
                     deadcrootc_xfer_to_deadcrootc_p(m) = 0.0_r8
                     livestemn_xfer_to_livestemn_p(m)   = 0.0_r8
                     deadstemn_xfer_to_deadstemn_p(m)   = 0.0_r8
                     livecrootn_xfer_to_livecrootn_p(m) = 0.0_r8
                     deadcrootn_xfer_to_deadcrootn_p(m) = 0.0_r8
                  end if
               ! set transfer pools to 0.0
                  leafc_xfer_p(m) = 0.0_r8
                  leafn_xfer_p(m) = 0.0_r8
                  frootc_xfer_p(m) = 0.0_r8
                  frootn_xfer_p(m) = 0.0_r8
                  if (woody(ivt) == 1) then
                     livestemc_xfer_p(m) = 0.0_r8
                     livestemn_xfer_p(m) = 0.0_r8
                     deadstemc_xfer_p(m) = 0.0_r8
                     deadstemn_xfer_p(m) = 0.0_r8
                     livecrootc_xfer_p(m) = 0.0_r8
                     livecrootn_xfer_p(m) = 0.0_r8
                     deadcrootc_xfer_p(m) = 0.0_r8
                     deadcrootn_xfer_p(m) = 0.0_r8
                  end if
               end if
            end if

         ! test for switching from dormant period to growth period
            if (dormant_flag_p(m) == 1.0_r8) then

            ! Test to turn on growing degree-day sum, if off.
            ! switch on the growing degree day sum on the winter solstice

               if (onset_gddflag_p(m) == 0._r8 .and. ws_flag == 1._r8) then
                  onset_gddflag_p(m) = 1._r8
                  onset_gdd_p(m) = 0._r8
               end if

            ! Test to turn off growing degree-day sum, if on.
            ! This test resets the growing degree day sum if it gets past
            ! the summer solstice without reaching the threshold value.
            ! In that case, it will take until the next winter solstice
            ! before the growing degree-day summation starts again.

               if (onset_gddflag_p(m) == 1._r8 .and. ws_flag == 0._r8) then
                  onset_gddflag_p(m) = 0._r8
                  onset_gdd_p(m) = 0._r8
               end if

            ! if the gdd flag is set, and if the soil is above freezing
            ! then accumulate growing degree days for onset trigger

               soilt = t_soisno(3,i)
               if (onset_gddflag_p(m) == 1.0_r8 .and. soilt > 273.15_r8) then
                  onset_gdd_p(m) = onset_gdd_p(m) + (soilt-273.15_r8)*(deltim/86400._r8)
               end if
   
            ! set onset_flag if critical growing degree-day sum is exceeded
               if (onset_gdd_p(m) > crit_onset_gdd) then
                  onset_flag_p(m) = 1.0_r8
                  dormant_flag_p(m) = 0.0_r8
                  onset_gddflag_p(m) = 0.0_r8
                  onset_gdd_p(m) = 0.0_r8
                  onset_counter_p(m) = ndays_on * 86400._r8

               ! move all the storage pools into transfer pools,
               ! where they will be transfered to displayed growth over the onset period.
               ! this code was originally handled with call cn_storage_to_xfer(p)
               ! inlined during vectorization

                  ! set carbon fluxes for shifting storage pools to transfer pools
!                  if(use_matrixcn)then
!                     matrix_phtransfer(p,ileafst_to_ileafxf_phc)   = matrix_phtransfer(p,ileafst_to_ileafxf_phc) + fstor2tran/dt
!                     matrix_phtransfer(p,ifrootst_to_ifrootxf_phc) = matrix_phtransfer(p,ifrootst_to_ifrootxf_phc) + fstor2tran/dt
!
!                     matrix_nphtransfer(p,ileafst_to_ileafxf_phn)   = matrix_nphtransfer(p,ileafst_to_ileafxf_phn) + fstor2tran/dt
!                     matrix_nphtransfer(p,ifrootst_to_ifrootxf_phn) = matrix_nphtransfer(p,ifrootst_to_ifrootxf_phn) + fstor2tran/dt
!                     if (woody(ivt) == 1) then
!                          matrix_phtransfer(p,ilivestemst_to_ilivestemxf_phc) = matrix_phtransfer(p,ilivestemst_to_ilivestemxf_phc) + fstor2tran/dt
!                          matrix_phtransfer(p,ideadstemst_to_ideadstemxf_phc) = matrix_phtransfer(p,ideadstemst_to_ideadstemxf_phc) + fstor2tran/dt
!                          matrix_phtransfer(p,ilivecrootst_to_ilivecrootxf_phc) = matrix_phtransfer(p,ilivecrootst_to_ilivecrootxf_phc) + fstor2tran/dt
!                          matrix_phtransfer(p,ideadcrootst_to_ideadcrootxf_phc) = matrix_phtransfer(p,ideadcrootst_to_ideadcrootxf_phc) + fstor2tran/dt
!                          matrix_nphtransfer(p,ilivestemst_to_ilivestemxf_phn) = matrix_nphtransfer(p,ilivestemst_to_ilivestemxf_phn) + fstor2tran/dt
!                          matrix_nphtransfer(p,ideadstemst_to_ideadstemxf_phn) = matrix_nphtransfer(p,ideadstemst_to_ideadstemxf_phn) + fstor2tran/dt
!                          matrix_nphtransfer(p,ilivecrootst_to_ilivecrootxf_phn) = matrix_nphtransfer(p,ilivecrootst_to_ilivecrootxf_phn) + fstor2tran/dt
!                          matrix_nphtransfer(p,ideadcrootst_to_ideadcrootxf_phn) = matrix_nphtransfer(p,ideadcrootst_to_ideadcrootxf_phn) + fstor2tran/dt
!                     end if
!                  else
                     ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
                     !                                        and CNNStateUpdate1::NStateUpdate1
!                  end if  ! use_matrixcn
                  leafc_storage_to_xfer_p(m)  = fstor2tran * leafc_storage_p(m)/deltim
                  frootc_storage_to_xfer_p(m) = fstor2tran * frootc_storage_p(m)/deltim
                  if (woody(ivt) == 1) then
                     livestemc_storage_to_xfer_p(m)  = fstor2tran * livestemc_storage_p(m)/deltim
                     deadstemc_storage_to_xfer_p(m)  = fstor2tran * deadstemc_storage_p(m)/deltim
                     livecrootc_storage_to_xfer_p(m) = fstor2tran * livecrootc_storage_p(m)/deltim
                     deadcrootc_storage_to_xfer_p(m) = fstor2tran * deadcrootc_storage_p(m)/deltim
                     gresp_storage_to_xfer_p(m)      = fstor2tran * gresp_storage_p(m)/deltim
                  end if

               ! set nitrogen fluxes for shifting storage pools to transfer pools
                  leafn_storage_to_xfer_p(m)  = fstor2tran * leafn_storage_p(m)/deltim
                  frootn_storage_to_xfer_p(m) = fstor2tran * frootn_storage_p(m)/deltim
                  if (woody(ivt) == 1) then
                     livestemn_storage_to_xfer_p(m)  = fstor2tran * livestemn_storage_p(m)/deltim
                     deadstemn_storage_to_xfer_p(m)  = fstor2tran * deadstemn_storage_p(m)/deltim
                     livecrootn_storage_to_xfer_p(m) = fstor2tran * livecrootn_storage_p(m)/deltim
                     deadcrootn_storage_to_xfer_p(m) = fstor2tran * deadcrootn_storage_p(m)/deltim
                  end if
               end if

            ! test for switching from growth period to offset period
            else if (offset_flag_p(m) == 0.0_r8) then
!             if (use_cndv) then
               ! If days_active > 355, then remove patch in
               ! CNDVEstablishment at the end of the year.
               ! days_active > 355 is a symptom of seasonal decid. patches occurring in
               ! gridcells where dayl never drops below crit_dayl.
               ! This results in TLAI>1e4 in a few gridcells.
!                days_active(p) = days_active(p) + fracday
!                if (days_active(p) > 355._r8) pftmayexist(p) = .false.
!             end if

                  ! only begin to test for offset daylength once past the summer sol
               if (ws_flag == 0._r8 .and. dayl(i) < crit_dayl) then
                  offset_flag_p(m) = 1._r8
                  offset_counter_p(m) = ndays_off * 86400._r8
                  prev_leafc_to_litter_p(m) = 0._r8
                  prev_frootc_to_litter_p(m) = 0._r8
               end if
            end if

         end if ! end if seasonal deciduous
      end do
   
   end subroutine CNSeasonDecidPhenology
   
  subroutine CNStressDecidPhenology(i,ps,pe,deltim,dayspyr)
   integer, intent(in) :: i
   integer ,intent(in) :: ps
   integer ,intent(in) :: pe
   real(r8),intent(in) :: deltim
   real(r8),intent(in) :: dayspyr         ! days per year

    ! !LOCAL VARIABLES:
   real(r8),parameter :: secspqtrday = 86400._r8 / 4  ! seconds per quarter day
   real(r8):: crit_onset_gdd  ! degree days for onset trigger
   real(r8):: soilt           ! temperature of top soil layer
   real(r8):: psi             ! soil water potential [MPa]
   real(r8):: rain_threshold  ! rain threshold for leaf on [mm]
   logical :: additional_onset_condition ! additional condition for leaf onset
   integer :: ivt, m
   !-----------------------------------------------------------------------


   ! specify rain threshold for leaf onset
   rain_threshold = 20._r8

   do m = ps , pe
      ivt = pftclass(m)
!   print*,'isstd',isstd(ivt),ivt
      if (isstd(ivt)) then
         soilt = t_soisno(3,i)
         psi = smp(3,i) * 1.e-5 ! mmH2O -> MPa

      ! onset gdd sum from Biome-BGC, v4.1.2
         crit_onset_gdd = exp(4.8_r8 + 0.13_r8*(annavg_tref_p(m) - 273.15_r8))

!      print*,'offset_flag',offset_flag_p(m)
!      print*,'onset_flag',onset_flag_p(m)
!      print*,'dormant_flag',dormant_flag_p(m)
!      print*,'annavg_tref',annavg_tref_p(m)
      ! update offset_counter and test for the end of the offset period
         if (offset_flag_p(m) == 1._r8) then
         ! decrement counter for offset period
            offset_counter_p(m) = offset_counter_p(m) - deltim

         ! if this is the end of the offset_period, reset phenology
         ! flags and indices
            if (offset_counter_p(m) == 0._r8) then
            ! this code block was originally handled by call cn_offset_cleanup(i)
            ! inlined during vectorization
               offset_flag_p(m) = 0._r8
               offset_counter_p(m) = 0._r8
               dormant_flag_p(m) = 1._r8
               days_active_p(m) = 0._r8

            ! reset the previous timestep litterfall flux memory
               prev_leafc_to_litter_p(m) = 0._r8
               prev_frootc_to_litter_p(m) = 0._r8
            end if
         end if

      ! update onset_counter and test for the end of the onset period
         if (onset_flag_p(m) == 1.0_r8) then
         ! decrement counter for onset period
            onset_counter_p(m) = onset_counter_p(m) - deltim
!         print*,'onset_flag_counting',onset_counter_p(m), deltim

         ! if this is the end of the onset period, reset phenology
         ! flags and indices
            if (onset_counter_p(m) == 0.0_r8) then
            ! this code block was originally handled by call cn_onset_cleanup(i)
            ! inlined during vectorization
               onset_flag_p(m) = 0._r8
               onset_counter_p(m) = 0._r8
            ! set all transfer growth rates to 0.0
               leafc_xfer_to_leafc_p(m)   = 0._r8
               frootc_xfer_to_frootc_p(m) = 0._r8
               leafn_xfer_to_leafn_p(m)   = 0._r8
               frootn_xfer_to_frootn_p(m) = 0._r8
               if (woody(ivt) == 1) then
                  livestemc_xfer_to_livestemc_p(m)   = 0._r8
                  deadstemc_xfer_to_deadstemc_p(m)   = 0._r8
                  livecrootc_xfer_to_livecrootc_p(m) = 0._r8
                  deadcrootc_xfer_to_deadcrootc_p(m) = 0._r8
                  livestemn_xfer_to_livestemn_p(m)   = 0._r8
                  deadstemn_xfer_to_deadstemn_p(m)   = 0._r8
                  livecrootn_xfer_to_livecrootn_p(m) = 0._r8
                  deadcrootn_xfer_to_deadcrootn_p(m) = 0._r8
               end if
               ! set transfer pools to 0.0
               leafc_xfer_p(m) = 0._r8
               leafn_xfer_p(m) = 0._r8
               frootc_xfer_p(m) = 0._r8
               frootn_xfer_p(m) = 0._r8
               if (woody(ivt) == 1) then
                  livestemc_xfer_p(m) = 0._r8
                  livestemn_xfer_p(m) = 0._r8
                  deadstemc_xfer_p(m) = 0._r8
                  deadstemn_xfer_p(m) = 0._r8
                  livecrootc_xfer_p(m) = 0._r8
                  livecrootn_xfer_p(m) = 0._r8
                  deadcrootc_xfer_p(m) = 0._r8
                  deadcrootn_xfer_p(m) = 0._r8
               end if
            end if
         end if

      ! test for switching from dormant period to growth period
         if (dormant_flag_p(m) == 1._r8) then

         ! keep track of the number of freezing degree days in this
         ! dormancy period (only if the freeze flag has not previously been set
         ! for this dormancy period

            if (onset_gddflag_p(m) == 0._r8 .and. soilt < 273.15_r8) onset_fdd_p(m) = onset_fdd_p(m) + deltim/86400._r8

         ! if the number of freezing degree days exceeds a critical value,
         ! then onset will require both wet soils and a critical soil
         ! temperature sum.  If this case is triggered, reset any previously
         ! accumulated value in onset_swi, so that onset now depends on
         ! the accumulated soil water index following the freeze trigger

            if (onset_fdd_p(m) > crit_onset_fdd) then
               onset_gddflag_p(m) = 1._r8
               onset_fdd_p(m) = 0._r8
               onset_swi_p(m) = 0._r8
            end if
   
         ! if the freeze flag is set, and if the soil is above freezing
         ! then accumulate growing degree days for onset trigger

            if (onset_gddflag_p(m) == 1._r8 .and. soilt > 273.15_r8) then
               onset_gdd_p(m) = onset_gdd_p(m) + (soilt-273.15_r8)*deltim/86400._r8
            end if

        ! if soils are wet, accumulate soil water index for onset trigger
            additional_onset_condition = .true.
!         if(CNParamsShareInst%constrain_stress_deciduous_onset) then
            ! if additional constraint condition not met,  set to false
               if ((prec10(i) * (3600.0_r8*10.0_r8*24.0_r8)) < rain_threshold) then
                  additional_onset_condition = .false.
               endif
!         endif

            if (psi >= soilpsi_on) then
               onset_swi_p(m) = onset_swi_p(m) + deltim/86400._r8
            endif

         ! if critical soil water index is exceeded, set onset_flag, and
         ! then test for soil temperature criteria

         ! Adding in Kyla's rainfall trigger when fun on. RF. prec10 (mm/s) needs to be higher than 8mm over 10 days.

!         print*,'onset_swi,crit_onset_swi',onset_swi_p(m),crit_onset_swi,&
!additional_onset_condition,psi,soilpsi_on,.true.,onset_gdd_p(m),crit_onset_gdd
            if (onset_swi_p(m) > crit_onset_swi.and. additional_onset_condition)  then
               onset_flag_p(m) = 1._r8

            ! only check soil temperature criteria if freeze flag set since
            ! beginning of last dormancy.  If freeze flag set and growing
            ! degree day sum (since freeze trigger) is lower than critical
            ! value, then override the onset_flag set from soil water.

               if (onset_gddflag_p(m) == 1._r8 .and. onset_gdd_p(m) < crit_onset_gdd) onset_flag_p(m) = 0._r8
            end if

         ! only allow onset if dayl > 6hrs
            if (onset_flag_p(m) == 1._r8 .and. dayl(i) <= secspqtrday) then
               onset_flag_p(m) = 0._r8
            end if

         ! if this is the beginning of the onset period
         ! then reset the phenology flags and indices
!         print*,'onset_flag2',onset_flag_p(m)

            if (onset_flag_p(m) == 1._r8) then
               dormant_flag_p(m) = 0._r8
               days_active_p(m) = 0._r8
               onset_gddflag_p(m) = 0._r8
               onset_fdd_p(m) = 0._r8
               onset_gdd_p(m) = 0._r8
               onset_swi_p(m) = 0._r8
               onset_counter_p(m) = ndays_on * 86400._r8

            ! call subroutine to move all the storage pools into transfer pools,
            ! where they will be transfered to displayed growth over the onset period.
            ! this code was originally handled with call cn_storage_to_xfer(i)
            ! inlined during vectorization

            ! set carbon fluxes for shifting storage pools to transfer pools
!            if (use_matrixcn) then
!               matrix_phtransfer(p,ileafst_to_ileafxf_phc)   = matrix_phtransfer(p,ileafst_to_ileafxf_phc) + fstor2tran/dt
!               matrix_phtransfer(p,ifrootst_to_ifrootxf_phc) = matrix_phtransfer(p,ifrootst_to_ifrootxf_phc) + fstor2tran/dt
!               matrix_nphtransfer(p,ileafst_to_ileafxf_phn)   = matrix_nphtransfer(p,ileafst_to_ileafxf_phn) + fstor2tran/dt
!               matrix_nphtransfer(p,ifrootst_to_ifrootxf_phn) = matrix_nphtransfer(p,ifrootst_to_ifrootxf_phn) + fstor2tran/dt
!               if (woody(ivt(p)) == 1.0_r8) then
!                  matrix_phtransfer(p,ilivestemst_to_ilivestemxf_phc) = matrix_phtransfer(p,ilivestemst_to_ilivestemxf_phc) + fstor2tran/dt
!                  matrix_phtransfer(p,ideadstemst_to_ideadstemxf_phc) = matrix_phtransfer(p,ideadstemst_to_ideadstemxf_phc) + fstor2tran/dt
!                  matrix_phtransfer(p,ilivecrootst_to_ilivecrootxf_phc) = matrix_phtransfer(p,ilivecrootst_to_ilivecrootxf_phc) + fstor2tran/dt
!                  matrix_phtransfer(p,ideadcrootst_to_ideadcrootxf_phc) = matrix_phtransfer(p,ideadcrootst_to_ideadcrootxf_phc) + fstor2tran/dt
!
!                  matrix_nphtransfer(p,ilivestemst_to_ilivestemxf_phn) = matrix_nphtransfer(p,ilivestemst_to_ilivestemxf_phn) + fstor2tran/dt
!                  matrix_nphtransfer(p,ideadstemst_to_ideadstemxf_phn) = matrix_nphtransfer(p,ideadstemst_to_ideadstemxf_phn) + fstor2tran/dt
!                  matrix_nphtransfer(p,ilivecrootst_to_ilivecrootxf_phn) = matrix_nphtransfer(p,ilivecrootst_to_ilivecrootxf_phn) + fstor2tran/dt
!                  matrix_nphtransfer(p,ideadcrootst_to_ideadcrootxf_phn) = matrix_nphtransfer(p,ideadcrootst_to_ideadcrootxf_phn) + fstor2tran/dt
!               end if
!            else
!               ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
!               !                                        and CNNStateUpdate1::NStateUpdate1
!            end if
!            print*,'leafc_storage_to_xfer in CNPHenology1',leafc_storage_to_xfer_p(m),fstor2tran,leafc_storage_p(m),deltim
!            print*,'frootc_storage_to_xfer in CNPHenology1',frootc_storage_to_xfer(i),fstor2tran,frootc_storage_p(m),deltim
               leafc_storage_to_xfer_p(m)  = fstor2tran * leafc_storage_p(m)/deltim
               frootc_storage_to_xfer_p(m) = fstor2tran * frootc_storage_p(m)/deltim
               if (woody(ivt) == 1) then
                  livestemc_storage_to_xfer_p(m)  = fstor2tran * livestemc_storage_p(m)/deltim
                  deadstemc_storage_to_xfer_p(m)  = fstor2tran * deadstemc_storage_p(m)/deltim
                  livecrootc_storage_to_xfer_p(m) = fstor2tran * livecrootc_storage_p(m)/deltim
                  deadcrootc_storage_to_xfer_p(m) = fstor2tran * deadcrootc_storage_p(m)/deltim
                  gresp_storage_to_xfer_p(m)      = fstor2tran * gresp_storage_p(m)/deltim
               end if
   
            ! set nitrogen fluxes for shifting storage pools to transfer pools
               leafn_storage_to_xfer_p(m)  = fstor2tran * leafn_storage_p(m)/deltim
               frootn_storage_to_xfer_p(m) = fstor2tran * frootn_storage_p(m)/deltim
               if (woody(ivt) == 1) then
                  livestemn_storage_to_xfer_p(m)  = fstor2tran * livestemn_storage_p(m)/deltim
                  deadstemn_storage_to_xfer_p(m)  = fstor2tran * deadstemn_storage_p(m)/deltim
                  livecrootn_storage_to_xfer_p(m) = fstor2tran * livecrootn_storage_p(m)/deltim
                  deadcrootn_storage_to_xfer_p(m) = fstor2tran * deadcrootn_storage_p(m)/deltim
               end if
            end if
   
         ! test for switching from growth period to offset period
         else if (offset_flag_p(m) == 0._r8) then

         ! if soil water potential lower than critical value, accumulate
         ! as stress in offset soil water index

            if (psi <= soilpsi_off) then
               offset_swi_p(m) = offset_swi_p(m) + deltim/86400._r8

            ! if the offset soil water index exceeds critical value, and
            ! if this is not the middle of a previously initiated onset period,
            ! then set flag to start the offset period and reset index variables

               if (offset_swi_p(m) >= crit_offset_swi .and. onset_flag_p(m) == 0._r8) offset_flag_p(m) = 1._r8

            ! if soil water potential higher than critical value, reduce the
            ! offset water stress index.  By this mechanism, there must be a
            ! sustained period of water stress to initiate offset.

            else if (psi >= soilpsi_on) then
               offset_swi_p(m) = offset_swi_p(m) - deltim/86400._r8
               offset_swi_p(m) = max(offset_swi_p(m),0._r8)
            end if

         ! decrease freezing day accumulator for warm soil
            if (offset_fdd_p(m) > 0._r8 .and. soilt > 273.15_r8) then
               offset_fdd_p(m) = offset_fdd_p(m) - deltim/86400._r8
               offset_fdd_p(m) = max(0._r8, offset_fdd_p(m))
            end if

         ! increase freezing day accumulator for cold soil
            if (soilt <= 273.15_r8) then
               offset_fdd_p(m) = offset_fdd_p(m) + deltim/86400._r8

            ! if freezing degree day sum is greater than critical value, initiate offset
               if (offset_fdd_p(m) > crit_offset_fdd .and. onset_flag_p(m) == 0._r8) offset_flag_p(m) = 1._r8
            end if

         ! force offset if daylength is < 6 hrs
            if (dayl(i) <= secspqtrday) then
               offset_flag_p(m) = 1._r8
            end if

         ! if this is the beginning of the offset period
         ! then reset flags and indices
            if (offset_flag_p(m) == 1._r8) then
               offset_fdd_p(m) = 0._r8
               offset_swi_p(m) = 0._r8
               offset_counter_p(m) = ndays_off * 86400._r8
               prev_leafc_to_litter_p(m) = 0._r8
               prev_frootc_to_litter_p(m) = 0._r8
            end if
         end if

      ! keep track of number of days since last dormancy for control on
      ! fraction of new growth to send to storage for next growing season

         if (dormant_flag_p(m) == 0.0_r8) then
            days_active_p(m) = days_active_p(m) + deltim/86400._r8
         end if

      ! calculate long growing season factor (lgsf)
      ! only begin to calculate a lgsf greater than 0.0 once the number
      ! of days active exceeds days/year.
         lgsf_p(m) = max(min(3.0_r8*(days_active_p(m)-leaf_long(ivt)*dayspyr )/dayspyr, 1._r8),0._r8)
      !print*,'lgsf',lgsf_p(m),days_active_p(m),leaf_long(ivt),dayspyr
      ! RosieF. 5 Nov 2015.  Changed this such that the increase in leaf turnover is faster after
      ! trees enter the 'fake evergreen' state. Otherwise, they have a whole year of
      ! cheating, with less litterfall than they should have, resulting in very high LAI.
      ! Further, the 'fake evergreen' state (where lgsf>0) is entered at the end of a single leaf lifespan
      ! and not a whole year. The '3' is arbitrary, given that this entire system is quite abstract.
      ! set background litterfall rate, when not in the phenological offset period
         if (offset_flag_p(m) == 1._r8) then
            bglfr_p(m) = 0._r8
         else
         ! calculate the background litterfall rate (bglfr)
         ! in units 1/s, based on leaf longevity (yrs) and correction for long growing season

            bglfr_p(m) = (1._r8/(leaf_long(ivt)*dayspyr*86400._r8))*lgsf_p(m)
         end if

      ! set background transfer rate when active but not in the phenological onset period
         if (onset_flag_p(m) == 1._r8) then
            bgtr_p(m) = 0._r8
         else
         ! the background transfer rate is calculated as the rate that would result
         ! in complete turnover of the storage pools in one year at steady state,
         ! once lgsf has reached 1.0 (after 730 days active).

            bgtr_p(m) = (1._r8/(dayspyr*86400._r8))*lgsf_p(m)

         ! set carbon fluxes for shifting storage pools to transfer pools

         ! reduced the amount of stored carbon flowing to display pool by only counting the delta
         ! between leafc and leafc_store in the flux. RosieF, Nov5 2015.
            leafc_storage_to_xfer_p(m)  = max(0.0_r8,(leafc_storage_p(m)-leafc_p(m))) * bgtr_p(m)
            frootc_storage_to_xfer_p(m) = max(0.0_r8,(frootc_storage_p(m)-frootc_p(m))) * bgtr_p(m)
         !print*,'frootc_storage_to_xfer in CNPhenology2',leafc_storage_to_xfer_p(m),frootc_storage_to_xfer_p(m),&
!         frootc_storage_p(m),bgtr_p(m)
!         if (use_matrixcn) then
!            if(leafc_storage(p) .gt. 0)then
!               matrix_phtransfer(p,ileafst_to_ileafxf_phc)   = matrix_phtransfer(p,ileafst_to_ileafxf_phc) &
!                                                      + leafc_storage_to_xfer(p) / leafc_storage(p)
!            else
!              leafc_storage_to_xfer(p) = 0
!            end if
!            if(frootc_storage(p) .gt. 0)then
!               matrix_phtransfer(p,ifrootst_to_ifrootxf_phc) = matrix_phtransfer(p,ifrootst_to_ifrootxf_phc) &
!                                                      + frootc_storage_to_xfer(p) / frootc_storage(p)
!            else
!               frootc_storage_to_xfer(p) = 0
!            end if
!          if (woody(ivt(p)) == 1.0_r8) then
!              matrix_phtransfer(p,ilivestemst_to_ilivestemxf_phc) = matrix_phtransfer(p,ilivestemst_to_ilivestemxf_phc) + bgtr(p)
!              matrix_phtransfer(p,ideadstemst_to_ideadstemxf_phc) = matrix_phtransfer(p,ideadstemst_to_ideadstemxf_phc) + bgtr(p)
!              matrix_phtransfer(p,ilivecrootst_to_ilivecrootxf_phc) = matrix_phtransfer(p,ilivecrootst_to_ilivecrootxf_phc) + bgtr(p)
!              matrix_phtransfer(p,ideadcrootst_to_ideadcrootxf_phc) = matrix_phtransfer(p,ideadcrootst_to_ideadcrootxf_phc) + bgtr(p)
!           end if
!        else
!           ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
!           !                                        and CNNStateUpdate1::NStateUpdate1
!        end if !use_matrixcn
            if (woody(ivt) == 1) then
               livestemc_storage_to_xfer_p(m)  = livestemc_storage_p(m) * bgtr_p(m)
               deadstemc_storage_to_xfer_p(m)  = deadstemc_storage_p(m) * bgtr_p(m)
               livecrootc_storage_to_xfer_p(m) = livecrootc_storage_p(m) * bgtr_p(m)
               deadcrootc_storage_to_xfer_p(m) = deadcrootc_storage_p(m) * bgtr_p(m)
               gresp_storage_to_xfer_p(m)      = gresp_storage_p(m) * bgtr_p(m)
            end if

         ! set nitrogen fluxes for shifting storage pools to transfer pools
            leafn_storage_to_xfer_p(m)  = leafn_storage_p(m) * bgtr_p(m)
            frootn_storage_to_xfer_p(m) = frootn_storage_p(m) * bgtr_p(m)
!         if (use_matrixcn) then
!            matrix_nphtransfer(p,ileafst_to_ileafxf_phn)   = matrix_nphtransfer(p,ileafst_to_ileafxf_phn) &
!                                                      + bgtr(p)
!            matrix_nphtransfer(p,ifrootst_to_ifrootxf_phn) = matrix_nphtransfer(p,ifrootst_to_ifrootxf_phn) &
!                                                      + bgtr(p)
!            if (woody(ivt(p)) == 1.0_r8) then
!               matrix_nphtransfer(p,ilivestemst_to_ilivestemxf_phn) = matrix_nphtransfer(p,ilivestemst_to_ilivestemxf_phn) + bgtr(p)
!               matrix_nphtransfer(p,ideadstemst_to_ideadstemxf_phn) = matrix_nphtransfer(p,ideadstemst_to_ideadstemxf_phn) + bgtr(p)
!               matrix_nphtransfer(p,ilivecrootst_to_ilivecrootxf_phn) = matrix_nphtransfer(p,ilivecrootst_to_ilivecrootxf_phn) + bgtr(p)
!               matrix_nphtransfer(p,ideadcrootst_to_ideadcrootxf_phn) = matrix_nphtransfer(p,ideadcrootst_to_ideadcrootxf_phn) + bgtr(p)
!            end if
!         else
           ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
           !                                        and CNNStateUpdate1::NStateUpdate1
!         end if !use_matrixcn
            if (woody(ivt) == 1) then
               livestemn_storage_to_xfer_p(m)  = livestemn_storage_p(m) * bgtr_p(m)
               deadstemn_storage_to_xfer_p(m)  = deadstemn_storage_p(m) * bgtr_p(m)
               livecrootn_storage_to_xfer_p(m) = livecrootn_storage_p(m) * bgtr_p(m)
               deadcrootn_storage_to_xfer_p(m) = deadcrootn_storage_p(m) * bgtr_p(m)
            end if
         end if
   
      end if ! end if stress deciduous
   end do !end pft loop

  end subroutine CNStressDecidPhenology


  subroutine CropPhenology(i,ps,pe,idate,h,deltim,dayspyr,npcropmin)

    ! !DESCRIPTION:
    ! Code from AgroIBIS to determine crop phenology and code from CN to
    ! handle CN fluxes during the phenological onset                       & offset periods.
    
    ! !USES:
    !

   integer, intent(in) :: i
   integer ,intent(in) :: ps
   integer ,intent(in) :: pe
   integer ,intent(in) :: idate(1:3)
   integer ,intent(in) :: h   ! 1 for north hemisphere; 2 for south hemisphere
   real(r8),intent(in) :: deltim
   real(r8),intent(in) :: dayspyr         ! days per year
   integer ,intent(in) :: npcropmin

    ! LOCAL VARAIBLES:
    integer kyr       ! current year
    integer kmo       ! month of year  (1, ..., 12)
    integer kda       ! day of month   (1, ..., 31)
    integer mcsec 
    integer jday
    integer fp,m      ! patch indices
    integer c         ! column indices
    integer g         ! gridcell indices
    integer idpp      ! number of days past planting
    real(r8) crmcorn  ! comparitive relative maturity for corn
    real(r8) ndays_on ! number of days to fertilize
    integer :: jdayyrstart(2)
    real(r8) :: initial_seed_at_planting = 3._r8 ! Initial seed at planting
    integer ivt

    !------------------------------------------------------------------------

    jdayyrstart(1) = 1
    jdayyrstart(2) = 182

    jday  = idate(2)
    mcsec = idate(3)
      ! get time info

!      if (use_fertilizer) then
       ndays_on = 20._r8 ! number of days to fertilize
!      else
!       ndays_on = 0._r8 ! number of days to fertilize
!      end if

         ! background litterfall and transfer rates; long growing season factor
      do m = ps, pe
         ivt = pftclass(m)
         if(ivt >= npcropmin)then
          bglfr_p(m) = 0._r8 ! this value changes later in a crop's life cycle
          bgtr_p(m)  = 0._r8
          lgsf_p(m)  = 0._r8

         ! ---------------------------------
         ! from AgroIBIS subroutine planting
         ! ---------------------------------

         ! in order to allow a crop to be planted only once each year
         ! initialize cropplant = .false., but hold it = .true. through the end of the year

         ! initialize other variables that are calculated for crops
         ! on an annual basis in cropresidue subroutine

          if ( jday == jdayyrstart(h) .and. mcsec == 1800 )then

            ! make sure variables aren't changed at beginning of the year
            ! for a crop that is currently planted, such as
            ! WINTER TEMPERATE CEREAL = winter (wheat + barley + rye)
            ! represented here by the winter wheat pft

            if (.not. croplive_p(m))  then
               cropplant_p(m) = .false.
               idop_p(m)      = NOT_Planted

               ! keep next for continuous, annual winter temperate cereal crop;
               ! if we removed elseif,
               ! winter cereal grown continuously would amount to a cereal/fallow
               ! rotation because cereal would only be planted every other year

            else if (croplive_p(m) .and. (ivt == nwwheat .or. ivt == nirrig_wwheat)) then
               cropplant_p(m) = .false.
               !           else ! not possible to have croplive and ivt==cornORsoy? (slevis)
            end if

          end if
 
          if ( (.not. croplive_p(m)) .and. (.not. cropplant_p(m)) ) then

            ! gdd needed for * chosen crop and a likely hybrid (for that region) *
            ! to reach full physiological maturity

            ! based on accumulated seasonal average growing degree days from
            ! April 1 - Sept 30 (inclusive)
            ! for corn and soybeans in the United States -
            ! decided upon by what the typical average growing season length is
            ! and the gdd needed to reach maturity in those regions

            ! first choice is used for spring temperate cereal and/or soybeans and maize

            ! slevis: ibis reads xinpdate in io.f from control.crops.nc variable name 'plantdate'
            !         According to Chris Kucharik, the dataset of
            !         xinpdate was generated from a previous model run at 0.5 deg resolution

            ! winter temperate cereal : use gdd0 as a limit to plant winter cereal
             if (ivt == nwwheat .or. ivt == nirrig_wwheat) then

               ! add check to only plant winter cereal after other crops (soybean, maize)
               ! have been harvested

               ! *** remember order of planting is crucial - in terms of which crops you want
               ! to be grown in what order ***

               ! in this case, corn or soybeans are assumed to be planted before
               ! cereal would be in any particular year that both patches are allowed
               ! to grow in the same grid cell (e.g., double-cropping)

               ! slevis: harvdate below needs cropplant(p) above to be cropplant(p,ivt(p))
               !         where ivt(p) has rotated to winter cereal because
               !         cropplant through the end of the year for a harvested crop.
               !         Also harvdate(p) should be harvdate(p,ivt(p)) and should be
               !         updated on Jan 1st instead of at harvest (slevis)
               if (a5tmin_p(m)             /= spval                  .and. &
                    a5tmin_p(m)             <= minplanttemp(ivt)   .and. &
                    jday                  >= minplantjday(ivt,h) .and. &
                    (gdd020_p(m)            /= spval                  .and. &
                    gdd020_p(m)             >= gddmin(ivt))) then

                  cumvd_p(m)       = 0._r8
                  hdidx_p(m)       = 0._r8
                  vf_p(m)          = 0._r8
                  croplive_p(m)    = .true.
                  cropplant_p(m)   = .true.
                  idop_p(m)        = jday
                  harvdate_p(m)    = NOT_Harvested
                  gddmaturity_p(m) = hybgdd(ivt)
                  leafc_xfer_p(m)  = initial_seed_at_planting
                  leafn_xfer_p(m)  = leafc_xfer_p(m) / leafcn(ivt) ! with onset
                  crop_seedc_to_leaf_p(m) = leafc_xfer_p(m)/deltim
                  crop_seedn_to_leaf_p(m) = leafn_xfer_p(m)/deltim

                  ! because leafc_xfer is set above rather than incremneted through the normal process, must also set its isotope
                  ! pools here.  use totvegc_patch as the closest analogue if nonzero, and use initial value otherwise
!                  if (use_c13) then
!                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
!                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
!                             c13_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
!                     else
!                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c13ratio
!                     endif
!                  endif
!                  if (use_c14) then
!                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
!                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
!                             c14_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
!                     else
!                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c14ratio
!                     endif
!                  endif

                  ! latest possible date to plant winter cereal and after all other 
                  ! crops were harvested for that year

               else if (jday       >=  maxplantjday(ivt,h) .and. &
                    gdd020_p(m)  /= spval                   .and. &
                    gdd020_p(m)  >= gddmin(ivt)) then

                  cumvd_p(m)       = 0._r8
                  hdidx_p(m)       = 0._r8
                  vf_p(m)          = 0._r8
                  croplive_p(m)    = .true.
                  cropplant_p(m)   = .true.
                  idop_p(m)        = jday
                  harvdate_p(m)    = NOT_Harvested
                  gddmaturity_p(m) = hybgdd(ivt)
                  leafc_xfer_p(m)  = initial_seed_at_planting
                  leafn_xfer_p(m)  = leafc_xfer_p(m) / leafcn(ivt) ! with onset
                  crop_seedc_to_leaf_p(m) = leafc_xfer_p(m)/deltim
                  crop_seedn_to_leaf_p(m) = leafn_xfer_p(m)/deltim

                  ! because leafc_xfer is set above rather than incremneted through the normal process, must also set its isotope
                  ! pools here.  use totvegc_patch as the closest analogue if nonzero, and use initial value otherwise
!                  if (use_c13) then
!                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
!                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
!                             c13_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
!                     else
!                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c13ratio
!                     endif
!                  endif
!                  if (use_c14) then
!                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
!                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
!                             c14_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
!                     else
!                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c14ratio
!                     endif
!                  endif
               else
                  gddmaturity_p(m) = 0._r8
               end if

             else ! not winter cereal... slevis: added distinction between NH and SH
               ! slevis: The idea is that jday will equal idop sooner or later in the year
               !         while the gdd part is either true or false for the year.
               if (t10_p(m) /= spval.and. a10tmin_p(m) /= spval   .and. &
                    t10_p(m)     > planttemp(ivt)             .and. &
                    a10tmin_p(m) > minplanttemp(ivt)          .and. &
                    jday       >= minplantjday(ivt,h)       .and. &
                    jday       <= maxplantjday(ivt,h)       .and. &
                    t10_p(m) /= spval .and. a10tmin_p(m) /= spval  .and. &
                    gdd820_p(m) /= spval                         .and. &
                    gdd820_p(m) >= gddmin(ivt)) then

                  ! impose limit on growing season length needed
                  ! for crop maturity - for cold weather constraints
                  croplive_p(m)  = .true.
                  cropplant_p(m) = .true.
                  idop_p(m)      = jday
                  harvdate_p(m)  = NOT_Harvested

                  ! go a specified amount of time before/after
                  ! climatological date
                  if ( ivt == ntmp_soybean .or. ivt == nirrig_tmp_soybean .or. &
                       ivt == ntrp_soybean .or. ivt == nirrig_trp_soybean) then
                     gddmaturity_p(m) = min(gdd1020_p(m), hybgdd(ivt))
                  end if
                  
                  if (ivt == ntmp_corn .or. ivt == nirrig_tmp_corn .or. &
                      ivt == ntrp_corn .or. ivt == nirrig_trp_corn .or. &
                      ivt == nsugarcane .or. ivt == nirrig_sugarcane .or. &
                      ivt == nmiscanthus .or. ivt == nirrig_miscanthus .or. &
                      ivt == nswitchgrass .or. ivt == nirrig_switchgrass) then
                     gddmaturity_p(m) = max(950._r8, min(gdd820_p(m)*0.85_r8, hybgdd(ivt)))
                     gddmaturity_p(m) = max(950._r8, min(gddmaturity_p(m)+150._r8, 1850._r8))
                  end if
                  if (ivt == nswheat .or. ivt == nirrig_swheat .or. &
                      ivt == ncotton .or. ivt == nirrig_cotton .or. &
                      ivt == nrice   .or. ivt == nirrig_rice) then
                     gddmaturity_p(m) = min(gdd020_p(m), hybgdd(ivt))
                  end if

                  leafc_xfer_p(m)  = initial_seed_at_planting
                  leafn_xfer_p(m) = leafc_xfer_p(m) / leafcn(ivt) ! with onset
                  crop_seedc_to_leaf_p(m) = leafc_xfer_p(m)/deltim
                  crop_seedn_to_leaf_p(m) = leafn_xfer_p(m)/deltim

                  ! because leafc_xfer is set above rather than incremneted through the normal process, must also set its isotope
                  ! pools here.  use totvegc_patch as the closest analogue if nonzero, and use initial value otherwise
!                  if (use_c13) then
!                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
!                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
!                             c13_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
!                     else
!                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c13ratio
!                     endif
!                  endif
!                  if (use_c14) then
!                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
!                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
!                             c14_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
!                     else
!                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c14ratio
!                     endif
!                  endif


                  ! If hit the max planting julian day -- go ahead and plant
               else if (jday == maxplantjday(ivt,h) .and. gdd820_p(m) > 0._r8 .and. &
                    gdd820_p(m) /= spval ) then
                  croplive_p(m)  = .true.
                  cropplant_p(m) = .true.
                  idop_p(m)      = jday
                  harvdate_p(m)  = NOT_Harvested

                  if (ivt == ntmp_soybean .or. ivt == nirrig_tmp_soybean .or. &
                      ivt == ntrp_soybean .or. ivt == nirrig_trp_soybean) then
                     gddmaturity_p(m) = min(gdd1020_p(m), hybgdd(ivt))
                  end if
                  
                  if (ivt == ntmp_corn .or. ivt == nirrig_tmp_corn .or. &
                      ivt == ntrp_corn .or. ivt == nirrig_trp_corn .or. &
                      ivt == nsugarcane .or. ivt == nirrig_sugarcane .or. &
                      ivt == nmiscanthus .or. ivt == nirrig_miscanthus .or. &
                      ivt == nswitchgrass .or. ivt == nirrig_switchgrass) then
                     gddmaturity_p(m) = max(950._r8, min(gdd820_p(m)*0.85_r8, hybgdd(ivt)))
                  end if
                  if (ivt == nswheat .or. ivt == nirrig_swheat .or. &
                      ivt == ncotton .or. ivt == nirrig_cotton .or. &
                      ivt == nrice   .or. ivt == nirrig_rice) then
                     gddmaturity_p(m) = min(gdd020_p(m), hybgdd(ivt))
                  end if

                  leafc_xfer_p(m)  = initial_seed_at_planting
                  leafn_xfer_p(m) = leafc_xfer_p(m) / leafcn(ivt) ! with onset
                  crop_seedc_to_leaf_p(m) = leafc_xfer_p(m)/deltim
                  crop_seedn_to_leaf_p(m) = leafn_xfer_p(m)/deltim

                  ! because leafc_xfer is set above rather than incremneted through the normal process, must also set its isotope
                  ! pools here.  use totvegc_patch as the closest analogue if nonzero, and use initial value otherwise
!                  if (use_c13) then
!                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
!                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
!                             c13_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
!                     else
!                        c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c13ratio
!                     endif
!                  endif
!                  if (use_c14) then
!                     if ( cnveg_carbonstate_inst%totvegc_patch(p) .gt. 0._r8) then
!                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * &
!                             c14_cnveg_carbonstate_inst%totvegc_patch(p) / cnveg_carbonstate_inst%totvegc_patch(p)
!                     else
!                        c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = leafc_xfer(p) * c14ratio
!                     endif
!                  endif

               else
                  gddmaturity_p(m) = 0._r8
               end if
            end if ! crop patch distinction
            ! crop phenology (gdd thresholds) controlled by gdd needed for
            ! maturity (physiological) which is based on the average gdd
            ! accumulation and hybrids in United States from April 1 - Sept 30

            ! calculate threshold from phase 1 to phase 2:
            ! threshold for attaining leaf emergence (based on fraction of
            ! gdd(i) -- climatological average)
            ! Hayhoe and Dwyer, 1990, Can. J. Soil Sci 70:493-497
            ! Carlson and Gage, 1989, Agric. For. Met., 45: 313-324
            ! J.T. Ritchie, 1991: Modeling Plant and Soil systems

            huileaf_p(m) = lfemerg(ivt) * gddmaturity_p(m) ! 3-7% in cereal

            ! calculate threshhold from phase 2 to phase 3:
            ! from leaf emergence to beginning of grain-fill period
            ! this hypothetically occurs at the end of tassling, not the beginning
            ! tassel initiation typically begins at 0.5-0.55 * gddmaturity

            ! calculate linear relationship between huigrain fraction and relative
            ! maturity rating for maize

            if (ivt == ntmp_corn .or. ivt == nirrig_tmp_corn .or. &
                ivt == ntrp_corn .or. ivt == nirrig_trp_corn .or. &
                ivt == nsugarcane .or. ivt == nirrig_sugarcane .or. &
                ivt == nmiscanthus .or. ivt == nirrig_miscanthus .or. &
                ivt == nswitchgrass .or. ivt == nirrig_switchgrass) then
               ! the following estimation of crmcorn from gddmaturity is based on a linear
               ! regression using data from Pioneer-brand corn hybrids (Kucharik, 2003,
               ! Earth Interactions 7:1-33: fig. 2)
               crmcorn = max(73._r8, min(135._r8, (gddmaturity_p(m)+ 53.683_r8)/13.882_r8))

               ! the following adjustment of grnfill based on crmcorn is based on a tuning
               ! of Agro-IBIS to give reasonable results for max LAI and the seasonal
               ! progression of LAI growth (pers. comm. C. Kucharik June 10, 2010)
               huigrain_p(m) = -0.002_r8  * (crmcorn - 73._r8) + grnfill(ivt)

               huigrain_p(m) = min(max(huigrain_p(m), grnfill(ivt)-0.1_r8), grnfill(ivt))
               huigrain_p(m) = huigrain_p(m) * gddmaturity_p(m)     ! Cabelguenne et
            else
               huigrain_p(m) = grnfill(ivt) * gddmaturity_p(m) ! al. 1999
            end if

          end if ! crop not live nor planted

         ! ----------------------------------
         ! from AgroIBIS subroutine phenocrop
         ! ----------------------------------

         ! all of the phenology changes are based on the total number of gdd needed
         ! to change to the next phase - based on fractions of the total gdd typical
         ! for  that region based on the April 1 - Sept 30 window of development

         ! crop phenology (gdd thresholds) controlled by gdd needed for
         ! maturity (physiological) which is based on the average gdd
         ! accumulation and hybrids in United States from April 1 - Sept 30

         ! Phase 1: Planting to leaf emergence (now in CNAllocation)
         ! Phase 2: Leaf emergence to beginning of grain fill (general LAI accumulation)
         ! Phase 3: Grain fill to physiological maturity and harvest (LAI decline)
         ! Harvest: if gdd past grain fill initiation exceeds limit
         ! or number of days past planting reaches a maximum, the crop has
         ! reached physiological maturity and plant is harvested;
         ! crop could be live or dead at this stage - these limits
         ! could lead to reaching physiological maturity or determining
         ! a harvest date for a crop killed by an early frost (see next comments)
         ! --- --- ---
         ! keeping comments without the code (slevis):
         ! if minimum temperature, tref_min_p <= freeze kill threshold, tkill
         ! for 3 consecutive days and lai is above a minimum,
         ! plant will be damaged/killed. This function is more for spring freeze events
         ! or for early fall freeze events

         ! spring temperate cereal is affected by this, winter cereal kill function
         ! is determined in crops.f - is a more elaborate function of
         ! cold hardening of the plant

         ! currently simulates too many grid cells killed by freezing temperatures

         ! removed on March 12 2002 - C. Kucharik
         ! until it can be a bit more refined, or used at a smaller scale.
         ! we really have no way of validating this routine
         ! too difficult to implement on 0.5 degree scale grid cells
         ! --- --- ---

          onset_flag_p(m)  = 0._r8 ! CN terminology to trigger certain
          offset_flag_p(m) = 0._r8 ! carbon and nitrogen transfers
 
          if (croplive_p(m)) then
            cphase_p(m) = 1._r8

            ! call vernalization if winter temperate cereal planted, living, and the
            ! vernalization factor is not 1;
            ! vf affects the calculation of gddtsoi & gddplant

            if (tref_min_p(m) < 1.e30_r8 .and. vf_p(m) /= 1._r8 .and. &
               (ivt == nwwheat .or. ivt == nirrig_wwheat)) then
!               call vernalization(i,m)  ! ignore winter wheat for now, add later
            end if

            ! days past planting may determine harvest

            if (jday >= idop_p(m)) then
               idpp = jday - idop_p(m)
            else
               idpp = int(dayspyr) + jday - idop_p(m)
            end if

            ! onset_counter initialized to zero when .not. croplive
            ! offset_counter relevant only at time step of harvest

            onset_counter_p(m) = onset_counter_p(m) - deltim

            ! enter phase 2 onset for one time step:
            ! transfer seed carbon to leaf emergence

            if (peaklai_p(m) >= 1) then
               gddplant_p(m) = max(gddplant_p(m),huigrain_p(m))
            endif

            if (gddtsoi_p(m) >= huileaf_p(m) .and. gddplant_p(m) < huigrain_p(m) .and. idpp < mxmat(ivt)) then
               cphase_p(m) = 2._r8
               if (abs(onset_counter_p(m)) > 1.e-6_r8) then
                  onset_flag_p(m)    = 1._r8
                  onset_counter_p(m) = deltim
                    fert_counter_p(m)  = ndays_on * 86400.
                    if (ndays_on .gt. 0) then
                       fert_p(m) = (manunitro(ivt) * 1000._r8 + fertnitro_p(m))/ fert_counter_p(m)
                    else
                       fert_p(m) = 0._r8
                    end if
               else
                  ! this ensures no re-entry to onset of phase2
                  ! b/c onset_counter(p) = onset_counter(p) - deltim
                  ! at every time step

                  onset_counter_p(m) = deltim
               end if

               ! enter harvest for one time step:
               ! - transfer live biomass to litter and to crop yield
               ! - send xsmrpool to the atmosphere
               ! if onset and harvest needed to last longer than one timestep
               ! the onset_counter would change from dt and you'd need to make
               ! changes to the offset subroutine below

            else if (gddplant_p(m) >= gddmaturity_p(m) .or. idpp >= mxmat(ivt)) then
               if (harvdate_p(m) >= NOT_Harvested) harvdate_p(m) = jday
               croplive_p(m) = .false.     ! no re-entry in greater if-block
               cphase_p(m) = 4._r8
               if (tlai_p(m) > 0._r8) then ! plant had emerged before harvest
                  offset_flag_p(m) = 1._r8
                  offset_counter_p(m) = deltim
               else                      ! plant never emerged from the ground
                  ! Revert planting transfers; this will replenish the crop seed deficit.
                  ! We subtract from any existing value in crop_seedc_to_leaf /
                  ! crop_seedn_to_leaf in the unlikely event that we enter this block of
                  ! code in the same time step where the planting transfer originally
                  ! occurred.
                  crop_seedc_to_leaf_p(m) = crop_seedc_to_leaf_p(m) - leafc_xfer_p(m)/deltim
                  crop_seedn_to_leaf_p(m) = crop_seedn_to_leaf_p(m) - leafn_xfer_p(m)/deltim
                  leafc_xfer_p(m) = 0._r8
                  leafn_xfer_p(m) = leafc_xfer_p(m) / leafcn(ivt)
!                  if(m .eq. 642821)print*,'here10,leafc_xfer',m,leafc_xfer_p(m)
!                  if (use_c13) then
!                     c13_cnveg_carbonstate_inst%leafc_xfer_patch_p(m) = 0._r8
!                  endif
!                  if (use_c14) then
!                     c14_cnveg_carbonstate_inst%leafc_xfer_patch_p(m) = 0._r8
!                  endif

               end if

               ! enter phase 3 while previous criteria fail and next is true;
               ! in terms of order, phase 3 occurs before harvest, but when
               ! harvest *can* occur, we want it to have first priority.
               ! AgroIBIS uses a complex formula for lai decline.
               ! Use CN's simple formula at least as a place holder (slevis)

            else if (gddplant_p(m) >= huigrain_p(m)) then
               cphase_p(m) = 3._r8
               bglfr_p(m) = 1._r8/(leaf_long(ivt)*dayspyr*86400.)
            end if

            ! continue fertilizer application while in phase 2;
            ! assumes that onset of phase 2 took one time step only

              if (fert_counter_p(m) <= 0._r8) then
                 fert_p(m) = 0._r8
              else ! continue same fert application every timestep
                 fert_counter_p(m) = fert_counter_p(m) - deltim
              end if

          else   ! crop not live
            ! next 2 lines conserve mass if leaf*_xfer > 0 due to interpinic.
            ! We subtract from any existing value in crop_seedc_to_leaf /
            ! crop_seedn_to_leaf in the unlikely event that we enter this block of
            ! code in the same time step where the planting transfer originally
            ! occurred.
            crop_seedc_to_leaf_p(m) = crop_seedc_to_leaf_p(m) - leafc_xfer_p(m)/deltim
            crop_seedn_to_leaf_p(m) = crop_seedn_to_leaf_p(m) - leafn_xfer_p(m)/deltim
            onset_counter_p(m) = 0._r8
            leafc_xfer_p(m) = 0._r8
            leafn_xfer_p(m) = leafc_xfer_p(m) / leafcn(ivt)
!            if (use_c13) then
!              c13_cnveg_carbonstate_inst%leafc_xfer_patch(p) = 0._r8
!            endif
!            if (use_c14) then
!               c14_cnveg_carbonstate_inst%leafc_xfer_patch(p) = 0._r8
!            endif
          end if ! croplive
         end if
      end do ! prognostic crops loop

  end subroutine CropPhenology

 !---------------------------------------------------------
  subroutine CNOnsetGrowth(i,ps,pe,deltim)
   integer, intent(in) :: i
   integer, intent(in) :: ps
   integer, intent(in) :: pe
   real(r8),intent(in) :: deltim

   ! !LOCAL VARIABLES:
   real(r8):: t1           ! temporary variable
   integer :: ivt, m

   ! only calculate these fluxes during onset period

   do m = ps, pe
      ivt = pftclass(m)
      if (onset_flag_p(m) == 1._r8) then

      ! The transfer rate is a linearly decreasing function of time,
      ! going to zero on the last timestep of the onset period

         if (onset_counter_p(m) == deltim) then
            t1 = 1.0_r8 / deltim
         else
            t1 = 2.0_r8 / (onset_counter_p(m))
         end if
         leafc_xfer_to_leafc_p(m)   = t1 * leafc_xfer_p(m)
         frootc_xfer_to_frootc_p(m) = t1 * frootc_xfer_p(m)
         leafn_xfer_to_leafn_p(m)   = t1 * leafn_xfer_p(m)
         frootn_xfer_to_frootn_p(m) = t1 * frootn_xfer_p(m)
         if (woody(ivt) == 1) then
            livestemc_xfer_to_livestemc_p(m)   = t1 * livestemc_xfer_p(m)
            deadstemc_xfer_to_deadstemc_p(m)   = t1 * deadstemc_xfer_p(m)
            livecrootc_xfer_to_livecrootc_p(m) = t1 * livecrootc_xfer_p(m)
            deadcrootc_xfer_to_deadcrootc_p(m) = t1 * deadcrootc_xfer_p(m)
            livestemn_xfer_to_livestemn_p(m)   = t1 * livestemn_xfer_p(m)
            deadstemn_xfer_to_deadstemn_p(m)   = t1 * deadstemn_xfer_p(m)
            livecrootn_xfer_to_livecrootn_p(m) = t1 * livecrootn_xfer_p(m)
            deadcrootn_xfer_to_deadcrootn_p(m) = t1 * deadcrootn_xfer_p(m)
         end if
   
      end if ! end if onset period

   ! calculate the background rate of transfer growth (used for stress
   ! deciduous algorithm). In this case, all of the mass in the transfer
   ! pools should be moved to displayed growth in each timestep.

      if (bgtr_p(m) > 0._r8) then
         leafc_xfer_to_leafc_p(m)   = leafc_xfer_p(m) / deltim
         frootc_xfer_to_frootc_p(m) = frootc_xfer_p(m) / deltim
         leafn_xfer_to_leafn_p(m)   = leafn_xfer_p(m) / deltim
         frootn_xfer_to_frootn_p(m) = frootn_xfer_p(m) / deltim
         if (woody(ivt) == 1) then
            livestemc_xfer_to_livestemc_p(m)   = livestemc_xfer_p(m) / deltim
            deadstemc_xfer_to_deadstemc_p(m)   = deadstemc_xfer_p(m) / deltim
            livecrootc_xfer_to_livecrootc_p(m) = livecrootc_xfer_p(m) / deltim
            deadcrootc_xfer_to_deadcrootc_p(m) = deadcrootc_xfer_p(m) / deltim
            livestemn_xfer_to_livestemn_p(m)   = livestemn_xfer_p(m) / deltim
            deadstemn_xfer_to_deadstemn_p(m)   = deadstemn_xfer_p(m) / deltim
            livecrootn_xfer_to_livecrootn_p(m) = livecrootn_xfer_p(m) / deltim
            deadcrootn_xfer_to_deadcrootn_p(m) = deadcrootn_xfer_p(m) / deltim
         end if
      end if ! end if bgtr
   end do

  end subroutine CNOnsetGrowth

  subroutine CNOffsetLitterfall(i,ps,pe,deltim,npcropmin)
    integer, intent(in) :: i
    integer, intent(in) :: ps
    integer, intent(in) :: pe
    real(r8),intent(in) :: deltim
    integer ,intent(in) :: npcropmin

    real(r8) :: t1           ! temporary variable
    real(r8) :: denom        ! temporary variable for divisor
    real(r8) :: ntovr_leaf
    real(r8) :: fr_leafn_to_litter ! fraction of the nitrogen turnover that goes to litter; remaining fraction is retranslocated
    integer  :: ivt, m

   do m = ps, pe
      ivt = pftclass(m)
   ! only calculate fluxes during offset period
      if (offset_flag_p(m) == 1._r8) then

         if (offset_counter_p(m) == deltim) then
            t1 = 1.0_r8 / deltim
            leafc_to_litter_p(m)  = t1 * leafc_p(m)  + cpool_to_leafc_p(m)
            frootc_to_litter_p(m) = t1 * frootc_p(m) + cpool_to_frootc_p(m)
          ! this assumes that offset_counter == dt for crops
          ! if this were ever changed, we'd need to add code to the "else"
            if (ivt >= npcropmin) then
             ! Replenish the seed deficits from grain, if there is enough
             ! available grain. (If there is not enough available grain, the seed
             ! deficits will accumulate until there is eventually enough grain to
             ! replenish them.)
               grainc_to_seed_p(m) = t1 * min(-cropseedc_deficit_p(m), grainc_p(m))
               grainn_to_seed_p(m) = t1 * min(-cropseedn_deficit_p(m), grainn_p(m))
             ! Send the remaining grain to the food product pool
               grainc_to_food_p(m) = t1 * grainc_p(m)  + cpool_to_grainc_p(m) - grainc_to_seed_p(m)
               grainn_to_food_p(m) = t1 * grainn_p(m)  + npool_to_grainn_p(m) - grainn_to_seed_p(m)

               livestemc_to_litter_p(m) = t1 * livestemc_p(m)  + cpool_to_livestemc_p(m)
            end if
         else
            t1 = deltim * 2.0_r8 / (offset_counter_p(m) * offset_counter_p(m))
            leafc_to_litter_p(m)  = prev_leafc_to_litter_p(m)  + t1*(leafc_p(m)  - prev_leafc_to_litter_p(m)*offset_counter_p(m))
            frootc_to_litter_p(m) = prev_frootc_to_litter_p(m) + t1*(frootc_p(m) - prev_frootc_to_litter_p(m)*offset_counter_p(m))

         end if

!       if ( use_fun ) then
!          if(leafc_to_litter(p)*dt.gt.leafc(p))then
!              leafc_to_litter(p) = leafc(p)/dt + cpool_to_leafc(p)
!             if (use_matrixcn) then
!                if(leafc(p) .gt. 0)then
!                   matrix_phtransfer(p,ileaf_to_iout_phc) = leafc_to_litter(p) / leafc(p)
!                else
!                   leafc_to_litter(p) = 0
!                end if
!             else
!                ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
!             end if !use_matrixcn
!          endif
!          if(frootc_to_litter(p)*dt.gt.frootc(p))then
!              frootc_to_litter(p) = frootc(p)/dt + cpool_to_frootc(p)
!             if (use_matrixcn) then
!                if(frootc(p) .gt. 0)then
!                   matrix_phtransfer(p,ifroot_to_iout_phc) = frootc_to_litter(p) / frootc(p)
!                else
!                   frootc_to_litter(p) = 0
!                end if
!             else
!                ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
!             end if !use_matrixcn
!          endif
!       end if


!       if ( use_fun ) then
!          leafc_to_litter_fun(p)      =  leafc_to_litter(p)
!          leafn_to_retransn(p)        =  paid_retransn_to_npool(p) + free_retransn_to_npool(p)
!          if (leafn(p).gt.0._r8) then
!             if (leafn(p)-leafn_to_retransn(p)*dt.gt.0._r8) then
!                 leafcn_offset(p)     =  leafc(p)/(leafn(p)-leafn_to_retransn(p)*dt)
!             else
!                 leafcn_offset(p)     =  leafc(p)/leafn(p)
!             end if
!          else
!             leafcn_offset(p)         =  leafcn(ivt(p))
!          end if
!          leafn_to_litter(p)          =  leafc_to_litter(p)/leafcn_offset(p) - leafn_to_retransn(p)
!          leafn_to_litter(p)          =  max(leafn_to_litter(p),0._r8)
!          if (use_matrixcn) then
!             if(leafn(p) .gt. 0)then
!                 matrix_nphtransfer(p,ileaf_to_iout_phn)       = (leafn_to_litter(p)) / leafn(p)
!                 matrix_nphtransfer(p,ileaf_to_iretransn_phn)  = (leafn_to_retransn(p)) / leafn(p)
!             else
!                 leafn_to_litter(p)   = 0
!                 leafn_to_retransn(p) = 0
!             end if
!          else
!             ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
!             !                                        and CNNStateUpdate1::NStateUpdate1
!          end if !use_matrixcn
!
!          denom = ( leafn_to_retransn(p) + leafn_to_litter(p) )
!          if ( denom /= 0.0_r8 ) then
!             fr_leafn_to_litter =  leafn_to_litter(p) / ( leafn_to_retransn(p) + leafn_to_litter(p) )
!          else if ( leafn_to_litter(p) == 0.0_r8 ) then
!             fr_leafn_to_litter =  0.0_r8
!          else
!             fr_leafn_to_litter =  1.0_r8
!          end if
!
!       else
!          if (CNratio_floating .eqv. .true.) then
!             fr_leafn_to_litter = 0.5_r8    ! assuming 50% of nitrogen turnover goes to litter
!          end if
          ! calculate the leaf N litterfall and retranslocation
            leafn_to_litter_p(m)   = leafc_to_litter_p(m)  / lflitcn(ivt)
            leafn_to_retransn_p(m) = (leafc_to_litter_p(m) / leafcn(ivt)) - leafn_to_litter_p(m)

!       end if

       ! calculate fine root N litterfall (no retranslocation of fine root N)
         frootn_to_litter_p(m) = frootc_to_litter_p(m) / frootcn(ivt)

!       if (CNratio_floating .eqv. .true.) then
!          if (leafc(p) == 0.0_r8) then
!             ntovr_leaf = 0.0_r8
!          else
!             ntovr_leaf = leafc_to_litter(p) * (leafn(p) / leafc(p))
!          end if

!          leafn_to_litter(p)   = fr_leafn_to_litter * ntovr_leaf
!          leafn_to_retransn(p) = ntovr_leaf - leafn_to_litter(p)
!!          if (use_matrixcn) then
!!             if(leafn(p) .gt. 0)then
!!                 matrix_nphtransfer(p,ileaf_to_iout_phn)  = (leafn_to_litter(p))/ leafn(p)
!!                 matrix_nphtransfer(p,ileaf_to_iretransn_phn)  = (leafn_to_retransn(p))/ leafn(p)
!!             else
!!                 leafn_to_litter(p)   = 0
!!                 leafn_to_retransn(p) = 0
!!             end if
!!          end if !use_matrixcn
!          if (frootc(p) == 0.0_r8) then
!              frootn_to_litter(p) = 0.0_r8
!           else
!              frootn_to_litter(p) = frootc_to_litter(p) * (frootn(p) / frootc(p))
!           end if
!!          if (use_matrixcn) then
!!             if(frootn(p) .gt. 0)then
!!                matrix_nphtransfer(p,ifroot_to_iout_phn)  = frootn_to_litter(p) / frootn(p)
!!             else
!!                frootn_to_litter(p) = 0
!!             end if
!!          else
!!             ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
!!          end if !use_matrixcn
!       end if

!       if ( use_fun ) then
!          if(frootn_to_litter(p)*dt.gt.frootn(p))then
!              frootn_to_litter(p) = frootn(p)/dt
!             if (use_matrixcn) then
!                matrix_nphtransfer(p,ifroot_to_iout_phn) = 1.0_r8 / dt
!             end if
!          endif
!       end if

         if (ivt >= npcropmin) then
          ! NOTE(slevis, 2014-12) results in -ve livestemn and -ve totpftn
          !X! livestemn_to_litter(p) = livestemc_to_litter(p) / livewdcn(ivt(p))
          ! NOTE(slevis, 2014-12) Beth Drewniak suggested this instead
            livestemn_to_litter_p(m) = livestemn_p(m) / deltim
!          if(use_matrixcn)then
!             matrix_nphtransfer(p,ilivestem_to_iout_phn) = 1.0_r8 / deltim
!          end if
         end if

       ! save the current litterfall fluxes
         prev_leafc_to_litter_p(m)  = leafc_to_litter_p(m)
         prev_frootc_to_litter_p(m) = frootc_to_litter_p(m)

      end if ! end if offset period
   end do

  end subroutine CNOffsetLitterfall

  subroutine CNBackgroundLitterfall(i,ps,pe)

    integer, intent(in) :: i
    integer, intent(in) :: ps
    integer, intent(in) :: pe

    ! !LOCAL VARIABLES:
    real(r8) :: fr_leafn_to_litter ! fraction of the nitrogen turnover that goes to litter; remaining fraction is retranslocated
    real(r8) :: ntovr_leaf  
    real(r8) :: denom       
    integer  :: ivt, m
    !-----------------------------------------------------------------------

      do m = ps , pe
         ! only calculate these fluxes if the background litterfall rate is non-zero
         ivt = pftclass(m)
         if (bglfr_p(m) > 0._r8) then
            ! units for bglfr are already 1/s
            leafc_to_litter_p(m)  = bglfr_p(m) * leafc_p(m)
            frootc_to_litter_p(m) = bglfr_p(m) * frootc_p(m)
!            if (use_matrixcn) then
!               matrix_phtransfer(p,ileaf_to_iout_phc)  = bglfr(p)
!               matrix_phtransfer(p,ifroot_to_iout_phc) = bglfr(p)
!            end if
!            if ( use_fun ) then
!               leafc_to_litter_fun(p)     = leafc_to_litter(p)
!               leafn_to_retransn(p)       = paid_retransn_to_npool(p) + free_retransn_to_npool(p)
!               if (leafn(p).gt.0._r8) then
!                  if (leafn(p)-leafn_to_retransn(p)*dt.gt.0._r8) then
!                     leafcn_offset(p)     = leafc(p)/(leafn(p)-leafn_to_retransn(p)*dt)
!                  else
!                     leafcn_offset(p)     = leafc(p)/leafn(p)
!                  end if
!               else
!                  leafcn_offset(p)        = leafcn(ivt(p))
!               end if
!               leafn_to_litter(p)         = leafc_to_litter(p)/leafcn_offset(p) - leafn_to_retransn(p)
!               leafn_to_litter(p)         = max(leafn_to_litter(p),0._r8)
!               if(use_matrixcn)then
!                  if(leafn(p) .ne. 0._r8)then
!                     matrix_nphtransfer(p,ileaf_to_iout_phn)      = (leafn_to_litter(p))/ leafn(p)
!                     matrix_nphtransfer(p,ileaf_to_iretransn_phn)  = (leafn_to_retransn(p))/ leafn(p)
!                  end if
!               else
!                  ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
!               end if !use_matrixcn

!               denom = ( leafn_to_retransn(p) + leafn_to_litter(p) )
!               if ( denom /= 0.0_r8 ) then
!                  fr_leafn_to_litter =  leafn_to_litter(p) / ( leafn_to_retransn(p) + leafn_to_litter(p) )
!               else if ( leafn_to_litter(p) == 0.0_r8 ) then
!                  fr_leafn_to_litter =  0.0_r8
!               else
!                  fr_leafn_to_litter =  1.0_r8
!               end if


!            else
!               if (CNratio_floating .eqv. .true.) then    
!                  fr_leafn_to_litter = 0.5_r8    ! assuming 50% of nitrogen turnover goes to litter
!               end if
               ! calculate the leaf N litterfall and retranslocation
               leafn_to_litter_p(m)   = leafc_to_litter_p(m)  / lflitcn(ivt)
               leafn_to_retransn_p(m) = (leafc_to_litter_p(m) / leafcn(ivt)) - leafn_to_litter_p(m)

!               if (use_matrixcn) then   
!                  if(leafn(p) .ne. 0)then
!                     matrix_nphtransfer(p,ileaf_to_iout_phn)      = (leafn_to_litter(p)) / leafn(p)
!                     matrix_nphtransfer(p,ileaf_to_iretransn_phn)  = (leafn_to_retransn(p))/ leafn(p)
!                  end if
!               else
!                  ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
!               end if !use_matrixcn
!            end if    

            ! calculate fine root N litterfall (no retranslocation of fine root N)
            frootn_to_litter_p(m) = frootc_to_litter_p(m) / frootcn(ivt)
            
!           if (CNratio_floating .eqv. .true.) then    
!              if (leafc(p) == 0.0_r8) then    
!                 ntovr_leaf = 0.0_r8    
!              else    
!                 ntovr_leaf = leafc_to_litter(p) * (leafn(p) / leafc(p))   
!              end if   
           
!              leafn_to_litter(p)   = fr_leafn_to_litter * ntovr_leaf
!              leafn_to_retransn(p) = ntovr_leaf - leafn_to_litter(p)
!               if (use_matrixcn) then   
!                 if(leafn(p) .gt. 0)then
!                    matrix_nphtransfer(p,ileaf_to_iout_phn)      = (leafn_to_litter(p))/ leafn(p)
!                    matrix_nphtransfer(p,ileaf_to_iretransn_phn)  = (leafn_to_retransn(p))/ leafn(p)
!                 else
!                    leafn_to_litter(p)   = 0
!                    leafn_to_retransn(p) = 0
!                 end if
!               else
!                  ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
!               end if !use_matrixcn
!               if (frootc(p) == 0.0_r8) then    
!                   frootn_to_litter(p) = 0.0_r8    
!               else    
!                   frootn_to_litter(p) = frootc_to_litter(p) * (frootn(p) / frootc(p))   
!               end if   
!            end if    

!            if ( use_fun ) then
!               if(frootn_to_litter(p)*dt.gt.frootn(p))then
!                    frootn_to_litter(p) = frootn(p)/dt
!               endif
!            end if

!            if (use_matrixcn) then   
!               if(frootn(p) .ne. 0)then
!                  matrix_nphtransfer(p,ifroot_to_iout_phn) = frootn_to_litter(p) / frootn(p)
!               end if
!            else
!               ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
!            end if !use_matrixcn
         end if
      end do

  end subroutine CNBackgroundLitterfall

  subroutine CNLivewoodTurnover(i,ps,pe)

    integer, intent(in) :: i
    integer, intent(in) :: ps
    integer, intent(in) :: pe

    ! !LOCAL VARIABLES:
    real(r8):: ctovr        ! temporary variable for carbon turnover
    real(r8):: ntovr        ! temporary variable for nitrogen turnover
    integer :: ivt, m
    !-----------------------------------------------------------------------
  
      do m = ps, pe
         ! only calculate these fluxes for woody types
         ivt = pftclass(m)
         if (woody(ivt) > 0._r8) then

            ! live stem to dead stem turnover

            ctovr = livestemc_p(m) * lwtop
            ntovr = ctovr / livewdcn(ivt)
            livestemc_to_deadstemc_p(m) = ctovr
            livestemn_to_deadstemn_p(m) = ctovr / deadwdcn(ivt)

!            if(use_matrixcn)then
!               matrix_phtransfer(p,ilivestem_to_ideadstem_phc)  = lwtop
!               matrix_nphtransfer(p,ilivestem_to_ideadstem_phn) = lwtop / deadwdcn(ivt(p))
!            else
!               ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
!               !                                        and CNNStateUpdate1::NStateUpdate1
!            end if
!            if (CNratio_floating .eqv. .true.) then    
!               if (livestemc(p) == 0.0_r8) then    
!                   ntovr = 0.0_r8    
!                else    
!                   ntovr = ctovr * (livestemn(p) / livestemc(p))   
!                end if   
!
!                livestemn_to_deadstemn(p) = 0.5_r8 * ntovr   ! assuming 50% goes to deadstemn 
!               if (use_matrixcn)then 
!                  if (livestemn(p) .gt. 0.0_r8) then
!                     matrix_nphtransfer(p,ilivestem_to_ideadstem_phn) = livestemn_to_deadstemn(p) / livestemn(p)
!                  else
!                     livestemn_to_deadstemn(p) = 0
!                  end if
!               else
!                  ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
!               end if
!            end if    
            
            livestemn_to_retransn_p(m)  = ntovr - livestemn_to_deadstemn_p(m)
            !matrix for livestemn_to_retransn will be added in allocation subroutine

            ! live coarse root to dead coarse root turnover

            ctovr = livecrootc_p(m) * lwtop
            ntovr = ctovr / livewdcn(ivt)
            livecrootc_to_deadcrootc_p(m) = ctovr
            livecrootn_to_deadcrootn_p(m) = ctovr / deadwdcn(ivt)
!            if(use_matrixcn)then
!               matrix_phtransfer(p,ilivecroot_to_ideadcroot_phc)  = lwtop
!               matrix_nphtransfer(p,ilivecroot_to_ideadcroot_phn) = lwtop / deadwdcn(ivt(p))
!            else
!               ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
!               !                                        and CNNStateUpdate1::NStateUpdate1
!            end if !use_matrixcn
            
!            if (CNratio_floating .eqv. .true.) then    
!              if (livecrootc(p) == 0.0_r8) then    
!                  ntovr = 0.0_r8    
!               else    
!                  ntovr = ctovr * (livecrootn(p) / livecrootc(p))   
!               end if   

!               livecrootn_to_deadcrootn(p) = 0.5_r8 * ntovr   ! assuming 50% goes to deadstemn 
!               if (use_matrixcn)then 
!                  if (livecrootn(p) .ne.0.0_r8 )then
!                     matrix_nphtransfer(p,ilivecroot_to_ideadcroot_phn) = livecrootn_to_deadcrootn(p) / livecrootn(p)
!                  end if
!               else
!                  ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
!               end if !use_matrixcn
!            end if    
            
            livecrootn_to_retransn_p(m)  = ntovr - livecrootn_to_deadcrootn_p(m)
!            if(use_fun)then
!               !TURNED OFF FLUXES TO CORRECT N ACCUMULATION ISSUE. RF. Oct 2015. 
!               livecrootn_to_retransn(p) = 0.0_r8
!               livestemn_to_retransn(p)  = 0.0_r8
!            endif
!            if(use_matrixcn)then
!               if(livecrootn(p) .gt. 0.0_r8) then
!                  matrix_nphtransfer(p,ilivecroot_to_iretransn_phn) = livecrootn_to_retransn(p) / livecrootn(p)
!               else
!                  livecrootn_to_retransn(p) = 0
!               end if
!               if(livestemn(p) .gt. 0.0_r8) then
!                  matrix_nphtransfer(p,ilivestem_to_iretransn_phn)  = livestemn_to_retransn(p) / livestemn(p)
!               else
!                  livestemn_to_retransn(p)  = 0
!               end if
!            end if !use_matrixcn

         end if
      end do ! end pft loop

  end subroutine CNLivewoodTurnover

  subroutine CNGrainToProductPools(i,ps,pe)

   integer ,intent(in) :: i
   integer ,intent(in) :: ps
   integer ,intent(in) :: pe
   
   integer m
   real(r8) wtcol

   do m = ps, pe
      wtcol = pftfrac(m)
      grainc_to_cropprodc(i) = grainc_to_cropprodc(i) + grainc_to_food_p(m) * wtcol
      grainn_to_cropprodn(i) = grainn_to_cropprodn(i) + grainn_to_food_p(m) * wtcol
   end do

  end subroutine CNGrainToProductPools

  subroutine CNLitterToColumn(i,ps,pe,nl_soil,npcropmin)

    integer ,intent(in) :: i
    integer ,intent(in) :: nl_soil
    integer ,intent(in) :: ps
    integer ,intent(in) :: pe
    integer ,intent(in) :: npcropmin

    integer j
    integer ivt,m
    real(r8):: wtcol
  
    do j = 1, nl_soil
       do m = ps,pe
          ivt   = pftclass(m)
          wtcol = pftfrac(m)
         ! leaf litter carbon fluxes
          phenology_to_met_c(j,i) = phenology_to_met_c(j,i) &
                 + leafc_to_litter_p(m) * lf_flab(ivt) * wtcol * leaf_prof_p(j,m)
          phenology_to_cel_c(j,i) = phenology_to_cel_c(j,i) &
                 + leafc_to_litter_p(m) * lf_fcel(ivt) * wtcol * leaf_prof_p(j,m)
          phenology_to_lig_c(j,i) = phenology_to_lig_c(j,i) &
                 + leafc_to_litter_p(m) * lf_flig(ivt) * wtcol * leaf_prof_p(j,m)
 !         print*,'CNLitterToProf,leaf2litter',leafc_to_litter(i),lf_flab(ivt),lf_fcel(ivt)

         ! leaf litter nitrogen fluxes
          phenology_to_met_n(j,i) = phenology_to_met_n(j,i) &
               + leafn_to_litter_p(m) * lf_flab(ivt) * wtcol * leaf_prof_p(j,m)
          phenology_to_cel_n(j,i) = phenology_to_cel_n(j,i) &
               + leafn_to_litter_p(m) * lf_fcel(ivt) * wtcol * leaf_prof_p(j,m)
          phenology_to_lig_n(j,i) = phenology_to_lig_n(j,i) &
               + leafn_to_litter_p(m) * lf_flig(ivt) * wtcol * leaf_prof_p(j,m)

         ! fine root litter carbon fluxes
          phenology_to_met_c(j,i) = phenology_to_met_c(j,i) &
               + frootc_to_litter_p(m) * fr_flab(ivt) * wtcol * froot_prof_p(j,m)
          phenology_to_cel_c(j,i) = phenology_to_cel_c(j,i) &
               + frootc_to_litter_p(m) * fr_fcel(ivt) * wtcol * froot_prof_p(j,m)
          phenology_to_lig_c(j,i) = phenology_to_lig_c(j,i) &
               + frootc_to_litter_p(m) * fr_flig(ivt) * wtcol * froot_prof_p(j,m)
!         print*,'CNLitterToProf,froot2litter',frootc_to_litter_p(m),fr_flab(ivt),fr_fcel(ivt)

         ! fine root litter nitrogen fluxes
          phenology_to_met_n(j,i) = phenology_to_met_n(j,i) &
               + frootn_to_litter_p(m) * fr_flab(ivt) * wtcol * froot_prof_p(j,m)
          phenology_to_cel_n(j,i) = phenology_to_cel_n(j,i) &
               + frootn_to_litter_p(m) * fr_fcel(ivt) * wtcol * froot_prof_p(j,m)
          phenology_to_lig_n(j,i) = phenology_to_lig_n(j,i) &
               + frootn_to_litter_p(m) * fr_flig(ivt) * wtcol * froot_prof_p(j,m)

         ! agroibis puts crop stem litter together with leaf litter
         ! so I've used the leaf lf_f* parameters instead of making
         ! new ones for now (slevis)
         ! also for simplicity I've put "food" into the litter pools

          if (ivt >= npcropmin) then ! add livestemc to litter
            ! stem litter carbon fluxes
             phenology_to_met_c(j,i) = phenology_to_met_c(j,i) &
                  + livestemc_to_litter_p(m) * lf_flab(ivt) * wtcol * leaf_prof_p(j,m)
             phenology_to_cel_c(j,i) = phenology_to_cel_c(j,i) &
                  + livestemc_to_litter_p(m) * lf_fcel(ivt) * wtcol * leaf_prof_p(j,m)
             phenology_to_lig_c(j,i) = phenology_to_lig_c(j,i) &
                  + livestemc_to_litter_p(m) * lf_flig(ivt) * wtcol * leaf_prof_p(j,m)

            ! stem litter nitrogen fluxes
             phenology_to_met_n(j,i) = phenology_to_met_n(j,i) &
                  + livestemn_to_litter_p(m) * lf_flab(ivt) * wtcol * leaf_prof_p(j,m)
             phenology_to_cel_n(j,i) = phenology_to_cel_n(j,i) &
                  + livestemn_to_litter_p(m) * lf_fcel(ivt) * wtcol * leaf_prof_p(j,m)
             phenology_to_lig_n(j,i) = phenology_to_lig_n(j,i) &
                  + livestemn_to_litter_p(m) * lf_flig(ivt) * wtcol * leaf_prof_p(j,m)

!            if (.not. use_grainproduct) then
             ! grain litter carbon fluxes
!                phenology_to_met_c(j,i) = phenology_to_met_c(j,i) &
!                   + grainc_to_food_p(m) * lf_flab(ivt) * wtcol * leaf_prof_p(j,m)
!                phenology_to_cel_c(j,i) = phenology_to_cel_c(j,i) &
!                   + grainc_to_food_p(m) * lf_fcel(ivt) * wtcol * leaf_prof_p(j,m)
!                phenology_to_lig_c(j,i) = phenology_to_lig_c(j,i) &
!                   + grainc_to_food_p(m) * lf_flig(ivt) * wtcol * leaf_prof_p(j,m)

             ! grain litter nitrogen fluxes
!                phenology_to_met_n(j,i) = phenology_to_met_n(j,i) &
!                   + grainn_to_food_p(m) * lf_flab(ivt) * wtcol * leaf_prof_p(j,m)
!                phenology_to_cel_n(j,i) = phenology_to_cel_n(j,i) &
!                   + grainn_to_food_p(m) * lf_fcel(ivt) * wtcol * leaf_prof_p(j,m)
!                phenology_to_lig_n(j,i) = phenology_to_lig_n(j,i) &
!                   + grainn_to_food_p(m) * lf_flig(ivt) * wtcol * leaf_prof_p(j,m)
!            end if
          end if
       end do  !end pft loop
    end do! end soil level loop

  end subroutine CNLitterToColumn

  subroutine vernalization(i,m)
   integer, intent(in) :: i
   integer, intent(in) :: m
  ! LOCAL VARAIBLES:
   real(r8) tcrown                     ! ?
   real(r8) vd, vd1, vd2               ! vernalization dependence
   real(r8) tkil                       ! Freeze kill threshold


       ! for all equations - temperatures must be in degrees (C)
       ! calculate temperature of crown of crop (e.g., 3 cm soil temperature)
       ! snow depth in centimeters

!       if (tref_p(m) < tfrz) then !slevis: tref inst of td=daily avg (K)
!          tcrown = 2._r8 + (tref_p(m) - tfrz) * (0.4_r8 + 0.0018_r8 * &
!               (min(snowdp(i)*100._r8, 15._r8) - 15._r8)**2)
!       else !slevis: snow_depth inst of adsnod=daily average (m)
!          tcrown = tref_p(m) - tfrz
!       end if

       ! vernalization factor calculation
       ! if vf(p) = 1.  then plant is fully vernalized - and thermal time
       ! accumulation in phase 1 will be unaffected
       ! refers to gddtsoi & gddplant, defined in the accumulation routines (slevis)
       ! reset vf, cumvd, and hdidx to 0 at planting of crop (slevis)

!       if (tref_max_p(p) > tfrz) then
!          if (tref_min_p(p) <= tfrz+15._r8) then
!             vd1      = 1.4_r8 - 0.0778_r8 * tcrown
!             vd2      = 0.5_r8 + 13.44_r8 / ((tref_max_p(p)-tref_min_p(p)+3._r8)**2) * tcrown
!             vd       = max(0._r8, min(1._r8, vd1, vd2))
!             cumvd(p) = cumvd(p) + vd
!          end if

!          if (cumvd(p) < 10._r8 .and. tref_max_p(p) > tfrz+30._r8) then
!             cumvd(p) = cumvd(p) - 0.5_r8 * (tref_max_p(p) - tfrz - 30._r8)
!          end if
!          cumvd(p) = max(0._r8, cumvd(p))       ! must be > 0

!          vf(p) = 1._r8 - p1v * (50._r8 - cumvd(p))
!          vf(p) = max(0._r8, min(vf(p), 1._r8)) ! must be between 0 - 1
!       end if

       ! calculate cold hardening of plant
       ! determines for winter cereal varieties whether the plant has completed
       ! a period of cold hardening to protect it from freezing temperatures. If
       ! not, then exposure could result in death or killing of plants.

       ! there are two distinct phases of hardening

!       if (tref_min_p(p) <= tfrz-3._r8 .or. hdidx(p) /= 0._r8) then
!          if (hdidx(p) >= hti) then   ! done with phase 1
!             hdidx(p) = hdidx(p) + 0.083_r8
!             hdidx(p) = min(hdidx(p), hti*2._r8)
!          end if

!          if (tref_max_p(p) >= tbase + tfrz + 10._r8) then
!             hdidx(p) = hdidx(p) - 0.02_r8 * (tref_max_p(p)-tbase-tfrz-10._r8)
!             if (hdidx(p) > hti) hdidx(p) = hdidx(p) - 0.02_r8 * (tref_max_p(p)-tbase-tfrz-10._r8)
!             hdidx(p) = max(0._r8, hdidx(p))
!          end if

!       else if (tcrown >= tbase-1._r8) then
!          if (tcrown <= tbase+8._r8) then
!             hdidx(p) = hdidx(p) + 0.1_r8 - (tcrown-tbase+3.5_r8)**2 / 506._r8
!             if (hdidx(p) >= hti .and. tcrown <= tbase + 0._r8) then
!                hdidx(p) = hdidx(p) + 0.083_r8
!                hdidx(p) = min(hdidx(p), hti*2._r8)
!             end if
!          end if

!          if (tref_max_p(p) >= tbase + tfrz + 10._r8) then
!             hdidx(p) = hdidx(p) - 0.02_r8 * (tref_max_p(p)-tbase-tfrz-10._r8)
!             if (hdidx(p) > hti) hdidx(p) = hdidx(p) - 0.02_r8 * (tref_max_p(p)-tbase-tfrz-10._r8)
!             hdidx(p) = max(0._r8, hdidx(p))
!          end if
!       end if

       ! calculate what the cereal killing temperature
       ! there is a linear inverse relationship between
       ! hardening of the plant and the killing temperature or
       ! threshold that the plant can withstand
       ! when plant is fully-hardened (hdidx = 2), the killing threshold is -18 C

       ! will have to develop some type of relationship that reduces LAI and
       ! biomass pools in response to cold damaged crop

!       if (tref_min_p(p) <= tfrz - 6._r8) then
!          tkil = (tbase - 6._r8) - 6._r8 * hdidx(p)
!          if (tkil >= tcrown) then
!             if ((0.95_r8 - 0.02_r8 * (tcrown - tkil)**2) >= 0.02_r8) then
!                write (*,*)  'crop damaged by cold temperatures at p,c =', p,c
!             else if (tlai(p) > 0._r8) then ! slevis: kill if past phase1
!                gddmaturity(p) = 0._r8      !         by forcing through
!                huigrain(p)    = 0._r8      !         harvest
!                write (*,*)  '95% of crop killed by cold temperatures at p,c =', p,c
!             end if
!          end if
!       end if

  end subroutine vernalization
  
end module bgc_veg_CNPhenologyMod
