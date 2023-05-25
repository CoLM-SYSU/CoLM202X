#include <define.h>
#ifdef BGC
module MOD_BGC_Veg_CNPhenology

  !--------------------
  ! !DESCRIPTION:
  ! This module holds all phenology related subroutines for natural vegetation and crop in the C and N cycle.
  ! CoLM Phenology controls the gain and loss of leaf carbon. LAI is then updated from leaf carbon.
  ! So, the seasonal variation in LAI can be simulated for different phenology types.
  !
  ! !ORIGINAL:
  ! The Community Land Model version 5.0 (CLM5)
  !
  ! !REVISION:
  ! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code sturcture.
  ! Fang Li, 2022, implemented GPAM crop model in this module.
  !
  ! !USES:
  use MOD_Const_PFT, only: &
      isevg  , issed  , isstd  , leaf_long, woody  , leafcn , frootcn, livewdcn, deadwdcn, &
      lflitcn, lf_flab, lf_fcel, lf_flig  , fr_flab, fr_fcel, fr_flig, &

! crop variables
       manunitro, lfemerg, mxmat, grnfill,  baset

  use MOD_BGC_Vars_TimeInvars, only: &
      ndays_on        , ndays_off      , fstor2tran, crit_dayl  , crit_onset_fdd, crit_onset_swi, &
      crit_offset_fdd , crit_offset_swi, soilpsi_on, soilpsi_off, lwtop, rice2pdt

  use MOD_Vars_Global, only: &
 !crop variables
      nswheat         , nirrig_swheat     , nsugarcane  , nirrig_sugarcane  , &
      nwwheat         , nirrig_wwheat     , ntmp_corn   , nirrig_tmp_corn   , &
      ntrp_corn       , nirrig_trp_corn   , nmiscanthus , nirrig_miscanthus , &
      nswitchgrass    , nirrig_switchgrass, ncotton     , nirrig_cotton     , &
      nrice           , nirrig_rice       , ntmp_soybean, nirrig_tmp_soybean, &
      ntrp_soybean    , nirrig_trp_soybean, &
      spval
  USE MOD_Const_Physical, only: tfrz

  use MOD_Vars_TimeVariables, only: &
      t_soisno, smp

  use MOD_BGC_Vars_TimeVars, only: &
      dayl, prev_dayl, prec10, prec60, prec365, prec_today, prec_daily, accumnstep

  use MOD_Vars_PFTimeVars, only: &
      tref_p       ,tlai_p

  use MOD_BGC_Vars_PFTimeVars, only: &
      tempavg_tref_p , annavg_tref_p  , gdd0_p        , gdd8_p            , &
      gdd10_p      , gdd020_p       , gdd820_p       , gdd1020_p     , nyrs_crop_active_p, &
      bglfr_p      , bgtr_p         , lgsf_p         , offset_flag_p , offset_counter_p  , &
      onset_flag_p , onset_counter_p, onset_gddflag_p, onset_gdd_p   , onset_fdd_p       , &
      onset_swi_p  , offset_fdd_p   , offset_swi_p   , dormant_flag_p, &

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

! crop variables
#ifdef CROP
      cropplant_p       , idop_p              , a5tmin_p            , a10tmin_p        , t10_p          , &
      cumvd_p           , vf_p                , cphase_p         , fert_counter_p , &
      croplive_p        , gddplant_p          , harvdate_p          , gddmaturity_p  , &
      hui_p             , peaklai_p           , &
      tref_min_p        , tref_max_p          , tref_min_inst_p  , tref_max_inst_p, &
      fertnitro_p       , plantdate_p         , &! input from files
#endif

      leaf_prof_p        , froot_prof_p        , &
      cropseedc_deficit_p, cropseedn_deficit_p


  use MOD_BGC_Vars_1DPFTFluxes, only: &
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
      cpool_to_livestemc_p           , fert_p

  use MOD_Vars_PFTimeInvars, only: pftclass, pftfrac

  use MOD_BGC_Vars_1DFluxes, only: &
      phenology_to_met_c , phenology_to_cel_c , phenology_to_lig_c, &
      phenology_to_met_n , phenology_to_cel_n , phenology_to_lig_n, &
      grainc_to_cropprodc, grainn_to_cropprodn

  use MOD_Vars_1DForcing, only: forc_prc, forc_prl

  use timemanager
  use precision
  use MOD_BGC_Daylength, only: daylength
  use spmd_task

  implicit none

  public CNPhenology

  integer, parameter :: NOT_Planted   = 999 ! If not planted   yet in year
  integer, parameter :: NOT_Harvested = 999 ! If not harvested yet in year

contains

  subroutine CNPhenology(i,ps,pe,nl_soil,idate,dz_soi,deltim,dlat,npcropmin,phase)

    ! !DESCRIPTION:
    ! The main driver of phenology model. Two phases are included:
    ! 1) phase==1: Calculates the phenology-related carbon and nitroge pool size changes,
    !              especially when specific phenology trigure is on (eg. leaf onset and offset).
    !              The pool size change rates is calculated from phase 2.
    ! 2) phase==2: Calculates phenology climatic diagnostics for onset and offset trigures
    !              Calculates the pool size change rates of all phenology processes.
    !
    ! !
    ! ORIGINAL:
    ! The Community Land Model version 5.0 (CLM5)
    !
    ! REVISION:
    ! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code sturcture.

    implicit none

    integer ,intent(in) :: i              ! patch index
    integer ,intent(in) :: ps             ! start pft index
    integer ,intent(in) :: pe             ! end pft index
    integer ,intent(in) :: nl_soil        ! number of total soil layers
    integer ,intent(in) :: idate(3)       ! current date (year, day of the year, seconds of the day)
    real(r8),intent(in) :: deltim         ! time step in seconds
    real(r8),intent(in) :: dlat           ! latitude (degrees)
    integer ,intent(in) :: npcropmin      ! first crop pft index
    real(r8),intent(in) :: dz_soi(nl_soil)! thicknesses of each soil layer
    integer ,intent(in) :: phase          ! indicator of the subroutine options (see DESCRIPTION above)

    real(r8) dayspyr
    integer  h    ! 1 for north hemsiphere; 2 for south hemisphere

    if(isleapyear(idate(1)))then
       dayspyr = 366
    else
       dayspyr = 365
    end if

    if ( phase == 1 ) then
       call CNPhenologyClimate    (i,ps,pe,idate(1:3),deltim,dayspyr,npcropmin,nl_soil,dz_soi,dlat)

       call CNEvergreenPhenology  (i,ps,pe,deltim,dayspyr)

       call CNSeasonDecidPhenology(i,ps,pe,idate(1:3),deltim,dayspyr,dlat)

       call CNStressDecidPhenology(i,ps,pe,deltim,dayspyr)

#ifdef CROP
       if(dlat >= 0)then
          h = 1
       else
          h = 2
       end if
       call CropPhenology(i,ps,pe,idate(1:3),h,deltim,dayspyr,npcropmin)
#endif
    else if ( phase == 2 ) then
    ! the same onset and offset routines are called regardless of
    ! phenology type - they depend only on onset_flag, offset_flag, bglfr, and bgtr

       call CNOnsetGrowth(i,ps,pe,deltim)

       call CNOffsetLitterfall(i,ps,pe,deltim,npcropmin)

       call CNBackgroundLitterfall(i,ps,pe)

       call CNLivewoodTurnover(i,ps,pe)

       call CNLitterToColumn(i,ps,pe,nl_soil,npcropmin)
    else
       write(*,*) 'bad phenology phase'
    end if

  end subroutine CNPhenology

  subroutine CNPhenologyClimate (i,ps,pe,idate,deltim,dayspyr,npcropmin,nl_soil,dz_soi,dlat)

    ! !DESCRIPTION:
    ! This subroutine summaries climate statistics, such as annual averaged temperature,
    ! maximum and minimum temperature, averaged precipitation over recent 10 days, 60 days and 365 days,
    ! growing degree days above 0, 8, and 10 degrees celsius. These climate statistics will be
    ! used in following phenology simulations.
    !
    ! ORIGINAL:
    ! The Community Land Model version 5.0 (CLM5)
    !
    ! REVISION:
    ! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code sturcture.

    !
    integer ,intent(in) :: i        ! patch index
    integer ,intent(in) :: ps       ! start pft index
    integer ,intent(in) :: pe       ! end pft index
    integer ,intent(in) :: idate(3) ! current date (year, days of the year, seconds of the day)
    real(r8),intent(in) :: deltim   ! time step in seconds
    real(r8),intent(in) :: dayspyr  ! days per year (days)
    integer ,intent(in) :: npcropmin! first crop pft index
    integer ,intent(in) :: nl_soil  ! number of total soil layers
    real(r8),intent(in) :: dz_soi(nl_soil) ! thicknesses of each soil layer
    real(r8),intent(in) :: dlat            ! latitude (degrees)

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

    nsteps = amin1(10._r8 * stepperday, accumnstep(i))
    prec10     (i)  = ( prec10  (i) * (nsteps - 1) + prec_today(i) ) / nsteps

    nsteps = amin1(60._r8 * stepperday, accumnstep(i))
    prec60     (i) = ( prec60 (i) * (nsteps - 1) + prec_today(i) ) / nsteps

    nsteps = amin1(365._r8 * stepperday, accumnstep(i))
    prec365    (i) = ( prec365(i) * (nsteps - 1) + prec_today(i) ) / nsteps

    call julian2monthday(idate(1),idate(2),month,mday)
    !calculate gdd0,gdd8,gdd10,gddplant for GPAM crop phenology F. Li
    do m = ps , pe
       ivt = pftclass(m)
       if(((month .ge. 4 .and. month .le. 9) .and. dlat .ge. 0) .or. &
          ((month .gt. 9 .or.  month .lt. 4) .and. dlat .lt. 0)) then
          gdd0_p (m) = gdd0_p (m) + max(0._r8, tref_p(m) - 273.15) * deltim / 86400._r8
          gdd8_p (m) = gdd8_p (m) + max(0._r8, tref_p(m) - 273.15 - 8) * deltim / 86400._r8
          gdd10_p(m) = gdd10_p(m) + max(0._r8, tref_p(m) - 273.15 - 10) * deltim / 86400._r8
       end if
#ifdef CROP
       if(croplive_p(m))then
          if((ivt == nwwheat .or. ivt == nirrig_wwheat).and.cphase_p(m) == 2._r8)then
             gddplant_p(m) = gddplant_p(m) +vf_p(m) *  max(0._r8, &
                               tref_p(m) - (273.15 + baset(ivt))) * deltim / 86400._r8
          else
             gddplant_p(m) = gddplant_p(m) + max(0._r8, &
                               tref_p(m) - (273.15 + baset(ivt))) * deltim / 86400._r8
          end if
       else
          gddplant_p(m) = 0._r8
       end if
#endif
    end do

    !calculate gdd020,gdd820,gdd1020 for gddmaturity in GPAM crop phenology F. Li
    do m = ps , pe
       ivt = pftclass(m)
       if (idate(2) == 1 .and. idate(3) ==1800)then
          if(nyrs_crop_active_p(m) == 0) then ! YR 1:
             gdd020_p(m)  = 0._r8                      ! set gdd..20 variables to 0
             gdd820_p(m)  = 0._r8                      ! and crops will not be planted
             gdd1020_p(m) = 0._r8
          else
             if (nyrs_crop_active_p(m) == 1) then                  ! <-- END of YR 1
                gdd020_p(m)  = gdd0_p(m)                           ! <-- END of YR 1
                gdd820_p(m)  = gdd8_p(m)                           ! <-- END of YR 1
                gdd1020_p(m) = gdd10_p(m)                          ! <-- END of YR 1
             else
                gdd020_p(m)  = (yravgm1* gdd020_p(m)  + gdd0_p(m))  / yravg  ! gdd..20 must be long term avgs
                gdd820_p(m)  = (yravgm1* gdd820_p(m)  + gdd8_p(m))  / yravg  ! so ignore results for yrs 1 & 2
                gdd1020_p(m) = (yravgm1* gdd1020_p(m) + gdd10_p(m)) / yravg
             end if
          end if                                                 ! <-- END of YR 1
          gdd0_p (m) = 0._r8
          gdd8_p (m) = 0._r8
          gdd10_p(m) = 0._r8
       end if
       if (isendofyear(idate,deltim)) then        ! <-- END of EVERY YR:
          nyrs_crop_active_p(m) = nyrs_crop_active_p(m) + 1
       end if
    end do

  end subroutine CNPhenologyClimate


  subroutine CNEvergreenPhenology (i,ps,pe,deltim,dayspyr)
    !
    ! !DESCRIPTION:
    ! Evergreen phenology assumes CN stock from vegetation storage pool go to transfer
    ! pool steadily with a constant rate 0.0002. All CN stock from transfer pool go to
    ! display pool immediately when it receives CN flow from storage pool. In recent version,
    ! Evergreen types only allocate NPP or N uptake to display pools. Storage and transfer
    ! pool stay 0 over whole simulation periods. Leaf litter fall simulation depends on a
    ! background turnover, which a constant parameter leaf_long was assigned (from MOD_Const_PFT.F90)
    ! to indicate the background turnover rates.
    !
    ! Allocation NPP -> DISPLAY pool -> litter
    !
    ! ORIGINAL:
    ! The Community Land Model version 5.0 (CLM5)
    !
    ! REVISION:
    ! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code sturcture.

    !
    ! !ARGUMENTS:
    integer ,intent(in) :: i        ! patch index
    integer ,intent(in) :: ps       ! start pft index
    integer ,intent(in) :: pe       ! end pft index
    real(r8),intent(in) :: deltim   ! time step in seconds
    real(r8),intent(in) :: dayspyr  ! Days per year
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

    do m = ps , pe
       ivt = pftclass(m)
       if (isevg(ivt)) then

          tranr=0.0002_r8
          ! set carbon fluxes for shifting storage pools to transfer pools
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
          leafn_storage_to_xfer_p(m)  = tranr * leafn_storage_p(m)/deltim
          frootn_storage_to_xfer_p(m) = tranr * frootn_storage_p(m)/deltim
          if (woody(ivt) == 1) then
             livestemn_storage_to_xfer_p(m)  = tranr * livestemn_storage_p(m)/deltim
             deadstemn_storage_to_xfer_p(m)  = tranr * deadstemn_storage_p(m)/deltim
             livecrootn_storage_to_xfer_p(m) = tranr * livecrootn_storage_p(m)/deltim
             deadcrootn_storage_to_xfer_p(m) = tranr * deadcrootn_storage_p(m)/deltim
          end if

          t1 = 1.0_r8 / deltim

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

  end subroutine CNEvergreenPhenology

  subroutine CNSeasonDecidPhenology(i,ps,pe,idate,deltim,dayspyr,dlat)

    !
    ! !DESCRIPTION:
    ! This routine handles the seasonal deciduous phenology code (temperate deciduous
    ! vegetation that has only one growing season per year). Seasonal deciduous phenology
    ! assumes 0 background turnover rates. NPP or N uptake only allocated to storage pool.
    ! Display pool size changes occur only in onset and offset period. All CN stock from
    ! transfer pool go to display pool immediately when it receives CN flow from storage
    ! pool.
    !
    ! Onset period:
    ! Allocation NPP -> STORAGE pool -> XFER pool -> DISPLAY pool
    ! Offset period:
    ! DISPLAY pool -> litter
    !
    ! ORIGINAL:
    ! The Community Land Model version 5.0 (CLM5)
    !
    ! REVISION:
    ! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code sturcture.

    ! !ARGUMENTS:
    integer ,intent(in) :: i       ! patch index
    integer ,intent(in) :: ps      ! start pft index
    integer ,intent(in) :: pe      ! end pft index
    integer ,intent(in) :: idate(3)! current date (year, day of the year, second of the day)
    real(r8),intent(in) :: deltim  ! time step in seconds
    real(r8),intent(in) :: dayspyr ! Days per year
    real(r8),intent(in) :: dlat    ! latitude (degree)

    !
    ! !LOCAL VARIABLES:
    real(r8):: ws_flag        !winter-summer solstice flag (0 or 1)
    real(r8):: crit_onset_gdd !critical onset growing degree-day sum
    real(r8):: soilt
    integer :: idate2_last
    integer :: ivt, m
    !-----------------------------------------------------------------------


    idate2_last = idate(2) - 1
    if(idate2_last .le. 0)idate2_last=idate2_last+365
    prev_dayl(i)=daylength(dlat,idate2_last)
    dayl(i)     =daylength(dlat,idate(2))

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

    !
    ! !DESCRIPTION:
    ! This routine handles the stress deciduous phenology code (deciduous vegetation with
    ! one or more growing season per year). NPP or N uptake only allocated to storage pool.
    !
    ! Onset period:
    ! Allocation NPP -> STORAGE pool -> XFER pool -> DISPLAY pool
    ! Offset period:
    ! DISPLAY pool -> litter
    !
    ! ORIGINAL:
    ! The Community Land Model version 5.0 (CLM5)
    !
    ! REVISION:
    ! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code sturcture.

    integer, intent(in) :: i       ! patch index
    integer ,intent(in) :: ps      ! start pft index
    integer ,intent(in) :: pe      ! end pft index
    real(r8),intent(in) :: deltim  ! time step in seconds
    real(r8),intent(in) :: dayspyr ! days per year

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
       if (isstd(ivt)) then
          soilt = t_soisno(3,i)
          psi = smp(3,i) * 1.e-5 ! mmH2O -> MPa

      ! onset gdd sum from Biome-BGC, v4.1.2
          crit_onset_gdd = exp(4.8_r8 + 0.13_r8*(annavg_tref_p(m) - 273.15_r8))

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
            ! if additional constraint condition not met,  set to false
             if ((prec10(i) * (3600.0_r8*10.0_r8*24.0_r8)) < rain_threshold) then
                additional_onset_condition = .false.
             endif

             if (psi >= soilpsi_on) then
                onset_swi_p(m) = onset_swi_p(m) + deltim/86400._r8
             endif

         ! if critical soil water index is exceeded, set onset_flag, and
         ! then test for soil temperature criteria

         ! Adding in Kyla's rainfall trigger when fun on. RF. prec10 (mm/s) needs to be higher than 8mm over 10 days.

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

#ifdef CROP
  subroutine CropPhenology(i,ps,pe,idate,h,deltim,dayspyr,npcropmin)

    ! !DESCRIPTION:
    ! GPAM crop phenology and code from CN to
    ! handle CN fluxes during the phenological onset & offset periods.
    !
    ! ORIGINAL: The Community Land Model version 5.0 (CLM5.0)
    !
    ! REVISION: F. Li, 2022, implemented GPAM in CoLM.

    integer, intent(in) :: i         ! patch index
    integer ,intent(in) :: ps        ! start pft index
    integer ,intent(in) :: pe        ! end pft index
    integer ,intent(in) :: idate(1:3)! current date (year, day of the year, second of the day)
    integer ,intent(in) :: h         ! hemisphere indicator: 1 for north hemisphere; 2 for south hemisphere
    real(r8),intent(in) :: deltim    ! timestep in seconds
    real(r8),intent(in) :: dayspyr   ! days per year
    integer ,intent(in) :: npcropmin ! first crop pft index

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

    ndays_on = 20._r8 ! number of days to fertilize

         ! background litterfall and transfer rates; long growing season factor
    do m = ps, pe
       ivt = pftclass(m)
       if(ivt >= npcropmin)then
          bglfr_p(m) = 0._r8 ! this value changes later in a crop's life cycle
          bgtr_p(m)  = 0._r8
          lgsf_p(m)  = 0._r8


         !  plantdate is read in
         !  determine if the cft is planted in this time step
          if ( (.not. croplive_p(m)) .and. (.not. cropplant_p(m)) ) then
             if (jday == int(plantdate_p(m))) then
                cumvd_p(m)       = 0._r8
                vf_p(m)          = 0._r8
                croplive_p(m)    = .true.
                cropplant_p(m)   = .true.
                idop_p(m)        = jday
                harvdate_p(m)    = NOT_Harvested
                leafc_xfer_p(m)  = initial_seed_at_planting
                leafn_xfer_p(m)  = leafc_xfer_p(m) / leafcn(ivt) ! with onset
                crop_seedc_to_leaf_p(m) = leafc_xfer_p(m)/deltim
                crop_seedn_to_leaf_p(m) = leafn_xfer_p(m)/deltim
             end if
          end if
         ! calculate gddmaturity
          if(croplive_p(m))then
             if (ivt == nwwheat .or. ivt == nirrig_wwheat)then
                gddmaturity_p(m) = 0.42_r8 * gdd1020_p(m) + 440._r8
             end if
             if ( ivt == ntmp_soybean .or. ivt == nirrig_tmp_soybean .or. &
                  ivt == ntrp_soybean .or. ivt == nirrig_trp_soybean) then
                gddmaturity_p(m) = 0.30_r8 * gdd1020_p(m) + 710._r8
             end if
             if (ivt == ntmp_corn    .or. ivt == nirrig_tmp_corn    .or. &
                 ivt == ntrp_corn    .or. ivt == nirrig_trp_corn    .or. &
                 ivt == nsugarcane   .or. ivt == nirrig_sugarcane   .or. &
                 ivt == nmiscanthus  .or. ivt == nirrig_miscanthus  .or. &
                 ivt == nswitchgrass .or. ivt == nirrig_switchgrass) then
                gddmaturity_p(m) = 0.30_r8 * gdd820_p(m) + 816._r8
             end if
             if (ivt == nswheat .or. ivt == nirrig_swheat .or. &
                 ivt == ncotton .or. ivt == nirrig_cotton)then
                gddmaturity_p(m) = 0.24_r8 * gdd020_p(m) + 1349._r8
             end if
             if (ivt == nrice   .or. ivt == nirrig_rice) then
                gddmaturity_p(m) = 0.35_r8 * gdd020_p(m) + 587._r8
             end if
             hui_p(m)=gddplant_p(m)/gddmaturity_p(m)
          end if

         ! all of the phenology changes are based on hui

         ! Phase 1: Planting to leaf emergence
         ! Phase 2: Leaf emergence to beginning of grain fill (LAI increase)
         ! Phase 3: Grain fill to physiological maturity and harvest (LAI decline)
         ! Harvest: if gdd past grain fill initiation exceeds limit
         ! or number of days past planting reaches a maximum, the crop has
         ! reached physiological maturity and plant is harvested;
         ! --- --- ---

          onset_flag_p(m)  = 0._r8 ! CN terminology to trigger certain
          offset_flag_p(m) = 0._r8 ! carbon and nitrogen transfers

          if (croplive_p(m)) then
             cphase_p(m) = 1._r8
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
                hui_p(m) = max(hui_p(m),grnfill(ivt))
             endif

             if (hui_p(m) >= lfemerg(ivt) .and. hui_p(m) < grnfill(ivt) .and. idpp < mxmat(ivt)) then
                cphase_p(m) = 2._r8
            ! call vernalization if winter temperate cereal planted, living, and the
            ! vernalization factor is not 1;
            ! vf affects the calculation of gddplant
                if ( vf_p(m) /= 1._r8 .and. (ivt == nwwheat .or. ivt == nirrig_wwheat) .and. hui_p(m) < 0.8_r8 * grnfill(ivt)) then
                   call vernalization(i,m,deltim)
                end if

           !fertilization

                if (abs(onset_counter_p(m)) > 1.e-6_r8) then
                   onset_flag_p(m)    = 1._r8
                   onset_counter_p(m) = deltim
                   fert_counter_p(m)  = ndays_on * 86400.
                   if (ndays_on .gt. 0) then
#ifdef FERT
                      fert_p(m) = (manunitro(ivt) * 1000._r8 + fertnitro_p(m))/ fert_counter_p(m)
#else
                      fert_p(m) = 0._r8
#endif
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

             else if (hui_p(m) >= 1._r8 .or. idpp >= mxmat(ivt)) then
                if (harvdate_p(m) >= NOT_Harvested) harvdate_p(m) = jday
                croplive_p(m) = .false.     ! no re-entry in greater if-block
                cropplant_p(m)=.false.
                cphase_p(m) = 4._r8
                hui_p(m)=0._r8
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
                end if

               ! enter phase 3 while previous criteria fail and next is true;
               ! in terms of order, phase 3 occurs before harvest, but when
               ! harvest *can* occur, we want it to have first priority.
               ! AgroIBIS uses a complex formula for lai decline.
               ! Use CN's simple formula at least as a place holder (slevis)

             else if (hui_p(m) >= grnfill(ivt)) then
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
          end if ! croplive
       end if
    end do ! prognostic crops loop

  end subroutine CropPhenology
#endif

 !---------------------------------------------------------
  subroutine CNOnsetGrowth(i,ps,pe,deltim)

    ! !DESCRIPTION:
    ! Calculates flux from transfer CN to display CN during onset period.
    ! Transfer CN -> DISPLAY CN

    ! ORIGINAL:
    ! The Community Land Model version 5.0 (CLM5)
    !
    ! REVISION:
    ! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code sturcture.

    integer, intent(in) :: i      ! patch index
    integer, intent(in) :: ps     ! start pft index
    integer, intent(in) :: pe     ! end pft index
    real(r8),intent(in) :: deltim ! time step in seconds

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
    ! !DESCRIPTION:
    ! Calculates flux from display CN to litter CN during offset period.
    ! DISPLAY CN -> litter CN

    ! ORIGINAL:
    ! The Community Land Model version 5.0 (CLM5)
    !
    ! REVISION:
    ! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code sturcture.

    integer, intent(in) :: i          ! patch index
    integer, intent(in) :: ps         ! start pft index
    integer, intent(in) :: pe         ! end pft index
    real(r8),intent(in) :: deltim     ! time step in seconds
    integer ,intent(in) :: npcropmin  ! first crop pft index

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

          leafn_to_litter_p(m)   = leafc_to_litter_p(m)  / lflitcn(ivt)
          leafn_to_retransn_p(m) = (leafc_to_litter_p(m) / leafcn(ivt)) - leafn_to_litter_p(m)


       ! calculate fine root N litterfall (no retranslocation of fine root N)
          frootn_to_litter_p(m) = frootc_to_litter_p(m) / frootcn(ivt)

          if (ivt >= npcropmin) then
          ! NOTE(slevis, 2014-12) results in -ve livestemn and -ve totpftn
          !X! livestemn_to_litter(p) = livestemc_to_litter(p) / livewdcn(ivt(p))
          ! NOTE(slevis, 2014-12) Beth Drewniak suggested this instead
             livestemn_to_litter_p(m) = livestemn_p(m) / deltim
          end if

       ! save the current litterfall fluxes
          prev_leafc_to_litter_p(m)  = leafc_to_litter_p(m)
          prev_frootc_to_litter_p(m) = frootc_to_litter_p(m)

       end if ! end if offset period
    end do

  end subroutine CNOffsetLitterfall

  subroutine CNBackgroundLitterfall(i,ps,pe)

    ! !DESCRIPTION:
    ! Calculate leaf and fine root background turnover.
    ! DISPLAY -> litter
    !
    ! ORIGINAL:
    ! The Community Land Model version 5.0 (CLM5)
    !
    ! REVISION:
    ! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code sturcture.

    integer, intent(in) :: i  ! patch index
    integer, intent(in) :: ps ! start pft index
    integer, intent(in) :: pe ! end pft index

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
               ! calculate the leaf N litterfall and retranslocation
          leafn_to_litter_p(m)   = leafc_to_litter_p(m)  / lflitcn(ivt)
          leafn_to_retransn_p(m) = (leafc_to_litter_p(m) / leafcn(ivt)) - leafn_to_litter_p(m)

          frootn_to_litter_p(m) = frootc_to_litter_p(m) / frootcn(ivt)
       end if
    end do

  end subroutine CNBackgroundLitterfall

  subroutine CNLivewoodTurnover(i,ps,pe)

    ! !DESCRIPTION:
    ! Livewood transfer to deadwood each year
    !
    ! !ORIGINAL:
    ! Community Land Model Version 5.0 (CLM5)
    !
    ! !REVISION:
    ! Xingjie Lu, 2022, Revise the code to be compatible with CoLM code structure.

    integer, intent(in) :: i  ! patch index
    integer, intent(in) :: ps ! start pft index
    integer, intent(in) :: pe ! end pft index

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

          livestemn_to_retransn_p(m)  = ntovr - livestemn_to_deadstemn_p(m)
            !matrix for livestemn_to_retransn will be added in allocation subroutine

            ! live coarse root to dead coarse root turnover

          ctovr = livecrootc_p(m) * lwtop
          ntovr = ctovr / livewdcn(ivt)
          livecrootc_to_deadcrootc_p(m) = ctovr
          livecrootn_to_deadcrootn_p(m) = ctovr / deadwdcn(ivt)

          livecrootn_to_retransn_p(m)  = ntovr - livecrootn_to_deadcrootn_p(m)
       end if
    end do ! end pft loop

  end subroutine CNLivewoodTurnover

  subroutine CNGrainToProductPools(i,ps,pe)
    ! !DESCRIPTION:
    ! Summary crop production C & N from pft to patch
    !
    ! !ORIGINAL:
    ! Community Land Model Version 5.0 (CLM5)
    !
    ! !REVISION:
    ! Xingjie Lu, 2022, revised the code to be compatible with CoLM code structure
    integer ,intent(in) :: i  ! patch index
    integer ,intent(in) :: ps ! start pft index
    integer ,intent(in) :: pe ! end pft index

    integer m
    real(r8) wtcol

    do m = ps, pe
       wtcol = pftfrac(m)
       grainc_to_cropprodc(i) = grainc_to_cropprodc(i) + grainc_to_food_p(m) * wtcol
       grainn_to_cropprodn(i) = grainn_to_cropprodn(i) + grainn_to_food_p(m) * wtcol
    end do

  end subroutine CNGrainToProductPools

  subroutine CNLitterToColumn(i,ps,pe,nl_soil,npcropmin)

    ! !DESCRIPTION:
    ! Calculate column level litterfall flux from pft level litterfall.

    ! ORIGINAL:
    ! The Community Land Model version 5.0 (CLM5)
    !
    ! REVISION:
    ! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code sturcture.

    integer ,intent(in) :: i           ! patch index
    integer ,intent(in) :: nl_soil     ! number of total soil layers
    integer ,intent(in) :: ps          ! start pft index
    integer ,intent(in) :: pe          ! end pft index
    integer ,intent(in) :: npcropmin   ! first crop pft index

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

          end if
       end do  !end pft loop
    end do! end soil level loop

  end subroutine CNLitterToColumn

#ifdef CROP
  subroutine vernalization(i,m,deltim)
    ! !DESCRIPTION
    ! vernalizatoin for winter wheat.

    ! ORITINAL: F. Li's GPAM model
   integer, intent(in) :: i      ! patch index
   integer, intent(in) :: m      ! pft index
   real(r8),intent(in) :: deltim ! time step in seconds
  ! LOCAL VARAIBLES:
    real(r8) vtmin,vtopt,vtmax          ! vernalization minimum, optimum, maximum temperature
    real(r8) alpha                      ! parameter in calculating vernalization rate
    real(r8) tc                         ! t_ref2m in degree C
    real(r8) dt          ! convert dtime from sec to hour


       ! for all equations - temperatures must be in degrees (C)
       ! calculate temperature of crown of crop (e.g., 3 cm soil temperature)
       ! snow depth in centimeters

        vtmin=-1.3_r8
        vtopt=4.9_r8
        vtmax=15.7_r8
        dt=deltim/3600.0_r8  !dt is the time step in hour
        alpha=log(2._r8)/log((vtmax-vtmin)/(vtopt-vtmin))

        tc = tref_p(m)-tfrz
        if(tc >=vtmin .and. tc <= vtmax) then
         cumvd_p(m)=cumvd_p(m) + (2._r8*((tc-vtmin)**alpha)*(vtopt-vtmin)**alpha &
                 - (tc-vtmin)**(2._r8*alpha))/(vtopt-vtmin)**(2._r8*alpha)*(dt/24._r8)
        end if

        vf_p(m)=(cumvd_p(m)**5._r8)/(22.5_r8**5._r8+cumvd_p(m)**5._r8)

  end subroutine vernalization
#endif

end module MOD_BGC_Veg_CNPhenology
#endif
