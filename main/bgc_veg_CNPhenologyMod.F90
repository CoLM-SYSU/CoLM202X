module bgc_veg_CNPhenologyMod

  use PFT_Const, only: &
      isevg  , issed  , isstd  , leaf_long, woody  , leafcn , frootcn, livewdcn, deadwdcn, &
      lflitcn, lf_flab, lf_fcel, lf_flig  , fr_flab, fr_fcel, fr_flig

  use MOD_TimeInvariants, only: &
      ndays_on        , ndays_off      , fstor2tran, crit_dayl  , crit_onset_fdd, crit_onset_swi, &
      crit_offset_fdd , crit_offset_swi, soilpsi_on, soilpsi_off, lwtop

  use MOD_TimeVariables, only: &
      t_soisno, smp, dayl, prev_dayl, prec10, prec60, prec365, prec_today, prec_daily

  use MOD_PFTimeVars, only: &
      tref_p       , tempavg_tref_p , annavg_tref_p  , gdd0_p        , gdd8_p            , &
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

      leaf_prof_p        , froot_prof_p        , &
      cropseedc_deficit_p, cropseedn_deficit_p
      
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

      grainc_to_seed_p               , grainn_to_seed_p               , &
      grainc_to_food_p               , grainn_to_food_p               , &
      cpool_to_grainc_p              , npool_to_grainn_p              , &
      livestemc_to_litter_p          , livestemn_to_litter_p          , &
      cpool_to_livestemc_p                

  use MOD_PFTimeInvars, only: pftclass, pftfrac

  use MOD_1D_Fluxes, only: &
      phenology_to_met_c , phenology_to_cel_c , phenology_to_lig_c, &
      phenology_to_met_n , phenology_to_cel_n , phenology_to_lig_n

  use MOD_1D_forcing, only: forc_prc, forc_prl

  use timemanager
  use precision
  use bgc_DaylengthMod, only: daylength

implicit none

public CNPhenology

contains

subroutine CNPhenology(i,ps,pe,nl_soil,idate,deltim,dlat,npcropmin,phase)

implicit none

integer ,intent(in) :: i
integer ,intent(in) :: ps
integer ,intent(in) :: pe
integer ,intent(in) :: nl_soil
integer ,intent(in) :: idate(3)
real(r8),intent(in) :: deltim
real(r8),intent(in) :: dlat
integer ,intent(in) :: npcropmin
integer ,intent(in) :: phase

real(r8) dayspyr

   if(isleapyear(idate(1)))then
      dayspyr = 366
   else
      dayspyr = 365
   end if

 if ( phase == 1 ) then
    call CNPhenologyClimate    (i,ps,pe,idate(1:3),deltim,dayspyr)

    call CNEvergreenPhenology  (i,ps,pe,deltim,dayspyr)

    call CNSeasonDecidPhenology(i,ps,pe,idate(1:3),deltim,dayspyr,dlat)

    call CNStressDecidPhenology(i,ps,pe,deltim,dayspyr)

!   if (doalb .and. num_pcropp > 0 ) then
!      call CropPhenology(num_pcropp, filter_pcropp, &
!           waterdiagnosticbulk_inst, temperature_inst, crop_inst, canopystate_inst, cnveg_state_inst, &
!           cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
!           c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst)
!   end if
 else if ( phase == 2 ) then
    ! the same onset and offset routines are called regardless of
    ! phenology type - they depend only on onset_flag, offset_flag, bglfr, and bgtr

    call CNOnsetGrowth(i,ps,pe,deltim)

    call CNOffsetLitterfall(i,ps,pe,deltim,npcropmin)

    call CNBackgroundLitterfall(i,ps,pe)

    call CNLivewoodTurnover(i,ps,pe)

!    call CNGrainToProductPools()

    call CNLitterToColumn(i,ps,pe,nl_soil,npcropmin)
 else
    write(*,*) 'bad phenology phase'
 end if

end subroutine CNPhenology

subroutine CNPhenologyClimate (i,ps,pe,idate,deltim,dayspyr)
   !
   ! !DESCRIPTION:
   integer ,intent(in) :: i
   integer ,intent(in) :: ps
   integer ,intent(in) :: pe
   integer ,intent(in) :: idate(3)
   real(r8),intent(in) :: deltim
   real(r8),intent(in) :: dayspyr ! days per year (days)
   ! For coupled carbon-nitrogen code (CN).
   !
   ! !LOCAL VARIABLES:
   real(r8), parameter :: yravg   = 20.0_r8      ! length of years to average for gdd
   real(r8), parameter :: yravgm1 = yravg-1.0_r8 ! minus 1 of above
   integer :: m
   !-----------------------------------------------------------------------

     ! set time steps

   prec_today(i) = prec_today(i) + (forc_prc(i) + forc_prl(i))*deltim
   if(isendofday(idate,deltim))then
      prec_daily (1:364,i) = prec_daily(2:365,i)
      prec_daily (365  ,i) = prec_today(i)
      prec_today (i) = 0._r8
      prec10     (i) = sum(prec_daily(356:365,i))/( 10._r8*86400._r8)
      prec60     (i) = sum(prec_daily(306:365,i))/( 60._r8*86400._r8)
      prec365    (i) = sum(prec_daily(  1:365,i))/(365._r8*86400._r8)
   end if

   do m = ps , pe
      tempavg_tref_p(m) = tempavg_tref_p(m) + tref_p(m) * (deltim/86400._r8/dayspyr)
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
      if (idate(2) == 1 .and. nyrs_crop_active_p(i) == 0) then ! YR 1:
         gdd020_p(m)  = 0._r8                             ! set gdd..20 variables to 0
         gdd820_p(m)  = 0._r8                             ! and crops will not be planted
         gdd1020_p(m) = 0._r8
      end if
      if (idate(2) == 1 .and. idate(3) == 0) then        ! <-- END of EVERY YR:
         if (nyrs_crop_active_p(m) == 1) then                     ! <-- END of YR 1
            gdd020_p(m)  = gdd0_p(m)                                ! <-- END of YR 1
            gdd820_p(m)  = gdd8_p(m)                                ! <-- END of YR 1
            gdd1020_p(m) = gdd10_p(m)                               ! <-- END of YR 1
         end if                                                 ! <-- END of YR 1
         gdd020_p(m)  = (yravgm1* gdd020_p(m)  + gdd0_p(m))  / yravg  ! gdd..20 must be long term avgs
         gdd820_p(m)  = (yravgm1* gdd820_p(m)  + gdd8_p(m))  / yravg  ! so ignore results for yrs 1 & 2
         gdd1020_p(m) = (yravgm1* gdd1020_p(m) + gdd10_p(m)) / yravg
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
            bglfr_p(i) = 1._r8/(leaf_long(ivt) * dayspyr * 86400._r8)
            bgtr_p(i)  = 0._r8
            lgsf_p(i)  = 0._r8
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
!      if (use_matrixcn)then
!         matrix_phtransfer(p,ileafxf_to_ileaf_phc)   = matrix_phtransfer(p,ileafxf_to_ileaf_phc) + t1
!         matrix_phtransfer(p,ifrootxf_to_ifroot_phc) = matrix_phtransfer(p,ifrootxf_to_ifroot_phc) + t1
!         matrix_nphtransfer(p,ileafxf_to_ileaf_phn)   = matrix_nphtransfer(p,ileafxf_to_ileaf_phn) + t1
!         matrix_nphtransfer(p,ifrootxf_to_ifroot_phn) = matrix_nphtransfer(p,ifrootxf_to_ifroot_phn) + t1
!         if (woody(ivt(p)) == 1.0_r8) then
!
!            matrix_phtransfer(p,ilivestemxf_to_ilivestem_phc)   = matrix_phtransfer(p,ilivestemxf_to_ilivestem_phc) + t1
!            matrix_phtransfer(p,ideadstemxf_to_ideadstem_phc)   = matrix_phtransfer(p,ideadstemxf_to_ideadstem_phc) + t1
!            matrix_phtransfer(p,ilivecrootxf_to_ilivecroot_phc) = matrix_phtransfer(p,ilivecrootxf_to_ilivecroot_phc) + t1
!            matrix_phtransfer(p,ideadcrootxf_to_ideadcroot_phc) = matrix_phtransfer(p,ideadcrootxf_to_ideadcroot_phc) + t1
!
!            matrix_nphtransfer(p,ilivestemxf_to_ilivestem_phn)   = matrix_nphtransfer(p,ilivestemxf_to_ilivestem_phn) + t1
!            matrix_nphtransfer(p,ideadstemxf_to_ideadstem_phn)   = matrix_nphtransfer(p,ideadstemxf_to_ideadstem_phn) + t1
!            matrix_nphtransfer(p,ilivecrootxf_to_ilivecroot_phn) = matrix_nphtransfer(p,ilivecrootxf_to_ilivecroot_phn) + t1
!            matrix_nphtransfer(p,ideadcrootxf_to_ideadcroot_phn) = matrix_nphtransfer(p,ideadcrootxf_to_ideadcroot_phn) + t1
!         end if
!      else
!         ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
!         !                                        and CNNStateUpdate1::NStateUpdate1
!      end if !use_matrixcn
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
!      if(use_matrixcn)then
!         matrix_phtransfer(p,ileafxf_to_ileaf_phc) = matrix_phtransfer(p,ileafxf_to_ileaf_phc) + 1.0_r8 / dt
!         matrix_phtransfer(p,ifrootxf_to_ifroot_phc) = matrix_phtransfer(p,ifrootxf_to_ifroot_phc) + 1.0_r8 / dt
!         matrix_nphtransfer(p,ileafxf_to_ileaf_phn) = matrix_nphtransfer(p,ileafxf_to_ileaf_phn) + 1.0_r8 / dt
!         matrix_nphtransfer(p,ifrootxf_to_ifroot_phn) = matrix_nphtransfer(p,ifrootxf_to_ifroot_phn) + 1.0_r8 / dt
!         if (woody(ivt(p)) == 1.0_r8) then
!
!            matrix_phtransfer(p,ilivestemxf_to_ilivestem_phc)   = matrix_phtransfer(p,ilivestemxf_to_ilivestem_phc) + 1.0_r8 / dt
!            matrix_phtransfer(p,ideadstemxf_to_ideadstem_phc)   = matrix_phtransfer(p,ideadstemxf_to_ideadstem_phc) + 1.0_r8 / dt
!            matrix_phtransfer(p,ilivecrootxf_to_ilivecroot_phc) = matrix_phtransfer(p,ilivecrootxf_to_ilivecroot_phc) + 1.0_r8 / dt
!            matrix_phtransfer(p,ideadcrootxf_to_ideadcroot_phc) = matrix_phtransfer(p,ideadcrootxf_to_ideadcroot_phc) + 1.0_r8 / dt
!
!            matrix_nphtransfer(p,ilivestemxf_to_ilivestem_phn)   = matrix_nphtransfer(p,ilivestemxf_to_ilivestem_phn) + 1.0_r8 / dt
!            matrix_nphtransfer(p,ideadstemxf_to_ideadstem_phn)   = matrix_nphtransfer(p,ideadstemxf_to_ideadstem_phn) + 1.0_r8 / dt
!            matrix_nphtransfer(p,ilivecrootxf_to_ilivecroot_phn) = matrix_nphtransfer(p,ilivecrootxf_to_ilivecroot_phn) + 1.0_r8 / dt
!            matrix_nphtransfer(p,ideadcrootxf_to_ideadcroot_phn) = matrix_nphtransfer(p,ideadcrootxf_to_ideadcroot_phn) + 1.0_r8 / dt
!         end if
!      else
!         ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
!         !                                        and CNNStateUpdate1::NStateUpdate1
!      end if !use_matrixcn
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
!          if (use_matrixcn) then
!             if(leafc(p) .gt. 0)then
!                matrix_phtransfer(p,ileaf_to_iout_phc)  = leafc_to_litter(p) / leafc(p)
!             else
!                leafc_to_litter(p) = 0
!             end if
!             if(frootc(p) .gt. 0)then
!                matrix_phtransfer(p,ifroot_to_iout_phc) = frootc_to_litter(p) / frootc(p)
!             else
!                frootc_to_litter(p) = 0
!             end if
!          else
!             ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
!             !                                        and CNNStateUpdate1::NStateUpdate1
!          end if ! use_matrixcn
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
!             if(use_matrixcn)then
!                if(grainc(p) .gt. 0)then
!                   matrix_phtransfer(p,igrain_to_iout_phc)  = (grainc_to_seed(p) + grainc_to_food(p)) / grainc(p)
!                else
!                   grainc_to_seed(p) = 0
!                   grainc_to_food(p) = 0
!                end if
!                if(grainn(p) .gt. 0)then
!                   matrix_nphtransfer(p,igrain_to_iout_phn) = (grainn_to_seed(p) + grainn_to_food(p)) / grainn(p)
!                else
!                   grainn_to_seed(p) = 0
!                   grainn_to_food(p) = 0
!                end if
!                if(livestemc(p) .gt. 0)then
!                   matrix_phtransfer(p,ilivestem_to_iout_phc) = livestemc_to_litter(p) / livestemc(p)
!                else
!                   livestemc_to_litter(p) = 0
!                end if
!             else
!                ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
!                !                                        and CNNStateUpdate1::NStateUpdate1
!             end if ! use_matrixcn
            end if
         else
            t1 = deltim * 2.0_r8 / (offset_counter_p(m) * offset_counter_p(m))
            leafc_to_litter_p(m)  = prev_leafc_to_litter_p(m)  + t1*(leafc_p(m)  - prev_leafc_to_litter_p(m)*offset_counter_p(m))
            frootc_to_litter_p(m) = prev_frootc_to_litter_p(m) + t1*(frootc_p(m) - prev_frootc_to_litter_p(m)*offset_counter_p(m))

!          if (use_matrixcn) then
!             if(leafc(p) .gt. 0)then
!                matrix_phtransfer(p,ileaf_to_iout_phc)  = leafc_to_litter(p) / leafc(p)
!             else
!                leafc_to_litter(p) = 0
!             end if
!             if(frootc(p) .gt. 0)then
!                matrix_phtransfer(p,ifroot_to_iout_phc) = frootc_to_litter(p) / frootc(p)
!             else
!                frootc_to_litter(p) = 0
!             end if
!          else
!             ! NOTE: The non matrix version of this is in CNCStateUpdate1::CStateUpdate1 EBK (11/26/2019)
!             !                                        and CNNStateUpdate1::NStateUpdate1
!          end if !use_matrixcn
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

!          if (use_matrixcn) then
!             if(leafn(p) .gt. 0)then
!                 matrix_nphtransfer(p,ileaf_to_iout_phn)       = (leafn_to_litter(p))/ leafn(p)
!                 matrix_nphtransfer(p,ileaf_to_iretransn_phn)   = (leafn_to_retransn(p))/ leafn(p)
!             else
!                 leafn_to_litter(p)   = 0
!                 leafn_to_retransn(p) = 0
!             end if
!          else
!             ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
!          end if !use_matrixcn
!       end if

       ! calculate fine root N litterfall (no retranslocation of fine root N)
         frootn_to_litter_p(m) = frootc_to_litter_p(m) / frootcn(ivt)
!       if (use_matrixcn) then
!          if(frootn(p) .gt. 0)then
!             matrix_nphtransfer(p,ifroot_to_iout_phn)  = frootn_to_litter(p) / frootn(p)
!          else
!             frootn_to_litter(p) = 0
!          end if
!       else
!          ! NOTE: The non matrix version of this is in CNNStateUpdate1::NStateUpdate1 EBK (11/26/2019)
!       end if !use_matrixcn

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

!  subroutine CNGrainToProductPools(i)
!
!   integer ,intent(in) :: i
!
!    grainc_to_cropprodc(i) = grainc_to_food(i)
!    grainn_to_cropprodn(i) = grainn_to_food(i)
!
!  end subroutine CNGrainToProductPools

  subroutine CNLitterToColumn(i,ps,pe,nl_soil,npcropmin)

    integer ,intent(in) :: i
    integer ,intent(in) :: nl_soil
    integer ,intent(in) :: ps
    integer ,intent(in) :: pe
    integer ,intent(in) :: npcropmin

    integer j
    integer ivt,m, wtcol
  
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
                phenology_to_met_c(j,i) = phenology_to_met_c(j,i) &
                   + grainc_to_food_p(m) * lf_flab(ivt) * wtcol * leaf_prof_p(j,m)
                phenology_to_cel_c(j,i) = phenology_to_cel_c(j,i) &
                   + grainc_to_food_p(m) * lf_fcel(ivt) * wtcol * leaf_prof_p(j,m)
                phenology_to_lig_c(j,i) = phenology_to_lig_c(j,i) &
                   + grainc_to_food_p(m) * lf_flig(ivt) * wtcol * leaf_prof_p(j,m)

             ! grain litter nitrogen fluxes
                phenology_to_met_n(j,i) = phenology_to_met_n(j,i) &
                   + grainn_to_food_p(m) * lf_flab(ivt) * wtcol * leaf_prof_p(j,m)
                phenology_to_cel_n(j,i) = phenology_to_cel_n(j,i) &
                   + grainn_to_food_p(m) * lf_fcel(ivt) * wtcol * leaf_prof_p(j,m)
                phenology_to_lig_n(j,i) = phenology_to_lig_n(j,i) &
                   + grainn_to_food_p(m) * lf_flig(ivt) * wtcol * leaf_prof_p(j,m)
!            end if
          end if
       end do  !end pft loop
    end do! end soil level loop

  end subroutine CNLitterToColumn

end module bgc_veg_CNPhenologyMod
