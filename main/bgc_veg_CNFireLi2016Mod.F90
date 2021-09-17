module bgc_veg_CNFireLi2016Mod

use precision
use timemanager
use PhysicalConstants, only: tfrz
use MOD_1D_Forcing, only: &
    forc_q, forc_t, forc_psrf, forc_us, forc_vs
use MOD_TimeInvariants, only: &
    i_cwd, occur_hi_gdp_tree, isshrub, isgrass, isbetr, isbdtr, isbare, iscrop, isnatveg, gdp_lf, abm_lf, peatf_lf, &
    lfuel, ufuel, cropfire_a1, borealat, troplat, non_boreal_peatfire_c, boreal_peatfire_c, rh_low, rh_hgh, &
    bt_min, bt_max, pot_hmn_ign_counts_alpha, g0, fd_pft, fsr_pft, psi0, porsl, bsw, rootfr
use MOD_TimeVariables, only: &
    leafc     , leafc_storage     , leafc_xfer     , frootc    , frootc_storage    , frootc_xfer    , &
    deadcrootc, deadcrootc_storage, deadcrootc_xfer, livecrootc, livecrootc_storage, livecrootc_xfer, &
    decomp_cpools_vr,  totlitc,  totvegc,  cropf   , lfwt      , fuelc             , fuelc_crop     , fsr            , &
    fd        , rootc             , lgdp           , lgdp1     , lpop              , wtlf           , &
    trotr1    , trotr2            , hdmlf          , lnfm      , baf_crop          , baf_peatf      , &
    burndate  , farea_burned      , nfire          , fsat      , prec60            , wf2            , &
    tsoi17    , rh30              , t_soisno       , wliq_soisno     


implicit none

public CNFireArea

contains

subroutine CNFireArea(i,dlat,nl_soil,ivt,idate,dz_soi)

integer ,intent(in) :: i
real(r8),intent(in) :: dlat
integer ,intent(in) :: nl_soil
integer ,intent(in) :: ivt
integer ,intent(in) :: idate(3)
real(r8),intent(in) :: dz_soi(1:nl_soil)

integer  :: g,l,c,p,j,fc,fp,kyr, kmo, kda, mcsec   ! index variables
real(r8) :: dayspyr  ! days per year
real(r8) :: fb       ! availability of fuel for regs A and C
real(r8) :: fhd      ! impact of hd on agricultural fire
real(r8) :: fgdp     ! impact of gdp on agricultural fire
real(r8) :: fire_m   ! combustability of fuel for fire occurrence
real(r8) :: spread_m ! combustability of fuel for fire spread
real(r8) :: Lb_lf    ! length-to-breadth ratio added by Lifang
real(r8) :: lh       ! anthro. ignitions (count/km2/hr)
real(r8) :: fs       ! hd-dependent fires suppression (0-1)
real(r8) :: ig       ! total ignitions (count/km2/hr)
real(r8) :: arh, arh30 !combustability of fuel related to RH and RH30
real(r8) :: afuel    !weight for arh and arh30
real(r8) :: eq
real(r8) :: deqdT
real(r8) :: qsatq
real(r8) :: qsatqdT
real(r8) :: forc_rh
real(r8) :: rootr(nl_soil)
real(r8) :: rresis(nl_soil)
real(r8) :: smp_node
real(r8) :: s_node
real(r8) :: tmp1d(nl_soil)
real(r8) :: tmp0d
real(r8) :: btran2

!logical  :: transient_landcover  ! whether this run has any prescribed transient landcover
!real(r8), target  :: prec60_col_target(bounds%begc:bounds%endc)
!real(r8), target  :: prec10_col_target(bounds%begc:bounds%endc)
!real(r8), target  :: rh30_col_target(bounds%begc:bounds%endc)
!real(r8) :: prec60(:)
!real(r8) :: rh30  (:)

real(r8),parameter :: secsphr = 3600._r8
real(r8),parameter :: secspday = 86400._r8
real(r8),parameter :: PI = 4.*atan(1.)

     tsoi17 = forc_t(i)  ! Temporarily use air temperature for tsoi17, need to revised later.
     prec60 = 0          ! Temporarily use 0 for prec60, need to revised later
     wf2    = 0.5        ! Temporarily set up, need to revise later.
     fsat   = 0          ! Temporarily set up, need to revise later.
     rh30   = 0          ! Temporarily set up, need to revise later.
     hdmlf  = 0._r8      ! Temporarily set up, need to revise later.
     gdp_lf = 0._r8      ! Temporarily set up, need to revise later.
     abm_lf = 8
     peatf_lf = 0._r8
     lnfm   = 0._r8      ! Temporarily set up, need to revise later.

     call julian2monthday(idate(1),idate(2),kmo,kda)

     call eroot(nl_soil,0._r8,porsl(1:,i),bsw(1:,i),psi0(1:,i),rootfr(1:,i),dz_soi(1:),&
          t_soisno(1:,i),wliq_soisno(1:,i),tmp1d,tmp0d,btran2)
     !
     ! Calculate fraction of crop (cropf_col) and non-crop and non-bare-soil 
     ! vegetation (lfwt) in vegetated column
     !
     cropf(i) = 0._r8 
     lfwt (i) = 0._r8   

     ! For crop veg types
     if( iscrop(ivt) )then
        cropf(i) = 1._r8
     end if
     ! For natural vegetation (non-crop and non-bare-soil)
     if( isnatveg(ivt))then
        lfwt (i) = 1._r8
     end if

     ! 
     ! Calculate crop fuel   
     !
     fuelc_crop(i)=0._r8

     ! For crop PFTs, fuel load includes leaf and litter; only
     ! column-level litter carbon
     ! is available, so we use leaf carbon to estimate the
     ! litter carbon for crop PFTs
     if( iscrop(ivt) .and. leafc(i) > 0._r8 )then
         fuelc_crop(i)=fuelc_crop(i) + leafc(i) + leafc_storage(i) + leafc_xfer(i) +  totlitc(i)
     end if
     !   
     ! Calculate noncrop column variables
     !
     fsr   (i) = 0._r8
     fd    (i) = 0._r8
     rootc (i) = 0._r8
     lgdp  (i) = 0._r8
     lgdp1 (i) = 0._r8
     lpop  (i) = 0._r8
!     btran (i) = 0._r8
     wtlf  (i) = 0._r8
     trotr1(i) = 0._r8
     trotr2(i) = 0._r8

!        if (transient_landcover) then
!     dtrotr_col(c)=0._r8
!        end if

     ! For non-crop -- natural vegetation and bare-soil
     if( isnatveg(ivt) .or. isbare(ivt) )then
!        if( .not. shr_infnan_isnan(btran2(p))) then
        if (btran2  <=  1._r8 ) then
           wtlf(i)      = 1._r8
        end if
!                 end if

                 ! NOTE(wjs, 2016-12-15) These calculations of the fraction of evergreen
                 ! and deciduous tropical trees (used to determine if a column is
                 ! tropical closed forest) use the current fractions. However, I think
                 ! they are used in code that applies to land cover change. Note that
                 ! land cover change is currently generated on the first time step of the
                 ! year (even though the fire code sees the annually-smoothed dwt). Thus,
                 ! I think that, for this to be totally consistent, this code should
                 ! consider the fractional coverage of each PFT prior to the relevant
                 ! land cover change event. (These fractions could be computed in the
                 ! code that handles land cover change, so that the fire code remains
                 ! agnostic to exactly how and when land cover change happens.)
                 !
                 ! For example, if a year started with fractional coverages of
                 ! nbrdlf_evr_trp_tree = 0.35 and nbrdlf_dcd_trp_tree = 0.35, but then
                 ! the start-of-year land cover change reduced both of these to 0.2: The
                 ! current code would consider the column to NOT be tropical closed
                 ! forest (because nbrdlf_evr_trp_tree+nbrdlf_dcd_trp_tree < 0.6),
                 ! whereas in fact the land cover change occurred when the column *was*
                 ! tropical closed forest.
        if( isbetr(ivt) )then
           trotr1(i)=1._r8
        end if
        if( isbdtr(ivt) .and. abs(dlat) .lt. troplat)then
           trotr2(i)=1._r8
        end if

!                 if (transient_landcover) then
!        if( ivt == nbrdlf_evr_trp_tree .or. ivt == nbrdlf_dcd_trp_tree )then
!           if(dwt_smoothed(p) < 0._r8)then
                          ! Land cover change in CLM happens all at once on the first time
                          ! step of the year. However, the fire code needs deforestation
                          ! rates throughout the year, in order to combine these
                          ! deforestation rates with the current season's climate. So we
                          ! use a smoothed version of dwt.
                          !
                          ! This isn't ideal, because the carbon stocks that the fire code
                          ! is operating on will have decreased by the full annual amount
                          ! before the fire code does anything. But the biggest effect of
                          ! these deforestation fires is as a trigger for other fires, and
                          ! the C fluxes are merely diagnostic so don't need to be
                          ! conservative, so this isn't a big issue.
                          !
                          ! (Actually, it would be even better if the fire code had a
                          ! realistic breakdown of annual deforestation into the
                          ! different seasons. But having deforestation spread evenly
                          ! throughout the year is much better than having it all
                          ! concentrated on January 1.)
!               dtrotr_col(c)=dtrotr_col(c)-dwt_smoothed(p)
!           end if
!        end if
!                end if
!                 if (spinup_state == 2) then         
!                    rootc_col(c) = rootc_col(c) + (frootc(p) + frootc_storage(p) + &
!                         frootc_xfer(p) + deadcrootc(p) * 10._r8 +       &
!                         deadcrootc_storage(p) + deadcrootc_xfer(p) +    &
!                         livecrootc(p)+livecrootc_storage(p) +           &
!                         livecrootc_xfer(p))*patch%wtcol(p)
!                 else
        rootc(i) = rootc(i) + frootc(i) + frootc_storage(i) + &
                         frootc_xfer(i) + deadcrootc(i) +                &
                         deadcrootc_storage(i) + deadcrootc_xfer(i) +    &
                         livecrootc(i)+livecrootc_storage(i) +           &
                         livecrootc_xfer(i)
!                 endif

        fsr(i) = fsr_pft(ivt)

!        hdmlf=this%forc_hdm(i)

        ! all these constants are in Li et al. BG (2012a,b;2013)

        if( hdmlf(i)  >  0.1_r8 )then           
           ! For NOT bare-soil
           if(.not. isbare(ivt) )then
              ! For shrub and grass (crop already excluded above)
              if( isshrub(ivt) .or. isgrass(ivt) )then      !for shurb and grass
                 lgdp(i)  = lgdp(i) + (0.1_r8 + 0.9_r8*    &
                                exp(-1._r8*PI* &
                                (gdp_lf(i)/8._r8)**0.5_r8))/(1.0_r8-cropf(i))
                 lgdp1(i) = lgdp1(i) + (0.2_r8 + 0.8_r8*   &
                                exp(-1._r8*PI* &
                                (gdp_lf(i)/7._r8)))/(1._r8-cropf(i))
                 lpop(i)  = lpop(i) + (0.2_r8 + 0.8_r8*    &
                                exp(-1._r8*PI* &
                                (hdmlf(i)/450._r8)**0.5_r8))/(1._r8-cropf(i))
              else   ! for trees
                 if( gdp_lf(i)  >  20._r8 )then
                    lgdp(i) =lgdp(i)+occur_hi_gdp_tree/(1._r8-cropf(i))
                    lgdp1(i) =lgdp1(i)+0.62_r8/(1._r8-cropf(i))
                 else
                    if( gdp_lf(i) > 8._r8 )then
                       lgdp(i)=lgdp(i)+0.79_r8/(1._r8-cropf(i))
                       lgdp1(i)=lgdp1(i)+0.83_r8/(1._r8-cropf(i))
                    else
                       lgdp(i) = lgdp(i)+1._r8/(1._r8-cropf(i))
                       lgdp1(i)=lgdp1(i)+1._r8/(1._r8-cropf(i))
                    end if
                 end if
                 lpop(i) = lpop(i) + (0.4_r8 + 0.6_r8*    &
                                  exp(-1._r8*PI* &
                                  (hdmlf(i)/125._r8)))/(1._r8-cropf(i))
              end if
           end if
        else
           lgdp(i)  = lgdp(i)  + 1._r8/(1._r8-cropf(i))
           lgdp1(i) = lgdp1(i) + 1._r8/(1._r8-cropf(i))
           lpop(i)  = lpop(i)  + 1._r8/(1._r8-cropf(i))
        end if

        fd(i) = fd_pft(ivt) * secsphr / (1.0_r8-cropf(i))         
     end if
     ! estimate annual decreased fractional coverage of BET+BDT
     ! land cover conversion in CLM4.5 is the same for each timestep except for the beginning  

!     if (transient_landcover) then
!     if( dtrotr_col(c)  >  0._r8 )then
!        if( kmo == 1 .and. kda == 1 .and. mcsec == 0)then
!           lfc(c) = 0._r8
!        end if
!        if( kmo == 1 .and. kda == 1 .and. mcsec == dt)then
!           lfc(c) = dtrotr_col(c)*dayspyr*secspday/dt
!        end if
!     else
!        lfc(c)=0._r8
!     end if
!     end if
     !
     ! calculate burned area fraction in cropland
     !
     baf_crop(i)=0._r8

     if( kmo == 1 .and. kda == 1 .and. idate(3) == 0 )then
        burndate(i) = 10000 ! init. value; actual range [0 365]
     end if

     ! For crop
     if( forc_t(i)  >=  tfrz .and. iscrop(ivt) .and.  &
        kmo == abm_lf(i) .and. burndate(i) >= 999)then ! catch  crop burn time

        ! calculate human density impact on ag. fire
        fhd = 0.04_r8+0.96_r8*exp(-1._r8*PI*(hdmlf(i)/350._r8)**0.5_r8)

        ! calculate impact of GDP on ag. fire
        fgdp = 0.01_r8+0.99_r8*exp(-1._r8*PI*(gdp_lf(i)/10._r8))

        ! calculate burned area
        fb   = max(0.0_r8,min(1.0_r8,(fuelc_crop(i)-lfuel)/(ufuel-lfuel)))

        ! crop fire only for generic crop types at this time
        ! managed crops are treated as grasses if crop model is turned on
        baf_crop(i) = baf_crop(i) + cropfire_a1/secsphr*fhd*fgdp
        if( fb*fhd*fgdp  >  0._r8)then
           burndate(i)  =  kda
        end if
     end if

     !
     ! calculate peatland fire
     !
!     do fc = 1, num_soilc
!        c = filter_soilc(fc)
!        g= col%gridcell(c)
     if(dlat < borealat )then
        baf_peatf(i) = non_boreal_peatfire_c/secsphr*max(0._r8, &
                min(1._r8,(4.0_r8-prec60(i)*secspday)/ &
                4.0_r8))**2*peatf_lf(i)*(1._r8-fsat(i))
     else
        baf_peatf(i) = boreal_peatfire_c/secsphr*exp(-PI*(max(wf2(i),0._r8)/0.3_r8))* &
                max(0._r8,min(1._r8,(tsoi17(i)-tfrz)/10._r8))*peatf_lf(i)* &
                (1._r8-fsat(i))
     end if
!     end do
     !
     ! calculate other fires
     !

     ! Set the number of timesteps for e-folding.
     ! When the simulation has run fewer than this number of steps,
     ! re-scale the e-folding time to get a stable early estimate.

     ! find which pool is the cwd pool
!     i_cwd = 0
!     do l = 1, ndecomp_pools
!        if ( is_cwd(l) ) then
!           i_cwd = l
!        endif
!     end do

     call qsadv(forc_t,forc_psrf,eq,deqdT,qsatq,qsatqdT)
     forc_rh = forc_q(i) / eq

     if( cropf(i)  <  1._r8 )then
        fuelc(i) = totlitc(i)+totvegc(i)-rootc(i)-fuelc_crop(i)*cropf(i)
!        if (spinup_state == 2) then
!           fuelc(c) = fuelc(c) + ((10._r8 - 1._r8)*deadstemc_col(c))
!           do j = 1, nlevdecomp  
!              fuelc(c) = fuelc(c)+decomp_cpools_vr(c,j,i_cwd) * dzsoi_decomp(j) * spinup_factor(i_cwd) &
!                         * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
!           end do
!        else
        do j = 1, nl_soil
           fuelc(i) = fuelc(i)+decomp_cpools_vr(j,i_cwd,i) * dz_soi(j)
        end do
!        end if
        fuelc(i) = fuelc(i)/(1._r8-cropf(i))
        fb       = max(0.0_r8,min(1.0_r8,(fuelc(i)-lfuel)/(ufuel-lfuel)))
        if (trotr1(i)+trotr2(i)<=0.6_r8) then  
           afuel  =min(1._r8,max(0._r8,(fuelc(i)-2500._r8)/(5000._r8-2500._r8)))
           arh=1._r8-max(0._r8, min(1._r8,(forc_rh-rh_low)/(rh_hgh-rh_low)))
           arh30=1._r8-max(0.7_r8, min(1._r8,rh30(i)/90._r8))
           if (forc_rh < rh_hgh.and. wtlf(i) > 0._r8 .and. tsoi17(i)> tfrz)then
              fire_m   = ((afuel*arh30+(1._r8-afuel)*arh)**1.5_r8)*((1._r8 -max(0._r8,&
                 min(1._r8,(btran2/wtlf(i)-bt_min)/(bt_max-bt_min))))**0.5_r8)
           else
              fire_m   = 0._r8
           end if
           lh       = pot_hmn_ign_counts_alpha*6.8_r8*hdmlf(i)**(0.43_r8)/30._r8/24._r8
           fs       = 1._r8-(0.01_r8+0.98_r8*exp(-0.025_r8*hdmlf(i)))
           ig       = (lh+lnfm(i)/(5.16_r8+2.16_r8*cos(PI/180._r8*3*min(60._r8,abs(dlat/PI*180))))*0.22_r8)  &
                      *(1._r8-fs)*(1._r8-cropf(i))
           nfire(i) = ig/secsphr*fb*fire_m*lgdp(i) !fire counts/km2/sec
           Lb_lf    = 1._r8+10._r8*(1._r8-EXP(-0.06_r8*sqrt(forc_us(i)*forc_us(i)+forc_vs(i)*forc_vs(i))))
           spread_m = fire_m**0.5_r8
           farea_burned(i) = min(1._r8,(g0*spread_m*fsr(i)* &
                   fd(i)/1000._r8)**2*lgdp1(i)* &
                   lpop(i)*nfire(i)*PI*Lb_lf+ &
                   baf_crop(i)+baf_peatf(i))  ! fraction (0-1) per sec
        else
           farea_burned(i)=min(1._r8,baf_crop(i)+baf_peatf(i))
        end if
           !
           ! if landuse change data is used, calculate deforestation fires and 
           ! add it in the total of burned area fraction
           !
!           if (transient_landcover) then
!        if( trotr1_col(c)+trotr2_col(c) > 0.6_r8 )then
!           if(( kmo == 1 .and. kda == 1 .and. mcsec == 0) .or. dtrotr_col(c) <=0._r8 )then
!              fbac1(c)        = 0._r8
!              farea_burned(c) = baf_crop(c)+baf_peatf(c)
!           else
!              cri = (4.0_r8*trotr1_col(c)+1.8_r8*trotr2_col(c))/(trotr1_col(c)+trotr2_col(c))
!              cli = (max(0._r8,min(1._r8,(cri-prec60_col(c)*secspday)/cri))**0.5)* &
!                    (max(0._r8,min(1._r8,(cri-prec10_col(c)*secspday)/cri))**0.5)* &
!                         max(0.0005_r8,min(1._r8,19._r8*dtrotr_col(c)*dayspyr*secspday/dt-0.001_r8))* &
!                         max(0._r8,min(1._r8,(0.25_r8-(forc_rain(c)+forc_snow(c))*secsphr)/0.25_r8))
!              farea_burned(c) = cli*(cli_scale/secspday)+baf_crop(c)+baf_peatf(c)
!              ! burned area out of conversion region due to land use fire
!              fbac1(c) = max(0._r8,fb*cli*(cli_scale/secspday) - 2.0_r8*lfc(c)/dt)   
!           end if
!           ! total burned area out of conversion 
!           fbac(c) = fbac1(c)+baf_crop(c)+baf_peatf(c) 
!        else
!           fbac(c) = farea_burned(c)
!        end if
     else
        farea_burned(i) = min(1._r8,baf_crop(i)+baf_peatf(i))
     end if


 end subroutine CNFireArea

end module bgc_veg_CNFireLi2016Mod
