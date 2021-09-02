!-----------------------------------------------------------------------
!BOP
!
!IROUTINE: ThreeDCanopy
!
! !INTERFACE:

SUBROUTINE ThreeDCanopy_wrap (ipatch, czen, albg, albv, ssun, ssha)

   USE precision
   USE GlobalVars
   USE PFT_Const
   USE MOD_PFTimeInvars
   USE MOD_PCTimeInvars
   USE MOD_PFTimeVars
   USE MOD_PCTimeVars

   IMPLICIT NONE

   INTEGER,  intent(in)  :: ipatch
   REAL(r8), Intent(in)  :: czen
   REAL(r8), Intent(in)  :: albg(2,2)
   REAL(r8), Intent(out) :: albv(2,2)
   REAL(r8), Intent(out) :: ssun(2,2)
   REAL(r8), Intent(out) :: ssha(2,2)

   ! local variables
   REAL(r8), dimension(0:N_PFT-1, 2) :: albd, albi, fabd, fabi, sun_fadd
   REAL(r8), dimension(0:N_PFT-1, 2) :: ftdd, ftid, ftii
   REAL(r8), dimension(0:N_PFT-1, 2) :: rho, tau
   REAL(r8), dimension(0:N_PFT-1)    :: csiz, chgt, lsai 
   REAL(r8), dimension(0:N_PFT-1)    :: fsun_id, fsun_ii, psun
   REAL(r8), dimension(0:N_PFT-1)    :: phi1, phi2, gdir

   INTEGER p, pc
  
   ! get PC patch index
   pc = patch2pc(ipatch)      

   ! initialization
   albd=1.; albi=1.; fabd=0.; fabi=0.; 
   ftdd=1.; ftid=0.; ftii=1.; sun_fadd=0.;
   csiz(:) = (htop_c(:,pc) - hbot_c(:,pc)) / 2
   chgt(:) = (htop_c(:,pc) + hbot_c(:,pc)) / 2
   lsai(:) = lai_c(:,pc) + sai_c(:,pc)
   
   ! calculate weighted plant optical properties
   ! loop for each PFT
   rho = 0.
   tau = 0.
   DO p = 0, N_PFT-1
      IF (lsai(p) > 0.) THEN
         rho(p,:) = rho_p(:,1,p)*lai_c(p,pc)/lsai(p) &
                  + rho_p(:,2,p)*sai_c(p,pc)/lsai(p) 
         tau(p,:) = tau_p(:,1,p)*lai_c(p,pc)/lsai(p) &
                  + tau_p(:,2,p)*sai_c(p,pc)/lsai(p) 
      ENDIF
   ENDDO 
   
   CALL ThreeDCanopy(N_PFT, canlay(:), pcfrac(:,pc), csiz, chgt, chil_p(:), czen, &
                     lsai, rho, tau, albg(:,1), albg(:,2), albd, albi, &
                     fabd, fabi, ftdd, ftid, ftii, sun_fadd, psun, &
                     thermk_c(:,pc), fshade_c(:,pc) )

   ! calculate extkb_c, extkd_c
   ! applied for 1D case
   extkd_c(:,pc) = 0.719 !used for scaling-up coefficients from leaf to canopy

   ! 11/07/2018: calculate gee FUNCTION consider LAD
   phi1 = 0.5 - 0.633 * chil_p(:) - 0.33 * chil_p(:) * chil_p(:)
   phi2 = 0.877 * ( 1. - 2. * phi1 )

   ! 11/07/2018: calculate gee FUNCTION consider LAD
   gdir = phi1 + phi2*czen
   extkb_c(:,pc) = gdir/czen

   fsun_id(:) = 0.
   fsun_ii(:) = 0.

   ! lai还是lsai? 应该是lsai. 吸收的辐射都算在叶子上? curr. no.
   ! par for lai or lai+sai. curr. lai+sai
   ! 如果lai比较小，sai比较大时，会不会造成虚假par吸收?
   ! 一般来讲，叶子在阳面的概率会更高，但并不能把所有阳面
   ! 的吸收都算在叶子上
   DO p = 1, N_PFT-1
      IF (lsai(p) > 0.) THEN
         fsun_id(p) = (1._r8 - exp(-2._r8*extkb_c(p,pc)*lsai(p))) / &
            (1._r8 - exp(-extkb_c(p,pc)*lsai(p))) / 2.0_r8 * psun(p)    

         fsun_ii(p) = (1._r8 - exp(-extkb_c(p,pc)*lsai(p)-0.5/0.5_r8*lsai(p))) / &
            (extkb_c(p,pc)+0.5/0.5_r8) / &
            (1._r8 - exp(-0.5/0.5_r8*lsai(p))) *  &
            (0.5/0.5_r8) * psun(p)
      ENDIF 
   ENDDO 

   ! calculate albv, ssun, ssha
   ! NOTE: CoLM (1/2,): vis/nir; (,1/2): dir/dif
   albv(1,1) = albd(1,1); albv(1,2) = albi(1,1)
   albv(2,1) = albd(1,2); albv(2,2) = albi(1,2)
      
   ! ssun(band, dir/dif, pft), fabd/sun_fadd(pft, band)
   ssun_c(1,1,:,pc) = sun_fadd(:,1) + (fabd(:,1)-sun_fadd(:,1))*fsun_id
   ssun_c(2,1,:,pc) = sun_fadd(:,2) + (fabd(:,2)-sun_fadd(:,2))*fsun_id
   ssha_c(1,1,:,pc) = (fabd(:,1)-sun_fadd(:,1)) * (1.-fsun_id)
   ssha_c(2,1,:,pc) = (fabd(:,2)-sun_fadd(:,2)) * (1.-fsun_id)
   ssun_c(1,2,:,pc) = fabi(:,1) * fsun_ii
   ssun_c(2,2,:,pc) = fabi(:,2) * fsun_ii
   ssha_c(1,2,:,pc) = fabi(:,1) * (1.-fsun_ii)
   ssha_c(2,2,:,pc) = fabi(:,2) * (1.-fsun_ii)

   ssun(1,1) = sum( ssun_c(1,1,:,pc) * pcfrac(:,pc) )
   ssun(2,1) = sum( ssun_c(2,1,:,pc) * pcfrac(:,pc) )
   ssun(1,2) = sum( ssun_c(1,2,:,pc) * pcfrac(:,pc) )
   ssun(2,2) = sum( ssun_c(2,2,:,pc) * pcfrac(:,pc) )

   ssha(1,1) = sum( ssha_c(1,1,:,pc) * pcfrac(:,pc) )
   ssha(2,1) = sum( ssha_c(2,1,:,pc) * pcfrac(:,pc) )
   ssha(1,2) = sum( ssha_c(1,2,:,pc) * pcfrac(:,pc) )
   ssha(2,2) = sum( ssha_c(2,2,:,pc) * pcfrac(:,pc) )

END SUBROUTINE ThreeDCanopy_wrap


SUBROUTINE ThreeDCanopy(npft, canlev, pwtcol, csiz, chgt, chil, coszen, &
                        lsai, rho, tau, albgrd, albgri, albd, albi, &
                        fabd, fabi, ftdd, ftid, ftii, sun_fadd, psun, &
                        thermk, fshade)

!
! !DESCRIPTION:
! ThreeDCanopy based on Dickinson (2008) using three canopy layer
! to calculate fluxes absorbed by vegetation, reflected by vegetation,
! and transmitted through vegetation for unit incoming direct or
! diffuse flux given an underlying surface with known albedo.

!
! !ARGUMENTS:
   IMPLICIT NONE

   INTEGER, parameter :: r8 = selected_real_kind(12)
   INTEGER, parameter :: numrad = 2

! !ARGUMENTS:
   INTEGER , intent(in) :: npft                       !number of vegetated pfts where coszen>0
   INTEGER , intent(in) :: canlev(0:npft-1)           !canopy level for current pft
   REAL(r8), intent(in) :: pwtcol(0:npft-1)           !weight of pft wrt corresponding column
   REAL(r8), intent(in) :: csiz  (0:npft-1)           !
   REAL(r8), intent(in) :: chgt  (0:npft-1)           !
   REAL(r8), intent(in) :: chil  (0:npft-1)           !leaf angle distribution parameter
   REAL(r8), intent(in) :: lsai  (0:npft-1)           !
   REAL(r8), intent(in) :: rho   (0:npft-1,numrad)    !leaf/stem refl weighted by fraction LAI and SAI
   REAL(r8), intent(in) :: tau   (0:npft-1,numrad)    !leaf/stem tran weighted by fraction LAI and SAI

   REAL(r8), intent(in) :: coszen                     !cosine solar zenith angle for next time step
   REAL(r8), intent(in) :: albgrd(numrad)             !ground albedo (direct) (column-level)
   REAL(r8), intent(in) :: albgri(numrad)             !ground albedo (diffuse)(column-level)

   REAL(r8), intent(out) :: albd(0:npft-1,numrad)     !surface albedo (direct)
   REAL(r8), intent(out) :: albi(0:npft-1,numrad)     !surface albedo (diffuse)
   REAL(r8), intent(out) :: fabd(0:npft-1,numrad)     !flux absorbed by veg per unit direct flux
   REAL(r8), intent(out) :: fabi(0:npft-1,numrad)     !flux absorbed by veg per unit diffuse flux
   REAL(r8), intent(out) :: ftdd(0:npft-1,numrad)     !down direct flux below veg per unit dir flx
   REAL(r8), intent(out) :: ftid(0:npft-1,numrad)     !down diffuse flux below veg per unit dir flx
   REAL(r8), intent(out) :: ftii(0:npft-1,numrad)     !down diffuse flux below veg per unit dif flx
   REAL(r8), intent(out) :: sun_fadd(0:npft-1,numrad) !
   REAL(r8), intent(out) :: psun  (0:npft-1)          !
   REAL(r8), intent(out) :: thermk(0:npft-1)          !
   REAL(r8), intent(out) :: fshade(0:npft-1)          !

! !LOCAL VARIABLES:
   REAL(selected_real_kind(12)), external::OverlapArea
   REAL(selected_real_kind(12)), external::tee

! !OTHER LOCAL VARIABLES:
   INTEGER,  parameter :: r16 = selected_real_kind(24) !16 byte REAL
   REAL(r8), parameter :: mpe = 1.0e-06_r8 !prevents overflow for division by zero
   INTEGER , parameter :: nlay=3           !number of canopy layers
   REAL(r8), parameter :: D0=0.0_r8        !double accuracy REAL number
   REAL(r8), parameter :: D1=1.0_r8        !double accuracy REAL number
   REAL(r8), parameter :: D2=2.0_r8        !double accuracy REAL number
   REAL(r8), parameter :: D3=3.0_r8        !double accuracy REAL number
   REAL(r8), parameter :: D4=4.0_r8        !double accuracy REAL number
   REAL(r8), parameter :: D6=6.0_r8        !double accuracy REAL number
   REAL(r8), parameter :: D7=7.0_r8        !double accuracy REAL number
   REAL(r8), parameter :: D8=8.0_r8        !double accuracy REAL number
   REAL(r8), parameter :: D9=9.0_r8        !double accuracy REAL number
   REAL(r8), parameter :: D10=10.0_r8      !double accuracy REAL number
   REAL(r8), parameter :: D16=16.0_r8      !double accuracy REAL number
   REAL(r8), parameter :: DH=0.5_r8        !quad accuracy REAL number
   REAL(r16),parameter :: DDH=0.5_r16      !quad accuracy REAL number
   REAL(r16),parameter :: DD0=0.0_r16      !quad accuracy REAL number
   REAL(r16),parameter :: DD1=1.0_r16      !quad accuracy REAL number
   REAL(r8) ,parameter :: pi=3.14159265358979323846_r8  !pi

   INTEGER  :: ib                   !band index 1:vis 2:nir
   INTEGER  :: ip,ic,ig,kband       !array indices for pft,column,grid 
   INTEGER  :: kfr                  !variable for layer radiation coming from
   INTEGER  :: klay                 !variable for layer absorbing radiation
   INTEGER  :: kto                  !variable for layer radiation is transmitted to
   INTEGER  :: lev                  !do loop variable
   INTEGER  :: nn                   !do loop variable
   INTEGER  :: nsoilveg             !number of pfts in gridcell with veg and cosz > 0
   INTEGER  :: nstep                !time step index
   INTEGER  :: clev                 !canopy level for current pft

   REAL(r8) :: albd_col(numrad)     !surface reflection (direct) for column
   REAL(r8) :: albi_col(numrad)     !surface reflection (diffuse) for column
   REAL(r8) :: bot_lay(nlay)        !avergae canopy bottom in layer
   REAL(r8) :: hgt_lay(nlay)        !average canopy height in layer
   REAL(r8) :: omg_lay(nlay,numrad) !average omega for all three layer
   REAL(r8) :: rho_lay(nlay,numrad) !average rho for all three layer
   REAL(r8) :: siz_lay(nlay)        !average canopy size in layer
   REAL(r8) :: tau_lay(nlay,numrad) !average tau for all three layer
   REAL(r8) :: lsai_lay(nlay)       !average lsai for each layer
   REAL(r8) :: cosz                 !0.001 <= coszen <= 1.000
   REAL(r8) :: cosd                 !0.001 <= coszen <= 1.000
   REAL(r8) :: delta                !variable for increment layer in loop
   REAL(r8) :: dif                  !diffuse radiation transmitted
   REAL(r8) :: dir                  !direct radiation transmitted
   REAL(r8) :: fabd_col(numrad)     !flux absorbed by veg per unit diffuse flux
   REAL(r8) :: fabd_lay(nlay,numrad)!layer absorption for direct beam
   REAL(r8) :: fabi_col(numrad)     !flux absorbed by veg per unit diffuse flux
   REAL(r8) :: fabi_lay(nlay,numrad)!layer absorption for diffuse beam
   REAL(r8) :: fabs_lay(0:4,numrad) !layer absorption for all five layers 
   REAL(r8) :: fabs_leq(0:4,numrad) !layer absorption for all five layers 
   REAL(r8) :: A(6,6)               !
   REAL(r8) :: B(6,2)               !
   REAL(r8) :: X(6,2)               !
   REAL(r8) :: fabsm                !pft absorption for multiple reflections
   REAL(r8) :: faid_lay(nlay)       !layer diffused absorption for direct beam
   REAL(r8) :: faid_p               !pft absorption direct beam
   REAL(r8) :: faii_lay(nlay)       !layer diffused absorption for diffuse beam
   REAL(r8) :: faii_p               !pft absorption diffuse beam
   REAL(r8) :: fc0(nlay)            !canopy fraction for layers
   REAL(r8) :: frid_lay(nlay)       !layer reflection for direct beam
   REAL(r8) :: frid_p               !pft reflection direct beam
   REAL(r8) :: frii_lay(nlay)       !layer reflection for indirect beam
   REAL(r8) :: ftdd_lay(nlay)       !unscattered layer transmission for direct beam
   REAL(r8) :: ftdi_lay(nlay)       !unscattered layer transmission for indirect beam
   REAL(r8) :: ftdd_lay_orig(nlay)  !unscattered layer transmission for direct beam without lad/crown_shape calibration
   REAL(r8) :: ftdi_lay_orig(nlay)  !unscattered layer transmission for indirect beam without lad/crown_shape calibratioin
   REAL(r8) :: ftid_lay(nlay)       !diffused layer transmission for direct beam
   REAL(r8) :: ftii_lay(nlay)       !diffused layer transmission for diffuse beam
   REAL(r8) :: ftran                !pft transmittance
   REAL(r8) :: gee=0.5_r8           !Ross factor geometric blocking
   REAL(r8) :: gdir(0:npft-1)       !Ross G factor considering LAD for incident direct radiation
   REAL(r8) :: gdif(0:npft-1)       !Ross G factor considering LAD for incident diffuse radiation
   REAL(r8) :: gdir_lay(nlay)       !Ross G factor considering LAD for incident direct radiation
   REAL(r8) :: gdif_lay(nlay)       !Ross G factor considering LAD for incident diffuse radiation
   REAL(r8) :: fcad(0:npft-1)       !calibration factor for LAD for direct radiation
   REAL(r8) :: fcai(0:npft-1)       !calibration factor for LAD for diffuse radiation
   REAL(r8) :: fcad_lay(nlay)       !calibration factor for LAD for direct radiation
   REAL(r8) :: fcai_lay(nlay)       !calibration factor for LAD for diffuse radiation
   REAL(r8) :: pad                  !probabilty function for absorption after two scat
   REAL(r8) :: pai                  !probabilty of asborption for diffuse incident beam
   REAL(r8) :: pfc                  !contribution of current pft in layer
   REAL(r8) :: probm                !prob photon reflect diffusly from grnd reach canopy
   REAL(r8) :: ref(0:nlay+1,0:nlay+1)!radiation reflected between five layers
   REAL(r8) :: fadd_lay(nlay,numrad)!
   REAL(r8) :: shad_oa(nlay,nlay)   !shadow overlaps (direct beam)
   REAL(r8) :: shadow_d(nlay)       !layer shadow for direct beam
   REAL(r8) :: shadow_i(nlay)       !layer shadow for diffuse beam
   REAL(r8) :: sum_fabd(3)          !sum of absorption for all pfts in grid (direct)
   REAL(r8) :: sum_fabi(3)          !sum of absorption for all pfts in grid (diffuse)
   REAL(r8) :: sum_fadd(nlay)       !
   REAL(r8) :: taud_lay(nlay)       !direct transmission for a layer 
   REAL(r8) :: taui_lay(nlay)       !diffuse transmission for a layer
   REAL(r8) :: trd(0:nlay+1,0:nlay+1)!direct radiation transmitted between five layers
   REAL(r8) :: tri(0:4,0:4)         !diffuse radiation transmitted between five layers
   REAL(r8) :: tt(0:4,0:4)          !unscattered direct radiation available at layer
   REAL(r8) :: wl                   !fraction of LAI+SAI that is LAI
   REAL(r8) :: ws                   !fraction of LAI+SAI that is SAI
   REAL(r8) :: zenith               !zenith angle
   REAL(r8) :: ftdd_col             !unscattered column transmission for direct beam

   REAL(r8) :: shadow_pd(0:npft-1)  !sky shadow area
   REAL(r8) :: shadow_pi(0:npft-1)  !sky shadow area
   REAL(r8) :: shadow_sky(0:npft-1) !sky shadow area
   REAL(r8) :: taud(0:npft-1)       !transmission to direct beam
   REAL(r8) :: taui(0:npft-1)       !transmission to diffuse beam
   REAL(r8) :: omega(0:npft-1,numrad)!leaf/stem transmitance weighted by frac veg
   REAL(r8) :: ftdi(0:npft-1,numrad)!leaf/stem transmitance weighted by frac veg
   REAL(r8) :: ftdd_orig(0:npft-1,numrad)!leaf/stem transmitance weighted by frac veg
   REAL(r8) :: ftdi_orig(0:npft-1,numrad)!leaf/stem transmitance weighted by frac veg
   LOGICAL  :: soilveg(0:npft-1)    !true if pft over soil with veg and cosz > 0

   REAL(r8) :: phi1(0:npft-1), phi2(0:npft-1)

   ! 11/07/2018: calculate gee FUNCTION consider LAD
   phi1 = 0.5 - 0.633 * chil - 0.33 * chil * chil
   phi2 = 0.877 * ( 1. - 2. * phi1 )
   
   cosz = coszen
   cosd = cos(60._r8/180._r8*pi)

   ! 11/07/2018: calculate gee FUNCTION consider LAD
   gdir = phi1 + phi2*cosz
   gdif = phi1 + phi2*cosd
   
   nsoilveg = 0

   fc0 = D0
   omg_lay  = D0
   rho_lay  = D0
   tau_lay  = D0
   hgt_lay  = D0
   siz_lay  = D0
   lsai_lay = D0
   gdir_lay = D0
   gdif_lay = D0

   DO ip=0, npft-1
      shadow_sky(ip) = D1

      ! check elai and pft weight are non-zero 
      IF( lsai(ip) > 1.E-6_r8 .and. pwtcol(ip) > D0 ) THEN

         soilveg(ip) = .true.
         nsoilveg = nsoilveg + 1

         clev      = canlev(ip)
         fc0(clev) = fc0(clev) + pwtcol(ip)

         siz_lay(clev)  = siz_lay(clev)  + pwtcol(ip)*csiz(ip)
         hgt_lay(clev)  = hgt_lay(clev)  + pwtcol(ip)*chgt(ip)
         lsai_lay(clev) = lsai_lay(clev) + pwtcol(ip)*lsai(ip)
         gdir_lay(clev) = gdir_lay(clev) + pwtcol(ip)*gdir(ip)
         gdif_lay(clev) = gdif_lay(clev) + pwtcol(ip)*gdif(ip)

         ! set optical properties
         DO ib = 1, numrad
            omega(ip,ib) = rho(ip,ib) + tau(ip,ib)

            ! sum of tau,rho and omega for pfts in a layer
            tau_lay(clev,ib) = tau_lay(clev,ib) + pwtcol(ip)*(tau(ip,ib))
            rho_lay(clev,ib) = rho_lay(clev,ib) + pwtcol(ip)*(rho(ip,ib))
            omg_lay(clev,ib) = omg_lay(clev,ib) + pwtcol(ip)*(omega(ip,ib))

         ENDDO ! ENDDO ib=1, numrad
      ELSE
         soilveg(ip) = .false.
      ENDIF 
   ENDDO ! ENDDO ip

!=============================================================
! layer average of lsai,tau,rho,omega...
!=============================================================

   DO lev=1,3
      IF( fc0(lev) > D0) THEN
         siz_lay(lev)  = max(siz_lay(lev)/fc0(lev),D0)
         hgt_lay(lev)  = max(hgt_lay(lev)/fc0(lev),D0)
         bot_lay(lev)  = hgt_lay(lev)-siz_lay(lev)
         lsai_lay(lev) = max(lsai_lay(lev)/fc0(lev),D0)
         DO ib = 1, numrad
            tau_lay(lev,ib) = max(tau_lay(lev,ib)/fc0(lev),D0)
            rho_lay(lev,ib) = max(rho_lay(lev,ib)/fc0(lev),D0)
            omg_lay(lev,ib) = max(omg_lay(lev,ib)/fc0(lev),D0)
         ENDDO
         gdir_lay(lev) = max(gdir_lay(lev)/fc0(lev),D0)
         gdif_lay(lev) = max(gdif_lay(lev)/fc0(lev),D0)
      ENDIF
   ENDDO ! ENDDO ib

!=============================================================
! layer shadows
!=============================================================

   shadow_d        = D0
   shadow_i        = D0
   DO lev =1,3
      IF( fc0(lev) >  D0 .and. cosz > D0 ) THEN
         shadow_d(lev) = (D1 - exp(-D1*fc0(lev)/cosz))/&
            (D1 - fc0(lev)*exp(-D1/cosz))
         shadow_d(lev) = max(fc0(lev), shadow_d(lev))
         shadow_i(lev) = (D1 - exp(-D1*fc0(lev)/cosd))/&
            (D1 - fc0(lev)*exp(-D1/cosd))
         shadow_i(lev) = max(fc0(lev), shadow_i(lev))
      ENDIF
   ENDDO

!=============================================================
! taud and ftdd for layers
!=============================================================

   taud_lay = D0
   taui_lay = D0
   ftdd_lay = D0
   ftdi_lay = D0
   fcad_lay = D1
   fcai_lay = D1
   ftdd_lay_orig = D0
   ftdi_lay_orig = D0

   DO lev=1,3
      IF( fc0(lev)>D0 .and. lsai_lay(lev)>D0 ) THEN

         taud_lay(lev) = D3/D4*gee*fc0(lev)*lsai_lay(lev)/&
            (cosz*shadow_d(lev))
         taui_lay(lev) = D3/D4*gee*fc0(lev)*lsai_lay(lev)/&
            (cosd*shadow_i(lev))

         ! 11/07/2018: LAD calibration
         ftdd_lay_orig(lev) = tee(DD1*taud_lay(lev))
         ftdi_lay_orig(lev) = tee(DD1*taui_lay(lev))
         
         ! 11/07/2018: gdir/gdif = FUNCTION(xl, cos)
         ftdd_lay(lev) = tee(DD1*taud_lay(lev)/gee*gdir_lay(lev))
         ftdi_lay(lev) = tee(DD1*taui_lay(lev)/gee*gdif_lay(lev))

         ! calibration for chil
         fcad_lay(lev) = (D1-ftdd_lay(lev)) / (D1-ftdd_lay_orig(lev))
         fcai_lay(lev) = (D1-ftdi_lay(lev)) / (D1-ftdi_lay_orig(lev))

      ENDIF
   ENDDO

!=============================================================
! initialize local variables for layers
!=============================================================

   albd_col = D0
   albi_col = D0
   fabd_col = D0
   fabd_lay = D0
   fabi_col = D0
   fabi_lay = D0
   frid_lay = D0
   frii_lay = D0
   tt       = D0

!=============================================================
! projection shadow overlapping fractions
!=============================================================

   zenith = acos(coszen)
   shad_oa(3,2) = fc0(3)*OverlapArea(siz_lay(3),hgt_lay(3)-bot_lay(2),&
      zenith)
   shad_oa(3,1) = fc0(3)*OverlapArea(siz_lay(3),hgt_lay(3)-bot_lay(1),&
      zenith)
   shad_oa(2,1) = fc0(2)*OverlapArea(siz_lay(2),hgt_lay(2)-bot_lay(1),&
      zenith)

   ! for test
   !shad_oa(3,2) = D0
   !shad_oa(3,1) = D0
   !shad_oa(2,1) = D0

!=============================================================
! unscattered direct sunlight available at  each layer
! 0:sky, 1:top 2:middle 3:bottom and 4:ground layer
!=============================================================

   ftdd_col=D0
   tt = D0
   tt(4,3) = shadow_d(3)
   tt(4,3) = min(D1, max(D0, tt(4,3)))
   tt(4,2) = shadow_d(2)*(D1-shadow_d(3)+shad_oa(3,2))
   tt(4,2) = min(1-tt(4,3), max(D0, tt(4,2)))
   tt(4,1) = shadow_d(1)*(D1-(shadow_d(2)-shad_oa(2,1)) &
      - (shadow_d(3)-shad_oa(3,1)) &
      + (shadow_d(2)-shad_oa(2,1))*(shadow_d(3)-shad_oa(3,2)))
   tt(4,1) = min(1-tt(4,3)-tt(4,2), max(D0, tt(4,1)))

   tt(4,0) = D1-(shadow_d(1)+shadow_d(2)+shadow_d(3) &
      - (shadow_d(2)-shad_oa(2,1))*shadow_d(1) &
      - (shadow_d(3)-shad_oa(3,2))*shadow_d(2) &
      - (shadow_d(3)-shad_oa(3,1))*shadow_d(1) &
      + (shadow_d(2)-shad_oa(2,1))*(shadow_d(3)-shad_oa(3,2))*shadow_d(1))
   tt(4,0) = min(1-tt(4,3)-tt(4,2)-tt(4,1), max(D0, tt(4,0)))

   IF (tt(4,0) < 0) THEN
      print *, abs(tt(4,0))
   ENDIF

   ! direct sunlight passing through top canopy layer
   IF (shadow_d(3) > 0) THEN
      tt(3,2) = shadow_d(2)*(shadow_d(3)-shad_oa(3,2))
      tt(3,2) = min(shadow_d(3), max(D0, tt(3,2)))
      tt(3,1) = shadow_d(1)*(shadow_d(3)-shad_oa(3,1)- &
         (shadow_d(3)-shad_oa(3,2))*(shadow_d(2)-shad_oa(2,1)))
      tt(3,1) = min(shadow_d(3)-tt(3,2), max(D0, tt(3,1)))
      tt(3,0) = shadow_d(3)-tt(3,2)-tt(3,1)

      tt(3,2) = tt(3,2)*ftdd_lay(3)
      tt(3,1) = tt(3,1)*ftdd_lay(3)
      tt(3,0) = tt(3,0)*ftdd_lay(3)
   ENDIF

   ! direct sunlight passing through middle canopy layer
   IF (shadow_d(2) > 0) THEN
      tt(2,1) = shadow_d(1)*(shadow_d(2)-shad_oa(2,1))
      tt(2,1) = min(shadow_d(2), max(D0, tt(2,1)))
      tt(2,0) = shadow_d(2)-tt(2,1)

      tt(2,1) = tt(2,1)*ftdd_lay(2)*(tt(4,2) + tt(3,2))/shadow_d(2)
      tt(2,0) = tt(2,0)*ftdd_lay(2)*(tt(4,2) + tt(3,2))/shadow_d(2)
   ENDIF

   ! direct sunlight passing through third canopy layer
   IF (shadow_d(1) > 0)  THEN
      tt(1,0) = ftdd_lay(1)*(tt(4,1) + tt(3,1) + tt(2,1))!*shadow_d(1)/shadow_d(1)
   ENDIF

!=============================================
! Aggregate direct radiation to layers
!=============================================

   tt(4,3) = tt(4,3)
   tt(3,2) = tt(4,2) + tt(3,2)
   tt(2,1) = tt(4,1) + tt(3,1) + tt(2,1)
   tt(1,0) = tt(4,0) + tt(3,0) + tt(2,0) + tt(1,0)
   ftdd_col = tt(1,0)

   tt(0:4,4) = D0
   tt(0:3,3) = D0
   tt(4:4,2) = D0; tt(0:2,2) = D0
   tt(3:4,1) = D0; tt(0:1,1) = D0
   tt(2:4,0) = D0; tt(0:0,0) = D0


!=======================================
! start radiation beam loop
! ib=1:visible band 2:nir band
!=======================================

   DO ib = 1, numrad

   !===============================
   ! get pft level tau and ftdd
   !===============================

      ! 10/12/2017
      ftdi(:,ib) = D1

      DO ip=0, npft-1

         taud(ip)=D0
         taui(ip)=D0
         shadow_pd(ip)=D0
         shadow_pi(ip)=D0

         IF( soilveg(ip) ) THEN
            clev = canlev(ip)

         !================================================
         ! fractional contribution of current pft in layer
         !================================================

            pfc = min( pwtcol(ip)/fc0(clev), D1)
            shadow_pd(ip)=pfc*shadow_d(clev)
            shadow_pi(ip)=pfc*shadow_i(clev)

         !=====================================
         ! get taud,taui at pft level
         !=====================================

            taud(ip)=D3/D4*gee*pwtcol(ip)*(lsai(ip))/&
               (cosz*shadow_pd(ip))

            taui(ip)=D3/D4*gee*pwtcol(ip)*(lsai(ip))/&
               (cosd*shadow_pi(ip))

         !====================================
         ! transmission at pft level
         !====================================

            ftdd_orig(ip,ib) = tee(DD1*taud(ip))
            ftdi_orig(ip,ib) = tee(DD1*taui(ip))
         
            ! 11/07/2018: gdir/gdif = FUNCTION(xl, cos)
            ftdd(ip,ib) = tee(DD1*taud(ip)/gee*gdir(ip))
            ftdi(ip,ib) = tee(DD1*taui(ip)/gee*gdif(ip))

            ! calibration for chil
            fcad(ip) = (D1-ftdd(ip,ib)) / (D1-ftdd_orig(ip,ib))
            fcai(ip) = (D1-ftdi(ip,ib)) / (D1-ftdi_orig(ip,ib))

         !========================================
         ! sum transmissions of pfts in layer
         !========================================

         ! Deleted by Yuan, 06/03/2012
         !sum_ftdd(clev) = sum_ftdd(clev) + pwtcol(ip)* &
         !     shadow_d(clev)*(D1-ftdd(ip,ib))
         ! 
         !sum_ftdi(clev) = sum_ftdi(clev) + pwtcol(ip)* &
         !     shadow_i(clev)*(D1-ftdi(ip,ib))

         ENDIF ! ENDIF soilveg
      ENDDO ! ENDDO ip

   !==========================================
   ! adjust ftdd using factor (1-sum(ftdd))
   !==========================================

      ! Deleted by Yuan, 06/03/2012
      !DO ip=pfti(ic),pftf(ic)
      !   IF( soilveg(ip)) THEN
      !      
      !      clev = canlev(ip)
      !      ftdd(ip,ib) = ftdd(ip,ib)*(D1 - sum_ftdd(clev))
      !      ftdi(ip,ib) = ftdi(ip,ib)*(D1 - sum_ftdi(clev))
      !      
      !   ENDIF ! ENDIF canlev
      !ENDDO ! ENDDO ip

   !===============================================================
   ! absorption, reflection and transmittance for three canopy layer
   ! using average optical properties of layers 
   ! subroutine CanopyRad calculates fluxes for unit input radiation
   !===============================================================

      ftid_lay=D0; ftii_lay=D1
      frid_lay=D0; frii_lay=D0
      faid_lay=D0; faii_lay=D0

      DO lev=1, 3
         IF( shadow_d(lev) > D0 ) THEN
            CALL CanopyRad(taud_lay(lev), taui_lay(lev), ftdd_lay_orig(lev),&
               ftdi_lay_orig(lev), cosz, cosd, shadow_d(lev), shadow_i(lev), &
               fc0(lev), omg_lay(lev,ib), lsai_lay(lev), &
               tau_lay(lev,ib), rho_lay(lev,ib), ftid_lay(lev), &
               ftii_lay(lev), frid_lay(lev), frii_lay(lev),&
               faid_lay(lev), faii_lay(lev))
         ENDIF
      ENDDO ! ENDDO lev

      ! 11/07/2018: calibration for LAD
      ftid_lay(:) = fcad_lay(:)*ftid_lay(:)
      ftii_lay(:) = fcai_lay(:)*(ftii_lay(:)-ftdi_lay_orig(:)) + ftdi_lay(:)
      frid_lay(:) = fcad_lay(:)*frid_lay(:)
      frii_lay(:) = fcai_lay(:)*frii_lay(:)
      faid_lay(:) = fcad_lay(:)*faid_lay(:)
      faii_lay(:) = fcai_lay(:)*faii_lay(:)

   !=============================================
   ! Calculate layer direct beam radiation absorbed
   ! in the sunlit canopy as direct
   !=============================================

      fadd_lay(:,ib) = D0

      DO lev = 1, nlay
         IF( fc0(lev)>D0 .and. lsai_lay(lev)>D0 ) THEN

            fadd_lay(lev,ib) = tt(lev+1,lev) * & 
               (D1-ftdd_lay(lev)) * (D1-omg_lay(lev,ib))

         ENDIF
      ENDDO

      A = D0; B = D0;
      fabs_leq = D0

      ! Calculate the coefficients matrix A
      A(1,1) = 1.0; A(1,3) = -shadow_i(3)*ftii_lay(3) + shadow_i(3) - 1.0;
      A(2,2) = 1.0; A(2,3) = -shadow_i(3)*frii_lay(3);
      A(3,3) = 1.0; A(3,2) = -shadow_i(2)*frii_lay(2);   A(3,5) = -shadow_i(2)*ftii_lay(2) + shadow_i(2) - 1.0;
      A(4,4) = 1.0; A(4,5) = -shadow_i(2)*frii_lay(2);   A(4,2) = -shadow_i(2)*ftii_lay(2) + shadow_i(2) - 1.0;
      A(5,5) = 1.0; A(5,4) = -shadow_i(1)*frii_lay(1);   A(5,6) =(-shadow_i(1)*ftii_lay(1) + shadow_i(1) - 1.0) * albgri(ib);
      A(6,6) = 1.0 - albgri(ib)*shadow_i(1)*frii_lay(1); A(6,4) = -shadow_i(1)*ftii_lay(1) + shadow_i(1) - 1.0;

      ! The constant vector B at right side
      B(1,1) = tt(4,3)*frid_lay(3); B(1,2) = shadow_i(3)*frii_lay(3);
      B(2,1) = tt(4,3)*ftid_lay(3); B(2,2) = shadow_i(3)*ftii_lay(3) - shadow_i(3) + 1.0;
      B(3,1) = tt(3,2)*frid_lay(2); B(3,2) = 0.0;
      B(4,1) = tt(3,2)*ftid_lay(2); B(4,2) = 0.0;
      B(5,1) = tt(2,1)*frid_lay(1) + tt(1,0)*albgrd(ib)*(shadow_i(1)*ftii_lay(1) - shadow_i(1) + 1.0); B(5,2) = 0.0;
      B(6,1) = tt(2,1)*ftid_lay(1) + tt(1,0)*albgrd(ib)*shadow_i(1)*frii_lay(1);                       B(6,2) = 0.0;

      ! Get the resolution
      CALL mGauss(A, B, X)

      ! ====================================================
      ! Set back to the absorption for each layer and albedo
      ! ====================================================

      ! Albedo
      fabs_leq(4,:) = X(1,:)

      ! Three layers' absorption
      fabs_leq(3,1) = tt(4,3)*faid_lay(3) + X(3,1)*shadow_i(3)*faii_lay(3)
      fabs_leq(3,2) = shadow_i(3)*faii_lay(3) + X(3,2)*shadow_i(3)*faii_lay(3)
      fabs_leq(2,1) = tt(3,2)*faid_lay(2) + (X(2,1) + X(5,1))*shadow_i(2)*faii_lay(2)
      fabs_leq(2,2) = (X(2,2) + X(5,2)) * shadow_i(2) * faii_lay(2)
      fabs_leq(1,1) = tt(2,1)*faid_lay(1) + (X(4,1) + X(6,1)*albgri(ib) + tt(1,0)*albgrd(ib))*shadow_i(1)*faii_lay(1)
      fabs_leq(1,2) = (X(4,2) + X(6,2)*albgri(ib)) * shadow_i(1) * faii_lay(1)

      ! Ground absorption
      fabs_leq(0,1) = X(6,1) * (1.0 - albgri(ib)) + tt(1,0) * (1.0 - albgrd(ib))
      fabs_leq(0,2) = X(6,2) * (1.0 - albgri(ib))

   !=============================================================
   ! main process for 3-layer  sunlight radiation
   !=============================================================
   ! available direct and diffuse energy continue to be absorbed
   ! reflect and transmit up and down direction till most of 
   ! enery is absorbed by ground, sky and three canopy layers
   !=============================================================
   ! layers numbering: 0:sky, 1,2,3=canopy layers 4:ground
   ! downward direction: from = 0,1,2,3  lay=1,2,3,4 to =2,3,4,5 
   ! upward direction:   from = 4,3,2    lay=3,2,1 to =2,1,0
   !=============================================================

      ! --- START delete ---
!      DO kband=1,2 !1:direct band  2:diffuse band
!
!         delta = -1 
!         klay  =  3
!         tri   = D0
!         trd   = D0
!         ref   = D0
!         fabs_lay(:,kband)=D0
!
!         !available radiation
!         IF (kband == 1) THEN     !direct sunlight case
!            trd(:,:)=tt
!            tri(:,:)=D0
!            ref(:,:)=D0
!         ELSE                     !diffuse sunlight case
!            tri(4,3)=D1
!            trd(:,:)=D0
!            ref(:,:)=D0
!         ENDIF
!
!         ! start layer integration
!         DO nn=1,40
!
!            ! terminate loop iteration IF total available energy is small
!            IF ( sum(trd+tri+ref) < 0.0001) EXIT
!            kfr = klay - delta
!            kto = klay + delta
!
!            ! direct and diffuse radiation at layer
!            dir = trd(kfr,klay)
!            dif = tri(kfr,klay) + ref(kfr,klay)
!
!            ! ground layer 
!            IF( kto == -1 ) THEN 
!               ref(klay,kfr) = (dir+dif)*albgrd(ib)
!               fabs_lay(klay,kband)=fabs_lay(klay,kband)+ &
!                  (dir+dif)*(D1-albgrd(ib))
!            ELSE            
!
!               tri(klay,kto)= dir*ftid_lay(klay) + &
!                  dif*(D1-shadow_i(klay)) + &
!                  dif*shadow_i(klay)*ftii_lay(klay)
!
!               ref(klay,kfr)= dir*frid_lay(klay) + &
!                  dif*shadow_i(klay)*frii_lay(klay)
!
!               fabs_lay(klay,kband)=fabs_lay(klay,kband) + &
!                  dir*faid_lay(klay) + &
!                  dif*shadow_i(klay)*faii_lay(klay)
!
!               ! reflection to sky --> albedo
!               IF( kfr == 4) THEN
!                  fabs_lay(kfr,kband) = fabs_lay(kfr,kband) + ref(klay,kfr)
!                  ref(klay,kfr)=D0
!               ENDIF
!
!               ! transmit to sky -> albedo
!               IF (kto == 4 ) THEN
!                  fabs_lay(kto,kband) = fabs_lay(kto,kband) + tri(klay,kto)
!                  tri(klay,kto) = 0.
!               ENDIF
!            ENDIF
!
!            ! clear the radiation source
!            trd(kfr,klay) = D0
!            tri(kfr,klay) = D0
!            ref(kfr,klay) = D0
!
!            ! change radiation direction: delta=1 down, 
!            ! delta=-1 upward radiation
!            IF (klay == 0) delta = 1
!            IF (klay == 3) delta =-1
!
!            ! Change lay
!            klay = klay + delta
!         ENDDO ! ENDDO nn loop
!
!         ! Need to assign the remain radiation
!         ! Simple method: totally absorbed by to
!         DO kfr=0,4
!            DO kto=0,4
!               IF (trd(kfr,kto) > D0)fabs_lay(kfr,kband)=fabs_lay(kfr,kband)+&
!                  trd(kfr,kto)
!               IF (tri(kfr,kto) > D0)fabs_lay(kfr,kband)=fabs_lay(kfr,kband)+&
!                  tri(kfr,kto)
!               IF (ref(kfr,kto) > D0)fabs_lay(kfr,kband)=fabs_lay(kfr,kband)+&
!                  ref(kfr,kto)
!            ENDDO
!         ENDDO
!
!         IF (maxval(abs(fabs_lay(:,kband)-fabs_leq(:,kband))) > 1.e-4) THEN
!            print *, maxval(abs(fabs_lay(:,kband)-fabs_leq(:,kband)))
!            print *, "linear equation solution error!"
!         ENDIF
!         !print *, maxval(abs(fabs_lay(:,kband)-fabs_leq(:,kband)))
!         
!         fabs_lay = fabs_leq
!
!         ! set column absorption and reflection
!         IF( kband==1) THEN
!            fabd_lay(1:3,ib) = fabs_lay(1:3,kband)
!            fabd_col(ib) = fabs_lay(1,kband)+fabs_lay(2,kband)+fabs_lay(3,kband)
!            albd_col(ib) = fabs_lay(4,kband)
!            IF (abs(fabd_col(ib)+albd_col(ib)+fabs_lay(0,kband)-1) > 1e-6) THEN
!               print *, "Imbalance kband=1"
!               print *, fabd_col(ib)+albd_col(ib)+fabs_lay(0,kband)-1
!            ENDIF
!         ELSE
!            fabi_lay(1:3,ib) = fabs_lay(1:3,kband)
!            fabi_col(ib) = fabs_lay(1,kband)+fabs_lay(2,kband)+fabs_lay(3,kband)
!            albi_col(ib) = fabs_lay(4,kband)
!            IF (abs(fabi_col(ib)+albi_col(ib)+fabs_lay(0,kband)-1) > 1e-6) THEN
!               print *, "Imbalance kband=2"
!               print *, fabi_col(ib)+albi_col(ib)+fabs_lay(0,kband)-1
!            ENDIF
!         ENDIF
!      ENDDO ! ENDDO kband (direct/diffuse beams)
!      ! --- END delete ---

      ! IF everything is ok, substitute fabs_lay for fabs_leq
      ! and delete the following line and the variables defined
      ! but not used anymore
      fabs_lay = fabs_leq

      ! set column absorption and reflection
      fabd_lay(1:3,ib) = fabs_lay(1:3,1)
      fabi_lay(1:3,ib) = fabs_lay(1:3,2)
      fabd_col(ib) = fabs_lay(1,1)+fabs_lay(2,1)+fabs_lay(3,1)
      fabi_col(ib) = fabs_lay(1,2)+fabs_lay(2,2)+fabs_lay(3,2)
      albd_col(ib) = fabs_lay(4,1)
      albi_col(ib) = fabs_lay(4,2)
      IF (abs(fabd_col(ib)+albd_col(ib)+fabs_lay(0,1)-1) > 1e-6) THEN
         print *, "Imbalance kband=1"
         print *, fabd_col(ib)+albd_col(ib)+fabs_lay(0,1)-1
      ENDIF
      IF (abs(fabi_col(ib)+albi_col(ib)+fabs_lay(0,2)-1) > 1e-6) THEN
         print *, "Imbalance kband=2"
         print *, fabi_col(ib)+albi_col(ib)+fabs_lay(0,2)-1
      ENDIF

   !====================================================
   ! Calculate individule PFT absorption
   !====================================================

      sum_fabd=D0
      sum_fabi=D0
      sum_fadd=D0

      DO ip = 0, npft-1
         clev = canlev(ip)
         IF(clev == D0) cycle
         IF(shadow_d(clev) > D0 .and. soilveg(ip)) THEN

         !=================================================
         ! fractional contribution of current pft in layer
         !=================================================

            pfc = min( pwtcol(ip)/fc0(clev), D1)

         !=========================================
         ! shadow contribution from ground to sky
         !=========================================

            !shadow_sky(ip) = pfc*shadow_i(clev)
            shadow_sky(ip) = shadow_pi(ip)

         !=======================================================
         ! absorption, reflection and transmittance fluxes for 
         ! unit incident radiation over pft.
         !=======================================================

            CALL CanopyRad(taud(ip), taui(ip), ftdd_orig(ip,ib), ftdi_orig(ip,ib), &
               cosz,cosd, shadow_pd(ip), shadow_pi(ip), pwtcol(ip),&
               omega(ip,ib), lsai(ip), tau(ip,ib),&
               rho(ip,ib), ftid(ip,ib), ftii(ip,ib), albd(ip,ib),&
               albi(ip,ib), faid_p, faii_p)

            ! calibration for LAD
            ! 11/07/2018: calibration for LAD
            ftid(ip,ib) = fcad(ip)*ftid(ip,ib)
            ftii(ip,ib) = fcai(ip)*(ftii(ip,ib)-ftdi_orig(ip,ib)) + ftdi(ip,ib)
            albd(ip,ib) = fcad(ip)*albd(ip,ib)
            albi(ip,ib) = fcai(ip)*albi(ip,ib)
            faid_p = fcad(ip)*faid_p
            faii_p = fcai(ip)*faii_p

            ! absorptions after multiple reflections for each pft
            probm = albi(ip,ib)*shadow_sky(ip)*albgri(ib)
            ftran = (D1-shadow_pd(ip)+shadow_pd(ip)*ftdd(ip,ib))*albgrd(ib) &
                  + shadow_pd(ip)*ftid(ip,ib)*albgri(ib)
            fabsm = ftran*faii_p*shadow_sky(ip)/(D1-probm)
            fabd(ip,ib) = shadow_pd(ip)*faid_p + fabsm

            probm = albi(ip,ib)*shadow_sky(ip)*albgri(ib)
            ftran = D1-shadow_pi(ip)*(D1 -ftii(ip,ib))
            fabsm = ftran*albgri(ib)*faii_p*shadow_sky(ip)/(D1-probm)
            fabi(ip,ib) = shadow_pi(ip)*faii_p + fabsm

            ! column albedo is assigned to each pft in column
            ! Deleted by Yuan, 06/03/2012 
            !albd(ip,ib) =albd_col(ib)
            !albi(ip,ib) =albi_col(ib)

            ! sum of pft absorptions in column
            sum_fabd(clev) = sum_fabd(clev) + fabd(ip,ib)
            sum_fabi(clev) = sum_fabi(clev) + fabi(ip,ib)

            ! pft absorption in sunlit as direct beam
            sun_fadd(ip,ib) = shadow_pd(ip) * (D1-ftdd(ip,ib)) * (D1-omega(ip,ib))

            ! sum of pft absorption in sunlit as direct beam
            sum_fadd(clev) = sum_fadd(clev) + sun_fadd(ip,ib)

         ENDIF ! ENDIF shadow & soilveg
      ENDDO ! ENDDO ip

      DO ip = 0, npft-1
         clev = canlev(ip)

         !IF(shadow_d(clev) > D0 .and. soilveg(ip)) THEN             
         ! Modified by Yuan, 06/03/2012
         !IF(soilveg(ip) .or. ivt(ip) == 0) THEN             

      !===========================================================
      ! adjust pft absorption for total column absorption per 
      ! unit column area
      !===========================================================

         IF (soilveg(ip)) THEN 
            fabd(ip,ib)=fabd(ip,ib)*fabd_lay(clev,ib)/&
               sum_fabd(clev)/pwtcol(ip)
            fabi(ip,ib)=fabi(ip,ib)*fabi_lay(clev,ib)/&
               sum_fabi(clev)/pwtcol(ip)

            sun_fadd(ip,ib) = sun_fadd(ip,ib)*fadd_lay(clev,ib)/&
               sum_fadd(clev)/pwtcol(ip)
            psun(ip) = tt(clev+1,clev)/shadow_d(clev)
         ELSE
            fabd(ip,ib) = D0
            fabi(ip,ib) = D0
            sun_fadd(ip,ib) = D0
            psun(ip) = D0
         ENDIF

         !fabd(ip,ib) = fabd_col(ib)
         !fabi(ip,ib) = fabi_col(ib)

         ! column albedo is assigned to each pft in column
         ! Added by Yuan, 06/03/2012 
         albd(ip,ib) =albd_col(ib)
         albi(ip,ib) =albi_col(ib)

         ! adjust ftdd and ftii for multi reflections between layers

! 03/06/2020, yuan: NOTE! there is no physical mean of ftdd,
! ftid, ftii anymore. they are the same for each PFT can only
! be used to calculate the ground absorption. 
         ftdd(ip,ib) = ftdd_col
         ftid(ip,ib)=(D1-albd(ip,ib)-fabd_col(ib)-&
            ftdd(ip,ib)*(D1-albgrd(ib)))/(D1-albgri(ib))

         !ftdd(ip,ib) = (D1 - albd(ip,ib) - fabd(ip,ib) - &
         !     ftid(ip,ib)*(D1-albgri(ib)))/(D1-albgrd(ib))
         ftii(ip,ib)=(D1-albi(ip,ib)-fabi_col(ib))/(D1-albgri(ib))

         !ftdd(ip,ib) = min(max(ftdd(ip,ib),D0),D1)
         !ftii(ip,ib) = min(max(ftii(ip,ib),D0),D1)
         !ftid(ip,ib) = min(max(ftid(ip,ib),D0),D1)

         ! check energy balance
         !fabd(ip,ib) = D1 - albd(ip,ib) - &
         !    ftdd(ip,ib)*(D1-albgrd(ib)) - &
         !    ftid(ip,ib)*(D1-albgri(ib))
         !abi(ip,ib) = D1 - albi(ip,ib) - &
         !    ftii(ip,ib)*(D1-albgri(ib))

      ENDDO ! ENDDO ip
   ENDDO !ENDDO ib

   ! set parameters for longwave calculation
   fshade(:) = shadow_pi(:)
   thermk(:) = ftdi(:,1)

END SUBROUTINE ThreeDCanopy

!=====================  
! FUNCTION tee
!=====================  

REAL(selected_real_kind(12)) FUNCTION tee(tau)

   IMPLICIT NONE
   INTEGER,parameter::r8 = selected_real_kind(12)  ! 64 bit REAL
   INTEGER,parameter::r16 = selected_real_kind(24) ! 128 bit REAL

   REAL(r16),parameter::DDH=0.50_r16 !128-bit accuracy REAL
   REAL(r16),parameter::DD1=1.0_r16  !128-bit accuracy REAL
   REAL(r16),parameter::DD2=2.0_r16  !128-bit accuracy REAL
   REAL(r16)::tau ! transmittance

   tee = DDH*(DD1/tau/tau-(DD1/tau/tau+DD2/tau)*exp(-DD2*tau))

END FUNCTION tee

!===========================================
! FUNCTION overlapArea
!===========================================

REAL(selected_real_kind(12)) FUNCTION OverlapArea(radius, hgt, zenith)

   IMPLICIT NONE
   INTEGER,parameter :: r8  = selected_real_kind(12) 
   INTEGER,parameter :: r16 = selected_real_kind(24) 

   REAL(r8),parameter :: rpi = 3.14159265358979323846_R8  ! pi
   REAL(r8),parameter::D0=0.0_r8  !128-bit accuracy REAL
   REAL(r8),parameter::D1=1.0_r8  !128-bit accuracy REAL

   REAL(r8)::radius !radius of bus
   REAL(r8)::hgt    !height of canopy
   REAL(r8)::zenith !zenith angle
   REAL(r8)::cost   !cosine of angle
   REAL(r8)::theta  !angle

   IF( radius == D0) THEN
      OverlapArea= D0
      RETURN
   ENDIF
   cost = hgt*tan(zenith)/radius/(D1+D1/cos(zenith))
   IF ( cost >= 1) THEN
      OverlapArea= D0
      RETURN
   ENDIF
   theta = acos(cost)
   OverlapArea = (theta-cost*sin(theta))*(D1+D1/cos(zenith))/rpi
   RETURN
END FUNCTION OverlapArea

!=========================================================
! FUNCTION to calculate scattering, absorption, reflection and 
! transmittance for unit input radiation
!=========================================================

SUBROUTINE CanopyRad(tau_d, tau_i, ftdd, ftdi, cosz,cosd, &
      shadow_d, shadow_i, fc, omg, lsai, &
      tau_p,  rho_p, ftid, ftii, frid, frii, faid, faii)
   IMPLICIT NONE

   INTEGER,parameter::r8 = selected_real_kind(12) ! 8 byte REAL
   INTEGER,parameter::r16 = selected_real_kind(24) ! 8 byte REAL

   ! input variables
   REAL(r8)::cosz      !0.001 <= coszen <= 1.000
   REAL(r8)::cosd      !0.001 <= coszen <= 1.000
   REAL(r8)::faid      !direct absorption
   REAL(r8)::faii      !diffuse absorption
   REAL(r8)::fc        !fraction of grid covered with canopy
   REAL(r8)::frid      !direct reflectance
   REAL(r8)::frii      !diffuse reflectance
   REAL(r8)::frio      !diffuse reflectance
   REAL(r8)::ftdd      !down direct flux below veg per unit dir flx
   REAL(r8)::ftdi      !down direct flux below veg per unit dif flux
   REAL(r8)::ftid      !direct transmittance 
   REAL(r8)::ftii      !diffuse transmittance
   REAL(r8)::omg       !frac of intercepted rad that is scattered
   REAL(r8)::rho_p     !leaf/stem reflectance weighted by fract of LAI and SAI
   REAL(r8)::shadow_d  !canopy shadow for direct solar
   REAL(r8)::shadow_i  !canopy shadow for diffuse solar
   REAL(r8)::tau_d     !radial optical depth for direct beam
   REAL(r8)::tau_i     !radial optical depth for indirect beam
   REAL(r8)::tau_p     !leaf/stem transmission weighted by frac of LAI & SAI
   REAL(r8)::lsai      !elai+esai

   ! output variables
   REAL(r8)::phi_dif_d !differnce of rad scattered forward-backward per direct beam
   REAL(r8)::phi_dif_i !difference of rad scattered forward-backward per direct beam
   REAL(r8)::phi_tot_d !total rad scattered in all direction per direct beam
   REAL(r8)::phi_tot_i !total rad scattered in all direction per diffuse beam
   REAL(r8)::phi_tot_o !total rad scattered in all direction per direct beam
   REAL(r8)::phi_dif_o !total rad scattered in all direction per diffuse beam
   REAL(r8)::pa2       !total rad scattered in all direction per direct beam

   ! local variables
   LOGICAL::runmode = .true.
   REAL(r8)::tau
   REAL(r8)::muv       !forward frac of 3D scat rad in all direction for diffuse
   REAL(r8)::ac        !forward frac of 3D scat rad in all direction for diffuse
   REAL(r8)::ald       !forward frac of 3D scat rad in all direction for diffuse
   REAL(r8)::ali       !forward frac of 3D scat rad in all direction for diffuse

   REAL(r8)::wb
   REAL(r8)::alpha
   REAL(r8)::nd
   REAL(r8)::ni
   REAL(r8)::gee=0.5_r8                !Ross factor geometric blocking

   REAL(r8),parameter::D0 = 0.0_r8     !64-bit REAL number
   REAL(r8),parameter::D1 = 1.0_r8     !64-bit REAL number
   REAL(r8),parameter::D2 = 2.0_r8     !64-bit REAL number
   REAL(r8),parameter::D3 = 3.0_r8     !64-bit REAL number
   REAL(r8),parameter::D4 = 4.0_r8     !64-bit REAL number
   REAL(r8),parameter::D6 = 6.0_r8     !64-bit REAL number
   REAL(r8),parameter::DH = 0.5_r8     !64-bit REAL number
   REAL(r16),parameter::DD1 = 1.0_r16  !128-bit REAL number

   REAL(r8),parameter :: pi = 3.14159265358979323846_R8  ! pi
   REAL(selected_real_kind(12)), external::tee

   tau = D3/D4*gee*lsai

   CALL phi(runmode, tau_d, omg, tau_p, rho_p, phi_tot_d, phi_dif_d, pa2)
   CALL phi(runmode, tau_i, omg, tau_p, rho_p, phi_tot_i, phi_dif_i, pa2)
   CALL phi(runmode, tau  , omg, tau_p, rho_p, phi_tot_o, phi_dif_o, pa2)

   IF (runmode) THEN
      ! NOTE: modified
      frio = DH*(phi_tot_o - DH*phi_dif_o)
      frio = max(min(frio,D1),D0)

      muv = D3*( D1 - sqrt(D1-sqrt(D3)*fc/(D2*pi)) ) + &
         D3*( D1 - sqrt(D1-sqrt(D3)*fc/(D6*pi)) )

      wb = D2/D3*rho_p + D1/D3*tau_p
      alpha = sqrt(D1-omg) * sqrt(D1-omg+D2*wb)
      nd = (D1 + D2*alpha) / (D1 + D2*alpha*cosz)
      ni = (D1 + D2*alpha) / (D1 + D2*alpha*cosd)

      ac  = phi_tot_o * muv * (D1-tee(DD1*tau)) * (D1-omg) / (D1-omg*pa2)
      ald = (nd-D1) * frio * fc * (D1/shadow_d - cosz/fc)
      ali = (ni-D1) * frio * fc * (D1/shadow_i - cosd/fc)
   ENDIF

!-----------------------------------------------------------------------
!frac indirect downward rad through canopy for black soil & direct solar
!-----------------------------------------------------------------------
   frid = DH*(phi_tot_d - DH*cosz*phi_dif_d)
   frii = DH*(phi_tot_i - DH*cosd*phi_dif_i)

   IF (runmode) THEN
      frid = frid + ald - DH*ac
      frii = frii + ali - DH*ac
   ENDIF

   frid = max(min(frid,D1),D0)
   frii = max(min(frii,D1),D0)

!---------------------------------------------------------------------
!downward diffuse fraction from direct and diffuse sun
!---------------------------------------------------------------------
   ftid = DH*(phi_tot_d + DH*cosz*phi_dif_d)
   ftii = DH*(phi_tot_i + DH*cosd*phi_dif_i)+ftdi 

   IF (runmode) THEN
      ftid = ftid - DH*ald - DH*ac
      ftii = ftii - DH*ali - DH*ac
   ENDIF

   ftid = max(min(ftid,D1),D0)
   ftii = max(min(ftii,D1),D0)

!---------------------------------------------------------------------
! canopy absorption for direct or diffuse beams
!---------------------------------------------------------------------
   IF (.not. runmode) THEN
      faid =  D1 - ftdd - phi_tot_d
      faii =  D1 - ftdi - phi_tot_i
   ELSE 
      faid =  D1 - ftdd - frid - ftid
      faii =  D1 - frii - ftii
   ENDIF

   faid = max(min(faid,D1),D0)
   faii = max(min(faii,D1),D0)

   IF (shadow_d == D0) THEN
      ! NOTE: corrected from D1 -> D0
      ftid = D0
      frid = D0
      faid = D0
   ENDIF
   IF ( shadow_i == D0) THEN
      ftii = D1
      frii = D0
      faii = D0
   ENDIF

END SUBROUTINE CanopyRad

SUBROUTINE phi(runmode, tau, omg, tau_p, rho_p, phi_tot, phi_dif, pa2)

   IMPLICIT NONE

   INTEGER,parameter::r8  = selected_real_kind(12) ! 8 byte REAL
   INTEGER,parameter::r16 = selected_real_kind(24) ! 8 byte REAL

   ! input variables
   LOGICAL::runmode
   REAL(r8)::omg       !frac of intercepted rad that is scattered
   REAL(r8)::rho_p     !leaf/stem reflectance weighted by fract of LAI and SAI
   REAL(r8)::tau       !radial optical depth for direct beam
   REAL(r8)::tau_p     !leaf/stem transmission weighted by frac of LAI & SAI

   ! output variables
   REAL(r8)::phi_dif   !differnce of rad scattered forward-backward 
   REAL(r8)::phi_tot   !total rad scattered in all direction
   REAL(r8)::pa2       !total rad scattered in all direction

   ! local variables
   REAL(r8)::pac       !probablity of absorption after two scatterings
   REAL(r8)::phi_1b    !backward single scattered radiation 
   REAL(r8)::phi_1f    !forward single scattered radiation 
   REAL(r8)::phi_2a    !average second-order scattered radiation
   REAL(r8)::phi_2b    !backward second-order scattered radiation
   REAL(r8)::phi_2f    !forward second-order scattered radiation
   REAL(r8)::phi_mb    !backward multiple scattered radiation
   REAL(r8)::phi_mf    !forward multiple scattered radiation
   REAL(r8)::phi_tb    !backward frac of 3D scat rad in all direction
   REAL(r8)::phi_tf    !forward frac of 3D scat rad in all direction
   REAL(r8)::aa,bb     !temporary constants

   REAL(r8),parameter::D0 = 0.0_r8     !64-bit REAL number
   REAL(r8),parameter::D1 = 1.0_r8     !64-bit REAL number

   REAL(r16),parameter::DD1 = 1.0_r16  !128-bit REAL number
   REAL(r16),parameter::DD2 = 2.0_r16  !128-bit REAL number
   REAL(r16),parameter::DD3 = 3.0_r16  !128-bit REAL number
   REAL(r16),parameter::DD4 = 4.0_r16  !128-bit REAL number
   REAL(r16),parameter::DD9 = 9.0_r16  !128-bit REAL number
   REAL(r16),parameter::DD10= 10.0_r16 !128-bit REAL number
   REAL(r16),parameter::DDH = 0.5_r16  !128-bit REAL number

   REAL(selected_real_kind(12)), external::tee

!---------------------------------------------------------------------
! single scattering terms for sphere with overlap corrections to path
! for direct and diffuse beams
!---------------------------------------------------------------------

   ! forward first order normalized scattering
   phi_1f =(DD1/tau/tau - (DD1/tau/tau + DD2/tau + DD2)*exp(-DD2*tau))

   ! backward first order normalized scattering
   phi_1b = DDH*(DD1 - tee(DD2*tau))

!---------------------------------------------------------------------
! sphere double scattering terms (RED 2008 Eq 19,20)
!---------------------------------------------------------------------

   IF (.not. runmode) THEN

      ! forward double scattering 
      phi_2f = DDH*(DD4*phi_1f/DD3 + tee(DD2*tau) + tee(DD4*tau)/DD9 &
         -DD10*tee(DD1*tau)/DD9)

      ! backward double scattering 
      phi_2b = DDH*(DD1/DD3 - tee(DD2*tau) + DD2*tee(DD3*tau)/DD3)

   ELSE 
      ! fitting FUNCTION for second order scattering
      aa = 0.70_r8
      bb = 1.74_r8

      phi_2b = aa*( DD1/(bb+DD1) -DD1/(bb-D1)*tee(DD2*tau) + &
         DD2/(bb+DD1)/(bb-DD1)*tee((DD1+bb)*tau) )

      phi_2f = aa*( DD2*bb/(bb*bb-DD1)*phi_1f - & 
         (DD1/(bb+DD1)/(bb+DD1) + DD1/(bb-DD1)/(bb-DD1))*tee(DD1*tau) + &
         DD1/(bb-DD1)/(bb-DD1)*tee(DD1*tau*bb) + &
         DD1/(bb+DD1)/(bb+DD1)*tee(DD1*(bb+DD2)*tau) )
   ENDIF

   ! second order avaerage scattering
   phi_2a = DDH*(phi_2b + phi_2f)

!---------------------------------------------------------------------
! probabilty of absorption after two scattering
!---------------------------------------------------------------------

   ! probabilty of absorption for diffuse beam
   ! corrected probabilty of absorption for direct beam
   pac = DD1-phi_2a /(DD1 - tee(DD1*tau) - (rho_p*phi_1b + &
      tau_p*phi_1f)/(tau_p+rho_p)) 

   ! NOTE: for test only
   pac = max(min(pac,D1),D0)
   pa2 = pac

!---------------------------------------------------
!third order and higher order scatterings
!---------------------------------------------------

   phi_mf = phi_2f + omg*pac*phi_2a/(DD1-omg*pac)
   phi_mb = phi_2b + omg*pac*phi_2a/(DD1-omg*pac)

!----------------------------------------------------------------------
! total sphere scattering,forward,backward, avg & diff for direct beam
!----------------------------------------------------------------------

   phi_tf  = tau_p*phi_1f + DDH*omg*omg*phi_mf
   phi_tb  = rho_p*phi_1b + DDH*omg*omg*phi_mb

   phi_tot = phi_tf + phi_tb
   phi_dif = phi_tf - phi_tb

END SUBROUTINE phi

SUBROUTINE mGauss(A, B, X)

   IMPLICIT NONE

   INTEGER,parameter :: r8  = selected_real_kind(12) 

   REAL(r8), intent(inout)  :: A(6,6)
   REAL(r8), intent(inout)  :: B(6,2)
   REAL(r8), intent(out)    :: X(6,2)

   INTEGER :: i, j
   INTEGER :: nstep(5) = (/0, 2, 1, 2, 1/)

   REAL(r8) :: f

   ! Elimination
   DO i = 1, 5
      DO j = i+1, i+nstep(i)
         IF (abs(A(i,i)) < 1.E-10) THEN
            print *, "Error in Gauss's solution"
            RETURN
         ENDIF
         f = - A(j,i)/A(i,i)
         A(j,:) = A(j,:) + f*A(i,:)
         B(j,:) = B(j,:) + f*B(i,:)
      ENDDO
   ENDDO

   ! Back substitution 
   X(6,:) = B(6,:)/A(6,6)
   DO i = 5, 1, -1
      X(i,1) = (B(i,1) - sum(A(i,i+1:6)*X(i+1:6,1))) / A(i,i)
      X(i,2) = (B(i,2) - sum(A(i,i+1:6)*X(i+1:6,2))) / A(i,i)
   ENDDO

END SUBROUTINE mGauss
