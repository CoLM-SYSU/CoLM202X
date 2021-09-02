#include <define.h>

MODULE ALBEDO

!-----------------------------------------------------------------------
 USE precision
 IMPLICIT NONE
 SAVE

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: albland
  PUBLIC :: albocean


! PRIVATE MEMBER FUNCTIONS:
  PRIVATE :: twostream
  PRIVATE :: twostream_mod
  PRIVATE :: twostream_wrap


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------



  SUBROUTINE albland (ipatch, patchtype,&
                      soil_s_v_alb,soil_d_v_alb,soil_s_n_alb,soil_d_n_alb,&
                      chil,rho,tau,fveg,green,lai,sai,coszen,&
                      wt,fsno,scv,sag,ssw,tg,&
                      alb,ssun,ssha,thermk,extkb,extkd) 

!=======================================================================
! Calculates fragmented albedos (direct and diffuse) in
! wavelength regions split at 0.7um.
! 
! (1) soil albedos: as in BATS formulations, which are the function of
!     soil color and moisture in the surface soil layer
! (2) snow albedos: as in BATS formulations, which are inferred from
!     the calculations of Wiscombe and Warren (1980) and the snow model
!     and data of Anderson(1976), and the function of snow age, grain size,
!     solar zenith angle, pollution, the amount of the fresh snow
! (3) canopy albedo: two-stream approximation model 
! (4) glacier albedos: as in BATS, which are set to constants (0.8 for visible beam,
!     0.55 for near-infrared)
! (5) lake and wetland albedos: as in BATS, which depend on cosine solar zenith angle,
!     based on data in Henderson-Sellers (1986). The frozen lake and wetland albedos
!     are set to constants (0.6 for visible beam, 0.4 for near-infrared)
! (6) over the snow covered tile, the surface albedo is estimated by a linear
!     combination of albedos for snow, canopy and bare soil (or lake, wetland, glacier).
!
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002, 03/2014
!=======================================================================

  USE precision
  USE PhysicalConstants, only: tfrz
  USE MOD_PFTimeInvars
  USE MOD_PFTimeVars
  USE MOD_PCTimeInvars
  USE MOD_PCTimeVars
  
  IMPLICIT NONE

!------------------------- Dummy Arguments -----------------------------
! ground cover index
 INTEGER, intent(in) :: &
      ipatch,   & ! patch index
      patchtype   ! land water type (0=soil, 1=urban or built-up, 2=wetland,
                  ! 3=land ice, 4=deep lake, 5=shallow lake)

                  ! parameters
 REAL(r8), intent(in) :: & 
      soil_s_v_alb, &! albedo of visible of the saturated soil
      soil_d_v_alb, &! albedo of visible of the dry soil
      soil_s_n_alb, &! albedo of near infrared of the saturated soil
      soil_d_n_alb, &! albedo of near infrared of the dry soil
      chil,      &! leaf angle distribution factor
      rho(2,2),  &! leaf reflectance (iw=iband, il=life and dead)
      tau(2,2),  &! leaf transmittance (iw=iband, il=life and dead)
      fveg,      &! fractional vegetation cover [-]
      green,     &! green leaf fraction
      lai,       &! leaf area index (LAI+SAI) [m2/m2]
      sai,       &! stem area index (LAI+SAI) [m2/m2]

                  ! variables
      coszen,    &! cosine of solar zenith angle [-]
      wt,        &! fraction of vegetation covered by snow [-]
      fsno,      &! fraction of soil covered by snow [-]
      ssw,       &! water volumetric content of soil surface layer [m3/m3]
      scv,       &! snow cover, water equivalent [mm]
      sag,       &! non dimensional snow age [-]
      tg          ! ground surface temperature [K]

 REAL(r8), intent(out) :: &
      alb(2,2),  &! averaged albedo [-]
      ssun(2,2), &! sunlit canopy absorption for solar radiation
      ssha(2,2), &! shaded canopy absorption for solar radiation,
                  ! normalized by the incident flux
      thermk,    &! canopy gap fraction for tir radiation
      extkb,     &! (k, g(mu)/mu) direct solar extinction coefficient
      extkd       ! diffuse and scattered diffuse PAR extinction coefficient

!-------------------------- Local variables ----------------------------
 INTEGER         &!
      iw,        &! wavelength (1=visible, 2=near-infrared)
      id,        &! 1=direct, 2=diffuse
      k           ! looping indx

 REAL(r8)        &!
      age,       &! factor to reduce visible snow alb due to snow age [-]
      albg0,     &! temporary varaiable [-]
      albsno(2,2),&! snow albedo [-]
      albg(2,2), &! albedo, ground
      albv(2,2), &! albedo, vegetation [-]
      alb_s_inc, &! decrease in soil albedo due to wetness [-]
      beta0,     &! upscattering parameter for direct beam [-]
      cff,       &! snow alb correction factor for zenith angle > 60 [-]
      conn,      &! constant (=0.5) for visible snow alb calculation [-]
      cons,      &! constant (=0.2) for nir snow albedo calculation [-]
      czen,      &! cosine of solar zenith angle > 0 [-]
      czf,       &! solar zenith correction for new snow albedo [-]
      dfalbl,    &! snow albedo for diffuse nir radiation [-]
      dfalbs,    &! snow albedo for diffuse visible solar radiation [-]
      dralbl,    &! snow albedo for visible radiation [-]
      dralbs,    &! snow albedo for near infrared radiation [-]
      fsol1,     &! solar flux fraction for wavelength < 0.7 micron [-]
      fsol2,     &! solar flux fraction for wavelength > 0.7 micron [-]
      lsai,      &! leaf and stem area index (LAI+SAI) [m2/m2]
      scat(2),   &! single scattering albedo for vir/nir beam [-]
      sl,        &! factor that helps control alb zenith dependence [-]
      snal0,     &! alb for visible,incident on new snow (zen ang<60) [-]
      snal1,     &! alb for NIR, incident on new snow (zen angle<60) [-]
      upscat,    &! upward scattered fraction for direct beam [-]
      tran(2,2)   ! canopy transmittances for solar radiation

   INTEGER ps, pe, pc
! ----------------------------------------------------------------------
! 1. Initial set
! ----------------------------------------------------------------------
! division of solar flux for wavelength less or greater than 0.7 micron
      fsol1 = 0.5      ! shortwave
      fsol2 = 0.5      ! longwave

! short and long wave albedo for new snow
      snal0 = 0.85     ! shortwave
      snal1 = 0.65     ! long wave

! set initial leaf scattering reflectance. Note: "scat" may use different
! value for different vegetation latter
      beta0 = 0.5
      scat(1) = 0.15
      scat(2) = 0.85

! ----------------------------------------------------------------------
! set default soil and vegetation albedos and solar absorption
      alb (:,:) = 0. ! averaged
      albg(:,:) = 0. ! ground
      albv(:,:) = 0. ! vegetation
      ssun(:,:) = 0.
      ssha(:,:) = 0.
      thermk    = 1.e-3
      extkb     = 1.
      extkd     = 0.718

! 08/25/2019, yuan:
IF (patchtype == 0) THEN
#ifdef PFT_CLASSIFICATION
      ps = patch_pft_s(ipatch)      
      pe = patch_pft_e(ipatch)
      ssun_p(:,:,ps:pe) = 0.
      ssha_p(:,:,ps:pe) = 0.
      thermk_p(ps:pe)   = 1.e-3
      extkb_p(ps:pe)    = 1.
      extkd_p(ps:pe)    = 0.718
#endif

#ifdef PC_CLASSIFICATION
      pc = patch2pc(ipatch)
      ssun_c(:,:,:,pc) = 0.
      ssha_c(:,:,:,pc) = 0.
      thermk_c(:,pc)   = 1.e-3
      extkb_c(:,pc)    = 1.
      extkd_c(:,pc)    = 0.718
#endif
ENDIF
      
      lsai=lai+sai
      IF(coszen<=0.) THEN
         RETURN  !only do albedo when coszen > 0
      ENDIF

      czen=max(coszen,0.001) 
      albsno(:,:)=0.         !set initial snow albedo

! ----------------------------------------------------------------------
! 2. albedo for snow cover.
!    snow albedo depends on snow-age, zenith angle, and thickness
!    of snow age gives reduction of visible radiation
! ----------------------------------------------------------------------
      IF(scv>0.)THEN
         cons = 0.2
         conn = 0.5
         sl  = 2.0           !sl helps control albedo zenith dependence

         ! correction for snow age
         age = 1.-1./(1.+sag) !correction for snow age
         dfalbs = snal0*(1.-cons*age)

         ! czf corrects albedo of new snow for solar zenith
         cff    = ((1.+1./sl)/(1.+czen*2.*sl )- 1./sl)
         cff    = max(cff,0.)
         czf    = 0.4*cff*(1.-dfalbs)
         dralbs = dfalbs+czf
         dfalbl = snal1*(1.-conn*age)
         czf    = 0.4*cff*(1.-dfalbl)
         dralbl = dfalbl+czf
   
         albsno(1,1) = dralbs
         albsno(2,1) = dralbl
         albsno(1,2) = dfalbs
         albsno(2,2) = dfalbl
      ENDIF

! ----------------------------------------------------------------------
! 3. get albedo over land
! ----------------------------------------------------------------------
! 3.1 bare soil albedos, depends on moisture
      IF(patchtype<=2)THEN    ! wetland, permanent ice and water
         alb_s_inc = max(0.11-0.40*ssw, 0.)
         albg(1,1) = min(soil_s_v_alb + alb_s_inc, soil_d_v_alb)
         albg(2,1) = min(soil_s_n_alb + alb_s_inc, soil_d_n_alb)
         albg(:,2) = albg(:,1)         !diffused albedos for bare soil

! 3.2 albedos for permanent ice sheet. 
      ELSE IF(patchtype==3) THEN       !permanent ice sheet
         albg(1,:) = 0.8
         albg(2,:) = 0.55

! 3.3 albedo for inland water (NOTE: wetland is removed)
      ELSE IF(patchtype>=4) THEN             
         albg0 = 0.05/(czen+0.15)
         albg(:,1) = albg0
         albg(:,2) = 0.1               !Subin (2012)

         IF(tg<tfrz)THEN               !frozen lake and wetland
            albg(1,:) = 0.6
            albg(2,:) = 0.4
         ENDIF
      ENDIF

! 3.4 correction due to snow cover
      albg(:,:) = (1.-fsno)*albg(:,:) + fsno*albsno(:,:)
      alb(:,:)  = albg(:,:)

! ----------------------------------------------------------------------
! 4. canopy albedos : two stream approximation  
! ----------------------------------------------------------------------
      IF (lai+sai > 1e-6) THEN

IF (patchtype == 0) THEN
#if(defined USGS_CLASSIFICATION || defined IGBP_CLASSIFICATION)
         CALL twostream (chil,rho,tau,green,lai,sai,&
                         czen,albg,albv,tran,thermk,extkb,extkd,ssun,ssha) 

         albv(:,:) = (1.-wt)*albv(:,:) + wt*albsno(:,:)
         alb(:,:)  = (1.-fveg)*albg(:,:) + fveg*albv(:,:)
#endif

#ifdef PFT_CLASSIFICATION
         CALL twostream_wrap (ipatch, czen, albg, &
                              albv, tran, thermk, extkb, extkd, ssun, ssha)
         alb(:,:) = albv(:,:)
#endif

#ifdef PC_CLASSIFICATION
         CALL ThreeDCanopy_wrap (ipatch, czen, albg, albv, ssun, ssha)
         alb(:,:) = albv(:,:)
#endif
      ELSE
         CALL twostream (chil,rho,tau,green,lai,sai,&
                         czen,albg,albv,tran,thermk,extkb,extkd,ssun,ssha) 

         albv(:,:) = (1.-wt)*albv(:,:) + wt*albsno(:,:)
         alb(:,:)  = (1.-fveg)*albg(:,:) + fveg*albv(:,:)
ENDIF
      ENDIF

!-----------------------------------------------------------------------

  END SUBROUTINE albland



  SUBROUTINE twostream ( chil, rho, tau, green, lai, sai, &
             coszen, albg, albv, tran, thermk, extkb, extkd, ssun, ssha )
                                                                              
!-----------------------------------------------------------------------
!                                                                               
!     calculation of canopy albedos via two stream approximation (direct
!     and diffuse ) and partition of incident solar
!                                                                               
! Original author: Yongjiu Dai, June 11, 2001
!
!-----------------------------------------------------------------------        
                                                                               
  USE precision
  IMPLICIT NONE

! parameters
  REAL(r8), intent(in) :: &
          ! static parameters associated with vegetation type
            chil,          &! leaf angle distribution factor
            rho(2,2),      &! leaf reflectance (iw=iband, il=life and dead)
            tau(2,2),      &! leaf transmittance (iw=iband, il=life and dead)

          ! time-space varying vegetation parameters
            green,         &! green leaf fraction
            lai,           &! leaf area index of exposed canopy (snow-free)
            sai             ! stem area index

! environmental variables
  REAL(r8), intent(in) :: &
            coszen,        &! consine of solar zenith angle
            albg(2,2)       ! albedos of ground

! output
  REAL(r8), intent(out) :: &
            albv(2,2),     &! albedo, vegetation [-]
            tran(2,2),     &! canopy transmittances for solar radiation
            thermk,        &! canopy gap fraction for tir radiation                    
            extkb,         &! (k, g(mu)/mu) direct solar extinction coefficient  
            extkd,         &! diffuse and scattered diffuse PAR extinction coefficient
            ssun(2,2),     &! sunlit canopy absorption for solar radiation
            ssha(2,2)       ! shaded canopy absorption for solar radiation,
                            ! normalized by the incident flux 
                                                                               
!-------------------------- local -----------------------------------           
  REAL(r8) :: &
           lsai,           &! lai+sai
           phi1,           &! (phi-1)
           phi2,           &! (phi-2)
           scat,           &! (omega)        
           proj,           &! (g(mu))        
           zmu,            &! (int(mu/g(mu)) 
           zmu2,           &! (zmu * zmu)    
           as,             &! (a-s(mu))      
           upscat,         &! (omega-beta)   
           betao,          &! (beta-0)       
           psi,            &! (h)            

           be,             &! (b)            
           ce,             &! (c)            
           de,             &! (d)            
           fe,             &! (f)            

           power1,         &! (h*lai)
           power2,         &! (k*lai)
           power3,         &!  

           sigma,          &! 
           s1,             &! 
           s2,             &! 
           p1,             &! 
           p2,             &! 
           p3,             &! 
           p4,             &! 
           f1,             &! 
           f2,             &!
           h1,             &!
           h4,             &!
           m1,             &!
           m2,             &!
           m3,             &!
           n1,             &!
           n2,             &!
           n3,             &!

           hh1,            &! (h1/sigma)     
           hh2,            &! (h2)           
           hh3,            &! (h3)           
           hh4,            &! (h4/sigma)     
           hh5,            &! (h5)           
           hh6,            &! (h6)           
           hh7,            &! (h7)           
           hh8,            &! (h8)           
           hh9,            &! (h9)           
           hh10,           &! (h10)          

           eup(2,2),       &! (integral of i_up*exp(-kx) )
           edown(2,2)       ! (integral of i_down*exp(-kx) )
  
  INTEGER iw                !
                                                                                
!-----------------------------------------------------------------------        
! projected area of phytoelements in direction of mu and 
! average inverse diffuse optical depth per unit leaf area

      phi1 = 0.5 - 0.633 * chil - 0.33 * chil * chil
      phi2 = 0.877 * ( 1. - 2. * phi1 )
         
      proj = phi1 + phi2 * coszen
      extkb = (phi1 + phi2 * coszen) / coszen

      extkd = 0.719

      IF (abs(phi1).gt.1.e-6 .and. abs(phi2).gt.1.e-6) THEN
         zmu = 1. / phi2 * ( 1. - phi1 / phi2 * log ( ( phi1 + phi2 ) / phi1 ) )
      ELSE IF (abs(phi1).le.1.e-6) THEN
         zmu = 1./0.877
      ELSE IF (abs(phi2).le.1.e-6) THEN
         zmu = 1./(2.*phi1)
      ENDIF
      zmu2 = zmu * zmu

      lsai   = lai + sai
      power3 = lsai / zmu
      power3 = min( 50., power3 )
      power3 = max( 1.e-5, power3 )
      thermk = exp(-power3)

      DO iw = 1, 2    ! WAVE_BAND_LOOP                                               

!-----------------------------------------------------------------------        
!     calculate average scattering coefficient, leaf projection and             
!     other coefficients for two-stream model.                                  
!-----------------------------------------------------------------------        

! + stem optical properties
      scat = lai/lsai * ( tau(iw,1) + rho(iw,1) ) + &
             sai/lsai * ( tau(iw,2) + rho(iw,2) ) 

      as = scat / 2. * proj / ( proj + coszen * phi2 )                               
      as = as * ( 1. - coszen * phi1 / ( proj + coszen * phi2 ) * &
               log ( ( proj + coszen * phi2 + coszen * phi1 ) / ( coszen * phi1 ) ) ) 

! + stem optical properties
      upscat = lai/lsai*tau(iw,1) + sai/lsai*tau(iw,2)                           
      upscat = 0.5 * ( scat + ( scat - 2. * upscat ) * &         
               (( 1. + chil ) / 2. ) ** 2 )                                     
      betao = ( 1. + zmu * extkb ) / ( scat * zmu * extkb ) * as            
                                                                                
!-----------------------------------------------------------------------        
!     intermediate variables identified in appendix of SE-85.                   
!-----------------------------------------------------------------------        
                                                                               
      be = 1. - scat + upscat                                                   
      ce = upscat                                                               
      de = scat * zmu * extkb * betao                                          
      fe = scat * zmu * extkb * ( 1. - betao )                                 

      psi = sqrt(be**2 - ce**2)/zmu                                            
      power1 = min( psi*lsai, 50. )                                            
      power2 = min( extkb*lsai, 50. )                                          
      s1 = exp( - power1 )                                                    
      s2 = exp( - power2 )                                                     

!-----------------------------------------------------------------------        
!     calculation of direct albedos and canopy transmittances.                  
!     albv(iw,1)     ( i-up ) 
!     tran(iw,irad)  ( i-down ) 
!-----------------------------------------------------------------------        

      p1 = be + zmu * psi
      p2 = be - zmu * psi
      p3 = be + zmu * extkb
      p4 = be - zmu * extkb

      f1 = 1. - albg(iw,2)*p1/ce
      f2 = 1. - albg(iw,2)*p2/ce
      
      h1 = - ( de * p4 + ce * fe )
      h4 = - ( fe * p3 + ce * de )

      sigma = ( zmu * extkb ) ** 2 + ( ce**2 - be**2 )                           

      IF (abs(sigma) .gt. 1.e-10) THEN     !<======

         hh1 = h1 / sigma
         hh4 = h4 / sigma
                                                                                
         m1 = f1 * s1
         m2 = f2 / s1
         m3 = ( albg(iw,1) - ( hh1 - albg(iw,2) * hh4 ) ) * s2
 
         n1 = p1 / ce
         n2 = p2 / ce
         n3 = - hh4

         hh2 = (m3*n2 - m2*n3) / (m1*n2 - m2*n1)
         hh3 = (m3*n1 - m1*n3) / (m2*n1 - m1*n2)

         hh5 = hh2 * p1 / ce
         hh6 = hh3 * p2 / ce

         albv(iw,1) = hh1 + hh2 + hh3                                   
         tran(iw,1) = hh4 * s2 + hh5 * s1 + hh6 / s1                    

         eup(iw,1) = hh1 * (1. - s2*s2) / (2.*extkb) &
                   + hh2 * (1. - s1*s2) / (extkb + psi) &
                   + hh3 * (1. - s2/s1) / (extkb - psi)

         edown(iw,1) = hh4 * (1. - s2*s2) / (2.*extkb) &
                     + hh5 * (1. - s1*s2) / (extkb + psi) &
                     + hh6 * (1. - s2/s1) / (extkb - psi)

      ELSE                               !<======

         m1 = f1 * s1
         m2 = f2 / s1
         m3 = h1 / zmu2 * ( lsai + 1. / (2.*extkb) ) * s2 &
            + albg(iw,2) / ce * ( - h1 / (2.*extkb) / zmu2 * &
              ( p3*lsai + p4 / (2.*extkb) ) - de ) * s2 &
            + albg(iw,1) * s2
 
         n1 = p1 / ce
         n2 = p2 / ce
         n3 = 1./ce * ( h1*p4 / (4.*extkb*extkb) / zmu2 + de)

         hh2 = (m3*n2 - m2*n3) / (m1*n2 - m2*n1)
         hh3 = (m3*n1 - m1*n3) / (m2*n1 - m1*n2)

         hh5 = hh2 * p1 / ce
         hh6 = hh3 * p2 / ce

         albv(iw,1) =  - h1 / (2.*extkb*zmu2) + hh2 + hh3
         tran(iw,1) = 1./ce * ( -h1 / (2.*extkb*zmu2) * &
                                ( p3*lsai + p4 / (2.*extkb) ) - de ) * s2 &
                     + hh5 * s1 + hh6 / s1

         eup(iw,1) = (hh2 - h1/(2.*extkb*zmu2)) * (1. - s2*s2) / (2.*extkb) &
                   + hh3 * (lsai - 0.) &
                   + h1/(2.*extkb*zmu2) * ( lsai*s2*s2 - (1. - s2*s2)/(2.*extkb) )

         edown(iw,1) = (hh5 - (h1*p4/(4.*extkb*extkb*zmu) + de)/ce) * &
                             (1. - s2*s2) / (2.*extkb) &
                     + hh6 * (lsai - 0.) &
                     + h1*p3/(ce*4.*extkb*extkb*zmu2) * &
                                         ( lsai*s2*s2 - (1. - s2*s2)/(2.*extkb) )
      
      ENDIF                              !<======

      ssun(iw,1) = (1.-scat) * ( 1.-s2 + 1. / zmu * (eup(iw,1) + edown(iw,1)) ) 
      ssha(iw,1) = scat * (1.-s2) &
               + ( albg(iw,2)*tran(iw,1) + albg(iw,1)*s2 - tran(iw,1) ) - albv(iw,1) &
               - ( 1. - scat ) / zmu * ( eup(iw,1) + edown(iw,1) ) 

!-----------------------------------------------------------------------        
!     calculation of diffuse albedos and canopy transmittances
!     albv(iw,2) ( i-up ) 
!     tran(iw,2) ( i-down ) 
!-----------------------------------------------------------------------        
                                                                               
      m1 = f1 * s1
      m2 = f2 / s1
      m3 = 0.

      n1 = p1 / ce
      n2 = p2 / ce
      n3 = 1.

      hh7 = -m2 / (m1*n2 - m2*n1)
      hh8 = -m1 / (m2*n1 - m1*n2)

      hh9 = hh7 * p1 / ce
      hh10 = hh8 * p2 / ce

      albv(iw,2) =  hh7 + hh8                                            
      tran(iw,2) = hh9 * s1 + hh10 / s1

      IF (abs(sigma) .gt. 1.e-10) THEN 
         eup(iw,2)   = hh7 * (1. - s1*s2) / (extkb + psi) &
                     + hh8 * (1. - s2/s1) / (extkb - psi)
         edown(iw,2) = hh9 * (1. - s1*s2) / (extkb + psi) &
                     + hh10 * (1. - s2/s1) / (extkb - psi)
      ELSE
         eup(iw,2)   = hh7 * (1. - s1*s2) / ( extkb + psi) + hh8 * (lsai - 0.)
         edown(iw,2) = hh9 * (1. - s1*s2) / ( extkb + psi) + hh10 * (lsai - 0.)
      ENDIF

      ssun(iw,2) = (1.-scat) / zmu * (eup(iw,2) + edown(iw,2))
      ssha(iw,2) = tran(iw,2) * ( albg(iw,2) -1. ) - ( albv(iw,2) - 1. ) &
                 - ( 1. - scat ) / zmu * ( eup(iw,2) + edown(iw,2) ) 

      ENDDO           ! WAVE_BAND_LOOP

  END SUBROUTINE twostream

  SUBROUTINE twostream_mod ( chil, rho, tau, green, lai, sai, &
             coszen, albg, albv, tran, thermk, extkb, extkd, ssun, ssha )
                                                                              
!-----------------------------------------------------------------------
!                                                                               
!     calculation of canopy albedos via two stream approximation (direct
!     and diffuse ) and partition of incident solar
!                                                                               
! Original author: Yongjiu Dai, June 11, 2001
!
!-----------------------------------------------------------------------        
                                                                               
  USE precision
  IMPLICIT NONE

! parameters
  REAL(r8), intent(in) :: &
          ! static parameters associated with vegetation TYPE
            chil,          &! leaf angle distribution factor
            rho(2,2),      &! leaf reflectance (iw=iband, il=life and dead)
            tau(2,2),      &! leaf transmittance (iw=iband, il=life and dead)

          ! time-space varying vegetation parameters
            green,         &! green leaf fraction
            lai,           &! leaf area index of exposed canopy (snow-free)
            sai             ! stem area index

! environmental variables
  REAL(r8), intent(in) :: &
            coszen,        &! consine of solar zenith angle
            albg(2,2)       ! albedos of ground

! output
  REAL(r8), intent(out) :: &
            albv(2,2),     &! albedo, vegetation [-]
            tran(2,2),     &! canopy transmittances for solar radiation
            thermk,        &! canopy gap fraction for tir radiation                    
            extkb,         &! (k, g(mu)/mu) direct solar extinction coefficient  
            extkd,         &! diffuse and scattered diffuse PAR extinction coefficient
            ssun(2,2),     &! sunlit canopy absorption for solar radiation
            ssha(2,2)       ! shaded canopy absorption for solar radiation,
                            ! normalized by the incident flux 
                                                                               
!-------------------------- local -----------------------------------           
  REAL(r8) :: &
           lsai,           &! lai+sai
           phi1,           &! (phi-1)
           phi2,           &! (phi-2)
           scat,           &! (omega)        
           proj,           &! (g(mu))        
           zmu,            &! (int(mu/g(mu)) 
           zmu2,           &! (zmu * zmu)    
           as,             &! (a-s(mu))      
           upscat,         &! (omega-beta)   
           betao,          &! (beta-0)       
           psi,            &! (h)            

           be,             &! (b)            
           ce,             &! (c)            
           de,             &! (d)            
           fe,             &! (f)            

           power1,         &! (h*lai)
           power2,         &! (k*lai)
           power3,         &!  

           sigma,          &! 
           s1,             &! 
           s2,             &! 
           p1,             &! 
           p2,             &! 
           p3,             &! 
           p4,             &! 
           f1,             &! 
           f2,             &!
           h1,             &!
           h4,             &!
           m1,             &!
           m2,             &!
           m3,             &!
           n1,             &!
           n2,             &!
           n3,             &!

           hh1,            &! (h1/sigma)     
           hh2,            &! (h2)           
           hh3,            &! (h3)           
           hh4,            &! (h4/sigma)     
           hh5,            &! (h5)           
           hh6,            &! (h6)           
           hh7,            &! (h7)           
           hh8,            &! (h8)           
           hh9,            &! (h9)           
           hh10,           &! (h10)          

           eup,            &! (integral of i_up*exp(-kx) )
           edw              ! (integral of i_down*exp(-kx) )
  
  INTEGER iw                ! band loop index
  INTEGER ic                ! direct/diffuse loop index

  ! variables for modified version
  REAL(r8) :: cosz, theta, cosdif, albgblk
  REAL(r8) :: tmptau, wrho, wtau
  REAL(r8) :: s2d, extkbd, sall(2,2), q, ssun_rev
                                                                                
!-----------------------------------------------------------------------        
! projected area of phytoelements in direction of mu and 
! average inverse diffuse optical depth per unit leaf area

      phi1 = 0.5 - 0.633 * chil - 0.33 * chil * chil
      phi2 = 0.877 * ( 1. - 2. * phi1 )
         
      extkd = 0.719

      IF (abs(phi1).gt.1.e-6 .and. abs(phi2).gt.1.e-6) THEN
         zmu = 1. / phi2 * ( 1. - phi1 / phi2 * log ( ( phi1 + phi2 ) / phi1 ) )
      ELSE IF (abs(phi1).le.1.e-6) THEN
         zmu = 1./0.877
      ELSE IF (abs(phi2).le.1.e-6) THEN
         zmu = 1./(2.*phi1)
      ENDIF
      zmu2 = zmu * zmu

      lsai   = lai + sai
      power3 = lsai / zmu
      power3 = min( 50., power3 )
      power3 = max( 1.e-5, power3 )
      thermk = exp(-power3)

      tmptau = 0.5_r8 * lsai
      cosdif = - tmptau / log(exp(-0.87_r8*tmptau) / (1+0.92_r8*tmptau))

      ! black ground case
      albgblk = 1.e-6_r8

      DO iw = 1, 2    ! WAVE_BAND_LOOP                                               

      ! ic 1: incident direct; 2: incident diffuse
      DO ic = 1, 2    
   
      IF (ic == 2) THEN 
         cosz = max(0.001_r8, cosdif)
         theta = acos(cosz)
         theta = theta/3.14159*180

         theta = theta + chil*5._r8
         cosz = cos(theta/180*3.14159)
      ELSE 
         cosz = coszen
      ENDIF
 
      proj = phi1 + phi2 * cosz
      extkb = (phi1 + phi2 * cosz) / cosz

!-----------------------------------------------------------------------        
!     calculate average scattering coefficient, leaf projection and             
!     other coefficients for two-stream model.                                  
!-----------------------------------------------------------------------        

! + stem optical properties
      wtau = lai/lsai*tau(iw,1) + sai/lsai*tau(iw,2)
      wrho = lai/lsai*rho(iw,1) + sai/lsai*rho(iw,2)
      
      scat = wtau + wrho
      !scat = lai/lsai * ( tau(iw,1) + rho(iw,1) ) + &
      !       sai/lsai * ( tau(iw,2) + rho(iw,2) ) 

      as = scat / 2. * proj / ( proj + cosz * phi2 )                               
      as = as * ( 1. - cosz * phi1 / ( proj + cosz * phi2 ) * &
               log ( ( proj + cosz * phi2 + cosz * phi1 ) / ( cosz * phi1 ) ) ) 

! + stem optical properties
      ! scat ~ omega
      ! upscat ~ betail*scat
      ! betao ~ betadl
      ! scat-2.*upscat ~ rho - tau
      upscat = lai/lsai*tau(iw,1) + sai/lsai*tau(iw,2)                           
      upscat = 0.5 * ( scat + ( scat - 2. * upscat ) * &         
               (( 1. + chil ) / 2. ) ** 2 )                                     
      betao = ( 1. + zmu * extkb ) / ( scat * zmu * extkb ) * as            
                                                                                
      ! [MODI 1]
      betao = 0.5_r8 * ( scat + 1._r8/extkb*(1._r8+chil)**2/4._r8*(wrho-wtau) )/scat

!-----------------------------------------------------------------------        
!     intermediate variables identified in appendix of SE-85.                   
!-----------------------------------------------------------------------        
                                                                               
      be = 1. - scat + upscat                                                   
      ce = upscat                                                               
      de = scat * zmu * extkb * betao                                          
      fe = scat * zmu * extkb * ( 1. - betao )                                 

      psi = sqrt(be**2 - ce**2)/zmu                                            
      power1 = min( psi*lsai, 50. )                                            
      power2 = min( extkb*lsai, 50. )                                          
      s1 = exp( - power1 )                                                    
      s2 = exp( - power2 )                                                     

!-----------------------------------------------------------------------        
!     calculation of direct albedos and canopy transmittances.                  
!     albv(iw,1)     ( i-up ) 
!     tran(iw,irad)  ( i-down ) 
!-----------------------------------------------------------------------        

      p1 = be + zmu * psi
      p2 = be - zmu * psi
      p3 = be + zmu * extkb
      p4 = be - zmu * extkb

      !f1 = 1. - albg(iw,2)*p1/ce
      !f2 = 1. - albg(iw,2)*p2/ce
      f1 = 1. - albgblk*p1/ce
      f2 = 1. - albgblk*p2/ce
      
      h1 = - ( de * p4 + ce * fe )
      h4 = - ( fe * p3 + ce * de )

      sigma = ( zmu * extkb ) ** 2 + ( ce**2 - be**2 )                           

      IF (ic == 1) THEN
         s2d = s2
         extkbd = extkb
      ENDIF

      IF (abs(sigma) .gt. 1.e-10) THEN     !<======

         hh1 = h1 / sigma
         hh4 = h4 / sigma
                                                                                
         m1 = f1 * s1
         m2 = f2 / s1
         !m3 = ( albg(iw,1) - ( hh1 - albg(iw,2) * hh4 ) ) * s2
         m3 = ( albgblk - ( hh1 - albgblk * hh4 ) ) * s2
 
         n1 = p1 / ce
         n2 = p2 / ce
         n3 = - hh4

         hh2 = (m3*n2 - m2*n3) / (m1*n2 - m2*n1)
         hh3 = (m3*n1 - m1*n3) / (m2*n1 - m1*n2)

         hh5 = hh2 * p1 / ce
         hh6 = hh3 * p2 / ce

         albv(iw,ic) = hh1 + hh2 + hh3                                   
         tran(iw,ic) = hh4 * s2 + hh5 * s1 + hh6 / s1                    

         eup = hh1 * (1. - s2*s2d) / (extkbd + extkb) &
             + hh2 * (1. - s2d*s1) / (extkbd + psi) &
             + hh3 * (1. - s2d/s1) / (extkbd - psi)

         edw = hh4 * (1. - s2*s2d) / (extkbd + extkb) &
             + hh5 * (1. - s2d*s1) / (extkbd + psi) &
             + hh6 * (1. - s2d/s1) / (extkbd - psi)

      ELSE                               !<======

         m1 = f1 * s1
         m2 = f2 / s1
         !m3 = h1 / zmu2 * ( lsai + 1. / (2.*extkb) ) * s2 &
         !   + albg(iw,2) / ce * ( - h1 / (2.*extkb) / zmu2 * &
         !     ( p3*lsai + p4 / (2.*extkb) ) - de ) * s2 &
         !   + albg(iw,1) * s2
         m3 = h1 / zmu2 * ( lsai + 1. / (extkb+extkbd) ) * s2 &
            + albgblk / ce * ( - h1 / (extkb+extkbd) / zmu2 * &
              ( p3*lsai + p4 / (extkb+extkbd) ) - de ) * s2 &
            + albgblk * s2
 
         n1 = p1 / ce
         n2 = p2 / ce
         n3 = 1./ce * ( h1*p4 / ((extkb+extkbd)*(extkb+extkbd)) / zmu2 + de)

         hh2 = (m3*n2 - m2*n3) / (m1*n2 - m2*n1)
         hh3 = (m3*n1 - m1*n3) / (m2*n1 - m1*n2)

         hh5 = hh2 * p1 / ce
         hh6 = hh3 * p2 / ce

         albv(iw,ic) =  - h1 / ((extkb+extkbd)*zmu2) + hh2 + hh3
         tran(iw,ic) = 1./ce * ( -h1 / ((extkb+extkbd)*zmu2) * &
            ( p3*lsai + p4 / (extkb+extkbd) ) - de ) * s2 &
            + hh5 * s1 + hh6 / s1

         eup = (hh2 - h1/((extkb+extkbd)*zmu2)) * (1. - s2*s2d)/(extkb+extkbd) &
            + hh3 * (lsai - 0.) &
            + h1/((extkb+extkbd)*zmu2) * ( lsai*s2*s2d - (1. - s2*s2d)/(extkb+extkbd) )

         edw = (hh5 - (h1*p4/((extkb+extkbd)*(extkb+extkbd)*zmu) + de)/ce) * &
            (1. - s2*s2d) / (extkb+extkbd) &
            + hh6 * (lsai - 0.) &
            + h1*p3/(ce*(extkb+extkbd)*(extkb+extkbd)*zmu2) * &
            ( lsai*s2*s2d - (1. - s2*s2d)/(extkb+extkbd) )
      
      ENDIF                              !<======

      sall(iw,ic) = 1. - albv(iw,ic) - (1.-albgblk)*(tran(iw,ic)+s2)

      IF (ic == 1) THEN
         ssun(iw,ic) = (1.-scat) * ( 1.-s2 + 1. / zmu * (eup + edw) ) 
         !ssha(iw,ic) = scat * (1.-s2) &
         !+ ( albg(iw,2)*tran(iw,ic) + albg(iw,1)*s2 - tran(iw,ic) ) - albv(iw,ic) &
         !- tran(iw,ic) - albv(iw,ic) &
         !- ( 1. - scat ) / zmu * ( eup + edw ) 
      ELSE 
         ssun(iw,ic) = (1.-scat) * ( extkb*(1.-s2*s2d)/(extkb+extkbd) + 1. / zmu * (eup + edw) ) 
      ENDIF
      
      ssha(iw,ic) = sall(iw,ic) - ssun(iw,ic)

      ENDDO ! ic

      ! for reversed diffuse radiation back from ground
      eup = hh1 * (1._r8 - s2/s2d) / (extkb - extkbd) &
          + hh2 * (1._r8 - s1/s2d) / (psi - extkbd) &
          + hh3 * (1._r8/s1/s2d - 1._r8) / (psi + extkbd)
  
      edw = hh4 * (1._r8 - s2/s2d) / (extkb - extkbd) &
          + hh5 * (1._r8 - s1/s2d) / (psi - extkbd) &
          + hh6 * (1._r8/s1/s2d - 1._r8) / (psi + extkbd)

      ssun_rev = s2d * (1._r8 - scat) * &
         ( extkb*(1._r8-s2/s2d)/(extkb-extkbd) + 1._r8 / zmu * (eup + edw ) )

      ! ++++++++++++++++++++++++++++++++++++
      ! consider the multiple reflectance between canopy and ground
      ! ++++++++++++++++++++++++++++++++++++

      ! common ratio for geometric series
      q = albg(iw,2) * albv(iw,2)

      ! black soil case
      !sall(iw,:) = ssun(iw,:) + ssha(iw,:)
   
      DO ic = 1, 2 ! from 1 to 2, cannot be reversed
         ! ++++++++++++++++++++++++++++++++++++
         ! re-calculate the absorption, transmission and albedo
         ! for direct radiation
         
! 03/06/2020, yuan: tran原来是表示漫射流，现在把直射透射率也算在内
! 03/14/2020, yuan: TODO-done, treat soil albedo in direct/diffuse cases
         IF (ic == 1) THEN
            tran(iw,ic) = (s2d*albg(iw,1)*albv(iw,2) + tran(iw,ic)) / (1.-q)

            sall(iw,ic) = sall(iw,ic) + &
               (tran(iw,ic)*albg(iw,2) + s2d*albg(iw,1)) * sall(iw,2)
            
            albv(iw,ic) = 1. - sall(iw,ic) - &
               (1.-albg(iw,2))*tran(iw,ic) - (1.-albg(iw,1))*s2d
            
            ssun(iw,ic) = ssun(iw,ic) + &
               (tran(iw,ic)*albg(iw,2) + s2d*albg(iw,1)) * ssun_rev

            ssha(iw,ic) = sall(iw,ic) - ssun(iw,ic)

         ELSE 
            tran(iw,ic) = (s2 + tran(iw,ic)) / (1.-q)

            sall(iw,ic) = sall(iw,ic) + tran(iw,ic)*albg(iw,2)*sall(iw,2)
            albv(iw,ic) = 1. - sall(iw,ic) - (1.-albg(iw,2))*tran(iw,ic)

            ssun(iw,ic) = ssun(iw,ic) + tran(iw,ic)*albg(iw,2)*ssun_rev
            ssha(iw,ic) = sall(iw,ic) - ssun(iw,ic)
         ENDIF 

      ENDDO !ic

      End DO !iw

      ! restore extkb
      extkb = extkbd

  END SUBROUTINE twostream_mod
  
  SUBROUTINE twostream_wrap ( ipatch, coszen, albg, &
             albv, tran, thermk, extkb, extkd, ssun, ssha )
                                                                              
      USE precision
      USE PFT_Const
      USE MOD_PFTimeInvars
      USE MOD_PFTimeVars
      IMPLICIT NONE

      ! parameters
      INTEGER, intent(in) :: &
            ipatch          ! patch index

      ! environmental variables
      REAL(r8), intent(in) ::  &
            coszen,        &! consine of solar zenith angle
            albg(2,2)       ! albedos of ground

      ! output
      REAL(r8), intent(out) :: &
            albv(2,2),     &! albedo, vegetation [-]
            tran(2,2),     &! canopy transmittances for solar radiation
            thermk,        &! canopy gap fraction for tir radiation                    
            extkb,         &! (k, g(mu)/mu) direct solar extinction coefficient  
            extkd,         &! diffuse and scattered diffuse PAR extinction coefficient
            ssun(2,2),     &! sunlit canopy absorption for solar radiation
            ssha(2,2)       ! shaded canopy absorption for solar radiation,
                            ! normalized by the incident flux 
 
      INTEGER :: i, p, ps, pe
      REAL(r8), allocatable :: tran_p(:,:,:)
      REAL(r8), allocatable :: albv_p (:,:,:)

      ps = patch_pft_s(ipatch)
      pe = patch_pft_e(ipatch)

      allocate ( tran_p (2,2,ps:pe) )
      allocate ( albv_p (2,2,ps:pe) ) 

      DO i = ps, pe
         p = pftclass(i)
         IF (lai_p(i)+sai_p(i) > 1.e-6) THEN
! 03/06/2020, yuan: TODO, not twostream_mod?
            CALL twostream (chil_p(p),rho_p(:,:,p),tau_p(:,:,p),1.,lai_p(i),sai_p(i),&
               coszen,albg,albv_p(:,:,i),tran_p(:,:,i),thermk_p(i),&
               extkb_p(i),extkd_p(i),ssun_p(:,:,i),ssha_p(:,:,i)) 
         ELSE
            albv_p(:,:,i) = albg(:,:)
            ssun_p(:,:,i) = 0.
            ssha_p(:,:,i) = 0.
            tran_p(:,1,i) = 1.
            tran_p(:,2,i) = 1.
         ENDIF
      ENDDO 

      albv(1,1) = sum( albv_p(1,1,ps:pe)*pftfrac(ps:pe) )
      albv(1,2) = sum( albv_p(1,2,ps:pe)*pftfrac(ps:pe) )
      albv(2,1) = sum( albv_p(2,1,ps:pe)*pftfrac(ps:pe) )
      albv(2,2) = sum( albv_p(2,2,ps:pe)*pftfrac(ps:pe) )
      
      ssun(1,1) = sum( ssun_p(1,1,ps:pe)*pftfrac(ps:pe) )
      ssun(1,2) = sum( ssun_p(1,2,ps:pe)*pftfrac(ps:pe) )
      ssun(2,1) = sum( ssun_p(2,1,ps:pe)*pftfrac(ps:pe) )
      ssun(2,2) = sum( ssun_p(2,2,ps:pe)*pftfrac(ps:pe) )

      ssha(1,1) = sum( ssha_p(1,1,ps:pe)*pftfrac(ps:pe) )
      ssha(1,2) = sum( ssha_p(1,2,ps:pe)*pftfrac(ps:pe) )
      ssha(2,1) = sum( ssha_p(2,1,ps:pe)*pftfrac(ps:pe) )
      ssha(2,2) = sum( ssha_p(2,2,ps:pe)*pftfrac(ps:pe) )

      tran(1,1) = sum( tran_p(1,1,ps:pe)*pftfrac(ps:pe) )
      tran(1,2) = sum( tran_p(1,2,ps:pe)*pftfrac(ps:pe) )
      tran(2,1) = sum( tran_p(2,1,ps:pe)*pftfrac(ps:pe) )
      tran(2,2) = sum( tran_p(2,2,ps:pe)*pftfrac(ps:pe) )

      ! fordebug ONLY below
      IF (ssun(1,1)<0 .or. ssun(1,2)<0 .or. ssun(2,1)<0 .or. ssun(2,2)<0) THEN
         print *, ipatch
         print *, ssun
      ENDIF

      deallocate ( tran_p )

  END SUBROUTINE twostream_wrap

  SUBROUTINE albocean (oro, scv, coszrs, alb)

!-----------------------------------------------------------------------
!
! Compute surface albedos
!
! Computes surface albedos for direct/diffuse incident radiation for
! two spectral intervals:
!   s = 0.2-0.7 micro-meters
!   l = 0.7-5.0 micro-meters
!
! Albedos specified as follows:
!
! Ocean           Uses solar zenith angle to compute albedo for direct
!                 radiation; diffuse radiation values constant; albedo
!                 independent of spectral interval and other physical
!                 factors such as ocean surface wind speed.
!
! Ocean with      Surface albs specified; combined with overlying snow
!   sea ice       
!
! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
!
! yongjiu dai and xin-zhong liang (08/01/2001)
!-----------------------------------------------------------------------

   USE precision
   IMPLICIT NONE

!------------------------------Arguments--------------------------------

  REAL(r8), intent(in) :: oro       ! /ocean(0)/seaice(2) flag
  REAL(r8), intent(in) :: scv       ! snow water equivalent) [mm]
  REAL(r8), intent(in) :: coszrs    ! Cosine solar zenith angle

  REAL(r8), intent(out) :: alb(2,2) ! srf alb for direct (diffuse) rad 0.2-0.7 micro-ms
                                    ! Srf alb for direct (diffuse) rad 0.7-5.0 micro-ms

!---------------------------Local variables-----------------------------

  REAL(r8) frsnow       ! horizontal fraction of snow cover
  REAL(r8) snwhgt       ! physical snow height
  REAL(r8) rghsnw       ! roughness for horizontal snow cover fractn

  REAL(r8) sasdir       ! snow alb for direct rad  0.2-0.7 micro-ms
  REAL(r8) saldir       ! snow alb for direct rad  0.7-5.0 micro-ms
  REAL(r8) sasdif       ! snow alb for diffuse rad  0.2-0.7 micro-ms
  REAL(r8) saldif       ! snow alb for diffuse rad  0.7-5.0 micro-ms

  REAL(r8), parameter :: asices = 0.70 ! sea ice albedo for 0.2-0.7 micro-meters [-]
  REAL(r8), parameter :: asicel = 0.50 ! sea ice albedo for 0.7-5.0 micro-meters [-]
  REAL(r8), parameter :: asnows = 0.95 ! snow    albedo for 0.2-0.7 micro-meters [-]
  REAL(r8), parameter :: asnowl = 0.70 ! snow    albedo for 0.7-5.0 micro-meters 

!-----------------------------------------------------------------------
! initialize all ocean/sea ice surface albedos to zero

      alb(:,:) = 0.
      IF(coszrs<=0.0) RETURN

      IF(nint(oro)==2)THEN
        alb(1,1) = asices
        alb(2,1) = asicel
        alb(1,2) = alb(1,1) 
        alb(2,2) = alb(2,1)
        sasdif = asnows
        saldif = asnowl

        IF(scv>0.)THEN
          IF (coszrs<0.5) THEN
          ! zenith angle regime 1 ( coszrs < 0.5 ).
          ! set direct snow albedos (limit to 0.98 max)
            sasdir = min(0.98,sasdif+(1.-sasdif)*0.5*(3./(1.+4.*coszrs)-1.))
            saldir = min(0.98,saldif+(1.-saldif)*0.5*(3./(1.+4.*coszrs)-1.))
          ELSE
          ! zenith angle regime 2 ( coszrs >= 0.5 )
            sasdir = asnows
            saldir = asnowl
          ENDIF

        ! compute both diffuse and direct total albedos
          snwhgt = 20.*scv / 1000.
          rghsnw = 0.25
          frsnow = snwhgt/(rghsnw+snwhgt)
          alb(1,1) = alb(1,1)*(1.-frsnow) + sasdir*frsnow
          alb(2,1) = alb(2,1)*(1.-frsnow) + saldir*frsnow
          alb(1,2) = alb(1,2)*(1.-frsnow) + sasdif*frsnow
          alb(2,2) = alb(2,2)*(1.-frsnow) + saldif*frsnow
        ENDIF
      ENDIF 

! ice-free ocean albedos function of solar zenith angle only, and
! independent of spectral interval:

      IF(nint(oro)==0)THEN
        alb(2,1) = .026/(coszrs**1.7+.065) & 
                 + .15*(coszrs-0.1)*(coszrs-0.5)*(coszrs-1.) 
        alb(1,1) = alb(2,1) 
        alb(1,2) = 0.06
        alb(2,2) = 0.06
      ENDIF
   
  END SUBROUTINE albocean

END MODULE ALBEDO
! --------- EOP ----------
