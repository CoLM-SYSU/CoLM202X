
MODULE ASSIM_STOMATA_conductance

!-----------------------------------------------------------------------
  use precision
  IMPLICIT NONE
  SAVE

! PUBLIC MEMBER FUNCTIONS:
  public :: stomata

! PRIVATE MEMBER FUNCTIONS:
  private :: sortin


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------



  subroutine stomata (vmax25,effcon,slti,hlti,shti, &
                      hhti,trda,trdm,trop,gradm,binter,tm, &
                      psrf,po2m,pco2m,pco2a,ea,ei,tlef,par, &
                      rb,ra,rstfac,cint,assim,respc,rst)  

!=======================================================================        
!                                                                               
!     calculation of canopy photosynthetic rate using the integrated            
!     model relating assimilation and stomatal conductance.                     
!
!     Original author: Yongjiu Dai, 08/11/2001
!
!     Reference: Dai et al., 2004: A two-big-leaf model for canopy temperature,
!         photosynthesis and stomatal conductance. J. Climate, 17: 2281-2299.
!
!
!     units are converted from mks to biological units in this routine.         
!                                                                               
!                          units                                                
!                         -------                                               
!                                                                               
!      pco2m, pco2a, pco2i, po2m                : pascals                       
!      co2a, co2s, co2i, h2oa, h2os, h2oa       : mol mol-1                     
!      vmax25, respcp, assim, gs, gb, ga        : mol m-2 s-1                   
!      effcon                                   : mol co2 mol quanta-1          
!      1/rb, 1/ra, 1/rst                        : m s-1                         
!                                                                               
!                       conversions                                             
!                      -------------                                            
!                                                                               
!      1 mol h2o           = 0.018 kg                                           
!      1 mol co2           = 0.044 kg                                           
!      h2o (mol mol-1)     = ea / psrf ( pa pa-1 )                              
!      h2o (mol mol-1)     = q*mm/(q*mm + 1)                                    
!      gs  (co2)           = gs (h2o) * 1./1.6                                  
!      gs  (mol m-2 s-1 )  = gs (m s-1) * 44.6*tf/t*p/po                        
!      par (mol m-2 s-1 )  = par(w m-2) * 4.6*1.e-6                             
!      mm  (molair/molh2o) = 1.611                                              
!                                                                               
!----------------------------------------------------------------------         

 use precision
 IMPLICIT NONE                                                                        

 real(r8),intent(in) :: &
      effcon,       &! quantum efficiency of RuBP regeneration (mol CO2 / mol quanta)
      vmax25,       &! maximum carboxylation rate at 25 C at canopy top 
                     ! the range : 30.e-6 <-> 100.e-6 (mol co2 m-2 s-1)

      slti,         &! slope of low temperature inhibition function      (0.2)
      hlti,         &! 1/2 point of low temperature inhibition function  (288.16)
      shti,         &! slope of high temperature inhibition function     (0.3)
      hhti,         &! 1/2 point of high temperature inhibition function (313.16)
      trda,         &! temperature coefficient in gs-a model             (1.3)
      trdm,         &! temperature coefficient in gs-a model             (328.16)
      trop,         &! temperature coefficient in gs-a model             (298.16)
      gradm,        &! conductance-photosynthesis slope parameter
      binter         ! conductance-photosynthesis intercept

 real(r8),intent(in) :: &
      tm,           &! atmospheric air temperature (K)
      psrf,         &! surface atmospheric pressure (pa)
      po2m,         &! O2 concentration in atmos. (pascals)
      pco2m,        &! CO2 concentration in atmos. (pascals)
      pco2a,        &! CO2 concentration in canopy air space (pa)
      ea,           &! canopy air space vapor pressure (pa)
      ei,           &! saturation h2o vapor pressure in leaf stomata (pa)
      tlef,         &! leaf temperature (K)
      par,          &! photosynthetic active radiation (W m-2)

      rb,           &! boundary resistance from canopy to cas (s m-1)
      ra,           &! aerodynamic resistance from cas to refence height (s m-1)
      rstfac         ! canopy resistance stress factors to soil moisture                         

 real(r8),intent(in), dimension(3) :: &
      cint           ! scaling up from leaf to canopy

 real(r8),intent(out) :: &! ATTENTION : all for canopy not leaf
      assim,        &! canopy assimilation rate (mol m-2 s-1)                       
      respc,        &! canopy respiration (mol m-2 s-1) 
      rst            ! canopy stomatal resistance (s m-1)

!-------------------- local --------------------------------------------       

 integer, parameter :: iterationtotal = 6   ! total iteration number in pco2i calculation
                                                                               
 real(r8) c3,       &! c3 vegetation : 1; 0 for c4
      c4,           &! c4 vegetation : 1; 0 for c3
      qt,           &! (tleaf - 298.16) / 10
      gammas,       &! co2 compensation point (pa)
      kc,           &! Michaelis-Menten constant for co2
      ko,           &! Michaelis-Menten constant for o2
      rrkk,         &! kc (1+o2/ko)

      templ,        &! intermediate value
      temph,        &! intermediate value
             
      vm,           &! maximum catalytic activity of Rubison (mol co2 m-2 s-1)
      jmax25,       &! potential rate of whole-chain electron transport at 25 C
      jmax,         &! potential rate of whole-chain electron transport (mol electron m-2 s-1)
      epar,         &! electron transport rate (mol electron m-2 s-1)
      respcp,       &! respiration fraction of vmax (mol co2 m-2 s-1)
      bintc,        &! residual stomatal conductance for co2 (mol co2 m-2 s-1)

      rgas,         &! universal gas contant (8.314 J mol-1 K-1)
      tprcor,       &! coefficient for unit transfer 
      gbh2o,        &! one side leaf boundary layer conductance (mol m-2 s-1)
      gah2o,        &! aerodynamic conductance between cas and reference height (mol m-2 s-1)
      gsh2o,        &! canopy conductance (mol m-2 s-1)                         

      atheta,       &! wc, we coupling parameter
      btheta,       &! wc & we, ws coupling parameter
      omss,         &! intermediate calcuation for oms
      omc,          &! rubisco limited assimilation (omega-c: mol m-2 s-1)               
      ome,          &! light limited assimilation (omega-e: mol m-2 s-1)                 
      oms,          &! sink limited assimilation (omega-s: mol m-2 s-1)                  
      omp,          &! intermediate calcuation for omc, ome

      co2m,         &! co2 concentration in atmos (mol mol-1)
      co2a,         &! co2 concentration at cas (mol mol-1)
      co2s,         &! co2 concentration at canopy surface (mol mol-1)             
      co2st,        &! co2 concentration at canopy surface (mol mol-1)
      co2i,         &! internal co2 concentration (mol mol-1)
      pco2in,       &! internal co2 concentration at the new iteration (pa)
      pco2i,        &! internal co2 concentration (pa) 
      es,           &! canopy surface h2o vapor pressure (pa)             

      sqrtin,       &! intermediate calculation for quadratic
      assmt,        &! net assimilation with a positive limitation (mol co2 m-2 s-1)
      assimn,       &! net assimilation (mol co2 m-2 s-1)
      hcdma,        &! a-1
      aquad,        &! a: ax^2 + bx + c = 0
      bquad,        &! b: ax^2 + bx + c = 0
      cquad          ! c: ax^2 + bx + c = 0

 real(r8) :: &
      eyy(iterationtotal),    &! differnce of pco2i at two iteration step
      pco2y(iterationtotal),  &! adjusted to total iteration number
      range                    !

 integer ic

!=======================================================================        

      c3 = 0.
      if( effcon .gt. 0.07 ) c3 = 1.
      c4 = 1. - c3 

!-----------------------------------------------------------------------        
! dependence on leaf temperature 
!     gammas - CO2 compensation point in the absence of day respiration 
!     ko     - Michaelis-Menton constant for carboxylation by Rubisco 
!     kc     - Michaelis-Menton constant for oxygenation by Rubisco
!-----------------------------------------------------------------------        

      qt = 0.1*( tlef - trop )

      kc = 30.     * 2.1**qt  
      ko = 30000.  * 1.2**qt 
      gammas = 0.5 * po2m / (2600. * 0.57**qt) * c3        ! = 0. for c4 plant ???

      rrkk = kc * ( 1. + po2m/ko ) * c3

!----------------------------------------------------------------------         
! maximun capacity 
! vm     - maximum catalytic activity of Rubisco in the presence of 
!          saturating level of RuP2 and CO2 (mol m-2s-1)
! jmax   - potential rate of whole-chain electron transport (mol m-2s-1)
! epar   - electron transport rate for a given absorbed photon radiation
! respc  - dark resipration (mol m-2s-1)
! omss   - capacity of the leaf to export or utilize the products of photosynthesis.
! binter - coefficient from observation, 0.01 for c3 plant, 0.04 for c4 plant 
!-----------------------------------------------------------------------        

      vm = vmax25 * 2.1**qt        ! (mol m-2 s-1) 
      templ = 1. + exp(slti*(hlti-tlef))
      temph = 1. + exp(shti*(tlef-hhti))
      vm = vm / temph * rstfac * c3 + vm / (templ*temph) * rstfac * c4
      vm = vm * cint(1)

      rgas = 8.314                 ! universal gas constant (J mol-1 K-1)
!---> jmax25 = 2.39 * vmax25 - 14.2e-6        ! (mol m-2 s-1)
!--->      jmax25 = 2.1 * vmax25        ! (mol m-2 s-1)
!/05/2014/
      jmax25 = 1.97 * vmax25       ! (mol m-2 s-1)   
      jmax = jmax25 * exp( 37.e3 * (tlef - trop) / (rgas*trop*tlef) ) * & 
             ( 1. + exp( (710.*trop-220.e3)/(rgas*trop) ) ) / &          
             ( 1. + exp( (710.*tlef-220.e3)/(rgas*tlef) ) ) 
                                   ! 37000  (J mol-1)
                                   ! 220000 (J mol-1)
                                   ! 710    (J K-1)

      jmax = jmax * rstfac
      jmax = jmax * cint(2)

!--->      epar = min(4.6e-6 * par * effcon, 0.25*jmax)
! /05/2014/
      epar = min(4.6e-6 * par * effcon, jmax)

      respcp = 0.015 * c3 + 0.025 * c4 
      respc = respcp * vmax25 * 2.0**qt / ( 1. + exp( trda*(tlef-trdm )) ) * rstfac
!     respc = 0.7e-6 * 2.0**qt / ( 1. + exp( trda*(tlef-trdm )) ) * rstfac
      respc = respc * cint(1)

      omss = ( vmax25/2. ) * (1.8**qt) / templ * rstfac * c3 &
           + ( vmax25/5. ) * (1.8**qt) * rstfac * c4 
      omss = omss * cint(1)

      bintc = binter * max( 0.1, rstfac )
      bintc = bintc * cint(3)

!-----------------------------------------------------------------------
      tprcor = 44.6*273.16*psrf/1.013e5

! one side leaf boundary layer conductance for water vapor [=1/(2*rb)]
! ATTENTION: rb in CLM is for one side leaf, but for SiB2 rb for 
! 2-side leaf, so the gbh2o shold be " 0.5/rb * tprcor/tlef "
!     gbh2o  = 0.5/rb * tprcor/tlef                   ! mol m-2 s-1
      gbh2o  = 1./rb * tprcor/tlef                    ! mol m-2 s-1

! rb is for single leaf, but here the flux is for canopy, thus
      gbh2o  = gbh2o * cint(3)

!  aerodynamic condutance between canopy and reference height atmosphere
      gah2o  = 1.0/ra * tprcor/tm                     ! mol m-2 s-1

!-----------------------------------------------------------------------        
!     first guess is midway between compensation point and maximum              
!     assimilation rate. ! pay attention on this iteration 

      co2m = pco2m/psrf                               ! mol mol-1
      co2a = pco2a/psrf

      range = pco2m * ( 1. - 1.6/gradm ) - gammas

      do ic = 1, iterationtotal    ! loop for total iteration number
      pco2y(ic) = 0.
      eyy(ic) = 0.
      enddo

      ITERATION_LOOP: do ic = 1, iterationtotal

      call sortin(eyy, pco2y, range, gammas, ic, iterationtotal)
      pco2i =  pco2y(ic)

!-----------------------------------------------------------------------
!                      NET ASSIMILATION
!     the leaf assimilation (or gross photosynthesis) rate is described 
!     as the minimum of three limiting rates:
!     omc: the efficiency of the photosynthetic enzyme system (Rubisco-limited);
!     ome: the amount of PAR captured by leaf chlorophyll;
!     oms: the capacity of the leaf to export or utilize the products of photosynthesis.
!     to aviod the abrupt transitions, two quadratic equations are used:
!             atheta*omp^2 - omp*(omc+ome) + omc*ome = 0
!         btheta*assim^2 - assim*(omp+oms) + omp*oms = 0
!-----------------------------------------------------------------------
 
      atheta = 0.877
      btheta = 0.95

      omc = vm   * ( pco2i-gammas ) / ( pco2i + rrkk ) * c3 + vm * c4 
      ome = epar * ( pco2i-gammas ) / ( pco2i+2.*gammas ) * c3 + epar * c4 
      oms   = omss * c3 + omss*pco2i * c4 

      sqrtin= max( 0., ( (ome+omc)**2 - 4.*atheta*ome*omc ) ) 
      omp   = ( ( ome+omc ) - sqrt( sqrtin ) ) / ( 2.*atheta )
      sqrtin= max( 0., ( (omp+oms)**2 - 4.*btheta*omp*oms ) )
      assim = ( ( oms+omp ) - sqrt( sqrtin ) ) / ( 2.*btheta )

      assimn= ( assim - respc)                         ! mol m-2 s-1

!-----------------------------------------------------------------------
!                      STOMATAL CONDUCTANCE
!
!  (1)   pathway for co2 flux
!                                                  co2m
!                                                   o
!                                                   |   
!                                                   |    
!                                                   <  |
!                                        1.37/gsh2o >  |  Ac-Rd-Rsoil
!                                                   <  v
!                                                   |
!                                     <--- Ac-Rd    |
!     o------/\/\/\/\/\------o------/\/\/\/\/\------o
!    co2i     1.6/gsh2o     co2s    1.37/gbh2o     co2a
!                                                   | ^                        
!                                                   | | Rsoil                        
!                                                   | |                        
!
!  (2)   pathway for water vapor flux 
!
!                                                  em
!                                                   o
!                                                   |   
!                                                   |    
!                                                   <  ^
!                                           1/gsh2o >  | Ea
!                                                   <  |
!                                                   |
!                                     ---> Ec       !
!     o------/\/\/\/\/\------o------/\/\/\/\/\------o
!     ei       1/gsh2o      es       1/gbh2o       ea
!                                                   | ^                        
!                                                   | | Eg
!                                                   | |                        
!
!  (3)   the relationship between net assimilation and tomatal conductance :
!        gsh2o = m * An * [es/ei] / [pco2s/p] + b
!        es = [gsh2o *ei + gbh2o * ea] / [gsh2o + gbh2o]
!        ===>
!        a*gsh2o^2 + b*gsh2o + c = 0
!
!-----------------------------------------------------------------------

!-->  co2a = co2m - 1.37/max(0.446,gah2o) * (assimn - 0.)     ! mol mol-1
      co2s = co2a - 1.37*assimn/gbh2o                  ! mol mol-1

      co2st = min( co2s, co2a )
      co2st = max( co2st,1.e-5 )
                               
      assmt = max( 1.e-12, assimn )
      hcdma = ei*co2st / ( gradm*assmt ) 

      aquad = hcdma                     
      bquad = gbh2o*hcdma - ei - bintc*hcdma
      cquad = -gbh2o*( ea + hcdma*bintc )  
                                             
      sqrtin= max( 0., ( bquad**2 - 4.*aquad*cquad ) )
      gsh2o = ( -bquad + sqrt ( sqrtin ) ) / (2.*aquad) 

      es  = ( gsh2o-bintc ) * hcdma                   ! pa 
      es  = min( es, ei )
      es  = max( es, 1.e-2)

      gsh2o = es/hcdma + bintc                        ! mol m-2 s-1                                      

      pco2in = ( co2s - 1.6 * assimn / gsh2o )*psrf   ! pa
      eyy(ic) = pco2i - pco2in                        ! pa

!-----------------------------------------------------------------------        

      if( abs(eyy(ic)) .lt. 0.1 ) exit

      enddo ITERATION_LOOP

! convert gsh2o (mol m-2 s-1) to resistance rst ( s m-1)
      rst   = min( 1.e6, 1./(gsh2o*tlef/tprcor) )     ! s m-1                  


  end subroutine stomata



  subroutine sortin( eyy, pco2y, range, gammas, ic, iterationtotal )

!-----------------------------------------------------------------------        
!     arranges successive pco2/error pairs in order of increasing pco2.         
!     estimates next guess for pco2 using combination of linear and             
!     quadratic fits.                                                           
!
!     original author: P. J. Sellers (SiB2)
!-----------------------------------------------------------------------        

      use precision
      IMPLICIT NONE
 
      integer, intent(in) :: ic,iterationtotal
      real(r8), INTENT(in) :: range      
      real(r8), INTENT(in) :: gammas     
      real(r8), INTENT(inout), dimension(iterationtotal) :: eyy, pco2y

!----- Local -----------------------------------------------------------        
      integer i, j, n, i1, i2, i3, is, isp, ix
      real(r8) a, b, pmin, emin, eyy_a
      real(r8) pco2b, pco2yl, pco2yq
      real(r8) ac1, ac2, bc1, bc2, cc1, cc2
      real(r8) bterm, aterm, cterm

!-----------------------------------------------------------------------        
                                                                               
      if( ic .ge. 4 ) go to 500                                                
      eyy_a = 1.0
      if(eyy(1).lt.0.) eyy_a = -1.0
      pco2y(1) = gammas + 0.5*range                                            
      pco2y(2) = gammas + range*( 0.5 - 0.3*eyy_a )                 
      pco2y(3) = pco2y(1) - (pco2y(1)-pco2y(2))/(eyy(1)-eyy(2)+1.e-10)*eyy(1)  
                                                                               
      pmin = min( pco2y(1), pco2y(2) )                                         
      emin = min(   eyy(1),   eyy(2) )                                         
      if ( emin .gt. 0. .and. pco2y(3) .gt. pmin ) pco2y(3) = gammas           
      go to 200                                                                
500   continue                                                                 
                                                                               
      n = ic - 1                                                               
      do 1000 j = 2, n                                                         
      a = eyy(j)                                                               
      b = pco2y(j)                                                             
      do 2000 i = j-1,1,-1                                                     
      if(eyy(i) .le. a ) go to 100                                             
      eyy(i+1) = eyy(i)                                                        
      pco2y(i+1) = pco2y(i)                                                    
2000  continue                                                                 
      i = 0                                                                    
100   eyy(i+1) = a                                                             
      pco2y(i+1) = b                                                           
1000  continue                                                                 
                                                                               
      pco2b = 0.                                                               
      is    = 1                                                                
      do 3000 ix = 1, n                                                        
      if( eyy(ix) .lt. 0. ) pco2b = pco2y(ix)                                  
      if( eyy(ix) .lt. 0. ) is = ix                                            
3000  continue                                                                 
      i1 = is-1                                                                
      i1 = max(1, i1)                                                          
      i1 = min(n-2, i1)                                                        
      i2 = i1 + 1                                                              
      i3 = i1 + 2                                                              
      isp   = is + 1                                                           
      isp = min( isp, n )                                                      
      is = isp - 1                                                             
                                                                               
      pco2yl=pco2y(is) - (pco2y(is)-pco2y(isp))/(eyy(is)-eyy(isp)+1.e-10)*eyy(is)     
                                                                               
!----------------------------------------------------------------------       
!   method using a quadratic fit                                                
!----------------------------------------------------------------------         
                                                                               
      ac1 = eyy(i1)*eyy(i1) - eyy(i2)*eyy(i2)                                  
      ac2 = eyy(i2)*eyy(i2) - eyy(i3)*eyy(i3)                                  
      bc1 = eyy(i1) - eyy(i2)                                                  
      bc2 = eyy(i2) - eyy(i3)                                                  
      cc1 = pco2y(i1) - pco2y(i2)                                              
      cc2 = pco2y(i2) - pco2y(i3)                                              
      bterm = (cc1*ac2-cc2*ac1)/(bc1*ac2-ac1*bc2+1.e-10)                              
      aterm = (cc1-bc1*bterm)/(ac1+1.e-10)                                              
      cterm = pco2y(i2) - aterm*eyy(i2)*eyy(i2) - bterm*eyy(i2)                
      pco2yq= cterm                                                            
      pco2yq= max( pco2yq, pco2b )                                             
      pco2y(ic) = ( pco2yl+pco2yq)/2.                                          
                                                                               
200   continue                                                                 
                                                                               
      pco2y(ic) = max ( pco2y(ic), 0.01 )                                      
                                                                               
  end subroutine sortin                                                    



END MODULE ASSIM_STOMATA_conductance
! -------------- EOP ---------------
