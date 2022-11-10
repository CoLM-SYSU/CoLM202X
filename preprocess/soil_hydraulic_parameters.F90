!=================================================================
SUBROUTINE soil_hydraulic_parameters( sand,clay,SOC,BD,soildepth, &
                          theta_s_l,psi_s_l,lambda_l,k_s_l )
!=================================================================
! (1) Dai et al., 2013, Development of a China Dataset of Soil 
!     Hydraulic Parameters Using Pedotransfer Functions for Land Surface Modeling. 
!     J. of Hydrometeorology, 14: 869-887. DOI: 10.1175/JHM-D-12-0149.1
! (2) Dai et al.,2014: Implementation of a New Global Soil Dataset in the Common Land Model.
!
! Created by Yongjiu Dai, 12/2013
! ----------------------------------------------------------------
use precision

IMPLICIT NONE

real(r8), intent(in) :: sand ! percent sand particle-size distribution (%, weight)
real(r8), intent(in) :: clay ! percent clay particle-size distribution (%, weight)
real(r8), intent(in) :: SOC  ! soil organic carbon concentration (%, weight)
real(r8), intent(in) :: BD   ! bulk density (g cm-3)
real(r8), intent(in) :: soildepth   ! soil depth (cm) 

real(r8), intent(out) :: theta_s_l  ! saturated water content (cm3/cm3)
real(r8), intent(out) :: psi_s_l    ! matric potential at saturation (cm)
real(r8), intent(out) :: lambda_l   ! pore size distribution index (dimensionless)
real(r8), intent(out) :: k_s_l      ! saturated hydraulic conductivity (cm/day)

real(r8) theta_s(18), psi_s(7), lambda(7), k_s(20)
real(r8) SOM, TOPSOIL  

! --------------------------------------------------------------------

      if(soildepth < 30.)then   ! cm
         TOPSOIL=1.  ! 1 if A or E horizon
      else
         TOPSOIL=0.  ! 0 if the subsoil
      endif

      SOM=1.72*SOC    

      CALL thetas(BD,SOM,SOC,sand,clay,TOPSOIL,theta_s)
      CALL debar(theta_s,18,theta_s_l)

      CALL ClappH(BD,sand,clay,SOM,SOC,theta_s_l,psi_s,lambda)
      CALL debar(psi_s,7,psi_s_l)
      CALL debar(lambda,7,lambda_l)

      CALL ksat(BD,SOM,SOC,sand,clay,topsoil,theta_s_l,psi_s_l,lambda_l,k_s)
      CALL debar(k_s,20,k_s_l)

END SUBROUTINE soil_hydraulic_parameters


!=====================================================
SUBROUTINE thetas(BD,SOM,SOC,sand,clay,topsoil,theta_s)
!=====================================================
! (1) Dai et al., 2013, Development of a China Dataset of Soil 
!     Hydraulic Parameters Using Pedotransfer Functions for Land Surface Modeling. 
!     J. of Hydrometeorology, 14: 869-887. DOI: 10.1175/JHM-D-12-0149.1
! (2) Dai et al.,2014: Implementation of a New Global Soil Dataset in the Common Land Model.
!
! Created by Yongjiu Dai, 12/2013
! ----------------------------------------------------
use precision

IMPLICIT NONE

real(r8), intent(in) :: BD   ! bulk density (g cm-3)
real(r8), intent(in) :: SOM  ! percent organic matter (%, weight)
real(r8), intent(in) :: SOC  ! soil organic carbon concentration (%, weight)
real(r8), intent(in) :: sand ! percent sand particle-size distribution (%)
real(r8), intent(in) :: clay ! percent clay particle-size distribution (%)
real(r8), intent(in) :: TOPSOIL  ! 1 if A or E horizon, 0 if the subsoil

real(r8), intent(out) :: theta_s(18)

real(r8) phi ! soil porosity (cm3 cm-3)
real(r8) silt     ! percent silt particle-size distribution (%)
real(r8) x,z1,z2
integer i,j,iclass
real(r8) params(7,10)
data params/0.025,0.403,0.0383,1.3774,0.2740, 1.2500,60.000, &
            0.010,0.439,0.0314,1.1804,0.1528,-2.3421,12.061, &
            0.010,0.430,0.0083,1.2539,0.2025,-0.5884, 2.272, &
            0.010,0.520,0.0367,1.1012,0.0919,-1.9772,24.800, &
            0.010,0.614,0.0265,1.1033,0.0936, 2.5000,15.000, &
            0.025,0.366,0.0430,1.5206,0.3424, 1.2500,70.000, &
            0.010,0.392,0.0249,1.1689,0.1445,-0.7437,10.755, &
            0.010,0.412,0.0082,1.2179,0.1789, 0.5000, 4.000, &
            0.010,0.481,0.0198,1.0861,0.0793,-3.7124, 8.500, &
            0.010,0.538,0.0168,1.0730,0.0680, 0.0001, 8.235/

      silt=max(1.,100.0-sand-clay)
      phi=1.0-BD/2.65

! 1) Saturated water content = total porosity:
      theta_s(1)=1.0-BD/2.65

! ---------------------------------------------------------
! 2) Williams et al. (1992):
! Williams, J., P. Ross, and K. Bristow. 1992. Prediction of the Campbell water
! retention function from texture, structure, and organic matter. p. 427-442.
! In M.Th . van Genuchten et al. (ed.) Proc. Int. Worksh. on Indirect Methods
! for Estimating the Hydraulic Properties of Unsaturated Soils, Riverside,
! CA. 11-13 Oct. 1989. U.S. Salinity Lab., Riverside, CA.
!
! based on large soil data sets from Australia, UK, and USA
! set the value of the residual water content equal to zero
! function # 1 all sand is fine sand and the structural index = 1
! ---------------------------------------------------------
      theta_s(2)=0.93*(1.0-BD/2.65)

! -------------------------------------------
! 3) Jauhiainen (2004):
! Jauhiainen, M., Relationships of particle size distribution curve, soil water retention curve 
! and unsaturated hydraulic conductivity and their implications on water balance of forested 
! and agricultural hillslopes, Helsinki University of Technology, Water Resources Publications, 
! TKK-VTR-12, Espoo, 2004, 165 pp.
! -------------------------------------------
      theta_s(3)=0.928*(1.0-BD/2.65)+0.021

! -------------------------------------------
! 4) Teepe et al. (2003):
! Teepe, R., H. Dilling, and F. Beese, 2003: Estimating water retention curves of forest soils
! from soil texture and bulk density. J. Plant Nutr. Soil Sci., 166, 11-119.
!
! the data of this study was a collective of 1850 water retention curve 
! from FOREST soils in German
! -------------------------------------------
      theta_s(4)=0.9786-0.36686*BD

! ------------------------------------------------------------
! 5) Cosby et al. (1984):
! B.J. COSBY, G.M. HORNBERGER, R.B. CLAPP, and T.R. GINN, 1984: A Statistical exploration
! of the Relationships of soil moisture characteristics to the physical properties of soils
! WATER RESOURCES RESEARCH, VOL. 20, NO. 6, PAGES 682-690, JUNE 1984.
!
! 1448 samples from 35 localities in 23 states in the United States.
! ------------------------------------------------------------
      theta_s(5)=0.489-0.00126*sand  ! (univariate regression)

! 6) Cosby et al. (1984):
      theta_s(6)=0.505-0.00142*sand-0.00037*clay ! (multi-variate regression)

! ------------------------------------------------------------
! 7) Saxton et al. (1986):
! Saxton, K.E., W.J. Rawls, J.S. Romberger, and R.I. Papendick. 1986.
! Estimating generalized soil water characteristics from texture.
! Soil Sci. Soc. Am. J. 50: 1031-1036
! ------------------------------------------------------------
      if((sand>=5.0) .and. (clay>=5.0 .and. clay<=60.))then
         theta_s(7)=0.332-0.0007251*sand+0.1276*log10(clay)
      else
         theta_s(7)=-1.e36
      endif

! ------------------------------------------------------------
! 8) Vereeken et al. (1989):
! Vereecken, H., J. Maes, J. Feyen, and P. Darius. 1989. Estimating the soil moisture
! retention characteristics from texture, bulk density and carbon content.
! Soil Sci. 148:389-403.
! ------------------------------------------------------------
      theta_s(8)=0.81-0.28*BD+0.001*clay

!--------------------------------------------------------------
! 9) Wosten et al. (1999) :
! Wosten, J.H.M., A. Lilly, A. Nemes, and C. Le Bas. 1999. Development and
! use of a database of hydraulic properties of European soils. Geoderma 90:169-185.
!
! Wosten et al. (1999) analyzed the all-Europe database and derived 
! the following PTFs to estimate van Genuchten parameters:
! where topsoil is an ordinal variable having the value of 1 or 
! of 0. %silt is the percentage of soil (2mm-50mm), BD= bulk density (g/cm3 = Mg/m3)
! the PTFs developd by using the European HYdraulic PRoperties of European Soils (HYPRES)
! soil database (5521 horizons)
! topsoil: 0.403 < theta_s < 0.614 (oganic: 0.766)
! subsoil: 0.366 < theta_s < 0.538
!--------------------------------------------------------------
      theta_s(9)=0.7919+0.001691*clay-0.29619*BD-0.000001491*(silt)**2 &
                +0.0000821*(SOM)**2+0.02427/(clay)+0.01113/(silt) &
                +0.01472*log(silt)-0.0000733*(SOM)*(clay)-0.000619*BD*(clay)&
                -0.001183*BD*SOM-0.0001664*(topsoil)*(silt)

! 10) Wosten et al. (1999):
! FAO class
      if(clay.GT.60.) then
         iclass=5
      else if(clay.GT.35.) then
         iclass=4
      else if(sand.LT.15.) then
         iclass=3
      else if(sand.GT.65.0.AND.clay.LT.18.0) then
         iclass=1
      else
         iclass=2
      endif

      j=iclass+nint(topsoil)*5
      theta_s(10)=params(2,j)

! -------------------------------------------------------------
! 11) Mayr and Jarvis (1999):
! T. Mayr, N.J. Jarvis, 1999: Pedotransfer functions to estimate soil water
! retention parameters for a modified Brooks-Corey type model. Geoderma 91:1-9
!
! 286 soil horizons of the soil physical properties database of England and Wales
! set the value of the residual water content equal to zero
! these function should NOT be applied to organic soils (organic carbon content >5%) 
! and/or soils of small dry bulk density (<0.9 g/cm3)
! -------------------------------------------------------------
      if(SOC < 5. .or. BD > 0.9)then
         theta_s(11)=0.2345971971 &
                 +0.0046614221*sand+0.0088163314*silt+0.0064338641*clay-0.3028160229*BD &
                 +1.79762e-05*(sand)**2-3.134631e-05*(silt)**2
      else
         theta_s(11)=-1.0e36
      endif

! -------------------------------------------
! 12) Merdun (2010) MLR PTFs:
! H. Merdun 2010: Alternative Methods in the Development of Pedotransfer Functions
! for Soil Hydraulic Characteristics. Eurasian Soil Science,Vol. 43, No. 1, 62-71.
!
! data from the UNSODA database
! 135 soil samples for the development and 45 for the validation
! 2.0  < sand <96.0   mean 47.2     s.d. 29.7
! 1.0  < silt <90.0   mean 35.2     s.d. 22.7
! 0.1  < clay <63.3   mean 0.176    s.d. 14.5
! 0.01 < SOM  <21.4;  mean 1.673    s.d. 2.478
! 0.59 < BD   <1.76   mean 1.42     s.d. 0.22
! 0.274 < theta_s < 0.837; mean: 0.443  
! -------------------------------------------
      if(((sand>=2.0) .and. (sand<=96.0)) .AND. & 
         ((silt>=1.0) .and. (silt<=90.0)) .AND. &
         ((clay>=0.1) .and. (clay<=63.3)) .AND. &
         ((SOM>=0.01) .and. (SOM<=21.4))  .AND. &
         ((BD>=0.59)  .and. (BD<=1.76)))then
         theta_s(12)=1.391-0.289*clay/100.-1.007*BD &
                    -0.026*SOM-0.096*(clay/100.)**2+0.22*BD**2-0.00039*SOM**2 &
                    +0.3*clay/100.*BD+0.0233*clay/100.*SOM+0.0229*BD*SOM

! 13) Merdun (2010) SUR PTFs:
         theta_s(13)=1.419680-0.37696*clay/100.-1.04082*BD &
                    -0.02362*SOM+0.08548*(clay/100.)**2+0.23091*BD**2+0.00002*SOM**2+0.32229*clay/100.*BD &
                    +0.004189*clay/100.*SOM +0.022525*BD*SOM
      else
         theta_s(12:13)=-1.0e36
      endif

! ----------------------------------------------------------
! 14) Saxton and Rawls (2006):
! K. E. Saxton and W. J. Rawls, 2006: Soil Water Characteristic Estimates by Texture and 
! Organic Matter for Hydrologic Solutions. Soil Sci. Soc. Am. J. 70:1569-1578.
!
! A-horizon samples from USDA/NRCS National Soil Characterization Dataset.
! Samples with "extreme" values were omitted from the data. 
! Excluded were those with bulk density < 1.0 and > 1.8 g cm-3,
! OM > 8 % (w) and clay > 60% (w). This reduced the A-horizon data set 
! from 2149 to 1722 samples. Set the value of the residual water content equal to zero
! the fitting samples with 
! 1.0  < BD <1.8 g/cm3
! SOM  < 8  %(w)
! clay < 60 %(w).
! ----------------------------------------------------------
      if((clay<=60.0) .and. (SOM<=8.0) .and. (BD>=1.0 .and. BD<=1.8))then
         x=-0.251*sand/100.+0.195*clay/100.+0.011*SOM &
           +0.006*(sand/100.*SOM)-0.027*(clay/100.*SOM) &
           +0.452*(sand/100.*clay/100.)+0.299
         z1=x+(1.283*x**2-0.374*x-0.015)
         z2=-0.107+1.636*(0.278*sand/100.+0.034*clay/100.+0.022*SOM &
           -0.018*(sand/100.*SOM)-0.027*(clay/100.*SOM) &
           -0.584*(sand/100.*clay/100.)+0.078)
         theta_s(14)=z1+z2-0.097*sand/100.+0.043
         if(theta_s(14) <= z2) theta_s(14) = -1.0e36  ! theta_s < field capacity
      else
         theta_s(14)=-1.0e36
      endif

! ----------------------------------------------------------
! 15) Tomasella and Hodnett (1998):
! Tomasella, J., and M.G. Hodnett. 1998. Estimating soil water retention characteristics
! from limited data in Brazilian Amazonia. Soil Sci. 163:190-202.
! ----------------------------------------------------------
      theta_s(15)=(37.937+2.24*SOC+0.298*silt+0.159*clay)/100.

! ------------------------------------------
! 16) Li et al. (2007):
! Y. Li, D. Chen, R.E. White, A. Zhu, J. Zhang, 2007: Estimating soil hydraulic 
! properties of Fengqiu County soils in the North China Plain using 
! pedo-transfer functions. Geoderma, 138:261-271.
!
! the residual water content ranged from 0.002-0.0001 cm3 cm-3
! when theta_r<0.001, set theta_r=0.
! 63 samples for SWRC, 
! 8.98< sand <93.03;    mean 53.28       s.d. 4.46
! 1.74< silt <79.5      mean 35.86       s.d. 3.83
! 0.54< clay <27.12;    mean 8.86        s.d. 1.01
! 0.12< SOM  <1.54;     mean 0.65        s.d. 0.07
! 1.2 < BD   <1.59      mean 1.42        s.d. 0.01
! all 63 samples: theta_s=0.477 +/- 0.006
! ------------------------------------------
      if(((sand>=8.98) .and. (sand<=93.03)) .AND. & 
         ((silt>=1.74) .and. (silt<=79.5))  .AND. &
         ((clay>=0.54) .and. (clay<=27.12)) .AND. &
         ((SOM>=0.12)  .and. (SOM<=1.54))   .AND. &
         ((BD>=1.2)    .and. (BD<=1.59)))    then
         theta_s(16)=exp(-1.531+0.212*log(sand)+0.006*silt-0.051*SOM-0.566*log(BD))
      else
         theta_s(16)=-1.0e36
      endif

! ----------------------------------------------------------
! 17) Al Majou et al (2007)
! Al Majou, H., A. Bruand, and O. Duval. 2008. The use of in situ volumetric
! water content at fi eld capacity to improve the prediction of soil water
! retention properties. Can. J. Soil Sci. 88:533-541.
!
! Class and continuous ptfs were developed using a set of 320 horizons, comprising 90 topsoils 
! (from 0 to 30 cm depth) and 230 subsoil horizons (> 30 cm depth)
! collected in Cambisols, Luvisols, Planosols, Albeluvisols,Podzols, and Fluvisols located 
! mainly in the Paris basin and secondarily in the western coastal marshlands 
! and Pyrenean piedmont plain. 
! 0.1< sand < 90.1       mean 24.9         s.d. 23.9
! 2.8< silt < 82.1       mean 46.2         s.d. 20.8
! 1.9< clay <92.9,       mean 28.9         s.d. 15.1
! 0.0< SOC <28.8 g kg-1, mean 5.7 g kg-1   s.d. 4.9
! 1.0< BD <1.84,         mean 1.53         s.d. 0.15 
! topsoil: 0.397 < theta_s < 0.587
! subsoil: 0.367 < theta_s < 0.742
! A set of 107 horizons comprising 39 topsoil and 68 subsoil horizons
! was constituted in order to test the ptfs established. These horizons were collected in Cambisols, 
! Luvisols and Fluvisols located in the South of the Paris basin
 
! the residual water content was fixed at 0.01 cm3 cm-3 
! except for texture coarse for which it was fixed at 0.025 cm3 cm-3
! ----------------------------------------------------------
      if(((sand>=0.1) .and. (sand<=90.1)) .AND. &
         ((silt>=2.8) .and. (sand<=82.1)) .AND. &
         ((clay>=1.9) .and. (clay<=92.9)) .AND. &
         ((SOC>=0.0)  .and. (SOC<=2.88 )) .AND. &
         ((BD>=1.0)   .and. (BD<=1.84)))   then 
         theta_s(17)=1.1658-0.0032*clay-0.4737*BD+2.0e-07*(silt)**2-0.0001*(SOC)**2 &
                    +0.0373/(SOC)+0.0131/silt-0.0072*log(silt)+0.00003*(SOC)*clay+0.0022*BD*clay &
                    -0.0002*BD*(SOC)-0.0001*silt
      else
         theta_s(17)=-1.0e36
      endif

! ------------------------------------------
! 18) Weynants et al. (2009):
! M. Weynants, H. Vereecken, and M. Javaux, 2009: Revisiting Vereecken Pedotransfer
! Functions: Introducing a Closed-Form Hydraulic Model. Vadose Zone J. 8:86-95.
!
! datasets: 136 hydraulic conductivity curve and 166 moisture retention curve in the middle of 
! each of the 182 horizons from northern Belgium. Training data value ranges:
! 5.6  < sand < 97.80
! 0    < clay < 54.50
! 0.1  < SOC  < 66 (g kg-1)
! 0.89 < BD   < 1.77 (g cm-3)
! the theta_r could be 0.
! ------------------------------------------
      if(((sand>=5.6) .and. (sand<=97.8)) .AND. &
         ((clay>=0.)  .and. (clay<=54.5)) .AND. &
         ((SOC>=0.01) .and. (SOC<=6.6))   .AND. &
         ((BD>=0.89)  .and. (BD<=1.77)))then 
         theta_s(18)=0.6355+0.0013*clay-0.1631*BD
      else
         theta_s(18)=-1.0e36
      endif

      do i = 1, 18
         if(abs(theta_s(i)) > 0.8) theta_s(i)=-1.0e36
      enddo

END SUBROUTINE thetas


!======================================================
SUBROUTINE ClappH(BD,sand,clay,SOM,SOC,phi,psi_s,lambda)
!======================================================
! (1) Dai et al., 2013, Development of a China Dataset of Soil 
!     Hydraulic Parameters Using Pedotransfer Functions for Land Surface Modeling. 
!     J. of Hydrometeorology, 14: 869-887. DOI: 10.1175/JHM-D-12-0149.1
! (2) Dai et al.,2014: Implementation of a New Global Soil Dataset in the Common Land Model.
!
! Created by Yongjiu Dai, 12/2013
! ----------------------------------------------------
use precision

IMPLICIT NONE
real(r8), intent(in) :: BD   ! bulk density (g cm-3)
real(r8), intent(in) :: sand ! percent sand particle-size distribution (%, weight)
real(r8), intent(in) :: clay ! percent clay particle-size distribution (%, weight)
real(r8), intent(in) :: SOM  ! percent organic matter (%, weight)
real(r8), intent(in) :: SOC  ! percent organic carbon (%, weight)
real(r8), intent(in) :: phi  ! soil porosity (cm3 cm-3)

REAL(r8), INTENT(OUT) :: PSI_S(7)  ! MATRIC POTENTIAL AT SATURATION
REAL(r8), INTENT(OUT) :: LAMBDA(7) ! 

real(r8) silt
real(r8) d_g, sigma_g, h_es, h_b   
real(r8) a, b, x, z, z33, zs33, zs, z1500, str
integer i

! h_b= bubbling capillary pressure in centimeters
! phi= pore size index

      silt=max(1.0,100.0-sand-clay)

! GROUP I: input soil physical parameters: Texture, Bulk Density, Porosity
! ------------------------------------------------------------
! 1) Cosby et al. (1984), uni-variate regression:
! B.J. COSBY, G.M. HORNBERGER, R.B. CLAPP, and T.R. GINN, 1984: A Statistical exploration
! of the Relationships of soil moisture characteristics to the physical properties of soils
! WATER RESOURCES RESEARCH, VOL. 20, NO. 6, PAGES 682-690, JUNE 1984. 
!
! 1448 samples from 35 localities in 23 states in the United States.
! set the value of the residual water content equal to zero
! ------------------------------------------------------------
      psi_s(1)=-10.0**(1.88-0.013*sand)                                                      ! cm
      lambda(1)=1./(2.91+0.159*clay)

! 2) Cosby et al. (1984), multi-variate regression:
! 1448 samples from 35 localities in 23 states in the United States.
! set the value of the residual water content equal to zero
      psi_s(2)=-10**(1.54-0.0095*sand+0.0063*silt)                                           ! cm
      lambda(2)=1.0/(3.10+0.157*clay-0.003*sand)

! ----------------------------------------------------------
! 3) Saxton (1986):
! Saxton, K.E., W.J. Rawls, J.S. Romberger, and R.I. Papendick. 1986.
! Estimating generalized soil water characteristics from texture.
! Soil Sci. Soc. Am. J. 50: 1031-1036
!
! 5 < clay < 60
! 5 < sand
! set the value of the residual water content equal to zero
! the equation is avaliable for a wide range of soil texture.
! ----------------------------------------------------------
      if((sand>=5.0) .and. (clay>=5.0 .and. clay<=60.))then
         lambda(3)=1./(3.14+0.00222*clay**2+0.00003484*sand**2*clay)
         psi_s(3)=-100.*exp(-4.396-0.0715*clay-0.000488*sand**2-0.00004285*sand**2*clay) &
                 /((0.332-0.0007251*sand+0.1276*log10(clay))**(1./lambda(3)))                ! kPa 
         psi_s(3)= 10.*psi_s(3)   
      else
         lambda(3)=-1.0e36
         psi_s(3)=-1.0e36
      endif

! ------------------------------------------------------------
! 4) Campbell and Shiozawa (1992):
! Campbell, G.S., and S. Shiozawa. 1992. Prediction of hydraulic properties of
! soils using particle size distribution and bulk density data. p. 317-328. In
! M.Th . van Genuchten et al. (ed.) Proc. Int. Worksh. on Indirect Methods
! for Estimating the Hydraulic Properties of Unsaturated Soils, Riverside,
! CA. 11-13 Oct. 1989. U.S. Salinity Lab., Riverside, CA.
!
! set the value of the residual water content equal to zero
! the arithmetic means 1.025mm for sand (2-0.05mm),0.026mm for silt (0.05-0.002mm)
! 0.001 for clay (<0.002mm)
! ------------------------------------------------------------
      d_g=exp(log(1.025)*sand/100.+log(0.026)*silt/100.+log(0.001)*clay/100.)
      sigma_g=exp(sqrt((log(1.025))**2*sand/100.+(log(0.026))**2*silt/100. &
             +(log(0.001))**2*clay/100.-(log(1.025)*sand/100.+log(0.026)*silt/100. &
             +log(0.001)*clay/100.)**2))

      lambda(4)=1.0/(1./sqrt(d_g)+0.2*sigma_g)
      h_es=0.05/sqrt(d_g)  
      h_b=h_es*(BD/1.3)**(0.67/lambda(4))
      psi_s(4)=-100.*h_b                                                                     ! cm 

! ---------------------------------------------------------
! 5) Rawls and Brakenssiek (1998):
! W. J. Rawls, D. Gimenez, R. Grossman, 1998: USE OF SOIL TEXTURE, BULK DENSITY, 
! AND SLOPE OF THEWATER RETENTION CURVE TO PREDICT SATURATED HYDRAULIC CONDUCTIVITY
! Transactions of the ASAE. VOL. 41(4): 983-988
!
! PTFs were developed from about 18,000 samples from the national soil characterization
! database that the clay activity ratio along with texture to characterized 
! the effect of soil minerals on soil hydraulic properties.
! the value of the residual water content NOT equal zero  
! the equation for theta_r was given in the subroutine thetar
! ------------------------------------------------------------ 
      h_b=exp(5.3396738+0.1845038*clay-2.48394546*phi-0.00213853*clay**2 &
         -0.04356349*sand*phi-0.61745089*clay*phi+0.00143598*(sand*phi)**2 &
         -0.00855375*(clay*phi)**2-0.00001282*sand**2*clay+0.00895359*clay**2*phi &
         -0.00072472*(sand)**2*phi+0.0000054*clay**2*sand+0.5002806*clay*phi**2)             ! cm
      psi_s(5)=-h_b                                                                          ! cm

      lambda(5)=exp(-0.7842831+0.0177544*sand-1.062498*phi-0.00005304*sand**2 &
               -0.00273493*clay**2+1.11134946*phi**2-0.03088295*sand*phi &
               +0.00026587*(sand*phi)**2-0.00610522*(clay*phi)**2 &
               -0.00000235*sand**2*clay+0.00798746*clay**2*phi-0.0067449*clay*phi**2)

! GROUP II: input soil physical parameters: Texture, Bulk Density and Organic Matter/Carbon
! -------------------------------------------------------------
! 6) Mayr and Jarvis (1999):
! T. Mayr, N.J. Jarvis, 1999: Pedotransfer functions to estimate soil water
! retention parameters for a modified Brooks-Corey type model. Geoderma 91:1-9.
!
! 286 soil horizons of the soil physical properties database of England and Wales
! set the value of the residual water content equal to zero
! these function should NOT be applied to organic soils (organic carbon content >5%) 
! and/or soils of small dry bulk density (<0.9 g/cm3)
! -------------------------------------------------------------
      if(SOC < 5. .or. BD > 0.9)then
      psi_s(6)=-10.0**(-4.9840297533+0.0509226283*sand+0.1575152771*silt &
              +0.1240901644*BD-0.1640033143*SOC &
              -0.0021767278*silt**2+0.00001438224*silt**3 &
              +0.0008040715*clay**2+0.0044067117*SOC**2)                                     ! cm

      lambda(6)=10.0**(-0.8466880654-0.0046806123*sand+0.0092463819*silt &
               -0.4542769707*BD-0.0497915563*SOC+3.294687e-04*sand**2-1.689056e-06*sand**3 &
               +0.0011225373*SOC**2)
      else
         lambda(6)=-1.0e36
         psi_s(6)=-1.0e36
      endif

! ----------------------------------------------------------
! 7) Saxton and Rawls(2006):
! K. E. Saxton and W. J. Rawls, 2006: Soil Water Characteristic Estimates by Texture and 
! Organic Matter for Hydrologic Solutions. Soil Sci. Soc. Am. J. 70:1569-1578.
!
! A-horizon samples from USDA/NRCS National Soil Characterization Dataset.
! samples with "extreme" values were omitted from the data. 
! excluded were those with bulk density < 1.0 and > 1.8 g cm-3,
! OM > 8 % (w) and clay > 60% (w). This reduced the A-horizon data set 
! from 2149 to 1722 samples. Set the value of the residual water content equal to zero
! the fitting samples with 
! 1.0  < BD <1.8 g/cm3
! SOM  < 8  %(w)
! clay < 60 %(w).
! ----------------------------------------------------------
      if((clay<=60.0) .and. (SOM<=8.0) .and. (BD>=1.0 .and. BD<=1.8))then
         x=-0.251*sand/100.+0.195*clay/100.+0.011*SOM &
          +0.006*(sand/100.*SOM)-0.027*(clay/100.*SOM) &
          +0.452*(sand/100.*clay/100.)+0.299
         z33=max(1.e-6, x+(1.283*x**2-0.374*x-0.015))

         zs33=-0.107+1.636*(0.278*sand/100.+0.034*clay/100.+0.022*SOM &
             -0.018*(sand/100.*SOM)-0.027*(clay/100.*SOM) &
             -0.584*(sand/100.*clay/100.)+0.078)
         zs33=max(0.001, zs33)

         zs=max(0.001, z33+zs33-0.097*sand/100.+0.043)

         z1500=-0.02+1.14*(0.024*sand/100.+0.487*clay/100.+0.006*SOM &
              +0.005*(sand/100.*SOM)-0.013*(clay/100.*SOM) &
              +0.068*(sand/100.*clay/100.)+0.031)
         z1500=max(0.001, z1500)

         lambda(7)=max(0.1, (log(z33)-log(z1500))/(log(1500.)-log(33.)))
         if(zs <= zs33 .or. zs <= z1500 .or. zs33 <= z1500)then 
             lambda(7)=-1.0e36
             psi_s(7)=-1.0e36
         else
             psi_s(7)=-10.*exp(log(33.)+1./lambda(7)*log(z33))/(zs**(1./lambda(7)))         ! cm
         endif
      else
         lambda(7)=-1.0e36
         psi_s(7)=-1.0e36
      endif

      do i = 1, 7
         if(lambda(i) > 1. .or. lambda(i) <= 0.) lambda(i)=-1.0e36
         if(psi_s(i) < -300. .or. psi_s(i) >= 0.) psi_s(i)=-1.0e36
      enddo
      

END SUBROUTINE clappH

!============================================================
SUBROUTINE ksat(BD,SOM,SOC,sand,clay,topsoil,phi,psi_l,lambda_l,k_s)
!============================================================
! (1) Dai et al., 2013, Development of a China Dataset of Soil 
!     Hydraulic Parameters Using Pedotransfer Functions for Land Surface Modeling. 
!     J. of Hydrometeorology, 14: 869-887. DOI: 10.1175/JHM-D-12-0149.1
! (2) Dai et al.,2014: Implementation of a New Global Soil Dataset in the Common Land Model.
!
! Created by Yongjiu Dai, 12/2013
! ----------------------------------------------------
use precision

IMPLICIT NONE
real(r8), intent(in) :: BD   ! bulk density (g cm-3)
real(r8), intent(in) :: SOM  ! soil organic matter (%, weight)
real(r8), intent(in) :: SOC  ! soil organic carbon concentration (%, weight)
real(r8), intent(in) :: sand ! percentage of sand particle-size distribution (%, weight)
real(r8), intent(in) :: clay ! percentage of clay particle-size distribution (%, weight)
real(r8), intent(in) :: phi  ! total porosity
real(r8), intent(in) :: psi_l      ! (cm)
real(r8), intent(in) :: lambda_l   !
real(r8), intent(in) :: TOPSOIL  ! 1 if A or E horizon, 0 if the subsoil

REAL(r8), INTENT(OUT) :: K_S (20)  ! SATURATED HYDRAULIC CONDUCTIVITY (m/d)

real(r8) :: silt              ! percentage of silt particle-size distribution
real(r8) :: theta_33          ! the water content at a potential of -33 kPa
real(r8) :: phi_e             ! phi minus theta_33
real(r8) x,B,zs,z33,zs33,z1500,lambda,d_g,sigma_g
integer i

      silt = max(1.,100.0-sand-clay)
      theta_33 = phi*(psi_l/(-33.*10.))**lambda_l     ! 1 kPa = 10 cm H2O

! GROUP I: input soil physical parameters: Texture, Bulk Density, Porosity
! ------------------------------------------------------------
! 1) Ahuja et al. (1989):
! Ahuja, L.R., D.K. Cassel, R.R. Bruce, and B.B. Barnes. 1989. Evaluation of
! spatial distribution of hydraulic conductivity using eff ective porosity data.
! Soil Sci. 148:404-411. 
! ------------------------------------------------------------
      phi_e=phi-theta_33
      k_s(1)=18348.*phi_e**3.295                                                             ! cm/d

! ------------------------------------------------------------
! 2) Suleiman et al. (2001):
! Suleiman, A. A., and J. T. Ritchie. 2000. Estimating Saturated Hydraulic Conductivity from Soil
! Porosity. American Society of Agricultural Engineers Vol. 44(2): 235-339. 
! ------------------------------------------------------------
      phi_e=phi-theta_33
      k_s(2)=12303.*phi_e**3.63                                                              ! cm/d

! ------------------------------------------------------------
! 3) Spychalski et al. (2007):
! Spychalski M., C. Kazmierowski, Z. Kaczmarek (2007): Estimation of saturated 
! hydraulic conductivity on the basis of drainage porosity. 
! Electronic J. Polish Agricultural Universities. 10(1), #04.
! for 0.003 < phi_e < 0.3, it provides a more accurate assessment.
! based on 35 samples.
! ------------------------------------------------------------
      phi_e=phi-theta_33
      k_s(3)=86400./10000.*(-2.52+581.598*phi_e**1.5-6966.14*phi_e**2.5+11693.78*phi_e**3)   ! cm/d

! ------------------------------------------------------------
! 4) Cosby et al. (1984):
! B.J. COSBY, G.M. HORNBERGER, R.B. CLAPP, and T.R. GINN, 1984: A Statistical exploration
! of the Relationships of soil moisture characteristics to the physical properties of soils
! WATER RESOURCES RESEARCH, VOL. 20, NO. 6, PAGES 682-690, JUNE 1984.
!
! 1448 samples from 35 localities in 23 states in the United States.
! ------------------------------------------------------------
! Multiple regression:
      k_s(4)= 60.96*10.0**(-0.6+0.0126*sand-0.0064*clay)                                     ! cm/d

! 5) Cosby et al. (1984):
! Univariate regression: 
      k_s(5)=60.96*10.0**(-0.884+0.0153*sand)                                                ! cm/d

! ----------------------------------------------------------
! 6) Saxton et al. (1986):
! Saxton, K.E., W.J. Rawls, J.S. Romberger, and R.I. Papendick. 1986.
! Estimating generalized soil water characteristics from texture.
! Soil Sci. Soc. Am. J. 50: 1031-1036
!
! 230 selected data points uniformly spaced on the 10 curves
! the equation is avaliable for a wide range of soil texture.
! ----------------------------------------------------------
      if((sand>=5.0) .and. (clay>=5.0 .and. clay<=60.))then
         k_s(6)=24.0*exp(12.012-0.0755*sand &
               +(-3.895+0.03671*sand-0.1103*clay+0.00087546*(clay)**2) &
               /(0.332-0.0007251*sand+0.1276*log10(clay)))                                   ! cm/d
      else
         k_s(6)=-1.0e36
      endif

! ----------------------------------------------------------
! 7) Dane and Puckett (1994): 
! Dane, J.H., W. Puckett. 1994. Field soil hydraulic properties based on physical and
! mineralogical information. In ¡®Proceedings of the International Workshop on Indirect
! Methods for Estimating the Hydraulic Properties of Unsaturated Soils. Eds MTh van
! Genuchten et al. pp. 389-403. University of California: Riverside, CA.
! ----------------------------------------------------------
      k_s(7)=729.16*exp(-0.144*clay)                                                         ! cm/d

! ----------------------------------------------------------
! 8) Jabro(1992):
! Jabro, J.D. 1992. Estimation of saturated hydraulic conductivity of soils from particle size
! distribution and bulk density data. Transactions of ASAE 35 (2), 557-560.
!
! 350 soil samples of varying types
! ----------------------------------------------------------
      k_s(8)=24.*10.**(9.56-0.82*log10(silt)-1.09*log10(clay)-4.64*BD)                       ! cm/d

! -----------------------------------------------------------
! 9) Brakensiek et al. (1984):
! Brakensiek, D. L., W. J. Rawls, G. R. Stephenson. 1984. Modifying SCS hydrologic soil groups
! and curve numbers for rangeland soils. ASAE paper no. PNR-84203, St. Joseph, MI.
!
! 1323 soil samples across the US.
! -----------------------------------------------------------
      x = 19.52348*phi-8.96847-0.028212*clay+0.00018107*(sand)**2 &
        - 0.0094125*(clay)**2-8.395215*phi**2+0.077718*phi*sand &
        - 0.00298*(phi*sand)**2-0.019492*(phi*clay)**2+0.0000173*(sand)**2*clay &
        + 0.02733*phi*(clay)**2+0.001434*phi*(sand)**2-0.0000035*(clay)**2*sand

      k_s(9)=24.0*exp(x)                                                                     ! cm/d

! ----------------------------------------------------------
! 10) Rawls et al. (1998):
! W. J. Rawls, D. Gimenez, R. Grossman, 1998: Use of soil texture, bulk density, 
! and slope of thewater retention curve to predict saturated hydraulic conductivity
! Transactions of the ASAE. VOL. 41(4): 983-988
!
! 900 saturated hydraulic conductivities for the soil matrix according to the USDA
! soil texture classes and two porosity classes. The
! Brooks Corey pore size distribution index (lambda) was obtained
! by fitting a log-log plot of water content vs pressure head
! using only the -33 and -1500 kPa water contents.
! ----------------------------------------------------------
      lambda=exp(-0.7842831+0.0177544*sand-1.062498*phi-0.00005304*sand**2 &
            -0.00273493*clay**2+1.11134946*phi**2-0.03088295*sand*phi &
            +0.00026587*(sand*phi)**2-0.00610522*(clay*phi)**2 &
            -0.00000235*sand**2*clay+0.00798746*clay**2*phi-0.0067449*clay*phi**2)
      k_s(10)=4632.*phi_e**(3.-lambda)                                                       ! cm/d

! ---------------------------------------------
! 11) Julia et al. (2004): 
! Julia, M. F., T. E Monreal, A. Sanchez del Corral Jimeneza, and E. Garc Melendez. 2004.
! Constructing a saturated hydraulic conductivity map of Spain using Pedotransfer
! functions and spatial prediction. Geoderma 123, 257-277.
! ---------------------------------------------
      k_s(11)=2.208*exp(0.0491*sand)                                                         ! cm/d

! ----------------------------------------------------------
! 12) Campbell (1985):
! Campbell GS (1985) Soil physics with basic development in soil science, 14. Elsevier, Amsterdam.
!
! d_g=geometric mean particle size (mm),
! sigma_g=geometric standard deviation of the particle size distribution,
! Note that this formula would lead to extremely inaccurate results if the bulk density 
! were low (BD << 1 ).  ! the arithmetic means 1.025mm for sand (2-0.05mm),
! 0.026mm for silt (0.05-0.002mm), 0.001 for clay (<0.002mm)
! ----------------------------------------------------------
      d_g=exp(log(1.025)*sand/100.+log(0.026)*silt/100.+log(0.001)*clay/100.)
      sigma_g=exp(sqrt((log(1.025))**2*sand/100.+(log(0.026))**2*silt/100. &
             +(log(0.001))**2*clay/100.-(log(1.025)*sand/100.+log(0.026)*silt/100. &
             +log(0.001)*clay/100.)**2))
      B=1./sqrt(d_g)+0.2*sigma_g
      k_s(12)=339.0*(1.3/BD)**(1.3*B)*exp(-6.88*clay/100.-3.63*silt/100.-0.025)              ! cm/d

! GROUP II: input soil physical parameters: Texture, Bulk Density and Organic Matter/Carbon
! ---------------------------------------------
! 13) Li et al. (2007): 
! Y. Li, D. Chen, R.E. White, A. Zhu, J. Zhang, 2007: Estimating soil hydraulic 
! properties of Fengqiu County soils in the North China Plain using 
! pedo-transfer functions. Geoderma, 138:261-271.
!
! 36 samples for k_s from Fengqiu county in the North China Plain
! 8.98< sand <93.03;    mean 53.28       s.d. 4.46
! 1.74< silt <79.5      mean 35.86       s.d. 3.83
! 0.54< clay <27.12;    mean 8.86        s.d. 1.01
! 0.12< SOM  <1.54;     mean 0.65        s.d. 0.07
! 1.2 < BD   <1.59      mean 1.42        s.d. 0.01
! ---------------------------------------------
      if(((sand>=8.98) .and. (sand<=93.03)) .AND. & 
         ((silt>=1.74) .and. (silt<=79.5))  .AND. &
         ((clay>=0.54) .and. (clay<=27.12)) .AND. &
         ((SOM>=0.12)  .and. (SOM<=1.54))   .AND. &
         ((BD>=1.2)    .and. (BD<=1.59)))    then
         k_s(13)=exp(13.262-1.914*log(sand)-0.974*log(silt)-0.058*clay &
                -1.709*log(SOM)+2.885*SOM-8.026*log(BD))                                     ! cm/d
      else
         k_s(13)=-1.e36
      endif

! ---------------------------------------------
! 14) Vereecken et al. (1990):
! Vereecken, H., J. Maes, and J. Feyen. 1990. Estimating unsaturated hydraulic
! conductivity from easily measured soil properties. Soil Sci. 149:1-12.
!
! using soils from Belgium.
! ---------------------------------------------
      if(SOM > 0.01)then
      k_s(14)=exp(20.62-0.96*log(clay)-0.66*log(sand)-0.46*log(SOM)-8.43*BD)                 ! cm/d
      else
         k_s(14)=-1.e36
      endif

! ---------------------------------------------
! 15) Wosten et al. (1999):
! Wosten, J.H.M., Lilly, A., Nemes, A., Le Bas, C., 1999. Development and use of a
! database of hydraulic properties of European soils. Geoderma 90,169-185.
!
! HYPRES (EU-wide) soil samples 5521 soil horizons
! where TOPSOIL which is a categorical variable, having a value of 1 
! if the soil sample comes from the topsoil (i.e., A or E horizon, 
! according to the FAO soil classification [FAO, 1990] or 0 if it is from the subsoil).
! topsoil: 2.272 < k_s < 60
! subsoil: 4.0 < k_s < 70.0
! oganic: 8.0
! the properties of EU soil datasets used to develop the PTF:
! 0.49<sand<100,   mean 30.24. s.d. 28.51
! 0   <silt<81.16, mean 41.93, s.d. 20.59
! 0   <clay<80,    mean 27.8,  s.d. 17.44
! 0.9 <BD  <1.9,   mean 1.46,  s.d. 0.19
! 0.01<k_s <2423,  mean 129,   s.d. 293 
! ---------------------------------------------
      if(sand>0.49 .or. silt<81.16 .or. clay<80. .or. (BD>0.9 .and. BD<1.9))then
      k_s(15)=exp(7.75+0.0352*silt+0.93*(TOPSOIL)-0.967*BD**2 &
             - 0.000484*(clay)**2-0.000322*(silt)**2 + 0.001/(silt) &
             - 0.0748/SOM-0.643*log(silt)-0.01398*BD*clay-0.1673*BD*SOM &
             + 0.02986*(TOPSOIL)*clay-0.03305*(TOPSOIL)*silt)                                ! cm/d
      else
         k_s(15)=-1.e36
      endif

! -------------------------------------------
! 16) Merdun (2010) MLR:
! H. Merdun 2010: Alternative Methods in the Development of Pedotransfer Functions
! for Soil Hydraulic Characteristics. Eurasian Soil Science,Vol. 43, No. 1, 62-71.
!
! data from the UNSODA database
! 135 soil samples for the development and 45 for the validation
! 2.0  < sand <96.0   mean 47.2     s.d. 29.7
! 1.0  < silt <90.0   mean 35.2     s.d. 22.7
! 0.1  < clay <63.3   mean 0.176    s.d. 14.5
! 0.01 < SOM  <21.4;  mean 1.673    s.d. 2.478
! 0.59 < BD   <1.76   mean 1.42     s.d. 0.22
! 0.1 < k_s < 3844.8 (cm/d); mean: 334.0 
! -------------------------------------------
      if(((sand>=2.0) .and. (sand<=96.0)) .AND. & 
         ((silt>=1.0) .and. (silt<=90.0)) .AND. &
         ((clay>=0.1) .and. (clay<=63.3)) .AND. &
         ((SOM>=0.01) .and. (SOM<=21.4))  .AND. &
         ((BD>=0.59)  .and. (BD<=1.76)))then
         k_s(16)=9509.-14437.*silt/100.-8169.*BD &
                -860.1*SOM+2332.*(silt/100.)**2+1620.*BD**2+9.113*SOM**2 &
                +7547.*silt/100.*BD+985.1*silt/100.*SOM+381.4*BD*SOM                         ! cm/d

! 17) Merdun (2010) SUR:
         k_s(17)=8777.937-12556.9*silt/100.-7849.05*BD &
                -728.977*SOM+2100.45*(silt/100.)**2+1666.81*BD**2 &
                +6.79971*SOM**2+6597.450*silt/100.*BD &
                +736.149*silt/100.*SOM+371.5434*BD*SOM                                       ! cm/d
      else
         k_s(16:17)=-1.0e36
      endif

! ----------------------------------------------------------
! 18) Aimrun and Amin (2009):
! W. Aimrun, M. S. M. Amin, 2009: Pedo-transfer function for saturated hydraulic conductivity
! of lowland paddy soils. Paddy Water Environ, 7:217-225
!
! Soil samples (n=404 samples) were collected randomly depending on 
! the soil series within the 2300 ha Sawah Sempadan rice cultivation area in Malaysia.
! 0.62  < BD   < 1.91       mean 1.09        s.d. 0.19.
! 19.84 < silt < 57.25      mean 42.67       s.d. 9.27
! 32.22 < clay < 76.80      mean 50.97       s.d. 10.43
! 0.07  < SOM  < 29.35      mean 8.57        s.d. 4.79
! d_g=geometric mean diameter
! arithmetic mean diameter of class: 1.025 mm for sand particle size, 
! 0.026 mm for silt particle size, 0.001 mm for silt particle size.
! the predicted value: 0.126<k_s<2.33; mean: 0.489 cm/d
! ----------------------------------------------------------
      if(((silt>=19.84) .and. (silt<=57.25)) .AND. &
         ((clay>=32.22) .and. (clay<=76.80)) .AND. &
         ((SOM>=0.07)   .and. (SOM <=29.35)) .AND. &
         ((BD>=0.62)    .and. (BD<=1.91   )))then
         d_g=exp(sand/100.*log(1.025)+silt/100.*log(0.026)+clay/100.*log(0.001))
         k_s(18)=100.0*exp(-2.368+3.846*BD+0.091*SOM-6.203*log(BD) &
                -0.343*log(SOM)-2.334*log(clay)-0.411*log(d_g))                              ! cm/d
      else
         k_s(18)=-1.0e36
      endif

! ----------------------------------------------------------
! 19) Saxton and Rawls (2006):
! K. E. Saxton and W. J. Rawls, 2006: Soil Water Characteristic Estimates by Texture and 
! Organic Matter for Hydrologic Solutions. Soil Sci. Soc. Am. J. 70:1569-1578.
!
! A-horizon samples from USDA/NRCS National Soil Characterization Dataset.
! Samples with "extreme" values were omitted from the data. 
! Excluded were those with bulk density < 1.0 and > 1.8 g cm-3,
! OM > 8 % (w) and clay > 60% (w). This reduced the A-horizon data set 
! from 2149 to 1722 samples. Set the value of the residual water content equal to zero
! the fitting samples with 
! 1.0  < BD <1.8 g/cm3
! SOM  < 8  %(w)
! clay < 60 %(w).
! ----------------------------------------------------------
      if((clay<=60.0) .and. (SOM<=8.0) .and. (BD>=1.0 .and. BD<=1.8))then
         x=-0.251*sand/100.+0.195*clay/100.+0.011*SOM &
          +0.006*(sand/100.*SOM)-0.027*(clay/100.*SOM) &
          +0.452*(sand/100.*clay/100.)+0.299
         z33=x+(1.283*x**2-0.374*x-0.015)

         zs33=-0.107+1.636*(0.278*sand/100.+0.034*clay/100.+0.022*SOM &
             -0.018*(sand/100.*SOM)-0.027*(clay/100.*SOM) &
             -0.584*(sand/100.*clay/100.)+0.078)

         zs=z33+zs33-0.097*sand/100.+0.043

         z1500=-0.02+1.14*(0.024*sand/100.+0.487*clay/100.+0.006*SOM &
              +0.005*(sand/100.*SOM)-0.013*(clay/100.*SOM) &
              +0.068*(sand/100.*clay/100.)+0.031)

         lambda=(log(z33)-log(z1500))/(log(1500.)-log(33.))
         k_s(19)=4632.*(zs-z33)**(3.-lambda)                                                    ! cm/d
      else
         k_s(19)=-1.0e36
      endif

! ------------------------------------------
! 20) Weynants et al. (2009):
! M. Weynants, H. Vereecken, and M. Javaux, 2009: Revisiting Vereecken Pedotransfer
! Functions: Introducing a Closed-Form Hydraulic Model. Vadose Zone J. 8:86-95.
!
! datasets: 136 hydraulic conductivity curve and 166 moisture retention curve in the middle of 
! each of the 182 horizons from northern Belgium. Training data value ranges:
! 5.6  < sand < 97.80
! 0    < clay < 54.50
! 0.1  < SOC  < 66 (g kg-1)
! 0.89 < BD   < 1.77 (g cm-3)
! the theta_r could be 0.
! ------------------------------------------
      if(((sand>=5.6) .and. (sand<=97.8)) .AND. &
         ((clay>=0.)  .and. (clay<=54.5)) .AND. &
         ((SOC>=0.01) .and. (SOC<=6.6))   .AND. &
         ((BD>=0.89)  .and. (BD<=1.77)))then 
         k_s(20)=exp(1.9582+0.0308*sand-0.6142*BD-0.1566*SOC/10.)
      else
         k_s(20)=-1.0e36
      endif

      do i = 1, 20
         if(k_s(i) >=900.) k_s(i)=-1.0e36
         if(k_s(i) < 0.1) k_s(i)=-1.0e36
      enddo

END SUBROUTINE ksat 



!======================
SUBROUTINE debar(a,n,x)
!======================
use precision

IMPLICIT NONE
integer n
real(r8) a(n), tmp(n)
real(r8) x
integer i, ii
real(r8), allocatable :: c(:)
real(r8), external :: median

      ii = 0
      do i = 1, n
         if(abs(a(i)) < 1.0e10)then
           ii = ii + 1
           tmp(ii) = a(i)
         endif
      enddo
      allocate (c(ii))
      c = tmp(1:ii)
      x =  median(c,ii)
      deallocate (c)

END SUBROUTINE debar
