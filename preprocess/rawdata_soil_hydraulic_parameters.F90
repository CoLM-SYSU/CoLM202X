SUBROUTINE soil_hydraulic_parameters(BD,SAND,CLAY,SOC,SOILDEPTH,&
           vf_gravels_ss,phi,CampBC_psi_s,CampBC_lambda_s,k_s,&
           VGM_theta_r,VGM_alpha,VGM_n,VGM_L,&
           VGM_theta_r_Rose,VGM_alpha_Rose,VGM_n_Rose,k_s_Rose)

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Calculate soil hydraulic parameters of soil water retension models
!  (Brooks and Corey, 1964 & van Genuchten, 1980) and soil saturated
!  hydraulic conductivity with multiple soil Pedotransfer functions by
!  using the rawdata soil properties.
!
! !REFERENCES:
!  (1) Dai et al.,2013: Development of a China Dataset of Soil
!      Hydraulic Parameters Using Pedotransfer Functions for Land
!      Surface Modeling.  J. of Hydrometeorology, 14: 869-887. DOI:
!      10.1175/JHM-D-12-0149.1
!  (2) Dai et al.,2019: A Global High-Resolution Data Set of Soil
!      Hydraulic and Thermal Properties for Land Surface Modeling. J. of
!      Advances in Modeling Earth Systems, DOI: 10.1029/2019MS001784
!
!  Original author: Yongjiu Dai, Wei Shangguan, 12/2013/
!
! !REVISIONS:
!  06/2018, Yongjiu Dai, Nan Wei and Yonggen Zhang:
!           add more highly cited or newly developed soil Pedotransfer
!           functions.
!
!  01/2019, Nan Wei: add algorithms for fitting soil hydraulic
!           parameters by multiple soil Pedotransfer functions.
!
!  06/2019, Yongjiu Dai and Nan Wei:
!           consider the gravel effects on soil hydraulic parameters
!-----------------------------------------------------------------------
USE MOD_Precision

IMPLICIT NONE

real(r8), intent(in) :: SAND             ! percent sand particle-size distribution (%, weight)
real(r8), intent(in) :: CLAY             ! percent clay particle-size distribution (%, weight)
real(r8), intent(in) :: SOC              ! soil organic carbon concentration (%, weight)
real(r8), intent(in) :: BD               ! bulk density (g cm-3)
real(r8), intent(in) :: SOILDEPTH        ! soil depth (cm)
real(r8), intent(in) :: vf_gravels_ss    ! volumetric fraction of gravels
real(r8), intent(in) :: phi              ! saturated water content (cm3/cm3)
real(r8), intent(in) :: VGM_theta_r_Rose ! residual moisture content by Rosetta H3
real(r8), intent(in) :: VGM_alpha_Rose   ! a parameter corresponding approximately to the inverse of the air-entry value by Rosetta H3
real(r8), intent(in) :: VGM_n_Rose       ! a shape parameter by Rosetta H3
real(r8), intent(in) :: k_s_Rose         ! saturated hydraulic conductivity (cm/day) by Rosetta H3

real(r8), intent(out) :: CampBC_psi_s    ! matric potential at saturation (cm)
real(r8), intent(out) :: CampBC_lambda_s ! pore size distribution index (dimensionless)

real(r8), intent(out) :: k_s             ! saturated hydraulic conductivity (cm/day)

real(r8), intent(out) :: VGM_theta_r     ! residual moisture content
real(r8), intent(out) :: VGM_alpha       ! a parameter corresponding approximately to the inverse of the air-entry value
real(r8), intent(out) :: VGM_n           ! a shape parameter
real(r8), intent(out) :: VGM_L           ! pore-connectivity parameter

real(r8) SOM, TOPSOIL, BD_om, BD_minerals,a,vf_gravels_s

! --------------------------------------------------------------------
      IF(SOILDEPTH < 30.)THEN   ! cm
         TOPSOIL=1.  ! 1 IF A or E horizon
      ELSE
         TOPSOIL=0.  ! 0 IF the subsoil
      ENDIF

      SOM=1.724*SOC
      vf_gravels_s = vf_gravels_ss/100.

! ---------------------------------------------------------------------
! SOIL  hydraulic parameters of the functions of Brooks and Corey (1964)
! simplified by Campbell (1974) with assumption of "zero residual saturation"
! Parameters of the retention curve of FINE EARTH
! ---------------------------------------------------------------------
      CALL CampBC(BD,SAND,CLAY,SOM,SOC,phi,CampBC_psi_s,CampBC_lambda_s)
      IF (vf_gravels_s > 0.2) CampBC_psi_s = 1.25 * CampBC_psi_s * (1 - vf_gravels_s)

! ---------------------------------------------------------------------
!*SOIL hydraulic parameters in the van Genuchten and Mualem functions
! ---------------------------------------------------------------------
      CALL VGM(BD,SAND,CLAY,SOM,SOC,TOPSOIL,phi,VGM_theta_r,VGM_alpha,VGM_n,VGM_L,&
               VGM_theta_r_Rose,VGM_alpha_Rose,VGM_n_Rose)

! ---------------------------------------------------------------------
! Saturated hydraulic conductivity of FINE EARTH
! ---------------------------------------------------------------------
      CALL ksat(BD,SOM,SOC,SAND,CLAY,TOPSOIL,phi,CampBC_psi_s,CampBC_lambda_s,k_s,k_s_Rose)

! ---------------------------------------------------------------------
! Hydraulic properties of the mixtures of FINE EARTH AND ROCK FRAGMENTS
! ---------------------------------------------------------------------
      a   = k_s
      k_s = a * (1.-1.1*vf_gravels_s)                     ! (Bouwer and Rice, 1984)
      IF(k_s <= 0.0) k_s = a * (1.-vf_gravels_s)
!      k_s = k_s * (1.-vf_gravels_s)/(1.0+vf_gravels_s)   ! (Peck and Watson, 1979)
!***  k_s = k_s * (1.-vf_gravels_s)                       ! (Brakensiek et al., 1986; Bagarello and Iovino, 2007)

! Soil water retention curves for stony soils were estimated under an assumption of zero retention capacity of the rock fragments


END SUBROUTINE soil_hydraulic_parameters


SUBROUTINE CampBC(BD,SAND,CLAY,SOM,SOC,phi,psi_s,lambda_s)
!-----------------------------------------------------------------------
! DESCRIPTION:
! SOIL water retention curve of Brooks and Corey (1964) simplified by Campbell (1974)
! with assumption of "zero residual saturation" parameters of the soil water retention curve of FINE EARTH
!
! REFERENCES:
! (1) Dai et al.,2013: Development of a China Dataset of Soil
!     Hydraulic Parameters Using Pedotransfer Functions for Land Surface Modeling.
!     J. of Hydrometeorology, 14: 869-887. DOI: 10.1175/JHM-D-12-0149.1
! (2) Dai et al.,2019: A Global High-Resolution Data Set of Soil Hydraulic and Thermal Properties
!     for Land Surface Modeling. J. of Advances in Modeling Earth Systems, DOI: 10.1029/2019MS001784
!
! Original author: Yongjiu Dai, Wei Shangguan, 12/2013/
!
! Revisions:
! Yongjiu Dai, Nan Wei and Yonggen Zhang,
!          06/2018: add more highly cited or newly developed soil Pedotransfer functions.
! Nan Wei, 01/2019: add algorithms for fitting soil hydraulic parameters by multiple soil Pedotransfer functions.
! ----------------------------------------------------
USE MOD_Precision
USE MOD_Utils

IMPLICIT NONE
real(r8), intent(in) :: BD   ! bulk density (g cm-3)
real(r8), intent(in) :: SAND ! percent sand particle-size distribution (%, weight)
real(r8), intent(in) :: CLAY ! percent clay particle-size distribution (%, weight)
real(r8), intent(in) :: SOM  ! percent organic matter (%, weight)
real(r8), intent(in) :: SOC  ! percent organic carbon (%, weight)
real(r8), intent(in) :: phi  ! saturated water content (cm3/cm3)

real(r8), intent(out) :: psi_s  ! MATRIC POTENTIAL AT SATURATION ! cm
real(r8), intent(out) :: lambda_s !

integer, parameter :: nc = 8  !  the number of PTFs in estimating SW retention parameters in the Campbell model
logical c(12) ! indicate whether a soil is in an class
! soil_classes=c('Sa','LoSa','SaLo','Lo','SaClLo','SaCl','ClLo','SiLo','Si','SiClLo','SiCl','Cl')

real(r8) CH_b(12),CH_ths(12),CH_psi_s(12)
data CH_b    /4.05 ,4.38,4.9  ,5.39 ,7.12,10.4 ,8.52  ,4.9  ,5.3  ,7.75 ,10.4 ,11.4/
data CH_ths  /0.395,0.41,0.435,0.451,0.42,0.426,0.476 ,0.485,0.485,0.477,0.492,0.482/
data CH_psi_s/12.10,9.0 ,21.8 ,47.80,29.9,15.3 ,63.0  ,78.6 ,78.6 ,35.6 ,49.0 ,40.5/  ! unit is cm

real(r8) Cosby0_b(12),Cosby_psi_s(12),Cosby_ths(12)
data Cosby0_b   /2.790,4.260,4.740,5.250,6.660,10.73,8.170,5.330,5.330,8.720,10.39,11.55/
data Cosby_psi_s/0.069,0.036,0.141,0.355,0.135,0.098,0.263,0.759,0.759,0.617,0.324,0.468/
data Cosby_ths  /0.339,0.421,0.434,0.439,0.404,0.406,0.465,0.476,0.476,0.464,0.468,0.468/

integer i, itype
real(r8) SILT
real(r8) d_g, sigma_g, h_es, h_b
real(r8) a, b, z, z33, zs33, zs, z1500, str
real(r8) psi(nc), lambda(nc)

! local variables for estimating the optimal parameters using the fitting methods
integer ,parameter :: npoint       = 24
real(r8),parameter :: xdat(npoint) = (/1.,5.,10.,20.,30.,40.,50.,60.,70.,90.,110.,130.,150.,170.,210.,300.,&
                                   345.,690.,1020.,5100.,15300.,20000.,100000.,1000000./)
                                    !  the points of soil pressure heads used for fitting SW retention curves
real(r8) ydatc(nc,npoint)           !  the estimated SW at the soil pressure heads using Campbell parameters

integer              :: m           ! the number of valid functions (PTFs)
integer,parameter    :: n = 2       ! the number of parameters to be fitted (psi and lambda)
integer              :: ldfjac,info,ipvt(n),maxfev,mode,nfev,njev,nprint
real(r8)             :: x(n),diag(n),factor,ftol,gtol,xtol,qtf(n)
real(r8),allocatable :: fjac(:,:),fvec(:)
integer isiter                      !  flags to tell whether the iteration is completed, 1=Yes, 0=No
external SW_CB_dist

! h_b= bubbling capillary pressure in centimeters
! phi= pore size index

      SILT=max(1.0,100.0-SAND-CLAY)

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
      psi(1)=-10.0**(1.88-0.013*SAND)                                                      ! cm
      lambda(1)=1./(2.91+0.159*CLAY)
      ydatc(1,:) = (-1.0*xdat/psi(1))**(-1.0*lambda(1)) * phi

! 2) Cosby et al. (1984), multi-variate regression:
! 1448 samples from 35 localities in 23 states in the United States.
! set the value of the residual water content equal to zero
      psi(2)=-10**(1.54-0.0095*SAND+0.0063*SILT)                                           ! cm
      lambda(2)=1.0/(3.10+0.157*CLAY-0.003*SAND)
      ydatc(2,:) = (-1.0*xdat/psi(2))**(-1.0*lambda(2)) * phi

! ------------------------------------------------------------
! 3) Campbell and Shiozawa (1992):
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
      d_g=exp(log(1.025)*SAND/100.+log(0.026)*SILT/100.+log(0.001)*CLAY/100.)
      sigma_g=exp(sqrt((log(1.025))**2*SAND/100.+(log(0.026))**2*SILT/100. &
             +(log(0.001))**2*CLAY/100.-(log(1.025)*SAND/100.+log(0.026)*SILT/100. &
             +log(0.001)*CLAY/100.)**2))

      lambda(3)=1.0/(1./sqrt(d_g)+0.2*sigma_g)
      h_es=0.05/sqrt(d_g)
      h_b=h_es*(BD/1.3)**(0.67/lambda(3))
      psi(3)=-100.*h_b                                                                     ! cm
      ydatc(3,:) = (-1.0*xdat/psi(3))**(-1.0*lambda(3)) * phi

! ---------------------------------------------------------
! 4) Rawls and Brakenssiek (1989):
! Rawls, W. J., & Brakensiek, D. L. (1989). Estimation of soil water retention and hydraulic properties.
! In H. J. Morel‐Seytoux (Ed.), Unsaturated flow in hydrologic modeling, NATO ASI Series
! (Series C: Mathematical and Physical Sciences) (Vol. 275, pp. 275-300).
!
! PTFs were developed from about 18,000 samples from the national soil characterization
! database that the clay activity ratio along with texture to characterized
! the effect of soil minerals on soil hydraulic properties.
! the value of the residual water content NOT equal zero
! the equation for theta_r was given in the subroutine thetar
! ------------------------------------------------------------
      h_b=exp(5.3396738+0.1845038*CLAY-2.48394546*phi-0.00213853*CLAY**2 &
         -0.04356349*SAND*phi-0.61745089*CLAY*phi+0.00143598*(SAND*phi)**2 &
         -0.00855375*(CLAY*phi)**2-0.00001282*SAND**2*CLAY+0.00895359*CLAY**2*phi &
         -0.00072472*(SAND)**2*phi+0.0000054*CLAY**2*SAND+0.5002806*CLAY*phi**2)             ! cm
      psi(4)=-h_b                                                                          ! cm

      lambda(4)=exp(-0.7842831+0.0177544*SAND-1.062498*phi-0.00005304*SAND**2 &
               -0.00273493*CLAY**2+1.11134946*phi**2-0.03088295*SAND*phi &
               +0.00026587*(SAND*phi)**2-0.00610522*(CLAY*phi)**2 &
               -0.00000235*SAND**2*CLAY+0.00798746*CLAY**2*phi-0.0067449*CLAY*phi**2)
      ydatc(4,:) = (-1.0*xdat/psi(4))**(-1.0*lambda(4)) * phi

! ---------------------------------------------------------
! GROUP II: input soil physical parameters: Texture, Bulk Density and Organic Matter/Carbon
! -------------------------------------------------------------
! 5) Mayr and Jarvis (1999):
! T. Mayr, N.J. Jarvis, 1999: Pedotransfer functions to estimate soil water
! retention parameters for a modified Brooks-Corey type model. Geoderma 91:1-9.
!
! 286 soil horizons of the soil physical properties database of England and Wales
! set the value of the residual water content equal to zero
! these function should NOT be applied to organic soils (organic carbon content >5%)
! and/or soils of small dry bulk density (<0.9 g/cm3)
! -------------------------------------------------------------
      IF(SOC < 5. .or. BD > 0.9)THEN
      psi(5)=-10.0**(-4.9840297533+0.0509226283*SAND+0.1575152771*SILT &
              +0.1240901644*BD-0.1640033143*SOC &
              -0.0021767278*SILT**2+0.00001438224*SILT**3 &
              +0.0008040715*CLAY**2+0.0044067117*SOC**2)                                     ! cm

      lambda(5)=10.0**(-0.8466880654-0.0046806123*SAND+0.0092463819*SILT &
               -0.4542769707*BD-0.0497915563*SOC+3.294687e-04*SAND**2-1.689056e-06*SAND**3 &
               +0.0011225373*SOC**2)
      ydatc(5,:) = (-1.0*xdat/psi(5))**(-1.0*lambda(5)) * phi
      ELSE
         lambda(5)=-1.0e36
         psi(5)=-1.0e36
         ydatc(5,:) = -1.0e36
      ENDIF

! ----------------------------------------------------------
! 6) Williams, J., P. Ross, and K. Bristow. 1992. Prediction of the Campbell water retention
! function from texture, structure, and organic matter.p. 427-442. In M.Th. van Genuchten et al (ed.)
! Proc. Int.Workshop on Indirect methods for Estimating the Hydraulic Properties of
! Unsaturated Soils. University of California, Riverside.
! ----------------------------------------------------------
      A=2.57+0.238*Alog(CLAY)-0.000192*SAND*SAND-0.0926*Alog(SOM)+0.0412*SOM
      B=-0.403+0.0871*Alog(CLAY)-0.00077*SAND
      lambda(6)=-B
      psi(6)=-1000.0*exp((alog(100.0*phi)-A)/B)  ! cm
      ydatc(6,:) = (-1.0*xdat/psi(6))**(-1.0*lambda(6)) * phi

! ----------------------------------------------------------
! 7) Clapp and Hornberger (1978) models, from Table 2
! 8) Cosby0 et al ., (1984)  Table 3: lookup table
! ----------------------------------------------------------
      CALL USDA_soil_classes(SILT,CLAY,c)

      IF(c(1))  itype = 12 ! ic = 'Cl'      clay
      IF(c(2))  itype = 11 ! ic = 'SiCl'    silty clay
      IF(c(3))  itype = 6  ! ic = 'SaCl'    sandy clay
      IF(c(4))  itype = 7  ! ic = 'ClLo'    clay loam
      IF(c(5))  itype = 10 ! ic = 'SiClLo'  silty clay loam
      IF(c(6))  itype = 5  ! ic = 'SaClLo'  sandy clay loam
      IF(c(7))  itype = 4  ! ic = 'Lo'      loam
      IF(c(8))  itype = 8  ! ic = 'SiLo'    silty loam
      IF(c(9))  itype = 3  ! ic = 'SaLo'    sandy loam
      IF(c(10)) itype = 9  ! ic = 'Si'      silt
      IF(c(11)) itype = 2  ! ic = 'LoSa'    loamy sand
      IF(c(12)) itype = 1  ! ic = 'Sa'      sand

      lambda(7)=1.0/CH_b(itype)
      psi(7)=-CH_psi_s(itype)  ! cm
      ydatc(7,:) = (-1.0*xdat/psi(7))**(-1.0*lambda(7)) * phi

      lambda(8)=1.0/Cosby0_b(itype)
      psi(8)=-100.0*Cosby_psi_s(itype)  ! m to cm
      ydatc(8,:) = (-1.0*xdat/psi(8))**(-1.0*lambda(8)) * phi

! ----------------------------------------------------------
      DO i = 1, nc
         IF(lambda(i) > 1. .or. lambda(i) <= 0.) THEN
            lambda(i)=-1.0e36
            ydatc(i,:)=-1.0e36
         ENDIF
         IF(psi(i) < -300. .or. psi(i) >= 0.) THEN
            psi(i)=-1.0e36
            ydatc(i,:)=-1.0e36
         ENDIF
      ENDDO

      m = 0
      DO i = 1, nc
         IF(abs(ydatc(i,1)) < 1.0e10) THEN
            m = m+1
            ydatc(m,:) = ydatc(i,:)
         ENDIF
      ENDDO

      ldfjac = npoint
      allocate(fjac(npoint,n))
      allocate(fvec(npoint))

! ------------------------------------------
! THE MEDIAN VALUE as the initial value of fitting
! ------------------------------------------
      CALL debar(psi,nc,psi_s)
      CALL debar(lambda,nc,lambda_s)

      x(1) = psi_s
      x(2) = lambda_s
      factor = 0.1
      maxfev = 100 * ( n + 1 )
      ftol = 1.0e-5
      xtol = 1.0e-4
      gtol = 0.0
      mode = 1
      nprint = 0
      isiter = 1

      CALL lmder ( SW_CB_dist, npoint, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
                   diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf,&
                   xdat, npoint, ydatc, m, phi, isiter)

      IF( x(1) >= -300. .and. x(1) < 0.0 .and. x(2) > 0.0 .and. x(2) <= 1.0 .and. isiter == 1)THEN
         psi_s    = x(1)
         lambda_s = x(2)
      ENDIF

      deallocate(fjac)
      deallocate(fvec)

END SUBROUTINE CampBC


SUBROUTINE VGM(BD,sand,clay,SOM,SOC,TOPSOIL,phi,theta_r_l,alpha_l,n_l,L_l,&
               VGM_theta_r_Rose,VGM_alpha_Rose,VGM_n_Rose)
!-----------------------------------------------------------------------
! DESCRIPTION:
! SOIL water retention curve of van Genuchten-Mualem (van Genuchten 1980).
!
! Note:
! The pore-connectivity parameter L in hydraulic conductivity was estimated by Mualem [1976a]
! to be 0.5 as an average for some 45 soils. Wosten and van Genuchten [1988],
! in an analysis of some 200 soil hydraulic data sets, found 4 to vary between -16 and more
! than 2. Fixing L at 0.5 in their study produced acceptable results for most coarse-textured soils,
! but not for many medium- and fine-textured soils. These results suggest that keeping L variable in the
! parameters optimization process will likely improve the analysis of individual soil hydraulic
! conductivity data sets, provided enough measured data points are available for the estimation process.
!
! REFERENCES:
! (1) Dai et al.,2013: Development of a China Dataset of Soil
!     Hydraulic Parameters Using Pedotransfer Functions for Land Surface Modeling.
!     J. of Hydrometeorology, 14: 869-887. DOI: 10.1175/JHM-D-12-0149.1
! (2) Dai et al.,2019: A Global High-Resolution Data Set of Soil Hydraulic and Thermal Properties
!     for Land Surface Modeling. J. of Advances in Modeling Earth Systems, DOI: 10.1029/2019MS001784
!
! Original author: Yongjiu Dai, Wei Shangguan, 12/2013/
!
! Revisions:
! Yongjiu Dai, Nan Wei and Yonggen Zhang,
!          06/2018: add more highly cited or newly developed soil Pedotransfer functions.
! Nan Wei, 01/2019: add algorithms for fitting soil hydraulic parameters by multiple soil Pedotransfer functions.
! ----------------------------------------------------
USE MOD_Precision
USE MOD_Utils

IMPLICIT NONE

real(r8), intent(in) :: BD   ! bulk density (g cm-3)
real(r8), intent(in) :: sand ! percent sand particle-size distribution (%)
real(r8), intent(in) :: clay ! percent clay particle-size distribution (%)
real(r8), intent(in) :: SOM  ! percent organic matter (%)
real(r8), intent(in) :: SOC  ! percent organic carbon (%)
real(r8), intent(in) :: TOPSOIL ! 1 for A-E soil horizon, 0 the subsoil
real(r8), intent(in) :: phi     ! saturated water content (cm3/cm3)
real(r8), intent(in) :: VGM_theta_r_Rose
real(r8), intent(in) :: VGM_alpha_Rose
real(r8), intent(in) :: VGM_n_Rose

real(r8), intent(out) :: theta_r_l
real(r8), intent(out) :: alpha_l
real(r8), intent(out) :: n_l
real(r8), intent(out) :: L_l

integer, parameter    :: nv = 11 !  the number of PTFs in estimating SW retention parameters in the VG model
real(r8) theta_r(nv),alpha(nv),n(nv),L(nv)

integer i,j,iclass,itype
real(r8) h_b ! bubbling capillary pressure in centimeters
real(r8) silt,x1,y
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

real PG(6)
real(r8) Theta(14),PRawls02(11),PRawls03(10),PRajkai(8),PGupta(13),PTomasella(9)
data PTomasella/0.0,10.0,30.0,60.0,100.0,330.0,1000.0,5000.0,15000.0/
data PGupta/0.0,40.0,70.0,100.0,200.0,330.0,600.0,1000.0,2000.0,4000.0,7000.0,10000.0,15000.0/
data PRajkai/0.0,3.0,10.0,32.0,500.0,2512.0,15849.0,1259000.0/
data PRawls02/0.0,100.0,200.0,330.0,600.0,1000.0,2000.0,4000.0,7000.0,10000.0,15000.0/
data PRawls03/0.0,200.0,300.0,600.0,1000.0,2000.0,4000.0,7000.0,10000.0,15000.0/

logical c(12) ! indicate whether a soil is in an class
real(r8) CP_thr(12),CP_alpha(12),CP_n(12)
data CP_thr  /0.045,0.057,0.065,0.078,0.095,0.1  ,0.095,0.067,0.034,0.089,0.07 ,0.068/
data CP_alpha/0.145,0.124,0.075,0.036,0.019,0.027,0.019,0.02 ,0.016,0.01 ,0.005,0.008/
data CP_n    /2.68 ,2.28 ,1.89 ,1.56 ,1.31 ,1.23 ,1.31 ,1.41 ,1.37 ,1.23 ,1.09 ,1.09/

real(r8) RoH1_thr(12),RoH1_alpha(12),RoH1_n(12)
data RoH1_thr  /0.055,0.058,0.061,0.09, 0.093,0.147,0.107,0.083,0.065,0.12, 0.124,0.131/
data RoH1_alpha/0.033,0.025,0.016,0.006,0.012,0.025,0.01, 0.003,0.006,0.006,0.01, 0.009/
data RoH1_n    /2.895,1.697,1.457,1.421,1.305,1.237,1.391,1.552,1.577,1.434,1.273,1.255/

! local variables for estimating the optimal parameters using the fitting methods
integer ,parameter :: npoint       = 24
real(r8),parameter :: xdat(npoint) = (/1.,5.,10.,20.,30.,40.,50.,60.,70.,90.,110.,130.,150.,170.,210.,300.,&
                                   345.,690.,1020.,5100.,15300.,20000.,100000.,1000000./)
                                    !  the points of soil pressure heads used for fitting SW retention curves
real(r8) ydatv(nv,npoint)           !  the estimated SW at the soil pressure heads using VG parameters
integer isiter                      !  flags to tell whether the iteration is completed, 1=Yes, 0=No
integer              :: m           !  the number of valid functions (PTFs)
integer,parameter    :: n1 = 3      !  the number of parameters to be fitted (namely theta_r, alpha and n)
integer              :: ldfjac,info,ipvt(n1),maxfev,mode,nfev,njev,nprint
real(r8)             :: x(n1),diag(n1),factor,ftol,gtol,xtol,qtf(n1)
real(r8),allocatable :: fjac(:,:),fvec(:)
external SW_VG_dist

silt=100.-sand-clay

! GROUP I: input soil physical parameters: Texture, Bulk Density, Porosity
! -------------------------------------------
! 1) Rawls and Brakenssiek (1989):
! Rawls, W. J., & Brakensiek, D. L. (1989). Estimation of soil water retention and hydraulic properties.
! In H. J. Morel‐Seytoux (Ed.), Unsaturated flow in hydrologic modeling, NATO ASI Series
! (Series C: Mathematical and Physical Sciences) (Vol. 275, pp. 275-300).
!
! 12,000 samples from the National Soil Characterization database
! the PTF is applicable for soil:
! 5 < sand < 70, 5 < clay < 60.
! -------------------------------------------
theta_r(1)=-0.0182482+0.00087269*sand+0.00513488*clay+0.02939286*phi &
-0.00015395*(clay)**2-0.0010827*sand*phi-0.00018233*(clay*phi)**2 &
+0.00030703*(clay)**2*phi-0.0023584*(clay)*(phi)**2

h_b=exp(5.3396738+0.1845038*clay-2.48394546*phi-0.00213853*clay**2 &
-0.04356349*sand*phi-0.61745089*clay*phi+0.00143598*(sand*phi)**2 &
-0.00855375*(clay*phi)**2-0.00001282*sand**2*clay+0.00895359*clay**2*phi &
-0.00072472*(sand)**2*phi+0.0000054*clay**2*sand+0.5002806*clay*phi**2)                ! cm

x1=exp(-0.7842831+0.0177544*sand-1.062498*phi-0.00005304*sand**2 &
-0.00273493*clay**2+1.11134946*phi**2-0.03088295*sand*phi &
+0.00026587*(sand*phi)**2-0.00610522*(clay*phi)**2 &
-0.00000235*sand**2*clay+0.00798746*clay**2*phi-0.0067449*clay*phi**2)

alpha(1)=1.0/h_b                                                                       ! 1/cm
n(1)=1.0+x1
L(1)=0.5
ydatv(1,:) = theta_r(1)+(phi - theta_r(1))*(1+(alpha(1)*xdat)**n(1))**(1.0/n(1)-1)

! -------------------------------------------
! 2) Wosten et al. (1999):
! Wosten, J.H.M., A. Lilly, A. Nemes, and C. Le Bas. 1999. Development and
! use of a database of hydraulic properties of European soils. Geoderma 90:169-185.
!
! Wosten et al. (1999) analyzed the all-Europe database and derived
! the following PTFs to estimate van Genuchten parameters:
! WHERE topsoil is an ordinal variable having the value of 1 or
! of 0. %silt is the percentage of soil (2mm-50mm), BD= bulk density (g/cm3 = Mg/m3)
! the PTFs developd by using the European HYdraulic PRoperties of European Soils (HYPRES)
! soil database (5521 horizons)
! topsoil: 0.0083 < alpha < 0.0383; 1.1033 < n < 1.3774
! subsoil: 0.0082 < alpha < 0.0430; 1.0730 < n < 1.5206
! -------------------------------------------
theta_r(2)=0.0

alpha(2)=exp(-14.96+0.03135*clay+0.0351*silt+0.646*SOM+15.29*BD-0.192*TOPSOIL &
-4.671*BD**2-0.000781*clay**2-0.00687*SOM**2+0.0449/SOM+0.0663*log(silt) &
+0.1482*log(SOM)-0.04546*BD*silt-0.4852*BD*SOM+0.00673*TOPSOIL*clay)                   ! 1/cm

n(2)=1.0+exp(-25.23-0.02195*clay+0.0074*silt-0.1940*SOM+45.5*BD-7.24*BD**2 &
+0.0003658*(clay)**2+0.002885*SOM**2-12.81/BD-0.1524/silt-0.01958/SOM &
-0.2876*log(silt)-0.0709*log(SOM)-44.6*log(BD)-0.02264*BD*clay &
+0.0896*BD*SOM+0.00718*TOPSOIL*clay)

x1=0.0202+0.0006193*clay**2-0.001136*SOM**2-0.2316*log(SOM) &
-0.03544*BD*clay+0.00283*BD*silt+0.0488*BD*SOM
L(2)=10.*(exp(x1)-1.0)/(exp(x1)+1.0)
ydatv(2,:) = theta_r(2)+(phi - theta_r(2))*(1+(alpha(2)*xdat)**n(2))**(1.0/n(2)-1)

! ------------------------------------------
! 3) Wosten et al. (1999):
! FAO class
! ------------------------------------------
IF(clay.gt.60.) THEN
   iclass=5
ELSEIF(clay.gt.35.) THEN
   iclass=4
ELSEIF(sand.lt.15.) THEN
   iclass=3
ELSEIF(sand.gt.65.0.and.clay.lt.18.0) THEN
   iclass=1
ELSE
   iclass=2
ENDIF

j=iclass+nint(TOPSOIL)*5

theta_r(3)=params(1,j)
alpha(3)=params(3,j)                                                                   ! 1/cm
n(3)=params(4,j)
L(3)=params(6,j)
ydatv(3,:) = theta_r(3)+(phi - theta_r(3))*(1+(alpha(3)*xdat)**n(3))**(1.0/n(3)-1)

! ------------------------------------------
! 4) Weynants et al. (2009):
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
theta_r(4)=theta_r(3)
alpha(4)=exp(-4.3003-0.0097*clay+0.0138*sand-0.0992*SOC/10.)
n(4)=1.0+exp(-1.0846-0.0236*clay-0.0085*sand+0.0001*(sand)**2)
L(4)=-1.8642-0.1317*clay+0.0067*sand
ydatv(4,:) = theta_r(4)+(phi - theta_r(4))*(1+(alpha(4)*xdat)**n(4))**(1.0/n(4)-1)

! ------------------------------------------
! 5) Rosetta H1
! ------------------------------------------
      CALL USDA_soil_classes(SILT,CLAY,c)

      IF(c(1))  itype = 12 ! ic = 'Cl'      clay
      IF(c(2))  itype = 11 ! ic = 'SiCl'    silty clay
      IF(c(3))  itype = 6  ! ic = 'SaCl'    sandy clay
      IF(c(4))  itype = 7  ! ic = 'ClLo'    clay loam
      IF(c(5))  itype = 10 ! ic = 'SiClLo'  silty clay loam
      IF(c(6))  itype = 5  ! ic = 'SaClLo'  sandy clay loam
      IF(c(7))  itype = 4  ! ic = 'Lo'      loam
      IF(c(8))  itype = 8  ! ic = 'SiLo'    silty loam
      IF(c(9))  itype = 3  ! ic = 'SaLo'    sandy loam
      IF(c(10)) itype = 9  ! ic = 'Si'      silt
      IF(c(11)) itype = 2  ! ic = 'LoSa'    loamy sand
      IF(c(12)) itype = 1  ! ic = 'Sa'      sand

      theta_r(5)=RoH1_thr(itype)
      alpha(5)=RoH1_alpha(itype)
      n(5)    =RoH1_n(itype)
      L(5)    =0.5
      ydatv(5,:) = theta_r(5)+(phi - theta_r(5))*(1+(alpha(5)*xdat)**n(5))**(1.0/n(5)-1)

! ------------------------------------------
! 6) Rosetta H3 provided by yonggen Zhang
! ------------------------------------------
      theta_r(6)=VGM_theta_r_Rose
      alpha(6)  =VGM_alpha_Rose
      n(6)      =VGM_n_Rose
      L(6)      =0.5
      IF (theta_r(6) .lt. 0.0 .or. theta_r(6)  .gt. phi .or. alpha(6) .lt. 0.00001 .or. &
          alpha(6)   .gt. 1.  .or. n(6)        .lt. 1.1 .or. n(6)     .gt. 10.    )THEN
          ydatv(6,:) = -1.0e36
      ELSE
          ydatv(6,:) = theta_r(6)+(phi - theta_r(6))*(1+(alpha(6)*xdat)**n(6))**(1.0/n(6)-1)
      ENDIF

! ------------------------------------------
! 7) Gupta, S.C., and W.E. Larson. 1979. Estimating soil water retention characteristics from
!    particle-size distribution, organic matter percent, and bulk density. Water Resour. Res. 15:1633-1635.
! ------------------------------------------
      Call GuLar(CLAY,SILT,SAND,SOC,BD,Theta)
      Theta(1)=phi
      CALL VGpar(Pgupta,Theta,13,PG,4,50)
      theta_r(7)=PG(1)
      alpha(7)=PG(3)
      n(7)=PG(4)
      L(7)=0.5
      ydatv(7,:) = theta_r(7)+(phi - theta_r(7))*(1+(alpha(7)*xdat)**n(7))**(1.0/n(7)-1)

! ------------------------------------------
! 8) Carsel and Parrish 1988. Developing joint probability distributions of soil water
!    retention characteristics. Water Resources Research, 24(5), 755-769.
! ------------------------------------------
      CALL USDA_soil_classes(SILT,CLAY,c)

      IF(c(1))  itype = 12 ! ic = 'Cl'      clay
      IF(c(2))  itype = 11 ! ic = 'SiCl'    silty clay
      IF(c(3))  itype = 6  ! ic = 'SaCl'    sandy clay
      IF(c(4))  itype = 7  ! ic = 'ClLo'    clay loam
      IF(c(5))  itype = 10 ! ic = 'SiClLo'  silty clay loam
      IF(c(6))  itype = 5  ! ic = 'SaClLo'  sandy clay loam
      IF(c(7))  itype = 4  ! ic = 'Lo'      loam
      IF(c(8))  itype = 8  ! ic = 'SiLo'    silty loam
      IF(c(9))  itype = 3  ! ic = 'SaLo'    sandy loam
      IF(c(10)) itype = 9  ! ic = 'Si'      silt
      IF(c(11)) itype = 2  ! ic = 'LoSa'    loamy sand
      IF(c(12)) itype = 1  ! ic = 'Sa'      sand

      theta_r(8)=CP_thr(itype)
      alpha(8)=CP_alpha(itype)
      n(8)    =CP_n(itype)
      L(8)    =0.5
      ydatv(8,:) = theta_r(8)+(phi - theta_r(8))*(1+(alpha(8)*xdat)**n(8))**(1.0/n(8)-1)

! ------------------------------------------
! 9) Rawls, W.J., D.L. Brakensiek, and K.E. Saxton. 1982. Estimation of soil water properties.
!    Trans. ASAE 25:1316-1320.
! ------------------------------------------
      Call Rawls82(CLAY,SILT,SAND,SOC,Theta)
      Theta(1)=phi
      CALL VGpar(PRawls02,Theta,11,PG,4,50)
      theta_r(9)=PG(1)
      alpha(9)=PG(3)
      n(9)=PG(4)
      L(9)=0.5
      ydatv(9,:) = theta_r(9)+(phi - theta_r(9))*(1+(alpha(9)*xdat)**n(9))**(1.0/n(9)-1)

! ------------------------------------------
! 10  Rawls, W.J., D.L. Brakensiek, and B. Soni. 1983. Agricultural management effects on soil
!     water processes. Part I. Soil water retention and Green-Ampt parameters. Trans. ASAE 26:1747-1752.
! ------------------------------------------
      CALL Rawls83(CLAY,SILT,SAND,BD,SOC,Theta)
      Theta(1)=phi
      CALL VGpar(PRawls03,Theta,10,PG,4,50)
      theta_r(10)=PG(1)
      alpha(10)=PG(3)
      n(10)=PG(4)
      L(10)=0.5
      ydatv(10,:) = theta_r(10)+(phi - theta_r(10))*(1+(alpha(10)*xdat)**n(10))**(1.0/n(10)-1)

! ------------------------------------------
! 11) Tomasella, J., and M.G. Hodnett. 1998. Estimating soil water retention characteristics from
!     limited data in Brazilian Amazonia. Soil Sci. 163:190-202.
! ------------------------------------------
      Call Tomasella(CLAY,SILT,SOC,Theta)
      Theta(1)=phi
      CALL VGpar(PTomasella,Theta,9,PG,4,50)
      theta_r(11)=PG(1)
      alpha(11)=PG(3)
      n(11)=PG(4)
      L(11)=0.5
      ydatv(11,:) = theta_r(11)+(phi - theta_r(11))*(1+(alpha(11)*xdat)**n(11))**(1.0/n(11)-1)

! ----------------------------------------------------------
      DO i = 1, nv
         IF(theta_r(i) > phi .or. theta_r(i) < 0.0) THEN
            theta_r(i)=-1.0e36
            ydatv(i,:)=-1.0e36
         ENDIF
         IF(alpha(i) < 1.0e-5 .or. alpha(i) > 1.0) THEN
            alpha(i)  =-1.0e36
            ydatv(i,:)=-1.0e36
         ENDIF
         IF(n(i) < 1.1 .or. n(i) > 10.0) THEN
            n(i)      =-1.0e36
            ydatv(i,:)=-1.0e36
         ENDIF
      ENDDO

      m = 0
      DO i = 1, nv
         IF(abs(ydatv(i,1)) < 1.0e10) THEN
            m = m+1
            ydatv(m,:) = ydatv(i,:)
         ENDIF
      ENDDO

      ldfjac = npoint
      allocate(fjac(npoint,n1))
      allocate(fvec(npoint))

! ------------------------------------------
! THE MEDIAN VALUE as the initial value of fitting
! ------------------------------------------
      CALL debar(theta_r,nv,theta_r_l)
      CALL debar(alpha,  nv,alpha_l)
      CALL debar(n,      nv,n_l)
      CALL debar(L,      nv,L_l)

      x(1) = theta_r_l
      x(2) = alpha_l
      x(3) = n_l
      factor = 0.1
      maxfev = 100 * ( n1 + 1 )
      ftol = 1.0e-5
      xtol = 1.0e-4
      gtol = 0.0
      mode = 1
      nprint = 0
      isiter = 1

      CALL lmder ( SW_VG_dist, npoint, n1, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
                   diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf, &
                   xdat, npoint, ydatv, m, phi, isiter )

      IF ( x(1) >= 0.0 .and. x(1) <= phi .and. x(2) >= 1.0e-5 .and. x(2) <= 1.0 .and. &
           x(3) >= 1.1 .and. x(3) <= 10.0 .and. isiter == 1) THEN
          theta_r_l = x(1)
          alpha_l   = x(2)
          n_l       = x(3)
      ENDIF

      deallocate(fjac)
      deallocate(fvec)

END SUBROUTINE VGM


SUBROUTINE ksat(BD,SOM,SOC,SAND,CLAY,TOPSOIL,phi,psi,lambda,k_s,k_s_Rose)
!-----------------------------------------------------------------------
! DESCRIPTION:
! SOIL saturated hydraulic conductivity.
!
! REFERENCES:
! (1) Dai et al.,2013: Development of a China Dataset of Soil
!     Hydraulic Parameters Using Pedotransfer Functions for Land Surface Modeling.
!     J. of Hydrometeorology, 14: 869-887. DOI: 10.1175/JHM-D-12-0149.1
! (2) Dai et al.,2019: A Global High-Resolution Data Set of Soil Hydraulic and Thermal Properties
!     for Land Surface Modeling. J. of Advances in Modeling Earth Systems, DOI: 10.1029/2019MS001784
!
! Original author: Yongjiu Dai, Wei Shangguan, 12/2013/
!
! Revisions:
! Yongjiu Dai, Nan Wei and Yonggen Zhang, 06/2018: add more highly cited or newly developed soil Pedotransfer functions.
! ----------------------------------------------------
USE MOD_Precision

IMPLICIT NONE
real(r8), intent(in) :: BD   ! bulk density (g cm-3)
real(r8), intent(in) :: SOM  ! soil organic matter (%, weight)
real(r8), intent(in) :: SOC  ! soil organic carbon concentration (%, weight)
real(r8), intent(in) :: SAND ! percentage of sand particle-size distribution (%, weight)
real(r8), intent(in) :: CLAY ! percentage of clay particle-size distribution (%, weight)
real(r8), intent(in) :: psi  ! (cm)
real(r8), intent(in) :: lambda  !
real(r8), intent(in) :: TOPSOIL ! 1 IF A or E horizon, 0 IF the subsoil
real(r8), intent(in) :: phi     ! saturated water content (cm3/cm3)
real(r8), intent(in) :: k_s_Rose

real(r8), intent(out) :: k_s  ! SATURATED HYDRAULIC CONDUCTIVITY (m/d)

integer i,j,iclass
real(r8) k(23)        ! SATURATED HYDRAULIC CONDUCTIVITY (m/d)
real(r8) :: SILT      ! percentage of silt particle-size distribution
real(r8) :: theta_33  ! the water content at a potential of -33 kPa
real(r8) :: phi_e     ! phi minus theta_33
real(r8) x,B,zs,z33,zs33,z1500,lam_g,d_g,sigma_g
logical c(12) ! indicate whether a soil is in an class
integer itype
real params(7,10)
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

real rahmati_ks_es(12), rahmati_ks_me(12)
data rahmati_ks_es/42.2,61.4,32.0,7.8,7.4,7.4,8.3,26.5,26.5,10.6,26.2,354.3/ ! cm/h
data rahmati_ks_me/43.6,24.6,41.2,4.9,5.4,5.4,7.6,29.0,29.0,12.3,44.8,148.8/ ! cm/h

real clapp_ks(12), cosby_ks(12)
data clapp_ks/1.056,0.938,0.208,0.0417,0.0378,0.0130,0.0147,0.0432,0.0432,0.0102,0.062,0.077/ !cm/min
data cosby_ks/0.82,0.3,-0.13,-0.32,-0.2,0.01,-0.46,-0.4,-0.4,-0.54,-0.72,-0.86/ !log10 inches/h

! ------------------------------------------------------------
      SILT = max(1.,100.0-SAND-CLAY)

! GROUP I: input soil physical parameters: Texture, Bulk Density, Porosity
! ------------------------------------------------------------
! 1) Ahuja et al. (1989):
! Ahuja, L.R., D.K. Cassel, R.R. Bruce, and B.B. Barnes. 1989. Evaluation of
! spatial distribution of hydraulic conductivity using eff ective porosity data.
! Soil Sci. 148:404-411.
! ------------------------------------------------------------
      theta_33 = phi*(psi/(-33.*10.))**lambda     ! 1 kPa = 10 cm H2O
      phi_e=phi-theta_33

      k(1)=18348.*phi_e**3.295                                                             ! cm/d

! ------------------------------------------------------------
! 2) Suleiman et al. (2001):
! Suleiman, A. A., and J. T. Ritchie. 2000. Estimating Saturated Hydraulic Conductivity from Soil
! Porosity. American Society of Agricultural Engineers Vol. 44(2): 235-339.
! ------------------------------------------------------------
      k(2)=12303.*phi_e**3.63                                                              ! cm/d

! ------------------------------------------------------------
! 3) Spychalski et al. (2007):
! Spychalski M., C. Kazmierowski, Z. Kaczmarek (2007): Estimation of saturated
! hydraulic conductivity on the basis of drainage porosity.
! Electronic J. Polish Agricultural Universities. 10(1), #04.
! for 0.003 < phi_e < 0.3, it provides a more accurate assessment.
! based on 35 samples.
! ------------------------------------------------------------
      k(3)=86400./10000.*(-2.52+581.598*phi_e**1.5-6966.14*phi_e**2.5+11693.78*phi_e**3)   ! cm/d

! ----------------------------------------------------------
! 4) Rawls et al. (1998):
! W. J. Rawls, D. Gimenez, R. Grossman, 1998: Use of soil texture, bulk density,
! and slope of thewater retention curve to predict saturated hydraulic conductivity
! Transactions of the ASAE. VOL. 41(4): 983-988
!
! 900 saturated hydraulic conductivities for the soil matrix according to the USDA
! soil texture classes and two porosity classes. The
! Brooks Corey pore size distribution index (lam_g) was obtained
! by fitting a log-log plot of water content vs pressure head
! using only the -33 and -1500 kPa water contents.
! ----------------------------------------------------------
      lam_g=exp(-0.7842831+0.0177544*SAND-1.062498*phi-0.00005304*SAND**2 &
            -0.00273493*CLAY**2+1.11134946*phi**2-0.03088295*SAND*phi &
            +0.00026587*(SAND*phi)**2-0.00610522*(CLAY*phi)**2 &
            -0.00000235*SAND**2*CLAY+0.00798746*CLAY**2*phi-0.0067449*CLAY*phi**2)
      k(4)=4632.*phi_e**(3.-lam_g)                                                         ! cm/d

! ------------------------------------------------------------
! 5) Cosby et al. (1984):
! B.J. COSBY, G.M. HORNBERGER, R.B. CLAPP, and T.R. GINN, 1984: A Statistical exploration
! of the Relationships of soil moisture characteristics to the physical properties of soils
! WATER RESOURCES RESEARCH, VOL. 20, NO. 6, PAGES 682-690, JUNE 1984.
!
! 1448 samples from 35 localities in 23 states in the United States.
! ------------------------------------------------------------
! Multiple regression:
      k(5)= 60.96*10.0**(-0.6+0.0126*SAND-0.0064*CLAY)                                     ! cm/d

! ------------------------------------------------------------
! 6) Cosby et al. (1984):
! Univariate regression:
      k(6)=60.96*10.0**(-0.884+0.0153*SAND)                                                ! cm/d

! ------------------------------------------------------------
! 7) Rosetta H3 by yonggen Zhang
      k(7)=k_s_Rose

! ----------------------------------------------------------
! 8) Dane and Puckett (1994):
! Dane, J.H., W. Puckett. 1994. Field soil hydraulic properties based on physical and
! mineralogical information. In ¡®Proceedings of the International Workshop on Indirect
! Methods for Estimating the Hydraulic Properties of Unsaturated Soils. Eds MTh van
! Genuchten et al. pp. 389-403. University of California: Riverside, CA.
! ----------------------------------------------------------
      k(8)=729.16*exp(-0.144*CLAY)                                                         ! cm/d

! ----------------------------------------------------------
! 9) Jabro(1992):
! Jabro, J.D. 1992. Estimation of saturated hydraulic conductivity of soils from particle size
! distribution and bulk density data. Transactions of ASAE 35 (2), 557-560.
!
! 350 soil samples of varying types
! ----------------------------------------------------------
      k(9)=24.*10.**(9.56-0.82*log10(SILT)-1.09*log10(CLAY)-4.64*BD)                       ! cm/d

! -----------------------------------------------------------
! 10) Brakensiek et al. (1984):
! Brakensiek, D. L., W. J. Rawls, G. R. Stephenson. 1984. Modifying SCS hydrologic soil groups
! and curve numbers for rangeland soils. ASAE paper no. PNR-84203, St. Joseph, MI.
!
! 1323 soil samples across the US.
! -----------------------------------------------------------
      x = 19.52348*phi-8.96847-0.028212*CLAY+0.00018107*(SAND)**2 &
        - 0.0094125*(CLAY)**2-8.395215*phi**2+0.077718*phi*SAND &
        - 0.00298*(phi*SAND)**2-0.019492*(phi*CLAY)**2+0.0000173*(SAND)**2*CLAY &
        + 0.02733*phi*(CLAY)**2+0.001434*phi*(SAND)**2-0.0000035*(CLAY)**2*SAND

      k(10)=24.0*exp(x)                                                                    ! cm/d

! ---------------------------------------------
! 11) Julia et al. (2004):
! Julia, M. F., T. E Monreal, A. Sanchez del Corral Jimeneza, and E. Garc Melendez. 2004.
! Constructing a saturated hydraulic conductivity map of Spain using Pedotransfer
! functions and spatial prediction. Geoderma 123, 257-277.
! ---------------------------------------------
      k(11)=2.208*exp(0.0491*SAND)                                                         ! cm/d

! ----------------------------------------------------------
! 12) Campbell (1985):
! Campbell GS (1985) Soil physics with basic development in soil science, 14. Elsevier, Amsterdam.
!
! d_g=geometric mean particle size (mm),
! sigma_g=geometric standard deviation of the particle size distribution,
! Note that this formula would lead to extremely inaccurate results IF the bulk density
! were low (BD << 1 ).  ! the arithmetic means 1.025mm for sand (2-0.05mm),
! 0.026mm for silt (0.05-0.002mm), 0.001 for clay (<0.002mm)
! ----------------------------------------------------------
      d_g=exp(log(1.025)*SAND/100.+log(0.026)*SILT/100.+log(0.001)*CLAY/100.)
      sigma_g=exp(sqrt((log(1.025))**2*SAND/100.+(log(0.026))**2*SILT/100. &
             +(log(0.001))**2*CLAY/100.-(log(1.025)*SAND/100.+log(0.026)*SILT/100. &
             +log(0.001)*CLAY/100.)**2))
      B=1./sqrt(d_g)+0.2*sigma_g
      k(12)=339.0*(1.3/BD)**(1.3*B)*exp(-6.88*CLAY/100.-3.63*SILT/100.-0.025)              ! cm/d

! GROUP II: input soil physical parameters: Texture, Bulk Density and Organic Matter/Carbon
! ---------------------------------------------
! 13) Vereecken et al. (1990):
! Vereecken, H., J. Maes, and J. Feyen. 1990. Estimating unsaturated hydraulic
! conductivity from easily measured soil properties. Soil Sci. 149:1-12.
!
! using soils from Belgium.
! ---------------------------------------------
      IF(SOM > 0.01)THEN
      k(13)=exp(20.62-0.96*log(CLAY)-0.66*log(SAND)-0.46*log(SOM)-8.43*BD)                 ! cm/d
      ELSE
         k(13)=-1.e36
      ENDIF

! ---------------------------------------------
! 14) Wosten et al. (1999):
! Wosten, J.H.M., Lilly, A., Nemes, A., Le Bas, C., 1999. Development and use of a
! database of hydraulic properties of European soils. Geoderma 90,169-185.
!
! HYPRES (EU-wide) soil samples 5521 soil horizons
! WHERE TOPSOIL which is a categorical variable, having a value of 1
! IF the soil sample comes from the topsoil (i.e., A or E horizon,
! according to the FAO soil classification [FAO, 1990] or 0 IF it is from the subsoil).
! topsoil: 2.272 < k < 60
! subsoil: 4.0 < k < 70.0
! oganic: 8.0
! the properties of EU soil datasets used to develop the PTF:
! 0.49<sand<100,   mean 30.24. s.d. 28.51
! 0   <silt<81.16, mean 41.93, s.d. 20.59
! 0   <clay<80,    mean 27.8,  s.d. 17.44
! 0.9 <BD  <1.9,   mean 1.46,  s.d. 0.19
! 0.01<k <2423,  mean 129,   s.d. 293
! ---------------------------------------------
      IF(SAND>0.49 .or. SILT<81.16 .or. CLAY<80. .or. (BD>0.9 .and. BD<1.9))THEN
      k(14)=exp(7.75+0.0352*SILT+0.93*(TOPSOIL)-0.967*BD**2 &
             - 0.000484*(CLAY)**2-0.000322*(SILT)**2 + 0.001/(SILT) &
             - 0.0748/SOM-0.643*log(SILT)-0.01398*BD*CLAY-0.1673*BD*SOM &
             + 0.02986*(TOPSOIL)*CLAY-0.03305*(TOPSOIL)*SILT)                              ! cm/d
      ELSE
         k(14)=-1.e36
      ENDIF

! -------------------------------------------
! 15) Merdun (2010) MLR:
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
! 0.1 < k < 3844.8 (cm/d); mean: 334.0
! -------------------------------------------
      IF(((SAND>=2.0) .and. (SAND<=96.0)) .and. &
         ((SILT>=1.0) .and. (SILT<=90.0)) .and. &
         ((CLAY>=0.1) .and. (CLAY<=63.3)) .and. &
         ((SOM>=0.01) .and. (SOM<=21.4))  .and. &
         ((BD>=0.59)  .and. (BD<=1.76)))THEN
         k(15)=9509.-14437.*SILT/100.-8169.*BD &
                -860.1*SOM+2332.*(SILT/100.)**2+1620.*BD**2+9.113*SOM**2 &
                +7547.*SILT/100.*BD+985.1*SILT/100.*SOM+381.4*BD*SOM                       ! cm/d

! 16) Merdun (2010) SUR:
         k(16)=8777.937-12556.9*SILT/100.-7849.05*BD &
                -728.977*SOM+2100.45*(SILT/100.)**2+1666.81*BD**2 &
                +6.79971*SOM**2+6597.450*SILT/100.*BD &
                +736.149*SILT/100.*SOM+371.5434*BD*SOM                                     ! cm/d
      ELSE
         k(15:16)=-1.0e36
      ENDIF

! ----------------------------------------------------------
! 17) Aimrun and Amin (2009):
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
! the predicted value: 0.126<k<2.33; mean: 0.489 cm/d
! ----------------------------------------------------------
      IF(((SILT>=19.84) .and. (SILT<=57.25)) .and. &
         ((CLAY>=32.22) .and. (CLAY<=76.80)) .and. &
         ((SOM>=0.07)   .and. (SOM <=29.35)) .and. &
         ((BD>=0.62)    .and. (BD<=1.91   )))THEN
         d_g=exp(SAND/100.*log(1.025)+SILT/100.*log(0.026)+CLAY/100.*log(0.001))
         k(17)=100.0*exp(-2.368+3.846*BD+0.091*SOM-6.203*log(BD) &
                -0.343*log(SOM)-2.334*log(CLAY)-0.411*log(d_g))                            ! cm/d
      ELSE
         k(17)=-1.0e36
      ENDIF

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
      IF(((SAND>=5.6) .and. (SAND<=97.8)) .and. &
         ((CLAY>=0.)  .and. (CLAY<=54.5)) .and. &
         ((SOC>=0.01) .and. (SOC<=6.6))   .and. &
         ((BD>=0.89)  .and. (BD<=1.77)))THEN
         k(18)=exp(1.9582+0.0308*SAND-0.6142*BD-0.1566*SOC/10.)
      ELSE
         k(18)=-1.0e36
      ENDIF

      DO i = 1, 18
!         IF(k(i) >=9000.) k(i)=-1.0e36
         IF(k(i) < 0.01) k(i)=-1.0e36
      ENDDO

! ------------------------------------------
! 19) Wosten et al. (1999):
! FAO class
IF(CLAY.gt.60.) THEN
   iclass=5
ELSEIF(CLAY.gt.35.) THEN
   iclass=4
ELSEIF(SAND.lt.15.) THEN
   iclass=3
ELSEIF(SAND.gt.65.0.and.CLAY.lt.18.0) THEN
   iclass=1
ELSE
   iclass=2
ENDIF

j=iclass+nint(TOPSOIL)*5

k(19)=params(7,j)

! ------------------------------------------
! soil_classes=c('Sa','LoSa','SaLo','Lo','SaClLo','SaCl','ClLo','SiLo','Si','SiClLo','SiCl','Cl')
      CALL USDA_soil_classes(SILT,CLAY,c)

      IF(c(1))  itype = 12 ! ic = 'Cl'      clay
      IF(c(2))  itype = 11 ! ic = 'SiCl'    silty clay
      IF(c(3))  itype = 6  ! ic = 'SaCl'    sandy clay
      IF(c(4))  itype = 7  ! ic = 'ClLo'    clay loam
      IF(c(5))  itype = 10 ! ic = 'SiClLo'  silty clay loam
      IF(c(6))  itype = 5  ! ic = 'SaClLo'  sandy clay loam
      IF(c(7))  itype = 4  ! ic = 'Lo'      loam
      IF(c(8))  itype = 8  ! ic = 'SiLo'    silty loam
      IF(c(9))  itype = 3  ! ic = 'SaLo'    sandy loam
      IF(c(10)) itype = 9  ! ic = 'Si'      silt
      IF(c(11)) itype = 2  ! ic = 'LoSa'    loamy sand
      IF(c(12)) itype = 1  ! ic = 'Sa'      sand

! Rahmati et al., 2018, Table 10
      k(20)=rahmati_ks_es(itype)*24.0         ! cm/h -> cm/d
      k(21)=rahmati_ks_me(itype)*24.0         ! cm/h -> cm/d

! Clapp et al., 1978, Table 2
      k(22)=clapp_ks(itype)*(60.0*24.0)       ! cm/min -> cm/d

! Cosby et al, 1984, Table 3
      k(23)=10.0**(cosby_ks(itype))*2.54*24.0 ! log10 inches/h -> cm/d

! ------------------------------------------
! THE MEDIAN VALUE
! ------------------------------------------
      CALL debar(k,23,k_s)

END SUBROUTINE ksat


!======================================================
!SUBROUTINE thetas(BD,SOC,sand,clay,soildepth,theta_s_l)
!======================================================
! (1) Dai et al., 2013, Development of a China Dataset of Soil
!     Hydraulic Parameters Using Pedotransfer Functions for Land Surface Modeling.
!     J. of Hydrometeorology, 14: 869-887. DOI: 10.1175/JHM-D-12-0149.1
! (2) Dai et al.,2019: A global high-resolution dataset of soil thermal
!                      and hydraulic properties for land surface modeling
!
! Yongjiu Dai, 12/2013, 06/2018
! ----------------------------------------------------
!USE MOD_Precision
!
!IMPLICIT NONE
!
!real(r8), intent(in) :: BD   ! bulk density (g cm-3)
!real(r8), intent(in) :: SOC  ! soil organic carbon concentration (%, weight)
!real(r8), intent(in) :: sand ! percent sand particle-size distribution (%)
!real(r8), intent(in) :: clay ! percent clay particle-size distribution (%)
!real(r8), intent(in) :: soildepth   ! soil depth (cm)
!real(r8), intent(out) :: theta_s_l  ! saturated water content (cm3/cm3)
!
!
!real(r8) theta_s(18)
!real(r8) SOM  ! percent organic matter (%, weight)
!real(r8) TOPSOIL  ! 1 IF A or E horizon, 0 IF the subsoil
!real(r8) phi ! soil porosity (cm3 cm-3)
!real(r8) silt     ! percent silt particle-size distribution (%)
!real(r8) x,z1,z2
!integer i,j,iclass
!real(r8) params(7,10)
!data params/0.025,0.403,0.0383,1.3774,0.2740, 1.2500,60.000, &
!            0.010,0.439,0.0314,1.1804,0.1528,-2.3421,12.061, &
!            0.010,0.430,0.0083,1.2539,0.2025,-0.5884, 2.272, &
!            0.010,0.520,0.0367,1.1012,0.0919,-1.9772,24.800, &
!            0.010,0.614,0.0265,1.1033,0.0936, 2.5000,15.000, &
!            0.025,0.366,0.0430,1.5206,0.3424, 1.2500,70.000, &
!            0.010,0.392,0.0249,1.1689,0.1445,-0.7437,10.755, &
!            0.010,0.412,0.0082,1.2179,0.1789, 0.5000, 4.000, &
!            0.010,0.481,0.0198,1.0861,0.0793,-3.7124, 8.500, &
!            0.010,0.538,0.0168,1.0730,0.0680, 0.0001, 8.235/
!
!! ---------------------------------------------------------
!      IF(soildepth < 30.)THEN   ! cm
!         TOPSOIL=1.  ! 1 IF A or E horizon
!      ELSE
!         TOPSOIL=0.  ! 0 IF the subsoil
!      ENDIF
!
!      SOM=1.724*SOC
!
!      silt=max(1.,100.0-sand-clay)
!      phi=1.0-BD/2.65
!
!! 1) Saturated water content = total porosity:
!      theta_s(1)=1.0-BD/2.65
!
!! ---------------------------------------------------------
!! 2) Williams et al. (1992):
!! Williams, J., P. Ross, and K. Bristow. 1992. Prediction of the Campbell water
!! retention function from texture, structure, and organic matter. p. 427-442.
!! In M.Th . van Genuchten et al. (ed.) Proc. Int. Worksh. on Indirect Methods
!! for Estimating the Hydraulic Properties of Unsaturated Soils, Riverside,
!! CA. 11-13 Oct. 1989. U.S. Salinity Lab., Riverside, CA.
!!
!! based on large soil data sets from Australia, UK, and USA
!! set the value of the residual water content equal to zero
!! function # 1 all sand is fine sand and the structural index = 1
!! ---------------------------------------------------------
!      theta_s(2)=0.93*(1.0-BD/2.65)
!
!! -------------------------------------------
!! 3) Jauhiainen (2004):
!! Jauhiainen, M., Relationships of particle size distribution curve, soil water retention curve
!! and unsaturated hydraulic conductivity and their implications on water balance of forested
!! and agricultural hillslopes, Helsinki University of Technology, Water Resources Publications,
!! TKK-VTR-12, Espoo, 2004, 165 pp.
!! -------------------------------------------
!      theta_s(3)=0.928*(1.0-BD/2.65)+0.021
!
!! -------------------------------------------
!! 4) Teepe et al. (2003):
!! Teepe, R., H. Dilling, and F. Beese, 2003: Estimating water retention curves of forest soils
!! from soil texture and bulk density. J. Plant Nutr. Soil Sci., 166, 11-119.
!!
!! the data of this study was a collective of 1850 water retention curve
!! from FOREST soils in German
!! -------------------------------------------
!      theta_s(4)=0.9786-0.36686*BD
!
!! ------------------------------------------------------------
!! 5) Cosby et al. (1984):
!! B.J. COSBY, G.M. HORNBERGER, R.B. CLAPP, and T.R. GINN, 1984: A Statistical exploration
!! of the Relationships of soil moisture characteristics to the physical properties of soils
!! WATER RESOURCES RESEARCH, VOL. 20, NO. 6, PAGES 682-690, JUNE 1984.
!!
!! 1448 samples from 35 localities in 23 states in the United States.
!! ------------------------------------------------------------
!      theta_s(5)=0.489-0.00126*sand  ! (univariate regression)
!
!! 6) Cosby et al. (1984):
!      theta_s(6)=0.505-0.00142*sand-0.00037*clay ! (multi-variate regression)
!
!! ------------------------------------------------------------
!! 7) Saxton et al. (1986):
!! Saxton, K.E., W.J. Rawls, J.S. Romberger, and R.I. Papendick. 1986.
!! Estimating generalized soil water characteristics from texture.
!! Soil Sci. Soc. Am. J. 50: 1031-1036
!! ------------------------------------------------------------
!      IF((sand>=5.0) .and. (clay>=5.0 .and. clay<=60.))THEN
!         theta_s(7)=0.332-0.0007251*sand+0.1276*log10(clay)
!      ELSE
!         theta_s(7)=-1.e36
!      ENDIF
!
!! ------------------------------------------------------------
!! 8) Vereeken et al. (1989):
!! Vereecken, H., J. Maes, J. Feyen, and P. Darius. 1989. Estimating the soil moisture
!! retention characteristics from texture, bulk density and carbon content.
!! Soil Sci. 148:389-403.
!! ------------------------------------------------------------
!      theta_s(8)=0.81-0.28*BD+0.001*clay
!
!!--------------------------------------------------------------
!! 9) Wosten et al. (1999) :
!! Wosten, J.H.M., A. Lilly, A. Nemes, and C. Le Bas. 1999. Development and
!! use of a database of hydraulic properties of European soils. Geoderma 90:169-185.
!!
!! Wosten et al. (1999) analyzed the all-Europe database and derived
!! the following PTFs to estimate van Genuchten parameters:
!! WHERE topsoil is an ordinal variable having the value of 1 or
!! of 0. %silt is the percentage of soil (2mm-50mm), BD= bulk density (g/cm3 = Mg/m3)
!! the PTFs developd by using the European HYdraulic PRoperties of European Soils (HYPRES)
!! soil database (5521 horizons)
!! topsoil: 0.403 < theta_s < 0.614 (oganic: 0.766)
!! subsoil: 0.366 < theta_s < 0.538
!!--------------------------------------------------------------
!      theta_s(9)=0.7919+0.001691*clay-0.29619*BD-0.000001491*(silt)**2 &
!                +0.0000821*(SOM)**2+0.02427/(clay)+0.01113/(silt) &
!                +0.01472*log(silt)-0.0000733*(SOM)*(clay)-0.000619*BD*(clay)&
!                -0.001183*BD*SOM-0.0001664*(TOPSOIL)*(silt)
!
!! 10) Wosten et al. (1999):
!! FAO class
!      IF(clay.gt.60.) THEN
!         iclass=5
!      ELSEIF(clay.gt.35.) THEN
!         iclass=4
!      ELSEIF(sand.lt.15.) THEN
!         iclass=3
!      ELSEIF(sand.gt.65.0.and.clay.lt.18.0) THEN
!         iclass=1
!      ELSE
!         iclass=2
!      ENDIF
!
!      j=iclass+nint(TOPSOIL)*5
!      theta_s(10)=params(2,j)
!
!! -------------------------------------------------------------
!! 11) Mayr and Jarvis (1999):
!! T. Mayr, N.J. Jarvis, 1999: Pedotransfer functions to estimate soil water
!! retention parameters for a modified Brooks-Corey type model. Geoderma 91:1-9
!!
!! 286 soil horizons of the soil physical properties database of England and Wales
!! set the value of the residual water content equal to zero
!! these function should not be applied to organic soils (organic carbon content >5%)
!! and/or soils of small dry bulk density (<0.9 g/cm3)
!! -------------------------------------------------------------
!      IF(SOC < 5. .or. BD > 0.9)THEN
!         theta_s(11)=0.2345971971 &
!                 +0.0046614221*sand+0.0088163314*silt+0.0064338641*clay-0.3028160229*BD &
!                 +1.79762e-05*(sand)**2-3.134631e-05*(silt)**2
!      ELSE
!         theta_s(11)=-1.0e36
!      ENDIF
!
!! -------------------------------------------
!! 12) Merdun (2010) MLR PTFs:
!! H. Merdun 2010: Alternative Methods in the Development of Pedotransfer Functions
!! for Soil Hydraulic Characteristics. Eurasian Soil Science,Vol. 43, No. 1, 62-71.
!!
!! data from the UNSODA database
!! 135 soil samples for the development and 45 for the validation
!! 2.0  < sand <96.0   mean 47.2     s.d. 29.7
!! 1.0  < silt <90.0   mean 35.2     s.d. 22.7
!! 0.1  < clay <63.3   mean 0.176    s.d. 14.5
!! 0.01 < SOM  <21.4;  mean 1.673    s.d. 2.478
!! 0.59 < BD   <1.76   mean 1.42     s.d. 0.22
!! 0.274 < theta_s < 0.837; mean: 0.443
!! -------------------------------------------
!      IF(((sand>=2.0) .and. (sand<=96.0)) .and. &
!         ((silt>=1.0) .and. (silt<=90.0)) .and. &
!         ((clay>=0.1) .and. (clay<=63.3)) .and. &
!         ((SOM>=0.01) .and. (SOM<=21.4))  .and. &
!         ((BD>=0.59)  .and. (BD<=1.76)))THEN
!         theta_s(12)=1.391-0.289*clay/100.-1.007*BD &
!                    -0.026*SOM-0.096*(clay/100.)**2+0.22*BD**2-0.00039*SOM**2 &
!                    +0.3*clay/100.*BD+0.0233*clay/100.*SOM+0.0229*BD*SOM
!
!! 13) Merdun (2010) SUR PTFs:
!         theta_s(13)=1.419680-0.37696*clay/100.-1.04082*BD &
!                    -0.02362*SOM+0.08548*(clay/100.)**2+0.23091*BD**2+0.00002*SOM**2+0.32229*clay/100.*BD &
!                    +0.004189*clay/100.*SOM +0.022525*BD*SOM
!      ELSE
!         theta_s(12:13)=-1.0e36
!      ENDIF
!
!! ----------------------------------------------------------
!! 14) Saxton and Rawls (2006):
!! K. E. Saxton and W. J. Rawls, 2006: Soil Water Characteristic Estimates by Texture and
!! Organic Matter for Hydrologic Solutions. Soil Sci. Soc. Am. J. 70:1569-1578.
!!
!! A-horizon samples from USDA/NRCS National Soil Characterization Dataset.
!! Samples with "extreme" values were omitted from the data.
!! Excluded were those with bulk density < 1.0 and > 1.8 g cm-3,
!! OM > 8 % (w) and clay > 60% (w). This reduced the A-horizon data set
!! from 2149 to 1722 samples. Set the value of the residual water content equal to zero
!! the fitting samples with
!! 1.0  < BD <1.8 g/cm3
!! SOM  < 8  %(w)
!! clay < 60 %(w).
!! ----------------------------------------------------------
!      IF((clay<=60.0) .and. (SOM<=8.0) .and. (BD>=1.0 .and. BD<=1.8))THEN
!         x=-0.251*sand/100.+0.195*clay/100.+0.011*SOM &
!           +0.006*(sand/100.*SOM)-0.027*(clay/100.*SOM) &
!           +0.452*(sand/100.*clay/100.)+0.299
!         z1=x+(1.283*x**2-0.374*x-0.015)
!         z2=-0.107+1.636*(0.278*sand/100.+0.034*clay/100.+0.022*SOM &
!           -0.018*(sand/100.*SOM)-0.027*(clay/100.*SOM) &
!           -0.584*(sand/100.*clay/100.)+0.078)
!         theta_s(14)=z1+z2-0.097*sand/100.+0.043
!         IF(theta_s(14) <= z2) theta_s(14) = -1.0e36  ! theta_s < field capacity
!      ELSE
!         theta_s(14)=-1.0e36
!      ENDIF
!
!! ----------------------------------------------------------
!! 15) Tomasella and Hodnett (1998):
!! Tomasella, J., and M.G. Hodnett. 1998. Estimating soil water retention characteristics
!! from limited data in Brazilian Amazonia. Soil Sci. 163:190-202.
!! ----------------------------------------------------------
!      theta_s(15)=(37.937+2.24*SOC+0.298*silt+0.159*clay)/100.
!
!! ------------------------------------------
!! 16) Li et al. (2007):
!! Y. Li, D. Chen, R.E. White, A. Zhu, J. Zhang, 2007: Estimating soil hydraulic
!! properties of Fengqiu County soils in the North China Plain using
!! pedo-transfer functions. Geoderma, 138:261-271.
!!
!! the residual water content ranged from 0.002-0.0001 cm3 cm-3
!! when theta_r<0.001, set theta_r=0.
!! 63 samples for SWRC,
!! 8.98< sand <93.03;    mean 53.28       s.d. 4.46
!! 1.74< silt <79.5      mean 35.86       s.d. 3.83
!! 0.54< clay <27.12;    mean 8.86        s.d. 1.01
!! 0.12< SOM  <1.54;     mean 0.65        s.d. 0.07
!! 1.2 < BD   <1.59      mean 1.42        s.d. 0.01
!! all 63 samples: theta_s=0.477 +/- 0.006
!! ------------------------------------------
!      IF(((sand>=8.98) .and. (sand<=93.03)) .and. &
!         ((silt>=1.74) .and. (silt<=79.5))  .and. &
!         ((clay>=0.54) .and. (clay<=27.12)) .and. &
!         ((SOM>=0.12)  .and. (SOM<=1.54))   .and. &
!         ((BD>=1.2)    .and. (BD<=1.59)))    THEN
!         theta_s(16)=exp(-1.531+0.212*log(sand)+0.006*silt-0.051*SOM-0.566*log(BD))
!      ELSE
!         theta_s(16)=-1.0e36
!      ENDIF
!
!! ----------------------------------------------------------
!! 17) Al Majou et al (2007)
!! Al Majou, H., A. Bruand, and O. Duval. 2008. The use of in situ volumetric
!! water content at fi eld capacity to improve the prediction of soil water
!! retention properties. Can. J. Soil Sci. 88:533-541.
!!
!! Class and continuous ptfs were developed using a set of 320 horizons, comprising 90 topsoils
!! (from 0 to 30 cm depth) and 230 subsoil horizons (> 30 cm depth)
!! collected in Cambisols, Luvisols, Planosols, Albeluvisols,Podzols, and Fluvisols located
!! mainly in the Paris basin and secondarily in the western coastal marshlands
!! and Pyrenean piedmont plain.
!! 0.1< sand < 90.1       mean 24.9         s.d. 23.9
!! 2.8< silt < 82.1       mean 46.2         s.d. 20.8
!! 1.9< clay <92.9,       mean 28.9         s.d. 15.1
!! 0.0< SOC <28.8 g kg-1, mean 5.7 g kg-1   s.d. 4.9
!! 1.0< BD <1.84,         mean 1.53         s.d. 0.15
!! topsoil: 0.397 < theta_s < 0.587
!! subsoil: 0.367 < theta_s < 0.742
!! A set of 107 horizons comprising 39 topsoil and 68 subsoil horizons
!! was constituted in order to test the ptfs established. These horizons were collected in Cambisols,
!! Luvisols and Fluvisols located in the South of the Paris basin
!
!! the residual water content was fixed at 0.01 cm3 cm-3
!! except for texture coarse for which it was fixed at 0.025 cm3 cm-3
!! ----------------------------------------------------------
!      IF(((sand>=0.1) .and. (sand<=90.1)) .and. &
!         ((silt>=2.8) .and. (sand<=82.1)) .and. &
!         ((clay>=1.9) .and. (clay<=92.9)) .and. &
!         ((SOC>=0.0)  .and. (SOC<=2.88 )) .and. &
!         ((BD>=1.0)   .and. (BD<=1.84)))   THEN
!         theta_s(17)=1.1658-0.0032*clay-0.4737*BD+2.0e-07*(silt)**2-0.0001*(SOC)**2 &
!                    +0.0373/(SOC)+0.0131/silt-0.0072*log(silt)+0.00003*(SOC)*clay+0.0022*BD*clay &
!                    -0.0002*BD*(SOC)-0.0001*silt
!      ELSE
!         theta_s(17)=-1.0e36
!      ENDIF
!
!! ------------------------------------------
!! 18) Weynants et al. (2009):
!! M. Weynants, H. Vereecken, and M. Javaux, 2009: Revisiting Vereecken Pedotransfer
!! Functions: Introducing a Closed-Form Hydraulic Model. Vadose Zone J. 8:86-95.
!!
!! datasets: 136 hydraulic conductivity curve and 166 moisture retention curve in the middle of
!! each of the 182 horizons from northern Belgium. Training data value ranges:
!! 5.6  < sand < 97.80
!! 0    < clay < 54.50
!! 0.1  < SOC  < 66 (g kg-1)
!! 0.89 < BD   < 1.77 (g cm-3)
!! the theta_r could be 0.
!! ------------------------------------------
!      IF(((sand>=5.6) .and. (sand<=97.8)) .and. &
!         ((clay>=0.)  .and. (clay<=54.5)) .and. &
!         ((SOC>=0.01) .and. (SOC<=6.6))   .and. &
!         ((BD>=0.89)  .and. (BD<=1.77)))THEN
!         theta_s(18)=0.6355+0.0013*clay-0.1631*BD
!      ELSE
!         theta_s(18)=-1.0e36
!      ENDIF
!
!      DO i = 1, 18
!         IF(abs(theta_s(i)) > 0.8) theta_s(i)=-1.0e36
!      ENDDO
!
!! ------------------------------------------
!! THE MEDIAN VALUE
!! ------------------------------------------
!      CALL debar(theta_s,18,theta_s_l)
!
!END SUBROUTINE thetas
!

!======================
SUBROUTINE debar(a,n,x)
!======================
USE MOD_Precision

IMPLICIT NONE
integer n
real(r8) a(n), tmp(n)
real(r8) x
integer i, ii
real(r8), allocatable :: c(:)
real(r8), external :: median

      ii = 0
      DO i = 1, n
         IF(abs(a(i)) < 1.0e10)THEN
           ii = ii + 1
           tmp(ii) = a(i)
         ENDIF
      ENDDO
      allocate (c(ii))
      c = tmp(1:ii)
      x =  median(c,ii)
      deallocate (c)

END SUBROUTINE debar

SUBROUTINE GuLar(Clay,Silt,Sand,OC,BD,Theta)
USE MOD_Precision
IMPLICIT NONE
real(r8) Clay,Silt,Sand,OC,BD
real(r8) a(5,12),Theta(14)
integer i
real(r8) OM
     data a/7.053,10.242,10.07,6.333,-321.2, &
            5.678,9.2280,9.135,6.103,-269.6, &
            5.018,8.5480,8.833,4.966,-242.3, &
            3.890,7.0660,8.408,2.817,-187.8, &
            3.075,5.8860,8.039,2.208,-143.4, &
            2.181,4.5570,7.557,2.191,-92.76, &
            1.563,3.6200,7.154,2.388,-57.59, &
            0.932,2.6430,6.636,2.717,-22.14, &
            0.483,1.9430,6.128,2.925,-2.040, &
            0.214,1.5380,5.908,2.855,15.300, &
            0.076,1.3340,5.802,2.653,21.450, &
            -0.059,1.142,5.766,2.228,26.710/
	OM=OC*1.724
      DO i=1,12
      Theta(i+1)=(Sand*a(1,i)+Silt*a(2,i)+Clay*a(3,i)+OM*a(4,i)+BD*a(5,i))/1000.
      ENDDO
END SUBROUTINE GuLar

SUBROUTINE Rajkai(Clay,Sand,BD,OC,Theta)
USE MOD_Precision
IMPLICIT NONE
real(r8) Clay,Sand,BD,OC
real(r8) b(6,8), Theta(14)
integer i
real(r8) Silt,X1,X2
	Data b/89.75,  -31.39,    0,    0.030,       0,       0, &
             85.05,  -27.17,      0,   -0.024,       0,       0, &
             78.58,  -23.94,      0,   -0.025,       0,       0, &
             69.78,  -21.74,      0,        0,       0,  0.0011, &
             20.87,    0.29,  -0.83,    0.030,       0,  0.0051, &
              2.19,    0.52,   3.93,    -0.07,       0,       0, &
              1.39,    0.36,      0,        0,       0,  0.220 , &
              0.73,       0,   0.32,        0,  0.0018,       0/
	Silt=100-Clay-Sand
        DO i=1,8
           select CASE (i)
		CASE (1)
			X1=BD
			X2=Silt
		CASE (2,3)
			X2=Sand
		CASE(4)
			X2=Clay+Silt
		CASE(5)
			X1=Clay+Silt
			X2=Sand/Silt
		CASE(6,7)
			X2=OC*1.724
		CASE(8)
			X1=Clay
	   End select
	   Theta(i)=(b(1,i)+b(2,i)*X1+b(3,i)*X2+b(4,i)*X1*X2+b(5,i)*X1*X1 &
      	            +b(6,i)*X2*X2)/100.
	ENDDO
END SUBROUTINE Rajkai

SUBROUTINE Rawls82(Clay,Silt,Sand,OC,Theta)
USE MOD_Precision
IMPLICIT NONE
real(r8) Clay,Silt,Sand,OC
real(r8) a(5,10), Theta(14)
integer i
	Data a/0.4118, -0.0030  ,       0   , 0.0023 ,  0.0317, &
             0.3121 ,  -0.0024  ,       0   , 0.0032 ,  0.0314, &
             0.2576 ,  -0.0020  ,       0   , 0.0036 ,  0.0299, &
             0.2065 ,  -0.0016  ,       0   , 0.0040 ,  0.0275, &
             0.0349 ,        0  ,   0.0014  , 0.0055 ,  0.0251, &
             0.0281 ,        0  ,   0.0011  , 0.0054 ,  0.0220, &
             0.0238 ,        0  ,   0.0008  , 0.0052 ,  0.0190, &
             0.0216 ,        0  ,   0.0006  , 0.0050 ,  0.0167, &
             0.0205 ,        0  ,   0.0005  , 0.0049 ,  0.0154, &
             0.0260 ,        0  ,       0   , 0.0050 ,  0.0158/
      DO i=1,10
      Theta(i+1)=a(1,i)+Sand*a(2,i)+Silt*a(3,i)+Clay*a(4,i)+OC*a(5,i)
      ENDDO
END SUBROUTINE Rawls82

SUBROUTINE Rawls83(Clay,Silt,Sand,BD,OC,Theta)
USE MOD_Precision
IMPLICIT NONE
real(r8) Clay,Silt,Sand,BD,OC
real(r8) a(5,9), Theta(14)
integer i
	Data a/0.4180, -0.0021  , 0.0035 ,  0.0232, -0.0859, &
             0.3486 ,  -0.0018  , 0.0039 ,  0.0228, -0.0738, &
             0.2819 ,  -0.0014  , 0.0042 ,  0.0216, -0.0612, &
             0.2352 ,  -0.0012  , 0.0043 ,  0.0202, -0.0517, &
             0.1837 ,  -0.0009  , 0.0044 ,  0.0181, -0.0407, &
             0.1426 ,  -0.0007  , 0.0045 ,  0.0160, -0.0315, &
             0.1155 ,  -0.0005  , 0.0045 ,  0.0143, -0.0253, &
             0.1005 ,  -0.0004  , 0.0044 ,  0.0133, -0.0218, &
             0.0854 ,  -0.0004  , 0.0044 ,  0.0122, -0.0182/
      DO i=1,9
      Theta(i+1)=a(1,i)+Sand*a(2,i)+Clay*a(3,i)+OC*a(4,i)+BD*a(5,i)
      ENDDO
END SUBROUTINE Rawls83

SUBROUTINE Tomasella(Clay,Silt,OC,Theta)
USE MOD_Precision
IMPLICIT NONE
real(r8) Clay,Silt,OC
real(r8) a(4,9),Theta(14)
integer i
	Data a/ 37.937,2.24,0.298,0.159, &
	 	23.839,0, 0.53,0.255, &
	 	18.495,0,0.552,0.262, &
	 	12.333,0,0.576,  0.3, &
	 	 9.806,0,0.543,0.321, &
	 	 4.046,0,0.426,0.404, &
	 	 3.198,0,0.369,0.351, &
	 	 1.567,0,0.258,0.361, &
	 	  0.91,0, 0.15,0.396/
      DO i=1,9
      Theta(i)=(a(1,i)+a(2,i)*OC+a(3,i)*Silt+a(4,i)*Clay)/100.
      ENDDO
END SUBROUTINE Tomasella


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Subprogram VGpar
!     This subprogram calculates parameters of the van Genuchten's
!     equation of soil water retention: ThR, ThS, Alpha,n, m
!     The van Genuchten's equation is used in its traditional form:
!       Theta = (ThS-ThR)/(1+(alpha*P)^n)^(1-1/n)
!     WHERE Theta stands for volumetric water content and P stands
!     for suction (or for absolute value of the matric potential)
!     The program uses raw water retention data
!
!     Input data:
!     y - water content (cm3/cm3)
!     x - suction (cm)
!     nob - number of observations
!     b - van Genuchten parameters
!     np - number of parameters
!     mit - maximun number of iterations
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +   The code is modification of van Genuchten, M.Th. 1980. Determining transport
! +   parameters from solute displacement experiments. Research report No.118,
! +   U.S. Salinity Laboratory, USDA-ARS-AR, Riverside, California.
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE VGpar(x,y,nob,b,np,mit)
USE MOD_Precision
IMPLICIT NONE
integer nob,np,mit
real(r8) y(nob),x(nob),f(50),r(50),st(5),b(5),e(5),&
         c(10), p(10), q(10), a(10,10), d(10,10),&
         delz(50,10),dz(50)
integer i,j,k,nit,IER
real(r8) ga,sumb,z,eps
real(r8) sdev,ssq,sum,sum1,sum2,sum3,step,temp,tmp,an,angle

      data eps /0.0005/

	b(1)=0.001
	b(2)=y(1)
	b(3)=0.01
	b(4)=1.5
      IER=0
      ga = 0.02
      sumb = 0.0
      CALL model(b,np,f,nob,x)
      DO 10 k = 1,nob
      z = y(k) - f(k)
      r(k) = z
      IF(abs(z) .gt. 1.0e-37) sumb = sumb + z * z
 10   continue

      DO 200 nit = 1,mit
      ssq = sumb
      ga = 0.1 * ga
      DO 30 j = 1,np
      temp = b(j)
      b(j) = 1.01 * b(j)
      CALL model(b,np,dz,nob,x)
      DO 15 i = 1,nob
      delz(i,j) = dz(i)
 15   continue
      sum = 0.0
      DO 20 k = 1,nob
      delz(k,j) = 100.0 * (delz(k,j) - f(k))
      tmp = delz(k,j) * r(k)
      sum = sum + tmp
 20   continue
      q(j) = sum / b(j)
      b(j) = temp
      c(j) = temp
 30   continue
      sum3 = 0.0
      DO 60 i = 1,np
      DO 50 j = 1,i
      sum = 0.0
      DO 40 k = 1,nob
      temp = delz(k,i) * delz(k,j)
      sum = sum + temp
 40   continue
      d(j,i) = sum / (b(j) * b(i))
      d(i,j) = d(j,i)
 50   continue
      e(i) = sqrt(d(i,i))
      If(e(i).le.0.) THEN
       IER=1
       Goto 500
      Endif
      q(i) = q(i) / e(i)
      IF(abs(q(i)) .gt. 1.0e-37) sum3 = sum3 + q(i) * q(i)
 60   continue
 70   DO 90 i = 1,np
      DO 80 j = 1,i
      a(j,i) = d(j,i) / e(j) / e(i)
      a(i,j) = a(j,i)
 80   continue
 90   continue
      DO 100 i = 1,np
      p(i) = q(i)
 100  a(i,i) = a(i,i) + ga
      CALL matinv(a,np,p)
      sum1 = 0.0
      sum2 = 0.0
      DO 110 i = 1,np
      temp = p(i) * q(i)
      sum1 = sum1 + temp
      temp = p(i) * p(i)
      sum2 = sum2 + temp
 110  continue
      an = sqrt((sum1/sum2)*(sum1/sum3))
      angle = 57.2958 * atan((sqrt(abs(1-an**2)))/an)
      step = 1.0
 120  DO 130 i = 1,np
 130  b(i) = p(i) * step / e(i) + c(i)
      DO 140 i = 1,np
      IF(c(i)*b(i) .le. 0.0) go to 160
 140  continue
      sumb = 0.0
      CALL model(b,np,f,nob,x)
      DO 150 k = 1,nob
      z = y(k) - f(k)
      r(k) = z
      IF(abs(z) .gt. 1.0e-37) sumb = sumb + z*z
 150  continue
      IF(sumb-ssq .lt. 1.0e-8) go to 180
 160  IF(angle .gt. 30.0) go to 170
      step = 0.5 * step
      go to 120
 170  ga = 10.0 * ga
      go to 70
 180  DO 190 i = 1,np
      IF(abs(c(i)-b(i)) .gt. eps*abs(b(i))) go to 200
 190  continue
      go to 210
 200  continue
 210  CALL matinv(d,np,p)
      sdev = sqrt(sumb/float(nob-np))
      DO 220 i = 1,np
      e(i) = sqrt(amax1(d(i,i),1.0e-20))
      st(i) = e(i) * sdev
 220  continue

 500  RETURN

END SUBROUTINE VGpar

SUBROUTINE matinv(a,np,b)
USE MOD_Precision
IMPLICIT NONE
integer np
real(r8) a(10,10),b(np)
integer indx1(10),indx2(10)
integer i,j,k,l,ir,ic
real(r8) p,amax

      DO 10 j = 1,np
 10   indx1(j) = 0
      i = 0
 20   amax = -1.0
      DO 40 j = 1,np
      IF(indx1(j) .ne. 0) go to 40
      DO 30 k = 1,np
      IF(indx1(k) .ne. 0) go to 30
      p = abs(a(j,k))
      IF(p .le. amax) go to 30
      ir = j
      ic = k
      amax = p
 30   continue
 40   continue
      IF(amax .le. 0.0) go to 120
      indx1(ic) = ir
      IF(ir .eq. ic) go to 60
      DO 50 l = 1,np
      p = a(ir,l)
      a(ir,l) = a(ic,l)
      a(ic,l) = p
 50   continue
      p = b(ir)
      b(ir) = b(ic)
      b(ic) = p
      i = i + 1
      indx2(i) = ic
 60   p = 1.0 / a(ic,ic)
      a(ic,ic) = 1.0
      DO 70 l = 1,np
      a(ic,l) = a(ic,l) * p
 70   continue
      b(ic) = b(ic) * p
      DO 90 k = 1,np
      IF(k .eq. ic) go to 90
      p = a(k,ic)
      a(k,ic) = 0.0
      DO 80 l = 1,np
      a(k,l) = a(k,l) - a(ic,l) * p
 80   continue
      b(k) = b(k) - b(ic) * p
 90   continue
      go to 20
 100  ic = indx2(i)
      ir = indx1(ic)
      DO 110 k = 1,np
      p = a(k,ir)
      a(k,ir) = a(k,ic)
      a(k,ic) = p
 110  continue
      i = i - 1
 120  IF(i .gt. 0) go to 100

END SUBROUTINE matinv


SUBROUTINE model(b,np,y,nob,x)
USE MOD_Precision
IMPLICIT NONE
integer, intent(in) :: np,nob
real(r8), intent(inout) :: b(np),x(nob)
real(r8), intent(out) :: y(nob)
integer i
      IF(b(3) .gt. 1.)      b(3) = 1.
      IF(b(3) .lt. 0.00001) b(3) = 0.00001
      IF(b(4) .gt. 10.) b(4) = 10.
      IF(b(4) .lt. 1.1) b(4) = 1.1
      DO i = 1,nob
         y(i)=b(1)+(b(2)-b(1))/(1+(b(3)*x(i))**b(4))**(1.-1/b(4))
      ENDDO
END SUBROUTINE model


SUBROUTINE SW_CB_dist ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydatc, nptf, phi, isiter)

!=================================================================
! DESCRIPTION:
! This is the subroutine for calculating the function/jacobian matrix
! of the distance between the fitted and prescribed SW retention curves
! for the Campbell model.
!
! Created by Nan Wei, 01/2019
! ----------------------------------------------------------------

      USE MOD_Precision
      IMPLICIT NONE

      integer m,n,ldfjac,iflag,i,nptf,isiter,npoint
      real(r8) x(n),fjac(ldfjac,n),fvec(m),xdat(npoint),ydatc(nptf,npoint),phi

      IF ( iflag == 0 ) THEN

         print*,x

      ELSEIF ( iflag == 1 ) THEN

         IF (x(1) >= 0.0) THEN
             isiter = 0
             RETURN
         ENDIF

         DO i = 1, m
            fvec(i) = sum(((-1.0*xdat(i)/x(1))**(-1.0*x(2)) * phi - ydatc(:,i))**2)
         ENDDO

      ELSEIF ( iflag == 2 ) THEN

         IF (x(1) >= 0.0) THEN
             isiter = 0
             RETURN
         ENDIF

         DO i = 1, m
            fjac(i,1) = sum(2.0*((-1.0*xdat(i)/x(1))**(-1.0*x(2)) * phi - ydatc(:,i))*&
                        phi * x(2) * (-1.0*xdat(i)/x(1))**(-1.0*x(2)) / x(1))
            fjac(i,2) = sum(-2.0*((-1.0*xdat(i)/x(1))**(-1.0*x(2)) * phi - ydatc(:,i))*&
                        phi * (-1.0*xdat(i)/x(1))**(-1.0*x(2)) * log(-1.0*xdat(i)/x(1)))
         ENDDO

      ENDIF

END SUBROUTINE SW_CB_dist

SUBROUTINE SW_VG_dist ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydatv, nptf, phi, isiter )

!=================================================================
! DESCRIPTION:
! This is the subroutine for calculating the function/jacobian matrix
! of the distance between the fitted and prescribed SW retention curves
! for the van Genuchten model.
!
! Created by Nan Wei, 01/2019
! ----------------------------------------------------------------

      USE MOD_Precision
      IMPLICIT NONE

      integer m,n,ldfjac,iflag,i,nptf,isiter,npoint
      real(r8) x(n),fjac(ldfjac,n),fvec(m),xdat(npoint),ydatv(nptf,npoint),phi

      IF ( iflag == 0 ) THEN

         print*,x

      ELSEIF ( iflag == 1 ) THEN

         IF (x(2) <= 0.0 .or. x(3) <= 0.1) THEN
             isiter = 0
             RETURN
         ENDIF

         DO i = 1, m
            fvec(i) = sum((x(1) + (phi - x(1))*(1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) - ydatv(:,i))**2)
         ENDDO

      ELSEIF ( iflag == 2 ) THEN

         IF (x(2) <= 0.0 .or. x(3) <= 0.1) THEN
             isiter = 0
             RETURN
         ENDIF

         DO i = 1, m
            fjac(i,1) = sum(2*(x(1) + (phi - x(1))*(1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) - ydatv(:,i))*&
                        (1 - (1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1)))
            fjac(i,2) = sum(2*(x(1) + (phi - x(1))*(1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) - ydatv(:,i))*&
      (phi - x(1)) * (1 - x(3)) * (1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-2) * x(2)**(x(3)-1) * xdat(i)**x(3))
            fjac(i,3) = sum(2*(x(1) + (phi - x(1))*(1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) - ydatv(:,i))*&
                        (phi - x(1)) * (1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) *&
                     ((1.0-x(3))*(x(2)*xdat(i))**x(3)*log(x(2)*xdat(i))/(x(3)*(1+(x(2)*xdat(i))**x(3))) &
                        - log(1+(x(2)*xdat(i))**x(3))/x(3)**2))
         ENDDO

      ENDIF

END SUBROUTINE SW_VG_dist
