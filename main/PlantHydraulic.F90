
MODULE PlantHydraulic

!-----------------------------------------------------------------------
  use precision
  IMPLICIT NONE
  SAVE

! PUBLIC MEMBER FUNCTIONS:
  public :: PlantHydraulicStress_oneleaf
  public :: PlantHydraulicStress_twoleaf
  public :: getvegwp_oneleaf
  public :: getvegwp_twoleaf

! PRIVATE MEMBER FUNCTIONS:
  private :: calcstress_oneleaf
  private :: calcstress_twoleaf


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------



  subroutine PlantHydraulicStress_oneleaf (nl_soil      ,nvegwcs    ,z_soi      ,&
                         dz_soi    ,rootfr    ,psrf       ,qsatl      ,qaf        ,tl         ,&
                         rb        ,ra        ,rd         ,rstfac     ,cint       ,lai        ,&
                         rhoair    ,fwet      ,sai        ,kmax_sun   ,kmax_sha   ,kmax_xyl   ,&
                         kmax_root ,psi50_sun ,psi50_sha  ,psi50_xyl  ,psi50_root ,htop       ,&
                         ck        ,smp       ,hk         ,hksati     ,vegwp      ,etr        ,&
                         rootr     ,sigf      ,qg         ,qm         ,gs0        ,k_soil_root,&
                         k_ax_root)  

!=======================================================================        
!                                                                               
!     calculation of plant hydraulic stress
!
!     Author: Xingjie Lu, 16/01/2019, modified from CLM5 plant_hydraulic_stress module
!
!----------------------------------------------------------------------         

 use precision
 IMPLICIT NONE                                                                        

 integer ,intent(in) :: nl_soil ! upper bound of array
 integer ,intent(in) :: nvegwcs ! upper bound of array
 real(r8),intent(in), dimension(nl_soil) :: &
      z_soi,        &! soil node depth (m)
      dz_soi         ! soil layer thicknesses (m)
 real(r8),intent(inout), dimension(nvegwcs) :: &
      vegwp          ! vegetation water potential
 real(r8),intent(inout):: &
      gs0            ! maximum stomata conductance

 real(r8),intent(in) :: &
      sigf           ! fraction of veg cover, excluding snow-covered veg [-]

 real(r8),intent(inout) :: &
      psrf,         &! surface atmospheric pressure (pa)
      qsatl,        &! specific humidity [kg/kg]
      qaf,          &! humidity of canopy air [kg/kg]
      qg,           &! specific humidity at ground surface [kg/kg]
      qm,           &! specific humidity at reference height [kg/kg]
      tl,           &! leaf temperature (K)

      rb,           &! boundary resistance from canopy to cas (s m-1)
      rd,           &! aerodynamical resistance between ground and canopy air
      ra             ! aerodynamic resistance from cas to refence height (s m-1)
 real(r8),intent(inout) :: &
      rstfac         ! canopy resistance stress factors to soil moisture 

 real(r8),intent(inout) :: &
      lai,          &! leaf area index, one-sided
      sai,          &! stem area index
      kmax_sun,     &
      kmax_sha,     &
      kmax_xyl,     &
      kmax_root,    &
      psi50_sun,    &! water potential at 50% loss of sunlit leaf tissue conductance (mmH2O)
      psi50_sha,    &! water potential at 50% loss of shaded leaf tissue conductance (mmH2O)
      psi50_xyl,    &! water potential at 50% loss of xylem tissue conductance (mmH2O)
      psi50_root,   &! water potential at 50% loss of root tissue conductance (mmH2O)
      htop,         &! canopy top [m]
      ck,           &! shape-fitting parameter for vulnerability curve (-)
      rhoair,       &! density [kg/m**3]
      fwet           ! fraction of foliage that is wet [-]

 real(r8),intent(in), dimension(3) :: &
      cint           ! scaling up from leaf to canopy

 real(r8),intent(in), dimension(nl_soil) :: &
      smp,      &    ! precipitation sensible heat from canopy
      rootfr,   &    ! root fraction
      hksati,   &    ! hydraulic conductivity at saturation [mm h2o/s] 
      hk             ! soil hydraulic conducatance [mm h2o/s]
      

 real(r8),intent(out) :: &! ATTENTION : all for canopy not leaf
      etr            ! transpiration (mm/s)
 
 real(r8),intent(out),dimension(nl_soil) :: &
      rootr          ! root water uptake from different layers

 real(r8),intent(inout),dimension(nl_soil) :: k_soil_root    ! radial root and soil conductance
 real(r8),intent(inout),dimension(nl_soil) :: k_ax_root      ! axial root conductance

!-------------------- local --------------------------------------------       

 integer, parameter :: iterationtotal = 6
                                                                               
 real(r8) c3,       &! c3 vegetation : 1; 0 for c4

      tprcor,       &! coefficient for unit transfer 
      gbh2o,        &! one side leaf boundary layer conductance of sunlit leaf (canopy scale:mol m-2 s-1)
      gb_mol         ! one side leaf boundary layer conductance of sunlit leaf (leaf scale:umol H2O m-2 s-1)

 real(r8), dimension(nl_soil) :: &
      fs     !root conductance scale factor (reduction in conductance due to decreasing (more negative) root water potential)
 real(r8), dimension(nl_soil) :: &
      rai    ! soil-root interface conductance [mm/s]

 real(r8)                :: soilflux     ! soil-root interface conductance [mm/s]
 real(r8)                :: soil_conductance ! soil conductance      
 real(r8)                :: root_conductance ! root conductance      
 real(r8)                :: r_soil ! root spacing [m]      
 real(r8)                :: root_biomass_density    ! root biomass density [g/m3]
 real(r8)                :: root_cross_sec_area     ! root cross sectional area [m2]
 real(r8)                :: root_length_density     ! root length density [m/m3]
 real(r8)                :: croot_average_length    ! average coarse root length [m]
 real(r8)                :: rs_resis                ! combined soil-root resistance [s]

 real(r8), parameter :: croot_lateral_length = 0.25_r8   ! specified lateral coarse root length [m]
 real(r8), parameter :: c_to_b               = 2.0_r8 !(g biomass /g C)
 real(r8), parameter :: rpi                  = 3.14159265358979_r8
 integer , parameter :: root                 = 4
 real(r8), parameter :: toldb                = 1.e-2_r8  ! tolerance for satisfactory bsun/bsha solution
 real(r8), parameter :: K_axs                = 2.0e-1

! temporary input
 real(r8), parameter :: froot_carbon = 288.392056287006_r8 
 real(r8), parameter :: root_radius  = 2.9e-4_r8
 real(r8), parameter :: root_density = 310000._r8
 real(r8), parameter :: froot_leaf   = 1.5_r8
 real(r8), parameter :: krmax        = 3.981071705534969e-009_r8

 real(r8),dimension(nvegwcs) :: x      ! vegetation water potential

 integer j

!----------------calculate root-soil interface conductance-----------------
do j = 1,nl_soil

! calculate conversion from conductivity to conductance
   root_biomass_density = c_to_b * froot_carbon * rootfr(j) / dz_soi(j)
! ensure minimum root biomass (using 1gC/m2)
   root_biomass_density = max(c_to_b*1._r8,root_biomass_density)

 ! Root length density: m root per m3 soil
   root_cross_sec_area = rpi*root_radius**2
   root_length_density = root_biomass_density / (root_density * root_cross_sec_area)

   ! Root-area index (RAI)
   rai(j) = (sai+lai) * froot_leaf * rootfr(j)

! fix coarse root_average_length to specified length
   croot_average_length = croot_lateral_length

! calculate r_soil using Gardner/spa equation (Bonan, GMD, 2014)
   r_soil = sqrt(1./(rpi*root_length_density))

   ! length scale approach
   soil_conductance = min(hksati(j),hk(j))/(1.e3*r_soil)

! use vegetation plc function to adjust root conductance
   fs(j)=  plc(smp(j),psi50_root,ck)

! krmax is root conductance per area per length
   root_conductance = (fs(j)*rai(j)*krmax)/(croot_average_length + z_soi(j))
   soil_conductance = max(soil_conductance, 1.e-16_r8)
   root_conductance = max(root_conductance, 1.e-16_r8)

! sum resistances in soil and root
   rs_resis = 1._r8/soil_conductance + 1._r8/root_conductance

! conductance is inverse resistance
! explicitly set conductance to zero for top soil layer
   if(rai(j)*rootfr(j) > 0._r8) then
      k_soil_root(j) =  1._r8/rs_resis
   else
      k_soil_root(j) =  0.
   end if
   k_ax_root(j) = (rootfr(j)/(dz_soi(j)*1000))*K_axs*0.6
end do
!=======================================================================        

      tprcor = 44.6*273.16*psrf/1.013e5

! one side leaf boundary layer conductance for water vapor [=1/(2*rb)]
! ATTENTION: rb in CLM is for one side leaf, but for SiB2 rb for 
! 2-side leaf, so the gbh2o shold be " 0.5/rb * tprcor/tl "
      gbh2o   = 1./rb * tprcor/tl                    ! mol m-2 s-1

      gb_mol = gbh2o / cint(3) * 1.e6  ! leaf to canopy

! rb is for single leaf, but here the flux is for canopy, thus
!      gbh2osun  = gbh2osun * cintsun(3)    ! debug by Xingjie Lu

      x = vegwp(1:nvegwcs)
      call calcstress_oneleaf(x, nvegwcs, rstfac, etr, rootr, gb_mol, gs0, &
             qsatl, qaf, qg, qm, rhoair, psrf, fwet, lai, sai, htop, tl, kmax_sun, &
             kmax_sha, kmax_xyl, kmax_root, psi50_sun, psi50_sha, psi50_xyl, psi50_root,&
             ck, nl_soil, z_soi, ra, rd, smp, k_soil_root, k_ax_root, sigf)
      vegwp(1:nvegwcs) = x

  end subroutine PlantHydraulicStress_oneleaf

  subroutine PlantHydraulicStress_twoleaf (nl_soil   ,nvegwcs   ,z_soi    ,&
                         dz_soi    ,rootfr    ,psrf       ,qsatlsun   ,qsatlsha   ,&
                         qaf       ,tlsun     ,tlsha      ,rbsun      ,rbsha      ,&
                         ra        ,rd        ,rstfacsun  ,rstfacsha  ,cintsun    ,&
                         cintsha   ,laisun    ,laisha     ,rhoair     ,fwet       ,&
                         sai       ,kmax_sun  ,kmax_sha   ,kmax_xyl   ,kmax_root  ,&
                         psi50_sun ,psi50_sha ,psi50_xyl  ,psi50_root ,htop       ,&
                         ck        ,smp       ,hk         ,hksati     ,vegwp      ,&
                         etrsun    ,etrsha    ,rootr      ,sigf       ,qg         ,&
                         qm        ,gs0sun    ,gs0sha     ,k_soil_root,k_ax_root  )  

!=======================================================================        
!                                                                               
!     calculation of plant hydraulic stress
!
!     Author: Xingjie Lu, 16/01/2019, modified from CLM5 plant_hydraulic_stress module
!
!----------------------------------------------------------------------         

 use precision
 IMPLICIT NONE                                                                        

 integer ,intent(in) :: nl_soil ! upper bound of array
 integer ,intent(in) :: nvegwcs ! upper bound of array
 real(r8),intent(in), dimension(nl_soil) :: &
      z_soi,      &! soil node depth (m)
      dz_soi       ! soil layer thicknesses (m)
 real(r8),intent(inout), dimension(nvegwcs) :: &
      vegwp       ! vegetation water potential
 real(r8),intent(inout):: &
      gs0sun,    & ! maximum stomata conductance of sunlit leaf
      gs0sha       ! maximum stomata conductance of shaded leaf

 real(r8),intent(in) :: &
      sigf,      & ! fraction of veg cover, excluding snow-covered veg [-]
      psrf,      & ! surface atmospheric pressure (pa)
      qg,           &! specific humidity at ground surface [kg/kg]
      qm             ! specific humidity at reference height [kg/kg]

 real(r8),intent(inout) :: &
      qsatlsun,     &! sunlit leaf specific humidity [kg/kg]
      qsatlsha,     &! shaded leaf specific humidity [kg/kg]
      qaf,          &! humidity of canopy air [kg/kg]
      tlsun,        &! sunlit leaf temperature (K)
      tlsha,        &! shaded leaf temperature (K)

      rbsun,        &! boundary resistance from sunlit canopy to cas (s m-1)
      rbsha,        &! boundary resistance from shaded canopy to cas (s m-1)
      rd,           &! aerodynamical resistance between ground and canopy air
      ra             ! aerodynamic resistance from cas to refence height (s m-1)
 real(r8),intent(inout) :: &
      rstfacsun,    &! canopy resistance stress factors to soil moisture for sunlit leaf                        
      rstfacsha      ! canopy resistance stress factors to soil moisture for shaded leaf                        

 real(r8),intent(in) :: &
      laisun,       &! sunlit leaf area index, one-sided
      laisha,       &! shaded leaf area index, one-sided
      sai,          &! stem area index
      kmax_sun,     &
      kmax_sha,     &
      kmax_xyl,     &
      kmax_root,    &
      psi50_sun,    &! water potential at 50% loss of sunlit leaf tissue conductance (mmH2O)
      psi50_sha,    &! water potential at 50% loss of shaded leaf tissue conductance (mmH2O)
      psi50_xyl,    &! water potential at 50% loss of xylem tissue conductance (mmH2O)
      psi50_root,   &! water potential at 50% loss of root tissue conductance (mmH2O)
      htop,         &! canopy top [m]
      ck,           &! shape-fitting parameter for vulnerability curve (-)
      rhoair,       &! density [kg/m**3]
      fwet           ! fraction of foliage that is wet [-]

 real(r8),intent(in), dimension(3) :: &
      cintsun,      &! scaling up from sunlit leaf to canopy
      cintsha        ! scaling up from shaded leaf to canopy

 real(r8),intent(in), dimension(nl_soil) :: &
      smp,      &    ! precipitation sensible heat from canopy
      rootfr,   &    ! root fraction
      hksati,   &    ! hydraulic conductivity at saturation [mm h2o/s] 
      hk             ! soil hydraulic conducatance [mm h2o/s]
      

 real(r8),intent(out) :: &! ATTENTION : all for canopy not leaf
      etrsun,       &! transpiration from sunlit leaf (mm/s)
      etrsha         ! transpiration from shaded leaf (mm/s)
 
 real(r8),intent(out),dimension(nl_soil) :: &
      rootr          ! root water uptake from different layers

 real(r8),intent(inout),dimension(nl_soil) :: k_soil_root    ! radial root and soil conductance
 real(r8),intent(inout),dimension(nl_soil) :: k_ax_root      ! axial root conductance


!-------------------- local --------------------------------------------       

 integer, parameter :: iterationtotal = 6
                                                                               
 real(r8) c3,       &! c3 vegetation : 1; 0 for c4

      tprcor,       &! coefficient for unit transfer 
      gbh2osun,     &! one side leaf boundary layer conductance of sunlit leaf (canopy scale:mol m-2 s-1)
      gbh2osha,     &! one side leaf boundary layer conductance of shaded leaf (canopy scale:mol m-2 s-1)
      gb_mol_sun,   &! one side leaf boundary layer conductance of sunlit leaf (leaf scale:umol H2O m-2 s-1)
      gb_mol_sha     ! one side leaf boundary layer conductance of shaded leaf (leaf scale:umol H2O m-2 s-1)

 real(r8), dimension(nl_soil) :: &
      fs     !root conductance scale factor (reduction in conductance due to decreasing (more negative) root water potential)
 real(r8), dimension(nl_soil) :: &
      rai    ! soil-root interface conductance [mm/s]

 real(r8)                :: soilflux     ! soil-root interface conductance [mm/s]
 real(r8)                :: soil_conductance ! soil conductance      
 real(r8)                :: root_conductance ! root conductance      
 real(r8)                :: r_soil ! root spacing [m]      
 real(r8)                :: root_biomass_density    ! root biomass density [g/m3]
 real(r8)                :: root_cross_sec_area     ! root cross sectional area [m2]
 real(r8)                :: root_length_density     ! root length density [m/m3]
 real(r8)                :: croot_average_length    ! average coarse root length [m]
 real(r8)                :: rs_resis                ! combined soil-root resistance [s]

 real(r8), parameter :: croot_lateral_length = 0.25_r8   ! specified lateral coarse root length [m]
 real(r8), parameter :: c_to_b               = 2.0_r8 !(g biomass /g C)
 real(r8), parameter :: rpi                  = 3.14159265358979_r8
 integer , parameter :: root                 = 4
 real(r8), parameter :: toldb                = 1.e-2_r8  ! tolerance for satisfactory bsun/bsha solution
 real(r8), parameter :: K_axs                = 2.0e-1

! temporary input
 real(r8), parameter :: froot_carbon = 288.392056287006_r8 
 real(r8), parameter :: root_radius  = 2.9e-4_r8
 real(r8), parameter :: root_density = 310000._r8
 real(r8), parameter :: froot_leaf   = 1.5_r8
 real(r8), parameter :: krmax        = 3.981071705534969e-009_r8

 real(r8),dimension(nvegwcs) :: x      ! vegetation water potential

 integer j

!----------------calculate root-soil interface conductance-----------------
do j = 1,nl_soil

! calculate conversion from conductivity to conductance
   root_biomass_density = c_to_b * froot_carbon * rootfr(j) / dz_soi(j)
! ensure minimum root biomass (using 1gC/m2)
   root_biomass_density = max(c_to_b*1._r8,root_biomass_density)

 ! Root length density: m root per m3 soil
   root_cross_sec_area = rpi*root_radius**2
   root_length_density = root_biomass_density / (root_density * root_cross_sec_area)

   ! Root-area index (RAI)
   rai(j) = (sai+laisun+laisha) * froot_leaf * rootfr(j)

! fix coarse root_average_length to specified length
   croot_average_length = croot_lateral_length

! calculate r_soil using Gardner/spa equation (Bonan, GMD, 2014)
   r_soil = sqrt(1./(rpi*root_length_density))

   ! length scale approach
   soil_conductance = min(hksati(j),hk(j))/(1.e3*r_soil)

! use vegetation plc function to adjust root conductance
   fs(j)=  plc(smp(j),psi50_root,ck)

! krmax is root conductance per area per length
   root_conductance = (fs(j)*rai(j)*krmax)/(croot_average_length + z_soi(j))
   soil_conductance = max(soil_conductance, 1.e-16_r8)
   root_conductance = max(root_conductance, 1.e-16_r8)

! sum resistances in soil and root
   rs_resis = 1._r8/soil_conductance + 1._r8/root_conductance

! conductance is inverse resistance
! explicitly set conductance to zero for top soil layer
   if(rai(j)*rootfr(j) > 0._r8) then
      k_soil_root(j) =  1._r8/rs_resis
   else
      k_soil_root(j) =  0.
   end if
   k_ax_root(j) = (rootfr(j)/(dz_soi(j)*1000))*K_axs*0.6
end do
!=======================================================================        

      tprcor = 44.6*273.16*psrf/1.013e5

! one side leaf boundary layer conductance for water vapor [=1/(2*rb)]
! ATTENTION: rb in CLM is for one side leaf, but for SiB2 rb for 
! 2-side leaf, so the gbh2o shold be " 0.5/rb * tprcor/tl "
      gbh2osun   = 1./rbsun * tprcor/tlsun                    ! mol m-2 s-1
      gbh2osha   = 1./rbsha * tprcor/tlsha                    ! mol m-2 s-1

      gb_mol_sun = gbh2osun / cintsun(3) * 1.e6  ! leaf to canopy
      gb_mol_sha = gbh2osha / cintsha(3) * 1.e6

! rb is for single leaf, but here the flux is for canopy, thus
!      gbh2osun  = gbh2osun * cintsun(3)    ! Commented by Xingjie Lu, 
!      gbh2osha  = gbh2osha * cintsha(3)

      x = vegwp(1:nvegwcs)

      call calcstress_twoleaf(x, nvegwcs, rstfacsun, rstfacsha, etrsun, etrsha, rootr,&
           gb_mol_sun, gb_mol_sha, gs0sun, gs0sha, qsatlsun, qsatlsha, qaf, qg, qm, rhoair, &
           psrf, fwet, laisun, laisha, sai, htop, tlsun, tlsha, kmax_sun, &
           kmax_sha, kmax_xyl, kmax_root, psi50_sun, psi50_sha, psi50_xyl, psi50_root, ck, &
           nl_soil, z_soi, ra, rd, smp, k_soil_root, k_ax_root, sigf)
   

      vegwp(1:nvegwcs) = x

  end subroutine PlantHydraulicStress_twoleaf

  subroutine calcstress_oneleaf(x,nvegwcs, rstfac, etr, rootr,&
             gb_mol, gs0, qsatl, qaf, qg, qm,rhoair, &
             psrf, fwet, lai, sai, htop, tl, kmax_sun, &
             kmax_sha, kmax_xyl, kmax_root, psi50_sun, psi50_sha, psi50_xyl, psi50_root, ck, & 
             nl_soil, z_soi, raw, rd, smp, k_soil_root, k_ax_root, sigf)
    !
    ! DESCRIPTIONS
    ! compute the transpiration stress using a plant hydraulics approach
    ! calls spacF, spacA, and getvegwp
    !
    ! !ARGUMENTS:
    integer                , intent(in)     :: nvegwcs
    real(r8)               , intent(inout)  :: x(nvegwcs)          ! working copy of vegwp(p,:)
    real(r8)               , intent(out)    :: rstfac              ! canopy transpiration wetness factor (0 to 1)
    real(r8)               , intent(out)    :: etr                 ! actual transpiration (mm/s)
    real(r8)               , intent(out)    :: rootr(nl_soil)        ! root water uptake from different layers
    integer                , intent(in)     :: nl_soil
    real(r8)               , intent(in)     :: z_soi(nl_soil)
    real(r8)               , intent(in)     :: gb_mol              ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)     :: gs0                 ! maximum shaded Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)     :: qsatl               ! leaf specific humidity [kg/kg]
    real(r8)               , intent(in)     :: qaf                 ! humidity of canopy air [kg/kg]
    real(r8)               , intent(in)     :: qg                  ! specific humidity at ground surface [kg/kg]
    real(r8)               , intent(in)     :: qm                  ! specific humidity at reference height [kg/kg]
    real(r8)               , intent(in)     :: rhoair              ! density [kg/m**3]
    real(r8)               , intent(in)     :: psrf                ! atmospheric pressure [Pa]
    real(r8)               , intent(in)     :: fwet                ! fraction of foliage that is green and dry [-]
    real(r8)               , intent(in)     :: raw                 ! moisture resistance [s/m]
    real(r8)               , intent(in)     :: rd                  ! aerodynamical resistance between ground and canopy air
    real(r8)               , intent(in)     :: lai                 ! leaf area index
    real(r8)               , intent(in)     :: sai                 ! stem area index
    real(r8)               , intent(in)     :: htop                ! canopy top [m] 
    real(r8)               , intent(in)     :: tl                  ! leaf temperature
    real(r8)               , intent(in)     :: kmax_sun
    real(r8)               , intent(in)     :: kmax_sha
    real(r8)               , intent(in)     :: kmax_xyl
    real(r8)               , intent(in)     :: kmax_root
    real(r8)               , intent(in)     :: psi50_sun              ! water potential at 50% loss of sunlit leaf tissue conductance (mmH2O)
    real(r8)               , intent(in)     :: psi50_sha              ! water potential at 50% loss of shaded leaf tissue conductance (mmH2O)
    real(r8)               , intent(in)     :: psi50_xyl              ! water potential at 50% loss of xylem tissue conductance (mmH2O)
    real(r8)               , intent(in)     :: psi50_root             ! water potential at 50% loss of root tissue conductance (mmH2O)
    real(r8)               , intent(in)     :: ck                     ! 
    real(r8)               , intent(in)     :: smp(nl_soil)           ! soil matrix potential
    real(r8)               , intent(in)     :: k_soil_root(nl_soil)   ! soil-root interface conductance [mm/s]
    real(r8)               , intent(in)     :: k_ax_root(nl_soil)     ! root axial-direction conductance [mm/s]
    real(r8)               , intent(in)     :: sigf                   ! fraction of veg cover, excluding snow-covered veg [-]

    real(r8) :: wtl                   ! heat conductance for leaf [m/s]
    real(r8) :: dx(nvegwcs)           ! change in vegwp from one iter to the next [mm]
    real(r8) :: efpot                 ! potential latent energy flux [kg/m2/s]
    real(r8) :: rppdry            ! fraction of potential evaporation through transp - shaded [-]
    real(r8) :: qflx              ! maximum transpiration without water stress [kg/m2/s]
    real(r8) :: gs         ! local gs_mol copies, actual stomata conductance
    real(r8) :: qeroot,dqeroot ! local gs_mol copies
    real(r8),dimension(nl_soil) :: xroot        ! local gs_mol copies
    integer  :: i,j                     ! index
    real(r8) :: cf                    ! s m**2/umol -> s/m
!    integer  :: iter,iterqflx         ! newton's method iteration number
!    logical  :: flag                  ! signal that matrix was not invertible
!    logical  :: night                 ! signal to store vegwp within this routine, b/c it is night-time and full suite won't be called
    integer, parameter :: itmax=50   ! exit newton's method if iters>itmax
!    real(r8),parameter :: toldx=1.e-9 !tolerances for a satisfactory solution
!    real(r8),parameter :: tolf         = 1.e-6_r8
!    real(r8),parameter :: tolf_leafxyl = 1.e-16_r8
!    real(r8),parameter :: tolf_root    = 1.e-14_r8 !tolerances for a satisfactory solution
!   logical  :: havegs                ! signals direction of calculation gs->qflx or qflx->gs 
!   logical  :: haroot                ! signals direction of calculation x_root_top->qeroot or qeroot->x_root_top
    real(r8) :: soilflux              ! total soil column transpiration [mm/s] 
    real(r8) :: x_root_top
    real(r8) :: maxscale
!    real(r8), parameter :: tol_lai=1.e-7_r8 ! minimum lai where transpiration is calc'd 
    integer, parameter :: leafsun=1
    integer, parameter :: leafsha=2
    integer, parameter :: xyl=3
    integer, parameter :: root=4
    !------------------------------------------------------------------------------
    

    !temporary flag for night time vegwp(sun)>0  

    gs=gs0
    call getqflx_gs2qflx_oneleaf(gb_mol,gs,qflx,qsatl,qaf,rhoair,psrf,lai,sai,fwet,tl,sigf,raw,rd,qg,qm)
    x_root_top  = x(root)
    if(qflx>0)then
       call getrootqflx_x2qe(nl_soil,smp,x_root_top ,z_soi,k_soil_root,k_ax_root,qeroot,dqeroot)

       call spacAF_oneleaf(x,nvegwcs,dx,nl_soil,qflx,lai,sai,htop,&
               qeroot,dqeroot,kmax_sun,kmax_sha,kmax_xyl,kmax_root,&
               psi50_sun,psi50_sha,psi50_xyl,psi50_root,ck)

       if ( maxval(abs(dx)) > 200000._r8) then
          maxscale = min(maxval(abs(dx)),maxval(abs(x))) / 2
          dx = maxscale * dx / maxval(abs(dx))! * log(maxval(abs(dx))/maxscale) !rescale step to max of 50000
       end if
    
       x=x+dx

    ! this is a catch to force spac gradient to atmosphere
       if ( x(xyl) > x(root) ) x(xyl) = x(root)
       if ( x(leafsun) > x(xyl) )  x(leafsun) = x(xyl)
       if ( x(leafsha) > x(xyl) )  x(leafsha) = x(xyl)

    ! compute attenuated flux
       etr=qflx*plc(x(leafsha),psi50_sha,ck)
    
    ! retrieve stressed stomatal conductance
       call getqflx_qflx2gs_oneleaf(gb_mol,gs,etr,qsatl,qaf,rhoair,psrf,lai,sai,fwet,tl,sigf,raw,rd,qg,qm)
    
    ! compute water stress
    ! .. generally -> B= gs_stressed / gs_unstressed
    ! .. when gs=0 -> B= plc( x )
       rstfac = amax1(gs/gs0,1.e-2_r8)
       call getrootqflx_qe2x(nl_soil,smp,z_soi,k_soil_root,k_ax_root,etr,xroot,x_root_top)

       x(root) = x_root_top
       do j = 1,nl_soil
          rootr(j) = k_soil_root(j)*(smp(j)-xroot(j))
       enddo
       soilflux = sum(rootr(:))
    else
       if ( x(xyl) > x(root) ) x(xyl) = x(root)
       if ( x(leafsun) > x(xyl) )  x(leafsun) = x(xyl)
       if ( x(leafsha) > x(xyl) )  x(leafsha) = x(xyl)
       etr = 0._r8
       rstfac = amax1(plc(x(leafsha),psi50_sha,ck),1.e-2_r8)
       gs = gs0 * rstfac
       rootr = 0._r8
    end if

    soilflux = sum(rootr(:))

  end subroutine calcstress_oneleaf

  subroutine calcstress_twoleaf(x,nvegwcs,rstfacsun, rstfacsha, etrsun, etrsha, rootr,&
             gb_mol_sun, gb_mol_sha, gs0sun, gs0sha, qsatlsun, qsatlsha, qaf, qg, qm,rhoair,&
             psrf, fwet, laisun, laisha, sai, htop, tlsun, tlsha, kmax_sun, kmax_sha, kmax_xyl, kmax_root, &
             psi50_sun, psi50_sha, psi50_xyl, psi50_root, ck, nl_soil, z_soi, raw, rd, smp, &
             k_soil_root, k_ax_root, sigf)
    !
    ! DESCRIPTIONS
    ! compute the transpiration stress using a plant hydraulics approach
    ! calls spacF, spacA, and getvegwp
    !
    ! !ARGUMENTS:
    integer                , intent(in)     :: nvegwcs
    real(r8)               , intent(inout)  :: x(nvegwcs)             ! working copy of vegwp(p,:)
    real(r8)               , intent(out)    :: rstfacsun              ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8)               , intent(out)    :: rstfacsha              ! shaded sunlit canopy transpiration wetness factor (0 to 1)
    real(r8)               , intent(out)    :: etrsun                 ! transpiration from sunlit leaf (mm/s)
    real(r8)               , intent(out)    :: etrsha                 ! transpiration from shaded leaf (mm/s)
    real(r8)               , intent(out)    :: rootr(nl_soil)        ! root water uptake from different layers
    integer                , intent(in)     :: nl_soil
    real(r8)               , intent(in)     :: z_soi(nl_soil)
    real(r8)               , intent(in)     :: gb_mol_sun             ! sunlit leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)     :: gb_mol_sha             ! shaded leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)     :: gs0sun                 ! sunlit Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)     :: gs0sha                 ! shaded Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(r8)               , intent(in)     :: qsatlsun               ! sunlit leaf specific humidity [kg/kg]
    real(r8)               , intent(in)     :: qsatlsha               ! shaded leaf specific humidity [kg/kg]
    real(r8)               , intent(in)     :: qaf                    ! humidity of canopy air [kg/kg]
    real(r8)               , intent(in)     :: qg                     ! specific humidity at ground surface [kg/kg]
    real(r8)               , intent(in)     :: qm                     ! specific humidity at reference height [kg/kg]
    real(r8)               , intent(in)     :: rhoair                 ! density [kg/m**3]
    real(r8)               , intent(in)     :: psrf                   ! atmospheric pressure [Pa]
    real(r8)               , intent(in)     :: fwet                   ! fraction of foliage that is green and dry [-]
    real(r8)               , intent(in)     :: raw                    ! moisture resistance [s/m]
    real(r8)               , intent(in)     :: rd                     ! aerodynamical resistance between ground and canopy air
    real(r8)               , intent(in)     :: laisun                 ! Sunlit leaf area index
    real(r8)               , intent(in)     :: laisha                 ! Shaded leaf area index
    real(r8)               , intent(in)     :: sai                    ! stem area index
    real(r8)               , intent(in)     :: htop                   ! canopy top [m] 
    real(r8)               , intent(in)     :: tlsun                ! sunlit leaf temperature
    real(r8)               , intent(in)     :: tlsha                ! shaded leaf temperature
    real(r8)               , intent(in)     :: kmax_sun
    real(r8)               , intent(in)     :: kmax_sha
    real(r8)               , intent(in)     :: kmax_xyl
    real(r8)               , intent(in)     :: kmax_root
    real(r8)               , intent(in)     :: psi50_sun              ! water potential at 50% loss of sunlit leaf tissue conductance (mmH2O)
    real(r8)               , intent(in)     :: psi50_sha              ! water potential at 50% loss of shaded leaf tissue conductance (mmH2O)
    real(r8)               , intent(in)     :: psi50_xyl              ! water potential at 50% loss of xylem tissue conductance (mmH2O)
    real(r8)               , intent(in)     :: psi50_root             ! water potential at 50% loss of root tissue conductance (mmH2O)
    real(r8)               , intent(in)     :: ck                     ! 
    real(r8)               , intent(in)     :: smp(nl_soil)           ! soil matrix potential
    real(r8)               , intent(in)     :: k_soil_root(nl_soil)   ! soil-root interface conductance [mm/s]
    real(r8)               , intent(in)     :: k_ax_root(nl_soil)     ! root axial-direction conductance [mm/s]
    real(r8)               , intent(in)     :: sigf                   ! fraction of veg cover, excluding snow-covered veg [-]
    real(r8) :: wtl                   ! heat conductance for leaf [m/s]
    real(r8) :: A(nvegwcs,nvegwcs)    ! matrix relating d(vegwp) and f: d(vegwp)=A*f 
    real(r8) :: f(nvegwcs)            ! flux divergence (mm/s)
    real(r8) :: dx(nvegwcs)           ! change in vegwp from one iter to the next [mm]
    real(r8) :: efpot                 ! potential latent energy flux [kg/m2/s]
    real(r8) :: rppdry_sun            ! fraction of potential evaporation through transp - sunlit [-]
    real(r8) :: rppdry_sha            ! fraction of potential evaporation through transp - shaded [-]
    real(r8) :: qflx_sun              ! [kg/m2/s]
    real(r8) :: qflx_sha              ! [kg/m2/s]
    real(r8) :: gssun,gssha         ! local gs_mol copies
    real(r8) :: qeroot,dqeroot
    real(r8),dimension(nl_soil) :: xroot        ! local gs_mol copies
    integer  :: i,j                     ! index
    real(r8) :: cf                    ! s m**2/umol -> s/m
    integer  :: iter,iterqflx         ! newton's method iteration number
    logical  :: flag                  ! signal that matrix was not invertible
    logical  :: night                 ! signal to store vegwp within this routine, b/c it is night-time and full suite won't be called
    integer, parameter :: itmax=50   ! exit newton's method if iters>itmax
    real(r8),parameter :: toldx=1.e-9 !tolerances for a satisfactory solution
    real(r8),parameter :: tolf         = 1.e-6_r8
    real(r8),parameter :: tolf_leafxyl = 1.e-16_r8
    real(r8),parameter :: tolf_root    = 1.e-14_r8 !tolerances for a satisfactory solution
    logical  :: havegs                ! signals direction of calculation gs->qflx or qflx->gs 
    logical  :: haroot                ! signals direction of calculation x_root_top->qeroot or qeroot->x_root_top
    real(r8) :: soilflux              ! total soil column transpiration [mm/s] 
    real(r8) :: x_root_top
    real(r8) :: x_root_top1
    real(r8) :: x_root_top2
    real(r8) :: dxsoiltop
    real(r8) :: maxscale
    real(r8), parameter :: tol_lai=1.e-7_r8 ! minimum lai where transpiration is calc'd 
    integer, parameter :: leafsun=1
    integer, parameter :: leafsha=2
    integer, parameter :: xyl=3
    integer, parameter :: root=4
    real(r8) fsto1,fsto2,fx,fr,grav1
    !------------------------------------------------------------------------------
    

    !temporary flag for night time vegwp(sun)>0  

    gssun=gs0sun
    gssha=gs0sha
    call getqflx_gs2qflx_twoleaf(gb_mol_sun,gb_mol_sha,gssun,gssha,qflx_sun,qflx_sha,qsatlsun,qsatlsha,qaf, &
                                 rhoair,psrf,laisun,laisha,sai,fwet,tlsun,tlsha,sigf,raw,rd,qg,qm)
    x_root_top  = x(root)

    if(qflx_sun .gt. 0 .or. qflx_sha .gt. 0)then
       call getrootqflx_x2qe(nl_soil,smp,x_root_top ,z_soi,k_soil_root,k_ax_root,qeroot,dqeroot)

       call spacAF_twoleaf(x,nvegwcs,dx,nl_soil,qflx_sun,qflx_sha,laisun,laisha,sai,htop,&
               qeroot,dqeroot,kmax_sun,kmax_sha,kmax_xyl,kmax_root,&
               psi50_sun,psi50_sha,psi50_xyl,psi50_root,ck)

       if ( maxval(abs(dx)) > 200000._r8) then
          maxscale = min(maxval(abs(dx)),maxval(abs(x))) / 2
          dx = maxscale * dx / maxval(abs(dx))! * log(maxval(abs(dx))/maxscale) !rescale step to max of 50000
       end if

       x=x+dx

    ! this is a catch to force spac gradient to atmosphere
       if ( x(xyl) > x(root) ) x(xyl) = x(root)
       if ( x(leafsun) > x(xyl) )  x(leafsun) = x(xyl)
       if ( x(leafsha) > x(xyl) )  x(leafsha) = x(xyl)

    ! compute attenuated flux; the actual transpiration
       etrsun=qflx_sun*plc(x(leafsun),psi50_sun,ck)
       etrsha=qflx_sha*plc(x(leafsha),psi50_sha,ck)
    
    ! retrieve stressed stomatal conductance
       call getqflx_qflx2gs_twoleaf(gb_mol_sun,gb_mol_sha,gssun,gssha,etrsun,etrsha,qsatlsun,qsatlsha,qaf, &
                                    rhoair,psrf,laisun,laisha,sai,fwet,tlsun,tlsha,sigf,raw,rd,qg,qm)
    
    ! compute water stress
    ! .. generally -> B= gs_stressed / gs_unstressed
    ! .. when gs=0 -> B= plc( x )
          rstfacsun = amax1(gssun/gs0sun,1.e-2_r8)
          rstfacsha = amax1(gssha/gs0sha,1.e-2_r8)
       qeroot = etrsun + etrsha
       call getrootqflx_qe2x(nl_soil,smp,z_soi,k_soil_root,k_ax_root,qeroot,xroot,x_root_top)
       x(root) = x_root_top
       do j = 1,nl_soil
          rootr(j) = k_soil_root(j)*(smp(j)-xroot(j))
       enddo
    else
       if ( x(xyl) > x(root) ) x(xyl) = x(root)
       if ( x(leafsun) > x(xyl) )  x(leafsun) = x(xyl)
       if ( x(leafsha) > x(xyl) )  x(leafsha) = x(xyl)
       etrsun = 0._r8
       etrsha = 0._r8
       rstfacsun = amax1(plc(x(leafsun),psi50_sun,ck),1.e-2_r8)
       rstfacsha = amax1(plc(x(leafsha),psi50_sha,ck),1.e-2_r8)
       gssun = gs0sun * rstfacsun
       gssha = gs0sha * rstfacsha
       rootr = 0._r8
    end if

    soilflux = sum(rootr(:))

  end subroutine calcstress_twoleaf
   
   !------------------------------------------------------------------------------
   
  subroutine spacAF_oneleaf(x,nvegwcs,dx,nl_soil,qflx,lai,sai,htop,&
                   qeroot,dqeroot,kmax_sun,kmax_sha,kmax_xyl,kmax_root,&
                   psi50_sun,psi50_sha,psi50_xyl,psi50_root,ck)
    !
    ! DESCRIPTION
    !  Returns invA, the inverse matrix relating delta(vegwp) to f
    !   d(vegwp)=invA*f
    !   evaluated at vegwp(p)
    !
    ! The methodology is currently hardcoded for linear algebra assuming the
    ! number of vegetation segments is four. Thus the matrix A and it's inverse
    ! invA are both 3x3 matrices. A more general method could be done using for
    ! example a LINPACK linear algebra solver.
    !
    ! !ARGUMENTS:
    integer                , intent(in)  :: nvegwcs
    real(r8)               , intent(in)  :: x(nvegwcs)      ! working copy of veg water potential for patch p [mm H2O] 
    real(r8)               , intent(out) :: dx(nvegwcs)   ! matrix relating d(vegwp) and f: d(vegwp)=invA*f
    integer                , intent(in)  :: nl_soil
    real(r8)               , intent(in)  :: qflx            ! leaf transpiration [kg/m2/s]
    real(r8)               , intent(in)  :: lai                    ! Shaded leaf area index
    real(r8)               , intent(in)  :: sai                    ! Stem area index
    real(r8)               , intent(in)  :: htop                   ! Canopy top [m]
    real(r8)               , intent(in)  :: qeroot  ! soil-root interface conductance [mm/s]
    real(r8)               , intent(in)  :: dqeroot  ! soil-root interface conductance [mm/s]
    real(r8)               , intent(in)  :: kmax_sun
    real(r8)               , intent(in)  :: kmax_sha
    real(r8)               , intent(in)  :: kmax_xyl
    real(r8)               , intent(in)  :: kmax_root
    real(r8)               , intent(in)  :: psi50_sun        ! water potential at 50% loss of sunlit leaf tissue conductance (mmH2O)
    real(r8)               , intent(in)  :: psi50_sha        ! water potential at 50% loss of shaded leaf tissue conductance (mmH2O)
    real(r8)               , intent(in)  :: psi50_xyl        ! water potential at 50% loss of xylem tissue conductance (mmH2O)
    real(r8)               , intent(in)  :: psi50_root       ! water potential at 50% loss of root tissue conductance (mmH2O)
    real(r8)               , intent(in)  :: ck
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtl                   ! heat conductance for leaf [m/s]
    real(r8) :: fsto                  ! transpiration reduction function [-] 
    real(r8) :: fx                    ! fraction of maximum conductance, xylem-to-leaf [-] 
    real(r8) :: fr                    ! fraction of maximum conductance, root-to-xylem [-] 
    real(r8) :: dfsto                 ! 1st derivative of fsto w.r.t. change in vegwp
    real(r8) :: dfx                   ! 1st derivative of fx w.r.t. change in vegwp
    real(r8) :: dfr                   ! 1st derivative of fr w.r.t. change in vegwp
    real(r8) :: A11, A12, A21, A22, A23, A32, A33  ! matrix relating vegwp to flux divergence f=A*d(vegwp)
    real(r8) :: leading               ! inverse of determiniant
    real(r8) :: determ                ! determinant of matrix
    real(r8) :: grav1                 ! gravitational potential surface to canopy top (mm H2O)
    real(r8) :: f(nvegwcs)
    real(r8), parameter :: tol_lai=1.e-7_r8 ! minimum lai where transpiration is calc'd
    integer, parameter :: leafsun=1
    integer, parameter :: leafsha=2
    integer, parameter :: xyl=3
    integer, parameter :: root=4
    integer  :: j                     ! index
    !------------------------------------------------------------------------------
    
    grav1 = htop*1000._r8

    !compute conductance attentuation for each segment
    fsto=   plc(x(leafsha),psi50_sha,ck)
    fx=     plc(x(xyl),psi50_xyl,ck)
    fr=     plc(x(root),psi50_root,ck)
    
    !compute 1st deriv of conductance attenuation for each segment
    dfsto=   d1plc(x(leafsha),psi50_sha,ck)
    dfx=     d1plc(x(xyl),psi50_xyl,ck)
    dfr=     d1plc(x(root),psi50_root,ck)

    A11 = - lai * kmax_sha * fx&
          - qflx * dfsto
    A12 =   lai * kmax_sha * dfx * (x(xyl)-x(leafsha))&
          + lai * kmax_sha * fx
    A21 =   lai * kmax_sha * fx
    A22 = - lai * kmax_sha * dfx * (x(xyl)-x(leafsha)) - lai * kmax_sha * fx&
          - sai * kmax_xyl / htop * fr
    A23 =   sai * kmax_xyl / htop * dfr * (x(root)-x(xyl)-grav1)&
          + sai * kmax_xyl / htop * fr
    A32 =   sai * kmax_xyl / htop * fr
    A33 = - sai * kmax_xyl / htop * fr&
          - sai * kmax_xyl / htop * dfr * (x(root)-x(xyl)-grav1)&
          + dqeroot

    f(leafsha)  = qflx * fsto - lai * kmax_sha * fx * (x(xyl)-x(leafsha))
    f(xyl)  = lai * kmax_sha * fx * (x(xyl)-x(leafsha)) &
            - sai * kmax_xyl / htop * fr * (x(root)-x(xyl)-grav1)
    f(root) = sai * kmax_xyl / htop * fr * (x(root)-x(xyl)-grav1) - qeroot
    determ=A11*A22*A33-A23*A11*A32-A12*A21*A33

    if(determ .ne. 0)then
       dx(leafsha) = (- A12*A33*f(xyl) + A12*A23*f(root) + (A22*A33 - A23*A32)*f(leafsha)) / determ
       dx(xyl)     = (  A11*A33*f(xyl) - A11*A23*f(root) - A21*A33*f(leafsha))             / determ
       dx(root)    = (- A11*A32*f(xyl) + (A11*A22 - A12*A21)*f(root) + A21*A32*f(leafsha)) / determ

       dx(leafsun) = x(leafsha) - x(leafsun) + dx(leafsha)
    else
       dx = 0._r8
    end if

  end subroutine spacAF_oneleaf
  
  !------------------------------------------------------------------------------
  subroutine spacAF_twoleaf(x,nvegwcs,dx,nl_soil,qflx_sun,qflx_sha,laisun,laisha,sai,htop,&
                   qeroot,dqeroot,kmax_sun,kmax_sha,kmax_xyl,kmax_root,&
                   psi50_sun,psi50_sha,psi50_xyl,psi50_root,ck)
    !
    ! DESCRIPTION
    !  Returns invA, the inverse matrix relating delta(vegwp) to f
    !   d(vegwp)=invA*f
    !   evaluated at vegwp(p)
    !
    ! The methodology is currently hardcoded for linear algebra assuming the
    ! number of vegetation segments is four. Thus the matrix A and it's inverse
    ! invA are both 4x4 matrices. A more general method could be done using for
    ! example a LINPACK linear algebra solver.
    !
    ! !ARGUMENTS:
    integer                , intent(in)  :: nvegwcs
    real(r8)               , intent(in)  :: x(nvegwcs)      ! working copy of veg water potential for patch p [mm H2O] 
    real(r8)               , intent(out) :: dx(nvegwcs)   ! matrix relating d(vegwp) and f: d(vegwp)=invA*f
    integer                , intent(in)  :: nl_soil
    real(r8)               , intent(in)  :: qflx_sun        ! Sunlit leaf transpiration [kg/m2/s] 
    real(r8)               , intent(in)  :: qflx_sha        ! Shaded leaf transpiration [kg/m2/s]
    real(r8)               , intent(in)  :: laisun                 ! Sunlit leaf area index
    real(r8)               , intent(in)  :: laisha                 ! Shaded leaf area index
    real(r8)               , intent(in)  :: sai                    ! Stem area index
    real(r8)               , intent(in)  :: htop                   ! Canopy top [m]
    real(r8)               , intent(in)  :: qeroot  ! soil-root interface conductance [mm/s]
    real(r8)               , intent(in)  :: dqeroot  ! soil-root interface conductance [mm/s]
    real(r8)               , intent(in)  :: kmax_sun
    real(r8)               , intent(in)  :: kmax_sha
    real(r8)               , intent(in)  :: kmax_xyl
    real(r8)               , intent(in)  :: kmax_root
    real(r8)               , intent(in)  :: psi50_sun        ! water potential at 50% loss of sunlit leaf tissue conductance (mmH2O)
    real(r8)               , intent(in)  :: psi50_sha        ! water potential at 50% loss of shaded leaf tissue conductance (mmH2O)
    real(r8)               , intent(in)  :: psi50_xyl        ! water potential at 50% loss of xylem tissue conductance (mmH2O)
    real(r8)               , intent(in)  :: psi50_root       ! water potential at 50% loss of root tissue conductance (mmH2O)
    real(r8)               , intent(in)  :: ck
    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtl                   ! heat conductance for leaf [m/s]
    real(r8) :: fsto1                 ! sunlit transpiration reduction function [-]
    real(r8) :: fsto2                 ! shaded transpiration reduction function [-] 
    real(r8) :: fx                    ! fraction of maximum conductance, xylem-to-leaf [-] 
    real(r8) :: fr                    ! fraction of maximum conductance, root-to-xylem [-] 
    real(r8) :: dfsto1                ! 1st derivative of fsto1 w.r.t. change in vegwp
    real(r8) :: dfsto2                ! 1st derivative of fsto2 w.r.t. change in vegwp
    real(r8) :: dfx                   ! 1st derivative of fx w.r.t. change in vegwp
    real(r8) :: dfr                   ! 1st derivative of fr w.r.t. change in vegwp
    real(r8) :: A11, A13, A22, A23, A31, A32, A33, A34, A43, A44   ! matrix relating vegwp to flux divergence f=A*d(vegwp)
    real(r8) :: leading               ! inverse of determiniant
    real(r8) :: determ                ! determinant of matrix
    real(r8) :: grav1                 ! gravitational potential surface to canopy top (mm H2O)
    real(r8) :: invfactor             ! 
    real(r8) :: f(nvegwcs)
    real(r8), parameter :: tol_lai=1.e-7_r8 ! minimum lai where transpiration is calc'd
    integer, parameter :: leafsun=1
    integer, parameter :: leafsha=2
    integer, parameter :: xyl=3
    integer, parameter :: root=4
    integer  :: j                     ! index
    !------------------------------------------------------------------------------
    
    grav1 = htop*1000._r8
    
    !compute conductance attentuation for each segment
    fsto1=  plc(x(leafsun),psi50_sun,ck)
    fsto2=  plc(x(leafsha),psi50_sha,ck)
    fx=     plc(x(xyl),psi50_xyl,ck)
    fr=     plc(x(root),psi50_root,ck)
    
    !compute 1st deriv of conductance attenuation for each segment
    dfsto1=  d1plc(x(leafsun),psi50_sun,ck)
    dfsto2=  d1plc(x(leafsha),psi50_sha,ck)
    dfx=     d1plc(x(xyl),psi50_xyl,ck)
    dfr=     d1plc(x(root),psi50_root,ck)

    
    A11= - laisun * kmax_sun * fx&
            - qflx_sun * dfsto1
    A13= laisun * kmax_sun * dfx * (x(xyl)-x(leafsun))&
            + laisun * kmax_sun * fx
    A22= - laisha * kmax_sha * fx&
         - qflx_sha * dfsto2
    A23= laisha * kmax_sha * dfx * (x(xyl)-x(leafsha))&
         + laisha * kmax_sha * fx
    A31= laisun * kmax_sun * fx
    A32= laisha * kmax_sha * fx
    A33= - laisun * kmax_sun * dfx * (x(xyl)-x(leafsun)) - laisun * kmax_sun * fx&
            - laisha * kmax_sha * dfx * (x(xyl)-x(leafsha)) - laisha * kmax_sha * fx&
            - sai * kmax_xyl / htop * fr
    A34= sai * kmax_xyl / htop * dfr * (x(root)-x(xyl)-grav1)&
            + sai * kmax_xyl / htop * fr
    A43= sai * kmax_xyl / htop * fr
    A44= - sai * kmax_xyl / htop * fr&
            - sai * kmax_xyl / htop * dfr * (x(root)-x(xyl)-grav1)&
            + dqeroot

    !compute flux divergence across each plant segment
    f(leafsun)  = qflx_sun * fsto1 - laisun * kmax_sun * fx * (x(xyl)-x(leafsun))
    f(leafsha)  = qflx_sha * fsto2 - laisha * kmax_sha * fx * (x(xyl)-x(leafsha))
    f(xyl)  = laisun * kmax_sun * fx * (x(xyl)-x(leafsun))&
            + laisha * kmax_sha * fx * (x(xyl)-x(leafsha)) &
            - sai * kmax_xyl / htop * fr * (x(root)-x(xyl)-grav1)
    f(root) = sai * kmax_xyl / htop * fr * (x(root)-x(xyl)-grav1) - qeroot

    if(qflx_sha > 0 )then
       determ=A44*A22*A33*A11-A44*A22*A31*A13-A44*A32*A23*A11-A43*A11*A22*A34

       if(determ .ne. 0)then
          dx(leafsun) = ((A22*A33*A44 - A22*A34*A43 - A23*A32*A44)*f(leafsun) + A13*A32*A44*f(leafsha) &
                        - A13*A22*A44*f(xyl) + A13*A22*A34*f(root)) / determ
          dx(leafsha) = ( A23*A31*A44*f(leafsun) + (A11*A33*A44 - A11*A34*A43 - A13*A31*A44)*f(leafsha) &
                        - A11*A23*A44*f(xyl) + A11*A23*A34*f(root)) / determ
          dx(xyl)     = (-A22*A31*A44*f(leafsun) - A11*A32*A44*f(leafsha) &
                        + A11*A22*A44*f(xyl) - A11*A22*A34*f(root)) / determ
          dx(root)    = ( A22*A31*A43*f(leafsun) + A11*A32*A43*f(leafsha) &
                        - A11*A22*A43*f(xyl) +(A11*A22*A33 - A11*A23*A32 - A13*A22*A31)*f(root)) / determ
       else
          dx = 0._r8
       end if
    else
       A33 = - laisun * kmax_sun * dfx * (x(xyl)-x(leafsun)) - laisun * kmax_sun * fx - sai * kmax_xyl / htop * fr
       f(xyl) = laisun * kmax_sun * fx * (x(xyl)-x(leafsun)) - sai * kmax_xyl / htop * fr * (x(root)-x(xyl)-grav1)
       determ=A11*A33*A44-A34*A11*A43-A13*A31*A44
       if(determ .ne. 0)then
          dx(leafsun) = (- A13*A44*f(xyl) + A13*A34*f(root) + (A33*A44 - A34*A43)*f(leafsun)) / determ
          dx(xyl)     = (  A11*A44*f(xyl) - A11*A34*f(root) - A31*A44*f(leafsun))             / determ
          dx(root)    = (- A11*A43*f(xyl) + (A11*A33 - A13*A31)*f(root) + A31*A43*f(leafsun)) / determ

          dx(leafsha) = x(leafsun) - x(leafsha) + dx(leafsun)
       else
          dx = 0._r8
       end if
    end if
     
  end subroutine spacAF_twoleaf
  
  !--------------------------------------------------------------------------------
  
  !--------------------------------------------------------------------------------
  subroutine getvegwp_oneleaf(x, nvegwcs, nl_soil, z_soi, gb_mol, gs_mol, &
             qsatl, qaf,qg,qm,rhoair, psrf, fwet, lai, htop, sai, tl, sigf,&
             raw, rd, smp, k_soil_root, k_ax_root, kmax_xyl, kmax_root, &
             psi50_sun, psi50_sha, psi50_xyl, psi50_root, ck, soilflux, rootr, etr)
    ! !DESCRIPTION:
    !  Calculates transpiration and returns corresponding vegwp in x
    !
    ! !USES:
    ! calls getqflx
  use PhysicalConstants, only : tfrz
    implicit none
    !
    ! !ARGUMENTS:
    integer       , intent(in)     :: nvegwcs
    real(r8)      , intent(out)    :: x(nvegwcs)             ! working copy of veg water potential for patch p
    integer       , intent(in)     :: nl_soil                ! number of soil layers
    real(r8)      , intent(in)     :: z_soi(nl_soil)         ! node depth [m]
    real(r8)      , intent(in)     :: gb_mol             ! Shaded leaf boundary layer conductance [umol H2O/m**2/s]
    real(r8)      , intent(inout)  :: gs_mol             ! Ball-Berry leaf conductance [umol H2O/m**2/s]
    real(r8)      , intent(in)     :: qsatl               ! Shalit leaf specific humidity [kg/kg]
    real(r8)      , intent(in)     :: qaf                    ! humidity of canopy air [kg/kg]
    real(r8)      , intent(in)     :: qg                     ! specific humidity at ground surface [kg/kg]
    real(r8)      , intent(in)     :: qm                     ! specific humidity at reference height [kg/kg]
    real(r8)      , intent(in)     :: rhoair                 ! density [kg/m**3]
    real(r8)      , intent(in)     :: psrf                   ! atmospheric pressure [Pa]
    real(r8)      , intent(in)     :: fwet                   ! fraction of foliage that is green and dry [-]
    real(r8)      , intent(in)     :: lai                 ! Shaded leaf area index
    real(r8)      , intent(in)     :: htop                   ! canopy top [m] 
    real(r8)      , intent(in)     :: sai                    ! stem area index
    real(r8)      , intent(in)     :: tl                ! shaded leaf temperature
    real(r8)      , intent(in)     :: sigf                   ! fraction of veg cover, excluding snow-covered veg [-]
    real(r8)      , intent(in)     :: kmax_xyl
    real(r8)      , intent(in)     :: kmax_root
    real(r8)      , intent(in)     :: psi50_sun              ! water potential at 50% loss of sunlit leaf tissue conductance (mmH2O)
    real(r8)      , intent(in)     :: psi50_sha              ! water potential at 50% loss of shaded leaf tissue conductance (mmH2O)
    real(r8)      , intent(in)     :: psi50_xyl              ! water potential at 50% loss of xylem tissue conductance (mmH2O)
    real(r8)      , intent(in)     :: psi50_root             ! water potential at 50% loss of root tissue conductance (mmH2O)
    real(r8)      , intent(in)     :: ck                     ! 
    real(r8)      , intent(in)     :: raw        ! moisture resistance [s/m]
    real(r8)      , intent(in)     :: rd         ! aerodynamical resistance between ground and canopy air
    real(r8)      , intent(in)     :: smp(nl_soil)           ! soil matrix potential
    real(r8)      , intent(in)     :: k_soil_root(nl_soil)   ! soil-root interface conductance [mm/s]
    real(r8)      , intent(in)     :: k_ax_root(nl_soil)     ! root axial-direction conductance [mm/s]
    real(r8)      , intent(out)    :: soilflux               ! total soil column transpiration [mm/s]
    real(r8)      , intent(out)    :: etr      ! transpiration from shaded leaf (mm/s)
    real(r8)      , intent(out)    :: rootr(nl_soil)        ! root water uptake from different layers
    !
    ! !LOCAL VARIABLES:
    real(r8) :: qeroot
    real(r8) :: dummy
    real(r8) :: fx                       ! fraction of maximum conductance, xylem-to-leaf [-]  
    real(r8) :: fr                       ! fraction of maximum conductance, root-to-xylem [-]  
    real(r8) :: x_root_top
    real(r8) :: xroot(nl_soil)
    real(r8) :: grav1                    ! gravitational potential surface to canopy top (mm H2O)
    real(r8) :: grav2(nl_soil)           ! soil layer gravitational potential relative to surface (mm H2O) 
    integer  :: j                        ! index
    logical  :: havegs                   ! signals direction of calculation gs->qflx or qflx->gs 
    logical  :: haroot                   ! signals direction of calculation x_root_top->qeroot or qeroot->x_root_top
    integer, parameter :: leafsun=1
    integer, parameter :: leafsha=2
    integer, parameter :: xyl=3
    integer, parameter :: root=4
    !----------------------------------------------------------------------
    grav1 = 1000._r8 * htop
    grav2(1:nl_soil) = 1000._r8 * z_soi(1:nl_soil)
    
    !compute transpiration demand
    havegs=.true.
    call getqflx_gs2qflx_oneleaf(gb_mol,gs_mol,etr,qsatl,qaf, &
                 rhoair,psrf,lai,sai,fwet,tl,sigf,raw,rd,qg,qm)

    !calculate root water potential
    call getrootqflx_qe2x(nl_soil,smp,z_soi,k_soil_root,k_ax_root,etr,xroot,x_root_top)

    !calculate xylem water potential
    fr = plc(x(root),psi50_root,ck)
    x(xyl) = x(root) - grav1 - etr/(fr*kmax_root/htop*sai)
    
    !calculate sun/sha leaf water potential
    fx = plc(x(xyl),psi50_xyl,ck)
    x(leafsha) = x(xyl) - etr/(fx*kmax_xyl*lai)
    x(leafsun) = x(xyl)

    !calculate soil flux
    do j = 1,nl_soil
       rootr(j) = k_soil_root(j)*(smp(j)-xroot(j))
    enddo
    soilflux = sum(rootr(:))

  end subroutine getvegwp_oneleaf
  
  !--------------------------------------------------------------------------------
  subroutine getvegwp_twoleaf(x, nvegwcs, nl_soil, z_soi, gb_mol_sun, gb_mol_sha, gs_mol_sun, gs_mol_sha, &
             qsatlsun, qsatlsha, qaf,qg,qm,rhoair, psrf, fwet, laisun, laisha, htop, sai, tlsun, tlsha, sigf,&
             raw, rd, smp, k_soil_root, k_ax_root, kmax_xyl, kmax_root, &
             psi50_sun, psi50_sha, psi50_xyl, psi50_root, ck, soilflux, rootr, etrsun, etrsha)
    ! !DESCRIPTION:
    !  Calculates transpiration and returns corresponding vegwp in x
    !
    ! !USES:
    ! calls getqflx
  use PhysicalConstants, only : tfrz
    implicit none
    !
    ! !ARGUMENTS:
    integer       , intent(in)     :: nvegwcs
    real(r8)      , intent(out)    :: x(nvegwcs)             ! working copy of veg water potential for patch p
    integer       , intent(in)     :: nl_soil                ! number of soil layers
    real(r8)      , intent(in)     :: z_soi(nl_soil)         ! node depth [m]
    real(r8)      , intent(in)     :: gb_mol_sun             ! Sunlit leaf boundary layer conductance [umol H2O/m**2/s]
    real(r8)      , intent(in)     :: gb_mol_sha             ! Shaded leaf boundary layer conductance [umol H2O/m**2/s]
    real(r8)      , intent(inout)  :: gs_mol_sun             ! Ball-Berry leaf conductance [umol H2O/m**2/s]
    real(r8)      , intent(inout)  :: gs_mol_sha             ! Ball-Berry leaf conductance [umol H2O/m**2/s]
    real(r8)      , intent(in)     :: qsatlsun               ! Sunlit leaf specific humidity [kg/kg]
    real(r8)      , intent(in)     :: qsatlsha               ! Shalit leaf specific humidity [kg/kg]
    real(r8)      , intent(in)     :: qaf                    ! humidity of canopy air [kg/kg]
    real(r8)      , intent(in)     :: qg                     ! specific humidity at ground surface [kg/kg]
    real(r8)      , intent(in)     :: qm                     ! specific humidity at reference height [kg/kg]
    real(r8)      , intent(in)     :: rhoair                 ! density [kg/m**3]
    real(r8)      , intent(in)     :: psrf                   ! atmospheric pressure [Pa]
    real(r8)      , intent(in)     :: fwet                   ! fraction of foliage that is green and dry [-]
    real(r8)      , intent(in)     :: laisun                 ! Sunlit leaf area index
    real(r8)      , intent(in)     :: laisha                 ! Shaded leaf area index
    real(r8)      , intent(in)     :: htop                   ! canopy top [m] 
    real(r8)      , intent(in)     :: sai                    ! stem area index
    real(r8)      , intent(in)     :: tlsun                ! sunlit leaf temperature
    real(r8)      , intent(in)     :: tlsha                ! shaded leaf temperature
    real(r8)      , intent(in)     :: sigf                   ! fraction of veg cover, excluding snow-covered veg [-]
    real(r8)      , intent(in)     :: kmax_xyl
    real(r8)      , intent(in)     :: kmax_root
    real(r8)      , intent(in)     :: psi50_sun              ! water potential at 50% loss of sunlit leaf tissue conductance (mmH2O)
    real(r8)      , intent(in)     :: psi50_sha              ! water potential at 50% loss of shaded leaf tissue conductance (mmH2O)
    real(r8)      , intent(in)     :: psi50_xyl              ! water potential at 50% loss of xylem tissue conductance (mmH2O)
    real(r8)      , intent(in)     :: psi50_root             ! water potential at 50% loss of root tissue conductance (mmH2O)
    real(r8)      , intent(in)     :: ck                     ! 
    real(r8)      , intent(in)     :: raw        ! moisture resistance [s/m]
    real(r8)      , intent(in)     :: rd         ! aerodynamical resistance between ground and canopy air
    real(r8)      , intent(in)     :: smp(nl_soil)           ! soil matrix potential
    real(r8)      , intent(in)     :: k_soil_root(nl_soil)   ! soil-root interface conductance [mm/s]
    real(r8)      , intent(in)     :: k_ax_root(nl_soil)     ! root axial-direction conductance [mm/s]
    real(r8)      , intent(out)    :: soilflux               ! total soil column transpiration [mm/s]
    real(r8)      , intent(out)    :: etrsun      ! transpiration from sunlit leaf (mm/s)
    real(r8)      , intent(out)    :: etrsha      ! transpiration from shaded leaf (mm/s)
    real(r8)      , intent(out)    :: rootr(nl_soil)        ! root water uptake from different layers
    !
    ! !LOCAL VARIABLES:
!    real(r8) :: qflx_sun                 ! Sunlit leaf transpiration [kg/m2/s]
!    real(r8) :: qflx_sha                 ! Shaded leaf transpiration [kg/m2/s] 
    real(r8) :: qeroot
    real(r8) :: dummy
    real(r8) :: fx                       ! fraction of maximum conductance, xylem-to-leaf [-]  
    real(r8) :: fr                       ! fraction of maximum conductance, root-to-xylem [-]  
    real(r8) :: x_root_top
    real(r8) :: xroot(nl_soil)
    real(r8) :: grav1                    ! gravitational potential surface to canopy top (mm H2O)
    real(r8) :: grav2(nl_soil)           ! soil layer gravitational potential relative to surface (mm H2O) 
    integer  :: j                        ! index
    logical  :: havegs                   ! signals direction of calculation gs->qflx or qflx->gs 
    logical  :: haroot                   ! signals direction of calculation x_root_top->qeroot or qeroot->x_root_top
    integer, parameter :: leafsun=1
    integer, parameter :: leafsha=2
    integer, parameter :: xyl=3
    integer, parameter :: root=4
    !----------------------------------------------------------------------
    grav1 = 1000._r8 * htop
    grav2(1:nl_soil) = 1000._r8 * z_soi(1:nl_soil)
    
    !compute transpiration demand
    havegs=.true.
    call getqflx_gs2qflx_twoleaf(gb_mol_sun,gb_mol_sha,gs_mol_sun,gs_mol_sha,etrsun,etrsha,qsatlsun,qsatlsha,qaf, &
                 rhoair,psrf,laisun,laisha,sai,fwet,tlsun,tlsha,sigf,raw,rd,qg,qm)
    
    !calculate root water potential
    qeroot = etrsun + etrsha

    call getrootqflx_qe2x(nl_soil,smp,z_soi,k_soil_root,k_ax_root,qeroot,xroot,x_root_top)

    !calculate xylem water potential
    fr = plc(x(root),psi50_root,ck)
    x(xyl) = x(root) - grav1 - (etrsun+etrsha)/(fr*kmax_root/htop*sai)
    
    !calculate sun/sha leaf water potential
    fx = plc(x(xyl),psi50_xyl,ck)
       x(leafsha) = x(xyl) - (etrsha/(fx*kmax_xyl*laisha))
       x(leafsun) = x(xyl) - (etrsun/(fx*kmax_xyl*laisun))

    !calculate soil flux
    do j = 1,nl_soil
       rootr(j) = k_soil_root(j)*(smp(j)-xroot(j))
    enddo
    soilflux = sum(rootr(:))

  end subroutine getvegwp_twoleaf
  
  !--------------------------------------------------------------------------------
  subroutine getqflx_gs2qflx_oneleaf(gb_mol,gs_mol,qflx,qsatl,qaf, &
                             rhoair,psrf,lai,sai,fwet,tl,sigf,raw,rd,qg,qm)
    ! !DESCRIPTION:
    !  calculate sunlit and shaded transpiration using gb_MOL and gs_MOL
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    real(r8) , intent(in)     :: gb_mol ! Shaded leaf boundary layer conductance (umol H2O/m**2/s), leaf scale
    real(r8) , intent(inout)  :: gs_mol ! Ball-Berry leaf conductance (umol H2O/m**2/s), leaf scale
    real(r8) , intent(inout)  :: qflx   ! Shaded leaf transpiration [kg/m2/s]
    real(r8) , intent(in)     :: qsatl   ! Shaded leaf specific humidity [kg/kg]
    real(r8) , intent(in)     :: qaf        ! humidity of canopy air [kg/kg]
    real(r8) , intent(in)     :: qg         ! specific humidity at ground surface [kg/kg]
    real(r8) , intent(in)     :: qm         ! specific humidity at reference height [kg/kg]
    real(r8) , intent(in)     :: rhoair   ! density (kg/m**3)
    real(r8) , intent(in)     :: psrf  ! atmospheric pressure (Pa)
    real(r8) , intent(in)     :: lai     ! shaded leaf area index (m2/m2)
    real(r8) , intent(in)     :: sai        ! stem area index (m2/m2)
    real(r8) , intent(in)     :: fwet       ! fraction of foliage that is green and dry [-]
    real(r8) , intent(in)     :: tl                ! shaded leaf temperature
    real(r8) , intent(in)     :: sigf       ! fraction of veg cover, excluding snow-covered veg [-]
    real(r8) , intent(in)     :: raw        ! moisture resistance [s/m]
    real(r8) , intent(in)     :: rd         ! aerodynamical resistance between ground and canopy air

    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtl                   ! heat conductance for shaded leaf boundary [m/s]
    real(r8) :: efpot                 ! potential latent energy flux for shaded leaf [kg/m2/s]
    real(r8) :: rppdry               ! fraction of potential evaporation through transp - shaded [-]
    real(r8) :: cf                       ! s m**2/umol -> s/m
    real(r8) :: tprcor                   !tf*psur*100./1.013e5

    real(r8) :: wtaq0                    ! normalized latent heat conductance for air [-]
    real(r8) :: wtgq0                    ! normalized latent heat conductance for ground [-]
    real(r8) :: wtlq0                 ! normalized latent heat cond. for air and shaded leaf [-]
    real(r8) :: wtsqi                    ! latent heat resistance for air, grd and leaf [-]

    real(r8) :: delta             
    real(r8) :: caw                      ! latent heat conductance for air [m/s]
    real(r8) :: cgw                      ! latent heat conductance for ground [m/s]
    real(r8) :: cfw                     ! latent heat conductance for sunlit leaf [m/s]

    !----------------------------------------------------------------------
    tprcor   = 44.6*273.16*psrf/1.013e5
    cf    = tprcor/tl * 1.e6_r8  ! gb->gbmol conversion factor
    wtl   = (lai+sai)*gb_mol

    delta = 0.0
    if(qsatl-qaf .gt. 0.) delta = 1.0

    caw = sigf / raw
    cgw = sigf / rd
    cfw = sigf * ( (1.-delta*(1.-fwet)) * (lai+sai) * gb_mol /cf &
                           + (1. - fwet) * delta * lai / (1._r8/gb_mol+1._r8/gs_mol)/cf)
    wtsqi = 1. / ( caw + cgw + cfw )

    wtaq0     = caw * wtsqi
    wtgq0     = cgw * wtsqi
    wtlq0     = cfw * wtsqi

    efpot     = sigf*rhoair*wtl*delta&
              * ((wtaq0+wtgq0)*qsatl-wtaq0*qm-wtgq0*qg)

    rppdry = (1.-fwet)/gb_mol*(lai/(1._r8/gb_mol+1._r8/gs_mol))/(lai+sai)
    qflx   = efpot*rppdry/cf
    if(qflx < 1.e-7_r8)then
       qflx   = 0._r8
    end if

  end subroutine getqflx_gs2qflx_oneleaf

  !--------------------------------------------------------------------------------
  subroutine getqflx_gs2qflx_twoleaf(gb_mol_sun,gb_mol_sha,gs_mol_sun,gs_mol_sha,qflx_sun,qflx_sha,qsatlsun,qsatlsha,qaf,&
                             rhoair,psrf,laisun,laisha,sai,fwet,tlsun,tlsha,sigf,raw,rd,qg,qm)
    ! !DESCRIPTION:
    !  calculate sunlit and shaded transpiration using gb_MOL and gs_MOL
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    real(r8) , intent(in)     :: gb_mol_sun ! Sunlit leaf boundary layer conductance (mol H2O/m**2/s), leaf scale
    real(r8) , intent(in)     :: gb_mol_sha ! Shaded leaf boundary layer conductance (mol H2O/m**2/s), leaf scale
    real(r8) , intent(inout)  :: gs_mol_sun ! Ball-Berry leaf conductance (mol H2O/m**2/s), leaf scale
    real(r8) , intent(inout)  :: gs_mol_sha ! Ball-Berry leaf conductance (mol H2O/m**2/s), leaf scale
    real(r8) , intent(inout)  :: qflx_sun   ! Sunlit leaf transpiration [kg/m2/s]
    real(r8) , intent(inout)  :: qflx_sha   ! Shaded leaf transpiration [kg/m2/s]
    real(r8) , intent(in)     :: qsatlsun   ! Sunlit leaf specific humidity [kg/kg]
    real(r8) , intent(in)     :: qsatlsha   ! Shaded leaf specific humidity [kg/kg]
    real(r8) , intent(in)     :: qaf        ! humidity of canopy air [kg/kg]
    real(r8) , intent(in)     :: qg         ! specific humidity at ground surface [kg/kg]
    real(r8) , intent(in)     :: qm         ! specific humidity at reference height [kg/kg]
    real(r8) , intent(in)     :: rhoair   ! density (kg/m**3)
    real(r8) , intent(in)     :: psrf  ! atmospheric pressure (Pa)
    real(r8) , intent(in)     :: laisun     ! sunlit leaf area index (m2/m2)
    real(r8) , intent(in)     :: laisha     ! shaded leaf area index (m2/m2)
    real(r8) , intent(in)     :: sai        ! stem area index (m2/m2)
    real(r8) , intent(in)     :: fwet       ! fraction of foliage that is green and dry [-]
    real(r8) , intent(in)     :: tlsun      ! shaded leaf temperature
    real(r8) , intent(in)     :: tlsha      ! shaded leaf temperature
    real(r8) , intent(in)     :: sigf       ! fraction of veg cover, excluding snow-covered veg [-]
    real(r8) , intent(in)     :: raw        ! moisture resistance [s/m]
    real(r8) , intent(in)     :: rd         ! aerodynamical resistance between ground and canopy air

    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtlsun                   ! heat conductance for sunlit leaf boundary [m/s]
    real(r8) :: wtlsha                   ! heat conductance for shaded leaf boundary [m/s]
    real(r8) :: efpotsun                 ! potential latent energy flux for sunlit leaf [kg/m2/s]
    real(r8) :: efpotsha                 ! potential latent energy flux for shaded leaf [kg/m2/s]
    real(r8) :: rppdry_sun               ! fraction of potential evaporation through transp - sunlit [-]
    real(r8) :: rppdry_sha               ! fraction of potential evaporation through transp - shaded [-]
    real(r8) :: cfsun                    ! s m**2/umol -> s/m
    real(r8) :: cfsha                    ! s m**2/umol -> s/m
    real(r8) :: tprcor                   !tf*psur*100./1.013e5

    real(r8) :: wtaq0                    ! normalized latent heat conductance for air [-]
    real(r8) :: wtgq0                    ! normalized latent heat conductance for ground [-]
    real(r8) :: wtlsunq0                 ! normalized latent heat cond. for air and sunlit leaf [-]
    real(r8) :: wtlshaq0                 ! normalized latent heat cond. for air and shaded leaf [-]
    real(r8) :: wtsqi                    ! latent heat resistance for air, grd and leaf [-]

    real(r8) :: delta1             
    real(r8) :: delta2             
    real(r8) :: caw                      ! latent heat conductance for air [m/s]
    real(r8) :: cgw                      ! latent heat conductance for ground [m/s]
    real(r8) :: cfsunw                   ! latent heat conductance for sunlit leaf [m/s]
    real(r8) :: cfshaw                   ! latent heat conductance for shaded leaf [m/s]

    !----------------------------------------------------------------------
    tprcor   = 44.6*273.16*psrf/1.013e5
    cfsun    = tprcor/tlsun * 1.e6_r8  ! gb->gbmol conversion factor
    cfsha    = tprcor/tlsha * 1.e6_r8  ! gb->gbmol conversion factor
    wtlsun   = (laisun+laisha+sai)*gb_mol_sun
    wtlsha   = (laisun+laisha+sai)*gb_mol_sha

    delta1 = 0.0
    delta2 = 0.0
    if(qsatlsun-qaf .gt. 0.) delta1 = 1.0
    if(qsatlsha-qaf .gt. 0.) delta2 = 1.0

    caw = sigf / raw
    cgw = sigf / rd
    cfsunw = sigf * ( (1.-delta1*(1.-fwet)) * laisun * gb_mol_sun /cfsun &
                         + (1. - fwet) * delta1 * laisun / (1._r8/gb_mol_sun+1._r8/gs_mol_sun)/cfsun)
    cfshaw = sigf * ( (1.-delta2*(1.-fwet)) * (laisha+sai) * gb_mol_sha /cfsha &
                         + (1. - fwet) * delta2 * laisha / (1._r8/gb_mol_sha+1._r8/gs_mol_sha)/cfsha)
    wtsqi = 1. / ( caw + cgw + cfsunw + cfshaw )

    wtaq0     = caw * wtsqi
    wtgq0     = cgw * wtsqi
    wtlsunq0  = cfsunw * wtsqi
    wtlshaq0  = cfshaw * wtsqi

    efpotsun = sigf*rhoair*wtlsun*delta1&
             * ((wtaq0+wtgq0+wtlshaq0)*qsatlsun-wtaq0*qm-wtgq0*qg-wtlshaq0*qsatlsha)
    efpotsha = sigf*rhoair*wtlsha*delta2&
             * ((wtaq0+wtgq0+wtlsunq0)*qsatlsha-wtaq0*qm-wtgq0*qg-wtlsunq0*qsatlsun)


    rppdry_sun = (1.-fwet)/gb_mol_sun*(laisun/(1._r8/gb_mol_sun+1._r8/gs_mol_sun))/(laisun+laisha+sai)
    qflx_sun   = efpotsun*rppdry_sun/cfsun
    if(qflx_sun < 1.e-7_r8)then
       qflx_sun   = 0._r8
    end if
    rppdry_sha = (1.-fwet)/gb_mol_sha*(laisha/(1._r8/gb_mol_sha+1._r8/gs_mol_sha))/(laisun+laisha+sai)
    qflx_sha   = efpotsha*rppdry_sha/cfsha
    if(qflx_sha < 1.e-7)then
       qflx_sha   = 0._r8
    end if

  end subroutine getqflx_gs2qflx_twoleaf

  subroutine getqflx_qflx2gs_oneleaf(gb_mol,gs_mol,qflx,qsatl,qaf,rhoair,psrf,lai,sai,&
                                     fwet,tl,sigf,raw,rd,qg,qm)
    ! !DESCRIPTION:
    !  calculate sunlit and shaded transpiration using gb_MOL and gs_MOL
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    real(r8) , intent(in)     :: gb_mol ! Shaded leaf boundary layer conductance (mol H2O/m**2/s), leaf scale
    real(r8) , intent(inout)  :: gs_mol ! Ball-Berry leaf conductance (mol H2O/m**2/s), leaf scale
    real(r8) , intent(in)     :: qflx   ! Shaded leaf transpiration [kg/m2/s]
    real(r8) , intent(in)     :: qsatl   ! Shaded leaf specific humidity [kg/kg]
    real(r8) , intent(in)     :: qaf        ! humidity of canopy air [kg/kg]
    real(r8) , intent(in)     :: qg         ! specific humidity at ground surface [kg/kg]
    real(r8) , intent(in)     :: qm         ! specific humidity at reference height [kg/kg]
    real(r8) , intent(in)     :: rhoair   ! density (kg/m**3)
    real(r8) , intent(in)     :: psrf  ! atmospheric pressure (Pa)
    real(r8) , intent(in)     :: lai     ! shaded leaf area index (m2/m2)
    real(r8) , intent(in)     :: sai        ! stem area index (m2/m2)
    real(r8) , intent(in)     :: fwet       ! fraction of foliage that is green and dry [-]
    real(r8) , intent(in)     :: tl                ! shaded leaf temperature
    real(r8) , intent(in)     :: sigf       ! fraction of veg cover, excluding snow-covered veg [-]
    real(r8) , intent(in)     :: raw        ! moisture resistance [s/m]
    real(r8) , intent(in)     :: rd         ! aerodynamical resistance between ground and canopy air

    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtl                   ! heat conductance for shaded leaf boundary [m/s]
    real(r8) :: efpot                 ! potential latent energy flux for shaded leaf [kg/m2/s]
    real(r8) :: rppdry               ! fraction of potential evaporation through transp - shaded [-]
    real(r8) :: cf                       ! s m**2/umol -> s/m
    real(r8) :: tprcor                   !tf*psur*100./1.013e5

    real(r8) :: wtaq0                    ! normalized latent heat conductance for air [-]
    real(r8) :: wtgq0                    ! normalized latent heat conductance for ground [-]
    real(r8) :: wtlq0                 ! normalized latent heat cond. for air and shaded leaf [-]
    real(r8) :: cqi_wet                  ! latent heat conductance for air, grd and wet leaf [-]

    real(r8) :: delta             
    real(r8) :: caw                      ! latent heat conductance for air [m/s]
    real(r8) :: cgw                      ! latent heat conductance for ground [m/s]
    real(r8) :: cfw_dry                  ! latent heat conductance for dry leaf [m/s]
    real(r8) :: cfw_wet                  ! latent heat conductance for wet leaf [m/s]

    !----------------------------------------------------------------------
    if (qflx > 0._r8)then
    tprcor   = 44.6*273.16*psrf/1.013e5
    cf    = tprcor/tl * 1.e6_r8  ! gb->gbmol conversion factor
    wtl   = (lai+sai)*gb_mol

    delta = 0.0
    if(qsatl-qaf .gt. 0.) delta = 1.0

    caw = sigf / raw
    cgw = sigf / rd
    cfw_wet = sigf * (1.-delta*(1.-fwet)) * (lai+sai) * gb_mol /cf 
    cqi_wet = caw + cgw + cfw_wet 

    cfw_dry = qflx/rhoair*cqi_wet / (caw*(qsatl-qm)+cgw*(qsatl-qg)-qflx/rhoair)
    gs_mol = 1._r8 / (sigf * (1. - fwet) * delta * lai / cfw_dry / cf - 1._r8 / gb_mol) 
    endif
    
  end subroutine getqflx_qflx2gs_oneleaf

  subroutine getqflx_qflx2gs_twoleaf(gb_mol_sun,gb_mol_sha,gs_mol_sun,gs_mol_sha,qflx_sun,qflx_sha,qsatlsun,qsatlsha,qaf, &
                     rhoair,psrf,laisun,laisha,sai,fwet,tlsun,tlsha,sigf,raw,rd,qg,qm)
    ! !DESCRIPTION:
    !  calculate sunlit and shaded transpiration using gb_MOL and gs_MOL
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    real(r8) , intent(in)     :: gb_mol_sun ! Sunlit leaf boundary layer conductance (mol H2O/m**2/s), leaf scale
    real(r8) , intent(in)     :: gb_mol_sha ! Shaded leaf boundary layer conductance (mol H2O/m**2/s), leaf scale
    real(r8) , intent(inout)  :: gs_mol_sun ! Ball-Berry leaf conductance (mol H2O/m**2/s), leaf scale
    real(r8) , intent(inout)  :: gs_mol_sha ! Ball-Berry leaf conductance (mol H2O/m**2/s), leaf scale
    real(r8) , intent(inout)  :: qflx_sun   ! Sunlit leaf transpiration [kg/m2/s]
    real(r8) , intent(inout)  :: qflx_sha   ! Shaded leaf transpiration [kg/m2/s]
    real(r8) , intent(in)     :: qsatlsun   ! Sunlit leaf specific humidity [kg/kg]
    real(r8) , intent(in)     :: qsatlsha   ! Shaded leaf specific humidity [kg/kg]
    real(r8) , intent(in)     :: qaf        ! humidity of canopy air [kg/kg]
    real(r8) , intent(in)     :: qg         ! specific humidity at ground surface [kg/kg]
    real(r8) , intent(in)     :: qm         ! specific humidity at reference height [kg/kg]
    real(r8) , intent(in)     :: rhoair   ! density (kg/m**3)
    real(r8) , intent(in)     :: psrf  ! atmospheric pressure (Pa)
    real(r8) , intent(in)     :: laisun     ! sunlit leaf area index (m2/m2)
    real(r8) , intent(in)     :: laisha     ! shaded leaf area index (m2/m2)
    real(r8) , intent(in)     :: sai        ! stem area index (m2/m2)
    real(r8) , intent(in)     :: fwet       ! fraction of foliage that is green and dry [-]
    real(r8) , intent(in)     :: tlsun                ! sunlit leaf temperature
    real(r8) , intent(in)     :: tlsha                ! shaded leaf temperature
    real(r8) , intent(in)     :: sigf       ! fraction of veg cover, excluding snow-covered veg [-]
    real(r8) , intent(in)     :: raw        ! moisture resistance [s/m]
    real(r8) , intent(in)     :: rd         ! aerodynamical resistance between ground and canopy air

    !
    ! !LOCAL VARIABLES:
    real(r8) :: wtlsun                   ! heat conductance for sunlit leaf boundary [m/s]
    real(r8) :: wtlsha                   ! heat conductance for shaded leaf boundary [m/s]
    real(r8) :: efpotsun                 ! potential latent energy flux for sunlit leaf [kg/m2/s]
    real(r8) :: efpotsha                 ! potential latent energy flux for shaded leaf [kg/m2/s]
    real(r8) :: rppdry_sun               ! fraction of potential evaporation through transp - sunlit [-]
    real(r8) :: rppdry_sha               ! fraction of potential evaporation through transp - shaded [-]
    real(r8) :: cfsun                       ! s m**2/umol -> s/m
    real(r8) :: cfsha                       ! s m**2/umol -> s/m
    real(r8) :: tprcor                   !tf*psur*100./1.013e5

    real(r8) :: wtaq0                    ! normalized latent heat conductance for air [-]
    real(r8) :: wtgq0                    ! normalized latent heat conductance for ground [-]
    real(r8) :: wtlsunq0                 ! normalized latent heat cond. for air and sunlit leaf [-]
    real(r8) :: wtlshaq0                 ! normalized latent heat cond. for air and shaded leaf [-]

    real(r8) :: delta1             
    real(r8) :: delta2             
    real(r8) :: caw                      ! latent heat conductance for air [m/s]
    real(r8) :: cgw                      ! latent heat conductance for ground [m/s]
    real(r8) :: cfsunw_dry               ! latent heat conductance for sunlit dry leaf [m/s]
    real(r8) :: cfsunw_wet               ! latent heat conductance for sunlit wet leaf [m/s]
    real(r8) :: cfshaw_dry               ! latent heat conductance for shaded dry leaf [m/s]
    real(r8) :: cfshaw_wet               ! latent heat conductance for shaded wet leaf [m/s]
    real(r8) :: cqi_wet                  ! latent heat conductance for air, grd and wet leaf [-]
    real(r8) :: Delta,deltaq,          & ! temporary variables to solve cfsunw_dry and cfshaw_dry 
                A1,B1,C1,A2,B2,C2        ! in binary quadratic equations

    !----------------------------------------------------------------------
    if(qflx_sun .gt. 0 .or. qflx_sha .gt. 0)then
    tprcor   = 44.6*273.16*psrf/1.013e5
    cfsun    = tprcor/tlsun * 1.e6_r8  ! gb->gbmol conversion factor
    cfsha    = tprcor/tlsha * 1.e6_r8  ! gb->gbmol conversion factor
    wtlsun   = (laisun+laisha+sai)*gb_mol_sun
    wtlsha   = (laisun+laisha+sai)*gb_mol_sha

    delta1 = 0.0
    delta2 = 0.0
    if(qsatlsun-qaf .gt. 0.) delta1 = 1.0
    if(qsatlsha-qaf .gt. 0.) delta2 = 1.0

    caw = sigf / raw
    cgw = sigf / rd
    cfsunw_wet = sigf * (1.-delta1*(1.-fwet)) * laisun * gb_mol_sun / cfsun
    cfshaw_wet = sigf * (1.-delta2*(1.-fwet)) * (laisha+sai) * gb_mol_sha / cfsha

    cqi_wet = ( caw + cgw + cfsunw_wet + cfshaw_wet )

    deltaq = qsatlsun - qsatlsha

       B1 = - qflx_sun / rhoair
       A1 = caw * (qsatlsun - qm) + cgw * (qsatlsun - qg) + cfshaw_wet * deltaq + B1
       C1 = qflx_sun / rhoair * cqi_wet
       B2 = - qflx_sha / rhoair
       A2 = caw * (qsatlsha - qm) + cgw * (qsatlsha - qg) - cfsunw_wet * deltaq + B2
       C2 = qflx_sha / rhoair * cqi_wet

    if(deltaq .ne. 0 .and. qflx_sun .gt. 1.e-20 .and. qflx_sha .gt. 1.e-20)then   !solve equations deltaq*cfsunw_dry*cfshaw_dry+A1*cfsunw_dry+B1*cfshaw_dry-C1 = 0
                                                                                 !             && -deltaq*cfsunw_dry*cfshaw_dry+A2*cfshaw_dry+B2*cfsunw_dry-C2 = 0
       Delta = A1**2*A2**2-2*A1*A2*B1*B2-2*A1*A2*C1*deltaq+2*A1*A2*C2*deltaq+4*A1*B1*C2*deltaq-4*A2*B2*C1*deltaq &
             + B1**2*B2**2-2*B1*B2*C1*deltaq+2*B1*B2*C2*deltaq+C1**2*deltaq**2+2*C1*C2*deltaq**2+C2**2*deltaq**2

       cfsunw_dry = (C1*deltaq+C2*deltaq-sqrt(Delta)+A1*A2-B1*B2)/(2*(A1+B2)*deltaq)
       cfshaw_dry = (C1*deltaq+C2*deltaq+sqrt(Delta)-A1*A2+B1*B2)/(2*(A2+B1)*deltaq)
    else ! solve equations A1*cfsunw_dry+B1*cfshaw_dry-C1 = 0
         !              && A2*cfshaw_dry+B2*cfsunw_dry-C2 = 0
       cfsunw_dry = (A2*C1-B1*C2)/(A1*A2-B1*B2) 
       cfshaw_dry = (A1*C2-B2*C1)/(A1*A2-B1*B2)
    end if

       if (qflx_sun > 0._r8) then
          gs_mol_sun = 1._r8 / (sigf * (1. - fwet) * delta1 * laisun / cfsunw_dry / cfsun - 1._r8 / gb_mol_sun) 
       endif
       if (qflx_sha > 0._r8) then
          gs_mol_sha = 1._r8 / (sigf * (1. - fwet) * delta2 * laisha / cfshaw_dry / cfsha - 1._r8 / gb_mol_sha)
       endif
    end if
    
  end subroutine getqflx_qflx2gs_twoleaf

  subroutine getrootqflx_x2qe(nl_soil,smp,x_root_top,z_soisno,krad,kax,qeroot,dqeroot)

  ! DESCRIPTION
  ! Return root water potential at top soil node. Return soil-root water flux.
  !

  integer ,intent(in)    :: nl_soil
  real(r8),intent(in)    :: smp      (nl_soil)
  real(r8),intent(in) :: x_root_top
  real(r8),intent(in)    :: z_soisno (nl_soil)
  real(r8),intent(in)    :: krad     (nl_soil)
  real(r8),intent(in)    :: kax      (nl_soil)
  real(r8),intent(out) :: qeroot
  real(r8),intent(out) :: dqeroot

! Local variables
  real(r8) :: den_AHR,den1,den2            ! used in calculating HR(Amenu model)
  real(r8),dimension(nl_soil-1) :: amx_hr             ! "a" left off diagonal of tridiagonal matrix
  real(r8),dimension(nl_soil-1) :: bmx_hr             ! "b" diagonal column for tridiagonal matrix
  real(r8),dimension(nl_soil-1) :: cmx_hr             ! "c" right off diagonal tridiagonal matrix
  real(r8),dimension(nl_soil-1) :: rmx_hr             ! "r" forcing term of tridiagonal matrix
  real(r8),dimension(nl_soil-1) :: drmx_hr            ! "dr" forcing term of tridiagonal matrix for d/dxroot(1)
  real(r8),dimension(nl_soil-1) :: x  ! root water potential from layer 2 to nl_soil
  real(r8),dimension(nl_soil-1) :: dx ! derivate of root water potential from layer 2 to nl_soil (dxroot(:)/dxroot(1))
  real(r8),dimension(nl_soil) :: xroot  ! root water potential from layer 2 to nl_soil
  real(r8) :: zmm(1:nl_soil)     ! layer depth [mm]
  real(r8) :: qeroot_nl(1:nl_soil) ! root water potential from layer 2 to nl_soil
  real(r8) :: dxroot2    ! dxroot(2)/dxroot(1)
  integer j

  ! Because the depths in this routine are in mm, use local
  ! variable arrays instead of pointers
  do j = 1, nl_soil
     zmm(j)  = z_soisno(j)*1000.
  end do

  xroot(1) = x_root_top + zmm(1)
  ! For the 2nd soil layer
  j            = 2
  den1         = zmm(j) - zmm(j-1)
  den2         = zmm(j+1) - zmm(j)
  amx_hr(j-1)  = 0
  bmx_hr(j-1)  = kax(j-1)/den1 + kax(j)/den2 + krad(j)
  cmx_hr(j-1)  = -kax(j)/den2
  rmx_hr(j-1)  = krad(j)*smp(j) + kax(j-1) - kax(j) + kax(j-1)/den1*xroot(1)
  drmx_hr(j-1) = kax(j-1)/den1

  ! For the middile soil layers
  do j = 3, nl_soil - 1
     den1   = zmm(j) - zmm(j-1)
     den2   = zmm(j+1) - zmm(j)
     amx_hr (j-1) = -kax(j-1)/den1
     bmx_hr (j-1) = kax(j-1)/den1 + kax(j)/den2 + krad(j)
     cmx_hr (j-1) = -kax(j)/den2
     rmx_hr (j-1) = krad(j)*smp(j) + kax(j-1) - kax(j)
     drmx_hr(j-1) = 0._r8
  end do

  ! For the bottom soil layer
  j           = nl_soil
  den_AHR     = zmm(j) - zmm(j-1)
  amx_hr (j-1) = -kax(j-1)/den_AHR
  bmx_hr (j-1) = kax(j-1)/den_AHR + krad(j)
  cmx_hr (j-1) = 0
  rmx_hr (j-1) = krad(j)*smp(j) + kax(j-1)
  drmx_hr(j-1) = 0._r8

  ! Solve for root pressure potential using tridiagonal matric solver x = A^-1 * r
  call tridia (nl_soil-1 ,amx_hr ,bmx_hr ,cmx_hr ,rmx_hr ,x)

  do j = 2,nl_soil
     xroot(j) = x(j-1)
  end do

     ! Solve the dx(:)/dxroot(1) = A^-1 * dr
  call tridia (nl_soil-1 ,amx_hr ,bmx_hr ,cmx_hr ,drmx_hr, dx)

  dxroot2 = dx(1)

  ! calculate the water flux
  j      = 1
  den2   = zmm(j+1) - zmm(j)
  qeroot = krad(j) * (smp(1) - xroot(1)) + (xroot(2) - xroot(1)) * kax(j)/den2 - kax(j)

  ! calculate the dqeroot/dx_root_top; 
  dqeroot = - krad(j) + (dxroot2 - 1) * kax(j)/den2
  do j = 1,nl_soil
     qeroot_nl(j) = krad(j)*(smp(j) - xroot(j))
  end do

  end subroutine getrootqflx_x2qe

  subroutine getrootqflx_qe2x(nl_soil,smp,z_soisno,krad,kax,qeroot,xroot,x_root_top)

  ! DESCRIPTION
  ! Return root water potential at top soil node. Return soil-root water flux.
  !

  integer ,intent(in)    :: nl_soil
  real(r8),intent(in)    :: smp      (nl_soil)
  real(r8),intent(in)    :: z_soisno (nl_soil)
  real(r8),intent(in)    :: krad     (nl_soil)
  real(r8),intent(in)    :: kax      (nl_soil)
  real(r8),intent(in)    :: qeroot
  real(r8),intent(out)   :: xroot    (nl_soil)
  real(r8),intent(out)   :: x_root_top

! Local variables
  real(r8) :: den_AHR,den1,den2            ! used in calculating HR(Amenu model)
  real(r8),dimension(nl_soil) :: amx_hr             ! "a" left off diagonal of tridiagonal matrix
  real(r8),dimension(nl_soil) :: bmx_hr             ! "b" diagonal column for tridiagonal matrix
  real(r8),dimension(nl_soil) :: cmx_hr             ! "c" right off diagonal tridiagonal matrix
  real(r8),dimension(nl_soil) :: rmx_hr             ! "r" forcing term of tridiagonal matrix
  real(r8),dimension(nl_soil) :: x  ! root water potential from layer 2 to nl_soil
  real(r8) :: zmm(1:nl_soil)     ! layer depth [mm]
  real(r8) :: qeroot_nl(1:nl_soil) ! root water potential from layer 2 to nl_soil
  integer j

  ! Because the depths in this routine are in mm, use local
  ! variable arrays instead of pointers
  do j = 1, nl_soil
     zmm(j)  = z_soisno(j)*1000.
  end do

  j           = 1 
  den2        = zmm(j+1) - zmm(j)
  amx_hr(j)   = 0
  bmx_hr(j) = kax(j)/den2 + krad(j)
  cmx_hr(j) = -kax(j)/den2
  rmx_hr(j) = krad(j)*smp(j) - qeroot - kax(j) 

  ! For the middile soil layers
  do j = 2, nl_soil - 1
      den1   = zmm(j) - zmm(j-1)
      den2   = zmm(j+1) - zmm(j)
      amx_hr(j) = -kax(j-1)/den1
      bmx_hr(j) = kax(j-1)/den1 + kax(j)/den2 + krad(j)
      cmx_hr(j) = -kax(j)/den2
      rmx_hr(j) = krad(j)*smp(j) + kax(j-1) - kax(j)
  end do

  ! For the bottom soil layer
  j      = nl_soil
  den_AHR    = zmm(j) - zmm(j-1)
  amx_hr(j) = -kax(j-1)/den_AHR
  bmx_hr(j) = kax(j-1)/den_AHR + krad(j)
  cmx_hr(j) = 0
  rmx_hr(j) = krad(j)*smp(j) + kax(j-1)

  ! Solve for root pressure potential using tridiagonal matric solver
  call tridia (nl_soil ,amx_hr ,bmx_hr ,cmx_hr ,rmx_hr ,x)

  xroot(1:nl_soil) = x(1:nl_soil)
  x_root_top = xroot(1) - zmm(1)

  end subroutine getrootqflx_qe2x

  !--------------------------------------------------------------------------------
  function plc(x,psi50,ck)
    ! !DESCRIPTION
    ! Return value of vulnerability curve at x
    !
    ! !ARGUMENTS
    real(r8) , intent(in)  :: x             ! water potential input
!    integer  , intent(in)  :: level         ! veg segment lvl (1:nvegwcs) 
!    integer  , intent(in)  :: plc_method    !
    real(r8) , intent(in)  :: psi50     ! water potential at 50% loss of sunlit leaf tissue conductance (mmH2O)
!    real(r8) , intent(in)  :: psi50_sun     ! water potential at 50% loss of sunlit leaf tissue conductance (mmH2O)
!    real(r8) , intent(in)  :: psi50_sha     ! water potential at 50% loss of shaded leaf tissue conductance (mmH2O)
!    real(r8) , intent(in)  :: psi50_xyl     ! water potential at 50% loss of xylem tissue conductance (mmH2O)
!    real(r8) , intent(in)  :: psi50_root    ! water potential at 50% loss of root tissue conductance (mmH2O)
    real(r8) , intent(in)  :: ck
    real(r8)               :: plc           ! attenuated conductance [0:1] 0=no flow
    !
    ! !PARAMETERS
!    integer , parameter :: vegetation_weibull=0  ! case number
!    integer , parameter :: leafsun = 1  ! index for sunlit leaf
!    integer , parameter :: leafsha = 2  ! index for shaded leaf
!    integer , parameter :: xyl     = 3  ! index for xylem
!    integer , parameter :: root    = 4  ! index for root

    ! !LOCAL VARIABLES
    !real(r8) psi50,tmp
    real(r8) tmp
    integer i
    
    !------------------------------------------------------------------------------
!    select case(level)
!    case (leafsun)
!       psi50 = psi50_sun
!    case (leafsha)
!       psi50 = psi50_sha
!    case (xyl)
!       psi50 = psi50_xyl
!    case (root)
!       psi50 = psi50_root
!    case default
!       write(*,*),'must choose level from 1 to 4 (sunlit leaf to root)'
!    end select

!    select case (plc_method)
       !possible to add other methods later
!    case (vegetation_weibull)
       tmp = amax1(-(x/psi50)**ck,-500._r8)
!       if(tmp .lt. -500._r8)then
!          plc = 0._r8
!       else
       plc=2._r8**tmp
!       end if
!       if ( plc < 0.00001_r8) plc = 0._r8
!    case default
!       write(*,*),'must choose plc method'
!    end select

  end function plc
  !--------------------------------------------------------------------------------
  
  !--------------------------------------------------------------------------------
  function d1plc(x,psi50,ck)
    ! !DESCRIPTION
    ! Return 1st derivative of vulnerability curve at x
    !
    ! !ARGUMENTS
    real(r8) , intent(in) :: x                ! water potential input
!    integer  , intent(in) :: level            ! veg segment lvl (1:nvegwcs)
!    integer  , intent(in) :: plc_method       ! 0 for vegetation, 1 for soil
    real(r8) , intent(in) :: psi50        ! water potential at 50% loss of sunlit leaf tissue conductance (mmH2O)
!    real(r8) , intent(in) :: psi50_sun        ! water potential at 50% loss of sunlit leaf tissue conductance (mmH2O)
!    real(r8) , intent(in) :: psi50_sha        ! water potential at 50% loss of shaded leaf tissue conductance (mmH2O)
!    real(r8) , intent(in) :: psi50_xyl        ! water potential at 50% loss of xylem tissue conductance (mmH2O)
!    real(r8) , intent(in) :: psi50_root       ! water potential at 50% loss of root tissue conductance (mmH2O)
    real(r8) , intent(in) :: ck
    real(r8)              :: d1plc            ! first deriv of plc curve at x
    !
    ! !PARAMETERS
!    integer , parameter :: vegetation_weibull=0  ! case number
!    integer , parameter :: leafsun = 1  ! index for sunlit leaf
!    integer , parameter :: leafsha = 2  ! index for shaded leaf
!    integer , parameter :: xyl     = 3  ! index for xylem
!    integer , parameter :: root    = 4  ! index for root

    ! !LOCAL VARIABLES
!    real(r8) psi50,tmp
    real(r8) tmp
    !------------------------------------------------------------------------------
!    select case(level)
!    case (leafsun)
!       psi50 = psi50_sun
!    case (leafsha)
!       psi50 = psi50_sha
!    case (xyl)
!       psi50 = psi50_xyl
!    case (root)
!       psi50 = psi50_root
!    case default
!       write(*,*),'must choose level from 1 to 4 (sunlit leaf to root)'
!    end select

!    select case (plc_method)
       !possible to add other methods later
!    case (vegetation_weibull)
       tmp = amax1(-(x/psi50)**ck,-500._r8)
!       if(tmp .lt. -500._r8)then
!          d1plc = 0._r8
!       else
          d1plc= ck * log(2._r8) * (2._r8**tmp) * tmp / x
!       end if
!    case default
!       write(*,*),'must choose plc method'
!    end select

  end function d1plc  
  

END MODULE PlantHydraulic
! -------------- EOP ---------------
