#include <define.h>

MODULE  MOD_Lake_Utils

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: LakeIni
   PUBLIC :: LakeRst
   PUBLIC :: LakeDisper
   PUBLIC :: LakeIceFrc
   PUBLIC :: LksoilLayer
   PUBLIC :: snowage
   PUBLIC :: snowdp2lev4lake
   PUBLIC :: reverse_array
    
!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------

   SUBROUTINE LakeIni( &
               ! "in" arguments
               ! -------------------
               nlake      , nsnow     , nsoil     , nlice     ,& 
               dplak      , tskin     ,&
               ! "out" arguments
               ! -------------------
               zlake      , zilak     , dzlak     , lktmp     ,& 
               rhosnw     , lkrho    , icefr     , stke1     ,&
               tmsno      , tmice     , tmmnw     , tmwml     ,&
               tmbot      , tmups     , mldp      , upsdp     ,&
               CTfrac     , icedp     , bicedp    , wicedp    ,&
               uwatv      , vwatv     , lksal     , tke       ,&
               eps        , etke      , num       , nuh       ,&
               ziarea     )
                  
! ---------------------------------- code history -------------------------------------
! Description:
!     +xum:  this one should be called after the cssptci, 
!            the function is to initilize the variables used in lake
!     +WMEJ: All schemes related to lake stratification should be retained here as options.
!
!
! Original author: 
!     Min Xu, 2012
!
! Revisions:
!     Xin-Zhong Liang,
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_Lake_Const, only: tfrz, denh2o, denice, tkwat, densnow, C_T_min, C_T_max,&
                                    ShallowDeepPartition, SHALLOW, DEEP 
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      nlake                   ,&! Maximum number of lake layer
      nsnow                   ,&! Maximum number of snow layer
      nsoil                   ,&! Maximum number of soil layer
      nlice                     ! Maximum number of ice layer
      
   real(r8), intent(in)     :: &
      dplak                   ,&! lake depth, for initialization [m]
      tskin                     ! ground surface temperature [k]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      dzlak(nlake)            ,&! lake layer thickness [m]
      zlake(nlake)            ,&! lake layer depth [m]
      zilak(nlake+1)          ,&! interface level below a "z" level [m]
      lktmp(nlake)            ,&! lake temperature (K)
      uwatv(nlake)            ,&! Water velocity in x-direction [m/s]
      vwatv(nlake)            ,&! Water velocity in y-direction [m/s]
      lksal(nlake)            ,&! Salinity [‰]
      lkrho(nlake)            ,&! Water density [kg/m^3]
      icefr(nlake)            ,&! lake ice fraction [-]
      ziarea(nlake+1)         ,&! lake layer interface area [m2]
      tke(nlake+1)            ,&! Turbulent kinetic energy (TKE) [J/kg]
      eps(nlake+1)            ,&! TKE dissipation rate [W/kg]
      num(nlake+1)            ,&! Turbulent viscosity (momentum)
      nuh(nlake+1)              ! Turbulent diffusivity (heat)
         
   real(r8), intent(out)    :: &
      rhosnw                  ,&! snow density (kg/m3)
      icedp                   ,&! ice depth [m]
      bicedp                  ,&! black snow depth (m)
      wicedp                  ,&! white snow depth (m)
      stke1                   ,&! top level eddy conductivity [W/m/K]
      etke                    ,&! Seiche energy [J]
      tmsno                   ,&! snow temperature [K]
      tmice                   ,&! Temperature at the snow-ice or air-ice interface [K]
      tmmnw                   ,&! Mean temperature of the water column [K]
      tmwml                   ,&! Mixed-layer temperature [K]
      tmbot                   ,&! Temperature at the water-bottom sediment interface [K]
      tmups                   ,&! Temperature at the bottom of the upper layer of the sediments [K]
      mldp                    ,&! mixed layer depth [m]
      upsdp                   ,&! bottom of the upper layer of the sediments [m]
      CTfrac                    ! Shape factor (thermocline)

!  ------------------------- local variables ---------------------------
   real(r8)                 :: &
      layerVolume             ,&! layer volume [m3]
      uparea                  ,&! upper layer interface area [m2]
      lowarea                 ,&! lower layer interface area [m2]
      dzlayer                   ! layer thickness [m]

   integer                  :: &
      scwat                   ,&! water body type
      b                       ,&! loop index
      i                         ! loop index

   real(r8)                 :: &
      wice_lake(nlake)        ,&! ice lens [kg/m2]
      wliq_lake(nlake)          ! liquid water [kg/m2]

   real(r8),parameter       :: &
      mlfr(1:2) = (/ 1./2., 1./3. /)                           ! fraction of mixed layer depth
   
   real(r8) :: step
   real(r8) :: xlist(5)                            
   real(r8) :: ylist(5)
!================================================================================
    
!****** Identify the water body type
      IF (dplak <= ShallowDeepPartition) THEN
         scwat = SHALLOW
      ELSE
         scwat = DEEP
      ENDIF


!****** Initialize the lake layer thickness and depth
      CALL LakeDisper( &  !-WMEJ Discretize the lake into stratified layers.
            ! "in" arguments
            ! -------------------
            nlake    , dplak    , scwat    ,&
            ! "inout" arguments
            ! -------------------
            zlake    , dzlak    , zilak    )


!****** Initialize the snow related variables
      tmsno = tskin


!****** Initialize the ice related variables
      icedp = 0.0     
      bicedp = 0.0    
      wicedp = 0.0    
      icefr = 0.0     
      tmice = tfrz    


!****** Initialize the lake temperature
      !+WMEJ TODO: IF possible, please initialize lake temperatures using
      !+           reanalysis data, as the current approach requires at least
      !+           one year of spin-up for the lake.
      lktmp(1) = tskin
      DO i = 2, nlake
         IF (zilak(i) > 50) THEN !+WMEJ: default setting in WRF-Lake
               lktmp(i) = lktmp(i-1)
         ELSE
               lktmp(i) = tskin + zilak(i)*(277. - tskin)/50. 
         ENDIF
         ! lktmp(i) = 10._r8 + 273.15_r8 !+WMEJ: default setting in CoLM mkinidata
      ENDDO


!****** Initialize the snow temperature
      !+WMEJ TODO: IF possible, please initialize snow temperatures using
      !+           reanalysis data, as the current approach requires at least
      !+           one year of spin-up for the lake botom soil.
      rhosnw = densnow
      lkrho = denh2o


!****** Initialize FLake special variables
      tmmnw = 0.0
      DO i = 1, nlake
         tmmnw = tmmnw + lktmp(i)*(dzlak(i)/dplak)    !+WMEJ Consider the proportion of thickness of each layer
      ENDDO
      tmwml = lktmp(1)
      tmbot = lktmp(nlake) 
      tmups = lktmp(nlake) 
      IF (abs(tmwml - tmbot) <= 1e-6) THEN !+WMEJ: Avoid division by zero, only influeces FLake
         tmbot = tmbot - 1.0
      ENDIF
      mldp  = mlfr(scwat) * dplak
      upsdp = mlfr(scwat) * dplak
      CTfrac = (tmwml-tmmnw)/(tmwml-tmbot) 
      CTfrac = MIN(C_T_max, MAX(CTfrac, C_T_min)) 
    

!****** Initialize Simstrat special variables
      uwatv = 1.0
      vwatv = 1.0
      lksal = 5.0
      tke = 3e-006
      eps = 3e-006
      etke = 0.
      num = 0.
      nuh = 0.
      !-WMEJ TODO: should be set using the bathymetric data,
      !            we dont consider the bathymetric.
      IF (scwat == SHALLOW) THEN
         ziarea = (/1.0000000000, 0.9999999999, 0.9999999998, 0.9999999997, 0.9999999996, 0.9999999995, 0.9999999994, 0.9999999993, 0.9999999992, 0.9999999991, 0.9999999990/)
      ELSE
         ziarea = (/1.0000000000, 0.9999999999, 0.9999999998, 0.9999999997, 0.9999999996, 0.9999999995, 0.9999999994, 0.9999999993, 0.9999999992, 0.9999999991, 0.9999999990/)
      ENDIF


!****** Initialize other variables
      stke1 = tkwat 

   END SUBROUTINE LakeIni



   SUBROUTINE LakeRst( &
               ! "in" arguments
               ! -------------------
               nlake     , nsnow     , nsoil     ,& 
               nlice     , scwat     ,&
               dplak     , snwdp     , icedp     ,&
               ! "inout" arguments
               ! -------------------
               zlake     , zilak     , dzlak     ,& 
               dzssb     , zssb      , zissb     ,& 
               snwcv     , snlay     )
   
! ---------------------------------- code history -------------------------------------
! Description:
!     +xum:  this one should be called after the cssptci, 
!            the function is to initilize the variables used in lake
!     +WMEJ: All schemes related to lake stratification should be retained here as options.
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_Lake_Const, only: tfrz, denh2o, denice, tkwat, densnow
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      nlake                   ,&! Maximum number of lake layer
      nsnow                   ,&! Maximum number of snow layer
      nsoil                   ,&! Maximum number of soil layer
      nlice                   ,&! Maximum number of ice layer
      scwat                     ! water body type
        
   real(r8), intent(in)     :: &
      dplak                   ,&! lake depth [m]
      icedp                   ,&! ice depth [m]
      snwdp                     ! snow depth [m]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      dzlak(nlake)            ,&! lake layer thickness [m]
      zlake(nlake)            ,&! lake layer depth [m]
      zilak(nlake+1)            ! lake layer interface level [m]

   real(r8), intent(out)    :: &
      dzssb(-nsnow+1:nsoil)   ,&! soil + snow layer thickness [m]
      zssb(-nsnow+1:nsoil)    ,&! soil + snow layer depth [m]
      zissb(-nsnow:nsoil+1)   ,&! soil + snow layer interface level [m]
      snwcv(-nsnow:nsoil+1)     ! snow cover fraction [-]

   integer, intent(out)     :: &
      snlay                     ! number of snow layers

!  ------------------------- local variables ---------------------------
   integer :: i                 ! loop index
!================================================================================
    
!****** Initialize the lake layer thickness and depth
      CALL LakeDisper( &  !-WMEJ Discretize the lake into stratified layers.
            ! "in" arguments
            ! -------------------
            nlake    , dplak    , scwat    ,&
            ! "inout" arguments
            ! -------------------
            zlake    , dzlak    , zilak    )


!****** Initialize the snow related variables
      CALL snowdp2lev4lake( & !-WMEJ Discretize the snow into stratified layers.
            ! "in" arguments
            ! -------------------
            nsnow            , snwdp           ,&
            ! "inout" arguments
            ! -------------------
            dzssb(-nsnow+1:0), zssb(-nsnow+1:0),& 
            zissb(-nsnow:0  ), snlay            ) 

!****** Initialize the lake bottom soil layer thickness and depth
      CALL LksoilLayer( &   !-WMEJ Discretize the lake bottom soil into stratified layers.
            ! "in" arguments
            ! -------------------
            nsoil     ,&
            ! "inout" arguments
            ! -------------------
            dzssb(1:nsoil), zssb(1:nsoil), zissb(1:nsoil+1) )

   END SUBROUTINE LakeRst



   SUBROUTINE LakeDisper ( &
               ! "in" arguments
               ! -------------------
               nlake    , dplak    , scwat    ,&
               ! "inout" arguments
               ! -------------------
               zlake    , dzlak    , zilak    )

! ------------------------ code history -------------------------------------
!
! Description:
!     Discretize the lake into stratified layers.
!     This SUBROUTINE is primarily designed to prepare for
!     future implementations of lake water level changes.
!
! Original author: 
!     Omarjan Obulkasim, 04/2024
!----------------------------------------------------------------------------
   USE MOD_Namelist, only: DEF_LAKE_LAYER_SCHEME
   USE MOD_Lake_Const, only: DZ_COLML, DZ_CSSPL, COLMLLayer,&
                                    CSSPLLayer, EqualLayer, XMlogLayer, WMlogLayer
!================================================================================
!  ------------------------- input variables ---------------------------
   integer , intent(in   ) :: nlake            ! Maximum number of lake layer
   integer , intent(in   ) :: scwat            ! water body type
   real(r8), intent(in   ) :: dplak            ! lake depth [m]
   real(r8), intent(inout) :: zlake(nlake)     ! lake layer depth [m]
   real(r8), intent(inout) :: dzlak(nlake)     ! lake layer thickness [m]
   real(r8), intent(inout) :: zilak(nlake+1)   ! lake layer interface level [m]

!  ------------------------- local variables ---------------------------
   integer :: i, layernum
   integer :: LayerScheme
   real(r8) :: rmdplak, sumrats
   real(r8) :: base, start_node, exp_value, strec
!==============================================================================
      !- check the lake layer scheme
      LayerScheme = DEF_LAKE_LAYER_SCHEME
      SELECT CASE (DEF_LAKE_LAYER_SCHEME)
         CASE(1)
            layernum = size(DZ_COLML(:,scwat), dim=1)
            IF (nlake > layernum .and. nlake > 0) THEN
               LayerScheme = WMlogLayer
            ELSE
               LayerScheme = 1
            ENDIF
         CASE(2)
            layernum = size(DZ_CSSPL(:,scwat), dim=1)
            IF (nlake > layernum .and. nlake > 0) THEN
               LayerScheme = WMlogLayer
            ELSE
               LayerScheme = 2
            ENDIF
      END SELECT


      SELECT CASE (LayerScheme)
         !------------------------------
         CASE (COLMLLayer)  !- CoLML layer scheme
         !------------------------------
            !--- Calculate the lake layer thickness 
            dzlak(1) = DZ_COLML(1,scwat)                       !- The top layer thickness is set to 0.1 m
            sumrats = dplak / sum(DZ_COLML(1:nlake, scwat))
            dzlak(2:nlake-1) = DZ_COLML(2:nlake-1, scwat)*sumrats
            dzlak(nlake) = DZ_COLML(nlake, scwat)*sumrats - (dzlak(1) - DZ_COLML(1, scwat)*sumrats)

            !--- Calculate the lake layer interface depth
            zilak(1) = 0.0_r8
            DO i = 2, nlake
               zilak(i) = zilak(i-1) + dzlak(i-1)
            ENDDO
            zilak(nlake+1) = dplak

            !--- Calculate the lake layer node depth
            DO i = 1, nlake
               zlake(i) = zilak(i) + 0.5_r8 * dzlak(i)
            ENDDO

         !------------------------------
         CASE (CSSPLLayer)  !- CSSPL layer scheme
         !------------------------------
            !--- Calculate the lake layer thickness 
            sumrats = sum(DZ_CSSPL(1:nlake,scwat))
            DO i = 1, nlake
               dzlak = dplak * DZ_CSSPL(i,scwat) / sumrats
            ENDDO

            !--- Calculate the lake layer interface depth
            zilak(1) = 0.0_r8
            DO i = 2, nlake
               zilak(i) = zilak(i-1) + dzlak(i-1)
            ENDDO
            zilak(nlake+1) = dplak

            !--- Calculate the lake layer node depth
            DO i = 1, nlake
               zlake(i) = zilak(i) + 0.5_r8 * dzlak(i)
            ENDDO

         !------------------------------
         CASE (EqualLayer)  !- equal layer scheme
         !------------------------------
            !--- Calculate the lake layer thickness 
            dzlak = dplak / nlake

            !--- Calculate the lake layer interface depth
            zilak(1) = 0.0_r8
            DO i = 2, nlake
               zilak(i) = zilak(i-1) + dzlak(i-1)
            ENDDO
            zilak(nlake+1) = dplak
            !--- Calculate the lake layer node depth
            DO i = 1, nlake
               zlake(i) = zilak(i) + 0.5_r8 * dzlak(i)
            ENDDO

         !------------------------------
         CASE (XMlogLayer)  !- MIN XU, log layer scheme
         !------------------------------
            !--- Calculate the lake layer interface depth
            strec = -4.
            DO i = 1, nlake
               zilak(i) = real(i) * (1.0 / real(nlake))
               zilak(i) = dplak / strec * log(1.0 - zilak(i) * (1.0 - exp(strec)))
            ENDDO
            zilak(nlake+1) = dplak

            !--- Calculate the lake layer thickness
            DO i = 1, nlake
               dzlak(i) = zilak(i+1) - zilak(i)
               zlake(i) = 0.5 * (zilak(i+1) + zilak(i))
            ENDDO

         !------------------------------
         CASE (WMlogLayer) !- WMEJ, log layer scheme
         !------------------------------
            base = 1.4_r8
            start_node = 0.05_r8
            !--- Calculate the lake layer node depth
            zlake(1) = start_node
            DO i = 2, nlake
               exp_value = exp((i-1) * log(dplak / start_node) / (real(nlake) - 1))
               zlake(i) = start_node * exp_value
            ENDDO

            !--- Calculate the lake layer interface depth
            zilak(1) = 0.0_r8  
            DO i = 1, nlake-1
               zilak(i+1) = 0.5_r8 * (zlake(i) + zlake(i+1))
            ENDDO
            zilak(nlake+1) = dplak  

            !--- Calculate the lake layer thickness
            dzlak(1) = zilak(2)
            DO i = 2, nlake
               dzlak(i) = zilak(i+1) - zilak(i)
            ENDDO

         !------------------------------
         CASE DEFAULT !- default equal layer scheme
         !------------------------------
            !--- Calculate the lake layer thickness 
            dzlak = dplak / nlake

            !--- Calculate the lake layer interface depth
            zilak(1) = 0.0_r8
            DO i = 2, nlake
               zilak(i) = zilak(i-1) + dzlak(i-1)
            ENDDO
            zilak(nlake+1) = dplak

            !--- Calculate the lake layer node depth
            DO i = 1, nlake
               zlake(i) = zilak(i) + 0.5_r8 * dzlak(i)
            ENDDO

      END Select

   END SUBROUTINE LakeDisper



   SUBROUTINE LakeIceFrc( &
               ! "in" arguments
               ! -------------------
               nlake     , dzlak     , zilak     , icedp     ,&
               ! "inout" arguments
               ! -------------------
               wice_lake , wliq_lake , icefr     )

! ------------------------ code history -------------------------------------
! Description:
!     Calculate the ice fraction in the lake
!
! Original author:
!     Omarjan Obulkasim, 04/2024
!----------------------------------------------------------------------------
   USE MOD_Lake_Const, only : denh2o, denice
!================================================================================
!  ------------------------- input variables ---------------------------
   integer , intent(in)     :: nlake
   real(r8), intent(in)     :: icedp
   real(r8), intent(in)     :: dzlak(nlake), zilak(nlake+1)
   ! ----- input/output variables -----
   real(r8), intent(inout)  :: icefr(nlake)
   real(r8)                 :: wice_lake(nlake), wliq_lake(nlake)
   integer                  :: i
!================================================================================

#ifdef WMEJICE 
      icefr = 0.
      IF (icedp < 0.00001) THEN 
         DO i = 1, nlake
         icefr(i) = 0.
         wice_lake(i) = 0.
         wliq_lake(i) = denh2o*dzlak(i)
         ENDDO        
      ELSE
         DO i = 1, nlake
            IF (icedp > zilak(i+1)) THEN
               icefr(i) = 1.
               wice_lake(i) = denice*dzlak(i)
               wliq_lake(i) = denh2o*dzlak(i)
            ELSE 
               wice_lake(i) = denice*(icedp - zilak(i))
               wliq_lake(i) = denh2o*dzlak(i)
               icefr(i) = wice_lake(i)/(wice_lake(i) + wliq_lake(i))
            ENDIF
         ENDDO
      ENDIF
#else
      !-GLERL Scheme, Maroni et al. (2010) also state 10 cm for whole ice cov
      ! Here, the ratio of ice depth to water layer depth is used to calculate the ice cover fraction.
      icefr = 0.
      IF (icedp < 0.00001) THEN 
         DO i = 1, nlake
            icefr(i) = 0.
         ENDDO        
      ELSE
         DO i = 1, nlake
            IF (icedp > zilak(i+1)) THEN
               icefr(i) = 1.
            ELSE 
               icefr(i) = (icedp - zilak(i))/dzlak(i)
            ENDIF
         ENDDO
      ENDIF
#endif

   END SUBROUTINE LakeIceFrc



   SUBROUTINE LksoilLayer( &
               ! "in" arguments
               ! -------------------
               nsoil     ,&
               ! "inout" arguments
               ! -------------------
               dzssb     , zssb      , zissb )

! ------------------------ code history -------------------------------------
!
! Description:
!     Discretize the lake bottom soil into stratified layers.
!     Soil thickness does not change over time.
!
! Original author: 
!     Omarjan Obulkasim, 04/2024
!----------------------------------------------------------------------------
!================================================================================
!  ------------------------- input variables ---------------------------
   integer , intent(in   )  :: nsoil
   ! ----- input/output variables -----
   real(r8), intent(inout)  :: dzssb(nsoil), zissb(nsoil+1), zssb(nsoil)
   integer                  :: j
!================================================================================

      !--- Calculate the soil layer node depth
      DO j = 1,nsoil
         zssb(j) = 0.025*(exp(0.5*(float(j)-0.5))-1.0)     ! node depths
      ENDDO

      !--- Calculate the soil layer thickness
      dzssb(1) = 0.5*(zssb(1) + zssb(2))
      j = nsoil
      dzssb(j) = zssb(j)-zssb(j-1)
      DO j = 2, nsoil-1
         dzssb(j) = 0.5*(zssb(j+1) - zssb(j-1))
      ENDDO

      !--- Calculate the soil layer interface depth
      zissb(1) = 0.
      j = nsoil+1
      zissb(j) = zssb(j) + 0.5*dzssb(j)
      DO j = 1, nsoil
         zissb(j) = 0.5*(zssb(j) + zssb(j+1))
      ENDDO

   END SUBROUTINE LksoilLayer



   SUBROUTINE snowage( &
               ! "in" arguments
               ! -------------------
               dtlak     , tskin    , scv    , scvold   ,&
               ! "inout" arguments
               ! -------------------
               snwag      ) 

! ------------------------ code history -------------------------------------
! Description:
!     Update snow cover and snow age, based on BATS code
!
! Original author:  
!     Robert Dickinson
!----------------------------------------------------------------------------
   USE MOD_Lake_Const, ONLY: tfrz
!================================================================================
!  ------------------------- input variables ---------------------------
   real(r8), intent(in) :: dtlak    ! time step [s]
   real(r8), intent(in) :: tskin    ! temperature of soil at surface [K]
   real(r8), intent(in) :: scv      ! snow cover, water equivalent [mm]
   real(r8), intent(in) :: scvold   ! snow cover for previous time step [mm]
   real(r8), intent(inout) :: snwag ! non dimensional snow age [-]
!  ------------------------- local variables ---------------------------
   real(r8) :: age1                 ! snow aging factor due to crystal growth [-]
   real(r8) :: age2                 ! snow aging factor due to surface growth [-]
   real(r8) :: age3                 ! snow aging factor due to accum of other particles [-]
   real(r8) :: arg                  ! temporary variable used in snow age calculation [-]
   real(r8) :: arg2                 ! temporary variable used in snow age calculation [-]
   real(r8) :: dela                 ! temporary variable used in snow age calculation [-]
   real(r8) :: dels                 ! temporary variable used in snow age calculation [-]
   real(r8) :: sge                  ! temporary variable used in snow age calculation [-]
!================================================================================

      IF (scv <= 0.) THEN
         snwag = 0.
      ! Over antarctica
      elseif (scv > 800.) THEN
         snwag = 0.
         ! Away from antarctica
      ELSE
         age3  = 0.3
         arg   = 5.e3*(1./tfrz-1./tskin)
         arg2  = min(0.,10.*arg)
         age2  = exp(arg2)
         age1  = exp(arg)
         dela  = 1.e-6*dtlak*(age1+age2+age3)
         dels  = 0.1*max(0.0,scv-scvold)
         sge   = (snwag+dela)*(1.0-dels)
         snwag   = max(0.0,sge)
      ENDIF

   END SUBROUTINE snowage



   SUBROUTINE snowdp2lev4lake( &
               ! "in" arguments
               ! -------------------
               nsnow     , snowdp    ,&
               ! "inout" arguments
               ! -------------------
               dzsno     , zsno      ,& 
               zisno     , snlay      )

! ---------------------------------- code history -------------------------------------
! DESCRIPTION:
!     Discretize the snow into stratified layers.
!
! Original author : 
!     Min Xu 
!
! Revisions:
!     xum: 2014/04/01: fixed a bug in CALL snowlayersdivide, 
!                  the lb should be set to -nsnow+1 instead of wrong one. snlay + 1 
!
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Updated input and output variables
!--------------------------------------------------------------------------------------
!================================================================================
!  ------------------------- input variables ---------------------------
   ! integer ,intent(in   ) :: nlake           ! number of lake layer
   integer ,intent(in   ) :: nsnow             ! number of snow layer  !positive
   real(r8),intent(in   ) :: snowdp            ! snow depth (m)
!  ------------------------- inout variables ---------------------------
   real(r8),intent(inout) :: dzsno(-nsnow+1:0) ! snow layer thichness (m)
   real(r8),intent(inout) :: zsno(-nsnow+1:0)  ! snow layer depth (m)
   real(r8),intent(inout) :: zisno(-nsnow:0)   ! interface level below a "z" level (m)
   integer ,intent(inout) :: snlay             ! number of snow layers
!  ------------------------- local variables ---------------------------
   integer  :: i, j                            ! loop counter
!================================================================================
      ! ---- Snow Layer ----
      snlay = 0
      zsno(-nsnow+1:0) = 0.
      zisno(-nsnow:0) = 0.
      dzsno(-nsnow+1:0) = 0.
      IF (snowdp .ge. 0.01) THEN
         IF ((snowdp .ge. 0.01) .and. (snowdp .le. 0.03)) THEN
            snlay = -1
            dzsno(0) = snowdp
         ELSE IF ((snowdp .gt. 0.03) .and. (snowdp .le. 0.04)) THEN
            snlay = -2
            dzsno(-1) = snowdp/2.
            dzsno(0)  = snowdp/2.
         ELSE IF ((snowdp .gt. 0.04) .and. (snowdp .le. 0.07)) THEN 
            snlay = -2
            dzsno(-1) = 0.02
            dzsno(0)  = snowdp - dzsno(-1)
         ELSE IF ((snowdp .gt. 0.07) .and. (snowdp .le. 0.12)) THEN
            snlay = -3
            dzsno(-2) = 0.02
            dzsno(-1) = (snowdp - 0.02) / 2.
            dzsno(0)  = dzsno(-1)
         ELSE IF ((snowdp .gt. 0.12) .and. (snowdp .le. 0.18)) THEN
            snlay = -3
            dzsno(-2) = 0.02
            dzsno(-1) = 0.05
            dzsno(0)  = snowdp - dzsno(-2) - dzsno(-1)
         ELSE IF ((snowdp .gt. 0.18) .and. (snowdp .le. 0.29)) THEN
            snlay = -4
            dzsno(-3) = 0.02
            dzsno(-2) = 0.05
            dzsno(-1) = (snowdp - dzsno(-3) - dzsno(-2))/2._r8
            dzsno( 0) = dzsno(-1)
         ELSE IF ((snowdp .gt. 0.29) .and. (snowdp .le. 0.41)) THEN
            snlay = -4
            dzsno(-3) = 0.02
            dzsno(-2) = 0.05
            dzsno(-1) = 0.11
            dzsno(0)  = snowdp - dzsno(-3) - dzsno(-2) - dzsno(-1)
         ELSE IF ((snowdp .gt. 0.41) .and. (snowdp .le. 0.64)) THEN
            snlay = -5 
            dzsno(-4) = 0.02
            dzsno(-3) = 0.05
            dzsno(-2) = 0.11
            dzsno(-1) = (snowdp - dzsno(-4) - dzsno(-3) - dzsno(-2))/2._r8
            dzsno(0)  = dzsno(-1)
         ELSE IF (snowdp .gt. 0.64) THEN
            snlay = -5
            dzsno(-4) = 0.02
            dzsno(-3) = 0.05
            dzsno(-2) = 0.11
            dzsno(-1) = 0.23
            dzsno(0)  = snowdp -dzsno(-4)-dzsno(-3)-dzsno(-2)-dzsno(-1)
         ENDIF
      ENDIF
      zisno(0) = 0.
      DO j = 0, snlay+1, -1
         zsno(j  ) = zisno(j) - 0.5_r8 * dzsno(j)
         zisno(j-1) = zisno(j) - dzsno(j)
      ENDDO

   END SUBROUTINE snowdp2lev4lake



   SUBROUTINE LinearInterp(numPoints, x1, y1, x2, y2, xlist, ylist)
! ------------------------ code history -------------------------------------
! Description:
!     Linear interpolation between two points
!
! Original author:  
!     Omarjan Obulkasim, 07/2024 : Following the GLM scheme (Hipsey et al., 2019)
!----------------------------------------------------------------------------
!================================================================================
!  ------------------------- input variables ---------------------------
   integer , intent(in)    :: numPoints          ! number of points  
   real(r8), intent(in)    :: x1, y1             ! x1
   real(r8), intent(in)    :: x2, y2             ! x2
   real(r8), intent(in)    :: xlist(numPoints)   ! the list of x values need to be interpolated
   real(r8), intent(out)   :: ylist(numPoints)   ! the list of y values interpolated
!  ------------------------- local variables ---------------------------
   real(r8) :: slope
   real(r8) :: intercept
   integer  :: i
!================================================================================

      slope = (y2 - y1) / (x2 - x1)
      intercept = y1 - slope * x1
      DO i = 1, numPoints
         ylist(i) = slope * xlist(i) + intercept
      ENDDO

   END SUBROUTINE LinearInterp



   ! Reverse an array without allocating a second array
   SUBROUTINE reverse_array(in_arr)
      real(r8), intent(inout) :: in_arr(:)
      real(r8) :: temp

      integer :: first, last, i, len

      first = lbound(in_arr, dim=1)
      last = ubound(in_arr, dim=1)
      len = size(in_arr)

      ! Works for even and odd sized arrays
      !(as len/2 is always integer and not rounded, but cutoff)
      DO i = last, first + int(len/2), -1
         temp = in_arr(i)
         in_arr(i) = in_arr(len + 1 - i)
         in_arr(len + 1 - i) = temp
      ENDDO

   END SUBROUTINE reverse_array


END MODULE  MOD_Lake_Utils