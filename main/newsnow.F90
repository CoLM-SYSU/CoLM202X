
 subroutine newsnow (patchtype,maxsnl,deltim,t_grnd,pg_rain,pg_snow,bifall,& 
                     t_precip,zi_soisno,z_soisno,dz_soisno,t_soisno,&
                     wliq_soisno,wice_soisno,fiold,snl,sag,scv,snowdp,fsno)

!=======================================================================
! add new snow nodes. 
! Original author : Yongjiu Dai, 09/15/1999; 08/31/2002, 07/2013, 04/2014
!=======================================================================
!
  use precision
  use PhysicalConstants, only : tfrz, cpliq, cpice
 
  implicit none
 
! ------------------------ Dummy Argument ------------------------------

  integer, INTENT(in) :: maxsnl     ! maximum number of snow layers
  integer, INTENT(in) :: patchtype  ! land water type (0=soil, 1=urban and built-up,
                                    ! 2=wetland, 3=land ice, 4=land water bodies, 99=ocean)
  real(r8), INTENT(in) :: deltim    ! model time step [second]
  real(r8), INTENT(in) :: t_grnd    ! ground surface temperature [k]
  real(r8), INTENT(in) :: pg_rain   ! rainfall onto ground including canopy runoff [kg/(m2 s)]
  real(r8), INTENT(in) :: pg_snow   ! snowfall onto ground including canopy runoff [kg/(m2 s)]
  real(r8), INTENT(in) :: bifall    ! bulk density of newly fallen dry snow [kg/m3]
  real(r8), INTENT(in) :: t_precip  ! snowfall/rainfall temperature [kelvin]

  real(r8), INTENT(inout) ::    zi_soisno(maxsnl:0)   ! interface level below a "z" level (m)
  real(r8), INTENT(inout) ::     z_soisno(maxsnl+1:0) ! layer depth (m)
  real(r8), INTENT(inout) ::    dz_soisno(maxsnl+1:0) ! layer thickness (m)
  real(r8), INTENT(inout) ::     t_soisno(maxsnl+1:0) ! soil + snow layer temperature [K]
  real(r8), INTENT(inout) ::  wliq_soisno(maxsnl+1:0) ! liquid water (kg/m2)
  real(r8), INTENT(inout) ::  wice_soisno(maxsnl+1:0) ! ice lens (kg/m2)
  real(r8), INTENT(inout) :: fiold(maxsnl+1:0) ! fraction of ice relative to the total water
   integer, INTENT(inout) :: snl               ! number of snow layers
  real(r8), INTENT(inout) :: sag               ! non dimensional snow age [-]
  real(r8), INTENT(inout) :: scv               ! snow mass (kg/m2)
  real(r8), INTENT(inout) :: snowdp            ! snow depth (m)
  real(r8), INTENT(inout) :: fsno              ! fraction of soil covered by snow [-]

! ----------------------- Local  Variables -----------------------------

  real(r8) dz_snowf  ! layer thickness rate change due to precipitation [mm/s]
  integer newnode    ! signification when new snow node is set, (1=yes, 0=no)
  integer lb

!-----------------------------------------------------------------------
      newnode = 0

      dz_snowf = pg_snow/bifall
      snowdp = snowdp + dz_snowf*deltim
      scv = scv + pg_snow*deltim            ! snow water equivalent (mm)

      if(patchtype==2 .AND. t_grnd>tfrz)then  ! snowfall on warmer wetland
         scv=0.; snowdp=0.; sag=0.; fsno = 0.
      endif

      zi_soisno(0) = 0.

! when the snow accumulation exceeds 10 mm, initialize a snow layer

      if(snl==0 .AND. pg_snow>0.0 .AND. snowdp>=0.01)then  
         snl = -1
         newnode = 1
         dz_soisno(0)  = snowdp             ! meter
         z_soisno (0)  = -0.5*dz_soisno(0)
         zi_soisno(-1) = -dz_soisno(0)

         sag = 0.                           ! snow age
! 03/16/2020, yuan: TODO, why not a combination of T 
         t_soisno (0) = min(tfrz, t_precip) ! K
         wice_soisno(0) = scv               ! kg/m2
         wliq_soisno(0) = 0.                ! kg/m2
         fiold(0) = 1.
! 03/16/2020, yuan: TODO, why don't apply the snowfall on snow pack case 
         fsno = min(1.,tanh(0.1*pg_snow*deltim))
         !fsno = 1. - (1. - tanh(0.1*pg_snow*deltim))*(1. - fsno)
         !fsno = min(1., fsno)
!**      write(6,*) 'snow layer is built'
      endif

      ! --------------------------------------------------
      ! snowfall on snow pack
      ! --------------------------------------------------
      ! the change of ice partial density of surface node due to precipitation
      ! only ice part of snowfall is added here, the liquid part will be added latter

      if(snl<0 .AND. newnode==0)then
         lb = snl + 1
         t_soisno(lb) = tfrz &
            + ( (wice_soisno(lb)*cpice+wliq_soisno(lb)*cpliq)*(t_soisno(lb)-tfrz) &
            +   (pg_rain*cpliq + pg_snow*cpice)*deltim*(t_precip-tfrz) ) &
            / ( wice_soisno(lb)*cpice + wliq_soisno(lb)*cpliq &
            +   pg_rain*deltim*cpliq + pg_snow*deltim*cpice )

         !print *, "new snow:", t_soisno(lb)    ! fordebug
         t_soisno(lb) = min(tfrz, t_soisno(lb))
         wice_soisno(lb) = wice_soisno(lb)+deltim*pg_snow
         dz_soisno(lb) = dz_soisno(lb)+dz_snowf*deltim
         z_soisno(lb) = zi_soisno(lb) - 0.5*dz_soisno(lb)
         zi_soisno(lb-1) = zi_soisno(lb) - dz_soisno(lb)

         ! update fsno by new snow event, add to previous fsno
         ! shape factor for accumulation of snow = 0.1
         fsno = 1. - (1. - tanh(0.1*pg_snow*deltim))*(1. - fsno)
         fsno = min(1., fsno)

      endif

 end subroutine newsnow
