#include <define.h>

MODULE SNOW_Layers_CombineDivide

!-----------------------------------------------------------------------
 use precision
 IMPLICIT NONE
 SAVE

! PUBLIC MEMBER FUNCTIONS:
  public :: snowcompaction
  public :: snowlayerscombine
  public :: snowlayersdivide


! PRIVATE MEMBER FUNCTIONS:
  private :: combo


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------



 subroutine snowcompaction (lb,deltim,imelt,fiold,t_soisno,wliq_soisno,wice_soisno,dz_soisno)
 
!=======================================================================
! Original author: Yongjiu Dai, September 15, 1999
!
! three of metamorphisms of changing snow characteristics are implemented,
! i.e., destructive, overburden, and melt. The treatments of the former two
! are from SNTHERM.89 and SNTHERM.99 (1991, 1999). The contribution due to
! melt metamorphism is simply taken as a ratio of snow ice fraction after
! the melting versus before the melting.
!
!=======================================================================

  use precision
  use PhysicalConstants, only : denice, denh2o, tfrz
  implicit none
!
!-------------------------- Dummy argument -----------------------------
!
  integer, INTENT(in) :: &
        lb             ! lower bound of array

  integer, INTENT(in) :: & 
        imelt(lb : 0)  ! signifies if node in melting (imelt = 1)

  real(r8), INTENT(in) :: &
        deltim,       &! seconds i a time step [second]
        fiold(lb : 0),&! fraction of ice relative to 
                       ! the total water content at the previous time step 
        t_soisno (lb : 0),   &! nodal temperature [K]
        wice_soisno(lb : 0), &! ice lens [kg/m2]
        wliq_soisno(lb : 0)   ! liquid water [kg/m2]

  real(r8), INTENT(inout) :: &
        dz_soisno  (lb : 0)   ! layer thickness [m]
!
!----------------------- local variables ------------------------------
  integer i            ! Numeber of doing loop
  
  real(r8) :: &
       c1,            &! = 2.777e-7 [m2/(kg s)]
       c2,            &! = 21e-3 [m3/kg]
       c3,            &! = 2.777e-6 [1/s]
       c4,            &! = 0.04 [1/K]
       c5,            &! = 2.0 
       c6,            &! = 5.15e-7.
       c7,            &! = 4.
       dm,            &! Upper Limit on Destructive Metamorphism Compaction [kg/m3]
       eta0            ! The Viscosity Coefficient Eta0 [kg-s/m2]

  real(r8) :: &
       burden,        &! pressure of overlying snow [kg/m2]
       ddz1,          &! Rate of settling of snowpack due to destructive metamorphism.
       ddz2,          &! Rate of compaction of snowpack due to overburden.
       ddz3,          &! Rate of compaction of snowpack due to melt [1/s]
       dexpf,         &! expf=exp(-c4*(273.15-t_soisno)).
       fi,            &! Fraction of ice relative to the total water content 
                       ! at the current time step
       td,            &! t_soisno - tfrz [K]
       pdzdtc,        &! Nodal rate of change in fractional-thickness 
                       ! due to compaction [fraction/s]
       void,          &! void (1 - vol_ice - vol_liq)
       wx,            &! water mass (ice+liquid) [kg/m2]
       bi              ! partitial density of ice [kg/m3]

  data c2,c3,c4,c5/23.e-3, 2.777e-6, 0.04, 2.0/
  data c6/5.15e-7/, c7/4./
  data dm/100./         
  data eta0/9.e5/      

!=======================================================================

      burden = 0.0

      do i = lb, 0
         wx = wice_soisno(i) + wliq_soisno(i)
         void = 1.- (wice_soisno(i)/denice + wliq_soisno(i)/denh2o)/dz_soisno(i)

! Disallow compaction for water saturated node and lower ice lens node.
         if(void <= 0.001 .or. wice_soisno(i) <= .1)then
            burden = burden+wx
            CYCLE
         endif

         bi = wice_soisno(i) / dz_soisno(i)
         fi = wice_soisno(i) / wx
         td = tfrz-t_soisno(i)
         
         dexpf = exp(-c4*td)

! Settling as a result of destructive metamorphism
         ddz1 = -c3*dexpf
         if(bi > dm) ddz1 = ddz1*exp(-46.0e-3*(bi-dm))

! Liquid water term
         if(wliq_soisno(i) > 0.01*dz_soisno(i)) ddz1=ddz1*c5

! Compaction due to overburden
         ddz2 = -burden*exp(-0.08*td - c2*bi)/eta0

! Compaction occurring during melt
         if(imelt(i) == 1)then
            ddz3 = - 1./deltim * max(0.,(fiold(i) - fi)/fiold(i))
         else
            ddz3 = 0.
         endif

! Time rate of fractional change in dz_soisno (units of s-1)
         pdzdtc = ddz1+ddz2+ddz3

! The change in dz_soisno due to compaction
         dz_soisno(i) = dz_soisno(i)*(1.+pdzdtc*deltim)

! Pressure of overlying snow
         burden = burden+wx

      end do

 end subroutine snowcompaction



 subroutine snowlayerscombine (lb,snl, &
            z_soisno,dz_soisno,zi_soisno,wliq_soisno,wice_soisno,t_soisno,scv,snowdp)

!=======================================================================
! Original author : Yongjiu Dai, September 15, 1999
!
! checks for elements which are below prescribed minimum for thickness or mass.
! If snow element thickness or mass is less than a prescribed minimum,
! it is combined with neighboring element to be best combine with,
! and executes the combination of mass and energy in clm_combo.f90
! 
!=======================================================================

  use precision
  implicit none

!-------------------------- Dummy argument -----------------------------
  integer, INTENT(in) :: lb               ! lower bound of array

! numbering from 1 (bottom) mss (surface)
  real(r8), INTENT(inout) :: wice_soisno(lb:1)   ! ice lens [kg/m2]
  real(r8), INTENT(inout) :: wliq_soisno(lb:1)   ! liquid water {kg/m2]
  real(r8), INTENT(inout) :: t_soisno (lb:1)     ! nodel temperature [K]
  real(r8), INTENT(inout) :: dz_soisno  (lb:1)   ! layer thickness [m]
  real(r8), INTENT(inout) :: z_soisno   (lb:1)   ! node depth [m]
  real(r8), INTENT(inout) :: zi_soisno  (lb-1:1) ! depth of layer interface [m]
  real(r8), INTENT(inout) :: snowdp       ! snow depth [m]
  real(r8), INTENT(inout) :: scv          ! snow mass - water equivalent [kg/m2]
  integer, INTENT(inout) :: snl           ! Number of snow

!----------------------- Local variables ------------------------------
  real(r8) :: drr           ! thickness of the combined [m]
  real(r8) :: dzmin(5)      ! minimum of snow layer 1 (top) to msn0 (bottom)
  real(r8) :: zwice         ! total ice mass in snow
  real(r8) :: zwliq         ! total liquid water in snow

  integer :: i              ! number of do looping
  integer :: j              ! node index
  integer :: k              ! number of do looping
  integer :: l              ! node index
  integer :: msn_old        ! number of snow layer 1 (top) to msn0 (bottom)
  integer :: mssi           ! node index
  integer :: neibor         ! adjacent node selected for combination

  data dzmin /0.010, 0.015, 0.025, 0.055, 0.115/

!-----------------------------------------------------------------------
! check the mass of ice lens of snow, when the total less than a small value,
! combine it with the underlying neighbor
      msn_old = snl
      do j = msn_old+1, 0
         if(wice_soisno(j) <= .1)then
            wliq_soisno(j+1) = wliq_soisno(j+1) + wliq_soisno(j)
            wice_soisno(j+1) = wice_soisno(j+1) + wice_soisno(j)

! shift all elements above this down one.
            if(j > snl+1 .AND. snl < -1)then
               do i =  j, snl+2, -1
                  t_soisno(i) = t_soisno(i-1)
                  wliq_soisno(i) = wliq_soisno(i-1)
                  wice_soisno(i) = wice_soisno(i-1)
                  dz_soisno(i) = dz_soisno(i-1)
               enddo
            endif

            snl = snl + 1
!*          write(6,*) 'one snow layer is gone'

         endif

      enddo

      if(snl == 0)then
         scv = 0.
         snowdp = 0.
!*       write(6,*) 'all snow has gone'
         return
      else
         scv = 0.
         snowdp = 0.
         zwice = 0.
         zwliq = 0.
         do j = snl + 1, 0
            scv = scv + wice_soisno(j) + wliq_soisno(j)
            snowdp = snowdp + dz_soisno(j)
            zwice = zwice + wice_soisno(j)
            zwliq = zwliq + wliq_soisno(j)
         enddo
      endif
!-----------------------------------------------------------------------
! check the snow depth

      if(snowdp < 0.01)then       !!! all snow gone 
          
         snl = 0
         scv = zwice
         if(scv <= 0.) snowdp = 0.

! the liquid water assumed ponding on soil surface
         wliq_soisno(1) = wliq_soisno(1) + zwliq
!*       write(6,'(17h all snow is gone)')
         return

      else                        !!! snow layers combined

! two or more layers 

        if(snl < -1)then
           msn_old = snl
           mssi = 1
           do i = msn_old+1, 0

! If top node is removed, combine with bottom neighbor
              if(dz_soisno(i) < dzmin(mssi))then
                 if(i == snl+1)then
                    neibor = i + 1

! If the bottom neighbor is not snow, combine with the top neighbor
                 else if(i == 0)then
                    neibor = i - 1

! If none of the above special cases apply, combine with the thinnest neighbor
                 else
                   neibor = i + 1
                   if((dz_soisno(i-1)+dz_soisno(i)) < (dz_soisno(i+1)+dz_soisno(i))) neibor = i-1
                 endif

! Node l and j are combined and stored as node j.

                 if(neibor > i)then
                    j = neibor
                    l = i
                 else
                    j = i
                    l = neibor
                 endif
                 call combo ( dz_soisno(j), wliq_soisno(j), wice_soisno(j), t_soisno(j),&
                              dz_soisno(l), wliq_soisno(l), wice_soisno(l), t_soisno(l) )

! Now shift all elements above this down one.

                 if(j-1 > snl+1) then
                    do k = j-1, snl+2, -1
                       t_soisno(k) = t_soisno(k-1)
                       wice_soisno(k) = wice_soisno(k-1)
                       wliq_soisno(k) = wliq_soisno(k-1)
                       dz_soisno(k) = dz_soisno(k-1)
                    enddo
                 endif

                 snl = snl + 1

!*    write(6,'(7h Nodes ,i4,4h and,i4,14h combined into,i4)') l,j,j

                 if(snl >= -1) EXIT

! The layer thickness great than the prescibed minimum value

              else
                 mssi = mssi + 1 
              endif
           enddo

        end if

! Reset the node depth and the depth of layer interface

        zi_soisno(0) = 0.
        do k = 0, snl+1, -1
           z_soisno(k) = zi_soisno(k) - 0.5*dz_soisno(k)
           zi_soisno(k-1) = zi_soisno(k) - dz_soisno(k)
        enddo

      endif                       !!! snow layers combined 

 end subroutine snowlayerscombine



 subroutine snowlayersdivide (lb,snl,z_soisno,dz_soisno,zi_soisno,wliq_soisno,wice_soisno,t_soisno)

!=======================================================================
! Original author : Yongjiu Dai, September 15, 1999
!
! subdivides snow layer when its thickness exceed the prescribed maximum
!=======================================================================

  use precision
  implicit none

!-------------------------- Dummy argument -----------------------------

   integer, INTENT(in) :: lb              ! lower bound of array
   integer, INTENT(inout) :: snl          ! Number of snow
  real(r8), INTENT(inout) :: wice_soisno(lb:0)   ! ice lens [kg/m2]
  real(r8), INTENT(inout) :: wliq_soisno(lb:0)   ! liquid water [kg/m2]
  real(r8), INTENT(inout) :: t_soisno   (lb:0)   ! Nodel temperature [K]
  real(r8), INTENT(inout) :: dz_soisno  (lb:0)   ! Layer thickness [m]
  real(r8), INTENT(inout) :: z_soisno   (lb:0)   ! Node depth [m]
  real(r8), INTENT(inout) :: zi_soisno  (lb-1:0) ! Depth of layer interface [m]

!----------------------- Local variables ------------------------------

! numbering from 1 (surface) msno (bottom)
  real(r8) :: drr      ! thickness of the combined [m]
  real(r8) :: dzsno(5) ! Snow layer thickness [m] 
  real(r8) :: swice(5) ! Partial volume of ice [m3/m3]
  real(r8) :: swliq(5) ! Partial volume of liquid water [m3/m3]
  real(r8) :: tsno(5)  ! Nodel temperature [K]

  integer k            ! number of do looping
  integer msno         ! number of snow layer 1 (top) to msno (bottom)

  real(r8) zwice,zwliq,propor

!-----------------------------------------------------------------------

      msno = abs(snl)
      do k = 1, msno
         dzsno(k) = dz_soisno  (k + snl)
         swice(k) = wice_soisno(k + snl)
         swliq(k) = wliq_soisno(k + snl)
         tsno(k)  = t_soisno (k + snl)
      enddo

      if(msno == 1)then
         if(dzsno(1) > 0.03)then
         msno = 2
! Specified a new snow layer
         dzsno(1) = dzsno(1)/2.
         swice(1) = swice(1)/2.
         swliq(1) = swliq(1)/2.

         dzsno(2) = dzsno(1)
         swice(2) = swice(1)
         swliq(2) = swliq(1)
         tsno(2)  = tsno(1)
!        write(6,*)'Subdivided Top Node into two layer (1/2)'
         endif
      endif

      if(msno > 1)then
         if(dzsno(1) > 0.02)then
         drr = dzsno(1) - 0.02
         propor = drr/dzsno(1)
         zwice = propor*swice(1)
         zwliq = propor*swliq(1)

         propor = 0.02/dzsno(1)
         swice(1) = propor*swice(1)
         swliq(1) = propor*swliq(1)
         dzsno(1) = 0.02

         call combo(dzsno(2),swliq(2),swice(2),tsno(2), &
                    drr,zwliq,zwice,tsno(1))

!        write(6,*) 'Subdivided Top Node &
!                    20 mm combined into underlying neighbor'

         if(msno <= 2 .AND. dzsno(2) > 0.07)then
! subdivided a new layer
            msno = 3
            dzsno(2) = dzsno(2)/2.
            swice(2) = swice(2)/2.
            swliq(2) = swliq(2)/2.

            dzsno(3) = dzsno(2)
            swice(3) = swice(2)
            swliq(3) = swliq(2)
            tsno(3)  = tsno(2)
         endif
         endif
      endif

      if(msno > 2)then
         if(dzsno(2) > 0.05)then
         drr = dzsno(2) - 0.05
         propor = drr/dzsno(2)
         zwice = propor*swice(2)
         zwliq = propor*swliq(2)

         propor = 0.05/dzsno(2)
         swice(2) = propor*swice(2)
         swliq(2) = propor*swliq(2)
         dzsno(2) = 0.05

         call combo(dzsno(3),swliq(3),swice(3),tsno(3), &
                    drr,     zwliq,   zwice,   tsno(2))

!        write(6,*)'Subdivided 50 mm from the subsface layer &
!                   &and combined into underlying neighbor'

         if(msno <= 3 .AND. dzsno(3) > 0.18)then
! subdivided a new layer
            msno =  4
            dzsno(3) = dzsno(3)/2.
            swice(3) = swice(3)/2.
            swliq(3) = swliq(3)/2.

            dzsno(4) = dzsno(3)
            swice(4) = swice(3)
            swliq(4) = swliq(3)
            tsno(4)  = tsno(3)
         endif
         endif
      endif

      if(msno > 3)then
         if(dzsno(3) > 0.11)then
         drr = dzsno(3) - 0.11
         propor = drr/dzsno(3)
         zwice = propor*swice(3)
         zwliq = propor*swliq(3)

         propor = 0.11/dzsno(3)
         swice(3) = propor*swice(3)
         swliq(3) = propor*swliq(3)
         dzsno(3) = 0.11

         call combo(dzsno(4),swliq(4),swice(4),tsno(4), &
                    drr,     zwliq,   zwice,   tsno(3))

!        write(6,*)'Subdivided 110 mm from the third Node &
!                   &and combined into underlying neighbor'

         if(msno <= 4 .AND. dzsno(4) > 0.41)then
! subdivided a new layer
            msno = 5
            dzsno(4) = dzsno(4)/2.
            swice(4) = swice(4)/2.
            swliq(4) = swliq(4)/2.

            dzsno(5) = dzsno(4)
            swice(5) = swice(4)
            swliq(5) = swliq(4)
            tsno(5)  = tsno(4)
         endif
         endif
      endif
         
      if(msno > 4)then
         if(dzsno(4) > 0.23)then
         drr = dzsno(4) - 0.23
         propor = drr/dzsno(4)
         zwice = propor*swice(4)
         zwliq = propor*swliq(4)

         propor = 0.23/dzsno(4)
         swice(4) = propor*swice(4)
         swliq(4) = propor*swliq(4)
         dzsno(4) = 0.23

         call combo(dzsno(5),swliq(5),swice(5),tsno(5), &
                    drr,     zwliq,   zwice,   tsno(4))

!        write(6,*)'Subdivided 230 mm from the fourth Node &
!                   'and combined into underlying neighbor'
         endif
      endif
         
      snl = - msno

      do k = snl+1, 0
         dz_soisno(k)   = dzsno(k - snl)
         wice_soisno(k) = swice(k - snl)
         wliq_soisno(k) = swliq(k - snl)
         t_soisno(k)  = tsno (k - snl)
      enddo

      zi_soisno(0) = 0.
      do k = 0, snl+1, -1
         z_soisno(k)    = zi_soisno(k) - 0.5*dz_soisno(k)
         zi_soisno(k-1) = zi_soisno(k) - dz_soisno(k)
      enddo

 end subroutine snowlayersdivide



 subroutine combo ( dz_soisno, wliq_soisno, wice_soisno, t, &
                   dz2, wliq2, wice2, t2 )

!=======================================================================
! Original author: Yongjiu Dai, September 15, 1999
!
! combines two elements and returns the following combined
! variabless: dz_soisno, t, wliq_soisno, wice_soisno.
! the combined temperature is based on the equation:
! the sum of the enthalpies of the two elements = that of the combined element.
!
!=======================================================================

  use precision
  use PhysicalConstants, only : cpice, cpliq, hfus, tfrz
  implicit none

!-------------------------- Dummy argument -----------------------------

  real(r8), INTENT(in) :: dz2     ! nodal thickness of 2 elements being combined [m]
  real(r8), INTENT(in) :: wliq2   ! liquid water of element 2 [kg/m2]
  real(r8), INTENT(in) :: wice2   ! ice of element 2 [kg/m2]
  real(r8), INTENT(in) :: t2      ! nodal temperature of element 2 [K]

  real(r8), INTENT(inout) :: dz_soisno   ! nodal thickness of 1 elements being combined [m]
  real(r8), INTENT(inout) :: wliq_soisno ! liquid water of element 1
  real(r8), INTENT(inout) :: wice_soisno ! ice of element 1 [kg/m2]
  real(r8), INTENT(inout) :: t    ! nodel temperature of elment 1 [K]

!----------------------- Local variables ------------------------------

  real(r8) dzc    ! Total thickness of nodes 1 and 2 (dzc=dz_soisno+dz2).
  real(r8) wliqc  ! Combined liquid water [kg/m2]
  real(r8) wicec  ! Combined ice [kg/m2]
  real(r8) tc     ! Combined node temperature [K]
  real(r8) h      ! enthalpy of element 1 [J/m2]
  real(r8) h2     ! enthalpy of element 2 [J/m2]
  real(r8) hc     ! temporary

!-----------------------------------------------------------------------

      dzc = dz_soisno+dz2
      wicec = (wice_soisno+wice2)
      wliqc = (wliq_soisno+wliq2)
      h   = (cpice*wice_soisno+cpliq*wliq_soisno)*(t-tfrz)+hfus*wliq_soisno
      h2  = (cpice*wice2+cpliq*wliq2)*(t2-tfrz)+hfus*wliq2

      hc = h + h2
      if(hc < 0.)then
         tc = tfrz + hc/(cpice*wicec+cpliq*wliqc)
      else if(hc.le.hfus*wliqc)then
         tc = tfrz
      else
         tc = tfrz + (hc - hfus*wliqc)/(cpice*wicec+cpliq*wliqc)
      endif

      dz_soisno = dzc
      wice_soisno = wicec 
      wliq_soisno = wliqc
      t = tc

 end subroutine combo


END MODULE SNOW_Layers_CombineDivide
