#include <define.h>

MODULE MOD_PhaseChange

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: meltf
   PUBLIC :: meltf_snicar
   PUBLIC :: meltf_urban

!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------


   subroutine meltf (patchtype,lb,nl_soil,deltim, &
                     fact,brr,hs,dhsdT, &
                     t_soisno_bef,t_soisno,wliq_soisno,wice_soisno,imelt, &
                     scv,snowdp,sm,xmf,porsl,psi0,&
#ifdef Campbell_SOIL_MODEL
                     bsw,&
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
                     theta_r,alpha_vgm,n_vgm,L_vgm,&
                     sc_vgm,fc_vgm,&
#endif
                     dz)

!-----------------------------------------------------------------------
! DESCRIPTION:
! calculation of the phase change within snow and soil layers:
! (1) check the conditions which the phase change may take place,
!     i.e., the layer temperature is great than the freezing point
!     and the ice mass is not equal to zero (i.e., melting),
!     or layer temperature is less than the freezing point
!     and the liquid water mass is not equal to zero (i.e., freezing);
! (2) assess the rate of phase change from the energy excess (or deficit)
!     after setting the layer temperature to freezing point;
! (3) re-adjust the ice and liquid mass, and the layer temperature
!
! Original author : Yongjiu Dai, /09/1999/, /03/2014/
!
! Revisions:
! Nan Wei, 04/2023: supercooled soil water is included IF supercool is defined.
!-----------------------------------------------------------------------

   use MOD_Precision
   USE MOD_Hydro_SoilFunction
   use MOD_Const_Physical, only : tfrz, hfus,grav
   USE MOD_Namelist
   IMPLICIT NONE

!-----------------------------------------------------------------------

    integer, INTENT(in) :: patchtype    !land patch type (0=soil,1=urban or built-up,2=wetland,
                                        !3=land ice, 4=deep lake, 5=shallow lake)
    integer, INTENT(in) :: nl_soil             ! upper bound of array (i.e., soil layers)
    integer, INTENT(in) :: lb                  ! lower bound of array (i.e., snl +1)
   real(r8), INTENT(in) :: deltim              ! time step [second]
   real(r8), INTENT(in) :: t_soisno_bef(lb:nl_soil)  ! temperature at previous time step [K]
   real(r8), INTENT(in) :: brr (lb:nl_soil)    !
   real(r8), INTENT(in) :: fact(lb:nl_soil)    ! temporary variables
   real(r8), INTENT(in) :: hs                  ! net ground heat flux into the surface
   real(r8), INTENT(in) :: dhsdT               ! temperature derivative of "hs"
   real(r8), INTENT(in) :: porsl(1:nl_soil)    ! soil porosity [-]
   real(r8), INTENT(in) :: psi0 (1:nl_soil)    ! soil water suction, negative potential [mm]
#ifdef Campbell_SOIL_MODEL
   real(r8), INTENT(in) :: bsw(1:nl_soil)      ! clapp and hornbereger "b" parameter [-]
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   real(r8), INTENT(in) :: theta_r  (1:nl_soil), &
                           alpha_vgm(1:nl_soil), &
                           n_vgm    (1:nl_soil), &
                           L_vgm    (1:nl_soil), &
                           sc_vgm   (1:nl_soil), &
                           fc_vgm   (1:nl_soil)
#endif
   real(r8), INTENT(in) :: dz(1:nl_soil)      ! soil layer thickiness [m]

   real(r8), INTENT(inout) :: t_soisno (lb:nl_soil) ! temperature at current time step [K]
   real(r8), INTENT(inout) :: wice_soisno(lb:nl_soil) ! ice lens [kg/m2]
   real(r8), INTENT(inout) :: wliq_soisno(lb:nl_soil) ! liquid water [kg/m2]
   real(r8), INTENT(inout) :: scv              ! snow mass [kg/m2]
   real(r8), INTENT(inout) :: snowdp           ! snow depth [m]

   real(r8), INTENT(out) :: sm                 ! rate of snowmelt [mm/s, kg/(m2 s)]
   real(r8), INTENT(out) :: xmf                ! total latent heat of phase change
    integer, INTENT(out) :: imelt(lb:nl_soil)  ! flag for melting or freezing [-]

! Local
   real(r8) :: hm(lb:nl_soil)                  ! energy residual [W/m2]
   real(r8) :: xm(lb:nl_soil)                  ! metling or freezing within a time step [kg/m2]
   real(r8) :: heatr                           ! energy residual or loss after melting or freezing
   real(r8) :: temp1                           ! temporary variables [kg/m2]
   real(r8) :: temp2                           ! temporary variables [kg/m2]
   REAL(r8) :: smp
   REAL(r8) :: supercool(1:nl_soil)            ! the maximum liquid water when the soil temperature is below the freezing point [mm3/mm3]
   real(r8), dimension(lb:nl_soil) :: wmass0, wice0, wliq0
   real(r8) :: propor, tinc, we, scvold
   integer j

!-----------------------------------------------------------------------
   sm = 0.
   xmf = 0.
   do j = lb, nl_soil
      imelt(j) = 0
      hm(j) = 0.
      xm(j) = 0.
      wice0(j) = wice_soisno(j)
      wliq0(j) = wliq_soisno(j)
      wmass0(j) = wice_soisno(j) + wliq_soisno(j)
   enddo

   scvold=scv
   we=0.
   if(lb<=0) we = sum(wice_soisno(lb:0)+wliq_soisno(lb:0))

! supercooling water
   IF (DEF_USE_SUPERCOOL_WATER) THEN
      DO j = 1, nl_soil
         supercool(j) = 0.0
         if(t_soisno(j) < tfrz .and. patchtype <=2 ) then
            smp = hfus * (t_soisno(j)-tfrz)/(grav*t_soisno(j)) * 1000.     ! mm
            if (porsl(j) > 0.) then
#ifdef Campbell_SOIL_MODEL
               supercool(j) = porsl(j)*(smp/psi0(j))**(-1.0/bsw(j))
#else
               supercool(j) = soil_vliq_from_psi(smp, porsl(j), theta_r(j), -10.0, 5, &
                  (/alpha_vgm(j), n_vgm(j), L_vgm(j), sc_vgm(j), fc_vgm(j)/))
#endif
            else
               supercool(j) = 0.
            end if
            supercool(j) = supercool(j)*dz(j)*1000.              ! mm
         end if
      END do
   ENDIF

   do j = lb, nl_soil
      ! Melting identification
      ! if ice exists above melt point, melt some to liquid.
      if(wice_soisno(j) > 0. .and. t_soisno(j) > tfrz)then
         imelt(j) = 1
         t_soisno(j) = tfrz
      endif

      ! Freezing identification
      ! if liquid exists below melt point, freeze some to ice.
      IF(j <= 0)then
         if(wliq_soisno(j) > 0. .and. t_soisno(j) < tfrz) then
            imelt(j) = 2
            t_soisno(j) = tfrz
         endif
      ELSE
         if (DEF_USE_SUPERCOOL_WATER) then
            if(wliq_soisno(j) > supercool(j) .and. t_soisno(j) < tfrz) then
               imelt(j) = 2
               t_soisno(j) = tfrz
            endif
         else
            if(wliq_soisno(j) > 0. .and. t_soisno(j) < tfrz) then
               imelt(j) = 2
               t_soisno(j) = tfrz
            endif
         endif
      END if
   enddo

! If snow exists, but its thickness less than the critical value (0.01 m)
   if(lb == 1 .and. scv > 0.)then
      if(t_soisno(1) > tfrz)then
         imelt(1) = 1
         t_soisno(1) = tfrz
      endif
   endif

! Calculate the energy surplus and loss for melting and freezing
   do j = lb, nl_soil
      if(imelt(j) > 0)then
         tinc = t_soisno(j)-t_soisno_bef(j)
         if(j > lb)then
            hm(j) = brr(j) - tinc/fact(j)
         else
            hm(j) = hs + dhsdT*tinc + brr(j) - tinc/fact(j)
         endif
      endif
   enddo

   do j = lb, nl_soil
      if(imelt(j) == 1 .and. hm(j) < 0.) then
        hm(j) = 0.
        imelt(j) = 0
      endif
! this error was checked carefully, it results from the the computed error
! of "Tridiagonal-Matrix" in subroutine "thermal".
      if(imelt(j) == 2 .and. hm(j) > 0.) then
        hm(j) = 0.
        imelt(j) = 0
      endif
   enddo

! The rate of melting and freezing
   do j = lb, nl_soil
      if(imelt(j) > 0 .and. abs(hm(j)) > .0) then
         xm(j) = hm(j)*deltim/hfus                    ! kg/m2

         ! if snow exists, but its thickness less than the critical value (1 cm)
         ! Note: more work is need on how to tune the snow depth at this case
         if(j == 1 .and. lb == 1 .and. scv > 0. .and. xm(j) > 0.)then
            temp1 = scv                               ! kg/m2
            scv = max(0.,temp1-xm(j))
            propor = scv/temp1
            snowdp = propor * snowdp
            heatr = hm(j) - hfus*(temp1-scv)/deltim   ! W/m2
            if(heatr > 0.) then
               xm(j) = heatr*deltim/hfus              ! kg/m2
               hm(j) = heatr                          ! W/m2
            else
               xm(j) = 0.
               hm(j) = 0.
            endif
            sm = max(0.,(temp1-scv))/deltim           ! kg/(m2 s)
            xmf = hfus*sm
         endif

         heatr = 0.
         if(xm(j) > 0.) then
            wice_soisno(j) = max(0., wice0(j)-xm(j))
            heatr = hm(j) - hfus*(wice0(j)-wice_soisno(j))/deltim
         else
            if(j <= 0) then  ! snow
               wice_soisno(j) = min(wmass0(j), wice0(j)-xm(j))
            else
               if (DEF_USE_SUPERCOOL_WATER) then
                  if(wmass0(j) < supercool(j)) then
                     wice_soisno(j) = 0.
                  else
                     wice_soisno(j) = min(wmass0(j)-supercool(j), wice0(j)-xm(j))
                  endif
               else
                  wice_soisno(j) = min(wmass0(j), wice0(j)-xm(j))
               endif
            endif
            heatr = hm(j) - hfus*(wice0(j)-wice_soisno(j))/deltim
         endif

         wliq_soisno(j) = max(0.,wmass0(j)-wice_soisno(j))

         if(abs(heatr) > 0.)then
            if(j > lb)then
               t_soisno(j) = t_soisno(j) + fact(j)*heatr
            else
               t_soisno(j) = t_soisno(j) + fact(j)*heatr/(1.-fact(j)*dhsdT)
            endif
            if (DEF_USE_SUPERCOOL_WATER) then
               IF(j <= 0 .or. patchtype == 3)THEN !snow
                  if(wliq_soisno(j)*wice_soisno(j) > 0.) t_soisno(j) = tfrz
               ENDIF
            ELSE
               if(wliq_soisno(j)*wice_soisno(j) > 0.) t_soisno(j) = tfrz
            ENDIF
         endif

         xmf = xmf + hfus * (wice0(j)-wice_soisno(j))/deltim

         if(imelt(j) == 1 .and. j < 1) &
         sm = sm + max(0.,(wice0(j)-wice_soisno(j)))/deltim

      endif
   enddo

   !scvold=scv
   if(lb<=0) then
   we = sum(wice_soisno(lb:0)+wliq_soisno(lb:0))-we
      if(abs(we)>1.e-6) then
         print*, 'meltf err : ', we
      endif
   endif

   end subroutine meltf


   subroutine meltf_snicar (patchtype,lb,nl_soil,deltim, &
                     fact,brr,hs,dhsdT,sabg_lyr, &
                     t_soisno_bef,t_soisno,wliq_soisno,wice_soisno,imelt, &
                     scv,snowdp,sm,xmf,porsl,psi0,&
#ifdef Campbell_SOIL_MODEL
                     bsw,&
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
                     theta_r,alpha_vgm,n_vgm,L_vgm,&
                     sc_vgm,fc_vgm,&
#endif
                     dz)

!-----------------------------------------------------------------------
! DESCRIPTION:
! calculation of the phase change within snow and soil layers:
! (1) check the conditions which the phase change may take place,
!     i.e., the layer temperature is great than the freezing point
!     and the ice mass is not equal to zero (i.e., melting),
!     or layer temperature is less than the freezing point
!     and the liquid water mass is not equal to zero (i.e., freezing);
! (2) assess the rate of phase change from the energy excess (or deficit)
!     after setting the layer temperature to freezing point;
! (3) re-adjust the ice and liquid mass, and the layer temperature
!
! Original author : Yongjiu Dai, /09/1999/, /03/2014/
!
! Revisions:
! Hua Yuan, 01/2023: added snow layer absorption in melting calculation
! Nan Wei , 04/2023: supercooled soil water is included IF supercool is defined.
!-----------------------------------------------------------------------

   use MOD_Precision
   USE MOD_Hydro_SoilFunction
   use MOD_Const_Physical, only : tfrz, hfus, grav
   USE MOD_Namelist
   IMPLICIT NONE

!-----------------------------------------------------------------------

    integer, INTENT(in) :: patchtype    !land patch type (0=soil,1=urban or built-up,2=wetland,
                                        !3=land ice, 4=deep lake, 5=shallow lake)
    integer, INTENT(in) :: nl_soil             ! upper bound of array (i.e., soil layers)
    integer, INTENT(in) :: lb                  ! lower bound of array (i.e., snl +1)
   real(r8), INTENT(in) :: deltim              ! time step [second]
   real(r8), INTENT(in) :: t_soisno_bef(lb:nl_soil)  ! temperature at previous time step [K]
   real(r8), INTENT(in) :: brr (lb:nl_soil)    !
   real(r8), INTENT(in) :: fact(lb:nl_soil)    ! temporary variables
   real(r8), INTENT(in) :: hs                  ! net ground heat flux into the surface
   real(r8), INTENT(in) :: dhsdT               ! temperature derivative of "hs"
   real(r8), INTENT(in) :: sabg_lyr (lb:1)     ! snow layer absorption [W/m-2]
   real(r8), INTENT(in) :: porsl(1:nl_soil)    ! soil porosity [-]
   real(r8), INTENT(in) :: psi0 (1:nl_soil)    ! soil water suction, negative potential [mm]
#ifdef Campbell_SOIL_MODEL
   real(r8), INTENT(in) :: bsw(1:nl_soil)      ! clapp and hornbereger "b" parameter [-]
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   real(r8), INTENT(in) :: theta_r  (1:nl_soil), &
                           alpha_vgm(1:nl_soil), &
                           n_vgm    (1:nl_soil), &
                           L_vgm    (1:nl_soil), &
                           sc_vgm   (1:nl_soil), &
                           fc_vgm   (1:nl_soil)
#endif
   real(r8), INTENT(in) :: dz(1:nl_soil)       ! soil layer thickiness [m]

   real(r8), INTENT(inout) :: t_soisno (lb:nl_soil) ! temperature at current time step [K]
   real(r8), INTENT(inout) :: wice_soisno(lb:nl_soil) ! ice lens [kg/m2]
   real(r8), INTENT(inout) :: wliq_soisno(lb:nl_soil) ! liquid water [kg/m2]
   real(r8), INTENT(inout) :: scv              ! snow mass [kg/m2]
   real(r8), INTENT(inout) :: snowdp           ! snow depth [m]

   real(r8), INTENT(out) :: sm                 ! rate of snowmelt [mm/s, kg/(m2 s)]
   real(r8), INTENT(out) :: xmf                ! total latent heat of phase change
    integer, INTENT(out) :: imelt(lb:nl_soil)  ! flag for melting or freezing [-]

! Local
   real(r8) :: hm(lb:nl_soil)                  ! energy residual [W/m2]
   real(r8) :: xm(lb:nl_soil)                  ! metling or freezing within a time step [kg/m2]
   real(r8) :: heatr                           ! energy residual or loss after melting or freezing
   real(r8) :: temp1                           ! temporary variables [kg/m2]
   real(r8) :: temp2                           ! temporary variables [kg/m2]
   REAL(r8) :: smp
   REAL(r8) :: supercool(1:nl_soil)            ! the maximum liquid water when the soil temperature is below the   freezing point [mm3/mm3]
   real(r8), dimension(lb:nl_soil) :: wmass0, wice0, wliq0
   real(r8) :: propor, tinc, we, scvold
   integer j

!-----------------------------------------------------------------------

   sm = 0.
   xmf = 0.
   do j = lb, nl_soil
      imelt(j) = 0
      hm(j) = 0.
      xm(j) = 0.
      wice0(j) = wice_soisno(j)
      wliq0(j) = wliq_soisno(j)
      wmass0(j) = wice_soisno(j) + wliq_soisno(j)
   enddo

   scvold=scv
   we=0.
   if(lb<=0) we = sum(wice_soisno(lb:0)+wliq_soisno(lb:0))

! supercooling water
   if (DEF_USE_SUPERCOOL_WATER) then
      DO j = 1, nl_soil
         supercool(j) = 0.0
         if(t_soisno(j) < tfrz .and. patchtype <= 2) then
            smp = hfus * (t_soisno(j)-tfrz)/(grav*t_soisno(j)) * 1000.     ! mm
            if (porsl(j) > 0.) then
#ifdef Campbell_SOIL_MODEL
               supercool(j) = porsl(j)*(smp/psi0(j))**(-1.0/bsw(j))
#else
               supercool(j) = soil_vliq_from_psi(smp, porsl(j), theta_r(j), -10.0, 5, &
                  (/alpha_vgm(j), n_vgm(j), L_vgm(j), sc_vgm(j), fc_vgm(j)/))
#endif
            else
               supercool(j) = 0.
            end if
            supercool(j) = supercool(j)*dz(j)*1000.              ! mm
         end if
      END do
   endif


   do j = lb, nl_soil
      ! Melting identification
      ! if ice exists above melt point, melt some to liquid.
      if(wice_soisno(j) > 0. .and. t_soisno(j) > tfrz)then
         imelt(j) = 1
         t_soisno(j) = tfrz
      endif

      ! Freezing identification
      ! if liquid exists below melt point, freeze some to ice.
      IF(j <= 0)then
         if(wliq_soisno(j) > 0. .and. t_soisno(j) < tfrz) then
            imelt(j) = 2
            t_soisno(j) = tfrz
         endif
      ELSE
         if (DEF_USE_SUPERCOOL_WATER) then
            if(wliq_soisno(j) > supercool(j) .and. t_soisno(j) < tfrz) then
               imelt(j) = 2
               t_soisno(j) = tfrz
            endif
         else
            if(wliq_soisno(j) > 0. .and. t_soisno(j) < tfrz) then
               imelt(j) = 2
               t_soisno(j) = tfrz
            endif
         endif
      END if
   enddo

! If snow exists, but its thickness less than the critical value (0.01 m)
   if(lb == 1 .and. scv > 0.)then
      if(t_soisno(1) > tfrz)then
         imelt(1) = 1
         t_soisno(1) = tfrz
      endif
   endif

! Calculate the energy surplus and loss for melting and freezing
   do j = lb, nl_soil
      if(imelt(j) > 0)then
         tinc = t_soisno(j)-t_soisno_bef(j)
         if(j > lb)then
            IF (j <= 1) THEN
               hm(j) = brr(j) - tinc/fact(j) + sabg_lyr(j)
            ELSE
               hm(j) = brr(j) - tinc/fact(j)
            ENDIF
         else
            hm(j) = hs + dhsdT*tinc + brr(j) - tinc/fact(j)
         endif
      endif
   enddo

   do j = lb, nl_soil
      if(imelt(j) == 1 .and. hm(j) < 0.) then
        hm(j) = 0.
        imelt(j) = 0
      endif
! this error was checked carefully, it results from the the computed error
! of "Tridiagonal-Matrix" in subroutine "thermal".
      if(imelt(j) == 2 .and. hm(j) > 0.) then
        hm(j) = 0.
        imelt(j) = 0
      endif
   enddo

! The rate of melting and freezing
   do j = lb, nl_soil
      if(imelt(j) > 0 .and. abs(hm(j)) > .0) then
         xm(j) = hm(j)*deltim/hfus                    ! kg/m2

         ! if snow exists, but its thickness less than the critical value (1 cm)
         ! Note: more work is need on how to tune the snow depth at this case
         if(j == 1 .and. lb == 1 .and. scv > 0. .and. xm(j) > 0.)then
            temp1 = scv                               ! kg/m2
            scv = max(0.,temp1-xm(j))
            propor = scv/temp1
            snowdp = propor * snowdp
            heatr = hm(j) - hfus*(temp1-scv)/deltim   ! W/m2
            if(heatr > 0.) then
               xm(j) = heatr*deltim/hfus              ! kg/m2
               hm(j) = heatr                          ! W/m2
            else
               xm(j) = 0.
               hm(j) = 0.
            endif
            sm = max(0.,(temp1-scv))/deltim           ! kg/(m2 s)
            xmf = hfus*sm
         endif

         heatr = 0.
         if(xm(j) > 0.) then
            wice_soisno(j) = max(0., wice0(j)-xm(j))
            heatr = hm(j) - hfus*(wice0(j)-wice_soisno(j))/deltim
         else
            IF(j <= 0) THEN ! snow
               wice_soisno(j) = min(wmass0(j), wice0(j)-xm(j))
            ELSE
               if (DEF_USE_SUPERCOOL_WATER) then
                  if(wmass0(j) < supercool(j)) then
                     wice_soisno(j) = 0.
                  else
                     wice_soisno(j) = min(wmass0(j)-supercool(j), wice0(j)-xm(j))
                  endif
               else
                  wice_soisno(j) = min(wmass0(j), wice0(j)-xm(j))
               endif
            endif
            heatr = hm(j) - hfus*(wice0(j)-wice_soisno(j))/deltim
         endif

         wliq_soisno(j) = max(0.,wmass0(j)-wice_soisno(j))

         if(abs(heatr) > 0.)then
            if(j > lb)then
               t_soisno(j) = t_soisno(j) + fact(j)*heatr
            else
               t_soisno(j) = t_soisno(j) + fact(j)*heatr/(1.-fact(j)*dhsdT)
            endif
            if (DEF_USE_SUPERCOOL_WATER) then
               IF(j <= 0 .or. patchtype == 3)THEN !snow
                  if(wliq_soisno(j)*wice_soisno(j) > 0.) t_soisno(j) = tfrz
               ENDIF
            ELSE
               if(wliq_soisno(j)*wice_soisno(j) > 0.) t_soisno(j) = tfrz
            ENDIF

         endif

         xmf = xmf + hfus * (wice0(j)-wice_soisno(j))/deltim

         if(imelt(j) == 1 .and. j < 1) &
         sm = sm + max(0.,(wice0(j)-wice_soisno(j)))/deltim

      endif
   enddo

   !scvold=scv
   if(lb<=0) then
   we = sum(wice_soisno(lb:0)+wliq_soisno(lb:0))-we
      if(abs(we)>1.e-6) then
         print*, 'meltf err : ', we
      endif
   endif

   end subroutine meltf_snicar

   subroutine meltf_urban (lb,nl_soil,deltim, &
                   fact,brr,hs,dhsdT, &
                   t_soisno_bef,t_soisno,wliq_soisno,wice_soisno,imelt, &
                   scv,snowdp,sm,xmf)

!-----------------------------------------------------------------------
! Original author : Yongjiu Dai, /09/1999/, /03/2014/
!
! calculation of the phase change within snow and soil layers:
!
! (1) check the conditions which the phase change may take place,
!     i.e., the layer temperature is great than the freezing point
!     and the ice mass is not equal to zero (i.e., melting),
!     or layer temperature is less than the freezing point
!     and the liquid water mass is not equal to zero (i.e., freezing);
! (2) assess the rate of phase change from the energy excess (or deficit)
!     after setting the layer temperature to freezing point;
! (3) re-adjust the ice and liquid mass, and the layer temperature
!
!-----------------------------------------------------------------------

     use MOD_Precision
     use MOD_Const_Physical, only : tfrz, hfus
     IMPLICIT NONE

!-----------------------------------------------------------------------

     integer, INTENT(in) :: nl_soil             ! upper bound of array (i.e., soil layers)
     integer, INTENT(in) :: lb                  ! lower bound of array (i.e., snl +1)
     real(r8), INTENT(in) :: deltim              ! time step [second]
     real(r8), INTENT(in) :: t_soisno_bef(lb:nl_soil)  ! temperature at previous time step [K]
     real(r8), INTENT(in) :: brr (lb:nl_soil)    !
     real(r8), INTENT(in) :: fact(lb:nl_soil)    ! temporary variables
     real(r8), INTENT(in) :: hs                  ! net ground heat flux into the surface
     real(r8), INTENT(in) :: dhsdT               ! temperature derivative of "hs"

     real(r8), INTENT(inout) :: t_soisno (lb:nl_soil) ! temperature at current time step [K]
     real(r8), INTENT(inout) :: wice_soisno(lb:nl_soil) ! ice lens [kg/m2]
     real(r8), INTENT(inout) :: wliq_soisno(lb:nl_soil) ! liquid water [kg/m2]
     real(r8), INTENT(inout) :: scv              ! snow mass [kg/m2]
     real(r8), INTENT(inout) :: snowdp           ! snow depth [m]

     real(r8), INTENT(out) :: sm                 ! rate of snowmelt [mm/s, kg/(m2 s)]
     real(r8), INTENT(out) :: xmf                ! total latent heat of phase change
     integer, INTENT(out) :: imelt(lb:nl_soil)  ! flag for melting or freezing [-]

! Local
     real(r8) :: hm(lb:nl_soil)                  ! energy residual [W/m2]
     real(r8) :: xm(lb:nl_soil)                  ! metling or freezing within a time step [kg/m2]
     real(r8) :: heatr                           ! energy residual or loss after melting or freezing
     real(r8) :: temp1                           ! temporary variables [kg/m2]
     real(r8) :: temp2                           ! temporary variables [kg/m2]

     real(r8), dimension(lb:nl_soil) :: wmass0, wice0, wliq0
     real(r8) :: propor, tinc, we, scvold
     integer j

!-----------------------------------------------------------------------

     sm = 0.
     xmf = 0.
     do j = lb, nl_soil
        imelt(j) = 0
        hm(j) = 0.
        xm(j) = 0.
        wice0(j) = wice_soisno(j)
        wliq0(j) = wliq_soisno(j)
        wmass0(j) = wice_soisno(j) + wliq_soisno(j)
     enddo

     scvold=scv
     we=0.
     if(lb<=0) we = sum(wice_soisno(lb:0)+wliq_soisno(lb:0))

     do j = lb, nl_soil
        ! Melting identification
        ! if ice exists above melt point, melt some to liquid.
        if(wice_soisno(j) > 0. .and. t_soisno(j) > tfrz)then
           imelt(j) = 1
           t_soisno(j) = tfrz
        endif

     ! Freezing identification
     ! if liquid exists below melt point, freeze some to ice.
        if(wliq_soisno(j) > 0. .and. t_soisno(j) < tfrz) then
           imelt(j) = 2
           t_soisno(j) = tfrz
        endif
     enddo

! If snow exists, but its thickness less than the critical value (0.01 m)
     if(lb == 1 .and. scv > 0.)then
        if(t_soisno(1) > tfrz)then
           imelt(1) = 1
           t_soisno(1) = tfrz
        endif
     endif

! Calculate the energy surplus and loss for melting and freezing
     do j = lb, nl_soil
        if(imelt(j) > 0)then
           tinc = t_soisno(j)-t_soisno_bef(j)
           if(j > lb)then
              hm(j) = brr(j) - tinc/fact(j)
           else
              hm(j) = hs + dhsdT*tinc + brr(j) - tinc/fact(j)
           endif
        endif
     enddo

     do j = lb, nl_soil
        if(imelt(j) == 1 .and. hm(j) < 0.) then
           hm(j) = 0.
           imelt(j) = 0
        endif
! this error was checked carefully, it results from the the computed error
! of "Tridiagonal-Matrix" in subroutine "thermal".
        if(imelt(j) == 2 .and. hm(j) > 0.) then
           hm(j) = 0.
           imelt(j) = 0
        endif
     enddo

! The rate of melting and freezing
     do j = lb, nl_soil
        if(imelt(j) > 0 .and. abs(hm(j)) > .0) then
           xm(j) = hm(j)*deltim/hfus                    ! kg/m2

        ! if snow exists, but its thickness less than the critical value (1 cm)
        ! Note: more work is need on how to tune the snow depth at this case
           if(j == 1 .and. lb == 1 .and. scv > 0. .and. xm(j) > 0.)then
              temp1 = scv                               ! kg/m2
              scv = max(0.,temp1-xm(j))
              propor = scv/temp1
              snowdp = propor * snowdp
              heatr = hm(j) - hfus*(temp1-scv)/deltim   ! W/m2
              if(heatr > 0.) then
                 xm(j) = heatr*deltim/hfus              ! kg/m2
                 hm(j) = heatr                          ! W/m2
              else
                 xm(j) = 0.
                 hm(j) = 0.
              endif
              sm = max(0.,(temp1-scv))/deltim           ! kg/(m2 s)
              xmf = hfus*sm
           endif

           heatr = 0.
           if(xm(j) > 0.) then
              wice_soisno(j) = max(0., wice0(j)-xm(j))
              heatr = hm(j) - hfus*(wice0(j)-wice_soisno(j))/deltim
           else
              wice_soisno(j) = min(wmass0(j), wice0(j)-xm(j))
              heatr = hm(j) - hfus*(wice0(j)-wice_soisno(j))/deltim
           endif

           wliq_soisno(j) = max(0.,wmass0(j)-wice_soisno(j))

           if(abs(heatr) > 0.)then
              if(j > lb)then
                 t_soisno(j) = t_soisno(j) + fact(j)*heatr
              else
                 t_soisno(j) = t_soisno(j) + fact(j)*heatr/(1.-fact(j)*dhsdT)
              endif
              if(wliq_soisno(j)*wice_soisno(j) > 0.) t_soisno(j) = tfrz
           endif

           xmf = xmf + hfus * (wice0(j)-wice_soisno(j))/deltim

           if(imelt(j) == 1 .and. j < 1) &
           sm = sm + max(0.,(wice0(j)-wice_soisno(j)))/deltim

        endif
     enddo

  !scvold=scv
     if(lb<=0) then
        we = sum(wice_soisno(lb:0)+wliq_soisno(lb:0))-we
        if(abs(we)>1.e-6) then
           print*, 'meltf err : ', we
        endif
     endif

  end subroutine meltf_urban

END MODULE MOD_PhaseChange
