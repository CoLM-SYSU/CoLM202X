#include <define.h>

#ifdef LULCC
MODULE MOD_Lulcc_MassEnergyConserve

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: LulccMassEnergyConserve


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------


  SUBROUTINE LulccMassEnergyConserve
! -------------------------------
! Created by Wanyi Lin and Hua Yuan, 07/2023
!
! !REVISONS:
!
! -------------------------------

!=======================================================================
! This is the main subroutine to execute the calculation of the restart
! variables for the begin of next year.
! There are mainly three ways to adjust restart variables:
! 1) variable related to mass: area weighted mean of the source patches,
!    e.g., ldew, wliq_soisno.
! 2) variable related to energy: keep energy conserve after the change
!    of temperature, e.g., t_soisno.
! 3) recalculate according to physical process, e.g., dz_sno, scv, fsno.
!=======================================================================

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_LandPatch
   USE MOD_LandElm
   USE MOD_Mesh
   USE MOD_SPMD_Task
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables
   USE MOD_Lulcc_Vars_TimeInvariants
   USE MOD_Lulcc_Vars_TimeVariables
   USE MOD_Lulcc_PatchTrace
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
   USE MOD_Vars_PFTimeInvariants
   USE MOD_Vars_PFTimeVariables
#endif
#ifdef URBAN_MODEL
   USE MOD_LandUrban
   USE MOD_Urban_Vars_TimeVariables
   USE MOD_Urban_Vars_TimeInvariants
#endif
   USE MOD_Const_Physical, only: cpice, cpliq, denh2o, denice
   USE MOD_GroundTemperature
   USE MOD_Namelist

   IMPLICIT NONE

   integer, allocatable, dimension(:) :: grid_patch_s , grid_patch_e
   integer, allocatable, dimension(:) :: grid_patch_s_, grid_patch_e_
   integer, allocatable, dimension(:) :: locpxl
   integer :: numpxl,ipxl

   integer, allocatable :: frnp_(:)     !index of source patches
   integer, allocatable :: gu_(:)       !index of urban patches in last year's grid
   real(r8),allocatable :: cvsoil_(:,:) !heat capacity [J/(m2 K)]

   integer :: k, ilc, num, inp_
   integer :: i, j, np, np_, selfnp_, l, ipft, ip, ip_, pc, pc_
   integer :: nsl      ! number of snow layer of the source patch with maximum area
   integer :: nsl_max  ! maximum number of snow layer considering all source patches
   integer :: u, u_, iu, selfu_, nurb, duclass
   integer :: nlc = N_land_classification
   real(r8), dimension(1:N_land_classification) :: lccpct_np
   real(r8):: sum_lccpct_np, wgt(maxsnl+1:nl_soil)
   real(r8):: zi_sno(maxsnl+1:0) ! local variable for snow node and depth calculation
   real(r8):: vf_water ! volumetric fraction liquid water within soil
   real(r8):: vf_ice   ! volumetric fraction ice len within soil
   real(r8):: hcap     ! J/(m3 K)
   real(r8):: c_water  ! Specific heat of water * density of liquid water
   real(r8):: c_ice    ! Specific heat of ice   * density of ice
   real(r8):: denice_np(maxsnl+1:0), denh2o_np(maxsnl+1:0)
   real(r8):: wbef,wpre! water before and water present for water calculation heck
   real(r8):: fmelt    ! dimensionless metling factor
   real(r8), parameter :: m = 1.0 ! the value of m used in CLM4.5 is 1.0.
   logical :: FROM_SOIL

   IF (p_is_worker) THEN
      ! allocate with numelm
      allocate(grid_patch_s (numelm ))
      allocate(grid_patch_e (numelm ))
      allocate(grid_patch_s_(numelm_))
      allocate(grid_patch_e_(numelm_))

      grid_patch_e (:) = -1
      grid_patch_s (:) = -1
      grid_patch_e_(:) = -1
      grid_patch_s_(:) = -1

      ! loop for numelm of next year, patches at the beginning and end of
      ! the element were recorded landpatch%eindex is arranged in order,
      ! and the not land element is skipped so, if element is missing, the
      ! recorder is -1.
      DO i=1, numelm
         ! how many patches in ith element in this worker
         numpxl = count(landpatch%eindex==landelm%eindex(i))

         IF (allocated(locpxl)) deallocate(locpxl)
         allocate(locpxl(numpxl))

         ! get all patches' index that eindex is equal the i element
         locpxl = pack([(ipxl, ipxl=1, numpatch)], &
                      landpatch%eindex==landelm%eindex(i))
         ! the min index is the start of patch's index
         grid_patch_s(i) = minval(locpxl)
         ! the max index is the end of patch's index
         grid_patch_e(i) = maxval(locpxl)
      ENDDO

      ! same as above, loop for numelm of previous year
      ! patches at the beginning and end of the element were recorded
      DO i=1, numelm_
         numpxl = count(landpatch_%eindex==landelm_%eindex(i))

         IF (allocated(locpxl)) deallocate(locpxl)
         allocate(locpxl(numpxl))

         locpxl = pack([(ipxl, ipxl=1, numpatch_)], &
                   landpatch_%eindex==landelm_%eindex(i))

         grid_patch_s_(i) = minval(locpxl)
         grid_patch_e_(i) = maxval(locpxl)
      ENDDO

      DO i=1, numelm
         DO j=1,numelm_
            IF (landelm%eindex(i) == landelm_%eindex(j)) THEN
               np = grid_patch_s (i)
               np_= grid_patch_s_(j)

               IF (np.le.0) CYCLE

               DO WHILE (np.le.grid_patch_e(i))

IF (patchtype(np) .ne. 3) THEN !not a glacier patch

IF (DEF_USE_PFT .or. DEF_FAST_PC) THEN
                  lccpct_np(:) = 0
                  lccpct_np(1) = sum(lccpct_patches(np,:), mask=patchtypes(:)==0)
                  lccpct_np(URBAN  )   = lccpct_patches(np,URBAN  )
                  lccpct_np(WETLAND)   = lccpct_patches(np,WETLAND)
                  lccpct_np(WATERBODY) = lccpct_patches(np,WATERBODY)
ELSE
                  lccpct_np(:) = lccpct_patches(np,1:nlc)
ENDIF

                  num = count(lccpct_np .gt. 0)
                  sum_lccpct_np = sum(lccpct_np)
                  allocate ( frnp_  (                 num))
                  allocate ( cvsoil_(maxsnl+1:nl_soil,num))

                  ! Source patch type which differs from np's type exists
                  IF ( (sum_lccpct_np - lccpct_np(patchclass(np))) .gt. 0 ) THEN

                     ! Get the index of source patches, and stored as frnp_
                     k = 0
                     DO ilc = 1, nlc
                        IF (lccpct_np(ilc) .gt. 0) THEN
                           k = k + 1
                           inp_ = np_
                           DO WHILE (inp_ .le. grid_patch_e_(j))

                              ! Get the index of source patch that has the same LC, and stored as selfnp_
                              IF (patchclass_(inp_) .eq. patchclass(np)) THEN
                                 selfnp_ = inp_
                              ENDIF

                              IF (patchclass_(inp_) .eq. ilc) THEN
                                 frnp_(k) = inp_
                                 EXIT
                              ENDIF
                              inp_ = inp_ + 1
                           ENDDO

                        ELSE
                          CYCLE
                        ENDIF
                     ENDDO

                     ! Initialize
                     wliq_soisno (:,np)                = 0  !liquid water in layers [kg/m2]
                     wice_soisno (:,np)                = 0  !ice lens in layers [kg/m2]
                     t_soisno    (:,np)                = 0  !soil + snow layer temperature [K]
                     z_sno       (:,np)                = 0  !node depth [m]
                     dz_sno      (:,np)                = 0  !interface depth [m]
                     t_grnd        (np)                = 0  !ground surface temperature [K]
                     ldew          (np)                = 0  !depth of water on foliage [mm]
                     sag           (np)                = 0  !non dimensional snow age [-]
                     scv           (np)                = 0  !snow cover, water equivalent [mm]
                     snowdp        (np)                = 0  !snow depth [meter]
                     fsno          (np)                = 0  !fraction of snow cover on ground
                     sigf          (np)                = 0  !fraction of veg cover, excluding snow-covered veg [-]
                     zwt           (np)                = 0  !the depth to water table [m]
                     wa            (np)                = 0  !water storage in aquifer [mm]

                     cvsoil_(:,:) = 0
                     ! Weight of temperature adjustment
                     c_water = cpliq * denh2o ! J/(m3 K) = 4188   [J/(kg K)]*1000(kg/m3)
                     c_ice   = cpice * denice ! J/(m3 K) = 2117.27[J/(kg K)]*917 (kg/m3)
                     wgt(maxsnl+1:nl_soil) = 0

                     DO k = 1, num
                        ! Soil ground and wetland heat capacity
                        DO l = 1, nl_soil
                           vf_water = wliq_soisno_(l,frnp_(k))/(dz_soi(l)*denh2o)
                           vf_ice   = wice_soisno_(l,frnp_(k))/(dz_soi(l)*denice)
                           hcap     = csol_(l,frnp_(k)) + vf_water*c_water + vf_ice*c_ice
                           cvsoil_(l,k) = hcap*dz_soi(l)
                        ENDDO

                        IF( dz_sno_(0,frnp_(k))>0 .and. scv_(frnp_(k))>0.) cvsoil_(1,k) = cvsoil_(1,k) + cpice*scv_(frnp_(k))

                        ! Snow heat capacity
                        IF( z_sno_(0,frnp_(k)) < 0 ) THEN
                           cvsoil_(:0,k) = cpliq*wliq_soisno_(:0,frnp_(k)) + cpice*wice_soisno_(:0,frnp_(k))
                        ENDIF
                        wgt(maxsnl+1:nl_soil) = wgt(maxsnl+1:nl_soil) + cvsoil_(maxsnl+1:nl_soil,k) * lccpct_np(patchclass_(frnp_(k)))
                     ENDDO

                     ! Get the maximum lccpct for snow layers assignment
                     inp_ = frnp_(1)
                     k = 2
                     DO WHILE (k .le. num)
                        IF ( lccpct_np(patchclass_(frnp_(k))) .gt. lccpct_np(patchclass_(inp_)) ) THEN
                           inp_ = frnp_(k)
                        ENDIF
                        k = k + 1
                     ENDDO

                     ! check if snow layer exist in patch inp_
                     nsl     = count(z_sno_(:,inp_) .lt. 0)
                     nsl_max = count(wgt(:0)        .gt. 0)
                     denh2o_np(maxsnl+1:0) = 0
                     denice_np(maxsnl+1:0) = 0
                     IF (nsl > 0) THEN
                        ! move wgt above nsl to nsl
                        IF ( nsl_max > nsl) THEN
                           DO l = nsl+1, nsl_max
                              wgt(-nsl+1) = wgt(-nsl+1) + wgt(-l+1)
                           ENDDO
                        ENDIF

                        DO k = 1, num
                           t_soisno (-nsl+1:0,np) = t_soisno (-nsl+1:0,np) + t_soisno_(-nsl+1:0,frnp_(k))*cvsoil_(-nsl+1:0,k)*lccpct_np(patchclass_(frnp_(k)))/wgt(-nsl+1:0)
                           wliq_soisno (-nsl+1:0,np) = wliq_soisno (-nsl+1:0,np) + wliq_soisno_(-nsl+1:0,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                           wice_soisno (-nsl+1:0,np) = wice_soisno (-nsl+1:0,np) + wice_soisno_(-nsl+1:0,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np

                           l = 1
                           DO WHILE ( (l .le. nsl) .and. (dz_sno_(-l+1,frnp_(k)) .gt. 0) )
                              denh2o_np (-l+1) = denh2o_np(-l+1) + wliq_soisno_(-l+1,frnp_(k))/dz_sno_(-l+1,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                              denice_np (-l+1) = denice_np(-l+1) + wice_soisno_(-l+1,frnp_(k))/dz_sno_(-l+1,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                              l = l + 1
                              IF (l .gt. -maxsnl) EXIT
                           ENDDO

                           ! if source patch has more snow layer than the main patch
                           IF (nsl .lt. -maxsnl) THEN
                              l = nsl+1
                              DO WHILE ( (l .le. -maxsnl) .and. (dz_sno_(-l+1,frnp_(k)) .gt. 0) )
                                 wliq_soisno(-nsl+1,np) = wliq_soisno (-nsl+1,np) + wliq_soisno_(-l+1,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                                 wice_soisno(-nsl+1,np) = wice_soisno (-nsl+1,np) + wice_soisno_(-l+1,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                                 t_soisno (-nsl+1,np) = t_soisno (-nsl+1,np) + t_soisno_(-l+1,frnp_(k))*cvsoil_(-l+1,k)*lccpct_np(patchclass_(frnp_(k)))/wgt(-nsl+1)
                                 denh2o_np (-nsl+1) = denh2o_np(-nsl+1) + wliq_soisno_(-l+1,frnp_(k))/dz_sno_(-l+1,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                                 denice_np (-nsl+1) = denice_np(-nsl+1) + wice_soisno_(-l+1,frnp_(k))/dz_sno_(-l+1,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                                 l = l + 1
                                 IF (l .gt. -maxsnl) EXIT
                              ENDDO
                           ENDIF
                        ENDDO

                        ! snow layer node and depth calculation according to new mass and density
                        zi_sno(0) = 0._r8
                        DO l = 0, -nsl+1, -1
                           IF (denice_np(l)>0 .and. denh2o_np(l)>0) THEN
                              dz_sno (l,np) = wice_soisno(l,np)/denice_np(l) + wliq_soisno(l,np)/denh2o_np(l)
                           ELSEIF (denice_np(l)==0 .and. denh2o_np(l)>0) THEN
                              dz_sno (l,np) = wliq_soisno(l,np)/denh2o_np(l)
                              ! print*, 'denice=0! np=',np,'igbp=',patchclass(np),'nsl=',nsl,'frnp_=',frnp_
                              ! DO k = 1,num
                              !    print*,'frnp_=',frnp_(k),'wice=',wice_soisno(:0,frnp_(k))
                              ! ENDDO
                           ELSEIF (denh2o_np(l)==0 .and. denice_np(l)>0) THEN
                              dz_sno (l,np) = wice_soisno(l,np)/denice_np(l)
                              ! print*, 'denh2o=0! np=',np,'igbp=',patchclass(np),'nsl=',nsl,'frnp_=',frnp_
                              ! DO k = 1,num
                              !    print*,'frnp_=',frnp_(k),'wliq=',wliq_soisno(:0,frnp_(k))
                              ! ENDDO
                           ELSE
                              print*, 'denh2o and denice == 0! np=',np,'igbp=',patchclass(np),'nsl=',nsl,'frnp_=',frnp_
                              DO k = 1,num
                                 print*,'frnp_=',frnp_(k),'wliq=',wliq_soisno(:0,frnp_(k))
                                 print*,'frnp_=',frnp_(k),'wice=',wice_soisno(:0,frnp_(k))
                              ENDDO
                              CALL CoLM_stop()
                           ENDIF
                           z_sno  (l,np) = zi_sno(l) - 0.5_r8*dz_sno(l,np)
                           IF (l-1 .lt. maxsnl+1) EXIT
                           zi_sno (l-1)  = zi_sno(l) - dz_sno(l,np)
                        ENDDO

                     ELSE
                        ! no snow layer exist in the main patch, add a layer
                        ! move wgt above soil to layer 0
                        IF ( nsl_max > nsl) THEN
                           DO l = nsl+1, nsl_max
                              wgt(0) = wgt(0) + wgt(-l+1)
                           ENDDO
                        ENDIF
                        l=0
                        DO WHILE (wgt(l) .gt. 0)
                           DO k = 1, num
                              wliq_soisno(0,np) = wliq_soisno(0,np) + wliq_soisno_(l,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                              wice_soisno(0,np) = wice_soisno(0,np) + wice_soisno_(l,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                              t_soisno   (0,np) = t_soisno   (0,np) + t_soisno_(l,frnp_(k))*cvsoil_(l,k)*lccpct_np(patchclass_(frnp_(k)))/wgt(0)
                              IF (dz_sno_(l,frnp_(k)) .gt. 0) THEN
                                 denh2o_np(0) = denh2o_np(0) + wliq_soisno_(l,frnp_(k))/dz_sno_(l,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                                 denice_np(0) = denice_np(0) + wice_soisno_(l,frnp_(k))/dz_sno_(l,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                              ENDIF
                           ENDDO
                           l=l-1
                           IF (l .lt. maxsnl+1) EXIT
                        ENDDO

                        IF (wgt(0) .gt. 0) THEN
                           ! snow layer node and depth calculation according to new mass and density
                           IF (denh2o_np(0)>0 .and. denh2o_np(0)>0) THEN
                              dz_sno (0,np) = wice_soisno(0,np)/denice_np(0) + wliq_soisno(0,np)/denh2o_np(0)
                           ELSEIF (denice_np(0)==0 .and. denh2o_np(0)>0) THEN
                              dz_sno (0,np) = wliq_soisno(0,np)/denh2o_np(0)
                              ! print*, 'denice=0! stop! np=',np,'igbp=',patchclass(np),'nsl=',nsl,'frnp_=',frnp_
                              ! DO k = 1,num
                              !    print*,'frnp_=',frnp_(k),'wice=',wice_soisno(:0,frnp_(k))
                              ! ENDDO
                           ELSEIF (denice_np(0)>0 .and. denh2o_np(0)==0) THEN
                              dz_sno (0,np) = wice_soisno(0,np)/denice_np(0)
                              ! print*, 'denh2o=0! stop! np=',np,'igbp=',patchclass(np),'nsl=',nsl,'frnp_=',frnp_
                              ! DO k = 1,num
                              !    print*,'frnp_=',frnp_(k),'wliq=',wliq_soisno(:0,frnp_(k))
                              ! ENDDO
                           ELSE
                              print*, 'denh2o and denice == 0! np=',np,'igbp=',patchclass(np),'nsl=',nsl,'frnp_=',frnp_
                              DO k = 1,num
                                 print*,'frnp_=',frnp_(k),'wliq=',wliq_soisno(:0,frnp_(k))
                                 print*,'frnp_=',frnp_(k),'wice=',wice_soisno(:0,frnp_(k))
                              ENDDO
                              CALL CoLM_stop()
                           ENDIF

                           zi_sno (0)    = 0._r8
                           z_sno  (0,np) = zi_sno(0) - 0.5_r8*dz_sno(0,np)
                        ENDIF

                     ENDIF

                     DO l = maxsnl+1, 0
                        IF ( z_sno(l,np) .lt. 0 ) THEN
                           scv(np) = scv(np) + wice_soisno(l,np) + wliq_soisno(l,np)
                           snowdp(np) = snowdp(np) + dz_sno(l,np)
                        ENDIF
                     ENDDO

                     ! Variable adjustment
                     DO k = 1, num
                        wliq_soisno (1:nl_soil,np) = wliq_soisno (1:nl_soil,np) + wliq_soisno_(1:nl_soil,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                        wice_soisno (1:nl_soil,np) = wice_soisno (1:nl_soil,np) + wice_soisno_(1:nl_soil,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                        t_soisno (1:nl_soil,np) = t_soisno (1:nl_soil,np) + t_soisno_(1:nl_soil,frnp_(k))*cvsoil_(1:nl_soil,k)*lccpct_np(patchclass_(frnp_(k)))/wgt(1:nl_soil)
                        ldew  (np) = ldew  (np) + ldew_   (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                        sag   (np) = sag   (np) + sag_    (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                        ! use MOD_SnowFraction.F90 to calculate sigf later
                        sigf  (np) = sigf  (np) + sigf_   (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                        wa    (np) = wa    (np) + wa_     (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                     ENDDO

                     ! Get the lowest zwt from source patches and assign to np suggested by Shupeng Zhang
                     zwt(np) = zwt_(frnp_(1))
                     k = 2
                     DO WHILE (k .le. num)
                        IF ( zwt_(frnp_(k)) .lt. zwt(np) ) zwt(np) = zwt_(frnp_(k))
                        k = k + 1
                     ENDDO

                     ! Use restart value from the same type of source patch or remain initialized
                     IF (lccpct_np(patchclass(np)) .gt. 0) THEN
                        tleaf         (np) = tleaf_         (selfnp_)
                        lake_icefrac(:,np) = lake_icefrac_(:,selfnp_)
                        t_lake      (:,np) = t_lake_      (:,selfnp_)
                     ENDIF

                     ! Fraction of soil covered by snow
                     zlnd     = 0.01    !Roughness length for soil [m]
                     fsno(np) = 0.0
                     IF (snowdp(np) > 0.) THEN
                        fmelt = (scv(np)/snowdp(np)/100.) ** m
                        fsno(np)  = tanh(snowdp(np)/(2.5 * zlnd * fmelt))
                     ENDIF

                     ! In case lai+sai come into existence this year, set sigf to 1
                     IF ( (lai(np) + sai(np)).gt.0 .and. sigf(np).eq.0 ) THEN
                        sigf(np) = 1
                     ENDIF

                     ! Set Groud temperature
                     IF ( sum( z_sno(:,np) ) .eq. 0 )  THEN
                        t_grnd(np) = t_soisno(1,np)
                     ELSE
                        DO k = maxsnl+1, 0
                           IF ( z_sno(k,np) .lt. 0 ) THEN
                              t_grnd(np) = t_soisno(k,np)
                              EXIT
                           ENDIF
                        ENDDO
                     ENDIF

                     ! check water balance
                     wbef = 0
                     wpre = 0
                     DO k = 1, num
                        wbef = wbef + ldew_(frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                        wbef = wbef + scv_ (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                        wbef = wbef + wa_  (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                        wbef = wbef + sum(wliq_soisno_(maxsnl+1:nl_soil,frnp_(k)))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                        wbef = wbef + sum(wice_soisno_(maxsnl+1:nl_soil,frnp_(k)))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                     ENDDO
                     wpre = ldew(np) + scv(np) + wa(np) + sum(wliq_soisno(maxsnl+1:nl_soil,np)) + sum(wice_soisno(maxsnl+1:nl_soil,np))
                     IF (wpre-wbef > 1.e-6) THEN
                        print*,'np=',np,'total err=',wpre-wbef
                     ENDIF

                     wbef = 0
                     wpre = 0
                     DO k = 1, num
                        wbef = wbef + sum(wliq_soisno_(maxsnl+1:nl_soil,frnp_(k)))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                        wbef = wbef + sum(wice_soisno_(maxsnl+1:nl_soil,frnp_(k)))*lccpct_np(patchclass_(frnp_(k)))/sum_lccpct_np
                     ENDDO
                     wpre = sum(wliq_soisno(maxsnl+1:nl_soil,np)) + sum(wice_soisno(maxsnl+1:nl_soil,np))
                     IF (wpre-wbef > 1.e-6) THEN
                        print*,'np=',np,'wice+wliq err=',wpre-wbef
                     ENDIF

                  ELSE
                  ! Patch area stay unchanged or decrease, use restart value or remain initialized

                     inp_ = np_
                     DO WHILE (inp_ .le. grid_patch_e_(j))
                        IF (patchclass_(inp_) .eq. patchclass(np)) THEN
                           selfnp_            = inp_
                           frnp_(1)           = inp_
                           wliq_soisno (:,np) = wliq_soisno_ (:,inp_)
                           wice_soisno (:,np) = wice_soisno_ (:,inp_)
                           t_soisno    (:,np) = t_soisno_    (:,inp_)
                           z_sno       (:,np) = z_sno_       (:,inp_)
                           dz_sno      (:,np) = dz_sno_      (:,inp_)
                           t_grnd        (np) = t_grnd_        (inp_)
                           tleaf         (np) = tleaf_         (inp_)
                           ldew          (np) = ldew_          (inp_)
                           sag           (np) = sag_           (inp_)
                           scv           (np) = scv_           (inp_)
                           snowdp        (np) = snowdp_        (inp_)
                           fsno          (np) = fsno_          (inp_)
                           sigf          (np) = sigf_          (inp_)
                           ! In case lai+sai come into existence this year, set sigf to 1
                           IF ( (lai(np) + sai(np)).gt.0 .and. sigf(np).eq.0 ) THEN
                              sigf(np) = 1
                           ENDIF
                           zwt           (np) = zwt_           (inp_)
                           wa            (np) = wa_            (inp_)
                           EXIT
                        ENDIF

                        inp_ = inp_ + 1
                     ENDDO

                  ENDIF

ELSEIF (patchtype(np)==3) THEN !glacier patch
                  ! Used restart value for GLACIERS patches if patchclass exists last year, or remain initialized
                  inp_ = np_
                  DO WHILE (inp_ .le. grid_patch_e_(j))
                     IF (patchclass_(inp_) .eq. patchclass(np)) THEN
                        wliq_soisno (:,np) = wliq_soisno_ (:,inp_)
                        wice_soisno (:,np) = wice_soisno_ (:,inp_)
                        t_soisno    (:,np) = t_soisno_    (:,inp_)
                        z_sno       (:,np) = z_sno_       (:,inp_)
                        dz_sno      (:,np) = dz_sno_      (:,inp_)
                        t_grnd        (np) = t_grnd_        (inp_)
                        tleaf         (np) = tleaf_         (inp_)
                        ldew          (np) = ldew_          (inp_)
                        sag           (np) = sag_           (inp_)
                        scv           (np) = scv_           (inp_)
                        snowdp        (np) = snowdp_        (inp_)
                        fsno          (np) = fsno_          (inp_)
                        sigf          (np) = sigf_          (inp_)
                        zwt           (np) = zwt_           (inp_)
                        wa            (np) = wa_            (inp_)
                        EXIT
                     ENDIF
                     inp_ = inp_ + 1
                  ENDDO
ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
                  IF (patchtype(np)==0) THEN
                     ! if totally come from other types,ldew set to zero since ldew_p(:)=0
                     ldew(np) = sum( ldew_p(patch_pft_s(np):patch_pft_e(np))*pftfrac(patch_pft_s(np):patch_pft_e(np)) )
                  ENDIF

                  IF (patchtype(np)==0 .and. lccpct_np(patchclass(np)) .gt. 0) THEN
                     ! Used restart value of the same pftclass for pft-specific variables
                     ! Note: For ip-specific variables, remain initialized value for new soil patch or pftclass
                     ip = patch_pft_s (np     )
                     ip_= patch_pft_s_(selfnp_)

                     IF (ip.le.0 .or. ip_.le.0) THEN
                        print *, "Error in LuLccMassEnergyConserve LULC_IGBP_PFT|LULC_IGBP_PC!"
                        STOP
                     ENDIF

                     DO WHILE (ip.le.patch_pft_e(np) .and. ip_.le.patch_pft_e_(np_))
                        ! if a PFT is missing, CYCLE
                        IF (pftclass(ip) > pftclass_(ip_)) THEN
                          ip_= ip_+ 1
                          CYCLE
                        ENDIF

                        ! if a PFT is added, CYCLE
                        IF (pftclass(ip) < pftclass_(ip_)) THEN
                          ip = ip + 1
                          CYCLE
                        ENDIF

                        ! for the same PFT, set PFT value
                        tleaf_p    (ip) = tleaf_p_    (ip_)
                        ldew_p     (ip) = ldew_p_     (ip_)
                        ! use MOD_SnowFraction.F90 later
                        sigf_p     (ip) = sigf_p_     (ip_)
                        IF ( (lai_p(ip) + sai_p(ip)).gt.0 .and. sigf_p(ip).eq.0 ) THEN
                           sigf_p(ip) = 1
                        ENDIF

                        ip = ip + 1
                        ip_= ip_+ 1
                     ENDDO
                     ldew(np) = sum( ldew_p(patch_pft_s(np):patch_pft_e(np))*pftfrac(patch_pft_s(np):patch_pft_e(np)) )
                  ENDIF
#endif

#ifdef URBAN_MODEL
                  IF (patchclass(np)==URBAN) THEN
                     ! If there isn't any urban patch in last year's grid,initialized value was remained. Though the first source soil patch would be used for pervious ground related variables.

                     u = patch2urban (np)
                     nurb = count( patchclass_(grid_patch_s_(j):grid_patch_e_(j)) == URBAN )
                     ! print*, 'URB patchclass_: ',patchclass_(grid_patch_s_(j):grid_patch_e_(j))

                     ! Get the index of urban patches in last year's grid, and index of urban patch with the same urbclass
                     IF (nurb > 0) THEN
                        allocate(gu_(nurb)) ! index of urban patches in last year's grid
                        selfu_   = -1       ! index of urban patch with the same urbclass in last year's grid
                        inp_     = np_
                        iu       = 0
                        DO WHILE (inp_ .le. grid_patch_e_(j))
                           IF (patchclass_(inp_) == URBAN) THEN
                              iu = iu + 1
                              gu_(iu) = patch2urban_(inp_)
                              IF (landurban%settyp(u) == urbclass_(gu_(iu))) THEN
                                 selfu_ = gu_(iu)
                              ENDIF
                           ENDIF
                           inp_ = inp_ + 1
                        ENDDO
                     ENDIF

                     ! Index of the same urbclass or the nearest class would be used for new year's assignment
                     IF (selfu_ > 0) THEN
                        u_ = selfu_
                     ELSE IF (nurb > 0) THEN
                        duclass = abs ( landurban%settyp(u) - urbclass_(gu_(1)) )
                        u_ = gu_(1)
                        iu = 2
                        DO WHILE (iu .le. nurb)
                           IF (duclass .gt. abs( landurban%settyp(u) - urbclass_(gu_(iu)) )) THEN
                              u_ = gu_(iu)
                           ENDIF
                        ENDDO
                     ENDIF


                     IF (u.le.0 .or. u_.le.0) THEN
                        print *, "Error in LuLccMassEnergyConserve URBAN_MODEL!"
                        STOP
                     ENDIF

                     fwsun          (u) = fwsun_          (u_)
                     dfwsun         (u) = dfwsun_         (u_)

                     sroof      (:,:,u) = sroof_      (:,:,u_)
                     swsun      (:,:,u) = swsun_      (:,:,u_)
                     swsha      (:,:,u) = swsha_      (:,:,u_)
                     sgimp      (:,:,u) = sgimp_      (:,:,u_)
                     sgper      (:,:,u) = sgper_      (:,:,u_)
                     slake      (:,:,u) = slake_      (:,:,u_)

                     z_sno_roof   (:,u) = z_sno_roof_   (:,u_)
                     z_sno_gimp   (:,u) = z_sno_gimp_   (:,u_)
                     z_sno_gper   (:,u) = z_sno_gper_   (:,u_)
                     z_sno_lake   (:,u) = z_sno_lake_   (:,u_)

                     dz_sno_roof  (:,u) = dz_sno_roof_  (:,u_)
                     dz_sno_gimp  (:,u) = dz_sno_gimp_  (:,u_)
                     dz_sno_gper  (:,u) = dz_sno_gper_  (:,u_)
                     dz_sno_lake  (:,u) = dz_sno_lake_  (:,u_)

                     lwsun          (u) = lwsun_          (u_)
                     lwsha          (u) = lwsha_          (u_)
                     lgimp          (u) = lgimp_          (u_)
                     lgper          (u) = lgper_          (u_)
                     lveg           (u) = lveg_           (u_)

                     t_roofsno    (:,u) = t_roofsno_    (:,u_)
                     t_wallsun    (:,u) = t_wallsun_    (:,u_)
                     t_wallsha    (:,u) = t_wallsha_    (:,u_)
                     t_gimpsno    (:,u) = t_gimpsno_    (:,u_)
                     t_gpersno    (:,u) = t_gpersno_    (:,u_)
                     t_lakesno    (:,u) = t_lakesno_    (:,u_)

                     troof_inner    (u) = troof_inner_    (u_)
                     twsun_inner    (u) = twsun_inner_    (u_)
                     twsha_inner    (u) = twsha_inner_    (u_)

                     wliq_roofsno (:,u) = wliq_roofsno_ (:,u_)
                     wice_roofsno (:,u) = wice_roofsno_ (:,u_)
                     wliq_gimpsno (:,u) = wliq_gimpsno_ (:,u_)
                     wice_gimpsno (:,u) = wice_gimpsno_ (:,u_)
                     wliq_gpersno (:,u) = wliq_gpersno_ (:,u_)
                     wice_gpersno (:,u) = wice_gpersno_ (:,u_)
                     wliq_lakesno (:,u) = wliq_lakesno_ (:,u_)
                     wice_lakesno (:,u) = wice_lakesno_ (:,u_)

                     sag_roof       (u) = sag_roof_       (u_)
                     sag_gimp       (u) = sag_gimp_       (u_)
                     sag_gper       (u) = sag_gper_       (u_)
                     sag_lake       (u) = sag_lake_       (u_)
                     scv_roof       (u) = scv_roof_       (u_)
                     scv_gimp       (u) = scv_gimp_       (u_)
                     scv_gper       (u) = scv_gper_       (u_)
                     scv_lake       (u) = scv_lake_       (u_)
                     fsno_roof      (u) = fsno_roof_      (u_)
                     fsno_gimp      (u) = fsno_gimp_      (u_)
                     fsno_gper      (u) = fsno_gper_      (u_)
                     fsno_lake      (u) = fsno_lake_      (u_)
                     snowdp_roof    (u) = snowdp_roof_    (u_)
                     snowdp_gimp    (u) = snowdp_gimp_    (u_)
                     snowdp_gper    (u) = snowdp_gper_    (u_)
                     snowdp_lake    (u) = snowdp_lake_    (u_)

                     t_room         (u) = t_room_         (u_)
                     tafu           (u) = tafu_           (u_)
                     Fhac           (u) = Fhac_           (u_)
                     Fwst           (u) = Fwst_           (u_)
                     Fach           (u) = Fach_           (u_)

                     ! used soil patch value for variable on pervious ground
                     FROM_SOIL = .false.
                     IF (selfu_ < 0) THEN
                        DO k = 1, num
                           IF (patchtype_(frnp_(k)) == 0) THEN
                              FROM_SOIL = .true.
                           ENDIF
                        ENDDO
                     ENDIF

                     ! Use the first source soil patch temporarily
                     IF (FROM_SOIL) THEN
                        z_sno_gper   (:,u) = z_sno_    (:,frnp_(1))
                        sag_gper       (u) = sag_        (frnp_(1))
                        scv_gper       (u) = scv_        (frnp_(1))
                        fsno_gper      (u) = fsno_       (frnp_(1))
                        snowdp_gper    (u) = snowdp_     (frnp_(1))
                     ENDIF
                  ENDIF
#endif

                  IF (allocated(frnp_  )) deallocate(frnp_  )
                  IF (allocated(gu_    )) deallocate(gu_    )
                  IF (allocated(cvsoil_)) deallocate(cvsoil_)
                  np = np + 1
               ENDDO
            ENDIF
         ENDDO
      ENDDO
   ENDIF

   IF (p_is_worker) THEN
      IF (allocated(grid_patch_s )) deallocate(grid_patch_s )
      IF (allocated(grid_patch_e )) deallocate(grid_patch_e )
      IF (allocated(grid_patch_s_)) deallocate(grid_patch_s_)
      IF (allocated(grid_patch_e_)) deallocate(grid_patch_e_)
      IF (allocated(locpxl       )) deallocate(locpxl       )
   ENDIF

  END SUBROUTINE LulccMassEnergyConserve

END MODULE MOD_Lulcc_MassEnergyConserve
#endif
! ---------- EOP ------------
