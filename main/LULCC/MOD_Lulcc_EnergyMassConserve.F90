#include <define.h>

  SUBROUTINE LulccEnergyConserve
! -------------------------------
! Created by Hua Yuan, 04/2022
! -------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   use MOD_LandPatch
   USE MOD_LandElm
   USE MOD_Mesh
   use MOD_SPMD_Task
   use MOD_Vars_TimeInvariants
   use MOD_Vars_TimeVariables
   USE MOD_Lulcc_Vars_TimeInvariants
   USE MOD_Lulcc_Vars_TimeVariables
   USE MOD_Lulcc_TMatrix
#ifdef LULC_IGBP_PFT
   USE MOD_LandPFT
   USE MOD_Vars_PFTimeInvariants
   USE MOD_Vars_PFTimeVariables
#endif
#ifdef LULC_IGBP_PC
   USE MOD_LandPC
   USE MOD_Vars_PCTimeInvariants
   USE MOD_Vars_PCTimeVariables
#endif
#ifdef URBAN_MODEL
   use MOD_LandUrban
   USE MOD_Urban_Vars_TimeVariables
   USE MOD_Urban_Vars_TimeInvariants
#endif

   IMPLICIT NONE
!TODO: need coding below...

   REAL(r8), allocatable, dimension(:) :: grid_patch_s , grid_patch_e
   REAL(r8), allocatable, dimension(:) :: grid_patch_s_, grid_patch_e_
   INTEGER , allocatable, dimension(:) :: locpxl
   INTEGER :: numpxl,ipxl
   INTEGER, allocatable :: frnp_(:) !index of source patches

   INTEGER  :: k,ilc,num,inp_
   INTEGER  :: i, j, np, np_, selfnp_
#ifdef LULC_IGBP
   INTEGER :: nlc = N_land_classification
   real(r8), dimension(1:N_land_classification) :: lccpct_np
#endif
   real(r8) :: wgt(1:nl_soil)


   IF (p_is_worker) THEN
      ! allocate with numelm
      allocate(grid_patch_s (numelm ))
      allocate(grid_patch_e (numelm ))
      allocate(grid_patch_s_(numelm_))
      allocate(grid_patch_e_(numelm_))

      grid_patch_e (:) = -1.
      grid_patch_s (:) = -1.
      grid_patch_e_(:) = -1.
      grid_patch_s_(:) = -1.

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

#ifdef LULC_IGBP
      DO i=1, numelm
         DO j=1,numelm_
            IF (landelm%eindex(i) == landelm_%eindex(j)) THEN
               np = grid_patch_s (i)
               np_= grid_patch_s_(j)

               IF (np.le.0) CYCLE

               DO WHILE (np.le.grid_patch_e(i))

                  IF ( patchclass(np)==WATERBODY .or. patchclass(np)==GLACIERS) THEN
                     
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
                           fveg          (np) = fveg_          (inp_)
                           fsno          (np) = fsno_          (inp_)
                           sigf          (np) = sigf_          (inp_)
                           green         (np) = green_         (inp_)
                           alb       (:,:,np) = alb_       (:,:,inp_)
                           ssun      (:,:,np) = ssun_      (:,:,inp_)
                           ssha      (:,:,np) = ssha_      (:,:,inp_)
                           thermk        (np) = thermk_        (inp_)
                           extkb         (np) = extkb_         (inp_)
                           extkd         (np) = extkd_         (inp_)
                           zwt           (np) = zwt_           (inp_)
                           wa            (np) = wa_            (inp_)

                           t_lake      (:,np) = t_lake_      (:,inp_)
                           lake_icefrac(:,np) = lake_icefrac_(:,inp_)

                           EXIT
                        ENDIF
                        inp_ = inp_ + 1
                     ENDDO
                     np = np + 1
                     CYCLE
                  ENDIF


                  ! ===============SOIL, URBAN, WETLAND patch==================
                  ! TODO：add case of def_urban_model
                  lccpct_np(:) = lccpct_patches(np,1:nlc)
                  num = count(lccpct_np .gt. 0)

                  allocate ( frnp_ (num) )

                  ! Source patch type which differs from np's type exists
                  IF ( (sum(lccpct_np) - lccpct_np(patchclass(np))) .gt. 0 ) THEN

                     !  Get the index of source patches, and stored as frnp_
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
                     t_soisno    (maxsnl+1:nl_soil,np) = 0  !soil + snow layer te ！mperature [K]
                     z_sno       (:,np)                = 0  !node depth [m]
                     dz_sno      (:,np)                = 0  !interface depth [m]
                     t_grnd        (np)                = 0  !ground surface temperature [K]
                     ldew          (np)                = 0  !depth of water on foliage [mm]
                     sag           (np)                = 0  !non dimensional snow age [-]
                     scv           (np)                = 0  !snow cover, water equivalent [mm]
                     snowdp        (np)                = 0  !snow depth [meter]
                     fveg          (np)                = 0  !fraction of vegetation cover
                     fsno          (np)                = 0  !fraction of snow cover on ground
                     sigf          (np)                = 0  !fraction of veg cover, excluding snow-covered veg [-]
                     green         (np)                = 0  !leaf greenness
                     zwt           (np)                = 0  !the depth to water table [m]
                     wa            (np)                = 0  !water storage in aquifer [mm]
                     tleaf         (np)                = 0  !leaf temperature [K]
                     alb       (:,:,np)                = 0  !averaged albedo [-]
                     ssun      (:,:,np)                = 0  !sunlit canopy absorption for solar radiation (0-1)
                     ssha      (:,:,np)                = 0  !shaded canopy absorption for solar radiation (0-1)
                     thermk        (np)                = 0  !canopy gap fraction for tir radiation
                     extkb         (np)                = 0  !(k, g(mu)/mu) direct solar extinction coefficient
                     extkd         (np)                = 0  !diffuse and scattered diffuse PAR extinction coefficient


                     ! Get the maximum snow layers within all the source patches
                     inp_ = frnp_(1)
                     IF ( num .ge. 2 ) THEN
                        DO k = 2, num
                           IF ( sum(z_sno_(:,frnp_(k))) .gt. sum(z_sno_(:,frnp_(k-1))) ) THEN
                              inp_ = frnp_(k)
                           ENDIF
                        ENDDO
                     ENDIF

                     t_soisno (:0,np) = t_soisno_(:0,inp_)
                     z_sno    (: ,np) = z_sno_   (: ,inp_)
                     dz_sno   (: ,np) = dz_sno_  (: ,inp_)


                     ! Weight of temperature adjustment
                     wgt(1:nl_soil) = 0
                     DO k = 1, num
                        wgt(1:nl_soil) = wgt(1:nl_soil) + cvsoil_(:,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))
                     ENDDO

                     ! Variable adjustment
                     DO k = 1, num
                        wliq_soisno (:,np) = wliq_soisno (:,np) + wliq_soisno_(:,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                        wice_soisno (:,np) = wice_soisno (:,np) + wice_soisno_(:,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                        t_soisno (1:nl_soil,np) = t_soisno (1:nl_soil,np) + t_soisno_(1:nl_soil,frnp_(k))*cvsoil_(:,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/wgt
                        ldew  (np) = ldew  (np) + ldew_   (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                        sag   (np) = sag   (np) + sag_    (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                        scv   (np) = scv   (np) + scv_    (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                        snowdp(np) = snowdp(np) + snowdp_ (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                        fsno  (np) = fsno  (np) + fsno_   (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                        sigf  (np) = sigf  (np) + sigf_   (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                        IF ( (lai(np) + sai(np)) .gt. 0) THEN
                           sigf(np) = 1 - fsno(np)
                        ENDIF
                        zwt   (np) = zwt   (np) + zwt_    (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                        wa    (np) = wa    (np) + wa_     (frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))

                        IF ( green_(frnp_(k)) .eq. 1 ) THEN
                           green (np) = green_(frnp_(k))
                        ENDIF

                     ENDDO


                     ! Use restart value from the same type of source patch
                     IF (lccpct_np(patchclass(np)) .gt. 0) THEN
                        tleaf    (np) = tleaf_    (selfnp_)
                        fveg     (np) = fveg_     (selfnp_)
                        alb  (:,:,np) = alb_  (:,:,selfnp_)
                        ssun (:,:,np) = ssun_ (:,:,selfnp_)
                        ssha (:,:,np) = ssha_ (:,:,selfnp_)
                        thermk   (np) = thermk_   (selfnp_)
                        extkb    (np) = extkb_    (selfnp_)
                        extkd    (np) = extkd_    (selfnp_)

                     ! Area weighted used
                     ELSE
                        ! print*, 'new patch, set area weighted average', 'i=',i, ',j=',j,';patchclass=',patchclass(np)
                        DO k = 1, num
                           tleaf    (np) = tleaf    (np) + tleaf_    (frnp_(k)) *lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                           fveg     (np) = fveg     (np) + fveg_     (frnp_(k)) *lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                           alb  (:,:,np) = alb  (:,:,np) + alb_  (:,:,frnp_(k)) *lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                           ssun (:,:,np) = ssun (:,:,np) + ssun_ (:,:,frnp_(k)) *lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                           ssha (:,:,np) = ssha (:,:,np) + ssha_ (:,:,frnp_(k)) *lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                           thermk   (np) = thermk   (np) + thermk_   (frnp_(k)) *lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                           extkb    (np) = extkb    (np) + extkb_    (frnp_(k)) *lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                           extkd    (np) = extkd    (np) + extkd_    (frnp_(k)) *lccpct_np(patchclass_(frnp_(k)))/sum(lccpct_np(1:nlc))
                        ENDDO

                     ENDIF

                     ! Set Groud temperature
                     IF ( sum( z_sno(:,np) ) .eq. 0 )  THEN
                        t_grnd(np) = t_soisno(1,np)
                     ELSE
                        DO k = maxsnl+1, 0
                           IF ( z_sno(k,np) .le. 0 ) THEN
                              t_grnd(np) = t_soisno(k,np)
                              EXIT
                           ENDIF
                        ENDDO
                     ENDIF

                     deallocate ( frnp_ )

                  ELSE
                  ! Patch area stay unchanged or decrease, use restart value

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
                           fveg          (np) = fveg_          (inp_)
                           fsno          (np) = fsno_          (inp_)
                           sigf          (np) = sigf_          (inp_)
                           green         (np) = green_         (inp_)
                           alb       (:,:,np) = alb_       (:,:,inp_)
                           ssun      (:,:,np) = ssun_      (:,:,inp_)
                           ssha      (:,:,np) = ssha_      (:,:,inp_)
                           thermk        (np) = thermk_        (inp_)
                           extkb         (np) = extkb_         (inp_)
                           extkd         (np) = extkd_         (inp_)
                           zwt           (np) = zwt_           (inp_)
                           wa            (np) = wa_            (inp_)
                           EXIT
                        ENDIF

                        inp_ = inp_ + 1
                     ENDDO

                     deallocate ( frnp_ )
                  ENDIF
                  np = np + 1
               ENDDO
            ENDIF
         ENDDO
      ENDDO
#endif  

   ENDIF

  END SUBROUTINE LulccEnergyConserve
! ---------- EOP ------------
