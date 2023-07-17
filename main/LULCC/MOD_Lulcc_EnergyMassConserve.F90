#include <define.h>

  SUBROUTINE LulccEnergyMassConserve
! -------------------------------
! Created by Hua Yuan, 04/2022
! -------------------------------

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
   USE MOD_LandUrban
   USE MOD_Urban_Vars_TimeVariables
   USE MOD_Urban_Vars_TimeInvariants
#endif

   IMPLICIT NONE
!TODO: need coding below...

   INTEGER, allocatable, dimension(:) :: grid_patch_s , grid_patch_e
   INTEGER, allocatable, dimension(:) :: grid_patch_s_, grid_patch_e_
   INTEGER, allocatable, dimension(:) :: locpxl
   INTEGER :: numpxl,ipxl

   INTEGER, allocatable :: frnp_(:) !index of source patches
   INTEGER, allocatable :: gu_(:)   !index of urban patches in last year's grid

   INTEGER :: k, ilc, num, inp_
   INTEGER :: i, j, np, np_, selfnp_, ipft, ip, ip_, pc, pc_
   INTEGER :: u, u_, iu, selfu_, nurb, duclass
   INTEGER :: nlc = N_land_classification
   REAL(r8), dimension(1:N_land_classification) :: lccpct_np
   REAL(r8):: wgt(1:nl_soil)
   LOGICAL :: FROM_SOIL

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

                  IF ( patchclass(np)==WATERBODY .or. patchclass(np)==GLACIERS) THEN
                     ! Used restart value for WATERBODY and GLACIERS patches if patchclass exists last year
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
#if (defined LULC_IGBP || defined LULC_IGBP_PC)
                  lccpct_np(:) = lccpct_patches(np,1:nlc)
#endif

#ifdef LULC_IGBP_PFT
                  lccpct_np(:) = 0
                  lccpct_np(1) = sum(lccpct_patches(np,1:10)) + lccpct_patches(np,14) &
                                 +   lccpct_patches(np,12)    + lccpct_patches(np,16)
                  lccpct_np(URBAN  ) = lccpct_patches(np,URBAN  )
                  lccpct_np(WETLAND) = lccpct_patches(np,WETLAND)
#endif

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
                        IF (any(wgt < 1.e-6)) THEN
                           ! e.g. WATERBODY -> SOIL, cvsoil_=0, wgt=0, need to adjust
                           ! Use the temperature of the previous soil patch
                           ! TODO: used the cvsoil underneath lake
                           t_soisno (1:nl_soil,np) = t_soisno (1:nl_soil,np-1)
                        ELSE
                           t_soisno (1:nl_soil,np) = t_soisno (1:nl_soil,np) + t_soisno_(1:nl_soil,frnp_(k))*cvsoil_(:,frnp_(k))*lccpct_np(patchclass_(frnp_(k)))/wgt
                        ENDIF
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

                  ELSE
                  ! Patch area stay unchanged or decrease, use restart value

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
                           fveg          (np) = fveg_          (inp_)
                           fsno          (np) = fsno_          (inp_)
                           sigf          (np) = sigf_          (inp_)
                           IF ( (lai(np) + sai(np)) .gt. 1e-6) THEN
                              sigf(np) = 1 - fsno(np)
                           ENDIF
                           green         (np) = green_         (inp_)
                           alb       (:,:,np) = alb_       (:,:,inp_)
                           ssun      (:,:,np) = ssun_      (:,:,inp_)
                           ssha      (:,:,np) = ssha_      (:,:,inp_)
                           extkb         (np) = extkb_         (inp_)
                           extkd         (np) = extkd_         (inp_)
                           zwt           (np) = zwt_           (inp_)
                           wa            (np) = wa_            (inp_)
                           EXIT
                        ENDIF

                        inp_ = inp_ + 1
                     ENDDO

                  ENDIF


#ifdef LULC_IGBP_PFT
                  IF (patchtype(np)==0 .and. lccpct_np(patchclass(np)) .gt. 0) THEN
                     ! Used restart value of the same pftclass for pft-specific variables
                     ! Note: For ip-specific variables, remain initialized value for new soil patch or pftclass
                     ip = patch_pft_s (np     )
                     ip_= patch_pft_s_(selfnp_)

                     IF (ip.le.0 .or. ip_.le.0) THEN
                        print *, "Error in LuLccEnergyMassConserve LULC_IGBP_PFT!"
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
                        sigf_p     (ip) = sigf_p_     (ip_)
                        lai_p      (ip) = lai_p_      (ip_)
                        sai_p      (ip) = sai_p_      (ip_)
                        IF ( (lai_p(ip) + sai_p(ip)) .gt. 0) THEN
                           sigf_p(ip) = 1
                        ENDIF
                        ssun_p (:,:,ip) = ssun_p_ (:,:,ip_)
                        ssha_p (:,:,ip) = ssha_p_ (:,:,ip_)
                        thermk_p   (ip) = thermk_p_   (ip_)
                        extkb_p    (ip) = extkb_p_    (ip_)
                        extkd_p    (ip) = extkd_p_    (ip_)

                        ip = ip + 1
                        ip_= ip_+ 1
                     ENDDO
                  ENDIF

#endif

#ifdef LULC_IGBP_PC
                  ! Note: For pc-specific variables, remain initialized value for newly appear soil patch
                  IF ( patchtype(np)==0 .and. lccpct_np(patchclass(np)) .gt. 0) THEN

                     pc = patch2pc (np )
                     pc_= patch2pc_(selfnp_)

                     IF (pc.le.0 .or. pc_.le.0) THEN
                        print *, "Error in LuLccEnergyMassConserve LULC_IGBP_PC!"
                        STOP
                     ENDIF

                     ! for the same patch TYPE
                     tleaf_c    (:,pc) = tleaf_c_    (:,pc_)
                     ldew_c     (:,pc) = ldew_c_     (:,pc_)
                     sigf_c     (:,pc) = sigf_c_     (:,pc_)
                     lai_c      (:,pc) = lai_c_      (:,pc_)
                     sai_c      (:,pc) = sai_c_      (:,pc_)
                     DO ipft = 0, N_PFT-1
                        IF ( (lai_c(ipft,pc) + sai_c(ipft,pc)) .gt. 0) THEN
                           sigf_c(ipft,pc) = 1
                        ENDIF
                     ENDDO
                     ssun_c (:,:,:,pc) = ssun_c_ (:,:,:,pc_)
                     ssha_c (:,:,:,pc) = ssha_c_ (:,:,:,pc_)
                     thermk_c   (:,pc) = thermk_c_   (:,pc_)
                     ! Note: should use initialize value because pcfrac(:,pc) could change next year
                     ! fshade_c   (:,pc) = fshade_c_   (:,pc_)
                     extkb_c    (:,pc) = extkb_c_    (:,pc_)
                     extkd_c    (:,pc) = extkd_c_    (:,pc_)

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
                        print *, "Error in LuLccEnergyMassConserve URBAN_MODEL!"
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
                     DO k = 1, num
                        IF (patchtype_(frnp_(k)) == 0) THEN
                           FROM_SOIL = .true.
                        ENDIF
                     ENDDO

                     ! Use the first source soil patch temporarily
                     IF (FROM_SOIL) THEN
                        sgper      (:,:,u) = ssun_   (:,:,frnp_(1)) + ssha_(:,:,frnp_(1)) !check if is the same meaning
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

  END SUBROUTINE LulccEnergyMassConserve
! ---------- EOP ------------
