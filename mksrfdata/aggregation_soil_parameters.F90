#include <define.h>

SUBROUTINE aggregation_soil_parameters ( &
      gland, dir_rawdata, dir_model_landdata)
   ! ----------------------------------------------------------------------
   !DESCRIPTION:
   !Create soil hydraulic and thermal parameters for the modeling reolustion
   !
   !ORIGINAL: 
   !Yongjiu Dai and Wei Shangguan, 02/2014
   !
   !REFERENCES:
   ! 1)Dai, Y., Q. Xin, N. Wei, Y. Zhang, W. Shangguan, H. Yuan, S. Zhang, S. Liu, X. Lu, 2019. A global high-resolution dataset of
   !          soil hydraulic and thermal properties for land surface modeling. Journal of Advances in Modeling Earth Systems,11, 2996-3023.
   ! 2)Dai, Y., N. Wei, H. Yuan, S. Zhang, W. Shangguan, S. Liu, and X. Lu, 2019. Evaluation of soil thermal conductivity schemes
   !          for use in land surface modelling, Journal of Advances in Modeling Earth Systems, 11, 3454-3473.
   ! 3)Dai, Y., W. Shangguan, Q. Duan, B. Liu, S. Fu, and G. Niu, 2013. Development of a China dataset of soil hydraulic parameters 
   !         using pedotransfer functions for land surface modeling. Journal of Hydrometeorology 14, 869–887
   !
   !REVISIONS:
   !Nan Wei, 06/2019, 02/2020: add three subroutines to  calculate the function/jacobian matrix
   !
   ! !USES:
   USE precision
   USE GlobalVars
   USE mod_namelist
   USE spmd_task
   USE mod_grid
   USE mod_landpatch
   USE ncio_block
   USE ncio_vector
   USE mod_aggregation
#ifdef CLMDEBUG 
   USE mod_colm_debug
#endif
   USE mod_utils
#ifdef SOILPAR_UPS_FIT
   USE par_fitting
#endif
#ifdef SinglePoint
   USE mod_single_srfdata
#endif

   IMPLICIT NONE

   ! arguments:
   TYPE(grid_type),  intent(in) :: gland
   CHARACTER(LEN=*), intent(in) :: dir_rawdata
   CHARACTER(LEN=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   CHARACTER(len=256) :: landdir, lndname
   CHARACTER(len=256) :: c
   INTEGER :: nsl, ipatch, L, np, LL, ipxstt, ipxend

   TYPE (block_data_real8_2d) :: vf_quartz_mineral_s_grid
   TYPE (block_data_real8_2d) :: vf_gravels_s_grid
   TYPE (block_data_real8_2d) :: vf_om_s_grid
   TYPE (block_data_real8_2d) :: vf_sand_s_grid
   TYPE (block_data_real8_2d) :: wf_gravels_s_grid
   TYPE (block_data_real8_2d) :: wf_sand_s_grid
   TYPE (block_data_real8_2d) :: OM_density_s_grid
   TYPE (block_data_real8_2d) :: BD_all_s_grid
   TYPE (block_data_real8_2d) :: theta_s_grid
   TYPE (block_data_real8_2d) :: psi_s_grid  
   TYPE (block_data_real8_2d) :: lambda_grid 
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   TYPE (block_data_real8_2d) :: theta_r_grid 
   TYPE (block_data_real8_2d) :: alpha_vgm_grid 
   TYPE (block_data_real8_2d) :: L_vgm_grid 
   TYPE (block_data_real8_2d) :: n_vgm_grid 
#endif
   TYPE (block_data_real8_2d) :: k_s_grid    
   TYPE (block_data_real8_2d) :: csol_grid   
   TYPE (block_data_real8_2d) :: tksatu_grid
   TYPE (block_data_real8_2d) :: tksatf_grid 
   TYPE (block_data_real8_2d) :: tkdry_grid 
   TYPE (block_data_real8_2d) :: k_solids_grid

   REAL(r8), allocatable :: vf_quartz_mineral_s_patches (:)
   REAL(r8), allocatable :: vf_gravels_s_patches (:)
   REAL(r8), allocatable :: vf_om_s_patches (:)
   REAL(r8), allocatable :: vf_sand_s_patches (:)
   REAL(r8), allocatable :: wf_gravels_s_patches (:)
   REAL(r8), allocatable :: wf_sand_s_patches (:)
   REAL(r8), allocatable :: OM_density_s_patches (:)
   REAL(r8), allocatable :: BD_all_s_patches (:)
   REAL(r8), allocatable :: theta_s_patches (:) 
   REAL(r8), allocatable :: psi_s_patches   (:) 
   REAL(r8), allocatable :: lambda_patches  (:) 
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   REAL(r8), allocatable :: theta_r_patches (:) 
   REAL(r8), allocatable :: alpha_vgm_patches  (:) 
   REAL(r8), allocatable :: L_vgm_patches  (:) 
   REAL(r8), allocatable :: n_vgm_patches  (:) 
#endif
   REAL(r8), allocatable :: k_s_patches     (:) 
   REAL(r8), allocatable :: csol_patches    (:) 
   REAL(r8), allocatable :: tksatu_patches  (:)
   REAL(r8), allocatable :: tksatf_patches  (:) 
   REAL(r8), allocatable :: tkdry_patches   (:) 
   REAL(r8), allocatable :: k_solids_patches  (:)
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
   REAL(r8), allocatable :: BA_alpha_patches  (:)
   REAL(r8), allocatable :: BA_beta_patches  (:)
#endif

   REAL(r8), allocatable :: vf_quartz_mineral_s_one (:)
   REAL(r8), allocatable :: vf_gravels_s_one (:)
   REAL(r8), allocatable :: vf_om_s_one (:)
   REAL(r8), allocatable :: vf_sand_s_one (:)
   REAL(r8), allocatable :: wf_gravels_s_one (:)
   REAL(r8), allocatable :: wf_sand_s_one (:) 
   REAL(r8), allocatable :: OM_density_s_one (:)
   REAL(r8), allocatable :: BD_all_s_one (:)
   REAL(r8), allocatable :: theta_s_one (:) 
   REAL(r8), allocatable :: psi_s_one   (:) 
   REAL(r8), allocatable :: lambda_one  (:) 
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   REAL(r8), allocatable :: theta_r_one  (:) 
   REAL(r8), allocatable :: alpha_vgm_one  (:) 
   REAL(r8), allocatable :: L_vgm_one  (:) 
   REAL(r8), allocatable :: n_vgm_one  (:) 
#endif
   REAL(r8), allocatable :: k_s_one     (:) 
   REAL(r8), allocatable :: csol_one    (:) 
   REAL(r8), allocatable :: tksatu_one  (:)
   REAL(r8), allocatable :: tksatf_one  (:) 
   REAL(r8), allocatable :: tkdry_one   (:) 
   REAL(r8), allocatable :: k_solids_one  (:)
   REAL(r8), allocatable :: area_one   (:)
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
   REAL(r8), allocatable :: BA_alpha_one  (:)
   REAL(r8), allocatable :: BA_beta_one  (:)
#endif

#ifdef SOILPAR_UPS_FIT
! local variables for estimating the upscaled soil parameters using the Levenberg–Marquardt fitting method
! ---------------------------------------------------------------
   integer, parameter   :: npointw  = 24      
   integer, parameter   :: npointb  = 20
   real(r8),parameter   :: xdat(npointw) = (/1.,5.,10.,20.,30.,40.,50.,60.,70.,90.,110.,130.,150.,&
                                170.,210.,300.,345.,690.,1020.,5100.,15300.,20000.,100000.,1000000./)
                           !  points of soil pressure heads used for fitting SW retention curves

   real(r8),parameter   :: xdatsr(npointb)=(/0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,&
                                             0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0/)
                           !  points of soil saturation levels (Sr) used for fitting Ke-Sr relationship

   real(r8),allocatable :: ydatc(:,:)     ! the Campbell SW retentions at fine grids 
   real(r8),allocatable :: ydatv(:,:)     ! the van Genuchten SW retentions at fine grids
   real(r8),allocatable :: ydatb(:,:)     ! the Balland and Arp (2005) Ke-Sr relationship at fine grids

   integer, parameter   :: nc = 2         ! the number of fitted parameters in Campbell SW retention curve (psi and lambda)
   integer, parameter   :: nv = 3         ! the number of fitted parameters in van Genuchten SW retention curve 
                                          ! (theta_r, alpha and n)
   integer, parameter   :: nb = 2         ! the number of fitted parameters in Ke-Sr relationship (alpha and beta)

! Variables needed for Levenberg–Marquardt algorithm in MINPACK library      
   real(r8),parameter   :: factor = 0.1
   real(r8),parameter   :: ftol = 1.0e-5
   real(r8),parameter   :: xtol = 1.0e-4
   real(r8),parameter   :: gtol = 0.0
   integer, parameter   :: mode = 1
   integer, parameter   :: nprint = 0
   integer              :: ldfjac,info,ipvtc(nc),ipvtv(nv),ipvtb(nb),maxfev,nfev,njev
   real(r8)             :: xc(nc),xv(nv),xb(nb),diagc(nc),diagv(nv),diagb(nb),qtfc(nc),qtfv(nv),qtfb(nb)
   real(r8),allocatable :: fjacc(:,:),fvecc(:),fjacv(:,:),fvecv(:),fjacb(:,:),fvecb(:)
   integer isiter                         ! flags to tell whether the iteration is completed, 1=Yes, 0=No

   external SW_CB_dist                    ! the objective function to be fitted for Campbell SW retention curve
   external SW_VG_dist                    ! the objective function to be fitted for van Genuchten SW retention curve
!   external Ke_Sr_dist                    ! the objective function to be fitted for Balland and Arp (2005) Ke-Sr relationship
#endif
      
   landdir = trim(dir_model_landdata) // '/soil'

   ! ........................................
   ! ... [2] aggregate the soil parameters from the resolution of raw data to modelling resolution
   ! ........................................
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A29)') 'Aggregate Soil Parameters ...'
      CALL system('mkdir -p ' // trim(adjustl(landdir)))
   ENDIF
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
      
#ifdef SinglePoint
   IF (USE_SITE_soilparameters) THEN
      RETURN
   ELSE
      allocate ( SITE_soil_vf_quartz_mineral (nl_soil) )
      allocate ( SITE_soil_vf_gravels        (nl_soil) )
      allocate ( SITE_soil_vf_om             (nl_soil) )
      allocate ( SITE_soil_vf_sand           (nl_soil) )
      allocate ( SITE_soil_wf_gravels        (nl_soil) )
      allocate ( SITE_soil_wf_sand           (nl_soil) )
      allocate ( SITE_soil_OM_density        (nl_soil) )
      allocate ( SITE_soil_BD_all            (nl_soil) )
      allocate ( SITE_soil_theta_s           (nl_soil) )
      allocate ( SITE_soil_psi_s             (nl_soil) )
      allocate ( SITE_soil_lambda            (nl_soil) )
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      allocate ( SITE_soil_theta_r   (nl_soil) )
      allocate ( SITE_soil_alpha_vgm (nl_soil) )
      allocate ( SITE_soil_L_vgm     (nl_soil) )
      allocate ( SITE_soil_n_vgm     (nl_soil) )
#endif
      allocate ( SITE_soil_k_s      (nl_soil) )
      allocate ( SITE_soil_csol     (nl_soil) )
      allocate ( SITE_soil_tksatu   (nl_soil) )
      allocate ( SITE_soil_tksatf   (nl_soil) )
      allocate ( SITE_soil_tkdry    (nl_soil) )
      allocate ( SITE_soil_k_solids (nl_soil) )
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
      allocate ( SITE_soil_BA_alpha (nl_soil) )
      allocate ( SITE_soil_BA_beta  (nl_soil) )
#endif
   ENDIF
#endif

   IF (p_is_worker) THEN

      allocate ( vf_quartz_mineral_s_patches(numpatch) )
      allocate ( vf_gravels_s_patches       (numpatch) )
      allocate ( vf_om_s_patches            (numpatch) )
      allocate ( vf_sand_s_patches          (numpatch) )
      allocate ( wf_gravels_s_patches       (numpatch) )
      allocate ( wf_sand_s_patches          (numpatch) )
      allocate ( OM_density_s_patches       (numpatch) )
      allocate ( BD_all_s_patches           (numpatch) )
      allocate ( theta_s_patches            (numpatch) )
      allocate ( psi_s_patches              (numpatch) )
      allocate ( lambda_patches             (numpatch) )
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      allocate ( theta_r_patches   (numpatch) )
      allocate ( alpha_vgm_patches (numpatch) )
      allocate ( L_vgm_patches     (numpatch) )
      allocate ( n_vgm_patches     (numpatch) )
#endif
      allocate ( k_s_patches       (numpatch) )
      allocate ( csol_patches      (numpatch) )
      allocate ( tksatu_patches    (numpatch) )
      allocate ( tksatf_patches    (numpatch) )
      allocate ( tkdry_patches     (numpatch) )
      allocate ( k_solids_patches  (numpatch) )
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
      allocate ( BA_alpha_patches  (numpatch) )
      allocate ( BA_beta_patches   (numpatch) )
#endif

   ENDIF

   DO nsl = 1, 8

      write(c,'(i1)') nsl

      ! (1) volumetric fraction of quartz within mineral soil
      IF (p_is_io) THEN
        
         CALL allocate_block_data (gland, vf_quartz_mineral_s_grid)
         lndname = trim(dir_rawdata)//'/soil/vf_quartz_mineral_s.nc'
         CALL ncio_read_block (lndname, 'vf_quartz_mineral_s_l'//trim(c), gland, vf_quartz_mineral_s_grid)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = vf_quartz_mineral_s_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                  data_r8_2d_in1 = vf_quartz_mineral_s_grid, data_r8_2d_out1 = vf_quartz_mineral_s_one)
               vf_quartz_mineral_s_patches (ipatch) = sum (vf_quartz_mineral_s_one * (area_one/sum(area_one)))
            ELSE
               vf_quartz_mineral_s_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(vf_quartz_mineral_s_patches(ipatch))) THEN
               write(*,*) "NAN appears in vf_quartz_mineral_s_patches."
               write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('vf_quartz_mineral_s lev '//trim(c), vf_quartz_mineral_s_patches)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/vf_quartz_mineral_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'vf_quartz_mineral_s_l'//trim(c)//'_patches', 'patch',& 
                              landpatch, vf_quartz_mineral_s_patches, 1)
#else
      SITE_soil_vf_quartz_mineral(nsl) = vf_quartz_mineral_s_patches(1)
#endif

      ! (2) volumetric fraction of gravels
      ! (3) volumetric fraction of sand
      ! (4) volumetric fraction of organic matter
      ! with the parameter alpha and beta in the Balland V. and P. A. Arp (2005) model if defined THERMAL_CONDUCTIVITY_SCHEME_4 
      IF (p_is_io) THEN
        
         CALL allocate_block_data (gland, vf_gravels_s_grid)
         lndname = trim(dir_rawdata)//'/soil/vf_gravels_s.nc'
         CALL ncio_read_block (lndname, 'vf_gravels_s_l'//trim(c), gland, vf_gravels_s_grid)

         CALL allocate_block_data (gland, vf_sand_s_grid)
         lndname = trim(dir_rawdata)//'/soil/vf_sand_s.nc'
         CALL ncio_read_block (lndname, 'vf_sand_s_l'//trim(c), gland, vf_sand_s_grid)

         CALL allocate_block_data (gland, vf_om_s_grid)
         lndname = trim(dir_rawdata)//'/soil/vf_om_s.nc'
         CALL ncio_read_block (lndname, 'vf_om_s_l'//trim(c), gland, vf_om_s_grid)

#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = vf_gravels_s_grid, &
            data_r8_2d_in2 = vf_sand_s_grid,  data_r8_2d_in3 = vf_om_s_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                  data_r8_2d_in1 = vf_gravels_s_grid, data_r8_2d_out1 = vf_gravels_s_one, &
                  data_r8_2d_in2 = vf_sand_s_grid,    data_r8_2d_out2 = vf_sand_s_one, &
                  data_r8_2d_in3 = vf_om_s_grid,      data_r8_2d_out3 = vf_om_s_one)

               vf_gravels_s_patches (ipatch) = sum (vf_gravels_s_one * (area_one/sum(area_one)))
               vf_sand_s_patches (ipatch) = sum (vf_sand_s_one * (area_one/sum(area_one)))
               vf_om_s_patches (ipatch) = sum (vf_om_s_one * (area_one/sum(area_one)))

#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
               ! the parameter values of Balland and Arp (2005) Ke-Sr relationship,
               ! modified by Barry-Macaulay et al.(2015), Evaluation of soil thermal conductivity models

               ipxstt = landpatch%ipxstt(ipatch)
               ipxend = landpatch%ipxend(ipatch)

               allocate(BA_alpha_one  (ipxstt:ipxend))
               allocate(BA_beta_one   (ipxstt:ipxend))
               where ((vf_gravels_s_one + vf_sand_s_one) > 0.4)
                       BA_alpha_one = 0.38
                       BA_beta_one = 35.0
               elsewhere ((vf_gravels_s_one + vf_sand_s_one) > 0.25)
                       BA_alpha_one = 0.24
                       BA_beta_one = 26.0
               elsewhere
                       BA_alpha_one = 0.2
                       BA_beta_one = 10.0
               end where

               BA_alpha_patches (ipatch) = median (BA_alpha_one, size(BA_alpha_one), spval)
               BA_beta_patches (ipatch) = median (BA_beta_one, size(BA_beta_one), spval)

#ifdef SOILPAR_UPS_FIT
!               np = size(BA_alpha_one)
!               IF( np > 1 ) then
!                  allocate ( ydatb  (ipxstt:ipxend,npointb) )
!! the jacobian matrix required in Levenberg–Marquardt fitting method
!                  allocate ( fjacb  (npointb,nb) )           ! calculated in Ke_Sr_dist
!! the values of objective functions to be fitted 
!                  allocate ( fvecb  (npointb)    )           ! calculated in Ke_Sr_dist 

!! Ke-Sr relationship at fine grids for each patch
!                  do LL = ipxstt,ipxend
!                     ydatb(LL,:) = xdatsr**(0.5*(1.0+vf_om_s_one(LL)-BA_alpha_one(LL)*vf_sand_s_one(LL) &
!                                 - vf_gravels_s_one(LL))) * ((1.0/(1.0+exp(-BA_beta_one(LL)*xdatsr)))**3 &
!                                 - ((1.0-xdatsr)/2.0)**3)**(1.0-vf_om_s_one(LL))
!                  end do                  
                  
!! Fitting the parameters in the Balland and Arp (2005) Ke-Sr relationship            
!                  ldfjac = npointb
!                  xb(1) = BA_alpha_patches (ipatch)  
!                  xb(2) = BA_beta_patches (ipatch)
!                  maxfev = 100 * ( nb + 1 )
!                  isiter = 1

!                  call lmder ( Ke_Sr_dist, npointb, nb, xb, fvecb, fjacb, ldfjac, ftol, xtol, gtol, maxfev, &
!                        diagb, mode, factor, nprint, info, nfev, njev, ipvtb, qtfb,&
!                        xdatsr,npointb,ydatb, np, vf_gravels_s_patches (ipatch), isiter,&
!                        vf_om_s_patches(ipatch),vf_sand_s_patches(ipatch),vf_gravels_s_patches(ipatch))

!                  ydatb(ipxstt,:) = xdatsr**(0.5*(1.0+vf_om_s_patches(ipatch)-xb(1)*vf_sand_s_patches(ipatch) &
!                             - vf_gravels_s_patches(ipatch))) * ((1.0/(1.0+exp(-xb(2)*xdatsr)))**3 &
!                             - ((1.0-xdatsr)/2.0)**3)**(1.0-vf_om_s_patches(ipatch))

!                  if ( all(ydatb(ipxstt,:) >= 0.) .and. all(ydatb(ipxstt,:) <= 1.) .and. isiter == 1 ) then
!                       BA_alpha_patches(ipatch) = xb(1)
!                       BA_beta_patches (ipatch) = xb(2)
!                  end if

!                  deallocate(ydatb)
!                  deallocate(fjacb)
!                  deallocate(fvecb)

!               ENDIF 
#endif
               deallocate(BA_alpha_one)
               deallocate(BA_beta_one)
#endif

            ELSE
               vf_gravels_s_patches (ipatch) = -1.0e36_r8
               vf_sand_s_patches (ipatch) = -1.0e36_r8
               vf_om_s_patches (ipatch) = -1.0e36_r8
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
               BA_alpha_patches (ipatch) = -1.0e36_r8
               BA_beta_patches (ipatch) = -1.0e36_r8
#endif
            ENDIF

            IF (isnan(vf_gravels_s_patches(ipatch))) THEN
               write(*,*) "NAN appears in vf_gravels_s_patches."
               write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

            IF (isnan(vf_sand_s_patches(ipatch))) THEN
               write(*,*) "NAN appears in vf_sand_s_patches."
               write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

            IF (isnan(vf_om_s_patches(ipatch))) THEN
               write(*,*) "NAN appears in vf_om_s_patches."
               write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('vf_gravels_s lev '//trim(c), vf_gravels_s_patches)
      CALL check_vector_data ('vf_sand_s lev '//trim(c), vf_sand_s_patches)
      CALL check_vector_data ('vf_om_s lev '//trim(c), vf_om_s_patches)
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
      CALL check_vector_data ('BA_alpha lev '//trim(c), BA_alpha_patches)
      CALL check_vector_data ('BA_beta lev '//trim(c), BA_beta_patches)
#endif
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/vf_gravels_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'vf_gravels_s_l'//trim(c)//'_patches', 'patch',& 
                              landpatch, vf_gravels_s_patches, 1)
#else
      SITE_soil_vf_gravels(nsl) = vf_gravels_s_patches(1)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/vf_sand_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'vf_sand_s_l'//trim(c)//'_patches', 'patch',&
                              landpatch, vf_sand_s_patches, 1)
#else
      SITE_soil_vf_sand(nsl) = vf_sand_s_patches(1)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/vf_om_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'vf_om_s_l'//trim(c)//'_patches', 'patch',&
                              landpatch, vf_om_s_patches, 1)
#else
      SITE_soil_vf_om(nsl) = vf_om_s_patches(1)
#endif

#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
#ifndef SinglePoint
      lndname = trim(landdir)//'/BA_alpha_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'BA_alpha_l'//trim(c)//'_patches', 'patch',& 
                              landpatch, BA_alpha_patches, 1)
#else
      SITE_soil_BA_alpha(nsl) = BA_alpha_patches(1)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/BA_beta_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'BA_beta_l'//trim(c)//'_patches', 'patch',&
                              landpatch, BA_beta_patches, 1)
#else
      SITE_soil_BA_beta(nsl) = BA_beta_patches(1)
#endif
#endif

      ! (5) gravimetric fraction of gravels
      IF (p_is_io) THEN
        
         CALL allocate_block_data (gland, wf_gravels_s_grid)
         lndname = trim(dir_rawdata)//'/soil/wf_gravels_s.nc'
         CALL ncio_read_block (lndname, 'wf_gravels_s_l'//trim(c), gland, wf_gravels_s_grid)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = wf_gravels_s_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                  data_r8_2d_in1 = wf_gravels_s_grid, data_r8_2d_out1 = wf_gravels_s_one)
               wf_gravels_s_patches (ipatch) = sum (wf_gravels_s_one * (area_one/sum(area_one)))
            ELSE
               wf_gravels_s_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(wf_gravels_s_patches(ipatch))) THEN
               write(*,*) "NAN appears in wf_gravels_s_patches."
               write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('wf_gravels_s lev '//trim(c), wf_gravels_s_patches)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/wf_gravels_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'wf_gravels_s_l'//trim(c)//'_patches', 'patch',& 
                              landpatch, wf_gravels_s_patches, 1)
#else
      SITE_soil_wf_gravels(nsl) = wf_gravels_s_patches(1)
#endif


      ! (6) gravimetric fraction of sand
      IF (p_is_io) THEN
        
         CALL allocate_block_data (gland, wf_sand_s_grid)
         lndname = trim(dir_rawdata)//'/soil/wf_sand_s.nc'
         CALL ncio_read_block (lndname, 'wf_sand_s_l'//trim(c), gland, wf_sand_s_grid)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = wf_sand_s_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                  data_r8_2d_in1 = wf_sand_s_grid, data_r8_2d_out1 = wf_sand_s_one)
               wf_sand_s_patches (ipatch) = sum (wf_sand_s_one * (area_one/sum(area_one)))
            ELSE
               wf_sand_s_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(wf_sand_s_patches(ipatch))) THEN
               write(*,*) "NAN appears in wf_sand_s_patches."
               write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('wf_sand_s lev '//trim(c), wf_sand_s_patches)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/wf_sand_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'wf_sand_s_l'//trim(c)//'_patches', 'patch',& 
                              landpatch, wf_sand_s_patches, 1)
#else
      SITE_soil_wf_sand(nsl) = wf_sand_s_patches(1)
#endif


#ifdef vanGenuchten_Mualem_SOIL_MODEL                           
      ! (7) VGM's pore-connectivity parameter (L)
      IF (p_is_io) THEN
        
         CALL allocate_block_data (gland, L_vgm_grid)
         lndname = trim(dir_rawdata)//'/soil/VGM_L.nc'
         CALL ncio_read_block (lndname, 'VGM_L_l'//trim(c), gland, L_vgm_grid)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = L_vgm_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, &
                  data_r8_2d_in1 = L_vgm_grid, data_r8_2d_out1 = L_vgm_one)
               L_vgm_patches (ipatch) = median (L_vgm_one, size(L_vgm_one), spval)
            ELSE
               L_vgm_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(L_vgm_patches(ipatch))) THEN
               write(*,*) "NAN appears in L_vgm_patches."
               write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('L VGM lev '//trim(c), L_vgm_patches)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/L_vgm_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'L_vgm_l'//trim(c)//'_patches', 'patch', landpatch, L_vgm_patches, 1)
#else
      SITE_soil_L_vgm(nsl) = L_vgm_patches(1)
#endif


      ! (8) VGM's residual water content (theta_r) [cm3/cm3]
      ! (9) VGM's parameter corresponding approximately to the inverse of the air-entry value (alpha)
      ! (10) VGM's shape parameter (n)
      ! (11) saturated water content [cm3/cm3]

      IF (p_is_io) THEN
        
         CALL allocate_block_data (gland, theta_r_grid)
         lndname = trim(dir_rawdata)//'/soil/VGM_theta_r.nc'
         CALL ncio_read_block (lndname, 'VGM_theta_r_l'//trim(c), gland, theta_r_grid)

         CALL allocate_block_data (gland, alpha_vgm_grid)
         lndname = trim(dir_rawdata)//'/soil/VGM_alpha.nc'
         CALL ncio_read_block (lndname, 'VGM_alpha_l'//trim(c), gland, alpha_vgm_grid)

         CALL allocate_block_data (gland, n_vgm_grid)
         lndname = trim(dir_rawdata)//'/soil/VGM_n.nc'
         CALL ncio_read_block (lndname, 'VGM_n_l'//trim(c), gland, n_vgm_grid)

         CALL allocate_block_data (gland, theta_s_grid)
         lndname = trim(dir_rawdata)//'/soil/theta_s.nc'
         CALL ncio_read_block (lndname, 'theta_s_l'//trim(c), gland, theta_s_grid)

#ifdef USEMPI
         CALL aggregation_data_daemon (gland, &
            data_r8_2d_in1 = theta_r_grid, data_r8_2d_in2 = alpha_vgm_grid, &
            data_r8_2d_in3 = n_vgm_grid,   data_r8_2d_in4 = theta_s_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                      data_r8_2d_in1 = theta_r_grid,   data_r8_2d_out1 = theta_r_one, &
                      data_r8_2d_in2 = alpha_vgm_grid, data_r8_2d_out2 = alpha_vgm_one, &
                      data_r8_2d_in3 = n_vgm_grid,     data_r8_2d_out3 = n_vgm_one, &
                      data_r8_2d_in4 = theta_s_grid,   data_r8_2d_out4 = theta_s_one)
               theta_r_patches (ipatch) = median (theta_r_one, size(theta_r_one), spval)
               alpha_vgm_patches (ipatch) = median (alpha_vgm_one, size(alpha_vgm_one), spval)
               n_vgm_patches (ipatch) = median (n_vgm_one, size(n_vgm_one), spval) 
               theta_s_patches (ipatch) = sum (theta_s_one * (area_one/sum(area_one)))

#ifdef SOILPAR_UPS_FIT
               np = size(theta_r_one)
               ipxstt = landpatch%ipxstt(ipatch)
               ipxend = landpatch%ipxend(ipatch)

               IF( np > 1 ) then
                  allocate ( ydatv  (ipxstt:ipxend,npointw) )
! the jacobian matrix required in Levenberg–Marquardt fitting method
                  allocate ( fjacv  (npointw,nv) )           ! calculated in SW_VG_dist
! the values of objective functions to be fitted 
                  allocate ( fvecv  (npointw)    )           ! calculated in SW_VG_dist 

! SW VG retentions at fine grids for each patch
                  do LL = ipxstt,ipxend
                     ydatv(LL,:) = theta_r_one(LL)+(theta_s_one(LL) - theta_r_one(LL)) &
                                 * (1+(alpha_vgm_one(LL)*xdat)**n_vgm_one(LL))**(1.0/n_vgm_one(LL)-1)
                  end do                  
                  
! Fitting the van Genuchten SW retention parameters             
                  ldfjac = npointw
                  xv(1) = theta_r_patches (ipatch)  
                  xv(2) = alpha_vgm_patches (ipatch) 
                  xv(3) = n_vgm_patches (ipatch)
                  maxfev = 100 * ( nv + 1 )
                  isiter = 1

                  call lmder ( SW_VG_dist, npointw, nv, xv, fvecv, fjacv, ldfjac, ftol, xtol, gtol, maxfev, &
                        diagv, mode, factor, nprint, info, nfev, njev, ipvtv, qtfv,&
                        xdat, npointw, ydatv, np, theta_s_patches(ipatch), isiter)

                  if ( xv(1) >= 0.0 .and. xv(1) <= theta_s_patches(ipatch) .and. xv(2) >= 1.0e-5 .and. xv(2) <= 1.0 .and. &
                       xv(3) >= 1.1 .and. xv(3) <= 10.0 .and. isiter == 1) then
                       theta_r_patches(ipatch)   = xv(1)
                       alpha_vgm_patches(ipatch) = xv(2)
                       n_vgm_patches(ipatch)     = xv(3)
                  end if

                  deallocate(ydatv)
                  deallocate(fjacv)
                  deallocate(fvecv)

               ENDIF 
#endif

            ELSE
               theta_r_patches (ipatch) = -1.0e36_r8
               alpha_vgm_patches (ipatch) = -1.0e36_r8
               n_vgm_patches (ipatch) = -1.0e36_r8
               theta_s_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(theta_r_patches(ipatch))) THEN
                write(*,*) "NAN appears in theta_r_patches."
                write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

            IF (isnan(alpha_vgm_patches(ipatch))) THEN
                write(*,*) "NAN appears in alpha_vgm_patches."
                write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

            IF (isnan(n_vgm_patches(ipatch))) THEN
                write(*,*) "NAN appears in n_vgm_patches."
                write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

            IF (isnan(theta_s_patches(ipatch))) THEN
                write(*,*) "NAN appears in theta_s_patches."
                write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('theta_r lev '//trim(c), theta_r_patches)
      CALL check_vector_data ('alpha VGM lev '//trim(c), alpha_vgm_patches)
      CALL check_vector_data ('n VGM lev '//trim(c), n_vgm_patches)
      CALL check_vector_data ('theta_s lev '//trim(c), theta_s_patches)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/theta_r_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'theta_r_l'//trim(c)//'_patches', 'patch', landpatch, theta_r_patches, 1)
#else
      SITE_soil_theta_r(nsl) = theta_r_patches(1)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/alpha_vgm_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'alpha_vgm_l'//trim(c)//'_patches', 'patch', landpatch, alpha_vgm_patches, 1)
#else
      SITE_soil_alpha_vgm(nsl) = alpha_vgm_patches(1)
#endif
      
#ifndef SinglePoint
      lndname = trim(landdir)//'/n_vgm_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'n_vgm_l'//trim(c)//'_patches', 'patch', landpatch, n_vgm_patches, 1)
#else
      SITE_soil_n_vgm(nsl) = n_vgm_patches(1)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/theta_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'theta_s_l'//trim(c)//'_patches', 'patch', landpatch, theta_s_patches, 1)
#else
      SITE_soil_theta_s(nsl) = theta_s_patches(1)
#endif
      
#endif

      ! (11) saturated water content [cm3/cm3]
      ! (12) matric potential at saturation (psi_s) [cm]
      ! (13) pore size distribution index [dimensionless] 

      IF (p_is_io) THEN

         CALL allocate_block_data (gland, theta_s_grid)
         lndname = trim(dir_rawdata)//'/soil/theta_s.nc'
         CALL ncio_read_block (lndname, 'theta_s_l'//trim(c), gland, theta_s_grid)

         CALL allocate_block_data (gland, psi_s_grid)
         lndname = trim(dir_rawdata)//'/soil/psi_s.nc'
         CALL ncio_read_block (lndname, 'psi_s_l'//trim(c), gland, psi_s_grid)

         CALL allocate_block_data (gland, lambda_grid)
         lndname = trim(dir_rawdata)//'/soil/lambda.nc'
         CALL ncio_read_block (lndname, 'lambda_l'//trim(c), gland, lambda_grid)

#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = theta_s_grid, &
                            data_r8_2d_in2 = psi_s_grid, data_r8_2d_in3 = lambda_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                            data_r8_2d_in1 = theta_s_grid, data_r8_2d_out1 = theta_s_one, &
                            data_r8_2d_in2 = psi_s_grid,   data_r8_2d_out2 = psi_s_one, &
                            data_r8_2d_in3 = lambda_grid,  data_r8_2d_out3 = lambda_one)
               theta_s_patches (ipatch) = sum (theta_s_one * (area_one/sum(area_one)))
               psi_s_patches (ipatch) = median (psi_s_one, size(psi_s_one), spval)
               lambda_patches (ipatch) = median (lambda_one, size(lambda_one), spval)

#ifdef SOILPAR_UPS_FIT
               np = size(psi_s_one)
               ipxstt = landpatch%ipxstt(ipatch)
               ipxend = landpatch%ipxend(ipatch)

               IF( np > 1 ) then
                  allocate ( ydatc  (ipxstt:ipxend,npointw) )
! the jacobian matrix required in Levenberg–Marquardt fitting method
                  allocate ( fjacc  (npointw,nc) )           ! calculated in SW_CB_dist
! the values of objective functions to be fitted 
                  allocate ( fvecc  (npointw)    )           ! calculated in SW_CB_dist 

! SW CB retentions at fine grids for each patch
                  do LL = ipxstt,ipxend
                     ydatc(LL,:) = (-1.0*xdat/psi_s_one(LL))**(-1.0*lambda_one(LL)) * theta_s_one(LL)
                  end do                  
                  
! Fitting the Campbell SW retention parameters             
                  ldfjac = npointw
                  xc(1) = psi_s_patches (ipatch)  
                  xc(2) = lambda_patches (ipatch) 
                  maxfev = 100 * ( nc + 1 )
                  isiter = 1

                  call lmder ( SW_CB_dist, npointw, nc, xc, fvecc, fjacc, ldfjac, ftol, xtol, gtol, maxfev, &
                        diagc, mode, factor, nprint, info, nfev, njev, ipvtc, qtfc,&
                        xdat, npointw, ydatc, np, theta_s_patches(ipatch), isiter)

                  if( xc(1) >= -300. .and. xc(1) < 0.0 .and. xc(2) > 0.0 .and. xc(2) <= 1.0 .and. isiter == 1)then
                        psi_s_patches (ipatch) = xc(1)
                        lambda_patches(ipatch) = xc(2)
                  end if

                  deallocate(ydatc)
                  deallocate(fjacc)
                  deallocate(fvecc)

               ENDIF 
#endif

            ELSE
               theta_s_patches (ipatch) = -1.0e36_r8
               psi_s_patches (ipatch) = -1.0e36_r8
               lambda_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(theta_s_patches(ipatch))) THEN
                write(*,*) "NAN appears in theta_s_patches."
                write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

            IF (isnan(psi_s_patches(ipatch))) THEN
                write(*,*) "NAN appears in psi_s_patches."
                write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

            IF (isnan(lambda_patches(ipatch))) THEN
                write(*,*) "NAN appears in lambda_patches."
                write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('theta_s lev '//trim(c), theta_s_patches)
      CALL check_vector_data ('psi_s lev '//trim(c), psi_s_patches)
      CALL check_vector_data ('lambda lev '//trim(c), lambda_patches)
#endif

#ifndef vanGenuchten_Mualem_SOIL_MODEL
#ifndef SinglePoint
      lndname = trim(landdir)//'/theta_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'theta_s_l'//trim(c)//'_patches', 'patch', landpatch, theta_s_patches, 1)
#else
      SITE_soil_theta_s(nsl) = theta_s_patches(1)
#endif
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/psi_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'psi_s_l'//trim(c)//'_patches', 'patch', landpatch, psi_s_patches, 1)
#else
      SITE_soil_psi_s(nsl) = psi_s_patches(1)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/lambda_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'lambda_l'//trim(c)//'_patches', 'patch', landpatch, lambda_patches, 1)
#else
      SITE_soil_lambda(nsl) = lambda_patches(1)
#endif

      ! (14) saturated hydraulic conductivity [cm/day]
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, k_s_grid)
         lndname = trim(dir_rawdata)//'/soil/k_s.nc'
         CALL ncio_read_block (lndname, 'k_s_l'//trim(c), gland, k_s_grid)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = k_s_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                  data_r8_2d_in1 = k_s_grid, data_r8_2d_out1 = k_s_one)
               k_s_patches (ipatch) = product(k_s_one**(area_one/sum(area_one)))
            ELSE
               k_s_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(k_s_patches(ipatch))) THEN
                write(*,*) "NAN appears in k_s_patches."
                write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('k_s lev '//trim(c), k_s_patches)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/k_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'k_s_l'//trim(c)//'_patches', 'patch', landpatch, k_s_patches, 1)
#else
      SITE_soil_k_s(nsl) = k_s_patches(1)
#endif

      ! (15) heat capacity of soil solids [J/(m3 K)]
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, csol_grid)
         lndname = trim(dir_rawdata)//'/soil/csol.nc'
         CALL ncio_read_block (lndname, 'csol_l'//trim(c), gland, csol_grid)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = csol_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                  data_r8_2d_in1 = csol_grid, data_r8_2d_out1 = csol_one)
               csol_patches (ipatch) = sum(csol_one*(area_one/sum(area_one)))
            ELSE
               csol_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(csol_patches(ipatch))) THEN
                write(*,*) "NAN appears in csol_patches."
                write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('csol lev '//trim(c), csol_patches)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/csol_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'csol_l'//trim(c)//'_patches', 'patch', landpatch, csol_patches, 1)
#else
      SITE_soil_csol(nsl) = csol_patches(1)
#endif

      ! (16) thermal conductivity of unfrozen saturated soil [W/m-K]
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, tksatu_grid)
         lndname = trim(dir_rawdata)//'/soil/tksatu.nc'
         CALL ncio_read_block (lndname, 'tksatu_l'//trim(c), gland, tksatu_grid)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = tksatu_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                  data_r8_2d_in1 = tksatu_grid, data_r8_2d_out1 = tksatu_one)
               tksatu_patches (ipatch) = product(tksatu_one**(area_one/sum(area_one)))
            ELSE
               tksatu_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(tksatu_patches(ipatch))) THEN
                write(*,*) "NAN appears in tksatu_patches."
                write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('tksatu lev '//trim(c), tksatu_patches)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/tksatu_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'tksatu_l'//trim(c)//'_patches', 'patch', landpatch, tksatu_patches, 1)
#else
      SITE_soil_tksatu(nsl) = tksatu_patches(1)
#endif

      ! (17) thermal conductivity of frozen saturated soil [W/m-K]
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, tksatf_grid)
         lndname = trim(dir_rawdata)//'/soil/tksatf.nc'
         CALL ncio_read_block (lndname, 'tksatf_l'//trim(c), gland, tksatf_grid)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = tksatf_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                  data_r8_2d_in1 = tksatf_grid, data_r8_2d_out1 = tksatf_one)
               tksatf_patches (ipatch) = product(tksatf_one**(area_one/sum(area_one)))
            ELSE
               tksatf_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(tksatf_patches(ipatch))) THEN
                write(*,*) "NAN appears in tksatf_patches."
                write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('tksatf lev '//trim(c), tksatf_patches)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/tksatf_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'tksatf_l'//trim(c)//'_patches', 'patch', landpatch, tksatf_patches, 1)
#else
      SITE_soil_tksatf(nsl) = tksatf_patches(1)
#endif

      ! (18) thermal conductivity for dry soil [W/(m-K)]
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, tkdry_grid)
         lndname = trim(dir_rawdata)//'/soil/tkdry.nc'
         CALL ncio_read_block (lndname, 'tkdry_l'//trim(c), gland, tkdry_grid)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = tkdry_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                  data_r8_2d_in1 = tkdry_grid, data_r8_2d_out1 = tkdry_one)
               tkdry_patches (ipatch) = product(tkdry_one**(area_one/sum(area_one)))
            ELSE
               tkdry_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(tkdry_patches(ipatch))) THEN
                write(*,*) "NAN appears in tkdry_patches."
                write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('tkdry lev '//trim(c), tkdry_patches)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/tkdry_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'tkdry_l'//trim(c)//'_patches', 'patch', landpatch, tkdry_patches, 1)
#else
      SITE_soil_tkdry(nsl) = tkdry_patches(1)
#endif

      ! (19) thermal conductivity of soil solids [W/m-K]
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, k_solids_grid)
         lndname = trim(dir_rawdata)//'/soil/k_solids.nc'
         CALL ncio_read_block (lndname, 'k_solids_l'//trim(c), gland, k_solids_grid)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = k_solids_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                  data_r8_2d_in1 = k_solids_grid, data_r8_2d_out1 = k_solids_one)
               k_solids_patches (ipatch) = product(k_solids_one**(area_one/sum(area_one)))
            ELSE
               k_solids_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(k_solids_patches(ipatch))) THEN
                write(*,*) "NAN appears in k_solids_patches."
                write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('k_solids lev '//trim(c), k_solids_patches)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/k_solids_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'k_solids_l'//trim(c)//'_patches', 'patch', landpatch, k_solids_patches, 1)
#else
      SITE_soil_k_solids(nsl) = k_solids_patches(1)
#endif

      ! (20) OM_density [kg/m3]
      IF (p_is_io) THEN
        
         CALL allocate_block_data (gland, OM_density_s_grid)
         lndname = trim(dir_rawdata)//'/soil/OM_density_s.nc'
         CALL ncio_read_block (lndname, 'OM_density_s_l'//trim(c), gland, OM_density_s_grid)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = OM_density_s_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                  data_r8_2d_in1 = OM_density_s_grid, data_r8_2d_out1 = OM_density_s_one)
               OM_density_s_patches (ipatch) = sum (OM_density_s_one * (area_one/sum(area_one)))
            ELSE
               OM_density_s_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(OM_density_s_patches(ipatch))) THEN
               write(*,*) "NAN appears in OM_density_s_patches."
               write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('OM_density_s lev '//trim(c), OM_density_s_patches)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/OM_density_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'OM_density_s_l'//trim(c)//'_patches', 'patch',& 
                              landpatch, OM_density_s_patches, 1)
#else
      SITE_soil_OM_density(nsl) = OM_density_s_patches(1)
#endif

      ! (21) bulk density of soil (GRAVELS + OM + Mineral Soils)
      IF (p_is_io) THEN
        
         CALL allocate_block_data (gland, BD_all_s_grid)
         lndname = trim(dir_rawdata)//'/soil/BD_all_s.nc'
         CALL ncio_read_block (lndname, 'BD_all_s_l'//trim(c), gland, BD_all_s_grid)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = BD_all_s_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
                  data_r8_2d_in1 = BD_all_s_grid, data_r8_2d_out1 = BD_all_s_one)
               BD_all_s_patches (ipatch) = sum (BD_all_s_one * (area_one/sum(area_one)))
            ELSE
               BD_all_s_patches (ipatch) = -1.0e36_r8
            ENDIF

            IF (isnan(BD_all_s_patches(ipatch))) THEN
               write(*,*) "NAN appears in BD_all_s_patches."
               write(*,*) landpatch%eindex(ipatch), landpatch%settyp(ipatch)
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('BD_all_s lev '//trim(c), BD_all_s_patches)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/BD_all_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'BD_all_s_l'//trim(c)//'_patches', 'patch',& 
                              landpatch, BD_all_s_patches, 1)
#else
      SITE_soil_BD_all(nsl) = BD_all_s_patches(1)
#endif


   ENDDO


   ! Deallocate the allocatable array
   ! --------------------------------

   IF (p_is_worker) THEN

      deallocate ( vf_quartz_mineral_s_patches )
      deallocate ( vf_gravels_s_patches )
      deallocate ( vf_om_s_patches )
      deallocate ( vf_sand_s_patches )
      deallocate ( wf_gravels_s_patches )
      deallocate ( wf_sand_s_patches )
      deallocate ( OM_density_s_patches )
      deallocate ( BD_all_s_patches )
      deallocate ( theta_s_patches )
      deallocate ( psi_s_patches   )
      deallocate ( lambda_patches  )
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      deallocate ( theta_r_patches  )
      deallocate ( alpha_vgm_patches)
      deallocate ( L_vgm_patches    )
      deallocate ( n_vgm_patches    )
#endif
      deallocate ( k_s_patches     )
      deallocate ( csol_patches    )
      deallocate ( tksatu_patches  )
      deallocate ( tkdry_patches   )
      deallocate ( tksatf_patches  )
      deallocate ( k_solids_patches)
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
      deallocate ( BA_alpha_patches)
      deallocate ( BA_beta_patches )
#endif

      IF (allocated(vf_quartz_mineral_s_one)) deallocate (vf_quartz_mineral_s_one)
      IF (allocated(vf_gravels_s_one))        deallocate (vf_gravels_s_one)
      IF (allocated(vf_om_s_one))             deallocate (vf_om_s_one)
      IF (allocated(vf_sand_s_one))           deallocate (vf_sand_s_one)
      IF (allocated(wf_gravels_s_one))        deallocate (wf_gravels_s_one)
      IF (allocated(wf_sand_s_one))           deallocate (wf_sand_s_one) 
      IF (allocated(OM_density_s_one))        deallocate (OM_density_s_one)
      IF (allocated(BD_all_s_one))            deallocate (BD_all_s_one)
      IF (allocated(theta_s_one))             deallocate (theta_s_one)
      IF (allocated(psi_s_one  ))             deallocate (psi_s_one  )
      IF (allocated(lambda_one ))     deallocate (lambda_one )
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      IF (allocated ( theta_r_one  )) deallocate ( theta_r_one  )
      IF (allocated ( alpha_vgm_one)) deallocate ( alpha_vgm_one)
      IF (allocated ( L_vgm_one    )) deallocate ( L_vgm_one    )
      IF (allocated ( n_vgm_one    )) deallocate ( n_vgm_one    )
#endif
      IF (allocated(k_s_one    ))     deallocate (k_s_one    )
      IF (allocated(csol_one   ))     deallocate (csol_one   )
      IF (allocated(tksatu_one ))     deallocate (tksatu_one )
      IF (allocated(tkdry_one  ))     deallocate (tkdry_one  )
      IF (allocated(tksatf_one ))     deallocate (tksatf_one)
      IF (allocated(k_solids_one))    deallocate (k_solids_one)
      IF (allocated(area_one))        deallocate (area_one)

   ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

END SUBROUTINE aggregation_soil_parameters

#ifdef SOILPAR_UPS_FIT

!SUBROUTINE Ke_Sr_dist ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydatb, nptf, phi, isiter, &
!                        vf_om_s, vf_sand_s, vf_gravels_s )

!!=================================================================
!! This is the subroutine for calculating the function/jacobian matrix 
!! of the distance between the fitted and prescribed Ke-Sr relationship 
!! for the Balland V. and P. A. Arp (2005) model.
!!
!! Nan Wei, 02/2020
!! ----------------------------------------------------------------

!      use precision
!      implicit none

!      integer m,n,ldfjac,iflag,i,nptf,isiter,npoint
!      real(r8) x(n),fjac(ldfjac,n),fvec(m),xdat(npoint),ydatb(nptf,npoint),phi
!      real(r8) vf_om_s, vf_sand_s, vf_gravels_s

!      if ( iflag == 0 ) then

!         print*,x

!      else if ( iflag == 1 ) then

!         if (x(2) <= 0.0) then
!             isiter = 0
!             return
!         end if

!         do i = 1, m
!             fvec(i) = sum((xdat(i)**(0.5*(1.0+vf_om_s-x(1)*vf_sand_s-vf_gravels_s)) * &
!                       ((1.0/(1.0+exp(-x(2)*xdat(i))))**3 - ((1.0-xdat(i))/2.0)**3)**(1.0-vf_om_s) - ydatb(:,i))**2)
!         end do

!      else if ( iflag == 2 ) then

!         if (x(2) <= 0.0) then
!             isiter = 0
!             return
!         end if

!         do i = 1, m
!             fjac(i,1) = sum(-2.0*(xdat(i)**(0.5*(1.0+vf_om_s-x(1)*vf_sand_s-vf_gravels_s)) * &
!                       ((1.0/(1.0+exp(-x(2)*xdat(i))))**3 - ((1.0-xdat(i))/2.0)**3)**(1.0-vf_om_s) - ydatb(:,i)) * &
!                       ((1.0/(1.0+exp(-x(2)*xdat(i))))**3 - ((1.0-xdat(i))/2.0)**3)**(1.0-vf_om_s) * &
!                       xdat(i)**(0.5*(1.0+vf_om_s-x(1)*vf_sand_s-vf_gravels_s)) * log(xdat(i)) * 0.5 * vf_sand_s)
!             fjac(i,2) = sum(2.0*(xdat(i)**(0.5*(1.0+vf_om_s-x(1)*vf_sand_s-vf_gravels_s)) * &
!                       ((1.0/(1.0+exp(-x(2)*xdat(i))))**3 - ((1.0-xdat(i))/2.0)**3)**(1.0-vf_om_s) - ydatb(:,i)) * &
!                       xdat(i)**(0.5*(1.0+vf_om_s-x(1)*vf_sand_s-vf_gravels_s)) * (1.0-vf_om_s) * &
!                       ((1.0/(1.0+exp(-x(2)*xdat(i))))**3 - ((1.0-xdat(i))/2.0)**3)**(-vf_om_s) * 3.0 * &
!                       (1.0/(1.0+exp(-x(2)*xdat(i))))**4 * exp(-x(2)*xdat(i)) * xdat(i))
!         end do

!      end if

!end subroutine Ke_Sr_dist

subroutine SW_CB_dist ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydatc, nptf, phi, isiter)

!=================================================================
! This is the subroutine for calculating the function/jacobian matrix 
! of the distance between the fitted and prescribed SW retention curves 
! for the Campbell model.
!
! Nan Wei, 06/2019
! ----------------------------------------------------------------

      use precision
      implicit none

      integer m,n,ldfjac,iflag,i,nptf,isiter,npoint
      real(r8) x(n),fjac(ldfjac,n),fvec(m),xdat(npoint),ydatc(nptf,npoint),phi

      if ( iflag == 0 ) then

         print*,x

      else if ( iflag == 1 ) then

         if (x(1) >= 0.0) then
             isiter = 0
             return
         end if

         do i = 1, m
            fvec(i) = sum(((-1.0*xdat(i)/x(1))**(-1.0*x(2)) * phi - ydatc(:,i))**2)
         end do

      else if ( iflag == 2 ) then

         if (x(1) >= 0.0) then
             isiter = 0
             return
         end if

         do i = 1, m
            fjac(i,1) = sum(2.0*((-1.0*xdat(i)/x(1))**(-1.0*x(2)) * phi - ydatc(:,i))*&
                        phi * x(2) * (-1.0*xdat(i)/x(1))**(-1.0*x(2)) / x(1))
            fjac(i,2) = sum(-2.0*((-1.0*xdat(i)/x(1))**(-1.0*x(2)) * phi - ydatc(:,i))*&
                        phi * (-1.0*xdat(i)/x(1))**(-1.0*x(2)) * log(-1.0*xdat(i)/x(1)))
         end do

      end if

end subroutine SW_CB_dist

subroutine SW_VG_dist ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydatv, nptf, phi, isiter )

!=================================================================
! This is the subroutine for calculating the function/jacobian matrix 
! of the distance between the fitted and prescribed SW retention curves 
! for the van Genuchten model.
!
! Nan Wei, 06/2019
! ----------------------------------------------------------------

      use precision
      implicit none

      integer m,n,ldfjac,iflag,i,nptf,isiter,npoint
      real(r8) x(n),fjac(ldfjac,n),fvec(m),xdat(npoint),ydatv(nptf,npoint),phi

      if ( iflag == 0 ) then

         print*,x

      else if ( iflag == 1 ) then

         if (x(2) <= 0.0 .or. x(3) <= 0.1) then
             isiter = 0
             return
         end if

         do i = 1, m
            fvec(i) = sum((x(1) + (phi - x(1))*(1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) - ydatv(:,i))**2)
         end do

      else if ( iflag == 2 ) then

         if (x(2) <= 0.0 .or. x(3) <= 0.1) then
             isiter = 0
             return
         end if

         do i = 1, m
            fjac(i,1) = sum(2*(x(1) + (phi - x(1))*(1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) - ydatv(:,i))*&
                        (1 - (1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1)))
            fjac(i,2) = sum(2*(x(1) + (phi - x(1))*(1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) - ydatv(:,i))*&
                        (phi - x(1)) * (1 - x(3)) * (1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-2) * x(2)**(x(3)-1) * xdat(i)**x(3))
            fjac(i,3) = sum(2*(x(1) + (phi - x(1))*(1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) - ydatv(:,i))*&
                        (phi - x(1)) * (1+(x(2)*xdat(i))**x(3))**(1.0/x(3)-1) *&
                        ((1.0-x(3))*(x(2)*xdat(i))**x(3)*log(x(2)*xdat(i))/(x(3)*(1+(x(2)*xdat(i))**x(3))) &
                        - log(1+(x(2)*xdat(i))**x(3))/x(3)**2))
         end do

      end if

end subroutine SW_VG_dist

#endif
!-----------------------------------------------------------------------
!EOP
