#include <define.h>

#ifdef SinglePoint
MODULE mod_single_srfdata
   
   USE precision, only: r8
   USE GlobalVars
   USE mod_namelist
   IMPLICIT NONE
   SAVE

   REAL(r8) :: SITE_htop = 1.
#if (defined PFT_CLASSIFICATION || defined PC_CLASSIFICATION)
   REAL(r8) :: SITE_htop_pfts (0:N_PFT-1) = 1. 
#endif
   
   INTEGER,  allocatable :: SITE_LAI_year (:)
   REAL(r8), allocatable :: SITE_LAI1 (:)
   REAL(r8), allocatable :: SITE_SAI1 (:)
   REAL(r8), allocatable :: SITE_LAI2 (:,:)
   REAL(r8), allocatable :: SITE_SAI2 (:,:)
#if (defined PFT_CLASSIFICATION || defined PC_CLASSIFICATION)
   REAL(r8), allocatable :: SITE_LAI_pfts2 (:,:)  
   REAL(r8), allocatable :: SITE_SAI_pfts2 (:,:)
   REAL(r8), allocatable :: SITE_LAI_pfts3 (:,:,:)  
   REAL(r8), allocatable :: SITE_SAI_pfts3 (:,:,:)
#endif

   REAL(r8) :: SITE_lakedepth = 1.

   REAL(r8), allocatable :: SITE_soil_s_v_alb (:)
   REAL(r8), allocatable :: SITE_soil_d_v_alb (:)
   REAL(r8), allocatable :: SITE_soil_s_n_alb (:)
   REAL(r8), allocatable :: SITE_soil_d_n_alb (:)

   REAL(r8), allocatable :: SITE_soil_vf_quartz_mineral (:)
   REAL(r8), allocatable :: SITE_soil_vf_gravels        (:) 
   REAL(r8), allocatable :: SITE_soil_vf_sand           (:) 
   REAL(r8), allocatable :: SITE_soil_vf_om             (:) 
   REAL(r8), allocatable :: SITE_soil_wf_gravels        (:) 
   REAL(r8), allocatable :: SITE_soil_wf_sand           (:) 
   REAL(r8), allocatable :: SITE_soil_theta_s           (:) 
   REAL(r8), allocatable :: SITE_soil_k_s               (:) 
   REAL(r8), allocatable :: SITE_soil_csol              (:) 
   REAL(r8), allocatable :: SITE_soil_tksatu            (:) 
   REAL(r8), allocatable :: SITE_soil_tksatf            (:) 
   REAL(r8), allocatable :: SITE_soil_tkdry             (:) 
   REAL(r8), allocatable :: SITE_soil_k_solids          (:) 
#ifdef Campbell_SOIL_MODEL
   REAL(r8), allocatable :: SITE_soil_psi_s             (:) 
   REAL(r8), allocatable :: SITE_soil_lambda            (:) 
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   REAL(r8), allocatable :: SITE_soil_theta_r           (:) 
   REAL(r8), allocatable :: SITE_soil_alpha_vgm         (:) 
   REAL(r8), allocatable :: SITE_soil_L_vgm             (:) 
   REAL(r8), allocatable :: SITE_soil_n_vgm             (:) 
#endif
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
   REAL(r8), allocatable :: SITE_soil_BA_alpha          (:) 
   REAL(r8), allocatable :: SITE_soil_BA_beta           (:) 
#endif

#ifdef USE_DEPTH_TO_BEDROCK
   REAL(r8) :: SITE_dbedrock = 1
#endif

CONTAINS

   ! -----
   SUBROUTINE read_surface_data_single (fsrfdata)

      USE ncio_serial
      USE mod_namelist
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: fsrfdata

      ! Local Variables
      INTEGER :: iyear, itime

      IF (USE_SITE_htop) THEN
#if (defined PFT_CLASSIFICATION || defined PC_CLASSIFICATION)
         CALL ncio_read_serial (fsrfdata, 'canopy_height_pfts', SITE_htop_pfts)
         SITE_htop(:) = sum(SITE_pct_pfts * SITE_htop_pfts(:))
#else
         CALL ncio_read_serial (fsrfdata, 'canopy_height', SITE_htop)
#endif
      ENDIF

      IF (USE_SITE_LAI) THEN
#if (defined PFT_CLASSIFICATION || defined PC_CLASSIFICATION)
         IF (DEF_LAI_CLIM) THEN
            CALL ncio_read_serial (fsrfdata, 'lai_pfts', SITE_LAI_pfts2)
            CALL ncio_read_serial (fsrfdata, 'sai_pfts', SITE_SAI_pfts2)

            allocate (SITE_LAI1 (size(SITE_LAI_pfts2,2))
            allocate (SITE_SAI1 (size(SITE_SAI_pfts2,2))
            DO itime = 1, 12
               SITE_LAI1(itime) = sum(SITE_pct_pfts * SITE_LAI_pfts2(:,itime))
               SITE_SAI1(itime) = sum(SITE_pct_pfts * SITE_SAI_pfts2(:,itime))
            ENDDO
         ELSE
            CALL ncio_read_serial (fsrfdata, 'lai_year', SITE_LAI_year)
            CALL ncio_read_serial (fsrfdata, 'lai_pfts', SITE_LAI_pfts3)
            CALL ncio_read_serial (fsrfdata, 'sai_pfts', SITE_SAI_pfts3)

            allocate (SITE_LAI2 (size(SITE_LAI_pfts3,2),size(SITE_LAI_pfts3,3))
            allocate (SITE_SAI2 (size(SITE_SAI_pfts3,2),size(SITE_SAI_pfts3,3))
            DO iyear = lbound(SITE_LAI_pfts3,3), ubound(SITE_LAI_pfts3,3)
               DO itime = 1, 46
                  SITE_LAI(:,itime,iyear) = sum(SITE_pct_pfts * SITE_LAI_pfts3(:,itime,iyear))
                  SITE_SAI(:,itime,iyear) = sum(SITE_pct_pfts * SITE_SAI_pfts3(:,itime,iyear))
               ENDDO
            ENDDO
         ENDIF
#else
         IF (DEF_LAI_CLIM) THEN
            CALL ncio_read_serial (fsrfdata, 'LAI_clim', SITE_LAI1)
         ELSE
            CALL ncio_read_serial (fsrfdata, 'LAI_year',  SITE_LAI_year)
            CALL ncio_read_serial (fsrfdata, 'LAI_modis', SITE_LAI2)
         ENDIF
#endif
      ENDIF

      IF (USE_SITE_lakedepth) THEN
         CALL ncio_read_serial (fsrfdata, 'lakedepth', SITE_lakedepth)
      ENDIF
      
      IF (USE_SITE_soilreflectance) THEN
         CALL ncio_read_serial (fsrfdata, 'soil_s_v_alb', SITE_soil_s_v_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_d_v_alb', SITE_soil_d_v_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_s_n_alb', SITE_soil_s_n_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_d_n_alb', SITE_soil_d_n_alb)
      ENDIF

      IF (USE_SITE_soilparameters) THEN
         CALL ncio_read_serial (fsrfdata, 'vf_quartz_mineral', SITE_soil_vf_quartz_mineral)
         CALL ncio_read_serial (fsrfdata, 'vf_gravels       ', SITE_soil_vf_gravels       ) 
         CALL ncio_read_serial (fsrfdata, 'vf_sand          ', SITE_soil_vf_sand          ) 
         CALL ncio_read_serial (fsrfdata, 'vf_om            ', SITE_soil_vf_om            ) 
         CALL ncio_read_serial (fsrfdata, 'wf_gravels       ', SITE_soil_wf_gravels       ) 
         CALL ncio_read_serial (fsrfdata, 'wf_sand          ', SITE_soil_wf_sand          ) 
         CALL ncio_read_serial (fsrfdata, 'theta_s          ', SITE_soil_theta_s          ) 
         CALL ncio_read_serial (fsrfdata, 'k_s              ', SITE_soil_k_s              ) 
         CALL ncio_read_serial (fsrfdata, 'csol             ', SITE_soil_csol             ) 
         CALL ncio_read_serial (fsrfdata, 'tksatu           ', SITE_soil_tksatu           ) 
         CALL ncio_read_serial (fsrfdata, 'tksatf           ', SITE_soil_tksatf           ) 
         CALL ncio_read_serial (fsrfdata, 'tkdry            ', SITE_soil_tkdry            ) 
         CALL ncio_read_serial (fsrfdata, 'k_solids         ', SITE_soil_k_solids         ) 
#ifdef Campbell_SOIL_MODEL
         CALL ncio_read_serial (fsrfdata, 'psi_s            ', SITE_soil_psi_s            ) 
         CALL ncio_read_serial (fsrfdata, 'lambda           ', SITE_soil_lambda           ) 
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         CALL ncio_read_serial (fsrfdata, 'theta_r          ', SITE_soil_theta_r          ) 
         CALL ncio_read_serial (fsrfdata, 'alpha_vgm        ', SITE_soil_alpha_vgm        ) 
         CALL ncio_read_serial (fsrfdata, 'L_vgm            ', SITE_soil_L_vgm            ) 
         CALL ncio_read_serial (fsrfdata, 'n_vgm            ', SITE_soil_n_vgm            ) 
#endif
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
         CALL ncio_read_serial (fsrfdata, 'BA_alpha         ', SITE_soil_BA_alpha         ) 
         CALL ncio_read_serial (fsrfdata, 'BA_beta          ', SITE_soil_BA_beta          ) 
#endif
      ENDIF

#ifdef USE_DEPTH_TO_BEDROCK
      IF (USE_SITE_dbedrock) THEN
         CALL ncio_read_serial (fsrfdata, 'depth_to_bedrock', SITE_dbedrock)
      ENDIF
#endif

   END SUBROUTINE read_surface_data_single

END MODULE mod_single_srfdata
#endif
