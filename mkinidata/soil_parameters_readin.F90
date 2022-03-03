#include <define.h>

SUBROUTINE soil_parameters_readin (dir_landdata)
   ! ======================================================================
   ! Read in soil parameters in (patches,lon_points,lat_points) and
   ! => 1d vector [numpatch]
   !
   ! Created by Yongjiu Dai, 03/2014
   !             
   ! ======================================================================
   use precision
   USE GlobalVars, only : nl_soil
   use spmd_task
   use ncio_vector
   use mod_landpatch
   use MOD_TimeInvariants
   use mod_colm_debug

   IMPLICIT NONE

   ! ----------------------------------------------------------------------
   
   character(LEN=*), INTENT(in) :: dir_landdata

   ! Local Variables
   real(r8), allocatable :: soil_theta_s_l (:)  ! saturated water content (cm3/cm3)
   real(r8), allocatable :: soil_psi_s_l   (:)  ! matric potential at saturation (cm)
#ifdef Campbell_SOIL_MODEL
   real(r8), allocatable :: soil_lambda_l  (:)  ! pore size distribution index (dimensionless)
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   real(r8), allocatable :: soil_theta_r_l   (:)  ! residual water content (cm3/cm3)
   real(r8), allocatable :: soil_alpha_vgm_l (:)  
   real(r8), allocatable :: soil_L_vgm_l     (:)  
   real(r8), allocatable :: soil_n_vgm_l     (:)  
#endif
   real(r8), allocatable :: soil_k_s_l     (:)  ! saturated hydraulic conductivity (cm/day)
   real(r8), allocatable :: soil_csol_l    (:)  ! heat capacity of soil solids [J/(m3 K)]
   real(r8), allocatable :: soil_tksatu_l  (:)  ! thermal conductivity of saturated unforzen soil [W/m-K]
   real(r8), allocatable :: soil_tkdry_l   (:)  ! thermal conductivity for dry soil  [W/(m-K)]

   integer  :: ipatch, m, nsl  ! indices

   character(len=256) :: c
   character(len=256) :: lndname

   ! ...............................................................

   if (p_is_worker) then

      if (numpatch > 0) then

         allocate ( soil_theta_s_l (numpatch) )
         allocate ( soil_psi_s_l   (numpatch) )
#ifdef Campbell_SOIL_MODEL
         allocate ( soil_lambda_l  (numpatch) )
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         allocate ( soil_theta_r_l   (numpatch) )  
         allocate ( soil_alpha_vgm_l (numpatch) )  
         allocate ( soil_L_vgm_l     (numpatch) )  
         allocate ( soil_n_vgm_l     (numpatch) )  
#endif
         allocate ( soil_k_s_l     (numpatch) )
         allocate ( soil_csol_l    (numpatch) )
         allocate ( soil_tksatu_l  (numpatch) )
         allocate ( soil_tkdry_l   (numpatch) )

      end if

   end if

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   DO nsl = 1, 8
      
      write(c,'(i1)') nsl 

      ! (1) read in the saturated water content [cm3/cm3]
      lndname = trim(dir_landdata)//'/theta_s_l'//trim(c)//'_patches.nc'
      call ncio_read_vector (lndname, 'theta_s_l'//trim(c)//'_patches', landpatch, soil_theta_s_l)

      ! (2) read in the matric potential at saturation [cm]
      lndname = trim(dir_landdata)//'/psi_s_l'//trim(c)//'_patches.nc'
      call ncio_read_vector (lndname, 'psi_s_l'//trim(c)//'_patches', landpatch, soil_psi_s_l)

#ifdef Campbell_SOIL_MODEL
      ! (3) read in the pore size distribution index [dimensionless]
      lndname = trim(dir_landdata)//'/lambda_l'//trim(c)//'_patches.nc'
      call ncio_read_vector (lndname, 'lambda_l'//trim(c)//'_patches', landpatch, soil_lambda_l)
#endif

#ifdef vanGenuchten_Mualem_SOIL_MODEL
      ! (3-1) read in saturated water content [cm3/cm3]
      lndname = trim(dir_landdata)//'/theta_r_l'//trim(c)//'_patches.nc'
      call ncio_read_vector (lndname, 'theta_r_l'//trim(c)//'_patches', landpatch, soil_theta_r_l)
      
      ! (3-2) read in alpha in VGM model
      lndname = trim(dir_landdata)//'/alpha_vgm_l'//trim(c)//'_patches.nc'
      call ncio_read_vector (lndname, 'alpha_vgm_l'//trim(c)//'_patches', landpatch, soil_alpha_vgm_l)
      
      ! (3-3) read in L in VGM model
      lndname = trim(dir_landdata)//'/L_vgm_l'//trim(c)//'_patches.nc'
      call ncio_read_vector (lndname, 'L_vgm_l'//trim(c)//'_patches', landpatch, soil_L_vgm_l)
      
      ! (3-4) read in n in VGM model
      lndname = trim(dir_landdata)//'/n_vgm_l'//trim(c)//'_patches.nc'
      call ncio_read_vector (lndname, 'n_vgm_l'//trim(c)//'_patches', landpatch, soil_n_vgm_l)
#endif


      ! (4) read in the saturated hydraulic conductivity [cm/day]
      lndname = trim(dir_landdata)//'/k_s_l'//trim(c)//'_patches.nc'

      call ncio_read_vector (lndname, 'k_s_l'//trim(c)//'_patches', landpatch, soil_k_s_l)

      ! (5) read in the heat capacity of soil solids [J/(m3 K)]
      lndname = trim(dir_landdata)//'/csol_l'//trim(c)//'_patches.nc'

      call ncio_read_vector (lndname, 'csol_l'//trim(c)//'_patches', landpatch, soil_csol_l)

      ! (6) read in the thermal conductivity of saturated soil [W/m-K]
      lndname = trim(dir_landdata)//'/tksatu_l'//trim(c)//'_patches.nc'

      call ncio_read_vector (lndname, 'tksatu_l'//trim(c)//'_patches', landpatch, soil_tksatu_l)

      ! (7) read in the thermal conductivity for dry soil [W/(m-K)]
      lndname = trim(dir_landdata)//'/tkdry_l'//trim(c)//'_patches.nc'

      call ncio_read_vector (lndname, 'tkdry_l'//trim(c)//'_patches', landpatch, soil_tkdry_l)

      if (p_is_worker) then

         do ipatch = 1, numpatch
            m = landpatch%ltyp(ipatch)
            if( m == 0 )then     ! ocean
               porsl (nsl,ipatch) = -1.e36
               psi0  (nsl,ipatch) = -1.e36
#ifdef Campbell_SOIL_MODEL
               bsw   (nsl,ipatch) = -1.e36
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
               theta_r  (nsl,ipatch) = -1.e36
               alpha_vgm(nsl,ipatch) = -1.e36
               L_vgm    (nsl,ipatch) = -1.e36
               n_vgm    (nsl,ipatch) = -1.e36
#endif
               hksati(nsl,ipatch) = -1.e36
               csol  (nsl,ipatch) = -1.e36
               dksatu(nsl,ipatch) = -1.e36
               dkdry (nsl,ipatch) = -1.e36
            else                 ! non ocean
               porsl (nsl,ipatch) =    soil_theta_s_l  (ipatch)               ! cm/cm
               psi0  (nsl,ipatch) =    soil_psi_s_l    (ipatch) * 10.         ! cm -> mm
#ifdef Campbell_SOIL_MODEL
               bsw   (nsl,ipatch) = 1./soil_lambda_l   (ipatch)               ! dimensionless
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
               theta_r  (nsl,ipatch) = soil_theta_r_l  (ipatch)
               alpha_vgm(nsl,ipatch) = soil_alpha_vgm_l(ipatch)
               L_vgm    (nsl,ipatch) = soil_L_vgm_l    (ipatch)
               n_vgm    (nsl,ipatch) = soil_n_vgm_l    (ipatch)
#endif
               hksati(nsl,ipatch) =    soil_k_s_l      (ipatch) * 10./86400.  ! cm/day -> mm/s
               csol  (nsl,ipatch) =    soil_csol_l     (ipatch)               ! J/(m2 K)
               dksatu(nsl,ipatch) =    soil_tksatu_l   (ipatch)               ! W/(m K)
               dkdry (nsl,ipatch) =    soil_tkdry_l    (ipatch)               ! W/(m K)
            endif
         end do

      end if

   ENDDO
         
   if (p_is_worker) then

      if (numpatch > 0) then
         deallocate ( soil_theta_s_l )
         deallocate ( soil_psi_s_l   )
#ifdef Campbell_SOIL_MODEL
         deallocate ( soil_lambda_l  )
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         deallocate ( soil_theta_r_l   )
         deallocate ( soil_alpha_vgm_l )
         deallocate ( soil_L_vgm_l     )
         deallocate ( soil_n_vgm_l     )
#endif
         deallocate ( soil_k_s_l     )
         deallocate ( soil_csol_l    )
         deallocate ( soil_tksatu_l  )
         deallocate ( soil_tkdry_l   )
      end if

   end if

   ! The parameters of the top NINTH soil layers were given by datasets
   ! [0-0.045 (LAYER 1-2), 0.045-0.091, 0.091-0.166, 0.166-0.289, 
   !  0.289-0.493, 0.493-0.829, 0.829-1.383 and 1.383-2.296 m].
   ! The NINTH layer's soil parameters will assigned to the bottom soil layer (2.296 - 3.8019m).

   if (p_is_worker) then

      do nsl = nl_soil, 9, -1
         porsl (nsl,:) = porsl (8,:)
         psi0  (nsl,:) = psi0  (8,:)
#ifdef Campbell_SOIL_MODEL
         bsw   (nsl,:) = bsw   (8,:)
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         theta_r  (nsl,:) = theta_r   (8,:)
         alpha_vgm(nsl,:) = alpha_vgm (8,:)
         L_vgm    (nsl,:) = L_vgm     (8,:)
         n_vgm    (nsl,:) = n_vgm     (8,:)
#endif
         hksati(nsl,:) = hksati(8,:)
         csol  (nsl,:) = csol  (8,:)
         dksatu(nsl,:) = dksatu(8,:)
         dkdry (nsl,:) = dkdry (8,:)
      end do

      do nsl = 8, 3, -1
         porsl (nsl,:) = porsl (nsl-1,:)
         psi0  (nsl,:) = psi0  (nsl-1,:)
#ifdef Campbell_SOIL_MODEL
         bsw   (nsl,:) = bsw   (nsl-1,:)
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         theta_r  (nsl,:) = theta_r   (nsl-1,:)
         alpha_vgm(nsl,:) = alpha_vgm (nsl-1,:)
         L_vgm    (nsl,:) = L_vgm     (nsl-1,:)
         n_vgm    (nsl,:) = n_vgm     (nsl-1,:)
#endif
         hksati(nsl,:) = hksati(nsl-1,:)
         csol  (nsl,:) = csol  (nsl-1,:)
         dksatu(nsl,:) = dksatu(nsl-1,:)
         dkdry (nsl,:) = dkdry (nsl-1,:)
      enddo

      porsl (2,:) = porsl (1,:)
      psi0  (2,:) = psi0  (1,:)
#ifdef Campbell_SOIL_MODEL
      bsw   (2,:) = bsw   (1,:)
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      theta_r  (2,:) = theta_r   (1,:)
      alpha_vgm(2,:) = alpha_vgm (1,:)
      L_vgm    (2,:) = L_vgm     (1,:)
      n_vgm    (2,:) = n_vgm     (1,:)
#endif
      hksati(2,:) = hksati(1,:)
      csol  (2,:) = csol  (1,:)
      dksatu(2,:) = dksatu(1,:)
      dkdry (2,:) = dkdry (1,:)

   end if

   ! Soil reflectance of broadband of visible(_v) and near-infrared(_n) of the sarurated(_s) and dry(_d) soil
#if(defined SOIL_REFL_GUESSED)
   if (p_is_worker) then
      do ipatch = 1, numpatch
         m = landpatch%ltyp(ipatch)
         CALL soil_color_refl(m,soil_s_v_alb(ipatch),soil_d_v_alb(ipatch),&
            soil_s_n_alb(ipatch),soil_d_n_alb(ipatch))
      enddo
   end if
#elif(defined SOIL_REFL_READ)

   ! (1) Read in the albedo of visible of the saturated soil
   lndname = trim(dir_landdata)//'/soil_s_v_alb_patches.nc'

   call ncio_read_vector (lndname, 'soil_s_v_alb', landpatch, soil_s_v_alb)

   ! (2) Read in the albedo of visible of the dry soil
   lndname = trim(dir_landdata)//'/soil_d_v_alb_patches.nc'

   call ncio_read_vector (lndname, 'soil_d_v_alb', landpatch, soil_d_v_alb)

   ! (3) Read in the albedo of near infrared of the saturated soil
   lndname = trim(dir_landdata)//'/soil_s_n_alb_patches.nc'

   call ncio_read_vector (lndname, 'soil_s_n_alb', landpatch, soil_s_n_alb)

   ! (4) Read in the albedo of near infrared of the dry soil
   lndname = trim(dir_landdata)//'/soil_d_n_alb_patches.nc'

   call ncio_read_vector (lndname, 'soil_d_n_alb', landpatch, soil_d_n_alb)

#endif

END SUBROUTINE soil_parameters_readin
! --------------------------------------------------
! EOP
