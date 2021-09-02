#include <define.h>

SUBROUTINE soil_parameters_readin (lon_points,lat_points,dir_model_landdata)
! ======================================================================
! Read in soil parameters in (patches,lon_points,lat_points) and
! => 1d vector [numpatch]
!
! Created by Yongjiu Dai, 03/2014
!             
! ======================================================================
   use precision
   use MOD_TimeInvariants
   USE GlobalVars

   IMPLICIT NONE

! ----------------------------------------------------------------------
   character(LEN=256), INTENT(in) :: dir_model_landdata
   integer, INTENT(in) :: lon_points ! number of longitude points on model grid
   integer, INTENT(in) :: lat_points ! number of latitude points on model grid

! ------------------------ local variables -----------------------------
  real(r8), allocatable :: soil_theta_s_l (:,:,:)  ! saturated water content (cm3/cm3)
  real(r8), allocatable :: soil_psi_s_l   (:,:,:)  ! matric potential at saturation (cm)
  real(r8), allocatable :: soil_lambda_l  (:,:,:)  ! pore size distribution index (dimensionless)
  real(r8), allocatable :: soil_k_s_l     (:,:,:)  ! saturated hydraulic conductivity (cm/day)
  real(r8), allocatable :: soil_csol_l    (:,:,:)  ! heat capacity of soil solids [J/(m3 K)]
  real(r8), allocatable :: soil_tksatu_l  (:,:,:)  ! thermal conductivity of saturated unforzen soil [W/m-K]
  real(r8), allocatable :: soil_tkdry_l   (:,:,:)  ! thermal conductivity for dry soil  [W/(m-K)]
#if(defined SOIL_REFL_READ)
  real(r8), allocatable :: s_v_alb        (:,:,:)  ! saturated visible soil reflectance
  real(r8), allocatable :: d_v_alb        (:,:,:)  ! dry visible soil reflectance
  real(r8), allocatable :: s_n_alb        (:,:,:)  ! saturated near infrared soil reflectance
  real(r8), allocatable :: d_n_alb        (:,:,:)  ! dry near infrared soil reflectance
#endif

  real(r8) :: a                        !
  integer  :: i,j,k,l,m,t,npatch,np,nsl  ! indices

  character(len=256) :: c
  CHARACTER(len=256) :: suffix
  character(len=256) :: lndname
  integer iunit
  integer MODEL_SOIL_LAYER

#ifdef USGS_CLASSIFICATION
  suffix        = ''
#else
  suffix        = '.igbp'
#endif

! ...............................................................

      iunit = 100
      DO nsl = 1, 8
         MODEL_SOIL_LAYER = nsl 
         write(c,'(i1)') MODEL_SOIL_LAYER

         allocate ( soil_theta_s_l (0:N_land_classification,1:lon_points,1:lat_points) )
         allocate ( soil_psi_s_l   (0:N_land_classification,1:lon_points,1:lat_points) )
         allocate ( soil_lambda_l  (0:N_land_classification,1:lon_points,1:lat_points) )
         allocate ( soil_k_s_l     (0:N_land_classification,1:lon_points,1:lat_points) )
         allocate ( soil_csol_l    (0:N_land_classification,1:lon_points,1:lat_points) )
         allocate ( soil_tksatu_l  (0:N_land_classification,1:lon_points,1:lat_points) )
         allocate ( soil_tkdry_l   (0:N_land_classification,1:lon_points,1:lat_points) )

! (1) read in the saturated water content [cm3/cm3]
         lndname = trim(dir_model_landdata)//'model_theta_s_l'//trim(c)//trim(suffix)//'.bin'
         print*,trim(lndname)
         OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
         READ(iunit,err=100) soil_theta_s_l
         close(iunit)

! (2) read in the matric potential at saturation [cm]
         lndname = trim(dir_model_landdata)//'model_psi_s_l'//trim(c)//trim(suffix)//'.bin'
         print*,trim(lndname)
         OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
         READ(iunit,err=100) soil_psi_s_l
         close(iunit)

! (3) read in the pore size distribution index [dimensionless]
         lndname = trim(dir_model_landdata)//'model_lambda_l'//trim(c)//trim(suffix)//'.bin'
         print*,trim(lndname)
         OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
         READ(iunit,err=100) soil_lambda_l
         close(iunit)

! (4) read in the saturated hydraulic conductivity [cm/day]
         lndname = trim(dir_model_landdata)//'model_k_s_l'//trim(c)//trim(suffix)//'.bin'
         print*,trim(lndname)
         OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
         READ(iunit,err=100) soil_k_s_l
         close(iunit)

! (5) read in the heat capacity of soil solids [J/(m3 K)]
         lndname = trim(dir_model_landdata)//'model_csol_l'//trim(c)//trim(suffix)//'.bin'
         print*,trim(lndname)
         OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
         READ(iunit,err=100) soil_csol_l
         close(iunit)

! (6) read in the thermal conductivity of saturated soil [W/m-K]
         lndname = trim(dir_model_landdata)//'model_tksatu_l'//trim(c)//trim(suffix)//'.bin'
         print*,trim(lndname)
         OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
         READ(iunit,err=100) soil_tksatu_l
         close(iunit)

! (7) read in the thermal conductivity for dry soil [W/(m-K)]
         lndname = trim(dir_model_landdata)//'model_tkdry_l'//trim(c)//trim(suffix)//'.bin'
         print*,trim(lndname)
         OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
         READ(iunit,err=100) soil_tkdry_l
         close(iunit)


         do npatch = 1, numpatch
            i = patch2lon(npatch)
            j = patch2lat(npatch)
! yuan, 12/28/2019: Bug, add patchtype
            t = patchtype(npatch)
            m = patchclass(npatch)
! for PFT structure, all the patches share the same soil data, storted at
! position 0 of array (0:N_land_classification,:,:)
#ifdef PFT_CLASSIFICATION
            IF( t == 0 ) THEN     
                porsl (nsl,npatch) =    soil_theta_s_l  (0,i,j)               ! cm/cm
                psi0  (nsl,npatch) =    soil_psi_s_l    (0,i,j) * 10.         ! cm -> mm
                bsw   (nsl,npatch) = 1./soil_lambda_l   (0,i,j)               ! dimensionless
                hksati(nsl,npatch) =    soil_k_s_l      (0,i,j) * 10./86400.  ! cm/day -> mm/s
                csol  (nsl,npatch) =    soil_csol_l     (0,i,j)               ! J/(m2 K)
                dksatu(nsl,npatch) =    soil_tksatu_l   (0,i,j)               ! W/(m K)
                dkdry (nsl,npatch) =    soil_tkdry_l    (0,i,j)               ! W/(m K)
            ELSE
                porsl (nsl,npatch) =    soil_theta_s_l  (m,i,j)               ! cm/cm
                psi0  (nsl,npatch) =    soil_psi_s_l    (m,i,j) * 10.         ! cm -> mm
                bsw   (nsl,npatch) = 1./soil_lambda_l   (m,i,j)               ! dimensionless
                hksati(nsl,npatch) =    soil_k_s_l      (m,i,j) * 10./86400.  ! cm/day -> mm/s
                csol  (nsl,npatch) =    soil_csol_l     (m,i,j)               ! J/(m2 K)
                dksatu(nsl,npatch) =    soil_tksatu_l   (m,i,j)               ! W/(m K)
                dkdry (nsl,npatch) =    soil_tkdry_l    (m,i,j)               ! W/(m K)
            END IF
#else
            if( m == 0 )then     
                ! ocean
                porsl (nsl,npatch) = -1.e36
                psi0  (nsl,npatch) = -1.e36
                bsw   (nsl,npatch) = -1.e36
                hksati(nsl,npatch) = -1.e36
                csol  (nsl,npatch) = -1.e36
                dksatu(nsl,npatch) = -1.e36
                dkdry (nsl,npatch) = -1.e36
            else
                porsl (nsl,npatch) =    soil_theta_s_l  (m,i,j)               ! cm/cm
                psi0  (nsl,npatch) =    soil_psi_s_l    (m,i,j) * 10.         ! cm -> mm
                bsw   (nsl,npatch) = 1./soil_lambda_l   (m,i,j)               ! dimensionless
                hksati(nsl,npatch) =    soil_k_s_l      (m,i,j) * 10./86400.  ! cm/day -> mm/s
                csol  (nsl,npatch) =    soil_csol_l     (m,i,j)               ! J/(m2 K)
                dksatu(nsl,npatch) =    soil_tksatu_l   (m,i,j)               ! W/(m K)
                dkdry (nsl,npatch) =    soil_tkdry_l    (m,i,j)               ! W/(m K)
            endif
#endif
         end do

         deallocate ( soil_theta_s_l )
         deallocate ( soil_psi_s_l   )
         deallocate ( soil_lambda_l  )
         deallocate ( soil_k_s_l     )
         deallocate ( soil_csol_l    )
         deallocate ( soil_tksatu_l  )
         deallocate ( soil_tkdry_l   )

      ENDDO

      ! The parameters of the top NINTH soil layers were given by datasets
      ! [0-0.045 (LAYER 1-2), 0.045-0.091, 0.091-0.166, 0.166-0.289, 
      !  0.289-0.493, 0.493-0.829, 0.829-1.383 and 1.383-2.296 m].
      ! The NINTH layer's soil parameters will assigned to the bottom soil layer (2.296 - 3.8019m).

      porsl (10,:) = porsl (8,:)
      psi0  (10,:) = psi0  (8,:)
      bsw   (10,:) = bsw   (8,:)
      hksati(10,:) = hksati(8,:)
      csol  (10,:) = csol  (8,:)
      dksatu(10,:) = dksatu(8,:)
      dkdry (10,:) = dkdry (8,:)

      porsl ( 9,:) = porsl (8,:)
      psi0  ( 9,:) = psi0  (8,:)
      bsw   ( 9,:) = bsw   (8,:)
      hksati( 9,:) = hksati(8,:)
      csol  ( 9,:) = csol  (8,:)
      dksatu( 9,:) = dksatu(8,:)
      dkdry ( 9,:) = dkdry (8,:)

      do nsl = 8, 3, -1
         porsl (nsl,:) = porsl (nsl-1,:)
         psi0  (nsl,:) = psi0  (nsl-1,:)
         bsw   (nsl,:) = bsw   (nsl-1,:)
         hksati(nsl,:) = hksati(nsl-1,:)
         csol  (nsl,:) = csol  (nsl-1,:)
         dksatu(nsl,:) = dksatu(nsl-1,:)
         dkdry (nsl,:) = dkdry (nsl-1,:)
      enddo

      porsl (2,:) = porsl (1,:)
      psi0  (2,:) = psi0  (1,:)
      bsw   (2,:) = bsw   (1,:)
      hksati(2,:) = hksati(1,:)
      csol  (2,:) = csol  (1,:)
      dksatu(2,:) = dksatu(1,:)
      dkdry (2,:) = dkdry (1,:)

#if(defined CLMDEBUG)
      print*,'porsl  =', minval(porsl , mask = porsl  .gt. -1.0e30), maxval(porsl , mask = porsl  .gt. -1.0e30)
      print*,'psi0   =', minval(psi0  , mask = psi0   .gt. -1.0e30), maxval(psi0  , mask = psi0   .gt. -1.0e30)
      print*,'bsw    =', minval(bsw   , mask = bsw    .gt.  0.0   ), maxval(bsw   , mask = bsw    .gt.  0.0   )
      print*,'hksati =', minval(hksati, mask = hksati .gt. -1.0e30), maxval(hksati, mask = hksati .gt. -1.0e30)
      print*,'csol   =', minval(csol  , mask = csol   .gt. -1.0e30), maxval(csol  , mask = csol   .gt. -1.0e30)
      print*,'dksatu =', minval(dksatu, mask = dksatu .gt. -1.0e30), maxval(dksatu, mask = dksatu .gt. -1.0e30)
      print*,'dkdry  =', minval(dkdry , mask = dkdry  .gt. -1.0e30), maxval(dkdry , mask = dkdry  .gt. -1.0e30)
#endif

! Soil reflectance of broadband of visible(_v) and near-infrared(_n) of the sarurated(_s) and dry(_d) soil
#if(defined SOIL_REFL_GUESSED)
      do i = 1, numpatch
         CALL soil_color_refl(mxy_patch(i),soil_s_v_alb(i),soil_d_v_alb(i),soil_s_n_alb(i),soil_d_n_alb(i))
      enddo
#elif(defined SOIL_REFL_READ)
      allocate ( s_v_alb (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( d_v_alb (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( s_n_alb (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( d_n_alb (0:N_land_classification,1:lon_points,1:lat_points) )

! (1) Read in the albedo of visible of the saturated soil
      lndname = trim(dir_model_landdata)//'soil_s_v_alb'//trim(suffix)//'.bin'
      print*,lndname
      OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
      read(iunit,err=100) s_v_alb
      close(iunit)

! (1) Read in the albedo of visible of the dry soil
      lndname = trim(dir_model_landdata)//'soil_d_v_alb'//trim(suffix)//'.bin'
      print*,lndname
      OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
      read(iunit,err=100) d_v_alb
      close(iunit)

! (3) Read in the albedo of near infrared of the saturated soil
      lndname = trim(dir_model_landdata)//'soil_s_n_alb'//trim(suffix)//'.bin'
      print*,lndname
      OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
      read(iunit,err=100) s_n_alb
      close(iunit)

! (4) Read in the albedo of near infrared of the dry soil
      lndname = trim(dir_model_landdata)//'soil_d_n_alb'//trim(suffix)//'.bin'
      print*,lndname
      OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
      read(iunit,err=100) d_n_alb
      close(iunit)


      do npatch = 1, numpatch
      !DO npatch = 275192,275192
         i = patch2lon(npatch)
         j = patch2lat(npatch)
! yuan, 12/28/2019: Bug, add patchtype
         t = patchtype(npatch)
         m = patchclass(npatch)
         !if( m == 0 )then 
! for PFT structure, all the patches share the same soil data, storted at
! position 0 of array (0:N_land_classification,:,:)
#ifdef PFT_CLASSIFICATION
         IF ( t == 0 ) THEN 
             ! ocean 
             soil_s_v_alb (npatch) = s_v_alb (0,i,j)
             soil_d_v_alb (npatch) = d_v_alb (0,i,j)
             soil_s_n_alb (npatch) = s_n_alb (0,i,j)
             soil_d_n_alb (npatch) = d_n_alb (0,i,j)
         ELSE
             soil_s_v_alb (npatch) = s_v_alb (m,i,j)
             soil_d_v_alb (npatch) = d_v_alb (m,i,j)
             soil_s_n_alb (npatch) = s_n_alb (m,i,j)
             soil_d_n_alb (npatch) = d_n_alb (m,i,j)
         END IF
#else
         if( m == 0 )then 
             soil_s_v_alb (npatch) = -1.e36
             soil_d_v_alb (npatch) = -1.e36
             soil_s_n_alb (npatch) = -1.e36
             soil_d_n_alb (npatch) = -1.e36
          else
             soil_s_v_alb (npatch) = s_v_alb (m,i,j)
             soil_d_v_alb (npatch) = d_v_alb (m,i,j)
             soil_s_n_alb (npatch) = s_n_alb (m,i,j)
             soil_d_n_alb (npatch) = d_n_alb (m,i,j)
          end if
#endif
      end do

      deallocate ( s_v_alb )
      deallocate ( d_v_alb )
      deallocate ( s_n_alb )
      deallocate ( d_n_alb )
#endif

#if(defined CLMDEBUG)
      print*,'soil_s_v_alb =', minval(soil_s_v_alb, mask = soil_s_v_alb .gt. -1.e30), &
                               maxval(soil_s_v_alb, mask = soil_s_v_alb .gt. -1.e30)
      print*,'soil_d_v_alb =', minval(soil_d_v_alb, mask = soil_d_v_alb .gt. -1.e30), &
                               maxval(soil_d_v_alb, mask = soil_d_v_alb .gt. -1.e30)
      print*,'soil_s_n_alb =', minval(soil_s_n_alb, mask = soil_s_n_alb .gt. -1.e30), &
                               maxval(soil_s_n_alb, mask = soil_s_n_alb .gt. -1.e30)
      print*,'soil_d_n_alb =', minval(soil_d_n_alb, mask = soil_d_n_alb .gt. -1.e30), &
                               maxval(soil_d_n_alb, mask = soil_d_n_alb .gt. -1.e30)
#endif

      go to 1000
100   print 101,lndname
101   format(' error occured on file: ',a50)
1000  continue


END SUBROUTINE soil_parameters_readin
! --------------------------------------------------
! EOP
