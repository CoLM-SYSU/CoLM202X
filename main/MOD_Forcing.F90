#include <define.h>

!-----------------------------------------------------------------------
MODULE MOD_Forcing

! DESCRIPTION:
! read in the atmospheric forcing using user defined interpolation method
! or downscaling forcing
!
! REVISIONS:
! Yongjiu Dai and Hua Yuan, 04/2014: initial code from CoLM2014 (metdata.F90,
!                                    GETMET.F90 and rd_forcing.F90
!
! Shupeng Zhang, 05/2023: 1) porting codes to MPI parallel version
!                         2) codes for dealing with missing forcing value
!                         3) interface for downscaling
!
! TODO...(need complement)

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Grid
   USE MOD_Mapping_Grid2Pset
   USE MOD_InterpBilinear
   USE MOD_UserSpecifiedForcing
   USE MOD_TimeManager
   USE MOD_SPMD_Task
   USE MOD_MonthlyinSituCO2MaunaLoa
   USE MOD_Vars_Global, only : pi
   USE MOD_OrbCoszen

   IMPLICIT NONE

   type (grid_type), PUBLIC :: gforc

   type (mapping_grid2pset_type) :: mg2p_forc   ! area weighted mapping from forcing to model unit
   type (interp_bilinear_type)   :: forc_interp ! bilinear interpolation from forcing to model unit

   logical, allocatable :: forcmask_pch (:)
   logical, allocatable :: forcmask_elm (:)

   ! for Forcing_Downscaling
   logical, allocatable :: glacierss    (:)
   
   ! local variables
   integer  :: deltim_int                ! model time step length
   ! real(r8) :: deltim_real               ! model time step length

   !  for SinglePoint
   type(timestamp), allocatable :: forctime (:)
   integer,  allocatable :: iforctime(:)

   logical :: forcing_read_ahead
   real(r8), allocatable :: forc_disk(:,:)

   type(timestamp), allocatable :: tstamp_LB(:)  ! time stamp of low boundary data
   type(timestamp), allocatable :: tstamp_UB(:)  ! time stamp of up boundary data

   type(block_data_real8_2d) :: avgcos   ! time-average of cos(zenith)
   type(block_data_real8_2d) :: metdata  ! forcing data
#ifdef URBAN_MODEL
   type(block_data_real8_2d) :: rainf
   type(block_data_real8_2d) :: snowf
#endif

   type(block_data_real8_2d), allocatable :: forcn    (:)  ! forcing data
   type(block_data_real8_2d), allocatable :: forcn_LB (:)  ! forcing data at lower bondary
   type(block_data_real8_2d), allocatable :: forcn_UB (:)  ! forcing data at upper bondary

   PUBLIC :: forcing_init
   PUBLIC :: read_forcing

CONTAINS

   !--------------------------------
   SUBROUTINE forcing_init (dir_forcing, deltatime, ststamp, lc_year, etstamp)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_LandElm
   USE MOD_LandPatch
#ifdef CROP
   USE MOD_LandCrop
#endif
   USE MOD_Mapping_Grid2Pset
   USE MOD_UserSpecifiedForcing
   USE MOD_NetCDFSerial
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_1DForcing
   IMPLICIT NONE

   character(len=*), intent(in) :: dir_forcing
   real(r8),         intent(in) :: deltatime  ! model time step
   type(timestamp),  intent(in) :: ststamp
   integer, intent(in) :: lc_year    ! which year of land cover data used
   type(timestamp),  intent(in), optional :: etstamp

   ! Local variables
   integer            :: idate(3)
   character(len=256) :: filename, lndname, cyear
   integer            :: ivar, year, month, day, time_i
   real(r8)           :: missing_value
   integer            :: ielm, istt, iend

      CALL init_user_specified_forcing

      ! CO2 data initialization
      CALL init_monthly_co2_mlo

      ! get value of fmetdat and deltim
      deltim_int  = int(deltatime)
      ! deltim_real = deltatime

      ! set initial values
      IF (allocated(tstamp_LB)) deallocate(tstamp_LB)
      IF (allocated(tstamp_UB)) deallocate(tstamp_UB)
      allocate (tstamp_LB(NVAR))
      allocate (tstamp_UB(NVAR))
      tstamp_LB(:) = timestamp(-1, -1, -1)
      tstamp_UB(:) = timestamp(-1, -1, -1)

      idate = (/ststamp%year, ststamp%day, ststamp%sec/)

      CALL metread_latlon (dir_forcing, idate)

      IF (p_is_io) THEN

         IF (allocated(forcn   )) deallocate(forcn   )
         IF (allocated(forcn_LB)) deallocate(forcn_LB)
         IF (allocated(forcn_UB)) deallocate(forcn_UB)
         allocate (forcn    (NVAR))
         allocate (forcn_LB (NVAR))
         allocate (forcn_UB (NVAR))

         DO ivar = 1, NVAR
            CALL allocate_block_data (gforc, forcn   (ivar))
            CALL allocate_block_data (gforc, forcn_LB(ivar))
            CALL allocate_block_data (gforc, forcn_UB(ivar))
         ENDDO

         ! allocate memory for forcing data
         CALL allocate_block_data (gforc, metdata)  ! forcing data
         CALL allocate_block_data (gforc, avgcos )  ! time-average of cos(zenith)
#if(defined URBAN_MODEL && defined SinglePoint)
         CALL allocate_block_data (gforc, rainf)
         CALL allocate_block_data (gforc, snowf)
#endif

      ENDIF

      IF (DEF_USE_Forcing_Downscaling) THEN
         IF (p_is_worker) THEN
#if (defined CROP)
            CALL elm_patch%build (landelm, landpatch, use_frac = .true., sharedfrac = pctshrpch)
#else
            CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
#endif 
         ENDIF
      ENDIF

      IF (DEF_forcing%has_missing_value) THEN

         CALL setstampLB(ststamp, 1, year, month, day, time_i)
         filename = trim(dir_forcing)//trim(metfilename(year, month, day, 1))
         tstamp_LB(1) = timestamp(-1, -1, -1)

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               allocate (forcmask_pch(numpatch));  forcmask_pch(:) = .true.
            ENDIF
            IF (DEF_USE_Forcing_Downscaling) THEN
               IF (numelm > 0) THEN
                  allocate (forcmask_elm(numelm)); forcmask_elm(:) = .true.
               ENDIF
            ENDIF
         ENDIF

         IF (p_is_master) THEN
            CALL ncio_get_attr (filename, vname(1), trim(DEF_forcing%missing_value_name), missing_value)
         ENDIF
#ifdef USEMPI
         CALL mpi_bcast (missing_value, 1, MPI_REAL8, p_root, p_comm_glb, p_err)
#endif

         CALL ncio_read_block_time (filename, vname(1), gforc, time_i, metdata)

      ENDIF

      IF (.not. DEF_USE_Forcing_Downscaling) THEN

         IF (.not. DEF_forcing%has_missing_value) THEN

            IF (trim(DEF_Forcing_Interp) == 'areaweight') THEN
               CALL mg2p_forc%build (gforc, landpatch)
            ELSEIF (trim(DEF_Forcing_Interp) == 'bilinear') THEN
               CALL forc_interp%build (gforc, landpatch)
            ENDIF

         ELSE

            IF (trim(DEF_Forcing_Interp) == 'areaweight') THEN
               CALL mg2p_forc%build (gforc, landpatch, metdata, missing_value, forcmask_pch)
            ELSEIF (trim(DEF_Forcing_Interp) == 'bilinear') THEN
               CALL forc_interp%build (gforc, landpatch, metdata, missing_value, forcmask_pch)
            ENDIF

         ENDIF

      ELSE

         IF (.not. DEF_forcing%has_missing_value) THEN

            IF (trim(DEF_Forcing_Interp) == 'areaweight') THEN
               CALL mg2p_forc%build (gforc, landelm)
            ELSEIF (trim(DEF_Forcing_Interp) == 'bilinear') THEN
               CALL forc_interp%build (gforc, landelm)
            ENDIF

         ELSE

            IF (trim(DEF_Forcing_Interp) == 'areaweight') THEN
               CALL mg2p_forc%build (gforc, landelm, metdata, missing_value, forcmask_elm)
            ELSEIF (trim(DEF_Forcing_Interp) == 'bilinear') THEN
               CALL forc_interp%build (gforc, landelm, metdata, missing_value, forcmask_elm)
            ENDIF
         
            IF (p_is_worker) THEN
               DO ielm = 1, numelm
                  istt = elm_patch%substt(ielm)
                  iend = elm_patch%subend(ielm)
                  forcmask_pch(istt:iend) = forcmask_elm(ielm)
               ENDDO
            ENDIF

         ENDIF
      ENDIF

      IF (DEF_USE_Forcing_Downscaling) THEN
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               
               forc_topo = topoelv

               DO ielm = 1, numelm
                  istt = elm_patch%substt(ielm)
                  iend = elm_patch%subend(ielm)
                  forc_topo_elm(ielm) = sum(forc_topo(istt:iend) * elm_patch%subfrc(istt:iend))
               ENDDO

               allocate (glacierss(numpatch))
               glacierss(:) = patchtype(:) == 3

            ENDIF
         ENDIF
      ENDIF

      forcing_read_ahead = .false.
      IF (trim(DEF_forcing%dataset) == 'POINT') THEN
         IF (USE_SITE_ForcingReadAhead .and. present(etstamp)) THEN
            forcing_read_ahead = .true.
            CALL metread_time (dir_forcing, ststamp, etstamp, deltatime)
         ELSE
            CALL metread_time (dir_forcing)
         ENDIF
         allocate (iforctime(NVAR))
      ENDIF

      IF (trim(DEF_forcing%dataset) == 'POINT') THEN

         filename = trim(dir_forcing)//trim(fprefix(1))

#ifndef URBAN_MODEL
         IF (ncio_var_exist(filename,'reference_height_v')) THEN
            CALL ncio_read_serial (filename, 'reference_height_v', Height_V)
         ENDIF

         IF (ncio_var_exist(filename,'reference_height_t')) THEN
            CALL ncio_read_serial (filename, 'reference_height_t', Height_T)
         ENDIF

         IF (ncio_var_exist(filename,'reference_height_q')) THEN
            CALL ncio_read_serial (filename, 'reference_height_q', Height_Q)
         ENDIF
#else
         IF (ncio_var_exist(filename,'measurement_height_above_ground')) THEN
            CALL ncio_read_serial (filename, 'measurement_height_above_ground', Height_V)
            CALL ncio_read_serial (filename, 'measurement_height_above_ground', Height_T)
            CALL ncio_read_serial (filename, 'measurement_height_above_ground', Height_Q)
         ENDIF
#endif

      ENDIF

   END SUBROUTINE forcing_init

   ! ---- forcing finalize ----
   SUBROUTINE forcing_final ()

   IMPLICIT NONE

      IF (allocated(forcmask_pch)) deallocate(forcmask_pch)
      IF (allocated(forcmask_elm)) deallocate(forcmask_elm)
      IF (allocated(glacierss   )) deallocate(glacierss   )
      IF (allocated(forctime    )) deallocate(forctime    )
      IF (allocated(iforctime   )) deallocate(iforctime   )
      IF (allocated(forc_disk   )) deallocate(forc_disk   )
      IF (allocated(tstamp_LB   )) deallocate(tstamp_LB   )
      IF (allocated(tstamp_UB   )) deallocate(tstamp_UB   )

   END SUBROUTINE forcing_final

   ! ------------
   SUBROUTINE forcing_reset ()

   IMPLICIT NONE

      tstamp_LB(:) = timestamp(-1, -1, -1)
      tstamp_UB(:) = timestamp(-1, -1, -1)

   END SUBROUTINE forcing_reset

   ! ------------
   SUBROUTINE forcing_xy2vec (f_xy, f_vec)

   IMPLICIT NONE
      
      type(block_data_real8_2d) :: f_xy
      real(r8) :: f_vec(:)

      IF (trim(DEF_Forcing_Interp) == 'areaweight') THEN
         CALL mg2p_forc%map_aweighted (f_xy, f_vec)
      ELSEIF (trim(DEF_Forcing_Interp) == 'bilinear') THEN
         CALL forc_interp%interp (f_xy, f_vec)
      ENDIF

   END SUBROUTINE forcing_xy2vec

   !--------------------------------
   SUBROUTINE read_forcing (idate, dir_forcing)

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Const_Physical, only: rgas, grav
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_1DForcing
   USE MOD_Vars_2DForcing
   USE MOD_Block
   USE MOD_SPMD_Task
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_LandPatch
   USE MOD_Mapping_Grid2Pset
   USE MOD_RangeCheck
   USE MOD_UserSpecifiedForcing
   USE MOD_ForcingDownscaling, only : rair, cpair, downscale_forcings

   IMPLICIT NONE
   integer, intent(in) :: idate(3)
   character(len=*), intent(in) :: dir_forcing

   ! local variables:
   integer  :: ivar, istt, iend
   integer  :: iblkme, ib, jb, i, j, ilon, ilat, np, ne
   real(r8) :: calday  ! Julian cal day (1.xx to 365.xx)
   real(r8) :: sunang, cloud, difrat, vnrat
   real(r8) :: a, hsolar, ratio_rvrf
   type(block_data_real8_2d) :: forc_xy_solarin

   type(timestamp) :: mtstamp
   integer  :: id(3)
   integer  :: dtLB, dtUB
   real(r8) :: cosz
   integer  :: year, month, mday
   logical  :: has_u,has_v

   real solar, frl, prcp, tm, us, vs, pres, qm
   real(r8) :: pco2m

      IF (p_is_io) THEN

         !------------------------------------------------------------
         ! READ in THE ATMOSPHERIC FORCING

         ! read lower and upper boundary forcing data
         CALL metreadLBUB(idate, dir_forcing)

         ! set model time stamp
         id(:) = idate(:)
         !CALL adj2end(id)
         mtstamp = id

         has_u = .true.
         has_v = .true.
         ! loop for variables
         DO ivar = 1, NVAR

            IF (ivar == 5 .and. trim(vname(ivar)) == 'NULL') has_u = .false.
            IF (ivar == 6 .and. trim(vname(ivar)) == 'NULL') has_v = .false.
            IF (trim(vname(ivar)) == 'NULL') CYCLE     ! no data, CYCLE
            IF (trim(tintalgo(ivar)) == 'NULL') CYCLE

            ! to make sure the forcing data calculated is in the range of time
            ! interval [LB, UB]
            IF ( (mtstamp < tstamp_LB(ivar)) .or. (tstamp_UB(ivar) < mtstamp) ) THEN
               write(6, *) "the data required is out of range! STOP!"; CALL CoLM_stop()
            ENDIF

            ! calcualte distance to lower/upper boundary
            dtLB = mtstamp - tstamp_LB(ivar)
            dtUB = tstamp_UB(ivar) - mtstamp

            ! nearest method, for precipitation
            IF (tintalgo(ivar) == 'nearest') THEN
               IF (dtLB <= dtUB) THEN
                  CALL block_data_copy (forcn_LB(ivar), forcn(ivar))
               ELSE
                  CALL block_data_copy (forcn_UB(ivar), forcn(ivar))
               ENDIF
            ENDIF

            ! linear method, for T, Pres, Q, W, LW
            IF (tintalgo(ivar) == 'linear') THEN
               IF ( (dtLB+dtUB) > 0 ) THEN
                  CALL block_data_linear_interp ( &
                     forcn_LB(ivar), real(dtUB,r8)/real(dtLB+dtUB,r8), &
                     forcn_UB(ivar), real(dtLB,r8)/real(dtLB+dtUB,r8), &
                     forcn(ivar))
               ELSE
                  CALL block_data_copy (forcn_LB(ivar), forcn(ivar))
               ENDIF
            ENDIF

            ! coszen method, for SW
            IF (tintalgo(ivar) == 'coszen') THEN
               DO iblkme = 1, gblock%nblkme
                  ib = gblock%xblkme(iblkme)
                  jb = gblock%yblkme(iblkme)

                  DO j = 1, gforc%ycnt(jb)
                     DO i = 1, gforc%xcnt(ib)

                        ilat = gforc%ydsp(jb) + j
                        ilon = gforc%xdsp(ib) + i
                        IF (ilon > gforc%nlon) ilon = ilon - gforc%nlon

                        calday = calendarday(mtstamp)
                        cosz = orb_coszen(calday, gforc%rlon(ilon), gforc%rlat(ilat))
                        cosz = max(0.001, cosz)
                        forcn(ivar)%blk(ib,jb)%val(i,j) = &
                           cosz / avgcos%blk(ib,jb)%val(i,j) * forcn_LB(ivar)%blk(ib,jb)%val(i,j)

                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

         ENDDO

         ! preprocess for forcing data, only for QIAN data right now?
         CALL metpreprocess (gforc, forcn)

         CALL allocate_block_data (gforc, forc_xy_solarin)

         CALL block_data_copy (forcn(1), forc_xy_t      )
         CALL block_data_copy (forcn(2), forc_xy_q      )
         CALL block_data_copy (forcn(3), forc_xy_psrf   )
         CALL block_data_copy (forcn(3), forc_xy_pbot   )
         CALL block_data_copy (forcn(4), forc_xy_prl, sca = 2/3._r8)
         CALL block_data_copy (forcn(4), forc_xy_prc, sca = 1/3._r8)
         CALL block_data_copy (forcn(7), forc_xy_solarin)
         CALL block_data_copy (forcn(8), forc_xy_frl    )
         IF (DEF_USE_CBL_HEIGHT) THEN
            CALL block_data_copy (forcn(9), forc_xy_hpbl    )
         ENDIF

         IF (has_u .and. has_v) THEN
            CALL block_data_copy (forcn(5), forc_xy_us )
            CALL block_data_copy (forcn(6), forc_xy_vs )
         ELSEif (has_u) THEN
            CALL block_data_copy (forcn(5), forc_xy_us , sca = 1/sqrt(2.0_r8))
            CALL block_data_copy (forcn(5), forc_xy_vs , sca = 1/sqrt(2.0_r8))
         ELSEif (has_v) THEN
            CALL block_data_copy (forcn(6), forc_xy_us , sca = 1/sqrt(2.0_r8))
            CALL block_data_copy (forcn(6), forc_xy_vs , sca = 1/sqrt(2.0_r8))
         ELSE
            IF (.not.trim(DEF_forcing%dataset) == 'CPL7') THEN
               write(6, *) "At least one of the wind components must be provided! STOP!";
               CALL CoLM_stop()
            ENDIF
         ENDIF

         CALL flush_block_data (forc_xy_hgt_u, real(HEIGHT_V,r8))
         CALL flush_block_data (forc_xy_hgt_t, real(HEIGHT_T,r8))
         CALL flush_block_data (forc_xy_hgt_q, real(HEIGHT_Q,r8))

         IF (solarin_all_band) THEN

            IF (trim(DEF_forcing%dataset) == 'QIAN') THEN
               !---------------------------------------------------------------
               ! 04/2014, yuan: NOTE! codes from CLM4.5-CESM1.2.0
               ! relationship between incoming NIR or VIS radiation and ratio of
               ! direct to diffuse radiation calculated based on one year's worth of
               ! hourly CAM output from CAM version cam3_5_55
               !---------------------------------------------------------------
               DO iblkme = 1, gblock%nblkme
                  ib = gblock%xblkme(iblkme)
                  jb = gblock%yblkme(iblkme)

                  DO j = 1, gforc%ycnt(jb)
                     DO i = 1, gforc%xcnt(ib)

                        hsolar = forc_xy_solarin%blk(ib,jb)%val(i,j)*0.5_R8

                        ! NIR (dir, diff)
                        ratio_rvrf = min(0.99_R8,max(0.29548_R8 + 0.00504_R8*hsolar  &
                           -1.4957e-05_R8*hsolar**2 + 1.4881e-08_R8*hsolar**3,0.01_R8))
                        forc_xy_soll %blk(ib,jb)%val(i,j) = ratio_rvrf*hsolar
                        forc_xy_solld%blk(ib,jb)%val(i,j) = (1._R8 - ratio_rvrf)*hsolar

                        ! VIS (dir, diff)
                        ratio_rvrf = min(0.99_R8,max(0.17639_R8 + 0.00380_R8*hsolar  &
                           -9.0039e-06_R8*hsolar**2 + 8.1351e-09_R8*hsolar**3,0.01_R8))
                        forc_xy_sols %blk(ib,jb)%val(i,j) = ratio_rvrf*hsolar
                        forc_xy_solsd%blk(ib,jb)%val(i,j) = (1._R8 - ratio_rvrf)*hsolar

                     ENDDO
                  ENDDO
               ENDDO

            ELSE
               !---------------------------------------------------------------
               ! as the downward solar is in full band, an empirical expression
               ! will be used to divide fractions of band and incident
               ! (visible, near-infrad, dirct, diffuse)
               ! Julian calday (1.xx to 365.xx)
               !---------------------------------------------------------------
               DO iblkme = 1, gblock%nblkme
                  ib = gblock%xblkme(iblkme)
                  jb = gblock%yblkme(iblkme)

                  DO j = 1, gforc%ycnt(jb)
                     DO i = 1, gforc%xcnt(ib)

                        ilat = gforc%ydsp(jb) + j
                        ilon = gforc%xdsp(ib) + i
                        IF (ilon > gforc%nlon) ilon = ilon - gforc%nlon

                        a = forc_xy_solarin%blk(ib,jb)%val(i,j)
                        calday = calendarday(idate)
                        sunang = orb_coszen (calday, gforc%rlon(ilon), gforc%rlat(ilat))

                        cloud = (1160.*sunang-a)/(963.*sunang)
                        cloud = max(cloud,0.)
                        cloud = min(cloud,1.)
                        cloud = max(0.58,cloud)

                        difrat = 0.0604/(sunang-0.0223)+0.0683
                        IF(difrat.lt.0.) difrat = 0.
                        IF(difrat.gt.1.) difrat = 1.

                        difrat = difrat+(1.0-difrat)*cloud
                        vnrat = (580.-cloud*464.)/((580.-cloud*499.)+(580.-cloud*464.))

                        forc_xy_sols %blk(ib,jb)%val(i,j) = a*(1.0-difrat)*vnrat
                        forc_xy_soll %blk(ib,jb)%val(i,j) = a*(1.0-difrat)*(1.0-vnrat)
                        forc_xy_solsd%blk(ib,jb)%val(i,j) = a*difrat*vnrat
                        forc_xy_solld%blk(ib,jb)%val(i,j) = a*difrat*(1.0-vnrat)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

         ENDIF

         ! [GET ATMOSPHERE CO2 CONCENTRATION DATA]
         year  = idate(1)
         CALL julian2monthday (idate(1), idate(2), month, mday)
         pco2m = get_monthly_co2_mlo(year, month)*1.e-6
         CALL block_data_copy (forc_xy_pbot, forc_xy_pco2m, sca = pco2m        )
         CALL block_data_copy (forc_xy_pbot, forc_xy_po2m , sca = 0.209_r8     )

      ENDIF


      IF (.not. DEF_USE_Forcing_Downscaling) THEN

         ! Mapping the 2d atmospheric fields [lon_points]x[lat_points]
         !     -> the 1d vector of subgrid points [numpatch]
         CALL forcing_xy2vec (forc_xy_pco2m,  forc_pco2m)
         CALL forcing_xy2vec (forc_xy_po2m ,  forc_po2m )
         CALL forcing_xy2vec (forc_xy_us   ,  forc_us   )
         CALL forcing_xy2vec (forc_xy_vs   ,  forc_vs   )

         CALL forcing_xy2vec (forc_xy_psrf ,  forc_psrf )

         CALL forcing_xy2vec (forc_xy_sols ,  forc_sols )
         CALL forcing_xy2vec (forc_xy_soll ,  forc_soll )
         CALL forcing_xy2vec (forc_xy_solsd,  forc_solsd)
         CALL forcing_xy2vec (forc_xy_solld,  forc_solld)

         CALL forcing_xy2vec (forc_xy_hgt_t,  forc_hgt_t)
         CALL forcing_xy2vec (forc_xy_hgt_u,  forc_hgt_u)
         CALL forcing_xy2vec (forc_xy_hgt_q,  forc_hgt_q)

         IF (DEF_USE_CBL_HEIGHT) THEN
            CALL forcing_xy2vec (forc_xy_hpbl, forc_hpbl)
         ENDIF

         CALL forcing_xy2vec (forc_xy_t    ,  forc_t    )
         CALL forcing_xy2vec (forc_xy_q    ,  forc_q    )
         CALL forcing_xy2vec (forc_xy_prc  ,  forc_prc  )
         CALL forcing_xy2vec (forc_xy_prl  ,  forc_prl  )
         CALL forcing_xy2vec (forc_xy_pbot ,  forc_pbot )
         CALL forcing_xy2vec (forc_xy_frl  ,  forc_frl  )

         IF (p_is_worker) THEN

            DO np = 1, numpatch
               IF (DEF_forcing%has_missing_value) THEN
                  IF (.not. forcmask_pch(np)) CYCLE
               ENDIF

               ! The standard measuring conditions for temperature are two meters above the ground
               ! Scientists have measured the most frigid temperature ever
               ! recorded on the continent's eastern highlands: about (180K) colder than dry ice.
               IF(forc_t(np) < 180.) forc_t(np) = 180.
               ! the highest air temp was found in Kuwait 326 K, Sulaibya 2012-07-31;
               ! Pakistan, Sindh 2010-05-26; Iraq, Nasiriyah 2011-08-03
               IF(forc_t(np) > 326.) forc_t(np) = 326.

               forc_rhoair(np) = (forc_pbot(np) &
                  - 0.378*forc_q(np)*forc_pbot(np)/(0.622+0.378*forc_q(np)))&
                  / (rgas*forc_t(np))

            ENDDO

         ENDIF

      ELSE

         ! Mapping the 2d atmospheric fields [lon_points]x[lat_points]
         !     -> the 1d vector of subgrid points [numelm]
         CALL forcing_xy2vec (forc_xy_pco2m,  forc_pco2m_elm)
         CALL forcing_xy2vec (forc_xy_po2m ,  forc_po2m_elm )
         CALL forcing_xy2vec (forc_xy_us   ,  forc_us_elm   )
         CALL forcing_xy2vec (forc_xy_vs   ,  forc_vs_elm   )

         CALL forcing_xy2vec (forc_xy_psrf ,  forc_psrf_elm )

         CALL forcing_xy2vec (forc_xy_sols ,  forc_sols_elm )
         CALL forcing_xy2vec (forc_xy_soll ,  forc_soll_elm )
         CALL forcing_xy2vec (forc_xy_solsd,  forc_solsd_elm)
         CALL forcing_xy2vec (forc_xy_solld,  forc_solld_elm)

         CALL forcing_xy2vec (forc_xy_hgt_t,  forc_hgt_t_elm)
         CALL forcing_xy2vec (forc_xy_hgt_u,  forc_hgt_u_elm)
         CALL forcing_xy2vec (forc_xy_hgt_q,  forc_hgt_q_elm)

         IF (DEF_USE_CBL_HEIGHT) THEN
            CALL forcing_xy2vec (forc_xy_hpbl, forc_hpbl_elm)
         ENDIF

         CALL forcing_xy2vec (forc_xy_t    ,  forc_t_elm    )
         CALL forcing_xy2vec (forc_xy_q    ,  forc_q_elm    )
         CALL forcing_xy2vec (forc_xy_prc  ,  forc_prc_elm  )
         CALL forcing_xy2vec (forc_xy_prl  ,  forc_prl_elm  )
         CALL forcing_xy2vec (forc_xy_pbot ,  forc_pbot_elm )
         CALL forcing_xy2vec (forc_xy_frl  ,  forc_lwrad_elm)
         CALL forcing_xy2vec (forc_xy_hgt_t,  forc_hgt_elm  )

         IF (p_is_worker) THEN

            DO ne = 1, numelm
               IF (DEF_forcing%has_missing_value) THEN
                  IF (.not. forcmask_elm(ne)) CYCLE
               ENDIF

               istt = elm_patch%substt(ne)
               iend = elm_patch%subend(ne)
         
               forc_pco2m(istt:iend) = forc_pco2m_elm (ne)
               forc_po2m (istt:iend) = forc_po2m_elm  (ne)
               forc_us   (istt:iend) = forc_us_elm    (ne)
               forc_vs   (istt:iend) = forc_vs_elm    (ne)
                             
               forc_psrf (istt:iend) = forc_psrf_elm  (ne)
                             
               forc_sols (istt:iend) = forc_sols_elm  (ne)
               forc_soll (istt:iend) = forc_soll_elm  (ne)
               forc_solsd(istt:iend) = forc_solsd_elm (ne)
               forc_solld(istt:iend) = forc_solld_elm (ne)
                             
               forc_hgt_t(istt:iend) = forc_hgt_t_elm (ne)
               forc_hgt_u(istt:iend) = forc_hgt_u_elm (ne)
               forc_hgt_q(istt:iend) = forc_hgt_q_elm (ne)

               IF (DEF_USE_CBL_HEIGHT) THEN
                  forc_hpbl(istt:iend) = forc_hpbl_elm(ne)
               ENDIF

               ! The standard measuring conditions for temperature are two meters above the ground
               ! Scientists have measured the most frigid temperature ever
               ! recorded on the continent's eastern highlands: about (180K) colder than dry ice.
               IF(forc_t_elm(ne) < 180.) forc_t_elm(ne) = 180.
               ! the highest air temp was found in Kuwait 326 K, Sulaibya 2012-07-31;
               ! Pakistan, Sindh 2010-05-26; Iraq, Nasiriyah 2011-08-03
               IF(forc_t_elm(ne) > 326.) forc_t_elm(ne) = 326.

               forc_rho_elm(ne) = (forc_pbot_elm(ne) &
                  - 0.378*forc_q_elm(ne)*forc_pbot_elm(ne)/(0.622+0.378*forc_q_elm(ne)))&
                  / (rgas*forc_t_elm(ne))

               forc_th_elm(ne) = forc_t_elm(ne) * (1.e5/forc_pbot_elm(ne)) ** (rair/cpair)

            ENDDO

            CALL downscale_forcings ( &
               numelm, numpatch, elm_patch%substt, elm_patch%subend, glacierss, elm_patch%subfrc,   &
               ! forcing in gridcells
               forc_topo_elm, forc_t_elm,   forc_th_elm,  forc_q_elm,     forc_pbot_elm, &
               forc_rho_elm,  forc_prc_elm, forc_prl_elm, forc_lwrad_elm, forc_hgt_elm,  &
               ! forcing in patches
               forc_topo,     forc_t,       forc_th,      forc_q,         forc_pbot,     &
               forc_rhoair,   forc_prc,     forc_prl,     forc_frl)

         ENDIF

      ENDIF

#ifdef RangeCheck
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) write(*,'(/, A20)') 'Checking forcing ...'

      CALL check_vector_data ('Forcing us    [m/s]   ', forc_us   )
      CALL check_vector_data ('Forcing vs    [m/s]   ', forc_vs   )
      CALL check_vector_data ('Forcing t     [kelvin]', forc_t    )
      CALL check_vector_data ('Forcing q     [kg/kg] ', forc_q    )
      CALL check_vector_data ('Forcing prc   [mm/s]  ', forc_prc  )
      CALL check_vector_data ('Forcing psrf  [pa]    ', forc_psrf )
      CALL check_vector_data ('Forcing prl   [mm/s]  ', forc_prl  )
      CALL check_vector_data ('Forcing sols  [W/m2]  ', forc_sols )
      CALL check_vector_data ('Forcing soll  [W/m2]  ', forc_soll )
      CALL check_vector_data ('Forcing solsd [W/m2]  ', forc_solsd)
      CALL check_vector_data ('Forcing solld [W/m2]  ', forc_solld)
      CALL check_vector_data ('Forcing frl   [W/m2]  ', forc_frl  )
      IF (DEF_USE_CBL_HEIGHT) THEN
         CALL check_vector_data ('Forcing hpbl  ', forc_hpbl )
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
#endif

   END SUBROUTINE read_forcing


   ! ------------------------------------------------------------
   !
   ! !DESCRIPTION:
   !    read lower and upper boundary forcing data, a major interface of this
   !    MODULE
   !
   ! REVISIONS:
   ! Hua Yuan, 04/2014: initial code
   ! ------------------------------------------------------------
   SUBROUTINE metreadLBUB (idate, dir_forcing)

   USE MOD_UserSpecifiedForcing
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_DataType
   USE MOD_Block
   USE MOD_NetCDFBlock
   USE MOD_RangeCheck
   IMPLICIT NONE

   integer, intent(in) :: idate(3)
   character(len=*), intent(in) :: dir_forcing

   ! Local variables
   integer         :: ivar, year, month, day, time_i
   integer         :: iblkme, ib, jb, i, j
   type(timestamp) :: mtstamp
   character(len=256) :: filename

      mtstamp = idate

      DO ivar = 1, NVAR

         IF (trim(vname(ivar)) == 'NULL') CYCLE     ! no data, CYCLE

         ! lower and upper boundary data already exist, CYCLE
         IF ( .not.(tstamp_LB(ivar)=='NULL') .and. .not.(tstamp_UB(ivar)=='NULL') .and. &
            tstamp_LB(ivar)<=mtstamp .and. mtstamp<tstamp_UB(ivar) ) THEN
            CYCLE
         ENDIF

         ! set lower boundary time stamp and get data
         IF (tstamp_LB(ivar) == 'NULL') THEN
            CALL setstampLB(mtstamp, ivar, year, month, day, time_i)

            ! read forcing data
            filename = trim(dir_forcing)//trim(metfilename(year, month, day, ivar))
            IF (trim(DEF_forcing%dataset) == 'POINT') THEN

#ifndef URBAN_MODEL
               IF (forcing_read_ahead) THEN
                  metdata%blk(gblock%xblkme(1),gblock%yblkme(1))%val = forc_disk(time_i,ivar)
               ELSE
                  CALL ncio_read_site_time (filename, vname(ivar), time_i, metdata)
               ENDIF
#else
               IF (trim(vname(ivar)) == 'Rainf') THEN
                  CALL ncio_read_site_time (filename, 'Rainf', time_i, rainf)
                  CALL ncio_read_site_time (filename, 'Snowf', time_i, snowf)

                  DO iblkme = 1, gblock%nblkme
                     ib = gblock%xblkme(iblkme)
                     jb = gblock%yblkme(iblkme)

                     metdata%blk(ib,jb)%val(1,1) = rainf%blk(ib,jb)%val(1,1) + snowf%blk(ib,jb)%val(1,1)
                     !DO j = 1, gforc%ycnt(jb)
                     !   DO i = 1, gforc%xcnt(ib)
                     !      metdata%blk(ib,jb)%val(i,j) = rainf%blk(ib,jb)%val(i,j) &
                     !                                    +  snowf%blk(ib,jb)%val(i,j)
                     !   ENDDO
                     !ENDDO
                  ENDDO
               ELSE
                  CALL ncio_read_site_time (filename, vname(ivar), time_i, metdata)
               ENDIF
#endif
            ELSE
               CALL ncio_read_block_time (filename, vname(ivar), gforc, time_i, metdata)
            ENDIF

            CALL block_data_copy (metdata, forcn_LB(ivar))
         ENDIF

         ! set upper boundary time stamp and get data
         IF (tstamp_UB(ivar) == 'NULL' .or. tstamp_UB(ivar) <= mtstamp) THEN
            IF ( .not. (tstamp_UB(ivar) == 'NULL') ) THEN
               CALL block_data_copy (forcn_UB(ivar), forcn_LB(ivar))
            ENDIF
            CALL setstampUB(ivar, year, month, day, time_i)
            ! when reaching the END of forcing data, always reuse the last time step data
            IF (year <= endyr) THEN
               ! read forcing data
               filename = trim(dir_forcing)//trim(metfilename(year, month, day, ivar))
               IF (trim(DEF_forcing%dataset) == 'POINT') THEN

#ifndef URBAN_MODEL
                  IF (forcing_read_ahead) THEN
                     metdata%blk(gblock%xblkme(1),gblock%yblkme(1))%val = forc_disk(time_i,ivar)
                  ELSE
                     CALL ncio_read_site_time (filename, vname(ivar), time_i, metdata)
                  ENDIF
#else
                  IF (trim(vname(ivar)) == 'Rainf') THEN
                     CALL ncio_read_site_time (filename, 'Rainf', time_i, rainf)
                     CALL ncio_read_site_time (filename, 'Snowf', time_i, snowf)

                     DO iblkme = 1, gblock%nblkme
                        ib = gblock%xblkme(iblkme)
                        jb = gblock%yblkme(iblkme)

                        metdata%blk(ib,jb)%val(1,1) = rainf%blk(ib,jb)%val(1,1) + snowf%blk(ib,jb)%val(1,1)
                        !DO j = 1, gforc%ycnt(jb)
                        !   DO i = 1, gforc%xcnt(ib)
                        !      metdata%blk(ib,jb)%val(i,j) = rainf%blk(ib,jb)%val(i,j) &
                        !                                    +  snowf%blk(ib,jb)%val(i,j)
                        !   ENDDO
                        !ENDDO
                     ENDDO
                  ELSE
                     CALL ncio_read_site_time (filename, vname(ivar), time_i, metdata)
                  ENDIF
#endif
               ELSE
                  CALL ncio_read_block_time (filename, vname(ivar), gforc, time_i, metdata)
               ENDIF

               CALL block_data_copy (metdata, forcn_UB(ivar))
            ELSE
               write(*,*) year, endyr
               print *, 'NOTE: reaching the END of forcing data, always reuse the last time step data!'
            ENDIF
            IF (ivar == 7) THEN  ! calculate time average coszen, for shortwave radiation
               CALL calavgcos(idate)
            ENDIF
         ENDIF

      ENDDO

   END SUBROUTINE metreadLBUB


   !-------------------------------------------------
   SUBROUTINE metread_latlon (dir_forcing, idate)

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_UserSpecifiedForcing
   USE MOD_Namelist
   IMPLICIT NONE

   character(len=*), intent(in) :: dir_forcing
   integer, intent(in) :: idate(3)

   ! Local variables
   character(len=256) :: filename
   integer         :: year, month, day, time_i
   type(timestamp) :: mtstamp
   real(r8), allocatable :: latxy (:,:)    ! latitude values in 2d
   real(r8), allocatable :: lonxy (:,:)    ! longitude values in 2d
   real(r8), allocatable :: lon_in(:)
   real(r8), allocatable :: lat_in(:)

      IF (trim(DEF_forcing%dataset) == 'POINT' .or. trim(DEF_forcing%dataset) == 'CPL7' ) THEN
         CALL gforc%define_by_ndims (360, 180)
      ELSE

         mtstamp = idate

         CALL setstampLB(mtstamp, 1, year, month, day, time_i)
         filename = trim(dir_forcing)//trim(metfilename(year, month, day, 1))
         tstamp_LB(1) = timestamp(-1, -1, -1)

         IF (dim2d) THEN
            CALL ncio_read_bcast_serial (filename, latname, latxy)
            CALL ncio_read_bcast_serial (filename, lonname, lonxy)

            allocate (lat_in (size(latxy,2)))
            allocate (lon_in (size(lonxy,1)))
            lat_in = latxy(1,:)
            lon_in = lonxy(:,1)

            deallocate (latxy)
            deallocate (lonxy)
         ELSE
            CALL ncio_read_bcast_serial (filename, latname, lat_in)
            CALL ncio_read_bcast_serial (filename, lonname, lon_in)
         ENDIF

         IF (.not. DEF_forcing%regional) THEN
            CALL gforc%define_by_center (lat_in, lon_in)
         ELSE
            CALL gforc%define_by_center (lat_in, lon_in, &
               south = DEF_forcing%regbnd(1), north = DEF_forcing%regbnd(2), &
               west  = DEF_forcing%regbnd(3), east  = DEF_forcing%regbnd(4))
         ENDIF

         deallocate (lat_in)
         deallocate (lon_in)
      ENDIF

      CALL gforc%set_rlon ()
      CALL gforc%set_rlat ()

   END SUBROUTINE metread_latlon

   !-------------------------------------------------
   SUBROUTINE metread_time (dir_forcing, ststamp, etstamp, deltime)

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_UserSpecifiedForcing
   USE MOD_Namelist
   IMPLICIT NONE

   character(len=*), intent(in) :: dir_forcing
   type(timestamp),  intent(in), optional :: ststamp, etstamp
   real(r8),         intent(in), optional :: deltime

   ! Local variables
   character(len=256) :: filename
   character(len=256) :: timeunit, timestr
   real(r8), allocatable :: forctime_sec(:), metcache(:,:,:)
   integer :: year, month, day, hour, minute, second
   integer :: itime, maxday, id(3)
   integer*8 :: sec_long
   integer :: ivar, ntime, its, ite, it

   type(timestamp) :: etstamp_f
   type(timestamp), allocatable :: forctime_ (:)


      filename = trim(dir_forcing)//trim(fprefix(1))

      CALL ncio_read_serial (filename, 'time', forctime_sec)
      CALL ncio_get_attr    (filename, 'time', 'units', timeunit)

      timestr = timeunit(15:18) // ' ' // timeunit(20:21) // ' ' // timeunit(23:24) &
         // ' ' // timeunit(26:27) // ' ' // timeunit(29:30) // ' ' // timeunit(32:33)
      read(timestr,*) year, month, day, hour, minute, second

      allocate (forctime (size(forctime_sec)))

      forctime(1)%year = year
      forctime(1)%day  = get_calday(month*100+day, isleapyear(year))
      forctime(1)%sec = hour*3600 + minute*60 + second + forctime_sec(1)

      id(:) = (/forctime(1)%year, forctime(1)%day, forctime(1)%sec/)
      CALL adj2end(id)
      forctime(1) = id

      ntime = size(forctime)

      DO itime = 2, ntime
         id(:) = (/forctime(itime-1)%year, forctime(itime-1)%day, forctime(itime-1)%sec/)
         CALL ticktime (forctime_sec(itime)-forctime_sec(itime-1), id)
         forctime(itime) = id
      ENDDO

      IF (forcing_read_ahead) THEN

         CALL ticktime (deltime, id)
         etstamp_f = id

         IF ((ststamp < forctime(1)) .or. (etstamp_f < etstamp)) THEN
            write(*,*) 'Error: Forcing does not cover simulation period!'
            write(*,*) 'Model start ', ststamp,     ' -> Model END ', etstamp
            write(*,*) 'Forc  start ', forctime(1), ' -> Forc END  ', etstamp_f
            CALL CoLM_stop ()
         ELSE
            its = 1
            DO WHILE (.not. (ststamp < forctime(its+1)))
               its = its + 1
               IF (its >= ntime) EXIT
            ENDDO

            ite = ntime
            DO WHILE (etstamp < forctime(ite-1))
               ite = ite - 1
               IF (ite <= 1) EXIT
            ENDDO

            ntime = ite-its+1

            allocate (forctime_(ntime))
            DO it = 1, ntime
               forctime_(it) = forctime(it+its-1)
            ENDDO

            deallocate (forctime)
            allocate (forctime (ntime))
            DO it = 1, ntime
               forctime(it) = forctime_(it)
            ENDDO

            deallocate(forctime_)
         ENDIF

         allocate (forc_disk (size(forctime),NVAR))

         filename = trim(dir_forcing)//trim(metfilename(-1,-1,-1,-1))
         DO ivar = 1, NVAR
            IF (trim(vname(ivar)) /= 'NULL') THEN
               CALL ncio_read_period_serial (filename, vname(ivar), its, ite, metcache)
               forc_disk(:,ivar) = metcache(1,1,:)
            ENDIF
         ENDDO

         IF (allocated(metcache)) deallocate(metcache)
      ENDIF

   END SUBROUTINE metread_time

! ------------------------------------------------------------
!
! !DESCRIPTION:
!    set the lower boundary time stamp and record information,
!    a KEY FUNCTION of this MODULE
!
! - for time stamp, set it regularly as the model time step.
! - for record information, account for:
!    o year alternation
!    o month alternation
!    o leap year
!    o required dada just beyond the first record
!
! REVISIONS:
! Hua Yuan, 04/2014: initial code
! ------------------------------------------------------------
   SUBROUTINE setstampLB(mtstamp, var_i, year, month, mday, time_i)

   IMPLICIT NONE
   type(timestamp), intent(in)  :: mtstamp
   integer,         intent(in)  :: var_i
   integer,         intent(out) :: year
   integer,         intent(out) :: month
   integer,         intent(out) :: mday
   integer,         intent(out) :: time_i

   integer :: i, day, sec, ntime
   integer :: months(0:12)

      year = mtstamp%year
      day  = mtstamp%day
      sec  = mtstamp%sec

      IF (trim(DEF_forcing%dataset) == 'POINT') THEN

         ntime  = size(forctime)
         time_i = 1

         IF ((mtstamp < forctime(1)) .or. (forctime(ntime) < mtstamp)) THEN
            write(*,*) 'Error: Forcing does not cover simulation period!'
            write(*,*) 'Need ', mtstamp, ', Forc start ', forctime(1), ', Forc END', forctime(ntime)
            CALL CoLM_stop ()
         ELSE
            DO WHILE (.not. (mtstamp < forctime(time_i+1)))
               time_i = time_i + 1
            ENDDO
            iforctime(var_i) = time_i
            tstamp_LB(var_i) = forctime(iforctime(var_i))
         ENDIF

         RETURN
      ENDIF

      tstamp_LB(var_i)%year = year
      tstamp_LB(var_i)%day  = day

      ! in the case of one year one file
      IF ( trim(groupby) == 'year' ) THEN

         ! calculate the intitial second
         sec    = 86400*(day-1) + sec
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
         sec    = (time_i-1)*dtime(var_i) + offset(var_i) - 86400*(day-1)
         tstamp_LB(var_i)%sec = sec

         ! set time stamp (ststamp_LB)
         IF (sec < 0) THEN
            tstamp_LB(var_i)%sec = 86400 + sec
            tstamp_LB(var_i)%day = day - 1
            IF (tstamp_LB(var_i)%day == 0) THEN
               tstamp_LB(var_i)%year = year - 1
               IF ( isleapyear(tstamp_LB(var_i)%year) ) THEN
                  tstamp_LB(var_i)%day = 366
               ELSE
                  tstamp_LB(var_i)%day = 365
               ENDIF
            ENDIF
         ENDIF

         ! set record info (year, time_i)
         IF ( sec<0 .or. (sec==0 .and. offset(var_i).NE.0) ) THEN

            ! IF the required dada just behind the first record
            ! -> set to the first record
            IF ( year==startyr .and. month==startmo .and. day==1 ) THEN
               sec = offset(var_i)

               ! ELSE, set to one record backward
            ELSE
               sec = 86400 + sec
               day = day - 1
               IF (day == 0) THEN
                  year = year - 1
                  IF ( isleapyear(year) .and. leapyear) THEN
                     day = 366
                  ELSE
                     day = 365
                  ENDIF
               ENDIF
            ENDIF
         ENDIF ! ENDIF (sec <= 0)

         ! in case of leapyear with a non-leayyear calendar
         ! USE the data 1 day before after FEB 28th (Julian day 59).
         IF ( .not. leapyear .and. isleapyear(year) .and. day>59 ) THEN
            day = day - 1
         ENDIF

         ! get record time index
         sec = 86400*(day-1) + sec
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
      ENDIF

      ! in the case of one month one file
      IF ( trim(groupby) == 'month' ) THEN

         IF ( isleapyear(year) ) THEN
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         ELSE
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         ENDIF

         ! calculate initial month and day values
         CALL julian2monthday(year, day, month, mday)

         ! calculate initial second value
         sec    = 86400*(mday-1) + sec
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
         sec    = (time_i-1)*dtime(var_i) + offset(var_i) - 86400*(mday-1)
         tstamp_LB(var_i)%sec  = sec

         ! set time stamp (ststamp_LB)
         IF (sec < 0) THEN
            tstamp_LB(var_i)%sec = 86400 + sec
            tstamp_LB(var_i)%day = day - 1
            IF (tstamp_LB(var_i)%day == 0) THEN
               tstamp_LB(var_i)%year = year - 1
               IF ( isleapyear(tstamp_LB(var_i)%year) ) THEN
                  tstamp_LB(var_i)%day = 366
               ELSE
                  tstamp_LB(var_i)%day = 365
               ENDIF
            ENDIF
         ENDIF

         ! set record info (year, month, time_i)
         IF ( sec<0 .or. (sec==0 .and. offset(var_i).NE.0) ) THEN

            ! IF just behind the first record -> set to first record
            IF ( year==startyr .and. month==startmo .and. mday==1 ) THEN
               sec = offset(var_i)

               ! set to one record backward
            ELSE
               sec = 86400 + sec
               mday = mday - 1
               IF (mday == 0) THEN
                  month = month - 1
                  ! bug found by Zhu Siguang & Zhang Xiangxiang, 05/19/2014
                  ! move the below line in the 'ELSE' statement
                  !mday = months(month) - months(month-1)
                  IF (month == 0) THEN
                     month = 12
                     year = year - 1
                     mday = 31
                  ELSE
                     mday = months(month) - months(month-1)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leayyear calendar
         ! USE the data 1 day before, i.e., FEB 28th.
         IF ( .not. leapyear .and. isleapyear(year) .and. month==2 .and. mday==29 ) THEN
            mday = 28
         ENDIF

         ! get record time index
         sec = 86400*(mday-1) + sec
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
      ENDIF

      ! in the case of one day one file
      IF ( trim(groupby) == 'day' ) THEN

         ! calculate initial month and day values
         CALL julian2monthday(year, day, month, mday)

         ! calculate initial second value
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
         sec    = (time_i-1)*dtime(var_i) + offset(var_i)
         tstamp_LB(var_i)%sec  = sec

         ! set time stamp (ststamp_LB)
         IF (sec < 0) THEN
            tstamp_LB(var_i)%sec = 86400 + sec
            tstamp_LB(var_i)%day = day - 1
            IF (tstamp_LB(var_i)%day == 0) THEN
               tstamp_LB(var_i)%year = year - 1
               IF ( isleapyear(tstamp_LB(var_i)%year) ) THEN
                  tstamp_LB(var_i)%day = 366
               ELSE
                  tstamp_LB(var_i)%day = 365
               ENDIF
            ENDIF

            IF ( year==startyr .and. month==startmo .and. mday==1 ) THEN
               sec = offset(var_i)
            ! set to one record backward
            ELSE
               sec = 86400 + sec
               year = tstamp_LB(var_i)%year
               CALL julian2monthday(tstamp_LB(var_i)%year, tstamp_LB(var_i)%day, month, mday)
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leayyear calendar
         ! USE the data 1 day before, i.e., FEB 28th.
         IF ( .not. leapyear .and. isleapyear(year) .and. month==2 .and. mday==29 ) THEN
            mday = 28
         ENDIF

         ! get record time index
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
      ENDIF

      IF (time_i <= 0) THEN
         write(6, *) "got the wrong time record of forcing! STOP!"; CALL CoLM_stop()
      ENDIF

      RETURN

   END SUBROUTINE setstampLB

! ------------------------------------------------------------
!
! !DESCRIPTION:
!    set the upper boundary time stamp and record information,
!    a KEY FUNCTION of this MODULE
!
! REVISIONS:
! Hua Yuan, 04/2014: initial code
! ------------------------------------------------------------
   SUBROUTINE setstampUB(var_i, year, month, mday, time_i)

   IMPLICIT NONE
   integer,         intent(in)  :: var_i
   integer,         intent(out) :: year
   integer,         intent(out) :: month
   integer,         intent(out) :: mday
   integer,         intent(out) :: time_i

   integer :: day, sec
   integer :: months(0:12)

      IF (trim(DEF_forcing%dataset) == 'POINT') THEN
         IF ( tstamp_UB(var_i) == 'NULL' ) THEN
            tstamp_UB(var_i) = forctime(iforctime(var_i)+1)
         ELSE
            iforctime(var_i) = iforctime(var_i) + 1
            tstamp_LB(var_i) = forctime(iforctime(var_i))
            tstamp_UB(var_i) = forctime(iforctime(var_i)+1)
         ENDIF

         time_i = iforctime(var_i)+1
         year = tstamp_UB(var_i)%year
         RETURN
      ENDIF

      ! calculate the time stamp
      IF ( tstamp_UB(var_i) == 'NULL' ) THEN
         tstamp_UB(var_i) = tstamp_LB(var_i) + dtime(var_i)
      ELSE
         tstamp_LB(var_i) = tstamp_UB(var_i)
         tstamp_UB(var_i) = tstamp_UB(var_i) + dtime(var_i)
      ENDIF

      ! calcualte initial year, day, and second values
      year = tstamp_UB(var_i)%year
      day  = tstamp_UB(var_i)%day
      sec  = tstamp_UB(var_i)%sec

      IF ( trim(groupby) == 'year' ) THEN

         ! adjust year value
         IF ( sec==86400 .and. offset(var_i).eq.0 ) THEN
            sec = 0
            day = day + 1
            IF( isleapyear(year) .and. day==367) THEN
               year = year + 1; day = 1
            ENDIF
            IF( .not. isleapyear(year) .and. day==366) THEN
               year = year + 1; day = 1
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leayyear calendar
         ! USE the data 1 day before after FEB 28th (Julian day 59).
         IF ( .not. leapyear .and. isleapyear(year) .and. day>59 ) THEN
            day = day - 1
         ENDIF

         ! set record index
         sec = 86400*(day-1) + sec
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
      ENDIF

      IF ( trim(groupby) == 'month' ) THEN

         IF ( isleapyear(year) ) THEN
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         ELSE
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         ENDIF

         ! calculate initial month and day values
         CALL julian2monthday(year, day, month, mday)

         ! record in the next day, adjust year, month and second values
         IF ( sec==86400 .and. offset(var_i).eq.0 ) THEN
            sec  = 0
            mday = mday + 1
            IF ( mday > (months(month)-months(month-1)) ) THEN
               mday = 1
               ! bug found by Zhu Siguang, 05/25/2014
               ! move the below line in the 'ELSE' statement
               !month = month + 1
               IF (month == 12) THEN
                  month = 1
                  year = year + 1
               ELSE
                  month = month + 1
               ENDIF
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leayyear calendar
         ! for day 29th Feb, USE the data 1 day before, i.e., 28th FEB.
         IF ( .not. leapyear .and. isleapyear(year) .and. month==2 .and. mday==29 ) THEN
            mday = 28
         ENDIF

         ! set record index
         sec    = 86400*(mday-1) + sec
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
      ENDIF

      IF ( trim(groupby) == 'day' ) THEN
         IF ( isleapyear(year) ) THEN
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         ELSE
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         ENDIF

         ! calculate initial month and day values
         CALL julian2monthday(year, day, month, mday)
         !mday = day

         ! record in the next day, adjust year, month and second values
         IF ( sec==86400 .and. offset(var_i).eq.0 ) THEN
            sec  = 0
            mday = mday + 1
            IF ( mday > (months(month)-months(month-1)) ) THEN
               mday = 1
               ! bug found by Zhu Siguang, 05/25/2014
               ! move the below line in the 'ELSE' statement
               !month = month + 1
               IF (month == 12) THEN
                  month = 1
                  year = year + 1
               ELSE
                  month = month + 1
               ENDIF
            ENDIF
         ENDIF

         ! in case of leapyear with a non-leayyear calendar
         ! for day 29th Feb, USE the data 1 day before, i.e., 28th FEB.
         IF ( .not. leapyear .and. isleapyear(year) .and. month==2 .and. mday==29 ) THEN
            mday = 28
         ENDIF

         ! set record index
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
      ENDIF

      IF (time_i < 0) THEN
         write(6, *) "got the wrong time record of forcing! STOP!"; CALL CoLM_stop()
      ENDIF

      RETURN

   END SUBROUTINE setstampUB

! ------------------------------------------------------------
! !DESCRIPTION:
! calculate time average coszen value bwteeen [LB, UB]
!
! REVISIONS:
! 04/2014, yuan: this method is adapted from CLM
! ------------------------------------------------------------
   SUBROUTINE calavgcos(idate)

   USE MOD_Block
   USE MOD_DataType
   IMPLICIT NONE

   integer, intent(in) :: idate(3)

   integer  :: ntime, iblkme, ib, jb, i, j, ilon, ilat
   real(r8) :: calday, cosz
   type(timestamp) :: tstamp

      tstamp = idate !tstamp_LB(7)
      ntime = 0
      DO WHILE (tstamp < tstamp_UB(7))
         ntime  = ntime + 1
         tstamp = tstamp + deltim_int
      ENDDO

      tstamp = idate !tstamp_LB(7)
      CALL flush_block_data (avgcos, 0._r8)

      DO WHILE (tstamp < tstamp_UB(7))

         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)
            DO j = 1, gforc%ycnt(jb)
               DO i = 1, gforc%xcnt(ib)

                  ilat = gforc%ydsp(jb) + j
                  ilon = gforc%xdsp(ib) + i
                  IF (ilon > gforc%nlon) ilon = ilon - gforc%nlon

                  calday = calendarday(tstamp)
                  cosz = orb_coszen(calday, gforc%rlon(ilon), gforc%rlat(ilat))
                  cosz = max(0.001, cosz)
                  avgcos%blk(ib,jb)%val(i,j) = avgcos%blk(ib,jb)%val(i,j) &
                     + cosz / real(ntime,r8) !  * deltim_real /real(tstamp_UB(7)-tstamp_LB(7))

               ENDDO
            ENDDO
         ENDDO

         tstamp = tstamp + deltim_int

      ENDDO

   END SUBROUTINE calavgcos

END MODULE MOD_Forcing
