#include <define.h>

!-----------------------------------------------------------------------
module MOD_Forcing

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

   use MOD_Precision
   USE MOD_Namelist
   use MOD_Grid
   use MOD_Mapping_Grid2Pset
   use MOD_UserSpecifiedForcing
   use MOD_TimeManager
   use MOD_SPMD_Task
   USE MOD_MonthlyinSituCO2MaunaLoa
   USE MOD_Vars_Global, only : pi
   USE MOD_OrbCoszen

   implicit none

   type (grid_type), public :: gforc
   type (mapping_grid2pset_type) :: mg2p_forc

   LOGICAL, allocatable :: forcmask (:)

   ! for Forcing_Downscaling
   type (mapping_grid2pset_type) :: mg2p_forc_elm
   LOGICAL, allocatable :: forcmask_elm (:)
   LOGICAL, allocatable :: glacierss    (:)

   ! local variables
   integer  :: deltim_int                ! model time step length
   real(r8) :: deltim_real               ! model time step length

   !  for SinglePoint
   TYPE(timestamp), allocatable :: forctime (:)
   INTEGER, allocatable :: iforctime(:)

   type(timestamp), allocatable :: tstamp_LB(:)  ! time stamp of low boundary data
   type(timestamp), allocatable :: tstamp_UB(:)  ! time stamp of up boundary data

   type(block_data_real8_2d) :: avgcos   ! time-average of cos(zenith)
   type(block_data_real8_2d) :: metdata  ! forcing data

   type(block_data_real8_2d), allocatable :: forcn    (:)  ! forcing data
   type(block_data_real8_2d), allocatable :: forcn_LB (:)  ! forcing data at lower bondary
   type(block_data_real8_2d), allocatable :: forcn_UB (:)  ! forcing data at upper bondary

   public :: forcing_init
   public :: read_forcing

contains

   !--------------------------------
   subroutine forcing_init (dir_forcing, deltatime, idate, lc_year)

      use MOD_SPMD_Task
      USE MOD_Namelist
      use MOD_DataType
      USE MOD_Mesh
      USE MOD_LandElm
      USE MOD_LandPatch
      use MOD_Mapping_Grid2Pset
      use MOD_UserSpecifiedForcing
      USE MOD_NetCDFSerial
      USE MOD_NetCDFVector
      USE MOD_NetCDFBlock
      USE MOD_Vars_TimeInvariants
      USE MOD_Vars_1DForcing
      implicit none

      character(len=*), intent(in) :: dir_forcing
      real(r8), intent(in) :: deltatime  ! model time step
      integer,  intent(in) :: idate(3)
      INTEGER, intent(in) :: lc_year    ! which year of land cover data used

      ! Local variables
      CHARACTER(len=256) :: filename, lndname, cyear
      type(timestamp)    :: mtstamp
      integer            :: ivar, year, month, day, time_i
      REAL(r8)           :: missing_value
      INTEGER            :: ielm, istt, iend

      call init_user_specified_forcing

      ! CO2 data initialization
      CALL init_monthly_co2_mlo

      ! get value of fmetdat and deltim
      deltim_int  = int(deltatime)
      deltim_real = deltatime

      ! set initial values
      IF (allocated(tstamp_LB)) deallocate(tstamp_LB)
      IF (allocated(tstamp_UB)) deallocate(tstamp_UB)
      allocate (tstamp_LB(NVAR))
      allocate (tstamp_UB(NVAR))
      tstamp_LB(:) = timestamp(-1, -1, -1)
      tstamp_UB(:) = timestamp(-1, -1, -1)

      call metread_latlon (dir_forcing, idate)

      if (p_is_io) then

         IF (allocated(forcn   )) deallocate(forcn   )
         IF (allocated(forcn_LB)) deallocate(forcn_LB)
         IF (allocated(forcn_UB)) deallocate(forcn_UB)
         allocate (forcn    (NVAR))
         allocate (forcn_LB (NVAR))
         allocate (forcn_UB (NVAR))

         do ivar = 1, NVAR
            call allocate_block_data (gforc, forcn   (ivar))
            call allocate_block_data (gforc, forcn_LB(ivar))
            call allocate_block_data (gforc, forcn_UB(ivar))
         end do

         ! allocate memory for forcing data
         call allocate_block_data (gforc, metdata)  ! forcing data
         call allocate_block_data (gforc, avgcos )  ! time-average of cos(zenith)

      end if

      IF (.not. DEF_forcing%has_missing_value) THEN
         call mg2p_forc%build (gforc, landpatch)
         IF (DEF_USE_Forcing_Downscaling) THEN
            call mg2p_forc_elm%build (gforc, landelm)
         ENDIF
      ELSE
         mtstamp = idate
         call setstampLB(mtstamp, 1, year, month, day, time_i)
         filename = trim(dir_forcing)//trim(metfilename(year, month, day, 1))
         tstamp_LB(1) = timestamp(-1, -1, -1)

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               allocate (forcmask(numpatch))
               forcmask(:) = .true.
            ENDIF
            IF (DEF_USE_Forcing_Downscaling) THEN
               IF (numelm > 0) THEN
                  allocate (forcmask_elm(numelm))
                  forcmask_elm(:) = .true.
               ENDIF
            ENDIF
         ENDIF

         IF (p_is_master) THEN
            CALL ncio_get_attr (filename, vname(1), 'missing_value', missing_value)
         ENDIF
#ifdef USEMPI
         CALL mpi_bcast (missing_value, 1, MPI_REAL8, p_root, p_comm_glb, p_err)
#endif

         call ncio_read_block_time (filename, vname(1), gforc, time_i, metdata)
         call mg2p_forc%build (gforc, landpatch, metdata, missing_value, forcmask)
         IF (DEF_USE_Forcing_Downscaling) THEN
            call mg2p_forc_elm%build (gforc, landelm, metdata, missing_value, forcmask_elm)
         ENDIF
      ENDIF

      IF (DEF_USE_Forcing_Downscaling) THEN

         write(cyear,'(i4.4)') lc_year
         lndname = trim(DEF_dir_landdata) // '/topography/'//trim(cyear)//'/topography_patches.nc'
         call ncio_read_vector (lndname, 'topography_patches', landpatch, forc_topo)

         IF (p_is_worker) THEN
#if (defined CROP)
            CALL elm_patch%build (landelm, landpatch, use_frac = .true., shadowfrac = pctcrop)
#else
            CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
#endif

            DO ielm = 1, numelm
               istt = elm_patch%substt(ielm)
               iend = elm_patch%subend(ielm)
               forc_topo_elm(ielm) = sum(forc_topo(istt:iend) * elm_patch%subfrc(istt:iend))
            ENDDO

            IF (numpatch > 0) THEN
               allocate (glacierss(numpatch))
               glacierss(:) = patchtype(:) == 3
            ENDIF
         ENDIF

      ENDIF

      IF (trim(DEF_forcing%dataset) == 'POINT') THEN
         CALL metread_time (dir_forcing)
         allocate (iforctime(NVAR))
      ENDIF

   end subroutine forcing_init

   ! ------------
   SUBROUTINE forcing_reset ()

      IMPLICIT NONE

      tstamp_LB(:) = timestamp(-1, -1, -1)
      tstamp_UB(:) = timestamp(-1, -1, -1)

   END SUBROUTINE forcing_reset

   !--------------------------------
   SUBROUTINE read_forcing (idate, dir_forcing)

      use MOD_Precision
      use MOD_Namelist
      use MOD_Const_Physical, only: rgas, grav
      use MOD_Vars_TimeInvariants
      use MOD_Vars_1DForcing
      use MOD_Vars_2DForcing
      use MOD_Block
      use MOD_SPMD_Task
      use MOD_DataType
      use MOD_Mesh
      use MOD_LandPatch
      use MOD_Mapping_Grid2Pset
      use MOD_RangeCheck
      use MOD_UserSpecifiedForcing
      USE MOD_ForcingDownscaling, only : rair, cpair, downscale_forcings

      IMPLICIT NONE
      integer, INTENT(in) :: idate(3)
      character(len=*), intent(in) :: dir_forcing

      ! local variables:
      integer  :: ivar
      integer  :: iblkme, ib, jb, i, j, ilon, ilat, np, ne
      real(r8) :: calday  ! Julian cal day (1.xx to 365.xx)
      real(r8) :: sunang, cloud, difrat, vnrat
      real(r8) :: a, hsolar, ratio_rvrf
      type(block_data_real8_2d) :: forc_xy_solarin

      type(timestamp) :: mtstamp
      integer  :: id(3)
      integer  :: dtLB, dtUB
      real(r8) :: cosz
      INTEGER  :: year, month, mday
      logical  :: has_u,has_v

      real solar, frl, prcp, tm, us, vs, pres, qm
      real(r8) :: pco2m

      if (p_is_io) then

         !------------------------------------------------------------
         ! READ IN THE ATMOSPHERIC FORCING

         ! read lower and upper boundary forcing data
         CALL metreadLBUB(idate, dir_forcing)

         ! set model time stamp
         id(:) = idate(:)
         call adj2end(id)
         mtstamp = id

         has_u = .true.
         has_v = .true.
         ! loop for variables
         do ivar = 1, NVAR

            IF (ivar == 5 .and. trim(vname(ivar)) == 'NULL') has_u = .false.
            IF (ivar == 6 .and. trim(vname(ivar)) == 'NULL') has_v = .false.
            if (trim(vname(ivar)) == 'NULL') cycle     ! no data, cycle
            if (trim(tintalgo(ivar)) == 'NULL') cycle

            ! to make sure the forcing data calculated is in the range of time
            ! interval [LB, UB]
            if ( (mtstamp < tstamp_LB(ivar)) .or. (tstamp_UB(ivar) < mtstamp) ) then
               write(6, *) "the data required is out of range! stop!"; stop
            end if

            ! calcualte distance to lower/upper boundary
            dtLB = mtstamp - tstamp_LB(ivar)
            dtUB = tstamp_UB(ivar) - mtstamp

            ! nearest method, for precipitation
            if (tintalgo(ivar) == 'nearest') then
               if (dtLB <= dtUB) then
                  call block_data_copy (forcn_LB(ivar), forcn(ivar))
               else
                  call block_data_copy (forcn_UB(ivar), forcn(ivar))
               end if
            end if

            ! linear method, for T, Pres, Q, W, LW
            if (tintalgo(ivar) == 'linear') then
               if ( (dtLB+dtUB) > 0 ) then
                  call block_data_linear_interp ( &
                     forcn_LB(ivar), real(dtUB,r8)/real(dtLB+dtUB,r8), &
                     forcn_UB(ivar), real(dtLB,r8)/real(dtLB+dtUB,r8), &
                     forcn(ivar))
               else
                  call block_data_copy (forcn_LB(ivar), forcn(ivar))
               end if
            end if

            ! coszen method, for SW
            if (tintalgo(ivar) == 'coszen') then
               DO iblkme = 1, gblock%nblkme
                  ib = gblock%xblkme(iblkme)
                  jb = gblock%yblkme(iblkme)

                  do j = 1, gforc%ycnt(jb)
                     do i = 1, gforc%xcnt(ib)

                        ilat = gforc%ydsp(jb) + j
                        ilon = gforc%xdsp(ib) + i
                        if (ilon > gforc%nlon) ilon = ilon - gforc%nlon

                        calday = calendarday(mtstamp)
                        cosz = orb_coszen(calday, gforc%rlon(ilon), gforc%rlat(ilat))
                        cosz = max(0.001, cosz)
                        forcn(ivar)%blk(ib,jb)%val(i,j) = &
                           cosz / avgcos%blk(ib,jb)%val(i,j) * forcn_LB(ivar)%blk(ib,jb)%val(i,j)

                     end do
                  end do
               end do
            end if

         end do

         ! preprocess for forcing data, only for QIAN data right now?
         CALL metpreprocess (gforc, forcn)

         call allocate_block_data (gforc, forc_xy_solarin)

         call block_data_copy (forcn(1), forc_xy_t      )
         call block_data_copy (forcn(2), forc_xy_q      )
         call block_data_copy (forcn(3), forc_xy_psrf   )
         call block_data_copy (forcn(3), forc_xy_pbot   )
         call block_data_copy (forcn(4), forc_xy_prl, sca = 2/3._r8)
         call block_data_copy (forcn(4), forc_xy_prc, sca = 1/3._r8)
         call block_data_copy (forcn(7), forc_xy_solarin)
         call block_data_copy (forcn(8), forc_xy_frl    )
         if (DEF_USE_CBL_HEIGHT) then
            call block_data_copy (forcn(9), forc_xy_hpbl    )
         endif

         if (has_u .and. has_v) then
            call block_data_copy (forcn(5), forc_xy_us )
            call block_data_copy (forcn(6), forc_xy_vs )
         ELSEif (has_u) then
            call block_data_copy (forcn(5), forc_xy_us , sca = 1/sqrt(2.0_r8))
            call block_data_copy (forcn(5), forc_xy_vs , sca = 1/sqrt(2.0_r8))
         ELSEif (has_v) then
            call block_data_copy (forcn(6), forc_xy_us , sca = 1/sqrt(2.0_r8))
            call block_data_copy (forcn(6), forc_xy_vs , sca = 1/sqrt(2.0_r8))
         ELSE
            write(6, *) "At least one of the wind components must be provided! stop!"; stop
         ENDIF

         call flush_block_data (forc_xy_hgt_u, real(HEIGHT_V,r8))
         call flush_block_data (forc_xy_hgt_t, real(HEIGHT_T,r8))
         call flush_block_data (forc_xy_hgt_q, real(HEIGHT_Q,r8))

         if (solarin_all_band) then

            if (trim(DEF_forcing%dataset) == 'QIAN') then
               !---------------------------------------------------------------
               ! 04/2014, yuan: NOTE! codes from CLM4.5-CESM1.2.0
               ! relationship between incoming NIR or VIS radiation and ratio of
               ! direct to diffuse radiation calculated based on one year's worth of
               ! hourly CAM output from CAM version cam3_5_55
               !---------------------------------------------------------------
               DO iblkme = 1, gblock%nblkme
                  ib = gblock%xblkme(iblkme)
                  jb = gblock%yblkme(iblkme)

                  do j = 1, gforc%ycnt(jb)
                     do i = 1, gforc%xcnt(ib)

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

                     end do
                  end do
               end do

            else
               !---------------------------------------------------------------
               ! as the downward solar is in full band, an empirical expression
               ! will be used to divide fractions of band and incident
               ! (visible, near-infrad, dirct, diffuse)
               ! Julian calday (1.xx to 365.xx)
               !---------------------------------------------------------------
               DO iblkme = 1, gblock%nblkme
                  ib = gblock%xblkme(iblkme)
                  jb = gblock%yblkme(iblkme)

                  do j = 1, gforc%ycnt(jb)
                     do i = 1, gforc%xcnt(ib)

                        ilat = gforc%ydsp(jb) + j
                        ilon = gforc%xdsp(ib) + i
                        if (ilon > gforc%nlon) ilon = ilon - gforc%nlon

                        a = forc_xy_solarin%blk(ib,jb)%val(i,j)
                        calday = calendarday(idate)
                        sunang = orb_coszen (calday, gforc%rlon(ilon), gforc%rlat(ilat))

                        cloud = (1160.*sunang-a)/(963.*sunang)
                        cloud = max(cloud,0.)
                        cloud = min(cloud,1.)
                        cloud = max(0.58,cloud)

                        difrat = 0.0604/(sunang-0.0223)+0.0683
                        if(difrat.lt.0.) difrat = 0.
                        if(difrat.gt.1.) difrat = 1.

                        difrat = difrat+(1.0-difrat)*cloud
                        vnrat = (580.-cloud*464.)/((580.-cloud*499.)+(580.-cloud*464.))

                        forc_xy_sols %blk(ib,jb)%val(i,j) = a*(1.0-difrat)*vnrat
                        forc_xy_soll %blk(ib,jb)%val(i,j) = a*(1.0-difrat)*(1.0-vnrat)
                        forc_xy_solsd%blk(ib,jb)%val(i,j) = a*difrat*vnrat
                        forc_xy_solld%blk(ib,jb)%val(i,j) = a*difrat*(1.0-vnrat)
                     end do
                  end do
               end do
            end if

         end if

         ! [GET ATMOSPHERE CO2 CONCENTRATION DATA]
         year  = idate(1)
         CALL julian2monthday (idate(1), idate(2), month, mday)
         pco2m = get_monthly_co2_mlo(year, month)*1.e-6
         call block_data_copy (forc_xy_pbot, forc_xy_pco2m, sca = pco2m        )
         call block_data_copy (forc_xy_pbot, forc_xy_po2m , sca = 0.209_r8     )

      end if

      ! Mapping the 2d atmospheric fields [lon_points]x[lat_points]
      !     -> the 1d vector of subgrid points [numpatch]
      call mg2p_forc%map_aweighted (forc_xy_pco2m,  forc_pco2m)
      call mg2p_forc%map_aweighted (forc_xy_po2m ,  forc_po2m )
      call mg2p_forc%map_aweighted (forc_xy_us   ,  forc_us   )
      call mg2p_forc%map_aweighted (forc_xy_vs   ,  forc_vs   )

      call mg2p_forc%map_aweighted (forc_xy_psrf ,  forc_psrf )

      call mg2p_forc%map_aweighted (forc_xy_sols ,  forc_sols )
      call mg2p_forc%map_aweighted (forc_xy_soll ,  forc_soll )
      call mg2p_forc%map_aweighted (forc_xy_solsd,  forc_solsd)
      call mg2p_forc%map_aweighted (forc_xy_solld,  forc_solld)

      call mg2p_forc%map_aweighted (forc_xy_hgt_t,  forc_hgt_t)
      call mg2p_forc%map_aweighted (forc_xy_hgt_u,  forc_hgt_u)
      call mg2p_forc%map_aweighted (forc_xy_hgt_q,  forc_hgt_q)
      if (DEF_USE_CBL_HEIGHT) then
	    call mg2p_forc%map_aweighted (forc_xy_hpbl,   forc_hpbl)
      endif

      IF (.not. DEF_USE_Forcing_Downscaling) THEN

         call mg2p_forc%map_aweighted (forc_xy_t    ,  forc_t    )
         call mg2p_forc%map_aweighted (forc_xy_q    ,  forc_q    )
         call mg2p_forc%map_aweighted (forc_xy_prc  ,  forc_prc  )
         call mg2p_forc%map_aweighted (forc_xy_prl  ,  forc_prl  )
         call mg2p_forc%map_aweighted (forc_xy_pbot ,  forc_pbot )
         call mg2p_forc%map_aweighted (forc_xy_frl  ,  forc_frl  )

         if (p_is_worker) then

            do np = 1, numpatch
               IF (DEF_forcing%has_missing_value) THEN
                  IF (.not. forcmask(np)) cycle
               ENDIF

               ! The standard measuring conditions for temperature are two meters above the ground
               ! Scientists have measured the most frigid temperature ever
               ! recorded on the continent's eastern highlands: about (180K) colder than dry ice.
               if(forc_t(np) < 180.) forc_t(np) = 180.
               ! the highest air temp was found in Kuwait 326 K, Sulaibya 2012-07-31;
               ! Pakistan, Sindh 2010-05-26; Iraq, Nasiriyah 2011-08-03
               if(forc_t(np) > 326.) forc_t(np) = 326.

               forc_rhoair(np) = (forc_pbot(np) &
                  - 0.378*forc_q(np)*forc_pbot(np)/(0.622+0.378*forc_q(np)))&
                  / (rgas*forc_t(np))

            end do

         end if

      ELSE

         call mg2p_forc_elm%map_aweighted (forc_xy_t    ,  forc_t_elm    )
         call mg2p_forc_elm%map_aweighted (forc_xy_q    ,  forc_q_elm    )
         call mg2p_forc_elm%map_aweighted (forc_xy_prc  ,  forc_prc_elm  )
         call mg2p_forc_elm%map_aweighted (forc_xy_prl  ,  forc_prl_elm  )
         call mg2p_forc_elm%map_aweighted (forc_xy_pbot ,  forc_pbot_elm )
         call mg2p_forc_elm%map_aweighted (forc_xy_frl  ,  forc_lwrad_elm)
         call mg2p_forc_elm%map_aweighted (forc_xy_hgt_t,  forc_hgt_elm  )

         if (p_is_worker) then

            do ne = 1, numelm
               IF (DEF_forcing%has_missing_value) THEN
                  IF (.not. forcmask_elm(ne)) cycle
               ENDIF

               ! The standard measuring conditions for temperature are two meters above the ground
               ! Scientists have measured the most frigid temperature ever
               ! recorded on the continent's eastern highlands: about (180K) colder than dry ice.
               if(forc_t_elm(ne) < 180.) forc_t_elm(ne) = 180.
               ! the highest air temp was found in Kuwait 326 K, Sulaibya 2012-07-31;
               ! Pakistan, Sindh 2010-05-26; Iraq, Nasiriyah 2011-08-03
               if(forc_t_elm(ne) > 326.) forc_t_elm(ne) = 326.

               forc_rho_elm(ne) = (forc_pbot_elm(ne) &
                  - 0.378*forc_q_elm(ne)*forc_pbot_elm(ne)/(0.622+0.378*forc_q_elm(ne)))&
                  / (rgas*forc_t_elm(ne))

               forc_th_elm(ne) = forc_t_elm(ne) * (1.e5/forc_pbot_elm(ne)) ** (rair/cpair)

            end do

            CALL downscale_forcings ( &
               numelm, numpatch, elm_patch%substt, elm_patch%subend, glacierss, elm_patch%subfrc,   &
               ! forcing in gridcells
            forc_topo_elm, forc_t_elm,   forc_th_elm,  forc_q_elm,     forc_pbot_elm, &
               forc_rho_elm,  forc_prc_elm, forc_prl_elm, forc_lwrad_elm, forc_hgt_elm,  &
               ! forcing in patches
            forc_topo,     forc_t,       forc_th,      forc_q,         forc_pbot,     &
               forc_rhoair,   forc_prc,     forc_prl,     forc_frl)

         end if

      ENDIF

#ifdef RangeCheck
#ifdef USEMPI
      call mpi_barrier (p_comm_glb, p_err)
#endif
      if (p_is_master) write(*,'(/, A20)') 'Checking forcing ...'

      call check_vector_data ('Forcing t     [kelvin]', forc_t    )
      call check_vector_data ('Forcing q     [kg/kg] ', forc_q    )
      call check_vector_data ('Forcing prc   [mm/s]  ', forc_prc  )
      call check_vector_data ('Forcing psrf  [pa]    ', forc_psrf )
      call check_vector_data ('Forcing prl   [mm/s]  ', forc_prl  )
      call check_vector_data ('Forcing sols  [W/m2]  ', forc_sols )
      call check_vector_data ('Forcing soll  [W/m2]  ', forc_soll )
      call check_vector_data ('Forcing solsd [W/m2]  ', forc_solsd)
      call check_vector_data ('Forcing solld [W/m2]  ', forc_solld)
      call check_vector_data ('Forcing frl   [W/m2]  ', forc_frl  )
      if (DEF_USE_CBL_HEIGHT) then
        call check_vector_data ('Forcing hpbl  ', forc_hpbl )
      endif

#ifdef USEMPI
      call mpi_barrier (p_comm_glb, p_err)
#endif
#endif

   END SUBROUTINE read_forcing


   ! ------------------------------------------------------------
   !
   ! !DESCRIPTION:
   !    read lower and upper boundary forcing data, a major interface of this
   !    module
   !
   ! REVISIONS:
   ! Hua Yuan, 04/2014: initial code
   ! ------------------------------------------------------------
   SUBROUTINE metreadLBUB (idate, dir_forcing)

      use MOD_UserSpecifiedForcing
      USE MOD_Namelist
      use MOD_DataType
      use MOD_NetCDFBlock
      use MOD_RangeCheck
      implicit none

      integer, intent(in) :: idate(3)
      character(len=*), intent(in) :: dir_forcing

      ! Local variables
      integer         :: ivar, year, month, day, time_i
      type(timestamp) :: mtstamp
      character(len=256) :: filename

      mtstamp = idate

      do ivar = 1, NVAR

         if (trim(vname(ivar)) == 'NULL') cycle     ! no data, cycle

         ! lower and upper boundary data already exist, cycle
         if ( .NOT.(tstamp_LB(ivar)=='NULL') .AND. .NOT.(tstamp_UB(ivar)=='NULL') .AND. &
            tstamp_LB(ivar)<=mtstamp .AND. mtstamp<=tstamp_UB(ivar) ) then
            cycle
         end if

         ! set lower boundary time stamp and get data
         if (tstamp_LB(ivar) == 'NULL') then
            call setstampLB(mtstamp, ivar, year, month, day, time_i)

            ! read forcing data
            filename = trim(dir_forcing)//trim(metfilename(year, month, day, ivar))
            IF (trim(DEF_forcing%dataset) == 'POINT') THEN
               CALL ncio_read_site_time (filename, vname(ivar), time_i, metdata)
            ELSE
               call ncio_read_block_time (filename, vname(ivar), gforc, time_i, metdata)
            ENDIF

            call block_data_copy (metdata, forcn_LB(ivar))
         end if

         ! set upper boundary time stamp and get data
         if (tstamp_UB(ivar) == 'NULL' .OR. tstamp_UB(ivar) < mtstamp) then
            if ( .NOT. (tstamp_UB(ivar) == 'NULL') ) then
               call block_data_copy (forcn_UB(ivar), forcn_LB(ivar))
            end if
            call setstampUB(ivar, year, month, day, time_i)
            ! when reaching the end of forcing data, always reuse the last time step data
            if (year <= endyr) then
               ! read forcing data
               filename = trim(dir_forcing)//trim(metfilename(year, month, day, ivar))
               IF (trim(DEF_forcing%dataset) == 'POINT') THEN
                  CALL ncio_read_site_time (filename, vname(ivar), time_i, metdata)
               ELSE
                  call ncio_read_block_time (filename, vname(ivar), gforc, time_i, metdata)
               ENDIF

               call block_data_copy (metdata, forcn_UB(ivar))
            else
               write(*,*) year, endyr
               print *, 'NOTE: reaching the end of forcing data, always reuse the last time step data!'
            end if
            if (ivar == 7) then  ! calculate time average coszen, for shortwave radiation
               call calavgcos()
            end if
         end if

      end do

   END SUBROUTINE metreadLBUB


   !-------------------------------------------------
   SUBROUTINE metread_latlon (dir_forcing, idate)

      use MOD_SPMD_Task
      use MOD_NetCDFSerial
      use MOD_UserSpecifiedForcing
      USE MOD_Namelist
      implicit none

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

      IF (trim(DEF_forcing%dataset) == 'POINT') THEN
         CALL gforc%define_by_ndims (360, 180)
      ELSE

         mtstamp = idate

         call setstampLB(mtstamp, 1, year, month, day, time_i)
         filename = trim(dir_forcing)//trim(metfilename(year, month, day, 1))
         tstamp_LB(1) = timestamp(-1, -1, -1)

         if (dim2d) then
            call ncio_read_bcast_serial (filename, latname, latxy)
            call ncio_read_bcast_serial (filename, lonname, lonxy)

            allocate (lat_in (size(latxy,2)))
            allocate (lon_in (size(lonxy,1)))
            lat_in = latxy(1,:)
            lon_in = lonxy(:,1)

            deallocate (latxy)
            deallocate (lonxy)
         else
            call ncio_read_bcast_serial (filename, latname, lat_in)
            call ncio_read_bcast_serial (filename, lonname, lon_in)
         ENDIF

         IF (.not. DEF_forcing%regional) THEN
            call gforc%define_by_center (lat_in, lon_in)
         ELSE
            call gforc%define_by_center (lat_in, lon_in, &
               south = DEF_forcing%regbnd(1), north = DEF_forcing%regbnd(2), &
               west  = DEF_forcing%regbnd(3), east  = DEF_forcing%regbnd(4))
         ENDIF

         deallocate (lat_in)
         deallocate (lon_in)
      ENDIF

      call gforc%set_rlon ()
      call gforc%set_rlat ()

   END SUBROUTINE metread_latlon

   !-------------------------------------------------
   SUBROUTINE metread_time (dir_forcing)

      use MOD_SPMD_Task
      use MOD_NetCDFSerial
      use MOD_UserSpecifiedForcing
      USE MOD_Namelist
      implicit none

      character(len=*), intent(in) :: dir_forcing

      ! Local variables
      character(len=256) :: filename
      character(len=256) :: timeunit, timestr
      REAL(r8), allocatable :: forctime_sec (:)
      INTEGER :: year, month, day, hour, minute, second
      INTEGER :: itime, maxday
      INTEGER*8 :: sec_long

      filename = trim(dir_forcing)//trim(fprefix(1))

      CALL ncio_read_serial (filename, 'time', forctime_sec)
      CALL ncio_get_attr    (filename, 'time', 'units', timeunit)

      timestr = timeunit(15:18) // ' ' // timeunit(20:21) // ' ' // timeunit(23:24) &
         // ' ' // timeunit(26:27) // ' ' // timeunit(29:30) // ' ' // timeunit(32:33)
      read(timestr,*) year, month, day, hour, minute, second

      allocate (forctime (size(forctime_sec)))

      forctime(1)%year = year
      forctime(1)%day  = get_calday(month*100+day, isleapyear(year))
      sec_long = hour*3600 + minute*60 + second + forctime_sec(1)

      DO itime = 1, size(forctime)
         IF (itime > 1) THEN
            forctime(itime) = forctime(itime-1)
            sec_long = sec_long + forctime_sec(itime) - forctime_sec(itime-1)
         ENDIF

         DO WHILE (sec_long > 86400)
            sec_long = sec_long - 86400
            IF( isleapyear(forctime(itime)%year) ) THEN
               maxday = 366
            ELSE
               maxday = 365
            ENDIF
            forctime(itime)%day = forctime(itime)%day + 1
            IF(forctime(itime)%day > maxday) THEN
               forctime(itime)%year = forctime(itime)%year + 1
               forctime(itime)%day = 1
            ENDIF
         ENDDO

         forctime(itime)%sec = sec_long
      ENDDO

   END SUBROUTINE metread_time

   ! ------------------------------------------------------------
   !
   ! !DESCRIPTION:
   !    set the lower boundary time stamp and record information,
   !    a KEY function of this module
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

      implicit none
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
            stop
         ELSE
            DO WHILE (.not. (mtstamp < forctime(time_i+1)))
               time_i = time_i + 1
            ENDDO
            iforctime(var_i) = time_i
            tstamp_LB(var_i) = forctime(iforctime(var_i))
         ENDIF
         write(*,*) mtstamp, forctime(time_i)

         RETURN
      ENDIF

      tstamp_LB(var_i)%year = year
      tstamp_LB(var_i)%day  = day

      ! in the case of one year one file
      if ( trim(groupby) == 'year' ) then

         ! calculate the intitial second
         sec    = 86400*(day-1) + sec
         time_i = floor( (sec-offset(var_i)-0.01) *1. / dtime(var_i) ) + 1
         sec    = (time_i-1)*dtime(var_i) + offset(var_i) - 86400*(day-1)
         tstamp_LB(var_i)%sec = sec

         ! set time stamp (ststamp_LB)
         if (sec <= 0) then
            tstamp_LB(var_i)%sec = 86400 + sec
            tstamp_LB(var_i)%day = day - 1
            if (tstamp_LB(var_i)%day == 0) then
               tstamp_LB(var_i)%year = year - 1
               if ( isleapyear(tstamp_LB(var_i)%year) ) then
                  tstamp_LB(var_i)%day = 366
               else
                  tstamp_LB(var_i)%day = 365
               end if
            end if
         end if

         ! set record info (year, time_i)
         if ( sec<0 .OR. (sec==0 .AND. offset(var_i).NE.0) ) then

            ! if the required dada just behind the first record
            ! -> set to the first record
            if ( year==startyr .AND. month==startmo .AND. day==1 ) then
               sec = offset(var_i)

               ! else, set to one record backward
            else
               sec = 86400 + sec
               day = day - 1
               if (day == 0) then
                  year = year - 1
                  if ( isleapyear(year) .AND. leapyear) then
                     day = 366
                  else
                     day = 365
                  end if
               end if
            end if
         end if ! end if (sec <= 0)

         ! in case of leapyear with a non-leayyear calendar
         ! use the data 1 day before after FEB 28th (Julian day 59).
         if ( .NOT. leapyear .AND. isleapyear(year) .AND. day>59 ) then
            day = day - 1
         end if

         ! get record time index
         sec = 86400*(day-1) + sec
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
      end if

      ! in the case of one month one file
      if ( trim(groupby) == 'month' ) then

         if ( isleapyear(year) ) then
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         else
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         end if

         ! calculate initial month and day values
         call julian2monthday(year, day, month, mday)

         ! calculate initial second value
         sec    = 86400*(mday-1) + sec
         time_i = floor( (sec-offset(var_i)-0.01) *1. / dtime(var_i) ) + 1
         sec    = (time_i-1)*dtime(var_i) + offset(var_i) - 86400*(mday-1)
         tstamp_LB(var_i)%sec  = sec

         ! set time stamp (ststamp_LB)
         if (sec <= 0) then
            tstamp_LB(var_i)%sec = 86400 + sec
            tstamp_LB(var_i)%day = day - 1
            if (tstamp_LB(var_i)%day == 0) then
               tstamp_LB(var_i)%year = year - 1
               if ( isleapyear(tstamp_LB(var_i)%year) ) then
                  tstamp_LB(var_i)%day = 366
               else
                  tstamp_LB(var_i)%day = 365
               end if
            end if
         end if

         ! set record info (year, month, time_i)
         if ( sec<0 .OR. (sec==0 .AND. offset(var_i).NE.0) ) then

            ! if just behind the first record -> set to first record
            if ( year==startyr .AND. month==startmo .AND. mday==1 ) then
               sec = offset(var_i)

               ! set to one record backward
            else
               sec = 86400 + sec
               mday = mday - 1
               if (mday == 0) then
                  month = month - 1
                  ! bug found by Zhu Siguang & Zhang Xiangxiang, 05/19/2014
                  ! move the below line in the 'else' statement
                  !mday = months(month) - months(month-1)
                  if (month == 0) then
                     month = 12
                     year = year - 1
                     mday = 31
                  else
                     mday = months(month) - months(month-1)
                  end if
               end if
            end if
         end if

         ! in case of leapyear with a non-leayyear calendar
         ! use the data 1 day before, i.e., FEB 28th.
         if ( .NOT. leapyear .AND. isleapyear(year) .AND. month==2 .AND. mday==29 ) then
            mday = 28
         end if

         ! get record time index
         sec = 86400*(mday-1) + sec
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
      end if

      ! in the case of one day one file
      if ( trim(groupby) == 'day' ) then

         ! calculate initial month and day values
         call julian2monthday(year, day, month, mday)

         ! calculate initial second value
         time_i = floor( (sec-offset(var_i)-0.01) *1. / dtime(var_i) ) + 1
         sec    = (time_i-1)*dtime(var_i) + offset(var_i)
         tstamp_LB(var_i)%sec  = sec

         ! set time stamp (ststamp_LB)
         if (sec <= 0) then
            tstamp_LB(var_i)%sec = 86400 + sec
            tstamp_LB(var_i)%day = day - 1
            if (tstamp_LB(var_i)%day == 0) then
               tstamp_LB(var_i)%year = year - 1
               if ( isleapyear(tstamp_LB(var_i)%year) ) then
                  tstamp_LB(var_i)%day = 366
               else
                  tstamp_LB(var_i)%day = 365
               end if
            end if

            if ( year==startyr .AND. month==startmo .AND. mday==1 ) then
                  sec = offset(var_i)
                  ! set to one record backward
               else
                  sec = 86400 + sec
                  year = tstamp_LB(var_i)%year
                  call julian2monthday(tstamp_LB(var_i)%year, tstamp_LB(var_i)%day, month, mday)
               end if
            end if

            ! in case of leapyear with a non-leayyear calendar
            ! use the data 1 day before, i.e., FEB 28th.
            if ( .NOT. leapyear .AND. isleapyear(year) .AND. month==2 .AND. mday==29 ) then
               mday = 28
            end if

            ! get record time index
            time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
         end if

         if (time_i <= 0) then
            write(6, *) "got the wrong time record of forcing! stop!"; stop
         end if

         return

      END SUBROUTINE setstampLB

      ! ------------------------------------------------------------
      !
      ! !DESCRIPTION:
      !    set the upper boundary time stamp and record information,
      !    a KEY function of this module
      !
      ! REVISIONS:
      ! Hua Yuan, 04/2014: initial code
      ! ------------------------------------------------------------
      SUBROUTINE setstampUB(var_i, year, month, mday, time_i)

         implicit none
         integer,         intent(in)  :: var_i
         integer,         intent(out) :: year
         integer,         intent(out) :: month
         integer,         intent(out) :: mday
         integer,         intent(out) :: time_i

         integer :: day, sec
         integer :: months(0:12)

      IF (trim(DEF_forcing%dataset) == 'POINT') THEN
         if ( tstamp_UB(var_i) == 'NULL' ) then
            tstamp_UB(var_i) = forctime(iforctime(var_i)+1)
         ELSE
            iforctime(var_i) = iforctime(var_i) + 1
            tstamp_LB(var_i) = forctime(iforctime(var_i))
            tstamp_UB(var_i) = forctime(iforctime(var_i) + 1)
         ENDIF

         time_i = iforctime(var_i)
         year = tstamp_UB(var_i)%year
         RETURN
      ENDIF

      ! calculate the time stamp
      if ( tstamp_UB(var_i) == 'NULL' ) then
         tstamp_UB(var_i) = tstamp_LB(var_i) + dtime(var_i)
      else
         tstamp_LB(var_i) = tstamp_UB(var_i)
         tstamp_UB(var_i) = tstamp_UB(var_i) + dtime(var_i)
      end if

         ! calcualte initial year, day, and second values
         year = tstamp_UB(var_i)%year
         day  = tstamp_UB(var_i)%day
         sec  = tstamp_UB(var_i)%sec

         if ( trim(groupby) == 'year' ) then

            ! adjust year value
            if ( sec==86400 .AND. offset(var_i).EQ.0 ) then
               sec = 0
               day = day + 1
               if( isleapyear(year) .AND. day==367) then
                  year = year + 1; day = 1
               end if
               if( .NOT. isleapyear(year) .AND. day==366) then
                  year = year + 1; day = 1
               end if
            end if

            ! in case of leapyear with a non-leayyear calendar
            ! use the data 1 day before after FEB 28th (Julian day 59).
            if ( .NOT. leapyear .AND. isleapyear(year) .AND. day>59 ) then
               day = day - 1
            end if

            ! set record index
            sec = 86400*(day-1) + sec
            time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
         end if

         if ( trim(groupby) == 'month' ) then

            if ( isleapyear(year) ) then
               months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
            else
               months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
            end if

            ! calculate initial month and day values
            call julian2monthday(year, day, month, mday)

            ! record in the next day, adjust year, month and second values
            if ( sec==86400 .AND. offset(var_i).EQ.0 ) then
               sec  = 0
               mday = mday + 1
               if ( mday > (months(month)-months(month-1)) ) then
                  mday = 1
                  ! bug found by Zhu Siguang, 05/25/2014
                  ! move the below line in the 'else' statement
                  !month = month + 1
                  if (month == 12) then
                     month = 1
                     year = year + 1
                  else
                     month = month + 1
                  end if
               end if
            end if

            ! in case of leapyear with a non-leayyear calendar
            ! for day 29th Feb, use the data 1 day before, i.e., 28th FEB.
            if ( .NOT. leapyear .AND. isleapyear(year) .AND. month==2 .AND. mday==29 ) then
               mday = 28
            end if

            ! set record index
            sec    = 86400*(mday-1) + sec
            time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
         end if

         if ( trim(groupby) == 'day' ) then
            if ( isleapyear(year) ) then
               months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
            else
               months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
            end if

            ! calculate initial month and day values
            call julian2monthday(year, day, month, mday)
            !mday = day

            ! record in the next day, adjust year, month and second values
            if ( sec==86400 .AND. offset(var_i).EQ.0 ) then
               sec  = 0
               mday = mday + 1
               if ( mday > (months(month)-months(month-1)) ) then
                  mday = 1
                  ! bug found by Zhu Siguang, 05/25/2014
                  ! move the below line in the 'else' statement
                  !month = month + 1
                  if (month == 12) then
                     month = 1
                     year = year + 1
                  else
                     month = month + 1
                  end if
               end if
            end if

            ! in case of leapyear with a non-leayyear calendar
            ! for day 29th Feb, use the data 1 day before, i.e., 28th FEB.
            if ( .NOT. leapyear .AND. isleapyear(year) .AND. month==2 .AND. mday==29 ) then
               mday = 28
            end if

            ! set record index
            time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
         end if

         if (time_i < 0) then
            write(6, *) "got the wrong time record of forcing! stop!"; stop
         end if

         return

      END SUBROUTINE setstampUB

      ! ------------------------------------------------------------
      ! !DESCRIPTION:
      ! calculate time average coszen value bwteeen [LB, UB]
      !
      ! REVISIONS:
      ! 04/2014, yuan: this method is adapted from CLM
      ! ------------------------------------------------------------
      SUBROUTINE calavgcos()

         use MOD_Block
         use MOD_DataType
         implicit none

         integer  :: iblkme, ib, jb, i, j, ilon, ilat
         real(r8) :: calday, cosz
         type(timestamp) :: tstamp

         tstamp = tstamp_LB(7)
         call flush_block_data (avgcos, 0._r8)

         do while (tstamp < tstamp_UB(7))

            tstamp = tstamp + deltim_int

            DO iblkme = 1, gblock%nblkme
               ib = gblock%xblkme(iblkme)
               jb = gblock%yblkme(iblkme)
               do j = 1, gforc%ycnt(jb)
                  do i = 1, gforc%xcnt(ib)

                     ilat = gforc%ydsp(jb) + j
                     ilon = gforc%xdsp(ib) + i
                     if (ilon > gforc%nlon) ilon = ilon - gforc%nlon

                     calday = calendarday(tstamp)
                     cosz = orb_coszen(calday, gforc%rlon(ilon), gforc%rlat(ilat))
                     cosz = max(0.001, cosz)
                     avgcos%blk(ib,jb)%val(i,j) = avgcos%blk(ib,jb)%val(i,j) &
                        + cosz*deltim_real /real(tstamp_UB(7)-tstamp_LB(7))

                  end do
               end do
            end do
         end do

      END SUBROUTINE calavgcos

   end module MOD_Forcing
