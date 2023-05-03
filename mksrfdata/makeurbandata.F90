#include <define.h>

! ======================================================
! Aggreate/screen high-resolution urban dataset
! to a lower resolutioin/subset data, suitable for running
! regional or point cases.
! ======================================================

SUBROUTINE makeurbandata (dir_rawdata, dir_srfdata, lc_year, &
      grid_urban_1km, grid_urban_5km, grid_urban_100m, grid_urban_500m)

   USE precision
   USE mod_namelist
   USE spmd_task
   USE mod_grid
   USE mod_landpatch
   USE ncio_serial
   USE ncio_vector
   USE ncio_block
   USE mod_colm_debug
   USE mod_aggregation
#ifdef SinglePoint
   USE mod_single_srfdata
#endif

   IMPLICIT NONE

   CHARACTER(len=256), intent(in) :: dir_rawdata
   CHARACTER(len=256), intent(in) :: dir_srfdata

   INTEGER , intent(in) :: lc_year

   TYPE(grid_type), intent(in) :: grid_urban_1km
   TYPE(grid_type), intent(in) :: grid_urban_5km
   TYPE(grid_type), intent(in) :: grid_urban_100m
   TYPE(grid_type), intent(in) :: grid_urban_500m

   ! Local Variables
   CHARACTER(len=*), parameter :: DATASRC = "MOD"
   CHARACTER(len=*), parameter :: Title   = "Land surface model input urban data"
   CHARACTER(len=*), parameter :: Authors = "Dai YJ group at Sun Yat-sen University"
   CHARACTER(len=*), parameter :: Address = "School of Atmospheric Sciences, Sun Yat-sen University, Zhuhai, China"
   CHARACTER(len=*), parameter :: Email   = "yuanh25@mail.sysu.edu.cn"

   INTEGER, parameter :: rid  = 33
   INTEGER, parameter :: ns   = 2
   INTEGER, parameter :: nr   = 2
   INTEGER, parameter :: ulev = 10

#ifndef USE_LCZ
   INTEGER , ALLOCATABLE, DIMENSION(:,:) :: urrgid
#endif

   ! input variables
   TYPE(block_data_int32_2d) :: LUCY_reg
   TYPE(block_data_real8_2d) :: pop
   TYPE(block_data_real8_2d) :: gfcc_tc
   TYPE(block_data_real8_2d) :: gedi_th
   TYPE(block_data_real8_2d) :: gl30_wt
   TYPE(block_data_real8_2d) :: wtrf
   TYPE(block_data_real8_2d) :: htrf
   TYPE(block_data_real8_2d) :: hlai
   TYPE(block_data_real8_2d) :: hsai

   ! output variables
   INTEGER , ALLOCATABLE, DIMENSION(:) :: LUCY_coun
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: pop_ur
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: pct_tc
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: htop_ur
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: pct_urwt
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: wt_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: ht_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: ur_lai
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: ur_sai

   REAL(r8), ALLOCATABLE, DIMENSION(:,:)   :: area
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: wgt_top
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: tc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: urwt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: htop

   INTEGER , allocatable, dimension(:) :: LUCY_reg_one
   REAL(r8), allocatable, dimension(:) :: area_one
   REAL(r8), allocatable, dimension(:) :: pop_one
   REAL(r8), allocatable, dimension(:) :: gfcc_tc_one
   REAL(r8), allocatable, dimension(:) :: gedi_th_one
   REAL(r8), allocatable, dimension(:) :: gl30_wt_one
   REAL(r8), allocatable, dimension(:) :: wt_rf_one
   REAL(r8), allocatable, dimension(:) :: ht_rf_one
   REAL(r8), allocatable, dimension(:) :: hlai_one
   REAL(r8), allocatable, dimension(:) :: slai_one
   INTEGER , allocatable, dimension(:) :: urrgid_one

#ifndef USE_LCZ
   REAL(r8), allocatable, DIMENSION(:,:)  :: hwrcan, wtrd, emrf, emwl, ncar_wt
   REAL(r8), allocatable, DIMENSION(:,:)  :: emimrd, emperd, ncar_ht
   REAL(r8), allocatable, DIMENSION(:,:)  :: thrf, thwl, tbmin, tbmax

   REAL(r8), allocatable, DIMENSION(:,:,:)  :: cvrf, cvwl, cvimrd, &
                                               tkrf, tkwl, tkimrd
   REAL(r8), allocatable, DIMENSION(:,:,:,:):: albrf, albwl, albimrd, albperd

   REAL(r8), ALLOCATABLE, DIMENSION(:) :: hwr_can
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: wt_rd
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: em_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: em_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: em_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: em_perd
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: th_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: th_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: tb_min
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: tb_max

   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: cv_wgt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: tk_wgt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: cv_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: cv_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: cv_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: tk_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: tk_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: tk_imrd

   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: alb_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: alb_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: alb_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: alb_perd
#endif

   CHARACTER(len=256) landsrfdir, landname, suffix
   CHARACTER(len=4)   cyear, c5year, cmonth

   write(cyear,'(i4.4)') lc_year
   landsrfdir = trim(dir_srfdata) // '/urban/' // trim(cyear)

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A)') 'Making urban data ('//trim(cyear)//') ...'
      CALL system('mkdir -p ' // trim(adjustl(landsrfdir)))
   ENDIF
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   write(c5year, '(i4.4)') (lc_year/5)*5

   ! ******* LUCY_id *******
   IF (p_is_io) THEN

      landname = TRIM(dir_rawdata)//'urban/LUCY_countryid.nc'
      CALL allocate_block_data (grid_urban_5km, LUCY_reg)
      CALL ncio_read_block (landname, 'Country_id', grid_urban_5km, LUCY_reg)

#ifdef USEMPI
      CALL aggregation_data_daemon (grid_urban_5km, data_i4_2d_in1 = LUCY_reg)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate ( LUCY_coun (numurban))

      DO iurban = 1, numpatch
         CALL aggregation_request_data (landurban, iurban, grid_urban_5km, &
            data_i4_2d_in1 = LUCY_reg, data_i4_2d_out1 = LUCY_reg_one)
         LUCY_coun(iurban) = num_max_frequency (LUCY_reg_one)
      ENDDO
#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

   landname = trim(landsrfdir)//'/urban/LUCY_country_id.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'LUCY_id', 'urban', landurban, LUCY_coun, 1)

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   ! ******* POP_DEN *******
   IF (p_is_io) THEN

      CALL allocate_block_data (grid_urban_1km, pop)

      landdir = TRIM(dir_rawdata)//'/urban/'
      suffix  = 'MOD'//trim(c5year)
      CALL read_5x5_data (landdir, suffix, grid_urban_1km, "POP", pop)

#ifdef USEMPI
      CALL aggregation_data_daemon (grid_urban_1km, data_r8_2d_in1 = pop)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate (pop_ur (numurban))

      DO iurban = 1, numurban
         CALL aggregation_request_data (landurban, iurban, grid_urban_1km, area = area_one, &
            data_r8_2d_in1 = pop, data_r8_2d_out1 = pop_one)
         pop_ur(iurban) = sum(pop_one * area_one) / sum(area_one)
      ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/urban.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'POP_DEN', 'urban', landurban, pop_ur, 1)

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   ! ******* Tree : PCT_Tree, HTOP *******
   IF (p_is_io) THEN

      landdir = TRIM(dir_rawdata)//'/urban/'
      suffix  = 'MOD'//trim(c5year)

      CALL allocate_block_data (grid_urban_500m, gfcc_tc)
      CALL read_5x5_data (landdir, suffix, grid_urban_500m, "PCT_Tree", gfcc_tc)

      CALL allocate_block_data (grid_urban_500m, gedi_th)
      CALL read_5x5_data (landdir, suffix, grid_urban_500m, "ETH_htop", gedi_th)

#ifdef USEMPI
      CALL aggregation_data_daemon (grid_urban_500m, &
         data_r8_2d_in1 = gfcc_tc, data_r8_2d_in2 = gedi_th)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate (pct_tc (numurban))
      allocate (htop_ur(numurban))

      DO iurban = 1, numurban
         CALL aggregation_request_data (landurban, iurban, grid_urban_500m, area = area_one, &
            data_r8_2d_in1 = gfcc_tc, data_r8_2d_out1 = gfcc_tc_one, &
            data_r8_2d_in2 = gedi_th, data_r8_2d_out2 = gedi_th_one)
         pct_tc (iurban) = sum(gfcc_tc_one * area_one) / sum(area_one)
         htop_ur(iurban) = sum(gedi_th_one * gfcc_tc_one * area_one) / sum(gfcc_tc_one * area_one)
      ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/PCT_Tree.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'PCT_Tree', 'urban', landurban, pct_tc, 1)

   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/htop_ur.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'URBAN_TREE_TOP', 'urban', landurban, htop_ur, 1)

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   ! ******* PCT_Water *******
   IF (p_is_io) THEN

      CALL allocate_block_data (grid_urban_500m, gl30_wt)

      landdir = TRIM(dir_rawdata)//'/urban/'
      suffix  = 'MOD'//trim(c5year)
      CALL read_5x5_data (landdir, suffix, grid_urban_500m, "PCT_Water", gl30_wt)

#ifdef USEMPI
      CALL aggregation_data_daemon (grid_urban_500m, gl30_wt)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate (pct_urwt (numurban))

      DO iurban = 1, numurban
         CALL aggregation_request_data (landurban, iurban, grid_urban_500m, area = area_one, &
            data_r8_2d_in1 = gl30_wt, data_r8_2d_out1 = gl30_wt_one)
         pct_urwt(iurban) = sum(gl30_wt_one * area_one) / sum(area_one)
      ENDDO

#ifdef USEMPI
      CALL aggregation_lc_worker_done ()
#endif
   ENDIF

   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/PCT_Water.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'PCT_Water', 'urban', landurban, pct_urwt, 1)

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   ! ******* Building : Weight, HTOP_Roof *******
   IF (p_is_io) THEN

      CALL allocate_block_data (grid_urban_500m, wtrf)
      CALL allocate_block_data (grid_urban_500m, htrf)

      landdir = TRIM(dir_rawdata)//'/urban/'
      suffix  = 'MOD'//trim(c5year)
      CALL read_5x5_data (landdir, suffix, grid_urban_500m, "WTLUNIT_ROOF", wtrf)

      landdir = TRIM(dir_rawdata)//'/urban/'
      suffix  = 'MOD'//trim(c5year)
      CALL read_5x5_data (landdir, suffix, grid_urban_500m, "HT_ROOF", htrf)

#ifdef USEMPI
      CALL aggregation_data_daemon (grid_urban_500m, &
         data_r8_2d_in1 = wtrf, data_r8_2d_in2 = htrf)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate (wt_rf (numurban))
      allocate (ht_rf (numurban))

      DO iurban = 1, numurban
         CALL aggregation_request_data (landurban, iurban, grid_urban_500m, area = area_one, &
            data_r8_2d_in1 = wtrf, data_r8_2d_out1 = wt_rf_one, &
            data_r8_2d_in2 = htrf, data_r8_2d_out2 = ht_rf_one)
         wt_rf(iurban) = sum(wt_rf_one * area_one) / sum(area_one)
         ht_rf(iurban) = sum(hr_rf_one * wt_rf_one * area_one) / sum(wt_rf_one * area_one)
      ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/WT_ROOF.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'WT_ROOF', 'urban', landurban, wt_rf, 1)

   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/HT_ROOF.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'HT_ROOF', 'urban', landurban, ht_rf, 1)

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   ! ******* LAI, SAI *******

   landdir = TRIM(dir_rawdata)//'/urban/'
   suffix  = 'MOD'//trim(c5year)

   CALL allocate_block_data (grid_urban_500m, hlai)
   CALL allocate_block_data (grid_urban_500m, slai)
   allocate (ur_lai (numurban))
   allocate (ur_sai (numurban))

   DO imonth = 1, 12

      write(cmonth, '(i2.2)') imonth

      IF (p_is_io) THEN

         CALL read_5x5_data_time (landdir, suffix, grid_urban_500m, "URBAN_TREE_LAI", imonth, hlai)
         CALL read_5x5_data_time (landdir, suffix, grid_urban_500m, "URBAN_TREE_SAI", imonth, slai)

#ifdef USEMPI
         CALL aggregation_data_daemon (grid, &
            data_r8_2d_in1 = gfcc_tc, data_r8_2d_in2 = hlai, data_r8_2d_in3 = slai)
#endif
      ENDIF

      IF (p_is_worker) THEN

         DO iurban = 1, numurban
            CALL aggregation_request_data (landurban, iurban, grid, area = area_one, &
               data_r8_2d_in1 = gfcc_tc, data_r8_2d_out1 = gfcc_tc_one, &
               data_r8_2d_in2 = hlai   , data_r8_2d_out2 = hlai_one   , &
               data_r8_2d_in3 = slai   , data_r8_2d_out3 = slai_one   )
            ur_lai(iurban) = sum(hlai_one * gfcc_tc_one * area_one) / sum(gfcc_tc_one * area_one)
            ur_sai(iurban) = sum(slai_one * gfcc_tc_one * area_one) / sum(gfcc_tc_one * area_one)
         ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
      ENDIF

      landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'urban.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'TREE_LAI'//'_'//trim(cmonth), 'urban', landurban, ur_lai, 1)

      landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'urban.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'TREE_SAI'//'_'//trim(cmonth), 'urban', landurban, ur_sai, 1)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
   ENDDO

#ifndef USE_LCZ
   ! look up table of NCAR urban properties 表数据
   landname = TRIM(dir_rawdata)//'urban/urban_properties.nc'

   CALL ncio_read_bcast_serial (landname,  "CANYON_HWR"         , hwrcan  )
   CALL ncio_read_bcast_serial (landname,  "WTLUNIT_ROOF"       , ncar_wt )
   CALL ncio_read_bcast_serial (landname,  "WTROAD_PERV"        , wtrd    )
   CALL ncio_read_bcast_serial (landname,  "EM_ROOF"            , emrf    )
   CALL ncio_read_bcast_serial (landname,  "EM_WALL"            , emwl    )
   CALL ncio_read_bcast_serial (landname,  "EM_IMPROAD"         , emimrd  )
   CALL ncio_read_bcast_serial (landname,  "EM_PERROAD"         , emperd  )
   CALL ncio_read_bcast_serial (landname,  "ALB_ROOF"           , albrf   )
   CALL ncio_read_bcast_serial (landname,  "ALB_WALL"           , albwl   )
   CALL ncio_read_bcast_serial (landname,  "ALB_IMPROAD"        , albimrd )
   CALL ncio_read_bcast_serial (landname,  "ALB_PERROAD"        , albperd )
   CALL ncio_read_bcast_serial (landname,  "HT_ROOF"            , ncar_ht )
   CALL ncio_read_bcast_serial (landname,  "TK_ROOF"            , tkrf    )
   CALL ncio_read_bcast_serial (landname,  "TK_WALL"            , tkwl    )
   CALL ncio_read_bcast_serial (landname,  "TK_IMPROAD"         , tkimrd  )
   CALL ncio_read_bcast_serial (landname,  "CV_ROOF"            , cvrf    )
   CALL ncio_read_bcast_serial (landname,  "CV_WALL"            , cvwl    )
   CALL ncio_read_bcast_serial (landname,  "CV_IMPROAD"         , cvimrd  )
   CALL ncio_read_bcast_serial (landname,  "THICK_ROOF"         , thrf    )
   CALL ncio_read_bcast_serial (landname,  "THICK_WALL"         , thwl    )
   CALL ncio_read_bcast_serial (landname,  "T_BUILDING_MIN"     , tbmin   )
   CALL ncio_read_bcast_serial (landname,  "T_BUILDING_MAX"     , tbmax   )

   IF (p_is_io) THEN

      CALL allocate_block_data (grid, urrgid)

      landdir = TRIM(dir_rawdata)//'urban_5x5/'
      suffix  = 'MOD'//trim(c5year)
      CALL read_5x5_data (landdir, suffix, grid_urban_500m, "REGION_ID", urrgid)

#ifdef USEMPI
      CALL aggregation_data_daemon (grid, data_i4_2d_in1 = urrgid)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate ( hwr_can  (numurban))
      allocate ( wt_rd    (numurban))
      allocate ( em_rf    (numurban))
      allocate ( em_wl    (numurban))
      allocate ( em_imrd  (numurban))
      allocate ( em_perd  (numurban))
      allocate ( th_rf    (numurban))
      allocate ( th_wl    (numurban))
      allocate ( tb_min   (numurban))
      allocate ( tb_max   (numurban))
      allocate ( tk_wgt   (ulev, numurban))
      allocate ( cv_wgt   (ulev, numurban))
      allocate ( cv_rf    (ulev, numurban))
      allocate ( cv_wl    (ulev, numurban))
      allocate ( cv_imrd  (ulev, numurban))
      allocate ( tk_rf    (ulev, numurban))
      allocate ( tk_wl    (ulev, numurban))
      allocate ( tk_imrd  (ulev, numurban))
      allocate ( alb_rf   (nr, ns, numurban))
      allocate ( alb_wl   (nr, ns, numurban))
      allocate ( alb_imrd (nr, ns, numurban))
      allocate ( alb_perd (nr, ns, numurban))

      hwr_can  (:) = 0.
      wt_rd    (:) = 0.
      em_rf    (:) = 0.
      em_wl    (:) = 0.
      em_imrd  (:) = 0.
      em_perd  (:) = 0.
      th_rf    (:) = 0.
      th_wl    (:) = 0.
      tb_min   (:) = 0.
      tb_max   (:) = 0.
      tk_wgt   (:,:) = 0.
      cv_wgt   (:,:) = 0.
      cv_rf    (:,:) = 0.
      cv_wl    (:,:) = 0.
      cv_imrd  (:,:) = 0.
      tk_rf    (:,:) = 0.
      tk_wl    (:,:) = 0.
      tk_imrd  (:,:) = 0.
      alb_rf   (:,:,:) = 0.
      alb_wl   (:,:,:) = 0.
      alb_imrd (:,:,:) = 0.
      alb_perd (:,:,:) = 0.

      DO iurban = 1, numurban
         CALL aggregation_request_data (landurban, iurban, grid, area = area_one, &
            data_i4_2d_in1 = urrgid, data_i4_2d_out1 = urrgid_one)

         inx = landurban%settyp(iurban)

         ipxstt = landurban%ipxstt(iurban)
         ipxend = landurban%ipxend(iurban)

         sumarea = sum(area_one)

         DO ipxl = ipxstt, ipxend

            uxid = urrgid_one(ipxl)

            hwr_can (iurban) = hwr_can (iurban) + hwrcan (inx,uxid) * area_one(ipxl)
            wt_rd   (iurban) = wt_rd   (iurban) + wtrd   (inx,uxid) * area_one(ipxl)
            em_rf   (iurban) = em_rf   (iurban) + emrf   (inx,uxid) * area_one(ipxl)
            em_wl   (iurban) = em_wl   (iurban) + emwl   (inx,uxid) * area_one(ipxl)
            em_imrd (iurban) = em_imrd (iurban) + emimrd (inx,uxid) * area_one(ipxl)
            em_perd (iurban) = em_perd (iurban) + emperd (inx,uxid) * area_one(ipxl)
            th_rf   (iurban) = th_rf   (iurban) + thrf   (inx,uxid) * area_one(ipxl)
            th_wl   (iurban) = th_wl   (iurban) + thwl   (inx,uxid) * area_one(ipxl)
            tb_min  (iurban) = tb_min  (iurban) + tbmin  (inx,uxid) * area_one(ipxl)
            tb_max  (iurban) = tb_max  (iurban) + tbmax  (inx,uxid) * area_one(ipxl)

            DO m = 1, 10
               IF (tkimrd(inx,uxid,m) .ne. -999.) THEN
                  tk_imrd(m,iurban) = tk_imrd(m,iurban) + tkimrd(inx,uxid,m) * area_one(ipxl)
                  tk_wgt (m,iurban) = tk_wgt (m,iurban) + area_one(ipxl)
               ENDIF
               IF (cvimrd(inx,uxid,m) .ne. -999.) THEN
                  cv_imrd(m,iurban) = cv_imrd(m,iurban) + cvimrd(inx,uxid,m) * area_one(ipxl)
                  cv_wgt (m,iurban) = cv_wgt (m,iurban) + area_one(ipxl)
               ENDIF
            ENDDO

            cv_rf (:,iurban) = cv_rf (:,iurban) + cvrf (inx,uxid,:) * area_one(ipxl)
            cv_wl (:,iurban) = cv_wl (:,iurban) + cvwl (inx,uxid,:) * area_one(ipxl)
            tk_rf (:,iurban) = tk_rf (:,iurban) + tkrf (inx,uxid,:) * area_one(ipxl)
            tk_wl (:,iurban) = tk_wl (:,iurban) + tkwl (inx,uxid,:) * area_one(ipxl)

            alb_rf  (:,:,iurban) = alb_rf  (:,:,iurban) + albrf  (inx,uxid,:,:) * area_one(ipxl)
            alb_wl  (:,:,iurban) = alb_wl  (:,:,iurban) + albwl  (inx,uxid,:,:) * area_one(ipxl)
            alb_imrd(:,:,iurban) = alb_imrd(:,:,iurban) + albimrd(inx,uxid,:,:) * area_one(ipxl)
            alb_perd(:,:,iurban) = alb_perd(:,:,iurban) + albperd(inx,uxid,:,:) * area_one(ipxl)

         ENDDO

         hwr_can (iurban) = hwr_can (iurban) / sumarea
         wt_rd   (iurban) = wt_rd   (iurban) / sumarea
         em_rf   (iurban) = em_rf   (iurban) / sumarea
         em_wl   (iurban) = em_wl   (iurban) / sumarea
         em_imrd (iurban) = em_imrd (iurban) / sumarea
         em_perd (iurban) = em_perd (iurban) / sumarea
         th_rf   (iurban) = th_rf   (iurban) / sumarea
         th_wl   (iurban) = th_wl   (iurban) / sumarea
         tb_min  (iurban) = tb_min  (iurban) / sumarea
         tb_max  (iurban) = tb_max  (iurban) / sumarea

         cv_rf (:,iurban) = cv_rf (:,iurban) / sumarea
         cv_wl (:,iurban) = cv_wl (:,iurban) / sumarea
         tk_rf (:,iurban) = tk_rf (:,iurban) / sumarea
         tk_wl (:,iurban) = tk_wl (:,iurban) / sumarea

         DO m = 1, 10
            tk_imrd(m,iurban) = tk_imrd(m,iurban) / tk_wgt(m,iurban)
            cv_imrd(m,iurban) = cv_imrd(m,iurban) / cv_wgt(m,iurban)
         ENDDO

         alb_rf  (:,:,iurban) = alb_rf  (:,:,iurban) / sumarea
         alb_wl  (:,:,iurban) = alb_wl  (:,:,iurban) / sumarea
         alb_imrd(:,:,iurban) = alb_imrd(:,:,iurban) / sumarea
         alb_perd(:,:,iurban) = alb_perd(:,:,iurban) / sumarea

      ENDDO
#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/urban.nc'
   CALL ncio_create_file_vector (landname, landurban)

   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_define_dimension_vector (landname, landurban, 'numsolar', ns)
   CALL ncio_define_dimension_vector (landname, landurban, 'numrad'  , nr)
   CALL ncio_define_dimension_vector (landname, landurban, 'ulev'    , ulev)

   CALL ncio_write_vector (landname, 'CANYON_HWR'    , 'urban', landurban, hwr_can, 1)
   CALL ncio_write_vector (landname, 'WTROAD_PERV'   , 'urban', landurban, wt_rd  , 1)
   CALL ncio_write_vector (landname, 'EM_ROOF'       , 'urban', landurban, em_rf  , 1)
   CALL ncio_write_vector (landname, 'EM_WALL'       , 'urban', landurban, em_wl  , 1)
   CALL ncio_write_vector (landname, 'EM_IMPROAD'    , 'urban', landurban, em_imrd, 1)
   CALL ncio_write_vector (landname, 'EM_PERROAD'    , 'urban', landurban, em_perd, 1)
   CALL ncio_write_vector (landname, 'THICK_ROOF'    , 'urban', landurban, th_rf  , 1)
   CALL ncio_write_vector (landname, 'THICK_WALL'    , 'urban', landurban, th_wl  , 1)
   CALL ncio_write_vector (landname, 'T_BUILDING_MIN', 'urban', landurban, tb_min , 1)
   CALL ncio_write_vector (landname, 'T_BUILDING_MAX', 'urban', landurban, tb_max , 1)

   CALL ncio_write_vector (landname, 'CV_ROOF'   , 'ulev', ulev, 'urban', landurban, cv_rf  , 1)
   CALL ncio_write_vector (landname, 'CV_WALL'   , 'ulev', ulev, 'urban', landurban, cv_wl  , 1)
   CALL ncio_write_vector (landname, 'TK_ROOF'   , 'ulev', ulev, 'urban', landurban, tk_rf  , 1)
   CALL ncio_write_vector (landname, 'TK_WALL'   , 'ulev', ulev, 'urban', landurban, tk_wl  , 1)
   CALL ncio_write_vector (landname, 'TK_IMPROAD', 'ulev', ulev, 'urban', landurban, tk_imrd, 1)
   CALL ncio_write_vector (landname, 'CV_IMPROAD', 'ulev', ulev, 'urban', landurban, cv_imrd, 1)

   CALL ncio_write_vector (landname, 'ALB_ROOF'   , 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_rf  , 1)
   CALL ncio_write_vector (landname, 'ALB_WALL'   , 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_wl  , 1)
   CALL ncio_write_vector (landname, 'ALB_IMPROAD', 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_imrd, 1)
   CALL ncio_write_vector (landname, 'ALB_PERROAD', 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_perd, 1)

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#endif

   IF (p_is_worker) THEN
      IF (allocated(LUCY_coun)) deallocate (LUCY_coun)
   ENDIF

END SUBROUTINE makeurbandata

#endif
