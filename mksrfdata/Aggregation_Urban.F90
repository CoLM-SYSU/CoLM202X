#include <define.h>

! ======================================================
! Aggreate/screen high-resolution urban dataset
! to a lower resolutioin/subset data, suitable for running
! regional or point cases.
! ======================================================

#ifdef URBAN_MODEL
SUBROUTINE Aggregation_Urban (dir_rawdata, dir_srfdata, lc_year, &
                              grid_urban_5km, grid_urban_500m)

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_NetCDFSerial
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
   USE MOD_AggregationRequestData
   USE MOD_5x5DataReadin
   USE MOD_DataType
   USE MOD_Namelist
   USE MOD_Utils, only: num_max_frequency
   USE MOD_LandUrban
   USE MOD_Vars_Global, only: N_URB
   USE MOD_Urban_Const_LCZ, only: wtroof_lcz, htroof_lcz
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif

   IMPLICIT NONE

   CHARACTER(len=256), intent(in) :: dir_rawdata
   CHARACTER(len=256), intent(in) :: dir_srfdata

   INTEGER , intent(in) :: lc_year

   TYPE(grid_type), intent(in) :: grid_urban_5km
   ! TYPE(grid_type), intent(in) :: grid_urban_100m
   TYPE(grid_type), intent(in) :: grid_urban_500m

   ! dimensions
   INTEGER, parameter :: rid  = 33
   INTEGER, parameter :: ns   = 2
   INTEGER, parameter :: nr   = 2
   INTEGER, parameter :: ulev = 10

   ! input variables
   TYPE(block_data_int32_2d) :: LUCY_reg
   TYPE(block_data_real8_2d) :: pop
   TYPE(block_data_real8_2d) :: gfcc_tc
   TYPE(block_data_real8_2d) :: gedi_th
   TYPE(block_data_real8_2d) :: gl30_wt
   TYPE(block_data_real8_2d) :: wtrf
   TYPE(block_data_real8_2d) :: htrf
   TYPE(block_data_real8_2d) :: ulai
   TYPE(block_data_real8_2d) :: usai
   TYPE(block_data_int32_2d) :: reg_typid

   ! output variables
   INTEGER , ALLOCATABLE, DIMENSION(:) :: LUCY_coun
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: pop_den
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: pct_tree
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: htop_urb
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: pct_urbwt
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: wt_roof
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: ht_roof
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: lai_urb
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: sai_urb

   ! delete variables not used
   INTEGER , allocatable, dimension(:) :: reg_typid_one
   INTEGER , allocatable, dimension(:) :: LUCY_reg_one
   REAL(r8), allocatable, dimension(:) :: area_one
   REAL(r8), allocatable, dimension(:) :: pop_one
   REAL(r8), allocatable, dimension(:) :: gfcc_tc_one
   REAL(r8), allocatable, dimension(:) :: gedi_th_one
   REAL(r8), allocatable, dimension(:) :: gl30_wt_one
   REAL(r8), allocatable, dimension(:) :: wt_roof_one
   REAL(r8), allocatable, dimension(:) :: ht_roof_one
   REAL(r8), allocatable, dimension(:) :: ulai_one
   REAL(r8), allocatable, dimension(:) :: slai_one

   ! urban morphological and thermal paras of NCAR data
   ! input variables, look-up-table data
   REAL(r8), allocatable, DIMENSION(:,:)  :: hwrcan, wtrd, emroof, emwall, ncar_wt
   REAL(r8), allocatable, DIMENSION(:,:)  :: emimrd, emperd, ncar_ht
   REAL(r8), allocatable, DIMENSION(:,:)  :: throof, thwall, tbmin, tbmax

   REAL(r8), allocatable, DIMENSION(:,:,:)  :: cvroof, cvwall, cvimrd, &
                                               tkroof, tkwall, tkimrd
   REAL(r8), allocatable, DIMENSION(:,:,:,:):: albroof, albwall, albimrd, albperd

   ! output variables, vector data
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: area_urb

   REAL(r8), ALLOCATABLE, DIMENSION(:) :: area_tb
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: area_hd
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: area_md
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: hwr_can
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: wt_rd
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: em_roof
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: em_wall
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: em_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: em_perd
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: th_roof
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: th_wall
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: tb_min
   REAL(r8), ALLOCATABLE, DIMENSION(:) :: tb_max

   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: cv_wgt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: tk_wgt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: cv_roof
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: cv_wall
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: cv_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: tk_roof
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: tk_wall
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: tk_imrd

   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: alb_roof
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: alb_wall
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: alb_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: alb_perd

   ! landfile variables
   CHARACTER(len=256) landsrfdir, landdir, landname, suffix
   CHARACTER(len=4)   cyear, c5year, cmonth, clay, c1, iyear

   ! local vars
   REAL(r8) :: sumarea

   ! index
   INTEGER :: iurban, urb_typidx, urb_regidx
   INTEGER :: pop_i, imonth, start_year, end_year
   INTEGER :: ipxstt, ipxend, ipxl, il, iy

   ! for surface data diag
#ifdef SrfdataDiag
   INTEGER  :: ityp
   INTEGER, allocatable, dimension(:) :: typindex

   allocate( typindex(N_URB) )
#endif

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

   write(c5year, '(i4.4)') int(lc_year/5)*5

   ! ******* LUCY_id *******
   ! allocate and read the LUCY id
   IF (p_is_io) THEN

      landname = TRIM(dir_rawdata)//'urban/LUCY_countryid.nc'
      CALL allocate_block_data (grid_urban_5km, LUCY_reg)
      CALL ncio_read_block (landname, 'LUCY_COUNTRY_ID', grid_urban_5km, LUCY_reg)

#ifdef USEMPI
      CALL aggregation_data_daemon (grid_urban_5km, data_i4_2d_in1 = LUCY_reg)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate ( LUCY_coun (numurban))

      LUCY_coun (:) = 0

      ! loop for each urban patch to get the LUCY id of all fine grid
      ! of iurban patch, then assign the most frequence id to this urban patch
      DO iurban = 1, numurban
         CALL aggregation_request_data (landurban, iurban, grid_urban_5km, &
            data_i4_2d_in1 = LUCY_reg, data_i4_2d_out1 = LUCY_reg_one)
         ! the most frequence id to this urban patch
         LUCY_coun(iurban) = num_max_frequency (LUCY_reg_one)
      ENDDO
#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

   ! output
   landname = trim(landsrfdir)//'/LUCY_country_id.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'LUCY_id', 'urban', landurban, LUCY_coun, 1)

#ifdef SrfdataDiag
   typindex = (/(ityp, ityp = 1, N_URB)/)
   landname  = trim(dir_srfdata) // '/diag/LUCY_country_id.nc'
   ! CALL srfdata_map_and_write (LUCY_coun*1.0, landurban%settyp, typindex, m_urb2diag, &
   !    -1.0e36_r8, landname, 'LUCY_id_'//trim(cyear), compress = 0, write_mode = 'one')
#endif

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
   CALL check_vector_data ('LUCY_ID ', LUCY_coun)
#endif

   ! ******* POP_DEN *******
   ! allocate and read the grided population raw data(500m)
   ! NOTE, the population is year-by-year
   IF (p_is_io) THEN

      CALL allocate_block_data (grid_urban_500m, pop)

      landdir = TRIM(dir_rawdata)//'/urban/'
      suffix  = 'URBSRF'//trim(c5year)

      ! populaiton data is year by year,
      ! so pop_i is calculated to determine the dimension of POP data reads
      IF (mod(lc_year,5) == 0) THEN
         pop_i = 1
      ELSE
         pop_i = 5 - (ceiling(lc_year*1./5.)*5 - lc_year) + 1
      ENDIF

      ! read the population data of total 5x5 region
      CALL read_5x5_data_time (landdir, suffix, grid_urban_500m, "POP_DEN", pop_i, pop)

#ifdef USEMPI
      CALL aggregation_data_daemon (grid_urban_500m, data_r8_2d_in1 = pop)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate (pop_den (numurban))

      pop_den (:) = 0.

      ! loop for urban patch to aggregate population data with area-weighted average
      DO iurban = 1, numurban
         ! request all fine grid data and area of the iurban urban patch
         ! a one dimension vector will be returned
         CALL aggregation_request_data (landurban, iurban, grid_urban_500m, area = area_one, &
            data_r8_2d_in1 = pop, data_r8_2d_out1 = pop_one)

         where (pop_one < 0)
            area_one = 0
         END where
         ! area-weighted average
         IF (sum(area_one) > 0._r8) THEN
            pop_den(iurban) = sum(pop_one * area_one) / sum(area_one)
         ENDIF
      ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

   ! output
   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/POP.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'POP_DEN', 'urban', landurban, pop_den, 1)

#ifdef SrfdataDiag
   typindex = (/(ityp, ityp = 1, N_URB)/)
   landname  = trim(dir_srfdata) // '/diag/POP.nc'
   CALL srfdata_map_and_write (pop_den, landurban%settyp, typindex, m_urb2diag, &
      -1.0e36_r8, landname, 'POP_DEN_'//trim(cyear), compress = 0, write_mode = 'one')
#endif

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
   CALL check_vector_data ('POP_DEN ', pop_den)
#endif

   ! ******* Tree : PCT_Tree, HTOP *******
   ! allocate and read the grided tree cover and tree height raw data(500m)
   ! NOTE, tree cover raw data is available every five years,
   ! tree height raw data is same from year to year
   IF (p_is_io) THEN

      landdir = TRIM(dir_rawdata)//'/urban/'
      suffix  = 'URBSRF'//trim(c5year)

      CALL allocate_block_data (grid_urban_500m, gfcc_tc)
      CALL read_5x5_data (landdir, suffix, grid_urban_500m, "PCT_Tree", gfcc_tc)

      CALL allocate_block_data (grid_urban_500m, gedi_th)
      CALL read_5x5_data (landdir, suffix, grid_urban_500m, "HTOP", gedi_th)

#ifdef USEMPI
      CALL aggregation_data_daemon (grid_urban_500m, &
         data_r8_2d_in1 = gfcc_tc, data_r8_2d_in2 = gedi_th)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate (pct_tree(numurban))
      allocate (htop_urb(numurban))

      pct_tree(:) = 0.
      htop_urb(:) = 0.

      ! loop for urban patch to aggregate tree cover and height data with area-weighted average
      DO iurban = 1, numurban
         CALL aggregation_request_data (landurban, iurban, grid_urban_500m, area = area_one, &
            data_r8_2d_in1 = gfcc_tc, data_r8_2d_out1 = gfcc_tc_one, &
            data_r8_2d_in2 = gedi_th, data_r8_2d_out2 = gedi_th_one)

         ! missing tree cover and tree height data (-999) were filtered
         where (gfcc_tc_one < 0)
            area_one = 0
         END where

         where (gedi_th_one < 0)
            area_one = 0
         END where

         ! area-weighted average
         IF (sum(area_one) > 0._r8) THEN
            ! print*, sum(area_one)
            pct_tree(iurban) = sum(gfcc_tc_one * area_one) / sum(area_one)
            htop_urb(iurban) = sum(gedi_th_one * area_one) / sum(area_one)
         ENDIF
      ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

   ! output
   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/PCT_Tree.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'PCT_Tree', 'urban', landurban, pct_tree, 1)

   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/htop_urb.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'URBAN_TREE_TOP', 'urban', landurban, htop_urb, 1)

#ifdef SrfdataDiag
   typindex = (/(ityp, ityp = 1, N_URB)/)
   landname  = trim(dir_srfdata) // '/diag/PCT_Tree.nc'
   CALL srfdata_map_and_write (pct_tree, landurban%settyp, typindex, m_urb2diag, &
      -1.0e36_r8, landname, 'pct_tree_'//trim(cyear), compress = 0, write_mode = 'one')

   typindex = (/(ityp, ityp = 1, N_URB)/)
   landname  = trim(dir_srfdata) // '/diag/htop_urb.nc'
   CALL srfdata_map_and_write (htop_urb, landurban%settyp, typindex, m_urb2diag, &
      -1.0e36_r8, landname, 'URBAN_TREE_TOP_'//trim(cyear), compress = 0, write_mode = 'one')
#endif

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
   CALL check_vector_data ('Urban Tree Cover ', pct_tree)
   CALL check_vector_data ('Urban Tree Top '  , htop_urb)
#endif

   ! ******* PCT_Water *******
   ! allocate and read grided water cover raw data
   IF (p_is_io) THEN

      CALL allocate_block_data (grid_urban_500m, gl30_wt)

      landdir = TRIM(dir_rawdata)//'/urban/'
      suffix  = 'URBSRF'//trim(c5year)
      CALL read_5x5_data (landdir, suffix, grid_urban_500m, "PCT_Water", gl30_wt)

#ifdef USEMPI
      CALL aggregation_data_daemon (grid_urban_500m, gl30_wt)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate (pct_urbwt (numurban))

      pct_urbwt (:) = 0.
      ! loop for urban patch to aggregate water cover data with area-weighted average
      DO iurban = 1, numurban
         CALL aggregation_request_data (landurban, iurban, grid_urban_500m, area = area_one, &
            data_r8_2d_in1 = gl30_wt, data_r8_2d_out1 = gl30_wt_one)

         where (gl30_wt_one < 0)
            area_one = 0
         END where
         ! only caculate when urban patch have water cover
         IF (sum(area_one) > 0) THEN
            pct_urbwt(iurban) = sum(gl30_wt_one * area_one) / sum(area_one)
         ENDIF
      ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

   ! output
   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/PCT_Water.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'PCT_Water', 'urban', landurban, pct_urbwt, 1)

#ifdef SrfdataDiag
   typindex = (/(ityp, ityp = 1, N_URB)/)
   landname  = trim(dir_srfdata) // '/diag/PCT_Water.nc'
   CALL srfdata_map_and_write (pct_urbwt, landurban%settyp, typindex, m_urb2diag, &
      -1.0e36_r8, landname, 'PCT_Water_'//trim(cyear), compress = 0, write_mode = 'one')
#endif

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
   CALL check_vector_data ('Urban Water Cover ', pct_urbwt)
#endif

   ! ******* Building : Weight, HTOP_Roof *******
   ! if building data is missing, how to look-up-table?
   ! a new arry with region id was used for look-up-table (urban_reg)
IF (DEF_URBAN_type_scheme == 1) THEN
   ! only used when urban patch have nan data of building height and fraction
   landname = TRIM(dir_rawdata)//'urban/NCAR_urban_properties.nc'

   CALL ncio_read_bcast_serial (landname,  "WTLUNIT_ROOF"       , ncar_wt )
   CALL ncio_read_bcast_serial (landname,  "HT_ROOF"            , ncar_ht )
ENDIF

   ! allocate and read grided building hegight and cover raw data
   IF (p_is_io) THEN
      CALL allocate_block_data (grid_urban_500m, reg_typid)
      CALL allocate_block_data (grid_urban_500m, wtrf)
      CALL allocate_block_data (grid_urban_500m, htrf)

      landdir = TRIM(dir_rawdata)//'urban_type/'
      suffix  = 'URBTYP'
      CALL read_5x5_data (landdir, suffix, grid_urban_500m, "REGION_ID", reg_typid)

      landdir = TRIM(dir_rawdata)//'/urban/'
      suffix  = 'URBSRF'//trim(c5year)
      CALL read_5x5_data (landdir, suffix, grid_urban_500m, "PCT_ROOF", wtrf)

      landdir = TRIM(dir_rawdata)//'/urban/'
      suffix  = 'URBSRF'//trim(c5year)
      CALL read_5x5_data (landdir, suffix, grid_urban_500m, "HT_ROOF", htrf)

#ifdef USEMPI
      CALL aggregation_data_daemon (grid_urban_500m, data_i4_2d_in1 = reg_typid, &
         data_r8_2d_in1 = wtrf, data_r8_2d_in2 = htrf)
#endif
   ENDIF

   IF (p_is_worker) THEN
      allocate (wt_roof (numurban))
      allocate (ht_roof (numurban))

      ! loop for urban patch to aggregate building height and fraction data with area-weighted average
      DO iurban = 1, numurban
         CALL aggregation_request_data (landurban, iurban, grid_urban_500m, area = area_one, &
            data_i4_2d_in1 = reg_typid, data_i4_2d_out1 = reg_typid_one, &
            data_r8_2d_in1 = wtrf, data_r8_2d_out1 = wt_roof_one, &
            data_r8_2d_in2 = htrf, data_r8_2d_out2 = ht_roof_one)

IF (DEF_URBAN_type_scheme == 1) THEN
         ! when urban patch has no data, use table data to fill gap
         ! urban type and region id for look-up-table
         urb_typidx = landurban%settyp(iurban)

         ! RG_-45_65_-50_70 of NCAR has no urban data,
         ! all urban patches of this area are assigned to region 30
         IF (all(reg_typid_one==0)) THEN
            reg_typid_one(:) = 30
         ENDIF

         IF (any(reg_typid_one==0)) THEN
            WHERE(reg_typid_one==0) reg_typid_one =  num_max_frequency(reg_typid_one)
         ENDIF

         where (wt_roof_one <= 0)
            wt_roof_one = ncar_wt(urb_typidx,reg_typid_one)
         END where

         where (ht_roof_one <= 0)
            ht_roof_one = ncar_ht(urb_typidx,reg_typid_one)
         END where
ELSE IF (DEF_URBAN_type_scheme == 2) THEN
         ! same for above, but for LCZ case
         ! LCZ type for look-up-table
         urb_typidx     = landurban%settyp(iurban)

         where (wt_roof_one <= 0)
            wt_roof_one = wtroof_lcz(urb_typidx)
         END where

         where (ht_roof_one <= 0)
            ht_roof_one = htroof_lcz(urb_typidx)
         END where
ENDIF

         ! area-weight average
         wt_roof(iurban) = sum(wt_roof_one * area_one) / sum(area_one)
         ht_roof(iurban) = sum(ht_roof_one * area_one) / sum(area_one)

      ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

   ! output
   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/WT_ROOF.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'WT_ROOF', 'urban', landurban, wt_roof, 1)

   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/HT_ROOF.nc'
   CALL ncio_create_file_vector (landname, landurban)
   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_write_vector (landname, 'HT_ROOF', 'urban', landurban, ht_roof, 1)

#ifdef SrfdataDiag
   typindex = (/(ityp, ityp = 1, N_URB)/)
   landname  = trim(dir_srfdata) // '/diag/ht_roof.nc'
   CALL srfdata_map_and_write (ht_roof, landurban%settyp, typindex, m_urb2diag, &
      -1.0e36_r8, landname, 'HT_ROOF_'//trim(cyear), compress = 0, write_mode = 'one')

   typindex = (/(ityp, ityp = 1, N_URB)/)
   landname  = trim(dir_srfdata) // '/diag/wt_roof.nc'
   CALL srfdata_map_and_write (wt_roof, landurban%settyp, typindex, m_urb2diag, &
      -1.0e36_r8, landname, 'WT_ROOF_'//trim(cyear), compress = 0, write_mode = 'one')
#endif

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
   CALL check_vector_data ('Urban Roof Fraction ', wt_roof)
   CALL check_vector_data ('Urban Roof Height '  , ht_roof)
#endif

   ! ******* LAI, SAI *******
#ifndef LULCC
   IF (DEF_LAI_CHANGE_YEARLY) THEN
      start_year = DEF_simulation_time%start_year
      end_year   = DEF_simulation_time%end_year
   ELSE
      start_year = lc_year
      end_year   = lc_year
   ENDIF
#else
   start_year = lc_year
   end_year   = lc_year
#endif



   IF (p_is_io) THEN
      CALL allocate_block_data (grid_urban_500m, ulai)
      CALL allocate_block_data (grid_urban_500m, usai)
   ENDIF

   IF (p_is_worker) THEN
      allocate (lai_urb (numurban))
      allocate (sai_urb (numurban))

      lai_urb(:) = 0.
      sai_urb(:) = 0.
   ENDIF

   DO iy = start_year, end_year
      write(iyear,'(i4.4)') iy
      landsrfdir = trim(dir_srfdata) // '/urban/' // trim(iyear) // '/LAI'
      CALL system('mkdir -p ' // trim(adjustl(landsrfdir)))

      ! allocate and read grided LSAI raw data
      landdir = TRIM(dir_rawdata)//'/urban_lai_5x5/'
      suffix  = 'UrbLAI_'//trim(iyear)
      ! loop for month
      DO imonth = 1, 12

         write(cmonth, '(i2.2)') imonth

         IF (p_is_master) THEN
            write(*,'(A,I4,A1,I3,A1,I3)') 'Aggregate LAI&SAI :', iy, ':', imonth, '/', 12
         ENDIF

         IF (p_is_io) THEN

            CALL read_5x5_data_time (landdir, suffix, grid_urban_500m, "URBAN_TREE_LAI", imonth, ulai)
            CALL read_5x5_data_time (landdir, suffix, grid_urban_500m, "URBAN_TREE_SAI", imonth, usai)

#ifdef USEMPI
            CALL aggregation_data_daemon (grid_urban_500m, &
               data_r8_2d_in1 = gfcc_tc, data_r8_2d_in2 = ulai, data_r8_2d_in3 = usai)
#endif
         ENDIF

         IF (p_is_worker) THEN

            ! loop for urban patch to aggregate LSAI data
            DO iurban = 1, numurban
               CALL aggregation_request_data (landurban, iurban, grid_urban_500m, area = area_one, &
                  data_r8_2d_in1 = gfcc_tc, data_r8_2d_out1 = gfcc_tc_one, &
                  data_r8_2d_in2 = ulai   , data_r8_2d_out2 = ulai_one   , &
                  data_r8_2d_in3 = usai   , data_r8_2d_out3 = slai_one   )

               WHERE (gfcc_tc_one < 0)
                  area_one = 0
               END WHERE

               ! area-weight average
               IF (sum(gfcc_tc_one * area_one) > 0) THEN
                  lai_urb(iurban) = sum(ulai_one * gfcc_tc_one * area_one) / &
                                    sum(gfcc_tc_one * area_one)
                  sai_urb(iurban) = sum(slai_one * gfcc_tc_one * area_one) / &
                                    sum(gfcc_tc_one * area_one)
               ENDIF
            ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
         ENDIF

      ! output
         landname = trim(dir_srfdata) // '/urban/'//trim(iyear)//'/LAI/urban_LAI_'//trim(cmonth)//'.nc'
         CALL ncio_create_file_vector (landname, landurban)
         CALL ncio_define_dimension_vector (landname, landurban, 'urban')
         CALL ncio_write_vector (landname, 'TREE_LAI', 'urban', landurban, lai_urb, 1)

         landname = trim(dir_srfdata) // '/urban/'//trim(iyear)//'/LAI/urban_SAI_'//trim(cmonth)//'.nc'
         CALL ncio_create_file_vector (landname, landurban)
         CALL ncio_define_dimension_vector (landname, landurban, 'urban')
         CALL ncio_write_vector (landname, 'TREE_SAI', 'urban', landurban, sai_urb, 1)

#ifdef SrfdataDiag
         typindex = (/(ityp, ityp = 1, N_URB)/)
         landname  = trim(dir_srfdata) // '/diag/Urban_Tree_LAI_' // trim(iyear) // '.nc'
         CALL srfdata_map_and_write (lai_urb, landurban%settyp, typindex, m_urb2diag, &
            -1.0e36_r8, landname, 'TREE_LAI_'//trim(cmonth), compress = 0, write_mode = 'one')

         typindex = (/(ityp, ityp = 1, N_URB)/)
         landname  = trim(dir_srfdata) // '/diag/Urban_Tree_SAI_' // trim(iyear) // '.nc'
         CALL srfdata_map_and_write (sai_urb, landurban%settyp, typindex, m_urb2diag, &
            -1.0e36_r8, landname, 'TREE_SAI_'//trim(cmonth), compress = 0, write_mode = 'one')
#endif

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

         write(c1,'(i2.2)') imonth

#ifdef RangeCheck
         CALL check_vector_data ('Urban Tree LAI '//trim(c1), lai_urb)
         CALL check_vector_data ('Urban Tree SAI '//trim(c1), sai_urb)
#endif
      ENDDO
   ENDDO

IF (DEF_URBAN_type_scheme == 1) THEN
   ! look up table of NCAR urban properties (using look-up tables)
   landname = TRIM(dir_rawdata)//'urban/NCAR_urban_properties.nc'

   CALL ncio_read_bcast_serial (landname,  "CANYON_HWR"    , hwrcan   )
   CALL ncio_read_bcast_serial (landname,  "WTROAD_PERV"   , wtrd     )
   CALL ncio_read_bcast_serial (landname,  "EM_ROOF"       , emroof   )
   CALL ncio_read_bcast_serial (landname,  "EM_WALL"       , emwall   )
   CALL ncio_read_bcast_serial (landname,  "EM_IMPROAD"    , emimrd   )
   CALL ncio_read_bcast_serial (landname,  "EM_PERROAD"    , emperd   )
   CALL ncio_read_bcast_serial (landname,  "ALB_ROOF"      , albroof  )
   CALL ncio_read_bcast_serial (landname,  "ALB_WALL"      , albwall  )
   CALL ncio_read_bcast_serial (landname,  "ALB_IMPROAD"   , albimrd  )
   CALL ncio_read_bcast_serial (landname,  "ALB_PERROAD"   , albperd  )
   CALL ncio_read_bcast_serial (landname,  "TK_ROOF"       , tkroof   )
   CALL ncio_read_bcast_serial (landname,  "TK_WALL"       , tkwall   )
   CALL ncio_read_bcast_serial (landname,  "TK_IMPROAD"    , tkimrd   )
   CALL ncio_read_bcast_serial (landname,  "CV_ROOF"       , cvroof   )
   CALL ncio_read_bcast_serial (landname,  "CV_WALL"       , cvwall   )
   CALL ncio_read_bcast_serial (landname,  "CV_IMPROAD"    , cvimrd   )
   CALL ncio_read_bcast_serial (landname,  "THICK_ROOF"    , throof   )
   CALL ncio_read_bcast_serial (landname,  "THICK_WALL"    , thwall   )
   CALL ncio_read_bcast_serial (landname,  "T_BUILDING_MIN", tbmin    )
   CALL ncio_read_bcast_serial (landname,  "T_BUILDING_MAX", tbmax    )

   IF (p_is_io) THEN

#ifdef USEMPI
      CALL aggregation_data_daemon (grid_urban_500m, data_i4_2d_in1 = reg_typid)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate (area_urb      (3, numurban))
      allocate (area_tb          (numurban))
      allocate (area_hd          (numurban))
      allocate (area_md          (numurban))
      allocate (hwr_can          (numurban))
      allocate (wt_rd            (numurban))
      allocate (em_roof          (numurban))
      allocate (em_wall          (numurban))
      allocate (em_imrd          (numurban))
      allocate (em_perd          (numurban))
      allocate (th_roof          (numurban))
      allocate (th_wall          (numurban))
      allocate (tb_min           (numurban))
      allocate (tb_max           (numurban))
      allocate (tk_wgt     (ulev, numurban))
      allocate (cv_wgt     (ulev, numurban))
      allocate (cv_roof    (ulev, numurban))
      allocate (cv_wall    (ulev, numurban))
      allocate (cv_imrd    (ulev, numurban))
      allocate (tk_roof    (ulev, numurban))
      allocate (tk_wall    (ulev, numurban))
      allocate (tk_imrd    (ulev, numurban))
      allocate (alb_roof (nr, ns, numurban))
      allocate (alb_wall (nr, ns, numurban))
      allocate (alb_imrd (nr, ns, numurban))
      allocate (alb_perd (nr, ns, numurban))

      ! initialization
      area_urb (:,:)   = 0.
      area_tb  (:)     = 0.
      area_hd  (:)     = 0.
      area_md  (:)     = 0.
      hwr_can  (:)     = 0.
      wt_rd    (:)     = 0.
      em_roof  (:)     = 0.
      em_wall  (:)     = 0.
      em_imrd  (:)     = 0.
      em_perd  (:)     = 0.
      th_roof  (:)     = 0.
      th_wall  (:)     = 0.
      tb_min   (:)     = 0.
      tb_max   (:)     = 0.
      tk_wgt   (:,:)   = 0.
      cv_wgt   (:,:)   = 0.
      cv_roof  (:,:)   = 0.
      cv_wall  (:,:)   = 0.
      cv_imrd  (:,:)   = 0.
      tk_roof  (:,:)   = 0.
      tk_wall  (:,:)   = 0.
      tk_imrd  (:,:)   = 0.
      alb_roof (:,:,:) = 0.
      alb_wall (:,:,:) = 0.
      alb_imrd (:,:,:) = 0.
      alb_perd (:,:,:) = 0.

      ! loop for each urban patch to aggregate NCAR urban morphological and thermal paras with area-weighted average
      DO iurban = 1, numurban
         CALL aggregation_request_data (landurban, iurban, grid_urban_500m, area = area_one, &
                                        data_i4_2d_in2 = reg_typid, data_i4_2d_out2 = reg_typid_one)

         ! urban region and type id for look-up-table
         urb_typidx = landurban%settyp(iurban)
         !urb_regidx = urban_reg(iurban)

         ipxstt = landurban%ipxstt(iurban)
         ipxend = landurban%ipxend(iurban)

         sumarea = sum(area_one)

         ! same for above, assign reg id for RG_-45_65_-50_70
         IF (all(reg_typid_one==0)) THEN
            reg_typid_one(:) = 30
         ENDIF

         IF (any(reg_typid_one==0)) THEN
            WHERE(reg_typid_one==0) reg_typid_one =  num_max_frequency(reg_typid_one)
         ENDIF

         ! loop for each finer grid to aggregate data
         DO ipxl = ipxstt, ipxend

            urb_regidx = reg_typid_one(ipxl)
            area_urb(urb_typidx,iurban) = area_urb(urb_typidx,iurban) + area_one(ipxl)
            hwr_can (iurban) = hwr_can (iurban) + hwrcan (urb_typidx,urb_regidx) * area_one(ipxl)
            wt_rd   (iurban) = wt_rd   (iurban) + wtrd   (urb_typidx,urb_regidx) * area_one(ipxl)
            em_roof (iurban) = em_roof (iurban) + emroof (urb_typidx,urb_regidx) * area_one(ipxl)
            em_wall (iurban) = em_wall (iurban) + emwall (urb_typidx,urb_regidx) * area_one(ipxl)
            em_imrd (iurban) = em_imrd (iurban) + emimrd (urb_typidx,urb_regidx) * area_one(ipxl)
            em_perd (iurban) = em_perd (iurban) + emperd (urb_typidx,urb_regidx) * area_one(ipxl)
            th_roof (iurban) = th_roof (iurban) + throof (urb_typidx,urb_regidx) * area_one(ipxl)
            th_wall (iurban) = th_wall (iurban) + thwall (urb_typidx,urb_regidx) * area_one(ipxl)
            tb_min  (iurban) = tb_min  (iurban) + tbmin  (urb_typidx,urb_regidx) * area_one(ipxl)
            tb_max  (iurban) = tb_max  (iurban) + tbmax  (urb_typidx,urb_regidx) * area_one(ipxl)

            ! tkimrd and cvimrd have nanvalues, and need to be calculated separately
            DO il = 1, 10

               IF (tkimrd(urb_typidx,urb_regidx,il) .ne. -999.) THEN
                  tk_imrd(il,iurban) = tk_imrd(il,iurban) + tkimrd(urb_typidx,urb_regidx,il) * area_one(ipxl)
                  tk_wgt (il,iurban) = tk_wgt (il,iurban) + area_one(ipxl)
               ENDIF

               IF (cvimrd(urb_typidx,urb_regidx,il) .ne. -999.) THEN
                  cv_imrd(il,iurban) = cv_imrd(il,iurban) + cvimrd(urb_typidx,urb_regidx,il) * area_one(ipxl)
                  cv_wgt (il,iurban) = cv_wgt (il,iurban) + area_one(ipxl)
               ENDIF
            ENDDO

            cv_roof (:,iurban) = cv_roof (:,iurban) + cvroof (urb_typidx,urb_regidx,:) * area_one(ipxl)
            cv_wall (:,iurban) = cv_wall (:,iurban) + cvwall (urb_typidx,urb_regidx,:) * area_one(ipxl)
            tk_roof (:,iurban) = tk_roof (:,iurban) + tkroof (urb_typidx,urb_regidx,:) * area_one(ipxl)
            tk_wall (:,iurban) = tk_wall (:,iurban) + tkwall (urb_typidx,urb_regidx,:) * area_one(ipxl)

            alb_roof(:,:,iurban) = alb_roof(:,:,iurban) + albroof(urb_typidx,urb_regidx,:,:) * area_one(ipxl)
            alb_wall(:,:,iurban) = alb_wall(:,:,iurban) + albwall(urb_typidx,urb_regidx,:,:) * area_one(ipxl)
            alb_imrd(:,:,iurban) = alb_imrd(:,:,iurban) + albimrd(urb_typidx,urb_regidx,:,:) * area_one(ipxl)
            alb_perd(:,:,iurban) = alb_perd(:,:,iurban) + albperd(urb_typidx,urb_regidx,:,:) * area_one(ipxl)

         ENDDO

         area_tb  (iurban) = area_urb(1,iurban) / sumarea
         area_hd  (iurban) = area_urb(2,iurban) / sumarea
         area_md  (iurban) = area_urb(3,iurban) / sumarea
         hwr_can  (iurban) = hwr_can   (iurban) / sumarea
         wt_rd    (iurban) = wt_rd     (iurban) / sumarea
         em_roof  (iurban) = em_roof   (iurban) / sumarea
         em_wall  (iurban) = em_wall   (iurban) / sumarea
         em_imrd  (iurban) = em_imrd   (iurban) / sumarea
         em_perd  (iurban) = em_perd   (iurban) / sumarea
         th_roof  (iurban) = th_roof   (iurban) / sumarea
         th_wall  (iurban) = th_wall   (iurban) / sumarea
         tb_min   (iurban) = tb_min    (iurban) / sumarea
         tb_max   (iurban) = tb_max    (iurban) / sumarea

         cv_roof(:,iurban) = cv_roof (:,iurban) / sumarea
         cv_wall(:,iurban) = cv_wall (:,iurban) / sumarea
         tk_roof(:,iurban) = tk_roof (:,iurban) / sumarea
         tk_wall(:,iurban) = tk_wall (:,iurban) / sumarea

         DO il = 1, 10
            IF (tk_wgt(il,iurban) > 0.) THEN
               tk_imrd(il,iurban) = tk_imrd(il,iurban) / tk_wgt(il,iurban)
            ENDIF
            IF (cv_wgt(il,iurban) > 0.) THEN
               cv_imrd(il,iurban) = cv_imrd(il,iurban) / cv_wgt(il,iurban)
            ENDIF
         ENDDO

         alb_roof(:,:,iurban) = alb_roof(:,:,iurban) / sumarea
         alb_wall(:,:,iurban) = alb_wall(:,:,iurban) / sumarea
         alb_imrd(:,:,iurban) = alb_imrd(:,:,iurban) / sumarea
         alb_perd(:,:,iurban) = alb_perd(:,:,iurban) / sumarea

      ENDDO
#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

   !output
   write(cyear,'(i4.4)') lc_year
   landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/urban.nc'
   CALL ncio_create_file_vector (landname, landurban)

   CALL ncio_define_dimension_vector (landname, landurban, 'urban')
   CALL ncio_define_dimension_vector (landname, landurban, 'numsolar', ns)
   CALL ncio_define_dimension_vector (landname, landurban, 'numrad'  , nr)
   CALL ncio_define_dimension_vector (landname, landurban, 'ulev'    , ulev)

   CALL ncio_write_vector (landname, 'PCT_TB'        , 'urban', landurban, area_tb, 1)
   CALL ncio_write_vector (landname, 'PCT_HD'        , 'urban', landurban, area_hd, 1)
   CALL ncio_write_vector (landname, 'PCT_MD'        , 'urban', landurban, area_md, 1)
   CALL ncio_write_vector (landname, 'CANYON_HWR'    , 'urban', landurban, hwr_can, 1)
   CALL ncio_write_vector (landname, 'WTROAD_PERV'   , 'urban', landurban, wt_rd  , 1)
   CALL ncio_write_vector (landname, 'EM_ROOF'       , 'urban', landurban, em_roof, 1)
   CALL ncio_write_vector (landname, 'EM_WALL'       , 'urban', landurban, em_wall, 1)
   CALL ncio_write_vector (landname, 'EM_IMPROAD'    , 'urban', landurban, em_imrd, 1)
   CALL ncio_write_vector (landname, 'EM_PERROAD'    , 'urban', landurban, em_perd, 1)
   CALL ncio_write_vector (landname, 'THICK_ROOF'    , 'urban', landurban, th_roof, 1)
   CALL ncio_write_vector (landname, 'THICK_WALL'    , 'urban', landurban, th_wall, 1)
   CALL ncio_write_vector (landname, 'T_BUILDING_MIN', 'urban', landurban, tb_min , 1)
   CALL ncio_write_vector (landname, 'T_BUILDING_MAX', 'urban', landurban, tb_max , 1)

   CALL ncio_write_vector (landname, 'CV_ROOF'   , 'ulev', ulev, 'urban', landurban, cv_roof, 1)
   CALL ncio_write_vector (landname, 'CV_WALL'   , 'ulev', ulev, 'urban', landurban, cv_wall, 1)
   CALL ncio_write_vector (landname, 'TK_ROOF'   , 'ulev', ulev, 'urban', landurban, tk_roof, 1)
   CALL ncio_write_vector (landname, 'TK_WALL'   , 'ulev', ulev, 'urban', landurban, tk_wall, 1)
   CALL ncio_write_vector (landname, 'TK_IMPROAD', 'ulev', ulev, 'urban', landurban, tk_imrd, 1)
   CALL ncio_write_vector (landname, 'CV_IMPROAD', 'ulev', ulev, 'urban', landurban, cv_imrd, 1)

   CALL ncio_write_vector (landname, 'ALB_ROOF'   , 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_roof, 1)
   CALL ncio_write_vector (landname, 'ALB_WALL'   , 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_wall, 1)
   CALL ncio_write_vector (landname, 'ALB_IMPROAD', 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_imrd, 1)
   CALL ncio_write_vector (landname, 'ALB_PERROAD', 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_perd, 1)

#ifdef SrfdataDiag
   typindex = (/(ityp, ityp = 1, N_URB)/)
   landname  = trim(dir_srfdata) // '/diag/hwr.nc'
   CALL srfdata_map_and_write (hwr_can, landurban%settyp, typindex, m_urb2diag, &
      -1.0e36_r8, landname, 'CANYON_HWR_'//trim(cyear), compress = 0, write_mode = 'one')

   typindex = (/(ityp, ityp = 1, N_URB)/)
   landname  = trim(dir_srfdata) // '/diag/cv_imrd' // trim(cyear) // '.nc'

   DO il = 1, 10
      write(clay, '(i2.2)') il
      CALL srfdata_map_and_write (cv_imrd(il,:), landurban%settyp, typindex, m_urb2diag, &
      -1.0e36_r8, landname, 'CV_IMPROAD', compress = 0, write_mode = 'one', lastdimname = 'layer', lastdimvalue = il)
      ! CALL srfdata_map_and_write (cv_imrd(il,:), landurban%settyp, typindex, m_urb2diag, &
      !    -1.0e36_r8, landname, 'CV_IMPROAD_'//trim(clay), compress = 0, write_mode = 'one')
   ENDDO
   deallocate(typindex)
#endif

#ifdef RangeCheck
   CALL check_vector_data ('CANYON_HWR '    , hwr_can )
   CALL check_vector_data ('WTROAD_PERV '   , wt_rd   )
   CALL check_vector_data ('EM_ROOF '       , em_roof )
   CALL check_vector_data ('EM_WALL '       , em_wall )
   CALL check_vector_data ('EM_IMPROAD '    , em_imrd )
   CALL check_vector_data ('EM_PERROAD '    , em_perd )
   CALL check_vector_data ('ALB_ROOF '      , alb_roof)
   CALL check_vector_data ('ALB_WALL '      , alb_wall)
   CALL check_vector_data ('ALB_IMPROAD '   , alb_imrd)
   CALL check_vector_data ('ALB_PERROAD '   , alb_perd)
   CALL check_vector_data ('TK_ROOF '       , tk_roof )
   CALL check_vector_data ('TK_WALL '       , tk_wall )
   CALL check_vector_data ('TK_IMPROAD '    , tk_imrd )
   CALL check_vector_data ('CV_ROOF '       , cv_roof )
   CALL check_vector_data ('CV_WALL '       , cv_wall )
   CALL check_vector_data ('CV_IMPROAD '    , cv_imrd )
   CALL check_vector_data ('THICK_ROOF '    , th_roof )
   CALL check_vector_data ('THICK_WALL '    , th_wall )
   CALL check_vector_data ('T_BUILDING_MIN ', tb_min  )
   CALL check_vector_data ('T_BUILDING_MAX ', tb_max  )
#endif

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

ENDIF

   IF (p_is_worker) THEN
      IF (allocated(LUCY_coun)) deallocate (LUCY_coun)
      IF (allocated(pop_den  )) deallocate (pop_den  )
      IF (allocated(pct_tree )) deallocate (pct_tree )
      IF (allocated(htop_urb )) deallocate (htop_urb )
      IF (allocated(pct_urbwt)) deallocate (pct_urbwt)
      IF (allocated(wt_roof  )) deallocate (wt_roof  )
      IF (allocated(ht_roof  )) deallocate (ht_roof  )
      IF (allocated(lai_urb  )) deallocate (lai_urb  )
      IF (allocated(sai_urb  )) deallocate (sai_urb  )
IF (DEF_URBAN_type_scheme == 1) THEN
      IF (allocated(ncar_ht  )) deallocate (ncar_ht  )
      IF (allocated(ncar_wt  )) deallocate (ncar_wt  )
      IF (allocated(area_urb )) deallocate (area_urb )
      IF (allocated(hwr_can  )) deallocate (hwr_can  )
      IF (allocated(wt_rd    )) deallocate (wt_rd    )
      IF (allocated(em_roof  )) deallocate (em_roof  )
      IF (allocated(em_wall  )) deallocate (em_wall  )
      IF (allocated(em_imrd  )) deallocate (em_imrd  )
      IF (allocated(em_perd  )) deallocate (em_perd  )
      IF (allocated(th_roof  )) deallocate (th_roof  )
      IF (allocated(th_wall  )) deallocate (th_wall  )
      IF (allocated(tb_min   )) deallocate (tb_min   )
      IF (allocated(tb_max   )) deallocate (tb_max   )
      IF (allocated(tk_wgt   )) deallocate (tk_wgt   )
      IF (allocated(cv_wgt   )) deallocate (cv_wgt   )
      IF (allocated(cv_roof  )) deallocate (cv_roof  )
      IF (allocated(cv_wall  )) deallocate (cv_wall  )
      IF (allocated(cv_imrd  )) deallocate (cv_imrd  )
      IF (allocated(tk_roof  )) deallocate (tk_roof  )
      IF (allocated(tk_wall  )) deallocate (tk_wall  )
      IF (allocated(tk_imrd  )) deallocate (tk_imrd  )
      IF (allocated(alb_roof )) deallocate (alb_roof )
      IF (allocated(alb_wall )) deallocate (alb_wall )
      IF (allocated(alb_imrd )) deallocate (alb_imrd )
      IF (allocated(alb_perd )) deallocate (alb_perd )
ENDIF
      IF (allocated(area_one    )) deallocate(area_one    )
      IF (allocated(LUCY_reg_one)) deallocate(LUCY_reg_one)
      IF (allocated(pop_one     )) deallocate(pop_one     )
      IF (allocated(gfcc_tc_one )) deallocate(gfcc_tc_one )
      IF (allocated(gedi_th_one )) deallocate(gedi_th_one )
      IF (allocated(gl30_wt_one )) deallocate(gl30_wt_one )
      IF (allocated(wt_roof_one )) deallocate(wt_roof_one )
      IF (allocated(ht_roof_one )) deallocate(ht_roof_one )
      IF (allocated(ulai_one    )) deallocate(ulai_one    )
      IF (allocated(slai_one    )) deallocate(slai_one    )
   ENDIF

END SUBROUTINE Aggregation_Urban
#endif
