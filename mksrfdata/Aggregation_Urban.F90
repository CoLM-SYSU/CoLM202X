#include <define.h>

!-----------------------------------------------------------------------
!
! !DESCRIPTION:
!  Aggregate/screen high-resolution urban dataset to a lower
!  resolution/subset data, suitable for running regional or point
!  cases.
!
!  Original authors: Hua Yuan and Wenzong Dong, 2021, OpenMP version.
!
! !REVISIONS:
!  05/2023, Wenzong Dong, Hua Yuan, Shupeng Zhang: porting codes to MPI
!           parallel version.
!
!-----------------------------------------------------------------------

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
   USE MOD_LandElm
   USE MOD_Mesh
   USE MOD_Vars_Global, only: N_URB
   USE MOD_Urban_Const_LCZ
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif

   IMPLICIT NONE

   character(len=256), intent(in) :: dir_rawdata
   character(len=256), intent(in) :: dir_srfdata

   integer , intent(in) :: lc_year

   type(grid_type), intent(in) :: grid_urban_5km
  !type(grid_type), intent(in) :: grid_urban_100m
   type(grid_type), intent(in) :: grid_urban_500m

   ! dimensions
   integer, parameter :: rid  = 33
   integer, parameter :: ns   = 2
   integer, parameter :: nr   = 2
   integer, parameter :: ulev = 10

   ! input variables
   type(block_data_int32_2d) :: LUCY_reg
   type(block_data_real8_2d) :: pop
   type(block_data_real8_2d) :: fvegu
   type(block_data_real8_2d) :: htopu
   type(block_data_real8_2d) :: flakeu
   type(block_data_real8_2d) :: wtrf
   type(block_data_real8_2d) :: htrf
   type(block_data_real8_2d) :: ulai
   type(block_data_real8_2d) :: usai
   type(block_data_int32_2d) :: reg_typid

   ! output variables
   integer , ALLOCATABLE, dimension(:) :: LUCY_rid
   real(r8), ALLOCATABLE, dimension(:) :: pop_den
   real(r8), ALLOCATABLE, dimension(:) :: pct_tree
   real(r8), ALLOCATABLE, dimension(:) :: htop_urb
   real(r8), ALLOCATABLE, dimension(:) :: pct_urbwt
   real(r8), ALLOCATABLE, dimension(:) :: wt_roof
   real(r8), ALLOCATABLE, dimension(:) :: ht_roof
   real(r8), ALLOCATABLE, dimension(:) :: lai_urb
   real(r8), ALLOCATABLE, dimension(:) :: sai_urb

   ! delete variables not used
   integer , allocatable, dimension(:) :: reg_typid_one
   integer , allocatable, dimension(:) :: LUCY_reg_one
   real(r8), allocatable, dimension(:) :: area_one
   real(r8), allocatable, dimension(:) :: pop_one
   real(r8), allocatable, dimension(:) :: fvegu_one
   real(r8), allocatable, dimension(:) :: htopu_one
   real(r8), allocatable, dimension(:) :: flakeu_one
   real(r8), allocatable, dimension(:) :: wt_roof_one
   real(r8), allocatable, dimension(:) :: ht_roof_one
   real(r8), allocatable, dimension(:) :: ulai_one
   real(r8), allocatable, dimension(:) :: slai_one

   ! urban morphological and thermal paras of NCAR data
   ! input variables, look-up-table data
   real(r8), allocatable, dimension(:,:)     :: hlrbld , wtrd   , ncar_ht, ncar_wt
   real(r8), allocatable, dimension(:,:)     :: emroof , emwall , emimrd , emperd
   real(r8), allocatable, dimension(:,:)     :: throof , thwall , tbmin  , tbmax
   real(r8), allocatable, dimension(:,:,:)   :: cvroof , cvwall , cvimrd , &
                                                tkroof , tkwall , tkimrd
   real(r8), allocatable, dimension(:,:,:,:) :: albroof, albwall, albimrd, albperd

   ! output variables, vector data
   real(r8), ALLOCATABLE, dimension(:)     :: area_urb
   real(r8), ALLOCATABLE, dimension(:)     :: sarea_urb
   real(r8), ALLOCATABLE, dimension(:)     :: urb_frc
   real(r8), ALLOCATABLE, dimension(:)     :: urb_pct

   real(r8), ALLOCATABLE, dimension(:)     :: hlr_bld
   real(r8), ALLOCATABLE, dimension(:)     :: wt_rd
   real(r8), ALLOCATABLE, dimension(:)     :: em_roof
   real(r8), ALLOCATABLE, dimension(:)     :: em_wall
   real(r8), ALLOCATABLE, dimension(:)     :: em_imrd
   real(r8), ALLOCATABLE, dimension(:)     :: em_perd
   real(r8), ALLOCATABLE, dimension(:)     :: th_roof
   real(r8), ALLOCATABLE, dimension(:)     :: th_wall
   real(r8), ALLOCATABLE, dimension(:)     :: tb_min
   real(r8), ALLOCATABLE, dimension(:)     :: tb_max

   real(r8), ALLOCATABLE, dimension(:,:)   :: cv_wgt
   real(r8), ALLOCATABLE, dimension(:,:)   :: tk_wgt
   real(r8), ALLOCATABLE, dimension(:,:)   :: cv_roof
   real(r8), ALLOCATABLE, dimension(:,:)   :: cv_wall
   real(r8), ALLOCATABLE, dimension(:,:)   :: cv_imrd
   real(r8), ALLOCATABLE, dimension(:,:)   :: tk_roof
   real(r8), ALLOCATABLE, dimension(:,:)   :: tk_wall
   real(r8), ALLOCATABLE, dimension(:,:)   :: tk_imrd

   real(r8), ALLOCATABLE, dimension(:,:,:) :: alb_roof
   real(r8), ALLOCATABLE, dimension(:,:,:) :: alb_wall
   real(r8), ALLOCATABLE, dimension(:,:,:) :: alb_imrd
   real(r8), ALLOCATABLE, dimension(:,:,:) :: alb_perd

   integer , allocatable, dimension(:)     :: locpxl

   ! landfile variables
   character(len=256) landsrfdir, landdir, landname, suffix
   character(len=4)   cyear, c5year, cmonth, clay, c1, iyear

   ! local vars
   real(r8) :: sumarea

   ! index
   integer :: iurban, urb_typidx, urb_regidx
   integer :: pop_i, imonth, start_year, end_year
   integer :: ipxstt, ipxend, ipxl, il, iy, i, numpxl, urb_s, urb_e, urb2p

   ! for surface data diag
#ifdef SrfdataDiag
   integer  :: ityp
   integer, allocatable, dimension(:) :: typindex
#endif

#ifdef SrfdataDiag
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

         landname = TRIM(dir_rawdata)//'urban/LUCY_regionid.nc'
         CALL allocate_block_data (grid_urban_5km, LUCY_reg)
         CALL ncio_read_block (landname, 'LUCY_REGION_ID', grid_urban_5km, LUCY_reg)

#ifdef USEMPI
         CALL aggregation_data_daemon (grid_urban_5km, data_i4_2d_in1 = LUCY_reg)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate ( LUCY_rid (numurban))

         LUCY_rid (:) = 0

         ! loop for each urban patch to get the LUCY id of all fine grid
         ! of iurban patch, then assign the most frequency id to this urban patch
         DO iurban = 1, numurban

            CALL aggregation_request_data (landurban, iurban, grid_urban_5km, zip = USE_zip_for_aggregation, &
               data_i4_2d_in1 = LUCY_reg, data_i4_2d_out1 = LUCY_reg_one)
            ! the most frequency id to this urban patch
            LUCY_rid(iurban) = num_max_frequency (LUCY_reg_one)
         ENDDO
#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

      ! output
      landname = trim(landsrfdir)//'/LUCY_region_id.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'LUCY_id', 'urban', landurban, LUCY_rid, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typindex = (/(ityp, ityp = 1, N_URB)/)
      landname  = trim(dir_srfdata) // '/diag/LUCY_region_id.nc'
      ! CALL srfdata_map_and_write (LUCY_rid*1.0, landurban%settyp, typindex, m_urb2diag, &
      !    -1.0e36_r8, landname, 'LUCY_id_'//trim(cyear), compress = 0, write_mode = 'one')
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('LUCY_ID ', LUCY_rid)
#endif

      ! ******* POP_DEN *******
      ! allocate and read the grided population raw data(500m)
      ! NOTE, the population is year-by-year
      IF (p_is_io) THEN

         CALL allocate_block_data (grid_urban_500m, pop)

         landdir = TRIM(dir_rawdata)//'/urban/'
         suffix  = 'URBSRF'//trim(c5year)

         ! population data is year by year,
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
            CALL aggregation_request_data (landurban, iurban, grid_urban_500m, zip = USE_zip_for_aggregation, &
               area = area_one, data_r8_2d_in1 = pop, data_r8_2d_out1 = pop_one)

            WHERE (pop_one < 0)
               area_one = 0
            END WHERE
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
      CALL ncio_write_vector (landname, 'POP_DEN', 'urban', landurban, pop_den, DEF_Srfdata_CompressLevel)

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

         CALL allocate_block_data (grid_urban_500m, fvegu)
         CALL read_5x5_data (landdir, suffix, grid_urban_500m, "PCT_Tree", fvegu)

         CALL allocate_block_data (grid_urban_500m, htopu)
         CALL read_5x5_data (landdir, suffix, grid_urban_500m, "HTOP", htopu)

#ifdef USEMPI
         CALL aggregation_data_daemon (grid_urban_500m, &
            data_r8_2d_in1 = fvegu, data_r8_2d_in2 = htopu)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (pct_tree(numurban))
         allocate (htop_urb(numurban))

         pct_tree(:) = 0.
         htop_urb(:) = 0.

         ! loop for urban patch to aggregate tree cover and height data with area-weighted average
         DO iurban = 1, numurban
            CALL aggregation_request_data (landurban, iurban, grid_urban_500m, zip = USE_zip_for_aggregation, area = area_one, &
               data_r8_2d_in1 = fvegu, data_r8_2d_out1 = fvegu_one, &
               data_r8_2d_in2 = htopu, data_r8_2d_out2 = htopu_one)

            ! missing tree cover and tree height data (-999) were filtered
            WHERE (fvegu_one < 0)
               area_one = 0
            END WHERE

            WHERE (htopu_one < 0)
               area_one = 0
            END WHERE

            ! area-weighted average
            IF (sum(area_one) > 0._r8) THEN
               ! print*, sum(area_one)
               pct_tree(iurban) = sum(fvegu_one * area_one) / sum(area_one)
               htop_urb(iurban) = sum(htopu_one * area_one) / sum(area_one)
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
      CALL ncio_write_vector (landname, 'PCT_Tree', 'urban', landurban, pct_tree, DEF_Srfdata_CompressLevel)

      landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/htop_urb.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'URBAN_TREE_TOP', 'urban', landurban, htop_urb, DEF_Srfdata_CompressLevel)

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

         CALL allocate_block_data (grid_urban_500m, flakeu)

         landdir = TRIM(dir_rawdata)//'/urban/'
         suffix  = 'URBSRF'//trim(c5year)
         CALL read_5x5_data (landdir, suffix, grid_urban_500m, "PCT_Water", flakeu)

#ifdef USEMPI
         CALL aggregation_data_daemon (grid_urban_500m, flakeu)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (pct_urbwt (numurban))

         pct_urbwt (:) = 0.
         ! loop for urban patch to aggregate water cover data with area-weighted average
         DO iurban = 1, numurban
            CALL aggregation_request_data (landurban, iurban, grid_urban_500m, zip = USE_zip_for_aggregation, area = area_one, &
               data_r8_2d_in1 = flakeu, data_r8_2d_out1 = flakeu_one)

            WHERE (flakeu_one < 0)
               area_one = 0
            END WHERE
            ! only calculate when urban patch have water cover
            IF (sum(area_one) > 0) THEN
               pct_urbwt(iurban) = sum(flakeu_one * area_one) / sum(area_one)
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
      CALL ncio_write_vector (landname, 'PCT_Water', 'urban', landurban, pct_urbwt, DEF_Srfdata_CompressLevel)

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
      ! if building data is missing, how to use look-up-table?
      ! a new array with region id was used for look-up-table (urban_reg)
      IF (DEF_URBAN_type_scheme == 1) THEN
         ! only used when urban patch have nan data of building height and fraction
         landname = TRIM(dir_rawdata)//'urban/NCAR_urban_properties.nc'

         CALL ncio_read_bcast_serial (landname,  "WTLUNIT_ROOF"       , ncar_wt )
         CALL ncio_read_bcast_serial (landname,  "HT_ROOF"            , ncar_ht )
      ENDIF

      ! allocate and read grided building height and cover raw data
      IF (p_is_io) THEN
         CALL allocate_block_data (grid_urban_500m, reg_typid)
         CALL allocate_block_data (grid_urban_500m, wtrf)
         CALL allocate_block_data (grid_urban_500m, htrf)

         landdir = TRIM(dir_rawdata)//'urban_type/'
         suffix  = 'URBTYP'
         CALL read_5x5_data (landdir, suffix, grid_urban_500m, "REGION_ID", reg_typid)

         landdir = TRIM(dir_rawdata)//'/urban/'
         suffix  = 'URBSRF'//trim(c5year)
IF (DEF_Urban_geom_data == 1) THEN
         CALL read_5x5_data (landdir, suffix, grid_urban_500m, "PCT_ROOF_GHSL", wtrf)
         CALL read_5x5_data (landdir, suffix, grid_urban_500m, "HT_ROOF_GHSL" , htrf)
ELSE
         CALL read_5x5_data (landdir, suffix, grid_urban_500m, "PCT_ROOF_Li", wtrf)
         CALL read_5x5_data (landdir, suffix, grid_urban_500m, "HT_ROOF_Li" , htrf)
ENDIF

#ifdef USEMPI
         CALL aggregation_data_daemon (grid_urban_500m, data_i4_2d_in1 = reg_typid, &
            data_r8_2d_in1 = wtrf, data_r8_2d_in2 = htrf)
#endif
      ENDIF

      IF (p_is_worker) THEN
         allocate (wt_roof  (numurban))
         allocate (ht_roof  (numurban))
         allocate (urb_pct  (numurban))
         allocate (urb_frc  (numurban))
         allocate (sarea_urb(numurban))
         allocate (area_urb (numurban))

         sarea_urb(:)     = 0.
         area_urb (:)     = 0.
         ! loop for urban patch to aggregate building height and fraction data with area-weighted average
         DO iurban = 1, numurban
            CALL aggregation_request_data (landurban, iurban, grid_urban_500m, zip = USE_zip_for_aggregation, area = area_one, &
               data_i4_2d_in1 = reg_typid, data_i4_2d_out1 = reg_typid_one, &
               data_r8_2d_in1 = wtrf, data_r8_2d_out1 = wt_roof_one, &
               data_r8_2d_in2 = htrf, data_r8_2d_out2 = ht_roof_one)

            area_urb(iurban) = sum(area_one)
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
                  WHERE(reg_typid_one==0) reg_typid_one = num_max_frequency(reg_typid_one)
               ENDIF

               WHERE (wt_roof_one <= 0)
                  wt_roof_one = ncar_wt(urb_typidx,reg_typid_one)
               END WHERE

               WHERE (ht_roof_one <= 0)
                  ht_roof_one = ncar_ht(urb_typidx,reg_typid_one)
               END WHERE
            ELSE IF (DEF_URBAN_type_scheme == 2) THEN
               ! same for above, but for LCZ case
               ! LCZ type for look-up-table
               urb_typidx     = landurban%settyp(iurban)

               WHERE (wt_roof_one <= 0)
                  wt_roof_one = wtroof_lcz(urb_typidx)
               END WHERE

               WHERE (ht_roof_one <= 0)
                  ht_roof_one = htroof_lcz(urb_typidx)
               END WHERE
            ENDIF

            ! area-weight average
            wt_roof(iurban) = sum(wt_roof_one * area_one) / sum(area_one)
            ht_roof(iurban) = sum(ht_roof_one * area_one) / sum(area_one)

         ENDDO

         DO i = 1, numelm
            numpxl = count(landurban%eindex==landelm%eindex(i))

            IF (allocated(locpxl)) deallocate(locpxl)
            allocate(locpxl(numpxl))

            locpxl = pack([(ipxl, ipxl=1, numurban)], &
                     landurban%eindex==landelm%eindex(i))

            urb_s = minval(locpxl)
            urb_e = maxval(locpxl)

            DO il = urb_s, urb_e
               sarea_urb(urb_s:urb_e) = sarea_urb(urb_s:urb_e) + area_urb(il)
            ENDDO
         ENDDO

         DO i = 1, numurban
            urb2p       = urban2patch(i)
            urb_frc (i) = elm_patch%subfrc(urb2p)
            urb_pct (i) = area_urb(i)/sarea_urb(i)
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

      ! output
      landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/WT_ROOF.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'WT_ROOF', 'urban', landurban, wt_roof, DEF_Srfdata_CompressLevel)

      landname = trim(dir_srfdata) // '/urban/'//trim(cyear)//'/HT_ROOF.nc'
      CALL ncio_create_file_vector (landname, landurban)
      CALL ncio_define_dimension_vector (landname, landurban, 'urban')
      CALL ncio_write_vector (landname, 'HT_ROOF', 'urban', landurban, ht_roof, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typindex = (/(ityp, ityp = 1, N_URB)/)
      landname  = trim(dir_srfdata) // '/diag/ht_roof.nc'
      CALL srfdata_map_and_write (ht_roof, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'HT_ROOF_'//trim(cyear), compress = 0, write_mode = 'one')

      typindex = (/(ityp, ityp = 1, N_URB)/)
      landname  = trim(dir_srfdata) // '/diag/wt_roof.nc'
      CALL srfdata_map_and_write (wt_roof, landurban%settyp, typindex, m_urb2diag, &
         -1.0e36_r8, landname, 'WT_ROOF_'//trim(cyear), compress = 0, write_mode = 'one')

      typindex = (/(ityp, ityp = 1, N_URB)/)
      landname  = trim(dir_srfdata) // '/diag/pct_urban' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (urb_pct(:), landurban%settyp, typindex, m_urb2diag, &
      -1.0e36_r8, landname, 'URBAN_PCT', compress = 0, write_mode = 'one')

      CALL srfdata_map_and_write (urb_frc(:), landurban%settyp, typindex, m_urb2diag, &
      -1.0e36_r8, landname, 'URBAN_PATCH_FRAC', compress = 0, write_mode = 'one')
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
         start_year = max(DEF_LAI_START_YEAR, DEF_simulation_time%start_year)
         end_year   = min(DEF_LAI_END_YEAR,   DEF_simulation_time%end_year  )
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

         IF (iy < 2000) THEN
            write(iyear,'(i4.4)') 2000
         ELSE
            write(iyear,'(i4.4)') iy
         ENDIF

         landsrfdir = trim(dir_srfdata) // '/urban/' // trim(iyear) // '/LAI'
         CALL system('mkdir -p ' // trim(adjustl(landsrfdir)))

         ! allocate and read grided LSAI raw data
         landdir = trim(dir_rawdata)//'/urban_lai_500m/'
         suffix  = 'URBLAI_'//trim(iyear)

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
                  data_r8_2d_in1 = fvegu, data_r8_2d_in2 = ulai, data_r8_2d_in3 = usai)
#endif
            ENDIF

            IF (p_is_worker) THEN

               ! loop for urban patch to aggregate LSAI data
               DO iurban = 1, numurban
                  CALL aggregation_request_data (landurban, iurban, grid_urban_500m, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_2d_in1 = fvegu, data_r8_2d_out1 = fvegu_one, &
                     data_r8_2d_in2 = ulai   , data_r8_2d_out2 = ulai_one   , &
                     data_r8_2d_in3 = usai   , data_r8_2d_out3 = slai_one   )

                  WHERE (fvegu_one < 0)
                     area_one = 0
                  END WHERE

                  ! area-weight average
                  IF (sum(fvegu_one * area_one) > 0) THEN
                     lai_urb(iurban) = sum(ulai_one * fvegu_one * area_one) / &
                                       sum(fvegu_one * area_one)
                     sai_urb(iurban) = sum(slai_one * fvegu_one * area_one) / &
                                       sum(fvegu_one * area_one)
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
            CALL ncio_write_vector (landname, 'TREE_LAI', 'urban', landurban, lai_urb, DEF_Srfdata_CompressLevel)

            landname = trim(dir_srfdata) // '/urban/'//trim(iyear)//'/LAI/urban_SAI_'//trim(cmonth)//'.nc'
            CALL ncio_create_file_vector (landname, landurban)
            CALL ncio_define_dimension_vector (landname, landurban, 'urban')
            CALL ncio_write_vector (landname, 'TREE_SAI', 'urban', landurban, sai_urb, DEF_Srfdata_CompressLevel)

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

         CALL ncio_read_bcast_serial (landname,  "CANYON_HWR"    , hlrbld   )
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

            allocate (hlr_bld          (numurban))
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
            hlr_bld  (:)     = 0.
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
               CALL aggregation_request_data (landurban, iurban, grid_urban_500m, zip = USE_zip_for_aggregation, area = area_one, &
                  data_i4_2d_in2 = reg_typid, data_i4_2d_out2 = reg_typid_one)

               ! urban region and type id for look-up-table
               urb_typidx = landurban%settyp(iurban)
               !urb_regidx = urban_reg(iurban)

               ! ipxstt = landurban%ipxstt(iurban)
               ! ipxend = landurban%ipxend(iurban)

               sumarea          = sum(area_one)
               area_urb(iurban) = sumarea

               ! same for above, assign reg id for RG_-45_65_-50_70
               IF (all(reg_typid_one==0)) THEN
                  reg_typid_one(:) = 30
               ENDIF

               IF (any(reg_typid_one==0)) THEN
                  WHERE(reg_typid_one==0) reg_typid_one =  num_max_frequency(reg_typid_one)
                  ENDIF

                  ! loop for each finer grid to aggregate data
                  DO ipxl = 1, size(area_one) ! ipxstt, ipxend

                     urb_regidx = reg_typid_one(ipxl)
                     hlr_bld (iurban) = hlr_bld (iurban) + hlrbld (urb_typidx,urb_regidx) * area_one(ipxl)
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

                  hlr_bld  (iurban) = hlr_bld   (iurban) / sumarea
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

         CALL ncio_write_vector (landname, 'URBAN_PCT'     , 'urban', landurban, urb_pct, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'URBAN_FRAC'    , 'urban', landurban, urb_frc, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'BUILDING_HLR'  , 'urban', landurban, hlr_bld, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'WTROAD_PERV'   , 'urban', landurban, wt_rd  , DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'EM_ROOF'       , 'urban', landurban, em_roof, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'EM_WALL'       , 'urban', landurban, em_wall, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'EM_IMPROAD'    , 'urban', landurban, em_imrd, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'EM_PERROAD'    , 'urban', landurban, em_perd, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'THICK_ROOF'    , 'urban', landurban, th_roof, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'THICK_WALL'    , 'urban', landurban, th_wall, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'T_BUILDING_MIN', 'urban', landurban, tb_min , DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'T_BUILDING_MAX', 'urban', landurban, tb_max , DEF_Srfdata_CompressLevel)

         CALL ncio_write_vector (landname, 'CV_ROOF'   , 'ulev', ulev, 'urban', landurban, cv_roof, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'CV_WALL'   , 'ulev', ulev, 'urban', landurban, cv_wall, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'TK_ROOF'   , 'ulev', ulev, 'urban', landurban, tk_roof, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'TK_WALL'   , 'ulev', ulev, 'urban', landurban, tk_wall, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'TK_IMPROAD', 'ulev', ulev, 'urban', landurban, tk_imrd, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'CV_IMPROAD', 'ulev', ulev, 'urban', landurban, cv_imrd, DEF_Srfdata_CompressLevel)

         CALL ncio_write_vector (landname, 'ALB_ROOF'   , 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_roof, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'ALB_WALL'   , 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_wall, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'ALB_IMPROAD', 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_imrd, DEF_Srfdata_CompressLevel)
         CALL ncio_write_vector (landname, 'ALB_PERROAD', 'numsolar', ns, 'numrad', nr, 'urban', landurban, alb_perd, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
         typindex = (/(ityp, ityp = 1, N_URB)/)
         landname  = trim(dir_srfdata) // '/diag/hlr.nc'
         CALL srfdata_map_and_write (hlr_bld, landurban%settyp, typindex, m_urb2diag, &
            -1.0e36_r8, landname, 'BUILDING_HLR_'//trim(cyear), compress = 0, write_mode = 'one')

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
         CALL check_vector_data ('BUILDING_HLR '  , hlr_bld  )
         CALL check_vector_data ('WTROAD_PERV '   , wt_rd    )
         CALL check_vector_data ('EM_ROOF '       , em_roof  )
         CALL check_vector_data ('EM_WALL '       , em_wall  )
         CALL check_vector_data ('EM_IMPROAD '    , em_imrd  )
         CALL check_vector_data ('EM_PERROAD '    , em_perd  )
         CALL check_vector_data ('ALB_ROOF '      , alb_roof )
         CALL check_vector_data ('ALB_WALL '      , alb_wall )
         CALL check_vector_data ('ALB_IMPROAD '   , alb_imrd )
         CALL check_vector_data ('ALB_PERROAD '   , alb_perd )
         CALL check_vector_data ('TK_ROOF '       , tk_roof  )
         CALL check_vector_data ('TK_WALL '       , tk_wall  )
         CALL check_vector_data ('TK_IMPROAD '    , tk_imrd  )
         CALL check_vector_data ('CV_ROOF '       , cv_roof  )
         CALL check_vector_data ('CV_WALL '       , cv_wall  )
         CALL check_vector_data ('CV_IMPROAD '    , cv_imrd  )
         CALL check_vector_data ('THICK_ROOF '    , th_roof  )
         CALL check_vector_data ('THICK_WALL '    , th_wall  )
         CALL check_vector_data ('T_BUILDING_MIN ', tb_min   )
         CALL check_vector_data ('T_BUILDING_MAX ', tb_max   )
#endif

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ENDIF

      IF (p_is_worker) THEN

         IF ( allocated (LUCY_rid ) ) deallocate (LUCY_rid  )
         IF ( allocated (pop_den  ) ) deallocate (pop_den   )
         IF ( allocated (pct_tree ) ) deallocate (pct_tree  )
         IF ( allocated (htop_urb ) ) deallocate (htop_urb  )
         IF ( allocated (pct_urbwt) ) deallocate (pct_urbwt )
         IF ( allocated (wt_roof  ) ) deallocate (wt_roof   )
         IF ( allocated (ht_roof  ) ) deallocate (ht_roof   )
         IF ( allocated (lai_urb  ) ) deallocate (lai_urb   )
         IF ( allocated (sai_urb  ) ) deallocate (sai_urb   )
         IF ( allocated (area_urb ) ) deallocate (area_urb  )
         IF ( allocated (sarea_urb) ) deallocate (sarea_urb )
         IF ( allocated (urb_frc  ) ) deallocate (urb_frc   )
         IF ( allocated (urb_pct  ) ) deallocate (urb_pct   )

         IF (DEF_URBAN_type_scheme == 1) THEN

            IF ( allocated (ncar_ht  ) ) deallocate (ncar_ht   )
            IF ( allocated (ncar_wt  ) ) deallocate (ncar_wt   )
            IF ( allocated (hlr_bld  ) ) deallocate (hlr_bld   )
            IF ( allocated (wt_rd    ) ) deallocate (wt_rd     )
            IF ( allocated (em_roof  ) ) deallocate (em_roof   )
            IF ( allocated (em_wall  ) ) deallocate (em_wall   )
            IF ( allocated (em_imrd  ) ) deallocate (em_imrd   )
            IF ( allocated (em_perd  ) ) deallocate (em_perd   )
            IF ( allocated (th_roof  ) ) deallocate (th_roof   )
            IF ( allocated (th_wall  ) ) deallocate (th_wall   )
            IF ( allocated (tb_min   ) ) deallocate (tb_min    )
            IF ( allocated (tb_max   ) ) deallocate (tb_max    )
            IF ( allocated (tk_wgt   ) ) deallocate (tk_wgt    )
            IF ( allocated (cv_wgt   ) ) deallocate (cv_wgt    )
            IF ( allocated (cv_roof  ) ) deallocate (cv_roof   )
            IF ( allocated (cv_wall  ) ) deallocate (cv_wall   )
            IF ( allocated (cv_imrd  ) ) deallocate (cv_imrd   )
            IF ( allocated (tk_roof  ) ) deallocate (tk_roof   )
            IF ( allocated (tk_wall  ) ) deallocate (tk_wall   )
            IF ( allocated (tk_imrd  ) ) deallocate (tk_imrd   )
            IF ( allocated (alb_roof ) ) deallocate (alb_roof  )
            IF ( allocated (alb_wall ) ) deallocate (alb_wall  )
            IF ( allocated (alb_imrd ) ) deallocate (alb_imrd  )
            IF ( allocated (alb_perd ) ) deallocate (alb_perd  )

         ENDIF

         IF ( allocated (area_one     ) ) deallocate (area_one     )
         IF ( allocated (LUCY_reg_one ) ) deallocate (LUCY_reg_one )
         IF ( allocated (pop_one      ) ) deallocate (pop_one      )
         IF ( allocated (fvegu_one    ) ) deallocate (fvegu_one    )
         IF ( allocated (htopu_one    ) ) deallocate (htopu_one    )
         IF ( allocated (flakeu_one   ) ) deallocate (flakeu_one   )
         IF ( allocated (wt_roof_one  ) ) deallocate (wt_roof_one  )
         IF ( allocated (ht_roof_one  ) ) deallocate (ht_roof_one  )
         IF ( allocated (ulai_one     ) ) deallocate (ulai_one     )
         IF ( allocated (slai_one     ) ) deallocate (slai_one     )

      ENDIF

END SUBROUTINE Aggregation_Urban
#endif
