#include <define.h>

SUBROUTINE initialize (casename, dir_landdata, dir_restart, &
      idate, greenwich)

   ! ======================================================================
   ! initialization routine for land surface model.
   !
   ! Created by Yongjiu Dai, 09/15/1999
   ! Revised by Yongjiu Dai, 08/30/2002
   ! Revised by Yongjiu Dai, 03/2014
   !             
   ! ======================================================================
   use precision
   USE GlobalVars
   use mod_namelist
   use spmd_task
   use mod_pixel
   use mod_landpatch
   use PhysicalConstants
   use MOD_TimeInvariants
   use MOD_TimeVariables
#ifdef PFT_CLASSIFICATION
   USE mod_landpft
   USE MOD_PFTimeInvars
   USE MOD_PFTimeVars
#endif
#ifdef PC_CLASSIFICATION
   USE mod_landpc
   USE MOD_PCTimeInvars
   USE MOD_PCTimeVars
#endif
   USE LC_Const
   USE PFT_Const
   use timemanager

   use mod_grid
   use mod_data_type
!   use mod_mapping_grid2pset
   use ncio_serial
   use ncio_block
   use mod_colm_debug
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   USE mod_soil_function
#endif

   IMPLICIT NONE

   ! ----------------------------------------------------------------------
   character(len=*), intent(in) :: casename      ! case name
   character(len=*), intent(in) :: dir_landdata
   character(len=*), intent(in) :: dir_restart
   integer, intent(inout) :: idate(3)   ! year, julian day, seconds of the starting time
   logical, intent(in)    :: greenwich  ! true: greenwich time, false: local time

   ! ------------------------ local variables -----------------------------
   real(r8) :: rlon, rlat
   INTEGER  :: month, mday, m

#if(defined SOILINI)
   character(len=256) :: fsoildat

   type(grid_type) :: gsoil
   type(mapping_grid2pset_type) :: ms2p

   integer :: nl_soil_ini
   
   real(r8), allocatable :: soil_z(:)
   
   type(block_data_real8_2d) :: snow_d_grid
   type(block_data_real8_3d) :: soil_t_grid
   type(block_data_real8_3d) :: soil_w_grid

   real(r8), allocatable :: snow_d(:)
   real(r8), allocatable :: soil_t(:,:)
   real(r8), allocatable :: soil_w(:,:)
#endif

   ! CLM soil layer thickiness and depths
   real(r8), allocatable :: z_soisno (:,:)
   real(r8), allocatable :: dz_soisno(:,:)

   real(r8) :: calday                   ! Julian cal day (1.xx to 365.xx)
   integer  :: year, jday, msec               ! Julian day and seconds
   integer  :: i,ipatch,nsl  ! indices

   integer :: Julian_8day
   integer :: ltyp

   real(r8), external :: orb_coszen     ! cosine of the solar zenith angle


   ! ----------------------------------------------------------------------
   ! [1] READ IN LAND INFORMATION
   ! read time-invariant boundary data on [lon_points] x [lat_points] grid.
   ! ----------------------------------------------------------------------

   ! ----------------------------------------------------------------------
   ! [2] MAPPING and ALLOCATE
   !-----------------------------------------------------------------------

   ! --------------------------------------------------------------------
   ! Allocates memory for CLM 1d [numpatch] variables
   ! --------------------------------------------------------------------

   CALL allocate_TimeInvariants 
   CALL allocate_TimeVariables  

   ! ---------------------------------------------------------------
   ! [3] INITIALIZE TIME INVARIANT VARIABLES
   ! ---------------------------------------------------------------

   if (p_is_worker) then
      
      patchclass = landpatch%ltyp

      DO ipatch = 1, numpatch
         patchtype(ipatch) = patchtypes(patchclass(ipatch))
      ENDDO

      call landpatch%get_lonlat_radian (patchlonr, patchlatr)

#ifdef PFT_CLASSIFICATION
      pftclass = landpft%ltyp
#endif
   
   ENDIF

#ifdef PFT_CLASSIFICATION
   CALL pct_readin (dir_landdata)
#endif
   
#ifdef PC_CLASSIFICATION
   CALL pct_readin (dir_landdata)
#endif

   ! ------------------------------------------
   ! Depth to bedrock
   ! ------------------------------------------
#ifdef USE_DEPTH_TO_BEDROCK
   CALL dbedrock_readin (dir_landdata)
#endif

   ! ------------------------------------------
   ! Lake depth and layers' thickness
   ! ------------------------------------------
   CALL lakedepth_readin (dir_landdata)

   ! ...............................................................
   ! 3.2 Read in the soil parameters of the patches of the gridcells
   ! ...............................................................

   CALL soil_parameters_readin (dir_landdata)

#ifdef vanGenuchten_Mualem_SOIL_MODEL
   IF (p_is_worker) THEN
      IF (numpatch > 0) THEN
         DO ipatch = 1, numpatch
            DO i = 1, nl_soil
               CALL get_derived_parameters_vGM ( &
                  porsl(i,ipatch), alpha_vgm(i,ipatch), n_vgm(i,ipatch), &
                  sc_vgm(i,ipatch), fc_vgm(i,ipatch))
            ENDDO 
         ENDDO 
      ENDIF 
   ENDIF  
#endif

   ! ...............................................................
   ! 3.3 Plant time-invariant variables (based on the look-up tables) 
   ! ...............................................................

   ! read global tree top height from nc file
   CALL HTOP_readin (dir_landdata)

   ! ................................
   ! 3.4 Initialize TUNABLE constants
   ! ................................
   zlnd   = 0.01    !Roughness length for soil [m]
   zsno   = 0.0024  !Roughness length for snow [m]
   csoilc = 0.004   !Drag coefficient for soil under canopy [-]
   dewmx  = 0.1     !maximum dew
   wtfact = 0.38    !Maximum saturated fraction (global mean; see Niu et al., 2005)
   capr   = 0.34    !Tuning factor to turn first layer T into surface T
   cnfac  = 0.5     !Crank Nicholson factor between 0 and 1
   ssi    = 0.033   !Irreducible water saturation of snow
   wimp   = 0.05    !Water impremeable if porosity less than wimp
   pondmx = 10.0    !Ponding depth (mm)
   smpmax = -1.5e5  !Wilting point potential in mm
   smpmin = -1.e8   !Restriction for min of soil poten. (mm)
   trsmx0 = 2.e-4   !Max transpiration for moist soil+100% veg. [mm/s]
   tcrit  = 2.5     !critical temp. to determine rain or snow

   ! ...............................................
   ! 3.5 Write out as a restart file [histTimeConst]
   ! ...............................................

#ifdef CLMDEBUG 
   call check_TimeInvariants ()
#endif

   CALL WRITE_TimeInvariants (casename, dir_restart)

#ifdef USEMPI
   call mpi_barrier (p_comm_glb, p_err)
#endif

   if (p_is_master) write (6,*) ('Successfully Initialize the Land Time-Invariants')

   ! ----------------------------------------------------------------------
   ! [4] INITIALIZE TIME-VARYING VARIABLES 
   ! as subgrid vectors of length [numpatch]
   ! initial run: create the time-varying variables based on :
   !              i) observation (NOT CODING CURRENTLY), or
   !             ii) some already-known information (NO CODING CURRENTLY), or
   !            iii) arbitrarily 
   ! continuation run: time-varying data read in from restart file 
   ! ----------------------------------------------------------------------

   ! 4.1 current time of model run
   ! ............................

   call initimetype(greenwich)
      
   IF (p_is_master) THEN
      IF(.not. greenwich)THEN
         print *, ".........greenwich false"
      ENDIF
   ENDIF

   year = idate(1)
   jday = idate(2)
   msec = idate(3)

   ! ................................
   ! 4.2 cosine of solar zenith angle 
   ! ................................
   calday = calendarday(idate)
   if (p_is_worker) then
      do i = 1, numpatch
         coszen(i) = orb_coszen(calday, patchlonr(i), patchlatr(i))
      enddo
   end if

#if(defined SOILINI)
      ! ...........................................
      !4.3 READ in or GUSSES land state information
      ! ...........................................
      !!! PLEASE CHANGE
      !!! PLEASE CHANGE
      !!! PLEASE CHANGE the root of directory when the soil T and W are ready !!!

      fsoildat = 'change/directory/and/file/path.nc'

      call gsoil%define_from_file (fsoildat)
      call ms2p%build (gsoil, ilev_patch)

      call ncio_read_bcast_serial (fsoildat, 'soil_z', soil_z)
      nl_soil_ini = size(soil_z)

      if (p_is_io) then

         call allocate_block_data (gsoil, snow_d_grid)
         call allocate_block_data (gsoil, soil_t_grid, nl_soil_ini)
         call allocate_block_data (gsoil, soil_w_grid, nl_soil_ini)

         call ncio_read_block (fsoildat, 'soil_t', gsoil, nl_soil_ini, soil_t_grid)  ! soil layer temperature (K)
         call ncio_read_block (fsoildat, 'soil_w', gsoil, nl_soil_ini, soil_w_grid)  ! soil layer wetness (-)
         call ncio_read_block (fsoildat, 'snow_d', gsoil, snow_d_grid)  ! snow depth (m)              

      end if

      if (p_is_worker) then

         allocate (snow_d(numpatch))
         allocate (soil_t(nl_soil_ini,numpatch))
         allocate (soil_w(nl_soil_ini,numpatch))

      end if

      call ms2p%map_aweighted (soil_t_grid, nl_soil_ini, soil_t)
      call ms2p%map_aweighted (soil_w_grid, nl_soil_ini, soil_w)
      call ms2p%map_aweighted (snow_d_grid, snow_d)

   end if
#endif

   ! ...................
   ! 4.4 LEAF area index
   ! ...................
#if(defined DYN_PHENOLOGY)
      ! CREAT fraction of vegetation cover, greenness, leaf area index, stem index
      if (p_is_worker) then

         tlai(:)=0.0; tsai(:)=0.0; green(:)=0.0; fveg(:)=0.0
         do i = 1, numpatch
#if(defined SOILINI)
               do nsl = 1, nl_soil
                  t_soisno(nsl,i) = soil_t(nl_soil_ini,i)
               enddo
#else
               t_soisno(1:,i) = 283.
#endif
            ! Call Ecological Model() 
            ltyp = patchtype(i)
            if(ltyp > 0) then
               call lai_empirical(ltyp, nl_soil,rootfr(1:,i), t_soisno(1:,i),tlai(i),tsai(i),fveg(i),green(i))
            endif
         enddo

      end if
#else

#ifdef USGS_CLASSIFICATION
      ! READ in Leaf area index and stem area index
      Julian_8day = int(calendarday(idate)-1)/8*8 + 1
      CALL LAI_readin (Julian_8day, dir_landdata)

      CALL check_vector_data ('LAI ', tlai)
      CALL check_vector_data ('SAI ', tsai)
#else
! yuan, 08/03/2019: read global LAI/SAI data
      CALL julian2monthday (year, jday, month, mday)
      CALL LAI_readin (month, dir_landdata)
#endif

#endif

CALL check_vector_data ('porsl', porsl)
   ! ..............................................................................
   ! 4.5 initialize time-varying variables, as subgrid vectors of length [numpatch]
   ! ..............................................................................
   if (p_is_worker) then
      
      allocate ( z_soisno (maxsnl+1:nl_soil,numpatch) )
      allocate ( dz_soisno(maxsnl+1:nl_soil,numpatch) )

      do i = 1, numpatch
         z_soisno (1:nl_soil ,i) = z_soi (1:nl_soil)
         dz_soisno(1:nl_soil ,i) = dz_soi(1:nl_soil)
      enddo

      do i = 1, numpatch
         m = patchclass(i)
         CALL iniTimeVar(i, patchtype(i)&
            ,porsl(1:,i)&
            ,soil_s_v_alb(i),soil_d_v_alb(i),soil_s_n_alb(i),soil_d_n_alb(i)&
            ,z0m(i),zlnd,chil(m),rho(1:,1:,m),tau(1:,1:,m)&
            ,z_soisno(maxsnl+1:,i),dz_soisno(maxsnl+1:,i)&
            ,t_soisno(maxsnl+1:,i),wliq_soisno(maxsnl+1:,i),wice_soisno(maxsnl+1:,i)&
            ,zwt(i),wa(i)&
            ,t_grnd(i),tleaf(i),ldew(i),sag(i),scv(i)&
            ,snowdp(i),fveg(i),fsno(i),sigf(i),green(i),lai(i),sai(i),coszen(i)&
            ,alb(1:,1:,i),ssun(1:,1:,i),ssha(1:,1:,i)&
            ,thermk(i),extkb(i),extkd(i)&
            ,trad(i),tref(i),qref(i),rst(i),emis(i),zol(i),rib(i)&
            ,ustar(i),qstar(i),tstar(i),fm(i),fh(i),fq(i)&
#if(defined SOILINI)
            ,nl_soil_ini,soil_z,soil_t(1:,i),soil_w(1:,i),snow_d(i))
#else
          )
#endif
      enddo

      do i = 1, numpatch
         z_sno (maxsnl+1:0,i) = z_soisno (maxsnl+1:0,i)
         dz_sno(maxsnl+1:0,i) = dz_soisno(maxsnl+1:0,i)
      end do

   end if

   ! ------------------------------------------
   ! PLEASE  
   ! PLEASE UPDATE
   ! PLEASE UPDATE when have the observed lake status
   if (p_is_worker) then

      t_lake      (:,:) = 285.
      lake_icefrac(:,:) = 0.

   end if
   ! ------------------------------------------

   ! ...............................................................
   ! 4.6 Write out the model variables for restart run [histTimeVar]
   ! ...............................................................
   
#ifdef CLMDEBUG 
   call check_TimeVariables ()
#endif

   CALL WRITE_TimeVariables (idate, casename, dir_restart)

#ifdef USEMPI
   call mpi_barrier (p_comm_glb, p_err)
#endif

   if (p_is_master) write (6,*) ('Successfully Initialize the Land Time-Vraying Variables')


   ! --------------------------------------------------
   ! Deallocates memory for CLM 1d [numpatch] variables
   ! --------------------------------------------------

   CALL deallocate_TimeInvariants
   CALL deallocate_TimeVariables 

   if (p_is_worker) then

      deallocate (z_soisno )
      deallocate (dz_soisno)

   end if

#if(defined SOILINI)
      deallocate (soil_z)
      deallocate (snow_d)
      deallocate (soil_t)
      deallocate (soil_w)
#endif

END SUBROUTINE initialize
! --------------------------------------------------
! EOP
