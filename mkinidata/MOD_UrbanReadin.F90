#include <define.h>

#ifdef URBAN_MODEL

MODULE MOD_UrbanReadin

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Urban_readin


!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE Urban_readin (dir_landdata, lc_year)!(dir_srfdata,dir_atmdata,nam_urbdata,nam_atmdata,lc_year)

! ===========================================================
! Read in the Urban dataset
! ===========================================================

      USE MOD_Precision
      USE MOD_SPMD_Task
      USE MOD_Vars_Global
      USE MOD_Namelist
      USE MOD_Const_LC
      USE MOD_Vars_TimeVariables
      USE MOD_Vars_TimeInvariants
      USE MOD_Urban_Vars_TimeInvariants
      USE MOD_NetCDFVector
      USE MOD_NetCDFSerial
      USE MOD_LandPatch
      USE MOD_LandUrban
      USE MOD_Urban_Const_LCZ

      IMPLICIT NONE

      INTEGER, intent(in) :: lc_year    ! which year of land cover data used
      CHARACTER(LEN=256), intent(in) :: dir_landdata

      CHARACTER(LEN=256) :: dir_rawdata
      CHARACTER(LEN=256) :: lndname
      CHARACTER(len=256) :: cyear

      INTEGER :: i, u, m, l, lucy_id, ns, nr, ulev

      REAL(r8) :: thick_roof, thick_wall

#ifdef USE_POINT_DATA
#ifdef USE_OBS_PARA
      REAL(r8) :: rfwt, rfht, tpct, wpct, hw_point, htop_point, prwt
#endif
#endif
      ! parameters for LUCY
      INTEGER , allocatable :: lucyid(:)          ! LUCY region id
      REAL(r8), allocatable :: popden(:)          ! population density [person/km2]

      INTEGER , allocatable :: lweek_holiday(:,:) ! week holidays
      REAL(r8), allocatable :: lwdh_prof    (:,:) ! Diurnal traffic flow profile [-]
      REAL(r8), allocatable :: lweh_prof    (:,:) ! Diurnal traffic flow profile [-]
      REAL(r8), allocatable :: lhum_prof    (:,:) ! Diurnal metabolic heat profile profile [W/person]
      REAL(r8), allocatable :: lfix_holiday (:,:) ! Fixed public holidays, holiday(0) or workday(1)
      REAL(r8), allocatable :: lvehicle     (:,:) ! vehicle numbers per thousand people

      ! thickness of roof and wall
      REAL(r8), allocatable :: thickroof     (:)  ! thickness of roof [m]
      REAL(r8), allocatable :: thickwall     (:)  ! thickness of wall [m]

      write(cyear,'(i4.4)') lc_year

      allocate (lucyid    (numurban))

IF (DEF_URBAN_type_scheme == 1) THEN
      allocate (thickroof (numurban))
      allocate (thickwall (numurban))

      ! READ in urban data
      lndname = trim(dir_landdata)//'/urban/'//trim(cyear)//'/urban.nc'
      print*,trim(lndname)
      CALL ncio_read_vector (lndname, 'CANYON_HWR  '  , landurban, hwr    ) ! average building height to their distance
      CALL ncio_read_vector (lndname, 'WTROAD_PERV'   , landurban, fgper  ) ! pervious fraction to ground area
      CALL ncio_read_vector (lndname, 'EM_ROOF'       , landurban, em_roof) ! emissivity of roof
      CALL ncio_read_vector (lndname, 'EM_WALL'       , landurban, em_wall) ! emissivity of wall
      CALL ncio_read_vector (lndname, 'EM_IMPROAD'    , landurban, em_gimp) ! emissivity of impervious road
      CALL ncio_read_vector (lndname, 'EM_PERROAD'    , landurban, em_gper) ! emissivity of pervious road

      CALL ncio_read_vector (lndname, 'T_BUILDING_MAX', landurban, t_roommax) ! maximum temperature of inner room [K]
      CALL ncio_read_vector (lndname, 'T_BUILDING_MIN', landurban, t_roommin) ! minimum temperature of inner room [K]
      CALL ncio_read_vector (lndname, 'THICK_ROOF'    , landurban, thickroof) ! thickness of roof [m]
      CALL ncio_read_vector (lndname, 'THICK_WALL'    , landurban, thickwall) ! thickness of wall [m]

      CALL ncio_read_vector (lndname, 'ALB_ROOF'      , 2, 2, landurban, alb_roof) ! albedo of roof
      CALL ncio_read_vector (lndname, 'ALB_WALL'      , 2, 2, landurban, alb_wall) ! albedo of wall
      CALL ncio_read_vector (lndname, 'ALB_IMPROAD'   , 2, 2, landurban, alb_gimp) ! albedo of impervious
      CALL ncio_read_vector (lndname, 'ALB_PERROAD'   , 2, 2, landurban, alb_gper) ! albedo of pervious road

      CALL ncio_read_vector (lndname, 'CV_ROOF'       , nl_roof, landurban, cv_roof) ! heat capacity of roof [J/(m2 K)]
      CALL ncio_read_vector (lndname, 'CV_WALL'       , nl_wall, landurban, cv_wall) ! heat capacity of wall [J/(m2 K)]
      CALL ncio_read_vector (lndname, 'CV_IMPROAD'    , nl_soil, landurban, cv_gimp) ! heat capacity of impervious road [J/(m2 K)]
      CALL ncio_read_vector (lndname, 'TK_ROOF'       , nl_roof, landurban, tk_roof) ! thermal conductivity of roof [W/m-K]
      CALL ncio_read_vector (lndname, 'TK_WALL'       , nl_wall, landurban, tk_wall) ! thermal conductivity of wall [W/m-K]
      CALL ncio_read_vector (lndname, 'TK_IMPROAD'    , nl_soil, landurban, tk_gimp) ! thermal conductivity of impervious road [W/m-K]
ENDIF

!TODO: add point case
#ifdef SinglePoint
      lndname = trim(dir_atmdata)//'/'//trim(nam_atmdata)
      print*, lndname
      CALL ncio_read_bcast_serial (landname, "impervious_area_fraction" , prwt    ) ! imperivous area fraciton to total surface
      CALL ncio_read_bcast_serial (landname, "tree_area_fraction"       , fveg_urb) ! urban tree percentage
      CALL ncio_read_bcast_serial (landname, "water_area_fraction"      , flake   ) ! urban lake precentage
      CALL ncio_read_bcast_serial (landname, "roof_area_fraction"       , froof   ) ! roof fractional cover
      CALL ncio_read_bcast_serial (landname, "building_mean_height"     , hroof   ) ! average building height
      CALL ncio_read_bcast_serial (landname, "tree_mean_height"         , htop_urb) ! urban tree crown top
      CALL ncio_read_bcast_serial (landname, "canyon_height_width_ratio", hwr     ) ! average building height to their distance

      wtperroad    (1,1,:) = 1 - (prwt-rfwt)/(1-rfwt-wpct) !1. - prwt
#endif

      !TODO: Variables distinguish between time-varying and time-invariant variables
      ! write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/urban/'//trim(cyear)//'/POP.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'POP_DEN'       , landurban, pop_den )
      ! write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/urban/'//trim(cyear)//'/LUCY_country_id.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'LUCY_id'       , landurban, lucyid  )
      ! write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/urban/'//trim(cyear)//'/WT_ROOF.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'WT_ROOF'       , landurban, froof   )
      ! write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/urban/'//trim(cyear)//'/HT_ROOF.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'HT_ROOF'       , landurban, hroof   )

      lndname = trim(dir_landdata)//'/urban/'//trim(cyear)//'/PCT_Water.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'PCT_Water'     , landurban, flake   )

      lndname = trim(dir_landdata)//'/urban/'//trim(cyear)//'/PCT_Tree.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'PCT_Tree'      , landurban, fveg_urb)

      lndname = trim(dir_landdata)//'/urban/'//trim(cyear)//'/htop_urb.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'URBAN_TREE_TOP', landurban, htop_urb)

      dir_rawdata = DEF_dir_rawdata
      lndname = trim(dir_rawdata)//'/urban/'//'/LUCY_rawdata.nc'
      print*, lndname
      CALL ncio_read_bcast_serial (lndname,  "NUMS_VEHC"             , lvehicle     )
      CALL ncio_read_bcast_serial (lndname,  "WEEKEND_DAY"           , lweek_holiday)
      CALL ncio_read_bcast_serial (lndname,  "TraffProf_24hr_holiday", lweh_prof    )
      CALL ncio_read_bcast_serial (lndname,  "TraffProf_24hr_work"   , lwdh_prof    )
      CALL ncio_read_bcast_serial (lndname,  "HumMetabolic_24hr"     , lhum_prof    )
      CALL ncio_read_bcast_serial (lndname,  "FIXED_HOLIDAY"         , lfix_holiday )

      IF (p_is_worker) THEN

         DO u = 1, numurban

            i       = urban2patch (u)
            lucy_id = lucyid      (u)

            IF (DEF_URBAN_LUCY) THEN
               IF (lucy_id > 0) THEN
                  vehicle      (:,u) = lvehicle      (lucy_id,:)
                  week_holiday (:,u) = lweek_holiday (lucy_id,:)
                  weh_prof     (:,u) = lweh_prof     (lucy_id,:)
                  wdh_prof     (:,u) = lwdh_prof     (lucy_id,:)
                  hum_prof     (:,u) = lhum_prof     (lucy_id,:)
                  fix_holiday  (:,u) = lfix_holiday  (lucy_id,:)
               ENDIF
            ELSE
               pop_den        (u) = 0.
               vehicle      (:,u) = 0.
               week_holiday (:,u) = 0.
               weh_prof     (:,u) = 0.
               wdh_prof     (:,u) = 0.
               hum_prof     (:,u) = 0.
               fix_holiday  (:,u) = 0.
            ENDIF

IF (DEF_URBAN_type_scheme == 1) THEN
            thick_roof = thickroof (u) !thickness of roof [m]
            thick_wall = thickwall (u) !thickness of wall [m]

            IF (all(cv_gimp(:,u)==0)) THEN
               fgper(u) = 1.
            ENDIF

            IF ( .not. DEF_URBAN_BEM) THEN
               t_roommax(u) = 373.16
               t_roommin(u) = 180.00
            ENDIF
ELSE IF (DEF_URBAN_type_scheme == 2) THEN
            ! read in LCZ constants
            hwr  (u) = canyonhwr_lcz (landurban%settyp(u)) !average building height to their distance
            fgper(u) = wtperroad_lcz (landurban%settyp(u)) !pervious fraction to ground area

            DO ns = 1,2
               DO nr = 1,2
                  alb_roof(ns,nr,u) = albroof_lcz    (landurban%settyp(u)) !albedo of roof
                  alb_wall(ns,nr,u) = albwall_lcz    (landurban%settyp(u)) !albedo of walls
                  alb_gimp(ns,nr,u) = albimproad_lcz (landurban%settyp(u)) !albedo of impervious
                  alb_gper(ns,nr,u) = albperroad_lcz (landurban%settyp(u)) !albedo of pervious road
               ENDDO
            ENDDO

            em_roof(u) = emroof_lcz    (landurban%settyp(u)) !emissivity of roof
            em_wall(u) = emwall_lcz    (landurban%settyp(u)) !emissiviry of wall
            em_gimp(u) = emimproad_lcz (landurban%settyp(u)) !emissivity of impervious
            em_gper(u) = emperroad_lcz (landurban%settyp(u)) !emissivity of pervious

            DO ulev = 1, nl_roof
               cv_roof(:,u) = cvroof_lcz (landurban%settyp(u)) !heat capacity of roof [J/(m2 K)]
               tk_roof(:,u) = tkroof_lcz (landurban%settyp(u)) !thermal conductivity of roof [W/m-K]
            ENDDO

            DO ulev = 1, nl_wall
               cv_wall(:,u) = cvwall_lcz (landurban%settyp(u)) !heat capacity of wall [J/(m2 K)]
               tk_wall(:,u) = tkwall_lcz (landurban%settyp(u)) !thermal conductivity of wall [W/m-K]
            ENDDO

            DO ulev = 1, nl_soil
               cv_gimp(:,u) = cvimproad_lcz (landurban%settyp(u)) !heat capacity of impervious [J/(m2 K)]
               tk_gimp(:,u) = tkimproad_lcz (landurban%settyp(u)) !thermal conductivity of impervious [W/m-K]
            ENDDO

            thick_roof = thickroof_lcz (landurban%settyp(u)) !thickness of roof [m]
            thick_wall = thickwall_lcz (landurban%settyp(u)) !thickness of wall [m]

            IF (DEF_URBAN_BEM) THEN
               t_roommax(u) = 297.65 !tbuildingmax  (landurban%settyp(u)) !maximum temperature of inner room [K]
               t_roommin(u) = 290.65 !tbuildingmin  (landurban%settyp(u)) !minimum temperature of inner room [K]
            ELSE
               t_roommax(u) = 373.16                !maximum temperature of inner room [K]
               t_roommin(u) = 180.00                !minimum temperature of inner room [K]
            ENDIF
ENDIF

            IF (DEF_URBAN_WATER) THEN
               flake(u) = flake(u)/100. !urban water fractional cover
            ELSE
               flake(u) = 0.
            ENDIF

            IF (DEF_URBAN_TREE) THEN
               ! set tree fractional cover (<= 1.-froof)
               fveg_urb(u) = fveg_urb(u)/100. !urban tree percent
               IF (flake(u) < 1.) THEN
                  fveg_urb(u) = fveg_urb(u)/(1.-flake(u))
               ELSE
                  fveg_urb(u) = 0.
               ENDIF

               ! Assuming that the tree coverage is less than or equal
               ! to the ground cover (assume no trees on the roof)
               fveg_urb(u) = min(fveg_urb(u), 1.-froof(u))
               fveg    (i) = fveg_urb(u)
            ELSE
               fveg_urb(u) = 0.
               fveg    (i) = fveg_urb(u)
            ENDIF

            ! set urban tree crown top and bottom [m]
            htop_urb(u) = min(hroof(u), htop_urb(u))
            htop_urb(u) = max(2., htop_urb(u))

            hbot_urb(u) = htop_urb(u)*hbot0(URBAN)/htop0(URBAN)
            hbot_urb(u) = max(1., hbot_urb(u))

            htop    (i) = htop_urb (u)
            hbot    (i) = hbot_urb (u)

            ! roof and wall layer depth
            DO l=1, nl_roof
               z_roof(l,u) = (l-0.5)*(thick_roof/nl_roof)
            ENDDO

            DO l=1, nl_wall
               z_wall(l,u) = (l-0.5)*(thick_wall/nl_wall)
            ENDDO

            dz_roof(1,u) = 0.5*(z_roof(1,u)+z_roof(2,u))
            DO l = 2, nl_roof-1
               dz_roof(l,u) = 0.5*(z_roof(l+1,u)-z_roof(l-1,u))
            ENDDO
            dz_roof(nl_roof,u) = z_roof(nl_roof,u)-z_roof(nl_roof-1,u)

            dz_wall(1,u) = 0.5*(z_wall(1,u)+z_wall(2,u))
            DO l = 2, nl_wall-1
               dz_wall(l,u) = 0.5*(z_wall(l+1,u)-z_wall(l-1,u))
            ENDDO
            dz_wall(nl_wall,u) = z_wall(nl_wall,u)-z_wall(nl_wall-1,u)

            ! lake depth and layer depth
            !NOTE: USE global lake depth right now, the below set to 1m
            !lakedepth(npatch) = 1.
            !dz_lake(:,npatch) = lakedepth(npatch) / nl_lake
         ENDDO
      ENDIF

      IF (p_is_worker) THEN
         IF (allocated(lvehicle     )) deallocate ( lvehicle      )
         IF (allocated(lwdh_prof    )) deallocate ( lwdh_prof     )
         IF (allocated(lweh_prof    )) deallocate ( lweh_prof     )
         IF (allocated(lhum_prof    )) deallocate ( lhum_prof     )
         IF (allocated(lweek_holiday)) deallocate ( lweek_holiday )
         IF (allocated(lfix_holiday )) deallocate ( lfix_holiday  )
         IF (allocated(thickroof    )) deallocate ( thickroof     )
         IF (allocated(thickwall    )) deallocate ( thickwall     )
         IF (allocated(lucyid       )) deallocate ( lucyid        )
      ENDIF

   END SUBROUTINE Urban_readin

END MODULE MOD_UrbanReadin

#endif
