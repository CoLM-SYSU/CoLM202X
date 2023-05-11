#include <define.h>

SUBROUTINE Urban_readin (dir_landdata)!(dir_srfdata,dir_atmdata,nam_urbdata,nam_atmdata,lc_year)

! ===========================================================
! Read in the Urban dataset
! ===========================================================

      USE precision
      USE spmd_task
      USE GlobalVars
      USE LC_Const
      USE MOD_TimeVariables
      USE MOD_TimeInvariants
      USE MOD_UrbanTimeInvars
      USE ncio_vector
      USE ncio_serial
      USE mod_landpatch
      USE mod_landurban
#ifdef USE_LCZ
      USE UrbanLCZ_Const
#endif

      IMPLICIT NONE

      ! INTEGER, intent(in) :: lc_year    ! which year of land cover data used
      ! CHARACTER(LEN=256), intent(in) :: dir_srfdata
      ! CHARACTER(LEN=256), intent(in) :: nam_urbdata
      ! CHARACTER(LEN=256), intent(in) :: dir_atmdata
      ! CHARACTER(LEN=256), intent(in) :: nam_atmdata
      CHARACTER(LEN=256), intent(in) :: dir_landdata

      CHARACTER(LEN=256) :: lndname
      CHARACTER(len=256) :: cyear

      INTEGER :: ncid
      INTEGER :: wtroof_vid, htroof_vid, canyonhwr_vid, wtperroad_vid, wtimproad_vid
      INTEGER :: urbanwaterpct_vid, urbantreepct_vid, urbantreetop_vid
      INTEGER :: albroof_vid, albwall_vid, albimproad_vid, albperroad_vid
      INTEGER :: emroof_vid, emwall_vid, emimproad_vid, emperroad_vid
      INTEGER :: cvroof_vid, cvwall_vid, cvimproad_vid
      INTEGER :: tkroof_vid, tkwall_vid, tkimproad_vid
      INTEGER :: tbuildingmax_vid, tbuildingmin_vid
      INTEGER :: thickroof_vid, thickwall_vid
      INTEGER :: reg_vid, pop_vid, veh_vid, wed_vid, weh_vid, wdh_vid, met_vid, hol_vid
      INTEGER :: urbanpop_vid, urbanlucy_vid
      INTEGER :: i, j, u, t, l, m, npatch, k, lucy_id, ns, nr, ulev, patch2urb

      REAL(r8) :: thick_roof, thick_wall

      ! define variables for reading in
      ! -------------------------------------------
#ifdef USE_POINT_DATA
#ifdef USE_OBS_PARA
      REAL(r8) :: rfwt, rfht, tpct, wpct, hw_point, htop_point, prwt
#endif
#endif
      ! parameters for LUCY
      INTEGER , allocatable :: urbanlucy(:)       ! LUCY region id

      REAL(r8), allocatable :: popden(:)          ! population density [person/km2]

      INTEGER , allocatable :: lweek_holiday(:,:) ! week holidays
      REAL(r8), allocatable :: lwdh_prof    (:,:) ! Diurnal traffic flow profile [-]
      REAL(r8), allocatable :: lweh_prof    (:,:) ! Diurnal traffic flow profile [-]
      REAL(r8), allocatable :: lhum_prof    (:,:) ! Diurnal metabolic heat profile profile [W/person]
      REAL(r8), allocatable :: lfix_holiday (:,:) ! Fixed public holidays, holiday(0) or workday(1)
      REAL(r8), allocatable :: lvehicle     (:,:) ! vehicle numbers per thousand people

      ! morphological parameter
      REAL(r8), allocatable :: wtroof        (:)  ! roof fraction [-]
      REAL(r8), allocatable :: htroof        (:)  ! height of roof [m]
      REAL(r8), allocatable :: canyonhwr     (:)  ! height/width ratio [-]
      REAL(r8), allocatable :: wtperroad     (:)  ! pervious road fraction [-]
      REAL(r8), allocatable :: urbanwaterpct (:)  ! urban water fraction [-]
      REAL(r8), allocatable :: urbantreepct  (:)  ! urban tree fraction [-]
      REAL(r8), allocatable :: urbantreetop  (:)  ! urban tree height [m]

      ! albedo
      REAL(r8), allocatable :: albroof   (:,:,:)  ! albedo of roof [-]
      REAL(r8), allocatable :: albwall   (:,:,:)  ! albedo of wall[-]
      REAL(r8), allocatable :: albimproad(:,:,:)  ! albedo of impervious road [-]
      REAL(r8), allocatable :: albperroad(:,:,:)  ! albedo of pervious road [-]

      ! emissivity
      REAL(r8), allocatable :: emroof        (:)  ! emissivity of roof [-]
      REAL(r8), allocatable :: emwall        (:)  ! emissivity of wall [-]
      REAL(r8), allocatable :: emimproad     (:)  ! emissivity of impervious road [-]
      REAL(r8), allocatable :: emperroad     (:)  ! emissivity of pervious road [-]

      ! thermal pars of roof, wall, impervious
      REAL(r8), allocatable :: cvroof      (:,:)  ! volumetric heat capacity of roof [J/m3*K]
      REAL(r8), allocatable :: cvwall      (:,:)  ! volumetric heat capacity of wall [J/m3*K]
      REAL(r8), allocatable :: cvimproad   (:,:)  ! volumetric heat capacity of impervious road [J/m3*K]

      ! thermal pars of roof, wall, impervious
      REAL(r8), allocatable :: tkroof      (:,:)  ! thermal conductivity of roof [W/m*K]
      REAL(r8), allocatable :: tkwall      (:,:)  ! thermal conductivity of wall [W/m*K]
      REAL(r8), allocatable :: tkimproad   (:,:)  ! thermal conductivity of impervious road [W/m*K]

      ! room maximum and minimum temperature
      REAL(r8), allocatable :: tbuildingmax  (:)  ! maximum interior building temperature [K]
      REAL(r8), allocatable :: tbuildingmin  (:)  ! minimum interior building temperature [K]

      ! thickness of roof and wall
      REAL(r8), allocatable :: thickroof     (:)  ! thickness of roof [m]
      REAL(r8), allocatable :: thickwall     (:)  ! thickness of wall [m]


#ifdef USE_LCZ

      ! write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/urban/2005/POP.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'POP_DEN'       , landurban, popden        )
      ! write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/urban/2005/LUCY_country_id.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'LUCY_id'       , landurban, urbanlucy     )
      ! write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/urban/2005/WT_ROOF.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'WT_ROOF'       , landurban, wtroof        )
      ! write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/urban/2005/HT_ROOF.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'HT_ROOF'       , landurban, htroof        )

      lndname = trim(dir_landdata)//'/urban/2005/PCT_Water.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'PCT_Water'     , landurban, urbanwaterpct )

      lndname = trim(dir_landdata)//'/urban/2005/PCT_Tree.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'PCT_Tree'      , landurban, urbantreepct  )

      lndname = trim(dir_landdata)//'/urban/2005/htop_urb.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'URBAN_TREE_TOP', landurban, urbantreetop  )

#else
      ! READ in urban data
      ! write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/urban/2005/urban.nc'
      print*,trim(lndname)
      CALL ncio_read_vector (lndname, 'CANYON_HWR  '  , landurban, canyonhwr    )
      CALL ncio_read_vector (lndname, 'WTROAD_PERV'   , landurban, wtperroad    )

      CALL ncio_read_vector (lndname, 'EM_ROOF'       , landurban, emroof       )
      CALL ncio_read_vector (lndname, 'EM_WALL'       , landurban, emwall       )
      CALL ncio_read_vector (lndname, 'EM_IMPROAD'    , landurban, emimproad    )
      CALL ncio_read_vector (lndname, 'EM_PERROAD'    , landurban, emperroad    )

      CALL ncio_read_vector (lndname, 'T_BUILDING_MAX', landurban, tbuildingmax )
      CALL ncio_read_vector (lndname, 'T_BUILDING_MIN', landurban, tbuildingmin )
      CALL ncio_read_vector (lndname, 'THICK_ROOF'    , landurban, thickroof    )
      CALL ncio_read_vector (lndname, 'THICK_WALL'    , landurban, thickwall    )

      CALL ncio_read_vector (lndname, 'ALB_ROOF'      , 2, 2, landurban, albroof   )
      CALL ncio_read_vector (lndname, 'ALB_WALL'      , 2, 2, landurban, albwall   )
      CALL ncio_read_vector (lndname, 'ALB_IMPROAD'   , 2, 2, landurban, albimproad)
      CALL ncio_read_vector (lndname, 'ALB_PERROAD'   , 2, 2, landurban, albperroad)

      CALL ncio_read_vector (lndname, 'CV_ROOF'       , nl_roof, landurban, cvroof   )
      CALL ncio_read_vector (lndname, 'CV_WALL'       , nl_wall, landurban, cvwall   )
      CALL ncio_read_vector (lndname, 'CV_IMPROAD'    , nl_soil, landurban, cvimproad)
      CALL ncio_read_vector (lndname, 'TK_ROOF'       , nl_roof, landurban, tkroof   )
      CALL ncio_read_vector (lndname, 'TK_WALL'       , nl_wall, landurban, tkwall   )
      CALL ncio_read_vector (lndname, 'TK_IMPROAD'    , nl_roof, landurban, tkimproad)

!TODO: add point case
!#ifdef USE_POINT_DATA
!#ifdef USE_OBS_PARA

      ! lndname = trim(dir_atmdata)//'/'//trim(nam_atmdata)
      ! print*, lndname
      ! CALL nccheck( nf90_open(lndname, nf90_nowrite, ncid) )

      ! CALL nccheck( nf90_inq_varid(ncid, "impervious_area_fraction" , wtimproad_vid    ) ) ! imperivous area fraciton
      ! CALL nccheck( nf90_inq_varid(ncid, "tree_area_fraction"       , urbantreepct_vid ) ) ! tree area fraction
      ! CALL nccheck( nf90_inq_varid(ncid, "water_area_fraction"      , urbanwaterpct_vid) ) ! water area fraction
      ! CALL nccheck( nf90_inq_varid(ncid, "roof_area_fraction"       , wtroof_vid       ) ) ! roof area fraction
      ! CALL nccheck( nf90_inq_varid(ncid, "building_mean_height"     , htroof_vid       ) ) ! building mean height
      ! CALL nccheck( nf90_inq_varid(ncid, "tree_mean_height"         , urbantreetop_vid ) ) ! tree mean height
      ! CALL nccheck( nf90_inq_varid(ncid, "canyon_height_width_ratio", canyonhwr_vid    ) ) ! H2W

      ! CALL nccheck( nf90_get_var(ncid, wtimproad_vid,     prwt      ) )
      ! CALL nccheck( nf90_get_var(ncid, urbantreepct_vid,  tpct      ) )
      ! CALL nccheck( nf90_get_var(ncid, urbanwaterpct_vid, wpct      ) )
      ! CALL nccheck( nf90_get_var(ncid, wtroof_vid,        rfwt      ) )
      ! CALL nccheck( nf90_get_var(ncid, htroof_vid,        rfht      ) )
      ! CALL nccheck( nf90_get_var(ncid, urbantreetop_vid,  htop_point) )
      ! CALL nccheck( nf90_get_var(ncid, canyonhwr_vid,     hw_point  ) )

      ! CALL nccheck( nf90_close(ncid) )

      ! wtperroad    (1,1,:) = 1 - (prwt-rfwt)/(1-rfwt-wpct) !1. - prwt
      ! urbantreepct (1,1,:) = tpct*100
      ! urbanwaterpct(1,1,:) = wpct*100
      ! wtroof       (1,1,:) = rfwt
      ! htroof       (1,1,:) = rfht
      ! urbantreetop (1,1,:) = htop_point
      ! canyonhwr    (1,1,:) = hw_point
      !albroof(:,:,:,:,:) = 0.4
!#endif
!#endif
#endif

      !TODO: duplication
      !TODO: 2005 -> year, 变量区分time变化和时间不变量
      ! write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/urban/2005/POP.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'POP_DEN'       , landurban, popden        )
      ! write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/urban/2005/LUCY_country_id.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'LUCY_id'       , landurban, urbanlucy     )
      ! write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/urban/2005/WT_ROOF.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'WT_ROOF'       , landurban, wtroof        )
      ! write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_landdata)//'/urban/2005/HT_ROOF.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'HT_ROOF'       , landurban, htroof        )

      lndname = trim(dir_landdata)//'/urban/2005/PCT_Water.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'PCT_Water'     , landurban, urbanwaterpct )

      lndname = trim(dir_landdata)//'/urban/2005/PCT_Tree.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'PCT_Tree'      , landurban, urbantreepct  )

      lndname = trim(dir_landdata)//'/urban/2005/htop_urb.nc'
      print*, lndname
      CALL ncio_read_vector (lndname, 'URBAN_TREE_TOP', landurban, urbantreetop  )

      lndname = trim("/stu01/dongwz/data/CLMrawdata/urban_5x5/LUCY_rawdata.nc")
      print*, lndname
      CALL ncio_read_bcast_serial (lndname,  "vehicle"    , lvehicle     )
      CALL ncio_read_bcast_serial (lndname,  "weekendday" , lweek_holiday)
      CALL ncio_read_bcast_serial (lndname,  "weekendhour", lweh_prof    )
      CALL ncio_read_bcast_serial (lndname,  "weekdayhour", lwdh_prof    )
      CALL ncio_read_bcast_serial (lndname,  "metabolism" , lhum_prof    )
      CALL ncio_read_bcast_serial (lndname,  "holiday"    , lfix_holiday )

      print*, numurban
      patch2urb = 1
      IF (p_is_worker) THEN
         DO i = 1, numpatch
            m = patchclass(i)
            IF (landpatch%settyp(i) == 13) THEN
                  lucy_id = urbanlucy (u)
                  u = patch2urb

                  IF (lucy_id > 0) THEN
                     vehicle     (:,u) = lvehicle     (lucy_id,:)
                     week_holiday(:,u) = lweek_holiday(lucy_id,:)
                     weh_prof    (:,u) = lweh_prof    (lucy_id,:)
                     wdh_prof    (:,u) = lwdh_prof    (lucy_id,:)
                     hum_prof    (:,u) = lhum_prof    (lucy_id,:)
                     fix_holiday (:,u) = lfix_holiday (lucy_id,:)
                  ENDIF

                  pop_den(u)      = popden        (u) !population density [person/km2]
                  froof(u)        = wtroof        (u) !roof fractional cover
                  hroof(u)        = htroof        (u) !average building height

#ifndef USE_LCZ
                  hwr(u)          = canyonhwr     (u) !average building height to their distance
                  fgper(u)        = wtperroad     (u) !pervious fraction to ground area

                  alb_roof(:,:,u) = albroof   (u,:,:) !albedo of roof
                  alb_wall(:,:,u) = albwall   (u,:,:) !albedo of walls
                  alb_gimp(:,:,u) = albimproad(u,:,:) !albedo of impervious
                  alb_gper(:,:,u) = albperroad(u,:,:) !albedo of pervious road

                  em_roof(u)      = emroof        (u) !emissivity of roof
                  em_wall(u)      = emwall        (u) !emissiviry of wall
                  em_gimp(u)      = emimproad     (u) !emissivity of impervious
                  em_gper(u)      = emperroad     (u) !emissivity of pervious

                  cv_roof(:,u)    = cvroof      (u,:) !heat capacity of roof [J/(m2 K)]
                  cv_wall(:,u)    = cvwall      (u,:) !heat capacity of wall [J/(m2 K)]
                  cv_gimp(:,u)    = cvimproad   (u,:) !heat capacity of impervious [J/(m2 K)]

                  tk_roof(:,u)    = tkroof      (u,:) !thermal conductivity of roof [W/m-K]
                  tk_wall(:,u)    = tkwall      (u,:) !thermal conductivity of wall [W/m-K]
                  tk_gimp(:,u)    = tkimproad   (u,:) !thermal conductivity of impervious [W/m-K]

                  thick_roof      = thickroof     (u) !thickness of roof [m]
                  thick_wall      = thickwall     (u) !thickness of wall [m]

#ifdef URBAN_BEM
                  t_roommax(u)    = tbuildingmax  (u) !maximum temperature of inner room [K]
                  t_roommin(u)    = tbuildingmin  (u) !minimum temperature of inner room [K]
#else
                  t_roommax(u)    = 373.16            !maximum temperature of inner room [K]
                  t_roommin(u)    = 180.00            !minimum temperature of inner room [K]
#endif
#else
                  ! read in LCZ constants
                  hwr  (u) = canyonhwr_lcz (landurban%settyp(u)) !average building height to their distance
                  fgper(u) = wtperroad_lcz (landurban%settyp(u)) !pervious fraction to ground area

                  DO ns = 1,2
                      DO nr = 1,2
                         alb_roof(ns,nr,u) = albroof_lcz   (landurban%settyp(u)) !albedo of roof
                         alb_wall(ns,nr,u) = albwall_lcz   (landurban%settyp(u)) !albedo of walls
                         alb_gimp(ns,nr,u) = albimproad_lcz(landurban%settyp(u)) !albedo of impervious
                         alb_gper(ns,nr,u) = albperroad_lcz(landurban%settyp(u)) !albedo of pervious road
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
                     tk_gimp(:,u) = tkperroad_lcz (landurban%settyp(u)) !thermal conductivity of impervious [W/m-K]
                  ENDDO

                  thick_roof = thickroof_lcz (landurban%settyp(u)) !thickness of roof [m]
                  thick_wall = thickwall_lcz (landurban%settyp(u)) !thickness of wall [m]

#ifdef URBAN_BEM
                  t_roommax(u) = 297.65 !tbuildingmax  (landurban%settyp(u)) !maximum temperature of inner room [K]
                  t_roommin(u) = 290.65 !tbuildingmin  (landurban%settyp(u)) !minimum temperature of inner room [K]
#else
                  t_roommax(u) = 373.16                !maximum temperature of inner room [K]
                  t_roommin(u) = 180.00                !minimum temperature of inner room [K]
#endif
#endif

#ifdef URBAN_WATER
                  flake(u) = urbanwaterpct(u)/100. !urban water fractional cover
#else
                  flake(u) = 0.
#endif

#ifdef URBAN_TREE
                  ! set tree fractional cover (<= 1.-froof)
                  fveg_urb(u) = urbantreepct(u)/100. !urban tree percent
                  IF (flake(u) < 1.) THEN
                     fveg_urb(u) = fveg_urb(u)/(1.-flake(u))
                  ELSE
                     fveg_urb(u) = 0.
                  ENDIF
                  ! Assuming that the tree coverage is less than or equal
                  ! to the ground cover (assume no trees on the roof)
                  fveg_urb(u) = min(fveg_urb(u), 1.-froof(u))
                  fveg    (i) = fveg_urb(u)
#else
                  fveg_urb(u) = 0.
#endif

                  ! set urban tree crown top and bottom [m]
                  htop_urb(u) = min(hroof(u), urbantreetop(u))
                  htop_urb(u) = max(2., htop_urb(u))

                  hbot_urb(u) = htop_urb(u)*hbot0(m)/htop0(m)
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
                  patch2urb = patch2urb+1
            ENDIF
         ENDDO
      ENDIF
      print*, patch2urb
      IF (p_is_worker) THEN
         IF (allocated(lvehicle     )) deallocate ( lvehicle      )
         IF (allocated(lwdh_prof    )) deallocate ( lwdh_prof     )
         IF (allocated(lweh_prof    )) deallocate ( lweh_prof     )
         IF (allocated(lhum_prof    )) deallocate ( lhum_prof     )
         IF (allocated(lweek_holiday)) deallocate ( lweek_holiday )
         IF (allocated(lfix_holiday )) deallocate ( lfix_holiday  )
         IF (allocated(wtroof       )) deallocate ( wtroof        )
         IF (allocated(htroof       )) deallocate ( htroof        )
         IF (allocated(canyonhwr    )) deallocate ( canyonhwr     )
         IF (allocated(wtperroad    )) deallocate ( wtperroad     )
         IF (allocated(urbanwaterpct)) deallocate ( urbanwaterpct )
         IF (allocated(urbantreetop )) deallocate ( urbantreetop  )
         IF (allocated(albroof      )) deallocate ( albroof       )
         IF (allocated(albwall      )) deallocate ( albwall       )
         IF (allocated(albimproad   )) deallocate ( albimproad    )
         IF (allocated(albperroad   )) deallocate ( albperroad    )
         IF (allocated(emroof       )) deallocate ( emroof        )
         IF (allocated(emwall       )) deallocate ( emwall        )
         IF (allocated(emimproad    )) deallocate ( emimproad     )
         IF (allocated(emperroad    )) deallocate ( emperroad     )
         IF (allocated(cvroof       )) deallocate ( cvroof        )
         IF (allocated(cvwall       )) deallocate ( cvwall        )
         IF (allocated(cvimproad    )) deallocate ( cvimproad     )
         IF (allocated(tkroof       )) deallocate ( tkroof        )
         IF (allocated(tkwall       )) deallocate ( tkwall        )
         IF (allocated(tkimproad    )) deallocate ( tkimproad     )
         IF (allocated(tbuildingmax )) deallocate ( tbuildingmax  )
         IF (allocated(tbuildingmin )) deallocate ( tbuildingmin  )
      ENDIF
END SUBROUTINE Urban_readin
