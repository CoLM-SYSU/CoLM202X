#include <define.h>

#if (defined URBAN_MODEL)
MODULE MOD_UrbanTimeVars

! -------------------------------
! Created by Hua Yuan, 12/2020
! -------------------------------

   USE precision
   IMPLICIT NONE
   SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run

   REAL(r8), allocatable :: fwsun          (:) !sunlit fraction of walls [-]
   REAL(r8), allocatable :: dfwsun         (:) !change of sunlit fraction of walls [-]

   ! shortwave absorption
   REAL(r8), allocatable :: sroof      (:,:,:) !roof aborption [-]
   REAL(r8), allocatable :: swsun      (:,:,:) !sunlit wall absorption [-]
   REAL(r8), allocatable :: swsha      (:,:,:) !shaded wall absorption [-]
   REAL(r8), allocatable :: sgimp      (:,:,:) !impervious absorptioin [-]
   REAL(r8), allocatable :: sgper      (:,:,:) !pervious absorptioin [-]
   REAL(r8), allocatable :: slake      (:,:,:) !urban lake absorptioin [-]

   ! net longwave radiation for last time temperature change
   REAL(r8), allocatable :: lwsun          (:) !net longwave of sunlit wall [W/m2]
   REAL(r8), allocatable :: lwsha          (:) !net longwave of shaded wall [W/m2]
   REAL(r8), allocatable :: lgimp          (:) !net longwave of impervious  [W/m2]
   REAL(r8), allocatable :: lgper          (:) !net longwave of pervious [W/m2]
   REAL(r8), allocatable :: lveg           (:) !net longwave of vegetation [W/m2]

   REAL(r8), allocatable :: z_sno_roof   (:,:) !node depth of roof [m]
   REAL(r8), allocatable :: z_sno_gimp   (:,:) !node depth of impervious [m]
   REAL(r8), allocatable :: z_sno_gper   (:,:) !node depth pervious [m]
   REAL(r8), allocatable :: z_sno_lake   (:,:) !node depth lake [m]

   REAL(r8), allocatable :: dz_sno_roof  (:,:) !interface depth of roof [m]
   REAL(r8), allocatable :: dz_sno_gimp  (:,:) !interface depth of impervious [m]
   REAL(r8), allocatable :: dz_sno_gper  (:,:) !interface depth pervious [m]
   REAL(r8), allocatable :: dz_sno_lake  (:,:) !interface depth lake [m]

   REAL(r8), allocatable :: troof_inner    (:) !temperature of roof [K]
   REAL(r8), allocatable :: twsun_inner    (:) !temperature of sunlit wall [K]
   REAL(r8), allocatable :: twsha_inner    (:) !temperature of shaded wall [K]

   REAL(r8), allocatable :: t_roofsno    (:,:) !temperature of roof [K]
   REAL(r8), allocatable :: t_wallsun    (:,:) !temperature of sunlit wall [K]
   REAL(r8), allocatable :: t_wallsha    (:,:) !temperature of shaded wall [K]
   REAL(r8), allocatable :: t_gimpsno    (:,:) !temperature of impervious [K]
   REAL(r8), allocatable :: t_gpersno    (:,:) !temperature of pervious [K]
   REAL(r8), allocatable :: t_lakesno    (:,:) !temperature of pervious [K]

   REAL(r8), allocatable :: wliq_roofsno (:,:) !liquid water in layers [kg/m2]
   REAL(r8), allocatable :: wliq_gimpsno (:,:) !liquid water in layers [kg/m2]
   REAL(r8), allocatable :: wliq_gpersno (:,:) !liquid water in layers [kg/m2]
   REAL(r8), allocatable :: wliq_lakesno (:,:) !liquid water in layers [kg/m2]
   REAL(r8), allocatable :: wice_roofsno (:,:) !ice lens in layers [kg/m2]
   REAL(r8), allocatable :: wice_gimpsno (:,:) !ice lens in layers [kg/m2]
   REAL(r8), allocatable :: wice_gpersno (:,:) !ice lens in layers [kg/m2]
   REAL(r8), allocatable :: wice_lakesno (:,:) !ice lens in layers [kg/m2]

   REAL(r8), allocatable :: sag_roof       (:) !roof snow age [-]
   REAL(r8), allocatable :: sag_gimp       (:) !impervious ground snow age [-]
   REAL(r8), allocatable :: sag_gper       (:) !pervious ground snow age [-]
   REAL(r8), allocatable :: sag_lake       (:) !urban lake snow age [-]

   REAL(r8), allocatable :: scv_roof       (:) !roof snow cover [-]
   REAL(r8), allocatable :: scv_gimp       (:) !impervious ground snow cover [-]
   REAL(r8), allocatable :: scv_gper       (:) !pervious ground snow cover [-]
   REAL(r8), allocatable :: scv_lake       (:) !urban lake snow cover [-]

   REAL(r8), allocatable :: fsno_roof      (:) !roof snow fraction [-]
   REAL(r8), allocatable :: fsno_gimp      (:) !impervious ground snow fraction [-]
   REAL(r8), allocatable :: fsno_gper      (:) !pervious ground snow fraction [-]
   REAL(r8), allocatable :: fsno_lake      (:) !urban lake snow fraction [-]

   REAL(r8), allocatable :: snowdp_roof    (:) !roof snow depth [m]
   REAL(r8), allocatable :: snowdp_gimp    (:) !impervious ground snow depth [m]
   REAL(r8), allocatable :: snowdp_gper    (:) !pervious ground snow depth [m]
   REAL(r8), allocatable :: snowdp_lake    (:) !urban lake snow depth [m]

   REAL(r8), allocatable :: t_room         (:) !temperature of inner building [K]
   !TODO: rename the below variables
   REAL(r8), allocatable :: tafu           (:) !temperature of outer building [K]
   REAL(r8), allocatable :: Fhac           (:) !sensible flux from heat or cool AC [W/m2]
   REAL(r8), allocatable :: Fwst           (:) !waste heat flux from heat or cool AC [W/m2]
   REAL(r8), allocatable :: Fach           (:) !flux from inner and outter air exchange [W/m2]
   REAL(r8), allocatable :: Fahe           (:) !flux from metabolism and vehicle [W/m2]
   REAL(r8), allocatable :: Fhah           (:) !sensible heat flux from heating [W/m2]
   REAL(r8), allocatable :: vehc           (:) !flux from vehicle [W/m2]
   REAL(r8), allocatable :: meta           (:) !flux from metabolism [W/m2]

   REAL(r8), allocatable :: fsen_roof      (:) !sensible heat flux from roof [W/m2]
   REAL(r8), allocatable :: fsen_wsun      (:) !sensible heat flux from sunlit wall [W/m2]
   REAL(r8), allocatable :: fsen_wsha      (:) !sensible heat flux from shaded wall [W/m2]
   REAL(r8), allocatable :: fsen_gimp      (:) !sensible heat flux from impervious road [W/m2]
   REAL(r8), allocatable :: fsen_gper      (:) !sensible heat flux from pervious road [W/m2]
   REAL(r8), allocatable :: fsen_urbl      (:) !sensible heat flux from urban vegetation [W/m2]

   REAL(r8), allocatable :: lfevp_roof     (:) !latent heat flux from roof [W/m2]
   REAL(r8), allocatable :: lfevp_gimp     (:) !latent heat flux from impervious road [W/m2]
   REAL(r8), allocatable :: lfevp_gper     (:) !latent heat flux from pervious road [W/m2]
   REAL(r8), allocatable :: lfevp_urbl     (:) !latent heat flux from urban vegetation [W/m2]

   REAL(r8), allocatable :: troof          (:) !temperature of roof [K]
   REAL(r8), allocatable :: twall          (:) !temperature of wall [K]

   REAL(r8), allocatable :: urb_green      (:) !fractional of green leaf in urban patch [-]
   REAL(r8), allocatable :: urb_lai        (:) !urban tree LAI [m2/m2]
   REAL(r8), allocatable :: urb_sai        (:) !urban tree SAI [m2/m2]


! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_UrbanTimeVars
   PUBLIC :: deallocate_UrbanTimeVars
   PUBLIC :: READ_UrbanTimeVars
   PUBLIC :: WRITE_UrbanTimeVars

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_UrbanTimeVars ()
! ------------------------------------------------------
! Allocates memory for CLM 1d [numurban] variables
! ------------------------------------------------------
      USE precision
      USE spmd_task
      USE mod_landurban
      USE GlobalVars
      IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numurban > 0) THEN
            allocate (fwsun                         (numurban))
            allocate (dfwsun                        (numurban))

            allocate (sroof                     (2,2,numurban))
            allocate (swsun                     (2,2,numurban))
            allocate (swsha                     (2,2,numurban))
            allocate (sgimp                     (2,2,numurban))
            allocate (sgper                     (2,2,numurban))
            allocate (slake                     (2,2,numurban))

            allocate (lwsun                         (numurban))
            allocate (lwsha                         (numurban))
            allocate (lgimp                         (numurban))
            allocate (lgper                         (numurban))
            allocate (lveg                          (numurban))

            allocate (z_sno_roof         (maxsnl+1:0,numurban))
            allocate (z_sno_gimp         (maxsnl+1:0,numurban))
            allocate (z_sno_gper         (maxsnl+1:0,numurban))
            allocate (z_sno_lake         (maxsnl+1:0,numurban))

            allocate (dz_sno_roof        (maxsnl+1:0,numurban))
            allocate (dz_sno_gimp        (maxsnl+1:0,numurban))
            allocate (dz_sno_gper        (maxsnl+1:0,numurban))
            allocate (dz_sno_lake        (maxsnl+1:0,numurban))

            allocate (t_roofsno    (maxsnl+1:nl_roof,numurban))
            allocate (t_wallsun    (maxsnl+1:nl_wall,numurban))
            allocate (t_wallsha    (maxsnl+1:nl_wall,numurban))
            allocate (t_gimpsno    (maxsnl+1:nl_soil,numurban))
            allocate (t_gpersno    (maxsnl+1:nl_soil,numurban))
            allocate (t_lakesno    (maxsnl+1:nl_soil,numurban))

            allocate (troof_inner                   (numurban))
            allocate (twsun_inner                   (numurban))
            allocate (twsha_inner                   (numurban))

            allocate (wliq_roofsno (maxsnl+1:nl_roof,numurban))
            allocate (wice_roofsno (maxsnl+1:nl_roof,numurban))
            allocate (wliq_gimpsno (maxsnl+1:nl_soil,numurban))
            allocate (wice_gimpsno (maxsnl+1:nl_soil,numurban))
            allocate (wliq_gpersno (maxsnl+1:nl_soil,numurban))
            allocate (wice_gpersno (maxsnl+1:nl_soil,numurban))
            allocate (wliq_lakesno (maxsnl+1:nl_soil,numurban))
            allocate (wice_lakesno (maxsnl+1:nl_soil,numurban))

            allocate (sag_roof                      (numurban))
            allocate (sag_gimp                      (numurban))
            allocate (sag_gper                      (numurban))
            allocate (sag_lake                      (numurban))
            allocate (scv_roof                      (numurban))
            allocate (scv_gimp                      (numurban))
            allocate (scv_gper                      (numurban))
            allocate (scv_lake                      (numurban))
            allocate (fsno_roof                     (numurban))
            allocate (fsno_gimp                     (numurban))
            allocate (fsno_gper                     (numurban))
            allocate (fsno_lake                     (numurban))
            allocate (snowdp_roof                   (numurban))
            allocate (snowdp_gimp                   (numurban))
            allocate (snowdp_gper                   (numurban))
            allocate (snowdp_lake                   (numurban))

            allocate (t_room                        (numurban))
            allocate (tafu                          (numurban))
            allocate (Fhac                          (numurban))
            allocate (Fwst                          (numurban))
            allocate (Fach                          (numurban))
            allocate (Fahe                          (numurban))
            allocate (Fhah                          (numurban))
            allocate (vehc                          (numurban))
            allocate (meta                          (numurban))

            allocate (fsen_roof                     (numurban))
            allocate (fsen_wsun                     (numurban))
            allocate (fsen_wsha                     (numurban))
            allocate (fsen_gimp                     (numurban))
            allocate (fsen_gper                     (numurban))
            allocate (fsen_urbl                     (numurban))

            allocate (lfevp_roof                    (numurban))
            allocate (lfevp_gimp                    (numurban))
            allocate (lfevp_gper                    (numurban))
            allocate (lfevp_urbl                    (numurban))

            allocate (troof                         (numurban))
            allocate (twall                         (numurban))

            allocate (urb_green                     (numurban))
            allocate (urb_lai                       (numurban))
            allocate (urb_sai                       (numurban))
         ENDIF
      ENDIF
   END SUBROUTINE allocate_UrbanTimeVars

   SUBROUTINE READ_UrbanTimeVars (file_restart)

      USE ncio_vector
      USE mod_landurban
      USE GlobalVars

      IMPLICIT NONE

      character(LEN=*), intent(in) :: file_restart

      call ncio_read_vector (file_restart, 'fwsun' , landurban, fwsun ) !
      call ncio_read_vector (file_restart, 'dfwsun', landurban, dfwsun) !

      call ncio_read_vector (file_restart, 'sroof', 2, 2, landurban, sroof) !
      call ncio_read_vector (file_restart, 'swsun', 2, 2, landurban, swsun) !
      call ncio_read_vector (file_restart, 'swsha', 2, 2, landurban, swsha) !
      call ncio_read_vector (file_restart, 'sgimp', 2, 2, landurban, sgimp) !
      call ncio_read_vector (file_restart, 'sgper', 2, 2, landurban, sgper) !
      call ncio_read_vector (file_restart, 'slake', 2, 2, landurban, slake) !

      call ncio_read_vector (file_restart, 'lwsun', landurban, lwsun) !
      call ncio_read_vector (file_restart, 'lwsha', landurban, lwsha) !
      call ncio_read_vector (file_restart, 'lgimp', landurban, lgimp) !
      call ncio_read_vector (file_restart, 'lgper', landurban, lgper) !

      call ncio_read_vector (file_restart, 'z_sno_roof' , -maxsnl, landurban, z_sno_roof ) !
      call ncio_read_vector (file_restart, 'z_sno_gimp' , -maxsnl, landurban, z_sno_gimp ) !
      call ncio_read_vector (file_restart, 'z_sno_gper' , -maxsnl, landurban, z_sno_gper ) !
      call ncio_read_vector (file_restart, 'z_sno_lake' , -maxsnl, landurban, z_sno_lake ) !

      call ncio_read_vector (file_restart, 'dz_sno_roof', -maxsnl, landurban, dz_sno_roof) !
      call ncio_read_vector (file_restart, 'dz_sno_gimp', -maxsnl, landurban, dz_sno_gimp) !
      call ncio_read_vector (file_restart, 'dz_sno_gper', -maxsnl, landurban, dz_sno_gper) !
      call ncio_read_vector (file_restart, 'dz_sno_lake', -maxsnl, landurban, dz_sno_lake) !

      call ncio_read_vector (file_restart, 't_roofsno', nl_roof-maxsnl, landurban, t_roofsno) !
      call ncio_read_vector (file_restart, 't_wallsun', nl_wall-maxsnl, landurban, t_wallsun) !
      call ncio_read_vector (file_restart, 't_wallsha', nl_wall-maxsnl, landurban, t_wallsha) !
      call ncio_read_vector (file_restart, 't_gimpsno', nl_soil-maxsnl, landurban, t_gimpsno) !
      call ncio_read_vector (file_restart, 't_gpersno', nl_soil-maxsnl, landurban, t_gpersno) !
      call ncio_read_vector (file_restart, 't_lakesno', nl_soil-maxsnl, landurban, t_lakesno) !

      call ncio_read_vector (file_restart, 'wliq_roofsno', nl_roof-maxsnl, landurban, wliq_roofsno) !
      call ncio_read_vector (file_restart, 'wliq_gimpsno', nl_soil-maxsnl, landurban, wliq_gimpsno) !
      call ncio_read_vector (file_restart, 'wliq_gpersno', nl_soil-maxsnl, landurban, wliq_gpersno) !
      call ncio_read_vector (file_restart, 'wliq_lakesno', nl_soil-maxsnl, landurban, wliq_lakesno) !
      call ncio_read_vector (file_restart, 'wice_roofsno', nl_roof-maxsnl, landurban, wice_roofsno) !
      call ncio_read_vector (file_restart, 'wice_gimpsno', nl_soil-maxsnl, landurban, wice_gimpsno) !
      call ncio_read_vector (file_restart, 'wice_gpersno', nl_soil-maxsnl, landurban, wice_gpersno) !
      call ncio_read_vector (file_restart, 'wice_lakesno', nl_soil-maxsnl, landurban, wice_lakesno) !

      call ncio_read_vector (file_restart, 'sag_roof'   , landurban, sag_roof   ) !
      call ncio_read_vector (file_restart, 'sag_gimp'   , landurban, sag_gimp   ) !
      call ncio_read_vector (file_restart, 'sag_gper'   , landurban, sag_gper   ) !
      call ncio_read_vector (file_restart, 'scv_roof'   , landurban, scv_roof   ) !
      call ncio_read_vector (file_restart, 'scv_gimp'   , landurban, scv_gimp   ) !
      call ncio_read_vector (file_restart, 'scv_gper'   , landurban, scv_gper   ) !
      call ncio_read_vector (file_restart, 'scv_lake'   , landurban, scv_lake   ) !
      call ncio_read_vector (file_restart, 'fsno_roof'  , landurban, fsno_roof  ) !
      call ncio_read_vector (file_restart, 'fsno_gimp'  , landurban, fsno_gimp  ) !
      call ncio_read_vector (file_restart, 'fsno_gper'  , landurban, fsno_gper  ) !
      call ncio_read_vector (file_restart, 'fsno_lake'  , landurban, fsno_lake  ) !
      call ncio_read_vector (file_restart, 'snowdp_roof', landurban, snowdp_roof) !
      call ncio_read_vector (file_restart, 'snowdp_gimp', landurban, snowdp_gimp) !
      call ncio_read_vector (file_restart, 'snowdp_gper', landurban, snowdp_gper) !
      call ncio_read_vector (file_restart, 'snowdp_lake', landurban, snowdp_lake) !
      call ncio_read_vector (file_restart, 't_room '    , landurban, t_room     ) !
      call ncio_read_vector (file_restart, 'tafu'       , landurban, tafu       ) !
      call ncio_read_vector (file_restart, 'Fhac'       , landurban, Fhac       ) !
      call ncio_read_vector (file_restart, 'Fwst'       , landurban, Fwst       ) !
      call ncio_read_vector (file_restart, 'Fach'       , landurban, Fach       ) !
      call ncio_read_vector (file_restart, 'Fahe'       , landurban, Fahe       ) !
      call ncio_read_vector (file_restart, 'Fhah'       , landurban, Fhah       ) !
      call ncio_read_vector (file_restart, 'vehc'       , landurban, vehc       ) !
      call ncio_read_vector (file_restart, 'meta'       , landurban, meta       ) !
      call ncio_read_vector (file_restart, 'fsen_roof'  , landurban, fsen_roof  ) !
      call ncio_read_vector (file_restart, 'fsen_wsun'  , landurban, fsen_wsun  ) !
      call ncio_read_vector (file_restart, 'fsen_wsha'  , landurban, fsen_wsha  ) !
      call ncio_read_vector (file_restart, 'fsen_gimp'  , landurban, fsen_gimp  ) !
      call ncio_read_vector (file_restart, 'fsen_gper'  , landurban, fsen_gper  ) !
      call ncio_read_vector (file_restart, 'fsen_urbl'  , landurban, fsen_urbl  ) !
      call ncio_read_vector (file_restart, 'lfevp_roof' , landurban, lfevp_roof ) !
      call ncio_read_vector (file_restart, 'lfevp_gimp' , landurban, lfevp_gimp ) !
      call ncio_read_vector (file_restart, 'lfevp_gper' , landurban, lfevp_gper ) !
      call ncio_read_vector (file_restart, 'lfevp_urbl' , landurban, lfevp_urbl ) !
      call ncio_read_vector (file_restart, 'troof'      , landurban, troof      ) !
      call ncio_read_vector (file_restart, 'twall'      , landurban, twall      ) !
      call ncio_read_vector (file_restart, 'urb_green'  , landurban, urb_green  ) !
      call ncio_read_vector (file_restart, 'tree_lai'   , landurban, urb_lai    ) !
      call ncio_read_vector (file_restart, 'tree_sai'   , landurban, urb_sai    ) !

   END SUBROUTINE READ_UrbanTimeVars

   SUBROUTINE WRITE_UrbanTimeVars (file_restart)

      USE mod_namelist, only : DEF_REST_COMPRESS_LEVEL
      USE mod_landurban
      USE ncio_vector
      USE GlobalVars
      IMPLICIT NONE

      character(LEN=*), intent(in) :: file_restart

      ! Local variables
      integer :: compress

      compress = DEF_REST_COMPRESS_LEVEL

      call ncio_create_file_vector (file_restart, landurban)
      CALL ncio_define_dimension_vector (file_restart, landurban, 'urban')

      CALL ncio_define_dimension_vector (file_restart, landurban, 'snow'    , -maxsnl       )
      CALL ncio_define_dimension_vector (file_restart, landurban, 'soil'    , nl_soil       )
      CALL ncio_define_dimension_vector (file_restart, landurban, 'roof'    , nl_roof       )
      CALL ncio_define_dimension_vector (file_restart, landurban, 'wall'    , nl_wall       )
      CALL ncio_define_dimension_vector (file_restart, landurban, 'soilsnow', nl_soil-maxsnl)
      CALL ncio_define_dimension_vector (file_restart, landurban, 'roofsnow', nl_roof-maxsnl)
      CALL ncio_define_dimension_vector (file_restart, landurban, 'wallsnow', nl_wall-maxsnl)


      CALL ncio_define_dimension_vector (file_restart, landurban, 'band', 2)
      CALL ncio_define_dimension_vector (file_restart, landurban, 'rtyp', 2)


      call ncio_write_vector (file_restart, 'fwsun' , 'urban', landurban, fwsun , compress) !
      call ncio_write_vector (file_restart, 'dfwsun', 'urban', landurban, dfwsun, compress) !

      call ncio_write_vector (file_restart, 'sroof', 'band', 2, 'rtyp', 2, 'urban', landurban, sroof, compress) !
      call ncio_write_vector (file_restart, 'swsun', 'band', 2, 'rtyp', 2, 'urban', landurban, swsun, compress) !
      call ncio_write_vector (file_restart, 'swsha', 'band', 2, 'rtyp', 2, 'urban', landurban, swsha, compress) !
      call ncio_write_vector (file_restart, 'sgimp', 'band', 2, 'rtyp', 2, 'urban', landurban, sgimp, compress) !
      call ncio_write_vector (file_restart, 'sgper', 'band', 2, 'rtyp', 2, 'urban', landurban, sgper, compress) !
      call ncio_write_vector (file_restart, 'slake', 'band', 2, 'rtyp', 2, 'urban', landurban, slake, compress) !

      call ncio_write_vector (file_restart, 'lwsun', 'urban', landurban, lwsun, compress) !
      call ncio_write_vector (file_restart, 'lwsha', 'urban', landurban, lwsha, compress) !
      call ncio_write_vector (file_restart, 'lgimp', 'urban', landurban, lgimp, compress) !
      call ncio_write_vector (file_restart, 'lgper', 'urban', landurban, lgper, compress) !

      call ncio_write_vector (file_restart, 'z_sno_roof' , 'snow', -maxsnl, 'urban', landurban, z_sno_roof , compress) !
      call ncio_write_vector (file_restart, 'z_sno_gimp' , 'snow', -maxsnl, 'urban', landurban, z_sno_gimp , compress) !
      call ncio_write_vector (file_restart, 'z_sno_gper' , 'snow', -maxsnl, 'urban', landurban, z_sno_gper , compress) !
      call ncio_write_vector (file_restart, 'z_sno_lake' , 'snow', -maxsnl, 'urban', landurban, z_sno_lake , compress) !

      call ncio_write_vector (file_restart, 'dz_sno_roof', 'snow', -maxsnl, 'urban', landurban, dz_sno_roof, compress) !
      call ncio_write_vector (file_restart, 'dz_sno_gimp', 'snow', -maxsnl, 'urban', landurban, dz_sno_gimp, compress) !
      call ncio_write_vector (file_restart, 'dz_sno_gper', 'snow', -maxsnl, 'urban', landurban, dz_sno_gper, compress) !
      call ncio_write_vector (file_restart, 'dz_sno_lake', 'snow', -maxsnl, 'urban', landurban, dz_sno_lake, compress) !

      call ncio_write_vector (file_restart, 't_roofsno', 'roofsnow', nl_roof-maxsnl, 'urban', landurban, t_roofsno, compress) !
      call ncio_write_vector (file_restart, 't_wallsun', 'wallsnow', nl_wall-maxsnl, 'urban', landurban, t_wallsun, compress) !
      call ncio_write_vector (file_restart, 't_wallsha', 'wallsnow', nl_wall-maxsnl, 'urban', landurban, t_wallsha, compress) !
      call ncio_write_vector (file_restart, 't_gimpsno', 'soilsnow', nl_soil-maxsnl, 'urban', landurban, t_gimpsno, compress) !
      call ncio_write_vector (file_restart, 't_gpersno', 'soilsnow', nl_soil-maxsnl, 'urban', landurban, t_gpersno, compress) !
      call ncio_write_vector (file_restart, 't_lakesno', 'soilsnow', nl_soil-maxsnl, 'urban', landurban, t_lakesno, compress) !

      call ncio_write_vector (file_restart, 'wliq_roofsno', 'roofsnow', nl_roof-maxsnl, 'urban', landurban, wliq_roofsno, compress) !
      call ncio_write_vector (file_restart, 'wliq_gimpsno', 'soilsnow', nl_soil-maxsnl, 'urban', landurban, wliq_gimpsno, compress) !
      call ncio_write_vector (file_restart, 'wliq_gpersno', 'soilsnow', nl_soil-maxsnl, 'urban', landurban, wliq_gpersno, compress) !
      call ncio_write_vector (file_restart, 'wliq_lakesno', 'soilsnow', nl_soil-maxsnl, 'urban', landurban, wliq_lakesno, compress) !
      call ncio_write_vector (file_restart, 'wice_roofsno', 'roofsnow', nl_roof-maxsnl, 'urban', landurban, wice_roofsno, compress) !
      call ncio_write_vector (file_restart, 'wice_gimpsno', 'soilsnow', nl_soil-maxsnl, 'urban', landurban, wice_gimpsno, compress) !
      call ncio_write_vector (file_restart, 'wice_gpersno', 'soilsnow', nl_soil-maxsnl, 'urban', landurban, wice_gpersno, compress) !
      call ncio_write_vector (file_restart, 'wice_lakesno', 'soilsnow', nl_soil-maxsnl, 'urban', landurban, wice_lakesno, compress) !

      call ncio_write_vector (file_restart, 'sag_roof'   , 'urban', landurban, sag_roof   , compress) !
      call ncio_write_vector (file_restart, 'sag_gimp'   , 'urban', landurban, sag_gimp   , compress) !
      call ncio_write_vector (file_restart, 'sag_gper'   , 'urban', landurban, sag_gper   , compress) !
      call ncio_write_vector (file_restart, 'scv_roof'   , 'urban', landurban, scv_roof   , compress) !
      call ncio_write_vector (file_restart, 'scv_gimp'   , 'urban', landurban, scv_gimp   , compress) !
      call ncio_write_vector (file_restart, 'scv_gper'   , 'urban', landurban, scv_gper   , compress) !
      call ncio_write_vector (file_restart, 'scv_lake'   , 'urban', landurban, scv_lake   , compress) !
      call ncio_write_vector (file_restart, 'fsno_roof'  , 'urban', landurban, fsno_roof  , compress) !
      call ncio_write_vector (file_restart, 'fsno_gimp'  , 'urban', landurban, fsno_gimp  , compress) !
      call ncio_write_vector (file_restart, 'fsno_gper'  , 'urban', landurban, fsno_gper  , compress) !
      call ncio_write_vector (file_restart, 'fsno_lake'  , 'urban', landurban, fsno_lake  , compress) !
      call ncio_write_vector (file_restart, 'snowdp_roof', 'urban', landurban, snowdp_roof, compress) !
      call ncio_write_vector (file_restart, 'snowdp_gimp', 'urban', landurban, snowdp_gimp, compress) !
      call ncio_write_vector (file_restart, 'snowdp_gper', 'urban', landurban, snowdp_gper, compress) !
      call ncio_write_vector (file_restart, 'snowdp_lake', 'urban', landurban, snowdp_lake, compress) !
      call ncio_write_vector (file_restart, 't_room '    , 'urban', landurban, t_room     , compress) !
      call ncio_write_vector (file_restart, 'tafu'       , 'urban', landurban, tafu       , compress) !
      call ncio_write_vector (file_restart, 'Fhac'       , 'urban', landurban, Fhac       , compress) !
      call ncio_write_vector (file_restart, 'Fwst'       , 'urban', landurban, Fwst       , compress) !
      call ncio_write_vector (file_restart, 'Fach'       , 'urban', landurban, Fach       , compress) !
      call ncio_write_vector (file_restart, 'Fahe'       , 'urban', landurban, Fahe       , compress) !
      call ncio_write_vector (file_restart, 'Fhah'       , 'urban', landurban, Fhah       , compress) !
      call ncio_write_vector (file_restart, 'vehc'       , 'urban', landurban, vehc       , compress) !
      call ncio_write_vector (file_restart, 'meta'       , 'urban', landurban, meta       , compress) !
      call ncio_write_vector (file_restart, 'fsen_roof'  , 'urban', landurban, fsen_roof  , compress) !
      call ncio_write_vector (file_restart, 'fsen_wsun'  , 'urban', landurban, fsen_wsun  , compress) !
      call ncio_write_vector (file_restart, 'fsen_wsha'  , 'urban', landurban, fsen_wsha  , compress) !
      call ncio_write_vector (file_restart, 'fsen_gimp'  , 'urban', landurban, fsen_gimp  , compress) !
      call ncio_write_vector (file_restart, 'fsen_gper'  , 'urban', landurban, fsen_gper  , compress) !
      call ncio_write_vector (file_restart, 'fsen_urbl'  , 'urban', landurban, fsen_urbl  , compress) !
      call ncio_write_vector (file_restart, 'lfevp_roof' , 'urban', landurban, lfevp_roof , compress) !
      call ncio_write_vector (file_restart, 'lfevp_gimp' , 'urban', landurban, lfevp_gimp , compress) !
      call ncio_write_vector (file_restart, 'lfevp_gper' , 'urban', landurban, lfevp_gper , compress) !
      call ncio_write_vector (file_restart, 'lfevp_urbl' , 'urban', landurban, lfevp_urbl , compress) !
      call ncio_write_vector (file_restart, 'troof'      , 'urban', landurban, troof      , compress) !
      call ncio_write_vector (file_restart, 'twall'      , 'urban', landurban, twall      , compress) !
      call ncio_write_vector (file_restart, 'urb_green'  , 'urban', landurban, urb_green  , compress) !
      call ncio_write_vector (file_restart, 'tree_lai'   , 'urban', landurban, urb_lai    , compress) !
      call ncio_write_vector (file_restart, 'tree_sai'   , 'urban', landurban, urb_sai    , compress) !

   END SUBROUTINE WRITE_UrbanTimeVars

   SUBROUTINE deallocate_UrbanTimeVars

      USE spmd_task
      USE mod_landurban

      IF (p_is_worker) THEN
         IF (numurban > 0) THEN
            deallocate (fwsun        )
            deallocate (dfwsun       )

            deallocate (sroof        )
            deallocate (swsun        )
            deallocate (swsha        )
            deallocate (sgimp        )
            deallocate (sgper        )
            deallocate (slake        )

            deallocate (lwsun        )
            deallocate (lwsha        )
            deallocate (lgimp        )
            deallocate (lgper        )
            deallocate (lveg         )

            deallocate (z_sno_roof   )
            deallocate (z_sno_gimp   )
            deallocate (z_sno_gper   )
            deallocate (z_sno_lake   )

            deallocate (dz_sno_roof  )
            deallocate (dz_sno_gimp  )
            deallocate (dz_sno_gper  )
            deallocate (dz_sno_lake  )

            deallocate (t_roofsno    )
            deallocate (t_wallsun    )
            deallocate (t_wallsha    )
            deallocate (t_gimpsno    )
            deallocate (t_gpersno    )
            deallocate (t_lakesno    )

            deallocate (troof_inner  )
            deallocate (twsun_inner  )
            deallocate (twsha_inner  )

            deallocate (wliq_roofsno )
            deallocate (wice_roofsno )
            deallocate (wliq_gimpsno )
            deallocate (wice_gimpsno )
            deallocate (wliq_gpersno )
            deallocate (wice_gpersno )
            deallocate (wliq_lakesno )
            deallocate (wice_lakesno )

            deallocate (sag_roof     )
            deallocate (sag_gimp     )
            deallocate (sag_gper     )
            deallocate (sag_lake     )
            deallocate (scv_roof     )
            deallocate (scv_gimp     )
            deallocate (scv_gper     )
            deallocate (scv_lake     )
            deallocate (fsno_roof    )
            deallocate (fsno_gimp    )
            deallocate (fsno_gper    )
            deallocate (fsno_lake    )
            deallocate (snowdp_roof  )
            deallocate (snowdp_gimp  )
            deallocate (snowdp_gper  )
            deallocate (snowdp_lake  )

            deallocate (t_room       )
            deallocate (tafu         )
            deallocate (Fhac         )
            deallocate (Fhah         )
            deallocate (vehc         )
            deallocate (meta         )
            deallocate (Fwst         )
            deallocate (Fach         )
            deallocate (Fahe         )

            deallocate (fsen_roof    )
            deallocate (fsen_wsun    )
            deallocate (fsen_wsha    )
            deallocate (fsen_gimp    )
            deallocate (fsen_gper    )
            deallocate (fsen_urbl    )

            deallocate (lfevp_roof   )
            deallocate (lfevp_gimp   )
            deallocate (lfevp_gper   )
            deallocate (lfevp_urbl   )

            deallocate (troof        )
            deallocate (twall        )

            deallocate (urb_green    )
            deallocate (urb_lai      )
            deallocate (urb_sai      )
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_UrbanTimeVars

END MODULE MOD_UrbanTimeVars
! ---------- EOP ------------
#endif
