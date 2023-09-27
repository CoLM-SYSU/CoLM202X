#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_GRACE

   USE MOD_DataType
   USE MOD_Mapping_Grid2Pset
   IMPLICIT NONE

   PUBLIC :: init_DA_GRACE
   PUBLIC :: do_DA_GRACE
   PUBLIC :: final_DA_GRACE 
   
   REAL(r8), allocatable, PUBLIC :: fslp_patch (:) ! slope factor of subsurface runoff

   PRIVATE
      
   CHARACTER(len=256) :: file_grace
   TYPE(grid_type)    :: grid_grace

   REAL(r8), allocatable :: longrace(:)
   REAL(r8), allocatable :: latgrace(:)

   INTEGER :: nobstime
   INTEGER,  allocatable :: obsyear  (:)
   INTEGER,  allocatable :: obsmonth (:)

   type (mapping_grid2pset_type) :: mg2p_grace
   
   REAL(r8), allocatable :: lwe_obs_this (:)
   REAL(r8), allocatable :: err_obs_this (:)
   
   REAL(r8), allocatable :: lwe_obs_prev (:)
   REAL(r8), allocatable :: err_obs_prev (:)
   
   REAL(r8), allocatable :: wat_prev_m (:)
   REAL(r8), allocatable :: wat_this_m (:)

   REAL(r8), allocatable :: rsub_acc_prev_m (:)
   REAL(r8), allocatable :: rsub_acc_this_m (:)
   
   REAL(r8), allocatable :: rsub_prev_m0 (:)
   REAL(r8), allocatable :: rsub_prev_m1 (:)
   REAL(r8), allocatable :: rsub_this_m  (:)

   LOGICAL :: has_prev_grace_obs
   INTEGER :: nac_grace

CONTAINS

   ! ----------
   SUBROUTINE init_DA_GRACE ()
      
      USE MOD_Spmd_Task
      USE MOD_Namelist, only : DEF_DA_obsdir
      USE MOD_Grid
      USE MOD_NetCDFSerial
      USE MOD_Mesh,     only : numelm
      USE MOD_LandElm,  only : landelm
      USE MOD_LandPatch
#ifdef CROP 
      USE MOD_LandCrop
#endif
      USE MOD_Pixelset
      USE MOD_Mapping_Grid2pset
      USE MOD_RangeCheck
      IMPLICIT NONE
      
      ! Local Variables

      REAL(r8), allocatable :: time_real8(:)
      INTEGER :: itime
      
      file_grace = trim(DEF_DA_obsdir) &
         // '/GRACE_JPL/GRCTellus.JPL.200204_202207.GLO.RL06M.MSCNv02CRI.nc'

      CALL ncio_read_bcast_serial (file_grace, 'time', time_real8)

      nobstime = size(time_real8)
      allocate (obsyear (nobstime))
      allocate (obsmonth(nobstime))

      DO itime = 1, nobstime
         CALL retrieve_yymm_from_days (time_real8(itime), obsyear(itime), obsmonth(itime))
      ENDDO

      IF (p_is_master) THEN
         write(*,*) 'Assimilate GRACE data at'
         DO itime = 1, nobstime
            write(*,*) obsyear(itime), obsmonth(itime)
         ENDDO
      ENDIF

      CALL ncio_read_bcast_serial (file_grace, 'lon', longrace)
      CALL ncio_read_bcast_serial (file_grace, 'lat', latgrace)

      CALL grid_grace%define_by_center (latgrace,longrace)

      call mg2p_grace%build (grid_grace, landelm)

      IF (p_is_worker) THEN
         IF (numelm > 0) THEN
            allocate (lwe_obs_this (numelm))
            allocate (err_obs_this (numelm))
            allocate (lwe_obs_prev (numelm))
            allocate (err_obs_prev (numelm))
         ENDIF

         IF (numpatch > 0) THEN
            allocate (wat_prev_m      (numpatch))
            allocate (wat_this_m      (numpatch))
            allocate (rsub_acc_prev_m (numpatch))
            allocate (rsub_acc_this_m (numpatch))
            allocate (rsub_prev_m0    (numpatch))
            allocate (rsub_prev_m1    (numpatch))
            allocate (rsub_this_m     (numpatch))
            
            allocate (fslp_patch      (numpatch))
         ENDIF
      ENDIF
      
      IF (p_is_worker) THEN
#ifdef CROP 
         CALL elm_patch%build (landelm, landpatch, use_frac = .true., sharedfrac = pctshrpch)
#else
         CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
#endif
      ENDIF
      
      has_prev_grace_obs = .false.
      
      IF (p_is_worker) THEN
         wat_this_m     (:) = 0.
         rsub_acc_this_m(:) = 0.
         rsub_this_m    (:) = 0.

         fslp_patch     (:) = 1.0
      ENDIF

      deallocate (time_real8)

   END SUBROUTINE init_DA_GRACE 

   ! ----------
   SUBROUTINE do_DA_GRACE (idate, deltim)
      
      USE MOD_Spmd_task
      USE MOD_TimeManager
      USE MOD_NetCDFBlock
      USE MOD_Mesh
      USE MOD_LandPatch
      USE MOD_Vars_1DFluxes,      only : rnof, rsur
      USE MOD_Vars_TimeVariables, only : wat, wa
      USE MOD_RangeCheck 
      IMPLICIT NONE
      
      INTEGER,  INTENT(in) :: idate(3)
      REAL(r8), INTENT(in) :: deltim

      ! Local Variables
      LOGICAL :: is_obs_time
      INTEGER :: month, mday, itime, ielm, istt, iend

      REAL(r8) :: w1, w0, r1, r0, var_o, var_m, var_c, dw_f, dw_o, dw_a, rr
      REAL(r8) :: fslp_rsub, fprev
   
      TYPE(block_data_real8_2d) :: f_grace_lwe ! unit: cm
      TYPE(block_data_real8_2d) :: f_grace_err ! unit: cm
      
      CALL julian2monthday (idate(1), idate(2), month, mday)
   
      is_obs_time = any((obsyear == idate(1)) .and. (obsmonth == month))

      IF (p_is_master) THEN
         IF (is_obs_time) THEN
            write(*,*) 'GRACE at this time.'
         ENDIF
      ENDIF

      IF (p_is_worker) THEN

         IF (has_prev_grace_obs) THEN
            rsub_acc_prev_m = rsub_acc_prev_m + (rnof - rsur) * deltim
         ENDIF

         IF (is_obs_time) THEN

            nac_grace = nac_grace + 1

            wat_this_m = wat_this_m + wat + wa

            rsub_acc_this_m = rsub_acc_this_m + (rnof - rsur) * deltim

            IF (has_prev_grace_obs) THEN
               rsub_prev_m1 = rsub_prev_m1 + rsub_acc_prev_m
            ENDIF

            rsub_this_m = rsub_this_m + rsub_acc_this_m
         ENDIF

      ENDIF

      IF (is_obs_time .and. (isendofmonth(idate, deltim))) then

         itime = findloc((obsyear == idate(1)) .and. (obsmonth == month), .true., dim=1)
      
         IF (p_is_io) THEN
            CALL allocate_block_data (grid_grace, f_grace_lwe)  
            CALL allocate_block_data (grid_grace, f_grace_err)  
            CALL ncio_read_block_time (file_grace, 'lwe_thickness', grid_grace, itime, f_grace_lwe)
            CALL ncio_read_block_time (file_grace, 'uncertainty'  , grid_grace, itime, f_grace_err)
         ENDIF

         CALL mg2p_grace%map_aweighted (f_grace_lwe, lwe_obs_this)
         CALL mg2p_grace%map_aweighted (f_grace_err, err_obs_this)

         IF (p_is_worker) THEN
            
            lwe_obs_this = lwe_obs_this * 10.0 ! from cm to mm
            err_obs_this = err_obs_this * 10.0 ! from cm to mm

            wat_this_m  = wat_this_m  / nac_grace

            IF (has_prev_grace_obs) THEN

               rsub_prev_m1 = rsub_prev_m1 / nac_grace

               DO ielm = 1, numelm
                  istt = elm_patch%substt(ielm)
                  iend = elm_patch%subend(ielm)

                  w1 = sum(wat_this_m  (istt:iend) * elm_patch%subfrc(istt:iend)) 
                  w0 = sum(wat_prev_m  (istt:iend) * elm_patch%subfrc(istt:iend)) 
                  r1 = sum(rsub_prev_m1(istt:iend) * elm_patch%subfrc(istt:iend)) 
                  r0 = sum(rsub_prev_m0(istt:iend) * elm_patch%subfrc(istt:iend)) 

                  var_o = err_obs_this(ielm)**2 + err_obs_prev(ielm)**2
                     
                  dw_f = w1 - w0
                  dw_o = lwe_obs_this(ielm) - lwe_obs_prev(ielm)
                  var_m = (dw_f-dw_o)**2 - var_o

                  fslp_rsub = -1.

                  IF (var_m > 0) THEN
                 
                     dw_a = (var_o * dw_f + var_m * dw_o) / (var_m+var_o)

                     ! var_a = var_o * var_m /(var_o+var_m)
                     var_c = 1.0**2

                     rr = r1 - r0
                     IF (rr > 0) THEN

                        fprev = fslp_patch(istt)
                        fslp_rsub = fprev + (dw_f - dw_a)/rr &
                           * (rr**2 * var_c) / (rr**2 * var_c + var_m) 
                        fslp_rsub = min(max(fslp_rsub, fprev*0.5), fprev*2.0)

                        fslp_patch(istt:iend) = fslp_rsub
                     ENDIF

                  ENDIF
                        
                  IF (fslp_rsub > 0) THEN
                     ! write(*,'(A,2I3,7ES11.2)') 'Check DA: ', p_iam_worker, ielm, &
                     !    dw_o, sqrt(var_o), dw_f, sign(sqrt(abs(var_m)),var_m), dw_a, r1-r0, fslp_rsub 
                  ELSE
                     ! write(*,'(A,2I3,6ES11.2)') 'Check DA: ', p_iam_worker, ielm, &
                     !    dw_o, sqrt(var_o), dw_f, sign(sqrt(abs(var_m)),var_m), r1-r0, fslp_rsub 
                  ENDIF

               ENDDO
            ENDIF

            lwe_obs_prev = lwe_obs_this
            err_obs_prev = err_obs_this

            wat_prev_m = wat_this_m
            wat_this_m = 0.

            rsub_acc_prev_m = rsub_acc_this_m
            rsub_acc_this_m = 0.

            rsub_prev_m0 = rsub_this_m / nac_grace
            rsub_prev_m1 = 0.
            rsub_this_m  = 0.

            nac_grace = 0

         ENDIF

         has_prev_grace_obs = .true.

      ENDIF

   END SUBROUTINE do_DA_GRACE

   ! ---------
   SUBROUTINE final_DA_GRACE ()

      IMPLICIT NONE

      IF (allocated(lwe_obs_this))    deallocate(lwe_obs_this)
      IF (allocated(err_obs_this))    deallocate(err_obs_this)
      IF (allocated(lwe_obs_prev))    deallocate(lwe_obs_prev)
      IF (allocated(err_obs_prev))    deallocate(err_obs_prev)
      IF (allocated(wat_prev_m     )) deallocate(wat_prev_m     )
      IF (allocated(wat_this_m     )) deallocate(wat_this_m     )
      IF (allocated(rsub_acc_prev_m)) deallocate(rsub_acc_prev_m)
      IF (allocated(rsub_acc_this_m)) deallocate(rsub_acc_this_m)
      IF (allocated(rsub_prev_m0   )) deallocate(rsub_prev_m0   )
      IF (allocated(rsub_prev_m1   )) deallocate(rsub_prev_m1   )
      IF (allocated(rsub_this_m    )) deallocate(rsub_this_m    )

      IF (allocated(fslp_patch)) deallocate(fslp_patch)

      IF (allocated(longrace)) deallocate(longrace)
      IF (allocated(latgrace)) deallocate(latgrace)

   END SUBROUTINE final_DA_GRACE 

   ! ---------
   SUBROUTINE retrieve_yymm_from_days (days, yy, mm)

      IMPLICIT NONE
      REAL(r8), intent(in)  :: days
      INTEGER,  intent(out) :: yy, mm

      ! Local Variables
      REAL(r8) :: resday
      INTEGER  :: mdays(12)

      yy = 2002
      mm = 1
      mdays = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      
      resday = days
      DO WHILE (resday > mdays(mm))
           
         resday = resday - mdays(mm)

         mm = mm + 1
         IF (mm > 12) THEN
            yy = yy + 1
            mm = 1
            IF( (mod(yy,4)==0 .AND. mod(yy,100)/=0) .OR. mod(yy,400)==0 ) THEN
               mdays(2) = 29
            ELSE
               mdays(2) = 28
            ENDIF
         ENDIF

      ENDDO

   END SUBROUTINE retrieve_yymm_from_days

END MODULE MOD_DA_GRACE
#endif
