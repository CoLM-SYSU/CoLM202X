#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_WriteParameters

CONTAINS

   SUBROUTINE write_catch_parameters 

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Catch_IO
   USE MOD_ElmVector
   USE MOD_Vars_TimeInvariants,  only : wf_sand, wf_clay, wf_om, wf_gravels, patchclass
   USE MOD_Catch_SubsurfaceFlow, only : hillslope_element, lake_id_elm
   USE MOD_HRUVector,   only : totalnumhru, hru_data_address, eindx_hru, htype_hru
   USE MOD_Mesh,        only : numelm
   USE MOD_LandHRU,     only : numhru
   USE MOD_LandPatch,   only : hru_patch
   USE MOD_Vars_Global, only : spval
   IMPLICIT NONE

   character(len=256) :: file_parameters
   character(len=1)   :: slev
   real(r8), allocatable :: hrupara(:)
   integer :: ilev, i, ps, pe, ielm, ihru, j

      file_parameters = trim(DEF_dir_output) // '/' // trim(DEF_CASE_NAME) // '/catch_parameters.nc'

      IF (p_is_master) THEN
         CALL ncio_create_file (trim(file_parameters))
         CALL ncio_define_dimension(file_parameters, 'hydrounit', totalnumhru)

         CALL ncio_write_serial (file_parameters, 'bsn_hru', eindx_hru, 'hydrounit')
         CALL ncio_put_attr (file_parameters, 'bsn_hru', &
            'long_name', 'basin index of hydrological units')

         CALL ncio_write_serial (file_parameters, 'hru_type' , htype_hru, 'hydrounit')
         CALL ncio_put_attr (file_parameters, 'hru_type' , &
            'long_name', 'index of hydrological units inside basin')
      ENDIF

      IF (p_is_worker) THEN
         IF (numhru > 0) allocate (hrupara (numhru))
      ENDIF

      DO ilev = 1, 5
         
         write(slev,'(i1.1)') ilev 
        
         ! sand
         IF (p_is_worker) THEN
            DO i = 1, numhru
               ps = hru_patch%substt(i)
               pe = hru_patch%subend(i)
               hrupara(i) = sum(wf_sand(ilev,ps:pe) * hru_patch%subfrc(ps:pe))
            ENDDO
         ENDIF

         CALL vector_write_basin (&
            file_parameters, hrupara, numhru, totalnumhru, 'wf_sand_l'//slev, 'hydrounit', hru_data_address)
         
         ! clay
         IF (p_is_worker) THEN
            DO i = 1, numhru
               ps = hru_patch%substt(i)
               pe = hru_patch%subend(i)
               hrupara(i) = sum(wf_clay(ilev,ps:pe) * hru_patch%subfrc(ps:pe))
            ENDDO
         ENDIF

         CALL vector_write_basin (&
            file_parameters, hrupara, numhru, totalnumhru, 'wf_clay_l'//slev, 'hydrounit', hru_data_address)
         
         ! organic matter 
         IF (p_is_worker) THEN
            DO i = 1, numhru
               ps = hru_patch%substt(i)
               pe = hru_patch%subend(i)
               hrupara(i) = sum(wf_om(ilev,ps:pe) * hru_patch%subfrc(ps:pe))
            ENDDO
         ENDIF

         CALL vector_write_basin (&
            file_parameters, hrupara, numhru, totalnumhru, 'wf_om_l'//slev, 'hydrounit', hru_data_address)
         
         ! silt 
         IF (p_is_worker) THEN
            DO i = 1, numhru
               ps = hru_patch%substt(i)
               pe = hru_patch%subend(i)
               hrupara(i) = sum((1-wf_sand(ilev,ps:pe)-wf_gravels(ilev,ps:pe)-wf_om(ilev,ps:pe)-wf_clay(ilev,ps:pe)) * hru_patch%subfrc(ps:pe))
            ENDDO
         ENDIF

         CALL vector_write_basin (&
            file_parameters, hrupara, numhru, totalnumhru, 'wf_silt_l'//slev, 'hydrounit', hru_data_address)
      ENDDO



      IF (p_is_worker) THEN
         DO i = 1, numhru
            ps = hru_patch%substt(i)
            pe = hru_patch%subend(i)
            hrupara(i) = patchclass(ps)
         ENDDO
      ENDIF

      CALL vector_write_basin (&
         file_parameters, hrupara, numhru, totalnumhru, 'lulc_igbp', 'hydrounit', hru_data_address)

      IF (p_is_worker) THEN
         IF (numhru > 0) hrupara(:) = spval
         DO ielm = 1, numelm
            IF (lake_id_elm(ielm) <= 0) THEN
               DO i = 1, hillslope_element(ielm)%nhru
                  IF (hillslope_element(ielm)%indx(i) /= 0) THEN
                     ihru = hillslope_element(ielm)%ihru(i)
                     hrupara(ihru) = hillslope_element(ielm)%plen(i) * 2.
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
         
      CALL vector_write_basin (&
         file_parameters, hrupara, numhru, totalnumhru, 'slope_length', 'hydrounit', hru_data_address)

      IF (p_is_master) THEN
         CALL ncio_put_attr (file_parameters, 'slope_length', 'units', 'm')
         CALL ncio_put_attr (file_parameters, 'slope_length', 'missing_value', spval)
      ENDIF

      IF (p_is_worker) THEN
         IF (numhru > 0) hrupara(:) = spval
         DO ielm = 1, numelm
            IF (lake_id_elm(ielm) <= 0) THEN
               DO i = 1, hillslope_element(ielm)%nhru
                  IF (hillslope_element(ielm)%indx(i) /= 0) THEN
                     ihru = hillslope_element(ielm)%ihru(i)
                     hrupara(ihru) = hillslope_element(ielm)%elva(i)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
         
      CALL vector_write_basin (&
         file_parameters, hrupara, numhru, totalnumhru, 'elevation', 'hydrounit', hru_data_address)
      
      IF (p_is_master) THEN
         CALL ncio_put_attr (file_parameters, 'elevation', 'units', 'm')
         CALL ncio_put_attr (file_parameters, 'elevation', 'missing_value', spval)
      ENDIF

      IF (p_is_worker) THEN
         IF (numhru > 0) hrupara(:) = spval
         DO ielm = 1, numelm
            IF (lake_id_elm(ielm) <= 0) THEN
               DO i = 1, hillslope_element(ielm)%nhru
                  IF (hillslope_element(ielm)%indx(i) /= 0) THEN
                     ihru = hillslope_element(ielm)%ihru(i)
                     j    = hillslope_element(ielm)%inext(i)
                     IF (j > 0) THEN
                        hrupara(ihru) = (hillslope_element(ielm)%hand(i) - hillslope_element(ielm)%hand(j)) &
                           / (hillslope_element(ielm)%plen(i) + hillslope_element(ielm)%plen(j))
                     ELSE
                        hrupara(ihru) = hillslope_element(ielm)%hand(i) / hillslope_element(ielm)%plen(i)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
         
      CALL vector_write_basin (&
         file_parameters, hrupara, numhru, totalnumhru, 'slope_ratio', 'hydrounit', hru_data_address)

      IF (p_is_master) THEN
         CALL ncio_put_attr (file_parameters, 'slope_ratio', 'units', '-')
         CALL ncio_put_attr (file_parameters, 'slope_ratio', 'missing_value', spval)
      ENDIF

      IF (allocated(hrupara)) deallocate(hrupara)

   END SUBROUTINE write_catch_parameters

END MODULE MOD_Catch_WriteParameters
#endif
