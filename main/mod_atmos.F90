#include <define.h>

module mod_atms

   use precision
   use mod_grid
   use mod_mapping_pset2grid
   implicit none

   type (grid_type) :: gatms
   type (mapping_pset2grid_type) :: mp2g_atms

   ! Required by atmospheric models initialization (such as GRAPES, RSM, ...)
   type(block_data_real8_2d) :: tg_xy    ! 
   type(block_data_real8_2d) :: albvb_xy ! 
   type(block_data_real8_2d) :: albvd_xy ! 
   type(block_data_real8_2d) :: albnb_xy ! 
   type(block_data_real8_2d) :: albnd_xy ! 
   type(block_data_real8_2d) :: trad_xy  ! 
   type(block_data_real8_2d) :: rib_xy   ! 
   type(block_data_real8_2d) :: fm_xy    ! 
   type(block_data_real8_2d) :: fh_xy    ! 
   type(block_data_real8_2d) :: fq_xy    ! 

contains

   subroutine atms_init 

      implicit none

      call gatm%define ()
      call mp2g_atms%build (pixl, patch, gatm)

      call allocate_block_data (gatms, tg_xy   )
      call allocate_block_data (gatms, albvb_xy)
      call allocate_block_data (gatms, albvd_xy)
      call allocate_block_data (gatms, albnb_xy)
      call allocate_block_data (gatms, albnd_xy)
      call allocate_block_data (gatms, trad_xy )
      call allocate_block_data (gatms, rib_xy  )
      call allocate_block_data (gatms, fm_xy   )
      call allocate_block_data (gatms, fh_xy   )
      call allocate_block_data (gatms, fq_xy   )

   end subroutine atms_init 

   subroutine patch_to_atm_model

      implicit none
      ! ----------------------------------------------------------------------
      ! average subgrid albedos, srf temperature, etc. for atmospheric model
      ! ----------------------------------------------------------------------
      tg_xy   (:,:) = 0.0
      albvb_xy(:,:) = 0.0
      albvd_xy(:,:) = 0.0
      albnb_xy(:,:) = 0.0
      albnd_xy(:,:) = 0.0
      trad_xy (:,:) = 0.0

      call mp2a%map (t_grnd,   tg_xy   )
      call mp2a%map (alb(1,1), albvb_xy)
      call mp2a%map (alb(1,2), albvd_xy)
      call mp2a%map (alb(2,1), albnb_xy)
      call mp2a%map (alb(2,2), albnd_xy)
      call mp2a%map (t_grnd,   trad_xy )

      rib_xy (:,:) = -0.1       
      fm_xy  (:,:) = alog(30.)  
      fh_xy  (:,:) = alog(30.)  
      fq_xy  (:,:) = alog(30.)  

   end subroutine patch_to_atm_model


   subroutine atms_final
      
      implicit none

      deallocate ( tg_xy    )
      deallocate ( albvb_xy )
      deallocate ( albvd_xy )
      deallocate ( albnb_xy )
      deallocate ( albnd_xy )
      deallocate ( trad_xy  )
      deallocate ( rib_xy   )
      deallocate ( fm_xy    )
      deallocate ( fh_xy    )
      deallocate ( fq_xy    )

   end subroutine atms_final

end module mod_atms
