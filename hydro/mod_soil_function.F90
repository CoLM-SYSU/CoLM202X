#include <define.h>

module mod_soil_function

   use precision

   implicit none

   ! Soil function type 1:
   ! Campbell model
   ! CAMPBELL, G. S. (1974), Soil Science, 117(6), 311-314.
   
   ! Soil function type 2:
   ! Modified van Genuchten & Mualem model by introducing an air-entry value 
   ! Ippisch et al. (2006), Advances in Water Resources, 29(12), 1780-1789. 

   PUBLIC :: get_derived_parameters_vGM 

   public :: soil_psi_from_vliq 
   public :: soil_hk_from_psi 
   public :: soil_vliq_from_psi 

contains

   !-------------------------------------
   subroutine get_derived_parameters_vGM ( &
         psi_s, alpha_vgm, n_vgm, sc_vgm, fc_vgm)

      real(r8), intent(in) :: psi_s
      real(r8), intent(in) :: alpha_vgm
      real(r8), intent(in) :: n_vgm
      
      real(r8), intent(out) :: sc_vgm
      real(r8), intent(out) :: fc_vgm

      ! Local variables
      real(r8) :: m_vgm

      m_vgm = 1.0_r8 - 1.0_r8 / n_vgm
      sc_vgm = (1.0_r8 + (- alpha_vgm * psi_s)**n_vgm) ** (-m_vgm)
      fc_vgm = 1.0_r8 - (1.0_r8 - sc_vgm ** (1.0_r8/m_vgm)) ** m_vgm

   end subroutine get_derived_parameters_vGM

   !------------------------------------------------------------------
   real(r8) function soil_hk_from_psi (psi, &
         psi_s, hksat, nprm, prms)

      implicit none

      real(r8), intent(in) :: psi

      real(r8), intent(in) :: psi_s
      real(r8), intent(in) :: hksat
      
      integer,  intent(in) :: nprm
      real(r8), intent(in) :: prms(nprm)
        
      ! Local variables
      real(r8) :: m_vgm, esat

      if (psi >= psi_s) then
         soil_hk_from_psi = hksat
         return
      end if

#ifdef Campbell_SOIL_MODEL
      ! bsw => prms(1)
      soil_hk_from_psi = hksat * (psi / psi_s)**(- 3.0_r8 / prms(1) - 2.0_r8)
#endif
      
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      ! alpha_vgm => prms(1), n_vgm => prms(2), L_vgm => prms(3), sc_vgm => prms(4), fc_vgm => prms(5)
      m_vgm = 1.0_r8 - 1.0_r8 / prms(2)
      esat = (1.0_r8 + (- prms(1) * psi)**(prms(2)))**(-m_vgm) / prms(4)
      soil_hk_from_psi = hksat * esat**prms(3) &
         * ((1.0_r8 - (1.0_r8 - (esat*prms(4))**(1.0_r8/m_vgm))**m_vgm) / prms(5))**2.0_r8 
#endif

   end function soil_hk_from_psi


   !-----------------------------------------------------------------
   real(r8) function soil_psi_from_vliq (vliq, &
         porsl, vl_r, psi_s, nprm, prms)
    
      implicit none
      
      real(r8), intent(in) :: vliq

      real(r8), intent(in) :: porsl
      real(r8), intent(in) :: vl_r
      real(r8), intent(in) :: psi_s
      
      integer,  intent(in) :: nprm
      real(r8), intent(in) :: prms(nprm)

      ! Local variables
      real(r8) :: esat, m_vgm

      if (vliq >= porsl) then
         soil_psi_from_vliq = psi_s
         return
      end if

#ifdef Campbell_SOIL_MODEL
      ! bsw => prms(1)
      soil_psi_from_vliq = psi_s * (vliq / porsl)**(-prms(1))
#endif
      
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      ! alpha_vgm => prms(1), n_vgm => prms(2), L_vgm => prms(3), sc_vgm => prms(4), fc_vgm => prms(5)
      m_vgm = 1.0_r8 - 1.0_r8 / prms(2)
      esat = (vliq - vl_r) / (porsl - vl_r)
      soil_psi_from_vliq = - ((esat*prms(4))**(- 1.0_r8/m_vgm) - 1.0_r8)**(1.0_r8/prms(2)) &
         / prms(1)
#endif

   end function soil_psi_from_vliq

   !------------------------------------------------------------------
   real(r8) function soil_vliq_from_psi (psi, &
         porsl, vl_r, psi_s, nprm, prms)

      implicit none
    
      real(r8), intent(in) :: psi

      real(r8), intent(in) :: porsl
      real(r8), intent(in) :: vl_r
      real(r8), intent(in) :: psi_s
      
      integer,  intent(in) :: nprm
      real(r8), intent(in) :: prms(nprm)
        
      ! Local variables
      real(r8) :: esat, m_vgm

      if (psi >= psi_s) then
         soil_vliq_from_psi = porsl
         return
      end if

#ifdef Campbell_SOIL_MODEL
      ! bsw => prms(1)
      soil_vliq_from_psi = porsl * (psi / psi_s)**(-1.0/prms(1))
#endif
      
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      ! alpha_vgm => prms(1), n_vgm => prms(2), L_vgm => prms(3), sc_vgm => prms(4), fc_vgm => prms(5)
      m_vgm = 1.0_r8 - 1.0_r8 / prms(2)
      esat = (1.0_r8 + (psi * (-prms(1)))**(prms(2))) ** (-m_vgm) / prms(4)
      soil_vliq_from_psi = (porsl - vl_r) * esat + vl_r  
#endif

   end function soil_vliq_from_psi


end module mod_soil_function
