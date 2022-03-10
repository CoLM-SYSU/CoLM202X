#include <define.h>

module MOD_1D_Acc_cama_Fluxes

   use precision

   real(r8) :: nacc              ! number of accumulation
   real(r8), allocatable :: a_rnof_cama   (:)
 

   public :: allocate_acc_cama_fluxes
   public :: deallocate_acc_cama_fluxes
   public :: flush_acc_cama_fluxes
   public :: accumulate_cama_fluxes

contains

   subroutine allocate_acc_cama_fluxes 

      use spmd_task
      USE GlobalVars
      use mod_landpatch, only : numpatch
      implicit none

      if (p_is_worker) then
         if (numpatch > 0) then
            allocate (a_rnof_cama      (numpatch))
         end if
      end if

   end subroutine allocate_acc_cama_fluxes

   subroutine deallocate_acc_cama_fluxes()

      use spmd_task
      use mod_landpatch, only : numpatch
      implicit none

      if (p_is_worker) then
         if (numpatch > 0) then
            deallocate (a_rnof_cama)
         end if
      end if

   end subroutine deallocate_acc_cama_fluxes

   !-----------------------
   SUBROUTINE FLUSH_acc_cama_fluxes ()

      use spmd_task
      use mod_landpatch, only : numpatch
      use GlobalVars,    only : spval 
      implicit none

      if (p_is_worker) then

         nacc = 0

         if (numpatch > 0) then
            ! flush the Fluxes for accumulation
            a_rnof_cama    (:) = spval
         end if
      end if

   END SUBROUTINE FLUSH_acc_cama_fluxes

   SUBROUTINE accumulate_cama_fluxes 
      ! ----------------------------------------------------------------------
      ! perfrom the grid average mapping: average a subgrid input 1d vector 
      ! of length numpatch to a output 2d array of length [ghist%xcnt,ghist%ycnt]
      !
      ! Created by Yongjiu Dai, 03/2014
      !---------------------------------------------------------------------

      use precision
      use spmd_task
      use mod_landpatch,     only : numpatch
      use MOD_TimeInvariants
      use MOD_TimeVariables
      use MOD_1D_Forcing
      use MOD_1D_Fluxes
      use FRICTION_VELOCITY
      use mod_colm_debug
      use GlobalVars

      IMPLICIT NONE

      if (p_is_worker) then
         if (numpatch > 0) then
            nacc = nacc + 1
            call acc1d_cama (rnof, a_rnof_cama   )  
         end if
      end if

   END SUBROUTINE accumulate_cama_fluxes


   !------
   SUBROUTINE acc1d_cama (var, s)

      use precision
      use GlobalVars, only: spval

      IMPLICIT NONE

      real(r8), intent(in)    :: var(:)
      real(r8), intent(inout) :: s  (:)
      ! Local variables
      integer :: i

      do i = lbound(var,1), ubound(var,1)
         if (var(i) /= spval) then
            if (s(i) /= spval) then
               s(i) = s(i) + var(i)
            else
               s(i) = var(i)
            end if
         end if
      end do
      
   END SUBROUTINE acc1d_cama

   
end module MOD_1D_Acc_cama_Fluxes
! ----- EOP ---------
