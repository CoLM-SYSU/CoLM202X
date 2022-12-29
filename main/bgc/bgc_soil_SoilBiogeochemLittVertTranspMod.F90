#include <define.h>
#ifdef BGC
module bgc_soil_SoilBiogeochemLittVertTranspMod

use precision
use MOD_BGCTimeInvars, only: &
    is_cwd, som_adv_flux, som_diffus, cryoturb_diffusion_k, max_altdepth_cryoturbation, max_depth_cryoturb
use MOD_BGCTimeVars, only: &
    altmax, altmax_lastyear, som_adv_coef, som_diffus_coef, &
    decomp_cpools_vr, decomp_npools_vr, &
    diagVX_c_vr_acc, upperVX_c_vr_acc, lowerVX_c_vr_acc, &
    diagVX_n_vr_acc, upperVX_n_vr_acc, lowerVX_n_vr_acc
use MOD_1D_BGCFluxes, only: &
    decomp_cpools_sourcesink, decomp_npools_sourcesink, &
    decomp_cpools_transport_tendency, decomp_npools_transport_tendency

implicit none

public SoilBiogeochemLittVertTransp

contains

subroutine SoilBiogeochemLittVertTransp(i,deltim,nl_soil,nl_soil_full,ndecomp_pools,nbedrock,z_soi,zi_soi,dz_soi)

integer ,intent(in) :: i
real(r8),intent(in) :: deltim                                                  ! land model time step (sec)
integer ,intent(in) :: nl_soil
integer ,intent(in) :: nl_soil_full
integer ,intent(in) :: ndecomp_pools
integer ,intent(in) :: nbedrock
real(r8),intent(in) :: z_soi (1:nl_soil_full)
real(r8),intent(in) :: zi_soi(0:nl_soil_full)
real(r8),intent(in) :: dz_soi(1:nl_soil_full)

! !LOCAL VARIABLES:
real(r8) :: diffus (1:nl_soil+1)                    ! diffusivity (m2/s)  (includes spinup correction, if any)
real(r8) :: adv_flux(1:nl_soil+1)                   ! advective flux (m/s)  (includes spinup correction, if any)
real(r8) :: aaa                                                                ! "A" function in Patankar
real(r8) :: pe                                                                 ! Pe for "A" function in Patankar
real(r8) :: w_m1, w_p1                                                         ! Weights for calculating harmonic mean of diffusivity
real(r8) :: d_m1, d_p1                                                         ! Harmonic mean of diffusivity
real(r8) :: a_tri(0:nl_soil+1)                      ! "a" vector for tridiagonal matrix
real(r8) :: b_tri(0:nl_soil+1)                      ! "b" vector for tridiagonal matrix
real(r8) :: c_tri(0:nl_soil+1)                      ! "c" vector for tridiagonal matrix
real(r8) :: r_tri_c(0:nl_soil+1)                    ! "r" vector for tridiagonal solution for soil C
real(r8) :: r_tri_n(0:nl_soil+1)                    ! "r" vector for tridiagonal solution for soil N
real(r8) :: d_p1_zp1(1:nl_soil+1)                   ! diffusivity/delta_z for next j  (set to zero for no diffusion)
real(r8) :: d_m1_zm1(1:nl_soil+1)                   ! diffusivity/delta_z for previous j (set to zero for no diffusion)
real(r8) :: f_p1(1:nl_soil+1)                       ! water flux for next j
real(r8) :: f_m1(1:nl_soil+1)                       ! water flux for previous j
real(r8) :: pe_p1(1:nl_soil+1)                      ! Peclet # for next j
real(r8) :: pe_m1(1:nl_soil+1)                      ! Peclet # for previous j
real(r8) :: dz_node(1:nl_soil+1)                                            ! difference between nodes
real(r8) :: conc_trcr_c(0:nl_soil+1)                  ! dummy term
real(r8) :: conc_trcr_n(0:nl_soil+1)                  ! dummy term
real(r8) :: a_p_0
integer  :: s,j,l                                                  ! indices
integer  :: jtop                                      ! top level at each column
real(r8) :: spinup_term                                                        ! spinup accelerated decomposition factor, used to accelerate transport as well
real(r8) :: epsilon                                                            ! small number


!real(r8) :: spinup_term = 1.0_r8
!real(r8) :: som_adv_flux =  0._r8
!real(r8) :: max_depth_cryoturb = 3._r8   ! (m) this is the maximum depth of cryoturbation
!real(r8) :: som_diffus                   ! [m^2/sec] = 1 cm^2 / yr
!real(r8) :: cryoturb_diffusion_k         ! [m^2/sec] = 5 cm^2 / yr = 1m^2 / 200 yr
!real(r8) :: max_altdepth_cryoturbation   ! (m) maximum active layer thickness for cryoturbation to occur

aaa (pe) = max (0._r8, (1._r8 - 0.1_r8 * abs(pe))**5)  ! A function from Patankar, Table 5.2, pg 95

   epsilon = 1.e-30
   spinup_term = 1._r8

   if  (( max(altmax(i), altmax_lastyear(i)) <= max_altdepth_cryoturbation ) .and. &
        ( max(altmax(i), altmax_lastyear(i)) > 0._r8) ) then
      ! use mixing profile modified slightly from Koven et al. (2009): constant through active layer, linear decrease from base of active layer to zero at a fixed depth
      do j = 1,nl_soil+1
         if ( j <= nbedrock+1 ) then
            if ( zi_soi(j) < max(altmax(i), altmax_lastyear(i)) ) then
               som_diffus_coef(j,i) = cryoturb_diffusion_k
               som_adv_coef(j,i) = 0._r8
            else
               som_diffus_coef(j,i) = max(cryoturb_diffusion_k * &
                 ( 1._r8 - ( zi_soi(j) - max(altmax(i), altmax_lastyear(i)) ) / &
                 ( min(max_depth_cryoturb, zi_soi(nbedrock+1)) - max(altmax(i), altmax_lastyear(i)) ) ), 0._r8)  ! go linearly to zero between ALT and max_depth_cryoturb
               som_adv_coef(j,i) = 0._r8
            endif
         else
            som_adv_coef(j,i) = 0._r8
            som_diffus_coef(j,i) = 0._r8
         endif
      end do
   elseif (  max(altmax(i), altmax_lastyear(i)) > 0._r8 ) then
      ! constant advection, constant diffusion
      do j = 1,nl_soil+1
         if ( j <= nbedrock+1 ) then
            som_adv_coef(j,i) = som_adv_flux
            som_diffus_coef(j,i) = som_diffus
         else
            som_adv_coef(j,i) = 0._r8
            som_diffus_coef(j,i) = 0._r8
         endif
      end do
   else
      ! completely frozen soils--no mixing
      do j = 1,nl_soil+1
         som_adv_coef(j,i) = 0._r8
         som_diffus_coef(j,i) = 0._r8
      end do
   endif

      ! Set the distance between the node and the one ABOVE it
   dz_node(1) = z_soi(1)
   do j = 2, nl_soil+1
      dz_node(j)= z_soi(j) - z_soi(j-1)
   enddo

   do s = 1, ndecomp_pools
      if ( .not. is_cwd(s) ) then
!         if(.not. use_soil_matrixcn .or. s .eq. 1)then
         do j = 1,nl_soil+1
               !
!               if ( spinup_state >= 1 ) then
!                  ! increase transport (both advection and diffusion) by the same factor as accelerated decomposition for a given pool
!                  spinup_term = spinup_factor(s)
!                        else
!                           spinup_term = 1._r8
!                        endif

!                        if (abs(spinup_term - 1._r8) > .000001_r8 ) then
!                           spinup_term = spinup_term * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
!                        endif

            if ( abs(som_adv_coef(j,i)) * spinup_term < epsilon ) then
               adv_flux(j) = epsilon
            else
                  adv_flux(j) = som_adv_coef(j,i) * spinup_term
            endif
            !
            if ( abs(som_diffus_coef(j,i)) * spinup_term < epsilon ) then
               diffus(j) = epsilon
            else
               diffus(j) = som_diffus_coef(j,i) * spinup_term
            endif
            !
         end do

                  ! Set Pe (Peclet #) and D/dz throughout column

         conc_trcr_c(0) = 0._r8
         conc_trcr_n(0) = 0._r8
         conc_trcr_c(nbedrock+1:nl_soil+1) = 0._r8
         conc_trcr_n(nbedrock+1:nl_soil+1) = 0._r8

         do j = 1,nl_soil+1
            conc_trcr_c(j) = decomp_cpools_vr(j,s,i)
            conc_trcr_n(j) = decomp_npools_vr(j,s,i)

            ! dz_tracer below is the difference between gridcell edges  (dz_soi)
            ! dz_node_tracer is difference between cell centers

            ! Calculate the D and F terms in the Patankar algorithm
            if (j == 1) then
               d_m1_zm1(j) = 0._r8
               w_p1 = (z_soi(j+1) - zi_soi(j)) / dz_node(j+1)
               if ( diffus(j+1) > 0._r8 .and. diffus(j) > 0._r8) then
                  d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus(j) + w_p1 / diffus(j+1)) ! Harmonic mean of diffus
               else
                  d_p1 = 0._r8
               endif
               d_p1_zp1(j) = d_p1 / dz_node(j+1)
               f_m1(j) = adv_flux(j)  ! Include infiltration here
               f_p1(j) = adv_flux(j+1)
               pe_m1(j) = 0._r8
               pe_p1(j) = f_p1(j) / d_p1_zp1(j) ! Peclet #
            elseif (j >= nbedrock+1) then
               ! At the bottom, assume no gradient in d_z (i.e., they're the same)
               w_m1 = (zi_soi(j-1) - z_soi(j-1)) / dz_node(j)
               if ( diffus(j) > 0._r8 .and. diffus(j-1) > 0._r8) then
                  d_m1 = 1._r8 / ((1._r8 - w_m1) / diffus(j) + w_m1 / diffus(j-1)) ! Harmonic mean of diffus
               else
                  d_m1 = 0._r8
               endif
               d_m1_zm1(j) = d_m1 / dz_node(j)
               d_p1_zp1(j) = d_m1_zm1(j) ! Set to be the same
               f_m1(j) = adv_flux(j)
               !f_p1(j) = adv_flux(j+1)
               f_p1(j) = 0._r8
               pe_m1(j) = f_m1(j) / d_m1_zm1(j) ! Peclet #
               pe_p1(j) = f_p1(j) / d_p1_zp1(j) ! Peclet #
            else
               ! Use distance from j-1 node to interface with j divided by distance between nodes
               w_m1 = (zi_soi(j-1) - z_soi(j-1)) / dz_node(j)
               if ( diffus(j-1) > 0._r8 .and. diffus(j) > 0._r8) then
                  d_m1 = 1._r8 / ((1._r8 - w_m1) / diffus(j) + w_m1 / diffus(j-1)) ! Harmonic mean of diffus
               else
                  d_m1 = 0._r8
               endif
               w_p1 = (z_soi(j+1) - zi_soi(j)) / dz_node(j+1)
               if ( diffus(j+1) > 0._r8 .and. diffus(j) > 0._r8) then
                  d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus(j) + w_p1 / diffus(j+1)) ! Harmonic mean of diffus
               else
                  d_p1 = (1._r8 - w_m1) * diffus(j) + w_p1 * diffus(j+1) ! Arithmetic mean of diffus
               endif
               d_m1_zm1(j) = d_m1 / dz_node(j)
               d_p1_zp1(j) = d_p1 / dz_node(j+1)
               f_m1(j) = adv_flux(j)
               f_p1(j) = adv_flux(j+1)
               pe_m1(j) = f_m1(j) / d_m1_zm1(j) ! Peclet #
               pe_p1(j) = f_p1(j) / d_p1_zp1(j) ! Peclet #
            end if
         enddo ! j; nl_soil

         ! Calculate the tridiagonal coefficients
         do j = 0,nl_soil +1

            if (j > 0 .and. j < nl_soil+1) then
               a_p_0 =  dz_soi(j) / deltim
            endif

            if (j == 0) then ! top layer (atmosphere)
               a_tri(j)   = 0._r8
               b_tri(j)   = 1._r8
               c_tri(j)   = -1._r8
               r_tri_c(j) = 0._r8
               r_tri_n(j) = 0._r8
            elseif (j == 1) then
               a_tri(j) = -(d_m1_zm1(j) * aaa(pe_m1(j)) + max( f_m1(j), 0._r8)) ! Eqn 5.47 Patankar
               c_tri(j) = -(d_p1_zp1(j) * aaa(pe_p1(j)) + max(-f_p1(j), 0._r8))
               b_tri(j) = - a_tri(j) - c_tri(j) + a_p_0
               r_tri_c(j) = decomp_cpools_sourcesink(j,s,i) * dz_soi(j) /deltim + (a_p_0 - adv_flux(j)) * conc_trcr_c(j)
               r_tri_n(j) = decomp_npools_sourcesink(j,s,i) * dz_soi(j) /deltim + (a_p_0 - adv_flux(j)) * conc_trcr_n(j)
#ifdef SASU
               upperVX_c_vr_acc(j,s,i) = upperVX_c_vr_acc(j,s,i) -  c_tri(j)          / dz_soi(j) * deltim * conc_trcr_c(j+1)! upwards transfer
               diagVX_c_vr_acc (j,s,i) = diagVX_c_vr_acc (j,s,i) + (b_tri(j) - a_p_0) / dz_soi(j) * deltim * conc_trcr_c(j)! exit flux
               upperVX_n_vr_acc(j,s,i) = upperVX_n_vr_acc(j,s,i) -  c_tri(j)          / dz_soi(j) * deltim * conc_trcr_n(j+1)! upwards transfer
               diagVX_n_vr_acc (j,s,i) = diagVX_n_vr_acc (j,s,i) + (b_tri(j) - a_p_0) / dz_soi(j) * deltim * conc_trcr_n(j)! exit flux
#endif
!               if(s .eq. 1 .and. i_type .eq. 1 .and. use_soil_matrixcn)then !vertical matrix are the same for all pools
!                  do i = 1,ndecomp_pools-1 !excluding cwd
!                     tri_ma_vr(c,1+(i-1)*(nl_soil*3-2)) = (b_tri(j,i) - a_p_0) / dz_soi(j) * (-deltim)
!                     tri_ma_vr(c,3+(i-1)*(nl_soil*3-2)) = c_tri(j,i) / dz_soi(j) * (-deltim)
!                  end do
!               end if
            elseif (j < nl_soil+1) then
               
               a_tri(j) = -(d_m1_zm1(j) * aaa(pe_m1(j)) + max( f_m1(j), 0._r8)) ! Eqn 5.47 Patankar
               c_tri(j) = -(d_p1_zp1(j) * aaa(pe_p1(j)) + max(-f_p1(j), 0._r8))
               b_tri(j) = - a_tri(j) - c_tri(j) + a_p_0
               r_tri_c(j) = decomp_cpools_sourcesink(j,s,i) * dz_soi(j) /deltim + a_p_0 * conc_trcr_c(j)
               r_tri_n(j) = decomp_npools_sourcesink(j,s,i) * dz_soi(j) /deltim + a_p_0 * conc_trcr_n(j)
               
#ifdef SASU
               if(j .le. nbedrock)then
                  lowerVX_c_vr_acc(j,s,i) = lowerVX_c_vr_acc(j,s,i) - a_tri(j) / dz_soi(j) * deltim * conc_trcr_c(j-1)
                  lowerVX_n_vr_acc(j,s,i) = lowerVX_n_vr_acc(j,s,i) - a_tri(j) / dz_soi(j) * deltim * conc_trcr_n(j-1)
                  if(j .ne. nl_soil)then
                     upperVX_c_vr_acc(j,s,i) = upperVX_c_vr_acc(j,s,i) - c_tri(j) / dz_soi(j) * deltim * conc_trcr_c(j+1)
                     upperVX_n_vr_acc(j,s,i) = upperVX_n_vr_acc(j,s,i) - c_tri(j) / dz_soi(j) * deltim * conc_trcr_n(j+1)
                  end if
                  diagVX_c_vr_acc(j,s,i) = diagVX_c_vr_acc(j,s,i) + (b_tri(j) - a_p_0) / dz_soi(j) * deltim * conc_trcr_c(j)
                  diagVX_n_vr_acc(j,s,i) = diagVX_n_vr_acc(j,s,i) + (b_tri(j) - a_p_0) / dz_soi(j) * deltim * conc_trcr_n(j)
               else
                  if(j .eq. nbedrock + 1 .and. j .ne. nl_soil .and. j .gt. 1)then
                     diagVX_c_vr_acc(j-1,s,i) = diagVX_c_vr_acc(j-1,s,i) + a_tri(j) / dz_soi(j-1) * deltim * conc_trcr_c(j-1)
                     diagVX_n_vr_acc(j-1,s,i) = diagVX_n_vr_acc(j-1,s,i) + a_tri(j) / dz_soi(j-1) * deltim * conc_trcr_n(j-1)
                  end if
               end if
#endif
!               if(s .eq. 1 .and. i_type .eq. 1 .and. use_soil_matrixcn)then
!                  if(j .le. col%nbedrock)then
!                     do i = 1,ndecomp_pools-1
!                        tri_ma_vr(c,j*3-4+(i-1)*(nl_soil*3-2)) = a_tri(c,j) / dz_soi(j) * (-deltim)
!                        if(j .ne. nl_soil)then
!                           tri_ma_vr(c,j*3  +(i-1)*(nl_soil*3-2)) = c_tri(c,j) / dz_soi(j) * (-deltim)
!                        end if
!                        tri_ma_vr(c,j*3-2+(i-1)*(nl_soil*3-2)) = (b_tri(c,j) - a_p_0) / dz_soi(j) * (-deltim)
!                     end do
!                  else
!                     if(j .eq. col%nbedrock + 1 .and. j .ne. nl_soil .and. j .gt. 1)then
!                        do i = 1,ndecomp_pools-1
!                           tri_ma_vr(c,(j-1)*3-2+(i-1)*(nl_soil*3-2)) = tri_ma_vr(c,(j-1)*3-2+(i-1)*(nl_soil*3-2)) &
!                                                                        + a_tri(c,j) / dz_soi(j-1)*(-deltim)
!                        end do
!                     end if
!                  end if
!               end if
            else ! j==nl_soil+1; 0 concentration gradient at bottom
               a_tri(j)   = -1._r8
               b_tri(j)   = 1._r8
               c_tri(j)   = 0._r8
               r_tri_c(j) = 0._r8
               r_tri_n(j) = 0._r8
            endif
         enddo ! j; nl_soil

         jtop = 0

         ! subtract initial concentration and source terms for tendency calculation
         do j = 1, nl_soil
!            if (.not. use_soil_matrixcn) then
            decomp_cpools_transport_tendency(j,s,i) = 0.-(conc_trcr_c(j) + decomp_cpools_sourcesink(j,s,i))
            decomp_npools_transport_tendency(j,s,i) = 0.-(conc_trcr_n(j) + decomp_npools_sourcesink(j,s,i))
!            else
!               decomp_cpools_transport_tendency(c,j,s) = 0.0_r8
!           end if !soil_matrix
!            if(i .eq. 79738)print*,'tendency in LittVertTransp',j,s,i,decomp_cpools_transport_tendency(j,s,i),&
!            conc_trcr_c(j), decomp_cpools_sourcesink(j,s,i)
         end do

!         if (.not. use_soil_matrixcn) then
         ! Solve for the concentration profile for this time step
         !print*,'before tridia,rc',s,r_tri_c
         !print*,'before tridia,rn',s,r_tri_n
         !print*,'conc_trcr_c',s,conc_trcr_c
         !print*,'conc_trcr_n',s,conc_trcr_n
         call tridia(nl_soil+2, a_tri  (:), b_tri(:), c_tri(:), r_tri_c(:), conc_trcr_c(0:nl_soil+1))
         call tridia(nl_soil+2, a_tri  (:), b_tri(:), c_tri(:), r_tri_n(:), conc_trcr_n(0:nl_soil+1))

!         if(i .eq. 79738)print*,'tendency in after tridia',j,s,i,decomp_cpools_transport_tendency(j,s,i),&
!         conc_trcr_c(j)
         ! add post-transport concentration to calculate tendency term
         do j = 1, nl_soil
            decomp_cpools_transport_tendency(j,s,i) = decomp_cpools_transport_tendency(j,s,i) + conc_trcr_c(j)
            decomp_cpools_transport_tendency(j,s,i) = decomp_cpools_transport_tendency(j,s,i) / deltim
            decomp_npools_transport_tendency(j,s,i) = decomp_npools_transport_tendency(j,s,i) + conc_trcr_n(j)
            decomp_npools_transport_tendency(j,s,i) = decomp_npools_transport_tendency(j,s,i) / deltim
         end do
!         else
!            do j = 1,nl_soil
!               do fc =1,num_soilc
!                  c = filter_soilc(fc)
!                  matrix_input(c,j+(s-1)*nl_soil) = matrix_input(c,j+(s-1)*nl_soil) + source(c,j,s)
!               end do
!            end do
!         end if  !soil_matrix
      else
         ! for CWD pools, just add
         do j = 1,nl_soil
 !           if(.not. use_soil_matrixcn)then
            conc_trcr_c(j) = decomp_cpools_vr(j,s,i) + decomp_cpools_sourcesink(j,s,i)
            conc_trcr_n(j) = decomp_npools_vr(j,s,i) + decomp_npools_sourcesink(j,s,i)
 !           else
 !              matrix_input(c,j+(s-1)*nl_soil) = matrix_input(c,j+(s-1)*nl_soil) + source(c,j,s)
 !           end if
            if (j > nbedrock .and. decomp_cpools_sourcesink(j,s,i) > 0._r8) then
               write(*,*) 'C source >0',i,j,s,decomp_cpools_sourcesink(j,s,i)
            end if
            if (j > nbedrock .and. decomp_cpools_vr(j,s,i) > 0._r8) then
               write(*,*) 'C conc_ptr >0',i,j,s,decomp_cpools_vr(j,s,i)
            end if
            if (j > nbedrock .and. decomp_npools_sourcesink(j,s,i) > 0._r8) then
               write(*,*) 'N source >0',i,j,s,decomp_npools_sourcesink(j,s,i)
            end if
            if (j > nbedrock .and. decomp_npools_vr(j,s,i) > 0._r8) then
               write(*,*) 'N conc_ptr >0',i,j,s,decomp_npools_vr(j,s,i)
            end if
         end do
      end if ! not CWD

 !    if (.not. use_soil_matrixcn) then
     do j = 1,nl_soil
        decomp_cpools_vr(j,s,i) = conc_trcr_c(j)
        decomp_npools_vr(j,s,i) = conc_trcr_n(j)
        !print*,'updating decomp_cpools_vr',j,s,decomp_cpools_vr(j,s,i)
        !print*,'updating decomp_npools_vr',j,s,decomp_npools_vr(j,s,i)
        ! Correct for small amounts of carbon that leak into bedrock
        if (j > nbedrock) then
           decomp_cpools_vr(nbedrock,s,i) = decomp_cpools_vr(nbedrock,s,i) + &
              conc_trcr_c(j) * (dz_soi(j) / dz_soi(nbedrock))
           decomp_cpools_vr(j,s,i) = 0._r8
           decomp_npools_vr(nbedrock,s,i) = decomp_npools_vr(nbedrock,s,i) + &
              conc_trcr_c(j) * (dz_soi(j) / dz_soi(nbedrock))
           decomp_npools_vr(j,s,i) = 0._r8
        end if
     end do
!               end if !not soil_matrix
  end do ! s (pool loop)


end subroutine SoilBiogeochemLittVertTransp

end module bgc_soil_SoilBiogeochemLittVertTranspMod
#endif
