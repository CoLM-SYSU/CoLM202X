#include <define.h>
#ifdef BGC
module MOD_BGC_Soil_BiogeochemLittVertTransp

  !----------------------------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Simulate the soil and litter CN veritical mixing (diffusion and advection) processes. Solve the dynamics
  ! of soil and litter vertical profile with a tridiagonal matrix.
  !
  ! !REFERENCE:
  ! Koven, C.D., Riley, W.J., Subin, Z.M., Tang, J.Y., Torn, M.S., Collins, W.D., Bonan, G.B., Lawrence,
  ! D.M. and Swenson, S.C., 2013. The effect of vertically resolved soil biogeochemistry and alternate
  ! soil C and N models on C dynamics of CLM4. Biogeosciences, 10(11), 7109-7131.
  ! Thornton, P.E., Law, B.E., Gholz, H.L., Clark, K.L., Falge, E., Ellsworth, D.S., Goldstein, A.H., Monson,
  ! R.K., Hollinger, D., Falk, M. and Chen, J., 2002. Modeling and measuring the effects of disturbance
  ! history and climate on carbon and water budgets in evergreen needleleaf forests.
  ! Agricultural and forest meteorology, 113(1-4), 185-222.
  !
  ! !ORIGINAL:
  ! The Community Land Model version 5.0 (CLM5)
  !
  ! !REVISION:
  ! Xingjie Lu, 2021, 1) Revised the CLM5 code to be compatible with CoLM code sturcture.
  !                   2) Record accumulated organic CN vertical transfer rates for semi-analytic spin-up.

  use MOD_Precision
  use MOD_Namelist, only : DEF_USE_SASU
  use MOD_BGC_Vars_TimeInvariants, only: &
      is_cwd, som_adv_flux, som_diffus, cryoturb_diffusion_k, max_altdepth_cryoturbation, max_depth_cryoturb
  use MOD_BGC_Vars_TimeVariables, only: &
      altmax, altmax_lastyear, som_adv_coef, som_diffus_coef, &
      decomp_cpools_vr, decomp_npools_vr, &
      diagVX_c_vr_acc, upperVX_c_vr_acc, lowerVX_c_vr_acc, &
      diagVX_n_vr_acc, upperVX_n_vr_acc, lowerVX_n_vr_acc
  use MOD_BGC_Vars_1DFluxes, only: &
      decomp_cpools_sourcesink, decomp_npools_sourcesink, &
      decomp_cpools_transport_tendency, decomp_npools_transport_tendency
  USE MOD_Utils, only: tridia

  implicit none

  public SoilBiogeochemLittVertTransp

contains

  subroutine SoilBiogeochemLittVertTransp(i,deltim,nl_soil,nl_soil_full,ndecomp_pools,nbedrock,z_soi,zi_soi,dz_soi)

    integer ,intent(in) :: i                        ! patch index
    real(r8),intent(in) :: deltim                   ! time step in seconds
    integer ,intent(in) :: nl_soil                  ! number of total soil layers
    integer ,intent(in) :: nl_soil_full             ! number of total soil layers plus bedrock layers
    integer ,intent(in) :: ndecomp_pools            ! number of total soil & litter pools in the decompositions
    integer ,intent(in) :: nbedrock                 ! where bedrock layer starts
    real(r8),intent(in) :: z_soi (1:nl_soil_full)   ! depth of each soil layer (m)
    real(r8),intent(in) :: zi_soi(0:nl_soil_full)   ! interface level below a zsoi level (m)
    real(r8),intent(in) :: dz_soi(1:nl_soil_full)   ! thicknesses of each soil layer (m)

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
          do j = 1,nl_soil+1
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
                if(DEF_USE_SASU)then
                   upperVX_c_vr_acc(j,s,i) = upperVX_c_vr_acc(j,s,i) -  c_tri(j)          / dz_soi(j) * deltim * conc_trcr_c(j+1)! upwards transfer
                   diagVX_c_vr_acc (j,s,i) = diagVX_c_vr_acc (j,s,i) + (b_tri(j) - a_p_0) / dz_soi(j) * deltim * conc_trcr_c(j)! exit flux
                   upperVX_n_vr_acc(j,s,i) = upperVX_n_vr_acc(j,s,i) -  c_tri(j)          / dz_soi(j) * deltim * conc_trcr_n(j+1)! upwards transfer
                   diagVX_n_vr_acc (j,s,i) = diagVX_n_vr_acc (j,s,i) + (b_tri(j) - a_p_0) / dz_soi(j) * deltim * conc_trcr_n(j)! exit flux
                end if
             elseif (j < nl_soil+1) then

                a_tri(j) = -(d_m1_zm1(j) * aaa(pe_m1(j)) + max( f_m1(j), 0._r8)) ! Eqn 5.47 Patankar
                c_tri(j) = -(d_p1_zp1(j) * aaa(pe_p1(j)) + max(-f_p1(j), 0._r8))
                b_tri(j) = - a_tri(j) - c_tri(j) + a_p_0
                r_tri_c(j) = decomp_cpools_sourcesink(j,s,i) * dz_soi(j) /deltim + a_p_0 * conc_trcr_c(j)
                r_tri_n(j) = decomp_npools_sourcesink(j,s,i) * dz_soi(j) /deltim + a_p_0 * conc_trcr_n(j)

                if(DEF_USE_SASU)then
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
                end if
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
             decomp_cpools_transport_tendency(j,s,i) = 0.-(conc_trcr_c(j) + decomp_cpools_sourcesink(j,s,i))
             decomp_npools_transport_tendency(j,s,i) = 0.-(conc_trcr_n(j) + decomp_npools_sourcesink(j,s,i))
          end do

          call tridia(nl_soil+2, a_tri  (:), b_tri(:), c_tri(:), r_tri_c(:), conc_trcr_c(0:nl_soil+1))
          call tridia(nl_soil+2, a_tri  (:), b_tri(:), c_tri(:), r_tri_n(:), conc_trcr_n(0:nl_soil+1))

          ! add post-transport concentration to calculate tendency term
          do j = 1, nl_soil
             decomp_cpools_transport_tendency(j,s,i) = decomp_cpools_transport_tendency(j,s,i) + conc_trcr_c(j)
             decomp_cpools_transport_tendency(j,s,i) = decomp_cpools_transport_tendency(j,s,i) / deltim
             decomp_npools_transport_tendency(j,s,i) = decomp_npools_transport_tendency(j,s,i) + conc_trcr_n(j)
             decomp_npools_transport_tendency(j,s,i) = decomp_npools_transport_tendency(j,s,i) / deltim
          end do
       else
          ! for CWD pools, just add
          do j = 1,nl_soil
             conc_trcr_c(j) = decomp_cpools_vr(j,s,i) + decomp_cpools_sourcesink(j,s,i)
             conc_trcr_n(j) = decomp_npools_vr(j,s,i) + decomp_npools_sourcesink(j,s,i)
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

       do j = 1,nl_soil
          decomp_cpools_vr(j,s,i) = conc_trcr_c(j)
          decomp_npools_vr(j,s,i) = conc_trcr_n(j)
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
    end do ! s (pool loop)

  end subroutine SoilBiogeochemLittVertTransp

end module MOD_BGC_Soil_BiogeochemLittVertTransp
#endif
