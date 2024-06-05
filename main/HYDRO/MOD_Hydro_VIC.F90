module MOD_Hydro_VIC
    use MOD_Hydro_VIC_Variables
    implicit none

    public :: compute_vic_runoff

    private :: compute_runoff_and_asat
    private :: calc_Q12
    private :: compute_zwt
    private :: wrap_compute_zwt

    contains

    ! ******************************************************************************
    subroutine Runoff_VIC(deltim, porsl, theta_r, hksati, bsw, &
                          wice_soisno, wliq_soisno, fevpg, rootflux, ppt, &
                          b_infilt, Dsmax, Ds, Ws, c, &
                          rsur,rsubst,wliq_soisno_tmp)

        USE MOD_Namelist
        USE MOD_Precision
        USE MOD_Vars_Global
        implicit none

        !-----------------------Arguments---------------------------------------
        type(soil_con_struct)  :: soil_con
        type(cell_data_struct) :: cell

        real(r8), intent(in) :: porsl(1:nl_soil), theta_r(1:nl_soil), hksati(1:nl_soil), bsw(1:nl_soil)
        real(r8), intent(in) :: wice_soisno(1:nl_soil)
        real(r8), intent(in) :: wliq_soisno(1:nl_soil)
        real(r8), intent(in) :: fevpg
        real(r8), intent(in) :: rootflux(1:nl_soil)
        real(r8), intent(in) :: ppt              ! /**< amount of liquid water coming to the surface   */
        real(r8), intent(in) :: deltim           ! int(DEF_simulation_time%timestep)

        real(r8), intent(in)    :: b_infilt, Dsmax, Ds, Ws, c

        real(r8), intent(inout) :: rsur, rsubst
        real(r8), intent(out)   :: wliq_soisno_tmp(1:nl_soil)

        !-----------------------Local Variables---------------------------------
        integer  :: ilay
        real(r8) :: vic_tmp(Nlayer), vic_tmp_(Nlayer)
        !-----------------------Arguments---------------------------------------

        call vic_para(porsl, theta_r, hksati, bsw, wice_soisno(1:nl_soil), wliq_soisno(1:nl_soil), fevpg, rootflux, &
                        b_infilt, Dsmax, Ds, Ws, c, &
                        soil_con, cell)

        call compute_vic_runoff(soil_con, ppt*deltim, soil_con%frost_fract, cell)

        DO ilay = 1, Nlayer
            vic_tmp(ilay) = cell%layer(ilay)%moist
        ENDDO
        wliq_soisno_tmp = 0.
        call VIC2CoLM(wliq_soisno_tmp, vic_tmp)

        DO ilay = 1, Nlayer
           vic_tmp_(ilay) = sum(cell%layer(ilay)%ice)
        ENDDO
        ! call VIC2CoLM(wice_soisno(1:nl_soil), vic_tmp_)

        if (ppt > 0.) rsur = cell%runoff/deltim
        rsubst = cell%baseflow/deltim

    end subroutine Runoff_VIC

    ! /******************************************************************************
    ! * @brief    Calculate infiltration and runoff from the surface, gravity driven
    ! *           drainage between all soil layers, and generates baseflow from the
    ! *           bottom layer.
    ! ******************************************************************************/
    subroutine compute_vic_runoff(soil_con, ppt, frost_fract, cell)
        use MOD_Hydro_VIC_Variables
        USE MOD_Namelist
        implicit none

        !-----------------------Arguments---------------------------------------
        type(soil_con_struct), intent(in) :: soil_con
        real(r8), intent(in) :: ppt              ! /**< amount of liquid water coming to the surface   */
        real(r8), intent(in) :: frost_fract(:)   ! /**< spatially distributed frost coverage fractions */
        type(cell_data_struct),intent(inout) :: cell

        !-----------------------Local Variables---------------------------------
        integer :: lindex, time_step
        integer :: last_index, tmplayer, fidx
        real(r8) :: A
        real(r8) :: frac
        real(r8) :: tmp_runoff
        real(r8) :: inflow
        real(r8) :: resid_moist(MAX_LAYERS) ! residual moisture (mm)
        real(r8) :: org_moist(MAX_LAYERS) ! total soil moisture (liquid and frozen) at beginning of this function (mm)
        real(r8) :: avail_liq(MAX_LAYERS, MAX_FROST_AREAS) ! liquid soil moisture available for evap/drainage (mm)
        real(r8) :: liq(MAX_LAYERS)
        real(r8) :: ice(MAX_LAYERS)
        real(r8) :: moist(MAX_LAYERS)
        real(r8) :: max_moist(MAX_LAYERS)
        real(r8) :: Ksat(MAX_LAYERS)
        real(r8) :: Q12(MAX_LAYERS - 1)
        real(r8) :: Dsmax
        real(r8) :: tmp_inflow
        real(r8) :: tmp_moist
        real(r8) :: tmp_moist_for_runoff(MAX_LAYERS)
        real(r8) :: tmp_liq
        real(r8) :: dt_inflow
        real(r8) :: dt_runoff
        real(r8) :: runoff(MAX_FROST_AREAS)
        real(r8) :: tmp_dt_runoff(MAX_FROST_AREAS)
        real(r8) :: baseflow(MAX_FROST_AREAS)
        real(r8) :: dt_baseflow
        real(r8) :: rel_moist
        real(r8) :: evap(MAX_LAYERS, MAX_FROST_AREAS)
        real(r8) :: sum_liq
        real(r8) :: evap_fraction
        real(r8) :: evap_sum
        type(layer_data_struct), dimension(MAX_LAYERS) :: layer

        real(r8) :: dltime                !/**< timestep in seconds */
        integer  :: runoff_steps_per_day  !/**< Number of runoff timesteps per day */
        integer  :: model_steps_per_day   !/**< Number of model timesteps per day */
        integer  :: runoff_steps_per_dt

        !-----------------------End Variable List-------------------------------

        dltime = DEF_simulation_time%timestep
        runoff_steps_per_day = 86400/dltime
        model_steps_per_day  = 86400/dltime

        ! /** Set Temporary Variables **/
        do lindex = 1, Nlayer
            resid_moist(lindex) = soil_con%resid_moist(lindex)
            max_moist(lindex) = soil_con%max_moist(lindex)
            Ksat(lindex) = soil_con%Ksat(lindex) / runoff_steps_per_day
        end do

        ! /** Allocate and Set Values for Soil Sublayers **/
        layer = cell%layer
        cell%runoff = 0
        cell%baseflow = 0
        cell%asat = 0

        runoff_steps_per_dt = runoff_steps_per_day / model_steps_per_day

        ! initialize baseflow
        do fidx = 1, Nfrost
            baseflow(fidx) = 0.0
        end do

        do lindex = 1, Nlayer
            evap(lindex, 1) = layer(lindex)%evap / real(runoff_steps_per_dt)
            org_moist(lindex) = layer(lindex)%moist
            layer(lindex)%moist = 0.0

            ! if there is positive evaporation
            if (evap(lindex, 1) > 0.0) then
                sum_liq = 0.0
                ! compute available soil moisture for each frost sub area
                do fidx = 1, Nfrost
                    avail_liq(lindex, fidx) = org_moist(lindex) - layer(lindex)%ice(fidx) - resid_moist(lindex)
                    !avail_liq(lindex, fidx) = org_moist(lindex) -  resid_moist(lindex)
                    if (avail_liq(lindex, fidx) < 0.0) then
                        avail_liq(lindex, fidx) = 0.0
                    end if
                    sum_liq = sum_liq + avail_liq(lindex, fidx) * frost_fract(fidx)
                end do

                ! compute fraction of available soil moisture that is evaporated
                if (sum_liq > 0.0) then
                    evap_fraction = evap(lindex, 1) / sum_liq
                else
                    evap_fraction = 1.0
                end if

                ! distribute evaporation between frost sub areas by percentage
                evap_sum = evap(lindex, 1)
                do fidx = Nfrost, 1, -1
                    evap(lindex, fidx) = avail_liq(lindex, fidx) * evap_fraction
                    avail_liq(lindex, fidx) = avail_liq(lindex, fidx) - evap(lindex, fidx)
                    evap_sum = evap_sum - evap(lindex, fidx) * frost_fract(fidx)
                end do
            else
                ! if no evaporation
                do fidx = Nfrost, 2, -1
                    evap(lindex, fidx) = evap(lindex, 1)
                end do
            end if
        end do


        do fidx = 1, Nfrost
            ! ppt = amount of liquid water coming to the surface
            inflow = ppt

            ! /**************************************************
            ! Initialize Variables
            ! **************************************************/
            do lindex = 1, Nlayer
                ! Set Layer Liquid Moisture Content
                liq(lindex) = org_moist(lindex) - layer(lindex)%ice(fidx)

                ! Set Layer Frozen Moisture Content
                ice(lindex) = layer(lindex)%ice(fidx)
            end do

            ! /******************************************************
            !    Runoff Based on Soil Moisture Level of Upper Layers
            ! ******************************************************/
            do lindex = 1, Nlayer
                tmp_moist_for_runoff(lindex) = liq(lindex) + ice(lindex)
            end do

            call compute_runoff_and_asat(soil_con, tmp_moist_for_runoff, inflow, A, runoff(fidx))

            ! Save dt_runoff based on initial runoff estimate
            tmp_dt_runoff(fidx) = runoff(fidx) / real(runoff_steps_per_dt, kind=r8)

            ! /**************************************************
            !     Compute Flow Between Soil Layers ()
            ! **************************************************/
            dt_inflow = inflow / real(runoff_steps_per_dt, kind=r8)

            Dsmax = soil_con%Dsmax / runoff_steps_per_day

            do time_step = 1, runoff_steps_per_dt
                inflow = dt_inflow

                ! /*************************************
                !     Compute Drainage between Sublayers
                ! *************************************/
                do lindex = 1, Nlayer - 1
                    ! Brooks & Corey relation for hydraulic conductivity
                    tmp_liq = liq(lindex) - evap(lindex, fidx)  ! Assume evap is a 2D array now, adjusted indexing

                    if (tmp_liq < resid_moist(lindex)) then
                        tmp_liq = resid_moist(lindex)
                    end if

                    if (tmp_liq > resid_moist(lindex)) then
                        call calc_Q12(Ksat(lindex), tmp_liq, resid_moist(lindex), max_moist(lindex), soil_con%expt(lindex),Q12(lindex))
                    else
                        Q12(lindex) = 0.0
                    end if
                end do

                ! /**************************************************
                !     Solve for Current Soil Layer Moisture, and
                !     Check Versus Maximum and Minimum Moisture Contents.
                ! **************************************************/
                last_index = 0
                do lindex = 1, Nlayer - 1
                    if (lindex == 1) then
                        dt_runoff = tmp_dt_runoff(fidx)
                    else
                        dt_runoff = 0.0
                    endif

                    ! transport moisture for all sublayers
                    tmp_inflow = 0.0

                    ! Update soil layer moisture content
                    liq(lindex) = liq(lindex) + (inflow - dt_runoff) - (Q12(lindex) + evap(lindex, fidx))

                    ! Verify that soil layer moisture is less than maximum
                    if ((liq(lindex) + ice(lindex)) > max_moist(lindex)) then
                        tmp_inflow = (liq(lindex) + ice(lindex)) - max_moist(lindex)
                        liq(lindex) = max_moist(lindex) - ice(lindex)

                        if (lindex == 1) then
                            Q12(lindex) = Q12(lindex) + tmp_inflow
                            tmp_inflow = 0.0
                        else
                            tmplayer = lindex
                            do while (tmp_inflow > 0)
                                tmplayer = tmplayer - 1
                                if (tmplayer < 1) then
                                    ! If top layer saturated, add to runoff
                                    runoff(fidx) = runoff(fidx) + tmp_inflow
                                    tmp_inflow = 0.0
                                else
                                    ! else add excess soil moisture to next higher layer
                                    liq(tmplayer) = liq(tmplayer) + tmp_inflow
                                    if ((liq(tmplayer) + ice(tmplayer)) > max_moist(tmplayer)) then
                                        tmp_inflow = (liq(tmplayer) + ice(tmplayer)) - max_moist(tmplayer)
                                        liq(tmplayer) = max_moist(tmplayer) - ice(tmplayer)
                                    else
                                        tmp_inflow = 0.0
                                    endif
                                endif
                            end do
                        endif ! /** end trapped excess moisture **/
                    endif ! /** end check if excess moisture in top layer **/

                    ! verify that current layer moisture is greater than minimum
                    if (liq(lindex) < 0.0) then
                        ! liquid cannot fall below 0
                        Q12(lindex) = Q12(lindex) + liq(lindex)
                        liq(lindex) = 0.0
                    endif

                    if ((liq(lindex) + ice(lindex)) < resid_moist(lindex)) then
                        ! moisture cannot fall below minimum
                        Q12(lindex) = Q12(lindex) + (liq(lindex) + ice(lindex)) - resid_moist(lindex)
                        liq(lindex) = resid_moist(lindex) - ice(lindex)
                    endif

                    inflow = Q12(lindex) + tmp_inflow
                    Q12(lindex) = Q12(lindex) + tmp_inflow

                    last_index = last_index + 1
                end do ! /* end loop through soil layers */

                ! /**************************************************
                !     Compute Baseflow
                ! **************************************************/
                ! ARNO model for the bottom soil layer (based on bottom
                ! soil layer moisture from previous time step)

                lindex = Nlayer

                ! Compute relative moisture
                rel_moist = (liq(lindex) - resid_moist(lindex)) / &
                            (max_moist(lindex) - resid_moist(lindex))

                ! Compute baseflow as function of relative moisture
                frac = Dsmax * soil_con%Ds / soil_con%Ws
                dt_baseflow = frac * rel_moist
                if (rel_moist > soil_con%Ws) then
                    frac = (rel_moist - soil_con%Ws) / (1 - soil_con%Ws)
                    dt_baseflow = dt_baseflow + Dsmax * (1 - soil_con%Ds / soil_con%Ws) * &
                                frac ** soil_con%c
                endif

                ! Make sure baseflow isn't negative
                if (dt_baseflow < 0) then
                    dt_baseflow = 0.0
                endif

                ! Extract baseflow from the bottom soil layer
                liq(lindex) = liq(lindex) + Q12(lindex - 1) - (evap(lindex, fidx) + dt_baseflow)

                ! Check Lower Sub-Layer Moistures
                tmp_moist = 0.0

                ! /* If soil moisture has gone below minimum, take water out
                !  * of baseflow and add back to soil to make up the difference
                !  * Note: this may lead to negative baseflow, in which case we will
                !  * reduce evap to make up for it */
                if ((liq(lindex) + ice(lindex)) < resid_moist(lindex)) then
                    dt_baseflow = dt_baseflow + &
                                (liq(lindex) + ice(lindex)) - resid_moist(lindex)
                    liq(lindex) = resid_moist(lindex) - ice(lindex)
                endif

                if ((liq(lindex) + ice(lindex)) > max_moist(lindex)) then
                    ! soil moisture above maximum
                    tmp_moist = (liq(lindex) + ice(lindex)) - max_moist(lindex)
                    liq(lindex) = max_moist(lindex) - ice(lindex)
                    tmplayer = lindex
                    do while (tmp_moist > 0)
                        tmplayer = tmplayer - 1
                        if (tmplayer < 1) then
                            ! If top layer saturated, add to runoff
                            runoff(fidx) = runoff(fidx) + tmp_moist
                            tmp_moist = 0.0
                        else
                            ! else if sublayer exists, add excess soil moisture
                            liq(tmplayer) = liq(tmplayer) + tmp_moist
                            if ((liq(tmplayer) + ice(tmplayer)) > max_moist(tmplayer)) then
                                tmp_moist = (liq(tmplayer) + ice(tmplayer)) - max_moist(tmplayer)
                                liq(tmplayer) = max_moist(tmplayer) - ice(tmplayer)
                            else
                                tmp_moist = 0.0
                            endif
                        endif
                    end do
                endif

                baseflow(fidx) = baseflow(fidx) + dt_baseflow
            end do ! /* end of sub-dt time step loop */

            ! If negative baseflow, reduce evap accordingly
            if (baseflow(fidx) < 0.0) then
                ! layer(lindex)%evap = layer(lindex)%evap + baseflow(fidx)   !!!! need check
                baseflow(fidx) = 0.0
            endif

            ! Recompute Asat based on final moisture level of upper layers
            do lindex = 1, Nlayer
                tmp_moist_for_runoff(lindex) = (liq(lindex) + ice(lindex))
            end do

            call compute_runoff_and_asat(soil_con, tmp_moist_for_runoff, real(0.0, kind=r8), A, tmp_runoff)

            ! Store tile-wide values
            do lindex = 1, Nlayer
                layer(lindex)%moist = layer(lindex)%moist + &
                                    ((liq(lindex) + ice(lindex)) * frost_fract(fidx))
            end do
            cell%asat = cell%asat + A * frost_fract(fidx)
            cell%runoff = cell%runoff + runoff(fidx) * frost_fract(fidx)
            cell%baseflow = cell%baseflow + baseflow(fidx) * frost_fract(fidx)

            ! ! /** Compute water table depth **/
            ! call wrap_compute_zwt(soil_con, cell)

        end do

    end subroutine compute_vic_runoff


    ! ******************************************************************************
    ! * @brief    Calculate the saturated area and runoff
    ! ******************************************************************************
    subroutine compute_runoff_and_asat(soil_con, moist, inflow, A, runoff)
        USE MOD_Hydro_VIC_Variables, only: soil_con_struct, Nlayer
        implicit none
        !-----------------------Arguments---------------------------------------
        type(soil_con_struct), intent(in) :: soil_con
        real(r8), intent(in) :: moist(Nlayer)
        real(r8), intent(in) :: inflow
        real(r8), intent(inout) :: A
        real(r8), intent(inout) :: runoff
        !-----------------------Local Variables---------------------------------
        real(r8) :: top_moist !! total moisture (liquid and frozen) in topmost soil layers (mm)
        real(r8) :: top_max_moist  !! maximum storable moisture (liquid and frozen) in topmost soil layers (mm)
        integer :: lindex
        real(r8) :: ex, max_infil, i_0, basis
        !-----------------------End Variable List-------------------------------

        top_moist = 0.0
        top_max_moist = 0.0
        do lindex = 1, Nlayer - 1
            top_moist = top_moist + moist(lindex)
            top_max_moist = top_max_moist + soil_con%max_moist(lindex)
        end do
        if (top_moist > top_max_moist) then
            top_moist = top_max_moist
        endif

        ! A as in Wood et al. in JGR 97, D3, 1992 equation (1)
        ex = soil_con%b_infilt / (1.0 + soil_con%b_infilt)
        A = 1.0 - (1.0 - top_moist / top_max_moist)**ex

        max_infil = (1.0 + soil_con%b_infilt) * top_max_moist
        i_0 = max_infil * (1.0 - (1.0 - A)**(1.0 / soil_con%b_infilt))

        ! equation (3a) Wood et al.
        if (inflow == 0.0) then
            runoff = 0.0
        else if (max_infil == 0.0) then
            runoff = inflow
        else if ((i_0 + inflow) > max_infil) then
            runoff = inflow - top_max_moist + top_moist
        ! equation (3b) Wood et al. (wrong in paper)
        else
            basis = 1.0 - (i_0 + inflow) / max_infil
            runoff = (inflow - top_max_moist + top_moist + &
                    top_max_moist * basis**(1.0 * (1.0 + soil_con%b_infilt)))
        endif
        if (runoff < 0.0) then
            runoff = 0.0
        endif
    end subroutine compute_runoff_and_asat


    ! ******************************************************************************
    ! * @brief    Calculate drainage between two layers
    ! ******************************************************************************
    subroutine calc_Q12(Ksat, init_moist, resid_moist, max_moist, expt, Q12)
        implicit none
        real(r8), intent(in) :: Ksat, init_moist, resid_moist, max_moist, expt
        real(r8), intent(out) :: Q12

        Q12 = init_moist - ((init_moist - resid_moist)**(1.0d0 - expt) - Ksat / &
            (max_moist - resid_moist)**expt * (1.0d0 - expt))**(1.0d0 / (1.0d0 - expt)) - resid_moist

    end subroutine calc_Q12


    ! /******************************************************************************
    !  * @brief    Compute spatial average water table position (zwt).  Water table
    !  *           position is measured in cm and is negative below the soil surface.
    !  *****************************************************************************/
    subroutine compute_zwt(soil_con,lindex, moist, zwt)
        use MOD_Hydro_VIC_Variables
        implicit none
        !-----------------------Arguments---------------------------------------
        type(soil_con_struct), intent(in) :: soil_con
        integer, intent(in) :: lindex
        real(r8), intent(in) :: moist
        real(r8), intent(out) :: zwt
        !-----------------------Local Variables---------------------------------
        integer :: i
        real(r8) :: MISSING =  -99999. !/**< missing value */
        !-----------------------End Variable List-------------------------------

        zwt = MISSING

        ! /** Compute zwt using soil moisture v zwt curve **/
        i = MAX_ZWTVMOIST - 1
        do while (i >= 1 .and. moist > soil_con%zwtvmoist_moist(lindex, i))
            i = i - 1
        end do

        if (i == MAX_ZWTVMOIST - 1) then
            if (moist < soil_con%zwtvmoist_moist(lindex, i)) then
                zwt = 999.0 ! 999 indicates water table not present in this layer
            else if (moist == soil_con%zwtvmoist_moist(lindex, i)) then
                zwt = soil_con%zwtvmoist_zwt(lindex, i) ! Just barely enough water for a water table
            end if
        else
            zwt = soil_con%zwtvmoist_zwt(lindex, i+1) + &
                    (soil_con%zwtvmoist_zwt(lindex, i) - soil_con%zwtvmoist_zwt(lindex, i+1)) * &
                    (moist - soil_con%zwtvmoist_moist(lindex, i+1)) / &
                    (soil_con%zwtvmoist_moist(lindex, i) - soil_con%zwtvmoist_moist(lindex, i+1))
        end if
    end subroutine compute_zwt


    ! /******************************************************************************
    !  * @brief    Function to compute spatial average water table position (zwt) for
    !  *           individual layers as well as various total-column versions of zwt.
    !  *           Water table position is measured in cm and is negative below the
    !  *           soil surface.
    !  *****************************************************************************/
    subroutine wrap_compute_zwt(soil_con, cell)
        use MOD_Hydro_VIC_Variables
        implicit none

        !-----------------------Arguments---------------------------------------
        type(soil_con_struct), intent(in) :: soil_con
        type(cell_data_struct), intent(inout) :: cell
        !-----------------------Local Variables---------------------------------
        integer :: lindex
        integer :: idx
        real(r8) :: total_depth
        real(r8) :: tmp_depth
        real(r8) :: tmp_moist
        integer, parameter :: CM_PER_M = 100  !/**< centimeters per meter */
        real(r8), parameter :: DBL_EPSILON = 2.2204460492503131E-16
        !-----------------------End Variable List-------------------------------

        ! /** Compute total soil column depth **/
        total_depth = 0.0
        do lindex = 1, Nlayer
            total_depth = total_depth + soil_con%depth(lindex)
        end do

        ! /** Compute each layer's zwt using soil moisture v zwt curve **/
        do lindex = 1, Nlayer
            call compute_zwt(soil_con, lindex, cell%layer(lindex)%moist, cell%layer(lindex)%zwt)
        end do
        if (cell%layer(Nlayer)%zwt == 999) then
            cell%layer(Nlayer)%zwt = -total_depth * CM_PER_M
        end if

        ! /** Compute total soil column's zwt; this will be the zwt of the lowest layer that isn't completely saturated **/
        idx = Nlayer
        tmp_depth = total_depth
        do while (idx >= 1 .and. soil_con%max_moist(idx) - cell%layer(idx)%moist <= DBL_EPSILON)
            tmp_depth = tmp_depth - soil_con%depth(idx)
            idx = idx - 1
        end do
        if (idx < 1) then
            cell%zwt = 0.0
        else if (idx < Nlayer) then
            if (cell%layer(idx)%zwt /= 999) then
                cell%zwt = cell%layer(idx)%zwt
            else
                cell%zwt = -tmp_depth * CM_PER_M
            end if
        else
            cell%zwt = cell%layer(idx)%zwt
        end if

        ! /** Compute total soil column's zwt_lumped; this will be the zwt of all N layers lumped together. **/
        tmp_moist = 0.0
        do lindex = 1, Nlayer
            tmp_moist = tmp_moist + cell%layer(lindex)%moist
        end do
        call compute_zwt(soil_con, Nlayer + 1, tmp_moist, cell%zwt_lumped)

        if (cell%zwt_lumped == 999) then
            cell%zwt_lumped = -total_depth * CM_PER_M   ! // in cm;
        end if
    end subroutine wrap_compute_zwt
end module MOD_Hydro_VIC
