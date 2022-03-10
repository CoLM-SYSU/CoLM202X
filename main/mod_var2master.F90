module mod_var2master
    use precision
    use mod_landpatch, only : numpatch
    use mod_grid
    use mod_mapping_pset2grid
    use mod_namelist
    use GlobalVars, only : spval
 



type :: segment_type
   integer :: blk
   integer :: cnt
   integer :: bdsp
   integer :: gdsp
end type segment_type

integer :: nxseg
type(segment_type), allocatable :: xsegment(:)

integer :: nyseg
type(segment_type), allocatable :: ysegment(:)

integer :: ndatablk



type(grid_type), target :: ghist
type(mapping_pset2grid_type) :: mp2g_hist1

integer :: hist_data_id_1=3000
public :: var_out

contains
    SUBROUTINE var_out (runoff_2d)

        !=======================================================================
        ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
        !=======================================================================

        use precision
        use mod_namelist
        use timemanager
        use spmd_task
        use mod_2d_fluxes
        use MOD_1D_Fluxes,only :rnof
        !use MOD_1D_Acc_Fluxes
        use mod_block
        use mod_data_type
       ! use mod_landset
        use mod_mapping_pset2grid
        use MOD_2D_Fluxes
        use mod_colm_debug
        use GlobalVars, only : spval
        use MOD_TimeInvariants, only: patchtype !, spval

        IMPLICIT NONE
        !integer,  INTENT(in)  :: ilev_patch
        real(r8), INTENT(out) ::  runoff_2d (DEF_nlon_hist, DEF_nlat_hist)
        ! Local variables
        !integer :: numpatch

        type(block_data_real8_2d) :: flux_xy
        type(block_data_real8_2d) :: sumwt
        real(r8), allocatable ::  vectmp(:)  
        logical,  allocatable ::  filter(:)
        integer :: xblk, yblk, xloc, yloc
        integer :: iblk, jblk, idata, ixseg, iyseg
        integer :: rmesg(3), smesg(3), isrc
        real(r8), allocatable :: rbuf(:,:), sbuf(:,:), vdata(:,:)
        integer :: xdsp, ydsp, xcnt, ycnt
        if(p_is_master)then 
            runoff_2d(:,:) = spval
        endif 
        if (p_is_worker) then
            ! numpatch = landset%levs(ilev_patch)%nset
             if (numpatch > 0) then
                 allocate (filter (numpatch))
                 allocate (vectmp (numpatch))
             end if
        print *, 'good 14'
         end if
        if (p_is_io) then
            call allocate_block_data (ghist1, sumwt)
        end if
 
 


    END SUBROUTINE var_out
end module mod_var2master
 
