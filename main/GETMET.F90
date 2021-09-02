#include <define.h>

MODULE GETMETMOD
   
   use precision
   use timemanager
   use omp_lib

   use METDAT, only: NVAR, tstamp_LB, tstamp_UB, tintalgo, avgcos, rlats, rlons, &
                     metinit, metreadLBUB, metreadpoint, metpreprocess, metfinal, sanity

   implicit none
   save
   private
   integer :: lat_i                                 ! start index of latitude
   integer :: lon_i                                 ! start index of longitude
   integer :: lat_n                                 ! number of latitudes
   integer :: lon_n                                 ! number of longitudes
   real(r8), allocatable :: forcn(:,:,:)            ! forcing data 
   real(r8), allocatable :: forcn_LB(:,:,:)         ! forcing data at lower bondary
   real(r8), allocatable :: forcn_UB(:,:,:)         ! forcing data at upper bondary   
   real(r8), external    :: orb_coszen              ! cosine of solar zenith angle

   INTERFACE GETMETINI
      MODULE PROCEDURE GETMETINI
   END INTERFACE

   public forcn
   public GETMET, GETMETINI, GETMETFINAL

CONTAINS
   
   SUBROUTINE GETMETINI(fmetdat, deltim, lat_points, lon_points, lat_start, lon_start)

      implicit none
      character(len=256), intent(in) :: fmetdat     ! directory of forcing data
      real(r8), intent(in) :: deltim                ! model time step
      integer,  intent(in) :: lat_points            ! model latitude points
      integer,  intent(in) :: lon_points            ! model longitude points
      integer,  optional   :: lat_start             ! model latitude start index
      integer,  optional   :: lon_start             ! model longitude start index
    
    ! initialization of forcing data
      call metinit(fmetdat, deltim)
      
    ! initialize the region
    ! number of points
      lat_n = lat_points
      lon_n = lon_points

    ! start points
      if (.NOT. present(lat_start)) then
         lat_i = 1 
      else
         lat_i = lat_start
      end if

      if (.NOT. present(lon_start)) then
         lon_i = 1 
      else
         lon_i = lon_start
      end if


    ! allocate memory
      allocate(forcn(lon_n, lat_n, NVAR))
      allocate(forcn_LB(lon_n, lat_n, NVAR))
      allocate(forcn_UB(lon_n, lat_n, NVAR))

   END SUBROUTINE GETMETINI
   
   SUBROUTINE GETMETFINAL

    ! finalization of forcing data
      call metfinal

    ! deallocate memory
      deallocate(forcn)
      deallocate(forcn_LB)
      deallocate(forcn_UB)

   END SUBROUTINE GETMETFINAL

 ! ------------------------------------------------------------
 ! FUNCTION NAME:
 !    GETMET
 ! 
 ! PURPOSE:
 !    major interface for getting forcing data
 ! ------------------------------------------------------------
   SUBROUTINE GETMET(idate)

      implicit none
      integer,  intent(in) :: idate(3)

#if(defined USE_POINT_DATA)
      CALL metreadpoint(forcn)
#else
      CALL GETMET_NC(idate)
#endif

   END SUBROUTINE GETMET

 ! ------------------------------------------------------------
 ! FUNCTION NAME:
 !    GETMET_NC
 !
 ! PURPOSE:
 !    read forcing data from PRINCETON/GSWP2/QIAN NetCDF data
 !    and interpolate data
 ! ------------------------------------------------------------
   SUBROUTINE GETMET_NC(idate)

      implicit none
      integer,  intent(in) :: idate(3)
   
      type(timestamp) :: mtstamp 
      integer  :: id(3)
      integer  :: i, j, k, dtLB, dtUB
      real(r8) :: calday, cosz
      
    ! read lower and upper boundary forcing data
      CALL metreadLBUB(idate, lat_i, lon_i, lat_n, lon_n, forcn_LB, forcn_UB)

    ! set model time stamp
      id(:) = idate(:)
      call adj2end(id)
      mtstamp = id
      
    ! loop for variables
      do i = 1, NVAR

         if (tintalgo(i) == 'NULL') cycle

       ! to make sure the forcing data calculated is in the range of time 
       ! interval [LB, UB]
         if ( .NOT. (tstamp_LB(i)<=mtstamp .AND. mtstamp<=tstamp_UB(i)) ) then
            write(6, *) "the data required is out of range! stop!"; stop
         end if
         
       ! calcualte distance to lower/upper boundary
         dtLB = mtstamp - tstamp_LB(i)
         dtUB = tstamp_UB(i) - mtstamp

       ! nearest method, for precipitation
         if (tintalgo(i) == 'nearest') then
            if (dtLB <= dtUB) then
               forcn(:,:,i) = forcn_LB(:,:,i)
            else
               forcn(:,:,i) = forcn_UB(:,:,i)
            end if
         end if
         
       ! linear method, for T, Pres, Q, W, LW
         if (tintalgo(i) == 'linear') then
            if ( (dtLB+dtUB) > 0 ) then
               forcn(:,:,i) = & 
                  real(dtUB)/real(dtLB+dtUB)*forcn_LB(:,:,i) + &
                  real(dtLB)/real(dtLB+dtUB)*forcn_UB(:,:,i) 
            else
               forcn(:,:,i) = forcn_LB(:,:,i)
            end if
         end if
         
       ! coszen method, for SW
         if (tintalgo(i) == 'coszen') then
            calday = calendarday(mtstamp)
! added by yuan, 06/20/2016
#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(k,cosz)
#endif
            do j = 1, lat_n
               do k = 1, lon_n
                  cosz = orb_coszen(calday, rlons(k+lon_i-1), rlats(j+lat_i-1))
                  cosz = max(0.001, cosz)
                  forcn(k,j,i) = cosz / avgcos(k+lon_i-1, j+lat_i-1) * forcn_LB(k,j,i)
               end do
            end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
         end if

      end do
    
    ! preprocess for forcing data, only for QIAN data right now?
      CALL metpreprocess(forcn)

   END SUBROUTINE GETMET_NC

END MODULE GETMETMOD
