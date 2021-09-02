
#include <define.h>

SUBROUTINE rd_forcing(idate,lon_points,lat_points,solarin_all_band,numpatch)

  use precision
  use PhysicalConstants, only: rgas, grav
  use MOD_TimeInvariants
  use MOD_1D_Forcing
  use MOD_2D_Forcing
  use GETMETMOD, only: forcn, GETMET
  use timemanager
  use omp_lib

      IMPLICIT NONE
      integer,  INTENT(in) :: idate(3)
      integer,  INTENT(in) :: lon_points
      integer,  INTENT(in) :: lat_points
      integer,  INTENT(in) :: numpatch

      logical,  INTENT(in) :: solarin_all_band

! local variables:
      integer  :: i, j, np
!added by yuan, 07/06/2016
      real(r8) :: pi      ! pie
      real(r8) :: coszen  ! cosine of solar zenith angle
      real(r8) :: calday  ! Julian cal day (1.xx to 365.xx)
      real(r8) :: sunang, cloud, difrat, vnrat
      real(r8) :: forc_xy_solarin(lon_points,lat_points)
      real(r8) :: a, hsolar, ratio_rvrf 

      real(r8), external :: orb_coszen

      real solar, frl, prcp, tm, us, vs, pres, qm

!added by yuan, 07/06/2016
      pi = 4.*atan(1.)        

!------------------------------------------------------------
    ! READ IN THE ATMOSPHERIC FORCING
      
      CALL GETMET(idate)

      forc_xy_t      (:,:) = forcn(:,:,1)
      forc_xy_q      (:,:) = forcn(:,:,2)
      forc_xy_psrf   (:,:) = forcn(:,:,3)
      forc_xy_pbot   (:,:) = forcn(:,:,3)
      forc_xy_solarin(:,:) = forcn(:,:,7)
      forc_xy_frl    (:,:) = forcn(:,:,8)

#if(defined USE_POINT_DATA)
      forc_xy_prl    (:,:) = forcn(:,:,4) * 2/3. ! ?
      forc_xy_prc    (:,:) = forcn(:,:,4) * 1/3. ! ?
      forc_xy_us     (:,:) = forcn(:,:,5)
      forc_xy_vs     (:,:) = forcn(:,:,6)
#else
      forc_xy_prl    (:,:) = forcn(:,:,4) * 2/3.
      forc_xy_prc    (:,:) = forcn(:,:,4) * 1/3.
      forc_xy_us     (:,:) = forcn(:,:,6)/sqrt(2.0)
      forc_xy_vs     (:,:) = forcn(:,:,6)/sqrt(2.0)
#endif


#if(defined HEIGHT_V)
      forc_xy_hgt_u  (:,:) = HEIGHT_V
#else
      forc_xy_hgt_u  (:,:) = 100.
#endif

#if(defined HEIGHT_T)
      forc_xy_hgt_t  (:,:) = HEIGHT_T
#else
      forc_xy_hgt_t  (:,:) = 50.
#endif

#if(defined HEIGHT_Q)
      forc_xy_hgt_q  (:,:) = HEIGHT_Q
#else
      forc_xy_hgt_q  (:,:) = 50.
#endif


      calday = calendarday(idate, gridlond(1)) 
      
      if(solarin_all_band)then

#if(defined USE_QIAN_DATA)
         !---------------------------------------------------------------
         ! NOTE: codes from CLM4.5-CESM1.2.0
         ! relationship between incoming NIR or VIS radiation and ratio of
         ! direct to diffuse radiation calculated based on one year's worth of
         ! hourly CAM output from CAM version cam3_5_55
         !---------------------------------------------------------------
! added by yuan, 06/20/2016
#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,hsolar,ratio_rvrf)
#endif
         do j = 1, lat_points
            do i = 1, lon_points

               hsolar = forc_xy_solarin(i,j)*0.5_R8

             ! NIR (dir, diff)
               ratio_rvrf = min(0.99_R8,max(0.29548_R8 + 0.00504_R8*hsolar  &
                  -1.4957e-05_R8*hsolar**2 + 1.4881e-08_R8*hsolar**3,0.01_R8))
               forc_xy_soll (i,j) = ratio_rvrf*hsolar
               forc_xy_solld(i,j) = (1._R8 - ratio_rvrf)*hsolar

             ! VIS (dir, diff)
               ratio_rvrf = min(0.99_R8,max(0.17639_R8 + 0.00380_R8*hsolar  &
                  -9.0039e-06_R8*hsolar**2 + 8.1351e-09_R8*hsolar**3,0.01_R8))
               forc_xy_sols (i,j) = ratio_rvrf*hsolar
               forc_xy_solsd(i,j) = (1._R8 - ratio_rvrf)*hsolar
            end do
         end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#else
         !---------------------------------------------------------------
         ! as the downward solar is in full band, an empirical expression
         ! will be used to divide fractions of band and incident
         ! (visible, near-infrad, dirct, diffuse)
         ! Julian calday (1.xx to 365.xx)
         !---------------------------------------------------------------
! added by yuan, 07/06/2016
#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,a,sunang,cloud,difrat,vnrat)
#endif
         do j = 1, lat_points
            do i = 1, lon_points
         !do np = 1, numpatch  
         !   i = ixy_patch(np) ! longitude index
         !   j = jxy_patch(np) ! latitude index

               a = forc_xy_solarin(i,j)
! modified by yuan, 07/06/2016
               !sunang = orb_coszen(calday,dlon(np),dlat(np))
               sunang = orb_coszen(calday,gridlond(i)*pi/180.,gridlatd(j)*pi/180.)

               cloud = (1160.*sunang-a)/(963.*sunang)
               cloud = max(cloud,0.)
               cloud = min(cloud,1.)
               cloud = max(0.58,cloud)

               difrat = 0.0604/(sunang-0.0223)+0.0683
               if(difrat.lt.0.) difrat = 0.
               if(difrat.gt.1.) difrat = 1.

               difrat = difrat+(1.0-difrat)*cloud
               vnrat = (580.-cloud*464.)/((580.-cloud*499.)+(580.-cloud*464.))

               forc_xy_sols  (i,j) = a*(1.0-difrat)*vnrat
               forc_xy_soll  (i,j) = a*(1.0-difrat)*(1.0-vnrat)
               forc_xy_solsd (i,j) = a*difrat*vnrat
               forc_xy_solld (i,j) = a*difrat*(1.0-vnrat)
            end do
         end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#endif
      endif

    ! Mapping the 2d atmospheric fields [lon_points]x[lat_points]
    !     -> the 1d vector of subgrid points [numpatch]
! added by yuan, 06/20/2016
#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j)
#endif
      do np = 1, numpatch  
         i = patch2lon(np) ! longitude index
         j = patch2lat(np) ! latitude index

       ! [CO2 concentration = 398.03ppm On March 12, 2014, NOAA MLO recorded]
         forc_pco2m (np) = forc_xy_pbot  (i,j)*398.03e-06 
         forc_po2m  (np) = forc_xy_pbot  (i,j)*0.209       
         forc_us    (np) = forc_xy_us    (i,j)            
         forc_vs    (np) = forc_xy_vs    (i,j)

         forc_t     (np) = forc_xy_t     (i,j)
! The standard measuring conditions for temperature are two meters above the ground
! Scientists have measured the most frigid temperature ever 
! recorded on the continent's eastern highlands: about (180K) colder than dry ice.
         if(forc_t(np) < 180.) forc_t(np) = 180.
! the highest air temp was found in Kuwait 326 K, Sulaibya 2012-07-31; 
! Pakistan, Sindh 2010-05-26; Iraq, Nasiriyah 2011-08-03
         if(forc_t(np) > 326.) forc_t(np) = 326.

         forc_q     (np) = forc_xy_q     (i,j)
         forc_prc   (np) = forc_xy_prc   (i,j)
         forc_prl   (np) = forc_xy_prl   (i,j)
         forc_psrf  (np) = forc_xy_psrf  (i,j)
         forc_pbot  (np) = forc_xy_pbot  (i,j)
         
         forc_sols  (np) = forc_xy_sols  (i,j)
         forc_soll  (np) = forc_xy_soll  (i,j)
         forc_solsd (np) = forc_xy_solsd (i,j)
         forc_solld (np) = forc_xy_solld (i,j)

         forc_frl   (np) = forc_xy_frl   (i,j)
         forc_hgt_u (np) = forc_xy_hgt_u (i,j)
         forc_hgt_t (np) = forc_xy_hgt_t (i,j)
         forc_hgt_q (np) = forc_xy_hgt_q (i,j)
         forc_rhoair(np) = (forc_pbot(np) &
                         - 0.378*forc_q(np)*forc_pbot(np)/(0.622+0.378*forc_q(np)))&
                         / (rgas*forc_t(np))  
      end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#if(defined CLMDEBUG)
      print*,'forc_t     =', minval(forc_t),     maxval(forc_t)
      print*,'forc_q     =', minval(forc_q),     maxval(forc_q)
      print*,'forc_prc   =', minval(forc_prc),   maxval(forc_prc)
      print*,'forc_psrf  =', minval(forc_psrf),  maxval(forc_psrf)
      print*,'forc_prl   =', minval(forc_prl),   maxval(forc_prl)
      print*,'forc_sols  =', minval(forc_sols),  maxval(forc_sols)
      print*,'forc_soll  =', minval(forc_soll),  maxval(forc_soll)
      print*,'forc_solsd =', minval(forc_solsd), maxval(forc_solsd)
      print*,'forc_solld =', minval(forc_solld), maxval(forc_solld)
      print*,'forc_frl   =', minval(forc_frl),   maxval(forc_frl)
#endif


END SUBROUTINE rd_forcing

