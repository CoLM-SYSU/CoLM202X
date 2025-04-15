#include <define.h>

MODULE MOD_PixelsetShared
!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    Shared pixelset refer to two or more pixelsets sharing the same geographic area.
!
!    For example, for patch of crops, multiple crops can be planted on a piece of land.
!    When planting these crops, different irrigation schemes may be used. Thus the water
!    and energy processes have difference in crops and should be modeled independently.
!    By using shared pixelset, crop patch is splitted to two or more shared patches.
!    Each shared patch is assigned with a percentage of area and has its own states.
!
!                Example of shared pixelsets
!        |<------------------- ELEMENT ------------------>| <-- level 1
!        |   subset 1  |       subset 2        | subset 3 | <-- level 2
!                      | subset 2 shared 1 50% |
!                      | subset 2 shared 2 20% |            <-- subset 2 shares
!                      | subset 2 shared 3 30% |
!
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   IMPLICIT NONE

CONTAINS

   SUBROUTINE pixelsetshared_build (pixelset, gshared, datashared, nmaxshared, typfilter, &
         fracout, sharedclass, fracin)

   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Pixel
   USE MOD_Pixelset
   USE MOD_Mesh
   USE MOD_Utils
   USE MOD_AggregationRequestData
   IMPLICIT NONE

   type(pixelset_type),       intent(inout) :: pixelset
   type(grid_type),           intent(in)    :: gshared
   type(block_data_real8_3d), intent(in)    :: datashared
   integer, intent(in) :: nmaxshared
   integer, intent(in) :: typfilter(:)

   real(r8), intent(out), allocatable :: fracout(:)
   integer,  intent(out), allocatable :: sharedclass(:)
   real(r8), intent(in),  optional    :: fracin (:)

   ! Local Variables
   real(r8), allocatable :: pctshared(:,:)
   real(r8), allocatable :: datashared1d(:,:), areapixel(:), rbuff(:,:)
   integer  :: nsetshared, ipset, jpset
   integer  :: ipxl, ie, ipxstt, ipxend, ishared
   integer*8,allocatable :: eindex1(:)
   integer,  allocatable :: ielm1(:), ipxstt1(:), ipxend1(:), settyp1(:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
      IF (p_is_io) THEN
         CALL aggregation_data_daemon (gshared, data_r8_3d_in1 = datashared, n1_r8_3d_in1 = nmaxshared)
      ENDIF
#endif

      IF (p_is_worker) THEN

         nsetshared = 0

         allocate (pctshared(nmaxshared,pixelset%nset))

         DO ipset = 1, pixelset%nset
            IF (any(typfilter(:) == pixelset%settyp(ipset))) THEN

               ie     = pixelset%ielm  (ipset)
               ipxstt = pixelset%ipxstt(ipset)
               ipxend = pixelset%ipxend(ipset)

               allocate (datashared1d (nmaxshared, ipxstt:ipxend))

               CALL aggregation_request_data (pixelset, ipset, gshared, zip = .false., &
                  data_r8_3d_in1 = datashared, data_r8_3d_out1 = rbuff, n1_r8_3d_in1 = nmaxshared)

               datashared1d = rbuff

               allocate (areapixel(ipxstt:ipxend))
               DO ipxl = ipxstt, ipxend
                  areapixel(ipxl) = areaquad (&
                     pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
                     pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)) )
               ENDDO

               DO ishared = 1, nmaxshared
                  pctshared(ishared,ipset) = sum(datashared1d(ishared,:) * areapixel)
               ENDDO

               IF (any(pctshared(:,ipset) > 0.)) THEN
                  nsetshared = nsetshared + count(pctshared(:,ipset) > 0.)
                  pctshared(:,ipset) = pctshared(:,ipset) / sum(pctshared(:,ipset))
               ENDIF

               deallocate (rbuff       )
               deallocate (areapixel   )
               deallocate (datashared1d)

            ELSE
               nsetshared = nsetshared + 1
            ENDIF

         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

      IF (p_is_worker) THEN

         IF (pixelset%nset > 0) THEN

            allocate (eindex1(pixelset%nset))
            allocate (ipxstt1(pixelset%nset))
            allocate (ipxend1(pixelset%nset))
            allocate (settyp1(pixelset%nset))
            allocate (ielm1  (pixelset%nset))

            eindex1 = pixelset%eindex
            ipxstt1 = pixelset%ipxstt
            ipxend1 = pixelset%ipxend
            settyp1 = pixelset%settyp
            ielm1   = pixelset%ielm

            deallocate (pixelset%eindex)
            deallocate (pixelset%ipxstt)
            deallocate (pixelset%ipxend)
            deallocate (pixelset%settyp)
            deallocate (pixelset%ielm  )

            allocate (pixelset%eindex(nsetshared))
            allocate (pixelset%ipxstt(nsetshared))
            allocate (pixelset%ipxend(nsetshared))
            allocate (pixelset%settyp(nsetshared))
            allocate (pixelset%ielm  (nsetshared))

            allocate (fracout    (nsetshared))
            allocate (sharedclass(nsetshared))

            fracout(:) = 1.0

            jpset = 0
            DO ipset = 1, pixelset%nset
               IF (any(typfilter(:) == settyp1(ipset))) THEN
                  IF (any(pctshared(:,ipset) > 0.)) THEN
                     DO ishared = 1, nmaxshared
                        IF (pctshared(ishared,ipset) > 0.) THEN
                           jpset = jpset + 1
                           pixelset%eindex(jpset) = eindex1(ipset)
                           pixelset%ipxstt(jpset) = ipxstt1(ipset)
                           pixelset%ipxend(jpset) = ipxend1(ipset)
                           pixelset%settyp(jpset) = settyp1(ipset)
                           pixelset%ielm  (jpset) = ielm1  (ipset)

                           IF (present(fracin)) THEN
                              fracout(jpset) = fracin(ipset) * pctshared(ishared,ipset)
                           ELSE
                              fracout(jpset) = pctshared(ishared,ipset)
                           ENDIF

                           sharedclass(jpset) = ishared
                        ENDIF
                     ENDDO
                  ENDIF
               ELSE
                  jpset = jpset + 1
                  pixelset%eindex(jpset) = eindex1(ipset)
                  pixelset%ipxstt(jpset) = ipxstt1(ipset)
                  pixelset%ipxend(jpset) = ipxend1(ipset)
                  pixelset%settyp(jpset) = settyp1(ipset)
                  pixelset%ielm  (jpset) = ielm1  (ipset)

                  IF (present(fracin)) THEN
                     fracout(jpset) = fracin(ipset)
                  ELSE
                     fracout(jpset) = 1.
                  ENDIF

                  sharedclass(jpset) = 0 ! no meaning
               ENDIF
            ENDDO

            pixelset%nset = nsetshared

            deallocate (eindex1)
            deallocate (ipxstt1)
            deallocate (ipxend1)
            deallocate (settyp1)
            deallocate (ielm1  )
            deallocate (pctshared)

         ENDIF

      ENDIF

      CALL pixelset%set_vecgs

   END SUBROUTINE pixelsetshared_build

END MODULE MOD_PixelsetShared
