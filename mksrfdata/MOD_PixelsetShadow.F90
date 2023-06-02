#include <define.h>

MODULE MOD_PixelsetShadow
   !----------------------------------------------------------------------------------------
   ! DESCRIPTION:
   !
   !    Shadows of pixelset refer to two or more pixelsets sharing the same geographic area.
   ! 
   !    For example, for patch of crops, multiple crops can be planted on a piece of land.
   !    When planting these crops, different irrigation schemes may be used. Thus the water 
   !    and energy processes have difference in crops and should be modeled independently.
   !    By using shadows, crop patch is splitted to two or more shadowed patches.
   !    Each shadow is assigned with a percentage of area and has its own states.
   !
   !                Example of shadowed pixelsets
   !        |<------------------- ELEMENT ------------------>| <-- level 1
   !        |   subset 1  |       subset 2        | subset 3 | <-- level 2
   !                      | subset 2 shadow 1 50% |            
   !                      | subset 2 shadow 2 20% |            <-- subset 2 shadows
   !                      | subset 2 shadow 3 30% |            
   !
   !
   ! Created by Shupeng Zhang, May 2023
   !----------------------------------------------------------------------------------------
   IMPLICIT NONE

CONTAINS

   SUBROUTINE pixelsetshadow_build (pixelset, gshadow, datashadow, nmaxshadow, typfilter, &
         fracout, shadowclass, fracin)

      USE MOD_SPMD_Task
      USE MOD_Grid
      USE MOD_DataType
      USE MOD_Pixel
      USE MOD_Pixelset
      USE MOD_Mesh
      USE MOD_Utils
      USE MOD_AggregationRequestData
      IMPLICIT NONE

      TYPE(pixelset_type),       intent(inout) :: pixelset
      TYPE(grid_type),           intent(in)    :: gshadow
      TYPE(block_data_real8_3d), intent(in)    :: datashadow
      INTEGER, intent(in) :: nmaxshadow
      INTEGER, intent(in) :: typfilter(:)

      REAL(r8), intent(out), allocatable :: fracout(:)
      INTEGER,  intent(out), allocatable :: shadowclass(:)
      REAL(r8), intent(in),  optional    :: fracin (:)

      ! Local Variables
      REAL(r8), allocatable :: pctshadow(:,:)
      REAL(r8), allocatable :: datashadow1d(:,:), areapixel(:), rbuff(:,:)
      INTEGER  :: nsetshadow, ipset, jpset
      INTEGER  :: ipxl, ie, ipxstt, ipxend, ishadow
      INTEGER,  allocatable :: eindex1(:), ielm1(:), ipxstt1(:), ipxend1(:), settyp1(:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
         
#ifdef USEMPI
      IF (p_is_io) THEN
         CALL aggregation_data_daemon (gshadow, data_r8_3d_in1 = datashadow, n1_r8_3d_in1 = nmaxshadow)
      ENDIF
#endif
         
      IF (p_is_worker) THEN

         nsetshadow = 0

         allocate (pctshadow(nmaxshadow,pixelset%nset))

         DO ipset = 1, pixelset%nset
            IF (any(typfilter(:) == pixelset%settyp(ipset))) THEN

               ie     = pixelset%ielm  (ipset)
               ipxstt = pixelset%ipxstt(ipset)
               ipxend = pixelset%ipxend(ipset)
      
               allocate (datashadow1d (nmaxshadow, ipxstt:ipxend))

               CALL aggregation_request_data (pixelset, ipset, gshadow, &
                  data_r8_3d_in1 = datashadow, data_r8_3d_out1 = rbuff, n1_r8_3d_in1 = nmaxshadow)

               datashadow1d = rbuff

               allocate (areapixel(ipxstt:ipxend))
               DO ipxl = ipxstt, ipxend
                  areapixel(ipxl) = areaquad (&
                     pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
                     pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)) )
               ENDDO

               DO ishadow = 1, nmaxshadow
                  pctshadow(ishadow,ipset) = sum(datashadow1d(ishadow,:) * areapixel)
               ENDDO

               IF (any(pctshadow(:,ipset) > 0.)) THEN
                  nsetshadow = nsetshadow + count(pctshadow(:,ipset) > 0.)
                  pctshadow(:,ipset) = pctshadow(:,ipset) / sum(pctshadow(:,ipset))
               ENDIF

               deallocate (rbuff       )
               deallocate (areapixel   )
               deallocate (datashadow1d)

            ELSE
               nsetshadow = nsetshadow + 1
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

            allocate (pixelset%eindex(nsetshadow))
            allocate (pixelset%ipxstt(nsetshadow))
            allocate (pixelset%ipxend(nsetshadow))
            allocate (pixelset%settyp(nsetshadow))
            allocate (pixelset%ielm  (nsetshadow))

            allocate (fracout    (nsetshadow))
            allocate (shadowclass(nsetshadow))

            jpset = 0
            DO ipset = 1, pixelset%nset
               IF (any(typfilter(:) == settyp1(ipset))) THEN
                  IF (any(pctshadow(:,ipset) > 0.)) THEN
                     DO ishadow = 1, nmaxshadow
                        IF (pctshadow(ishadow,ipset) > 0.) THEN
                           jpset = jpset + 1
                           pixelset%eindex(jpset) = eindex1(ipset)
                           pixelset%ipxstt(jpset) = ipxstt1(ipset)
                           pixelset%ipxend(jpset) = ipxend1(ipset)
                           pixelset%settyp(jpset) = settyp1(ipset)
                           pixelset%ielm  (jpset) = ielm1  (ipset)

                           IF (present(fracin)) THEN
                              fracout(jpset) = fracin(ipset) * pctshadow(ishadow,ipset)
                           ELSE
                              fracout(jpset) = pctshadow(ishadow,ipset)
                           ENDIF

                           shadowclass(jpset) = ishadow
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

                  shadowclass(jpset) = 0 ! no meaning
               ENDIF
            ENDDO

            pixelset%nset = nsetshadow

            deallocate (eindex1)
            deallocate (ipxstt1)
            deallocate (ipxend1)
            deallocate (settyp1)
            deallocate (ielm1  )
            deallocate (pctshadow)

         ENDIF

      ENDIF
         
      CALL pixelset%set_vecgs

   END SUBROUTINE pixelsetshadow_build

END MODULE MOD_PixelsetShadow
