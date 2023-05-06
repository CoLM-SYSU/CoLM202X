#include <define.h>

MODULE mod_pixelsetshadow

   IMPLICIT NONE

CONTAINS

   SUBROUTINE pixelsetshadow_build (pixelset, gshadow, datashadow, nshadow, ntypshadow, filter, &
         fracout, shadowclass, fracin)

      USE spmd_task
      USE mod_grid
      USE mod_data_type
      USE mod_pixel
      USE mod_pixelset
      USE mod_mesh
      USE mod_utils
      USE mod_aggregation
      IMPLICIT NONE

      TYPE(pixelset_type),       intent(inout) :: pixelset
      TYPE(grid_type),           intent(in)    :: gshadow
      TYPE(block_data_real8_3d), intent(in)    :: datashadow
      INTEGER, intent(in) :: nshadow, ntypshadow
      INTEGER, intent(in) :: filter(ntypshadow)

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
         CALL aggregation_data_daemon (gshadow, data_r8_3d_in1 = datashadow, n1_r8_3d_in1 = nshadow)
      ENDIF
#endif
         
      IF (p_is_worker) THEN

         nsetshadow = 0

         allocate (pctshadow(nshadow,pixelset%nset))

         DO ipset = 1, pixelset%nset
            IF (any(filter(:) == pixelset%settyp(ipset))) THEN

               ie     = pixelset%ielm  (ipset)
               ipxstt = pixelset%ipxstt(ipset)
               ipxend = pixelset%ipxend(ipset)
      
               allocate (datashadow1d (nshadow, ipxstt:ipxend))

               CALL aggregation_request_data (pixelset, ipset, gshadow, &
                  data_r8_3d_in1 = datashadow, data_r8_3d_out1 = rbuff, n1_r8_3d_in1 = nshadow)

               datashadow1d = rbuff

               allocate (areapixel(ipxstt:ipxend))
               DO ipxl = ipxstt, ipxend
                  areapixel(ipxl) = areaquad (&
                     pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
                     pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)) )
               ENDDO

               DO ishadow = 1, nshadow
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
               IF (any(filter(:) == settyp1(ipset))) THEN
                  IF (any(pctshadow(:,ipset) > 0.)) THEN
                     DO ishadow = 1, nshadow
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

END MODULE mod_pixelsetshadow
