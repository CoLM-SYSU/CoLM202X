#include <define.h>

MODULE mod_pixelsetshadow

   IMPLICIT NONE

CONTAINS

   SUBROUTINE pixelsetshadow_build (pixelset, gshadow, datashadow, nshadow, ntypshadow, filter, &
         fracout, shadowclass, fracin)

      USE spmd_task
      USE mod_block
      USE mod_grid
      USE mod_data_type
      USE mod_pixel
      USE mod_pixelset
      USE mod_landbasin
      USE mod_utils
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
      INTEGER :: nreq, ireq, rmesg(2), smesg(2), isrc, idest, iproc
      INTEGER :: ilon, ilat, xblk, yblk, xloc, yloc
      INTEGER,  allocatable :: ylist(:), xlist(:), ipt(:), ibuf(:)
      REAL(r8), allocatable :: pctshadow(:,:), sbuf(:,:), rbuf(:,:)
      REAL(r8), allocatable :: datashadow1d(:,:), areapixel(:)
      REAL(r8) :: areatotal
      INTEGER  :: nsetshadow, ipset, jpset
      INTEGER  :: ipxl, iu, ipxstt, ipxend, ishadow
      INTEGER,  allocatable :: bindex1(:), ibasin1(:), ipxstt1(:), ipxend1(:), ltyp1(:)
      LOGICAL,  allocatable :: msk(:), worker_done(:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
         
#ifdef USEMPI
      if (p_is_io) then
         
         allocate (worker_done (0:p_np_worker-1))

         worker_done(:) = .false.
         DO while (any(.not. worker_done))

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc = rmesg(1)
            nreq = rmesg(2)

            IF (nreq > 0) THEN

               allocate (xlist (nreq))
               allocate (ylist (nreq))

               CALL mpi_recv (xlist, nreq, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (ylist, nreq, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               
               allocate (sbuf (nshadow,nreq))

               DO ireq = 1, nreq
                  xblk = gshadow%xblk(xlist(ireq))
                  yblk = gshadow%yblk(ylist(ireq))
                  xloc = gshadow%xloc(xlist(ireq))
                  yloc = gshadow%yloc(ylist(ireq))

                  sbuf(:,ireq) = datashadow%blk(xblk,yblk)%val(:,xloc,yloc)
               ENDDO

               idest = isrc
               CALL mpi_send (sbuf, nreq*nshadow, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
               
               deallocate (ylist)
               deallocate (xlist)
               deallocate (sbuf )

            ELSE
               worker_done(p_itis_worker(isrc)) = .true.
            ENDIF

         ENDDO

         deallocate (worker_done)

      ENDIF 
#endif

      IF (p_is_worker) THEN

         nsetshadow = 0

         allocate (pctshadow(nshadow,pixelset%nset))

         DO ipset = 1, pixelset%nset
            IF (any(filter(:) == pixelset%ltyp(ipset))) THEN

               iu   = pixelset%ibasin(ipset)
               ipxstt = pixelset%ipxstt(ipset)
               ipxend = pixelset%ipxend(ipset)
      
               allocate (datashadow1d (nshadow, ipxstt:ipxend))

#ifdef USEMPI
               allocate (xlist (ipxstt:ipxend))
               allocate (ylist (ipxstt:ipxend))
               allocate (ipt   (ipxstt:ipxend))
               allocate (msk   (ipxstt:ipxend))

               DO ipxl = ipxstt, ipxend
                  xlist(ipxl) = gshadow%xgrd(landbasin(iu)%ilon(ipxl))
                  ylist(ipxl) = gshadow%ygrd(landbasin(iu)%ilat(ipxl))

                  xblk = gshadow%xblk(xlist(ipxl))
                  yblk = gshadow%yblk(ylist(ipxl))
                  ipt(ipxl) = gblock%pio(xblk,yblk)
               ENDDO

               DO iproc = 0, p_np_io-1
                  msk = (ipt == p_address_io(iproc))
                  nreq = count(msk)

                  IF (nreq > 0) THEN

                     smesg = (/p_iam_glb, nreq/)
                     idest = p_address_io(iproc)
                     CALL mpi_send (smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

                     allocate (ibuf (nreq))
                     allocate (rbuf (nshadow,nreq))

                     ibuf = pack(xlist, msk)
                     CALL mpi_send (ibuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

                     ibuf = pack(ylist, msk)
                     CALL mpi_send (ibuf, nreq, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

                     isrc = idest

                     CALL mpi_recv (rbuf, nreq*nshadow, MPI_REAL8, &
                        isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
                     CALL unpack_inplace (rbuf, msk, datashadow1d)

                     deallocate (ibuf)
                     deallocate (rbuf)
                  ENDIF
               ENDDO

               deallocate (xlist)
               deallocate (ylist)
               deallocate (ipt  )
               deallocate (msk  )

#else
               DO ipxl = ipxstt, ipxend
                  ilon = gshadow%xgrd(landbasin(iu)%ilon(ipxl))
                  ilat = gshadow%ygrd(landbasin(iu)%ilat(ipxl))
                  xblk = gshadow%xblk(ilon)
                  yblk = gshadow%yblk(ilat)
                  xloc = gshadow%xloc(ilon)
                  yloc = gshadow%yloc(ilat)

                  datashadow1d(:,ipxl) = datashadow%blk(xblk,yblk)%val(:,xloc,yloc)
               ENDDO
#endif

               allocate (areapixel(ipxstt:ipxend))
               DO ipxl = ipxstt, ipxend
                  areapixel(ipxl) = areaquad (&
                     pixel%lat_s(landbasin(iu)%ilat(ipxl)), pixel%lat_n(landbasin(iu)%ilat(ipxl)), &
                     pixel%lon_w(landbasin(iu)%ilon(ipxl)), pixel%lon_e(landbasin(iu)%ilon(ipxl)) )
               ENDDO

               areatotal = sum(areapixel)

               DO ishadow = 1, nshadow
                  pctshadow(ishadow,ipset) = sum(datashadow1d(ishadow,:) * areapixel) / areatotal 
               ENDDO

               IF (any(pctshadow(:,ipset) > 0.)) THEN
                  nsetshadow = nsetshadow + count(pctshadow(:,ipset) > 0.)
                  pctshadow(:,ipset) = pctshadow(:,ipset) / sum(pctshadow(:,ipset))
               ENDIF

               deallocate (areapixel   )
               deallocate (datashadow1d)

            ELSE
               nsetshadow = nsetshadow + 1
            ENDIF
             
         ENDDO
         
#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            smesg = (/p_iam_glb, -1/)
            idest = p_address_io(iproc)
            CALL mpi_send (smesg, 2, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
         ENDDO
#endif

         allocate (bindex1(pixelset%nset))
         allocate (ibasin1(pixelset%nset))
         allocate (ipxstt1(pixelset%nset))
         allocate (ipxend1(pixelset%nset))
         allocate (ltyp1(pixelset%nset))

         bindex1 = pixelset%bindex
         ibasin1 = pixelset%ibasin
         ipxstt1 = pixelset%ipxstt
         ipxend1 = pixelset%ipxend
         ltyp1 = pixelset%ltyp

         deallocate (pixelset%bindex)
         deallocate (pixelset%ibasin)
         deallocate (pixelset%ipxstt)
         deallocate (pixelset%ipxend)
         deallocate (pixelset%ltyp)

         allocate (pixelset%bindex(nsetshadow))
         allocate (pixelset%ibasin(nsetshadow))
         allocate (pixelset%ipxstt(nsetshadow))
         allocate (pixelset%ipxend(nsetshadow))
         allocate (pixelset%ltyp(nsetshadow))
         
         allocate (fracout    (nsetshadow))
         allocate (shadowclass(nsetshadow))

         jpset = 0
         DO ipset = 1, pixelset%nset
            IF (any(filter(:) == ltyp1(ipset))) THEN
               IF (any(pctshadow(:,ipset) > 0.)) THEN
                  DO ishadow = 1, nshadow
                     IF (pctshadow(ishadow,ipset) > 0.) THEN
                        jpset = jpset + 1
                        pixelset%bindex(jpset) = bindex1(ipset)
                        pixelset%ibasin(jpset) = ibasin1(ipset)
                        pixelset%ipxstt(jpset) = ipxstt1(ipset)
                        pixelset%ipxend(jpset) = ipxend1(ipset)
                        pixelset%ltyp(jpset) = ltyp1(ipset)

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
               pixelset%bindex(jpset) = bindex1(ipset)
               pixelset%ibasin(jpset) = ibasin1(ipset)
               pixelset%ipxstt(jpset) = ipxstt1(ipset)
               pixelset%ipxend(jpset) = ipxend1(ipset)
               pixelset%ltyp(jpset) = ltyp1(ipset)
                     
               IF (present(fracin)) THEN
                  fracout(jpset) = fracin(ipset)
               ELSE
                  fracout(jpset) = 1.
               ENDIF

               shadowclass(jpset) = 0 ! no meaning
            ENDIF
         ENDDO

         pixelset%nset = nsetshadow

         deallocate (bindex1)
         deallocate (ibasin1)
         deallocate (ipxstt1)
         deallocate (ipxend1)
         deallocate (ltyp1)
         deallocate (pctshadow)

      ENDIF
         
      CALL pixelset%set_vecgs

   END SUBROUTINE pixelsetshadow_build

END MODULE mod_pixelsetshadow
