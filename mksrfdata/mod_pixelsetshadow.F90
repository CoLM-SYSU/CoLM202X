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
      USE mod_landunit
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
      INTEGER  :: ipxl, iu, istt, iend, ishadow
      INTEGER,  allocatable :: unum1(:), iunt1(:), istt1(:), iend1(:), ltyp1(:)
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

               iu   = pixelset%iunt(ipset)
               istt = pixelset%istt(ipset)
               iend = pixelset%iend(ipset)
      
               allocate (datashadow1d (nshadow, istt:iend))

#ifdef USEMPI
               allocate (xlist (istt:iend))
               allocate (ylist (istt:iend))
               allocate (ipt   (istt:iend))
               allocate (msk   (istt:iend))

               DO ipxl = istt, iend
                  xlist(ipxl) = gshadow%xgrd(landunit(iu)%ilon(ipxl))
                  ylist(ipxl) = gshadow%ygrd(landunit(iu)%ilat(ipxl))

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
               DO ipxl = istt, iend
                  ilon = gshadow%xgrd(landunit(iu)%ilon(ipxl))
                  ilat = gshadow%ygrd(landunit(iu)%ilat(ipxl))
                  xblk = gshadow%xblk(ilon)
                  yblk = gshadow%yblk(ilat)
                  xloc = gshadow%xloc(ilon)
                  yloc = gshadow%yloc(ilat)

                  datashadow1d(:,ipxl) = datashadow%blk(xblk,yblk)%val(:,xloc,yloc)
               ENDDO
#endif

               allocate (areapixel(istt:iend))
               DO ipxl = istt, iend
                  areapixel(ipxl) = areaquad (&
                     pixel%lat_s(landunit(iu)%ilat(ipxl)), pixel%lat_n(landunit(iu)%ilat(ipxl)), &
                     pixel%lon_w(landunit(iu)%ilon(ipxl)), pixel%lon_e(landunit(iu)%ilon(ipxl)) )
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

         allocate (unum1(pixelset%nset))
         allocate (iunt1(pixelset%nset))
         allocate (istt1(pixelset%nset))
         allocate (iend1(pixelset%nset))
         allocate (ltyp1(pixelset%nset))

         unum1 = pixelset%unum
         iunt1 = pixelset%iunt
         istt1 = pixelset%istt
         iend1 = pixelset%iend
         ltyp1 = pixelset%ltyp

         deallocate (pixelset%unum)
         deallocate (pixelset%iunt)
         deallocate (pixelset%istt)
         deallocate (pixelset%iend)
         deallocate (pixelset%ltyp)

         allocate (pixelset%unum(nsetshadow))
         allocate (pixelset%iunt(nsetshadow))
         allocate (pixelset%istt(nsetshadow))
         allocate (pixelset%iend(nsetshadow))
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
                        pixelset%unum(jpset) = unum1(ipset)
                        pixelset%iunt(jpset) = iunt1(ipset)
                        pixelset%istt(jpset) = istt1(ipset)
                        pixelset%iend(jpset) = iend1(ipset)
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
               pixelset%unum(jpset) = unum1(ipset)
               pixelset%iunt(jpset) = iunt1(ipset)
               pixelset%istt(jpset) = istt1(ipset)
               pixelset%iend(jpset) = iend1(ipset)
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

         deallocate (unum1)
         deallocate (iunt1)
         deallocate (istt1)
         deallocate (iend1)
         deallocate (ltyp1)
         deallocate (pctshadow)

      ENDIF
         
      CALL pixelset%set_vecgs

   END SUBROUTINE pixelsetshadow_build

END MODULE mod_pixelsetshadow
