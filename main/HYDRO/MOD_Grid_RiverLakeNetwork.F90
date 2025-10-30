#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeNetwork
!--------------------------------------------------------------------------------
! DESCRIPTION:
!--------------------------------------------------------------------------------

   USE MOD_WorkerPushData
   IMPLICIT NONE

   ! ----- River Lake network -----

   integer :: totalnumucat
   integer :: numucat

   integer, allocatable :: numucat_wrk (:)
   type(pointer_int32_1d), allocatable :: ucat_data_address (:)

   ! ----- Part 1: between elements and unit catchments -----
   integer, allocatable :: ucat_elid (:)   ! index in element numbering
   integer, allocatable :: ucat_ucid (:)   ! index in unit catchment numbering

   type(worker_pushdata_type) :: push_ucat2elm

   ! ----- Part 2: between runoff input elements and unit catchments -----
   integer :: inpn
   integer,  allocatable :: inpmat_el2uc    (:,:)
   real(r8), allocatable :: inpmat_area_e2u (:,:)

   integer :: nucpart
   integer,  allocatable :: inpmat_uc2el    (:,:)
   real(r8), allocatable :: inpmat_area_u2e (:,:)

   type(worker_pushdata_type) :: push_inpmat2ucat
   type(worker_pushdata_type) :: push_ucat2inpmat

   ! ----- Part 3: between upstream and downstream unit catchments -----
   integer, allocatable :: ucat_next (:)   ! next unit catchment
   integer :: upnmax
   integer, allocatable :: ucat_ups  (:,:) ! upstream unit catchments

   type(worker_pushdata_type) :: push_next2ucat
   type(worker_pushdata_type) :: push_ups2ucat

   ! ----- Part 4: river systems -----
   integer :: numrivsys
   logical :: rivsys_by_multiple_procs
   integer, allocatable :: irivsys (:)
#ifdef USEMPI
   integer :: p_comm_rivsys
#endif


   ! ----- Parameters for River and Lake -----

   integer,  allocatable :: lake_type      (:)   ! 0: river; 2: reservoir.

   real(r8), allocatable :: topo_rivelv    (:)   ! river bed elevation [m]
   real(r8), allocatable :: topo_rivhgt    (:)   ! river channel depth [m]
   real(r8), allocatable :: topo_rivlen    (:)   ! river channel length [m]
   real(r8), allocatable :: topo_rivman    (:)   ! river manning coefficient [m]
   real(r8), allocatable :: topo_rivwth    (:)   ! river channel width [m]
   real(r8), allocatable :: topo_rivare    (:)   ! river channel area [m^2]
   real(r8), allocatable :: topo_rivstomax (:)   ! max river channel storage [m^3]

   real(r8), allocatable :: topo_area      (:)   ! floodplain area [m^2]
   real(r8), allocatable :: topo_fldhgt    (:,:) ! floodplain height profile [m]

   real(r8), allocatable :: bedelv_next    (:)   ! downstream river bed elevation [m]
   real(r8), allocatable :: outletwth      (:)   ! river outlet width [m]

   type :: vol_dep_curve_type
      integer  :: nlfp
      real(r8) :: rivhgt
      real(r8) :: rivare
      real(r8) :: rivstomax
      real(r8), allocatable :: flphgt    (:) ! floodplain height profile [m]
      real(r8), allocatable :: flparea   (:) ! flood plain area [m^2]
      real(r8), allocatable :: flpaccare (:) ! flood plain accumulated area [m^2]
      real(r8), allocatable :: flpstomax (:) ! max flood plain storage [m^3]
   CONTAINS
      procedure, PUBLIC :: depth     => retrieve_depth_from_volume
      procedure, PUBLIC :: volume    => retrieve_volume_from_depth
      procedure, PUBLIC :: floodarea => retrieve_area_from_depth
      final :: vol_depth_curve_free_mem
   END type vol_dep_curve_type

   type(vol_dep_curve_type), allocatable :: floodplain_curve (:)

CONTAINS

   ! ----------
   SUBROUTINE build_riverlake_network ()

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Mesh
   USE MOD_LandElm
   USE MOD_Utils
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   ! Local Variables
   character(len=256)    :: parafile
   integer,  allocatable :: seq_x(:), seq_y(:)
   real(r8), allocatable :: x(:)

   integer,  allocatable :: inpmat_x(:,:), inpmat_y(:,:)

   integer :: numrivmth
   integer,  allocatable :: rivermouth(:)

   integer,  allocatable :: nups_nst  (:), iups_nst  (:), nups_all(:)
   integer,  allocatable :: uc_up2down(:), order_ucat(:)
   integer,  allocatable :: addr_ucat (:)

   integer , allocatable :: nuc_rs(:), iwrk_rs(:), nwrk_rs(:), nave_rs(:)
   real(r8), allocatable :: wt_uc (:), wt_rs  (:), wt_wrk (:), nuc_wrk(:)

   integer,  allocatable :: elmindex(:)

   integer,  allocatable :: idata1d(:), idata2d(:,:)
   real(r8), allocatable :: rdata2d(:,:)

   integer,  allocatable :: allelm_in_inp (:), nucat_elm(:), iucat_elm(:)

   integer,  allocatable :: inpmat_uc2el_all    (:,:)
   real(r8), allocatable :: inpmat_area_u2e_all (:,:)

   integer  :: nucat, iriv, nelmall, nelm
   integer  :: p_np_rivsys, color
   integer  :: iworker, iwrkdsp
   integer  :: iloc, i, j, ithis
   real(r8) :: sumwt
   logical  :: is_new


#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ! read in parameters from file.
      IF (p_is_master) THEN

         parafile = DEF_UnitCatchment_file

         CALL ncio_read_serial (parafile, 'seq_next', ucat_next)

         CALL ncio_read_serial (parafile, 'seq_x', seq_x)
         CALL ncio_read_serial (parafile, 'seq_y', seq_y)
         CALL ncio_read_serial (parafile, 'lon', x)

         CALL ncio_read_serial (parafile, 'inpmat_x', inpmat_x)
         CALL ncio_read_serial (parafile, 'inpmat_y', inpmat_y)
         CALL ncio_read_serial (parafile, 'inpmat_area', inpmat_area_e2u)

      ENDIF

      IF (p_is_master) THEN

         totalnumucat = size(ucat_next)

         allocate (nups_nst (totalnumucat))
         allocate (iups_nst (totalnumucat))

         nups_nst(:) = 0
         DO i = 1, totalnumucat
            j = ucat_next(i)
            IF (j > 0) THEN
               nups_nst(j) = nups_nst(j) + 1
            ENDIF
         ENDDO

#ifdef USEMPI
         ! divide unit catchments into groups and assign to workers

         allocate (wt_uc (totalnumucat));  wt_uc(:) = 1.

         ! sort unit catchment from upstream to downstream, recorded by "uc_up2down"

         allocate (uc_up2down (totalnumucat))

         ithis = 0
         iups_nst(:) = 0
         DO i = 1, totalnumucat
            IF (iups_nst(i) == nups_nst(i)) THEN

               ithis = ithis + 1
               uc_up2down(ithis) = i
               iups_nst(i) = -1

               j = ucat_next(i)
               DO WHILE (j > 0)

                  iups_nst(j) = iups_nst(j) + 1

                  IF (iups_nst(j) == nups_nst(j)) THEN
                     ithis = ithis + 1
                     uc_up2down(ithis) = j
                     iups_nst(j) = -1

                     j = ucat_next(j)
                  ELSE
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         allocate (rivermouth (totalnumucat))
         numrivmth = 0
         DO i = totalnumucat, 1, -1
            j = ucat_next(uc_up2down(i))
            IF (j <= 0) THEN
               numrivmth = numrivmth + 1
               rivermouth(uc_up2down(i)) = numrivmth
            ELSE
               rivermouth(uc_up2down(i)) = rivermouth(j)
            ENDIF
         ENDDO

         allocate (nuc_rs (numrivmth)); nuc_rs(:) = 0
         allocate (wt_rs  (numrivmth)); wt_rs (:) = 0.
         DO i = 1, totalnumucat
            nuc_rs(rivermouth(i)) = nuc_rs(rivermouth(i)) + 1
            wt_rs (rivermouth(i)) = wt_rs (rivermouth(i)) + wt_uc(i)
         ENDDO

         sumwt = sum(wt_rs)

         allocate (iwrk_rs (numrivmth))
         allocate (nwrk_rs (numrivmth))
         allocate (nave_rs (numrivmth))

         iwrkdsp = -1
         DO i = 1, numrivmth
            nwrk_rs(i) = floor(wt_rs(i)/sumwt * p_np_worker)
            IF (nwrk_rs(i) > 1) THEN

               nave_rs(i) = nuc_rs(i) / nwrk_rs(i)
               IF (mod(nuc_rs(i), nwrk_rs(i)) /= 0) THEN
                  nave_rs(i) = nave_rs(i) + 1
               ENDIF

               iwrk_rs(i) = iwrkdsp + 1
               iwrkdsp = iwrkdsp + nwrk_rs(i)
            ENDIF
         ENDDO

         allocate (nups_all (totalnumucat));  nups_all(:) = 1

         DO i = 1, totalnumucat
            j = ucat_next(uc_up2down(i))
            IF (j > 0) THEN
               nups_all(j) = nups_all(j) + nups_all(uc_up2down(i))
            ENDIF
         ENDDO

         allocate (addr_ucat (totalnumucat));  addr_ucat(:) = -1

         allocate (wt_wrk (0:p_np_worker-1));  wt_wrk (:) = 0
         allocate (nuc_wrk(0:p_np_worker-1));  nuc_wrk(:) = 0

         allocate (order_ucat (totalnumucat))
         order_ucat(uc_up2down) = (/(i, i = 1, totalnumucat)/)

         ithis = totalnumucat
         DO WHILE (ithis > 0)

            i = uc_up2down(ithis)

            IF (addr_ucat(i) >= 0) THEN
               ithis = ithis - 1
               CYCLE
            ENDIF

            j = ucat_next(i)
            IF (j > 0) THEN
               IF (addr_ucat(j) >= 0) THEN
                  addr_ucat(i) = addr_ucat(j)
                  ithis = ithis - 1
                  CYCLE
               ENDIF
            ENDIF

            iriv = rivermouth(i)
            IF (nwrk_rs(iriv) > 1) THEN
               iworker = iwrk_rs(iriv)
               IF (nups_all(i) <= nave_rs(iriv)-nuc_wrk(iworker)) THEN

                  addr_ucat(i) = p_address_worker(iworker)

                  nuc_wrk(iworker) = nuc_wrk(iworker) + nups_all(i)
                  IF (nuc_wrk(iworker) == nave_rs(iriv)) THEN
                     iwrk_rs(iriv) = iwrk_rs(iriv) + 1
                  ENDIF

                  j = ucat_next(i)
                  IF (j > 0) THEN
                     DO WHILE (j > 0)
                        nups_all(j) = nups_all(j) - nups_all(i)
                        ithis = order_ucat(j)
                        j = ucat_next(j)
                     ENDDO
                  ELSE
                     ithis = ithis - 1
                  ENDIF
               ELSE
                  ithis = ithis - 1
               ENDIF
            ELSE
               iworker = minloc(wt_wrk(iwrkdsp+1:p_np_worker-1), dim=1) + iwrkdsp

               addr_ucat(i) = p_address_worker(iworker)

               wt_wrk(iworker) = wt_wrk(iworker) + wt_rs(iriv)
               ithis = ithis - 1
            ENDIF

         ENDDO

         deallocate (uc_up2down)
         deallocate (order_ucat)
         deallocate (nups_all  )
         deallocate (nuc_rs    )
         deallocate (iwrk_rs   )
         deallocate (nwrk_rs   )
         deallocate (nave_rs   )
         deallocate (wt_uc     )
         deallocate (wt_rs     )
         deallocate (wt_wrk    )
         deallocate (nuc_wrk   )

      ENDIF

      IF (p_is_master) THEN

         allocate(ucat_ucid (totalnumucat))
         ucat_ucid = (/(i, i = 1, totalnumucat)/)

         allocate (numucat_wrk       (0:p_np_worker-1))
         allocate (ucat_data_address (0:p_np_worker-1))

         DO iworker = 0, p_np_worker-1
            nucat = count(addr_ucat == p_address_worker(iworker))
            numucat_wrk(iworker) = nucat
            IF (nucat > 0) THEN
               allocate (ucat_data_address(iworker)%val (nucat))
               ucat_data_address(iworker)%val = &
                  pack(ucat_ucid, mask = (addr_ucat == p_address_worker(iworker)))
            ENDIF
         ENDDO

         deallocate (ucat_ucid)
      ENDIF

      CALL mpi_bcast (totalnumucat, 1, mpi_integer, p_address_master, p_comm_glb, p_err)

      ! send unit catchment index to workers
      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1

            CALL mpi_send (numucat_wrk(iworker), 1, MPI_INTEGER, p_address_worker(iworker), &
               mpi_tag_mesg, p_comm_glb, p_err)

            nucat = numucat_wrk(iworker)
            IF (nucat > 0) THEN
               CALL mpi_send (ucat_data_address(iworker)%val, nucat, MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO

      ELSEIF (p_is_worker) THEN

         CALL mpi_recv (numucat, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

         IF (numucat > 0) THEN
            allocate (ucat_ucid (numucat))
            CALL mpi_recv (ucat_ucid, numucat, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      numucat = totalnumcat

      allocate(ucat_ucid (totalnumucat))
      ucat_ucid = (/(i, i = 1, totalnumcat)/)

      allocate (numucat_wrk (0:0))
      numucat_wrk(0) = numucat

      allocate (ucat_data_address (0:0))
      allocate (ucat_data_address(0)%val (nucat))
      ucat_data_address(0)%val = ucat_ucid
#endif

      IF (allocated(addr_ucat)) deallocate(addr_ucat)

      ! ----- Part 1: between elements and unit catchments -----

      IF (p_is_master) THEN
         allocate(ucat_elid (totalnumucat))
         ucat_elid = (seq_y-1)*size(x) + seq_x
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      ! send unit catchment index to workers
      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1
            nucat = numucat_wrk(iworker)
            IF (nucat > 0) THEN
               allocate (idata1d (nucat))

               idata1d = ucat_elid(ucat_data_address(iworker)%val)
               CALL mpi_send (idata1d, nucat, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)

               deallocate (idata1d)
            ENDIF
         ENDDO

         deallocate (ucat_elid)

      ELSEIF (p_is_worker) THEN

         IF (numucat > 0) THEN
            allocate (ucat_elid (numucat))
            CALL mpi_recv (ucat_elid, numucat, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_worker) THEN
         IF (numelm > 0) THEN
            allocate (elmindex (numelm))
            elmindex = landelm%eindex
         ENDIF
      ENDIF

      CALL build_worker_pushdata (numucat, ucat_elid, numelm,  elmindex,  push_ucat2elm)


      IF (allocated(seq_x)) deallocate(seq_x)
      IF (allocated(seq_y)) deallocate(seq_y)

      ! ----- Part 2: between runoff input elements and unit catchments -----

      IF (p_is_master) THEN

         inpn = size(inpmat_area_e2u,1)

         allocate(inpmat_el2uc (inpn,totalnumucat))

         inpmat_el2uc = (inpmat_y-1)*size(x) + inpmat_x

         WHERE ((inpmat_area_e2u <= 0) .or. (inpmat_el2uc <= 0))
            inpmat_el2uc    = 0
            inpmat_area_e2u = 0.
         END WHERE

         allocate (allelm_in_inp (inpn*totalnumucat))
         allocate (nucat_elm     (inpn*totalnumucat))

         nelmall = 0
         nucat_elm(:) = 0
         DO i = 1, totalnumucat
            DO j = 1, inpn
               IF (inpmat_el2uc(j,i) > 0) THEN

                  CALL insert_into_sorted_list1 (inpmat_el2uc(j,i), nelmall, allelm_in_inp, iloc, is_new)

                  IF (is_new) THEN
                     IF (iloc < nelmall) THEN
                        nucat_elm(iloc+1:nelmall) = nucat_elm(iloc:nelmall-1)
                     ENDIF
                     nucat_elm(iloc) = 1
                  ELSE
                     nucat_elm(iloc) = nucat_elm(iloc) + 1
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

         nucpart = maxval(nucat_elm (1:nelmall))

         allocate (inpmat_uc2el_all    (nucpart, nelmall));  inpmat_uc2el_all   (:,:) = 0
         allocate (inpmat_area_u2e_all (nucpart, nelmall));  inpmat_area_u2e_all(:,:) = 0.

         allocate (iucat_elm (nelmall)); iucat_elm(:) = 0

         DO i = 1, totalnumucat
            DO j = 1, inpn
               IF (inpmat_el2uc(j,i) > 0) THEN
                  iloc = find_in_sorted_list1 (inpmat_el2uc(j,i), nelmall, allelm_in_inp(1:nelmall))
                  iucat_elm(iloc) = iucat_elm(iloc) + 1
                  inpmat_uc2el_all   (iucat_elm(iloc),iloc) = i
                  inpmat_area_u2e_all(iucat_elm(iloc),iloc) = inpmat_area_e2u(j,i)
               ENDIF
            ENDDO
         ENDDO

      ENDIF

#ifdef USEMPI
      CALL mpi_bcast (inpn, 1, mpi_integer, p_address_master, p_comm_glb, p_err)

      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1

            nucat = numucat_wrk(iworker)

            IF (nucat > 0) THEN
               allocate (idata2d (inpn, nucat))
               DO i = 1, nucat
                  idata2d(:,i) = inpmat_el2uc(:,ucat_data_address(iworker)%val(i))
               ENDDO

               CALL mpi_send (idata2d, inpn*nucat, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)

               allocate (rdata2d (inpn, nucat))
               DO i = 1, nucat
                  rdata2d(:,i) = inpmat_area_e2u(:,ucat_data_address(iworker)%val(i))
               ENDDO

               CALL mpi_send (rdata2d, inpn*nucat, MPI_REAL8, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)

               deallocate (idata2d)
               deallocate (rdata2d)
            ENDIF
         ENDDO

         deallocate (inpmat_el2uc   )
         deallocate (inpmat_area_e2u)

      ELSEIF (p_is_worker) THEN

         IF (numucat > 0) THEN

            allocate (inpmat_el2uc (inpn, numucat))
            CALL mpi_recv (inpmat_el2uc, inpn*numucat, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (inpmat_area_e2u (inpn, numucat))
            CALL mpi_recv (inpmat_area_e2u, inpn*numucat, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

         ENDIF

      ENDIF

      CALL mpi_bcast (nucpart, 1, mpi_integer, p_address_master, p_comm_glb, p_err)

      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (nelm, 1, MPI_INTEGER, &
               p_address_worker(iworker), mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            IF (nelm > 0) THEN

               allocate (elmindex (nelm))
               allocate (idata2d  (nucpart, nelm));  idata2d(:,:) = 0
               allocate (rdata2d  (nucpart, nelm));  rdata2d(:,:) = 0.

               CALL mpi_recv (elmindex, nelm, MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO i = 1, nelm
                  iloc = find_in_sorted_list1 (elmindex(i), nelmall, allelm_in_inp(1:nelmall))
                  IF (iloc > 0) THEN
                     idata2d(:,i) = inpmat_uc2el_all   (:,iloc)
                     rdata2d(:,i) = inpmat_area_u2e_all(:,iloc)
                  ENDIF
               ENDDO

               CALL mpi_send (idata2d, nucpart*nelm, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)
               CALL mpi_send (rdata2d, nucpart*nelm, MPI_REAL8,   p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)

               deallocate (elmindex)
               deallocate (idata2d )
               deallocate (rdata2d )
            ENDIF
         ENDDO

      ELSEIF (p_is_worker) THEN

         CALL mpi_send (numelm, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (numelm > 0) THEN

            CALL mpi_send (elmindex, numelm, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_err)

            allocate (inpmat_uc2el (nucpart, numelm))
            CALL mpi_recv (inpmat_uc2el, nucpart*numelm, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (inpmat_area_u2e (nucpart, numelm))
            CALL mpi_recv (inpmat_area_u2e, nucpart*numelm, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      allocate (inpmat_uc2el    (nucpart, numelm))
      allocate (inpmat_area_u2e (nucpart, numelm))
      inpmat_uc2el    = 0
      inpmat_area_u2e = 0.

      DO = 1, numelm
         iloc = find_in_sorted_list1 (elmindex(i), nelmall, allelm_in_inp(1:nelmall))
         IF (iloc > 0) THEN
            inpmat_uc2el    (:,i) = inpmat_uc2el_all   (:,iloc)
            inpmat_area_u2e (:,i) = inpmat_area_u2e_all(:,iloc)
         ENDIF
      ENDDO
#endif


      CALL build_worker_pushdata (numelm,  elmindex,  numucat, inpmat_el2uc, push_inpmat2ucat)
      CALL build_worker_pushdata (numucat, ucat_ucid, numelm,  inpmat_uc2el, push_ucat2inpmat)


      IF (p_is_master) THEN
         deallocate (inpmat_x           )
         deallocate (inpmat_y           )
         deallocate (allelm_in_inp      )
         deallocate (nucat_elm          )
         deallocate (iucat_elm          )
         deallocate (inpmat_uc2el_all   )
         deallocate (inpmat_area_u2e_all)
         deallocate (x                  )
      ENDIF

      IF (allocated(elmindex)) deallocate(elmindex)

      ! ----- Part 3: between upstream and downstream unit catchments -----

      IF (p_is_master) THEN

         upnmax = maxval(nups_nst)
         allocate (ucat_ups (upnmax,totalnumucat))
         ucat_ups(:,:) = 0

         iups_nst(:) = 0
         DO i = 1, totalnumucat
            j = ucat_next(i)
            IF (j > 0) THEN
               iups_nst(j) = iups_nst(j) + 1
               ucat_ups(iups_nst(j),j) = i
            ENDIF
         ENDDO

      ENDIF


#ifdef USEMPI
      CALL mpi_bcast (upnmax, 1, mpi_integer, p_address_master, p_comm_glb, p_err)

      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1

            nucat = numucat_wrk(iworker)

            IF (nucat > 0) THEN
               allocate (idata1d (nucat))
               idata1d = ucat_next(ucat_data_address(iworker)%val)
               CALL mpi_send (idata1d, nucat, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)

               allocate (idata2d (upnmax, nucat))
               DO i = 1, nucat
                  idata2d(:,i) = ucat_ups(:,ucat_data_address(iworker)%val(i))
               ENDDO
               CALL mpi_send (idata2d, upnmax*nucat, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)

               deallocate (idata1d)
               deallocate (idata2d)
            ENDIF
         ENDDO

         deallocate (ucat_ups )

      ELSEIF (p_is_worker) THEN

         IF (numucat > 0) THEN

            allocate (ucat_next (numucat))
            CALL mpi_recv (ucat_next, numucat, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (ucat_ups (upnmax, numucat))
            CALL mpi_recv (ucat_ups, upnmax*numucat, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#endif


      CALL build_worker_pushdata (numucat, ucat_ucid, numucat, ucat_next, push_next2ucat)
      CALL build_worker_pushdata (numucat, ucat_ucid, numucat, ucat_ups,  push_ups2ucat )

#ifdef CoLMDEBUG
      ! IF (p_is_worker) THEN
      !    write(*,'(A,I0,A,I0,A,I0,A)') 'worker ', p_iam_worker, ' has ', numucat, &
      !       ' unit catchment with ', sum(push_next2ucat%n_from_other), ' downstream to other workers'
      ! ENDIF
#endif

      IF (allocated(nups_nst  )) deallocate(nups_nst  )
      IF (allocated(iups_nst  )) deallocate(iups_nst  )

      ! ----- Part 4: river systems -----

#ifdef USEMPI
      IF (p_is_master) THEN
         DO iworker = 0, p_np_worker-1
            nucat = numucat_wrk(iworker)
            IF (nucat > 0) THEN
               allocate (idata1d (nucat))
               idata1d = rivermouth(ucat_data_address(iworker)%val)
               CALL mpi_send (idata1d, nucat, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)
               deallocate (idata1d)
            ENDIF
         ENDDO
      ELSEIF (p_is_worker) THEN
         IF (numucat > 0) THEN
            allocate (rivermouth (numucat))
            CALL mpi_recv (rivermouth, numucat, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
      ENDIF

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            color = maxval(rivermouth)
            CALL mpi_comm_split (p_comm_worker, color, p_iam_worker, p_comm_rivsys, p_err)
         ELSE
            CALL mpi_comm_split (p_comm_worker, MPI_UNDEFINED, p_iam_worker, p_comm_rivsys, p_err)
         ENDIF

         rivsys_by_multiple_procs = .false.
         IF (p_comm_rivsys /= MPI_COMM_NULL) THEN
            CALL mpi_comm_size (p_comm_rivsys, p_np_rivsys, p_err)
            IF (p_np_rivsys > 1) THEN
               rivsys_by_multiple_procs = .true.
            ENDIF
         ENDIF
      ENDIF
#else
      rivsys_by_multiple_procs = .false.
#endif

      IF (p_is_worker) THEN

         IF (numucat > 0) allocate (irivsys (numucat))

         IF (.not. rivsys_by_multiple_procs) THEN
            IF (numucat > 0) THEN

               allocate (order_ucat (numucat))
               order_ucat = (/(i, i = 1, numucat)/)

               CALL quicksort (numucat, rivermouth, order_ucat)

               numrivsys = 1
               irivsys(order_ucat(1)) = numrivsys
               DO i = 2, numucat
                  IF (rivermouth(i) /= rivermouth(i-1)) THEN
                     numrivsys = numrivsys + 1
                  ENDIF
                  irivsys(order_ucat(i)) = numrivsys
               ENDDO

            ENDIF
         ELSE
            numrivsys  = 1
            irivsys(:) = 1
         ENDIF

      ENDIF


      IF (allocated(rivermouth)) deallocate(rivermouth)
      IF (allocated(order_ucat)) deallocate(order_ucat)

      ! ----- Parameters for River and Lake -----

      CALL readin_riverlake_parameter (parafile, 'topo_rivelv',    rdata1d = topo_rivelv   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivhgt',    rdata1d = topo_rivhgt   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivlen',    rdata1d = topo_rivlen   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivman',    rdata1d = topo_rivman   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivwth',    rdata1d = topo_rivwth   )
      CALL readin_riverlake_parameter (parafile, 'topo_rivstomax', rdata1d = topo_rivstomax)
      CALL readin_riverlake_parameter (parafile, 'topo_area',      rdata1d = topo_area     )
      CALL readin_riverlake_parameter (parafile, 'topo_fldhgt',    rdata2d = topo_fldhgt   )

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN

            allocate (lake_type (numucat))
            lake_type(:) = 0

            allocate (topo_rivare (numucat))
            topo_rivare = topo_rivstomax / topo_rivhgt

            allocate (floodplain_curve (numucat))

            DO i = 1, numucat
               floodplain_curve(i)%nlfp      = size(topo_fldhgt,1)
               floodplain_curve(i)%rivhgt    = topo_rivhgt(i)
               floodplain_curve(i)%rivstomax = topo_rivstomax(i)
               floodplain_curve(i)%rivare    = topo_rivare(i)

               allocate (floodplain_curve(i)%flphgt    (0:floodplain_curve(i)%nlfp))
               allocate (floodplain_curve(i)%flparea   (0:floodplain_curve(i)%nlfp))
               allocate (floodplain_curve(i)%flpaccare (0:floodplain_curve(i)%nlfp))
               allocate (floodplain_curve(i)%flpstomax (0:floodplain_curve(i)%nlfp))

               floodplain_curve(i)%flphgt(0)  = 0.
               floodplain_curve(i)%flphgt(1:) = topo_fldhgt(:,i)

               floodplain_curve(i)%flparea(0)  = 0.
               floodplain_curve(i)%flparea(1:) = topo_area(i) / floodplain_curve(i)%nlfp

               floodplain_curve(i)%flpaccare(0) = 0.
               DO j = 1, floodplain_curve(i)%nlfp
                  floodplain_curve(i)%flpaccare(j) = &
                     floodplain_curve(i)%flpaccare(j-1) + floodplain_curve(i)%flparea(j)
               ENDDO

               floodplain_curve(i)%flpstomax(0) = 0.
               DO j = 1, floodplain_curve(i)%nlfp
                  floodplain_curve(i)%flpstomax(j) = floodplain_curve(i)%flpstomax(j-1)        &
                     + 0.5 * (floodplain_curve(i)%flparea(j) + floodplain_curve(i)%flparea(j-1)) &
                           * (floodplain_curve(i)%flphgt(j)  - floodplain_curve(i)%flphgt(j-1))
               ENDDO
            ENDDO

            allocate (bedelv_next (numucat))
            allocate (outletwth   (numucat))

         ENDIF
      ENDIF

      CALL worker_push_data (push_next2ucat, topo_rivelv, bedelv_next, fillvalue = spval)
      CALL worker_push_data (push_next2ucat, topo_rivwth, outletwth  , fillvalue = spval)

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            WHERE (ucat_next > 0)
               outletwth = (outletwth + topo_rivwth) * 0.5
            ELSEWHERE
               outletwth = topo_rivwth
            END WHERE
         ENDIF
      ENDIF

   END SUBROUTINE build_riverlake_network

   ! ---------
   SUBROUTINE readin_riverlake_parameter (parafile, varname, rdata1d, rdata2d, idata1d)

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   character(len=*), intent(in) :: parafile
   character(len=*), intent(in) :: varname

   real(r8), allocatable, intent(inout), optional :: rdata1d (:)
   real(r8), allocatable, intent(inout), optional :: rdata2d (:,:)
   integer,  allocatable, intent(inout), optional :: idata1d (:)

   ! Local Variables
   integer :: iworker, nucat, ndim1, i
   real(r8), allocatable :: rsend1d (:)
   real(r8), allocatable :: rsend2d (:,:)
   integer,  allocatable :: isend1d (:)

      IF (p_is_master) THEN
         IF (present(rdata1d))  CALL ncio_read_serial (parafile, varname, rdata1d)
         IF (present(rdata2d))  CALL ncio_read_serial (parafile, varname, rdata2d)
         IF (present(idata1d))  CALL ncio_read_serial (parafile, varname, idata1d)
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (present(rdata2d)) THEN
         IF (p_is_master) ndim1 = size(rdata2d,1)
         CALL mpi_bcast (ndim1, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      ENDIF

      ! send unit catchment index to workers
      IF (p_is_master) THEN

         DO iworker = 0, p_np_worker-1

            nucat = numucat_wrk(iworker)

            IF (nucat > 0) THEN
               IF (present(rdata1d)) THEN
                  allocate (rsend1d (nucat))

                  rsend1d = rdata1d(ucat_data_address(iworker)%val)
                  CALL mpi_send (rsend1d, nucat, MPI_REAL8, p_address_worker(iworker), &
                     mpi_tag_data, p_comm_glb, p_err)

                  deallocate (rsend1d)
               ENDIF

               IF (present(rdata2d)) THEN
                  allocate (rsend2d (ndim1,nucat))

                  DO i = 1, nucat
                     rsend2d(:,i) = rdata2d(:,ucat_data_address(iworker)%val(i))
                  ENDDO
                  CALL mpi_send (rsend2d, ndim1*nucat, MPI_REAL8, p_address_worker(iworker), &
                     mpi_tag_data, p_comm_glb, p_err)

                  deallocate (rsend2d)
               ENDIF

               IF (present(idata1d)) THEN
                  allocate (isend1d (nucat))

                  isend1d = idata1d(ucat_data_address(iworker)%val)
                  CALL mpi_send (isend1d, nucat, MPI_INTEGER, p_address_worker(iworker), &
                     mpi_tag_data, p_comm_glb, p_err)

                  deallocate (rsend1d)
               ENDIF
            ENDIF
         ENDDO

         IF (present(rdata1d))  deallocate (rdata1d)
         IF (present(rdata2d))  deallocate (rdata2d)
         IF (present(idata1d))  deallocate (idata1d)

      ELSEIF (p_is_worker) THEN

         IF (numucat > 0) THEN
            IF (present(rdata1d)) THEN
               allocate (rdata1d (numucat))
               CALL mpi_recv (rdata1d, numucat, MPI_REAL8, p_address_master, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF

            IF (present(rdata2d)) THEN
               allocate (rdata2d (ndim1,numucat))
               CALL mpi_recv (rdata2d, ndim1*numucat, MPI_REAL8, p_address_master, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF

            IF (present(idata1d)) THEN
               allocate (idata1d (numucat))
               CALL mpi_recv (idata1d, numucat, MPI_INTEGER, p_address_master, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF
         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE readin_riverlake_parameter

   !
   FUNCTION retrieve_depth_from_volume (this, volume) result(depth)

   IMPLICIT NONE

   class(vol_dep_curve_type) :: this
   real(r8), intent(in)      :: volume
   real(r8) :: depth

   ! Local Variables
   real(r8) :: v0, g
   integer  :: i

      v0 = volume - this%rivstomax
      IF (v0 <= 0) THEN
         depth = volume / this%rivare
      ELSE
         i = 1
         DO WHILE (i <= this%nlfp)
            IF (v0 > this%flpstomax(i)) THEN
               i = i + 1
            ELSE
               EXIT
            ENDIF
         ENDDO

         IF (i == this%nlfp+1) THEN
            depth = this%rivhgt + this%flphgt(this%nlfp) &
               + (v0-this%flpstomax(this%nlfp)) / this%flpaccare(this%nlfp)
         ELSE
            g = (this%flphgt(i)-this%flphgt(i-1))/this%flparea(i)
            depth = this%rivhgt + this%flphgt(i-1) &
               + g * (-this%flpaccare(i-1)+sqrt((this%flpaccare(i-1))**2+2*(v0-this%flpstomax(i-1))/g))
         ENDIF
      ENDIF

   END FUNCTION retrieve_depth_from_volume

   !
   FUNCTION retrieve_volume_from_depth (this, depth) result(volume)

   IMPLICIT NONE

   class(vol_dep_curve_type) :: this
   real(r8), intent(in)  :: depth
   real(r8) :: volume

   ! Local Variables
   real(r8) :: h, d
   integer  :: i

      IF (depth <= this%rivhgt) THEN
         volume = this%rivare * depth
      ELSE
         i = 1
         DO WHILE (i <= this%nlfp)
            IF (depth > this%rivhgt+this%flphgt(i)) THEN
               i = i + 1
            ELSE
               EXIT
            ENDIF
         ENDDO

         d = depth - this%rivhgt - this%flphgt(i-1)
         IF (i == this%nlfp+1) THEN
            volume = this%rivstomax + this%flpstomax(this%nlfp) + d * this%flpaccare(this%nlfp)
         ELSE
            h = this%flphgt(i)-this%flphgt(i-1)
            volume = this%rivstomax + this%flpstomax(i-1) &
               + (d/h*this%flparea(i)+2*this%flpaccare(i-1))*d*0.5
         ENDIF
      ENDIF

   END FUNCTION retrieve_volume_from_depth

   !
   FUNCTION retrieve_area_from_depth (this, depth) result(area)

   IMPLICIT NONE

   class(vol_dep_curve_type) :: this
   real(r8), intent(in)  :: depth
   real(r8) :: area

   ! Local Variables
   real(r8) :: h, d
   integer  :: i

      IF (depth <= this%rivhgt) THEN
         area = 0.
      ELSE
         i = 1
         DO WHILE (i <= this%nlfp)
            IF (depth > this%rivhgt+this%flphgt(i)) THEN
               i = i + 1
            ELSE
               EXIT
            ENDIF
         ENDDO

         IF (i == this%nlfp+1) THEN
            area = this%flpaccare(this%nlfp)
         ELSE
            h = this%flphgt(i)-this%flphgt(i-1)
            d = depth - this%rivhgt - this%flphgt(i-1)
            area = this%flpaccare(i-1) + d/h * this%flparea(i)
         ENDIF
      ENDIF

   END FUNCTION retrieve_area_from_depth

   ! ---
   SUBROUTINE vol_depth_curve_free_mem (this)

   IMPLICIT NONE
   type(vol_dep_curve_type) :: this

      IF (allocated(this%flphgt   )) deallocate (this%flphgt   )
      IF (allocated(this%flparea  )) deallocate (this%flparea  )
      IF (allocated(this%flpaccare)) deallocate (this%flpaccare)
      IF (allocated(this%flpstomax)) deallocate (this%flpstomax)

   END SUBROUTINE vol_depth_curve_free_mem

   ! ---------
   SUBROUTINE riverlake_network_final ()

   IMPLICIT NONE

      IF (allocated(numucat_wrk      )) deallocate(numucat_wrk      )
      IF (allocated(ucat_data_address)) deallocate(ucat_data_address)

      IF (allocated(ucat_elid        )) deallocate(ucat_elid        )
      IF (allocated(ucat_ucid        )) deallocate(ucat_ucid        )
      IF (allocated(inpmat_el2uc     )) deallocate(inpmat_el2uc     )
      IF (allocated(inpmat_area_e2u  )) deallocate(inpmat_area_e2u  )
      IF (allocated(inpmat_uc2el     )) deallocate(inpmat_uc2el     )
      IF (allocated(inpmat_area_u2e  )) deallocate(inpmat_area_u2e  )
      IF (allocated(ucat_next        )) deallocate(ucat_next        )
      IF (allocated(ucat_ups         )) deallocate(ucat_ups         )
      IF (allocated(irivsys          )) deallocate(irivsys          )

      IF (allocated(topo_rivelv      )) deallocate(topo_rivelv      )
      IF (allocated(topo_rivhgt      )) deallocate(topo_rivhgt      )
      IF (allocated(topo_rivlen      )) deallocate(topo_rivlen      )
      IF (allocated(topo_rivman      )) deallocate(topo_rivman      )
      IF (allocated(topo_rivwth      )) deallocate(topo_rivwth      )
      IF (allocated(topo_rivare      )) deallocate(topo_rivare      )
      IF (allocated(topo_rivstomax   )) deallocate(topo_rivstomax   )
      IF (allocated(topo_area        )) deallocate(topo_area        )
      IF (allocated(topo_fldhgt      )) deallocate(topo_fldhgt      )
      IF (allocated(bedelv_next      )) deallocate(bedelv_next      )
      IF (allocated(outletwth        )) deallocate(outletwth        )

      IF (allocated(floodplain_curve )) deallocate(floodplain_curve )

   END SUBROUTINE riverlake_network_final

END MODULE MOD_Grid_RiverLakeNetwork
#endif
