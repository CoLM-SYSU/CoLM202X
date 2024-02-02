MODULE CMF_CTRL_SED_MOD
!==========================================================
!* PURPOSE: physics for sediment transport
! (C) M.Hatono  (Hiroshima-U)  May 2021
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not USE this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
#ifdef UseMPI_CMF
   USE MPI
#endif
   USE PARKIND1,                only: JPIM, JPRB, JPRM
   USE YOS_CMF_INPUT,           only: LOGNAM, NX, NY
   USE YOS_CMF_MAP,             only: MPI_COMM_CAMA, NSEQALL, NSEQMAX, REGIONTHIS
   USE CMF_UTILS_MOD,           only: INQUIRE_FID
   USE yos_cmf_sed,             only: lyrdph, nsed, totlyrnum, sDiam

   IMPLICIT NONE
   !*** namelist/sediment_map/
   character(len=256)              :: crocdph
   character(len=256)              :: csedfrc
   character(len=256)              :: sedD
   namelist/sediment_map/ crocdph,sedD,csedfrc
CONTAINS
   !####################################################################
   ! -- cmf_sed_nmlist
   ! -- cmf_sed_init
   !####################################################################
   SUBROUTINE cmf_sed_nmlist
   USE yos_cmf_sed,             only: lambda, psedD, pset, pwatD, sedDT, &
                                       revEgia, visKin, vonKar

   IMPLICIT NONE
   integer(kind=JPIM)              :: nsetfile

   namelist/sediment_param/  lambda, lyrdph, nsed, sedDT, psedD, &
                              pset, pwatD, revEgia, totlyrnum, &
                              visKin, vonKar

      nsetfile = INQUIRE_FID()
      open(nsetfile,file='input_sed.nam',status='OLD')
      
      lambda = 0.4d0
      lyrdph = 0.00005d0
      nsed = 3
      sedDT = 3600
      psedD = 2.65d0
      pset = 1.d0
      pwatD = 1.d0
      revEgia = .true.
      totlyrnum = 5
      visKin = 1.d-6
      vonKar = 0.4d0

      rewind(nsetfile)
      read(nsetfile,nml=sediment_param)
      !defaults
      write(LOGNAM,*) 'nml sediment_param'
      write(LOGNAM,*) 'lambda    :', lambda
      write(LOGNAM,*) 'lyrdph    :', lyrdph
      write(LOGNAM,*) 'sedDT     :', sedDT
      write(LOGNAM,*) 'psedD     :', psedD
      write(LOGNAM,*) 'pset      :', pset
      write(LOGNAM,*) 'pwatD     :', pwatD
      write(LOGNAM,*) 'revEgia   :', revEgia
      write(LOGNAM,*) 'totlyrnum :', totlyrnum
      write(LOGNAM,*) 'visKin    :', visKin
      write(LOGNAM,*) 'vonKar    :', vonKar

      rewind(nsetfile)
      read(nsetfile,nml=sediment_map)
      !defaults
      write(LOGNAM,*) 'nml sediment_map'
      write(LOGNAM,*) 'crocdph  :', trim(crocdph)
      write(LOGNAM,*) 'sDiam   :', sedD
      write(LOGNAM,*) 'csedfrc  :', csedfrc

      close(nsetfile)
   END SUBROUTINE cmf_sed_nmlist
   !==========================================================
   !+
   !==========================================================
   SUBROUTINE cmf_sed_init
   USE YOS_CMF_INPUT,           only: LOUTPUT
   USE cmf_ctrl_sedinp_mod,     only: sediment_input_init
   USE cmf_ctrl_sedout_mod,     only: sediment_output_init
   USE cmf_ctrl_sedrest_mod,    only: sediment_restart_init
      
   IMPLICIT NONE

      CALL sediment_vars_init

      CALL sediment_map_init

      CALL sediment_input_init

      IF ( LOUTPUT ) THEN
         CALL sediment_output_init
      ENDIF

      CALL sediment_restart_init

   CONTAINS
      !==================================
      SUBROUTINE sediment_map_init
      USE YOS_CMF_INPUT,           only: NLFP, PGRV
      USE CMF_UTILS_MOD,           only: mapR2vecD
      USE yos_cmf_sed,             only: d2sedfrc, psedD, pset, pwatD, setVel, visKin
      USE cmf_calc_sedpar_mod,     only: calc_settingVelocity
      USE sed_utils_mod,           only: splitchar

      IMPLICIT NONE
      integer(kind=JPIM)              :: i, ierr, ised, iseq, tmpnam
      real(kind=JPRM)                 :: r2temp(NX,NY), sTmp1(nsed)
      character(len=256)              :: ctmp(20)

         !------------------------!
         ! get sediment diameters !
         !------------------------!
         ctmp(:) = '-999'
         CALL splitchar(sedD,ctmp)
         ised = 0
         allocate(sDiam(nsed))
         DO i = 1, nsed
            IF ( ctmp(i) /= '-999' ) THEN
               ised = ised + 1
               read(ctmp(i),*) sDiam(ised)
            ENDIF
         ENDDO
         IF ( ised /= nsed ) THEN
            write(LOGNAM,*) 'nsed and sedD DO not match',ised,nsed
            STOP
         ENDIF
         write(LOGNAM,*) ised,' grain sizes: ',sDiam(:)

         !----------------------------!
         ! calculate setting velocity !
         !----------------------------!
         allocate(setVel(nsed))
         setVel(:) = calc_settingVelocity()

         !-----------------------------!
         ! read sediment fraction file !
         !-----------------------------!
         allocate(d2sedfrc(NSEQMAX,nsed))
         IF ( REGIONTHIS == 1 ) THEN
            tmpnam = INQUIRE_FID()
            open(tmpnam,file=csedfrc,form='unformatted',access='direct',recl=4*NX*NY)
         ENDIF
         DO ised = 1, nsed
            IF ( REGIONTHIS == 1 ) read(tmpnam,rec=ised) r2temp
#ifdef UseMPI_CMF
            CALL MPI_Bcast(r2temp(1,1),NX*NY,mpi_real4,0,MPI_COMM_CAMA,ierr)
#endif
            CALL mapR2vecD(r2temp,d2sedfrc(:,ised))
         ENDDO
         IF ( REGIONTHIS == 1 ) close(tmpnam)
      
         ! adjust if any fractions are negative or if sum is not equal to 1
         IF ( nsed == 1 ) THEN
            d2sedfrc(:,:) = 1.d0
         ELSE
!$omp parallel DO
            DO iseq = 1, NSEQALL
               IF ( minval(d2sedfrc(iseq,:)) < 0.d0 .or. sum(d2sedfrc(iseq,:)) == 0.d0 ) THEN
                  d2sedfrc(iseq,:) = 1.d0 / dble(nsed)
               ELSE IF ( sum(d2sedfrc(iseq,:)) /= 1.d0 ) THEN
                  d2sedfrc(iseq,:) = d2sedfrc(iseq,:) / sum(d2sedfrc(iseq,:))
               ENDIF
            ENDDO
!$omp END parallel do
         ENDIF
      END SUBROUTINE sediment_map_init
   !==================================
      SUBROUTINE sediment_vars_init
      USE yos_cmf_sed,             only: d2bedout, d2netflw, d2seddep, &
                                          d2bedout_avg, d2netflw_avg,   &
                                          d2sedout, d2sedcon, d2sedinp, &
                                          d2sedout_avg, d2sedinp_avg, d2layer, &
                                          d2sedv, d2sedv_avg, d2depv,   &
                                          sedDT, step_sed
      USE YOS_CMF_INPUT,           only: DT
      IMPLICIT NONE

         IF ( mod(sedDT,DT) /= 0 ) THEN
            write(lognam,*) 'sedDT ',sedDT,'is not a multiple of DT',DT
            STOP
         ENDIF
         step_sed = int(sedDT/DT)

         allocate(d2sedv(NSEQMAX,nsed,6))
         d2sedv(:,:,:) = 0._JPRB
         d2sedout => d2sedv(:,:,1)
         d2sedcon => d2sedv(:,:,2)
         d2sedinp => d2sedv(:,:,3)
         d2bedout => d2sedv(:,:,4)
         d2netflw => d2sedv(:,:,5)
         d2layer => d2sedv(:,:,6)

         allocate(d2depv(NSEQMAX,totlyrnum,nsed))
         d2depv(:,:,:) = 0._JPRB
         d2seddep => d2depv

         allocate(d2sedv_avg(NSEQMAX,nsed,4))
         d2sedv_avg(:,:,:) = 0._JPRB
         d2sedout_avg => d2sedv_avg(:,:,1)
         d2sedinp_avg => d2sedv_avg(:,:,2)
         d2bedout_avg => d2sedv_avg(:,:,3)
         d2netflw_avg => d2sedv_avg(:,:,4)
      END SUBROUTINE sediment_vars_init
   END SUBROUTINE cmf_sed_init
!####################################################################
END MODULE CMF_CTRL_SED_MOD
