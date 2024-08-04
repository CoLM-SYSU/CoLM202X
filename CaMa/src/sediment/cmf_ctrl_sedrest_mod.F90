MODULE cmf_ctrl_sedrest_mod
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
   USE PARKIND1,                only: JPIM, JPRM
   USE YOS_CMF_INPUT,           only: LOGNAM, NX, NY
   USE YOS_CMF_MAP,             only: MPI_COMM_CAMA, NSEQALL, NSEQMAX, REGIONTHIS
   USE CMF_UTILS_MOD,           only: INQUIRE_FID
   USE yos_cmf_sed,             only: nsed, totlyrnum, d2layer, d2sedcon, d2seddep

   IMPLICIT NONE
   SAVE
   integer(kind=JPIM)              :: ifrq_rst_sed
   character(len=256)              :: sedrest_infile, sedrest_outpre

   namelist/sediment_restart/  sedrest_infile, sedrest_outpre, ifrq_rst_sed

CONTAINS
   !####################################################################
   ! -- sediment_restart_init
   ! -- sediment_restart_write
   ! --
   !####################################################################
   SUBROUTINE sediment_restart_init
   USE YOS_CMF_MAP,             only: D2RIVLEN, D2RIVWTH
   USE YOS_CMF_PROG,            only: P2RIVSTO
   USE yos_cmf_sed,             only: d2sedfrc, lyrdph, d2rivsto_pre, &
                                       d2rivout_sed, d2rivvel_sed, sadd_riv, sadd_out
   USE CMF_UTILS_MOD,           only: mapR2vecD, inquire_fid

   IMPLICIT NONE
   SAVE
#ifdef UseMPI_CMF
   integer(kind=JPIM)             :: ierr
#endif
   integer(kind=JPIM)             :: ilyr, irec, ised, iseq, tmpnam, nsetfile
   real(kind=JPRM)                :: r2temp(NX,NY)

      nsetfile = inquire_fid()
      open(nsetfile, file='input_sed.nam', status='old')
      rewind(nsetfile)
      read(nsetfile,nml=sediment_restart)
      close(nsetfile)

      IF ( sedrest_infile == "" ) THEN  ! set layer/bedload if no restart file
!$omp parallel DO
         DO iseq = 1, NSEQALL
            d2layer(iseq,:) = lyrdph * D2RIVWTH(iseq,1) * D2RIVLEN(iseq,1) * d2sedfrc(iseq,:)
            DO ilyr = 1, totlyrnum-1
               d2seddep(iseq,ilyr,:) = d2layer(iseq,:)
            ENDDO
            d2seddep(iseq,totlyrnum,:) = ( max(10.d0-lyrdph*totlyrnum,0.d0) ) * D2RIVWTH(iseq,1) * D2RIVLEN(iseq,1) * d2sedfrc(iseq,:)
         ENDDO
!$omp END parallel DO

      ELSE
         IF ( REGIONTHIS == 1 ) THEN
            tmpnam = INQUIRE_FID()
            open(tmpnam,file=sedrest_infile,form='unformatted',access='direct',recl=4*NX*NY)
         ENDIF
         DO irec = 1, 2
            DO ised = 1, nsed
               IF ( REGIONTHIS == 1 ) read(tmpnam,rec=(irec-1)*nsed+ised) r2temp
#ifdef UseMPI_CMF
               CALL MPI_Bcast(r2temp(1,1),NX*NY,mpi_real4,0,MPI_COMM_CAMA,ierr)
#endif
               select CASE(irec)
                  CASE (1)
                     CALL mapR2vecD(r2temp,d2layer(:,ised))
                  CASE (2)
                     CALL mapR2vecD(r2temp,d2sedcon(:,ised))
               END select
            ENDDO
         ENDDO

         DO irec = 1, totlyrnum
            DO ised = 1, nsed
               IF ( REGIONTHIS == 1 ) read(tmpnam,rec=(irec+1)*nsed+ised) r2temp
#ifdef UseMPI_CMF
               CALL MPI_Bcast(r2temp(1,1),NX*NY,mpi_real4,0,MPI_COMM_CAMA,ierr)
#endif
               CALL mapR2vecD(r2temp,d2seddep(:,irec,ised))
            ENDDO
         ENDDO
         IF ( REGIONTHIS == 1 ) close(tmpnam)
         write(LOGNAM,*) 'read restart sediment',maxval(d2seddep(:,totlyrnum,:))
      ENDIF

      allocate(d2rivsto_pre(NSEQMAX), d2rivout_sed(NSEQMAX), d2rivvel_sed(NSEQMAX))
      sadd_riv = 0.d0
      sadd_out = 0.d0
      d2rivsto_pre(:) = P2RIVSTO(:,1)
      d2rivout_sed(:) = 0.d0
      d2rivvel_sed(:) = 0.d0
   END SUBROUTINE sediment_restart_init
!==========================================================
!+
!==========================================================
   SUBROUTINE sediment_restart_write
   USE YOS_CMF_TIME,            only: KSTEP, NSTEPS, JDD, JHHMM, JHOUR, JMIN, JYYYYMMDD
   USE YOS_CMF_INPUT,           only: CSUFBIN, RMIS
   USE CMF_CTRL_RESTART_MOD,    only: CRESTDIR
   USE CMF_UTILS_MOD,           only: vecD2mapR
#ifdef UseMPI_CMF
   USE CMF_CTRL_MPI_MOD,        only: CMF_MPI_AllReduce_R2MAP
#endif
  
   IMPLICIT NONE
   SAVE
   integer(kind=JPIM)              :: irec, irest, ised, tmpnam
   real(kind=JPRM)                 :: r3final(NX,NY,nsed), r2temp(NX,NY)
   character(len=256)              :: cdate, cfile

      irest = 0

      IF ( ifrq_rst_sed>=0 .and. KSTEP==NSTEPS ) THEN  !! END of run
         irest = 1
      ENDIF

      IF ( ifrq_rst_sed>=1 .and. ifrq_rst_sed<=24 ) THEN  !! at selected hour
         IF ( mod(JHOUR,ifrq_rst_sed)==0 .and. JMIN==0 ) THEN
            irest = 1
         ENDIF
      ENDIF

      IF ( ifrq_rst_sed==30 ) THEN  !! at END of month
         IF ( JDD==1 .and. JHOUR==0 .and. JMIN==0 ) THEN
            irest = 1
         ENDIF
      ENDIF

      IF ( irest==1 ) THEN
         write(LOGNAM,*) ""
         write(LOGNAM,*) "!---------------------!"
         write(LOGNAM,*) 'cmf::sediment_restart_write: write time: ' , JYYYYMMDD, JHHMM

         write(cdate,'(I8.8,I2.2)') JYYYYMMDD,JHOUR
         cfile=trim(CRESTDIR)//TRIM(sedrest_outpre)//TRIM(cdate)//TRIM(CSUFBIN)
         write(LOGNAM,*) 'wrte_rest_bin: restart file:',cfile

         !*** write restart data (2D map)
         IF ( REGIONTHIS == 1 ) THEN
            tmpnam = INQUIRE_FID()
            open(TMPNAM,file=cfile,form='unformatted',access='direct',recl=4*NX*NY*nsed)
         ENDIF
         DO irec = 1, 2
            r3final(:,:,:) = RMIS
            DO ised = 1, nsed
               select CASE(irec)
                  CASE (1)
                     CALL vecD2mapR(d2layer(:,ised),r2temp)
                  CASE (2)
                     CALL vecD2mapR(d2sedcon(:,ised),r2temp)
               END select
#ifdef UseMPI_CMF
               CALL CMF_MPI_AllReduce_R2MAP(r2temp)
#endif
               r3final(:,:,ised) = r2temp(:,:)
            ENDDO
            IF ( REGIONTHIS == 1 ) write(tmpnam,rec=irec) r3final
         ENDDO

         DO irec = 1, totlyrnum
            r3final(:,:,:) = RMIS
            DO ised = 1, nsed
            CALL vecD2mapR(d2seddep(:,irec,ised),r2temp)
#ifdef UseMPI_CMF
            CALL CMF_MPI_AllReduce_R2MAP(r2temp)
#endif
            r3final(:,:,ised) = r2temp
            ENDDO
            IF ( REGIONTHIS == 1 ) write(tmpnam,rec=irec+2) r3final
         ENDDO

         IF ( REGIONTHIS == 1 ) close(tmpnam)

      ENDIF
   END SUBROUTINE sediment_restart_write
   !####################################################################

END MODULE cmf_ctrl_sedrest_mod
