MODULE cmf_ctrl_sedinp_mod
!==========================================================
!* PURPOSE: Manage sediment input
! (C) M.Hatono  (Hiroshima-U)  Oct 2022
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
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
   USE YOS_CMF_INPUT,           only: LOGNAM
   USE YOS_CMF_MAP,             only: NSEQALL, NSEQMAX, D2GRAREA
   USE CMF_CTRL_FORCING_MOD,    only: INPX, INPY, INPA

   IMPLICIT NONE
   SAVE
   character(len=256)              :: sedinput_dir, sedinput_pre, sedinput_suf

   real(kind=JPRB),allocatable     :: d2slope(:,:)     ! floodplain slope [deg]
   integer(kind=JPIM)              :: iseq

   real(kind=JPRB)                 :: dsylunit         ! unit conversion for sediment [m3/km2] -> [m3/m2]
   real(kind=JPRB)                 :: pyld, pyldc, pyldpc  ! parameters for sediment erosion calculation

CONTAINS
   !####################################################################
   !-- sediment_input_init
   !-- cmf_sed_forcing
   !-- calc_sedyld
   !-- sedinp_interp
   !####################################################################
   SUBROUTINE sediment_input_init
#ifdef UseMPI_CMF
   USE YOS_CMF_MAP,             only: MPI_COMM_CAMA
#endif
   USE CMF_UTILS_MOD,           only: INQUIRE_FID, mapR2vecD

   IMPLICIT NONE
   SAVE
   character(len=256)              :: cslope       ! slope file
   character(len=256)              :: cinpmat_sed  ! input matrix for sediment

      CALL read_sedinp_nmlist

      CALL read_slope
   CONTAINS
      !================================
      SUBROUTINE read_sedinp_nmlist
      IMPLICIT NONE
      integer(kind=JPIM)            :: nsetfile
      
      namelist/sediment_input/ sedinput_dir, sedinput_pre, sedinput_suf, &
                                 cslope, dsylunit, pyld, pyldc, pyldpc,    &
                                 cinpmat_sed

         nsetfile = INQUIRE_FID()
         open(nsetfile,file='input_sed.nam',status='OLD')

         sedinput_dir='./'
         sedinput_pre='./'
         sedinput_suf='./'
         cslope='./slope.bin'
         dsylunit = 1.d-6
         pyld = 0.01d0
         pyldc = 2.d0
         pyldpc = 2.d0
         cinpmat_sed = './inpmat.bin'

         rewind(nsetfile)
         read(nsetfile,nml=sediment_input)
         !defaults
         write(LOGNAM,*) 'nml sediment_input'
         write(LOGNAM,*) 'cslope    :', trim(cslope)
         write(LOGNAM,*) 'dsylunit  :', dsylunit
         write(LOGNAM,*) 'pyld      :', pyld
         write(LOGNAM,*) 'pyldc     :', pyldc
         write(LOGNAM,*) 'pyldpc    :', pyldpc
         write(LOGNAM,*) 'cinpmat_sed:', trim(cinpmat_sed)
      END SUBROUTINE

      SUBROUTINE read_slope
      USE YOS_CMF_INPUT,         only: NX,NY, NLFP
      USE YOS_CMF_MAP,           only: REGIONTHIS

      IMPLICIT NONE
      integer                       :: ierr, tmpnam, i
      real(kind=jprm)               :: r2temp(nx,ny)
         allocate(d2slope(NSEQMAX,NLFP))
         IF ( REGIONTHIS == 1 ) THEN
            tmpnam = INQUIRE_FID()
            open(tmpnam,file=cslope,form='unformatted',access='direct',recl=4*NX*NY)
         ENDIF
         DO i = 1, NLFP
            IF ( REGIONTHIS == 1 ) read(tmpnam,rec=i) r2temp
#ifdef UseMPI_CMF
            CALL MPI_Bcast(r2temp(1,1),NX*NY,mpi_real4,0,MPI_COMM_CAMA,ierr)
#endif
            CALL mapR2vecD(r2temp,d2slope(:,i))
         ENDDO
         IF ( REGIONTHIS == 1 ) close(tmpnam)
      END SUBROUTINE read_slope
  
   END SUBROUTINE sediment_input_init
   !==========================================================
   !+
   !==========================================================
   SUBROUTINE cmf_sed_forcing
   ! read forcing from file
   USE YOS_CMF_INPUT,           only: TMPNAM,NXIN,NYIN,DTIN
   USE YOS_CMF_TIME,            only: IYYYY, IMM, IDD, IHOUR, IMIN
   USE CMF_UTILS_MOD,           only: CONV_END,INQUIRE_FID
  
   IMPLICIT NONE
   SAVE
   real(kind=JPRB)                 :: d2temp(nseqmax)
   !* local variables
   integer(kind=jpim)              :: irecinp
   integer(kind=jpim)              :: isec
   character(len=256)              :: cifname             !! input file
   character(len=256)              :: cdate               !!
   real(kind=jprm)                 :: r2tmp(nxin,nyin)
      
      !*** 1. calculate irec for sub-daily precipitation
      isec    = ihour*60*60+imin*60   !! current second in a day
      irecinp = int( isec/dtin ) +1   !! precipitation irec (sub-daily precipitation)

      !*** 2. set file name
      write(cdate,'(i4.4,i2.2,i2.2)') iyyyy,imm,idd
      cifname=trim(sedinput_dir)//'/'//trim(sedinput_pre)//trim(cdate)//trim(sedinput_suf)
      write(LOGNAM,*) "cmf::sed_forcing_get_bin:",trim(cifname)

      !*** 3. open & read forcing data
      tmpnam=inquire_fid()
      open(tmpnam,file=cifname,form='unformatted',access='direct',recl=4*nxin*nyin)
      read(tmpnam,rec=irecinp) r2tmp
      close(tmpnam)

      !*** 5. conduct necessary conversion
      CALL sedinp_interp(r2tmp,d2temp)   ! interpolate forcing grid to model grid
      CALL calc_sedyld(d2temp)           ! calculate sediment yield into rivers
   END SUBROUTINE cmf_sed_forcing
   !==========================================================
   !+
   !==========================================================
   SUBROUTINE calc_sedyld(pbuffin)
   USE PARKIND1,                only: JPIM, JPRB
   USE YOS_CMF_INPUT,           only: DTIN
   USE yos_cmf_sed,             only: d2sedinp, d2sedinp_avg, d2sedfrc
   
   IMPLICIT NONE
   SAVE
   real(kind=JPRB), intent(in)     :: pbuffin(:)
   real(kind=JPRB)                 :: sbuff(NSEQMAX)
  !================================================

      
      CALL prcp_convert_sed(pbuffin, sbuff) ! convert precipitation to sediment yield based on Sunada&Hasegawa(1993)
 
  !$omp parallel DO 
      DO iseq = 1, NSEQALL
         d2sedinp(iseq,:) = sbuff(iseq) * d2sedfrc(iseq,:)  ! distribute sediment yield to proportionate to sediment grain fraction
         d2sedinp_avg(iseq,:) = d2sedinp_avg(iseq,:) + d2sedinp(iseq,:) * DTIN
      ENDDO
  !$omp end parallel DO

   CONTAINS
      !=============================
      !+ prcp_convert_sed
      !=============================
      SUBROUTINE prcp_convert_sed(pbuffin,pbuffout)
      USE YOS_CMF_DIAG,          only: D2FLDFRC
      USE YOS_CMF_INPUT,         only: NLFP

      IMPLICIT NONE
      SAVE
      real(kind=JPRB), intent(in)   :: pbuffin(:)     !! kg/m2/s
      real(kind=JPRB), intent(out)  :: pbuffout(:)  !! m3/s
      integer(kind=JPIM)            :: i, iseq

!$omp parallel DO
         DO iseq = 1, NSEQALL
            pbuffout(iseq) = 0.d0
            IF ( pbuffin(iseq) * 86400.d0 <= 10.d0 ) CYCLE

            DO i = 1, NLFP
               IF ( D2FLDFRC(iseq,1) * NLFP > dble(i) ) CYCLE  ! no erosion if submerged
               pbuffout(iseq) = pbuffout(iseq) + pyld * (pbuffin(iseq)*3600.d0)**pyldpc * d2slope(iseq,i)**pyldc / 3600.d0 & 
                  & * D2GRAREA(iseq,1) * min(dble(i)/dble(NLFP)-D2FLDFRC(iseq,1), 1.d0/dble(NLFP)) * dsylunit
            ENDDO
         ENDDO
!$omp end parallel do
      END SUBROUTINE prcp_convert_sed
    
   END SUBROUTINE calc_sedyld
   !==========================================================
   !+
   !==========================================================
   SUBROUTINE sedinp_interp(pbuffin,pbuffout)
   ! interporlate sediment forcing data using "input matrix"
   USE YOS_CMF_INPUT,           only: NXIN, NYIN, INPN, RMIS
  
   IMPLICIT NONE
   real(kind=JPRM),intent(in)      :: pbuffin(:,:)     !! default for prcp[kg/m2/s]
   real(kind=JPRB),intent(out)     :: pbuffout(:)    !! kg/m2/s
   ! save for omp
   integer(kind=jpim),SAVE  ::  iseq, ixin, iyin, inpi  !! for output
!$omp threadprivate    (ixin, iyin)
  !============================
!$omp parallel do
      DO iseq=1, NSEQALL
         pbuffout(iseq)=0._JPRB
         DO inpi=1, INPN
            ixin=INPX(iseq,inpi)
            iyin=INPY(iseq,inpi)
            IF( ixin>0 )THEN
               IF( ixin > NXIN .or. iyin > NYIN ) THEN
                  write(LOGNAM,*)  "error"
                  write(LOGNAM,*)  'xxx',iseq,inpi,ixin,iyin
                  CYCLE
               ENDIF
               IF( pbuffin(ixin,iyin).ne.RMIS )THEN
                  pbuffout(iseq) = pbuffout(iseq) + pbuffin(ixin,iyin) * INPA(iseq,inpi) / D2GRAREA(iseq,1)
               ENDIF
            ENDIF
         ENDDO
         pbuffout(iseq)=max(pbuffout(iseq), 0._JPRB)
      ENDDO
!$omp END parallel do
   END SUBROUTINE sedinp_interp
!####################################################################

END MODULE cmf_ctrl_sedinp_mod
