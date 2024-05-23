MODULE cmf_ctrl_sedout_mod
!==========================================================
!* PURPOSE: Output MODULE for CaMa-Flood sediment scheme
! (C) M. Hatono (Hiroshima Univ)  Jan 2023
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not USE this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
   USE PARKIND1,                only: JPIM, JPRB, JPRM
   USE YOS_CMF_INPUT,           only: LOGNAM, NX, NY
   USE YOS_CMF_MAP,             only: REGIONTHIS
   USE CMF_CTRL_OUTPUT_MOD,     only: COUTDIR, IRECOUT, LOUTCDF, LOUTVEC, TVAROUT
   USE yos_cmf_sed,             only: nsed

   IMPLICIT NONE
   SAVE
   type(TVAROUT),allocatable       :: varout(:)          ! output variable type set

   integer(kind=JPIM)              :: nvarsout

   !*** namelist/sediment_output
   character(len=256)         :: csedsout
   namelist/sediment_output/ csedsout

CONTAINS
   !####################################################################
   !-- sediment_output_init
   !-- cmf_sed_output
   !-- sediment_output_end
   !####################################################################
   SUBROUTINE sediment_output_init
   USE CMF_CTRL_OUTPUT_MOD,     only: COUTTAG
   USE CMF_UTILS_MOD,           only: INQUIRE_FID
   USE sed_utils_mod,           only: splitchar
   
   IMPLICIT NONE
   SAVE
   integer(kind=JPIM)              :: jf, j
   integer(kind=JPIM)              :: nvars, nsetfile
   parameter                         (nvars=30)
   character(len=256)              :: cvnames(nvars), fName
      
      nsetfile = INQUIRE_FID()
      open(nsetfile,file='input_sed.nam',status='OLD')
      rewind(nsetfile)
      read(nsetfile,nml=sediment_output)
      close(nsetfile)    

      !---------------------------!
      ! get output variable names !
      !---------------------------!
      cvnames(:) = 'NONE'
      nvarsout = 0
      CALL splitchar(csedsout,cvnames)
      DO j = 1, nvars
         IF ( cvnames(j) /= 'NONE' ) THEN
            nvarsout = nvarsout + 1
         ENDIF
      ENDDO

      IF ( nvarsout == 0 ) THEN
         write(LOGNAM,*) "cmf::sed_output_init: no output files will be produced!"
         RETURN
      ENDIF

      allocate(varout(nvarsout))

      !* loop on variables and create files
      DO jf=1,nvarsout
         write(LOGNAM,*) "creating output for variable:", trim( cvnames(jf) )
         select CASE (cvnames(jf))
            CASE ('sedout')
            varout(jf)%cvname=cvnames(jf)
            varout(jf)%cvlname='suspended sediment flow'
            varout(jf)%cvunits='m3/s'
            CASE ('sedcon')
            varout(jf)%cvname=cvnames(jf)
            varout(jf)%cvlname='suspended sediment concentration'
            varout(jf)%cvunits='m3/m3'
            CASE ('sedinp')
            varout(jf)%cvname=cvnames(jf)
            varout(jf)%cvlname='sediment inflow from land'
            varout(jf)%cvunits='m3/s'
            CASE ('bedout')
            varout(jf)%cvname=cvnames(jf)
            varout(jf)%cvlname='bedload'
            varout(jf)%cvunits='m3/s'
            CASE ('netflw')
            varout(jf)%cvname=cvnames(jf)
            varout(jf)%cvlname='net entrainment flow'
            varout(jf)%cvunits='m3/s'
            CASE ('layer')
            varout(jf)%cvname=cvnames(jf)
            varout(jf)%cvlname='exchange layer volume'
            varout(jf)%cvunits='m3'
            CASE default  ! should only be seddep
            IF ( cvnames(jf)(:6) == 'deplyr' ) THEN
               varout(jf)%cvname=cvnames(jf)
               varout(jf)%cvlname='river bed volume (vertical layer)'
               varout(jf)%cvunits='m3'
            ELSE
               write(LOGNAM,*) trim(cvnames(jf)), 'Not defined in sediment output init'
            ENDIF
         END select
         varout(jf)%binid=INQUIRE_FID()

         IF ( trim(varout(jf)%cvname(:6)) == 'deplyr' ) THEN
            fName = trim(varout(jf)%cvname)//'_'//trim(COUTTAG)
         ELSE
            fName = trim(varout(jf)%cvname)//trim(COUTTAG)
         ENDIF

         IF ( LOUTCDF ) THEN
            IF ( REGIONTHIS==1 ) THEN
            CALL create_outcdf
            ENDIF
         ELSE
            CALL create_outbin
         ENDIF
      ENDDO

   CONTAINS
      SUBROUTINE create_outcdf
#ifdef UseCDF_CMF
      USE YOS_CMF_INPUT,           only: RMIS, CSUFCDF
      USE YOS_CMF_TIME,            only: ISYYYY, ISMM,   ISDD,   ISHOUR, ISMIN
      USE YOS_CMF_MAP,             only: D1LON, D1LAT
      USE CMF_UTILS_MOD,           only: NCERROR
      USE CMF_CTRL_OUTPUT_MOD,     only: NDLEVEL
      USE yos_cmf_sed,             only: sDiam
      USE NETCDF
    
      IMPLICIT NONE
      SAVE
      integer(kind=JPIM)              :: timeid, varid, latid, lonid, sedid
      character(len=256)              :: ctime
 
         varout(jf)%irecnc = 1

         varout(jf)%cfile = trim(COUTDIR)//trim(fName)//trim(CSUFCDF)
         CALL NCERROR( nf90_create(varout(jf)%cfile,nf90_netcdf4,varout(jf)%ncid),&
                     'creating file:'//trim(varout(jf)%cfile) )
         !=== set dimension ===
         CALL NCERROR( nf90_def_dim(varout(jf)%ncid, 'time', nf90_unlimited, timeid) )
         CALL NCERROR( nf90_def_dim(varout(jf)%ncid, 'lat', NY, latid) )
         CALL NCERROR( nf90_def_dim(varout(jf)%ncid, 'lon', NX, lonid) )
         CALL NCERROR( nf90_def_dim(varout(jf)%ncid, 'sedD', nsed, sedid) )   
   
         !=== define variables ===
         CALL NCERROR( nf90_def_var(varout(jf)%ncid, 'sedD', nf90_double, (/sedid/), varid) )
         CALL NCERROR( nf90_put_att(varout(jf)%ncid, varid, 'long_name','sediment grain size') )
         CALL NCERROR( nf90_put_att(varout(jf)%ncid, varid, 'units','meters') )

         CALL NCERROR( nf90_def_var(varout(jf)%ncid, 'lat', nf90_float, (/latid/), varid) )
         CALL NCERROR( nf90_put_att(varout(jf)%ncid, varid, 'long_name','latitude') )
         CALL NCERROR( nf90_put_att(varout(jf)%ncid, varid, 'units','degrees_north') )
      
         CALL NCERROR( nf90_def_var(varout(jf)%ncid, 'lon', nf90_float, (/lonid/), varid) )
         CALL NCERROR( nf90_put_att(varout(jf)%ncid, varid, 'long_name','longitude') )
         CALL NCERROR( nf90_put_att(varout(jf)%ncid, varid, 'units','degrees_east') )
      
         write(ctime,'(a14,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') 'seconds since ',ISYYYY,'-',ISMM,'-',ISDD,' ',ISHOUR,":",ISMIN
         CALL NCERROR( nf90_def_var(varout(jf)%ncid, 'time', nf90_double, (/timeid/), varout(jf)%timid) )
         CALL NCERROR( nf90_put_att(varout(jf)%ncid, varout(jf)%timid, 'long_name','time') )
         CALL NCERROR( nf90_put_att(varout(jf)%ncid, varout(jf)%timid, 'units',ctime) )
      
         !===
         CALL NCERROR( nf90_def_var(varout(jf)%ncid, varout(jf)%cvname, nf90_float, &
                     (/lonid,latid,sedid,timeid/), varout(jf)%varid,deflate_level=ndlevel),     &
                     'creating variable')
      
         CALL NCERROR( nf90_put_att(varout(jf)%ncid, varout(jf)%varid, 'long_name', trim(varout(jf)%cvlname)) )
         CALL NCERROR( nf90_put_att(varout(jf)%ncid, varout(jf)%varid, 'units',     trim(varout(jf)%cvunits)) )
         CALL NCERROR( nf90_put_att(varout(jf)%ncid, varout(jf)%varid, '_fillvalue',rmis) )
      
         CALL NCERROR( nf90_enddef(varout(jf)%ncid) )
      
         !=== put nsed lon lat info ===
         CALL NCERROR ( nf90_inq_varid(varout(jf)%ncid,'sedD',varid),'getting id' )
         CALL NCERROR( nf90_put_var(varout(jf)%ncid,varid,sDiam))

         CALL NCERROR ( nf90_inq_varid(varout(jf)%ncid,'lon',varid),'getting id' )
         CALL NCERROR( nf90_put_var(varout(jf)%ncid,varid,D1LON))
      
         CALL NCERROR ( nf90_inq_varid(varout(jf)%ncid,'lat',varid),'getting id' )
         CALL NCERROR( nf90_put_var(varout(jf)%ncid,varid,D1LAT))
      
         write(LOGNAM,*) 'cfile: ',trim(varout(jf)%cfile),' cvar:',trim(varout(jf)%cvname),&
                     ' clname: ',trim(varout(jf)%cvlname),' cunits: ',trim(varout(jf)%cvunits)
         write(LOGNAM,*) 'open in unit: ',varout(jf)%ncid
#endif
      END SUBROUTINE create_outcdf

      SUBROUTINE create_outbin
      USE YOS_CMF_INPUT,           only: CSUFBIN, CSUFVEC
      USE YOS_CMF_MAP,             only: NSEQMAX, REGIONALL
      
      IMPLICIT NONE
      
         IF ( LOUTVEC ) THEN
            varout(jf)%cfile=trim(coutdir)//trim(fName)//trim(CSUFVEC)
            open(varout(jf)%binid,file=varout(jf)%cfile,form='unformatted',access='direct',recl=4*NSEQMAX*nsed)
         ELSE
            IF ( REGIONTHIS==1 ) THEN
               varout(jf)%cfile=trim(coutdir)//trim(fName)//trim(CSUFBIN)
               open(varout(jf)%binid,file=varout(jf)%cfile,form='unformatted',access='direct',recl=4*NX*NY*nsed)
            ENDIF
         ENDIF
         write(LOGNAM,*) "output file opened in unit: ", TRIM(VAROUT(JF)%CFILE), VAROUT(JF)%BINID
      END SUBROUTINE create_outbin

   END SUBROUTINE sediment_output_init
   !==========================================================
   !+
   !==========================================================
   SUBROUTINE cmf_sed_output
   USE CMF_UTILS_MOD,           only: vecD2mapR
   USE YOS_CMF_INPUT,           only: IFRQ_OUT, RMIS
   USE YOS_CMF_MAP,             only: NSEQMAX
   USE YOS_CMF_TIME,            only: JHOUR, JMIN
   USE yos_cmf_sed,             only: d2layer, d2sedcon, d2seddep, d2bedout_avg, d2netflw_avg, &
                                       d2sedout_avg, d2sedinp_avg, d2sedv_avg, sadd_out
   USE cmf_ctrl_sedrest_mod,    only: sediment_restart_write
#ifdef UseMPI_CMF
   USE CMF_CTRL_MPI_MOD,        only: CMF_MPI_AllReduce_R2MAP
#endif
  
   IMPLICIT NONE
   SAVE
   integer(kind=JPIM)              :: ilyr, ised
   integer(kind=JPIM)              :: jf
   real(kind=JPRB),pointer         :: d2vec(:,:) ! point data location to output
   !*** local
   real(kind=JPRM)                 :: r3out(NX,NY,nsed)
      !================================================
      CALL sediment_restart_write

      d2sedv_avg(:,:,:) = d2sedv_avg(:,:,:) / dble(sadd_out)
      write(LOGNAM,*) 'cmf_sed_output: average ',sadd_out,' seconds'

      !*** 0. check date:hour with output frequency
      IF ( mod(JHOUR,IFRQ_OUT)==0 .and. JMIN==0 ) THEN             ! JHOUR: END of time step , nfpph: output frequency (hour)

         !*** 1. calc average variable
         write(LOGNAM,*) 'cmf::sediment_output_write: write irec: ', IRECOUT

         !*** 2. check variable name & allocate data to pointer dvec
         DO jf=1,nvarsout
            select CASE (varout(jf)%cvname)
            CASE ('sedout')
               d2vec => d2sedout_avg
            CASE ('sedcon')
               d2vec => d2sedcon
            CASE ('sedinp')
               d2vec => d2sedinp_avg
            CASE ('bedout')
               d2vec => d2bedout_avg
            CASE ('netflw')
               d2vec => d2netflw_avg
            CASE ('layer')
               d2vec => d2layer
            CASE default
               IF ( varout(jf)%cvname(:6) == 'deplyr' ) THEN
                  read(varout(jf)%cvname(7:8),*) ilyr
                  d2vec => d2seddep(:,ilyr,:)
               ELSE
                  write(LOGNAM,*) varout(jf)%cvname, ' not defined in cmf_output_mod'
               ENDIF
            END select   !! variable name select

            !! convert 1dvector to 3dmap
            r3out(:,:,:) = RMIS
         
            IF ( .not. LOUTVEC ) THEN
               DO ised = 1, nsed
                  CALL vecD2mapR(d2vec(:,ised),r3out(:,:,ised))             !! mpi node data is gathered by vec2map
#ifdef UseMPI_CMF
                  CALL CMF_MPI_AllReduce_R2MAP(r3out(:,:,ised))
#endif
               ENDDO
        
               IF ( REGIONTHIS==1 ) THEN
                  IF ( LOUTCDF ) THEN
                     CALL wrte_outcdf
                  ELSE
                     CALL wrte_outbin(varout(jf)%binid,IRECOUT,r3out)
                  ENDIF
               ENDIF
            ELSE 
               CALL wrte_outvec(varout(jf)%binid,IRECOUT,d2vec)
             ENDIF
         ENDDO

         write(LOGNAM,*) 'cmf::sediment_output_write: END'
      ENDIF

      d2sedv_avg(:,:,:) = 0._JPRB
      sadd_out = 0._JPRB

   CONTAINS
      SUBROUTINE wrte_outcdf
#ifdef UseCDF_CMF
      USE NETCDF
      USE YOS_CMF_TIME,            only: KMINSTART, KMINNEXT
      USE CMF_UTILS_MOD,           only: NCERROR
    
      IMPLICIT NONE
      SAVE
      real(kind=JPRB)                 :: xtime
         xtime = real( (KMINNEXT-KMINSTART), JPRB) *60._JPRB
         CALL NCERROR( nf90_put_var(varout(jf)%ncid,varout(jf)%timid,xtime,(/varout(jf)%irecnc/)) )

         CALL NCERROR( nf90_put_var(varout(jf)%ncid,varout(jf)%varid,r3out(1:NX,1:NY,1:nsed),&
                     (/1,1,1,varout(jf)%irecnc/),(/NX,NY,nsed,1/)) )
      
         ! update irec
         varout(jf)%irecnc=varout(jf)%irecnc+1
#endif
      END SUBROUTINE wrte_outcdf
      !==========================================================
      SUBROUTINE wrte_outbin(ifn,irec,r2outdat)

      IMPLICIT NONE
      !*** input
      SAVE
      integer(kind=JPIM),intent(in)   :: ifn                 !! file number
      integer(kind=JPIM),intent(in)   :: irec                !! record
      real(kind=JPRM)                 :: r2outdat(NX,NY,nsed)
      !================================================
         write(ifn,rec=irec) r2outdat
      END SUBROUTINE wrte_outbin
      !==========================================================
      SUBROUTINE wrte_outvec(ifn,irec,d2outdat)
    
      IMPLICIT NONE
      !*** input
      SAVE
      integer(kind=JPIM),intent(in)   :: ifn                 !! file number
      integer(kind=JPIM),intent(in)   :: irec                !! record
      real(kind=JPRB),intent(in)      :: d2outdat(NSEQMAX,nsed) !! output data
      !*** local
      real(kind=JPRM)                 :: r2outdat(NSEQMAX,nsed)
      !================================================
         r2outdat(:,:)=real(d2outdat(:,:))
         write(ifn,rec=irec) r2outdat
      END SUBROUTINE wrte_outvec
   !==========================================================
   END SUBROUTINE cmf_sed_output
   !==========================================================
   !+
   !==========================================================
   SUBROUTINE sediment_output_end
#ifdef UseCDF_CMF
   USE NETCDF
   USE CMF_UTILS_MOD,           only: NCERROR
#endif
   USE YOS_CMF_MAP,             only: REGIONTHIS
  
   IMPLICIT NONE
   SAVE
   integer(kind=JPIM)              :: jf

      write(LOGNAM,*) ""
      write(LOGNAM,*) "!---------------------!"
      write(LOGNAM,*) "sediment_output_end: finalize output MODULE"

      IF ( LOUTVEC ) THEN
         DO jf = 1, nvarsout
            close(varout(jf)%binid)
         ENDDO
      ELSE IF ( REGIONTHIS==1 ) THEN
         DO jf = 1, nvarsout
            IF ( LOUTCDF ) THEN
#ifdef UseCDF_CMF
               CALL NCERROR( nf90_close(varout(jf)%ncid) )
#endif
            ELSE
               close(varout(jf)%binid)
            ENDIF
         ENDDO
      ENDIF 
  
      write(LOGNAM,*) 'sediment_output_end: END'
   END SUBROUTINE sediment_output_end
  
END MODULE cmf_ctrl_sedout_mod
