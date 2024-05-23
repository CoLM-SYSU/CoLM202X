MODULE YOS_CMF_TIME
!==========================================================
!* PURPOSE: Shared time-related variables
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
   !==========================================================
   USE PARKIND1,                only: JPIM, JPRB, JPRM
   IMPLICIT NONE
   !======================================
   SAVE
   ! simulation time step
   integer(KIND=JPIM)              :: KSTEP              !! time step since start
   integer(KIND=JPIM)              :: NSTEPS             !! total time step (from start to end), given in CMF_TIME_INIT
   ! elapsed minute from base date (YYYY0,MM0,DD0)
   integer(KIND=JPIM)              :: KMIN               !! KMIN at the start of time step
   integer(KIND=JPIM)              :: KMINNEXT           !! KMIN at the end   of time step
   !
   integer(KIND=JPIM)              :: KMINSTART          !! KMIN at the start of simulation
   integer(KIND=JPIM)              :: KMINEND            !! KMIN at the end   of simulation
   !
   integer(KIND=JPIM)              :: KMINSTAIN          !! KMIN at the start of forcing runoff  data (netCDF)
   integer(KIND=JPIM)              :: KMINSTASL          !! KMIN at the start of boundary sealev data (netCDF)
   ! simulation start date:hour (KMINSTART)
   integer(KIND=JPIM)              :: ISYYYYMMDD         !! date     at simulation start time
   integer(KIND=JPIM)              :: ISHHMM             !! hour+min at simulation start time
   integer(KIND=JPIM)              :: ISYYYY
   integer(KIND=JPIM)              :: ISMM
   integer(KIND=JPIM)              :: ISDD
   integer(KIND=JPIM)              :: ISHOUR
   integer(KIND=JPIM)              :: ISMIN
   ! simulation end   date:hour (KMINEND)
   integer(KIND=JPIM)              :: IEYYYYMMDD         !! date     of simulation end time
   integer(KIND=JPIM)              :: IEHHMM             !! hour+min of simulation end time
   integer(KIND=JPIM)              :: IEYYYY
   integer(KIND=JPIM)              :: IEMM
   integer(KIND=JPIM)              :: IEDD
   integer(KIND=JPIM)              :: IEHOUR
   integer(KIND=JPIM)              :: IEMIN
   !*** date:hour at START of time steop (KMIN)
   integer(KIND=JPIM)              :: IYYYYMMDD          !! date     at the start of time-step
   integer(KIND=JPIM)              :: IYYYY              !! year     at the start of time-step
   integer(KIND=JPIM)              :: IMM                !! month    at the start of time-step
   integer(KIND=JPIM)              :: IDD                !! day      at the start of time-step
   integer(KIND=JPIM)              :: IHHMM              !! hour+min at the start of time-step
   integer(KIND=JPIM)              :: IHOUR              !! hour     at the start of time-step
   integer(KIND=JPIM)              :: IMIN               !! min      at the start of time-step
   !*** date:hour at END   of time steop (KMINNEXT)
   integer(KIND=JPIM)              :: JYYYYMMDD          !! date     at the end   of time-step
   integer(KIND=JPIM)              :: JYYYY              !! year     at the end   of time-step
   integer(KIND=JPIM)              :: JMM                !! month    at the end   of time-step
   integer(KIND=JPIM)              :: JDD                !! day      at the end   of time-step
   integer(KIND=JPIM)              :: JHHMM              !! hour+min at the end   of time-step
   integer(KIND=JPIM)              :: JHOUR              !! hour     at the end   of time-step
   integer(KIND=JPIM)              :: JMIN               !! min      at the end   of time-step

   !*** base time to define kmin
   integer(KIND=JPIM)              :: YYYY0              !! base year
   integer(KIND=JPIM)              :: MM0                !! base month
   integer(KIND=JPIM)              :: DD0                !! base day
   !==========================================================
END MODULE YOS_CMF_TIME
