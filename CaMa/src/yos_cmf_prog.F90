MODULE YOS_CMF_PROG
!==========================================================
!* PURPOSE: Shared prognostic variables for CaMa-Flood
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  2Aug 2019

! Modified by Zhongwang Wei @ SYSU 2022.11.20: add water re-infiltration calculation 

! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
   USE PARKIND1, only: JPIM, JPRB, JPRM, JPRD
   IMPLICIT NONE
   SAVE
   !================================================
   ! Pointer was removed in v4.08 in order to keep simple codes when activating Single Precision Mode
   !*** prognostics / state variables initial conditions

   ! Dammy variable for input/output
   real(KIND=JPRB),allocatable,target     :: D2DAMMY(:,:)       !! Dammy Array for unused variables
   real(KIND=JPRB),allocatable,target     :: D2COPY(:,:)        !! Dammy Array for Float64/32 switch

   !================================================
   !*** input runoff (interporlated)
   real(KIND=JPRB),allocatable,target     :: D2RUNOFF(:,:)         !! input runoff             [m3/s]
   real(KIND=JPRB),allocatable,target     :: D2ROFSUB(:,:)         !! input sub-surface runoff [m3/s]

   !TODO: check d2wevap and d2winfilt units!
   real(KIND=JPRB),allocatable,target     :: D2WEVAP(:,:)          !! input Evaporation [m3/s]
   real(KIND=JPRB),allocatable,target     :: D2WINFILT(:,:)        !! input Infiltration [m3/s]

   !================================================
   !*** river & floodpain
   ! storage variables are always in double precision
   real(KIND=JPRD),allocatable,target     :: P2RIVSTO(:,:)         !! river      storage [m3]
   real(KIND=JPRD),allocatable,target     :: P2FLDSTO(:,:)         !! floodplain storage [m3]

   real(KIND=JPRB),allocatable,target     :: D2RIVOUT(:,:)         !! river      outflow [m3/s]
   real(KIND=JPRB),allocatable,target     :: D2FLDOUT(:,:)         !! floodplain outflow [m3/s]

   !================================================
   !*** for implicit schemes of the local inertial equation
   real(KIND=JPRB),allocatable,target     :: D2RIVOUT_PRE(:,:)     !! river      outflow [m3/s] (prev t-step)
   real(KIND=JPRB),allocatable,target     :: D2RIVDPH_PRE(:,:)     !! river      depth   [m]    (prev t-step)
   real(KIND=JPRB),allocatable,target     :: D2FLDOUT_PRE(:,:)     !! floodplain outflow [m3/s] (prev t-step)
   real(KIND=JPRB),allocatable,target     :: D2FLDSTO_PRE(:,:)     !! floodplain storage [m3]   (prev t-step)

   !================================================
   !*** Groundwater Delay
   real(KIND=JPRD),allocatable,target     :: P2GDWSTO(:,:)         !! ground water storage  [m3]
   real(KIND=JPRB),allocatable,target     :: D2GDWRTN(:,:)         !! Ground water return flow [m3/s]

   !================================================
   !*** These have a different share, not part of the D2PROG array
   real(KIND=JPRB),allocatable,target     :: D1PTHFLW(:,:)         !! flood path outflow [m3/s]
   real(KIND=JPRB),allocatable,target     :: D1PTHFLW_PRE(:,:)     !! flood path outflow [m3/s] (prev t-step)

   !================================================
   !!!*** dam variables
   real(KIND=JPRD),allocatable,target     :: P2DAMSTO(:,:)         !! reservoir storage [m3]
   real(KIND=JPRD),allocatable,target     :: P2DAMINF(:,:)         !! reservoir inflow [m3/s]; discharge before operation

   !================================================
   !!!*** levee variables
   real(KIND=JPRD),allocatable,target     :: P2LEVSTO(:,:)         !! flood storage in protected side (storage betwen river & levee)

END MODULE YOS_CMF_PROG
