MODULE YOS_CMF_MAP
!==========================================================
!* PURPOSE: Shared map/topography variables for CaMa-Flood
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
   USE PARKIND1, ONLY: JPIM, JPRM, JPRB,JPRD
   IMPLICIT NONE
   SAVE 
   !================================================
   !*** river network
   integer(KIND=JPIM),allocatable           ::  I2NEXTX(:,:)       !! POINT DOWNSTREAM HORIZONTAL
   integer(KIND=JPIM),allocatable           ::  I2NEXTY(:,:)       !! POINT DOWNSTREAM VERTICAL

   integer(KIND=JPIM),allocatable           ::  I1SEQX(:)          !! 1D SEQUENCE HORIZONTAL
   integer(KIND=JPIM),allocatable           ::  I1SEQY(:)          !! 1D SEQUENCE VERTICAL
   integer(KIND=JPIM),allocatable           ::  I1NEXT(:)          !! 1D DOWNSTREAM
   integer(KIND=JPIM)                       ::  NSEQRIV            !! LENGTH OF 1D SEQUNECE FOR RIVER
   integer(KIND=JPIM)                       ::  NSEQALL            !! LENGTH OF 1D SEQUNECE FOR RIVER AND MOUTH
   integer(KIND=JPIM)                       ::  NSEQMAX            !! MAX OF NSEQALL (PARALLEL)

   integer(KIND=JPIM),allocatable           ::  I2VECTOR(:,:)      !! VECTOR INDEX
   integer(KIND=JPIM),allocatable           ::  I2REGION(:,:)      !! REGION INDEX
   integer(KIND=JPIM)                       ::  REGIONALL          !! REGION TOTAL
   integer(KIND=JPIM)                       ::  REGIONTHIS         !! REGION THIS CPU
   integer(KIND=JPIM)                       ::  MPI_COMM_CAMA      !! MPI COMMUNICATOR

   !================================================
   !*** lat, lon
   real(KIND=JPRB),allocatable              ::  D1LON(:)           !! longitude [degree_east]
   real(KIND=JPRB),allocatable              ::  D1LAT(:)           !! latitude  [degree_north]

   !================================================
   !*** River + Floodplain topography (map)
   real(KIND=JPRB),allocatable              ::  D2GRAREA(:,:)      !! GRID AREA [M2]
   real(KIND=JPRB),allocatable              ::  D2ELEVTN(:,:)      !! ELEVATION [M]
   real(KIND=JPRB),allocatable              ::  D2NXTDST(:,:)      !! DISTANCE TO THE NEXT GRID [M]
   real(KIND=JPRB),allocatable              ::  D2RIVLEN(:,:)      !! RIVER LENGTH [M]
   real(KIND=JPRB),allocatable              ::  D2RIVWTH(:,:)      !! RIVER WIDTH [M]
   real(KIND=JPRB),allocatable              ::  D2RIVMAN(:,:)      !! RIVER MANNING COEFFICIENT
   real(KIND=JPRB),allocatable              ::  D2RIVHGT(:,:)      !! RIVER HEIGHT [M]
   real(KIND=JPRB),allocatable              ::  D2FLDHGT(:,:,:)    !! FLOODPLAIN HEIGHT [M]

   real(KIND=JPRB),allocatable              ::  D2GDWDLY(:,:)      !! Ground water delay
   real(KIND=JPRB),allocatable              ::  D2ELEVSLOPE(:,:)   !! River bed slope
   integer(KIND=JPIM),allocatable           ::  I2MASK(:,:)        !! Mask 

   !================================================
   !*** Floodplain Topography (diagnosed)
   real(KIND=JPRD),allocatable              ::  P2RIVSTOMAX(:,:)   !! maximum river storage [m3]
   real(KIND=JPRD),allocatable              ::  P2FLDSTOMAX(:,:,:) !! MAXIMUM FLOODPLAIN STORAGE [M3]
   real(KIND=JPRB),allocatable              ::  D2RIVELV(:,:)      !! elevation of river bed [m3]

   real(KIND=JPRB),allocatable              ::  D2FLDGRD(:,:,:)    !! FLOODPLAIN GRADIENT
   real(KIND=JPRB)                          ::  DFRCINC            !! FLOODPLAIN FRACTION INCREMENT [-] (1/NLFP)

   !================================================
   !*** Downstream boundary
   real(KIND=JPRB),allocatable              ::  D2MEANSL(:,:)      !! MEAN SEA LEVEL [M]
   real(KIND=JPRB),allocatable              ::  D2SEALEV(:,:)        !! sea level variation [m]
   real(KIND=JPRB),allocatable              ::  D2DWNELV(:,:)        !! downstream boundary elevation [m]

   !================================================
   !*** bifurcation channel
   integer(KIND=JPIM)                       ::  NPTHOUT            !! NUMBER OF FLOODPLAIN PATH
   integer(KIND=JPIM)                       ::  NPTHLEV            !! NUMBER OF FLOODPLAIN PATH LAYER
   integer(KIND=JPIM),allocatable           ::  PTH_UPST(:)        !! FLOOD PATHWAY UPSTREAM   ISEQ
   integer(KIND=JPIM),allocatable           ::  PTH_DOWN(:)        !! FLOOD PATHWAY DOWNSTREAM JSEQ
   real(KIND=JPRB),allocatable              ::  PTH_DST(:)         !! FLOOD PATHWAY DISTANCE [m]
   real(KIND=JPRB),allocatable              ::  PTH_ELV(:,:)       !! FLOOD PATHWAY ELEVATION [m]
   real(KIND=JPRB),allocatable              ::  PTH_WTH(:,:)       !! FLOOD PATHWAY WIDTH [m]
   real(KIND=JPRB),allocatable              ::  PTH_MAN(:)         !! FLOOD PATHWAY Manning

   DATA REGIONALL  /1/
   DATA REGIONTHIS /1/

END MODULE YOS_CMF_MAP
