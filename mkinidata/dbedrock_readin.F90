#include <define.h>

#ifdef USE_DEPTH_TO_BEDROCK

subroutine dbedrock_readin (dir_landdata)

   use precision
   use spmd_task
   use mod_landpatch
   use ncio_vector
   USE GlobalVars, only : nl_soil, dz_soi
   use MOD_TimeInvariants, only : dbedrock, ibedrock

   IMPLICIT NONE

   character(LEN=256), INTENT(in) :: dir_landdata

   ! Local Variables
   character(len=256) :: lndname
   integer  :: ipatch, L, ibd
   real(r8) :: dres

   ! Read bedrock
  
   lndname = trim(dir_landdata)//'/dbedrock_patches.nc'

   call ncio_read_vector (lndname, 'dbedrock_patches', landpatch, dbedrock) 

   if (p_is_worker) then

      do ipatch = 1, numpatch

         L = landpatch%ltyp(ipatch)

         if (L == 0) then
            ibedrock(ipatch) = 0
         else
            
            dbedrock(ipatch) = dbedrock(ipatch) / 100.0 ! from cm to meter
            dbedrock(ipatch) = max(dbedrock(ipatch), dz_soi(1))

            ibd = 1
            dres = dbedrock(ipatch)
            do while (ibd <= nl_soil)
               if (dres > dz_soi(ibd)) then
                  dres = dres - dz_soi(ibd)
                  ibd = ibd + 1
               else
                  exit
               end if
            end do

            ibedrock(ipatch) = ibd

         end if

      end do

   end if

end subroutine dbedrock_readin

#endif
