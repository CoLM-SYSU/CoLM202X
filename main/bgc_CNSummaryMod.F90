module bgc_CNSummaryMod

use precision
implicit none

public CNDriverSummarizeStates
public CNDriverSummarizeFluxes

private soilbiogeochem_carbonstate_summary
private soilbiogeochem_nitrogenstate_summary
private cnveg_carbonstate_summary
private cnveg_nitrogenstate_summary
private soilbiogeochem_carbonflux_summary
private soilbiogeochem_nitrogenflux_summary
private cnveg_carbonflux_summary
private cnveg_nitrogenflux_summary

contains

subroutine CNDriverSummarizeStates()

call soilbiogeochem_carbonstate_summary()
call soilbiogeochem_nitrogenstate_summary()

call cnveg_carbonstate_summary()
call cnveg_nitrogenstate_summary()

end subroutine CNDriverSummarizeStates

subroutine CNDriverSummarizeFluxes(i,ps,pe,nl_soil,dz_soi)

integer, intent(IN) :: i
integer, intent(IN) :: ps
integer, intent(IN) :: pe
integer, intent(IN) :: nl_soil
real(r8),intent(IN) :: dz_soi(1:nl_soil)

call soilbiogeochem_carbonflux_summary()
call soilbiogeochem_nitrogenflux_summary()

call cnveg_carbonflux_summary()
call cnveg_nitrogenflux_summary()

end subroutine CNDriverSummarizeFluxes

subroutine soilbiogeochem_carbonstate_summary()

end subroutine soilbiogeochem_carbonstate_summary

subroutine soilbiogeochem_nitrogenstate_summary()

end subroutine soilbiogeochem_nitrogenstate_summary

subroutine cnveg_carbonstate_summary()

end subroutine cnveg_carbonstate_summary

subroutine cnveg_nitrogenstate_summary()

end subroutine cnveg_nitrogenstate_summary

subroutine soilbiogeochem_carbonflux_summary()

end subroutine soilbiogeochem_carbonflux_summary

subroutine soilbiogeochem_nitrogenflux_summary()

end subroutine soilbiogeochem_nitrogenflux_summary

subroutine cnveg_carbonflux_summary()

end subroutine cnveg_carbonflux_summary

subroutine cnveg_nitrogenflux_summary()

end subroutine cnveg_nitrogenflux_summary

end module bgc_CNSummaryMod


