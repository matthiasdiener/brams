module sib_vars

  use ModNamelistFile, only: namelistFile

  use ref_sounding, only : maxsndg

  integer :: n_co2 ! from RAMSIN
  real    :: co2_init(maxsndg) ! from RAMSIN


contains
    subroutine StoreNamelistFileAtSib_vars(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    co2_init = oneNamelistFile%co2_init
    n_co2 = oneNamelistFile%n_co2
  end subroutine StoreNamelistFileAtSib_vars

end module sib_vars
