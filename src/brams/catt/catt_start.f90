! Active or not the CATT module
module catt_start

  use ModNamelistFile, only: namelistFile

  implicit none

  integer :: catt ! from RAMSIN

contains
  subroutine StoreNamelistFileAtCatt_start(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    catt = oneNamelistFile%catt
  end subroutine StoreNamelistFileAtCatt_start

end module catt_start
