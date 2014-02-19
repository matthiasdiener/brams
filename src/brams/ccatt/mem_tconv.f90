module mem_tconv
  use chem1_list, only: nspecies_chem=> nspecies

  use aer1_list, only:  nspecies_aer=> nspecies,nmodes

  implicit none
  integer, parameter :: ntotal=nspecies_chem+nspecies_aer*nmodes
  integer :: nspecies
  integer :: nchem_a,nchem_z
  integer, dimension(ntotal) :: ind_chem
  integer :: naer_a,naer_z
  integer, dimension(ntotal)   :: ind_aer
  integer, dimension(ntotal)   :: ind_mode
 
  logical           :: trans_conv_alloc = .false.
  real, allocatable ::&
       se     (:,:),    & ! environment scalar profile Z levels
       se_cup (:,:),    & ! environment scalar profile Z_cup levels
       sc_up  (:,:),    & ! updraft   gas-phase  scalar profile
       sc_dn  (:,:),    & ! downdraft gas-phase  scalar profile
       stcum1d(:,:),    & ! 1d convective tendency
       sc_up_c(:,:),    & ! updraft   aqueous-phase scalar profile
       sc_dn_c(:,:),    & ! DOwndraft aqueous-phase scalar profile
       pw_up  (:,:),    & ! updraft precitable gas/aer
       pw_dn  (:,:),	& ! downdraft precitable gas/aer
       henry_coef(:,:), & ! Henry's constant 
       dn01d  (:)       ! 1d air density

contains

  subroutine alloc_trans_conv(NSPECIES_TRANSPORTED,mgmzp)
    
    use chem1_list, only: nspecies_chem=> nspecies &
                         ,spc_alloc_chem=> spc_alloc &
		         ,transport,on,off,spc_name
    use aer1_list, only:  nspecies_aer=> nspecies &
                         ,spc_alloc_aer=> spc_alloc &
			 ,mode_alloc, nucle, accum, coarse 
  USE memMatrix, ONLY:  &
                        aerosol !use aerosol Matrix?
    
    implicit none
    integer, intent (IN)::NSPECIES_TRANSPORTED,mgmzp
    integer :: ispc,nspecies,imode


    if(trans_conv_alloc) then
       print *,'ERROR: trans_conv already allocated'
       print *,'Routine: trans_conv File: trans_conv.f90'
       stop
    end if

    nspecies = 0
    
    ! chemistry section : mapping

    nchem_a  = 1
    nchem_z  = 0
    do ispc = 1,nspecies_chem

     if(spc_alloc_chem(transport,ispc) == off) cycle
     nchem_z = nchem_z + 1
     ind_chem(nchem_z) = ispc
     !if(spc_name(ispc)=='CO') print*, 'CO INDEX= ',nchem_z, ind_chem(nchem_z) 
    enddo	
    nspecies = nchem_z ! number of chemical species with conv transport

    ! aerosol section : mapping  
    naer_a = 1
    naer_z = 0
    if(AEROSOL == on) then

      naer_a = nchem_z + 1
      naer_z = nchem_z 
      
      do ispc = 1,nspecies_aer
         do imode=1,nmodes
	 
           if(mode_alloc   (imode,ispc          ) == on .and. & 
              spc_alloc_aer(transport,imode,ispc) == on) then 
	      
	      naer_z = naer_z + 1
	      ind_aer (naer_z) = ispc
	      ind_mode(naer_z) = imode
	       
            endif
         enddo
      enddo 
     
      !- total number of species (chem + aer) to be transported 
      nspecies  =  naer_z
    endif
    
    !print*,'xx=',nspecies,NSPECIES_TRANSPORTED,naer_z
    if(nspecies /=NSPECIES_TRANSPORTED) stop 'nspecies must not be different of nspecies_transported' 
    if(nspecies == 0) return

    allocate( se    (nspecies,mgmzp), &
              se_cup(nspecies,mgmzp), & 
	      sc_up (nspecies,mgmzp), &
	      sc_dn (nspecies,mgmzp), &
             stcum1d(nspecies,mgmzp), &
             sc_up_c(nspecies,mgmzp), &
	     sc_dn_c(nspecies,mgmzp), &
	     pw_up  (nspecies,mgmzp), &
	     pw_dn  (nspecies,mgmzp), &
             henry_coef(nspecies,mgmzp))	    
	     

    allocate(dn01d(mgmzp) )

    trans_conv_alloc=.true.

  end subroutine alloc_trans_conv

  subroutine zero_tconv()

    implicit none

!!$    print *,'LFR->Zero tconv !'
    se=.0
    se_cup=.0
    sc_up=.0
    sc_dn=.0
    stcum1d=.0
    dn01d=.0
    sc_up_c=.0
    sc_dn_c=.0
    henry_coef=.0
    pw_up=.0
    pw_dn=0.

  end subroutine zero_tconv

end module mem_tconv
