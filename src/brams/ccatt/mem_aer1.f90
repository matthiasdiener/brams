!###########################################################################
!  B - Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_aer1
  use grid_dims, only : maxgrds,maxsclr
  use aer1_list, only : nspecies_aer=> nspecies,nmodes
  use mem_chem1, only : chem_assim,RECYCLE_TRACERS

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  use ModNamelistFile, only: namelistFile

  include "i8.h"
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

  type aer1_vars   
!--- All families
     real, pointer, dimension(:,:,:)  :: sc_p,sc_pp,sc_pf,sc_src
     real, pointer, dimension(:,:  )  :: sc_dd,sc_wd
     real, pointer, dimension(:    )  :: sc_t
!-----------
  end type aer1_vars
  type (aer1_vars)    , allocatable :: aer1_g(:,:,:), aer1m_g(:,:,:)


!  integer, parameter :: nsrc=3  !number_sources
  
 ! type aer1_src_vars   
 !    real, pointer, dimension(:,:,:)  :: sc_src
 ! end type aer1_src_vars

  !type (aer1_src_vars), allocatable :: aer1_src_g(:,:,:),aer1m_src_g(:,:,:)

  !- dimension of sources arrays (=1 for 2dim, =m1 for 3dim)
  integer :: aer1_src_z_dim_g(nspecies_aer,maxgrds)



!  character(LEN=20),dimension(nsrc),parameter :: src_name= &
!! '12345678901234567890'
!(/                      &
!  'antro               '&    
!, 'bburn               '&
!, 'bioge               '/)!

!  integer     :: &
!   antro   = 01, & ! anthropogenic sources
!   bburn   = 02, & ! biomass burning sources 
!   bioge   = 03    ! biogenic sources ! must be equal to "nsrc"
  
    integer :: AEROSOL
    integer :: ind_mode_sedim(maxsclr+nspecies_aer*nmodes)
    integer :: ind_aer_sedim (maxsclr+nspecies_aer*nmodes)
    integer :: num_scalar_aer_1st
contains
  !---------------------------------------------------------------

  subroutine alloc_aer1(aer1,nvert_src,n1,n2,n3,nmodes,nspecies)

    use aer1_list, only : spc_alloc,spc_name, src, ddp, wdp, fdda, offline, on ,off &
                          ,mode_alloc,    nucle,      accum, coarse , transport
    !use mem_chem1, only : chem_assim,RECYCLE_TRACERS
    implicit none

    integer,intent(in) :: n1,n2,n3,nspecies,nmodes
    integer :: ispc,isrc,imode
    
    type (aer1_vars)    ,dimension(     nmodes,nspecies) :: aer1
!    type (aer1_src_vars),dimension(nsrc,nmodes,nspecies) :: aer1_src
    integer,dimension(nspecies)    :: nvert_src
    
    !print*,'----------------------------------------------------------------'
    !print*,' memory allocation for aerosol species:'
    
    do ispc=1,nspecies
    
      do imode=1,nmodes
      
	 !- for aerosols all allocated must be transported
	 if(spc_alloc(transport,imode,ispc) /= 1 .and. mode_alloc(imode,ispc) == 1) then
	   print*,'aerosol=',ispc,'mode=', imode
	   print*,'all allocated must be transported, change aer1-list'	   
	   print*,'transport =',spc_alloc(transport,imode,ispc)
	   print*,'allocation=',mode_alloc(imode,ispc)
	   stop 3334
	 endif
	  
	 !print*,'spc=',spc_name(ispc),'mode=',imode
	 
	 !1st test: if the mode does not exist, cycle
	 if(mode_alloc(imode,ispc) /= 1 ) cycle
         
         !- allocate memory for the past time tracer mixing ratio
         allocate (aer1(imode,ispc)%sc_p  (n1,n2,n3)) ;aer1(imode,ispc)%sc_p = 0.
         

         !- allocate memory for the dry deposition of the tracer mixing ratio
         if(spc_alloc(ddp,imode,ispc) == on) then 
               allocate (aer1(imode,ispc)%sc_dd    (n2,n3)) ;aer1(imode,ispc)%sc_dd= 0.
         endif

         !- allocate memory for the wet deposition of the tracer mixing ratio
         if(spc_alloc(wdp,imode,ispc) == on) then
               allocate (aer1(imode,ispc)%sc_wd    (n2,n3)) ;aer1(imode,ispc)%sc_wd= 0.
         endif
         
         !- allocate memory for the past/future time tracer mixing ratio 
         !- (4D data assimilation)
         if(chem_assim == on) then
            if(spc_alloc(fdda,imode,ispc) == on) then 
               allocate (aer1(imode,ispc)%sc_pp (n1,n2,n3),aer1(imode,ispc)%sc_pf (n1,n2,n3))
               aer1(imode,ispc)%sc_pp =0.; aer1(imode,ispc)%sc_pf = 0.
            endif
         endif

         !- allocate memory for the sources (3d=bburn and 2d= urban/bioge/marin/sdust)
         if(spc_alloc(src,imode,ispc) == on) then 
        	 allocate (aer1(imode,ispc)%sc_src  (nvert_src(ispc),n2,n3))
        	           aer1(imode,ispc)%sc_src = 0.
         endif

     enddo
    enddo
    return
  end subroutine alloc_aer1

  !--------------------------------------------------------------------------

  subroutine dealloc_aer1(aer1,nmodes,nspecies)

   implicit none

   integer,intent(in) :: nspecies,nmodes
   type (aer1_vars)    ,dimension(nmodes,nspecies) :: aer1
   integer :: ispc,imode
  
    !  Deallocate arrays
    do ispc=1,nspecies
      do imode=1,nmodes

        if (associated(aer1(imode,ispc)%sc_p ) )    deallocate(aer1(imode,ispc)%sc_p  )
        if (associated(aer1(imode,ispc)%sc_src))    deallocate(aer1(imode,ispc)%sc_src)
        if (associated(aer1(imode,ispc)%sc_dd) )    deallocate(aer1(imode,ispc)%sc_dd )
        if (associated(aer1(imode,ispc)%sc_wd) )    deallocate(aer1(imode,ispc)%sc_wd )
        if (associated(aer1(imode,ispc)%sc_pp) )    deallocate(aer1(imode,ispc)%sc_pp )
        if (associated(aer1(imode,ispc)%sc_pf) )    deallocate(aer1(imode,ispc)%sc_pf )

       enddo
    enddo

    return
  end subroutine dealloc_aer1

  !---------------------------------------------------------------
  !---------------------------------------------------------------

  subroutine nullify_aer1(aer1,nmodes,nspecies)

    implicit none

    integer,intent(in) ::nspecies,nmodes
    type (aer1_vars),dimension(nmodes,nspecies) :: aer1
    integer :: ispc,imode

    do ispc=1,nspecies
      do imode=1,nmodes

         if (associated(aer1(imode,ispc)%sc_p ) )    nullify (aer1(imode,ispc)%sc_p  )
         if (associated(aer1(imode,ispc)%sc_src))    nullify (aer1(imode,ispc)%sc_src)
         if (associated(aer1(imode,ispc)%sc_dd) )    nullify (aer1(imode,ispc)%sc_dd )
         if (associated(aer1(imode,ispc)%sc_wd) )    nullify (aer1(imode,ispc)%sc_wd )
         if (associated(aer1(imode,ispc)%sc_pp) )    nullify (aer1(imode,ispc)%sc_pp )
         if (associated(aer1(imode,ispc)%sc_pf) )    nullify (aer1(imode,ispc)%sc_pf )
        
      enddo
    enddo

    return
  end subroutine nullify_aer1

  !---------------------------------------------------------------

  subroutine filltab_aer1(aer1,aer1m,imean,nvert_src,n1,n2,n3,nmodes,nspecies,ng)

!    use var_tables
    use aer1_list, only : spc_alloc,spc_name, src, ddp, wdp, fdda, offline, on ,off &
                         ,mode_alloc, mode_name,   nucle,      accum, coarse 

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    use var_tables, only: InsertVTab
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

    implicit none

    integer, intent(in) :: imean,n1,n2,n3,nmodes,nspecies,ng
    type (aer1_vars)    ,dimension(nmodes,nspecies) :: aer1,aer1m
    integer,dimension(nspecies)    :: nvert_src

    integer :: ispc,imode  

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    integer(kind=i8) :: npts
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!   integer :: npts
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------
    
    character(len=8) :: str_recycle
    character(len=1) :: str_src_dim
   
    str_recycle = ''; str_src_dim = ''
    if (RECYCLE_TRACERS == 1) then
       str_recycle = ':recycle'
    endif

    !- Fill pointers to arrays into variable tables
    do ispc=1,nspecies
      do imode=1,nmodes

  	 if (associated(aer1(imode,ispc)%sc_p)) then
!---- tracer mixing ratio (dimension 3d)
  	   npts = n1 * n2 * n3

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  	   call InsertVTab(aer1(imode,ispc)%sc_p,aer1m(imode,ispc)%sc_p,  &
  		ng, npts, imean,                                          &
                trim(spc_name(ispc))//trim(mode_name(imode))              &
		//'P :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle))

!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!  	   call vtables2 (aer1(imode,ispc)%sc_p(1,1,1),aer1m(imode,ispc)%sc_p(1,1,1)  &
!  		,ng, npts, imean, trim(spc_name(ispc))//trim(mode_name(imode)) &
!		//'P :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle))
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------


!---- sources (3 and 2 dimension)
  	  !- old way
  	  ! npts=n1*n2*n3
  	  ! if(spc_alloc(1,ispc) == 1) &
  	  ! call vtables2 (aer1(imode,ispc)%sc_s(1,1,1),aer1m(imode,ispc)%sc_s(1,1,1)  &
  	  !	 ,ng, npts, imean, trim(spc_name(ispc))//'S :3:hist:anal:mpti:mpt3:mpt1')
  	  !- new way
  	   if(spc_alloc(src,imode,ispc) == on) then
  		
		
  		    npts= nvert_src(ispc)  * n2 * n3
  		    if(nvert_src(ispc) ==  1) str_src_dim = '2'  ! for 2d sources
  		    if(nvert_src(ispc) == n1) str_src_dim = '3'  ! for 3d sources

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  		    call InsertVTab(aer1(imode,ispc)%sc_src,                       &
  				    aer1m(imode,ispc)%sc_src,                      &
  				    ng, npts, imean,                               &
                                    trim(spc_name(ispc))//trim(mode_name(imode))// &
                                    '_SRC :'//trim(str_src_dim)//':hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!  		    call vtables2 ( aer1(imode,ispc)%sc_src(1,1,1)  &
!  				  ,aer1m(imode,ispc)%sc_src(1,1,1)  &
!  				  ,ng, npts, imean, trim(spc_name(ispc))//trim(mode_name(imode))//&
!   				  '_SRC :'//trim(str_src_dim)//':hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

  		!print*,'src alloc=', trim(spc_name(ispc))//'_'//trim(src_name(isrc)),' npts=',npts
  	   endif
  	   
!---- dry and wet deposition (dimension 2d)
  	   npts = n2 * n3
  	   if(spc_alloc(ddp,imode,ispc) == on) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  	   call InsertVTab(aer1(imode,ispc)%sc_dd,aer1m(imode,ispc)%sc_dd, &
  		           ng, npts, imean,                                &
                           trim(spc_name(ispc))//trim(mode_name(imode))//'DD :2:hist:anal:mpti:mpt3')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!  	   call vtables2 (aer1(imode,ispc)%sc_dd(1,1),aer1m(imode,ispc)%sc_dd (1,1)  &
!  		,ng, npts, imean, trim(spc_name(ispc))//trim(mode_name(imode))//'DD :2:hist:anal:mpti:mpt3')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

  	   npts = n2 * n3
  	   if(spc_alloc(wdp,imode,ispc) == on) &

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  	   call InsertVTab(aer1(imode,ispc)%sc_wd,aer1m(imode,ispc)%sc_wd,  &
  	                   ng, npts, imean,                                 &
                           trim(spc_name(ispc))//trim(mode_name(imode))//'WD :2:hist:anal:mpti:mpt3')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!  	   call vtables2 (aer1(imode,ispc)%sc_wd(1,1),aer1m(imode,ispc)%sc_wd(1,1)  &
!  	       ,ng, npts, imean, trim(spc_name(ispc))//trim(mode_name(imode))//'WD :2:hist:anal:mpti:mpt3')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------


!----  data assimilation (dimension 3d)
  	   if(chem_assim == on) then
  	      npts = n1 * n2 * n3
  	      if(spc_alloc(fdda,imode,ispc) == on) then

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  		 call InsertVTab(aer1(imode,ispc)%sc_pp,aer1m(imode,ispc)%sc_pp, &
  		                 ng, npts, imean,                                &
                                 trim(spc_name(ispc))//trim(mode_name(imode))//'PP :3:mpti')
  		 call InsertVTab(aer1(imode,ispc)%sc_pf,aer1m(imode,ispc)%sc_pf, &
  		                 ng, npts, imean,                                &
                                 trim(spc_name(ispc))//trim(mode_name(imode))//'PF :3:mpti')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!  		 call vtables2 (aer1(imode,ispc)%sc_pp(1,1,1),aer1m(imode,ispc)%sc_pp(1,1,1)  &
!  		  ,ng, npts, imean, trim(spc_name(ispc))//trim(mode_name(imode))//'PP :3:mpti')
!  		 call vtables2 (aer1(imode,ispc)%sc_pf(1,1,1),aer1m(imode,ispc)%sc_pf(1,1,1)  &
!  		  ,ng, npts, imean, trim(spc_name(ispc))//trim(mode_name(imode))//'PF :3:mpti')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

  	      endif
  	    endif

  	 endif

      enddo
    enddo
  end subroutine filltab_aer1


  !---------------------------------------------------------------

  subroutine alloc_tend_aer1(nmzp,nmxp,nmyp,ngrs,nmodes,nspecies,proc_type)
   use  aer1_list, only :  mode_alloc
   implicit none
   integer,intent(in)                   :: ngrs,proc_type,nspecies,nmodes
   integer,intent(in), dimension (ngrs) :: nmzp,nmxp,nmyp
   integer :: ng,ntpts,ispc,imode

!         Find the maximum number of grid points needed for any grid.

   if(proc_type==1) then
      ntpts=1
   else
      ntpts=0
      do ng=1,ngrs
         ntpts=max( nmxp(ng)*nmyp(ng)*nmzp(ng),ntpts )
      enddo
   endif

!!!!!  WE ARE ONLY CHECKING GRID 1 !!!!!!!!!
!!!!!    All grids must have same scalars defined !!!!!!!
    do ispc=1,nspecies
      do imode=1,nmodes

	 !1st test: if the mode does not exist, cycle
	 if(mode_alloc(imode,ispc) /= 1 ) cycle
         
	 if (associated(aer1_g(imode,ispc,1)%sc_p)) allocate (aer1_g(imode,ispc,1)%sc_t(ntpts))
         
         do ng=2,ngrs
            aer1_g(imode,ispc,ng)%sc_t => aer1_g(imode,ispc,1)%sc_t
         enddo
            
       enddo
    enddo

  end subroutine alloc_tend_aer1
  !---------------------------------------------------------------

  subroutine nullify_tend_aer1(nmodes,nspecies)

    implicit none
    integer,intent(in) :: nspecies,nmodes
    integer ::ispc,imode

    do ispc=1,nspecies
      do imode=1,nmodes
        if (associated(aer1_g(imode,ispc,1)%sc_t)) nullify (aer1_g(imode,ispc,1)%sc_t)
      enddo
    enddo


  end subroutine nullify_tend_aer1

  !---------------------------------------------------------------
  subroutine dealloc_tend_aer1(nmodes,nspecies)
    implicit none
    integer,intent(in) ::  nspecies,nmodes
    integer ::ispc,imode

    do ispc=1,nspecies
      do imode=1,nmodes
       if (associated(aer1_g(imode,ispc,1)%sc_t)) deallocate (aer1_g(imode,ispc,1)%sc_t)
      enddo
    enddo

  end  subroutine dealloc_tend_aer1
   
  !---------------------------------------------------------------

  subroutine filltab_tend_aer1(nmodes,nspecies,ng)
   use var_tables, only : num_scalar
   use aer1_list, only:spc_name,mode_name
   use mem_chem1, only: nspecies_transported ! this is first calculated at chemistry 
                                             ! "filltab_tend_chem1" routine
   implicit none

    integer,intent(in) :: nmodes,nspecies,ng
    integer ::ispc,imode
    integer :: elements
    
    num_scalar_aer_1st = 0

    do ispc=1,nspecies
      do imode=1,nmodes

! Fill pointers to scalar arrays into scalar tables

         if (associated(aer1_g(imode,ispc,ng)%sc_t)) then
           call vtables_scalar (aer1_g(imode,ispc,ng)%sc_p(1,1,1),aer1_g(imode,ispc,ng)%sc_t(1),&
        		      ng,trim(spc_name(ispc))//trim(mode_name(imode))//'P')
           elements = size(aer1_g(imode,ispc,ng)%sc_t)
           call vtables_scalar_new (aer1_g(imode,ispc,ng)%sc_p(1,1,1),aer1_g(imode,ispc,ng)%sc_t(1),&
        		      ng,trim(spc_name(ispc))//trim(mode_name(imode))//'P',elements)
	   !- total number of transported species (CHEM + AER)
	   nspecies_transported = nspecies_transported + 1		   

           !-save the aerosol identity for sedimentation/advection routine
           !ind_mode_sedim(num_scalar(ng)) = imode
           !ind_aer_sedim (num_scalar(ng)) = ispc
	   if(num_scalar_aer_1st == 0) num_scalar_aer_1st = num_scalar(ng)  



         endif  	      

      enddo
    enddo

  end subroutine filltab_tend_aer1
  

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  subroutine StoreNamelistFileAtMem_aer1(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile
    aerosol = oneNamelistFile%aerosol
  end subroutine StoreNamelistFileAtMem_aer1
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

end module mem_aer1

!--------------------------------------------------------------------------

subroutine define_aer1_src_zdim(aer1_src_z_dim,n1)

  use aer1_list, only: nspecies,sdust,bburn,urban,bioge,marin,v_ash
  
  implicit none
  integer,intent(in) :: n1
  
  integer,dimension(nspecies)    :: aer1_src_z_dim

   !- determination of the dimension of Z-dir of source field array
   aer1_src_z_dim(sdust) = n1	     ! 3d 
   aer1_src_z_dim(bburn) = n1	     ! 3d
   aer1_src_z_dim(urban) = n1	     ! 2d
   aer1_src_z_dim(bioge) = n1	     ! 3d
   aer1_src_z_dim(marin) = n1	     ! 3d
   aer1_src_z_dim(v_ash) = n1	     ! 3d
  return
end subroutine define_aer1_src_zdim

  !--------------------------------------------------------------------------
