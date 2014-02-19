!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_micro

  type micro_vars

     ! Variables to be dimensioned by (nzp,nxp,nyp)
     real, pointer, dimension(:,:,:) :: &
          rcp,rrp,rpp,rsp,rap,rgp,rhp &
          ,ccp,crp,cpp,csp,cap,cgp,chp &
          ,cccnp,cifnp,q2,q6,q7

     ! Variables to be dimensioned by (nnxp,nyp)
     real, pointer, dimension(:,:) :: &
          accpr,accpp,accps,accpa,accpg,accph &
          ,pcprr,pcprp,pcprs,pcpra,pcprg,pcprh &
          ,pcpg,qpcpg,dpcpg 

  end type micro_vars

  type (micro_vars), allocatable :: micro_g(:), microm_g(:)

contains                  

  subroutine alloc_micro(micro,n1,n2,n3,ng)

    use micphys

    implicit none          
    type (micro_vars) :: micro
    integer, intent(in) :: n1,n2,n3,ng

    ! Allocate arrays based on options (if necessary)

    if (level >= 2 ) then
       allocate (micro%rcp(n1,n2,n3))

       !--(DMK-LFR NEC-SX6)----------------------------------------------	         
       micro%rcp = 0.
       !--(DMK-LFR NEC-SX6)----------------------------------------------	  

    endif
    if (level >= 3) then
       if(irain >= 1)  then
          allocate (micro%rrp(n1,n2,n3))
          allocate (micro%accpr(n2,n3))
          allocate (micro%pcprr(n2,n3))
          allocate (micro%q2(n1,n2,n3))

          !--(DMK-LFR NEC-SX6)----------------------------------------------	  
          micro%rrp   = 0.
	  micro%accpr = 0.
	  micro%pcprr = 0.
	  micro%q2    = 0.
   !--(DMK-LFR NEC-SX6)----------------------------------------------	  

       endif
       if(ipris >= 1)  then
          allocate (micro%rpp(n1,n2,n3))
          allocate (micro%accpp(n2,n3))
          allocate (micro%pcprp(n2,n3))

          !--(DMK-LFR NEC-SX6)----------------------------------------------	  
          micro%rpp   = 0.
	  micro%accpp = 0.
	  micro%pcprp = 0.
   !--(DMK-LFR NEC-SX6)----------------------------------------------	  

       endif
       if(isnow >= 1)  then
          allocate (micro%rsp(n1,n2,n3))
          allocate (micro%accps(n2,n3))
          allocate (micro%pcprs(n2,n3))

          !--(DMK-LFR NEC-SX6)----------------------------------------------	  
          micro%rsp   = 0.
	  micro%accps = 0.
	  micro%pcprs = 0.
   !--(DMK-LFR NEC-SX6)----------------------------------------------	  

       endif
       if(iaggr >= 1)  then
          allocate (micro%rap(n1,n2,n3))
          allocate (micro%accpa(n2,n3))
          allocate (micro%pcpra(n2,n3))

          !--(DMK-LFR NEC-SX6)----------------------------------------------	  
          micro%rap   = 0.
	  micro%accpa = 0.
	  micro%pcpra = 0.
   !--(DMK-LFR NEC-SX6)----------------------------------------------	  

       endif
       if(igraup >= 1) then
          allocate (micro%rgp(n1,n2,n3))
          allocate (micro%accpg(n2,n3))
          allocate (micro%pcprg(n2,n3))
          allocate (micro%q6(n1,n2,n3))

          !--(DMK-LFR NEC-SX6)----------------------------------------------	  
	  micro%rgp   = 0.
	  micro%accpg = 0.
	  micro%pcprg = 0.
	  micro%q6    = 0.
   !--(DMK-LFR NEC-SX6)----------------------------------------------

       endif
       if(ihail >= 1)  then
          allocate (micro%rhp(n1,n2,n3))
          allocate (micro%accph(n2,n3))
          allocate (micro%pcprh(n2,n3))
          allocate (micro%q7(n1,n2,n3))

          !--(DMK-LFR NEC-SX6)----------------------------------------------	  
          micro%rhp   = 0.
	  micro%accph = 0.
	  micro%pcprh = 0.
	  micro%q7    = 0.
   !--(DMK-LFR NEC-SX6)----------------------------------------------	  

       endif
       !--(DMK-LFR NEC-SX6)----------------------------------------------
       !       if(icloud == 5) allocate (micro%ccp(n1,n2,n3))
       !       if(irain == 5)  allocate (micro%crp(n1,n2,n3))
       !       if(ipris == 5)  allocate (micro%cpp(n1,n2,n3))
       !       if(isnow == 5)  allocate (micro%csp(n1,n2,n3))
       !       if(iaggr == 5)  allocate (micro%cap(n1,n2,n3))
       !       if(igraup == 5) allocate (micro%cgp(n1,n2,n3))
       !       if(ihail == 5)  allocate (micro%chp(n1,n2,n3))

!--(DMK-CARRIO-INI)--------------------------------------------------
       !Carrio 2012 ---------------------------------
       ! ccp and cccnp are needed for icloud 5, 6, and 7 !!!
       if(icloud >= 5) then 
            allocate (micro%ccp(n1,n2,n3))
	    micro%ccp = 0.
	    allocate (micro%cccnp(n1,n2,n3))
	    micro%cccnp = 0.
       endif     
!--(DMK-CARRIO-OLD)--------------------------------------------------
!       if(icloud == 5) then
!          allocate (micro%ccp(n1,n2,n3))
!          micro%cpp = 0.
!       endif
!--(DMK-CARRIO-FIM)--------------------------------------------------

       if(irain == 5)  then
          allocate (micro%crp(n1,n2,n3))
          micro%crp = 0.
       endif
       if(ipris == 5)  then
          allocate (micro%cpp(n1,n2,n3))
          micro%cpp = 0.
       endif
       if(isnow == 5)  then
          allocate (micro%csp(n1,n2,n3))
          micro%csp = 0.
       endif
       if(iaggr == 5)  then
          allocate (micro%cap(n1,n2,n3))
          micro%cap = 0.
       endif
       if(igraup == 5) then
          allocate (micro%cgp(n1,n2,n3))
          micro%cgp = 0.
       endif
       if(ihail == 5)  then
          allocate (micro%chp(n1,n2,n3))
          micro%chp = 0.
       endif
       !--(DMK-LFR NEC-SX6)----------------------------------------------

       ! *************Prob.in NEC SX-6
       ! * 240 Subscript error  array=cccnp size=1 subscript=4 eln=155
       ! PROG=initlz ELN=155(400940320)
       !IF(icloud == 7) ALLOCATE (micro%cccnp(n1,n2,n3))
       allocate (micro%cccnp(n1,n2,n3))
       ! * 240 Subscript error  array=cifnp size=1 subscript=4 eln=155
       ! PROG=initlz ELN=155(40094081c)
       !IF(ipris == 7)  ALLOCATE (micro%cifnp(n1,n2,n3))
       allocate (micro%cifnp(n1,n2,n3))

       allocate (micro%pcpg(n2,n3))
       allocate (micro%qpcpg(n2,n3))
       allocate (micro%dpcpg(n2,n3))

       !--(DMK-LFR NEC-SX6)----------------------------------------------       
       micro%cccnp = 0.
       micro%cifnp = 0.
       micro%pcpg  = 0.
       micro%qpcpg = 0.
       micro%dpcpg = 0.
       !--(DMK-LFR NEC-SX6)----------------------------------------------

    endif

    return
  end subroutine alloc_micro


  subroutine nullify_micro(micro)

    implicit none
    type (micro_vars) :: micro

    if (associated(micro%rcp))     nullify (micro%rcp)
    if (associated(micro%rrp))     nullify (micro%rrp)
    if (associated(micro%rpp))     nullify (micro%rpp)
    if (associated(micro%rsp))     nullify (micro%rsp)
    if (associated(micro%rap))     nullify (micro%rap)
    if (associated(micro%rgp))     nullify (micro%rgp)
    if (associated(micro%rhp))     nullify (micro%rhp)
    if (associated(micro%ccp))     nullify (micro%ccp)
    if (associated(micro%crp))     nullify (micro%crp)
    if (associated(micro%cpp))     nullify (micro%cpp)
    if (associated(micro%csp))     nullify (micro%csp)
    if (associated(micro%cap))     nullify (micro%cap)
    if (associated(micro%cgp))     nullify (micro%cgp)
    if (associated(micro%chp))     nullify (micro%chp)
    if (associated(micro%cccnp))   nullify (micro%cccnp)
    if (associated(micro%cifnp))   nullify (micro%cifnp)
    if (associated(micro%q2))      nullify (micro%q2)
    if (associated(micro%q6))      nullify (micro%q6)
    if (associated(micro%q7))      nullify (micro%q7)

    if (associated(micro%accpr))   nullify (micro%accpr)
    if (associated(micro%accpp))   nullify (micro%accpp)
    if (associated(micro%accps))   nullify (micro%accps)
    if (associated(micro%accpa))   nullify (micro%accpa)
    if (associated(micro%accpg))   nullify (micro%accpg)
    if (associated(micro%accph))   nullify (micro%accph)
    if (associated(micro%pcprr))   nullify (micro%pcprr)
    if (associated(micro%pcprp))   nullify (micro%pcprp)
    if (associated(micro%pcprs))   nullify (micro%pcprs)
    if (associated(micro%pcpra))   nullify (micro%pcpra)
    if (associated(micro%pcprg))   nullify (micro%pcprg)
    if (associated(micro%pcprh))   nullify (micro%pcprh)
    if (associated(micro%pcpg))    nullify (micro%pcpg)
    if (associated(micro%qpcpg))   nullify (micro%qpcpg)
    if (associated(micro%dpcpg))   nullify (micro%dpcpg)

    return
  end subroutine nullify_micro

  subroutine dealloc_micro(micro)

    implicit none
    type (micro_vars) :: micro

    if (associated(micro%rcp))     deallocate (micro%rcp)
    if (associated(micro%rrp))     deallocate (micro%rrp)
    if (associated(micro%rpp))     deallocate (micro%rpp)
    if (associated(micro%rsp))     deallocate (micro%rsp)
    if (associated(micro%rap))     deallocate (micro%rap)
    if (associated(micro%rgp))     deallocate (micro%rgp)
    if (associated(micro%rhp))     deallocate (micro%rhp)
    if (associated(micro%ccp))     deallocate (micro%ccp)
    if (associated(micro%crp))     deallocate (micro%crp)
    if (associated(micro%cpp))     deallocate (micro%cpp)
    if (associated(micro%csp))     deallocate (micro%csp)
    if (associated(micro%cap))     deallocate (micro%cap)
    if (associated(micro%cgp))     deallocate (micro%cgp)
    if (associated(micro%chp))     deallocate (micro%chp)
    if (associated(micro%cccnp))   deallocate (micro%cccnp)
    if (associated(micro%cifnp))   deallocate (micro%cifnp)
    if (associated(micro%q2))      deallocate (micro%q2)
    if (associated(micro%q6))      deallocate (micro%q6)
    if (associated(micro%q7))      deallocate (micro%q7)

    if (associated(micro%accpr))   deallocate (micro%accpr)
    if (associated(micro%accpp))   deallocate (micro%accpp)
    if (associated(micro%accps))   deallocate (micro%accps)
    if (associated(micro%accpa))   deallocate (micro%accpa)
    if (associated(micro%accpg))   deallocate (micro%accpg)
    if (associated(micro%accph))   deallocate (micro%accph)
    if (associated(micro%pcprr))   deallocate (micro%pcprr)
    if (associated(micro%pcprp))   deallocate (micro%pcprp)
    if (associated(micro%pcprs))   deallocate (micro%pcprs)
    if (associated(micro%pcpra))   deallocate (micro%pcpra)
    if (associated(micro%pcprg))   deallocate (micro%pcprg)
    if (associated(micro%pcprh))   deallocate (micro%pcprh)
    if (associated(micro%pcpg))    deallocate (micro%pcpg)
    if (associated(micro%qpcpg))   deallocate (micro%qpcpg)
    if (associated(micro%dpcpg))   deallocate (micro%dpcpg)

    return
  end subroutine dealloc_micro


  subroutine filltab_micro(micro,microm,imean,n1,n2,n3,ng)
    use var_tables, only: InsertVTab
    implicit none
    include "i8.h"
    type (micro_vars) :: micro,microm
    integer, intent(in) :: imean,n1,n2,n3,ng
    integer(kind=i8) :: npts
    real, pointer :: var,varm

    ! Fill pointers to arrays into variable tables

    npts=n1*n2*n3
    if (associated(micro%rcp))   &
         call InsertVTab (micro%rcp,microm%rcp  &
         ,ng, npts, imean,  &
         'RCP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rrp))   &
         call InsertVTab (micro%rrp,microm%rrp  &
         ,ng, npts, imean,  &
         'RRP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rpp))   &
         call InsertVTab (micro%rpp,microm%rpp  &
         ,ng, npts, imean,  &
         'RPP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rsp))   &
         call InsertVTab (micro%rsp,microm%rsp  &
         ,ng, npts, imean,  &
         'RSP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rap))   &
         call InsertVTab (micro%rap,microm%rap  &
         ,ng, npts, imean,  &
         'RAP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rgp))   &
         call InsertVTab (micro%rgp,microm%rgp  &
         ,ng, npts, imean,  &
         'RGP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%rhp))   &
         call InsertVTab (micro%rhp,microm%rhp  &
         ,ng, npts, imean,  &
         'RHP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%ccp))   &
         call InsertVTab (micro%ccp,microm%ccp  &
         ,ng, npts, imean,  &
         'CCP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%crp))   &
         call InsertVTab (micro%crp,microm%crp  &
         ,ng, npts, imean,  &
         'CRP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cpp))   &
         call InsertVTab (micro%cpp,microm%cpp  &
         ,ng, npts, imean,  &
         'CPP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%csp))   &
         call InsertVTab (micro%csp,microm%csp  &
         ,ng, npts, imean,  &
         'CSP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cap))   &
         call InsertVTab (micro%cap,microm%cap  &
         ,ng, npts, imean,  &
         'CAP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cgp))   &
         call InsertVTab (micro%cgp,microm%cgp  &
         ,ng, npts, imean,  &
         'CGP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%chp))   &
         call InsertVTab (micro%chp,microm%chp  &
         ,ng, npts, imean,  &
         'CHP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cccnp)) &
         call InsertVTab (micro%cccnp,microm%cccnp  &
         ,ng, npts, imean,  &
         'CCCNP :3:hist:anal:mpti:mpt3:mpt1')
    if (associated(micro%cifnp)) &
         call InsertVTab (micro%cifnp,microm%cifnp  &
         ,ng, npts, imean,  &
         'CIFNP :3:hist:anal:mpti:mpt3:mpt1')

    if (associated(micro%q2))   &
         call InsertVTab (micro%q2,microm%q2  &
         ,ng, npts, imean,  &
         'Q2 :3:hist:anal:mpti:mpt3')
    if (associated(micro%q6)) &
         call InsertVTab (micro%q6,microm%q6  &
         ,ng, npts, imean,  &
         'Q6 :3:hist:anal:mpti:mpt3')
    if (associated(micro%q7)) &
         call InsertVTab (micro%q7,microm%q7  &
         ,ng, npts, imean,  &
         'Q7 :3:hist:anal:mpti:mpt3')

    npts=n2*n3
    if (associated(micro%accpr)) &
         call InsertVTab (micro%accpr,microm%accpr  &
         ,ng, npts, imean,  &
         'ACCPR :2:hist:anal:mpti:mpt3')
    if (associated(micro%accpp)) &
         call InsertVTab (micro%accpp,microm%accpp  &
         ,ng, npts, imean,  &
         'ACCPP :2:hist:anal:mpti:mpt3')
    if (associated(micro%accps)) &
         call InsertVTab (micro%accps,microm%accps  &
         ,ng, npts, imean,  &
         'ACCPS :2:hist:anal:mpti:mpt3')
    if (associated(micro%accpa)) &
         call InsertVTab (micro%accpa,microm%accpa  &
         ,ng, npts, imean,  &
         'ACCPA :2:hist:anal:mpti:mpt3')
    if (associated(micro%accpg)) &
         call InsertVTab (micro%accpg,microm%accpg  &
         ,ng, npts, imean,  &
         'ACCPG :2:hist:anal:mpti:mpt3')
    if (associated(micro%accph)) &
         call InsertVTab (micro%accph,microm%accph  &
         ,ng, npts, imean,  &
         'ACCPH :2:hist:anal:mpti:mpt3')
    if (associated(micro%pcprr)) &
         call InsertVTab (micro%pcprr,microm%pcprr  &
         ,ng, npts, imean,  &
         'PCPRR :2:anal:mpt3')
    if (associated(micro%pcprp)) &
         call InsertVTab (micro%pcprp,microm%pcprp  &
         ,ng, npts, imean,  &
         'PCPRP :2:anal:mpt3')
    if (associated(micro%pcprs)) &
         call InsertVTab (micro%pcprs,microm%pcprs  &
         ,ng, npts, imean,  &
         'PCPRS :2:anal:mpt3')
    if (associated(micro%pcpra)) &
         call InsertVTab (micro%pcpra,microm%pcpra  &
         ,ng, npts, imean,  &
         'PCPRA :2:anal:mpt3')
    if (associated(micro%pcprg)) &
         call InsertVTab (micro%pcprg,microm%pcprg  &
         ,ng, npts, imean,  &
         'PCPRG :2:anal:mpt3')
    if (associated(micro%pcprh)) &
         call InsertVTab (micro%pcprh,microm%pcprh  &
         ,ng, npts, imean,  &
         'PCPRH :2:anal:mpt3')
    if (associated(micro%pcpg)) &
         call InsertVTab (micro%pcpg,microm%pcpg  &
         ,ng, npts, imean,  &
         'PCPG :2:hist:mpti:mpt3')
    if (associated(micro%qpcpg)) &
         call InsertVTab (micro%qpcpg,microm%qpcpg  &
         ,ng, npts, imean,  &
         'QPCPG :2:hist:mpti:mpt3')
    if (associated(micro%dpcpg)) &
         call InsertVTab (micro%dpcpg,microm%dpcpg  &
         ,ng, npts, imean,  &
         'DPCPG :2:hist:mpti:mpt3')

    return
  end subroutine filltab_micro

end module mem_micro
