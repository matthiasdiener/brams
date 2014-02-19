!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_carma

  use grid_dims, only: nzpmax ! INTENT(IN)

  implicit none

  type carma_v

     real, pointer :: aot(:,:,:)

  end type carma_v

  type(carma_v), allocatable :: carma(:), carma_m(:)

  !Used in radcomp_carma (radcom) - radriv

!--(DMK-CCATT-INI)----------------------------------------------------------  
  integer, parameter:: ntotal_aer = 2 !ner & rmf 
!--(DMK-CCATT-FIM)----------------------------------------------------------  

  real :: tdec
  real :: sdec
  real :: cdec
  real :: declin
  real :: rvr(nzpmax)
  real :: rtr(nzpmax)
  real :: dn0r(nzpmax)
  real :: pird(nzpmax)
  real :: prd(nzpmax)
  real :: temprd(nzpmax+1)  
  real :: fthrl(nzpmax)
  real :: dzmr(nzpmax)
  real :: dztr(nzpmax)
  real :: fthrs(nzpmax) 
  
  !Used for vetorization
  integer :: p_isize
  integer :: p_jsize

  real, allocatable :: p_surf(:)
  real, allocatable :: p_top(:)
  real, allocatable :: t_surf(:)
  real, allocatable :: tabove_aerad(:)
  real, allocatable :: solnet(:)
  real, allocatable :: xirdown(:)
    
  real, allocatable :: p(:,:)
  real, allocatable :: t(:,:)
  real, allocatable :: t_aerad(:,:)
  real, allocatable :: p_aerad(:,:)
  real, allocatable :: qv_aerad(:,:)
!kmlnew
  real, allocatable :: LWL_aerad(:,:)
  real, allocatable :: IWL_aerad(:,:)
  real, allocatable :: LWP_aerad(:,:)
  real, allocatable :: IWP_aerad(:,:)
  real, allocatable :: xland_aerad(:)
!srf  REAL, ALLOCATABLE, DIMENSION(:,:)     :: RAIN_aerad
  real, allocatable :: RAIN_aerad(:)
!kmlnew   
  real, allocatable :: press(:,:)
  real, allocatable :: dpg(:,:)
  real, allocatable :: tt(:,:)
  real, allocatable :: rhoa(:,:)
  real, allocatable :: rsfx(:,:)
  real, allocatable :: heats_aerad(:,:)
  real, allocatable :: heati_aerad(:,:)
  real, allocatable :: rdh2o(:,:)
  real, allocatable :: ptempg(:,:)
  real, allocatable :: ptempt(:,:)
  real, allocatable :: u1i(:,:)

  real, allocatable :: gc(:,:,:)
  real, allocatable :: pah2o(:,:,:)
  real, allocatable :: paco2(:,:,:)
  real, allocatable :: pao2(:,:,:)
  real, allocatable :: pao3(:,:,:)
  real, allocatable :: taugas(:,:,:)
  real, allocatable :: paray(:,:,:)
  real, allocatable :: tauaer(:,:,:)
  real, allocatable :: taul(:,:,:)
  real, allocatable :: taucld(:,:,:)
  real, allocatable :: wcld(:,:,:)
  real, allocatable :: gcld(:,:,:)
  real, allocatable :: w0(:,:,:)
  real, allocatable :: g0(:,:,:)
!kmlnew
  real, allocatable :: taucldlw(:,:,:)
  real, allocatable :: wolc(:,:,:)
  real, allocatable :: gl(:,:,:)
  real, allocatable :: taucldice(:,:,:)
  real, allocatable :: woice(:,:,:)
  real, allocatable :: gice(:,:,:)

!--(DMK-CCATT-INI)----------------------------------------------------------  
  real, allocatable, dimension(:,:,:,:) :: tauaer_x 
  real, allocatable, dimension(:,:,:,:) :: wol_x   
  real, allocatable, dimension(:,:,:,:) :: gol_x
!--(DMK-CCATT-FIM)----------------------------------------------------------  
    
!kmlnew    
  real, allocatable :: opd(:,:,:)
  real, allocatable :: uopd(:,:,:)
  real, allocatable :: ptemp(:,:,:)
  real, allocatable :: slope(:,:,:)
  real, allocatable :: b1(:,:,:)
  real, allocatable :: b2(:,:,:)
  real, allocatable :: b3(:,:,:)
  real, allocatable :: ak(:,:,:)
  real, allocatable :: gami(:,:,:)
  real, allocatable :: ee1(:,:,:)
  real, allocatable :: el1(:,:,:)
  real, allocatable :: em1(:,:,:)
  real, allocatable :: el2(:,:,:)
  real, allocatable :: em2(:,:,:)
  real, allocatable :: af(:,:,:)
  real, allocatable :: bf(:,:,:)
  real, allocatable :: ef(:,:,:)
  real, allocatable :: cp(:,:,:)
  real, allocatable :: cpb(:,:,:)
  real, allocatable :: cmb(:,:,:)
  real, allocatable :: ck1(:,:,:)
  real, allocatable :: ck2(:,:,:)
  real, allocatable :: fnet(:,:,:)
  real, allocatable :: tmi(:,:,:)
  real, allocatable :: tmid(:,:,:)
  real, allocatable :: tmiu(:,:,:)
  real, allocatable :: direc(:,:,:)
  real, allocatable :: directu(:,:,:)
    
  real, allocatable :: pc(:,:,:,:)    
  real, allocatable :: pc_aerad(:,:,:,:)
  real, allocatable :: caer(:,:,:,:)
  real, allocatable :: y3(:,:,:,:)

  integer, allocatable :: isl_aerad(:)
  integer, allocatable :: lla(:)
  integer, allocatable :: lls(:)
  
  real, parameter :: emisir_aerad = 1.0
  real, parameter :: h2ocol_aerad = 0.01

  integer :: ir_aerad

!--(DMK-CCATT-INI)----------------------------------------------------------  
  !for tuv
  INTEGER, PARAMETER :: na = 4  ! number of aerosols
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: g_tauaer
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: g_wol
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: g_gol
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: g_taucld 
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: g_wcld
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: g_gcld
!--(DMK-CCATT-FIM)----------------------------------------------------------  

contains

  subroutine alloc_carma(car, ng, nmxp, nmyp, nw)

    implicit none
    
    type(carma_v), intent(inout) :: car(:)
    integer, intent(in)          :: ng, nw, nmxp, nmyp
    
    allocate(car(ng)%aot(nmxp,nmyp,nw))
    
  end subroutine alloc_carma

  !---------------------------------------------------------------

  subroutine nullify_carma(car, ng)
    
    implicit none
    
    type(carma_v), intent(inout) :: car(:)
    integer, intent(in)          :: ng
    
    if (associated(car(ng)%aot))  nullify (car(ng)%aot)

  end subroutine nullify_carma

  !---------------------------------------------------------------

  subroutine dealloc_carma(car, ng)

    implicit none

    type(carma_v), intent(inout) :: car(:)
    integer, intent(in)          :: ng

    if (associated(car(ng)%aot)) deallocate (car(ng)%aot)

  end subroutine dealloc_carma

  !-----------------------------------------------------------------

  subroutine zero_carma(car, ng)

    implicit none

    type(carma_v), intent(inout) :: car(:)
    integer, intent(in)          :: ng

    car(ng)%aot(:,:,:) = 0.

  end subroutine zero_carma

  !---------------------------------------------------------------

  subroutine filltab_carma(cv, cvm, ng, imean, n1, n2, n3)
    use var_tables, only: InsertVTab
    use mem_scalar, only: RECYCLE_TRACERS ! INTENT(IN)
    use io_params, only : ipastin         ! INTENT(IN)

    implicit none
    
!--(DMK-LFR NEC-SX6)----------------------------------------------
    include 'i8.h'
!--(DMK-LFR NEC-SX6)----------------------------------------------
    
    ! Arguments:
    integer, intent(in) :: ng, n1, n2, n3, imean
    type(carma_v), intent(in) :: cv, cvm
    ! Local Variables:

!--(DMK-LFR NEC-SX6)----------------------------------------------
!    integer          :: npts
    integer(kind=i8) :: npts
!--(DMK-LFR NEC-SX6)----------------------------------------------

    character(len=7) :: sname
    ! ALF
    character(len=8) :: str_recycle

    ! ALF
    str_recycle = ''
    if (RECYCLE_TRACERS==1 .or. ipastin==1) then
       str_recycle = ':recycle'
    endif

    if (associated(cv%aot)) then
       npts = n1*n2*n3
       write(sname,'(a4)') 'AOT'
       call InsertVTab (cv%aot, cvm%aot, ng, &
            npts, imean,sname//' :7:hist:anal:mpti:mpt3'//trim(str_recycle))
       ! Not necessary MPT1 - Comunication NODE to NODE on DTLONG
       ! Radiation is a Column oriented process
    endif

  end subroutine filltab_carma

  !---------------------------------------------------------------

  subroutine init_carma(ia, iz, ja, jz, m1, m2, m3)

    use mem_globrad, only: ntotal, nlayer, ndbl, ngauss, nrad
    use mem_aerad, only: nz_rad, ngas, nelem, ngroup, nbin
  
    IMPLICIT NONE
    ! Arguments:
    integer,intent(IN) :: ia, iz, ja, jz, m1, m2, m3
    ! Local Variables:
    integer :: iend
  
    iend = (iz-ia+1)*(jz-ja+1)
    !2d
    allocate(p_surf(iend)) 
    allocate(p_top(iend))  
    allocate(t_surf(iend)) 
    allocate(tabove_aerad(iend))  
    allocate(solnet(iend))
    allocate(xirdown(iend))
    !3d  
    allocate(p(iend,m1))      
    allocate(t(iend,m1))      
    allocate(t_aerad(iend,m1)) 
    allocate(p_aerad(iend,m1)) 
    allocate(qv_aerad(iend,m1)) 
    !kmlnew
    allocate(LWL_aerad(iend,m1)) 
    allocate(IWL_aerad(iend,m1)) 
    allocate(LWP_aerad(iend,m1)) 
    allocate(IWP_aerad(iend,m1)) 
    allocate(xland_aerad(iend)) 
    !srf  ALLOCATE(RAIN_aerad(iend,m1))
    allocate(RAIN_aerad(iend))
    !kmlnew

    allocate(press(iend,nlayer))   
    allocate(dpg(iend,nlayer))    
    allocate(tt(iend,nlayer))     
    allocate(rdh2o(iend,nlayer))     
    allocate(rhoa(iend,m1))   
    allocate(rsfx(iend,ntotal))  
    allocate(heats_aerad(iend,nz_rad)) 
    allocate(heati_aerad(iend,nz_rad)) 
    allocate(ptempg(iend,ntotal))  
    allocate(ptempt(iend,ntotal)) 
    allocate(u1i(iend,ntotal)) 
    !4d
    allocate(gc(iend,m1,ngas))     
    allocate(pah2o(iend,ntotal,nlayer)) 
    allocate(paco2(iend,ntotal,nlayer)) 
    allocate(pao2(iend,ntotal,nlayer))  
    allocate(pao3(iend,ntotal,nlayer))  
    allocate(taugas(iend,ntotal,nlayer))
    allocate(paray(iend,ntotal,nlayer)) 
    allocate(tauaer(iend,ntotal,nlayer)) 
    allocate(taul(iend,ntotal,nlayer))  
    allocate(taucld(iend,ntotal,nlayer)) 
    allocate(wcld(iend,ntotal,nlayer))   
    allocate(gcld(iend,ntotal,nlayer)) 
    !kmlnew
    allocate(taucldlw (iend,ntotal,nlayer)) 
    allocate(wolc     (iend,ntotal,nlayer))   
    allocate(gl       (iend,ntotal,nlayer)) 
    allocate(taucldice(iend,ntotal,nlayer)) 
    allocate(woice    (iend,ntotal,nlayer))   
    allocate(gice     (iend,ntotal,nlayer)) 

!--(DMK-CCATT-INI)----------------------------------------------------------  
    allocate(tauaer_x (iend,ntotal,nlayer,ntotal_aer)) 
    allocate(wol_x    (iend,ntotal,nlayer,ntotal_aer))   
    allocate(gol_x    (iend,ntotal,nlayer,ntotal_aer)) 
!--(DMK-CCATT-FIM)----------------------------------------------------------  

    !kmlnew
    allocate(w0(iend,ntotal,nlayer))   
    allocate(g0(iend,ntotal,nlayer))   
    allocate(opd(iend,ntotal,nlayer))  
    allocate(uopd(iend,ntotal,nlayer)) 
    allocate(ptemp(iend,ntotal,nlayer))
    allocate(slope(iend,ntotal,nlayer))  
    allocate(b1(iend,ntotal,nlayer))
    allocate(b2(iend,ntotal,nlayer))
    allocate(b3(iend,ntotal,nlayer))
    allocate(ak(iend,ntotal,nlayer))
    allocate(gami(iend,ntotal,nlayer))
    allocate(ee1(iend,ntotal,nlayer))
    allocate(el1(iend,ntotal,nlayer))
    allocate(em1(iend,ntotal,nlayer))
    allocate(el2(iend,ntotal,nlayer))
    allocate(em2(iend,ntotal,nlayer))
    allocate(af(iend,ntotal,ndbl))
    allocate(bf(iend,ntotal,ndbl))
    allocate(ef(iend,ntotal,ndbl))
    allocate(cp(iend,ntotal,nlayer))
    allocate(cpb(iend,ntotal,nlayer))
    allocate(cmb(iend,ntotal,nlayer))
    allocate(ck1(iend,ntotal,nlayer))
    allocate(ck2(iend,ntotal,nlayer))
    allocate(fnet(iend,ntotal,nlayer))
    allocate(tmi(iend,ntotal,nlayer ))
    allocate(tmid(iend,ntotal,nlayer))
    allocate(tmiu(iend,ntotal,nlayer))
    allocate(direc(iend,ntotal,nlayer))
    allocate(directu(iend,ntotal,nlayer))
    !5d  
    allocate(pc(iend,m1,nbin,nelem))  
    allocate(pc_aerad(iend,nz_rad,nbin,ngroup)) 
    allocate(caer(iend,nlayer,nrad,ngroup))    
    allocate(y3(iend,ntotal,ngauss,nlayer))    
    !2d int
    allocate(isl_aerad(iend)) 
    allocate(lla(iend))
    allocate(lls(iend))

    !Zero all variables
    p_surf=0.0 
    p_top=0.0  
    t_surf=0.0 
    tabove_aerad=0.0  
    lla=0
    lls=0
    solnet=0.0
    xirdown=0.0
    !3d  
    p=0.0      
    t=0.0      
    t_aerad=0.0 
    p_aerad=0.0 
    qv_aerad=0.0 
    !kmlnew
    LWL_aerad=0.0 
    IWL_aerad=0.0 
    LWP_aerad=0.0 
    IWP_aerad=0.0
    xland_aerad=0.0
    RAIN_aerad=0.0 
    !kmlnew
    press=0.0   
    dpg=0.0    
    tt=0.0     
    rdh2o=0.0     
    rhoa=0.0   
    rsfx=0.0  
    heats_aerad=0.0 
    heati_aerad=0.0 
    ptempg=0.0 
    ptempt=0.0 
    u1i=0.0
    !4d
    gc=0.0     
    pah2o=0.0 
    paco2=0.0 
    pao2=0.0  
    pao3=0.0  
    taugas=0.0
    paray=0.0 
    tauaer=0.0 
    taul=0.0  
    taucld=0.0
    !kmlnew
    taucldlw=0.0
    wolc=0.0 
    gl=0.0 
    taucldice=0.0
    woice=0.0	
    gice=0.0 

!--(DMK-CCATT-INI)----------------------------------------------------------  
    tauaer_x=0.0
    wol_x=0.0
    gol_x=0.0 
!--(DMK-CCATT-FIM)----------------------------------------------------------     
    
    !kmlnew  
    wcld=0.0   
    gcld=0.0   
    w0=0.0   
    g0=0.0   
    opd=0.0  
    uopd=0.0 
    ptemp=0.0
    slope=0.0  
    b1=0.0
    b2=0.0
    b3=0.0
    ak=0.0
    gami=0.0
    ee1=0.0
    el1=0.0
    em1=0.0
    el2=0.0
    em2=0.0
    af=0.0
    bf=0.0
    ef=0.0
    cp=0.0
    cpb=0.0
    cmb=0.0
    ck1=0.0
    ck2=0.0
    fnet=0.0
    tmi=0.0
    tmid=0.0
    tmiu=0.0
    direc=0.0
    directu=0.0
    !5d  
    pc=0.0  
    pc_aerad=0.0 
    caer=0.0    
    y3=0.0    
    !2d int
    isl_aerad=0    

  end subroutine init_carma

  !---------------------------------------------------------------

  subroutine end_carma()

    implicit none

    !2d
    deallocate(p_surf) 
    deallocate(p_top)  
    deallocate(t_surf) 
    deallocate(tabove_aerad)  
    deallocate(lla)
    deallocate(lls)
    deallocate(solnet)
    deallocate(xirdown)
    !3d  
    deallocate(p)      
    deallocate(t)      
    deallocate(t_aerad) 
    deallocate(p_aerad) 
    deallocate(qv_aerad) 
    !kmlnew
    deallocate(LWL_aerad) 
    deallocate(IWL_aerad) 
    deallocate(LWP_aerad) 
    deallocate(IWP_aerad) 
    deallocate(xland_aerad)
    deallocate(RAIN_aerad)
    !kmlnew
    deallocate(press)   
    deallocate(dpg)    
    deallocate(tt)     
    deallocate(rdh2o)     
    deallocate(rhoa)   
    deallocate(rsfx)  
    deallocate(heats_aerad) 
    deallocate(heati_aerad) 
    deallocate(ptempg)  
    deallocate(ptempt) 
    deallocate(u1i) 
    !4d
    deallocate(gc)     
    deallocate(pah2o) 
    deallocate(paco2) 
    deallocate(pao2)  
    deallocate(pao3)  
    deallocate(taugas)
    deallocate(paray) 
    deallocate(tauaer) 
    deallocate(taul)  
    deallocate(taucld) 
    deallocate(wcld)   
    deallocate(gcld)
    !kmlnew
    deallocate(taucldlw) 
    deallocate(wolc)   
    deallocate(gl) 
    deallocate(taucldice) 
    deallocate(woice)   
    deallocate(gice) 

!--(DMK-CCATT-INI)----------------------------------------------------------  
    deallocate(tauaer_x) 
    deallocate(wol_x) 
    deallocate(gol_x) 
!--(DMK-CCATT-FIM)----------------------------------------------------------  

    !kmlnew

    deallocate(w0)   
    deallocate(g0)   
    deallocate(opd)  
    deallocate(uopd) 
    deallocate(ptemp)
    deallocate(slope)  
    deallocate(b1)
    deallocate(b2)
    deallocate(b3)
    deallocate(ak)
    deallocate(gami)
    deallocate(ee1)
    deallocate(el1)
    deallocate(em1)
    deallocate(el2)
    deallocate(em2)
    deallocate(af)
    deallocate(bf)
    deallocate(ef)
    deallocate(cp)
    deallocate(cpb)
    deallocate(cmb)
    deallocate(ck1)
    deallocate(ck2)
    deallocate(fnet)
    deallocate(tmi)
    deallocate(tmid)
    deallocate(tmiu)
    deallocate(direc)
    deallocate(directu)
    !5d  
    deallocate(pc)  
    deallocate(pc_aerad) 
    deallocate(caer)    
    deallocate(y3)    
    !2d int
    deallocate(isl_aerad)

  end subroutine end_carma

end module mem_carma
