!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module micphys

  use ModNamelistFile, only: namelistFile

  use grid_dims, only: &
       NZPMAX, maxgrds

  implicit none

  include "files.h"

  !--------------------------------------------------------------------------
  !     The product [(nthz-1)  * dthz ] must equal 25.0.
  !     The product [(nrhhz-1) * drhhz] must equal 0.18.
  !     The product [(ntc-1)   * dtc  ] must equal 20.0.
  !     The product [(ndnc-1)  * ddnc ] must equal 20.e-6.

  integer, parameter :: nthz=26
  integer, parameter :: nrhhz=10
  integer, parameter :: ngam=5000
  integer, parameter :: ninc=201
  integer, parameter :: ndns=15
  integer, parameter :: ntc=21
  integer, parameter :: ndnc=11
  integer, parameter :: nd1cc=30
  integer, parameter :: nd1cr=15
  integer, parameter :: nr2cr=10
  integer, parameter :: nd2cr=30
  integer, parameter :: nr2rr=20
  integer, parameter :: nd2rr=20
  integer, parameter :: ncat=7
  integer, parameter :: nhcat=15
  integer, parameter :: npairc=93
  integer, parameter :: npairr=131
  integer, parameter :: nembc=20

  real, parameter    :: dtc=1.
  real, parameter    :: ddnc=2.e-6
  real, parameter    :: dthz=1.
  real, parameter    :: drhhz=.02
  !--------------------------------------------------------------------------

  integer            :: level ! from RAMSIN
  integer            :: icloud ! from RAMSIN
  integer            :: irain ! from RAMSIN
  integer            :: ipris ! from RAMSIN
  integer            :: isnow ! from RAMSIN
  integer            :: iaggr ! from RAMSIN
  integer            :: igraup ! from RAMSIN
  integer            :: ihail ! from RAMSIN
  integer            :: mkcoltab ! from RAMSIN
  integer            :: jnmb(ncat)
  integer            :: ipairc(nhcat,nhcat)
  integer            :: ipairr(nhcat,nhcat)
  integer            :: jhabtab(31,100,2)
  integer            :: jhcat(nzpmax,ncat)
  integer            :: ict1(nzpmax,ncat)
  integer            :: ict2(nzpmax,ncat)

  real               :: cparm ! from RAMSIN
  real               :: rparm ! from RAMSIN
  real               :: pparm ! from RAMSIN
  real               :: sparm ! from RAMSIN
  real               :: aparm ! from RAMSIN
  real               :: gparm ! from RAMSIN
  real               :: hparm ! from RAMSIN
  real               :: rictmin
  real               :: rictmax
  real               :: dps
  real               :: dps2
  real               :: d1min
  real               :: r2min
  real               :: d2min
  real               :: d1max
  real               :: r2max
  real               :: d2max
  real               :: d1ecc
  real               :: d1ecr
  real               :: r2ecr
  real               :: r2err
  real               :: colf
  real               :: pi4dt
  real               :: sedtime0
  real               :: sedtime1

  real               :: emb0(ncat)
  real               :: emb1(ncat)
  real               :: gnu(ncat) ! from RAMSIN
  real               :: parm(ncat)
  real               :: emb0log(ncat)
  real               :: emb1log(ncat)
  real               :: dict(ncat)


  ! Changing name of variabel SHAPE to VAR_SHAPE
  ! To avoid confusing with SHAPE FUNCTION (intrinsic)
  ! ALF

  real               :: var_shape(nhcat)
  real               :: cfmas(nhcat)
  real               :: pwmas(nhcat)
  real               :: cfvt(nhcat)
  real               :: pwvt(nhcat)
  real               :: dpsmi(nhcat)
  real               :: cfden(nhcat)
  real               :: pwden(nhcat)
  real               :: cfemb0(nhcat)
  real               :: cfen0(nhcat)
  real               :: pwemb0(nhcat)
  real               :: pwen0(nhcat)
  real               :: vtfac(nhcat)
  real               :: frefac1(nhcat)
  real               :: frefac2(nhcat)
  real               :: cfmasft(nhcat)
  real               :: dnfac(nhcat)
  real               :: sipfac(nhcat)
  real               :: pwmasi(nhcat)
  real               :: ch1(nhcat)
  real               :: ch3(nhcat)
  real               :: cdp1(nhcat)
  real               :: pwvtmasi(nhcat)

!--(DMK-CARRIO-INI)-----------------------------------------------------
!change_MP for chem
  real        :: xcoll(nzpmax)
  real        :: rsedim(nzpmax)
!end change_MP
!--(DMK-CARRIO-FIM)-----------------------------------------------------

  real         :: tair(nzpmax)
  real         :: tairc(nzpmax)
  real         :: tairstrc(nzpmax)
  real         :: til(nzpmax)
  real         :: rvstr(nzpmax)
  real         :: press(nzpmax)
  real         :: pitot(nzpmax)
  real         :: rliq(nzpmax)
  real         :: rice(nzpmax)
  real         :: qhydm(nzpmax)
  real         :: rvlsair(nzpmax)
  real         :: rvisair(nzpmax)
  real         :: rvs0(nzpmax)
  real         :: thrmcon(nzpmax)
  real         :: vapdif(nzpmax)
  real         :: dynvisc(nzpmax)
  real         :: rdynvsci(nzpmax)
  real         :: denfac(nzpmax)
  real         :: dn0i(nzpmax)
  real         :: colfacr(nzpmax)
  real         :: colfacr2(nzpmax)
  real         :: colfacc(nzpmax)
  real         :: colfacc2(nzpmax)
  real         :: sumuy(nzpmax)
  real         :: sumuz(nzpmax)
  real         :: sumvr(nzpmax)
  real         :: scrmic1(nzpmax)
  real         :: scrmic2(nzpmax)
  real         :: scrmic3(nzpmax)
  real         :: cccnx(nzpmax)
  real         :: cifnx(nzpmax)

  real         :: rx(nzpmax,ncat)
  real         :: cx(nzpmax,ncat)
  real         :: qr(nzpmax,ncat)
  real         :: qx(nzpmax,ncat)
  real         :: tx(nzpmax,ncat)
  real         :: emb(nzpmax,ncat)
  real         :: vterm(nzpmax,ncat)
  real         :: vap(nzpmax,ncat)
  real         :: ttest(nzpmax,ncat)
  real         :: wct1(nzpmax,ncat)
  real         :: wct2(nzpmax,ncat)
  real         :: sb(nzpmax,ncat)
  real         :: sd(nzpmax,ncat)
  real         :: se(nzpmax,ncat)
  real         :: sf(nzpmax,ncat)
  real         :: sg(nzpmax,ncat)
  real         :: sh(nzpmax,ncat)
  real         :: sm(nzpmax,ncat)
  real         :: ss(nzpmax,ncat)
  real         :: su(nzpmax,ncat)
  real         :: sw(nzpmax,ncat)
  real         :: sy(nzpmax,ncat)
  real         :: sz(nzpmax,ncat)
  real         :: tref(nzpmax,2)
  real         :: rvsref(nzpmax,2)
  real         :: rvsrefp(nzpmax,2)
  real         :: sa(nzpmax,9)
  real         :: eff(nzpmax,10)

  real         :: rxfer(nzpmax,ncat,ncat)
  real         :: qrxfer(nzpmax,ncat,ncat)
  real         :: enxfer(nzpmax,ncat,ncat)
  real         :: dispemb0(nhcat,maxgrds)
  real         :: dispemb1(nhcat,maxgrds)
  real         :: ch2(nhcat,maxgrds)

  real         :: coltabc(nembc,nembc,npairc)
  real         :: coltabr(nembc,nembc,npairr)

  real         :: frachz(nrhhz,nthz)
  real         :: fracc(ndnc,ntc,maxgrds)
  real         :: gamm(4)
  real         :: gamn1(4)
  real         :: gam(ngam,3)
  real         :: gaminc(ngam,2)
  real         :: gamsip13(ngam)
  real         :: gamsip24(ngam)
  real         :: rmlttab(ninc)
  real         :: enmlttab(ninc,nhcat)
  real         :: shedtab(ninc,ndns)
  real         :: sc(2)
  real         :: sk(2)
  real         :: sl(2)
  real         :: sj(7)
  real         :: pcprx(7)
  real         :: accpx(7)

  real         :: r1tabcc(nd1cc)
  real         :: c1tabcc(nd1cc)
  real         :: c2tabcc(nd1cc)
  real         :: r1tabcr(nd1cr,nr2cr,nd2cr)
  real         :: c1tabcr(nd1cr,nr2cr,nd2cr)
  real         :: c2tabrr(nr2rr,nd2rr)

  character(len=f_name_length) :: coltabfn ! from RAMSIN
  ! Modif. by ALF

  !---------------------------------------------------------------------------

contains
  subroutine StoreNamelistFileAtMicphys(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    aparm = oneNamelistFile%aparm
    coltabfn = oneNamelistFile%coltabfn
    cparm = oneNamelistFile%cparm
    gnu = oneNamelistFile%gnu
    gparm = oneNamelistFile%gparm
    hparm = oneNamelistFile%hparm
    iaggr = oneNamelistFile%iaggr
    icloud = oneNamelistFile%icloud
    igraup = oneNamelistFile%igraup
    ihail = oneNamelistFile%ihail
    ipris = oneNamelistFile%ipris
    irain = oneNamelistFile%irain
    isnow = oneNamelistFile%isnow
    level = oneNamelistFile%level
    mkcoltab = oneNamelistFile%mkcoltab
    pparm = oneNamelistFile%pparm
    rparm = oneNamelistFile%rparm
    sparm = oneNamelistFile%sparm
  end subroutine StoreNamelistFileAtMicphys
end Module micphys
