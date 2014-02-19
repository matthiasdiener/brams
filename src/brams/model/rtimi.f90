!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine tend0()

  use mem_grid, only: ngrid
  use mem_tend, only: tend
  use var_tables, only: num_scalar, scalar_tab
  use node_mod, only: mxp, myp, mzp

  implicit none
  include "i8.h"
  integer :: n
  integer(kind=i8) :: mxyzp

  !     This routine simply sets all tendency arrays to zero.

  !     First u,v tendencies

  mxyzp = mxp*myp*mzp
!!$call azero(mxyzp,tend%ut(1))
!!$call azero(mxyzp,tend%vt(1))
!!$call azero(mxyzp,tend%wt(1))
!!$call azero(mxyzp,tend%pt(1))
  tend%ut = 0.
  tend%vt = 0.
  tend%wt = 0.
  tend%pt = 0.
  
  !     Now sclrr tendencies

  do n = 1,num_scalar(ngrid)
     call azero_l(mxyzp, scalar_tab(n,ngrid)%var_t)
  enddo

end subroutine tend0

!**************************************************************************

subroutine hadvance(iac)

  use mem_grid, only: eps, ngrid, dtlv, icorflg, jdim
  use mem_tend, only: tend
  use mem_basic, only: basic_g
  use mem_scratch, only: scratch
  use node_mod, only: mxp, myp, mzp

  implicit none
  include "i8.h"
  integer :: iac

  integer(kind=i8) :: mxyzp

  !     It is here that the Asselin filter is applied.  For the velocities
  !     and pressure, this must be done in two stages, the first when
  !     IAC=1 and the second when IAC=2.

  mxyzp = mxp * myp * mzp
  eps = .2

  !     For both IAC=1 and IAC=2, call PREDICT for U, V, W, and P.

  call predict(mxyzp,basic_g(ngrid)%uc(1,1,1)   &
       ,basic_g(ngrid)%up(1,1,1),tend%ut(1),scratch%vt3da(1),iac,dtlv)

  if (icorflg .eq. 1 .or. jdim .eq. 1) then
     call predict(mxyzp,basic_g(ngrid)%vc(1,1,1)  &
          ,basic_g(ngrid)%vp(1,1,1),tend%vt(1),scratch%vt3da(1),iac,dtlv)
  endif

  call predict(mxyzp,basic_g(ngrid)%wc(1,1,1),basic_g(ngrid)%wp(1,1,1)  &
       ,tend%wt(1),scratch%vt3da(1),iac,dtlv)
  call predict(mxyzp,basic_g(ngrid)%pc(1,1,1),basic_g(ngrid)%pp(1,1,1)  &
       ,tend%pt(1),scratch%vt3da(1),iac,dtlv)

  return
end subroutine hadvance


!**************************************************************************

subroutine predict(npts,ac,ap,fa,af,iac,dtlp)

  use mem_grid, only: eps, ngbegun, ngrid, nzp, nxp, nyp
  use node_mod, only: nmachs

  implicit none
  include "i8.h"
  integer :: iac   !,npts, m
  integer(kind=i8) :: m
  integer(kind=i8), intent(in) :: npts
  real :: epsu,dtlp
  real, dimension(npts) :: ac,ap,fa,af !(*)

  !     For IAC=3, this routine moves the arrays AC and AP forward by
  !     1 time level by adding in the prescribed tendency. It also
  !     applies the Asselin filter given by:

  !              {AC} = AC + EPS * (AP - 2 * AC + AF)

  !     where AP,AC,AF are the past, current and future time levels of A.
  !     All IAC=1 does is to perform the {AC} calculation without the AF
  !     term present.  IAC=2 completes the calculation of {AC} by adding
  !     the AF term only, and advances AC by filling it with input AP
  !     values which were already updated in ACOUSTC.
  !
  epsu = eps
  if (ngbegun(ngrid) .eq. 0) epsu = 0.5

  if (iac .eq. 1) then
     ac = ac + epsu * (ap -2. * ac)
     return
  elseif (iac .eq. 2) then
     af = ap
     ap = ac + epsu * af
  elseif (iac .eq. 3) then
     af = ap + dtlp * fa
     ap = ac + epsu * (ap - 2. * ac + af)
  endif

  ac = af

  return
end subroutine predict

!**************************************************************************

subroutine predtr()

  use mem_grid, only: ngrid, dtlt
  use var_tables, only: num_scalar, scalar_tab
  use node_mod, only: mxp, myp, mzp

  implicit none
  include "i8.h"
  integer :: n !mxyzp
  integer(kind=i8) :: mxyzp

  !   -  Step thermodynamic variables from  t  to  t+1.
  !   -  Set top, lateral and bottom boundary conditions on some variables
  !        if needed.
  !   -  Call adjustment to assure all positive definite quantities
  !        remain positive.
  !   -  Rediagnose some thermodynamic quantities for use on the small
  !        timestep.

  !     Update the scalars and apply lateral, top, and bottom boundary
  !     conditions.

  mxyzp = mxp * myp * mzp

  do n = 1,num_scalar(ngrid)
     call update_long(mxyzp, scalar_tab(n,ngrid)%var_p,  &
          scalar_tab(n,ngrid)%var_t, dtlt)
  enddo

end subroutine predtr






