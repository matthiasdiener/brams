!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


!*** read_nl: master reads namelist file and propagate to remaining processes ***


subroutine read_nl(file)

  ! read_nl:
  !   reads namelist and checks consistency of a few options

  use node_mod, only: &
       mchnum,        &
       master_num

  use mem_grid, only: &
       timeunit,      &
       timmax,        &
       dtlong                 !DSM

  use io_params, only: &
       timstr

  implicit none

  include "files.h"

!--(DMK-LFR NEC-SX6)----------------------------------------------
!  character(len=f_name_length), intent(in) :: file
  character(len=*), intent(in) :: file
!--(DMK-LFR NEC-SX6)----------------------------------------------

  character(len=*), parameter :: h="**(read_nl)**"

  real :: tfact

  ! i/o process reads namelist,
  ! estabilishes default options and checks option consistency

!--(DMK-LFR NEC-SX6)----------------------------------------------
!  call ReadNamelist(file(1:len_trim(file)))
  call ReadNamelist(file)
!--(DMK-LFR NEC-SX6)----------------------------------------------

  ! broadcast namelist data to all processes

  ! change some input time specifications into seconds

  if(timeunit == 'd'.or.timeunit == 'D') tfact=86400.
  if(timeunit == 'h'.or.timeunit == 'H') tfact=3600.
  if(timeunit == 'm'.or.timeunit == 'M') tfact=60.
  if(timeunit == 's'.or.timeunit == 'S') tfact=1.

  timmax=timmax*tfact
  timstr=timstr*tfact
end subroutine read_nl



subroutine ReadNamelist(fileName)

  ! ReadNamelist:
  !    open, reads and close namelist file
  !    implements defaults for namelist variables
  !    check input options consistency

  use io_params, only: frqboth, &
       afilout,                 &
       avgtim,                  &
       frqanl,                  &
       frqhis,                  &
       frqlite,                 &
       frqmean,                 &
       frqprt,                  &
       hfilin,                  &
       hfilout,                 &
       iclobber,                &
       ihistdel,                &
       initfld,                 &
       prtcputime,              &
       ioutput,                 &
       ipastin,                 &
       iplfld,                  &
       isbval,                  &
       isoilflg,                &
       isoilfn,                 &
       isstflg,                 &
       isstfn,                  &
       itopsflg,                &
       itoptflg,                &
       itoptfn,                 &
       iupdndvi,                &
       iupdsst,                 &
       ivegtflg,                &
       ivegtfn,                 &
       ixsctn,                  &
       iz0flg,                  &
       kwrite,                  &
       lite_vars,               &
       ndviflg,                 &
       ndvifn,                  &
       ndvifpfx,                &
       nlite_vars,              &
       nofilflg,                &
       nplt,                    &
       pastfn,                  &
       sfcfiles,                &
       sstfpfx,                 &
       timstr,                  &
       topfiles,                &
       toptenh,                 &
       toptwvl,                 &
       xlite,                   &
       ylite,                   &
       z0fact,                  &
       z0max,                   &
       zlite,                   &
                                ! TEB
       ifusflg,                 &
       ifusfn,                  &
       fusfiles
  use isan_coms, only: gobrad, &
       gobsep, &
       gridwt, &
       guess1st, &
       hybbot, &
       hybtop, &
       i1st_flg, &
       iapr, &
       iarawi, &
       iasrfce, &
       igridfl, &
       iobswin, &
       ioflgisz, &
       ioflgvar, &
       isan_inc, &
       isfc_flg, &
       iszstage, &
       iupa_flg, &
       ivrstage, &
       levth, &
       maxsfc, &
       maxsta, &
       nfeedvar, &
       nigrids, &
       nisn, &
       notid, &
       notsta, &
       respon, &
       sfcinf, &
       sigzwt, &
       stasep, &
       swvlnth, &
       topsigz, &
       varpfx, &
       wvlnth
  use mem_cuparm, only: confrq, &
       cu_prefix, &
       cu_tel, &
       cu_til, &
       if_cuinv, &
       nnqparm, &
       tcu_beg, &
       tcu_end, &
       tnudcu, &
       wcldbs, &
       wt_cu_grid
  use cuparm_grell3, only: &
      g3d_spread, & 
      g3d_smoothh, &
      g3d_smoothv
  use mem_globrad, only: raddatfn !, & !CARMA
!!$       rdatfnum                      !CARMA
  use mem_grell_param, only: closure_type
  use mem_grid, only: centlat, &
       centlon, &
       cphas, &
       deltax, &
       deltay, &
       deltaz, &
       distim, &
       dtlong, &
       dzmax, &
       dzrat, &
       expnme, &
       gridu, &
       gridv, &
       ibnd, &
       icorflg, &
       idate1, &
       ideltat, &
       if_adap, &
       ihtran, &
       imonth1, &
       initial, &
       itime1, &
       iyear1, &
       jbnd, &
       lsflg, &
       nacoust, &
       naddsc, &
       nestz1, &
       nestz2, &
       nfpt, &
       ngrids, &
       ninest, &
       njnest, &
       nknest, &
       nndtrat, &
       nnstbot, &
       nnsttop, &
       nnxp, &
       nnyp, &
       nnzp, &
       npatch, &
       nstratx, &
       nstraty, &
       nstratz1, &
       nstratz2, &
       nxtnest, &
       nzg, &
       nzs, &
       polelat, &
       polelon, &
       runtype, &
       timeunit, &
       timmax, &
       zz,time
  use mem_leaf, only: albedo, &
       drtcon, &
       dthcon, &
       isfcl, &
       nslcon, &
       nvegpat, &
       nvgcon, &
       pctlcon, &
       seatmp, &
       slmstr, &
       slz, &
       stgoff, &
       zrough
  use mem_oda, only: frqoda, &
       if_oda, &
       oda_sfc_tel, &
       oda_sfc_til, &
       oda_sfcprefix, &
       oda_upa_tel, &
       oda_upa_til, &
       oda_upaprefix, &
       roda_hgt, &
       roda_sfc0, &
       roda_sfce, &
       roda_upa0, &
       roda_upae, &
       roda_zfact, &
       tnudoda, &
       todabeg, &
       todaend, &
       wt_oda_grid, &
       wt_oda_pi, &
       wt_oda_rt, &
       wt_oda_th, &
       wt_oda_uv
  use mem_radiate, only: ilwrtyp, &
       iswrtyp, &
       lonrad, &
       radfrq
  use soilMoisture, only: soil_moist, &
       soil_moist_fail, &
       usdata_in, &
       usmodel_in
  use mem_turb, only: akmin, &
       csx, &
       csz, &
       idiffk, &
       if_urban_canopy, &
       ihorgrad, &
       xkhkm, &
       zkhkm
  use mem_varinit, only: cond_hfile, &
       nud_cond, &
       nud_hfile, &
       nud_type, &
       nudlat, &
       t_nudge_rc, &
       tcond_beg, &
       tcond_end, &
       tnudcent, &
       tnudlat, &
       tnudtop, &
       varfpfx, &
       vwait1, &
       vwaittot, &
       wt_nudge_grid, &
       wt_nudge_pi, &
       wt_nudge_rt, &
       wt_nudge_th, &
       wt_nudge_uv, &
       wt_nudgec_grid, &
       znudtop
  use micphys, only: &
       aparm, &
       coltabfn, &
       cparm, &
       gnu, &
       gparm, &
       hparm, &
       iaggr, &
       icloud, &
       igraup, &
       ihail, &
       ipris, &
       irain, &
       isnow, &
       level, &
       mkcoltab, &
       pparm, &
       rparm, &
       sparm
  use node_mod, only: &
       load_bal
  use ref_sounding, only: &
       hs, &
       ipsflg, &
       irtsflg, &
       itsflg, &
       iusflg, &
       ps, &
       rts, &
       ts, &
       us, &
       vs
  use shcu_vars_const, only: &
       nnshcu, &
       shcufrq
  use sib_vars, only: &
       co2_init, &
       n_co2

  use catt_start, only: &
       CATT

  use emission_source_map, only: &
       firemapfn, &
       tracersfn,                           &
       plumerise,                           &
       define_proc

  use plume_utils, only: &
       prfrq

  use mem_scalar, only: &
       recycle_tracers

  use teb_spm_start, only: &
       teb_spm

  use mem_emiss, only : &
       ichemi,          &
       ichemi_in,       &
       chemdata_in,     &
       isource,         &
       weekdayin,       &
       efsat,           &
       efsun,           &
       eindno,          &
       eindno2,         &
       eindpm,          &
       eindco,          &
       eindso2,         &
       eindvoc,         &
       eveino,          &
       eveino2,         &
       eveipm,          &
       eveico,          &
       eveiso2,         &
       eveivoc

  use teb_vars_const, only : &
       rushh1,               &
       rushh2,               &
       daylight,             &
       iteb,                 &
       tminbld,              &
       nteb,                 &
       hc_roof,              &
       tc_roof,              &
       d_roof,               &
       hc_road,              &
       d_road,               &
       tc_road,              &
       d_wall,               &
       tc_wall,              &
       hc_wall,              &
       nurbtype,             &
       ileafcod,             &
       z0_town,              &
       bld,                  &
       bld_height,           &
       bld_hl_ratio,         &
       aroof,                &
       eroof,                &
       aroad,                &
       eroad,                &
       awall,                &
       ewall,                &
       htraf,                &
       hindu,                &
       pletraf,              &
       pleindu

  ! Explicit domain decomposition
  use domain_decomp, only: &
       domain_fname

  implicit none

  include "files.h"
  
!--(DMK-LFR NEC-SX6)----------------------------------------------
!  character(len=f_name_length), intent(in) :: fileName  ! file name with namelists
  character(len=*), intent(in) :: fileName  ! file name with namelists
!--(DMK-LFR NEC-SX6)----------------------------------------------

  integer :: tam    !DSM
  integer :: i                        ! loop count
  integer :: iunit                    ! io unit number
  integer, parameter :: firstUnit=20  ! lowest io unit number available
  integer, parameter :: lastUnit=99   ! highest io unit number available
  logical :: op                       ! io unit number opened or not
  logical :: ex                       ! namelist file exists?
  integer :: err                      ! return code on iostat
  character(len=10) :: c0             ! scratch
  character(len=*), parameter :: h="**(ReadNamelist)**"  ! program unit name

  namelist /MODEL_GRIDS/                                               &
       expnme, runtype, timeunit, timmax, load_bal, imonth1, idate1,   &
       iyear1, itime1, ngrids, nnxp, nnyp, nnzp, nzg, nzs, nxtnest,    &
       domain_fname,                                                   &
       if_adap, ihtran, deltax, deltay, deltaz, dzrat, dzmax, zz,      &
       dtlong, nacoust, ideltat, nstratx, nstraty, nndtrat, nestz1,    &
       nstratz1, nestz2, nstratz2, polelat, polelon, centlat, centlon, &
       ninest, njnest, nknest, nnsttop, nnstbot, gridu, gridv, advmnt

  namelist /CATT_INFO/                                                 &
       catt,                                                           &
       firemapfn, recycle_tracers,                                     &
       plumerise, define_proc, prfrq

  namelist /TEB_SPM_INFO/                                              &
       teb_spm,                                                        &
       fusfiles, ifusflg, ifusfn,                                      &
       ichemi, ichemi_in, chemdata_in, isource, weekdayin, rushh1,     &
       rushh2, daylight, efsat, efsun, eindno, eindno2, eindpm,        &
       eindco, eindso2, eindvoc, eveino, eveino2, eveipm, eveico,      &
       eveiso2, eveivoc, iteb, tminbld, nteb, hc_roof, tc_roof,        &
       d_roof, hc_road, tc_road, d_road, hc_wall, tc_wall, d_wall,     &
       nurbtype, ileafcod, z0_town, bld, bld_height, bld_hl_ratio,     &
       aroof, eroof, aroad, eroad, awall, ewall, htraf, hindu,         &
       pletraf, pleindu

  namelist /MODEL_FILE_INFO/                                           &
       initial, nud_type, varfpfx, vwait1, vwaittot, nud_hfile, nudlat,&
       tnudlat, tnudcent, tnudtop, znudtop, wt_nudge_grid, wt_nudge_uv,&
       wt_nudge_th, wt_nudge_pi, wt_nudge_rt, nud_cond, cond_hfile,    &
       tcond_beg, tcond_end, t_nudge_rc, wt_nudgec_grid, if_oda,       &
       oda_upaprefix,oda_sfcprefix, frqoda, todabeg, todaend, tnudoda, &
       wt_oda_grid, wt_oda_uv, wt_oda_th, wt_oda_pi, wt_oda_rt,        &
       roda_sfce, roda_sfc0, roda_upae,roda_upa0, roda_hgt,            &
       roda_zfact, oda_sfc_til, oda_sfc_tel, oda_upa_til, oda_upa_tel, &
       if_cuinv, cu_prefix, tnudcu, wt_cu_grid, tcu_beg, tcu_end,      &
       cu_tel, cu_til, timstr, hfilin, ipastin, pastfn, ioutput,       &
       hfilout, afilout, iclobber, ihistdel, frqhis, frqanl, frqlite,  &
       xlite, ylite, zlite, nlite_vars, lite_vars, avgtim, frqmean,    &
       frqboth, kwrite, frqprt, initfld, prtcputime, topfiles,         &
       sfcfiles, sstfpfx, ndvifpfx, itoptflg, isstflg, ivegtflg,       &
       isoilflg, ndviflg, nofilflg, iupdndvi, iupdsst, itoptfn, isstfn,&
       ivegtfn, isoilfn, ndvifn, itopsflg, toptenh, toptwvl, iz0flg,   &
       z0max, z0fact, mkcoltab, coltabfn

  namelist /MODEL_OPTIONS/ &
       naddsc, icorflg, ibnd, jbnd, cphas, lsflg, nfpt, distim,          &
       iswrtyp, ilwrtyp,                                                 &
       raddatfn,                                                         &
       radfrq, lonrad, nnqparm, closure_type, nnshcu, confrq,            &
       shcufrq, wcldbs, g3d_spread, g3d_smoothh,             &
       g3d_smoothv, npatch, nvegpat, isfcl, n_co2, co2_init, nvgcon,    &
       pctlcon, nslcon, drtcon, zrough, albedo, seatmp, dthcon,          &
       soil_moist, soil_moist_fail, usdata_in, usmodel_in, slz, slmstr,  &
       stgoff, if_urban_canopy, idiffk, ihorgrad, csx, csz, xkhkm, zkhkm,&
       akmin, level, icloud, irain, ipris, isnow, iaggr, igraup, ihail,  &
       cparm, rparm, pparm, sparm, aparm, gparm, hparm,                  &
       gnu

  namelist /MODEL_SOUND/ &
       ipsflg, itsflg, irtsflg, iusflg, hs, ps, ts, rts, us, vs

  namelist /MODEL_PRINT/ &
       nplt, iplfld, ixsctn, isbval

  namelist /ISAN_CONTROL/ &
       iszstage, ivrstage, isan_inc, guess1st, i1st_flg, iupa_flg,       &
       isfc_flg, iapr, iarawi, iasrfce, varpfx, ioflgisz, ioflgvar

  namelist /ISAN_ISENTROPIC/ &
       nisn, levth, nigrids, topsigz, hybbot, hybtop, sfcinf, sigzwt,    &
       nfeedvar, maxsta, maxsfc, notsta, notid, iobswin, stasep, igridfl,&
       gridwt, gobsep, gobrad, wvlnth, swvlnth, respon


  ! defaults (just for arrays); should be revised to accomodate
  ! precise defaults
  !PPL  nnxp = 0
  !PPL  nnyp = 0
  !PPL  nnzp = 0
  !PPL  centlat = 0.0
  !PPL  centlon = 0.0
  !PPL  nxtnest = 0
  !PPL  zz = 0.0
  !PPL  nstratx = 0
  !PPL  nstraty = 0
  !PPL  nndtrat = 0
  !PPL  nstratz1 = 0
  !PPL  nstratz2 = 0
  !PPL  ninest = 0
  !PPL  njnest = 0
  !PPL  nknest = 0
  !PPL  nnsttop = 0
  !PPL  nnstbot = 0
  !PPL  gridu = 0.0
  !PPL  gridv = 0.0
  !PPL  wt_nudge_grid = 0.0
  !PPL  wt_nudgec_grid = 0.0
  !PPL  wt_oda_grid = 0.0
  !PPL  roda_sfce = 0.0
  !PPL  roda_sfc0 = 0.0
  !PPL  roda_upae = 0.0
  !PPL  roda_upa0 = 0.0
  !PPL  roda_hgt  = 0.0
  !PPL  roda_zfact = 0.0
  !PPL  wt_cu_grid = 0.0
  !PPL  lite_vars = " "
  !PPL  itoptflg = 0
  !PPL  isstflg = 0
  !PPL  ivegtflg = 0
  !PPL  isoilflg = 0
  !PPL  ndviflg = 0
  !PPL  nofilflg = 0
  !PPL  itoptfn = " "
  !PPL  isstfn = " "
  !PPL  ivegtfn = " "
  !PPL  isoilfn = " "
  !PPL  ndvifn = " "
  !PPL  itopsflg = 0
  !PPL  toptenh = 0.0
  !PPL  toptwvl = 0.0
  !PPL  iz0flg = 0
  !PPL  z0max = 0.0
  !PPL  nnqparm = 0
  !PPL  gnu = 0.0
  !PPL  iplfld = " "
  !PPL  ixsctn = 0
  !PPL  isbval = 0
  !PPL  us = 0.0
  !PPL  vs = 0.0
  !PPL  ts = 0.0
  !PPL  ps = 0.0
  !PPL  hs = 0.0
  !PPL  rts = 0.0
  !PPL  nnshcu = 0
  !PPL  co2_init = 0.0
  !PPL  slz = 0.0
  !PPL  slmstr = 0.0
  !PPL  stgoff = 0.0
  !PPL  idiffk = 0
  !PPL  csx = 0.0
  !PPL  csz = 0.0
  !PPL  xkhkm = 0.0
  !PPL  zkhkm = 0.0
  !PPL  akmin = 0.0
  !PPL  levth = 0
  !PPL  notid = " "
  !PPL  gridwt = 0.0
  !PPL  wvlnth = 0.0
  !PPL  swvlnth = 0.0
  !PPL  respon = 0.0
  !PPL  domain_fname = ''


  ! PPL - begin defaults
  ! CATT_INFO
  catt                = 0
  firemapfn           = ''
  recycle_tracers     = 0
  plumerise           = 0
  define_proc         = 0
  prfrq               = 3600. ! Initial Value for PlumeRise Frequency - CATT

  ! ISAN_CONTROL
  iszstage            = 1
  ivrstage            = 1
  isan_inc            = 0600
  guess1st	      = 'PRESS'
  i1st_flg	      = 1
  iupa_flg	      = 3
  isfc_flg	      = 3
  iapr		      = './dprep/dp' ! 2
  iarawi	      = ''
  iasrfce	      = ''
  varpfx	      = './ivar/iv-brams4' ! 2
  ioflgisz 	      = 0
  ioflgvar 	      = 1

  ! ISAN_ISENTROPIC
  nisn                = 43
  levth               = 800
  levth(1:nisn)       = (/280,282,284,286,288,290,292,294,296,298,300,303,306,309,&
       312,315,318,321,324,327,330,335,340,345,350,355,360,380,400,420,440, &
       460,480,500,520,540,570,600,630,670,700,750,800/)
  nigrids             = 1
  topsigz	      = 20000.
  hybbot	      = 4000.
  hybtop	      = 6000.
  sfcinf	      = 1000.
  sigzwt	      = 1.
  nfeedvar            = 1.
  maxsta 	      = 150
  maxsfc 	      = 1000
  notsta 	      = 0
  notid               = ''
  iobswin	      = 1800
  stasep	      = .1
  igridfl	      = 3
  gridwt	      = .01
  gobsep	      = 5.
  gobrad	      = 5.
  wvlnth	      = 1000
  swvlnth	      = 500.
  respon	      = .90

  ! MODEL_FILE_INFO
  initial	      = 2 ! 2
  nud_type	      = 2 ! 2
  varfpfx	      = varpfx ! will be rewrited on the namelist read, and after too
  vwait1	      = 0.
  vwaittot	      = 0.
  nud_hfile	      = '' ! 2
  nudlat	      = 5 ! 2
  tnudlat	      = 1800. ! 2
  tnudcent	      = 0. ! 2
  tnudtop	      = 10800. ! 2
  znudtop	      = 16000. ! 2
  wt_nudge_grid       = 1. ! 2
  wt_nudge_grid(1:4)  = (/1.,0.8,0.7,0.5/) ! 2
  wt_nudge_uv	      = 1.
  wt_nudge_th	      = 1.
  wt_nudge_pi	      = 1.
  wt_nudge_rt	      = 1.
  nud_cond	      = 0
  cond_hfile	      = ''
  tcond_beg	      = 0.
  tcond_end	      = 21600.
  t_nudge_rc	      = 3600.
  wt_nudgec_grid      = 0.5
  wt_nudgec_grid(1:4) = (/1.,0.8,0.7,0.5/)
  if_oda	      = 0
  oda_upaprefix       = ''
  oda_sfcprefix       = ''
  frqoda 	      = 300.
  todabeg	      = 0.
  todaend	      = 9999999.
  tnudoda	      = 900.
  wt_oda_grid         = 1.
  wt_oda_grid(1:4)    = (/1.,0.8,0.7,0.5/)
  wt_oda_uv	      = 1.
  wt_oda_th	      = 1.
  wt_oda_pi	      = 1.
  wt_oda_rt	      = 1.
  roda_sfce	      = 0.
  roda_sfce(1:4)      = (/50000.,100.,100.,100./)
  roda_sfc0           = 0.
  roda_sfc0(1:4)      = (/100000.,100000.,100000.,100000./)
  roda_upae           = 0.
  roda_upae(1:4)      = (/100000.,200.,200.,200./)
  roda_upa0           = 0.
  roda_upa0(1:4)      = (/200000.,2000.,2000.,2000./)
  roda_hgt            = 0.
  roda_hgt(1:4)       = (/3000.,3000.,3000.,3000./)
  roda_zfact          = 0.
  roda_zfact(1:4)     = (/100.,100.,100.,100./)
  oda_sfc_til 	      = 21600.
  oda_sfc_tel 	      = 900.
  oda_upa_til 	      = 43200.
  oda_upa_tel 	      = 21600.
  oda_upa_tel 	      = 21600.
  if_cuinv            = 0
  cu_prefix           = ''
  tnudcu              = 900.
  wt_cu_grid          = 1.
  wt_cu_grid(1:4)     = (/1.,1.,0.5,0.5/)
  tcu_beg	      = 0.
  tcu_end	      = 7200.
  cu_tel 	      = 3600.
  cu_til 	      = 21600.
  timstr 	      = 0.
  hfilin 	      = ''
  ipastin	      = 0 ! 2
  pastfn 	      = '' ! 2
  ioutput	      = 2 ! 2
  hfilout	      = './H/H-brams4' ! 2
  afilout	      = './A/A-brams4' ! 2
  iclobber 	      = 0 ! 2
  ihistdel 	      = 0 ! 2
  frqhis	      = 21600. ! 2
  frqanl	      = 10800. ! 2
  frqlite             = 0.
  xlite 	      = '/0:0/'
  ylite 	      = '/0:0/'
  zlite 	      = '/0:0/'
  nlite_vars	      = 4
  lite_vars 	      = ''
  lite_vars(1)        = 'UP'
  lite_vars(2)        = 'VP'
  lite_vars(3)        = 'WP'
  lite_vars(4)        = 'THETA'
  avgtim 	      = 0.
  frqmean	      = 0.
  frqboth	      = 0.
  kwrite 	      = 0
  frqprt 	      = 21600.
  initfld	      = 1
  topfiles	      = './data/toph-'
  sfcfiles	      = './data/sfc-'
  sstfpfx 	      = './data/sst-'
  ndvifpfx	      = './data/ndvi-'
  itoptflg	      = 2 ! 2
  itoptflg(1:4)       = (/2,2,2,2/) ! 2
  isstflg             = 2 ! 2
  isstflg(1:4)        = (/2,2,2,2/) ! 2
  ivegtflg            = 2 ! 2
  ivegtflg(1:4)       = (/2,2,2,2/) ! 2
  isoilflg            = 2 ! 2
  isoilflg(1:4)       = (/2,2,2,2/) ! 2
  ndviflg             = 2 ! 2
  ndviflg(1:4)        = (/2,2,2,2/) ! 2
  nofilflg            = 2
  nofilflg(1:4)       = (/2,2,2,2/)
  iupdndvi            = 1
  iupdsst	      = 1
  itoptfn	      = ''
  itoptfn(1)          = './topo10km/H'
  itoptfn(2:4)        = (/'./topo/EL','./topo/EL','./topo/EL'/)
  isstfn              = ''
  isstfn(1:4)         = (/'./sst/S','./sst/S','./sst/S','./sst/S'/)
  ivegtfn             = ''
  ivegtfn(1:4)        = (/'./soil-fao/FAO','./soil-fao/FAO','./soil-fao/FAO','./soil-fao/FAO'/)
  isoilfn             = ''
  isoilfn(1:4)        = (/'./oge_brams4/OGE','./oge_brams4/OGE','./oge_brams4/OGE','./oge_brams4/OGE'/)
  ndvifn              = ''
  ndvifn(1:4)         = (/'./ndvi-modis/N','./ndvi-modis/N','./ndvi-modis/N','./ndvi-modis/N'/)
  itopsflg            = 0
  itopsflg(1:4)       = (/0,0,0,0/)
  toptenh             = 1
  toptenh(1:4)        = (/1,1,1,1/)
  toptwvl             = 1
  toptwvl(1:4)        = (/1,1,1,1/)
  iz0flg              = 0
  iz0flg(1:4)         = (/0,0,0,0/)
  z0max               = 5.
  z0max(1:4)          = (/5.,5.,5.,5./)
  z0fact              = 0.005
  mkcoltab            = 0 ! duvida!
  coltabfn            = './micro/ct2.0' ! duvida!
  prtcputime          = 0

  ! MODEL_GRIDS
  expnme              = 'BRAMS 41' ! 2
  timeunit 	      = 'h' ! 2
  load_bal 	      = 0
  ngrids 	      = 1 ! 2
  nzg    	      = 9 ! 2
  nzs    	      = 32
  nxtnest	      = 1 ! 2
  nxtnest(1:4)        = (/0,1,2,3/) ! 2
  if_adap	      = 0 ! 2
  ihtran 	      = 1
  deltaz 	      = 250 ! 2
  dzrat  	      = 1.2 ! 2
  dzmax  	      = 1000.
  zz     	      = 19700.
  zz(1:41)            = (/0.0,20.0,46.0,80.0,120.0,165.0,220.0,290.0,380.0,480.0,590.0, &
       720.0,870.0,1030.0,1200.0,1380.0,1595.0,1850.0,2120.0,2410.0,2715.0,  &
       3030.0,3400.0,3840.0,4380.0,5020.0,5800.0,6730.0,7700.0,8700.0,9700.0,&
       10700.,11700.,12700.,13700.,14700.,15700.,16700.,17700.,18700.,19700./)
  dtlong  	      = 30. ! 2
  nacoust 	      = 3
  ideltat 	      = 1 ! 2
  nstratx 	      = 3
  nstratx(1:4)        = (/1,3,3,3/)
  nstraty             = 3
  nstraty(1:4)        = (/1,3,3,3/)
  nndtrat             = 3
  nndtrat(1:4)        = (/1,3,3,3/)
  nestz1              = 0
  nstratz1            = 1
  nstratz1(1:4)       = (/2,2,2,1/)
  nestz2              = 0
  nstratz2            = 1
  nstratz2(1:4)       = (/3,3,2,1/)
  ninest              = 3 ! 2
  ninest(1:4)         = (/1,1,2,3/) ! 2
  njnest              = 3 ! 2
  njnest(1:4)         = (/1,1,2,3/) ! 2
  nknest              = 1
  nknest(1:4)         = (/1,1,1,1/)
  nnsttop             = 1
  nnsttop(1:4)        = (/1,1,1,1/)
  nnstbot             = 1
  nnstbot(1:4)        = (/1,1,1,1/)
  gridu               = 0.
  gridu(1:4)          = (/0.,0.,0.,0./)
  gridv               = 0.
  gridv(1:4)          = (/0.,0.,0.,0./)
  advmnt              = 0
  domain_fname        = ''

  ! MODEL_OPTIONS
  naddsc  	      = 0
  icorflg 	      = 1
  ibnd    	      = 1
  jbnd    	      = 1
  cphas   	      = 20.
  lsflg   	      = 0
  nfpt    	      = 0
  distim  	      = 400.
  iswrtyp 	      = 1
  ilwrtyp 	      = 1
  raddatfn	      = "./carma/rad_param.data" ! 2
  radfrq  	      = 900. !2
  lonrad  	      = 1
  nnqparm 	      = 2
  closure_type        = "EN" ! 2
  nnshcu      	      = 1 ! 2
  nnshcu(1:4) 	      = (/1,1,1,1/) !2
  confrq      	      = 600. ! 2
  shcufrq     	      = 600. ! 2
  wcldbs      	      = .0005
  g3d_spread          = 0
  g3d_smoothh         = 0
  g3d_smoothv	      = 0
  npatch      	      = 2 ! 2
  nvegpat     	      = 1
  n_co2       	      = 1
  co2_init    	      = 340.
  co2_init(1:9)       = (/360.,360.,360.,355.,355.,355.,350.,350.,340./)
  isfcl               = 1 ! 2
  nvgcon              = 6 ! 2
  pctlcon             = 1.
  nslcon 	      = 6 ! 2
  zrough 	      = .05
  albedo 	      = .2
  seatmp 	      = 298. ! 2
  dthcon 	      = 0.
  drtcon 	      = 0.
  soil_moist          = "i" ! 2
  soil_moist_fail     = "l" ! 2
  usdata_in  	      = "./soil-moisture/GL_SM.GMNR." ! 2
  usmodel_in 	      = "./data/SM-" ! 2
  slz        	      = -0.1 ! 2
  slz(1:9)   	      = (/-2.0,-1.75,-1.50,-1.25,-1.00,-0.75,-0.50,-0.25,-0.1/)
  slmstr     	      = 0.28 ! 2
  slmstr(1:9)	      = (/0.40,0.37,0.35,0.33,0.32,0.31,0.30,0.29,0.28/) ! 2
  stgoff     	      = 0.0
  stgoff(1:9)	      = (/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
  if_urban_canopy     = 0
  idiffk              = 1 ! 2
  idiffk(1:4)         = (/1,1,1,1/) ! 2
  ihorgrad            = 2
  csx                 = .2
  csx(1:4) 	      = (/.2,.2,.2,.2/)
  csz      	      = .2
  csz(1:4) 	      = (/.35,.35,.35,.35/)
  xkhkm    	      = 3.
  xkhkm(1:4) 	      = (/3.,3.,3.,3./)
  zkhkm      	      = 3.
  zkhkm(1:4) 	      = (/3.,3.,3.,3./)
  akmin      	      = 1.
  akmin(1:4) 	      = (/1.,1.,1.,1./)
  level  	      = 3 ! 2
  icloud 	      = 4 ! 2
  irain  	      = 2 ! 2
  ipris  	      = 5
  isnow  	      = 2 ! 2
  iaggr  	      = 2 ! 2
  igraup 	      = 2 ! 2
  ihail  	      = 2 ! 2
  cparm  	      = .1e9
  rparm  	      = 1e-3
  pparm  	      = 0.
  sparm  	      = 1e-3
  aparm  	      = 1e-3
  gparm  	      = 1e-3
  hparm  	      = 3e-3
  gnu    	      = 2.
  gnu(1:4)	      = (/2.,2.,2.,2./)

  ! MODEL_PRINT
  nplt   	      = 0
  iplfld 	      = ""
  ixsctn 	      = 3
  ixsctn(1:4)         = (/3,3,3,3/)
  isbval              = 2
  isbval(1:4)         = (/2,2,2,2/)

  ! MODEL_SOUND
  ipsflg  	      = 1 ! 2
  itsflg  	      = 0 ! 2
  irtsflg 	      = 3 ! 2
  iusflg  	      = 0 ! 2
  hs      	      = 0. ! 2
  ps(1:11)	      = (/1010.,1000.,2000.,3000.,4000.,6000.,8000.,11000.,15000.,20000.,25000./) ! 2
  ts(1:11)	      = (/25.,18.5,12.,4.5,-11.,-24.,-37.,-56.5,-56.5,-56.5,-56.5/) ! 2
  rts(1:11) 	      = (/70.,70.,70.,70.,20.,20.,20.,20.,10.,10.,10./) ! 2
  us(1:11)  	      = (/10.,10.,10.,10.,10.,10.,10.,10.,10.,10.,10./) ! 2
  vs(1:11)  	      = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./) ! 2

  ! TEB
  teb_spm  	      = 0
  fusfiles 	      = ''
  ifusflg  	      = 0
  ifusflg(1:4)        = (/1,1,1,1/)
  ifusfn 	      = ''
  ifusfn(1:4)	      = (/'./fusos/fuso','./fusos/fuso','./fusos/fuso','./fusos/fuso'/) ! 2
  ichemi     	      = 0    !Photochemical module activation - 1=on, 0=off
  ichemi_in  	      = 1    !Use initial values from previous run (1=yes,0=no) ! 2
  chemdata_in	      = ''
  isource    	      = 1    !Emission module activation - 1=on, 0=off ! 2
  weekdayin  	      = 'SUN'  !Initial weeakday of the simulation
  rushh1     	      = 7.81  !Morning Rush Hour (Local Time in Hours)
  rushh2     	      = 16.0  !Afternoon/Evening Rush Hour (Local Time)
  daylight   	      = 0.    !Daylight saving time (horario de verao)
  ! Emission factor (fraction of weekdays) for Saturdays and Sundays
  ! They are used in the emission module and TEB. - EDF
  efsat     	      = 0.8
  efsun     	      = 0.5
  ! Input GMT difference time variable (To define local time)
  !Industrial emissions (kg/s/m2)
  eindno     	      = 2.6636227e-10
  eindno2    	      = 2.9595805e-11
  eindpm     	      = 4.3421278e-10
  eindco     	      = 8.1599860e-10
  eindso2    	      = 3.6149164e-10
  eindvoc    	      = 2.5367833e-10 
  !Veicular emissions (kg/day/m2)
  eveino     	      = 4.3196708e-04
  eveino2    	      = 6.8566209e-05
  eveipm     	      = 6.2648396e-06
  eveico     	      = 7.5433785e-03
  eveiso2    	      = 4.0730592e-05
  eveivoc    	      = 1.1892237e-03
  !----- Urban canopy parameterization using TEB (Masson, 2000)-------------
  iteb      	      = 0     !1=on, 0=off
  tminbld   	      = 12.   !Minimum internal building temperature (degrees Celsius)
  nteb      	      = 3     !Number of roof,road and wall layers used in TEB, Max.3
  ! ROOF layers properties
  ! heat capacity
  hc_roof             = 0.
  hc_roof(1:3)        = (/2110000.,280000.,290000./)
  ! thermal conductivity
  tc_roof             = 0.
  tc_roof(1:3)        = (/0.41,0.05,0.03/)
  ! depth
  d_roof              = 0.
  d_roof(1:3)         = (/0.05,0.4,0.05/)
  ! ROAD layers properties
  ! heat capacity
  hc_road             = 0.
  hc_road(1:3)        = (/1240000.,1280000.,1280000./)
  ! thermal conductivity 1.01
  tc_road             = 1.0103
  ! depth
  d_road              = 0.
  d_road(1:3)         = (/0.05,0.1,1.0/)
  ! WALL layers properties
  ! heat capacity J/m3/K 10e6
  hc_wall             = 1000000.
  ! thermal conductivity 0.81 W/m/K
  tc_wall             = 0.81
  ! depth
  d_wall              = 0.
  d_wall(1:3)         = (/0.02,0.125,0.02/)
  nurbtype            = 2	  !Number of urban types (maximum of 3)

  !Leaf class code to identify each urban type
  ileafcod            = 19
  ileafcod(1:2)       = (21,19)
  !Urban type properties
  !Urban type roughness length 5 e 1
  z0_town             = 0.0
  z0_town(1:2)        = (/3.0,0.5/)
  !Fraction occupied by buildings in the grid cell
  bld      	      = 0.0
  bld(1:2) 	      = (/0.5,0.7/)
  !Building Height
  bld_height          = 0
  bld_height(1:2)     = (/50.,5.0/)
  !Vertical/Horizontal rate 3 e 0.5
  bld_hl_ratio        = 2.4
  bld_hl_ratio(1:2)   = (/4.4,2.4/)
  !Roof albedo
  aroof               = 0.15
  !Roof emissivitiy
  eroof               = 0.9
  !Road albedo
  aroad               = 0.1
  !Road emissivity 90% masson
  eroad               = 0.9
  !Wall albedo
  awall               = 0.25
  !Wall emissivity
  ewall               = 0.85
  !Maximum value of sensible heat
  htraf      	      = 0.
  htraf(1:2) 	      = (/90.0,60.0/)
  !Maximum value of sensible heat
  hindu      	      = 14.
  hindu(1:2) 	      = (/10.0,14.0/)
  !released by Industry (W/m2)
  !Maximum value of latent heat
  pletraf             = 5.
  pletraf(1:2)        = (/10.0,5.0/)
  !released by Traffic (W/m2)
  !Maximum value of latent heat
  pleindu             = 50.
  pleindu(1:2)        = (/30.0,50.0/)

  !PPL defaults for necessary variables
  ! (values from the first-time users)
  ! ATENTION: those variables below should be declared on the simplest RAMSIN
  runtype 	      = "initial"
  timmax  	      = 24
  imonth1 	      = 01
  idate1  	      = 25
  iyear1  	      = 2005
  itime1  	      = 0000
  nnxp    	      = 35
  nnyp    	      = 34
  nnzp    	      = 32
  deltax  	      = 112000.
  deltay  	      = 112000.
  polelat 	      = -23.
  polelon 	      = -52.5
  centlat 	      = -22.
  centlon             = -56.

  ! PPL end defaults

  ! select unused i/o unit

  do iunit = firstUnit, lastUnit
     inquire(iunit,opened=op)
     if (.not. op) exit
  end do

  if (iunit > lastUnit) then
     call fatal_error(h//" all i/o units in use")
  end if

  ! if namelist file exists, open, read each section and close

  inquire(file=fileName(1:len_trim(fileName)), exist=ex)
  if (.not. ex) then
     call fatal_error(h//" namelist file "//trim(fileName)//&
          " does not exist")
  end if

  open(iunit, file=fileName(1:len_trim(fileName)), status="old", action="read",&
       iostat=err)
  if (err /= 0) then
     write(c0,"(i10)") err
     call fatal_error(h//" open namelist file "//trim(fileName)//&
          " returned iostat="//trim(adjustl(c0)))
  end if

  read (iunit, iostat=err, NML=MODEL_GRIDS)
  if (err /= 0) then
     write(*,"(a)") h//"**(ERROR)** reading section MODEL_GRIDS "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") h//" compare values read with file contents:"
     write(*,*) "expnme=",trim(expnme)
     write(*,*) "runtype=",trim(runtype)
     write(*,*) "timeunit=",trim(timeunit)
     write(*,*) "timmax=",timmax
     write(*,*) "load_bal=",load_bal
     write(*,*) "imonth1=",imonth1
     write(*,*) "idate1=",idate1
     write(*,*) "iyear1=",iyear1
     write(*,*) "itime1=",itime1
     write(*,*) "ngrids=",ngrids
     write(*,*) "nnxp=",nnxp
     write(*,*) "nnyp=",nnyp
     write(*,*) "nnzp=",nnzp
     write(*,*) "nzg=",nzg
     write(*,*) "nzs=",nzs
     write(*,*) "nxtnest=",nxtnest
     write(*,*) "DOMAIN_FNAME=", DOMAIN_FNAME
     write(*,*) "if_adap=",if_adap
     write(*,*) "ihtran=",ihtran
     write(*,*) "deltax=",deltax
     write(*,*) "deltay=",deltay
     write(*,*) "deltaz=",deltaz
     write(*,*) "dzrat=",dzrat
     write(*,*) "dzmax=",dzmax
     write(*,*) "zz=",zz
     write(*,*) "dtlong=",dtlong
     write(*,*) "nacoust=",nacoust
     write(*,*) "ideltat=",ideltat
     write(*,*) "nstratx=",nstratx
     write(*,*) "nstraty=",nstraty
     write(*,*) "nndtrat=",nndtrat
     write(*,*) "nestz1=",nestz1
     write(*,*) "nstratz1=",nstratz1
     write(*,*) "nestz2=",nestz2
     write(*,*) "nstratz2=",nstratz2
     write(*,*) "polelat=",polelat
     write(*,*) "polelon=",polelon
     write(*,*) "centlat=",centlat
     write(*,*) "centlon=",centlon
     write(*,*) "ninest=",ninest
     write(*,*) "njnest=",njnest
     write(*,*) "nknest=",nknest
     write(*,*) "nnsttop=",nnsttop
     write(*,*) "nnstbot=",nnstbot
     write(*,*) "gridu=",gridu
     write(*,*) "gridv=",gridv
     write(*,*) "advmnt=",advmnt
     call fatal_error(h//" reading namelist")
  end if

  read (iunit, iostat=err, NML=CATT_INFO)
  if (err /= 0) then
     write(*,"(a)") h//"**(ERROR)** reading section CATT_INFO "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") h//" compare values read with file contents:"
     print *, "CATT=", CATT
     write(*,*) "FIREMAPFN=", FIREMAPFN
     !!write(*,*) "TRACERSFN=", TRACERSFN
     write(*,*) "RECYCLE_TRACERS=", RECYCLE_TRACERS
     write(*,*) "PLUMERISE=", PLUMERISE
     write(*,*) "DEFINE_PROC=", define_proc
     write(*,*) "PRFRQ=", PRFRQ
     call fatal_error(h//" reading namelist")
  end if

  read (iunit, iostat=err, NML=TEB_SPM_INFO)
  if (err /= 0) then
     write(*,"(a)") h//"**(ERROR)** reading section TEB_SPM_INFO "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") h//" compare values read with file contents:"
     print *, "TEB_SPM=", TEB_SPM
     print *, "ifusflg=", ifusflg
     print *, "ifusfn=", ifusfn
     print *, "fusfiles=", trim(fusfiles)
     print *, "ICHEMI=", ICHEMI
     print *, "ICHEMI_IN=", ICHEMI_IN
     print *, "CHEMDATA_IN=", CHEMDATA_IN
     print *, "ISOURCE=", ISOURCE
     print *, "WEEKDAYIN=", trim(WEEKDAYIN)
     print *, "RUSHH1=", RUSHH1
     print *, "RUSHH2=", RUSHH2
     print *, "DAYLIGHT=", DAYLIGHT
     print *, "EFSAT=", EFSAT
     print *, "EFSUN=", EFSUN
     print *, "EINDNO=", EINDNO
     print *, "EINDNO2=", EINDNO2
     print *, "EINDPM=", EINDPM
     print *, "EINDCO=", EINDCO
     print *, "EINDSO2=", EINDSO2
     print *, "EINDVOC=", EINDVOC
     print *, "EVEINO=", EVEINO
     print *, "EVEINO2=", EVEINO2
     print *, "EVEIPM=", EVEIPM
     print *, "EVEICO=", EVEICO
     print *, "EVEISO2=", EVEISO2
     print *, "EVEIVOC=", EVEIVOC
     print *, "ITEB=", ITEB
     print *, "TMINBLD=", TMINBLD
     print *, "NTEB=", NTEB
     print *, "HC_ROOF=", HC_ROOF
     print *, "TC_ROOF=", TC_ROOF
     print *, "D_ROOF=", D_ROOF
     print *, "HC_ROAD=", HC_ROAD
     print *, "TC_ROAD=", TC_ROAD
     print *, "D_ROAD=", D_ROAD
     print *, "HC_WALL=", HC_WALL
     print *, "TC_WALL=", TC_WALL
     print *, "D_WALL=", D_WALL
     print *, "NURBTYPE=", NURBTYPE
     print *, "ILEAFCOD=", ILEAFCOD
     print *, "Z0_TOWN=", Z0_TOWN
     print *, "BLD=", BLD
     print *, "BLD_HEIGHT=", BLD_HEIGHT
     print *, "BLD_HL_RATIO=", BLD_HL_RATIO
     print *, "AROOF=", AROOF
     print *, "EROOF=", EROOF
     print *, "AROAD=", AROAD
     print *, "EROAD=", EROAD
     print *, "AWALL=", AWALL
     print *, "EWALL=", EWALL
     print *, "HTRAF=", HTRAF
     print *, "HINDU=", HINDU
     print *, "PLETRAF=", PLETRAF
     print *, "PLEINDU=", PLEINDU
     call fatal_error(h//" reading namelist")
  end if

  read (iunit, iostat=err, NML=MODEL_FILE_INFO)
  if (err /= 0) then
     write(*,"(a)") h//"**(ERROR)** reading section MODEL_FILE_INFO "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") h//" compare values read with file contents:"
     write (*, "(a)") " namelist MODEL_FILE_INFO: "
     write (*,*) "initial=", initial
     write (*,*) "nud_type=", nud_type
     write (*,*) "varfpfx=", trim(varfpfx)
     write (*,*) "vwait1=", vwait1
     write (*,*) "vwaittot=", vwaittot
     write (*,*) "nud_hfile=", trim(nud_hfile)
     write (*,*) "nudlat=", nudlat
     write (*,*) "tnudlat=", tnudlat
     write (*,*) "tnudcent=", tnudcent
     write (*,*) "tnudtop=", tnudtop
     write (*,*) "znudtop=", znudtop
     write (*,*) "wt_nudge_grid=", wt_nudge_grid
     write (*,*) "wt_nudge_uv=", wt_nudge_uv
     write (*,*) "wt_nudge_th=", wt_nudge_th
     write (*,*) "wt_nudge_pi=", wt_nudge_pi
     write (*,*) "wt_nudge_rt=", wt_nudge_rt
     write (*,*) "nud_cond=", nud_cond
     write (*,*) "cond_hfile=", trim(cond_hfile)
     write (*,*) "tcond_beg=", tcond_beg
     write (*,*) "tcond_end=", tcond_end
     write (*,*) "t_nudge_rc=", t_nudge_rc
     write (*,*) "wt_nudgec_grid=", wt_nudgec_grid
     write (*,*) "if_oda=", if_oda
     write (*,*) "oda_upaprefix=", trim(oda_upaprefix)
     write (*,*) "oda_sfcprefix=", trim(oda_sfcprefix)
     write (*,*) "frqoda=", frqoda
     write (*,*) "todabeg=", todabeg
     write (*,*) "todaend=", todaend
     write (*,*) "tnudoda=", tnudoda
     write (*,*) "wt_oda_grid=", wt_oda_grid
     write (*,*) "wt_oda_uv=", wt_oda_uv
     write (*,*) "wt_oda_th=", wt_oda_th
     write (*,*) "wt_oda_pi=", wt_oda_pi
     write (*,*) "wt_oda_rt=", wt_oda_rt
     write (*,*) "roda_sfce=", roda_sfce
     write (*,*) "roda_sfc0=", roda_sfc0
     write (*,*) "roda_upae=", roda_upae
     write (*,*) "roda_upa0=", roda_upa0
     write (*,*) "roda_hgt=", roda_hgt
     write (*,*) "roda_zfact=", roda_zfact
     write (*,*) "oda_sfc_til=", oda_sfc_til
     write (*,*) "oda_sfc_tel=", oda_sfc_tel
     write (*,*) "oda_upa_til=", oda_upa_til
     write (*,*) "oda_upa_tel=", oda_upa_tel
     write (*,*) "if_cuinv=", if_cuinv
     write (*,*) "cu_prefix=", trim(cu_prefix)
     write (*,*) "tnudcu=", tnudcu
     write (*,*) "wt_cu_grid=", wt_cu_grid
     write (*,*) "tcu_beg=", tcu_beg
     write (*,*) "tcu_end=", tcu_end
     write (*,*) "cu_tel=", cu_tel
     write (*,*) "cu_til=", cu_til
     write (*,*) "timstr=", timstr
     write (*,*) "hfilin=", trim(hfilin)
     write (*,*) "ipastin=", ipastin
     write (*,*) "pastfn=", trim(pastfn)
     write (*,*) "ioutput=", ioutput
     write (*,*) "hfilout=", trim(hfilout)
     write (*,*) "afilout=", trim(afilout)
     write (*,*) "iclobber=", iclobber
     write (*,*) "ihistdel=", ihistdel
     write (*,*) "frqhis=", frqhis
     write (*,*) "frqanl=", frqanl
     write (*,*) "frqlite=", frqlite
     write (*,*) "xlite=", xlite
     write (*,*) "ylite=", ylite
     write (*,*) "zlite=", zlite
     write (*,*) "nlite_vars=", nlite_vars
     write (*,*) "lite_vars=", (trim(lite_vars(i))//";", i=1,size(lite_vars))
     write (*,*) "avgtim=", avgtim
     write (*,*) "frqmean=", frqmean
     write (*,*) "frqboth=", frqboth
     write (*,*) "kwrite=", kwrite
     write (*,*) "frqprt=", frqprt
     write (*,*) "initfld=", initfld
     write (*,*) "prtcputime", prtcputime
     write (*,*) "topfiles=", trim(topfiles)
     write (*,*) "sfcfiles=", trim(sfcfiles)
     write (*,*) "sstfpfx=", trim(sstfpfx)
     write (*,*) "ndvifpfx=", trim(ndvifpfx)
     write (*,*) "itoptflg=", itoptflg
     write (*,*) "isstflg=", isstflg
     write (*,*) "ivegtflg=", ivegtflg
     write (*,*) "isoilflg=", isoilflg
     write (*,*) "ndviflg=", ndviflg
     write (*,*) "nofilflg=", nofilflg
     write (*,*) "iupdndvi=", iupdndvi
     write (*,*) "iupdsst=", iupdsst
     write (*,*) "itoptfn=", (trim(itoptfn(i))//";", i =1,size(itoptfn))
     write (*,*) "isstfn=", (trim(isstfn(i))//";", i=1,size(isstfn))
     write (*,*) "ivegtfn=", (trim(ivegtfn(i))//";", i = 1, size(ivegtfn))
     write (*,*) "isoilfn=", (trim(isoilfn(i))//";", i = 1, size(isoilfn))
     write (*,*) "ndvifn=", (trim(ndvifn(i))//";", i=1,size(ndvifn))
     write (*,*) "itopsflg=", itopsflg
     write (*,*) "toptenh=", toptenh
     write (*,*) "toptwvl=", toptwvl
     write (*,*) "iz0flg=", iz0flg
     write (*,*) "z0max=", z0max
     write (*,*) "z0fact=", z0fact
     write (*,*) "mkcoltab=", mkcoltab
     write (*,*) "coltabfn=", trim(coltabfn)
     call fatal_error(h//" reading namelist")
  end if

  read (iunit, iostat=err, NML=MODEL_OPTIONS)
  if (err /= 0) then
     write(*,"(a)") h//"**(ERROR)** reading section MODEL_OPTIONS "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") h//" compare values read with file contents:"
     write (*,*) "naddsc=",naddsc
     write (*, *) "icorflg=",icorflg
     write (*, *) "ibnd=",ibnd
     write (*, *) "jbnd=",jbnd
     write (*, *) "cphas=",cphas
     write (*, *) "lsflg=",lsflg
     write (*, *) "nfpt=",nfpt
     write (*, *) "distim=",distim
     write (*, *) "iswrtyp=",iswrtyp
     write (*, *) "ilwrtyp=",ilwrtyp
     write (*, *) "raddatfn=", RADDATFN
     write (*, *) "radfrq=",radfrq
     write (*, *) "lonrad=",lonrad
     write (*, *) "nnqparm=",nnqparm
     write (*, *) "closure_type=",closure_type
     write (*, *) "nnshcu=",nnshcu
     write (*, *) "confrq=",confrq
     write (*, *) "shcufrq=",shcufrq
     write (*, *) "wcldbs=",wcldbs
     
     write (*, *) "g3d_spread=", g3d_spread
     write (*, *) "g3d_smoothh=", g3d_smoothh
     write (*, *) "g3d_smoothv=", g3d_smoothv
     
     write (*, *) "npatch=",npatch
     write (*, *) "nvegpat=",nvegpat
     write (*, *) "isfcl=",isfcl
     write (*, *) "n_co2=",n_co2
     write (*, *) "co2_init=",co2_init
     write (*, *) "nvgcon=",nvgcon
     write (*, *) "pctlcon=",pctlcon
     write (*, *) "nslcon=",nslcon
     write (*, *) "drtcon=",drtcon
     write (*, *) "zrough=",zrough
     write (*, *) "albedo=",albedo
     write (*, *) "seatmp=",seatmp
     write (*, *) "dthcon=",dthcon
     write (*, *) "soil_moist=",soil_moist
     write (*, *) "soil_moist_fail=",soil_moist_fail
     write (*, *) "usdata_in=",trim(usdata_in)
     write (*, *) "usmodel_in=",trim(usmodel_in)
     write (*, *) "slz=",slz
     write (*, *) "slmstr=",slmstr
     write (*, *) "stgoff=",stgoff
     write (*, *) "if_urban_canopy=",if_urban_canopy
     write (*, *) "idiffk=",idiffk
     write (*, *) "ihorgrad=",ihorgrad
     write (*, *) "csx=",csx
     write (*, *) "csz=",csz
     write (*, *) "xkhkm=",xkhkm
     write (*, *) "zkhkm=",zkhkm
     write (*, *) "akmin=",akmin
     write (*, *) "level=",level
     write (*, *) "icloud=",icloud
     write (*, *) "irain=",irain
     write (*, *) "ipris=",ipris
     write (*, *) "isnow=",isnow
     write (*, *) "iaggr=",iaggr
     write (*, *) "igraup=",igraup
     write (*, *) "ihail=",ihail
     write (*, *) "cparm=",cparm
     write (*, *) "rparm=",rparm
     write (*, *) "pparm=",pparm
     write (*, *) "sparm=",sparm
     write (*, *) "aparm=",aparm
     write (*, *) "gparm=",gparm
     write (*, *) "hparm=",hparm
     write (*, *) "gnu=",gnu
     call fatal_error(h//" reading namelist")
  end if

  read (iunit, iostat=err, NML=MODEL_SOUND)
  if (err /= 0) then
     write(*,"(a)") h//"**(ERROR)** reading section MODEL_SOUND "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") h//" compare values read with file contents:"
     write (*, *) "ipsflg=",ipsflg
     write (*, *) "itsflg=",itsflg
     write (*, *) "irtsflg=",irtsflg
     write (*, *) "iusflg=",iusflg
     write (*, *) "hs=",hs
     write (*, *) "ps=",ps
     write (*, *) "ts=",ts
     write (*, *) "rts=",rts
     write (*, *) "us=",us
     write (*, *) "vs=",vs
     call fatal_error(h//" reading namelist")
  end if

  read (iunit, iostat=err, NML=MODEL_PRINT)
  if (err /= 0) then
     write(*,"(a)") h//"**(ERROR)** reading section MODEL_PRINT "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") h//" compare values read with file contents:"
     write (*, *) "nplt=",nplt
     write (*, *) "iplfld=",(trim(iplfld(i))//";", i=1,size(iplfld))
     write (*, *) "ixsctn=",ixsctn
     write (*, *) "isbval=",isbval
     call fatal_error(h//" reading namelist")
  end if


  read (iunit, iostat=err, NML=ISAN_CONTROL)
  if (err /= 0) then
     write(*,"(a)") h//"**(ERROR)** reading section ISAN_CONTROL "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") h//" compare values read with file contents:"
     write (*, *) "iszstage=",iszstage
     write (*, *) "ivrstage=",ivrstage
     write (*, *) "isan_inc=",isan_inc
     write (*, *) "guess1st=",guess1st
     write (*, *) "i1st_flg=",i1st_flg
     write (*, *) "iupa_flg=",iupa_flg
     write (*, *) "isfc_flg=",isfc_flg
     write (*, *) "iapr=",trim(iapr)
     write (*, *) "iarawi=",trim(iarawi)
     write (*, *) "iasrfce=",trim(iasrfce)
     write (*, *) "varpfx=",trim(varpfx)
     write (*, *) "ioflgisz=",ioflgisz
     write (*, *) "ioflgvar=",ioflgvar
     call fatal_error(h//" reading namelist")
  end if

  read (iunit, iostat=err, NML=ISAN_ISENTROPIC)
  if (err /= 0) then
     write(*,"(a)") h//"**(ERROR)** reading section ISAN_ISENTROPIC "//&
          &"of namelist file "//trim(fileName)
     write(*,"(a)") h//" compare values read with file contents:"
     write (*, *) "nisn=",nisn
     write (*, *) "levth=",levth
     write (*, *) "nigrids=",nigrids
     write (*, *) "topsigz=",topsigz
     write (*, *) "hybbot=",hybbot
     write (*, *) "hybtop=",hybtop
     write (*, *) "sfcinf=",sfcinf
     write (*, *) "sigzwt=",sigzwt
     write (*, *) "nfeedvar=",nfeedvar
     write (*, *) "maxsta=",maxsta
     write (*, *) "maxsfc=",maxsfc
     write (*, *) "notsta=",notsta
     write (*, *) "notid=",(trim(notid(i))//";", i=1,size(notid))
     write (*, *) "iobswin=",iobswin
     write (*, *) "stasep=",stasep
     write (*, *) "igridfl=",igridfl
     write (*, *) "gridwt=",gridwt
     write (*, *) "gobsep=",gobsep
     write (*, *) "gobrad=",gobrad
     write (*, *) "wvlnth=",wvlnth
     write (*, *) "swvlnth=",swvlnth
     write (*, *) "respon=",respon
     call fatal_error(h//" reading namelist")
  end if

  close(iunit, iostat=err)
  if (err /= 0) then
     write(c0,"(i10)") err
     call fatal_error(h//" closing file "//&
          trim(fileName)//" returned iostat="//&
          trim(adjustl(c0)))
  end if

  ! PPL - VARFPFX is the same variable of VARPFX
  VARFPFX = VARPFX

  !DSM{
  if (trim(RUNTYPE)=='HISTORY') then
      tam=len_trim(hfilin)
      read(hfilin(tam-20:tam-19),'(i2.2)') IMONTH1
      read(hfilin(tam-17:tam-16),'(i2.2)') IDATE1
      read(hfilin(tam-25:tam-22),'(i4.4)') IYEAR1
      read(hfilin(tam-14:tam-11),'(i4.4)') ITIME1
   endif
  !DSM}

end subroutine ReadNamelist


!-------------------------------------------------------------

subroutine get_akmin2d(ngr, n2, n3, akmin2d)

  use mem_grid, only: &
       platn, plonn, xmn, ymn ! INTENT(IN)

  use node_mod, only: &
       mynum,         & ! INTENT(IN)
       nodei0,        & ! INTENT(IN)
       nodej0           ! INTENT(IN)

!!$  !tmp-srf para diferente difusao numerica >
!!$  use extras, only: extra2d, NA_EXTRA2D

  implicit none
  ! Arguments:
  integer, intent(IN) :: n2, n3, ngr
  real, intent(OUT)   :: akmin2d(n2,n3)
  ! Local Variables:
  real    :: rlat(n2,n3), rlon(n2,n3)
  integer :: i, j

!!$  print *, "DEBUG-ALF:get_akmin2d"

  !srf-define diferentes AKMINs para melhorar estabilidade 
  !srf-sobre os Andes
!!$  !testa numero de extras:
!!$  if (NA_EXTRA2D<5) call fatal_error('NA_EXTRA2d deve ser no minimo 5')
  !default 
!!$  extra2d(5,ngr)%d2(:,:) = 1.
  akmin2d = 1.
  !----

  !print*,n2,n3,platn(ngr),plonn(ngr),xmn(1,ngr),ymn(1,ngr)
  !-calculate lat, lon of each grid box T-points
  do j=1,n3
     do i=1,n2
        call xy_ll(rlat(i,j), rlon(i,j), platn(ngr), plonn(ngr), &
             xmn(i+nodei0(mynum,ngr),ngr), ymn(j+nodej0(mynum,ngr),ngr))
     enddo
  enddo

  do j=1,n3
     do i=1,n2

        if (rlat(i,j)<-15.) then
           if (rlon(i,j)<-60.) then 
              !       TOPTWVL_2d(i,j)=15.
!!$              extra2d(5,ngr)%d2(i,j) = 3.
              akmin2d(i,j) = 3.
           endif
        endif

        if (rlat(i,j)>=-15. .and. rlat(i,j)<-9.) then
           if (rlon(i,j)<-62.) then
              !       TOPTWVL_2d(i,j)=15.
!!$              extra2d(5,ngr)%d2(i,j) = 3.
              akmin2d(i,j) = 3.
           endif
        endif

        if (rlat(i,j)>=-9. .and. rlat(i,j)<-1.) then
           if (rlon(i,j)<-70.) then 
              !      TOPTWVL_2d(i,j)=15.
!!$              extra2d(5,ngr)%d2(i,j) = 3.
              akmin2d(i,j) = 3.
           endif
        endif

        if (rlat(i,j)>=-1.) then
           if (rlon(i,j)<=-57 .and. rlon(i,j)>-67.) then
              !      TOPTWVL_2d(i,j)=11.
!!$              extra2d(5,ngr)%d2(i,j) = 1.5
              akmin2d(i,j) = 1.5              
           endif
        endif

        if (rlat(i,j)>=-1.) then
           if (rlon(i,j)<=-67.) then
              !       TOPTWVL_2d(i,j)=15.
!!$              extra2d(5,ngr)%d2(i,j) = 3.
              akmin2d(i,j) = 3.
           endif
        endif

        if (rlat(i,j)<=-10.) then
           if (rlon(i,j)>=40.) then 
              !      TOPTWVL_2d(i,j)=11.
!!$              extra2d(5,ngr)%d2(i,j) = 1.5
              akmin2d(i,j) = 1.5
           endif
        endif
     enddo
  enddo

!!$  print *, "DEBUG-ALF:get_akmin2d:sum(akmin2d),media=", &
!!$       sum(akmin2d), sum(akmin2d)/(n2*n3)
!!$  call flush(6)

end subroutine get_akmin2d
!**********************************************************************

!-------------------------------------------------------------

subroutine get_traning_grell(ngr, n2, n3, train)

  use mem_grid, only: &
       platn, plonn, xmn, ymn ! INTENT(IN)

  use node_mod, only: &
       mynum,         & ! INTENT(IN)
       nodei0,        & ! INTENT(IN)
       nodej0           ! INTENT(IN)

!!$  !tmp-srf para diferente difusao numerica >
!!$  use extras, only: extra2d, NA_EXTRA2D

  implicit none
  ! Arguments:
  integer, intent(IN) :: n2, n3, ngr
  real, intent(OUT)   :: train(n2,n3)
  ! Local Variables:
  real    :: rlat(n2,n3), rlon(n2,n3), rm
  integer :: i, j, np,ii,jj

!!$  print *, "DEBUG-ALF:get_train"

  !srf-define diferentes AKMINs para melhorar estabilidade 
  !srf-sobre os Andes
!!$  !testa numero de extras:
!!$  if (NA_EXTRA2D<5) call fatal_error('NA_EXTRA2d deve ser no minimo 5')
  !default 
!!$  extra2d(5,ngr)%d2(:,:) = 1.
  !train = 1.17 !valor do La Plata.
  !----

  !print*,n2,n3,platn(ngr),plonn(ngr),xmn(1,ngr),ymn(1,ngr)
  !-calculate lat, lon of each grid box T-points
  do j=1,n3
     do i=1,n2
        call xy_ll(rlat(i,j), rlon(i,j), platn(ngr), plonn(ngr), &
             xmn(i+nodei0(mynum,ngr),ngr), ymn(j+nodej0(mynum,ngr),ngr))
     enddo
  enddo

  do j=1,n3
     do i=1,n2

!  train(i,j) = 1.17 !valor do La Plata.
  train(i,j) = 1.3 !valor do La Plata.! 30102011



!regiao norte
        if (rlat(i,j) > -14.) then
           if (rlon(i,j)<-45.5) then 
 !             train(i,j) = 0.43! 30102011
              train(i,j) = 0.3
           endif
        endif

!regiao NORDESTE
        if (rlat(i,j) > -18.5) then
           if (rlon(i,j) .GE. -45.5) then 
!              train(i,j) = 0.25! 30102011
              train(i,j) = 0.3
           endif
        endif
!regiao CENTRO OESTE
        if (rlat(i,j) > -25. .AND. RLAT(I,J) < -12.) then
           if (rlon(i,j) < -45.5) then 
              train(i,j) = 0.75
           endif
        endif

!regiao SULDESTE
        if (rlat(i,j) > -26.5 .AND. RLAT(I,J) < -14.) then
           if (rlon(i,j) .GE. 45.5) then 
              train(i,j) = 0.32
           endif
        endif

!REGIAO SUL
        if (rlat(i,j) > -34. .AND. RLAT(I,J) .LE. -26.5) then
           if (rlon(i,j) .GE. 45.5) then 
!              train(i,j) = 0.93 ! 30102011
              train(i,j) = 1.3  
           endif
        endif

     enddo
  enddo

!media movel ( 2 viz)
  do j=3,n3-2
     do i=3,n2-2
      rm=0.
      np=0
       do jj=j-1,j+1
         do ii=i-1,i+1
	  np=np+1
	  rm = rm + train(ii,jj)
         enddo
        enddo
        train(I,J)=rm/np
     enddo
  enddo

!!$  print *, "DEBUG-ALF:get_train:sum(train),media=", &
!!$       sum(train), sum(train)/(n2*n3)
!!$  call flush(6)

end subroutine get_traning_grell
!**********************************************************************

