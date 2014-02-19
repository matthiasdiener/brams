!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


subroutine rams_mem_alloc(proc_type)

  use grid_dims, only: &
       maxmach, &
       maxsclr, &
       maxgrds

  use mem_tend, only: &
       nullify_tend, &
       alloc_tend, &
       filltab_tend

  use mem_scratch1, only: &
       alloc_scratch1, &
       nullify_scratch1

  use mem_nestb, only: &
       alloc_nestb

  use mem_scratch, only: &
       alloc_scratch,    &
       nullify_scratch,  &
       createvctr

  use io_params, only: &
       avgtim

  use var_tables, only: &
       num_scalar, &
       maxvars, &
       num_var, &
       nvgrids, &
       vtab_r,  &
       scalar_tab

  use mem_grid, only: &
       nxtnest, &
       oneGlobalGridData, &
       grid_g, &
       gridm_g, &
       naddsc, &
       maxnxp, &
       maxnyp, &
       maxnzp, &
       if_adap, &
       npatch, &
       nnxp, &
       nnyp, &
       nnzp, &
       nzg,  &
       nzs,  &
       ngrids, &
       nullify_grid, &
       alloc_grid, &
       nullify_GlobalGridData, &
       alloc_GlobalGridData, &

!--(DMK-CCATT-INI)-----------------------------------------------------
       dtlong, &
       nndtrat, &
!--(DMK-CCATT-FIM)-----------------------------------------------------
      
       filltab_grid

  use mem_oda, only: &
       oda_g, &
       odam_g, &
       nullify_oda, &
       alloc_oda, &
       filltab_oda

  use mem_cuparm, only: &
       nnqparm, &
       cuparm_g_sh, &
       cuparmm_g_sh, &
       cuparm_g, &
       cuparmm_g, &
       nullify_cuparm, &
       alloc_cuparm, &
       alloc_cuparm_sh, &
       filltab_cuparm, &
       filltab_cuparm_sh

  use mem_scalar, only: &
       scalar_g, &
       scalarm_g, &
       alloc_scalar, &
       nullify_scalar, &
       filltab_scalar

  use mem_radiate, only: &
       ilwrtyp, &
       iswrtyp, &
       radiate_g, &
       radiatem_g, &
       nullify_radiate, &
       alloc_radiate, &
       filltab_radiate

  use mem_micro, only: &
       micro_g, &
       microm_g, &
       nullify_micro, &
       alloc_micro, &
       filltab_micro

  use mem_leaf, only: &
       isfcl,  &
       leaf_g, &
       leafm_g, &
       nullify_leaf, &
       alloc_leaf, &
       filltab_leaf

  use mem_jules, only: &
       jules_g, &
       julesm_g, &
       nullify_jules, &
       alloc_jules, &
       filltab_jules

  use mem_varinit, only: &
       varinit_g, &
       varinitm_g, &
       nullify_varinit, &
       alloc_varinit, &
       filltab_varinit

  use mem_turb, only: &
       ihorgrad, &
       turb_g, &
       turbm_g, &
       nullify_turb, &
       alloc_turb, &
       filltab_turb

  use mem_basic, only: &
       basic_g,        &
       basicm_g,       &
       nullify_basic,  &
       alloc_basic,    &
       filltab_basic

  use node_mod, only: &
       alloc_paths,   & !Subroutine
       nmachs,        &
       ipaths,        &
       iget_paths,    &
       nodemxp,       &
       nodemyp,       &
       nodebounds,    &
       mchnum,         &
       master_num,     &
       mynum

  use mem_shcu, only: &
       shcu_g,        &
       shcum_g,       &
       nullify_shcu,  &
       alloc_shcu,    &
       filltab_shcu

  use shcu_vars_const, only : &
       nnshcu                           ! INTENT(IN)

  use mem_grell_param, only: &
       ngrids_cp,            &
       flag_grell,           &
       closure_type,         &
       icoic,                &
       icoic_sh,             &
       define_memory

  use mem_grell, only: &
       grell_g,        &
       grellm_g,       &
       grell_g_sh,     &
       grellm_g_sh,    &
       alloc_grell,    &
       nullify_grell,  &
       filltab_grell,  &
       filltab_grell_sh, &
       cuforc_g, &
       cuforc_sh_g, &
       cuforcm_g, &
       cuforcm_sh_g, &
       nullify_cuforc, &
       alloc_cu_forcings, &
       filltab_cuforc_sh, &
       filltab_cuforc, &
       alloc_grell_sh

  use mem_scratch1_grell, only: &
       alloc_scratch1_grell

  use mem_scratch2_grell, only: &
       alloc_scratch2_grell

  use mem_scratch3_grell, only: &
       alloc_scratch3_grell

  use sib_vars, only : &
       N_CO2 ! INTENT(IN)

  use mem_sib_co2, only: &
       sib_g,            &
       sibm_g,           &
       alloc_sib_co2,    &
       zero_sib_co2,     &
       filltab_sib_co2,  &
       nullify_sib_co2

  use mem_sib, only :   &
       sib_brams_g,     &
       sib_bramsm_g,    &
       alloc_sib_brams, &
       nullify_sib_brams, &
       dealloc_sib_brams, &
       filltab_sib_brams, &
       zero_sib_brams

  use mem_opt, only: &
       nullify_opt_scratch, &
       alloc_opt_scratch

!--(DMK-CCATT-INI)-----------------------------------------------------
!  use catt_start, only: &
!       CATT           ! intent(in)
!--(DMK-CCATT-FIM)-----------------------------------------------------

  use mem_carma, only: &
       carma,          &
       carma_m,        &
       nullify_carma,  &
       filltab_carma,  &
       alloc_carma,    &
       zero_carma

  use mem_aerad, only: &
       nwave,          &           !INTENT(IN)
       initial_definitions_aerad !Subroutine

  use mem_globaer, only: &
       initial_definitions_globaer !Subroutine

  use mem_globrad, only: &
       initial_definitions_globrad !Subroutine

!--(DMK-CCATT-INI)-----------------------------------------------------
  use ccattExtras, only: &
!--(DMK-CCATT-OLD)-----------------------------------------------------
!  use extras, only: &
!--(DMK-CCATT-FIM)-----------------------------------------------------
       extra2d,     &
       extra3d,     &
       extra2dm,    &
       extra3dm,    &
       na_extra2d,  &
       na_extra3d,  &
       nullify_extra2d, &
       alloc_extra2d,   &
       zero_extra2d,    &
       filltab_extra2d, &
       nullify_extra3d, &
       alloc_extra3d,   &
       filltab_extra3d, &
       zero_extra3d

  use mem_turb_scalar, only: &
       turb_s,         &
       nullify_turb_s, &
       alloc_turb_s,   &
       filltab_turb_s

  use mem_scratch2_grell_sh, only: &
       alloc_scratch2_grell_sh

  use mem_scratch3_grell_sh, only: &
       alloc_scratch3_grell_sh


  ! Data for Optimization for vector machines
  use mem_micro_opt, only: &
       alloc_micro_opt

  ! For specific optimization depending the type of machine
  use machine_arq, only: machine ! INTENT(IN)

  ! Global grid dimension definitions
  use mem_grid_dim_defs, only: define_grid_dim_pointer ! subroutine

  ! TEB_SPM
  use teb_spm_start, only: TEB_SPM ! INTENT(IN)

  use mem_emiss, only: isource ! INTENT(IN)

  use mem_gaspart, only: &
       gaspart_g,        & ! INTENT(IN)
       gaspartm_g,       & ! INTENT(IN)
       gaspart_vars,     & ! Type
       nullify_gaspart,  & ! Subroutine
       alloc_gaspart,    & ! Subroutine
       zero_gaspart,     & ! Subroutine
       filltab_gaspart     ! Subroutine

  use mem_teb_common, only: &
       tebc_g,              & ! INTENT(IN)
       tebcm_g,             & ! INTENT(IN)
       nullify_tebc,        & ! Subroutine
       alloc_tebc,          & ! Subroutine
       filltab_tebc           ! Subroutine

  use teb_vars_const, only: iteb ! INTENT(IN)

  use mem_teb, only: &
       teb_g,        & ! INTENT(IN)
       tebm_g,       & ! INTENT(IN)
       nullify_teb,  & ! Subroutine
       alloc_teb,    & ! Subroutine
       filltab_teb     ! Subroutine
       
  use cuparm_grell3

!--(DMK-CCATT-INI)-----------------------------------------------------
  use chem1_list, only:  &
       PhotojMethod,    &       ! INTENT(IN)
       nspecies_chem=> nspecies ! INTENT(IN)

  use aer1_list, only: &        ! INTENT(IN)
       on,             &        ! INTENT(IN)
       mode_alloc,     &        ! INTENT(IN)
       nmodes,         &
       nspecies_aer=> nspecies  ! INTENT(IN)

  use chem1aq_list, only: &
       nspeciesaq_chem=> nspeciesaq ! INTENT(IN)

  use ccatt_start, only: &
       ccatt                ! CHEM_RAMSIN

  use mem_chem1, only:     &
       nullify_chem1,      & ! Subroutine
       alloc_chem1,        & ! Subroutine
       filltab_chem1,      & ! Subroutine
       define_n_dyn_chem,  & ! Subroutine
       nullify_tend_chem1, & ! Subroutine
       alloc_tend_chem1,   & ! Subroutine
       filltab_tend_chem1, & ! Subroutine
       nsrc,               & ! INTENT(IN)
       max_ntimes_src,     & ! INTENT(IN)
       chemistry,          & ! CHEM_RAMSIN
       chem1_g,            &
       chem1m_g,           &
       chem1_src_g,        &
       chem1m_src_g,       &
       chem1_src_z_dim_g

  use mem_chemic, only: &
       nullify_chemic,  & ! Subroutine
       alloc_chemic,    & ! Subroutine
       chemic_g

  use mem_chem1aq, only:     &
       nullify_chem1aq,      & ! Subroutine
       alloc_chem1aq,        & ! Subroutine
       filltab_chem1aq,      & ! Subroutine
       nullify_tend_chem1aq, & ! Subroutine
       alloc_tend_chem1aq,   & ! Subroutine
       filltab_tend_chem1aq, & ! Subroutine
       chemistry_aq,         & ! CHEM_RAMSIN 
       chem1aq_g,            &
       chem1maq_g

  use mem_aer1, only:        &
       nullify_aer1,         & ! Subroutine
       alloc_aer1,           & ! Subroutine
       filltab_aer1,         & ! Subroutine     
       nullify_tend_aer1,    & ! Subroutine
       alloc_tend_aer1,      & ! Subroutine
       filltab_tend_aer1,    & ! Subroutine
       aer1_g,               &
       aer1m_g,              &
       aer1_src_z_dim_g

  use mem_plume_chem1, only: &
       nullify_plume_chem1,  & ! Subroutine
       alloc_plume_chem1,    & ! Subroutine
       filltab_plume_chem1,  & ! Subroutine
       nveg_agreg,           & ! INTENT(IN)
       plumerise,            & ! CHEM_RAMSIN
       plume_g,              &
       plumem_g,             &
       plume_mean_g,         &
       plume_meanm_g

  use mem_volc_chem1, only: &
       nullify_volc_chem1,  & ! Subroutine
       alloc_volc_chem1,    & ! Subroutine
       filltab_volc_chem1,  & ! Subroutine
       volcanoes,           & ! CHEM_RAMSIN
       volc_mean_g,         &
       volc_meanm_g

  use chem_sources, only:       &
       oneGlobalEmissData,      &
       alloc_GlobalEmissData,   &
       nullify_GlobalEmissData
      
  use module_dry_dep, only: &
       alloc_aer_sedim        ! Subroutine

  use mem_stilt, only: &
       nullify_stilt,  & ! Subroutine
       alloc_stilt,    & ! Subroutine
       filltab_stilt,  & ! Subroutine
       stilt_g,        &
       stiltm_g

  use carma_fastjx, only: &
       alloc_carma_fastjx   ! Subroutine

  use mem_tuv, only : &
  	alloc_carma_tuv     ! Subroutine
!--(DMK-CCATT-FIM)-----------------------------------------------------
       
  use digitalFilter, only: 	&
  initDigitalFilter, 		& ! subroutine
  dfVars, &			  ! intent(out) - initializing
  applyDF
  
  !matrix
   USE mem_aerosol, only: &
      alloc_aerosol, &
      fill_aerosol, &
      filltab_aerosol
      
   USE memMatrix, only: &
   	aerosol
  !
  
  implicit none

  ! Argumenst:
  integer, intent(IN) :: proc_type

  ! Local Variables:
  integer, pointer :: nmzp(:), nmxp(:), nmyp(:)
  integer :: ng, nv, imean, na
  ! Local variable to control Shallow Cumulus memory allocation
  integer :: Alloc_ShCu_Flag
  ! Local variable to control Grell memory allocation
  integer :: Alloc_Grell_Flag
  integer :: ng_cp
  ! Flag to control new Grell MEmory allocation
  integer :: Alloc_Grell3_Flag
  ! Local variables because of TEB_SPM
  type(gaspart_vars), pointer :: gaspart_p
  character(len=*), parameter :: h="**(rams_mem_alloc)**"
  integer :: ierr
  integer :: ne2d, ne3d


  ! Initiating values:
  Alloc_ShCu_Flag  = 0
  Alloc_Grell_Flag = 0
  Alloc_Grell3_Flag =0 

  call alloc_paths(ngrids, nmachs)

  ! First, depending on type of process, define grid point pointers correctly..

  ! TEB_SPM
  nullify(gaspart_p)
  !

  if (proc_type==0 .or. proc_type==1) then
     !  This is the call for either a single processor run or
     !    for the master process
     nmzp => nnzp
     nmxp => nnxp
     nmyp => nnyp
  elseif (proc_type==2) then
     !  This is the call for a initial compute node process
     nmzp => nnzp
     nmxp => nodemxp(mynum,:)
     nmyp => nodemyp(mynum,:)
  elseif (proc_type==3) then
     !  This is the call for a dynamic balance compute node process
     nmzp => nnzp
     nmxp => nodemxp(mynum,:)
     nmyp => nodemyp(mynum,:)
     call dealloc_all()
  endif

  ! Call global grid dimension definitions 
  !**(JP)** To be eliminated
  call define_grid_dim_pointer(proc_type, ngrids, maxgrds, &
       nnzp, nnxp, nnyp, nnzp, nodemxp(mynum,:), nodemyp(mynum,:))

  !  If we are doing time-averaging for output, set flag ...
  imean = 0
  if (avgtim/=0.) imean = 1

  ! Allocate universal variable tables

  allocate (num_var(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating num_var")
  allocate (vtab_r(maxvars,ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating vtab_r")
  num_var(:) = 0
  nvgrids    = ngrids

  ! Allocate scalar table

  allocate(num_scalar(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating num_scalar")
  allocate(scalar_tab(maxsclr,ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating scalar_tab")
  num_scalar(:) = 0

  ! Allocate Basic variables data type
  allocate(basic_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating basic_g")
  allocate(basicm_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating basicm_g")
  do ng=1,ngrids
     call nullify_basic(basic_g(ng)); call nullify_basic(basicm_g(ng))
     call alloc_basic(basic_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
     if (imean==1) then
        call alloc_basic(basicm_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
     elseif (imean==0) then
        call alloc_basic(basicm_g(ng),        1,        1,        1, ng)
     endif

     call filltab_basic(basic_g(ng), basicm_g(ng), imean,  &
          nmzp(ng), nmxp(ng), nmyp(ng), ng)
  enddo

  ! Allocate Cuparm variables data type
  allocate(cuparm_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating cuparm_g")
  allocate(cuparmm_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating cuparmm_g")

  allocate(cuparm_g_sh(ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating cuparm_g_sh")
     allocate(cuparmm_g_sh(ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating cuparmm_g_sh")


  do ng=1,ngrids
     call nullify_cuparm(cuparm_g(ng))
     call nullify_cuparm(cuparmm_g(ng))
     call alloc_cuparm(cuparm_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)

     !-srf-feb2012: for shallow cumulus 
     if (nnshcu(ng) > 1) then
        call nullify_cuparm(cuparm_g_sh(ng))
        call nullify_cuparm(cuparmm_g_sh(ng))
        call alloc_cuparm_sh(cuparm_g_sh(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
     endif

     if (imean==1) then
     	call alloc_cuparm(cuparmm_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
        !-srf-feb2012: for shallow cumulus 
        if (nnshcu(ng) > 1) call alloc_cuparm_sh(cuparmm_g_sh(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
     elseif (imean==0) then
        call alloc_cuparm(cuparmm_g(ng), 1, 1, 1, ng)
        if (nnshcu(ng) > 1) call alloc_cuparm_sh(cuparmm_g_sh(ng), 1, 1, 1, ng)
     endif

     call filltab_cuparm(cuparm_g(ng),cuparmm_g(ng),imean,nmzp(ng),nmxp(ng),nmyp(ng),ng)

     !-srf-feb2012: for shallow cumulus 
     if (nnshcu(ng) > 1) call filltab_cuparm_sh(cuparm_g_sh(ng), cuparmm_g_sh(ng), imean, nmzp(ng), nmxp(ng), nmyp(ng), ng)

  enddo

  ! Allocate Leaf type
  allocate(leaf_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating leaf_g")
  allocate(leafm_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating leafm_g")
  do ng=1,ngrids
     call nullify_leaf(leaf_g(ng)); call nullify_leaf(leafm_g(ng))
     call alloc_leaf(leaf_g(ng), nmzp(ng), nmxp(ng), nmyp(ng),  &
          nzg, nzs, npatch, ng)
     if (imean==1) then
        call alloc_leaf(leafm_g(ng), nmzp(ng), nmxp(ng), nmyp(ng),  &
             nzg, nzs, npatch, ng)
     elseif (imean==0) then
        call alloc_leaf(leafm_g(ng),        1,        1,        1, 1, 1, 1, 1)
     endif

     call filltab_leaf(leaf_g(ng), leafm_g(ng), imean,  &
          nmzp(ng), nmxp(ng), nmyp(ng), nzg, nzs, npatch, ng)
  enddo
  ! Bob (1/10/2002) added the following line.  Is this the right place for
  ! the long term??
  call alloc_leafcol(nzg, nzs)


  ! Allocate Jules type
  allocate(jules_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating jules_g")
  allocate(julesm_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating julesm_g")
  do ng=1,ngrids
     call nullify_jules(jules_g(ng)); call nullify_jules(julesm_g(ng))
     call alloc_jules(jules_g(ng), nmzp(ng), nmxp(ng), nmyp(ng),  &
          nzg, nzs, npatch, ng)
     if (imean==1) then
        call alloc_jules(julesm_g(ng), nmzp(ng), nmxp(ng), nmyp(ng),  &
             nzg, nzs, npatch, ng)
     elseif (imean==0) then
        call alloc_jules(julesm_g(ng),        1,        1,        1, 1, 1, 1, 1)
     endif

     call filltab_jules(jules_g(ng), julesm_g(ng), imean,  &
          nmzp(ng), nmxp(ng), nmyp(ng), nzg, nzs, npatch, ng)
  enddo


  ! Allocate Micro variables data type
  allocate(micro_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating micro_g")
  allocate(microm_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating microm_g")
  do ng=1,ngrids
     call nullify_micro(micro_g(ng)); call nullify_micro(microm_g(ng))
     call alloc_micro(micro_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
     if (imean==1) then
        call alloc_micro(microm_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
     elseif (imean==0) then
        call alloc_micro(microm_g(ng),        1,        1,        1, ng)
     endif

     call filltab_micro(micro_g(ng), microm_g(ng), imean,  &
          nmzp(ng), nmxp(ng), nmyp(ng), ng)
  enddo


  ! Allocate Optimized Micro variables
  ! Only for use with SX-6 specific optimization
  if (machine==1) then
     call alloc_micro_opt(nmzp, nmxp, nmyp)
  endif

  ! Allocate radiate variables data type
  allocate(radiate_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating radiate_g")
  allocate(radiatem_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating radiatem_g")
  do ng=1,ngrids
     call nullify_radiate(radiate_g(ng)); call nullify_radiate(radiatem_g(ng))
     call alloc_radiate(radiate_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
     if (imean==1) then
        call alloc_radiate(radiatem_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
     elseif (imean==0) then
        call alloc_radiate(radiatem_g(ng),        1,        1,        1, ng)
     endif

     call filltab_radiate(radiate_g(ng), radiatem_g(ng), imean,  &
          nmzp(ng), nmxp(ng), nmyp(ng), ng)
  enddo

  ! Allocate turb variables data type
  allocate(turb_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating turb_g")
  allocate(turbm_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating turbm_g")
  do ng=1,ngrids
     call nullify_turb(turb_g(ng)); call nullify_turb(turbm_g(ng))
     call alloc_turb(turb_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
     if (imean==1) then
        call alloc_turb(turbm_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
     elseif (imean==0) then
        call alloc_turb(turbm_g(ng),        1,        1,        1, ng)
     endif

     call filltab_turb(turb_g(ng), turbm_g(ng), imean,  &
          nmzp(ng), nmxp(ng), nmyp(ng), ng)
  enddo

  ! Allocate varinit variables data type.
  !    These do not need "mean" type ever.
  allocate(varinit_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating varinit_g")
  allocate(varinitm_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating varinitm_g")
  do ng=1,ngrids
     call nullify_varinit(varinit_g(ng)); call nullify_varinit(varinitm_g(ng))
     call alloc_varinit(varinit_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
     call alloc_varinit(varinitm_g(ng),       1,        1,        1, ng)

     call filltab_varinit(varinit_g(ng), varinitm_g(ng), 0,  &
          nmzp(ng), nmxp(ng), nmyp(ng), ng)
  enddo

  ! Allocate oda variables data type.
  !    These do not need "mean" type ever.
  allocate(oda_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating oda_g")
  allocate(odam_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating odam_g")
  do ng=1,ngrids

     call nullify_oda(oda_g(ng)); call nullify_oda(odam_g(ng))

     call alloc_oda(oda_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng, proc_type)
     call alloc_oda(odam_g(ng),       1,        1,        1, ng, proc_type)

     call filltab_oda(oda_g(ng), odam_g(ng), 0,  &
          nmzp(ng), nmxp(ng), nmyp(ng), ng)

  enddo


  ! Allocate grid variables data type.
  allocate(grid_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating grid_g")
  allocate(gridm_g(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating gridm_g")
  do ng=1,ngrids
     call nullify_grid(grid_g(ng)); call nullify_grid(gridm_g(ng))
     call alloc_grid(grid_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng, if_adap)
     call alloc_grid(gridm_g(ng),       1,        1,        1, ng, if_adap)

     call filltab_grid(grid_g(ng), gridm_g(ng), 0,  &
          nmzp(ng), nmxp(ng), nmyp(ng), ng)
  enddo

  ! Allocate global grid variables data type.
  allocate(oneGlobalGridData(ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating oneGlobalGridData")
  do ng=1,ngrids
     call nullify_GlobalGridData(oneGlobalGridData(ng))
     call alloc_GlobalGridData(oneGlobalGridData(ng), nnxp(ng), nnyp(ng))
  end do

  ! Allocate any added Scalar types
  ! NOT ALLOWING DIFFERENT NUMBERS OF SCALARS ON DIFFERENT NESTS

  !   Allocate length 1 of these datatypes by default
  allocate(scalar_g(1,ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating scalar_g")
  allocate(scalarm_g(1,ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating scalarm_g")

  if (naddsc>0) then
     ! deallocate datatypes, then re-alloc to correct length
     deallocate(scalar_g, STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Deallocating scalar_g")
     deallocate(scalarm_g, STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Deallocating scalarm_g")
     allocate(scalar_g(naddsc,ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating scalar_g")
     allocate(scalarm_g(naddsc,ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating scalarm_g")
     do ng=1,ngrids
        call nullify_scalar(scalar_g(:,ng),  naddsc)
        call nullify_scalar(scalarm_g(:,ng), naddsc)
        call alloc_scalar(scalar_g(:,ng), nmzp(ng), nmxp(ng), nmyp(ng), naddsc)
        if (imean==1) then
           call alloc_scalar(scalarm_g(:,ng), nmzp(ng), nmxp(ng), nmyp(ng), &
                naddsc)
        elseif (imean==0) then
           call alloc_scalar(scalarm_g(:,ng),        1,        1,        1, &
                naddsc)
        endif

     enddo
     ! For CATT
     do ng=1,ngrids
        do na=1,naddsc ! For CATT
           call filltab_scalar(scalar_g(na,ng), scalarm_g(na,ng), imean,  &
                nmzp(ng), nmxp(ng), nmyp(ng), ng, na)
        end do
     enddo

  elseif (ilwrtyp==4 .or. iswrtyp==4) then
     ! necessary for CARMA Radiation, assume naddsc==3
     ! deallocate datatypes, then re-alloc to correct length
     ! deallocate datatypes, then re-alloc to correct length
     deallocate(scalar_g, STAT=ierr)
     if (ierr/=0) call fatal_error("Deallocating scalar_g")
     deallocate(scalarm_g, STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Deallocating scalarm_g")
     allocate(scalar_g(3,ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating scalar_g for CARMA")
     allocate(scalarm_g(3,ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating scalarm_g for CARMA")
     do ng=1,ngrids
        call nullify_scalar(scalar_g(:,ng), 3)
        call nullify_scalar(scalarm_g(:,ng), 3)
        ! assume grid 1x1xz
        call alloc_scalar(scalar_g(:,ng),  1, 1, 1, 3)
        call alloc_scalar(scalarm_g(:,ng), 1, 1, 1, 3)
        !do na=1, 3
        !call filltab_scalar(scalar_g(na,ng), scalarm_g(na,ng), imean  &
        !     , 1, 1, 1, ng, na)
        !enddo
     enddo
  endif


  ! TEB_SPM
  if (TEB_SPM==1) then
     if (isource==1) then
        !-----------------------------------------------------------------------
        ! Allocate  gaspart vars for emission
        !
        ! Defining pointers
        allocate(gaspart_g(ngrids), STAT=ierr)
        if (ierr/=0) call fatal_error(h//"Allocating gaspart_g")
        allocate(gaspartm_g(ngrids), STAT=ierr)
        if (ierr/=0) call fatal_error(h//"Allocating gaspartm_g")
        do ng=1,ngrids
           call nullify_gaspart(gaspart_g(ng))
           call nullify_gaspart(gaspartm_g(ng))
           call alloc_gaspart(gaspart_g(ng), nmzp(ng), nmxp(ng), nmyp(ng))
           if (imean==1) then
              call alloc_gaspart(gaspartm_g(ng), nmzp(ng), nmxp(ng), nmyp(ng))
           elseif (imean==0) then
              call alloc_gaspart(gaspartm_g(ng), 1,        1,        1)
           endif
           call zero_gaspart(gaspart_g(ng), nmzp(ng), nmxp(ng), nmyp(ng))
           call filltab_gaspart(gaspart_g(ng), gaspartm_g(ng), imean,  &
                nmzp(ng), nmxp(ng), nmyp(ng), ng)
        enddo
     endif
  endif

  ! Allocate Tendency data type,  filltab_tendency is responsible
  !   for filling the main atmospheric model variables in the scalar table,
  !   so make sure to call any routines that define scalar variables first.

  ! Assuming same scalars on all grids!!!!!

  call nullify_tend(naddsc)

  call alloc_tend(nmzp, nmxp, nmyp, ngrids, naddsc, proc_type)

  !Aerosol Matrix
  IF (aerosol>1) THEN
    CALL alloc_aerosol(ngrids,nmzp,nmxp,nmyp)
    CALL fill_aerosol(ngrids,nmzp,nmxp,nmyp)
    call filltab_aerosol(ngrids,imean,nmzp,nmxp,nmyp)
  END IF

  
!--(DMK-CCATT-INI)-----------------------------------------------------
  if (ccatt == 1  .and. chemistry >= 0) then
    
     ! Allocate global grid variables data type.
     allocate(oneGlobalEmissData(ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating oneGlobalEmissData")
     do ng=1,ngrids
        call nullify_GlobalEmissData(oneGlobalEmissData(ng))
        call alloc_GlobalEmissData(oneGlobalEmissData(ng), nnxp(ng), nnyp(ng))
     end do


     call define_n_dyn_chem(ngrids,dtlong,nndtrat,mynum)
    
     allocate(chem1_g    (     nspecies_chem,ngrids),chem1m_g    (     nspecies_chem,ngrids))
     allocate(chem1_src_g (max_ntimes_src,nsrc,nspecies_chem,ngrids),&
              chem1m_src_g(max_ntimes_src,nsrc,nspecies_chem,ngrids))
  
     do ng=1,ngrids

        call define_chem1_src_zdim(chem1_src_z_dim_g(:,ng),nmzp(ng))   
        call nullify_chem1(chem1_g (:,ng),chem1_src_g (:,:,:,ng), nspecies_chem) 
        call nullify_chem1(chem1m_g(:,ng),chem1m_src_g(:,:,:,ng), nspecies_chem)
     
        call alloc_chem1(chem1_g(:,ng),chem1_src_g(:,:,:,ng),chem1_src_z_dim_g(:,ng) &
                     ,nmzp(ng),nmxp(ng),nmyp(ng),nspecies_chem, ng,volcanoes) 
     
        if (imean == 1) then  
           call alloc_chem1(chem1m_g(:,ng),chem1m_src_g(:,:,:,ng),chem1_src_z_dim_g(:,ng) &
           		   ,nmzp(ng),nmxp(ng),nmyp(ng),nspecies_chem, ng,volcanoes)
        
        elseif (imean == 0) then
           call alloc_chem1(chem1m_g(:,ng),chem1m_src_g(:,:,:,ng),chem1_src_z_dim_g(:,ng) &
        		      ,1,1,1,nspecies_chem, ng,volcanoes)
        endif
        
        call filltab_chem1(chem1_g     (:,ng)	 ,chem1m_g    (:,ng)	  &
        		  ,chem1_src_g (:,:,:,ng),chem1m_src_g(:,:,:,ng)  &
        		  ,imean ,chem1_src_z_dim_g(:,ng) &
        		  ,nmzp(ng),nmxp(ng),nmyp(ng),nspecies_chem,ng,volcanoes) 
     end do

     call nullify_tend_chem1(nspecies_chem)

     call alloc_tend_chem1(nmzp,nmxp,nmyp,ngrids,nspecies_chem,proc_type)
 
     if( plumerise == 1) then
     	!-- plumerise  section 
     	allocate(plume_g     (nveg_agreg,nspecies_chem,ngrids),plumem_g     (nveg_agreg,nspecies_chem,ngrids))
     	allocate(plume_mean_g(nveg_agreg,ngrids)	      ,plume_meanm_g(nveg_agreg,ngrids))
     	
     	do ng=1,ngrids
     	
     	   call nullify_plume_chem1(plume_g  (:,:,ng),plume_mean_g  (:,ng),nspecies_chem) 
     	   call nullify_plume_chem1(plumem_g (:,:,ng),plume_meanm_g (:,ng),nspecies_chem)	
     	   call alloc_plume_chem1(plume_g(:,:,ng),plume_mean_g  (:,ng)&
     		,nmzp(ng),nmxp(ng),nmyp(ng),nspecies_chem) 
     	   if (imean == 1) then  
     	      
     	      call alloc_plume_chem1(plumem_g(:,:,ng),plume_meanm_g (:,ng)&
     		   ,nmzp(ng),nmxp(ng),nmyp(ng),nspecies_chem)	  
     	   elseif (imean == 0) then
     	   
     	      call alloc_plume_chem1(plumem_g(:,:,ng),plume_meanm_g (:,ng),1,1,1,nspecies_chem)
     	   endif
     	 
     	   call filltab_plume_chem1(plume_g   (:,:,ng),plumem_g    (:,:,ng)&
     		,plume_mean_g(:,ng)  ,plume_meanm_g (:,ng)  &
     		,imean,nmzp(ng),nmxp(ng),nmyp(ng),nspecies_chem,ng) 
     	enddo
     ENDIF
     
     !-- aerosol section
     if( aerosol == 1) then

     	   allocate(aer1_g(nmodes,nspecies_aer,ngrids),aer1m_g(nmodes,nspecies_aer,ngrids))
 	   
	   ! (DMK) dry_dep limpo-----------------------------------------
           CALL alloc_aer_sedim(npatch,ngrids, &
                              nmodes,nspecies_aer,mode_alloc, &
                              on,nmzp,nmxp,nmyp) 
           ! (DMK) dry_dep limpo-----------------------------------------

	   
     	   do ng=1,ngrids
           call define_aer1_src_zdim(aer1_src_z_dim_g(:,ng),nmzp(ng))
    	 
           call nullify_aer1(aer1_g (:,:,ng),nmodes, nspecies_aer) 
           call nullify_aer1(aer1m_g(:,:,ng),nmodes, nspecies_aer)
    	 
           call alloc_aer1(aer1_g(:,:,ng),aer1_src_z_dim_g(:,ng) &
                ,nmzp(ng),nmxp(ng),nmyp(ng),nmodes,nspecies_aer) 
    	 
           if (imean == 1) then  
              call alloc_aer1(aer1m_g(:,:,ng),aer1_src_z_dim_g(:,ng) &
    			    ,nmzp(ng),nmxp(ng),nmyp(ng),nmodes,nspecies_aer)
    	 
           elseif (imean == 0) then
              call alloc_aer1(aer1m_g(:,:,ng),aer1_src_z_dim_g(:,ng) &
    			    ,1,1,1,nmodes,nspecies_aer)
           endif
    	 
           call filltab_aer1(aer1_g(:,:,ng)	,aer1m_g(:,:,ng)       &
                ,imean ,aer1_src_z_dim_g(:,ng) &
                ,nmzp(ng),nmxp(ng),nmyp(ng),nmodes,nspecies_aer,ng) 
  
        enddo

        call nullify_tend_aer1(nmodes,nspecies_aer)

        call alloc_tend_aer1(nmzp,nmxp,nmyp,ngrids,nmodes,nspecies_aer,proc_type)
     endif
     
     !-- volcanoes section
     allocate(volc_mean_g(ngrids),volc_meanm_g(ngrids))
     if( volcanoes == 1) then
!         allocate(volc_mean_g(ngrids),volc_meanm_g(ngrids))

         do ng=1,ngrids

            call nullify_volc_chem1(volc_mean_g(ng))
            call nullify_volc_chem1(volc_meanm_g(ng))
             
            call alloc_volc_chem1(volc_mean_g(ng),nmzp(ng),nmxp(ng),nmyp(ng))

            if (imean == 1) then  
             call alloc_volc_chem1(volc_meanm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng))
            elseif (imean == 0) then
              call alloc_volc_chem1(volc_meanm_g(ng),1,1,1)
            endif

            call filltab_volc_chem1(volc_mean_g(ng),volc_meanm_g(ng)&
        			    ,imean,nmzp(ng),nmxp(ng),nmyp(ng),ng) 
         enddo
     endif
     !- end of volcanoes section ---
     
     
     !-- CHEM1AQ section
     !change MP---------------------------------------------- for CHEM1aq 
     if (chemistry_aq > 0) then
        allocate(chem1aq_g  (  nspeciesaq_chem,ngrids),chem1maq_g    (  nspeciesaq_chem,ngrids))
        allocate(chemic_g(ngrids)) 

        do ng=1,ngrids
   	 
           call nullify_chem1aq(chem1aq_g (:,ng), nspeciesaq_chem) 
           call nullify_chem1aq(chem1maq_g(:,ng), nspeciesaq_chem)
   	   call nullify_chemic(chemic_g(ng))
   	 
           call alloc_chem1aq(chem1aq_g(:,ng) &
                ,nmzp(ng),nmxp(ng),nmyp(ng),nspeciesaq_chem) 
          
	   call alloc_chemic(chemic_g(ng),nmzp(ng),nmxp(ng),nmyp(ng))
  	 
           if (imean == 1) then  
              call alloc_chem1aq(chem1maq_g(:,ng) &
                   ,nmzp(ng),nmxp(ng),nmyp(ng),nspeciesaq_chem)
 	 
           elseif (imean == 0) then
              call alloc_chem1aq(chem1maq_g(:,ng) &
                   ,1,1,1,nspeciesaq_chem)
           endif
   	 
           call filltab_chem1aq(chem1aq_g(:,ng)	    ,chem1maq_g(:,ng)	    &
                ,imean,nmzp(ng),nmxp(ng),nmyp(ng),nspeciesaq_chem,ng) 
        enddo

        call nullify_tend_chem1aq(nspeciesaq_chem)

        call alloc_tend_chem1aq(nmzp,nmxp,nmyp,ngrids,nspeciesaq_chem,proc_type) 
     endif
     !end change MP--------------------------------------------fo CHEM1 aq - end
  
  endif
!--(DMK-CCATT-FIM)-----------------------------------------------------


  do ng=1,ngrids
     ! TEB_SPM
     if (TEB_SPM==1) then
        nullify(gaspart_p)
        gaspart_p => gaspart_g(ng)
     endif

     call filltab_tend(basic_g(ng), micro_g(ng), turb_g(ng),  &
          scalar_g(:,ng),                                     &
                                ! TEB_SPM
                                !!gaspart_g(ng),                                      &
          gaspart_p,                                          &
                                !
          naddsc, ng)
	  
!--(DMK-CCATT-INI)-----------------------------------------------------
     if (ccatt == 1  .and. chemistry >= 0)  then 
         
	 call filltab_tend_chem1(nspecies_chem,ng)
	 
	 !change MP ---chem1aq
         if(chemistry_aq >= 1) call filltab_tend_chem1aq(nspeciesaq_chem,ng)
         !end change MP --chem1aq- END     

	 if (aerosol == 1)     call filltab_tend_aer1(nmodes,nspecies_aer,ng)
         	 
     endif
!--(DMK-CCATT-FIM)-----------------------------------------------------
	  
  enddo

  ! Allocate Scratch data type, This also fills the max's that are needed
  !    by nesting stuff.
  call createVctr(ngrids, nodemxp(mynum,:), nodemyp(mynum,:), nnzp)

  call nullify_scratch()

  call alloc_scratch(nmzp, nmxp, nmyp, nnzp, nnxp, nnyp, maxgrds, ngrids,  &
       nzg, nzs, npatch, proc_type, maxnxp, maxnyp, maxnzp)

  ! For CATT - LFR
!--(DMK-CCATT-INI)-----------------------------------------------------
   if (CCATT==1) then
!--(DMK-CCATT-OLD)-----------------------------------------------------
!  if (CATT==1) then
!--(DMK-CCATT-FIM)-----------------------------------------------------

     allocate(turb_s(ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating turb_s")

     do ng=1,ngrids
        call nullify_turb_s(turb_s(ng))
        call alloc_turb_s(turb_s(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
        call filltab_turb_s(turb_s(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
     enddo

  endif

  ! Reproducibility - Saulo Barros
  call nullify_scratch1()
  call alloc_scratch1(nodebounds, maxgrds, ngrids, nnzp, mynum)

  ! For optmization - ALF
  call nullify_opt_scratch()
  if ((if_adap==0) .and. (ihorgrad==2)) &
       call alloc_opt_scratch(proc_type, ngrids, nnzp, nnxp, nnyp, 1000, 1000)

  ! Allocate nested boundary interpolation arrays. All grids will be allocated.
  ! Changed by Alvaro L.Fazenda
  ! To correct a problem when running in a NEC SX-6
  ! Master process needs allocation for nesting in a parallel run
  if (proc_type==0 .or. proc_type==2 .or. proc_type==1) then
     do ng=1,ngrids
        if (nxtnest(ng)==0 ) then
           call alloc_nestb(ng,        1,        1,        1)
        else
           call alloc_nestb(ng, nnxp(ng), nnyp(ng), nnzp(ng))
        endif
     enddo
  endif

  !--------------------------------------------------------------------------
  ! Allocate data for Shallow Cumulus
  ! Verifying if the allocation is necessary in any grids
  do ng=1, ngrids
     if (NNSHCU(ng)==1) Alloc_ShCu_Flag = 1
  enddo
  !Alloc_ShCu_Flag = 1  !!Provisorio
  if (Alloc_ShCu_Flag==1) then
     allocate(shcu_g(ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating shcu_g")
     allocate(shcum_g(ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating shcum_g")
     do ng=1,ngrids
        call nullify_shcu(shcu_g(ng))
        call nullify_shcu(shcum_g(ng))
        call alloc_shcu(shcu_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
        if (imean==1) then
           call alloc_shcu(shcum_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
        elseif (imean==0) then
           call alloc_shcu(shcum_g(ng),        1,        1,        1, ng)
        endif

        call filltab_shcu(shcu_g(ng), shcum_g(ng), imean,  &
             nmzp(ng), nmxp(ng), nmyp(ng), ng)
     enddo
  endif

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !-
  !-srf-fev-2012 - reformulating cumulus (deep/shallow) convection options 2 and above
  !--------------------------------------------------------------------------
  ! Allocate data for Shallow/Deep Cumulus 
  !
  Alloc_SHCU_Flag   =0
  Alloc_Grell_Flag  =0
  Alloc_Grell3_Flag =0 
  do ng=1, ngrids
     if (NNSHCU (ng) == 2) Alloc_SHCU_Flag   = 1
     if (NNQPARM(ng) == 2) Alloc_Grell_Flag  = 1
     if (NNQPARM(ng) >  2) Alloc_Grell3_Flag = 1
  enddo

  if (Alloc_Grell_Flag == 1 .or. Alloc_SHCU_Flag == 1 .or. Alloc_Grell3_Flag == 1) then
     ! Calculating the necessary space for scratch data
     call define_memory(nmxp, nmyp, nmzp, ngrids, nnqparm, nnshcu)
     ! Allocating data for scratch data
     call alloc_scratch1_grell()
    !remove call zero_scratch1_grell() 
    !call alloc_scratch2_grell()
    !call alloc_scratch3_grell()
     
      allocate(cuforc_g   (ngrids_cp),cuforcm_g   (ngrids_cp)) !usar ierr do allocate
      allocate(cuforc_sh_g(ngrids_cp),cuforcm_sh_g(ngrids_cp))
      do ng=1,ngrids_cp
           call nullify_cuforc(cuforc_g    (ng))
           call nullify_cuforc(cuforcm_g   (ng))
           call nullify_cuforc(cuforc_sh_g (ng))
           call nullify_cuforc(cuforcm_sh_g(ng))
	   call alloc_cu_forcings(cuforc_g   (ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
	   call alloc_cu_forcings(cuforc_sh_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
	   if (imean == 1) then
              call alloc_cu_forcings(cuforcm_g   (ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
	      call alloc_cu_forcings(cuforcm_sh_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)

	    elseif (imean == 0) then
              call alloc_cu_forcings(cuforcm_g   (ng),1,1,1,ng)
 	      call alloc_cu_forcings(cuforcm_sh_g(ng),1,1,1,ng)

           endif
	   call filltab_cuforc_sh(cuforc_sh_g(ng),cuforcm_sh_g(ng),imean  &
                                 ,nmzp(ng),nmxp(ng),nmyp(ng),ng)
	   
	   call filltab_cuforc   (cuforc_g(ng),cuforcm_g(ng),imean  &
                                 ,nmzp(ng),nmxp(ng),nmyp(ng),ng)
      enddo
  endif
     
     
  !- shallow convection version 2
  if  (Alloc_SHCU_Flag  == 1) then
        call alloc_scratch2_grell_sh()
        call alloc_scratch3_grell_sh()
	
	allocate(grell_g_sh(ngrids_cp), STAT=ierr)
        if (ierr/=0) call fatal_error(h//"Allocating grell_g_sh_g")
        allocate(grellm_g_sh(ngrids_cp), STAT=ierr)
        if (ierr/=0) call fatal_error(h//"Allocating grell_g_sh_g")
	
	do ng=1,ngrids
        if (NNSHCU(ng) == 2) then
           call nullify_grell(grell_g_sh (ng))
           call nullify_grell(grellm_g_sh(ng))
          
	   call alloc_grell_sh(grell_g_sh(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
          
	   if (imean == 1) then
  	      call alloc_grell_sh(grellm_g_sh(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)

           elseif (imean == 0) then
              call alloc_grell_sh(grellm_g_sh(ng),1,1,1,ng)

           endif
           
	   call filltab_grell_sh(grell_g_sh(ng),grellm_g_sh(ng),imean  &
                   ,nmzp(ng),nmxp(ng),nmyp(ng),ng)
         endif
      enddo
   endif

!--------------------------------------------------------------------------
  ! Allocate data for Grell's deep cumulus - old GD
  !
  if (Alloc_Grell_Flag == 1) then
     
     call alloc_scratch2_grell()
     call alloc_scratch3_grell()

     ! Allocating data for main Grell data
     ! Allocate to a quantity set by ngrids_cp not ngrids
     allocate(grell_g(ngrids_cp), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating grell_g")
     allocate(grellm_g(ngrids_cp), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating grellm_g")
     
     ng_cp = 1
     do ng=1,ngrids
        if (NNQPARM(ng) == 2) then
           call nullify_grell(grell_g (ng))
           call nullify_grell(grellm_g(ng))

           call alloc_grell(grell_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
	
           if (imean == 1) then
              call alloc_grell(grellm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
           elseif (imean == 0) then
              call alloc_grell(grellm_g(ng),1,1,1,ng)
           endif

           call filltab_grell(grell_g(ng),grellm_g(ng),imean  &
                ,nmzp(ng),nmxp(ng),nmyp(ng),ng)
           ng_cp = ng_cp + 1
        endif
     enddo

    ! call zero_scratch3_grell() 

     Flag_Grell = 2
     ! For New Grell Param.
     if     (CLOSURE_TYPE == 'EN') then
        icoic = 0
     elseif (CLOSURE_TYPE == 'GR') then
        icoic = 1
     elseif (CLOSURE_TYPE == 'LO') then
        icoic = 4
     elseif (CLOSURE_TYPE == 'MC') then
        icoic = 7
     elseif (CLOSURE_TYPE == 'SC') then
        icoic = 10
     elseif (CLOSURE_TYPE == 'AS') then
        icoic = 13
     else
        print *, "****Grell Closure type ERROR!!"
        ! the subroutine opspec3 stop the program before this point.
     endif

     ! For CATT

!--(DMK-CCATT-INI)-----------------------------------------------------
     if (CCATT == 1) then
!--(DMK-CCATT-OLD)-----------------------------------------------------
!     if (CATT == 1) then
!--(DMK-CCATT-FIM)-----------------------------------------------------

        icoic_sh=icoic !LFR
     endif

  endif

 !--------------------------------------------------------------------------
  ! Allocate data for Grell Cumulus version 3d and FIM

  if (Alloc_Grell3_Flag == 1) then
     ! Calculating the necessary space for scratch data
     !call define_memory(nmxp, nmyp, nmzp, ngrids, nnqparm)
     ! Allocating data for scratch data
     !call alloc_scratch1_grell
     !call zero_scratch1_grell
     
     !- for shallow convection
     !if (Alloc_SHCU_Flag==1) then
     !   call alloc_scratch2_grell_sh 
     !   call alloc_scratch3_grell_sh 
     !endif

     ! Allocating data for main Grell data
     ! Allocate to a quantity set by ngrids_cp not ngrids
     !allocate(grell_g(ngrids_cp),grellm_g(ngrids_cp))
     ! These arrays must be allocated anyway for the forcing
     !allocate(grell_g_sh(ngrids_cp),grellm_g_sh(ngrids_cp))

     allocate(g3d_ens_g(train_dim,ngrids_cp),g3d_ensm_g(train_dim,ngrids_cp)) ! arrumar iErr aqui
     allocate(g3d_g(ngrids_cp),g3dm_g(ngrids_cp))

     ng_cp = 1
     do ng=1,ngrids
     	if (NNQPARM(ng) > 2) then
           !call nullify_grell(grell_g (ng_cp))
           !call nullify_grell(grellm_g(ng_cp))

           !call nullify_grell(grell_g_sh (ng_cp))
           !call nullify_grell(grellm_g_sh(ng_cp))


           !call alloc_grell   (grell_g   (ng_cp),nmzp(ng),nmxp(ng),nmyp(ng),ng)
           !call alloc_grell_sh(grell_g_sh(ng_cp),nmzp(ng),nmxp(ng),nmyp(ng),ng)

           !if (imean == 1) then
           !   call alloc_grell   (grellm_g   (ng_cp),nmzp(ng),nmxp(ng),nmyp(ng),ng)
	   !   call alloc_grell_sh(grellm_g_sh(ng_cp),nmzp(ng),nmxp(ng),nmyp(ng),ng)

           !elseif (imean == 0) then
           !   call alloc_grell   (grellm_g   (ng_cp),1,1,1,ng)
           !   call alloc_grell_sh(grellm_g_sh(ng_cp),1,1,1,ng)

           !endif
           !call filltab_grell   (grell_g(ng_cp),grellm_g(ng_cp),imean  &
           !                     ,nmzp(ng),nmxp(ng),nmyp(ng),ng)
           !print*,'call filltab_grell_sh'; call flush(6)
	   !call filltab_grell_sh(grell_g_sh(ng_cp),grellm_g(ng_cp),imean  &
           !                     ,nmzp(ng),nmxp(ng),nmyp(ng),ng)

          !-- arrays needed for G3d and GD-FIM

           call nullify_grell3(g3d_ens_g (:,ng_cp) , g3d_g(ng_cp),train_dim)
           call nullify_grell3(g3d_ensm_g(:,ng_cp) ,g3dm_g(ng_cp),train_dim)

           call alloc_grell3(g3d_ens_g(:,ng_cp),g3d_g(ng_cp),nmzp(ng),nmxp(ng),nmyp(ng),ng,train_dim)
           
           if (imean == 1) then
           	call alloc_grell3(g3d_ensm_g(:,ng_cp),g3dm_g(ng_cp),nmzp(ng),nmxp(ng),nmyp(ng),ng,train_dim)
	   
          elseif (imean == 0) then
 	   
          	call alloc_grell3(g3d_ensm_g(:,ng_cp),g3dm_g(ng_cp),1,1,1,ng,train_dim)
          endif
	   	   
	  call filltab_grell3(g3d_ens_g(:,ng_cp),g3d_g(ng_cp),g3d_ensm_g(:,ng_cp), &
                              g3dm_g(ng_cp),imean,nmzp(ng),nmxp(ng),nmyp(ng),ng,train_dim)
			      
          ng_cp = ng_cp + 1
      	endif
    enddo
 endif

!-srf-g3d-end--------------------------------------------------------
!-srf-g3d-end--------------------------------------------------------
!-srf-g3d-end--------------------------------------------------------











     
     
     
     
     
     
     
     
     





  ! Allocation for SiB
  if (ISFCL==3) then
     allocate(sib_g(ngrids, N_CO2), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating sib_g")
     allocate(sibm_g(ngrids, N_CO2), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating sibm_g")
     allocate(sib_brams_g(ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating sib_brams_g")
     allocate(sib_bramsm_g(ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating sib_brams_g")
     do ng=1,ngrids
        ! SiB CO2
        call nullify_sib_co2(sib_g(ng, 1))
        call nullify_sib_co2(sibm_g(ng, 1))
        call alloc_sib_co2(sib_g(ng, 1), nmxp(ng), nmyp(ng))
        ! Putting zero on all values
        call zero_sib_co2(sib_g(ng, 1), nmxp(ng), nmyp(ng))
        if (imean==1) then
           call alloc_sib_co2(sibm_g(ng, 1), nmxp(ng), nmyp(ng))
           ! Putting zero on all values
           call zero_sib_co2(sibm_g(ng, 1),  nmxp(ng), nmyp(ng))
        elseif (imean==0) then
           call alloc_sib_co2(sibm_g(ng, 1), 1, 1)
           ! Putting zero on all values
           call zero_sib_co2(sibm_g(ng, 1),  1, 1)
        endif
        call filltab_sib_co2(sib_g(ng, 1), sibm_g(ng, 1), imean,  &
             nmxp(ng), nmyp(ng), ng)
        ! Putting zero on all values
        ! SiB BRAMS types
        call nullify_sib_brams(sib_brams_g(ng))
        call nullify_sib_brams(sib_bramsm_g(ng))
        call alloc_sib_brams(sib_brams_g(ng), nmxp(ng), nmyp(ng))
        ! Putting zero on all values
        call zero_sib_brams(sib_brams_g(ng), nmxp(ng), nmyp(ng))
        if (imean==1) then
           call alloc_sib_brams(sib_bramsm_g(ng), nmxp(ng), nmyp(ng))
           ! Putting zero on all values
           call zero_sib_brams(sib_bramsm_g(ng),  nmxp(ng), nmyp(ng))
        elseif (imean==0) then
           call alloc_sib_brams(sib_bramsm_g(ng), 1, 1)
           ! Putting zero on all values
           call zero_sib_brams(sib_bramsm_g(ng),  1, 1)
        endif
        call filltab_sib_brams(sib_brams_g(ng), sib_bramsm_g(ng), imean,  &
             nmxp(ng), nmyp(ng), ng)
     enddo
  endif
!--(DMK-CCATT-INI)-----------------------------------------------------
  
  !srf- definir a configura��o que isto � necess�rio
  ![ML
  !Allocate stilt variables data type
  allocate(stilt_g(ngrids), stiltm_g(ngrids))
  do ng=1, ngrids
    call nullify_stilt(stilt_g(ng)); call nullify_stilt(stiltm_g(ng))
    call alloc_stilt(stilt_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
    if (imean == 1) then
      call alloc_stilt(stiltm_g(ng),nmzp(ng),nmxp(ng),nmyp(ng),ng)
    else if (imean == 0) then
      call alloc_stilt(stiltm_g(ng),1,1,1,ng)
    end if  
    call filltab_stilt(stilt_g(ng),stiltm_g(ng),imean, &
                       nmzp(ng),nmxp(ng),nmyp(ng),ng)
  end do
  !ML]
!--(DMK-CCATT-END)-----------------------------------------------------

  ! CATT - Allocation for Transporte
  ! If not used just allocate the array but not the item of the structure
  allocate(extra2d(na_extra2d,ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating extra2d")
  allocate(extra3d(na_extra3d,ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating extra3d")
  allocate(extra2dm(na_extra2d,ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating extra2dm")
  allocate(extra3dm(na_extra3d,ngrids), STAT=ierr)
  if (ierr/=0) call fatal_error(h//"Allocating extra3dm")
  call nullify_extra2d(extra2d,  na_extra2d, ngrids)
  call nullify_extra2d(extra2dm, na_extra2d, ngrids)
  call nullify_extra3d(extra3d,  na_extra3d, ngrids)
  call nullify_extra3d(extra3dm, na_extra3d, ngrids)
  ! Allocating itens of the structure
  do ng=1,ngrids
     call alloc_extra2d(extra2d, nmxp(ng), nmyp(ng), na_extra2d, ng)
     call zero_extra2d(extra2d, na_extra2d, ng)
     if (imean==1) then
        call alloc_extra2d(extra2dm, nmxp(ng), nmyp(ng), na_extra2d, ng)
        call zero_extra2d(extra2dm, na_extra2d, ng)
     elseif (imean==0) then
        call alloc_extra2d(extra2dm, 1, 1, na_extra2d, ng)
        call zero_extra2d(extra2dm, na_extra2d, ng)
     end if
     call alloc_extra3d(extra3d, nmzp(ng), nmxp(ng), nmyp(ng), na_extra3d, ng)
     call zero_extra3d(extra3d, na_extra3d, ng)
     if (imean==1) then
        call alloc_extra3d(extra3dm,  &
             nmzp(ng), nmxp(ng), nmyp(ng), na_extra3d, ng)
        call zero_extra3d(extra3dm, na_extra3d, ng)
     elseif (imean==0) then
        call alloc_extra3d(extra3dm, 1, 1, 1, na_extra3d, ng)
        call zero_extra3d(extra3dm, na_extra3d, ng)
     end if
  end do
  do ng=1,ngrids
     do na=1,na_extra2d
        call filltab_extra2d(extra2d(na,ng), extra2dm(na,ng), imean, &
             nmxp(ng), nmyp(ng), ng, na)
     end do
     do na=1,na_extra3d
        call filltab_extra3d(extra3d(na,ng), extra3dm(na,ng), imean, &
             nmzp(ng), nmxp(ng), nmyp(ng), ng, na)
     end do
  end do

  ! Carma Radiation
  if (ilwrtyp==4 .or. iswrtyp==4) then
     call initial_definitions_aerad()
     call initial_definitions_globrad()
     call initial_definitions_globaer()

!--(DMK-CCATT-INI)-----------------------------------------------------
     if (trim(PhotojMethod) == 'FAST-JX') &
	  CALL alloc_carma_fastjx(maxval(nmzp(1:ngrids))&
          ,maxval(nmxp(1:ngrids))&
          ,maxval(nmyp(1:ngrids)))
     if (trim(PhotojMethod) == 'FAST-TUV' .and. chemistry > 0) &
          call alloc_carma_tuv()
!--(DMK-CCATT-END)-----------------------------------------------------

     allocate(carma(ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating carma")
     allocate(carma_m(ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating carma_m")
     do ng=1,ngrids
        call nullify_carma(carma,ng)
        call alloc_carma(carma, ng, nmxp(ng), nmyp(ng), nwave)
        call zero_carma(carma, ng)
        call nullify_carma(carma_m, ng)
        if(imean==1) then
           call alloc_carma(carma_m, ng, nmxp(ng), nmyp(ng), nwave)
           call zero_carma(carma_m, ng)
        elseif (imean==0) then
           call alloc_carma(carma_m, ng,        1,        1, nwave)
           call zero_carma(carma_m, ng)
        end if

        call filltab_carma(carma(ng), carma_m(ng), ng, imean,  &
             nmxp(ng), nmyp(ng), nwave)
     end do
  endif

  if (TEB_SPM==1) then
     !---------------------------------------------------------------------
     ! Allocate common use variables (TEB_SPM,LEAF)
     allocate(tebc_g(ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating tebc_g")
     allocate(tebcm_g(ngrids), STAT=ierr)
     if (ierr/=0) call fatal_error(h//"Allocating tebcm_g")
     do ng=1,ngrids
        call nullify_tebc(tebc_g(ng))
        call nullify_tebc(tebcm_g(ng))
        call alloc_tebc(tebc_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
        if (imean==1) then
           call alloc_tebc(tebcm_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
        elseif (imean==0) then
           call alloc_tebc(tebcm_g(ng),        1,        1,        1, ng)
        endif

        call filltab_tebc(tebc_g(ng), tebcm_g(ng), imean,  &
             nmzp(ng), nmxp(ng), nmyp(ng), ng)
     enddo

     !---------------------------------------------------------------------
     ! Allocate data for urban canopy parameterization
     if (iteb==1) then
        allocate(teb_g(ngrids), STAT=ierr)
        if (ierr/=0) call fatal_error(h//"Allocating teb_g")
        allocate(tebm_g(ngrids), STAT=ierr)
        if (ierr/=0) call fatal_error(h//"Allocating tebm_g")
        do ng=1,ngrids
           call nullify_teb(teb_g(ng))
           call nullify_teb(tebm_g(ng))
           call alloc_teb(teb_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
           if (imean==1) then
              call alloc_teb(tebm_g(ng), nmzp(ng), nmxp(ng), nmyp(ng), ng)
           elseif (imean==0) then
              call alloc_teb(tebm_g(ng),        1,        1,        1, ng)
           endif

           call filltab_teb(teb_g(ng), tebm_g(ng), imean,  &
                nmzp(ng), nmxp(ng), nmyp(ng), ng)
        enddo
     endif
  endif
  
  !---------------------------------------------------------------------------
  ! allocate date for digital filter - rmf

	if(applyDF) call initDigitalFilter(dfVars, ngrids, nmzp, nmxp, nmyp)

  !--------------------------------------------------------------------------

  ! Set "Lite" variable flags according to namelist input LITE_VARS.
  if (proc_type==0 .or. proc_type==2 .or. proc_type==1) then
     call lite_varset(proc_type)
  endif

  ! Set ALL variables in the vtab_r variable table to zero by default. These
  !  are variables processed in the filltab_* routines with a call to vtables2.
  !  This does NOT include scratch arrays, tendencies, or mean arrays.
  do ng=1,ngrids
     do nv=1,num_var(ng)
        call azero_l(vtab_r(nv,ng)%npts, vtab_r(nv,ng)%var_p)
     enddo
  enddo

end subroutine rams_mem_alloc
