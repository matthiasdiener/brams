#Makefile objects_model.mk

# Define main source.

MAIN = $(MODEL)/rammain.f90
MAINOBJ = rammain.o


# Define objects.

OBJ_MODEL = \
        node_mod.o    \
	ModVarfFile.o \
	ModBuffering.o \
	ModMessageData.o \
	ModMessagePassing.o \
	ModMessageSet.o \
	ModGridDims.o \
	ModGrid.o \
	ModGridTree.o \
	ModNeighbourNodes.o \
	ModFieldSectionList.o \
	ModDomainDecomp.o \
	ModDateUtils.o \
	ReadBcst.o \
	mpi_io_engine-5d.o \
	hdf5_parallel_engine.o \
	ModUkmoAdapt.o \
	Rad_UKMO.o \
	Phys_const.o \
	sfclyr_jules.o \
	alloc.o   \
	an_header.o \
	aobj.o \
	asgen.o \
	asti.o \
	asti2.o \
	astp.o \
	avarf.o \
	catt_start.o \
	ccatt_start.o \
	emission_source_map.o \
	machine_arq.o \
	teb_spm_start.o \
	mem_grid_dim_defs.o \
	cond_read.o \
	cond_update.o \
	conv_coms.o \
	coriolis.o \
	cu_read.o \
	cup_dn.o \
	cup_env.o \
	cup_env_catt.o \
	cup_grell_catt_deep.o \
	cup_grell_catt_shallow.o \
	cup_output_vars.o \
	cup_up.o \
	diffsclr.o \
	diffuse.o \
	extra.o  \
	file_inv.o \
	first_rams.o \
	geodat.o \
	grid_dims.o \
	grid_struct.o \
	inithis.o  \
	io_params.o \
	isan_coms.o \
	isan_io.o \
	ke_coms.o \
	landuse_input.o \
	leaf3.o \
	leaf3_hyd.o \
	leaf3_init.o \
	leaf_coms.o \
	leaf3_teb.o \
	mem_aerad.o \
	mem_all.o \
	mem_basic.o \
	mem_carma.o \
	mem_cuparm.o \
	mem_cutrans.o \
	mem_globaer.o \
	mem_globrad.o \
	mem_grell.o \
	mem_grell_param2.o \
	mem_grid.o \
	mem_leaf.o \
	mem_jules.o \
	mem_micro.o \
	mem_mksfc.o \
	mem_nestb.o \
	mem_oda.o \
	mem_opt_scratch.o \
	mem_precision.o \
	mem_radiate.o \
	mem_scalar.o \
	mem_scratch.o \
	mem_scratch1_brams.o \
	mem_scratch1_grell.o \
	mem_scratch2_grell.o \
	mem_scratch2_grell_sh.o \
	mem_scratch3_grell.o \
	mem_scratch3_grell_sh.o \
	mem_shcu.o \
	mem_sib.o \
	mem_sib_co2.o \
	mem_tconv.o \
	mem_tend.o \
	mem_turb.o \
	mem_turb_scalar.o \
	mem_varinit.o \
	mem_micro_optij.o \
	mic_coll.o \
	mic_driv.o \
	mic_driv_new.o \
	mic_gamma.o \
	mic_init.o \
	mic_misc.o \
	mic_nuc.o \
	mic_tabs.o \
	mic_vap.o \
	micphys.o \
	mksfc_driver.o \
	mksfc_ndvi.o \
	mksfc_sfc.o \
	mksfc_sst.o \
	mksfc_top.o \
	ModTimeStamp.o \
	modsched.o \
	local_proc.o \
	mpass_dtl.o \
	mpass_feed.o \
	mpass_full.o      \
	mpass_nest.o \
	ndvi_read.o \
	nest_drivers.o \
	nest_feed.o \
	nest_filldens.o \
	nest_geosst.o \
	nest_intrp.o \
	nud_analysis.o \
	nud_read.o \
	nud_update.o \
	obs_input.o \
	oda_krig.o \
	oda_nudge.o \
	oda_proc_obs.o \
	oda_read.o \
	oda_sta_count.o \
	oda_sta_input.o \
	opspec.o \
	domain_decomp.o \
	para_init.o \
	Phys_const.o \
	raco.o \
	raco_adap.o \
	rad_carma.o \
	dry_dep.o \
	rad_ccmp.o \
	rad_driv.o \
	rad_mclat.o \
	rad_stable.o \
	radvc.o \
	radvc_adap.o \
        radvc_mnt.o \
	mod_GhostBlock.o \
	mod_GhostBlockPartition.o \
	mod_advect_kit.o \
	radvc_new.o \
	rams_grid.o \
	gridset.o \
	adap_init.o \
	ModOneProc.o \
	rams_mem_alloc.o   \
	rams_read_header.o \
	ranlavg.o \
	rbnd.o \
	rbnd_adap.o \
	rcio.o \
	rconstants.o \
	rconv.o \
	rconv_grell_catt.o \
        chem_conv_transp.o \
	ModNamelistFile.o  \
	ModParallelEnvironment.o  \
	read_ralph.o \
	recycle.o  \
	ref_sounding.o \
	refstate.o \
	rgrad.o \
	rhhi.o  \
	rinit.o  \
	rio.o \
	rnest_par.o \
	rnode.o \
	rrad2.o \
	rrad3.o \
	rshcupar.o \
	rthrm.o \
	rtimh.o \
	rtimi.o \
	ruser.o \
	shcu_vars_const.o \
	sib2_co2.o \
	sib2_driver.o \
	sib2_init.o \
	sib_vars.o \
	memSoilMoisture.o \
	soilMoisture.o \
	sst_read.o \
	turb_diff.o \
	turb_diff_adap.o \
	turb_k.o \
	turb_k_adap.o \
	turb_ke.o \
	upcase.o \
	urban_canopy.o \
	v_interps.o \
	ccatt_extra.o \
	mem_stilt.o \
	aer1_list.o \
	mem_aer1.o \
	chem1_list.o \
	chem1aq_list.o \
	mem_chem1aq.o \
	mem_chem1.o \
	mem_chemic.o \
	mem_plume_chem1.o \
	mem_volc_chem1.o \
	chem_sources.o \
	chem_plumerise_scalar.o \
	ChemSourcesDriver.o \
	chem_dry_dep.o \
	ChemDryDepDriver.o \
	mem_tuv.o \
	tuvParameter.o \
	ModTuv2.7.o \
	ModTuvDriver2.7.o \
	var_tables.o \
	varf_update.o \
	vtab_fill.o \
	plumerise_vector.o \
	mksfc_fuso.o \
	mem_teb.o \
	mem_teb_common.o \
	mem_teb_vars_const.o \
	mem_gaspart.o \
	mem_emiss.o \
	urban.o \
	gaspart.o \
	ozone.o \
	mod_ozone.o \
	chem_isan_coms.o \
	chem_aobj.o \
	chem_asgen.o \
	chem_asti2.o \
	chem_asti.o \
	chem_astp.o \
	chem_avarf.o \
	chem_file_inv.o \
	chem_first_rams.o \
	chem_isan_io.o \
	chem_refstate.o \
	chem_v_interps.o \
	carma_fastjx.o \
	chem_fastjx57.o \
	chem_fastjx_data.o \
	chem_fastjx_driv.o \
	chem_spack_utils.o \
	mem_spack.o \
	chem_uv_att.o\
	chem_dry_dep.o \
	chem_spack_solve_sparse.o \
	chem_spack_lu.o \
	chem_spack_jacdchemdc.o \
	chem_spack_dratedc.o \
	chem_spack_kinetic.o \
	chem_spack_fexprod.o \
	chem_spack_fexloss.o \
	chem_spack_rates.o \
	chem_spack_fexchem.o \
	chem_spack_ros.o \
	chem_spack_ros_dyndt.o \
	chem_spack_rodas3_dyndt.o\
	chem_spack_qssa.o \
	chem_trans_gasaq.o \
	chem_trans_liq.o \
	chem_orage.o \
	chemistry.o \
	ModPostProcess.o \
	ModPostOneField.o \
	ModPostOneFieldUtils.o \
	ModPostUtils.o \
	ModPostGrid.o \
	ModBramsGrid.o \
	ModOutputUtils.o \
	kbcon_ecmwf.o \
	module_cu_g3.o \
        module_cu_gf.o \
	module_cu_gd_fim.o \
	cup_grell3.o \
	rexev.o \
	rstilt.o \
	turb_constants.o \
	tkenn.o \
	digitalFilter.o \
	GridMod.o \
	MapMod.o \
	ProcessorMod.o \
	BoundaryMod.o \
	errorMod.o \
	advSendMod.o \
	InitAdvect.o \
	isrpia.o \
	quad.o \
	actv.o \
	coag.o \
	depv.o \
	diam.o \
	dicrete.o \
	init.o \
        solut.o \
	issoropia.o \
	isofwd.o \
	isorev.o \
	matrix.o \
	npf.o \
	quad.o \
	setup.o \
	subs.o \
	thermo_isorr.o \
	isrpia.o \
	MatrixDriver.o \
        ModParticle.o \
        mem_aerosol.o \
        memMatrix.o \
	seasalt.o  \
	meteogram.o \
	meteogramType.o
