module ModPostOneFieldUtils

  use ModOutputUtils, only: GetVarFromMemToOutput

  use ModPostGrid, only: OutputGradsField
  use ModBramsGrid, only: BramsGrid
  use ModPostGrid, only: PostGrid

  use ModPostUtils, only: rams_comp_tempk
  use ModPostUtils, only: rams_comp_tempc
  use ModPostUtils, only: rams_comp_dewk
  use ModPostUtils, only: rams_comp_thetv
  use ModPostUtils, only: rams_get_surface
  use ModPostUtils, only: rams_comp_slpmm5
  use ModPostUtils, only: rams_comp_rh
  use ModPostUtils, only: rams_comp_press
  use ModPostUtils, only: cape_cine
  use ModPostUtils, only: rams_comp_dn0
  use ModPostUtils, only: rams_comp_z
  use ModPostUtils, only: rams_comp_rotate
  use ModPostUtils, only: rams_comp_avgu
  use ModPostUtils, only: rams_comp_avgv
  use ModPostUtils, only: rams_comp_speed
  use ModPostUtils, only: rams_reduced_temp
  use ModPostUtils, only: rams_fill_sst
  use ModPostUtils, only: get_leaf_soil
  use ModPostUtils, only: rams_comp_pbl
  use ModPostUtils, only: calc_omeg
  use ModPostUtils, only: comp_vertint
   use ModPostUtils, only: comp_vertint_press
  use ModPostUtils, only: comp_slp_metar
  use ModPostUtils, only: rams_reduced_rv
  use ModPostUtils, only: rams_comp_dewk_2m
  use ModPostUtils, only: rams_reduced_wind
  use ModPostUtils, only: rams_comp_dir
  use ModPostUtils, only: calc_u10m
  use ModPostUtils, only: calc_v10m
  use ModPostUtils, only: rams_comp_vegclass
  
  implicit none

  private

!--(DMK-CCATT-INI)-------------------------------------------------------
  real, parameter    :: PMAR  = 28.96
!--(DMK-CCATT-FIM)-------------------------------------------------------

  public :: Brams2Post_u 
  public :: Brams2Post_v 
  public :: Brams2Post_tempk 
  public :: Brams2Post_tveg 
  public :: Brams2Post_totpcp 
  public :: Brams2Post_acccon 
  public :: Brams2Post_dewptc 
  public :: Brams2Post_rshort 
  public :: Brams2Post_rlong 
  public :: Brams2Post_sea_press 
  public :: Brams2Post_cape 
  public :: Brams2Post_cine 
  public :: Brams2Post_topo 
  public :: Brams2Post_precip 
  public :: Brams2Post_le 
  public :: Brams2Post_h 
  public :: Brams2Post_rv 
  public :: Brams2Post_rlongup 
  public :: Brams2Post_albedt 
  public :: Brams2Post_geo 
  public :: Brams2Post_ue_avg 
  public :: Brams2Post_ve_avg 
  public :: Brams2Post_tempc2m 
  public :: Brams2Post_tempc 
  public :: Brams2Post_rh 
  public :: Brams2Post_w 
  public :: Brams2Post_sst 
  public :: Brams2Post_land 
  public :: Brams2Post_smoist 
  public :: Brams2Post_zi 
  public :: Brams2Post_tke 
  public :: Brams2Post_cloud
  public :: Brams2Post_omeg
  public :: Brams2Post_pwt
  public :: Brams2Post_slp_metar
  public :: Brams2Post_td2m
  public :: Brams2Post_u10m
  public :: Brams2Post_v10m
  public :: Brams2Post_vtype 
  public :: Brams2Post_sltex_p
  public :: Brams2Post_ndvi
  public :: Brams2Post_lai
  public :: Brams2Post_u10mj
  public :: Brams2Post_v10mj
  public :: Brams2Post_t2mj
  public :: Brams2Post_csj
  public :: Brams2Post_rv2mj
  public :: Brams2Post_td2mj
  public :: Brams2Post_co2
  public :: Brams2Post_co2_antro
  public :: Brams2Post_co2_burn
  public :: Brams2Post_co2_burn2D
  public :: Brams2Post_co2_bioge
  public :: Brams2Post_co2_total
  public :: Brams2Post_theta

!--(DMK-CCATT-INI)-------------------------------------------------------
  public :: Brams2Post_CO
  public :: Brams2Post_NO
  public :: Brams2Post_HNO3
  public :: Brams2Post_O3
  public :: Brams2Post_NO2
  public :: Brams2Post_PM25
  public :: Brams2Post_PM25wd
  public :: Brams2Post_AOT550
  public :: Brams2Post_PMINT
!--(DMK-CCATT-FIM)-------------------------------------------------------
  public :: Brams2Post_cuthdp
  public :: Brams2Post_curtdp
  public :: Brams2Post_cucldp
     
contains

  subroutine Brams2Post_u (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'UP' from var_table

    call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'u'
    onePostGrid%fieldUnits = 'm/s'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_u



  subroutine Brams2Post_v (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'VP' from var_table

    call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'v'
    onePostGrid%fieldUnits = 'm/s'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_v



  subroutine Brams2Post_tempk (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'THETA' from var_table
    ! ScrT3N01 <- field 'PI' from var_table
    ! OutputField <- rams_comp_tempk (OutputField, ScrT3N01)

    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, OutputField)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
    call rams_comp_tempk (OutputField, ScrT3N01)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'temperature'
    onePostGrid%fieldUnits = 'K'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_tempk



  subroutine Brams2Post_tveg (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)

    ! OutputField <- field 'VEG_TEMP' from var_table
    ! OutputField <- rams_comp_tempc (OutputField)

    call GetVarFromMemToOutput ('VEG_TEMP', oneBramsGrid%currGrid, OutputField)
    call rams_comp_tempc (OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 7
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'vegetation temperature'
    onePostGrid%fieldUnits = 'C'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_tveg



  subroutine Brams2Post_totpcp (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)

    ! OutputField <- 0.0
    ! ScrT2N01 <- field 'ACCPR' from var_table
    ! OutputField <- OutputField + ScrT2N01
    ! ScrT2N01 <- field 'ACCPP' from var_table
    ! OutputField <- OutputField + ScrT2N01
    ! ScrT2N01 <- field 'ACCPS' from var_table
    ! OutputField <- OutputField + ScrT2N01
    ! ScrT2N01 <- field 'ACCPA' from var_table
    ! OutputField <- OutputField + ScrT2N01
    ! ScrT2N01 <- field 'ACCPG' from var_table
    ! OutputField <- OutputField + ScrT2N01
    ! ScrT2N01 <- field 'ACCPH' from var_table
    ! OutputField <- OutputField + ScrT2N01
    ! OutputField <- max(OutputField, 0.0)

    OutputField = 0.0
    call GetVarFromMemToOutput ('ACCPR', oneBramsGrid%currGrid, ScrT2N01)
    OutputField = OutputField + ScrT2N01
    call GetVarFromMemToOutput ('ACCPP', oneBramsGrid%currGrid, ScrT2N01)
    OutputField = OutputField + ScrT2N01
    call GetVarFromMemToOutput ('ACCPS', oneBramsGrid%currGrid, ScrT2N01)
    OutputField = OutputField + ScrT2N01
    call GetVarFromMemToOutput ('ACCPA', oneBramsGrid%currGrid, ScrT2N01)
    OutputField = OutputField + ScrT2N01
    call GetVarFromMemToOutput ('ACCPG', oneBramsGrid%currGrid, ScrT2N01)
    OutputField = OutputField + ScrT2N01
    call GetVarFromMemToOutput ('ACCPH', oneBramsGrid%currGrid, ScrT2N01)
    OutputField = OutputField + ScrT2N01
    OutputField = max(OutputField, 0.0)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'total resolved precip'
    onePostGrid%fieldUnits = 'mm liq'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_totpcp



  subroutine Brams2Post_acccon (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    ! OutputField <- field 'ACONPR' from var_table
    ! OutputField <- max(OutputField, 0.0)

    call GetVarFromMemToOutput ('ACONPR', oneBramsGrid%currGrid, OutputField)
    OutputField = max(OutputField, 0.0)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'accum convective pcp'
    onePostGrid%fieldUnits = 'mm'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_acccon



  subroutine Brams2Post_dewptc (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'RV' from var_table
    ! ScrT3N01 <- field 'PI' from var_table
    ! ScrT3N02 <- field 'THETA' from var_table
    ! OutputField <- rams_comp_dewk (OutputField, ScrT3N01, ScrT3N02)
    ! OutputField <- rams_comp_tempc (OutputField)

    call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, OutputField)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
    call rams_comp_dewk (OutputField, ScrT3N01, ScrT3N02)
    call rams_comp_tempc (OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'dewpoint temp'
    onePostGrid%fieldUnits = 'C'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_dewptc



  subroutine Brams2Post_rshort (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    ! OutputField <- field 'RSHORT' from var_table

    call GetVarFromMemToOutput ('RSHORT', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'rshort'
    onePostGrid%fieldUnits = 'W/m2'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_rshort



  subroutine Brams2Post_rlong (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    ! OutputField <- field 'RLONG' from var_table

    call GetVarFromMemToOutput ('RLONG', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'rlong'
    onePostGrid%fieldUnits = 'W/m2'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_rlong



  subroutine Brams2Post_sea_press (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    ! ScrT2N01 <- field 'TOPT' from var_table
    ! ScrT3N01 <- field 'PI' from var_table
    ! ScrT3N02 <- field 'THETA' from var_table
    ! ScrT3N03 <- field 'RV' from var_table
    ! ScrT3N02 <- rams_comp_thetv (ScrT3N02, ScrT3N03)
    ! OutputField <- ScrT3N03(:,:,1)
    ! OutputField <- rams_comp_slpmm5 (ScrT3N02, ScrT3N01, ScrT2N01)

    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
    call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, ScrT3N03)
    call rams_comp_thetv (ScrT3N02, ScrT3N03)
    OutputField = ScrT3N03(:,:,1)
    call rams_comp_slpmm5 (ScrT3N02, ScrT3N01, ScrT2N01, OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'sea level pressure;'
    onePostGrid%fieldUnits = 'mb;'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_sea_press



  subroutine Brams2Post_cape (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    !- rel hum (e)
    ! ScrT3N01 <- field 'RV' from var_table
    ! ScrT3N02 <- field 'PI' from var_table
    ! ScrT3N03 <- field 'THETA' from var_table
    ! ScrT3N01 <- rams_comp_rh (ScrT3N01, ScrT3N02, ScrT3N03)
    ! ScrT3N01 <- max(ScrT3N01, 0.0)
    !- tempk (d)
    ! ScrT3N03 <- field 'THETA' from var_table
    ! ScrT3N02 <- field 'PI' from var_table
    ! ScrT3N03 <- rams_comp_tempk (ScrT3N03, ScrT3N02)
    !- press (c)
    ! ScrT3N02 <- field 'PI' from var_table
    ! ScrT3N02 <- rams_comp_press (ScrT3N02)
    !- cape
    ! OutputField <- cape_cine (ScrT3N02, ScrT3N03, ScrT3N01, OutputField, 'cape')

    call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, ScrT3N01)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N03)
    call rams_comp_rh (ScrT3N01, ScrT3N02, ScrT3N03)
    ScrT3N01 = max(ScrT3N01, 0.0)
    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N03)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
    call rams_comp_tempk (ScrT3N03, ScrT3N02)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
    call rams_comp_press(ScrT3N02)
    call cape_cine (ScrT3N02, ScrT3N03, ScrT3N01, OutputField, 'cape')

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'cape'
    onePostGrid%fieldUnits = 'J/kg'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_cape



  subroutine Brams2Post_cine (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    !- rel hum (e)
    ! ScrT3N01 <- field 'RV' from var_table
    ! ScrT3N02 <- field 'PI' from var_table
    ! ScrT3N03 <- field 'THETA' from var_table
    ! ScrT3N01 <- rams_comp_rh (ScrT3N01, ScrT3N02, ScrT3N03)
    ! ScrT3N01 <- max(ScrT3N01, 0.0)
    !- tempk (d)
    ! ScrT3N03 <- field 'THETA' from var_table
    ! ScrT3N02 <- field 'PI' from var_table
    ! ScrT3N03 <- rams_comp_tempk (ScrT3N03, ScrT3N02)
    !- press (c)
    ! ScrT3N02 <- field 'PI' from var_table
    ! ScrT3N02 <- rams_comp_press (ScrT3N02)
    !- cine
    ! OutputField <- cape_cine (ScrT3N02, ScrT3N03, ScrT3N01, OutputField, 'cine')

    call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, ScrT3N01)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N03)
    call rams_comp_rh (ScrT3N01, ScrT3N02, ScrT3N03)
    ScrT3N01 = max(ScrT3N01, 0.0)
    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N03)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
    call rams_comp_tempk (ScrT3N03, ScrT3N02)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N02)
    call rams_comp_press(ScrT3N02)
    call cape_cine (ScrT3N02, ScrT3N03, ScrT3N01, OutputField, 'cine')

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'cine'
    onePostGrid%fieldUnits = 'J/kg'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_cine



  subroutine Brams2Post_topo (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    ! OutputField <- field 'TOPT' from var_table

    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'topo'
    onePostGrid%fieldUnits = 'm'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_topo



  subroutine Brams2Post_precip (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)

    ! OutputField <- 0.0
    ! ScrT2N01 <- field 'ACCPR' from var_table
    ! OutputField <- OutputField + ScrT2N01
    ! ScrT2N01 <- field 'ACCPP' from var_table
    ! OutputField <- OutputField + ScrT2N01
    ! ScrT2N01 <- field 'ACCPS' from var_table
    ! OutputField <- OutputField + ScrT2N01
    ! ScrT2N01 <- field 'ACCPA' from var_table
    ! OutputField <- OutputField + ScrT2N01
    ! ScrT2N01 <- field 'ACCPG' from var_table
    ! OutputField <- OutputField + ScrT2N01
    ! ScrT2N01 <- field 'ACCPH' from var_table
    ! OutputField <- OutputField + ScrT2N01
    ! ScrT2N01 <- field 'ACONPR' from var_table
    ! OutputField <- OutputField + ScrT2N01
    ! OutputField <- max(OutputField, 0.0)

    OutputField = 0.0
    call GetVarFromMemToOutput ('ACCPR', oneBramsGrid%currGrid, ScrT2N01)
    OutputField = OutputField + ScrT2N01
    call GetVarFromMemToOutput ('ACCPP', oneBramsGrid%currGrid, ScrT2N01)
    OutputField = OutputField + ScrT2N01
    call GetVarFromMemToOutput ('ACCPS', oneBramsGrid%currGrid, ScrT2N01)
    OutputField = OutputField + ScrT2N01
    call GetVarFromMemToOutput ('ACCPA', oneBramsGrid%currGrid, ScrT2N01)
    OutputField = OutputField + ScrT2N01
    call GetVarFromMemToOutput ('ACCPG', oneBramsGrid%currGrid, ScrT2N01)
    OutputField = OutputField + ScrT2N01
    call GetVarFromMemToOutput ('ACCPH', oneBramsGrid%currGrid, ScrT2N01)
    OutputField = OutputField + ScrT2N01
    call GetVarFromMemToOutput ('ACONPR', oneBramsGrid%currGrid, ScrT2N01)
    OutputField = OutputField + ScrT2N01
    OutputField = max(OutputField, 0.0)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'total accum precip'
    onePostGrid%fieldUnits = 'mm liq'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_precip



  subroutine Brams2Post_le (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'SFLUX_R' from var_table
    ! ScrT2N01 <- field 'TOPT' from var_table
    ! ScrT3N01, ScrT3N02, ScrT3N03,  <- rams_comp_dn0 (ScrT2N01, pi01dn, th01dn, ztn, ztop, dzmn)
    ! ScrT2N01 <- ScrT3N03(:,:,1)
    ! OutputField <- OutputField * ScrT2N01
    ! OutputField <- OutputField * 2.5e6

    call GetVarFromMemToOutput ('SFLUX_R', oneBramsGrid%currGrid, OutputField)
    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
    call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%pi01dn, &
         oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)
    call rams_get_surface(ScrT2N01, ScrT3N03)
    OutputField = OutputField * ScrT2N01
    OutputField = OutputField * 2.5e6

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'sfc lat heat flx'
    onePostGrid%fieldUnits = 'W/m2'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_le



  subroutine Brams2Post_h (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'SFLUX_T' from var_table
    ! ScrT2N01 <- field 'TOPT' from var_table
    ! ScrT3N01, ScrT3N02, ScrT3N03,  <- rams_comp_dn0 (ScrT2N01, pi01dn, th01dn, ztn, ztop, dzmn)
    ! ScrT2N01 <- ScrT3N03(:,:,1)
    ! OutputField <- OutputField * ScrT2N01
    ! OutputField <- OutputField * 1004.

    call GetVarFromMemToOutput ('SFLUX_T', oneBramsGrid%currGrid, OutputField)
    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
    call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%pi01dn, &
         oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)
    call rams_get_surface(ScrT2N01, ScrT3N03)
    OutputField = OutputField * ScrT2N01
    OutputField = OutputField * 1004.

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'sfc sens heat flx'
    onePostGrid%fieldUnits = 'W/m2'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_h



  subroutine Brams2Post_rv (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'RV' from var_table
    ! OutputField <- OutputField * 1.e3
    ! OutputField <- max(OutputField, 0.0)

    call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, OutputField)
    OutputField = OutputField * 1.e3
    OutputField = max(OutputField, 0.0)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'vapor mix ratio'
    onePostGrid%fieldUnits = 'g/kg'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_rv



  subroutine Brams2Post_rlongup (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    ! OutputField <- field 'RLONGUP' from var_table

    call GetVarFromMemToOutput ('RLONGUP', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'rlongup'
    onePostGrid%fieldUnits = 'W/m2'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_rlongup



  subroutine Brams2Post_albedt (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    ! OutputField <- field 'ALBEDT' from var_table

    call GetVarFromMemToOutput ('ALBEDT', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'albedt'
    onePostGrid%fieldUnits = ' '

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_albedt



  subroutine Brams2Post_geo (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! ScrT2N01 <- field 'TOPT' from var_table
    ! OutputField <- rams_comp_z (OutputField, ScrT2N01, oneBramsGrid%ztn, oneBramsGrid%ztop)

    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
    call rams_comp_z (OutputField, ScrT2N01, oneBramsGrid%ztn, oneBramsGrid%ztop)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'geopotential height'
    onePostGrid%fieldUnits = 'm'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_geo



  subroutine Brams2Post_ue_avg (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    integer :: firstX, lastX, firstY, lastY
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'UP' from var_table
    ! ScrT3N01 <- field 'VP' from var_table
    ! OutputField, ScrT3N01,  <- rams_comp_rotate (OutputField, ScrT3N01, xtn, ytn, polelat, polelon)
    ! OutputField <- rams_comp_avgu (OutputField)

    firstX = oneBramsGrid%nodei0(oneBramsGrid%mynum)+1
    lastX = oneBramsGrid%nodei0(oneBramsGrid%mynum)+oneBramsGrid%mxp
    firstY = oneBramsGrid%nodej0(oneBramsGrid%mynum)+1
    lastY = oneBramsGrid%nodej0(oneBramsGrid%mynum)+oneBramsGrid%myp
    call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, OutputField)
    call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N01)
    call rams_comp_rotate (OutputField, ScrT3N01, &
         oneBramsGrid%xtn(firstX:lastX), oneBramsGrid%ytn(firstY:lastY), &
         oneBramsGrid%polelat, oneBramsGrid%polelon)
    call rams_comp_avgu (OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'ue_avg'
    onePostGrid%fieldUnits = 'm/s'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_ue_avg



  subroutine Brams2Post_ve_avg (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    integer :: firstX, lastX, firstY, lastY
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'VP' from var_table
    ! ScrT3N01 <- field 'UP' from var_table
    ! ScrT3N01, OutputField,  <- rams_comp_rotate (ScrT3N01, OutputField, xtn, ytn, polelat, polelon)
    ! OutputField <- rams_comp_avgv (OutputField)

    firstX = oneBramsGrid%nodei0(oneBramsGrid%mynum)+1
    lastX = oneBramsGrid%nodei0(oneBramsGrid%mynum)+oneBramsGrid%mxp
    firstY = oneBramsGrid%nodej0(oneBramsGrid%mynum)+1
    lastY = oneBramsGrid%nodej0(oneBramsGrid%mynum)+oneBramsGrid%myp
    call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, OutputField)
    call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
    call rams_comp_rotate (ScrT3N01, OutputField, &
         oneBramsGrid%xtn(firstX:lastX), oneBramsGrid%ytn(firstY:lastY), &
         oneBramsGrid%polelat, oneBramsGrid%polelon)
    call rams_comp_avgv (OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 've_avg'
    onePostGrid%fieldUnits = 'm/s'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_ve_avg



  subroutine Brams2Post_tempc2m (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT6N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N04(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N05(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT1N01
    real :: ScrT1N02

    ! ScrT3N01 <- field 'UP' from var_table
    ! ScrT3N02 <- field 'VP' from var_table
    ! ScrT3N01 <- rams_comp_speed (ScrT3N01, ScrT3N02)
    ! ScrT3N02 <- field 'THETA' from var_table
    ! ScrT3N03 <- field 'PI' from var_table
    ! ScrT2N01 <- field 'TOPT' from var_table
    ! ScrT6N01 <- field 'USTAR' from var_table
    ! ScrT6N02 <- field 'PATCH_ROUGH' from var_table
    ! ScrT6N03 <- field 'CAN_TEMP' from var_table
    ! ScrT6N04 <- field 'PATCH_AREA' from var_table
    ! ScrT6N05 <- field 'TSTAR' from var_table
    ! OutputField <- rams_reduced_temp (ScrT3N01, ScrT6N01, ScrT6N05, 2., ztn, ScrT6N02, ScrT6N04, ScrT6N03, ScrT3N02, ScrT3N03, ScrT2N01, ztop)
    ! ScrT2N01 <- ScrT3N03(:,:,1)
    ! OutputField <- rams_comp_tempk (OutputField, ScrT2N01)
    ! OutputField <- rams_comp_tempc (OutputField)

    call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
    call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
    call rams_comp_speed (ScrT3N01, ScrT3N02)
    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N03)
    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
    call GetVarFromMemToOutput ('USTAR', oneBramsGrid%currGrid, ScrT6N01)
    call GetVarFromMemToOutput ('PATCH_ROUGH', oneBramsGrid%currGrid, ScrT6N02)
    call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N03)
    call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N04)
    call GetVarFromMemToOutput ('TSTAR', oneBramsGrid%currGrid, ScrT6N05)
    call rams_reduced_temp (OutputField, ScrT3N01, ScrT6N01, ScrT6N05, 2., oneBramsGrid%ztn(2), &
         ScrT6N02, ScrT6N04, ScrT6N03, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%ztop)
    call rams_get_surface(ScrT2N01, ScrT3N03)
    call rams_comp_tempk (OutputField, ScrT2N01)
    call rams_comp_tempc (OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'temp - 2m AGL;'
    onePostGrid%fieldUnits = 'C'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_tempc2m



  subroutine Brams2Post_tempc (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer     :: oneBramsGrid
    type(PostGrid), pointer      :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'THETA' from var_table
    ! ScrT3N01 <- field 'PI' from var_table
    ! OutputField <- rams_comp_tempk (OutputField, ScrT3N01)
    ! OutputField <- rams_comp_tempc (OutputField)

    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, OutputField)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
    call rams_comp_tempk (OutputField, ScrT3N01)
    call rams_comp_tempc (OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'temperature'
    onePostGrid%fieldUnits = 'C'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_tempc



  subroutine Brams2Post_rh (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'RV' from var_table
    ! ScrT3N01 <- field 'PI' from var_table
    ! ScrT3N02 <- field 'THETA' from var_table
    ! OutputField <- rams_comp_rh (OutputField, ScrT3N01, ScrT3N02)
    ! OutputField <- max(OutputField, 0.0)

    call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, OutputField)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
    call rams_comp_rh (OutputField, ScrT3N01, ScrT3N02)
    OutputField = max(OutputField, 0.0)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'relative humidity'
    onePostGrid%fieldUnits = 'pct'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_rh



  subroutine Brams2Post_w (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'WP' from var_table

    call GetVarFromMemToOutput ('WP', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'w'
    onePostGrid%fieldUnits = 'm/s'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_w
  
  

  subroutine Brams2Post_omeg (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'WP' from var_table
    ! OutputField <- field 'TOPT' from var_table

    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
    call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%pi01dn, &
         oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)
    call GetVarFromMemToOutput ('WP', oneBramsGrid%currGrid, ScrT3N01)
    call calc_omeg (OutputField,ScrT3N01,ScrT3N03)
    
    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'omega'
    onePostGrid%fieldUnits = 'Pa/s'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_omeg



  subroutine Brams2Post_sst (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: ScrT4N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%nzg, oneBramsGrid%npatch)
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    ! ScrT4N01 <- field 'SOIL_ENERGY' from var_table
    !kp        = nzg
    ! OutputField <- rams_fill_sst (oneBramsGrid%nzg, OutputField, ScrT4N01)

    call GetVarFromMemToOutput ('SOIL_ENERGY', oneBramsGrid%currGrid, ScrT4N01)
    call rams_fill_sst (oneBramsGrid%nzg, OutputField, ScrT4N01)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'water temperature'
    onePostGrid%fieldUnits = 'C'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_sst



  subroutine Brams2Post_land (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: ScrT6N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    ! ScrT6N01 <- field 'PATCH_AREA' from var_table
    ! OutputField <- ScrT6N01(:,:,1)
    ! OutputField <- 1.0 - OutputField

    call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N01)
    call rams_get_surface(OutputField, ScrT6N01)
    !DSM OutputField = 1.0 - OutputField

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'land frac area'
    onePostGrid%fieldUnits = ''

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_land



  subroutine Brams2Post_smoist (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: ScrT4N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%nzg, oneBramsGrid%npatch)
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%nzg, oneBramsGrid%npatch)

    ! ScrT4N01 <- field 'SOIL_WATER' from var_table
    ! OutputField <- get_leaf_soil (ScrT4N01, OutputField)

    call GetVarFromMemToOutput ('SOIL_WATER', oneBramsGrid%currGrid, ScrT4N01)
    call get_leaf_soil (ScrT4N01, OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 8
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'soil moisture'
    onePostGrid%fieldUnits = 'm3/m3'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_smoist



  subroutine Brams2Post_zi (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)

    ! OutputField <- field 'TKEP' from var_table
    ! OutputField <- max(OutputField, 0.0)
    ! ScrT2N01 <- field 'TOPT' from var_table
    ! OutputField <- rams_comp_pbl (OutputField, ScrT2N01, oneBramsGrid%ztn, oneBramsGrid%ztop)

    call GetVarFromMemToOutput ('TKEP', oneBramsGrid%currGrid, OutputField)
    OutputField = max(OutputField, 0.0)
    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
    call rams_comp_pbl (OutputField, ScrT2N01, oneBramsGrid%ztn, oneBramsGrid%ztop)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'Zi'
    onePostGrid%fieldUnits = 'm'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_zi






  subroutine Brams2Post_tke (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'TKEP' from var_table
    ! OutputField <- max(OutputField, 0.0)

    call GetVarFromMemToOutput ('TKEP', oneBramsGrid%currGrid, OutputField)
    OutputField = max(OutputField, 0.0)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'turb kinetic energy'
    onePostGrid%fieldUnits = 'm2/s2'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_tke




  subroutine Brams2Post_cloud (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'RCP' from var_table
    ! OutputField <- OutputField * 1.e3
    ! OutputField <- max(OutputField, 0.0)

    call GetVarFromMemToOutput ('RCP', oneBramsGrid%currGrid, OutputField)
    OutputField = OutputField * 1.e3
    OutputField = max(OutputField, 0.0)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'cloud mix ratio'
    onePostGrid%fieldUnits = 'g/kg'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_cloud
  
  
  
  subroutine Brams2Post_pwt (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid
    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N04(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N05(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
    
    ScrT3N05 = 0.0
    call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid, ScrT3N04)
    ScrT3N05 = ScrT3N05 + ScrT3N04
    call GetVarFromMemToOutput ('RCP', oneBramsGrid%currGrid, ScrT3N04)
    ScrT3N05 = ScrT3N05 + ScrT3N04
    call GetVarFromMemToOutput ('RRP', oneBramsGrid%currGrid, ScrT3N04)
    ScrT3N05 = ScrT3N05 + ScrT3N04
    call GetVarFromMemToOutput ('RPP', oneBramsGrid%currGrid, ScrT3N04)
    ScrT3N05 = ScrT3N05 + ScrT3N04
    call GetVarFromMemToOutput ('RSP', oneBramsGrid%currGrid, ScrT3N04)
    ScrT3N05 = ScrT3N05 + ScrT3N04
    call GetVarFromMemToOutput ('RAP', oneBramsGrid%currGrid, ScrT3N04)
    ScrT3N05 = ScrT3N05 + ScrT3N04
    call GetVarFromMemToOutput ('RGP', oneBramsGrid%currGrid, ScrT3N04)
    ScrT3N05 = ScrT3N05 + ScrT3N04
    call GetVarFromMemToOutput ('RHP', oneBramsGrid%currGrid, ScrT3N04)
    ScrT3N05 = ScrT3N05 + ScrT3N04
    
    ScrT3N05 = max(ScrT3N05, 0.0)
    
    
    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
    call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, ScrT2N01, oneBramsGrid%pi01dn, &
         oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)
    ScrT3N05=ScrT3N05*ScrT3N03
    call comp_vertint (ScrT3N05,ScrT3N05,ScrT2N01,oneBramsGrid%ztop,oneBramsGrid%zmn)
    
    call rams_get_surface(OutputField, ScrT3N05)

    OutputField=OutputField*0.1
    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'precipitable total water'
    onePostGrid%fieldUnits = 'cm'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_pwt
  
  
  
  
  subroutine Brams2Post_slp_metar(varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid
    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid,ScrT3N01)
    call RAMS_comp_press(ScrT3N01)
    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
    call rams_comp_tempk(ScrT3N02,ScrT3N01)
    call comp_slp_metar(OutputField,ScrT3N01,ScrT2N01,ScrT3N02,oneBramsGrid%ztn)
    

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'sea level pressure -  METAR formulation'
    onePostGrid%fieldUnits = 'mb'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_slp_metar
  
  
  
  
subroutine Brams2Post_td2m(varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid
    
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT2N02(oneBramsGrid%mxp, oneBramsGrid%myp)    
    real :: ScrT2N03(oneBramsGrid%mxp, oneBramsGrid%myp)    
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N04(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT6N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N04(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N05(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N06(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)

    
    call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
    call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
    call rams_comp_speed (ScrT3N01,ScrT3N02)
    call GetVarFromMemToOutput ('RV', oneBramsGrid%currGrid,ScrT3N02 )
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid,ScrT3N03)
    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N04)
    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
    call GetVarFromMemToOutput ('USTAR', oneBramsGrid%currGrid, ScrT6N01)
    call GetVarFromMemToOutput ('PATCH_ROUGH', oneBramsGrid%currGrid, ScrT6N02)
    call GetVarFromMemToOutput ('CAN_RVAP', oneBramsGrid%currGrid, ScrT6N03)
    call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N04)
    call GetVarFromMemToOutput ('RSTAR', oneBramsGrid%currGrid, ScrT6N05)
    call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N06)
    
    call rams_reduced_rv (OutputField, ScrT3N01,ScrT6N01,ScrT6N05,2.,oneBramsGrid%ztn(2), &
                           ScrT6N02,ScrT6N04,ScrT6N03,ScrT3N02,ScrT3N03,ScrT2N01, &
                           oneBramsGrid%ztop,ScrT6N06,ScrT3N04)
    
    OutputField=OutputField*1.e3
    OutputField = max(OutputField, 0.0)
    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
    call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N03)
    call GetVarFromMemToOutput ('TSTAR', oneBramsGrid%currGrid, ScrT6N05)
    
    call rams_reduced_temp (ScrT2N02, ScrT3N01, ScrT6N01, ScrT6N05, 2., oneBramsGrid%ztn(2), &
                             ScrT6N02, ScrT6N04, ScrT6N03, ScrT3N02, ScrT3N03, ScrT2N01, &
                             oneBramsGrid%ztop)
                             
    ScrT2N03=(ScrT3N03(:,:,1)+ScrT3N03(:,:,2))*0.5
    call rams_comp_tempk (ScrT2N02, ScrT2N03)
    call RAMS_comp_dewK_2m (OutputField, ScrT2N03, ScrT2N02)
    call RAMS_comp_tempC (OutputField)
    
    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'dewpoint temp in 2m'
    onePostGrid%fieldUnits = 'C'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_td2m
  
  
  
  
  subroutine Brams2Post_u10m(varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid
    
    integer :: firstX, lastX, firstY, lastY
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT6N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N04(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    
    call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
    call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
    call rams_comp_speed (ScrT3N01,ScrT3N02)
    
    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid,ScrT3N03)
    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
    call GetVarFromMemToOutput ('USTAR', oneBramsGrid%currGrid, ScrT6N01)
    call GetVarFromMemToOutput ('PATCH_ROUGH', oneBramsGrid%currGrid, ScrT6N02)
    call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N03)
    call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N04)
    
    call rams_reduced_wind (OutputField,ScrT3N01,ScrT6N01,10.,oneBramsGrid%ztn(2),ScrT6N02,ScrT6N04, &
                             ScrT6N03,ScrT3N02,ScrT3N03,ScrT2N01,oneBramsGrid%ztop)
    call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
    call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
    
    firstX = oneBramsGrid%nodei0(oneBramsGrid%mynum)+1
    lastX = oneBramsGrid%nodei0(oneBramsGrid%mynum)+oneBramsGrid%mxp
    firstY = oneBramsGrid%nodej0(oneBramsGrid%mynum)+1
    lastY = oneBramsGrid%nodej0(oneBramsGrid%mynum)+oneBramsGrid%myp
    call RAMS_comp_dir (ScrT3N01, ScrT3N02, oneBramsGrid%xtn(firstX:lastX), oneBramsGrid%ytn(firstY:lastY), &
                        oneBramsGrid%polelat, oneBramsGrid%polelon)
    call calc_u10m (OutputField,ScrT3N01)
    
    !--- Indefinindo campos "zerados" - ex: analise ---
    if (maxval(OutputField)==0.0 .and. minval(OutputField)==0.0) OutputField=-9.99e+33
    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'Zonal Wind - 10m'
    onePostGrid%fieldUnits = 'm/s'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_u10m
  
  
  
  
  subroutine Brams2Post_v10m(varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid
    
    integer :: firstX, lastX, firstY, lastY
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT6N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    real :: ScrT6N04(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)
    
    call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
    call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
    call rams_comp_speed (ScrT3N01,ScrT3N02)
    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, ScrT3N02)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid,ScrT3N03)
    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, ScrT2N01)
    call GetVarFromMemToOutput ('USTAR', oneBramsGrid%currGrid, ScrT6N01)
    call GetVarFromMemToOutput ('PATCH_ROUGH', oneBramsGrid%currGrid, ScrT6N02)
    call GetVarFromMemToOutput ('CAN_TEMP', oneBramsGrid%currGrid, ScrT6N03)
    call GetVarFromMemToOutput ('PATCH_AREA', oneBramsGrid%currGrid, ScrT6N04)
    
    call rams_reduced_wind (OutputField,ScrT3N01,ScrT6N01,10.,oneBramsGrid%ztn(2),ScrT6N02,ScrT6N04, &
                             ScrT6N03,ScrT3N02,ScrT3N03,ScrT2N01,oneBramsGrid%ztop)
    call GetVarFromMemToOutput ('UP', oneBramsGrid%currGrid, ScrT3N01)
    call GetVarFromMemToOutput ('VP', oneBramsGrid%currGrid, ScrT3N02)
    
    firstX = oneBramsGrid%nodei0(oneBramsGrid%mynum)+1
    lastX = oneBramsGrid%nodei0(oneBramsGrid%mynum)+oneBramsGrid%mxp
    firstY = oneBramsGrid%nodej0(oneBramsGrid%mynum)+1
    lastY = oneBramsGrid%nodej0(oneBramsGrid%mynum)+oneBramsGrid%myp
    call RAMS_comp_dir (ScrT3N01, ScrT3N02, oneBramsGrid%xtn(firstX:lastX), oneBramsGrid%ytn(firstY:lastY), &
                        oneBramsGrid%polelat, oneBramsGrid%polelon)
    call calc_v10m (OutputField,ScrT3N01)
    
    !--- Indefinindo campos "zerados" - ex: analise ---
    if (maxval(OutputField)==0.0 .and. minval(OutputField)==0.0) OutputField=-9.99e+33
    
    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'Meridional Wind - 10m'
    onePostGrid%fieldUnits = 'm/s'

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_v10m
  
  subroutine Brams2Post_vtype (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)


    call GetVarFromMemToOutput ('LEAF_CLASS', oneBramsGrid%currGrid, OutputField)
    call rams_comp_vegclass(OutputField)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 7
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'vegetation class'
    onePostGrid%fieldUnits = '#'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_vtype
  
 
  
  subroutine Brams2Post_sltex_p (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%nzg, oneBramsGrid%npatch)

    call GetVarFromMemToOutput ('SOIL_TEXT', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 8
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'soil textural class'
    onePostGrid%fieldUnits = '#'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_sltex_p
    
 
  
  
  subroutine Brams2Post_ndvi (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)

    call GetVarFromMemToOutput ('VEG_NDVIC', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 7
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'ndvi'
    onePostGrid%fieldUnits = '#'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_ndvi
    
     
 
  
  subroutine Brams2Post_lai (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%npatch)

    call GetVarFromMemToOutput ('VEG_LAI', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 7
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'green leaf area index'
    onePostGrid%fieldUnits = ''

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_lai


    
  subroutine Brams2Post_u10mj (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    call GetVarFromMemToOutput ('U10MJ', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'Zonal Wind at 10m - from JULES'
    onePostGrid%fieldUnits = 'm/s'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_u10mj


    
  subroutine Brams2Post_v10mj (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    call GetVarFromMemToOutput ('V10MJ', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'Meridional Wind at 10m - from JULES'
    onePostGrid%fieldUnits = 'm/s'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_v10mj

    
  subroutine Brams2Post_t2mj (varName, oneBramsGrid, onePostGrid)
    real, parameter :: undef=-9.99e+33
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    call GetVarFromMemToOutput ('T2MJ', oneBramsGrid%currGrid, OutputField)
    OutputField=OutputField - 273.16
    where (OutputField<-70.) OutputField=undef


    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'Temperature at 2m - from JULES'
    onePostGrid%fieldUnits = 'C'
    
    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_t2mj


  subroutine Brams2Post_csj (varName, oneBramsGrid, onePostGrid)
    real, parameter :: undef=-9.99e+33
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    call GetVarFromMemToOutput ('CSJ', oneBramsGrid%currGrid, OutputField)


    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'Soil Carbon - from JULES'
    onePostGrid%fieldUnits = 'kg/m2'
    
    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_csj

!---RB
  subroutine Brams2Post_theta (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    call GetVarFromMemToOutput ('THETA', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'Potential temp'
    onePostGrid%fieldUnits = 'K'
    
    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_theta

  subroutine Brams2Post_rv2mj (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    call GetVarFromMemToOutput ('RV2MJ', oneBramsGrid%currGrid, OutputField)
    OutputField=OutputField*1000.
    
    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'Mixing rate at 2m - from JULES'
    onePostGrid%fieldUnits = 'g/kg'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_rv2mj

!----
  subroutine Brams2Post_td2mj (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT2N01(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT2N03(oneBramsGrid%mxp, oneBramsGrid%myp)    
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    call GetVarFromMemToOutput ('RV2MJ', oneBramsGrid%currGrid, OutputField)
    OutputField=OutputField*1.e3
    call GetVarFromMemToOutput ('T2MJ', oneBramsGrid%currGrid, ScrT2N01)
    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)

    ScrT2N03=(ScrT3N01(:,:,1)+ScrT3N01(:,:,2))*0.5
    call RAMS_comp_dewK_2m(OutputField,ScrT2N03,ScrT2N01)

    OutputField=OutputField - 273.16
    
    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'Dewpoint temp at 2m - from JULES'
    onePostGrid%fieldUnits = 'C'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_td2mj


!----
  subroutine Brams2Post_co2 (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    call GetVarFromMemToOutput ('CO2P', oneBramsGrid%currGrid, OutputField)

    OutputField = OutputField*(28.96/44.)*1.E-3  ! TRANSFORMACAO DE kg[CO2]/kg[AR] para ppm (PARTE POR MILHAO)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'CO2 Concentration'
    onePostGrid%fieldUnits = 'ppmv'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_co2



!----RB
  subroutine Brams2Post_co2_total (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: OutputField2(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: OutputField3(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N05(oneBramsGrid%mxp, oneBramsGrid%myp)

   
    call GetVarFromMemToOutput ('CO2P', oneBramsGrid%currGrid, OutputField3)
    ! call RAMS_comp_noneg(n1,n2,n3,a)
    OutputField3 = max(OutputField3, 0.0)
    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, OutputField2)
    call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, OutputField2, oneBramsGrid%pi01dn, &
         oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)
    OutputField3 = OutputField3 * ScrT3N03

    call GetVarFromMemToOutput ('PI', oneBramsGrid%currGrid, ScrT3N01)
    call rams_comp_press(ScrT3N01)
!Ate aqui acho que esta OK
    
    call comp_vertint_press (ScrT3N05,OutputField3,ScrT3N01,OutputField2,oneBramsGrid%ztop,oneBramsGrid%zmn)

    OutputField = ScrT3N05 * 1e-9*1000.*1./44.*1e-4*6.022e23

    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'CO2_total column'
    onePostGrid%fieldUnits = 'molec/cm^2'
    
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_co2_total

!----
  subroutine Brams2Post_co2_antro (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: OutputField2(oneBramsGrid%mxp, oneBramsGrid%myp)

    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, OutputField2)
    call GetVarFromMemToOutput ('CO2_antro_SRC', oneBramsGrid%currGrid, OutputField)

    OutputField(:,:,2) = OutputField(:,:,1)*(1.-OutputField2(:,:)/oneBramsGrid%ztop) &
                          /oneBramsGrid%dztn(2)/86400.*1.e-3  ! convertendo de kg/m3/dia para mg/m2/s

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'CO2_antro_SRC flux'
    onePostGrid%fieldUnits = 'mg/m2/s'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_co2_antro

!----
  subroutine Brams2Post_co2_bioge (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: OutputField2(oneBramsGrid%mxp, oneBramsGrid%myp)
  
    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, OutputField2)
    call GetVarFromMemToOutput ('CO2_bioge_SRC', oneBramsGrid%currGrid, OutputField)

    OutputField(:,:,2) = OutputField(:,:,1)*(1.-OutputField2(:,:)/oneBramsGrid%ztop) &
                          /oneBramsGrid%dztn(2)*1.e-3  ! convertendo de kg/m3/dia para mg/m2/s

    OutputField = OutputField

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'CO2_bioge_SRC flux'
    onePostGrid%fieldUnits = 'mg/m2/s'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_co2_bioge

!----
  subroutine Brams2Post_co2_burn2D (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: OutputField2(oneBramsGrid%mxp, oneBramsGrid%myp)

    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, OutputField2)
    call GetVarFromMemToOutput ('CO2_bburn_SRC', oneBramsGrid%currGrid, OutputField)

    OutputField(:,:,2) = OutputField(:,:,1)*(1.-OutputField2(:,:)/oneBramsGrid%ztop) &
                          /oneBramsGrid%dztn(2)/86400.*1.e-3  ! convertendo de kg/m3/dia para mg/m2/s

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'CO2_bburn_SRC flux in first level'
    onePostGrid%fieldUnits = 'mg/m2/s'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_co2_burn2D

!----
  subroutine Brams2Post_co2_burn (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: OutputField2(oneBramsGrid%mxp, oneBramsGrid%myp)
    integer :: k

    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, OutputField2)
    call GetVarFromMemToOutput ('CO2_bburn_SRC', oneBramsGrid%currGrid, OutputField)

    do k=1,oneBramsGrid%mzp
       OutputField(:,:,k) = OutputField(:,:,k)*(1.-OutputField2(:,:)/oneBramsGrid%ztop) &
                            /oneBramsGrid%dztn(2)/86400.*1.e-3  ! convertendo de kg/m3/dia para mg/m2/s
    enddo

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'CO2_bburn_SRC flux'
    onePostGrid%fieldUnits = 'mg/m2/s'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_co2_burn


!--(DMK-CCATT-INI)-------------------------------------------------------
  subroutine Brams2Post_CO (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'COP' from var_table

    ! ierr= RAMS_getvar('COP',idim_type,ngrd,a,b,flnm)
    call GetVarFromMemToOutput ('COP', oneBramsGrid%currGrid, OutputField)
    
    ! call RAMS_comp_noneg(n1,n2,n3,a)
    OutputField = max(OutputField, 0.0)
    
    ! call RAMS_transf_ppb(n1,n2,n3,a,28.)
    OutputField = OutputField*(PMAR/28.)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'COP Concentration'
    onePostGrid%fieldUnits = 'ppbv'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_CO




  subroutine Brams2Post_NO (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'NOP' from var_table

    ! ierr= RAMS_getvar('NO2P',idim_type,ngrd,a,b,flnm)
    call GetVarFromMemToOutput ('NOP', oneBramsGrid%currGrid, OutputField)

    ! call RAMS_comp_noneg(n1,n2,n3,a)
    OutputField = max(OutputField, 0.0)

    ! call RAMS_transf_ppb(n1,n2,n3,a,30.)
    OutputField = OutputField*(PMAR/30.)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'NOP Concentration'
    onePostGrid%fieldUnits = 'ppbv'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_NO




  subroutine Brams2Post_HNO3 (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'HNO3P' from var_table

    ! ierr= RAMS_getvar('HNO3P',idim_type,ngrd,a,b,flnm)
    call GetVarFromMemToOutput ('HNO3P', oneBramsGrid%currGrid, OutputField)

    ! call RAMS_comp_noneg(n1,n2,n3,a)
    OutputField = max(OutputField, 0.0)

    ! call RAMS_transf_ppb(n1,n2,n3,a,63.)
    OutputField = OutputField*(PMAR/63.)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'HNO3 Concentration'
    onePostGrid%fieldUnits = 'ppbv'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_HNO3




  subroutine Brams2Post_O3 (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'O3P' from var_table

    ! ierr= RAMS_getvar('O3P',idim_type,ngrd,a,b,flnm)
    call GetVarFromMemToOutput ('O3P', oneBramsGrid%currGrid, OutputField)
    
    ! call RAMS_comp_noneg(n1,n2,n3,a)
    OutputField = max(OutputField, 0.0)
    
    ! call RAMS_transf_ppb(n1,n2,n3,a,48.)
    OutputField = OutputField*(PMAR/48.)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'O3P Concentration'
    onePostGrid%fieldUnits = 'ppbv'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_O3




  subroutine Brams2Post_NO2 (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'NO2P' from var_table

    ! ierr= RAMS_getvar('NO2P',idim_type,ngrd,a,b,flnm)
    call GetVarFromMemToOutput ('NO2P', oneBramsGrid%currGrid, OutputField)
    
    ! call RAMS_comp_noneg(n1,n2,n3,a)
    OutputField = max(OutputField, 0.0)
    
    ! call RAMS_transf_ppb(n1,n2,n3,a,46.)
    OutputField = OutputField*(PMAR/46.)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'NO2 mixing ratio'
    onePostGrid%fieldUnits = 'ppbv'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_NO2




  subroutine Brams2Post_PM25 (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    real :: OutputField2(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'bburn2P' from var_table

    ! ierr= RAMS_getvar('bburn2P',idim_type,ngrd,a,b,flnm)
    call GetVarFromMemToOutput ('bburn2P', oneBramsGrid%currGrid, OutputField)

    ! call RAMS_comp_noneg(n1,n2,n3,a)
    OutputField = max(OutputField, 0.0)

    ! call RAMS_comp_mults(n1,n2,n3,a,1.e-9)  ! converte de 1e-9 kg/kg para 1 kg/kg   
    OutputField = OutputField*1.e-9

    ! ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, OutputField2)

    ! call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
    call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, OutputField2, oneBramsGrid%pi01dn, &
         oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)

    ! call RAMS_transf_ugm3(n1,n2,n3,a,d)
    OutputField = OutputField*scrT3N03*1.e+9

    ! call RAMS_comp_noneg(n1,n2,n3,a)
    OutputField = max(OutputField, 0.0)

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'PM25 Concentration'
    onePostGrid%fieldUnits = 'ug/m3'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_PM25




  subroutine Brams2Post_PM25wd (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)


    ! OutputField <- field 'bburn2WD' from var_table

    ! ierr= RAMS_getvar('bburn2WD',idim_type,ngrd,a,b,flnm)
    call GetVarFromMemToOutput ('bburn2WD', oneBramsGrid%currGrid, OutputField)
    
    ! call RAMS_comp_mults(n1,n2,n3,a,1.e-9)  ! converte de 1e-9 kg/kg para 1 kg/kg
    OutputField = OutputField*1.e-9

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'WET deposition mass tracer PM25'
    onePostGrid%fieldUnits = 'kg/m2'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_PM25wd




  subroutine Brams2Post_AOT550 (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: AOT(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp)

    ! OutputField <- field 'AOT' from var_table

    ! ierr= RAMS_getvar('AOT',idim_type,ngrd,aot,b,flnm)
    call GetVarFromMemToOutput ('AOT', oneBramsGrid%currGrid, AOT)

    ! call D3toD2(n1,n2,nwave,12,a,aot)
    Outputfield(:,:) = AOT(:,:,12)
  
    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'AOT 550nm'
    onePostGrid%fieldUnits = ' '

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_AOT550




  subroutine Brams2Post_PMINT (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    real :: OutputField2(oneBramsGrid%mxp, oneBramsGrid%myp)
    real :: ScrT3N01(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N02(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)
    real :: ScrT3N03(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'bburn2P' from var_table

    ! ierr= RAMS_getvar('bburn2P',idim_type,ngrd,a,b,flnm)
    call GetVarFromMemToOutput ('bburn2P', oneBramsGrid%currGrid, OutputField)

    ! call RAMS_comp_noneg(n1,n2,n3,a)
    OutputField = max(OutputField, 0.0)

    ! ierr= RAMS_getvar('TOPT',idim_type,ngrd,e,b,flnm)
    call GetVarFromMemToOutput ('TOPT', oneBramsGrid%currGrid, OutputField2)

    ! call RAMS_comp_dn0(n1,n2,n3,b,c,d,e,ngrd)
    call rams_comp_dn0 (ScrT3N01, ScrT3N02, ScrT3N03, OutputField2, oneBramsGrid%pi01dn, &
         oneBramsGrid%th01dn, oneBramsGrid%ztn, oneBramsGrid%ztop, oneBramsGrid%dzmn)

    ! call RAMS_comp_mult(n1,n2,n3,a,d) !Unit: kg[pm25]/m3
    OutputField = OutputField * ScrT3N03

    !call RAMS_comp_vertint(n1,n2,n3,a,e,ngrd) ! Unit: kg[pm25]/m2
    call comp_vertint (OutputField,OutputField,OutputField2,oneBramsGrid%ztop,oneBramsGrid%zmn)
    
    !call RAMS_comp_mults(n1,n2,n3,a,1.e+6)  ! converte de kg/m2 para mg/m2
    OutputField = OutputField * 1.e+6

    ! register Grads fields on PostGrid
    onePostGrid%ivar_type = 2
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'PM25 vert int'
    onePostGrid%fieldUnits = 'mg/m2'

    ! output binary field
    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_PMINT
!--(DMK-CCATT-FIM)-------------------------------------------------------

  subroutine Brams2Post_cuthdp (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'WP' from var_table

    call GetVarFromMemToOutput ('THSRC', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'deep conv heat rate'
    onePostGrid%fieldUnits = 'K/day'
    OutputField = OutputField*86400.0

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_cuthdp
  
  
  subroutine Brams2Post_curtdp (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'WP' from var_table

    call GetVarFromMemToOutput ('RTSRC', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'deep conv moist rate'
    onePostGrid%fieldUnits = 'g/kg/day'
    OutputField = OutputField*86400.0*1000.0

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_curtdp
  
    subroutine Brams2Post_cucldp (varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    real :: OutputField(oneBramsGrid%mxp, oneBramsGrid%myp, oneBramsGrid%mzp)

    ! OutputField <- field 'WP' from var_table

    call GetVarFromMemToOutput ('CLSRC', oneBramsGrid%currGrid, OutputField)

    ! register Grads fields on PostGrid

    onePostGrid%ivar_type = 3
    onePostGrid%fieldName = varName
    onePostGrid%fieldDescription = 'deep conv cloud/ice rate'
    onePostGrid%fieldUnits = 'g/kg/day'
    OutputField = OutputField*86400.0*1000.0

    ! output binary field

    call OutputGradsField (oneBramsGrid, onePostGrid, OutputField)
  end subroutine Brams2Post_cucldp




end module ModPostOneFieldUtils





