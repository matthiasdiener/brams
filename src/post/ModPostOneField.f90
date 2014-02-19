module ModPostOneField

  USE mem_grid ,     ONLY : time   ! INTENT(IN)  !DSM

  use ModBramsGrid, only: BramsGrid
  use ModPostGrid, only: PostGrid

  use ModPostOneFieldUtils, only: Brams2Post_u
  use ModPostOneFieldUtils, only: Brams2Post_v
  use ModPostOneFieldUtils, only: Brams2Post_tempk
  use ModPostOneFieldUtils, only: Brams2Post_tveg
  use ModPostOneFieldUtils, only: Brams2Post_totpcp
  use ModPostOneFieldUtils, only: Brams2Post_acccon
  use ModPostOneFieldUtils, only: Brams2Post_dewptc
  use ModPostOneFieldUtils, only: Brams2Post_rshort
  use ModPostOneFieldUtils, only: Brams2Post_rlong
  use ModPostOneFieldUtils, only: Brams2Post_sea_press
  use ModPostOneFieldUtils, only: Brams2Post_cape
  use ModPostOneFieldUtils, only: Brams2Post_cine
  use ModPostOneFieldUtils, only: Brams2Post_topo
  use ModPostOneFieldUtils, only: Brams2Post_precip
  use ModPostOneFieldUtils, only: Brams2Post_le
  use ModPostOneFieldUtils, only: Brams2Post_h
  use ModPostOneFieldUtils, only: Brams2Post_rv
  use ModPostOneFieldUtils, only: Brams2Post_rlongup
  use ModPostOneFieldUtils, only: Brams2Post_albedt
  use ModPostOneFieldUtils, only: Brams2Post_geo
  use ModPostOneFieldUtils, only: Brams2Post_ue_avg
  use ModPostOneFieldUtils, only: Brams2Post_ve_avg
  use ModPostOneFieldUtils, only: Brams2Post_tempc2m
  use ModPostOneFieldUtils, only: Brams2Post_tempc
  use ModPostOneFieldUtils, only: Brams2Post_rh
  use ModPostOneFieldUtils, only: Brams2Post_w
  use ModPostOneFieldUtils, only: Brams2Post_sst
  use ModPostOneFieldUtils, only: Brams2Post_land
  use ModPostOneFieldUtils, only: Brams2Post_smoist
  use ModPostOneFieldUtils, only: Brams2Post_zi
  use ModPostOneFieldUtils, only: Brams2Post_tke
  use ModPostOneFieldUtils, only: Brams2Post_cloud
  use ModPostOneFieldUtils, only: Brams2Post_omeg
  use ModPostOneFieldUtils, only: Brams2Post_pwt
  use ModPostOneFieldUtils, only: Brams2Post_slp_metar
  use ModPostOneFieldUtils, only: Brams2Post_td2m 
  use ModPostOneFieldUtils, only: Brams2Post_u10m
  use ModPostOneFieldUtils, only: Brams2Post_v10m 
  use ModPostOneFieldUtils, only: Brams2Post_vtype
  use ModPostOneFieldUtils, only: Brams2Post_sltex_p
  use ModPostOneFieldUtils, only: Brams2Post_ndvi
  use ModPostOneFieldUtils, only: Brams2Post_lai
  use ModPostOneFieldUtils, only: Brams2Post_u10mj
  use ModPostOneFieldUtils, only: Brams2Post_v10mj
  use ModPostOneFieldUtils, only: Brams2Post_t2mj
  use ModPostOneFieldUtils, only: Brams2Post_csj
  use ModPostOneFieldUtils, only: Brams2Post_rv2mj
  use ModPostOneFieldUtils, only: Brams2Post_td2mj
  use ModPostOneFieldUtils, only: Brams2Post_co2
  use ModPostOneFieldUtils, only: Brams2Post_co2_antro
  use ModPostOneFieldUtils, only: Brams2Post_co2_burn
  use ModPostOneFieldUtils, only: Brams2Post_co2_burn2D
  use ModPostOneFieldUtils, only: Brams2Post_co2_bioge
  use ModPostOneFieldUtils, only: Brams2Post_co2_total
  use ModPostOneFieldUtils, only: Brams2Post_theta

 !--(DMK-CCATT-INI)-------------------------------------------------------
  use ModPostOneFieldUtils, only: Brams2Post_CO
  use ModPostOneFieldUtils, only: Brams2Post_NO 
  use ModPostOneFieldUtils, only: Brams2Post_HNO3 
  use ModPostOneFieldUtils, only: Brams2Post_O3 
  use ModPostOneFieldUtils, only: Brams2Post_NO2 
  use ModPostOneFieldUtils, only: Brams2Post_PM25 
  use ModPostOneFieldUtils, only: Brams2Post_PM25wd
  use ModPostOneFieldUtils, only: Brams2Post_AOT550
  use ModPostOneFieldUtils, only: Brams2Post_PMINT
 !--(DMK-CCATT-FIM)-------------------------------------------------------
  use ModPostOneFieldUtils, only: Brams2Post_curtdp
  use ModPostOneFieldUtils, only: Brams2Post_cuthdp
  use ModPostOneFieldUtils, only: Brams2Post_cucldp
   
   
  implicit none

  private
  public :: PostOneField

contains


  subroutine PostOneField(varName, oneBramsGrid, onePostGrid)
    character(len=*), intent(in) :: varName
    type(BramsGrid), pointer :: oneBramsGrid
    type(PostGrid), pointer :: onePostGrid

    select case (varName)

    case ('u')
       call Brams2Post_u (varName, oneBramsGrid, onePostGrid)

    case ('v')
       call Brams2Post_v (varName, oneBramsGrid, onePostGrid)

    case ('tempk')
       call Brams2Post_tempk (varName, oneBramsGrid, onePostGrid)

    case ('tveg')
       call Brams2Post_tveg (varName, oneBramsGrid, onePostGrid)

    case ('tke')
       call Brams2Post_tke (varName, oneBramsGrid, onePostGrid)

    case ('totpcp')
       call Brams2Post_totpcp (varName, oneBramsGrid, onePostGrid)

    case ('acccon')
       call Brams2Post_acccon (varName, oneBramsGrid, onePostGrid)

    case ('dewptc')
       call Brams2Post_dewptc (varName, oneBramsGrid, onePostGrid)

    case ('cloud')
       call Brams2Post_cloud (varName, oneBramsGrid, onePostGrid)

    case ('rshort')
       call Brams2Post_rshort (varName, oneBramsGrid, onePostGrid)

    case ('rlong')
       call Brams2Post_rlong (varName, oneBramsGrid, onePostGrid)

    case ('sea_press')
       call Brams2Post_sea_press (varName, oneBramsGrid, onePostGrid)

    case ('cape')
       call Brams2Post_cape (varName, oneBramsGrid, onePostGrid)

    case ('cine')
       call Brams2Post_cine (varName, oneBramsGrid, onePostGrid)

    case ('topo')
       call Brams2Post_topo (varName, oneBramsGrid, onePostGrid)

    case ('precip')
       call Brams2Post_precip (varName, oneBramsGrid, onePostGrid)

    case ('le')
       call Brams2Post_le (varName, oneBramsGrid, onePostGrid)

    case ('h')
       call Brams2Post_h (varName, oneBramsGrid, onePostGrid)

    case ('rv')
       call Brams2Post_rv (varName, oneBramsGrid, onePostGrid)

    case ('rlongup')
       call Brams2Post_rlongup (varName, oneBramsGrid, onePostGrid)

    case ('albedt')
       call Brams2Post_albedt (varName, oneBramsGrid, onePostGrid)

    case ('geo')
       call Brams2Post_geo (varName, oneBramsGrid, onePostGrid)

    case ('ue_avg')
       call Brams2Post_ue_avg (varName, oneBramsGrid, onePostGrid)

    case ('ve_avg')
       call Brams2Post_ve_avg (varName, oneBramsGrid, onePostGrid)

    case ('tempc2m')
       call Brams2Post_tempc2m (varName, oneBramsGrid, onePostGrid)

    case ('tempc')
       call Brams2Post_tempc (varName, oneBramsGrid, onePostGrid)

    case ('rh')
       call Brams2Post_rh (varName, oneBramsGrid, onePostGrid)

    case ('w')
       call Brams2Post_w (varName, oneBramsGrid, onePostGrid)

    case ('omeg')
       call Brams2Post_omeg (varName, oneBramsGrid, onePostGrid)

    case ('sst')
       call Brams2Post_sst (varName, oneBramsGrid, onePostGrid)

    case ('land')
       call Brams2Post_land (varName, oneBramsGrid, onePostGrid)

    case ('smoist')
       call Brams2Post_smoist (varName, oneBramsGrid, onePostGrid)

    case ('zi')
       call Brams2Post_zi (varName, oneBramsGrid, onePostGrid)

    case ('pwt')
       call Brams2Post_pwt (varName, oneBramsGrid, onePostGrid)

    case ('slp_metar')
       call Brams2Post_slp_metar (varName, oneBramsGrid, onePostGrid)

    case ('td2m')
       call Brams2Post_td2m (varName, oneBramsGrid, onePostGrid)

    case ('u10m')
       call Brams2Post_u10m (varName, oneBramsGrid, onePostGrid)

    case ('v10m')
       call Brams2Post_v10m (varName, oneBramsGrid, onePostGrid)


    case ('vtype')
       call Brams2Post_vtype (varName, oneBramsGrid, onePostGrid)

    case ('sltex_p')
       call Brams2Post_sltex_p (varName, oneBramsGrid, onePostGrid)

    case ('ndvi')
       call Brams2Post_ndvi (varName, oneBramsGrid, onePostGrid)

    case ('lai')
       call Brams2Post_lai (varName, oneBramsGrid, onePostGrid)

    case ('u10mj')
       if (time==0.) then
         call Brams2Post_u10m (varName, oneBramsGrid, onePostGrid)
       else
         call Brams2Post_u10mj (varName, oneBramsGrid, onePostGrid)
       endif

    case ('v10mj')
       if (time==0.) then
         call Brams2Post_v10m (varName, oneBramsGrid, onePostGrid)
       else
         call Brams2Post_v10mj (varName, oneBramsGrid, onePostGrid)
       endif

    case ('t2mj')
       if (time==0.) then
         call Brams2Post_tempc2m (varName, oneBramsGrid, onePostGrid)
       else
         call Brams2Post_t2mj (varName, oneBramsGrid, onePostGrid)
       endif

    case ('csj')
       if (time==0.) then
         call Brams2Post_tempc2m (varName, oneBramsGrid, onePostGrid)
       else
         call Brams2Post_csj (varName, oneBramsGrid, onePostGrid)
       endif

    case ('rv2mj')
       if (time==0.) then
         call Brams2Post_rv (varName, oneBramsGrid, onePostGrid)
       else
         call Brams2Post_rv2mj (varName, oneBramsGrid, onePostGrid)
       endif

    case ('td2mj')
       if (time==0.) then
         call Brams2Post_td2m (varName, oneBramsGrid, onePostGrid)
       else
         call Brams2Post_td2mj (varName, oneBramsGrid, onePostGrid)
       endif

    case ('CO2')
       call Brams2Post_co2 (varName, oneBramsGrid, onePostGrid)

    case ('CO2_antro')
       call Brams2Post_co2_antro (varName, oneBramsGrid, onePostGrid)

    case ('CO2_burn')
       call Brams2Post_co2_burn (varName, oneBramsGrid, onePostGrid)

    case ('CO2_burn2D')
       call Brams2Post_co2_burn2D (varName, oneBramsGrid, onePostGrid)

    case ('CO2_bioge')
       call Brams2Post_co2_bioge (varName, oneBramsGrid, onePostGrid)

!---RB       
    case ('co2_total')
       call Brams2Post_co2_total (varName, oneBramsGrid, onePostGrid)
    
    case ('theta')
       call Brams2Post_theta (varName, oneBramsGrid, onePostGrid)
    

!--(DMK-CCATT-INI)-------------------------------------------------------
    case ('CO')
       call Brams2Post_CO (varName, oneBramsGrid, onePostGrid)

    case ('NO')   
       call Brams2Post_NO (varName, oneBramsGrid, onePostGrid)

    case ('HNO3')   
       call Brams2Post_HNO3 (varName, oneBramsGrid, onePostGrid)

    case ('O3')       
       call Brams2Post_O3 (varName, oneBramsGrid, onePostGrid)

    case ('NO2')   
       call Brams2Post_NO2(varName, oneBramsGrid, onePostGrid)

    case ('PM25')   
       call Brams2Post_PM25(varName, oneBramsGrid, onePostGrid)

    case ('PM25wd')   
       call Brams2Post_PM25wd(varName, oneBramsGrid, onePostGrid)

    case ('aot550')   
       call Brams2Post_AOT550(varName, oneBramsGrid, onePostGrid)

    case ('PMINT')   
       call Brams2Post_PMINT(varName, oneBramsGrid, onePostGrid)
!--(DMK-CCATT-FIM)-------------------------------------------------------

    case ('cuthdp')   
       call Brams2Post_cuthdp(varName, oneBramsGrid, onePostGrid)
    case ('curtdp')   
       call Brams2Post_curtdp(varName, oneBramsGrid, onePostGrid)
    case ('cucldp')   
       call Brams2Post_cucldp(varName, oneBramsGrid, onePostGrid)

    case default
       write(*,"(a)") "**(OnePostField)** Post field "//trim(varName)//" unknown to this procedure"
    end select

  end subroutine PostOneField
end module ModPostOneField
