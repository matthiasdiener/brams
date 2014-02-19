MODULE UkmoAdapt

  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers

  TYPE t_ukmo
     REAL(KIND=r8), DIMENSION(:), POINTER :: &
          ! Downward Surface SW flux visible beam (clear)    
          SfcVisBeamC &
          ! Downward Surface SW flux visible diffuse (clear)
          , SfcVisDiffC &
          ! Downward Surface SW flux Near-IR beam (clear)
          , SfcNirBeamC &
          ! Downward Surface SW flux Near-IR diffuse (clear)
          , SfcNirDiffC &
          ! Downward Surface SW flux visible beam (cloudy)
          , SfcVisBeam  &
          ! Downward Surface SW flux visible diffuse (cloudy)
          , SfcVisDiff  &
          ! Downward Surface SW flux Near-IR beam (cloudy)
          , SfcNirBeam  &
          ! Downward Surface SW flux Near-IR diffuse (cloudy)
          , SfcNirDiff  &
          ! Net Solar flux at top of atmosphere (cloudy)
          , ToaNet      &
          ! Net Solar flux at top of atmosphere (clear)
          , ToaNetC     &
          ! Net Solar flux at surface (cloudy)
          , SfcNet      &
          ! Net Solar flux at surface (clear)
          , SfcNetC     &
          ! Solar flux at top of atmosphere
          , ToaDown     &
          ! Upward TOA longwave flux (clear)
          , lw_toa_up_clr &
          ! Upward TOA longwave flux
          , lw_toa_up &
          ! Net Surface longwave flux (clear)
          , lw_sfc_net_clr &
          ! Net Surface longwave flux
          , lw_sfc_net &
          ! Downward Surface longwave flux (clear)
          , lw_sfc_down_clr &
          ! Downward Surface longwave flux
          , lw_sfc_down &
	  ! Cos of Zenital Angle
          , coszLocal &
	  ! Ground surface temperature (K)
          , gtg &
	  ! visible diffuse surface albedo
	  , AlbVisDiff &
	  ! near-ir diffuse surface albedo
          , AlbNirDiff &
	  ! visible beam surface albedo
          , AlbVisBeam &
	  ! near-ir beam surface albedo
          , AlbNirBeam
	  
     REAL(KIND=r8), DIMENSION(:,:), POINTER :: &
          ! Heating rate (clear case) (K/s)
          HeatRateC &
          ! Heating rate (cloudy case) (K/s)
          , HeatRate &
          ! Cooling rate (K/s, clear)
          , lw_cool_clr &
          ! Cooling rate (K/s)
          , lw_cool &
	  ! Temperature at middle of Layer (K)
          , FlipTe &
          ! Specific Humidity at middle of layer (g/g)
	  , FlipQe &
          ! Pressure at bottom of layers (mb)
	  , FlipPbot &
          ! Pressure difference bettween levels (mb)
	  , FlipDP &
          ! Pressure at Middle of Layer(mb)
	  , FlipPMid    &  
	  ! Ozone Mixing ratio at middle of layer (g/g)
          , FlipO3      &
	  ! Ice particle Effective Radius (microns)
          , FlipRei     &
	  ! Liquid particle Effective Radius (microns)
          , FlipRel     &
	  ! fractional amount of cloud that is ice
          , FlipFice    &
	  ! ice/water mixing ratio
          , FlipLMixR  &
	  ! Large scale cloud amount in layers
	  , cld &
	  ! Cumulus cloud amount in layers
	  , clu
	  
     REAL(KIND=r8) :: co2val, solar     
     INTEGER, DIMENSION(:), POINTER :: mapI,mapJ
     INTEGER(KIND=i8), DIMENSION(:), POINTER :: iMask
     REAL(KIND=r8), POINTER, dimension(:) :: sigmid
     REAL(KIND=r8), POINTER, dimension(:) :: sig 
     integer, POINTER, dimension(:) :: flip
     LOGICAL :: radUkmoNotInitiated=.true.
     REAL(KIND=r8),POINTER,DIMENSION(:) :: delsig
     INTEGER :: nls
  END TYPE t_ukmo
  TYPE(t_ukmo),ALLOCATABLE,DIMENSION(:) :: mem_UKMO

  LOGICAL,DIMENSION(:),ALLOCATABLE :: mem_Rad_UKMOCreated
  REAL(KIND=r8), allocatable, dimension(:,:,:) :: ozone
  INTEGER, PARAMETER :: nlm_getoz=18
  LOGICAL :: first_getoz,inter_getoz
  INTEGER :: year_getoz,mon_getoz
  INTEGER, PARAMETER :: nl=37
  INTEGER, PARAMETER :: ns=4
  REAL(KIND=r8) :: ozsig(18)
  REAL(KIND=r8)  :: yrl=365.2500_r8
  LOGICAL :: ukmo_created=.false.

CONTAINS

  SUBROUTINE Create_mem_ukmo(ngrids)
    INTEGER,INTENT(IN) :: ngrids

    IF(.not. allocated(mem_UKMO)) THEN
      allocate(mem_UKMO(ngrids),mem_Rad_UKMOCreated(ngrids))
      !Rad_UKMOCreated=.false.
     END IF
    
  END SUBROUTINE Create_mem_ukmo

  SUBROUTINE Destroy_mem_ukmo()
    
    IF(allocated(mem_UKMO)) deallocate(mem_UKMO)

  END SUBROUTINE Destroy_mem_ukmo

  SUBROUTINE Alloc_mem_Rad_UKMO(ngrid,ncols,kmax)
    INTEGER,INTENT(in) :: ncols,kmax,ngrid

    IF(.not. mem_Rad_UKMOCreated(ngrid)) THEN
       allocate(&
                mem_UKMO(ngrid)%sigmid(kmax+1)    &
              , mem_UKMO(ngrid)%sig   (kmax+1)    &
              , mem_UKMO(ngrid)%flip  (kmax+1)    &
              , mem_UKMO(ngrid)%delsig(kMax+1)    &
              )

       allocate(&
            mem_UKMO(ngrid)%SfcVisBeamC      (ncols) &
            , mem_UKMO(ngrid)%SfcVisDiffC    (ncols) &
            , mem_UKMO(ngrid)%SfcNirBeamC    (ncols) &
            , mem_UKMO(ngrid)%SfcNirDiffC    (ncols) &
            , mem_UKMO(ngrid)%SfcVisBeam     (ncols) &
            , mem_UKMO(ngrid)%SfcVisDiff     (ncols) &
            , mem_UKMO(ngrid)%SfcNirBeam     (ncols) &
            , mem_UKMO(ngrid)%SfcNirDiff     (ncols) &
            , mem_UKMO(ngrid)%ToaNet         (ncols) &
            , mem_UKMO(ngrid)%ToaNetC        (ncols) &
            , mem_UKMO(ngrid)%SfcNet         (ncols) &
            , mem_UKMO(ngrid)%SfcNetC        (ncols) &
            , mem_UKMO(ngrid)%ToaDown        (ncols) &
            , mem_UKMO(ngrid)%lw_toa_up_clr  (ncols) &
            , mem_UKMO(ngrid)%lw_toa_up      (ncols) &
            , mem_UKMO(ngrid)%lw_sfc_net_clr (ncols) &
            , mem_UKMO(ngrid)%lw_sfc_net     (ncols) &
            , mem_UKMO(ngrid)%lw_sfc_down_clr(ncols) &
            , mem_UKMO(ngrid)%lw_sfc_down    (ncols) &
            , mem_UKMO(ngrid)%MapI           (ncols) &
            , mem_UKMO(ngrid)%MapJ           (ncols) &
            , mem_UKMO(ngrid)%cosZLocal      (ncols) &
            , mem_UKMO(ngrid)%gtg            (ncols) &
            , mem_UKMO(ngrid)%iMask          (ncols) &
            , mem_UKMO(ngrid)%AlbVisDiff     (ncols) &
            , mem_UKMO(ngrid)%AlbNirDiff     (ncols) &
            , mem_UKMO(ngrid)%AlbVisBeam     (ncols) &
            , mem_UKMO(ngrid)%AlbNirBeam     (ncols) &
            )

       allocate( &
              mem_UKMO(ngrid)%HeatRateC  (ncols,kmax) &
            , mem_UKMO(ngrid)%HeatRate   (ncols,kmax) &
            , mem_UKMO(ngrid)%lw_cool_clr(ncols,kmax) &
            , mem_UKMO(ngrid)%lw_cool    (ncols,kmax) &
            , mem_UKMO(ngrid)%FlipTe     (ncols,kmax) &
	    , mem_UKMO(ngrid)%FlipQe     (ncols,kmax) &
	    , mem_UKMO(ngrid)%FlipPbot   (ncols,kmax) &
	    , mem_UKMO(ngrid)%FlipDP     (ncols,kmax) &
	    , mem_UKMO(ngrid)%FlipPMid   (ncols,kmax) &  
            , mem_UKMO(ngrid)%FlipO3     (ncols,kmax) &
            , mem_UKMO(ngrid)%FlipRei    (ncols,kmax) &
            , mem_UKMO(ngrid)%FlipRel    (ncols,kmax) &
            , mem_UKMO(ngrid)%FlipFice   (ncols,kmax) &
            , mem_UKMO(ngrid)%FlipLMixR  (ncols,kmax) &
            , mem_UKMO(ngrid)%cld        (ncols,kmax) &
            , mem_UKMO(ngrid)%clu        (ncols,kmax) &
            )

       mem_UKMO(ngrid)%SfcVisBeamC=0.0_r8
       mem_UKMO(ngrid)%SfcVisDiffC=0.0_r8
       mem_UKMO(ngrid)%SfcNirBeamC=0.0_r8
       mem_UKMO(ngrid)%SfcNirDiffC=0.0_r8
       mem_UKMO(ngrid)%SfcVisBeam =0.0_r8
       mem_UKMO(ngrid)%SfcVisDiff =0.0_r8
       mem_UKMO(ngrid)%SfcNirBeam =0.0_r8
       mem_UKMO(ngrid)%SfcNirDiff =0.0_r8
       mem_UKMO(ngrid)%ToaNet     =0.0_r8
       mem_UKMO(ngrid)%ToaNetC    =0.0_r8
       mem_UKMO(ngrid)%SfcNet     =0.0_r8
       mem_UKMO(ngrid)%SfcNetC    =0.0_r8
       mem_UKMO(ngrid)%ToaDown    =0.0_r8
       mem_UKMO(ngrid)%HeatRateC  =0.0_r8
       mem_UKMO(ngrid)%HeatRate   =0.0_r8
       !Long wave
       mem_UKMO(ngrid)%lw_toa_up_clr  =0.0_r8
       mem_UKMO(ngrid)%lw_toa_up      =0.0_r8
       mem_UKMO(ngrid)%lw_sfc_net_clr =0.0_r8
       mem_UKMO(ngrid)%lw_sfc_net     =0.0_r8
       mem_UKMO(ngrid)%lw_sfc_down_clr=0.0_r8
       mem_UKMO(ngrid)%lw_sfc_down    =0.0_r8
       mem_UKMO(ngrid)%lw_cool_clr    =0.0_r8
       mem_UKMO(ngrid)%lw_cool        =0.0_r8
       mem_UKMO(ngrid)%cosZLocal      =0.0_r8
       mem_UKMO(ngrid)%FlipTe	      =0.0_r8
       mem_UKMO(ngrid)%FlipQe	      =0.0_r8
       mem_UKMO(ngrid)%FlipPbot       =0.0_r8
       mem_UKMO(ngrid)%FlipDP	      =0.0_r8
       mem_UKMO(ngrid)%FlipPMid       =0.0_r8
       mem_UKMO(ngrid)%AlbVisDiff     =0.0_r8
       mem_UKMO(ngrid)%AlbNirDiff     =0.0_r8
       mem_UKMO(ngrid)%AlbVisBeam     =0.0_r8
       mem_UKMO(ngrid)%AlbNirBeam     =0.0_r8  
       mem_UKMO(ngrid)%FlipO3         =0.0_r8
       mem_UKMO(ngrid)%FlipRei        =0.0_r8
       mem_UKMO(ngrid)%FlipRel        =0.0_r8
       mem_UKMO(ngrid)%FlipFice       =0.0_r8
       mem_UKMO(ngrid)%FlipLMixR      =0.0_r8 
       mem_UKMO(ngrid)%cld            =0.0_r8
       mem_UKMO(ngrid)%clu            =0.0_r8 
       mem_UKMO(ngrid)%gtg  	  =0.0_r8
       mem_UKMO(ngrid)%MapI  =0
       mem_UKMO(ngrid)%MapJ =0
       mem_UKMO(ngrid)%iMask =0

       mem_Rad_UKMOCreated(ngrid)=.true.	  
       mem_UKMO(ngrid)%sigmid=0.0
       mem_UKMO(ngrid)%sig=0.0
       
    END IF
  END SUBROUTINE Alloc_mem_Rad_UKMO

  SUBROUTINE Dealloc_mem_Rad_UKMO(ngrid)
    INTEGER,INTENT(IN) :: ngrid
    IF(mem_Rad_UKMOCreated(ngrid)) THEN
       deallocate(&
            mem_UKMO(ngrid)%SfcVisBeamC       &
            , mem_UKMO(ngrid)%SfcVisDiffC     &
            , mem_UKMO(ngrid)%SfcNirBeamC     &
            , mem_UKMO(ngrid)%SfcNirDiffC     &
            , mem_UKMO(ngrid)%SfcVisBeam      &
            , mem_UKMO(ngrid)%SfcVisDiff      &
            , mem_UKMO(ngrid)%SfcNirBeam      &
            , mem_UKMO(ngrid)%SfcNirDiff      &
            , mem_UKMO(ngrid)%ToaNet          &
            , mem_UKMO(ngrid)%ToaNetC         &
            , mem_UKMO(ngrid)%SfcNet          &
            , mem_UKMO(ngrid)%SfcNetC         &
            , mem_UKMO(ngrid)%ToaDown         &
            , mem_UKMO(ngrid)%lw_toa_up_clr   &
            , mem_UKMO(ngrid)%lw_toa_up       &
            , mem_UKMO(ngrid)%lw_sfc_net_clr  &
            , mem_UKMO(ngrid)%lw_sfc_net      &
            , mem_UKMO(ngrid)%lw_sfc_down_clr &
            , mem_UKMO(ngrid)%lw_sfc_down     &
            , mem_UKMO(ngrid)%MapI            &
            , mem_UKMO(ngrid)%MapJ            &
            , mem_UKMO(ngrid)%cosZLocal       &
            , mem_UKMO(ngrid)%gtg             &
            , mem_UKMO(ngrid)%iMask           &
            , mem_UKMO(ngrid)%AlbVisDiff      &
            , mem_UKMO(ngrid)%AlbNirDiff      &
            , mem_UKMO(ngrid)%AlbVisBeam      &
            , mem_UKMO(ngrid)%AlbNirBeam      &
            )

       deallocate( &
              mem_UKMO(ngrid)%HeatRateC   &
            , mem_UKMO(ngrid)%HeatRate    &
            , mem_UKMO(ngrid)%lw_cool_clr &
            , mem_UKMO(ngrid)%lw_cool     &
            , mem_UKMO(ngrid)%FlipTe      &
	    , mem_UKMO(ngrid)%FlipQe      &
	    , mem_UKMO(ngrid)%FlipPbot    &
	    , mem_UKMO(ngrid)%FlipDP      &
	    , mem_UKMO(ngrid)%FlipPMid    &  
            , mem_UKMO(ngrid)%FlipO3      &
            , mem_UKMO(ngrid)%FlipRei     &
            , mem_UKMO(ngrid)%FlipRel     &
            , mem_UKMO(ngrid)%FlipFice    &
            , mem_UKMO(ngrid)%FlipLMixR   &
            , mem_UKMO(ngrid)%cld         &
            , mem_UKMO(ngrid)%clu         &
            )

       mem_Rad_UKMOCreated(ngrid)=.false.

    END IF

  END SUBROUTINE Dealloc_mem_Rad_UKMO

  SUBROUTINE InitGetoz(yrl,kmax,sl)
    IMPLICIT NONE
    REAL(KIND=r8),    INTENT(IN   ) :: yrl
    INTEGER, INTENT(IN   ) :: kmax
    REAL(KIND=r8),    INTENT(in   ) :: sl (kmax)
    INTEGER                :: l, ll


    ALLOCATE(ozone(nlm_getoz,nl,ns))

    !
    !     four season climatological ozone data in nmc sigma layers
    !
    !     for seasonal variation
    !     season=1 - winter          season=2 - spring
    !     season=3 - summer          season=4 - fall
    !     unit of ozone mixing ratio is in ( 10**-4 g/g ).  the data is
    !     in 18 sigma layers from top to bottom.  for every layer, there
    !     are 37 latitudes at 5 degree interval from north pole to south
    !     pole.
    !     mrf86 18 layers
    !
    !
    !     1. winter
    !
    !     wint1(18,6)
    !
    ozone(1:18, 1:6, 1) = RESHAPE( (/ &
         .068467e0_r8,.052815e0_r8,.035175e0_r8,.022334e0_r8,.013676e0_r8,.007363e0_r8, &
         .003633e0_r8,.001582e0_r8,.001111e0_r8,.000713e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .069523e0_r8,.052249e0_r8,.034255e0_r8,.021379e0_r8,.012306e0_r8,.006727e0_r8, &
         .003415e0_r8,.001578e0_r8,.001072e0_r8,.000681e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .070579e0_r8,.051684e0_r8,.033335e0_r8,.020423e0_r8,.010935e0_r8,.006091e0_r8, &
         .003197e0_r8,.001573e0_r8,.001034e0_r8,.000650e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .074885e0_r8,.049987e0_r8,.030140e0_r8,.017894e0_r8,.009881e0_r8,.005543e0_r8, &
         .002907e0_r8,.001379e0_r8,.000961e0_r8,.000644e0_r8,.000512e0_r8,.000463e0_r8, &
         .000451e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8, &
         .079190e0_r8,.048290e0_r8,.026945e0_r8,.015366e0_r8,.008826e0_r8,.004995e0_r8, &
         .002616e0_r8,.001184e0_r8,.000887e0_r8,.000637e0_r8,.000508e0_r8,.000486e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .082443e0_r8,.047591e0_r8,.025358e0_r8,.014294e0_r8,.008233e0_r8,.004664e0_r8, &
         .002430e0_r8,.001068e0_r8,.000851e0_r8,.000644e0_r8,.000508e0_r8,.000474e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8/), &
         (/18,6/))
    !
    !     wint2(18,6)
    !
    ozone(1:18, 7:12, 1) = RESHAPE( (/ &
         .085695e0_r8,.046892e0_r8,.023772e0_r8,.013223e0_r8,.007640e0_r8,.004333e0_r8, &
         .002244e0_r8,.000951e0_r8,.000815e0_r8,.000650e0_r8,.000508e0_r8,.000463e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .089618e0_r8,.042869e0_r8,.019963e0_r8,.010502e0_r8,.005966e0_r8,.003525e0_r8, &
         .001936e0_r8,.000906e0_r8,.000769e0_r8,.000625e0_r8,.000508e0_r8,.000452e0_r8, &
         .000451e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8, &
         .093540e0_r8,.038846e0_r8,.016155e0_r8,.007781e0_r8,.004292e0_r8,.002716e0_r8, &
         .001628e0_r8,.000862e0_r8,.000724e0_r8,.000600e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .097097e0_r8,.034916e0_r8,.012983e0_r8,.006240e0_r8,.003666e0_r8,.002259e0_r8, &
         .001336e0_r8,.000730e0_r8,.000629e0_r8,.000549e0_r8,.000499e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .100654e0_r8,.030986e0_r8,.009812e0_r8,.004698e0_r8,.003041e0_r8,.001803e0_r8, &
         .001044e0_r8,.000599e0_r8,.000533e0_r8,.000499e0_r8,.000491e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .101724e0_r8,.026500e0_r8,.007228e0_r8,.003391e0_r8,.002058e0_r8,.001285e0_r8, &
         .000811e0_r8,.000531e0_r8,.000478e0_r8,.000449e0_r8,.000440e0_r8,.000421e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint3(18,6)
    !
    ozone(1:18, 13:18, 1) = RESHAPE( (/ &
         .102794e0_r8,.022015e0_r8,.004645e0_r8,.002084e0_r8,.001076e0_r8,.000767e0_r8, &
         .000577e0_r8,.000463e0_r8,.000423e0_r8,.000399e0_r8,.000389e0_r8,.000401e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .103456e0_r8,.018235e0_r8,.003195e0_r8,.001379e0_r8,.000771e0_r8,.000585e0_r8, &
         .000474e0_r8,.000411e0_r8,.000380e0_r8,.000362e0_r8,.000343e0_r8,.000348e0_r8, &
         .000346e0_r8,.000328e0_r8,.000317e0_r8,.000305e0_r8,.000302e0_r8,.000302e0_r8, &
         .104118e0_r8,.014455e0_r8,.001745e0_r8,.000674e0_r8,.000467e0_r8,.000403e0_r8, &
         .000370e0_r8,.000359e0_r8,.000337e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104106e0_r8,.012997e0_r8,.001479e0_r8,.000639e0_r8,.000468e0_r8,.000422e0_r8, &
         .000392e0_r8,.000372e0_r8,.000342e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104093e0_r8,.011539e0_r8,.001213e0_r8,.000604e0_r8,.000468e0_r8,.000442e0_r8, &
         .000414e0_r8,.000385e0_r8,.000347e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104087e0_r8,.010726e0_r8,.000971e0_r8,.000538e0_r8,.000440e0_r8,.000434e0_r8, &
         .000418e0_r8,.000397e0_r8,.000375e0_r8,.000343e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint4(18,6)
    !
    ozone(1:18, 19:24, 1) = RESHAPE( (/ &
         .102665e0_r8,.010977e0_r8,.001237e0_r8,.000590e0_r8,.000498e0_r8,.000479e0_r8, &
         .000458e0_r8,.000436e0_r8,.000421e0_r8,.000387e0_r8,.000326e0_r8,.000298e0_r8, &
         .000246e0_r8,.000227e0_r8,.000211e0_r8,.000200e0_r8,.000194e0_r8,.000186e0_r8, &
         .100892e0_r8,.012873e0_r8,.001886e0_r8,.000785e0_r8,.000643e0_r8,.000568e0_r8, &
         .000519e0_r8,.000487e0_r8,.000471e0_r8,.000437e0_r8,.000368e0_r8,.000305e0_r8, &
         .000201e0_r8,.000151e0_r8,.000117e0_r8,.000098e0_r8,.000090e0_r8,.000093e0_r8, &
         .100534e0_r8,.013704e0_r8,.002028e0_r8,.000861e0_r8,.000701e0_r8,.000604e0_r8, &
         .000546e0_r8,.000513e0_r8,.000504e0_r8,.000462e0_r8,.000381e0_r8,.000307e0_r8, &
         .000201e0_r8,.000151e0_r8,.000117e0_r8,.000098e0_r8,.000090e0_r8,.000093e0_r8, &
         .100218e0_r8,.015035e0_r8,.002537e0_r8,.001037e0_r8,.000790e0_r8,.000726e0_r8, &
         .000673e0_r8,.000628e0_r8,.000579e0_r8,.000512e0_r8,.000440e0_r8,.000374e0_r8, &
         .000307e0_r8,.000253e0_r8,.000227e0_r8,.000208e0_r8,.000194e0_r8,.000186e0_r8, &
         .099903e0_r8,.016365e0_r8,.003045e0_r8,.001214e0_r8,.000879e0_r8,.000848e0_r8, &
         .000801e0_r8,.000744e0_r8,.000654e0_r8,.000562e0_r8,.000499e0_r8,.000441e0_r8, &
         .000410e0_r8,.000358e0_r8,.000342e0_r8,.000322e0_r8,.000302e0_r8,.000302e0_r8, &
         .099547e0_r8,.017725e0_r8,.003693e0_r8,.001578e0_r8,.001125e0_r8,.000985e0_r8, &
         .000879e0_r8,.000795e0_r8,.000712e0_r8,.000643e0_r8,.000584e0_r8,.000521e0_r8, &
         .000482e0_r8,.000384e0_r8,.000351e0_r8,.000322e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint5(18,6)
    !
    ozone(1:18, 25:30, 1) = RESHAPE( (/ &
         .099191e0_r8,.019085e0_r8,.004340e0_r8,.001943e0_r8,.001371e0_r8,.001122e0_r8, &
         .000957e0_r8,.000847e0_r8,.000770e0_r8,.000724e0_r8,.000669e0_r8,.000601e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .098107e0_r8,.020617e0_r8,.004758e0_r8,.002137e0_r8,.001516e0_r8,.001211e0_r8, &
         .000999e0_r8,.000848e0_r8,.000778e0_r8,.000730e0_r8,.000677e0_r8,.000603e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .097023e0_r8,.022148e0_r8,.005177e0_r8,.002332e0_r8,.001660e0_r8,.001300e0_r8, &
         .001041e0_r8,.000849e0_r8,.000786e0_r8,.000737e0_r8,.000686e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .093464e0_r8,.026177e0_r8,.008525e0_r8,.003892e0_r8,.002452e0_r8,.001609e0_r8, &
         .001116e0_r8,.000851e0_r8,.000809e0_r8,.000762e0_r8,.000690e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .089906e0_r8,.030206e0_r8,.011873e0_r8,.005453e0_r8,.003244e0_r8,.001918e0_r8, &
         .001192e0_r8,.000852e0_r8,.000832e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .080939e0_r8,.032414e0_r8,.014163e0_r8,.007241e0_r8,.004328e0_r8,.002522e0_r8, &
         .001481e0_r8,.000934e0_r8,.000861e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint6(18,6)
    !
    ozone(1:18, 31:36, 1) = RESHAPE( (/ &
         .071972e0_r8,.034622e0_r8,.016453e0_r8,.009029e0_r8,.005413e0_r8,.003127e0_r8, &
         .001770e0_r8,.001015e0_r8,.000890e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .069820e0_r8,.035028e0_r8,.016929e0_r8,.009389e0_r8,.005645e0_r8,.003260e0_r8, &
         .001843e0_r8,.001055e0_r8,.000905e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .067669e0_r8,.035434e0_r8,.017406e0_r8,.009749e0_r8,.005876e0_r8,.003393e0_r8, &
         .001916e0_r8,.001094e0_r8,.000920e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .065518e0_r8,.035975e0_r8,.017854e0_r8,.010100e0_r8,.006534e0_r8,.003985e0_r8, &
         .002321e0_r8,.001240e0_r8,.000966e0_r8,.000774e0_r8,.000640e0_r8,.000548e0_r8, &
         .000479e0_r8,.000384e0_r8,.000346e0_r8,.000316e0_r8,.000302e0_r8,.000302e0_r8, &
         .063367e0_r8,.036516e0_r8,.018302e0_r8,.010452e0_r8,.007192e0_r8,.004577e0_r8, &
         .002727e0_r8,.001387e0_r8,.001012e0_r8,.000762e0_r8,.000585e0_r8,.000490e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .061216e0_r8,.037359e0_r8,.019151e0_r8,.010633e0_r8,.006845e0_r8,.004382e0_r8, &
         .002691e0_r8,.001511e0_r8,.001061e0_r8,.000749e0_r8,.000568e0_r8,.000465e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    !
    !     wint7(18)
    !
    ozone(1:18, 37, 1) = (/ &
         .059066e0_r8,.038201e0_r8,.019999e0_r8,.010813e0_r8,.006498e0_r8,.004188e0_r8, &
         .002656e0_r8,.001636e0_r8,.001110e0_r8,.000737e0_r8,.000551e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/)
    !
    !     2. spring
    !
    ozone(1:18, 1:6, 2) = RESHAPE( (/ &
         .074229e0_r8,.050084e0_r8,.030930e0_r8,.018676e0_r8,.011965e0_r8,.008165e0_r8, &
         .005428e0_r8,.003399e0_r8,.002098e0_r8,.001138e0_r8,.000780e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .074927e0_r8,.049459e0_r8,.029215e0_r8,.018025e0_r8,.011754e0_r8,.007786e0_r8, &
         .004972e0_r8,.002926e0_r8,.001817e0_r8,.001025e0_r8,.000758e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .075625e0_r8,.048835e0_r8,.027500e0_r8,.017375e0_r8,.011544e0_r8,.007407e0_r8, &
         .004516e0_r8,.002453e0_r8,.001536e0_r8,.000912e0_r8,.000737e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .077409e0_r8,.048159e0_r8,.026661e0_r8,.016596e0_r8,.010962e0_r8,.006972e0_r8, &
         .004160e0_r8,.002132e0_r8,.001391e0_r8,.000868e0_r8,.000686e0_r8,.000601e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .079194e0_r8,.047483e0_r8,.025822e0_r8,.015818e0_r8,.010380e0_r8,.006537e0_r8, &
         .003804e0_r8,.001811e0_r8,.001245e0_r8,.000825e0_r8,.000635e0_r8,.000570e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .084591e0_r8,.046553e0_r8,.025037e0_r8,.015156e0_r8,.009841e0_r8,.006124e0_r8, &
         .003534e0_r8,.001693e0_r8,.001170e0_r8,.000793e0_r8,.000631e0_r8,.000537e0_r8, &
         .000551e0_r8,.000509e0_r8,.000486e0_r8,.000516e0_r8,.000548e0_r8,.000446e0_r8/), &
         (/18,6/))
    ozone(1:18, 7:12, 2) = RESHAPE( (/ &
         .089988e0_r8,.045622e0_r8,.024253e0_r8,.014495e0_r8,.009303e0_r8,.005711e0_r8, &
         .003264e0_r8,.001574e0_r8,.001096e0_r8,.000762e0_r8,.000627e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .092863e0_r8,.042419e0_r8,.020704e0_r8,.012034e0_r8,.007417e0_r8,.004504e0_r8, &
         .002590e0_r8,.001334e0_r8,.000977e0_r8,.000731e0_r8,.000622e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .095737e0_r8,.039215e0_r8,.017155e0_r8,.009572e0_r8,.005532e0_r8,.003296e0_r8, &
         .001916e0_r8,.001094e0_r8,.000858e0_r8,.000699e0_r8,.000618e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .097501e0_r8,.035382e0_r8,.014856e0_r8,.008207e0_r8,.004619e0_r8,.002720e0_r8, &
         .001610e0_r8,.001012e0_r8,.000829e0_r8,.000687e0_r8,.000610e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .099264e0_r8,.031548e0_r8,.012557e0_r8,.006841e0_r8,.003705e0_r8,.002144e0_r8, &
         .001304e0_r8,.000930e0_r8,.000799e0_r8,.000675e0_r8,.000601e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .101718e0_r8,.026523e0_r8,.008473e0_r8,.004382e0_r8,.002392e0_r8,.001505e0_r8, &
         .001036e0_r8,.000836e0_r8,.000727e0_r8,.000618e0_r8,.000550e0_r8,.000494e0_r8, &
         .000501e0_r8,.000479e0_r8,.000473e0_r8,.000509e0_r8,.000541e0_r8,.000445e0_r8/), &
         (/18,6/))
    ozone(1:18, 13:18, 2) = RESHAPE( (/ &
         .104172e0_r8,.021499e0_r8,.004389e0_r8,.001922e0_r8,.001078e0_r8,.000865e0_r8, &
         .000767e0_r8,.000743e0_r8,.000654e0_r8,.000562e0_r8,.000499e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .104145e0_r8,.018082e0_r8,.003274e0_r8,.001493e0_r8,.000919e0_r8,.000762e0_r8, &
         .000678e0_r8,.000641e0_r8,.000584e0_r8,.000531e0_r8,.000495e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .104118e0_r8,.014665e0_r8,.002159e0_r8,.001063e0_r8,.000759e0_r8,.000659e0_r8, &
         .000589e0_r8,.000539e0_r8,.000514e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .107719e0_r8,.013052e0_r8,.001822e0_r8,.000953e0_r8,.000701e0_r8,.000604e0_r8, &
         .000551e0_r8,.000525e0_r8,.000509e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .111320e0_r8,.011439e0_r8,.001485e0_r8,.000843e0_r8,.000642e0_r8,.000549e0_r8, &
         .000512e0_r8,.000512e0_r8,.000504e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .112375e0_r8,.011255e0_r8,.001357e0_r8,.000744e0_r8,.000585e0_r8,.000533e0_r8, &
         .000512e0_r8,.000512e0_r8,.000504e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8/), &
         (/18,6/))
    ozone(1:18, 19:24, 2) = RESHAPE( (/ &
         .109850e0_r8,.010424e0_r8,.001079e0_r8,.000567e0_r8,.000498e0_r8,.000479e0_r8, &
         .000463e0_r8,.000448e0_r8,.000418e0_r8,.000399e0_r8,.000389e0_r8,.000367e0_r8, &
         .000351e0_r8,.000328e0_r8,.000320e0_r8,.000337e0_r8,.000355e0_r8,.000304e0_r8, &
         .107002e0_r8,.009961e0_r8,.001025e0_r8,.000533e0_r8,.000497e0_r8,.000460e0_r8, &
         .000422e0_r8,.000385e0_r8,.000332e0_r8,.000300e0_r8,.000288e0_r8,.000249e0_r8, &
         .000202e0_r8,.000158e0_r8,.000132e0_r8,.000114e0_r8,.000104e0_r8,.000093e0_r8, &
         .107735e0_r8,.010146e0_r8,.001120e0_r8,.000576e0_r8,.000526e0_r8,.000477e0_r8, &
         .000430e0_r8,.000385e0_r8,.000332e0_r8,.000300e0_r8,.000288e0_r8,.000249e0_r8, &
         .000202e0_r8,.000158e0_r8,.000132e0_r8,.000114e0_r8,.000104e0_r8,.000093e0_r8, &
         .107021e0_r8,.012233e0_r8,.001533e0_r8,.000643e0_r8,.000556e0_r8,.000505e0_r8, &
         .000471e0_r8,.000448e0_r8,.000403e0_r8,.000362e0_r8,.000355e0_r8,.000296e0_r8, &
         .000251e0_r8,.000207e0_r8,.000180e0_r8,.000161e0_r8,.000152e0_r8,.000140e0_r8, &
         .106308e0_r8,.014320e0_r8,.001946e0_r8,.000709e0_r8,.000585e0_r8,.000533e0_r8, &
         .000512e0_r8,.000512e0_r8,.000473e0_r8,.000425e0_r8,.000423e0_r8,.000342e0_r8, &
         .000301e0_r8,.000257e0_r8,.000232e0_r8,.000212e0_r8,.000205e0_r8,.000209e0_r8, &
         .100592e0_r8,.015718e0_r8,.002411e0_r8,.001007e0_r8,.000802e0_r8,.000642e0_r8, &
         .000559e0_r8,.000526e0_r8,.000501e0_r8,.000474e0_r8,.000470e0_r8,.000439e0_r8, &
         .000430e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 25:30, 2) = RESHAPE( (/ &
         .094877e0_r8,.017116e0_r8,.002877e0_r8,.001305e0_r8,.001018e0_r8,.000751e0_r8, &
         .000606e0_r8,.000539e0_r8,.000529e0_r8,.000524e0_r8,.000516e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .094163e0_r8,.020198e0_r8,.004594e0_r8,.001772e0_r8,.001077e0_r8,.000806e0_r8, &
         .000649e0_r8,.000565e0_r8,.000547e0_r8,.000537e0_r8,.000521e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .093449e0_r8,.023279e0_r8,.006312e0_r8,.002240e0_r8,.001135e0_r8,.000862e0_r8, &
         .000692e0_r8,.000591e0_r8,.000564e0_r8,.000549e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .089886e0_r8,.026029e0_r8,.008558e0_r8,.003312e0_r8,.001655e0_r8,.001124e0_r8, &
         .000807e0_r8,.000631e0_r8,.000602e0_r8,.000568e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .086323e0_r8,.028778e0_r8,.010805e0_r8,.004383e0_r8,.002175e0_r8,.001386e0_r8, &
         .000923e0_r8,.000671e0_r8,.000640e0_r8,.000587e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .082715e0_r8,.031096e0_r8,.013350e0_r8,.006131e0_r8,.003205e0_r8,.002043e0_r8, &
         .001304e0_r8,.000842e0_r8,.000734e0_r8,.000631e0_r8,.000555e0_r8,.000494e0_r8, &
         .000480e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8/), &
         (/18,6/))
    ozone(1:18, 31:36, 2) = RESHAPE( (/ &
         .079108e0_r8,.033415e0_r8,.015895e0_r8,.007878e0_r8,.004234e0_r8,.002700e0_r8, &
         .001686e0_r8,.001014e0_r8,.000829e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .074807e0_r8,.034651e0_r8,.017056e0_r8,.008574e0_r8,.004769e0_r8,.002986e0_r8, &
         .001827e0_r8,.001079e0_r8,.000853e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .070506e0_r8,.035887e0_r8,.018218e0_r8,.009270e0_r8,.005304e0_r8,.003271e0_r8, &
         .001969e0_r8,.001145e0_r8,.000878e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .067669e0_r8,.037799e0_r8,.019680e0_r8,.009612e0_r8,.005481e0_r8,.003476e0_r8, &
         .002093e0_r8,.001123e0_r8,.000837e0_r8,.000631e0_r8,.000546e0_r8,.000447e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .064832e0_r8,.039712e0_r8,.021142e0_r8,.009954e0_r8,.005658e0_r8,.003681e0_r8, &
         .002218e0_r8,.001100e0_r8,.000796e0_r8,.000587e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .063734e0_r8,.039842e0_r8,.022004e0_r8,.010859e0_r8,.005712e0_r8,.003589e0_r8, &
         .002155e0_r8,.001174e0_r8,.000856e0_r8,.000612e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 37, 2) = (/ &
         .062636e0_r8,.039972e0_r8,.022867e0_r8,.011765e0_r8,.005766e0_r8,.003498e0_r8, &
         .002092e0_r8,.001248e0_r8,.000917e0_r8,.000637e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/)
    !
    !     3. summer
    !
    ozone(1:18, 1:6, 3) = RESHAPE( (/ &
         .059066e0_r8,.038201e0_r8,.019999e0_r8,.010813e0_r8,.006498e0_r8,.004188e0_r8, &
         .002656e0_r8,.001636e0_r8,.001110e0_r8,.000737e0_r8,.000551e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .061216e0_r8,.037359e0_r8,.019151e0_r8,.010633e0_r8,.006845e0_r8,.004382e0_r8, &
         .002691e0_r8,.001511e0_r8,.001061e0_r8,.000749e0_r8,.000568e0_r8,.000465e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .063367e0_r8,.036516e0_r8,.018302e0_r8,.010452e0_r8,.007192e0_r8,.004577e0_r8, &
         .002727e0_r8,.001387e0_r8,.001012e0_r8,.000762e0_r8,.000585e0_r8,.000490e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .065518e0_r8,.035975e0_r8,.017854e0_r8,.010100e0_r8,.006534e0_r8,.003985e0_r8, &
         .002321e0_r8,.001240e0_r8,.000966e0_r8,.000774e0_r8,.000640e0_r8,.000548e0_r8, &
         .000479e0_r8,.000384e0_r8,.000346e0_r8,.000316e0_r8,.000302e0_r8,.000302e0_r8, &
         .067669e0_r8,.035434e0_r8,.017406e0_r8,.009749e0_r8,.005876e0_r8,.003393e0_r8, &
         .001916e0_r8,.001094e0_r8,.000920e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .069820e0_r8,.035028e0_r8,.016929e0_r8,.009389e0_r8,.005645e0_r8,.003260e0_r8, &
         .001843e0_r8,.001055e0_r8,.000905e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 7:12, 3) = RESHAPE( (/ &
         .071972e0_r8,.034622e0_r8,.016453e0_r8,.009029e0_r8,.005413e0_r8,.003127e0_r8, &
         .001770e0_r8,.001015e0_r8,.000890e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .080939e0_r8,.032414e0_r8,.014163e0_r8,.007241e0_r8,.004328e0_r8,.002522e0_r8, &
         .001481e0_r8,.000934e0_r8,.000861e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .089906e0_r8,.030206e0_r8,.011873e0_r8,.005453e0_r8,.003244e0_r8,.001918e0_r8, &
         .001192e0_r8,.000852e0_r8,.000832e0_r8,.000787e0_r8,.000694e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .093464e0_r8,.026177e0_r8,.008525e0_r8,.003892e0_r8,.002452e0_r8,.001609e0_r8, &
         .001116e0_r8,.000851e0_r8,.000809e0_r8,.000762e0_r8,.000690e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .097023e0_r8,.022148e0_r8,.005177e0_r8,.002332e0_r8,.001660e0_r8,.001300e0_r8, &
         .001041e0_r8,.000849e0_r8,.000786e0_r8,.000737e0_r8,.000686e0_r8,.000606e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .098107e0_r8,.020617e0_r8,.004758e0_r8,.002137e0_r8,.001516e0_r8,.001211e0_r8, &
         .000999e0_r8,.000848e0_r8,.000778e0_r8,.000730e0_r8,.000677e0_r8,.000603e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 13:18, 3) = RESHAPE( (/ &
         .099191e0_r8,.019085e0_r8,.004340e0_r8,.001943e0_r8,.001371e0_r8,.001122e0_r8, &
         .000957e0_r8,.000847e0_r8,.000770e0_r8,.000724e0_r8,.000669e0_r8,.000601e0_r8, &
         .000557e0_r8,.000412e0_r8,.000362e0_r8,.000326e0_r8,.000309e0_r8,.000302e0_r8, &
         .099547e0_r8,.017725e0_r8,.003693e0_r8,.001578e0_r8,.001125e0_r8,.000985e0_r8, &
         .000879e0_r8,.000795e0_r8,.000712e0_r8,.000643e0_r8,.000584e0_r8,.000521e0_r8, &
         .000482e0_r8,.000384e0_r8,.000351e0_r8,.000322e0_r8,.000302e0_r8,.000302e0_r8, &
         .099903e0_r8,.016365e0_r8,.003045e0_r8,.001214e0_r8,.000879e0_r8,.000848e0_r8, &
         .000801e0_r8,.000744e0_r8,.000654e0_r8,.000562e0_r8,.000499e0_r8,.000441e0_r8, &
         .000410e0_r8,.000358e0_r8,.000342e0_r8,.000322e0_r8,.000302e0_r8,.000302e0_r8, &
         .100218e0_r8,.015035e0_r8,.002537e0_r8,.001037e0_r8,.000790e0_r8,.000726e0_r8, &
         .000673e0_r8,.000628e0_r8,.000579e0_r8,.000512e0_r8,.000440e0_r8,.000374e0_r8, &
         .000307e0_r8,.000253e0_r8,.000227e0_r8,.000208e0_r8,.000194e0_r8,.000186e0_r8, &
         .100534e0_r8,.013704e0_r8,.002028e0_r8,.000861e0_r8,.000701e0_r8,.000604e0_r8, &
         .000546e0_r8,.000513e0_r8,.000504e0_r8,.000462e0_r8,.000381e0_r8,.000307e0_r8, &
         .000201e0_r8,.000151e0_r8,.000117e0_r8,.000098e0_r8,.000090e0_r8,.000093e0_r8, &
         .100892e0_r8,.012873e0_r8,.001886e0_r8,.000785e0_r8,.000643e0_r8,.000568e0_r8, &
         .000519e0_r8,.000487e0_r8,.000471e0_r8,.000437e0_r8,.000368e0_r8,.000305e0_r8, &
         .000201e0_r8,.000151e0_r8,.000117e0_r8,.000098e0_r8,.000090e0_r8,.000093e0_r8/), &
         (/18,6/))
    ozone(1:18, 19:24, 3) = RESHAPE( (/ &
         .102665e0_r8,.010977e0_r8,.001237e0_r8,.000590e0_r8,.000498e0_r8,.000479e0_r8, &
         .000458e0_r8,.000436e0_r8,.000421e0_r8,.000387e0_r8,.000326e0_r8,.000298e0_r8, &
         .000246e0_r8,.000227e0_r8,.000211e0_r8,.000200e0_r8,.000194e0_r8,.000186e0_r8, &
         .104087e0_r8,.010726e0_r8,.000971e0_r8,.000538e0_r8,.000440e0_r8,.000434e0_r8, &
         .000418e0_r8,.000397e0_r8,.000375e0_r8,.000343e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104093e0_r8,.011539e0_r8,.001213e0_r8,.000604e0_r8,.000468e0_r8,.000442e0_r8, &
         .000414e0_r8,.000385e0_r8,.000347e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104106e0_r8,.012997e0_r8,.001479e0_r8,.000639e0_r8,.000468e0_r8,.000422e0_r8, &
         .000392e0_r8,.000372e0_r8,.000342e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .104118e0_r8,.014455e0_r8,.001745e0_r8,.000674e0_r8,.000467e0_r8,.000403e0_r8, &
         .000370e0_r8,.000359e0_r8,.000337e0_r8,.000325e0_r8,.000296e0_r8,.000294e0_r8, &
         .000293e0_r8,.000302e0_r8,.000306e0_r8,.000302e0_r8,.000302e0_r8,.000302e0_r8, &
         .103456e0_r8,.018235e0_r8,.003195e0_r8,.001379e0_r8,.000771e0_r8,.000585e0_r8, &
         .000474e0_r8,.000411e0_r8,.000380e0_r8,.000362e0_r8,.000343e0_r8,.000348e0_r8, &
         .000346e0_r8,.000328e0_r8,.000317e0_r8,.000305e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 25:30, 3) = RESHAPE( (/ &
         .102794e0_r8,.022015e0_r8,.004645e0_r8,.002084e0_r8,.001076e0_r8,.000767e0_r8, &
         .000577e0_r8,.000463e0_r8,.000423e0_r8,.000399e0_r8,.000389e0_r8,.000401e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .101724e0_r8,.026500e0_r8,.007228e0_r8,.003391e0_r8,.002058e0_r8,.001285e0_r8, &
         .000811e0_r8,.000531e0_r8,.000478e0_r8,.000449e0_r8,.000440e0_r8,.000421e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .100654e0_r8,.030986e0_r8,.009812e0_r8,.004698e0_r8,.003041e0_r8,.001803e0_r8, &
         .001044e0_r8,.000599e0_r8,.000533e0_r8,.000499e0_r8,.000491e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .097097e0_r8,.034916e0_r8,.012983e0_r8,.006240e0_r8,.003666e0_r8,.002259e0_r8, &
         .001336e0_r8,.000730e0_r8,.000629e0_r8,.000549e0_r8,.000499e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .093540e0_r8,.038846e0_r8,.016155e0_r8,.007781e0_r8,.004292e0_r8,.002716e0_r8, &
         .001628e0_r8,.000862e0_r8,.000724e0_r8,.000600e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .089618e0_r8,.042869e0_r8,.019963e0_r8,.010502e0_r8,.005966e0_r8,.003525e0_r8, &
         .001936e0_r8,.000906e0_r8,.000769e0_r8,.000625e0_r8,.000508e0_r8,.000452e0_r8, &
         .000451e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8/), &
         (/18,6/))
    ozone(1:18, 31:36, 3) = RESHAPE( (/ &
         .085695e0_r8,.046892e0_r8,.023772e0_r8,.013223e0_r8,.007640e0_r8,.004333e0_r8, &
         .002244e0_r8,.000951e0_r8,.000815e0_r8,.000650e0_r8,.000508e0_r8,.000463e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .082443e0_r8,.047591e0_r8,.025358e0_r8,.014294e0_r8,.008233e0_r8,.004664e0_r8, &
         .002430e0_r8,.001068e0_r8,.000851e0_r8,.000644e0_r8,.000508e0_r8,.000474e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .079190e0_r8,.048290e0_r8,.026945e0_r8,.015366e0_r8,.008826e0_r8,.004995e0_r8, &
         .002616e0_r8,.001184e0_r8,.000887e0_r8,.000637e0_r8,.000508e0_r8,.000486e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .074885e0_r8,.049987e0_r8,.030140e0_r8,.017894e0_r8,.009881e0_r8,.005543e0_r8, &
         .002907e0_r8,.001379e0_r8,.000961e0_r8,.000644e0_r8,.000512e0_r8,.000463e0_r8, &
         .000451e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8, &
         .070579e0_r8,.051684e0_r8,.033335e0_r8,.020423e0_r8,.010935e0_r8,.006091e0_r8, &
         .003197e0_r8,.001573e0_r8,.001034e0_r8,.000650e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .069523e0_r8,.052249e0_r8,.034255e0_r8,.021379e0_r8,.012306e0_r8,.006727e0_r8, &
         .003415e0_r8,.001578e0_r8,.001072e0_r8,.000681e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 37, 3) = (/ &
         .068467e0_r8,.052815e0_r8,.035175e0_r8,.022334e0_r8,.013676e0_r8,.007363e0_r8, &
         .003633e0_r8,.001582e0_r8,.001111e0_r8,.000713e0_r8,.000517e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/)
    !
    !     4. fall
    !
    ozone(1:18, 1:6, 4) = RESHAPE( (/ &
         .062636e0_r8,.039972e0_r8,.022867e0_r8,.011765e0_r8,.005766e0_r8,.003498e0_r8, &
         .002092e0_r8,.001248e0_r8,.000917e0_r8,.000637e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .063734e0_r8,.039842e0_r8,.022004e0_r8,.010859e0_r8,.005712e0_r8,.003589e0_r8, &
         .002155e0_r8,.001174e0_r8,.000856e0_r8,.000612e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .064832e0_r8,.039712e0_r8,.021142e0_r8,.009954e0_r8,.005658e0_r8,.003681e0_r8, &
         .002218e0_r8,.001100e0_r8,.000796e0_r8,.000587e0_r8,.000508e0_r8,.000441e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .067669e0_r8,.037799e0_r8,.019680e0_r8,.009612e0_r8,.005481e0_r8,.003476e0_r8, &
         .002093e0_r8,.001123e0_r8,.000837e0_r8,.000631e0_r8,.000546e0_r8,.000447e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .070506e0_r8,.035887e0_r8,.018218e0_r8,.009270e0_r8,.005304e0_r8,.003271e0_r8, &
         .001969e0_r8,.001145e0_r8,.000878e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .074807e0_r8,.034651e0_r8,.017056e0_r8,.008574e0_r8,.004769e0_r8,.002986e0_r8, &
         .001827e0_r8,.001079e0_r8,.000853e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8/), &
         (/18,6/))
    ozone(1:18, 7:12, 4) = RESHAPE( (/ &
         .079108e0_r8,.033415e0_r8,.015895e0_r8,.007878e0_r8,.004234e0_r8,.002700e0_r8, &
         .001686e0_r8,.001014e0_r8,.000829e0_r8,.000675e0_r8,.000584e0_r8,.000454e0_r8, &
         .000401e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .082715e0_r8,.031096e0_r8,.013350e0_r8,.006131e0_r8,.003205e0_r8,.002043e0_r8, &
         .001304e0_r8,.000842e0_r8,.000734e0_r8,.000631e0_r8,.000555e0_r8,.000494e0_r8, &
         .000480e0_r8,.000408e0_r8,.000385e0_r8,.000361e0_r8,.000351e0_r8,.000349e0_r8, &
         .086323e0_r8,.028778e0_r8,.010805e0_r8,.004383e0_r8,.002175e0_r8,.001386e0_r8, &
         .000923e0_r8,.000671e0_r8,.000640e0_r8,.000587e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .089886e0_r8,.026029e0_r8,.008558e0_r8,.003312e0_r8,.001655e0_r8,.001124e0_r8, &
         .000807e0_r8,.000631e0_r8,.000602e0_r8,.000568e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .093449e0_r8,.023279e0_r8,.006312e0_r8,.002240e0_r8,.001135e0_r8,.000862e0_r8, &
         .000692e0_r8,.000591e0_r8,.000564e0_r8,.000549e0_r8,.000525e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .094163e0_r8,.020198e0_r8,.004594e0_r8,.001772e0_r8,.001077e0_r8,.000806e0_r8, &
         .000649e0_r8,.000565e0_r8,.000547e0_r8,.000537e0_r8,.000521e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8/), &
         (/18,6/))
    ozone(1:18, 13:18, 4) = RESHAPE( (/ &
         .094877e0_r8,.017116e0_r8,.002877e0_r8,.001305e0_r8,.001018e0_r8,.000751e0_r8, &
         .000606e0_r8,.000539e0_r8,.000529e0_r8,.000524e0_r8,.000516e0_r8,.000535e0_r8, &
         .000558e0_r8,.000459e0_r8,.000436e0_r8,.000416e0_r8,.000406e0_r8,.000395e0_r8, &
         .100592e0_r8,.015718e0_r8,.002411e0_r8,.001007e0_r8,.000802e0_r8,.000642e0_r8, &
         .000559e0_r8,.000526e0_r8,.000501e0_r8,.000474e0_r8,.000470e0_r8,.000439e0_r8, &
         .000430e0_r8,.000358e0_r8,.000333e0_r8,.000311e0_r8,.000302e0_r8,.000302e0_r8, &
         .106308e0_r8,.014320e0_r8,.001946e0_r8,.000709e0_r8,.000585e0_r8,.000533e0_r8, &
         .000512e0_r8,.000512e0_r8,.000473e0_r8,.000425e0_r8,.000423e0_r8,.000342e0_r8, &
         .000301e0_r8,.000257e0_r8,.000232e0_r8,.000212e0_r8,.000205e0_r8,.000209e0_r8, &
         .107021e0_r8,.012233e0_r8,.001533e0_r8,.000643e0_r8,.000556e0_r8,.000505e0_r8, &
         .000471e0_r8,.000448e0_r8,.000403e0_r8,.000362e0_r8,.000355e0_r8,.000296e0_r8, &
         .000251e0_r8,.000207e0_r8,.000180e0_r8,.000161e0_r8,.000152e0_r8,.000140e0_r8, &
         .107735e0_r8,.010146e0_r8,.001120e0_r8,.000576e0_r8,.000526e0_r8,.000477e0_r8, &
         .000430e0_r8,.000385e0_r8,.000332e0_r8,.000300e0_r8,.000288e0_r8,.000249e0_r8, &
         .000202e0_r8,.000158e0_r8,.000132e0_r8,.000114e0_r8,.000104e0_r8,.000093e0_r8, &
         .107002e0_r8,.009961e0_r8,.001025e0_r8,.000533e0_r8,.000497e0_r8,.000460e0_r8, &
         .000422e0_r8,.000385e0_r8,.000332e0_r8,.000300e0_r8,.000288e0_r8,.000249e0_r8, &
         .000202e0_r8,.000158e0_r8,.000132e0_r8,.000114e0_r8,.000104e0_r8,.000093e0_r8/), &
         (/18,6/))
    ozone(1:18, 19:24, 4) = RESHAPE( (/ &
         .109850e0_r8,.010424e0_r8,.001079e0_r8,.000567e0_r8,.000498e0_r8,.000479e0_r8, &
         .000463e0_r8,.000448e0_r8,.000418e0_r8,.000399e0_r8,.000389e0_r8,.000367e0_r8, &
         .000351e0_r8,.000328e0_r8,.000320e0_r8,.000337e0_r8,.000355e0_r8,.000304e0_r8, &
         .112375e0_r8,.011255e0_r8,.001357e0_r8,.000744e0_r8,.000585e0_r8,.000533e0_r8, &
         .000512e0_r8,.000512e0_r8,.000504e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .111320e0_r8,.011439e0_r8,.001485e0_r8,.000843e0_r8,.000642e0_r8,.000549e0_r8, &
         .000512e0_r8,.000512e0_r8,.000504e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .107719e0_r8,.013052e0_r8,.001822e0_r8,.000953e0_r8,.000701e0_r8,.000604e0_r8, &
         .000551e0_r8,.000525e0_r8,.000509e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .104118e0_r8,.014665e0_r8,.002159e0_r8,.001063e0_r8,.000759e0_r8,.000659e0_r8, &
         .000589e0_r8,.000539e0_r8,.000514e0_r8,.000499e0_r8,.000491e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .104145e0_r8,.018082e0_r8,.003274e0_r8,.001493e0_r8,.000919e0_r8,.000762e0_r8, &
         .000678e0_r8,.000641e0_r8,.000584e0_r8,.000531e0_r8,.000495e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8/), &
         (/18,6/))
    ozone(1:18, 25:30, 4) = RESHAPE( (/ &
         .104172e0_r8,.021499e0_r8,.004389e0_r8,.001922e0_r8,.001078e0_r8,.000865e0_r8, &
         .000767e0_r8,.000743e0_r8,.000654e0_r8,.000562e0_r8,.000499e0_r8,.000486e0_r8, &
         .000501e0_r8,.000502e0_r8,.000509e0_r8,.000561e0_r8,.000607e0_r8,.000515e0_r8, &
         .101718e0_r8,.026523e0_r8,.008473e0_r8,.004382e0_r8,.002392e0_r8,.001505e0_r8, &
         .001036e0_r8,.000836e0_r8,.000727e0_r8,.000618e0_r8,.000550e0_r8,.000494e0_r8, &
         .000501e0_r8,.000479e0_r8,.000473e0_r8,.000509e0_r8,.000541e0_r8,.000445e0_r8, &
         .099264e0_r8,.031548e0_r8,.012557e0_r8,.006841e0_r8,.003705e0_r8,.002144e0_r8, &
         .001304e0_r8,.000930e0_r8,.000799e0_r8,.000675e0_r8,.000601e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .097501e0_r8,.035382e0_r8,.014856e0_r8,.008207e0_r8,.004619e0_r8,.002720e0_r8, &
         .001610e0_r8,.001012e0_r8,.000829e0_r8,.000687e0_r8,.000610e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .095737e0_r8,.039215e0_r8,.017155e0_r8,.009572e0_r8,.005532e0_r8,.003296e0_r8, &
         .001916e0_r8,.001094e0_r8,.000858e0_r8,.000699e0_r8,.000618e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .092863e0_r8,.042419e0_r8,.020704e0_r8,.012034e0_r8,.007417e0_r8,.004504e0_r8, &
         .002590e0_r8,.001334e0_r8,.000977e0_r8,.000731e0_r8,.000622e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8/), &
         (/18,6/))
    ozone(1:18, 31:36, 4) = RESHAPE( (/ &
         .089988e0_r8,.045622e0_r8,.024253e0_r8,.014495e0_r8,.009303e0_r8,.005711e0_r8, &
         .003264e0_r8,.001574e0_r8,.001096e0_r8,.000762e0_r8,.000627e0_r8,.000503e0_r8, &
         .000501e0_r8,.000459e0_r8,.000436e0_r8,.000460e0_r8,.000486e0_r8,.000398e0_r8, &
         .084591e0_r8,.046553e0_r8,.025037e0_r8,.015156e0_r8,.009841e0_r8,.006124e0_r8, &
         .003534e0_r8,.001693e0_r8,.001170e0_r8,.000793e0_r8,.000631e0_r8,.000537e0_r8, &
         .000551e0_r8,.000509e0_r8,.000486e0_r8,.000516e0_r8,.000548e0_r8,.000446e0_r8, &
         .079194e0_r8,.047483e0_r8,.025822e0_r8,.015818e0_r8,.010380e0_r8,.006537e0_r8, &
         .003804e0_r8,.001811e0_r8,.001245e0_r8,.000825e0_r8,.000635e0_r8,.000570e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .077409e0_r8,.048159e0_r8,.026661e0_r8,.016596e0_r8,.010962e0_r8,.006972e0_r8, &
         .004160e0_r8,.002132e0_r8,.001391e0_r8,.000868e0_r8,.000686e0_r8,.000601e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .075625e0_r8,.048835e0_r8,.027500e0_r8,.017375e0_r8,.011544e0_r8,.007407e0_r8, &
         .004516e0_r8,.002453e0_r8,.001536e0_r8,.000912e0_r8,.000737e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8, &
         .074927e0_r8,.049459e0_r8,.029215e0_r8,.018025e0_r8,.011754e0_r8,.007786e0_r8, &
         .004972e0_r8,.002926e0_r8,.001817e0_r8,.001025e0_r8,.000758e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8/), &
         (/18,6/))
    ozone(1:18, 37, 4) = (/ &
         .074229e0_r8,.050084e0_r8,.030930e0_r8,.018676e0_r8,.011965e0_r8,.008165e0_r8, &
         .005428e0_r8,.003399e0_r8,.002098e0_r8,.001138e0_r8,.000780e0_r8,.000632e0_r8, &
         .000603e0_r8,.000559e0_r8,.000538e0_r8,.000574e0_r8,.000614e0_r8,.000515e0_r8/)

    ozsig(:) = (/ &
         .020747_r8,.073986_r8,.124402_r8,.174576_r8,.224668_r8,.274735_r8, &
         .324767_r8,.374806_r8,.424818_r8,.497450_r8,.593540_r8,.688125_r8, &
         .777224_r8,.856317_r8,.920400_r8,.960480_r8,.981488_r8,.995004_r8/)


    first_getoz=.TRUE.

    IF(first_getoz)THEN
       mon_getoz=yrl/12.0_r8
       year_getoz=yrl
       IF(nlm_getoz.NE.kmax)THEN
          inter_getoz=.TRUE.
       ELSE
          inter_getoz=.FALSE.
          DO l=1,nlm_getoz
             ll=nlm_getoz-l+1
             IF(ABS(ozsig(l)-sl(ll)).GT.0.0001_r8) inter_getoz=.TRUE.
          END DO
       ENDIF
       first_getoz=.FALSE.
    ENDIF
  END SUBROUTINE InitGetoz

  SUBROUTINE getoz (ncols,adjNCols, kmax,sigmid,colrad,date,o3l)
    !
    ! input parameters and variables:
    !     ncols  =  number of atmospheric columns
    !     kmax   =  number of atmospheric layers
    !     sigmid =  sigma coordinate at middle of layer
    !     colrad =  colatitude of each column (0-3.14 from np to sp in radians)
    !     date   =  model julian date
    !
    ! tabulated data
    !     ozone  =  climatological ozone mixing ratio in 18 sigma layers 
    !               and in 5 degree latitude interval
    !
    ! output variables:
    !     o3l   =  18 layers ozone mixing ratio in given lat and date
    !
    !==========================================================================
    ! :: kmax.....Number of grid points at vertical
    ! :: sigmid.......sigma coordinate at middle of layer
    ! :: pai......constant pi=3.1415926
    ! :: yrl......length of year in days
    !==========================================================================
    
    ! Input variables
    INTEGER      ,    INTENT(IN   ) :: ncols,adjNCols
    INTEGER      ,    INTENT(IN   ) :: kmax
    REAL(KIND=r8),    INTENT(IN   ) :: sigmid (kmax)
    REAL(KIND=r8),    INTENT(IN   ) :: colrad(ncols)
    REAL(KIND=r8),    INTENT(INOUT) :: date
    REAL(KIND=r8),PARAMETER :: pai=3.1415926e00_r8

    ! Output variable
    REAL(KIND=r8),    INTENT(OUT) :: o3l(ncols,kmax)

    ! Local variables
    REAL(KIND=r8) :: a1   (nlm_getoz)
    REAL(KIND=r8) :: a2   (nlm_getoz)
    REAL(KIND=r8) :: a3   (nlm_getoz)
    REAL(KIND=r8) :: a4   (nlm_getoz)
    REAL(KIND=r8) :: b1   (nlm_getoz)
    REAL(KIND=r8) :: b2   (nlm_getoz)
    REAL(KIND=r8) :: b3   (nlm_getoz)
    REAL(KIND=r8) :: b4   (nlm_getoz)
    REAL(KIND=r8) :: do3a (nlm_getoz)
    REAL(KIND=r8) :: do3b (nlm_getoz)
    REAL(KIND=r8) :: ozo3l(ncols,nlm_getoz)

    REAL(KIND=r8), PARAMETER :: rlag = 14.8125e0_r8

    INTEGER :: l
    INTEGER :: la   (ncols)
    INTEGER :: ll   (ncols)
    INTEGER :: kmx
    INTEGER :: imon
    INTEGER :: isea
    INTEGER :: k
    INTEGER :: i
    INTEGER :: kk    (ncols,kmax)
    REAL(KIND=r8)    :: theta
    REAL(KIND=r8)    :: flat
    REAL(KIND=r8)    :: rang
    REAL(KIND=r8)    :: rsin1(ncols)
    REAL(KIND=r8)    :: rcos1(ncols)
    REAL(KIND=r8)    :: rcos2(ncols)
    REAL(KIND=r8)    :: rate (ncols)
    REAL(KIND=r8)    :: aa
    REAL(KIND=r8)    :: bb
    LOGICAL :: notfound(kMax)

    kmx=nlm_getoz
    !
    !     find closest place in the data according to input slat.
    !
    IF(date.GT.year_getoz) date=date-year_getoz

    imon=date/mon_getoz + 1

    IF(imon.LT.1)imon=1

    isea=imon/3 + 1

    IF(isea.EQ.5) isea=1
    IF(isea.GT.5) THEN
       WRITE(nfprt,"('0 ERROR IN ISEA - TERMINATION IN SUBROUTINE GETOZ')")
       WRITE(nferr,"('0 ERROR IN ISEA - TERMINATION IN SUBROUTINE GETOZ')")
       STOP 9954
    END IF
    DO i=1,adjNCols
       theta = 90.0_r8-(180.0_r8/pai)*colrad(i) ! colatitude -> latitude
       ! the 180 degrees are divided into 37 bands with 5deg each
       ! except for the first and last, which have 2.5 deg
       ! The centers of the bands are located at:
       !   90, 85, 80, ..., 5, 0, -5, ..., -85, -90 (37 latitudes)
       flat  = 0.2_r8*theta ! indexing the latitudes: goes from -18. to +18.
       ! find the latitude index before and after each latitude
       la(i)    = 19.501e0_r8-flat !
       ll(i)    = 19.001e0_r8-flat

       !
       !     find sin and cos coefficients for time interpolation.
       !
       rang=2.0e0_r8*pai*(date-rlag)/year_getoz
       rsin1(i)=SIN(rang)
       rcos1(i)=COS(rang)
       rcos2(i)=COS(2.0e0_r8*rang)
       rate(i)=REAL(19-ll(i),r8)-flat
       !
       !     ozone interpolation in latitude and time
       !
    END DO
    DO k=1,kmx
       DO i=1,adjNCols
          a1(k) =2.5e-1_r8*(ozone(k,la(i),1)+ozone(k,la(i),2)+ &
               ozone(k,la(i),3)+ozone(k,la(i),4))
          a2(k) =0.5e0_r8*(ozone(k,la(i),2)-ozone(k,la(i),4))
          a3(k) =0.5e0_r8*(ozone(k,la(i),1)-ozone(k,la(i),3))
          a4(k) =2.5e-1_r8*(ozone(k,la(i),1)+ozone(k,la(i),3)- &
               ozone(k,la(i),2)-ozone(k,la(i),4))
          b1(k) =2.5e-1_r8*(ozone(k,ll(i),1)+ozone(k,ll(i),2)+ &
               ozone(k,ll(i),3)+ozone(k,ll(i),4))
          b2(k) =0.5e0_r8*(ozone(k,ll(i),2)-ozone(k,ll(i),4))
          b3(k) =0.5e0_r8*(ozone(k,ll(i),1)-ozone(k,ll(i),3))
          b4(k) =2.5e-1_r8*(ozone(k,ll(i),1)+ozone(k,ll(i),3)- &
               ozone(k,ll(i),2)-ozone(k,ll(i),4))
          do3a(k)=a1(k)+rsin1(i)*a2(k)+rcos1(i)*a3(k)+rcos2(i)*a4(k)
          do3b(k)=b1(k)+rsin1(i)*b2(k)+rcos1(i)*b3(k)+rcos2(i)*b4(k)
          ozo3l(i,k)=do3a(k)+rate(i)*(do3b(k)-do3a(k))
          ozo3l(i,k)=1.0e-04_r8*ozo3l(i,k)
       END DO
    END DO
    !print *,'LFR->Oz: ', (ozo3l(1,k),k=1,kmx)
    IF(inter_getoz)THEN
       DO l=1,kmax
!          print *,'LFR->MQB: ',sigmid(l),ozsig(1),sigmid(l) > ozsig(1)
          notfound(l) = sigmid(l) > ozsig(1)
          IF (notfound(l)) THEN
             kk(1,l)=kmx
          ELSE
             kk(1,l)=2
          END IF
!	  print *,'lfr->1.kk(1,l) :',l,kk(1,l)
       END DO
       DO l=1,kmax
          IF (notfound(l)) THEN
             DO k=2,kmx
                IF(sigmid(l).GT.ozsig(k-1).AND.sigmid(l).LE.ozsig(k))THEN
                   kk(1,l)=k
                   EXIT
                END IF
             END DO
          END IF
!	  print *,'lfr->2.kk(1,l) :',l,kk(1,l)
	  
       END DO
       DO l = 1, kmax
          DO i= 2, adjNCols
             kk(i,l) = kk(1,l)
          END DO
!	  print *,'lfr->3.kk(1,l) :',l,kk(1,l)
       END DO
    END IF
!    print *,'LFR->inter_getoz=',inter_getoz
    IF(inter_getoz)THEN
       DO l=1,kmax
          DO i=1,adjNCols
             aa=(ozo3l(i,kk(i,l))-ozo3l(i,kk(i,l)-1))/(ozsig(kk(i,l))-ozsig(kk(i,l)-1))
             bb=ozo3l(i,kk(i,l)-1)-aa*ozsig(kk(i,l)-1)
             o3l(i,kmax+1-l)=bb+aa*sigmid(l)
!	     print *,'LFR->',l,i,kmax+1-l,aa,bb,o3l(i,kmax+1-l)
          END DO
       END DO
    END IF
    IF(.NOT.inter_getoz)THEN
       DO l=1,nlm_getoz
          DO i=1,adjNCols
             o3l(i,l)=ozo3l(i,l)
          END DO
       END DO
    ENDIF
  END SUBROUTINE getoz

  SUBROUTINE getco2(time,co2val)
    !==========================================================================
    ! getco2: Interpolates Mauna Loa data for a given time
    !
    ! *** Atmospheric CO2 concentrations (ppmv) derived from in situ  ***
    ! *** air samples collected at Mauna Loa Observatory, Hawaii      ***
    !
    ! Data:
    !
    !   http://cdiac.ornl.gov/trends/co2/contents.htm
    !   http://cdiac.ornl.gov/ftp/trends/co2/maunaloa.co2
    !
    ! Parabolic fitting by hbarbosa@cptec.inpe.br, 17 Jan 2007:
    !
    !   co2val = a*(time-2000)^2 + b*(time-2000) + c
    !
    !       a  = 0.0116696   +/- 0.0005706    (4.89%)
    !       b  = 1.79984     +/- 0.022        (1.222%)
    !       c  = 369         +/- 0.1794       (0.04863%)
    !
    !==========================================================================
    !     time.......date of current data
    !     time(1)....hour(00/12)
    !     time(2)....month
    !     time(3)....day of month
    !     time(4)....year
    !
    !    co2val....co2val is wgne standard value in ppm "co2val = /345.0/
    !==========================================================================

    IMPLICIT NONE
    INTEGER,PARAMETER :: r8 = SELECTED_REAL_KIND(15)
    REAL(KIND=r8), PARAMETER :: A = 0.0116696
    REAL(KIND=r8), PARAMETER :: B = 1.79984
    REAL(KIND=r8), PARAMETER :: C = 369.0

    INTEGER,       INTENT(IN ) :: time(4)
    REAL(KIND=r8), INTENT(OUT) :: co2val

    REAL(KIND=r8) :: TDIF

    tdif=time(4) + (time(2)-1.)/12. + (time(3)-1.+ time(1)/24.)/365. - 2000.

    co2val = A*tdif**2 + B*tdif + C

    !    WRITE(*,123) time,tdif+2000.,co2val
    !123 format('hmjb co2val date=',3(I2,1x),I4,' fyear=',F10.5,' val=',F7.3)

    RETURN
  END SUBROUTINE getco2

 SUBROUTINE CreateFlip(k,ngrid)
     INTEGER, INTENT(in) :: k,ngrid
     INTEGER :: i,j          
          
     j=k
     DO i=1,k
        mem_UKMO(ngrid)%flip(i)=j
	j=j-1
     END DO

  END SUBROUTINE Createflip


  SUBROUTINE Cloud_Micro_CCM3(&
       ! Model info
       ncols,kmax , sigmid, sigbot, delsig, imask   , &
       ! Atmospheric Fields
       Ps   , Te   , Qe    , tsea  , FlipPbot ,pptop,        &
       ! Cloud properties
       clwp , lmixr, fice  , rei   , rel   , taud       )
    IMPLICIT NONE

    ! As in the CCM2, cloud optical properties in the CCM3 are accounted for using
    ! the Slingo (1989) parameterization for liquid water droplet clouds. This
    ! scheme relates the extinction optical depth, the single-scattering albedo,
    ! and the asymmetry parameter to the cloud liquid water path and cloud drop
    ! effective radius. The latter two microphysical cloud properties were
    ! statically specified in the CCM2. In particular, in-cloud liquid water paths
    ! were evaluated from a prescribed, meridionally and height varying, but
    ! time independent, cloud liquid water density profile, rho_l(z), which
    ! was analytically determined on the basis of a meridionally specified
    ! liquid water scale height (e.g. see Kiehl et al., 1994; Kiehl, 1991).
    ! The cloud drop effective radius was simplly specified to be 10microns
    ! for all clouds. The CCM3 continues to diagnose cloud optical properties,
    ! but relaxes the rigid CCM2 framework. CCM3 employs the same exponentially
    ! decaying vertical profile for in-cloud water concentration
    !
    !             rho_l=rho_l^0*exp(-z/h_l)               eq 4.a.11
    !
    ! , where rho_l^0=0.21g/m3. Instead of specifying a zonally symmetric meridional
    ! dependence for the cloud water scale heigh, h_l, it is locally diagnosed
    ! as a function of the vertically integrated water vapor (precipitable water) 
    !
    !          h_l=700 ln [1+\frac{1}{g} \int_pT^ps q dp]  eq 4.a.12
    !
    ! hmjb> It is not explained, but the units of h_l must be meters, the same 
    ! hmjb> of the height, z.
    !
    ! The cloud water path (CWP) is determined by integrating the liquid
    ! water concentration using
    !
    !                 cwp = int rho_l dz     eq. 4.a.13
    ! 
    ! Which can be analytically evaluated for an arbitrary layer k as
    !
    !  rho_l^0 h_l [exp(-z_bot(k)/h_l) - exp(-z_top(k)/h_l)]   eq. 4.a.14
    !
    ! Where z_bot and z_top are the heights of the k'th layer interfaces.
    !
    ! hmjb> It is not explained, but the units of clwp must be g/m2
    ! hmjb> since it is the integral of rho_l*dz (eq.4.a.13)
    !
    ! CCM3 Documentation, pg 50
    ! Observational studies have shown a distinct difference between
    ! maritime and continental effective cloud drop size, r_e, for warm
    ! clouds. For this reason, the CCM3 differentiates between the cloud
    ! drop effective radius for clouds diagnosed over maritime and
    ! continental regimes (Kiehl, 1994). Over the ocean, the cloud drop
    ! effective radius for liquid water clouds, r_el, is specified to be
    ! 10microns, as in the CCM3. Over land masses r_el is determinedusing
    !
    ! r_el = 5 microns             T > -10oC
    !      = 5-5(t+10)/20 microns  -30oC <= T <= -10oC     eq. 4.a.14.1
    !      = r_ei                  T < -30oC
    !
    ! An ice particle effective radius, r_ei, is also diagnosed by CCM3,
    ! which at the moment amounts to a specification of ice radius as a
    ! function of normalized pressure
    !
    ! r_ei = 10 microns                                 p/ps > p_I^high   
    !      = r_ei^max - (r_ei^max - r_ei^min)           p/ps <= p_I^high       eq. 4.a.15.1
    !            *[(p/ps)-p_I^high/(p_I^high-p_I^low)]
    ! where r_ei^max=30microns, r_ei^min=10microns, p_I^high=0.4 and p_I^low=0.0
    !
    ! hmjb>> I think there is a typo in the equation, otherwise the 
    ! hmjb>> expression for r_ei is not a continuous funcion of p/ps.
    ! hmjb>> For p/ps=p_I^high, r_ei should be r_ei^min and not r_ei^max.
    ! hmjb>> The correct equation is:
    ! hmjb>> r_ei = 10 microns                                 p/ps > p_I^high   
    ! hmjb>>      = R_EI^MIN - (r_ei^max - r_ei^min)           p/ps <= p_I^high
    ! hmjb>>            *[(p/ps)-p_I^high/(p_I^high-p_I^low)]
    !
    !--------------------------------------------------------------------------------- 
    ! Input/Output Variables
    !--------------------------------------------------------------------------------- 
     INTEGER,PARAMETER :: r8 = SELECTED_REAL_KIND(15)
     INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
     REAL (KIND=r8), PARAMETER   :: gasr  =  287.05_r8! gas constant of dry air        (j/kg/k)
     REAL (KIND=r8), PARAMETER   :: grav  =   9.8e0_r8! gravity constant               (m/s**2)

    ! Model info
    INTEGER         , INTENT(IN) :: ncols 
    INTEGER         , INTENT(IN) :: kmax   
    REAL(KIND=r8)   , INTENT(IN) :: sigmid(kmax)  ! Sigma cordinate at middle of layer
    REAL(KIND=r8)   , INTENT(IN) :: sigbot(kmax)  ! Sigma cordinate at bottom of layer
    REAL(KIND=r8)   , INTENT(IN) :: delsig(kmax)  ! Layer thickness (sigma)
    INTEGER(KIND=i8), INTENT(IN) :: imask (ncols) ! Ocean/Land mask

    ! Atmospheric Fields
    REAL(KIND=r8), INTENT(IN) :: Ps  (ncols)      ! Surface pressure (mb)
    REAL(KIND=r8), INTENT(IN) :: Te  (ncols,kmax) ! Temperature (K)
    REAL(KIND=r8), INTENT(IN) :: Qe  (ncols,kmax) ! Specific Humidity (g/g)
    REAL(KIND=r8), INTENT(IN) :: tsea(ncols)
    REAL(KIND=r8), INTENT(IN) :: FlipPbot (ncols,kmax)  ! Pressure at bottom of layer (mb)
    REAL(KIND=r8), INTENT(IN) :: pptop  ! Model-top presure    

    ! Cloud properties
    REAL(KIND=r8), INTENT(OUT) :: clwp (ncols,kmax) ! Cloud Liquid Water Path
    REAL(KIND=r8), INTENT(OUT) :: lmixr(ncols,kmax) ! Ice/Water mixing ratio
    REAL(KIND=r8), INTENT(OUT) :: fice (ncols,kmax) ! Fractional amount of cloud that is ice
    REAL(KIND=r8), INTENT(OUT) :: rei  (ncols,kmax) ! Ice particle Effective Radius (microns)
    REAL(KIND=r8), INTENT(OUT) :: rel  (ncols,kmax) ! Liquid particle Effective Radius (microns)
    REAL(KIND=r8), INTENT(OUT) :: taud (ncols,kmax) ! Shortwave cloud optical depth

    !--------------------------------------------------------------------------------- 
    ! Parameters
    !--------------------------------------------------------------------------------- 

    REAL(KIND=r8), PARAMETER :: abarl=2.261e-2_r8
    REAL(KIND=r8), PARAMETER :: bbarl=1.4365_r8
    REAL(KIND=r8), PARAMETER :: abari=3.448e-3_r8
    REAL(KIND=r8), PARAMETER :: bbari=2.431_r8

    REAL(KIND=r8), PARAMETER :: clwc0   = 0.21_r8 ! Reference liquid water concentration (g/m3)        
    REAL(KIND=r8), PARAMETER :: reimin  = 10.0_r8 ! Minimum of Ice particle efective radius (microns)  
    REAL(KIND=r8), PARAMETER :: reirnge = 20.0_r8 ! Range of Ice particle efective radius (microns)   
    REAL(KIND=r8), PARAMETER :: sigrnge = 0.4_r8  ! Normalized pressure range                         
    REAL(KIND=r8), PARAMETER :: sigmax  = 0.4_r8  ! Normalized pressure maximum                       

!    REAL(KIND=r8), PARAMETER :: pptop = 0.5       ! Model-top presure                                 

    !--------------------------------------------------------------------------------- 
    ! Local Variables
    !--------------------------------------------------------------------------------- 

    REAL(KIND=r8) :: hl     (ncols)        ! cloud water scale heigh (m)
    REAL(KIND=r8) :: rhl    (ncols)        ! cloud water scale heigh (m)
    REAL(KIND=r8) :: pw     (ncols)        ! precipitable water (kg/m2)
    REAL(KIND=r8) :: Zibot  (ncols,kmax+1) ! Height at middle of layer (m)
    REAL(KIND=r8) :: emziohl(ncols,kmax+1) ! exponential of Minus zi Over hl (no dim)

    REAL(KIND=r8) :: tauxcl(ncols,kmax)    ! extinction optical depth of liquid phase
    REAL(KIND=r8) :: tauxci(ncols,kmax)    ! extinction optical depth of ice phase

    !-- Aux variables

    INTEGER :: i,k
    REAL(KIND=r8) :: weight

    !--------------------------------------------------------------------------------- 
    !--------------------------------------------------------------------------------- 

    clwp=0.0_r8
    lmixr=0.0_r8
    fice=0.0_r8
    rei=0.0_r8
    rel=0.0_r8
    pw=0.0_r8

    ! Heights corresponding to sigma at middle of layer: sig(k)
    ! Assuming isothermal atmosphere within each layer
    DO i=1,nCols
       Zibot(i,1) = 0.0_r8
       DO k=2,kMax
          Zibot(i,k) = Zibot(i,k-1) + (gasr/grav)*Te(i,k-1)* &
!               LOG(sigbot(k-1)/sigbot(k))
               LOG(FlipPbot(i,kmax+2-k)/FlipPbot(i,kmax+1-k))
          
       END DO
    END DO

    DO i=1,nCols
       Zibot(i,kmax+1)=Zibot(i,kmax)+(gasr/grav)*Te(i,kmax)* &
            LOG(FlipPbot(i,1)/pptop)
    END DO

! precitable water, pw = sum_k { delsig(k) . Qe(k) } . Ps . 100 / g
!                   pw = sum_k { Dp(k) . Qe(k) } / g
!
! 100 is to change from mbar to pascal
! Dp(k) is the difference of pressure (N/m2) between bottom and top of layer
! Qe(k) is specific humidity in (g/g)
! gravity is m/s2 => so pw is in Kg/m2
    DO k=1,kmax
       DO i = 1,nCols
          pw(i) = pw(i) + delsig(k)*Qe(i,k)
       END DO
    END DO
    DO i = 1,nCols
       pw(i)=100.0_r8*pw(i)*Ps(i)/grav !LFR
!       pw(i)=pw(i)*Ps(i)/grav
    END DO
    !
    ! diagnose liquid water scale height from precipitable water
    DO i=1,nCols
       hl(i)  = 700.0_r8*LOG(MAX(pw(i)+1.0_r8,1.0_r8))
       rhl(i) = 1.0_r8/hl(i)
    END DO
    !hmjb> emziohl stands for Exponential of Minus ZI Over HL
    DO k=1,kmax+1
       DO i=1,nCols
!          emziohl(i,k) = EXP(-zibot(i,k)/hl(i))
          emziohl(i,k) = EXP(-zibot(i,k)*rhl(i))
       END DO
    END DO
!    DO i=1,ncols
!       emziohl(i,kmax+1) = 0.0_r8
!    END DO

    ! The units are g/m2.
    DO k=1,kmax
       DO i=1,nCols
          clwp(i,k) = clwc0*hl(i)*(emziohl(i,k) - emziohl(i,k+1))
       END DO
    END DO

! If we want to calculate the 'droplets/cristals' mixing ratio, we need
! to find the amount of dry air in each layer. 
!
!             dry_air_path = int rho_air dz  
!
! This can be simply done using the hydrostatic equation:
!
!              dp/dz = -rho grav
!            dp/grav = -rho dz
!
!
! The units are g/m2. The factor 1e5 accounts for the change
!  mbar to Pa and kg/m2 to g/m2. 
    DO k=1,kmax
       DO i=1,nCols
          lmixr(i,k)=clwp(i,k)*grav*1e-5/delsig(k)/Ps(i)
       END DO
    END DO

    ! determine Ice particle Effective Radius (rei)
    ! as function of normalized pressure 
    ! docs CCM3, eq 4.a.15.1
    DO k=1,kmax
       weight   = MIN((sigmid(k)-sigmax)/sigrnge,0.0_r8)
       DO i=1,nCols
          rei(i,k) = reimin - reirnge*weight
       END DO
    END DO

    ! define fractional amount of cloud that is ice
    ! if warmer than -10 degrees c then water phase
    ! docs CCM3, eq 4.a.16.1     
    fice=MAX(MIN((263.16_r8-Te)*0.05_r8,1.0_r8),0.0_r8)

    ! determine Liquid particle Effective Radius (rel) 
    ! as function of normalized pressure
    ! docs CCM3, eq 4.a.15.1

    DO i=1,nCols
       IF (imask(i) .lt. 1) THEN
          rel(i,:) = 10.0_r8
       ELSE
          rel(i,:) = 5.0_r8+5.0_r8*fice(i,:)
          WHERE(fice(i,:).eq.1) rel(i,:) = rei(i,:)
       END IF
    END DO    

    ! Compute optical depth from liquid water
    DO k=1,kmax
       DO i=1,nCols
          ! ccm3 manual, page 53, eqs 4.b.3 and 4.b.7
          tauxcl(i,k) = clwp(i,k)*(abarl + bbarl/rel(i,k))*(1.0_r8-fice(i,k))
          tauxci(i,k) = clwp(i,k)*(abari + bbari/rei(i,k))*fice(i,k)
          IF (tsea(i) > 0.0_r8) THEN
             taud(i,k)=0.70_r8*(tauxcl(i,k)+tauxci(i,k))
          ELSE
             taud(i,k)=1.00_r8*(tauxcl(i,k)+tauxci(i,k))
          END IF
       END DO
    END DO


  END SUBROUTINE Cloud_Micro_CCM3

END MODULE UkmoAdapt
