MODULE aer1_list
  IMPLICIT NONE
  
  
  INTEGER,PARAMETER :: maxnspecies= 200
  INTEGER,PARAMETER :: nspecies=5
  INTEGER,PARAMETER :: nmodes=18
  
  
  !Name of species 
  CHARACTER(LEN=8),PARAMETER,DIMENSION(nspecies) :: spc_name=(/ &
      'sulf   ' & !
     ,'bcar   ' & !
     ,'ocar   ' & !
     ,'dust   ' & !
     ,'seas   ' & !
   /)
  
  
  !Number of each specie   
  INTEGER,PARAMETER :: sulf=001
  INTEGER,PARAMETER :: bcar=002
  INTEGER,PARAMETER :: ocar=003
  INTEGER,PARAMETER :: dust=004
  INTEGER,PARAMETER :: seas=005
  
  
  !for memory allocattion: 
  INTEGER,PARAMETER :: on = 1
  INTEGER,PARAMETER :: off = 0
  
  
  !Name of species 
  CHARACTER(LEN=8),PARAMETER,DIMENSION(nmodes) :: mode_name=(/ &
      'akk' & !
     ,'acc' & !
     ,'dd1' & !
     ,'ds1' & !
     ,'dd2' & !
     ,'ds2' & !
     ,'ssa' & !
     ,'ssc' & !
     ,'sss' & !
     ,'occ' &
      'bc1' & !
     ,'bc2' & !
     ,'bc3' & !
     ,'ocs' & !
     ,'dbc' & !
     ,'boc' & !
     ,'bcs' & !
     ,'mxx' & !
           /)
 
  INTEGER,PARAMETER :: nucle = 1 ! nucleation mode
  INTEGER,PARAMETER :: accum = 2 ! accumulation mode
  INTEGER,PARAMETER :: coarse = 3 ! coarse mode 
  INTEGER, PARAMETER :: X = 0
  
  !define if a specific mode will exist (=1), not (=0), No permission by mechanism (=X)
  !for  bins 1 to 10  
  !Mechanism = 8
  INTEGER,PARAMETER,DIMENSION(nmodes,nspecies) :: mode_alloc=RESHAPE((/ &
!----------------------------------------------------------
!        akk acc dd1 ds1 dd2 ds2 ssa ssc sss occ bc1 bc2 bc3 ocs dbc boc bcs mxx
!         0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1
!-mode    1   2   3   4   5   6   7   8   9   0   1   2   3   4   5   6   7   8
!----------------------------------------------------------
          X , 1 , 1 , 1 , X , X , X , X , 1 , 1 , 1 , 1 , X , X , X , X , X , 1 , & ! sulf
          X , 1 , 1 , 0 , X , X , X , X , 0 , 0 , 1 , 1 , X , X , X , X , X , 1 , & ! bcar
          X , 1 , 1 , 0 , X , X , X , X , 0 , 1 , 0 , 0 , X , X , X , X , X , 1 , & ! ocar
          X , 0 , 1 , 1 , X , X , X , X , 0 , 0 , 0 , 0 , X , X , X , X , X , 1 , & ! dust
          X , 0 , 0 , 0 , X , X , X , X , 1 , 0 , 0 , 0 , X , X , X , X , X , 1   & ! seas
!----------------------------------------------------------
          /),(/nmodes,nspecies/))
    
  CHARACTER(LEN=8),PARAMETER,DIMENSION(nmodes,nspecies) :: aer_name=RESHAPE((/ &
      'sulf_akk','sulf_acc','sulf_dd1','sulf_ds1','sulf_dd2','sulf_ds2','sulf_ssa','sulf_ssc','sulf_sss','sulf_occ','sulf_bc1','sulf_bc2','sulf_bc3','sulf_ocs','sulf_dbc','sulf_boc','sulf_bcs','sulf_mxx',
     ,'bcar_akk','bcar_acc','bcar_dd1','bcar_ds1','bcar_dd2','bcar_ds2','bcar_ssa','bcar_ssc','bcar_sss','bcar_occ','bcar_bc1','bcar_bc2','bcar_bc3','bcar_ocs','bcar_dbc','bcar_boc','bcar_bcs','bcar_mxx',
     ,'ocar_akk','ocar_acc','ocar_dd1','ocar_ds1','ocar_dd2','ocar_ds2','ocar_ssa','ocar_ssc','ocar_sss','ocar_occ','ocar_bc1','ocar_bc2','ocar_bc3','ocar_ocs','ocar_dbc','ocar_boc','ocar_bcs','ocar_mxx',
     ,'dust_akk','dust_acc','dust_dd1','dust_ds1','dust_dd2','dust_ds2','dust_ssa','dust_ssc','dust_sss','dust_occ','dust_bc1','dust_bc2','dust_bc3','dust_ocs','dust_dbc','dust_boc','dust_bcs','dust_mxx',
     ,'seas_akk','seas_acc','seas_dd1','seas_ds1','seas_dd2','seas_ds2','seas_ssa','seas_ssc','seas_sss','seas_occ','seas_bc1','seas_bc2','seas_bc3','seas_ocs','seas_dbc','seas_boc','seas_bcs','seas_mxx',
    /),(/nmodes,nspecies/))

  real :: mass_bin_dist(nmodes) ! only for ash

  
  !section for aer type 1: dus
  !This parameters are use for documentation only. 
  !Use them in a program in substitution of numerical terms.
  INTEGER,PARAMETER :: src     = 1 ! source term 
  INTEGER,PARAMETER :: ddp     = 2 ! dry deposition 
  INTEGER,PARAMETER :: wdp     = 3 ! wet deposition 
  INTEGER,PARAMETER :: fdda    = 4 ! four-dim assimilation 
  INTEGER,PARAMETER :: offline = 5 ! off-line emissions: 
                                   !=1, emission will be read from file
				   !=0, emission will be calculated during the model simulation (on-line emission)
  INTEGER,PARAMETER :: transport = 6 ! transported species 
                                   !=1, yes
				   !=0, not

  ! spaction(specie,[1=source,2=drydep,3=wetdep,4=fdda, 5=offline emission, 6=transport]) 
  ! attention : for aerosols,  mode_alloc(ispc) = spc_alloc(transport,imode,ispc)
  INTEGER,PARAMETER,DIMENSION(6,nmodes,nspecies) :: spc_alloc=RESHAPE((/ &
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf akk
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf acc
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf dd1
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf ds1
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf dd2
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf ds2
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf ssa
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf ssc
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf sss
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf occ
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf bc1
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf bc2
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf bc3
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf ocs
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf dbc
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf boc
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf bcs
    0 , 0 , 0 , 0 , 0 ,  & ! aer sulf mxx

!
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar akk
    1 , 1 , 1 , 0 , 0 ,  & ! aer bcar acc
    1 , 1 , 1 , 0 , 0 ,  & ! aer bcar dd1
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar ds1
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar dd2
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar ds2
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar ssa
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar ssc
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar sss
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar occ
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar bc1
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar bc2
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar bc3
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar ocs
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar dbc
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar boc
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar bcs
    0 , 0 , 0 , 0 , 0 ,  & ! aer bcar mxx
!
!
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar akk
    1 , 1 , 1 , 0 , 0 ,  & ! aer ocar acc
    0 , 1 , 1 , 0 , 0 ,  & ! aer ocar dd1
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar ds1
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar dd2
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar ds2
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar ssa
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar ssc
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar sss
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar occ
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar bc1
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar bc2
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar bc3
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar ocs
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar dbc
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar boc
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar bcs
    0 , 0 , 0 , 0 , 0 ,  & ! aer ocar mxx
!
!
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust akk
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust acc
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust dd1
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust ds1
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust dd2
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust ds2
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust ssa
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust ssc
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust sss
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust occ
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust bc1
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust bc2
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust bc3
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust ocs
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust dbc
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust boc
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust bcs
    0 , 0 , 0 , 0 , 0 ,  & ! aer dust mxx
!
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas akk
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas acc
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas dd1
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas ds1
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas dd2
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas ds2
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas ssa
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas ssc
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas sss
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas occ
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas bc1
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas bc2
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas bc3
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas ocs
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas dbc
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas boc
    0 , 0 , 0 , 0 , 0 ,  & ! aer seas bcs
    0 , 0 , 0 , 0 , 0    & ! aer seas mxx
!
    /),(/6,nmodes,nspecies/))
  
  ! effective particle radius (meter)
  REAL,PARAMETER,DIMENSION(nmodes,nspecies) :: part_radius=RESHAPE((/ &
!------------------------------------------------------------------------------------------------------------------
!-bin size  1          2        3        4          5       6       7       8     9    10
!------------------------------------------------------------------------------------------------------------------
     	  1.95e-7 , 1.95e-7 , 1.95e-7 , 999., 999., 999., 999., 999., 999., 999.,   & ! sdust
     	  1.95e-7 , 1.95e-7 , 1.00e-5 , 999., 999., 999., 999., 999., 999., 999.,   & ! bburn (0, pm25, pm10) meters
     	  1.95e-7 , 1.95e-7 , 1.95e-7 , 999., 999., 999., 999., 999., 999., 999.,   & ! urban
     	  1.95e-7 , 1.95e-7 , 1.95e-7 , 999., 999., 999., 999., 999., 999., 999.,   & ! bioge
     	  1.95e-7 , 1.95e-7 , 1.95e-7 , 999., 999., 999., 999., 999., 999., 999.,   & ! marin
     	  0.98e-6,  2.93e-6,  5.89e-6, 11.72e-6, 23.44e-6, 46.88e-6, 93.75e-6, 0.1875e-3, 0.375e-3, 0.750e-3   & ! v_ash
!---------------------------------------------------------- ---------------------------- ----------------------------
    /),(/nmodes,nspecies/))




  ! particle density kg/m^3
  REAL,PARAMETER,DIMENSION(nmodes,nspecies) :: part_dens=RESHAPE((/ &
!------------------------------------------------------------------------------------------------------------------
!-bin size 1   2   3  4   5   6   7   8   9  10
!------------------------------------------------------------------------------------------------------------------
    2.65e+3 , 2.65e+3 , 2.65e+3 , 999., 999., 999., 999., 999., 999., 999.,   & ! sdust
    1.35e+3 , 1.35e+3 , 1.35e+3 , 999., 999., 999., 999., 999., 999., 999.,   & ! bburn (0, pm25, pm10) kg/m^3
    1.35e+3 , 1.35e+3 , 1.35e+3 , 999., 999., 999., 999., 999., 999., 999.,   & ! urban
    1.35e+3 , 1.35e+3 , 1.35e+3 , 999., 999., 999., 999., 999., 999., 999.,   & ! bioge
    2.20e+3 , 2.20e+3 , 2.20e+3 , 999., 999., 999., 999., 999., 999., 999.,   & ! marin
    2.50e+3 , 2.50e+3 , 2.50e+3 , 2.50e+3, 2.50e+3, 2.50e+3, 2.50e+3, 2.50e+3, 2.50e+3, 2.50e+3    & ! v_ash
!------------------------------------------------------------------------------------------------------------------
    /),(/nmodes,nspecies/))
 
END MODULE aer1_list
