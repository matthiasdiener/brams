!> @brief This module adapts the data from/to CCATT-BRAMS to MAtrix
!! it uses blocks of athmosferic points to prepare the data. 
!! @author Luiz Flavio
!! @date Oct/2011
MODULE DriverMatrix   
      USE mem_aerosol, ONLY: &
                              matrix2mode, &
                              matrix2specie
      USE memMatrix
      USE mem_aer1, ONLY: aer1_g
      USE mem_chem1, ONLY:               &
         CHEMISTRY
      USE aer1_list, ONLY: mode_alloc
      USE mem_aerosol, ONLY: aero
         
      PUblic ::  MatrixDrive
      
      LOGICAL :: Test_firstTime=.true.
      LOGICAL, PARAMETER :: test=.false.
      INTEGER, PARAMETER :: test_type=1 !1 - normal , 2= fixed
      LOGICAL :: first_call=.true.
      !Correlacao entre naer e aerosois normais
      !                                01,02,03,04,05,06,07,08,09,10
      INTEGER, DIMENSION(29) :: nspc=(/00,00,00,00,00,00,01,00,00,01, &
                                       00,00,05,00,00,02,00,00,04,00, &
                                       00,03,00,00,00,00,00,00,00/)
      !                                01,02,03,04,05,06,07,08,09,10
      INTEGER, DIMENSION(29) :: nmde=(/00,00,00,02,00,00,00,00,00,00, &
                                       00,02,02,00,02,02,00,03,03,00, &
                                       03,03,00,00,00,00,00,00,00/)
      !                                   01     02       03      04     05      06      07      08      09      10
      logical, DIMENSION(29) :: doit=(/.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
                                       .false.,.false.,.true. ,.false.,.false.,.true. ,.false.,.false.,.true. ,.false., &
                                       .false.,.true. ,.false.,.false.,.false.,.false.,.false.,.false.,.false./)
                      

CONTAINS

   !> @brief To prepare and adapt data from brams to Matrix
   !! @author Luiz Flavio
   !! @date Oct/2011
   !! @todo verify sc_t of aerosols
   !! @todo Verify the wupdraft because I use w wind
   !! @todo aqso4rate is fixed using the box value, must be changed
   !! @todo emis_map is fixed to zero. Must be changed to real values
   subroutine MatrixDrive(ia,iz,ja,jz,m1,m2,m3)
      USE mem_basic, ONLY: basic_g
      use mem_grid, ONLY: dtlt,ngrid,time,nzg
      USE rconstants, ONLY: cp, cpor, p00, pi180, cpi ! INTENT(IN)
      USE ModParticle, ONLY: copyDensRadius
      USE mem_radiate,  ONLY: radfrq
      USE mem_leaf, ONLY: &
                        leaf_g
      USE mem_turb         ,  ONLY :  &
                                 turb_g
      USE mem_micro        ,  ONLY :  &
                                 micro_g
      USE mem_grid         ,  ONLY :  &
                                 grid_g, &
                                 zt


      implicit none

      !ngases     = 3      ! number of gas-phase species
      !nmass_spcs = 5      ! total number of mass species
      !gas_h2so4  = 1      !\
      !gas_hno3   = 2      !-indices in the gas array
      !gas_nh3    = 3      !/
 
      INTEGER, INTENT(IN) :: ia !< i initial point
      INTEGER, INTENT(IN) :: iz !< i final point
      INTEGER, INTENT(IN) :: ja !< j initial point
      INTEGER, INTENT(IN) :: jz !< j final point
      INTEGER, INTENT(IN) :: m1 !< z size               
      INTEGER, INTENT(IN) :: m2 !< i size               
      INTEGER, INTENT(IN) :: m3 !< j size               
       
      INTEGER :: i,j,k,noc,point,na,ng,gasOfChem,md,sp,nm,nsp,nmd
      double precision :: picpi,tstep,pblht
      DOUBLE PRECision,DIMENSION(m1,ia:iz,ja:jz,nmodes) :: ac  
      DOUBLE PRECision,DIMENSION(m1,ia:iz,ja:jz,nmodes) :: dens_mode_dry 
      REAL,PARAMETER :: tkethrsh=0.001       !   tke threshold for PBL height in m2/s2
                                             !   tkmin    = 5.e-4   minimum TKE in RAMS 
      REAL,PARAMETER :: rcmin=1.e-4          !   min liq water = 0.1 g/kg
      INTEGER :: toci,tocj,fcoi,fcoj
 
      IF (.not. (mod(time + .001,aerfrq) .lt. dtlt .or. time .lt. 0.001)) RETURN

      CALL allocateMatrix(ia,iz,ja,jz,m1,m2,m3,nzg) !Allocating the matrix with right size
      CALL MakeCrossMap(ia,iz,ja,jz,m1) !Doing point to point cross reference
      
      !Pressure, air temperature and Moisture
      DO i=ia,iz
         do j=ja,jz
            noc=nColumn(i,j)
            matrixVar(noc)%ustar(1)=leaf_g(ngrid)%ustar(i,j,1) !verificar npatch
            matrixVar(noc)%tstar(1)=leaf_g(ngrid)%tstar(i,j,1) !   "        "
            DO k = 1,m1
               !calculating and moving to right position
               picpi = (basic_g(ngrid)%pi0(k,i,j) + basic_g(ngrid)%pp(k,i,j)) * cpi
               !Pressure
               matrixVar(noc)%pres(k)= dble(p00 * picpi ** cpor)
               !Temperature
               matrixVar(noc)%tk(k) = dble(basic_g(ngrid)%theta(k,i,j) * picpi)
               !Humidity
               matrixVar(noc)%rh(k)=dble(basic_g(ngrid)%rv(k,i,j)) !*100.0)
               !Vertical Wind Component
               matrixVar(noc)%wupdraft(k)=dble( max(1.0,basic_g(ngrid)%wp(k,i,j)) )
               matrixVar(noc)%dn0(k)=basic_g(ngrid)%dn0(k,i,j)  
               matrixVar(noc)%pi0(k)=basic_g(ngrid)%pi0(k,i,j)  
               matrixVar(noc)%pp(k)=basic_g(ngrid)%pp(k,i,j) 
            END DO
            matrixVar(noc)%Zi = 0.
            !- convective layer
            if(turb_g(ngrid)%sflux_t(i,j) >= 1.e-8) then 
              pblht=0.
              do k=2,m1-1
                 pblht=zt(k)*grid_g(ngrid)%rtgt(i,j) 
                 if( micro_g(ngrid)%rcp(k,i,j) .gt. rcmin     ) EXIT ! dry convective layer
                 if( turb_g(ngrid)%tkep(k,i,j) .le. tkethrsh  ) EXIT 
              enddo
              matrixVar(noc)%Zi=pblht
            endif
         END DO
      END DO

      !Must be changed
      DO noc=1,NumberOfColunms
         matrixVar(noc)%aqso4rate(:)=0.5D-4 !Fixed
         matrixVar(noc)%emis_mas(:,:)=0.0 !No emission of any aerossols
         matrixVar(noc)%aero=0.0
      END DO
      
      !Gas input
      DO i=ia,iz
         do j=ja,jz
            noc=nColumn(i,j)
            DO k=1,m1
               matrixVar(noc)%gas(k,1)=0.0
               matrixVar(noc)%gas(k,2)=0.0
               matrixVar(noc)%gas(k,3)=0.0            
            END DO
         END DO
      END DO
      CALL fillAero(ia,iz,ja,jz,m1) !For the first time copy aerosol data to Matrix
            
      !Fill the matrixVar type with the data from aerosol ammount from model
      DO i=ia,iz
         do j=ja,jz
            noc=nColumn(i,j)
            DO k = 1,m1
               DO nmd=1,nmodes
                  DO nsp=0,nmass_spcs
                     If (.not. aero(1,nsp,nmd)%specieIsPresent) CYCLE
                     matrixVar(noc)%aero(k,aero(1,nsp,nmd)%MapInMatrix)= 0.0!&
                           !aero(ngrid,nsp,nmd)%ammount(k,iPos(noc),jPos(noc))
                  END DO
               END DO
            END DO
         END DO
      END DO
      
      tstep=aerfrq !Fazendo tstep=aerfrq está correto???????????????????????????????????????????                 
      !Calling the matrix adapter (will be changed in future to blocking and vectorizing)

      DO noc=1,NumberOfColunms
         CALL Matrix(matrixVar(noc),tstep,m1,noc)
      END DO
      IF(test) THEN
         !CALL compareTest(testData,matrixVar(1),1,m1)
      END IF
      
      DO noc=1,NumberOfColunms
         DO k=1,m1
            DO nmd=1,nmodes
               DO nsp=0,nmass_spcs
                  If (.not. aero(1,nsp,nmd)%specieIsPresent) CYCLE
                     aero(ngrid,nsp,nmd)%ammount(k,iPos(noc),jPos(noc))= &
                        matrixVar(noc)%aero(k,aero(ngrid,nsp,nmd)%MapInMatrix)
               END DO
            END DO
         END DO
      END DO
      
      !CALL DumpAero(ia,iz,ja,jz,m1)
      !***************************************************************!
      !TESTES TESTES TESTES TESTES TESTES TESTES TESTES TESTES        !
      !***************************************************************!
                                !
      !CALL deAllocateMatrix !Deallocating the matrix before return    !
      RETURN                                                          !
                                                                      !
      !***************************************************************!
      
   
   END subroutine MatrixDrive
   
   !> @brief create a xreference map between ia,iz,ja,jz,m1 and pointOfBlock  \n
   !! @author Luiz Flavio
   !! @date Oct/2011
   subroutine makeCrossMap(ia,iz,ja,jz,m1)
      implicit none
      
      INTEGER, INTENT(IN) :: ia,iz,ja,jz,m1
      INTEGER :: i,j,k,nb,p,column
      column=0
      do i=ia,iz
         do j=ja,jz
               column=column+1
               iPos(column)=i
               jPos(column)=j
               nColumn(i,j)=column
         END DO
      END DO

   END subroutine makeCrossMap
   
   
   SUBROUTINE fillAero(ia,iz,ja,jz,m1)
      use node_mod, only: &
             mynum
   
      INTEGER, INTENT(IN) :: ia,iz,ja,jz,m1
  
      INTEGER :: i,j,k,nmd,nsp
   
      !IF (mynum==1) THEN
      !   DO nmd=1,nmodes
      !      DO nsp=0,nmass_spcs
      !          WRITE(*,fmt='("MAP: ",2(I2.2,1X),L1,1X,I2.2)') &
      !             nmd,nsp,(aero(1,nsp,nmd)%specieIsPresent),(aero(1,nsp,nmd)%MapInMatrix)
      !      END DO
      !   END DO
      !END IF
      
   
      !Test copy loop - depende de tabela nova
      IF (CHEMISTRY>=1 .and. first_call ) then
      DO i=ia,iz
         do j=ja,jz
            DO k = 1,m1
               DO nmd=1,nmodes
                  DO nsp=0,nmass_spcs
                     If (aero(1,nsp,nmd)%specieIsPresent) THEN
                        IF(.not. doit(aero(1,nsp,nmd)%MapInMatrix)) CYCLE
                        IF(mode_alloc(nmde(aero(1,nsp,nmd)%MapInMatrix),nspc(aero(1,nsp,nmd)%MapInMatrix))==1) THEN

                           !IF (mynum==1) print *, 'LFR-DBG: ',k,i,j,aero(1,nsp,nmd)%MapInMatrix, &
                           !nspc(aero(1,nsp,nmd)%MapInMatrix),nmde(aero(1,nsp,nmd)%MapInMatrix), &
                           !mode_alloc(nmde(aero(1,nsp,nmd)%MapInMatrix),nspc(aero(1,nsp,nmd)%MapInMatrix)), &
                           !associated(aer1_g(nmde(aero(1,nsp,nmd)%MapInMatrix),nspc(aero(1,nsp,nmd)%MapInMatrix),1)%sc_p)
  

                           aero(ngrid,nsp,nmd)%ammount(k,iPos(noc),jPos(noc))= &
                           aer1_g(nmde(aero(1,nsp,nmd)%MapInMatrix),nspc(aero(1,nsp,nmd)%MapInMatrix),1)%sc_p(k,i,j)
                        END IF
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO

      END IF
      
      first_call=.false.
      
   END SUBROUTINE fillAero
   
   SUBROUTINE DumpAero(ia,iz,ja,jz,m1)

   INTEGer, intent(in) :: ia,iz,ja,jz,m1
   integer :: nsp,nmd,i,j,k

   DO nsp=0,nmass_spc
      DO nmd=1,nmodes
         If (aero(1,nsp,nmd)%specieIsPresent) THEN
            IF(.not. doit(aero(1,nsp,nmd)%MapInMatrix)) CYCLE
            IF(mode_alloc(nmde(aero(1,nsp,nmd)%MapInMatrix),nspc(aero(1,nsp,nmd)%MapInMatrix))==1) THEN
               write(77,FMT='(A)') aero(1,nsp,nmd)%aername
               DO i=ia,iz
                  DO j=ja,jz
                     WRite (77,FMT='(2(I2.2,1X),38(E18.8,1X))') i,j,(aero(1,nsp,nmd)%ammount(k,i,j),k=1,38)
                  END DO
               END DO
               write(77,FMT='(A)') '-----------------------------------'
            END IF
         END IF
      END DO
   END DO
   
   
   END SUBROUTINE DumpAero

   
END MODULE DriverMatrix  

