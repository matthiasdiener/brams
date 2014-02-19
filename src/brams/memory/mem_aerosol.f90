module mem_aerosol
   
   TYPE aer_var
      REAL, POINTER  :: ammount(:,:,:) !total of number or mass nzp,nxp,nyp
      CHARACTER(LEN=13) :: aerName                 !Name for specie of aerosol
      LOGICAL :: SpecieIsPresent                   !If exist the specie/mode in mechanism
      INTEGER :: MapInMatrix                       !Map from aerosol (mode,specie) to matrix position
   END TYPE aer_var
   TYPE(aer_var),ALLOCATABLE :: aero(:,:,:), aero_m(:,:,:) !ngrids,nspec,nmodes
   integer,ALLOCATABLE :: matrix2mode(:) !Map back from matrix aerosol position to mode
   INTEGER,ALLOCATABLE :: matrix2specie(:) !Map back from matrix aerosol position to specie
   
   contains
   
   SUBROUTINE alloc_aerosol(ngrids,n1,n2,n3)
      USE memMatrix, ONLY : &
                            nmodes, & !Number of modes
                            nmass_spcs, &  !number of mass species
                            naerobox, & !The total number of aerosols (number and mass)
                            aero_spcs, &       !Aerosol species names in Matrix
                            mech      , &!Selected mechanism
                            naerovars, &
                            nextra, &
	                    nWeights, &
                            nPoints , &
                            allocateMatrixSetup !Subroutine to allocate setup
   
      use node_mod, only: mynum
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ngrids
      INTEGER, INTENT(IN) :: n1(ngrids),n2(ngrids),n3(ngrids)
      integer :: ngrid,nsp,nmd

      !Fill the aerosol modes from the TRAMP_param and for the
      !selected mechanism. In each mechanism there are a specific
      !number of different modes of aerosol
      !WRITE (*,*) 'Matrix: Ajusting aerosol setup and sizes'
      SELECT CASE(mech)
         CASE (1)
            naerovars=51;nextra=3;nmodes=16
	    nAeroBox=nAeroVars+nExtra
	    nWeights=nModes*nPoints
         CASE (2)
            naerovars=51;nextra=3;nmodes=16
	    nAeroBox=nAeroVars+nExtra
	    nWeights=nModes*nPoints
          CASE (3)
            naerovars=41;nextra=3;nmodes=13
	    nAeroBox=nAeroVars+nExtra
	    nWeights=nModes*nPoints
         CASE (4)
            naerovars=34;nextra=1;nmodes=10
	    nAeroBox=nAeroVars+nExtra
	    nWeights=nModes*nPoints
         CASE (5)
            naerovars=45;nextra=3;nmodes=14
	    nAeroBox=nAeroVars+nExtra
	    nWeights=nModes*nPoints
         CASE (6)
            naerovars=45;nextra=3;nmodes=14
	    nAeroBox=nAeroVars+nExtra
	    nWeights=nModes*nPoints
         CASE (7)
            naerovars=35;nextra=3;nmodes=11
	    nAeroBox=nAeroVars+nExtra
	    nWeights=nModes*nPoints
          CASE (8)
            naerovars=28;nextra=1;nmodes= 8
	    nAeroBox=nAeroVars+nExtra
	    nWeights=nModes*nPoints
      CASE DEFAULT
         PRINT *, 'Error: invalid mechanism for aerosol. Must be 1 until 8'
	 STOP
      END SELECT
      IF (mynum==1) THEN
         write (*,*) '==== Matrix Ajusted ===='
         write (*,FMT='(A,I2.2)') 'Mechanism: ',mech
         write (*,FMT='(A,I2.2)') 'naerovars: ',naerovars
         write (*,FMT='(A,I2.2)') 'nextra   : ',nextra   
         write (*,FMT='(A,I2.2)') 'nmodes   : ',nmodes   
         write (*,FMT='(A,I2.2)') 'nAeroBox : ',nAeroBox 
         write (*,FMT='(A,I2.2)') 'nWeights : ',nWeights 
         write (*,FMT='(A,I2.2)') 'nPoints  : ',nPoints  
         write (*,*) '========================'
      END IF


      ALLOCATE(aero(ngrids,0:nmass_spcs,nmodes)) !0 indicates number, >0 mass
      ALLOCATE(aero_m(ngrids,0:nmass_spcs,nmodes)) !0 indicates number, >0 mass

      DO ngrid=1,ngrids
         DO nsp=0,nmass_spcs
            DO nmd=1,nmodes
               ALLOCATE(aero(ngrid,nsp,nmd)%ammount(n1(ngrid) &
                                         ,n2(ngrid),n3(ngrid)))
                aero(ngrid,nsp,nmd)%ammount=0.0
                aero(ngrid,nsp,nmd)%SpecieIsPresent=.false.
                aero(ngrid,nsp,nmd)%MapInMatrix=0
               ALLOCATE(aero_m(ngrid,nsp,nmd)%ammount(n1(ngrid) &
                                         ,n2(ngrid),n3(ngrid)))
                aero_m(ngrid,nsp,nmd)%ammount=0.0
                aero_m(ngrid,nsp,nmd)%SpecieIsPresent=.false.
                aero_m(ngrid,nsp,nmd)%MapInMatrix=0
            END DO
         END DO
      END DO
      ALLOCATE(matrix2mode(naerobox)) 
      ALLOCATE(matrix2specie(naerobox))

      CALL allocateMatrixSetup()
      
   END SUBROUTINE alloc_aerosol
   
   SUBroutine fill_aerosol(ngrids,n1,n2,n3)
      USE memMatrix, ONLY: &
                            nmass_spcs, & !number of mass species
                            chem_spc_name, & !Name of species
                            ngases, & !Number of gas species
                            gas_name, & !Name of gases
                            nmodes_max, &  !maximum number os modes
                            modes1,modes2,modes3,modes4,modes5,modes6,modes7,modes8,& !Map of modes for each mech
                            mname , &!Name of each mode (complete)
                            mspcs, & !Relation between aerosol species and modes
                            nmodes, & !Number of modes
                            naerobox, & !The total number of aerosols (number and mass)
                            aero_spcs, &       !Aerosol species names in Matrix
                            mech      , &!Selected mechanism
                            naerovars, &
                            nextra, &
	                    nWeights, &
                            nPoints

      USE Aero_setup, ONLY: &
                            Setup_Config, & !Create and setup species in Matrix
                            Setup_Species_Maps, &
                            Setup_Aero_Mass_Map, & 
                            Setup_Coag_Tensors, &
                            Setup_Dp0, &
                            Setup_Emis !, &  
                            !Setup_Kci   
                            
      USE Aero_coag,   ONLY: &
                            setup_kij_diameters,       & !subroutine
                            setup_kij_tables,          & !subroutine
                            get_kbarnij                  !subroutine
      
      USE Aero_subs,   ONLY: &
                            no_microphysics            !intent()                            

      USE Aero_npf,    ONLY: &
                            setup_npfmass              !subroutine
      use node_mod, only: mynum
      
      INTEGER, INTENT(IN) :: ngrids,n1(ngrids),n2(ngrids),n3(ngrids)
      INTEGER :: ngrid,nsp,nmd,i,j,k
      INTEGER :: modSel(nmodes_max)
      CHARACTER(LEN=4) :: prefix,sname
      CHARACTER(LEN=3) :: cmode
      LOGICAL, PARAMETER :: test=.true.
      
      !Fill the aerosol modes from the TRAMP_param and for the
      !selected mechanism. In each mechanism there are a specific
      !number of different modes of aerosol
      !WRITE (*,*) 'Matrix: Ajusting aerosol setup and sizes'
      SELECT CASE(mech)
         CASE (1)
            modSel(1:size(modes1))=modes1
         CASE (2)
            modSel(1:size(modes2))=modes2
          CASE (3)
            modSel(1:size(modes3))=modes3
         CASE (4)
            modSel(1:size(modes4))=modes4
         CASE (5)
            modSel(1:size(modes5))=modes5
         CASE (6)
            modSel(1:size(modes6))=modes6
         CASE (7)
            modSel(1:size(modes7))=modes7
          CASE (8)
            modSel(1:size(modes8))=modes8
      CASE DEFAULT
         PRINT *, 'Error: invalid mechanism for aerosol. Must be 1 until 8'
	 STOP
      END SELECT
      
      !Creating aerosol species name. The name is select from all the
      !names in mname but just for selected modes.
      !In this case we used just NUM because is the sum of number 
      !of particles in each mode
      DO nmd=1,nmodes
         aero(:,0,nmd)%aerName='NUMB_'//mname(modSel(nmd))
         aero(:,0,nmd)%specieIsPresent=.true.
         !print *,'Number: ', aero(:,0,nmd)%aerName; CALL flush(6)
      END DO
         
      DO nmd=1,nmodes
         DO nsp=1,nmass_spcs
            !Check if the specie nsp has the mode nmd
            IF(mspcs(nsp,modSel(nmd)) == 0) THEN
               aero(:,nsp,nmd)%aerName='Unused_in_mec'
               CYCLE
            END IF            
            aero(:,nsp,nmd)%aerName='MASS_'//mname(modSel(nmd))// &
                                    '_'//chem_spc_name(nsp)
            aero(:,nsp,nmd)%specieIsPresent=.true.
            !print *,'Mass: ', aero(:,nsp,nmd)%aerName; CALL flush(6)
         END DO
      END DO
      IF (mynum==1) write (*,*) 'Calling matrix setup routines'
      CALL Setup_Config()
      call Setup_Species_Maps
      call Setup_Aero_Mass_Map 
      if ( .not. no_microphysics ) then
         call Setup_Coag_Tensors
         call Setup_Dp0     
         call Setup_Kij_Diameters
         call Setup_Kij_Tables    
         call Setup_Emis  
         !call Setup_Kci   
         call Setup_Npfmass
      END IF  
      IF (mynum==1) write (*,*) 'Giving Aerosol Names'
      Do naer=1,naerobox
         !Cut the prefix (NUMB or MASS) and the Mode
         prefix=aero_spcs(naer)(1:4)
         cmode= aero_spcs(naer)(6:8)
         sname= aero_spcs(naer)(10:13)
         IF(prefix=='NUMB') THEN !Is number os particles
            DO nmd=1,nmodes
               If (.not. aero(1,0,nmd)%specieIsPresent) CYCLE
               IF( aero(1,0,nmd)%aerName(6:8)==cmode) THEN
                  aero(:,0,nmd)%MapInMatrix=naer
                  matrix2mode(naer)=nmd
                  matrix2specie(naer)=0
                  exit
               END IF
            END DO
         ELSE !Is mass of particles
            DO nmd=1,nmodes
               DO nsp=1,nmass_spcs
                  If (.not. aero(1,nsp,nmd)%specieIsPresent) CYCLE
                  IF(aero(1,nsp,nmd)%aerName(6:8)==cmode .and. aero(1,nsp,nmd)%aerName(10:13)==sname) THEN
                     aero(:,nsp,nmd)%MapInMatrix=naer
                     matrix2mode(naer)=nmd
                     matrix2specie(naer)=nsp                                         
                  END IF
               END DO
            END DO
         END IF
!         write(*,fmt='(A)')  '====== Aerosol naer,nsp,nmd,Name ======'      
!         if (matrix2specie(naer)/=0 .and. matrix2mode(naer)/=0) &
!                 write(*,fmt='(3(I2.2,1X),A)') naer,matrix2specie(naer),matrix2mode(naer),trim(aero(1,matrix2specie(naer),matrix2mode(naer))%aerName)
      END DO
      
   END SUBROUTINE fill_aerosol
   
  subroutine filltab_aerosol(ngrids,imean,n1,n2,n3)
    USE memMatrix, ONLY : &
                            nmodes, & !Number of modes
                            nmass_spcs  !number of mass species

    use var_tables, only: &
                           InsertVTab
    implicit none
    include "i8.h"
    integer, intent(in) :: ngrids,imean
    INTEGER, INTENT(IN) :: n1(ngrids),n2(ngrids),n3(ngrids)
    integer(kind=i8) :: npts
    INTEGER :: ngrid,nsp,nmd,i
    
    ! Fill pointers to arrays into variable tables    
    DO ngrid=1,ngrids
       npts = n1(ngrid)*n2(ngrid)*n3(ngrid)
       DO nmd=1,nmodes
          DO nsp=0,nmass_spcs
             if (associated(aero(ngrid,nsp,nmd)%ammount) .and. aero(ngrid,nsp,nmd)%SpecieIsPresent)  &
             call InsertVTab (aero(ngrid,nsp,nmd)%ammount,aero_m(ngrid,nsp,nmd)%ammount  &
            ,ngrid, npts, imean,  &
            aero(ngrid,nsp,nmd)%aername//' :3:anal:mpti:mpt3')
          END DO
       END DO
    END DO
    
  end subroutine filltab_aerosol

end module mem_aerosol
