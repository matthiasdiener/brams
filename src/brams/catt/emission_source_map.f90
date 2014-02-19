module emission_source_map
  use ModNamelistFile, only: namelistFile

  integer            :: nvar
  integer            :: nlat
  integer            :: nlon
  real               :: lonW
  real               :: latS
  real               :: dlon
  real               :: dlat

  character(len=256) :: firemapfn ! from RAMSIN
  character(len=256) :: tracersfn ! from RAMSIN

  integer, parameter :: nsources=3
  integer, parameter :: ngases=4
  integer, parameter :: nplumesource=12
  integer, parameter :: nz_orig=50
  character(len=10), parameter :: gas(ngases)=(/&
       'CO        ', &
       'COstc     ', &
       'PM25      ', &
       'COANT     '/)

  type qsc_type
     real, pointer :: q(:,:,:)
  end type qsc_type

  type(qsc_type), allocatable :: qsc(:)

  integer :: iplume(nplumesource)
  data (iplume(i),i=1,nplumesource) / &
       01, 02, 03, & !arrays with mean area for bioma
       04, 05, 06, & !arrays with STD of mean area for bioma
       11, 12, 13, & !CO flaming emission for bioma
       21, 22, 23  / !PM25 flaming emission for bioma

  integer :: plumerise ! from RAMSIN ( = 1 on, /= 1 off)

  integer, parameter :: k_CO_total   =  9 ! stores CO total emission
  integer, parameter :: k_CO_smold   = 10 ! stores CO emis. in smoldering phase
  integer, parameter :: k_PM25_total = 19 ! stores PM25 total emisision
  integer, parameter :: k_PM25_smold = 20 ! stores PM25 emis. smoldering phase

  integer :: define_proc ! from RAMSIN

  contains

    subroutine read_emission_sources_map()

      use mem_basic, only: basic_g ! intent(in)
      use mem_grid, only: ngrids, time, iyear1, imonth1, idate1, & ! intent(in)
           grid_g, dzt, nnxp, nnyp, nnzp, & ! intent(in)
           xt, yt, platn, plonn ! intent(in)
      use mem_scalar, only: scalar_g ! intent(inout)
      use extras, only: &
           na_extra3d, & ! intent(in)
           extra3d       ! intent(inout)
      use node_mod, only: mchnum, nmachs, master_num, & ! intent(in)
           ixb, ixe, iyb, iye, & ! intent(in)
           mynum, nodeia, nodeja, nodeiz, nodejz, & ! intent(in)
           nodei0, nodej0, nodemxp, nodemyp, nodemzp ! intent(in)
      use ParLib, only: parf_bcast ! Subroutine

      implicit none
      include "i8.h"
      include "files.h"
      ! Local Variables:
      character(len=f_name_length) :: fname
      character(len=2)             :: cgrid
      character(len=10)            :: tracer
      logical                      :: there
      integer                      :: ifname, kgas, len1
      real, pointer                :: qsc_orig_map(:,:,:)
      real, pointer                :: extra3d_orig(:,:,:)
      real, allocatable            :: lat_orig(:,:), lon_orig(:,:)
      integer                      :: ng
      integer                      :: i, j, ie, je, isrc
      real                         :: dist_ij, dist_iej, dist_iwj, &
           dist_ijn, dist_ijs
      integer                      :: ext, ierr, k, local_i, local_j !, countpt
      real                         :: minlat(ngrids), maxlat(ngrids)
      real                         :: minlon(ngrids), maxlon(ngrids)

      type latlon
         real, pointer :: lat(:,:)
         real, pointer :: lon(:,:)
      end type latlon

      type(latlon), allocatable :: globalGrid(:)

      ! scalar_g(IG,GRID)%srcsc    
      ! IG =>  Sources
      !        1 => CO 
      !        2 => COstc 
      !        3 => PM25 
      !        4 => COANT
      !        5 => COTOT

      ! Allocating globalGrid
      allocate(globalGrid(ngrids),STAT=ierr)
      if (ierr/=0) then
         call fatal_error("ERROR:impossible to allocate globalGrid")
      endif
      ! Nullifing pointers to receive a global Latitude and Longitude
      do ng = 1, ngrids
         nullify(globalGrid(ng)%lat)
         nullify(globalGrid(ng)%lon)
      enddo
      ! Allocating global Latitude and Longitude
      do ng = 1, ngrids
         allocate(globalGrid(ng)%lat(nnxp(ng),nnyp(ng)),STAT=ierr)
         if (ierr/=0) then
            call fatal_error&
                 ("ERROR:impossible to allocate globalGrid(ng)%lat")
         endif
         allocate(globalGrid(ng)%lon(nnxp(ng),nnyp(ng)),STAT=ierr)
         if (ierr/=0) then
            call fatal_error&
                 ("ERROR:impossible to allocate globalGrid(ng)%lon")
         endif
      enddo

      ! Calculating global Latitude and Longitude
      do ng = 1, ngrids
         do j = 1, nnyp(ng)
            do i = 1, nnxp(ng)
               call xy_ll(globalGrid(ng)%lat(i,j),globalGrid(ng)%lon(i,j),&
                    platn(ng),plonn(ng),xt(i),yt(j))
            enddo
         enddo
      enddo

      ! Nullifing pointer variables refered to original source map
      nullify(qsc_orig_map)
      nullify(extra3d_orig)
      ! The 2 arrays will be allocated inside subroutine read_sources_plume

      ! Allocating qsc_type with total of grids used
      if(.not. allocated(qsc)) then
         allocate(qsc(ngrids))
         ! Allocating qsc with dimension corresponding to the current grid
         do ng=1, ngrids
            allocate(qsc(ng)%q(nnzp(ng),nodemxp(mynum,ng),nodemyp(mynum,ng)), &
                 STAT=ierr)
            if (ierr/=0) then
               call fatal_error&
                    ("ERROR:impossible to allocate qsc(ng)%q")
            endif
            qsc(ng)%q(:,:,:)=0
         enddo
      endif

      ! Putting Zeros in Sources
      do ng=1, ngrids
         do kgas=1, ngases
            scalar_g(kgas,ng)%srcsc(:,:,:) = 0.
         enddo
      enddo

      if (mchnum==master_num) then !io-process only

         !Name of file with emissions map
         write(cgrid,'(a1,i1)') 'g', 1
         call makefnam(fname, firemapfn//' ', &
              time, iyear1, imonth1, idate1, 0*100, 'T', cgrid, 'vfm')

         ifname=len_trim(fname)
  
         !Open the file
         print *, "Opening EmissionMap File=", trim(fname)
         call flush(6)
         inquire(file=fname(1:ifname), exist=there)

      endif !io-process

      call parf_bcast(there, master_num)

      if(.not.there) then

         print*,'-----------------------------------------------------'
         print*,'FILE= ',fname(1:IFname),' does not exist.'
         
         !case 1 : stop the program
         if(define_proc == 1 .or. define_proc == 2) then
            print*, 'STOP at subroutine read_emission_sources'
            print*,'-----------------------------------------------------'
            call fatal_error(" Reading emission sources")
            
            !case 2 : Keeping previous values
         elseif(define_proc == 0)then
            print*,'Keeping previous values.'
            print*,'-----------------------------------------------------'
            return
            
         endif
         
      else

         ! EmissionMap File exists
         if(define_proc == 1 ) define_proc = 0

         ! Read the emission and plumerise data
         call read_sources_plume(plumerise, fname, &
              extra3d_orig, qsc_orig_map, mchnum, master_num)

         ! Allocating Latitude and Longitude for input grid
         if((.not. allocated(lat_orig)) .and. (.not.allocated(lon_orig))) &
              allocate(lat_orig(nlon,nlat), lon_orig(nlon,nlat))
         
         ! plumerise 
         if(plumerise == 1) then 
            
            ! Split the emission in two phases: flaming/smoldering
            call emis_flam_smold(nlon, nlat,&
                 extra3d_orig, qsc_orig_map) !extra3d_orig_smold, 
            
         endif

         ! Interpolating the input grid with local grid
         ! Creating LAT/LON data for emission input grid
         call set_latlon(nlon, nlat, lon_orig, lat_orig, &
              dlon, dlat, lonW, latS)

         ! Initiating arrays for each grid
         do ng = 1, ngrids
            !!qsc(ng)%q(:,:,:)        = 0.
            extra3d(6,ng)%d3(:,:,:) = 0.
         enddo

         do ng=1, ngrids
            minlat(ng) = minval(globalGrid(ng)%lat)
            maxlat(ng) = maxval(globalGrid(ng)%lat)
            minlon(ng) = minval(globalGrid(ng)%lon)
            maxlon(ng) = maxval(globalGrid(ng)%lon)
         enddo

!!$         countpt = 0 !DEBUG

         ! Loop in the input emission grid
         do je = 1, nlat
            
            do ie = 1, nlon
               
               ! Looking for points in local grids to set as emission point
               do ng=1, ngrids
                  
                  ! Discard points with null values
                  if (sum(qsc_orig_map(1:nsources,ie,je))==0.0) then
                     if (plumerise==1) then
                        if (sum(extra3d_orig(1:maxval(iplume),ie,je))==0.0)&
                             exit
                     else
                        exit
                     endif
                  endif
                  
                  ! Avoiding point outside local domain
                  if(  lat_orig(ie,je)<minlat(ng) .or. &
                       lat_orig(ie,je)>maxlat(ng) .or. &
                       lon_orig(ie,je)<minlon(ng) .or. &
                       lon_orig(ie,je)>maxlon(ng) )    &
                       cycle
                  
                  call newgrid(ng) ! Setting some parameters

                  ! Looking for local grid point in the neigborhood 
                  ! of emission point
                  
                  ! Local grid Initial indexes
                  i = 1
                  j = 1
                  
                  ! Calculating Euclidian distances from emission point
                  ! relationed to X axis
                  dist_ij  = sqrt( &
                       abs(lat_orig(ie,je)-globalGrid(ng)%lat(i,j))**2 + &
                       abs(lon_orig(ie,je)-globalGrid(ng)%lon(i,j))**2 )
                  dist_iwj = dist_ij
                  dist_iej = sqrt( &
                       abs(lat_orig(ie,je)-globalGrid(ng)%lat(i+1,j))**2 + &
                       abs(lon_orig(ie,je)-globalGrid(ng)%lon(i+1,j))**2 )
                  
                  ! Loop to looking for the closer point in local grid
                  do
                     ! look in the neighborhood in X axis
                     if (dist_iej<dist_ij) then
                        i = i + 1
                     elseif (dist_iwj<dist_ij) then
                        i = i - 1
                     endif
                     
                     ! Calculating Euclidian distances from emission point
                     ! relationed to Y axis
                     dist_ij = sqrt( &
                          abs(lat_orig(ie,je)-globalGrid(ng)%lat(i,j))**2 +&
                          abs(lon_orig(ie,je)-globalGrid(ng)%lon(i,j))**2 )
                     if (j>1) then
                        dist_ijs = sqrt( &
                             abs(&
                             lat_orig(ie,je)-globalGrid(ng)%lat(i,j-1))**2+&
                             abs(&
                             lon_orig(ie,je)-globalGrid(ng)%lon(i,j-1))**2 )
                     else
                        dist_ijs = dist_ij
                     endif
                     if (j<nnyp(ng)) then
                        dist_ijn = sqrt( &
                             abs(&
                             lat_orig(ie,je)-globalGrid(ng)%lat(i,j+1))**2+&
                             abs(&
                             lon_orig(ie,je)-globalGrid(ng)%lon(i,j+1))**2 )
                     else
                        dist_ijn = dist_ij
                     endif
                     
                     ! look in the neighborhood in Y axis
                     if (dist_ijn<dist_ij) then
                        j = j + 1
                     elseif (dist_ijs<dist_ij) then
                        j = j - 1
                     endif
                     
                     ! Calculating Euclidian distances from emission point
                     ! relationed to X axis
                     ! To be used in the next iteration
                     dist_ij = sqrt( &
                          abs(lat_orig(ie,je)-globalGrid(ng)%lat(i,j))**2 +&
                          abs(lon_orig(ie,je)-globalGrid(ng)%lon(i,j))**2 )
                     
                     if (i>1) then
                        dist_iwj = sqrt( &
                             abs(&
                             lat_orig(ie,je)-globalGrid(ng)%lat(i-1,j))**2+&
                             abs(&
                             lon_orig(ie,je)-globalGrid(ng)%lon(i-1,j))**2 )
                     else
                        dist_iwj = dist_ij
                     endif
                     if (i<nnxp(ng)) then
                        dist_iej = sqrt( &
                             abs(&
                             lat_orig(ie,je)-globalGrid(ng)%lat(i+1,j))**2+&
                             abs(&
                             lon_orig(ie,je)-globalGrid(ng)%lon(i+1,j))**2 )
                     else
                        dist_iej = dist_ij
                     endif
                     
                     ! Calculating Euclidian distances from emission point
                     ! relationed to Y axis
                     ! To be used in the next iteration
                     if (j>1) then
                        dist_ijs = sqrt( &
                             abs(&
                             lat_orig(ie,je)-globalGrid(ng)%lat(i,j-1))**2+&
                             abs(&
                             lon_orig(ie,je)-globalGrid(ng)%lon(i,j-1))**2 )
                     else
                        dist_ijs = dist_ij
                     endif
                     if (j<nnyp(ng)) then
                        dist_ijn = sqrt( &
                             abs(&
                             lat_orig(ie,je)-globalGrid(ng)%lat(i,j+1))**2+&
                             abs(&
                             lon_orig(ie,je)-globalGrid(ng)%lon(i,j+1))**2 )
                     else
                        dist_ijn = dist_ij
                     endif
                     
                     ! Stop Criteria
                     if ( dist_ij<=dist_iwj .and. dist_ij<=dist_iej .and. &
                          dist_ij<=dist_ijs .and. dist_ij<=dist_ijn ) exit
                     
                  enddo
                  
                  if ( lat_orig(ie,je)<globalGrid(ng)%lat(i,1) .or. &
                       lat_orig(ie,je)>globalGrid(ng)%lat(i,nnyp(ng)) &
                       .or. &
                       lon_orig(ie,je)<globalGrid(ng)%lon(1,j) .or. &
                       lon_orig(ie,je)>globalGrid(ng)%lon(nnxp(ng),j) ) &
                       cycle

                  ! Setting parameters relationed to the actual grid
                  call newgrid(ng)

                  local_i = i-nodei0(mynum,ng)
                  local_j = j-nodej0(mynum,ng)

                  ! Extract information from Global to Local domain
                  if ( i>=ixb(mynum,ng) .and. &
                       i<=ixe(mynum,ng) .and. &
                       j>=iyb(mynum,ng) .and. &
                       j<=iye(mynum,ng) ) then

!!$                     countpt = countpt + 1 !DEBUG
                  
                     ! Sum emission value to local grid point
                     do isrc=1, nsources !ngases
                        qsc(ng)%q(isrc,local_i,local_j) =      &
                             qsc(ng)%q(isrc,local_i,local_j) + &
                             qsc_orig_map(isrc,ie,je)
                     enddo
                  
                     if (plumerise==1) then
                        do isrc=1, maxval(iplume)
                           extra3d(6,ng)%d3(isrc,local_i,local_j) =      &
                                extra3d(6,ng)%d3(isrc,local_i,local_j) + &
                                extra3d_orig(isrc,ie,je)
                        enddo
                     endif

                  endif
                  
               enddo !ngrids
               
            enddo !nlon
               
         enddo !nlat

         ! Loop throw the gases

         print ('(A,I5,A)'), "Process: ", mynum, &
              " - gases and tracers summary:"
            
         do isrc=1, ngases

            print *,'---- > tracer: ',isrc,gas(isrc)
            call flush(6)

            tracer=gas(isrc)
            len1 = len_trim(tracer)+ 1
               
            do ng = 1, ngrids
                  
               ! Setting parameters relationed to the actual grid
               call newgrid(ng)
                  
               if(plumerise == 1) then
                  call plume2qsc(ng, isrc, extra3d(6,ng)%d3)!qsc(ng)%plume_smold
               endif

               ! Building sources in SCALAR_G
               scalar_g(isrc,ng)%srcsc(1,:,:) =      &
                    scalar_g(isrc,ng)%srcsc(1,:,:) + &
                    qsc(ng)%q(isrc,:,:)

               if(plumerise == 1) then
                  !- Putting emission smoldering in level 2 of array
                  scalar_g(isrc,ng)%srcsc(2,:,:) = &
                       scalar_g(isrc,ng)%srcsc(1,:,:)
               else
                  call reorganize_sources(nodemxp(mynum,ng), &
                       nodemyp(mynum,ng), &
                       scalar_g(isrc,ng)%srcsc)
               endif

               call convert_to_misture_ratio(nodemxp(mynum,ng), &
                    nodemyp(mynum,ng), basic_g(ng)%dn0, &
                    scalar_g(isrc,ng)%srcsc, grid_g(ng)%rtgt, dzt)

            enddo
  
         enddo

      endif !Found file

      ! Deallocating global Latitude and Longitude
      do ng = 1, ngrids
         deallocate(globalGrid(ng)%lat)
         deallocate(globalGrid(ng)%lon)
      enddo
      ! Deallocating globalGrid
      deallocate(globalGrid)

    end subroutine read_emission_sources_map
    
    !----------------------------------------------------------------------

    subroutine extGlob2Loc2d3d(global, nxp, nyp, i0, j0, local, qualproc, mynum)

      implicit none
      real,    intent(in   ) :: global(:,:,:)
      integer, intent(in   ) :: nxp
      integer, intent(in   ) :: nyp
      integer, intent(in   ) :: i0
      integer, intent(in   ) :: j0
      real,    intent(inout) :: local(:,:,:)
      integer, intent(inout) :: qualproc(:,:) !DEBUG
      integer, intent(in   ) :: mynum         !DEBUG
      
      integer :: i, ig
      integer :: j, jg
      
      do j= 1, nyp
         jg = j+j0
         do i = 1, nxp
            ig = i+i0
            local(:,i,j) = global(:,ig,jg)
            ! DEBUG
            qualproc(ig,jg) = mynum
         enddo
      enddo

    end subroutine extGlob2Loc2d3d

    !----------------------------------------------------------------------

    subroutine read_sources_plume(plumerise, fname, plume, qsc_orig, &
         mchnum, master_num) 
      
      use ParLib, only: parf_bcast ! Subroutine
      
      implicit none
      include "i8.h"
      include "files.h"

      ! Arguments:
      integer, intent(in)          :: plumerise
      character(len=f_name_length), intent(in) :: fname
      real, pointer                :: plume(:,:,:)
      real, pointer                :: qsc_orig(:,:,:)
      integer, intent(in)          :: mchnum
      integer, intent(in)          :: master_num
      ! Local Variables:
      integer            :: iunit, ifname, i
      character(len=128) :: cmessage
      character(len=9)   :: date
      
      iunit=2
      
      if (mchnum==master_num) then !io-process only
         
         ifname=len_trim(fname)
         print*,'-----------------------------------------------------'
         print*,'Opening emission file= ',fname(1:IFname)
         
         open(UNIT=iunit, FILE=fname(1:len_trim(fname)), FORM='formatted', STATUS='old')
         ! Reading header of source map file:
         ! Extrating information about size of grid map, geografic localization
         ! and resolution.
         do i=1,11
            read(iunit,'(A80)') cmessage
            !print*,cmessage(1:lastchar(cmessage))
         enddo
         read(iunit,*) cmessage,nlon,lonW,dlon
         !print*,cmessage(1:lastchar(cmessage)),nlon,lonW,dlon
         read(iunit,*) cmessage,nlat,latS,dlat
         !print*,cmessage(1:lastchar(cmessage)),nlat,latS,dlat
         read(iunit,*) cmessage,date
         !print*,cmessage(1:lastchar(cmessage)),date(1:lastchar(date))
         read(iunit,*) cmessage,nvar
         !print*,cmessage(1:lastchar(cmessage)),nvar
         
      endif !io-process
      
      ! Broadcasting array bounds
      call parf_bcast(nlon, master_num)
      call parf_bcast(lonW, master_num)
      call parf_bcast(dlon, master_num)
      call parf_bcast(nlat, master_num)
      call parf_bcast(latS, master_num)
      call parf_bcast(dlat, master_num)
      call parf_bcast(nvar, master_num)

      ! Allocating pointer array for source of gases and 
      ! particulate material
      if (.not.associated(qsc_orig)) allocate(qsc_orig(nsources,nlon,nlat))
      ! Allocating pointer array for emission sources data for Plume Rise
      if (.not.associated(plume))    allocate(plume(nz_orig,nlon,nlat))
      
      ! Initiating arrays for emission input grid
      ! OBS.: Using exclusive extra3D like array to plumerise model
      if(plumerise == 1) then
         qsc_orig(:,:,:) = 0.
         plume(:,:,:)    = 0.
      endif
      
      if (mchnum==master_num) then !io-process only

         ! Reading sources of gases and particulate material
         do i=1, nsources
            ! Reading header:
            read(iunit,'(A80)') cmessage
            print *, 'Reading emission file source: ', i, ':', &
                 trim(cmessage)
            read(iunit,'(A80)') cmessage
            print *, trim(cmessage)
            
            ! Reading data in VFM format
            CALL vfirec(iunit, qsc_orig(i,:,:),nlon*nlat,'LIN')
            ! Eliminating negative values
            qsc_orig(i,:,:)=max(0.,qsc_orig(i,:,:))
         enddo
         
         ! Reading sources of emission data for Plume Rise
         if (plumerise == 1) then
            !- plumerise start
            do i=1, nplumesource
               ! Reading header:
               read(iunit,'(A80)') cmessage
               print *, 'Reading PLUMERISE emission file source: ', iplume(i) !i
               PRINT *, trim(cmessage)
               read(iunit,'(A80)') cmessage
               print *, trim(cmessage)
               
               ! Reading data in VFM format
               CALL vfirec(iunit,plume(iplume(i),:,:),nlon*nlat,'LIN')
               ! Eliminating negative values
               plume(iplume(i),:,:)=max(0.,plume(iplume(i),:,:))
            end do
         endif
         
         print*,'-----------------------------------------------------'
         close (iunit)
         
      endif ! io-process

      ! Broadcasting arrays
      call parf_bcast(qsc_orig, &
           int(nsources,i8), int(nlon,i8), int(nlat,i8), master_num)
      call parf_bcast(plume, &
           int(nz_orig,i8),  int(nlon,i8), int(nlat,i8), master_num)
      
    end subroutine READ_sources_plume

  !----------------------------------------------------------------------

  subroutine  emis_flam_smold(n2, n3, plume, qsc_orig)

    implicit none
    ! Arguments:
    integer, intent(in) :: n2, n3 ! nsources, ngases
    real, intent(inout) :: plume(:,:,:) !(nsources, n2, n3)
    real, intent(in)    :: qsc_orig(:,:,:) !(ngases, n2, n3)

    ! Local Variables:
    real    :: frac_smold(n2,n3)
    integer :: k, kk, ksource, k_smold
    integer :: k_total

    !----------------------------------------------------------------------
    !               indexing to the array "plume(k,i,j)" 
    ! k 
    ! 1 => average area (m^2) for biome points in forest inside gribox i,j
    ! 2 => average area (m^2) for biome points in savana inside gribox i,j
    ! 3 => average area (m^2) for biome points in grass  inside gribox i,j
    ! 4 => Standard deviation for average area (m^2) for points in forest
    ! 5 => Standard deviation for average area (m^2) for points in savana
    ! 6 => Standard deviation for average area (m^2) for points in grass
    ! 7 a 9 => not used
    !10 => part of CO total emission corresponding to smoldering phase
    !11, 12 e 13 => Initialy stores the ratio between CO emission in 
    !               flaming phase and total emission, for forest, savana and
    !               grass respectivally
    !               qCO(flaming, forest) = plume(11,i,j)*qscCO_total.
    !               After this procedure this array stores the relation betwen:
    !               qCO(flaming, forest) and the amount of total emited for:
    !               qCO(flaming, floresta) =  plume(11,i,j)*plume(10,i,j)
    !               qCO(flaming, savana  ) =  plume(12,i,j)*plume(10,i,j)
    !               qCO(flaming, pastagem) =  plume(13,i,j)*plume(10,i,j)
    !20,21,22 e 23 the same for PM25               
    !
    !24-NZP => not used
    !----------------------------------------------------------------------
    !
    ! Calculate the smoldering emission and factors to obtain flaming fraction 
    ! as function of smoldering emission

    frac_smold(:,:) = 0.
  
    do k=1, 2 ! loop in burns sources for CO and PM25

       ksource = k
       ! assuming that qsc_orig(1,:,:) corresponding to CO, 
       ! and           qsc_orig(2,:,:) corresponding to PM25
       
       if (k==1) then
          k_smold = k_CO_smold
       elseif (k==2) then
          k_smold = k_PM25_smold
       endif
  
       ! Smoldering fraction
       frac_smold(:,:) = (1. - ( &
            plume(k_smold+1,:,:) + &
            plume(k_smold+2,:,:) + &
            plume(k_smold+3,:,:) ) )
  
       ! smoldering emission = total emission * smoldering fraction
       plume(k_smold,:,:) = qsc_orig(ksource,:,:) * frac_smold(:,:) ! kg/m^2/day
      
       ! conversion of flaming fraction to emission in smoldering phase
       do kk=1, 3
          plume(k_smold+kk,:,:) = plume(k_smold+kk,:,:)/(1.e-8+frac_smold(:,:))
       enddo

       ! saves the total emission value
       if(k==1) then
          k_total = k_CO_total   ! k =  9 to CO
       elseif(k==2) then
          k_total = k_PM25_total ! k = 19 to PM25
       endif
       plume(k_total,:,:) = qsc_orig(ksource,:,:)

   enddo

  end subroutine  emis_flam_smold

  !----------------------------------------------------------------

  ! Sub-routine that sets the latitudes and longitudes for input grid points
  subroutine set_latlon(nlon, nlat, prlon, prlat, dlon, dlat, lonni, latni)
    implicit none
    ! Arguments:
    integer, intent(in) :: nlon, nlat
    real, intent(in)    :: dlon, dlat, lonni, latni
    real, intent(out)   :: prlon(:,:), prlat(:,:)
    ! Local Variables:
    integer :: i, j

    do j=1,nlat
       do i=1,nlon
          prlon(i,j)=lonni+(i-1)*dlon
          prlat(i,j)=latni+(j-1)*dlat
       enddo
    enddo
  end subroutine set_latlon
  
  !----------------------------------------------------------------------

  subroutine plume2qsc(ng, ig, plume)

    implicit none

    integer,intent(IN) :: ng, ig !n1,n2,n3
    real,intent(IN)    :: plume(:,:,:) !(n1,n2,n3)
    
    ! OBS: QSC refers to emission source and IG refers to gas
    !      IG=1=2 (CO and COstc) uses the same array qsc(1,:,:)
 
    ! CO from burning with plumerise
    if(ig == 1) qsc(ng)%q(1,:,:) = plume(k_CO_smold,:,:) !plume_smold(1,:,:)
     
    ! CO from burning without plumerise
    if(ig == 2) qsc(ng)%q(1,:,:) = plume(k_CO_total,:,:)

    ! PM25 form burning with plumerise
    if(ig == 3) qsc(ng)%q(3,:,:) = plume(k_PM25_smold,:,:) !plume_smold(2,:,:)

  end subroutine plume2qsc

  !----------------------------------------------------------------------

  subroutine reorganize_sources(n2, n3, srcsc)
    
    implicit none

    ! Arguments
    integer,intent(IN)  :: n2,n3 !n1
    real, intent(inout) :: srcsc(:,:,:)

    ! Local variables
    integer :: i,j,ii,jj,k,kk
    real :: f,si,sf

    ! Perform a distribution in the sources to avoid strong located values
    ! (concentrated in 1 grid point)
    ! Distribution factor in 10% for each 18 neighbors
    ! (including the i,j central grid point)
    ! f=0.2
    ! f=0.1
    f=0.05 ! 5% for neighbors
    kk=1 ! First Level
  
    srcsc(2:,:,:)=0. ! all zeros, except k=1

!!$    do j=3,n3-2
!!$       do i=3,n2-2
    do j=1,n3
       do i=1,n2
          srcsc(2,i,j) = srcsc(2,i,j) + (1.-f)*srcsc(1,i,j)

          ! Distributes in the 18 grid points around i,j !  j+1  .   .   .
!!$          do jj = j-1,j+1             !   j   .   .   .  K=2
!!$             do ii = i-1,i+1                             !  j-1  .   .   .
          do jj = max(j-1,1),min(j+1,n3)
             !   j   .   .   .  K=2
             do ii = max(i-1,1),min(i+1,n2)
                !  j-1  .   .   .
                !      i-1  i  i+1
                do k=2,3
                   srcsc(k,ii,jj)=srcsc(k,ii,jj) + &
                        (1./18.) * f * srcsc(1,i,j)
                enddo
             enddo
          enddo
       enddo
    enddo

    ! Checks mass conservation

    si=0.
    sf=0.
!!$    do j=3,n3-2
!!$       do i=3,n2-2
    do j=1,n3
       do i=1,n2
          si   = si + srcsc(1,i,j)
          !	  print*,'si=',i,j,srcsc(1,i,j),si
          do k=2,3
             sf   = sf + srcsc(k,i,j)
             ! 	  print*,'sf=',k,srcsc(k,i,j),sf
          enddo

       enddo
    enddo

    if(si > 0.) then
       if ((abs(si-sf)/si) > 0.01) then          
          print*,'----------------------------------------------------'
          print*,'Mass conservation does not apply to REORGANIZE_SOURCES subroutine.'
          print*,'si e sf=',si,sf,si-sf,100*(si-sf)/si
          print*,'STOP'
          !	  stop 4455
          print*,'----------------------------------------------------'
       endif
    endif
    !   stop 442
  end subroutine reorganize_sources

  !-------------------------------------------------------------

  subroutine convert_to_misture_ratio(n2, n3, dn0, srcsc, rtgt, dzt)

    implicit none

    ! Arguments
    integer,intent(IN)  :: n2, n3
    real, intent(in)    :: dn0(:,:,:)   !n1,n2,n3
    real, intent(inout) :: srcsc(:,:,:) !n1,n2,n3
    real, intent(in)    :: rtgt(:,:)    !   n2,n3
    real, intent(in)    :: dzt(:)       !n1

    ! Local variables:
    ! Unit conversion factors:
    !!fcu=1.        !=> kg [gas/part] /kg [ar]
    !!fcu =1.e+12   !=> ng [gas/part] /kg [ar]
    real, parameter     :: fcu =1.e+6 !=> mg [gas/part] /kg [ar]
    real                :: dz
    integer             :: i,j,k

    do j=1,n3
       do i=1,n2
          do k=2,3

             ! All units in kg/m2/day => using 'dz' instead of 'vol'
             !    vol = 1./(dxt(i,j)*dyt(i,j)*dzt(k))*rtgt(i,j)
             dz = rtgt(i,j)/dzt(k) ! dzt=1/(z(k)-z(k-1))
         
             !=> units = mg[gas]/kg[ar]/day
	     srcsc(k,i,j)=fcu*srcsc(k,i,j)/(dz*dn0(k,i,j))

          enddo
       enddo
    enddo

    ! Doing source(k=1) = source(k=2)
    srcsc(1,:,:) = srcsc(2,:,:)

    ! Multipling qsc(kgas) by fcu to compatibilize units to output
    !do j=1,n3
    !   do i=1,n2
    !      srcsc(1,i,j) = srcsc(1,i,j)*fcu
    !   enddo
    !enddo

  end subroutine convert_to_misture_ratio

  !-----------------------------------------------------------------------

  subroutine burns(ng, mzp, mxp, myp, ia, iz, ja, jz, &
       scalar_g, time, iyear1, imonth1, idate1)

    use mem_scalar, only: scalar_vars ! type

    implicit none

    ! Arguments:
    integer, intent(in)              :: ng, mzp, mxp, myp, ia, iz, ja, jz
    type(scalar_vars), intent(inout) :: scalar_g(:,:)
    real, intent(in)                 :: time
    integer, intent(in)              :: iyear1, imonth1, idate1

    ! Local variables:
    integer :: npts
    real :: sclt1(mzp,mxp,myp), sclt2(mzp,mxp,myp), &
         sclt3(mzp,mxp,myp), sclt4(mzp,mxp,myp)
    
    if( plumerise == 1 ) &
         call plumerise_driver(mzp, mxp, myp, ia, iz, ja, jz, &
         k_CO_smold, k_PM25_smold)

    ! Coping sclt# information to a local array to call 'sources' and 'sink'
    npts = mzp*mxp*myp
    sclt1(:,:,:) = reshape(scalar_g(1,ng)%sclt(1:(npts)), (/mzp,mxp,myp/))
    sclt2(:,:,:) = reshape(scalar_g(2,ng)%sclt(1:(npts)), (/mzp,mxp,myp/))
    sclt3(:,:,:) = reshape(scalar_g(3,ng)%sclt(1:(npts)), (/mzp,mxp,myp/))
    sclt4(:,:,:) = reshape(scalar_g(4,ng)%sclt(1:(npts)), (/mzp,mxp,myp/))

    call sources(mzp, ia, iz, ja, jz, time, imonth1, idate1, iyear1, &
         sclt1,             &
         sclt2,             &
         sclt3,             &
         sclt4,             &
         scalar_g(1,ng)%srcsc, &
         scalar_g(2,ng)%srcsc, &
         scalar_g(3,ng)%srcsc, &
         scalar_g(4,ng)%srcsc  )

    call sink(mzp, ia, iz, ja, jz, &
         sclt1,  	&
         sclt2,  	&
         !!sclt3,  	&
         sclt4,  	&
         scalar_g(1,ng)%sclp, &
         scalar_g(2,ng)%sclp, &
         !!scalar_g(3,ng)%sclp, &
         scalar_g(4,ng)%sclp  )

    !Returning local values from sclt# to global scalar_g(#,ngrid)%sclt(1:(npts))
    scalar_g(1,ng)%sclt(1:(npts)) = reshape(sclt1(:,:,:), (/npts/))
    scalar_g(2,ng)%sclt(1:(npts)) = reshape(sclt2(:,:,:), (/npts/))
    scalar_g(3,ng)%sclt(1:(npts)) = reshape(sclt3(:,:,:), (/npts/))
    scalar_g(4,ng)%sclt(1:(npts)) = reshape(sclt4(:,:,:), (/npts/))

  end subroutine burns

  !------------------------------------------------------------

  subroutine sources(m1, ia, iz, ja, jz, time, imonth1, idate1, iyear1, &
       s1t, s2t, s3t, s4t, src1, src2, src3, src4)
    use ModDateUtils
    implicit none

    ! Arguments
    integer, intent(in) :: m1, ia, iz, ja, jz
    real, intent(in)    :: time
    integer, intent(in) :: imonth1, idate1, iyear1
    real, intent(inout) :: s1t(:,:,:) !(m1,m2,m3)
    real, intent(inout) :: s2t(:,:,:) !(m1,m2,m3)
    real, intent(inout) :: s3t(:,:,:) !(m1,m2,m3)
    real, intent(inout) :: s4t(:,:,:) !(m1,m2,m3)
    real, intent(in)    :: src1(:,:,:) !(m1,m2,m3)
    real, intent(in)    :: src2(:,:,:) !(m1,m2,m3)
    real, intent(in)    :: src3(:,:,:) !(m1,m2,m3)
    real, intent(in)    :: src4(:,:,:) !(m1,m2,m3)

    ! Local Variables:
!    integer, external :: julday
    integer :: iweek,idays,j,i,k
    real :: tign,strtim,ax,bx,cx,timeq,rinti,r_q,r_antro
    real, dimension(7) :: week_CYCLE
    !             weak day     :  DOM   SEG   TER   QUA   QUI   SEX   SAB
    !             iweak        :   1     2     3     4     5     6     7
    data (week_CYCLE(iweek),iweek=1,7) /0.25, 1.25, 1.25, 1.25, 1.25, 1.25, 0.50/

    !-------------Emission of biomass buning--------------------
    !Time
    tign = 0. * 3600. !UTC time of ignition
    strtim=0.
    ! number of days of simulation
    idays = int((strtim+time/3600. - tign/3600. )/24.+.00001)
    tign = tign + real(idays)*24.*3600.

    ! Burning mean modulation during day cycle (1/s)
    ! with int(r_q dt) (0 - 24h)= 1.
    ax   = 2000.6038
    ! bx   = 15.041288 * 3600. ! Local time peak
    bx   = 18.041288 * 3600.  ! Peak in 18 UTC
    cx   =  2.184936 * 3600.
    timeq= ( time - tign ) - bx
    rinti=2.1813936e-8
    r_q  = rinti*( ax * exp( -timeq**2/(2.*cx**2) ) + 100. -  &
         5.6712963e-4*( time - tign ))

    !-------------Antropogenic Emission (industrial, home, ...)
    ! Weak cycle
    ! Weak day
    iweek= int(((float(julday(imonth1,idate1,iyear1))/7. - &
         int(julday(imonth1,idate1,iyear1)/7))*7.)) + 3
    if(iweek>7) iweek = iweek-7
    ! Day cycle
    bx   = 15.041288 * 3600. ! Peak in 15 UTC
    timeq= ( time - tign ) - bx
    r_antro  =1.4041297e-05*(exp(-((timeq)**2)/(   43200.    **2))+0.1)

    ! Day cycle + weak
    r_antro = r_antro * week_CYCLE(iweek)
    
    do j=ja,jz
       do i=ia,iz
          do k=2,m1-1
             ! All units are in kg[n]/(kg[ar] s

             ! Burnign emission use r_q
             ! CO with plume rise
             s1t(k,i,j) = s1t(k,i,j) + src1(k,i,j)*r_q
	     	     
	     ! CO without plume rise 
             s2t(k,i,j) = s2t(k,i,j) + src2(k,i,j)*r_q

	     ! PM25 with plume rise
             s3t(k,i,j) = s3t(k,i,j) + src3(k,i,j)*r_q

             ! Antropogenic emission using emission cycle r_antro
             s4t(k,i,j) = s4t(k,i,j) + src4(k,i,j)*r_antro

          enddo
       enddo
    enddo

  end subroutine sources

  !------------------------------------------------------------

  subroutine sink(m1, ia, iz, ja, jz, s1t, s2t, s4t, s1p, s2p, s4p)
    implicit none

    ! Arguments:
    integer,intent(in)  :: m1,ia,iz,ja,jz
    real, intent(inout) :: s1t(:,:,:) !(m1,m2,m3)
    real, intent(inout) :: s2t(:,:,:) !(m1,m2,m3)
    !!real, intent(inout) :: s3t(:,:,:) !(m1,m2,m3)
    real, intent(inout) :: s4t(:,:,:) !(m1,m2,m3)
    real, intent(in)    :: s1p(:,:,:) !(m1,m2,m3)
    real, intent(in)    :: s2p(:,:,:) !(m1,m2,m3)
    !!real, intent(in)    :: s3p(:,:,:) !(m1,m2,m3)
    real, intent(in)    :: s4p(:,:,:) !(m1,m2,m3)

    ! Local variables:
    ! sink for CO and PM25 mean life
    real, parameter :: vmco = 30.*24.*3600., vmpm25 = 6.*24.*3600.
    integer :: j,i,k

    do j = ja,jz
       do i = ia,iz
          do k = 2,m1-1
             s1t(k,i,j) = s1t(k,i,j) - s1p(k,i,j)/vmco
             s2t(k,i,j) = s2t(k,i,j) - s2p(k,i,j)/vmco
             !srf- using param dry deposition
             !s3t(k,i,j) = s3t(k,i,j) - s3p(k,i,j)/vmpm25
             s4t(k,i,j) = s4t(k,i,j) - s4p(k,i,j)/vmco

             !do ig=1,6
             !if( s1t(k,i,j) > 0.) then
             !    print*,'kgas kij 1=',k,i,j,s1p(k,i,j)
             !    print*,'kgas kij 2=',k,i,j,s2p(k,i,j)
             !    print*,'kgas kij 3=',k,i,j,s3p(k,i,j)
             !    print*,'kgas kij 4=',k,i,j,s4p(k,i,j)
             !    print*,'kgas kij 5=',k,i,j,s5p(k,i,j)
             !    print*,'kgas kij 6=',k,i,j,s6p(k,i,j)
             !endif

          enddo
       enddo
    enddo

  end subroutine sink


  subroutine StoreNamelistFileAtEmission_source_map(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    firemapfn = oneNamelistFile%firemapfn
    plumerise = oneNamelistFile%plumerise
    define_proc = oneNamelistFile%define_proc
  end subroutine StoreNamelistFileAtEmission_source_map


end module emission_source_map


!--(DMK-BRAMS-5.0-FIM)---------------------------------------------------------
!!$!-----------------------------------------------------------------------
!!$
!!$subroutine latset_tracer(m1, m2, m3, ia, iz, ja, jz, ibcon, ap, uc, vc, &
!!$     dxu, dxm, dyv, dym)
!!$
!!$  use mem_grid, only: dtlt, ibnd, jdim, jbnd, nz
!!$  use mem_scratch, only: vctr17, vctr18
!!$
!!$  implicit none
!!$
!!$  ! Arguments:
!!$  integer, intent(in) :: m1,m2,m3,ia,iz,ja,jz,ibcon
!!$  real, intent(inout) :: ap(m1,m2,m3)
!!$  real, intent(in)    :: uc(m1,m2,m3), vc(m1,m2,m3)
!!$  real, intent(in)    :: dxu(m2,m3), dxm(m2,m3), dyv(m2,m3), dym(m2,m3)
!!$
!!$  ! Local variables:
!!$  integer :: i,j,k,lbw,lbe,lbs,lbn
!!$  real :: thresh,dtlx,c1,dxr,dyr
!!$
!!$  if (iand(ibcon,1) .gt. 0) lbw = ia - 1
!!$  if (iand(ibcon,2) .gt. 0) lbe = iz + 1
!!$  if (iand(ibcon,4) .gt. 0) lbs = ja - 1
!!$  if (iand(ibcon,8) .gt. 0) lbn = jz + 1
!!$
!!$  thresh = 0.
!!$  dtlx = dtlt
!!$
!!$  if (ibnd .ne. 4) then
!!$
!!$     ! Western boundary for lsflg = 2  == constant inflow, radiative b.c. outflow
!!$     ! Notice that the field ap(k,i,j) in the boundary is modified only in case
!!$     ! of outflow. So in inflow regions the values are constant and equal to
!!$     ! initial value (t=0).
!!$     ! To introduce some kind of nudging in boundary condition, change the array
!!$     ! ap(k,lbw,j) to western boundary. The same procedure is required to others
!!$     ! boundarys
!!$
!!$     if (iand(ibcon,1) .gt. 0) then
!!$        do j = 1,m3
!!$           dxr = dxu(ia,j) / dxu(lbw,j)
!!$           c1 = dtlx * dxu(lbw,j)
!!$           do k = 1,m1
!!$              vctr17(k) = -c1 * uc(k,lbw,j)
!!$              vctr18(k) = ap(k,ia,j) + dxr * (ap(k,ia,j) - ap(k,ia+1,j))
!!$           enddo
!!$           do k = 1,m1
!!$              if (vctr17(k) .ge. thresh) ap(k,lbw,j) = vctr18(k)
!!$           enddo
!!$        enddo
!!$     endif
!!$
!!$     !     Eastern Boundary for LSFLG =  2
!!$
!!$     if (iand(ibcon,2) .gt. 0) then
!!$        do j = 1,m3
!!$           dxr = dxu(iz-1,j) / dxu(iz,j)
!!$           c1 = dtlx * dxu(iz,j)
!!$           do k = 1,m1
!!$              vctr17(k) = c1 * uc(k,iz,j)
!!$              vctr18(k) = ap(k,iz,j) + dxr * (ap(k,iz,j) - ap(k,iz-1,j))
!!$           enddo
!!$           do k = 1,m1
!!$              if (vctr17(k) .ge. thresh)  ap(k,lbe,j) = vctr18(k)
!!$           enddo
!!$        enddo
!!$     endif
!!$  endif
!!$
!!$  if(jdim.eq.1.and.jbnd.ne.4)then
!!$
!!$     !     Southern boundary for LSFLG  2
!!$
!!$     if (iand(ibcon,4) .gt. 0) then
!!$        do i = 1,m2
!!$           dyr = dyv(i,ja) / dyv(i,lbs)
!!$           c1 = dtlx * dyv(i,lbs)
!!$           do k = 1,nz
!!$              vctr17(k) = -c1 * vc(k,i,lbs)
!!$              vctr18(k) = ap(k,i,ja) + dyr * (ap(k,i,ja) - ap(k,i,ja+1))
!!$           enddo
!!$           do k = 1,m1
!!$              if (vctr17(k) .ge. thresh) ap(k,i,lbs) = vctr18(k)
!!$           enddo
!!$        enddo
!!$     endif
!!$
!!$     !     Northern Boundary for LSFLG =  2
!!$
!!$     if (iand(ibcon,8) .gt. 0) then
!!$        do i = 1,m2
!!$           dyr = dyv(i,jz-1) / dyv(i,jz)
!!$           c1 = dtlx * dyv(i,jz)
!!$           do k = 1,m1
!!$              vctr17(k) = c1 * vc(k,i,jz)
!!$              vctr18(k) = ap(k,i,jz) + dyr * (ap(k,i,jz) - ap(k,i,jz-1))
!!$           enddo
!!$           do k = 1,m1
!!$              if (vctr17(k) .ge. thresh) ap(k,i,lbn) = vctr18(k)
!!$           enddo
!!$        enddo
!!$     endif
!!$  endif
!!$end subroutine latset_tracer
!--(DMK-BRAMS-5.0-FIM)---------------------------------------------------------
