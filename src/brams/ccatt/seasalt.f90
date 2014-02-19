MODULE seaSalt

   INTEGER,PARAMETER :: nbeans=3
   real,parameter,dimension(nbeans,0:4) :: c=RESHAPE((/   &
      -2.881e+06,-6.743e+06, 2.181e+06,   &
      -3.003e+13, 1.183e+14,-4.165e+12,   &
      -2.867e+21,-8.148e+20, 3.132e+18,   &
       5.932e+28, 2.404e+27,-9.841e+23,   &
      -2.576e+35,-2.452e+33, 1.085e+29/), &
       (/3,5/))
   real,parameter,dimension(nbeans,0:4) :: d=RESHAPE((/   &
       7.609e+08, 2.279e+09,-5.800e+08,   &
       1.829e+16,-3.787e+16, 1.105e+15,   &
       6.791e+23, 2.528e+23,-8.297e+20,   &
      -1.616e+31,-7.310e+29, 2.601e+26,   &
       7.188e+37, 7.368e+35,-2.859e+31/), &
       (/3,5/))
   REAL, PARAMETER,DIMENSION(nBeans) :: minDiam=(/0.020e-06,0.145e-06,0.419e-06/)
   REAL, PARAMETER,DIMENSION(nBeans) :: maxDiam=(/0.145e-06,0.419e-06,2.800e-06/)
   REAL, PARAMETER,DIMENSION(nBeans) :: Diam=(minDiam+maxDiam)/2
   !      ak(k)=c(k,4)*diam(k)**4+ &
   !            c(k,3)*diam(k)**3+ &
   !            c(k,2)*diam(k)**2+ &
   !            c(k,1)*diam(k)   + &
   !            c(k,0)
   !      bk(k)=d(k,4)*diam(k)**4+ &
   !            d(k,3)*diam(k)**3+ &
   !            d(k,2)*diam(k)**2+ &
   !            d(k,1)*diam(k)   + &
   !            d(k,0)
   REAL, PARAMETER,DIMENSION(nBeans) :: ak=(/-3.496e+06, 2.264e+05, 1.433e+06/)
   REAL, PARAMETER,DIMENSION(nBeans) :: bk=(/ 1.148e+09,-3.034e+07,-3.815e+08/)
   REAL, PARAMETER :: DrySeaSaltDendity=2.17 !g/cm3
   REAL, PARAMETER :: pi=3.1415926
   REAL, PARAMETER :: mp1=(4/3)*Pi*(diam(1)/2)**3*DrySeaSaltDendity*1.0e-06
   REAL, PARAMETER :: mp2=(4/3)*Pi*(diam(2)/2)**3*DrySeaSaltDendity*1.0e-06
   REAL, PARAMETER :: mp3=(4/3)*Pi*(diam(3)/2)**3*DrySeaSaltDendity*1.0e-06
    
   REAL, PARAMETER,DIMENSION(nBeans) :: massSSaer=(/mp1,mp2,mp3/)
   LOGICAL, PARAMETER :: dump=.false. ! WARNING: Must be used just when running with one thread
   
contains


   subroutine aersaltFlux(flux,w10,sst)
      implicit none
      
      real, intent(in) :: w10
      real, intent(in) :: sst
      real, intent(out) :: flux(3)
      INTEGER :: k
   
      DO k=1,nBeans

         flux(k)=3.84e-06*(ak(k)*sst+bk(k))*w10**3.41
         !print *, massSSaer(k)
      END DO
      
   end subroutine aersaltFlux


   Subroutine massSeaSalt(ia,iz,ja,jz,w10,sst,landMask,seaSaltMass)
   
      INTEGER, INTENT(IN) :: ia,iz,ja,jz           ! Dimensions
      REAL,INTENT(IN) :: w10(ia:iz,ja:jz)          ! wind at 10m [m/s]
      REAL,INTENT(IN) :: sst(ia:iz,ja:jz)          ! Sea Surface Temperature [k]
      LOGICAL, INTENT(IN) :: landMask(ia:iz,ja:jz) !Mask land surface (.false.=ocean)
      
      REAL,intent(out) :: seaSaltMass(nBeans,ia:iz,ja:jz) !Mass flux of SeaSalt for each grid cell [kg/m2s]
      
      INTEGER :: i,j,k
      real,dimension(nBeans) :: flux
      
      seaSaltMass=0.0
      
      do i=ia,iz
         do j=ja,jz
            !WRITE(*,FMT='("LFR: ",2(I2.2,1X),L1)'), i,j,landmask(i,j)
            if(landMask(i,j)) CYCLE
            CALL aersaltFlux(flux,w10(i,j),sst(i,j))
            DO k=1,nBeans
               seaSaltMass(k,i,j)=flux(k)*massSSaer(k)*1.0e-3 !Kg/m2S 
            END DO
         END DO
      END DO
   
   end subroutine massSeaSalt

end module seaSalt

subroutine seaSaltDriver(ia,iz,ja,jz,ngrid,m2,m3)
   USE mem_leaf, ONLY: &
         leaf_g           ! INTENT(IN)
   USE mem_grid, ONLY: &
         npatch,       &  ! INTENT(IN)
         deltaxn,      &  ! INTENT(IN)
         deltayn,      &  ! INTENT(IN)
         time,         &
         grid_g,       &
         dzt
!   use leaf_coms, only: &
!         timefac_sst      ! INTENT(IN)
   use io_params, only: &  
         iupdsst,       & ! INTENT(IN)
         ssttime1,      & ! INTENT(IN)
         ssttime2         ! INTENT(IN)
  USE mem_basic, ONLY:  &
         basic_g  !INTENT(IN) 
  USE seaSalt, ONLY:  &
         massSeaSalt, &
         dump
  use node_mod, only: &
       mynum
  USE mem_aer1, ONLY: &
       aer1_g
  USE aer1_list, ONLY: &
      marin
      
   INTEGER, INTENT(IN) :: ia,iz,ja,jz,ngrid,m2,m3
   
   REAL :: seaSaltMass(3,ia:iz,ja:jz)
   
   INTEGER :: np,i,j,nmd
   LOGICAL :: landMask(ia:iz,ja:jz) ! .false. is ocean
   REAL :: sst(ia:iz,ja:jz)
   REAL :: w10(ia:iz,ja:jz)
   REAL :: timefac_sst
   CHARACTER(LEN=8) :: ctime
   CHARACTER(LEN=6) :: clat,clon
   
   DO i=ia,iz
      DO j=ja,jz
         IF(leaf_g(ngrid)%patch_area(i,j,1)==1) THEN
            landMask(i,j)=.false.
         ELSE
            landMask(i,j)=.true.
         END IF
         !WRITE(*,FMT='("LFR: ",2(I2.2,1X),2(F6.2,1X),L1)'), i,j,grid_g(ngrid)%glat(i,j),grid_g(ngrid)%glon(i,j),landmask(i,j)
         ! Time interpolation factor for updating SST
         if (iupdsst == 0) then
            timefac_sst = 0.
         else
            timefac_sst = (time - ssttime1(ngrid)) / (ssttime2(ngrid) - ssttime1(ngrid))
         endif
         sst(i,j) =   leaf_g(ngrid)%seatp(i,j) + &
                    ( leaf_g(ngrid)%seatf(i,j) - leaf_g(ngrid)%seatp(i,j) ) * timefac_sst
         w10(i,j)=(sqrt(basic_g(ngrid)%up(1,i,j)**2+basic_g(ngrid)%vp(1,i,j)**2)+ &
                   sqrt(basic_g(ngrid)%up(2,i,j)**2+basic_g(ngrid)%vp(2,i,j)**2))/2
      END DO
   END DO

   CALL massSeaSalt(ia,iz,ja,jz,w10,sst,landMask,seaSaltMass)
  
   DO i=ia,iz
      DO j=ja,jz
         DO nmd=1,3
            dz = grid_g(ngrid)%rtgt(i,j)/dzt(2)
            aer1_g(nmd,marin,ngrid)%sc_src(1,i,j)=seaSaltMass(nmd,i,j)/dz!/basic_g(ngrid)%dn0(2,i,j)/dz
         END DO
      END DO
   END DO
   
   !If dump is ON the code will create a grads file with output of sea salt aerosols
   !WARNING: Use it only if the model run with one processor (thread).
   IF(dump) CALL writeSeaSalt(m2,m3,ia,iz,ja,jz,              &
                     aer1_g(1,marin,ngrid)%sc_src(2,:,:),     &
                     aer1_g(2,marin,ngrid)%sc_src(2,:,:),     &
                     aer1_g(3,marin,ngrid)%sc_src(2,:,:),     &
                     sst,                                     &
                     w10,                                     &
                     landMask,                                &
                     grid_g(ngrid)%glat(1,1),                 &
                     grid_g(ngrid)%glon(1,1),                 &
                     time                                     &
                     )
   !STOP 777
end subroutine seaSaltDriver

SUBROUTINE writeSeaSalt(nx,ny,ia,iz,ja,jz,var1,var2,var3,sst,w10,mask,lat,lon,time)
   use mem_grid, ONLY: dtlt
   
   INTEGER, INTENT(IN) :: nx,ny
   real,INTENT(IN) :: var1(nx,ny), var2(nx,ny),var3(nx,ny),lat,lon,time
   REAL,INTENT(IN) :: sst(ia:iz,ja:jz),w10(ia:iz,ja:jz)
   LOGICAL, INTENT(IN) :: mask(ia:iz,ja:jz)
   
   integer,parameter :: nvar =6 
   INTEGER,save :: nt = 1
   integer,SAVE :: it=0
   ! loop counter
   INTEGER,SAVE :: irec=0  ! Record number
   REAL, PARAMETER :: timewrite=1200.0
   REAL :: var4(nx,ny),var5(nx,ny),var6(nx,ny)


   IF (.not. (mod(time + .001,timewrite) .lt. dtlt .or. time .lt. 0.001)) RETURN
   
   var4=0.0;var5=0.0;var6=0.0
   
   DO i=ia,iz
      Do j=ja,jz
         var4(i,j)=sst(i,j)
         var5(i,j)=w10(i,j)
         IF(.not. mask(i,j)) var6(i,j)=1.0
      END DO
   END DO
   

   IF(it==0) THEN
      OPEN(30,file='Outfile.gra',action='write',form='unformatted',access='direct',&
           recl=4*nVar*nx*ny,status='replace')
      CALL writeSeaSaltCtl(nx,ny,lat,lon)
   END IF
   
   !loop over timesteps
   it  = it+1
   irec = irec + 1
   WRITE(30,rec=irec) Var1(:,:),Var2(:,:),Var3(:,:),Var4(:,:),Var5(:,:),Var6(:,:)
   
END SUBROUTINE writeSeaSalt

SUBROUTINE writeSeaSaltCtl(nx,ny,lat,lon)
   INTEGER, INTENT(IN) :: nx,ny
   REAL, INTENT(IN) :: lat,lon
   OPEN(31,FILE='Outfile.ctl',STATUS='replace')
   write (31,fmt='(A)') 'dset ./Outfile.gra'
   write (31,fmt='(A)') 'undef  -0.9990000E+34'
   write (31,fmt='(A)') 'title Sea Salt Aerosol'
   write (31,fmt='(A)') 'options little_endian'
   write (31,fmt='(A,I2.2,A,F8.4,1X,F8.4)') 'xdef ',nx,' linear ',lon,0.95
   write (31,fmt='(A,I2.2,A,F8.4,1X,F8.4)') 'ydef ',ny,' linear ',lat,0.89
   write (31,fmt='(A)') 'zdef 1 linear 1.0 1.0'
   write (31,fmt='(A)') 'tdef 5 linear 00:00Z01jan2012 1hr'
   write (31,fmt='(A)') 'vars 6'
   write (31,fmt='(A)') 'seas1       0 99    - Sea Salt 0.020e-06 a 0.145e-06  [ ]'
   write (31,fmt='(A)') 'seas2       0 99    - Sea Salt 0.145e-06 a 0.419e-06  [ ]'
   write (31,fmt='(A)') 'seas3       0 99    - Sea Salt 0.419e-06 a 2.800e-06  [ ]'
   write (31,fmt='(A)') 'sst         0 99    - Sea Surface Temperature         [K]'
   write (31,fmt='(A)') 'w10         0 99    - Wind 10m                        [m/s]'
   write (31,fmt='(A)') 'mask        0 99    - Mask land Ocean=1.0             [ ]'   
   write (31,fmt='(A)') 'endvars'
END SUBROUTINE writeSeaSaltCtl
