!------------------------------------------------------------------------------------------!
! Subroutine prep_advflx_to_stilt                                                          !
! Developed by Saulo R. Freitas (CPTEC/INPE)                                               !
!                                                                                          !
!   This subroutine prepares the advective fluxes to be used in STILT (or other            !
! Lagrangian models).                                                                      !
!------------------------------------------------------------------------------------------!
subroutine prep_advflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz,ng)
   
   use mem_grid   ,  only: dtlt
   use mem_scratch,  only: scratch
   use mem_stilt   ,  only: stilt_g                 & ! structure
                         , frqmassave             & ! intent(in)
                         , etime_adve             & ! intent(inout)
                         , zero_average_mass_adve ! ! subroutine

   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: mzp
   integer, intent(in) :: mxp
   integer, intent(in) :: myp
   integer, intent(in) :: ia
   integer, intent(in) :: iz
   integer, intent(in) :: ja
   integer, intent(in) :: jz
   integer, intent(in) :: ng
   !----- Local variables. ----------------------------------------------------------------!
   real                :: dtlti
   real                :: frqmassi
   !---------------------------------------------------------------------------------------!


   !----- Find the inverse of the time step and the average time. -------------------------!
   dtlti    = 1./dtlt
   frqmassi = 1./frqmassave
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Update the step counter for advective fluxes.  If this exceeds the maximum number !
   ! of steps, it means it is time to reset it and start integrating again.                ! 
   !---------------------------------------------------------------------------------------!
   etime_adve(ng) = etime_adve(ng) + dtlt
   if (etime_adve(ng) > frqmassave + 0.1 * dtlt) then
      etime_adve(ng) = dtlt
      call zero_average_mass_adve(stilt_g(ng))
   end if
   !---------------------------------------------------------------------------------------!

   call compute_mass_flux(mzp,mxp,myp,ia,iz,ja,jz,dtlti,frqmassi                           &
           ,scratch%vt3da            , scratch%vt3db           , scratch%vt3dc             &
           ,stilt_g(ng)%afxu          , stilt_g(ng)%afxv         , stilt_g(ng)%afxw           &
           ,stilt_g(ng)%afxub         , stilt_g(ng)%afxvb        , stilt_g(ng)%afxwb          )

   return

end subroutine prep_advflx_to_stilt
!------------------------------------------------------------------------------------------!

!==========================================================================================!
!==========================================================================================!
! Subroutine compute_mass_flux                                                             !
! Based on original Saulo R. Freitas (CPTEC/INPE) subroutine                               !
!                                                                                          !
! This subroutine compute the integrated mass flux from advection.                         !
!------------------------------------------------------------------------------------------!
subroutine compute_mass_flux(mzp,mxp,myp,ia,iz,ja,jz,dtlti,frqmassi,vt3da,vt3db,vt3dc,afxu &
                            ,afxv,afxw,afxub,afxvb,afxwb)
   implicit none
   integer                        , intent(in)    :: mzp,mxp,myp
   integer                        , intent(in)    :: ia,iz,ja,jz
   real                           , intent(in)    :: dtlti,frqmassi
   real   , dimension(mzp,mxp,myp), intent(in)    :: vt3da,vt3db,vt3dc
   real   , dimension(mzp,mxp,myp), intent(out)   :: afxu,afxv,afxw
   real   , dimension(mzp,mxp,myp), intent(inout) :: afxub,afxvb,afxwb
   integer                                        :: i,j,k

   do k=1,mzp
      do i=ia,iz
         do j=ja,jz
            afxu(k,i,j)  =                vt3da(k,i,j) * dtlti
            afxv(k,i,j)  =                vt3db(k,i,j) * dtlti
            afxw(k,i,j)  =                vt3dc(k,i,j) * dtlti
            afxub(k,i,j) = afxub(k,i,j) + vt3da(k,i,j) * frqmassi
            afxvb(k,i,j) = afxvb(k,i,j) + vt3db(k,i,j) * frqmassi
            afxwb(k,i,j) = afxwb(k,i,j) + vt3dc(k,i,j) * frqmassi
         end do
      end do
   end do

   return
end subroutine compute_mass_flux
!==========================================================================================!
!==========================================================================================!

!------------------------------------------------------------------------------------------!
! Subroutine prep_convflx_to_stilt                                                         !
! Developed by Saulo R. Freitas (CPTEC/INPE)                                               !
!                                                                                          !
!   This subroutine prepares the convective fluxes to be used in STILT (or other           !
! Lagrangian models).                                                                      !
!------------------------------------------------------------------------------------------!
subroutine prep_convflx_to_stilt(m1,m2,m3,ia,iz,ja,jz,mgmxp,mgmyp,mgmzp,maxiens,ngrid      &
                                ,ngrids_cp,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d       &
                                ,kpbl4d,kstabi4d,kstabm4d,xmb4d,edt4d,zcup5d,pcup5d,enup5d &
                                ,endn5d,deup5d,dedn5d,zup5d,zdn5d,iens)

use mem_stilt

implicit none

integer, intent(in) :: mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,iens, maxiens,m1,m2,m3,ia,iz,ja,jz

integer, intent(in),dimension(mgmxp,mgmyp,maxiens,ngrids_cp) ::                            &
                         ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d,kstabi4d,kstabm4d
                       
real, intent(in), dimension(mgmxp,mgmyp,maxiens,ngrids_cp) :: xmb4d,edt4d
                        
real, intent(in), dimension(mgmzp,mgmxp,mgmyp,maxiens,ngrids_cp) ::                        &
                                      enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d,zcup5d,pcup5d

integer :: i, j


!zero out all convflx
if (iens == 1) then
!Deep conv

!--(DMK-BRAMS-5.0)---------------------------------------------
!  call azero(m1*m2*m3,stilt_g(ngrid)%cfxup1(1,1,1))
!  call azero(m1*m2*m3,stilt_g(ngrid)%cfxdn1(1,1,1))
!  call azero(m1*m2*m3,stilt_g(ngrid)%dfxup1(1,1,1))
!  call azero(m1*m2*m3,stilt_g(ngrid)%efxup1(1,1,1))
!  call azero(m1*m2*m3,stilt_g(ngrid)%dfxdn1(1,1,1))
!  call azero(m1*m2*m3,stilt_g(ngrid)%efxdn1(1,1,1))

  stilt_g(ngrid)%cfxup1 = 0.
  stilt_g(ngrid)%cfxdn1 = 0.
  stilt_g(ngrid)%dfxup1 = 0.
  stilt_g(ngrid)%efxup1 = 0.
  stilt_g(ngrid)%dfxdn1 = 0.
  stilt_g(ngrid)%efxdn1 = 0.
!--(DMK-BRAMS-5.0)---------------------------------------------

!shallow conv
elseif (iens == 2) then

!--(DMK-BRAMS-5.0)---------------------------------------------
! call azero(m1*m2*m3,stilt_g(ngrid)%cfxup2(1,1,1))
! call azero(m1*m2*m3,stilt_g(ngrid)%dfxup2(1,1,1))       
! call azero(m1*m2*m3,stilt_g(ngrid)%efxup2(1,1,1))

 stilt_g(ngrid)%cfxup2 = 0.
 stilt_g(ngrid)%dfxup2 = 0.
 stilt_g(ngrid)%efxup2 = 0.
!--(DMK-BRAMS-5.0)---------------------------------------------

end if

do j=ja,jz
  do i=ia,iz
!    if((iens == 1 .and. ierr4d(i,j,iens,ngrid) == 0) .or. iens == 2) then
    if(ierr4d(i,j,iens,ngrid) == 0) then
      call get_convflx(iens,i,j,mgmzp,m1,m2,m3                                             &
                  ,   xmb4d(i,j,iens,ngrid),   edt4d(i,j,iens,ngrid)                       &
                  ,  jmin4d(i,j,iens,ngrid),  kdet4d(i,j,iens,ngrid)                       &
                  ,   k224d(i,j,iens,ngrid), kbcon4d(i,j,iens,ngrid)                       &
                  ,  ktop4d(i,j,iens,ngrid),  kpbl4d(i,j,iens,ngrid)                       &
                  ,kstabi4d(i,j,iens,ngrid),kstabm4d(i,j,iens,ngrid)                       &
                  ,zcup5d(1:mgmzp,i,j,iens,ngrid),pcup5d(1:mgmzp,i,j,iens,ngrid)           &
                  ,deup5d(1:mgmzp,i,j,iens,ngrid),enup5d(1:mgmzp,i,j,iens,ngrid)           &
                  ,dedn5d(1:mgmzp,i,j,iens,ngrid),endn5d(1:mgmzp,i,j,iens,ngrid)           &
                  , zup5d(1:mgmzp,i,j,iens,ngrid), zdn5d(1:mgmzp,i,j,iens,ngrid)           &
                  ,stilt_g(ngrid)%cfxup1(1,1,1),stilt_g(ngrid)%cfxdn1(1,1,1)               &
                  ,stilt_g(ngrid)%dfxup1(1,1,1),stilt_g(ngrid)%efxup1(1,1,1)               &
                  ,stilt_g(ngrid)%dfxdn1(1,1,1),stilt_g(ngrid)%efxdn1(1,1,1)               &
                  ,stilt_g(ngrid)%cfxup2(1,1,1),stilt_g(ngrid)%dfxup2(1,1,1)               &
                  ,stilt_g(ngrid)%efxup2(1,1,1))
    endif 
  enddo 
enddo  

return
end subroutine prep_convflx_to_stilt
!------------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------------!
! Subroutine get_convflx                                                                   !
! Developed by Saulo R. Freitas (CPTEC/INPE)                                               !
!                                                                                          !
!   This subroutine aims at getting the convective fluxes from Grell shallow and           !
! deep convective parameterizations.                                                       !
!------------------------------------------------------------------------------------------!
subroutine get_convflx(iens,i,j,mgmzp,m1,m2,m3,xmb,edt,jmin,kdet,k22,kbcon,ktop,kpbl       &
                      ,kstabi,kstabm,z_cup,p_cup,cd,entr,cdd,entrd,zu,zd,cfxup1,cfxdn1     &
                      ,dfxup1,efxup1,dfxdn1,efxdn1,cfxup2,dfxup2,efxup2)

implicit none

integer, intent(in)                      ::  iens,i,j,mgmzp,m1,m2,m3,jmin,kdet,k22,kbcon   &
                                            ,ktop,kpbl,kstabi,kstabm
                                            
real, intent(in)                         ::  xmb,edt

real, intent(in),    dimension(mgmzp)    ::  z_cup,p_cup,cd,entr,cdd,entrd,zu,zd

real, intent(inout), dimension(m1,m2,m3) ::  cfxup1,cfxdn1,dfxup1,efxup1,dfxdn1,efxdn1     &
			                    ,cfxup2,dfxup2,efxup2  
                                             
integer                                  ::  k,kr
real                                     ::  dz,totmas,entup,detup,entdoj,entupk,detupk    &
                                            ,detdo,entdo,subdown,subin,detdo1,detdo2

do k=2,ktop+1
        
  kr= k + 1   ! level K of conv grid  corresponds to level K + 1 of RAMS grid
  dz =  z_cup(kr) - z_cup(k)
   
  entup   = 0.
  detup   = 0.
  entdoj  = 0.
  entupk  = 0.
  detupk  = 0.
  detdo   = edt*cdd(k)*dz*zd(kr)
  entdo   = edt*entrd(k)*dz*zd(kr)
  subdown = zu(k ) - edt*zd(k )
  subin   = zu(kr) - edt*zd(kr)

  if (k >= kbcon .and. k < ktop) then
    entup  = entr(k)   *dz*zu(k)
    detup  =   cd(kr) *dz*zu(k)
  end if

  if(k == jmin)  entdoj = edt*zd(k)
  if(k == k22-1) entupk = zu(kpbl)
  if(k == ktop)  detupk = zu(ktop)
  if(k > kdet)   detdo  = 0.
  if(k == ktop)  subin  = 0.
  if(k < kbcon)  detup  = 0.
  
  if(iens == 1) then ! Deep convection
      cfxup1(k ,i,j) =     xmb* zu(k)
      cfxdn1(k ,i,j) =-edt*xmb* zd(k)
      dfxup1(kr,i,j) =     xmb*(detup + detupk)
      efxup1(kr,i,j) =    -xmb*(entup + entupk)
      dfxdn1(kr,i,j) =     xmb*(detdo         ) !edt already is at detdo
      efxdn1(kr,i,j) =    -xmb*(entdo + entdoj) !edt already is at entdo,entdoj
  elseif(iens == 2)  then ! Shallow convection
      cfxup2(k ,i,j) =     xmb* zu(k)
      dfxup2(kr,i,j) =     xmb*(detup + detupk)
      efxup2(kr,i,j) =    -xmb*(entup + entupk)
  end if
!------------------------------------------------------------------------------------------!
! Checking the mass conservation                                                           !
!------------------------------------------------------------------------------------------!
  totmas=subin-subdown+detup-entup-entdo+detdo-entupk-entdoj+detupk
  if(abs(totmas) > 1.e-6) then
    write (unit=*,fmt='(a)')                 '----------- Subroutine Get_convflx ----------'
    write(unit=*, fmt='(4(a,1x,i3,1x))')     '  K= ',k,'   I=',i,'   J=',j,'   IENS=',iens
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  subin=  ',    subin,'subdown=',subdown
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  detup=  ',    detup,'entup=  ',entup
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entdo=  ',    entdo,'detdo=  ',detdo
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entupk= ',   entupk,'detupk= ',detupk
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  zu(k)=  ',    zu(k),'zd(k)=  ',zd(k)
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  zu(kr)= ',   zu(kr),'zd(kr)= ',zd(kr)
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  entdoj= ',   entdoj,'edt=    ',edt
    write(unit=*, fmt='(1(a,1x,es10.3,1x))') '  totmas= ',   totmas
    write(unit=*, fmt='(a)')                 '---------------------------------------------'
    stop 'The model will stop since it is not conserving mass...' 
  end if

end do

!------------------------------------------------------------------------------------------!
! Bottom layer                                                                             !
!------------------------------------------------------------------------------------------!
if(iens == 1) then  ! Deep convection
  k = 1
  kr= k + 1  ! the K-level of Grell is equivalent to the BRAMS K+1-level

  dz        =  z_cup(2)-z_cup(1)

  detdo1    = edt*zd(2)*  cdd(1)*dz
  detdo2    = edt*zd(1)
  entdo     = edt*zd(2)*entrd(1)*dz
  subin     =-edt*zd(2)           

  cfxup1(kr,i,j) = 0.
  cfxdn1(k,i,j) =-edt*xmb* zd(1)
  dfxup1(kr,i,j) = 0.
  efxup1(kr,i,j) = 0.
  dfxdn1(kr,i,j) = xmb*(detdo1+detdo2) !edt already is at detdo1,2
  efxdn1(kr,i,j) =-xmb* entdo          !edt already is at entdo
 

!------------------------------------------------------------------------------------------!
! Checking the mass conservation                                                           !
!------------------------------------------------------------------------------------------!
  totmas = detdo1+detdo2-entdo+subin
  if(abs(totmas) > 1.e-6) then
    write (unit=*,fmt='(a)')                 '----------- Subroutine Get_convflx ----------'
    write(unit=*, fmt='(4(a,1x,i3,1x))')     '  K= ',k,'   I=',i,'   J=',j,'   IENS=',iens
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  subin=  ',    subin,'entdo=  ',entdo
    write(unit=*, fmt='(2(a,1x,es10.3,1x))') '  detdo1= ',   detdo1,'detdo2= ',detdo2
    write(unit=*, fmt='(1(a,1x,es10.3,1x))') '  totmas= ',   totmas
    write(unit=*, fmt='(a)')                 '---------------------------------------------'
    stop 'The model will stop since it is not conserving mass...' 
  end if
end if

return
end subroutine get_convflx
!------------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------------!
! Subroutine output_mean_at_inst_analysis_advec                                            !
! Developed by Saulo Freitas (CPTEC/INPE)                                                  !
! Modified by Marcos Longo (EPS/Harvard U.) - Turbulence needs to be called separately,    !
!                                             because parameters vary accordinly to idiffk !
!                                                                                          !
!   The aim of this subroutine is simply save the mean values of some advection fluxes as  !
! at the regular and lite analysis.                                                        !
!------------------------------------------------------------------------------------------!
subroutine output_mean_at_inst_analysis_advec(m1,m2,m3,afxub,afxvb,afxwb,afxum,afxvm,afxwm)
implicit none
integer, intent(in)                    :: m1,m2,m3
real, dimension(m1,m2,m3), intent(in)  :: afxum, afxvm, afxwm
real, dimension(m1,m2,m3), intent(out) :: afxub, afxvb, afxwb
integer                                :: i, j, k

do k=1, m1
  do i=1, m2
    do j=1, m3
      afxub(k,i,j)=afxum(k,i,j)
      afxvb(k,i,j)=afxvm(k,i,j)
      afxwb(k,i,j)=afxwm(k,i,j)
    end do
  end do
end do

return
end subroutine output_mean_at_inst_analysis_advec
!------------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------------!
! Subroutine output_mean_at_inst_analysis_tke                                              !
! Developed by Saulo Freitas (CPTEC/INPE)                                                  !
! Modified by Marcos Longo (EPS/Harvard U.) - Turbulence needs to be called separately,    !
!                                             because just IDIFFK=1,4,5,6,7 produces TKE   !
!                                                                                          !
!   The aim of this subroutine is simply save the mean TKE value at the regular and lite   !
! analysis.                                                                                !
!------------------------------------------------------------------------------------------!
subroutine output_mean_at_inst_analysis_tke(m1,m2,m3,tkepb,tkepm)
implicit none
integer, intent(in)                    :: m1,m2,m3
real, dimension(m1,m2,m3), intent(in)  :: tkepm
real, dimension(m1,m2,m3), intent(out) :: tkepb
integer                                :: i, j, k

do k=1, m1
  do i=1, m2
    do j=1, m3
      tkepb(k,i,j)=tkepm(k,i,j)
    end do
  end do
end do

return
end subroutine output_mean_at_inst_analysis_tke
!------------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------------!
! Subroutine output_mean_at_inst_analysis_turb                                             !
! Developed by Saulo Freitas (CPTEC/INPE)                                                  !
! Modified by Marcos Longo (EPS/Harvard U.) - Turbulence needs to be called separately,    !
!                                             because just IDIFFK=7 produces sigmaw and Tl !
!                                                                                          !
!   The aim of this subroutine is simply save the mean values of Lagrangian time scale and !
! Sigma-w at the regular and lite analysis.                                                !
!------------------------------------------------------------------------------------------!
subroutine output_mean_at_inst_analysis_turb(m1,m2,m3,sigwb,tlb,sigwm,tlm)
implicit none
integer, intent(in)                    :: m1,m2,m3
real, dimension(m1,m2,m3), intent(in)  :: sigwm, tlm
real, dimension(m1,m2,m3), intent(out) :: sigwb, tlb
integer                                :: i, j, k

do k=1, m1
  do i=1, m2
    do j=1, m3
      sigwb(k,i,j)=sigwm(k,i,j)
      tlb  (k,i,j)=tlm  (k,i,j)
    end do
  end do
end do

return
end subroutine output_mean_at_inst_analysis_turb
!------------------------------------------------------------------------------------------!


