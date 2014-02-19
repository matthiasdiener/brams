!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

MODULE RADIATION

CONTAINS

  SUBROUTINE radiate(mzp, mxp, myp, ia, iz, ja, jz, mynum)

    USE mem_tend, ONLY: tend ! INTENT(INOUT)

    USE mem_grid, ONLY: &
         ngrid, time, dtlt, itime1, nzg, nzs, npatch, grid_g, & ! INTENT(IN)
         nnzp, if_adap, zm, zt, naddsc, nzpmax, imonth1,      & ! INTENT(IN)
         idate1, iyear1, centlat, centlon, ztop,dzm ,ngrids     ! INTENT(IN)

    USE mem_leaf, ONLY: leaf_g ! INTENT(IN)

    USE mem_radiate, ONLY: &
         ilwrtyp, iswrtyp, &       ! INTENT(IN)
         radiate_g,        &       ! INTENT(INOUT)
         radfrq, &                 ! INTENT(IN)
         ncall_i, prsnz, prsnzp, & ! INTENT(INOUT)
         lonrad                    ! INTENT(IN)

    USE mem_basic, ONLY: basic_g ! INTENT(INOUT)

    USE mem_scratch, ONLY: &
         vctr1, vctr2, vctr3, vctr4, vctr5, vctr6, vctr7, & ! INTENT(OUT)
         vctr8, vctr9, vctr10, vctr11, vctr12, scratch      ! INTENT(OUT)

    USE mem_micro, ONLY: micro_g ! INTENT(IN)

    USE rconstants, ONLY: cp, cpor, p00, pi180 ! INTENT(IN)

    USE rrad3, ONLY: jday, solfac, ng, nb, nsolb, npsb, nuum, prf, alpha, &
         trf, beta, xp, wght, wlenlo, wlenhi, solar0, ralcs, a0, a1, a2, a3, &
         exptabc, ulim, npartob, npartg, ncog, ncb, ocoef, bcoef, gcoef

    USE ref_sounding, ONLY: pi01dn ! INTENT(IN)

    ! CATT
    !kmlnew
    USE micphys, ONLY: &
         gnu, level, icloud, irain, ipris, & ! INTENT(IN)
         isnow, iaggr, igraup, ihail         ! INTENT(IN)
    USE mem_cuparm, ONLY: cuparm_g, nnqparm  ! INTENT(IN)
    !kmlnew

    USE rad_carma, ONLY: radcomp_carma ! Subrotuine
!srf- this obsolete with ccatt-brams
!    USE catt_start, ONLY: catt ! INTENT(IN)
!    USE mem_scalar, ONLY: scalar_g ! INTENT(INOUT)
!srf-
    ! TEB_SPM
    use teb_spm_start, only: TEB_SPM ! INTENT(IN)
    use mem_teb_common, only: tebc_g ! INTENT(INOUT)
    !
    USE Rad_UKMO, ONLY: ukmo_swintf, &
                        ukmo_lwintf, &
                        InitRadUKMO
    USE UkmoAdapt,ONLY: r8, &
                        Alloc_mem_Rad_UKMO, &
                        Create_mem_ukmo, &
                        initGetOz
                            
    USE rconstants , ONLY : solar

    IMPLICIT NONE
    ! Arguments:
    INTEGER, INTENT(IN) :: mzp, mxp, myp, ia, iz, ja, jz, mynum
    ! Local Variables:
    INTEGER :: koff
    REAL, pointer :: rc_ptr
    
!--(DMK-CCATT-INI)-------------------------------------------------------------------
    REAL :: hrAngleLocal
!--(DMK-CCATT-FIM)-------------------------------------------------------------------
    
    REAL :: solc
    !kmlnew
    REAL :: lwl(mzp,mxp,myp), iwl(mzp,mxp,myp), fracao_liq(mzp,mxp,myp)
    REAL :: rain(mxp,myp)
    REAL :: dummy
!!$  INTEGER :: ncall
    INTEGER :: i, j, k, ncols
    !kmlnew
    ! TEB_SPM
    real, pointer :: emis_town(:,:), alb_town(:,:), ts_town(:,:), g_urban(:,:,:)
    !
    real          :: cdec
    !
!!$    real, allocatable, target :: pm(:,:,:)
!!$    real, pointer     :: p_pm(:,:,:)

    REAL(KIND=r8), PARAMETER :: Sigma_Estrat=0.099999_r8
   
    character(len=*),parameter :: swFile= &
          '/home/lflavio/run/ftu/tables/sp_sw_hadgem1_3r'
    character(len=*),parameter :: lwFile= &
          '/home/lflavio/run/ftu/tables/sp_lw_hadgem1_3' 
    
    !srf -03/03/2013 cloud fraction
    real, parameter :: fxx = 1. ! DMK ORI => 0.2                   

!!$  ! Initial value
!!$  ncall = 0

    IF ((ilwrtyp + iswrtyp)==0) RETURN

!!$    ! Avoiding problemns with scalar_g not allocated
!!$    nullify(p_pm)
!!$    IF (naddsc<3) THEN
!!$       allocate(pm(mzp, mxp, myp), STAT=ierr)
!!$       if (ierr/=0) call fatal_error("ERROR allocating pm (radiate)")
!!$       pm = 0.
!!$       p_pm => pm
!!$    ELSE
!!$       p_pm => scalar_g(3,ngrid)%sclp
!!$    ENDIF

    ! TEB_SPM
    if (TEB_SPM==1) then
       EMIS_TOWN => tebc_g(ngrid)%EMIS_TOWN
       ALB_TOWN  => tebc_g(ngrid)%ALB_TOWN
       TS_TOWN   => tebc_g(ngrid)%TS_TOWN
       G_URBAN   => leaf_g(ngrid)%G_URBAN
    else
       nullify(EMIS_TOWN)
       nullify(ALB_TOWN)
       nullify(TS_TOWN)
       nullify(G_URBAN)
    endif
    !

    CALL tend_accum(mzp, mxp, myp, ia, iz, ja, jz, &
         tend%tht(1), radiate_g(ngrid)%fthrd(1,1,1))

    IF (mod(time+.001, radfrq)<dtlt .or. time<0.001) THEN

       ! Compute solar zenith angle, multiplier for solar constant, sfc albeDO,
       ! and surface upward longwave radiation.
 
       if (TEB_SPM==1) then 
          CALL radprep(TEB_SPM, imonth1, idate1, iyear1, time, itime1, &
               centlat, centlon, lonrad, pi180,                        &
               nzg, nzs, npatch, ia, iz, ja, jz, jday,       &
               leaf_g(ngrid)%soil_water,                               &
               leaf_g(ngrid)%soil_energy,                              &
               leaf_g(ngrid)%soil_text,                                &
               leaf_g(ngrid)%sfcwater_energy,                          &
               leaf_g(ngrid)%sfcwater_depth,                           &
               leaf_g(ngrid)%leaf_class,                               &
               leaf_g(ngrid)%veg_fracarea,                             &
               leaf_g(ngrid)%veg_height,                               &
               leaf_g(ngrid)%veg_albedo,                               &
               leaf_g(ngrid)%patch_area,                               &
               leaf_g(ngrid)%sfcwater_nlev,                            &
               leaf_g(ngrid)%veg_temp,                                 &
               leaf_g(ngrid)%can_temp,                                 &
               solfac,                                                 &
               grid_g(ngrid)%glat,                                     &
               grid_g(ngrid)%glon,                                     &
               radiate_g(ngrid)%rshort,                                &
               radiate_g(ngrid)%rlong,                                 &
               radiate_g(ngrid)%rlongup,                               &
               radiate_g(ngrid)%albedt,                                &
               radiate_g(ngrid)%cosz,                                  &
!--(DMK-CCATT-INI)-------------------------------------------------------------------       
               hrAngleLocal,                                           &	       
!--(DMK-CCATT-FIM)-------------------------------------------------------------------
               cdec,                                                   &
               ! TEB_SPM
               EMIS_TOWN, ALB_TOWN, TS_TOWN, G_URBAN                   )
       else
          
           CALL radprep(TEB_SPM, imonth1, idate1, iyear1, time, itime1, &
               centlat, centlon, lonrad, pi180,                        &
               nzg, nzs, npatch, ia, iz, ja, jz, jday,       &
               leaf_g(ngrid)%soil_water,                               &
               leaf_g(ngrid)%soil_energy,                              &
               leaf_g(ngrid)%soil_text,                                &
               leaf_g(ngrid)%sfcwater_energy,                          &
               leaf_g(ngrid)%sfcwater_depth,                           &
               leaf_g(ngrid)%leaf_class,                               &
               leaf_g(ngrid)%veg_fracarea,                             &
               leaf_g(ngrid)%veg_height,                               &
               leaf_g(ngrid)%veg_albedo,                               &
               leaf_g(ngrid)%patch_area,                               &
               leaf_g(ngrid)%sfcwater_nlev,                            &
               leaf_g(ngrid)%veg_temp,                                 &
               leaf_g(ngrid)%can_temp,                                 &
               solfac,                                                 &
               grid_g(ngrid)%glat,                                     &
               grid_g(ngrid)%glon,                                     &
               radiate_g(ngrid)%rshort,                                &
               radiate_g(ngrid)%rlong,                                 &
               radiate_g(ngrid)%rlongup,                               &
               radiate_g(ngrid)%albedt,                                &
               radiate_g(ngrid)%cosz,                                  &
!--(DMK-CCATT-INI)-------------------------------------------------------------------       
               hrAngleLocal,                                           &	       
!--(DMK-CCATT-FIM)-------------------------------------------------------------------
               cdec                                                    )
       endif
!LFR>        WRITE(88,fmt='(A,I1.1)') 'voltou radprep',ilwrtyp; CALL flush(88)
!LFR> 
!LFR>        i=2
!LFR>        j=15
!LFR>        WRITE(88,fmt='(A,8(E12.5,1X))') '2DX ' &
!LFR> 		      ,radiate_g(ngrid)%albedt(i,j) &
!LFR>                       ,radiate_g(ngrid)%rshort(i,j) &
!LFR> 		      ,radiate_g(ngrid)%rlong(i,j)  &
!LFR>                       ,grid_g(ngrid)%rtgt(i,j)  &
!LFR>                       ,grid_g(ngrid)%f13t(i,j)  &
!LFR>                       ,grid_g(ngrid)%f23t(i,j)  &
!LFR>                       ,grid_g(ngrid)%glat(i,j)  &
!LFR>                       ,grid_g(ngrid)%glon(i,j)  
!LFR>        WRITE(88,fmt='(A,32(E12.5,1X))') '3D1 ' &
!LFR>                       ,(basic_g(ngrid)%theta(i,j,k),k=1,32)   
!LFR>       WRITE(88,fmt='(A,32(E12.5,1X))') '3D2 ' &
!LFR>                       ,(basic_g(ngrid)%pi0(i,j,k),k=1,32)     
!LFR>       WRITE(88,fmt='(A,32(E12.5,1X))') '3D3 ' &
!LFR>                       ,(basic_g(ngrid)%pp(i,j,k),k=1,32)      
!LFR>       WRITE(88,fmt='(A,32(E12.5,1X))') '3D4 ' &
!LFR>                       ,(basic_g(ngrid)%rv(i,j,k),k=1,32)      
!LFR>       WRITE(88,fmt='(A,32(E12.5,1X))') '3D5 ' &
!LFR>                       ,(basic_g(ngrid)%dn0(i,j,k),k=1,32)     
!LFR>       WRITE(88,fmt='(A,32(E12.5,1X))') '3D6 ' &
!LFR>                       ,(basic_g(ngrid)%rtp(i,j,k),k=1,32)     
!LFR>       WRITE(88,fmt='(A,32(E12.5,1X))') '3D7 ' &
!LFR> 	              ,(radiate_g(ngrid)%fthrd(i,j,k),k=1,32) 


       
 
!!$     CALL azero(mzp*mxp*myp,radiate_g(ngrid)%fthrd(1,1,1))
       radiate_g(ngrid)%fthrd = 0.

       !UKMet radiation
       IF (ilwrtyp==5 .or. iswrtyp==5) THEN
          CALL radCompUkmo(ngrids,ia,iz,ja,jz,mxp,myp,mzp,ngrid,swFile,lwFile)
       END IF
       
       ! Carma radiation
       IF (ilwrtyp==4 .or. iswrtyp==4) THEN

!--(DMK-LFR NEC-SX6)----------------------------------------------
          RAIN = 0.
          LWL = 0.
          IWL = 0.
          fracao_liq = 0.
!--(DMK-LFR NEC-SX6)----------------------------------------------

          IF (level==1) then

             ! Not used with LEVEL 1 - Putting zeros
             LWL  = 0.
             IWL  = 0.
             RAIN = 0.

          ELSEIF (level==2) THEN

!!$           call atob(mxp*myp*mzp, micro_g(ngrid)%rcp(:,:,:), LWL)
             LWL = micro_g(ngrid)%rcp
!!$           CALL azero(mzp*mxp*myp,IWL) ! Not used with LEVEL 2 - Putting zeros
             IWL = 0. ! Not used with LEVEL 2 - Putting zeros
             IF (nnqparm(ngrid)/=0) THEN
!!$              DO i=ia,iz
!!$                 DO j=ja,jz
!!$                    ! conv  precip at mm/h
!!$                    RAIN(i,j) = cuparm_g(ngrid)%conprr(i,j)
!!$                 ENDDO
!!$              ENDDO

!--(DMK-CCATT-INI)-----------------------------------------------------
                RAIN(ia:iz,ja:jz)= cuparm_g(ngrid)%conprr(ia:iz,ja:jz)* 3600.
!--(DMK-CCATT-OLD)-----------------------------------------------------
!                rain(ia:iz,ja:jz) = cuparm_g(ngrid)%conprr(ia:iz,ja:jz)
!--(DMK-CCATT-FIM)-----------------------------------------------------

             ENDIF

          ELSEIF (level>=3) THEN

             IF (nnqparm(ngrid)/=0) THEN
!!$              DO i=ia,iz
!!$                 DO j=ja,jz
!!$                    ! conv + explic  precip at mm/h
!!$                    RAIN(i,j) = (cuparm_g(ngrid)%conprr(i,j) + &
!!$                         micro_g(ngrid)%pcpg(i,j))*3600.
!!$                 ENDDO
!!$              ENDDO
                rain(ia:iz,ja:jz) = cuparm_g(ngrid)%conprr(ia:iz,ja:jz) + &
                     micro_g(ngrid)%pcpg(ia:iz,ja:jz)
                rain(ia:iz,ja:jz) = rain(ia:iz,ja:jz)*3600.
             ENDIF

!!$           CALL azero2(mzp*mxp*myp,LWL,IWL)
             LWL = 0.
             IWL = 0.
             IF (icloud>0) &
!!$                call accum (mxp*myp*mzp, LWL, micro_g(ngrid)%rcp(:,:,:))
                  LWL = LWL + micro_g(ngrid)%rcp

             !kml           IF(irain>0)         call accum (mxp * myp * mzp,LWL,micro_g(ngrid)%rrp(1,1,1))

             IF (igraup>0) THEN
                DO k=1,mzp
                   DO i=ia,iz
                      DO j=ja,jz
                         call qtc(micro_g(ngrid)%Q6(k,i,j), dummy, &
                              fracao_liq(k,i,j))
                      ENDDO
                   ENDDO
                ENDDO

!!$              call ae1t1p1(mxp*myp*mzp, LWL, fracao_liq,  &
!!$                   micro_g(ngrid)%rgp(:,:,:), LWL)
                LWL = fracao_liq*micro_g(ngrid)%rgp + LWL

                fracao_liq(:,:,:) = 1. - fracao_liq(:,:,:)

!!$              call ae1t1p1(mxp*myp*mzp, IWL, fracao_liq, &
!!$                      micro_g(ngrid)%rgp(:,:,:), IWL)
                IWL = fracao_liq*micro_g(ngrid)%rgp + IWL

             ENDIF

             IF (ihail>0) THEN
                DO k=1,mzp
                   DO i=ia,iz
                      DO j=ja,jz
                         call qtc(micro_g(ngrid)%Q7(k,i,j), dummy, &
                              fracao_liq(k,i,j))
                      ENDDO
                   ENDDO
                ENDDO

!!$              call ae1t1p1(mxp*myp*mzp, LWL, fracao_liq, &
!!$                   micro_g(ngrid)%rhp(:,:,:), LWL)
                LWL = fracao_liq*micro_g(ngrid)%rhp + LWL

                fracao_liq(:,:,:) = 1. - fracao_liq(:,:,:)

!!$              call ae1t1p1(mxp*myp*mzp, IWL, fracao_liq, &
!!$                   micro_g(ngrid)%rhp(:,:,:), IWL)
                IWL = fracao_liq*micro_g(ngrid)%rhp + IWL

             ENDIF

             IF (iaggr>0) &
!!$                call accum (mxp*myp*mzp, IWL, micro_g(ngrid)%rap(:,:,:))
                  IWL = IWL + micro_g(ngrid)%rap

             IF (isnow>0) &
!!$                call accum (mxp*myp*mzp, IWL, micro_g(ngrid)%rsp(:,:,:))
                  IWL = IWL + micro_g(ngrid)%rsp

             IF (ipris>0) &
!!$                call accum (mxp*myp*mzp, IWL, micro_g(ngrid)%rpp(:,:,:))
                  IWL = IWL + micro_g(ngrid)%rpp

          ENDIF

          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>SRF
          !             DO k=2,mzp-1
          !              DO i=ia,iz
          !               DO j=ja,jz
          !                LWL(k,i,j) = max(LWL(k,i,j),0.)
          !                IWL(k,i,j) = max(IWL(k,i,j),0.)
          !                !if(LWL(k,i,j) .gt. 1. .or. IWL(k,i,j).gt.1.) then
          !                !  print*,'WL=',i,j,k,IWL(k,i,j), LWL(k,i,j)        
          !                !  stop 444
          !                !endif
          !               ENDDO
          !              ENDDO
          !             ENDDO
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>SRF

          !kmln
	  
!--(DMK-CCATT-INI)-------------------------------------------------------------------       
          !-srf - including a 'cloud fraction" :
          !- included 01.03.2013 to improve near surf temp and precip
          RAIN= fxx * RAIN
          LWL = fxx * LWL
          IWL = fxx * IWL

          CALL radcomp_carma(mzp,mxp,myp,ia,iz,ja,jz,solfac  &
                  ,basic_g(ngrid)%theta       &
                  ,basic_g(ngrid)%pi0         &
                  ,basic_g(ngrid)%pp          &
                  ,basic_g(ngrid)%rv          &
                  !kmlnew
                  ,RAIN,LWL,IWL               &
                  !kmlnew
                  ,basic_g(ngrid)%dn0         &
                  ,basic_g(ngrid)%rtp         &
                  ,radiate_g(ngrid)%fthrd     &
                  ,grid_g(ngrid)%rtgt         &
                  ,grid_g(ngrid)%f13t         &
                  ,grid_g(ngrid)%f23t         &
                  ,grid_g(ngrid)%glat         &
                  ,grid_g(ngrid)%glon         &
                  ,radiate_g(ngrid)%rshort    &
                  ,radiate_g(ngrid)%rlong     &
                  ,radiate_g(ngrid)%albedt    &
                  ,radiate_g(ngrid)%cosz      &
                  ,radiate_g(ngrid)%rlongup   &
                  ,mynum                      &
                  !srf - Carma arrays
                  ,grid_g(ngrid)%fmapt        &
                  !kmlnew
                  ,leaf_g(ngrid)%patch_area   &
                  ,npatch                     &
                  ,hrAngleLocal               &
                  )
!--(DMK-CCATT-OLD)-------------------------------------------------------------------       
!- this is now obsolete with CCATT-BRAMS
!          ! Avoiding problemns with scalar_g not allocated
!	   IF (naddsc>=3) THEN !SCALAR_G PRESENT
!	      CALL radcomp_carma(mzp,mxp,myp,ia,iz,ja,jz,solfac  &
!		   ,basic_g(ngrid)%theta       &
!		   ,basic_g(ngrid)%pi0         &
!		   ,basic_g(ngrid)%pp	       &
!		   ,basic_g(ngrid)%rv	       &
!		   !kmlnew
!		   ,RAIN,LWL,IWL	       &
!		   !kmlnew
!		   ,basic_g(ngrid)%dn0         &
!		   ,basic_g(ngrid)%rtp         &
!		   ,radiate_g(ngrid)%fthrd     &
!		   ,grid_g(ngrid)%rtgt         &
!		   ,grid_g(ngrid)%f13t         &
!		   ,grid_g(ngrid)%f23t         &
!		   ,grid_g(ngrid)%glat         &
!		   ,grid_g(ngrid)%glon         &
!		   ,radiate_g(ngrid)%rshort    &
!		   ,radiate_g(ngrid)%rlong     &
!		   ,radiate_g(ngrid)%albedt    &
!		   ,radiate_g(ngrid)%cosz      &
!		   ,radiate_g(ngrid)%rlongup   &
!		   ,mynum		       &
!		   !srf - Carma arrays
!		   ,grid_g(ngrid)%fmapt        &
!		   !kmlnew
!		   ,leaf_g(ngrid)%patch_area   &
!		   ,npatch		       &
!		   !klmlnew
!		   ,scalar_g(3,ngrid)%sclp     & ! Optional argument for particulate material
!		   )
!	   ELSE
!	      !-srf - including a 'cloud fraction" :
!	      !- included 01.03.2013 to improve near surf temp and precip
!	      RAIN= fxx * RAIN
!	      LWL = fxx * LWL
!	      IWL = fxx * IWL
!	      !-------------------------  
!	      CALL radcomp_carma(mzp,mxp,myp,ia,iz,ja,jz,solfac  &
!		  ,basic_g(ngrid)%theta       &
!                  ,basic_g(ngrid)%pi0         &
!                  ,basic_g(ngrid)%pp          &
!                  ,basic_g(ngrid)%rv          &
!                  !kmlnew
!                  ,RAIN,LWL,IWL               &
!                  !kmlnew
!                  ,basic_g(ngrid)%dn0         &
!                  ,basic_g(ngrid)%rtp         &
!                  ,radiate_g(ngrid)%fthrd     &
!                  ,grid_g(ngrid)%rtgt         &
!                  ,grid_g(ngrid)%f13t         &
!                  ,grid_g(ngrid)%f23t         &
!                  ,grid_g(ngrid)%glat         &
!                  ,grid_g(ngrid)%glon         &
!                  ,radiate_g(ngrid)%rshort    &
!                  ,radiate_g(ngrid)%rlong     &
!                  ,radiate_g(ngrid)%albedt    &
!                  ,radiate_g(ngrid)%cosz      &
!                  ,radiate_g(ngrid)%rlongup   &
!                  ,mynum                      &
!                  !srf - Carma arrays
!                  ,grid_g(ngrid)%fmapt        &
!                  !kmlnew
!                  ,leaf_g(ngrid)%patch_area   &
!                  ,npatch                     &
!                  !klmlnew
!                  !,scalar_g(3,ngrid)%sclp     & ! Optional argument for particulate material
!                  )
!!          ENDIF
!srf-
!--(DMK-CCATT-FIM)-------------------------------------------------------------------       


          !xxxxxxxxxxxxxxxxxxxxxx                
          !do i=ia,iz
          !   do j=ja,jz
          !      if (radiate_g(ngrid)%rlong(i,j).lt.10.0)&
          !         print*,'apos radcomp_carma',i,j,radiate_g(ngrid)%rlong(i,j)
          !      do k=2,mzp-1
          !         if ( (radiate_g(ngrid)%fthrd(k,i,j))*86400 .lt. -20.0)then
          !            print*,'apos radcomp_carma I J K MYNUM',i,j,k,mynum
          !            print*,'apos radcomp_carma RLONG FTHRD',  &
          !                   radiate_g(ngrid)%rlong(i,j), &
          !                   (radiate_g(ngrid)%fthrd(k,i,j))*86400
          !         endif
          !
          !      enddo
          !   enddo
          !enddo
          !xxxxxxxxxxxxxxxxxxxxxx

       END IF
       !END IF

       IF (ilwrtyp<=2 .or. iswrtyp<=2) THEN

          ! IF using Mahrer-Pielke and/or Chen-Cotton radiation, CALL radcomp.
          CALL radcomp(mzp, mxp, myp, ia, iz, ja, jz, solfac  &
               ,basic_g(ngrid)%theta     (1,1,1)  &
               ,basic_g(ngrid)%pi0       (1,1,1)  &
               ,basic_g(ngrid)%pp        (1,1,1)  &
               ,basic_g(ngrid)%rv        (1,1,1)  &
               ,basic_g(ngrid)%dn0       (1,1,1)  &
               ,basic_g(ngrid)%rtp       (1,1,1)  &
               ,radiate_g(ngrid)%fthrd   (1,1,1)  &
               ,grid_g(ngrid)%rtgt       (1,1)    &
               ,grid_g(ngrid)%f13t       (1,1)    &
               ,grid_g(ngrid)%f23t       (1,1)    &
               ,grid_g(ngrid)%glat       (1,1)    &
               ,grid_g(ngrid)%glon       (1,1)    &
               ,radiate_g(ngrid)%rshort  (1,1)    &
               ,radiate_g(ngrid)%rlong   (1,1)    &
               ,radiate_g(ngrid)%albedt  (1,1)    &
               ,radiate_g(ngrid)%cosz    (1,1)    &
               ,radiate_g(ngrid)%rlongup (1,1)    &
               ,cdec                              &
               ,mynum                             )

       END IF

       IF (iswrtyp==3 .or. ilwrtyp==3) THEN

          ! Using Harrington radiation

          IF (nCALL_i==0) THEN

             ! IF first CALL for this node, initialize several quantities & 
             ! Mclatchy sounding DATA.

             CALL radinit(ng, nb, nsolb, npsb, nuum, prf, alpha, trf, beta  &
                  ,xp, wght, wlenlo, wlenhi, solar0, ralcs, a0, a1, a2, a3, solc&
                  ,exptabc, ulim, npartob, npartg, ncog, ncb  &
                  ,ocoef, bcoef, gcoef, gnu)

             prsnz  = (pi01dn(nnzp(1)-1,1)/cp)**cpor*p00
             prsnzp = (pi01dn(nnzp(1)  ,1)/cp)**cpor*p00

             CALL mclatchy(1,mzp,IF_adap,koff  &
                  ,prsnz,prsnzp  &
                  ,grid_g(ngrid)%glat       (1,1)  &
                  ,grid_g(ngrid)%rtgt       (1,1)  &
                  ,grid_g(ngrid)%topt       (1,1)  &
                  ,radiate_g(ngrid)%rlongup (1,1)  &
                  ,zm,zt,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7  &
                  ,vctr8,vctr9,vctr10,vctr11,vctr12)

             nCALL_i = nCALL_i + 1
          END IF

          ! For any CALL, interpolate the mclatchy sounding DATA by latitude and
          ! season.

          CALL mclatchy(2,mzp,IF_adap,koff  &
               ,prsnz,prsnzp  &
               ,grid_g(ngrid)%glat       (1,1)  &
               ,grid_g(ngrid)%rtgt       (1,1)  &
               ,grid_g(ngrid)%topt       (1,1)  &
               ,radiate_g(ngrid)%rlongup (1,1)  &
               ,zm,zt,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7  &
               ,vctr8,vctr9,vctr10,vctr11,vctr12)

          ! IF using Harrington radiation with moisture complexity LEVEL < 3,
          ! CALL radcomp3 which is a substitute driving structure to the bulk
          ! microphysics.

          IF (level<=2) THEN
             IF (level==2) THEN
                rc_ptr => micro_g(ngrid)%rcp(1,1,1)
             ELSE
!!$              CALL azero(mzp*mxp*myp, scratch%vt3dp(1))
                scratch%vt3dp = 0.
                rc_ptr => scratch%vt3dp(1)
             END IF

             CALL radcomp3(mzp,mxp,myp,ia,iz,ja,jz  &
                  ,grid_g(ngrid)%lpw        (1,1)    &
                  ,grid_g(ngrid)%glat       (1,1)    &
                  ,grid_g(ngrid)%rtgt       (1,1)    &
                  ,grid_g(ngrid)%topt       (1,1)    &
                  ,radiate_g(ngrid)%albedt  (1,1)    &
                  ,radiate_g(ngrid)%cosz    (1,1)    &
                  ,radiate_g(ngrid)%rlongup (1,1)    &
                  ,radiate_g(ngrid)%rshort  (1,1)    &
                  ,radiate_g(ngrid)%rlong   (1,1)    &
                  ,basic_g(ngrid)%rv        (1,1,1)  &
                  ,basic_g(ngrid)%dn0       (1,1,1)  &
                  ,radiate_g(ngrid)%fthrd   (1,1,1)  &
                  ,basic_g(ngrid)%pi0       (1,1,1)  &
                  ,basic_g(ngrid)%pp        (1,1,1)  &
                  ,basic_g(ngrid)%theta     (1,1,1)  &
                  ,rc_ptr  )
          END IF

       END IF
      !LFR  WRITE(88,fmt='(A,I2.2)') 'voltou radiacao',ilwrtyp; CALL flush(88)
              !LFR  DO i=1,mxp
                 !LFR  DO j=1,myp
                    !LFR  WRITE(88,fmt='(2(I2.2,1X),2(E12.5,1X))') i,j, &
 		          !LFR  radiate_g(ngrid)%rshort(i,j),radiate_g(ngrid)%rlong(i,j)
                 !LFR  END DO
              !LFR  END DO
        !LFR  WRITE(88,fmt='(A)') '-----------------------...-----------------------'; call flush(88)
    END IF

  
  END SUBROUTINE radiate

! ****************************************************************************

  SUBROUTINE radprep(TEB_SPM, imonth1, idate1, iyear1, time, itime1,        &
       centlat, centlon, lonrad, pi180,                                     &
       mzg, mzs, np, ia, iz, ja, jz, jday,                          &
       soil_water, soil_energy, soil_text, sfcwater_energy, sfcwater_depth, &
       leaf_class, veg_fracarea, veg_height, veg_albedo, patch_area,        &
       sfcwater_nlev, veg_temp, can_temp,                                   &
       solfac, glat, glon, rshort, rlong, rlongup, albedt, cosz, cdec,      &
!--(DMK-CCATT-INI)-------------------------------------------------------------------       
       hrAngleLocal,                                                        &
!--(DMK-CCATT-FIM)-------------------------------------------------------------------       
       ! TEB_SPM
       EMIS_TOWN, ALB_TOWN, TS_TOWN, G_URBAN                                &
       !
       )

    ! TEB_SPM
!!$    use teb_spm_start, only: TEB_SPM !INTENT(IN)

    use mem_leaf, only: isfcl  !DSM
  
    IMPLICIT NONE

    ! Arguments:
    INTEGER, INTENT(IN)      :: TEB_SPM
    INTEGER, INTENT(IN)      :: imonth1, idate1, iyear1, itime1
    REAL, INTENT(IN)         :: time
    REAL, INTENT(IN)         :: centlat(:), centlon(:)
    INTEGER, INTENT(IN)      :: lonrad
    REAL, INTENT(IN)         :: pi180
    INTEGER, INTENT(IN)      :: mzg, mzs, np, ia, iz, ja, jz
    INTEGER, INTENT(OUT)     :: jday
    !DIMENSION(mzg,m2,m3,np)
    REAL, INTENT(IN)         :: soil_water(:,:,:,:), soil_energy(:,:,:,:), &
         soil_text(:,:,:,:)
    !DIMENSION(mzs,m2,m3,np)
    REAL, INTENT(IN)         :: sfcwater_energy(:,:,:,:), &
         sfcwater_depth(:,:,:,:)
    !DIMENSION(m2,m3,np)
    REAL, INTENT(IN)         :: leaf_class(:,:,:), veg_fracarea(:,:,:), &
         veg_height(:,:,:), veg_albedo(:,:,:), patch_area(:,:,:),       &
         sfcwater_nlev(:,:,:), veg_temp(:,:,:), can_temp(:,:,:)
    !TEB_SPM
    !DIMENSION(m2,m3)
    REAL, pointer, optional :: EMIS_TOWN(:,:), ALB_TOWN(:,:), TS_TOWN(:,:)
    !DIMENSION(m2,m3,np)
    real, pointer, optional :: G_URBAN(:,:,:)
    !
    REAL, INTENT(OUT)        :: solfac
    !DIMENSION(m2,m3)
    REAL, INTENT(IN)         :: glat(:,:), glon(:,:), rshort(:,:), rlong(:,:)
    REAL, INTENT(OUT)        :: rlongup(:,:), albedt(:,:)
    REAL, INTENT(INOUT)      :: cosz(:,:)
    REAL, INTENT(OUT)        :: cdec
    
!--(DMK-CCATT-INI)-------------------------------------------------------------------       
    REAL, INTENT(OUT)        :: hrAngleLocal
!--(DMK-CCATT-FIM)-------------------------------------------------------------------       
        
    ! Local Variables
    ! Needed by TEB_SPM
    real, pointer :: L_EMIS_TOWN, L_ALB_TOWN, L_TS_TOWN, L_G_URBAN
    !
    INTEGER :: ip, i, j
    ! REAL :: c1,c2

    INTERFACE

       SUBROUTINE sfcrad(mzg, mzs, ip,                  &
            soil_energy, soil_water, soil_text,         &
            sfcwater_energy, sfcwater_depth,            &
            patch_area, can_temp, veg_temp, leaf_class, &
            veg_height, veg_fracarea, veg_albedo,       &
            sfcwater_nlev, rshort, rlong, albedt,       &
            rlongup, cosz,                              &
            G_URBAN, ETOWN, ALBTOWN, TSTOWN             )
         integer, intent(IN) :: mzg, mzs, ip
         real, intent(IN)    :: soil_energy(mzg)
         real, intent(IN)    :: soil_water(mzg)
         real, intent(IN)    :: soil_text(mzg)
         real, intent(IN)    :: sfcwater_energy(mzs)
         real, intent(IN)    :: sfcwater_depth(mzs)
         real, intent(IN)    :: patch_area
         real, intent(IN)    :: can_temp
         real, intent(IN)    :: veg_temp
         real, intent(IN)    :: leaf_class
         real, intent(IN)    :: veg_height
         real, intent(IN)    :: veg_fracarea
         real, intent(IN)    :: veg_albedo
         real, intent(IN)    :: sfcwater_nlev
         real, intent(IN)    :: rshort
         real, intent(IN)    :: rlong
         real, intent(INOUT) :: albedt
         real, intent(INOUT) :: rlongup
         real, intent(IN)    :: cosz
         real, pointer, optional :: G_URBAN, ETOWN, ALBTOWN, TSTOWN
       END SUBROUTINE sfcrad
         
    END INTERFACE

    ! TEB_SPM
    nullify(L_EMIS_TOWN)
    nullify(L_ALB_TOWN)
    nullify(L_TS_TOWN)
    nullify(L_G_URBAN)

    ! Compute solar zenith angle [cosz(i,j)] & solar constant factr [solfac].

    CALL zen(imonth1, idate1, iyear1, time, itime1, centlat, centlon, &
         lonrad, pi180, ia, iz, ja, jz, jday, glat, glon, cosz, &
!--(DMK-CCATT-INI)-------------------------------------------------------------------       
         solfac, hrAngleLocal, cdec)
!--(DMK-CCATT-OLD)-------------------------------------------------------------------       
!         solfac, cdec)
!--(DMK-CCATT-FIM)-------------------------------------------------------------------       

    ! Compute patch-averaged surface albeDO [albedt(i,j)] and up longwave
    ! radiative flux [rlongup(i,j)].

  if (isfcl /= 5 .or. time==0.) then !DSM ---- o rlongup e o albedt deve vir do proprio jules
!!$  CALL azero2(m2*m3,albedt,rlongup)
    albedt  = 0.
    rlongup = 0.

    DO ip = 1,np
       DO j = 1,jz
          DO i = 1,iz

             ! TEB_SPM
             if (TEB_SPM==1 .and. present(G_URBAN) .and.           &
                  present(EMIS_TOWN) .and. present(ALB_TOWN) .and. &
                  present(TS_TOWN))                                then
                L_G_URBAN   => G_URBAN(i,j,ip)
                L_EMIS_TOWN => EMIS_TOWN(i,j)
                L_ALB_TOWN  => ALB_TOWN(i,j)
                L_TS_TOWN   => TS_TOWN(i,j)
                CALL sfcrad(mzg, mzs, ip,                                     &
                     soil_energy(1:mzg,i,j,ip), soil_water(1:mzg,i,j,ip),     &
                     soil_text(1:mzg,i,j,ip),   sfcwater_energy(1:mzs,i,j,ip),&
                     sfcwater_depth(1:mzs,i,j,ip), patch_area(i,j,ip),        &
                     can_temp(i,j,ip),         veg_temp(i,j,ip),              &
                     leaf_class(i,j,ip),       veg_height(i,j,ip),            &
                     veg_fracarea(i,j,ip),     veg_albedo(i,j,ip),            &
                     sfcwater_nlev(i,j,ip),                                   &
                     rshort(i,j), rlong(i,j), albedt(i,j),                    &
                     rlongup(i,j), cosz(i,j),                                 &
                     ! TEB_SPM
                     L_G_URBAN, L_EMIS_TOWN, L_ALB_TOWN, L_TS_TOWN            &
                     !
                     )
             else
                CALL sfcrad(mzg, mzs, ip,                                     &
                     soil_energy(1:mzg,i,j,ip), soil_water(1:mzg,i,j,ip),     &
                     soil_text(1:mzg,i,j,ip),   sfcwater_energy(1:mzs,i,j,ip),&
                     sfcwater_depth(1:mzs,i,j,ip), patch_area(i,j,ip),        &
                     can_temp(i,j,ip),         veg_temp(i,j,ip),              &
                     leaf_class(i,j,ip),       veg_height(i,j,ip),            &
                     veg_fracarea(i,j,ip),     veg_albedo(i,j,ip),            &
                     sfcwater_nlev(i,j,ip),                                   &
                     rshort(i,j), rlong(i,j), albedt(i,j),                    &
                     rlongup(i,j), cosz(i,j)                                  &
                     )
             endif

!!$             CALL sfcrad(mzg, mzs, ip,                                        &
!!$                  soil_energy(1:mzg,i,j,ip), soil_water(1:mzg,i,j,ip),        &
!!$                  soil_text(1:mzg,i,j,ip),   sfcwater_energy(1:mzs,i,j,ip),   &
!!$                  sfcwater_depth(1:mzs,i,j,ip), patch_area(i,j,ip),           &
!!$                  can_temp(i,j,ip),         veg_temp(i,j,ip),                 &
!!$                  leaf_class(i,j,ip),       veg_height(i,j,ip),               &
!!$                  veg_fracarea(i,j,ip),     veg_albedo(i,j,ip),               &
!!$                  sfcwater_nlev(i,j,ip),                                      &
!!$                  rshort(i,j), rlong(i,j), albedt(i,j),rlongup(i,j),cosz(i,j),&
!!$                  ! TEB_SPM
!!$                  L_G_URBAN, L_EMIS_TOWN, L_ALB_TOWN, L_TS_TOWN               &
!!$                  !
!!$                  )

          END DO
       END DO
    END DO

  endif  !DSM
  
  END SUBROUTINE radprep

  ! ***************************************************************************

  SUBROUTINE zen(imonth1, idate1, iyear1, time, itime1, centlat, centlon, &
       lonrad, pi180, &
!--(DMK-CCATT-INI)-------------------------------------------------------------------       
       ia, iz, ja, jz, jday, glat, glon, cosz, solfac, hrangle, cdec)
!--(DMK-CCATT-OLD)-------------------------------------------------------------------       
!       ia, iz, ja, jz, jday, glat, glon, cosz, solfac, cdec)
!--(DMK-CCATT-FIM)-------------------------------------------------------------------       

    USE ModDateUtils
!!$    USE mem_grid, ONLY: imonth1, idate1, iyear1, time, itime1, &
!!$         centlat, centlon
!!$    USE mem_radiate, ONLY: lonrad
!!$    USE rconstants , ONLY: pi180

!--(DMK-CCATT-INI)-----------------------------------------------------
    USE mem_radiate, only: radfrq
!--(DMK-CCATT-FIM)-----------------------------------------------------

    IMPLICIT NONE
    ! Arguments:
    INTEGER, INTENT(IN)  :: imonth1, idate1, iyear1, itime1
    REAL, INTENT(IN)     :: time
    REAL, INTENT(IN)     :: centlat(:), centlon(:)
    INTEGER, INTENT(IN)  :: lonrad
    REAL, INTENT(IN)     :: pi180
    INTEGER, INTENT(IN)  :: ia, iz, ja, jz
    INTEGER, INTENT(OUT) :: jday
    REAL, INTENT(IN)     :: glat(:,:) !(m2,m3)
    REAL, INTENT(IN)     :: glon(:,:) !(m2,m3)
    REAL, INTENT(INOUT)  :: cosz(:,:) !(m2,m3)
    REAL, INTENT(OUT)    :: solfac

!--(DMK-CCATT-INI)-----------------------------------------------------
    REAL, INTENT(OUT)    :: hrangle
!--(DMK-CCATT-FIM)-----------------------------------------------------

    REAL, INTENT(OUT)    :: cdec
    ! Local Variables:
    INTEGER :: i, j! , julday
    REAL    :: sdec, declin, d0, d02, dayhr, radlat, cslcsd, snlsnd, gglon, &
         dayhrr !,tdec
 
!--(DMK-CCATT-OLD)-----------------------------------------------------
!    REAL    :: hrangl
!--(DMK-CCATT-INI)-----------------------------------------------------

!--(DMK-CCATT-INI)-----------------------------------------------------
    REAL :: eqt
!--(DMK-CCATT-FIM)-----------------------------------------------------

    !common /radcom/ tdec,sdec,cdec,declin,rvr,rtr,dn0r,pird,prd,temprd  &
    !     ,fthrl,dzmr,dztr,fthrs

    jday   = julday(imonth1, idate1, iyear1)
    jday   = jday + nint(time/86400.)
    !      sdec - sine of declination, cdec - cosine of declination
    declin = -23.5*cos(6.283/365.*(jday + 9))*pi180
    sdec   = sin(declin)
    cdec   = cos(declin)

    ! Find the factor, solfac, to multiply the solar constant to correct
    ! for Earth's varying distance to the sun.

    d0     = 6.2831853*float(jday-1)/365.
    d02    = d0*2.
    solfac = 1.000110 + 0.034221*cos(d0) + 0.001280*sin(d0) + &
         0.000719*cos(d02) + 0.000077*sin(d02)

    ! Find the hour angle, THEN get cosine of zenith angle.

!--(DMK-CCATT-INI)-----------------------------------------------------
    !NER_i - including solar time equation ("eqt" must be defined, it is a new variable)
    eqt = (0.000075 + 0.001868*cos(d0) - 0.032077*sin(d0) - 0.014615*cos(d02) &
 	  - 0.040849*sin(d02))*1440/(2*3.141593)	   
    !NER_f - including solar time equation	
  
    dayhr = (time / 3600. + float(itime1/100) + float(mod(itime1,100)) / 60.)+ (radfrq/(2*3600))
    !NER (radfrq/(2*3600)) - rad transfer shift half of radfrq(improving rad tendency representativity)
!--(DMK-CCATT-OLD)-----------------------------------------------------
!    dayhr  = time/3600. + float(itime1/100) + float(mod(itime1,100))/60.
!--(DMK-CCATT-FIM)-----------------------------------------------------

    DO j = ja,jz
       DO i = ia,iz
          radlat = glat(i,j)*pi180
          IF (lonrad==0)      radlat = centlat(1)*pi180
          IF (radlat==declin) radlat = radlat + 1.e-5
          cslcsd = cos(radlat)*cdec
          snlsnd = sin(radlat)*sdec
          gglon  = glon(i,j)
          IF (lonrad==0)      gglon = centlon(1)

!--(DMK-CCATT-INI)-----------------------------------------------------
	!NER_i new hour angle calculation		
	hrangle=((dayhr+(gglon/15.)-(12.-eqt/60.))*15./1.)*3.141593/180.	
!--(DMK-CCATT-OLD)-----------------------------------------------------
!          dayhrr    = mod(dayhr+gglon/15.+24., 24.)
!          hrangl    = 15.*(dayhrr - 12.)*pi180
!--(DMK-CCATT-FIM)-----------------------------------------------------

          cosz(i,j) = snlsnd + cslcsd*cos(hrangle)
          !-srf - cosz > 1 no SX6

!--(DMK-CCATT-OLD)-----------------------------------------------------
!          cosz(i,j) = min(cosz(i,j)+1.0E-10, 1.0) !LFR: Prevent 90 degrees z angle
!--(DMK-CCATT-FIM)-----------------------------------------------------

          !cosz(i,j) = min(cosz(i,j), 1.0)
          !cosz(i,j) = max(cosz(i,j),-1.0)

!--(DMK-CCATT-INI)-----------------------------------------------------
          cosz(i,j) = max(cosz(i,j),1.0E-10)
!--(DMK-CCATT-FIM)-----------------------------------------------------

          !-srf        
       END DO
    END DO

  END SUBROUTINE zen

END MODULE RADIATION

!*****************************************************************************

SUBROUTINE tend_accum(m1,m2,m3,ia,iz,ja,jz,at,at2)

  IMPLICIT NONE
  ! Arguments:
  INTEGER, INTENT(IN) :: m1, m2, m3, ia, iz, ja, jz
  REAL, INTENT(INOUT) :: at(m1,m2,m3)
  REAL, INTENT(IN)    :: at2(m1,m2,m3)
  ! Local Variables:
  INTEGER :: i, j, k

  DO j=ja,jz
     DO i=ia,iz
        DO k=1,m1
           at(k,i,j) = at(k,i,j) + at2(k,i,j)
        END DO
     END DO
  END DO

END SUBROUTINE tend_accum

!*****************************************************************************




!*****************************************************************************

SUBROUTINE radcomp(m1,m2,m3,ia,iz,ja,jz,solfac  &
     ,theta,pi0,pp,rv,dn0,rtp,fthrd  &
     ,rtgt,f13t,f23t,glat,glon,rshort,rlong,albedt,cosz,rlongup  &
     ,cdec,mynum)

  USE mem_grid   , ONLY : dzm,dzt,nzp,itopo,plonn,ngrid,time,itime1,centlon, &

!--(DMK-CCATT-INI)-----------------------------------------------------
                          grid_g

  USE mem_radiate, ONLY : radfrq			  
!--(DMK-CCATT-FIM)-----------------------------------------------------

  USE grid_dims  , ONLY : nzpmax
  USE mem_scratch, ONLY : vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7, &
                          vctr8,vctr9,vctr10,vctr11,vctr12,vctr13,vctr14, &
                          vctr15,scratch
  USE mem_radiate, ONLY : ilwrtyp,iswrtyp,lonrad
  USE rconstants , ONLY : cp,cpor,p00,stefan,solar,pi180

  IMPLICIT NONE

  INTEGER :: m1,m2,m3,ia,iz,ja,jz,i,j,k,mynum,l

  REAL :: solfac,cdec,dzsdx,dzsdy,dlon,a1,a2,dayhr,gglon  &
       ,dayhrr,hrangl,sinz,sazmut,slazim,slangl,cosi !tdec, sdec, declin
  REAL, DIMENSION(nzpmax) :: rvr,rtr,dn0r,pird,prd,fthrl,dzmr,dztr,fthrs
  REAL, DIMENSION(nzpmax+1) :: temprd
  REAL, DIMENSION(m1,m2,m3) :: theta,pi0,pp,rv,dn0,rtp,fthrd
  REAL, DIMENSION(m2,m3) :: rtgt,f13t,f23t,glat,glon,rshort,rlong,cosz  &
       ,albedt,rlongup

!--(DMK-CCATT-INI)-----------------------------------------------------
  integer :: k1,k2,nlev  !-srf-rams60
!--(DMK-CCATT-FIM)-----------------------------------------------------

  !common /radcom/ tdec,sdec,cdec,declin,rvr,rtr,dn0r,pird,prd,temprd  &
  !     ,fthrl,dzmr,dztr,fthrs

  DO j = ja,jz
     DO i = ia,iz

!--(DMK-CCATT-INI)-----------------------------------------------------
 !-srf-rams60
      k2=grid_g(ngrid)%lpw(i,j)
      k1=k2-1
      nlev=m1-k1+1
 !-srf-rams60
!--(DMK-CCATT-FIM)-----------------------------------------------------

        DO k = 1,m1

!--(DMK-CCATT-INI)-----------------------------------------------------
!-srf-rams60
          pird(k-k1+1) = (pp(k,i,j) + pi0(k,i,j)) / cp
          temprd(k-k1+1) = theta(k,i,j) * pird(k-k1+1)
          rvr(k-k1+1) = max(0.,rv(k,i,j))
          rtr(k-k1+1) = max(rvr(k-k1+1),rtp(k,i,j))
          ! Convert the next 7 variables to cgs for now.
          prd(k-k1+1) = pird(k-k1+1) ** cpor * p00 * 10.
          dn0r(k-k1+1) = dn0(k,i,j) * 1.e-3
          dzmr(k-k1+1) = dzm(k) / rtgt(i,j) * 1.e-2
          dztr(k-k1+1) = dzt(k) / rtgt(i,j) * 1.e-2
!-srf-rams60
!--(DMK-CCATT-OLD)-----------------------------------------------------
!!-srf-rams60
!           pird(k) = (pp(k,i,j) + pi0(k,i,j)) / cp
!           temprd(k) = theta(k,i,j) * pird(k)
!           rvr(k) = max(0.,rv(k,i,j))
!           rtr(k) = max(rvr(k),rtp(k,i,j))
!           ! Convert the next 7 variables to cgs for now.
!           prd(k) = pird(k) ** cpor * p00 * 10.
!           dn0r(k) = dn0(k,i,j) * 1.e-3
!           dzmr(k) = dzm(k) / rtgt(i,j) * 1.e-2
!           dztr(k) = dzt(k) / rtgt(i,j) * 1.e-2
!!-srf-rams60
!--(DMK-CCATT-FIM)-----------------------------------------------------

        END DO
        temprd(1) = (rlongup(i,j) / stefan) ** 0.25

!--(DMK-CCATT-INI)-----------------------------------------------------
!-srf-rams60        temprd(nzp+1) = temprd(nzp)
        temprd(m1+1-k1+1) = temprd(m1-k1+1)
!--(DMK-CCATT-OLD)-----------------------------------------------------
!        temprd(nzp+1) = temprd(nzp)
!--(DMK-CCATT-FIM)-----------------------------------------------------

        ! CALL the longwave parameterizations.
	!LFR  IF(i==4 .and. j==13) THEN
	   !LFR  WRITE(*,fmt='(33(F12.6,1X))') cp,(pird(k),k=1,m1)
	   !LFR  WRITE(*,fmt='(32(F12.6,1X))') (pp(k,4,13),k=1,m1)
	   !LFR  WRITE(*,fmt='(32(F12.6,1X))') (pi0(k,4,13),k=1,m1)	   
	!LFR  END IF
        IF (ilwrtyp .eq. 2) THEN

!--(DMK-CCATT-INI)-----------------------------------------------------
           CALL lwradp(nlev,temprd,rvr,dn0r,dztr,pird,vctr1,fthrl,rlong(i,j))
!--(DMK-CCATT-OLD)-----------------------------------------------------
!           CALL lwradp(nzp,temprd,rvr,dn0r,dztr,pird,vctr1,fthrl,rlong(i,j))
!--(DMK-CCATT-FIM)-----------------------------------------------------

       ELSEIF (ilwrtyp .eq. 1) THEN

!--(DMK-CCATT-INI)-----------------------------------------------------
             CALL lwradc(nlev+1,rvr,rtr,dn0r,temprd,prd,dztr,fthrl,rlong(i,j)  &
                  ,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7,vctr8,vctr9,vctr10  &
                  ,vctr11,vctr12,vctr13,vctr14,vctr15)
!--(DMK-CCATT-OLD)-----------------------------------------------------
!           CALL lwradc(nzp+1,rvr,rtr,dn0r,temprd,prd,dztr,fthrl,rlong(i,j)  &
!                ,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7,vctr8,vctr9,vctr10  &
!                ,vctr11,vctr12,vctr13,vctr14,vctr15)
!--(DMK-CCATT-FIM)-----------------------------------------------------

        END IF

        ! The shortwave parameterizations are ONLY valid IF the cosine
        !    of the zenith angle is greater than .03 .

        IF (cosz(i,j) .gt. .03) THEN

           IF (iswrtyp .eq. 2) THEN

!--(DMK-CCATT-INI)-----------------------------------------------------
              CALL shradp(nlev,rvr,dn0r,dzmr,vctr1,pird,cosz(i,j)  &
                   ,albedt(i,j),solar*1e3*solfac,fthrs,rshort(i,j))
!--(DMK-CCATT-OLD)-----------------------------------------------------
!              CALL shradp(nzp,rvr,dn0r,dzmr,vctr1,pird,cosz(i,j)  &
!                   ,albedt(i,j),solar*1e3*solfac,fthrs,rshort(i,j))
!--(DMK-CCATT-FIM)-----------------------------------------------------

           ELSEIF (iswrtyp .eq. 1) THEN

!--(DMK-CCATT-INI)-----------------------------------------------------
              CALL shradc(nlev+1,rvr,rtr,dn0r,dztr,prd,scratch%scr1(:)  &
                   ,albedt(i,j),solar*1.e3*solfac,cosz(i,j),fthrs,rshort(i,j))
!--(DMK-CCATT-OLD)-----------------------------------------------------
!              CALL shradc(nzp+1,rvr,rtr,dn0r,dztr,prd,scratch%scr1(:)  &
!                   ,albedt(i,j),solar*1.e3*solfac,cosz(i,j),fthrs,rshort(i,j))
!--(DMK-CCATT-FIM)-----------------------------------------------------

           END IF

           ! ModIFy the DOwnward surface shortwave flux by considering
           !    the slope of the topography.

           IF (itopo .eq. 1) THEN
              dzsdx = f13t(i,j) * rtgt(i,j)
              dzsdy = f23t(i,j) * rtgt(i,j)

              ! The y- and x-directions must be true north and east for
              ! this correction. the following rotates the model y/x
              ! to the true north/east.

              ! The following rotation seems to be incorrect, so CALL this instead:
              !      SUBROUTINE uvtoueve(u,v,ue,ve,qlat,qlon,platn(ngrid),plonn(ngrid))

              dlon = (plonn(ngrid) - glon(i,j)) * pi180
              a1 = dzsdx*cos(dlon) + dzsdy * sin(dlon)
              a2 = -dzsdx*sin(dlon) + dzsdy * cos(dlon)
              dzsdx = a1
              dzsdy = a2

!--(DMK-CCATT-INI)-----------------------------------------------------
	      !NER (radfrq/(2*3600)) - rad transfer shift half of radfrq(improving rad tendency representativity)   
              dayhr = (time / 3600. + float(itime1/100) + float(mod(itime1,100)) / 60.) &
	              + (radfrq/(2*3600))
!--(DMK-CCATT-OLD)-----------------------------------------------------
!              dayhr = time / 3600. + float(itime1/100)  &
!                   + float(mod(itime1,100)) / 60.
!--(DMK-CCATT-FIM)-----------------------------------------------------
		   
              gglon = glon(i,j)
              IF (lonrad .eq. 0) gglon = centlon(1)
              dayhrr = mod(dayhr+gglon/15.+24.,24.)
              hrangl = 15. * (dayhrr - 12.) * pi180
              ! To avoid floating divide by zero: max(..., 1.e-20)
              sinz = max(sqrt(1. - cosz(i,j)**2), 1.e-20)
!!$              sinz = sqrt(1. - cosz(i,j) ** 2)
              sazmut = asin(max(-1.,min(1.,cdec*sin(hrangl)/sinz)))
              IF (abs(dzsdx) .lt. 1e-20) dzsdx = 1.e-20
              IF (abs(dzsdy) .lt. 1e-20) dzsdy = 1.e-20
              slazim = 1.571 - atan2(dzsdy,dzsdx)
              slangl = atan(sqrt(dzsdx*dzsdx+dzsdy*dzsdy))
              cosi = cos(slangl) * cosz(i,j) + sin(slangl) * sinz  &
                   * cos(sazmut-slazim)

!--(DMK-CCATT-INI)-----------------------------------------------------
              !srf- mod rams60
	      if (cosi > 0.) then
                 !rshort(i,j) = rshort(i,j) * cosi / cosz(i,j) !Previous
		 rshort(i,j) = rshort(i,j)                     !NER 
              else
                ! rshort(i,j) =  0.           !Previous
		 rshort(i,j) = rshort(i,j)    !NER
              endif
              !rshort(i,j) = rshort(i,j) * cosi / cosz(i,j)
	      !srf- mod rams60
!--(DMK-CCATT-OLD)-----------------------------------------------------
!              rshort(i,j) = rshort(i,j) * cosi / cosz(i,j)
!--(DMK-CCATT-FIM)-----------------------------------------------------

           END IF

        ELSE
           DO k = 1,nzp
              fthrs(k) = 0.
           END DO
           rshort(i,j) = 0.
        END IF

!--(DMK-CCATT-INI)-----------------------------------------------------
      ! Add fluxes, adjusting for adap if necessary.
        do k = k2,m1-1
            fthrd(k,i,j) = fthrl(k-k1+1) + fthrs(k-k1+1)
        enddo
!--(DMK-CCATT-OLD)-----------------------------------------------------
!        DO k = 2,m1-1
!           fthrd(k,i,j) = fthrl(k) + fthrs(k)
!        END DO
!--(DMK-CCATT-FIM)-----------------------------------------------------

        ! Convert the DOwnward flux at the ground to SI.

        rshort(i,j) = rshort(i,j) * 1.e-3 / (1. - albedt(i,j))
        rlong(i,j) = rlong(i,j) * 1.e-3

!--(DMK-CCATT-INI)-----------------------------------------------------
      fthrd(k1,i,j) = fthrd(k2,i,j)
!--(DMK-CCATT-OLD)-----------------------------------------------------
!        fthrd(1,i,j) = fthrd(2,i,j)
!--(DMK-CCATT-FIM)-----------------------------------------------------

     END DO
  END DO
  RETURN
END SUBROUTINE radcomp

!******************************************************************************

SUBROUTINE radcomp3(m1,m2,m3,ia,iz,ja,jz,lpw  &
     ,glat,rtgt,topt,albedt,cosz,rlongup,rshort,rlong  &
     ,rv,dn0,fthrd,pi0,pp,theta,rcp)

  USE mem_grid   , ONLY: maxnzp,if_adap,zm,zt,time,ngrid
  USE mem_radiate, ONLY: iswrtyp,ilwrtyp
  USE rconstants , ONLY: cpi,p00,cpor
  USE micphys    , ONLY: level,pwmas,pwmasi,cparm,emb0,emb1,gnu,dnfac,jhcat, &
                         press,tair,emb,cx

  IMPLICIT NONE

  INTEGER :: m1,m2,m3,ia,iz,ja,jz,mcat,i,j,k
  INTEGER, DIMENSION(m2,m3) :: lpw

  REAL :: cfmasi,cparmi,glg,glgm,gammln,picpi
  REAL, DIMENSION(m2,m3) :: glat,rtgt,topt,cosz,albedt,rlongup,rshort,rlong
  REAL, DIMENSION(m1,m2,m3) :: dn0,rv,fthrd,pi0,pp,theta,rcp

  IF (level .le. 1) THEN
     mcat = 0
  ELSE
     mcat = 1
     pwmas(1) = 3.
     pwmasi(1) = 1. / pwmas(1)
     cfmasi = 1. / 524.
     cparmi = 1. / cparm
     emb0(1) = 5.24e-16
     emb1(1) = 3.35e-11
     glg = gammln(gnu(1))
     glgm = gammln(gnu(1) + pwmas(1))
     dnfac(1) = (cfmasi * exp(glg - glgm)) ** pwmasi(1)
     DO k = 2,m1-1
        jhcat(k,1) = 1
     END DO
  END IF

  DO j = ja,jz
     DO i = ia,iz

        DO k = 2,m1-1
           picpi = (pi0(k,i,j) + pp(k,i,j)) * cpi
           press(k) = p00 * picpi ** cpor
           tair(k) = theta(k,i,j) * picpi
	   
!--(DMK-CCATT-INI)-----------------------------------------------------
	   if (level >= 2) then
!--(DMK-CCATT-)-----------------------------------------------------

               emb(k,1) = max(emb0(1),min(emb1(1),rcp(k,i,j) * cparmi))
               cx(k,1) = rcp(k,i,j) / emb(k,1)

!--(DMK-CCATT-INI)-----------------------------------------------------
	   endif
!--(DMK-CCATT-FIM)-----------------------------------------------------
	   
        END DO

        ! ModIF. ALF-Sfreitas

        CALL radcalc3(m1,maxnzp,mcat,iswrtyp,ilwrtyp,IF_adap,lpw(i,j)  &
             ,glat(i,j),rtgt(i,j),topt(i,j),albedt(i,j),cosz(i,j)  &
             ,rlongup(i,j),rshort(i,j),rlong(i,j)  &
             ,zm,zt,rv(1,i,j),dn0(1,i,j),fthrd(1,i,j),i,j,time,ngrid)

     END DO
  END DO
  RETURN
END SUBROUTINE radcomp3

!******************************************************************************



!****************************************************************************

! SUBROUTINE radcalc3: column driver for twostream radiation code

! variables USEd within SUBROUTINE radcalc3:
! ==================================================================

! Variables in rrad3 module parameter statement

!  mb               : maximum allowed number of bands [=8]
!  mg                  : maximum allowed number of gases [=3]
!  mk               : maximum number of pseuDObands allowed for any gas [=7]
!  ncog             : number of fit coefficients (omega and asym) [=5]
!  ncb              : number of fit coefficients (extinction) [=2]
!  npartob          : number of hydrometeor categories (including dIFferent habits)
!  npartg           : number of hydrometeor categories USEd for gc coefficients [=7]
!  nrad                  : total number of radiation levels USEd (m1 - 1 + narad)
!  narad            : number of radiation levels added above model
!  nsolb            : active number of solar bands
!  nb               : active number of bands
!  ng                : active number of gases
!  jday             : julian day
!  solfac           : solar constant multiplier for variable E-S distance
!  ralcs (mb)       : rayleigh scattering integration constants
!  solar1 (mb)      : solar fluxes at top of atmosphere - corrected for ES distance
!  solar0 (mb)      : solar fluxes at top of atmosphere - uncorrected for ES distance
!  nuum (mb)        :    continuum flags
!  a0,a1,a2,a3 (mb) : Planck function fit coefficients
!  npsb (mg,mb)     : number of pseuDO bands
!  trf (mg,mb)      : reference temperature for xp and wght coefficients
!  prf (mg,mb)      : reference pressure for xp and wght coefficients
!  ulim (mg,mb)     : upper bound on pathlength for gases
!  xp (mg,mk,mb)    : coefficient USEd in computing gaseous optical depth
!  alpha (mg,mk,mb) : pressure correction factor exponent
!  beta (mg,mk,mb)  : temperature correction factor exponent
!  wght (mg,mk,mb)  : pseuDO band weight
!  exptabc (150)    : table of exponential function values
!  ocoef(ncog,mb,npartob)  : fit coefficients for hyd. single scatter.
!  bcoef(ncb,mb ,npartob)  : fit coefficients for hyd. extinction coefficient.
!  gcoef(ncog,mb,npartg)   : fit coefficients for hyd. asymmetry parameter.

! Input variables from model

!  m1               : number of vertical levels in model grid
!  ncat             : max number of hydrometeor categories [=7]
!  mcat             : actual number of hydrometeor categories [= 0, 1, or 7]
!  nhcat            : number of hydrometeor categories including ice habits [=15]
!  iswrtyp          : shortwave radiation parameterization selection flag
!  ilwrtyp          : longwave radiation parameterization selection flag
!  glat             : latitude
!  rtgt             : terrain-following coordinate metric factor
!  topt             : topography height
!  albedt          : surface albeDO
!  cosz             : solar zenith angle
!  rlongup          : upward longwave radiation at surface (W/m^2)
!  rshort           : DOwnward shortwave radiation at surface (W/m^2)
!  rlong            : DOwnward longwave radiation at surface (W/m^2)
!  jnmb (ncat)      : microphysics category flag
!  dnfac (nhcat)    : factor for computing dn from emb
!  pwmasi (nhcat)   : inverse of mass power law exponent for hydrometeors
!  zm (m1)          : model physical heights of W points (m)
!  zt (m1)          : model physical heights of T points (m)
!  press (nzpmax)   : model pressure (Pa)
!  tair (nzpmax)    : model temperature (K)
!  rv (m1)          : model vapor mixing ratio (kg/kg)
!  dn0 (m1)         : model air density (kg/m^3)
!  fthrd (m1)       : theta_il tENDency from radiation
!  jhcat (nzpmax,ncat)  : microphysics category array
!  cx (nzpmax,ncat) : hydrometeor number concentration (#/kg)
!  emb (nzpmax,ncat): hydrometeor mean mass (kg)

! Variables input from model scratch space (redefined internally on each CALL)

!  zml (nrad)       : physical heights of W points of all radiation levels (m)
!  ztl (nrad)       : physical heights of T points of all radiation levels (m)
!  dzl (nrad)       : delta-z (m) of all radiation levels
!  pl (nrad)        : pressure (Pa)
!  tl (nrad)        : temperature (K)
!  dl (nrad)        : air density of all radiation levels (kg/m^3)
!  rl (nrad)        : vapor density of all radiation levels (kg/m^3)
!  vp (nrad)        : vapor pressure (Pa)
!  o3l (nrad)       : stores the calculated ozone profile (g/m^3)
!  flxu (nrad)      : Total upwelling flux (W/m^2)
!  flxd (nrad)      : Total DOwnwelling flux (W/m^2)
!  t (nrad)         : layer transmission function
!  r (nrad)         : layer reflection function
!  tc (nrad)        : cumulative optical depth
!  sigu (nrad)      : upwelling layer source function
!  sigd (nrad)      : DOwnwelling layer source function
!  re (nrad)        : cumulative reflection function
!  vd (nrad)        : multi-scat. dIFfUSE DOwnwelling contributions
!                         from source functions
!  td (nrad)        : inverse of cumulative transmission fnct
!  vu (nrad)        : multi-scat. dIFfUSE upwelling contributions
!                         from source functions
!  tg (nrad)        : gaseous optical depth
!  tcr (nrad)       : continuum/Rayleigh optical depth
!  src (nrad)       : Planck function source for each band
!  fu (nrad*6)      : upwelling fluxes for pseuDO-bands (W/m^2)
!  fd (nrad*6)      : DOwnwelling fluxes for pseuDO-bands (W/m^2)
!  u (nrad*mg)      : path-length for gases (H_2O, CO_2, O_3)  (Pa)
!  tp (nrad*mb)     : optical depth of hydrometeors (m^-1)
!  omgp (nrad*mb)   : Single scatter albeDO of hydrometeors
!  gp (nrad*mb)     : Asymmetry factor of hydrometeors

! LoCALLy-defined variables

!  ngass (mg)       : flags indicating IF H20, CO2, O3 are active for solar wavelengths
!  ngast (mg)       : flags indicating IF H20, CO2, O3 are active for long wavelengths
!  prsnz,prsnzp     : pressure in top two model reference state levels

! Additional variables USEd ONLY within SUBROUTINE mclatchy:
! ==================================================================

!  namax            : maximum allowed number of added rad levels above model top[=10]
!                       USEd for oc and bc coefficients [=13]
! mcdat (33,9,6)    : Mclatchy sounding DATA (33 levels, 9 soundings, 6 vars)
! mclat (33,9,6)    : mcdat interpolated by season to latitude bands
! mcol (33,6)       : mclat interpolated to lat-lon of grid column

! Additional variables USEd ONLY within SUBROUTINE cloud_opt:
! ==================================================================

!  ib .......... band number
!  dn .......... characteristic diameter (m)
!  oc .......... scattering albeDO fit coefficients
!  bc .......... extinction fit coefficients
!  gc .......... asymmetery fit coefficients
!  kradcat ..... cross-reference table giving Jerry's 13 hydrometeor category
!                   numbers as a function of 15 microphysics category numbers

! Particle Numbers describe the following particles:

!     Harrington radiation code             RAMS microphysics
! ----------------------------------------------------------------
!  1:   cloud drops                 1.  cloud drops
!  2:   rain                        2.  rain
!  3:   pristine ice columns        3.  pristine ice columns
!  4:   pristine ice rosettes       4.  snow columns
!  5:   pristine ice plates         5.  aggregates
!  6:   snow columns                6.  graupel
!  7:   snow rosettes               7.  hail
!  8:   snow plates                 8.  pristine ice hexagonal plates
!  9:   aggregates columns          9.  pristine ice dENDrites
!  10:  aggregates rosettes        10.  pristine ice needles
!  11:  aggregates plates          11.  pristine ice rosettes
!  12:  graupel                    12.  snow hexagonal plates
!  13:  hail                       13.  snow dENDrites
!                                  14.  snow needles
!                                  15.  snow rosettes

! for the asymmetery parameter, since we ONLY have spherical
! particles, there are ONLY 7 particle types...

!  1:   cloud drops
!  2:   rain
!  3:   pristine ice
!  4:   snow
!  5:   aggregates
!  6:   graupel
!  7:   hail

!******************************************************************************

SUBROUTINE radcalc3(m1,maxnzp,mcat,iswrtyp,ilwrtyp,IF_adap,lpw  &
     ,glat,rtgt,topt,albedt,cosz,rlongup,rshort,rlong  &
     ,zm,zt,rv,dn0,fthrd,i,j,time,ngrid)

  USE rconstants, ONLY: cp
  USE rrad3     , ONLY: mg,namax,mb,nrad,narad,nb,ocoef,bcoef,gcoef,ncog,ncb, &
                        npartob,npartg,ng,nsolb,npsb,xp,alpha,beta,wght,prf, &
                        trf,ralcs,solar1,ulim,nuum,a0,a1,a2,a3,exptabc
  USE micphys   , ONLY: press,tair

  IMPLICIT NONE

  INTEGER :: m1,maxnzp,mcat,ilwrtyp,iswrtyp,IF_adap,lpw,ngrid
  INTEGER :: i,j,k,kk,koff,nradmax  !,ik,ib,ig

  ! ORIGINAL
  !INTEGER, SAVE :: nCALL = 0
  ! ModIFied by Alvaro L.FazENDa:
  ! SAVE in line 734 automatiCALLy SAVE all variables
  INTEGER :: ncall = 0

  INTEGER ngass(mg),ngast(mg)
  REAL prsnz,prsnzp

  REAL glat,rtgt,topt,cosz,albedt,rlongup,rshort,rlong,slr,time  !,rmix
  REAL eps   !,dzl9,rvk1,rvk0

  REAL zm(m1),zt(m1),dn0(m1),rv(m1),fthrd(m1)

  REAL, allocatable, DIMENSION(:)     :: zml,ztl,dzl,pl,tl,dl,rl,o3l  &
       ,vp,flxu,flxd,tg,tcr,src,t,r,tc  &
       ,sigu,sigd,re,vd,td,vu  &
       ,u,fu,fd,tp,omgp,gp

  DATA eps/1.e-15/

  !     one can choose the gases of importance here,
  !       ngas = 1    gas active
  !            = 0    gas not active
  !
  !       ngas(1) = H2O
  !       ngas(2) = CO2
  !       ngas(3) =  O3

  DATA ngass/1, 1, 1/,ngast/1, 1, 1/
  SAVE

  IF (nCALL == 0) THEN
     nCALL = 1
     nradmax = maxnzp + namax
     ALLOCATE(zml(nradmax),ztl (nradmax),dzl (nradmax),pl (nradmax)  &
          ,tl (nradmax),dl  (nradmax),rl  (nradmax),o3l(nradmax)  &
          ,vp (nradmax),flxu(nradmax),flxd(nradmax),tg (nradmax)  &
          ,tcr(nradmax),src (nradmax),t   (nradmax),r  (nradmax)  &
          ,tc (nradmax),sigu(nradmax),sigd(nradmax),re (nradmax)  &
          ,vd (nradmax),td  (nradmax),vu  (nradmax)               )
     ALLOCATE(u(nradmax*mg),fu(nradmax*6),fd(nradmax*6))
     ALLOCATE(tp(nradmax*mb),omgp(nradmax*mb),gp(nradmax*mb))

     !**ALF: Putting zero on all variables ALLOCATEd:
!!$     CALL azero5(nradmax, zml, ztl, dzl, pl, tl)
!!$     CALL azero5(nradmax, dl, rl, o3l, vp, flxu)
!!$     CALL azero5(nradmax, flxd, tg, tcr, src, t)
!!$     CALL azero5(nradmax, r, tc, sigu, sigd, re)
!!$     CALL azero3(nradmax, vd, td, vu)
!!$     CALL azero(nradmax*mg, u)
!!$     CALL azero2(nradmax*6, fu, fd)
!!$     CALL azero3(nradmax*mb, tp, omgp, gp)
     zml  = 0.
     ztl  = 0.
     dzl  = 0.
     pl   = 0.
     tl   = 0.
     dl   = 0.
     rl   = 0.
     o3l  = 0.
     vp   = 0.
     flxu = 0.
     flxd = 0.
     tg   = 0.
     tcr  = 0.
     src  = 0.
     t    = 0.
     r    = 0.
     tc   = 0.
     sigu = 0.
     sigd = 0.
     re   = 0.
     vd   = 0.
     td   = 0.
     vu   = 0.
     u    = 0.
     fu   = 0.
     fd   = 0.
     tp   = 0.
     omgp = 0.
     gp   = 0.

  END IF


  koff = lpw - 2
  nrad = m1 - 1 + narad - koff

  CALL mclatchy(3,m1,IF_adap,koff  &
       ,prsnz,prsnzp,glat,rtgt,topt,rlongup  &
       ,zm,zt,press,tair,dn0,rv,zml,ztl,pl,tl,dl,rl,o3l,dzl)

  ! zero out scratch arrays
!!$  CALL azero(nrad*mg,u)
!!$  CALL azero(nrad*6,fu)
!!$  CALL azero(nrad*6,fd)
!!$  CALL azero(nrad*mb,tp)
!!$  CALL azero(nrad*mb,omgp)
!!$  CALL azero(nrad*mb,gp)
  u    = 0.
  fu   = 0.
  fd   = 0.
  tp   = 0.
  omgp = 0.
  gp   = 0.

  ! Compute hydrometeor optical properties

  CALL cloud_opt(mb,nb,nrad,m1,koff,mcat,dzl  &
       ,dn0,tp,omgp,gp &
       ,ocoef,bcoef,gcoef,ncog,ncb,npartob,npartg,i,j,time)

  ! Get the path lengths for the various gases...

  CALL path_lengths(nrad,u,rl,dzl,dl,o3l,vp,pl,eps)

  DO k = 1,nrad
     IF (rl(k) .lt.   0. .or.  &
          dl(k) .lt.   0. .or.  &
          pl(k) .lt.   0. .or.  &
          o3l(k) .lt.   0. .or.  &
          tl(k) .lt. 183.) THEN

        PRINT*, 'Temperature too low or negative value of'
        PRINT*, 'density, vapor, pressure, or ozone'
        PRINT*, 'when CALLing Harrington radiation'
        PRINT*, 'at k,i,j = ',k,i,j,'   ngrid=',ngrid
        PRINT*, 'STOPping model'
        PRINT*, 'rad: rl(k), dl(k), pl(k), o3l(k), tl(k)'
        DO kk=1,nrad
           WRITE (*,'(5g15.6)') rl(kk), dl(kk), pl(kk), o3l(kk), tl(kk)
        END DO
        STOP 'STOP: radiation CALL'
     END IF
  END DO

  ! CALL shortwave and longwave schemes...

  IF (iswrtyp .eq. 3 .and. cosz .gt. 0.03) THEN
!!$     CALL azero2(nrad,flxu,flxd)
     flxu = 0.
     flxd = 0.

     CALL swrad(nrad,ng,nb,nsolb,npsb,   &         !  counters
          u,pl,tl,dzl,vp,                  &      !  model variables
          xp,alpha,beta,wght,prf,trf,ralcs,  &   !  band specIFics
          solar1,ngass,                      &    !        "
          albedt,slr,cosz,               & !  boundaries
          tg,tcr,tp,omgp,gp,   &!  loCALLy defined
          t,r,tc,                              & !  for fluxes
          sigu,sigd,re,vd,td,vu,fu,fd,         & !       "
          flxu,flxd,ulim,i,time)                     !  sw fluxes

     rshort = flxd(1)

!--(DMK-CCATT-INI)-----------------------------------------------------
     do k = lpw,m1-1
       kk = k - koff
       fthrd(k) = fthrd(k)  &
          + (flxd(kk) - flxd(kk-1) + flxu(kk-1) - flxu(kk)) &
             / (dl(kk) * dzl(kk) * cp)
     enddo
!--(DMK-CCATT-OLD)-----------------------------------------------------
!     DO k = lpw,m1-1
!        fthrd(k+koff) = fthrd(k+koff)  &
!             + (flxd(k) - flxd(k-1) + flxu(k-1) - flxu(k)) / (dl(k) * dzl(k) * cp)
!     END DO
!--(DMK-CCATT-FIM)-----------------------------------------------------

  ELSE
     rshort = 0.
  END IF

  IF (ilwrtyp .eq. 3) THEN
!!$     CALL azero2(nrad,flxu,flxd)
     flxu = 0.
     flxd = 0.

     CALL lwrad(nrad,ng,nb,nsolb,npsb,nuum,   &    !  counters
          u,pl,tl,dzl,vp,                       & !  model variables
          xp,alpha,beta,wght,prf,trf,ralcs,     &!  band specIFics
          a0,a1,a2,a3,                          &!        "
          exptabc,ngast,                        &!  boundaries
          tg,tcr,tp,omgp,gp,src,   &  !  loCALLy defined
          t,r,tc,                               &!  for fluxes
          sigu,sigd,re,vd,td,vu,fu,fd,          &!       "
          flxu,flxd,ulim)                          !  fluxes, output

     rlong = flxd(1)

!--(DMK-CCATT-INI)-----------------------------------------------------
   do k = lpw,m1-1
      kk = k - koff
      fthrd(k) = fthrd(k)  &
         + (flxd(kk) - flxd(kk-1) + flxu(kk-1) - flxu(kk)) &
            / (dl(kk) * dzl(kk) * cp)
   enddo
!--(DMK-CCATT-OLD)-----------------------------------------------------
!     DO k = lpw,m1-1
!        fthrd(k+koff) = fthrd(k+koff)  &
!             + (flxd(k) - flxd(k-1) + flxu(k-1) - flxu(k)) / (dl(k) * dzl(k) * cp)
!     END DO
!--(DMK-CCATT-FIM)-----------------------------------------------------

  END IF

  RETURN
END SUBROUTINE radcalc3

!******************************************************************************

SUBROUTINE cloud_opt(mb,nb,nrad,m1,koff,mcat,dzl  &
     ,dn0,tp,omgp,gp,oc,bc,gc,ncog,ncb,npartob,npartg,i,j,time)

  ! computing properties of spherical liquid water and irregular ice
  ! using fits to adt theory
  !
  ! ib .......... band number
  ! mb .......... maximum number of bands
  ! nb .......... total number of bands
  ! m1 .......... number of vertical levels
  ! dzl .......... delta z in each level (m)
  ! dn .......... characteristic diameter (m)
  ! emb ......... mean hydrometeor mass (kg)
  ! cx .......... hydrometeor concentration (#/kg)
  ! tp .......... optical depth
  ! omgp ........ scattering albeDO
  ! gp .......... asymmetry parameter
  ! oc .......... scattering albeDO fit coefficients
  ! bc .......... extinction fit coefficients
  ! gc .......... asymmetry fit coefficients
  ! ncog ........ number of fit coefficients (omega and asym)
  ! ncb ......... number of fit coefficients (extinction)
  ! kradcat ..... cross-reference table giving Jerry's 13 hydrometeor category
  !                 numbers as a function of 15 microphysics category numbers
  ! dn0 ......... model air density (kg/m^3)
  ! dnfac ....... factor for computing dn from emb
  ! pwmasi ...... inverse of power USEd in mass power law
  ! npartob ..... number of hydrometeor categories (including dIFferent habits)
  !                 USEd for oc and bc coefficients
  ! npartg ...... number of hydrometeor categories USEd for gc coefficients
  !

  USE micphys, ONLY: jnmb,cx,jhcat,dnfac,emb,pwmasi

  IMPLICIT NONE

  INTEGER mb,nb,ib,nrad,m1,ncog,ncb,krc,npartob,npartg,koff  !,iz
  INTEGER icat,mcat,k,i,j,ihcat

  INTEGER kradcat(15)
  REAL dzl(nrad),tp(nrad,mb),omgp(nrad,mb),gp(nrad,mb),dn0(m1)
  REAL oc(ncog,mb,npartob),bc(ncb,mb,npartob),gc(ncog,mb,npartg)
  REAL ext,om,gg,dn,time

  REAL dnmin(7),dnmax(7)
  DATA dnmin /   1.,   10.,   1.,  125.,   10.,   10.,   10./
  DATA dnmax /1000.,10000., 125.,10000.,10000.,10000.,10000./

  DATA kradcat/1,2,3,6,10,12,13,5,5,3,4,8,8,6,7/

  DO icat = 1,mcat
     IF (jnmb(icat) .gt. 0) THEN

        DO k = 2,m1-1-koff

           IF (cx(k+koff,icat) .gt. 1.e-9) THEN
              ihcat = jhcat(k+koff,icat)
              krc = kradcat(ihcat)
              dn = dnfac(ihcat) * emb(k+koff,icat) ** pwmasi(ihcat) * 1.e6
              dn = max(dnmin(icat),min(dnmax(icat),dn))

              DO ib = 1,nb

                 ext = cx(k+koff,icat) * dn0(k+koff) * dzl(k)  &
                      * bc(1,ib,krc) * dn ** bc(2,ib,krc)
                 om = oc(1,ib,krc)  &
                      + oc(2,ib,krc) * exp(oc(3,ib,krc) * dn)  &
                      + oc(4,ib,krc) * exp(oc(5,ib,krc) * dn)
                 gg = gc(1,ib,icat)  &
                      + gc(2,ib,icat) * exp(gc(3,ib,icat) * dn)  &
                      + gc(4,ib,icat) * exp(gc(5,ib,icat) * dn)


                 tp(k,ib) = tp(k,ib) + ext

                 omgp(k,ib) = omgp(k,ib) + om * ext
                 gp(k,ib) = gp(k,ib) + gg * om * ext

              END DO

           END IF
        END DO

     END IF
  END DO

  ! Combine the optical properties....

  DO ib = 1,nb
     DO k = 2,m1-1-koff
        IF (tp(k,ib) .gt. 0.0) THEN
           gp(k,ib) = gp(k,ib) / omgp(k,ib)
           omgp(k,ib) = omgp(k,ib) / tp(k,ib)
        ELSE
           omgp(k,ib) = 0.0
           gp(k,ib) = 0.0
        END IF

        ! Check for validity of opt values before CALLing radiation
        !      IF (tp(k,ib) .lt. 0) THEN
        !         PRINT*, 'tp(k,ib) less than zero for k,ib = ',k,ib
        !         PRINT*, 'tp(k,ib) = ',tp(k,ib)
        !         STOP 'opt1'
        !      END IF
        !      IF (omgp(k,ib) .lt. 0. .or. omgp(k,ib) .gt. 1.) THEN
        !         PRINT*, 'omgp(k,ib) out of range [0,1] for k,ib = ',k,ib
        !         PRINT*, 'omgp(k,ib) = ',omgp(k,ib)
        !         STOP 'opt2'
        !      END IF
        !      IF (gp(k,ib) .lt. 0. .or. gp(k,ib) .gt. 1.) THEN
        !         PRINT*, 'gp(k,ib) out of range [0,1] for k,ib = ',k,ib
        !         PRINT*, 'gp(k,ib) = ',gp(k,ib)
        !         STOP 'opt3'
        !      END IF

     END DO
  END DO

  RETURN
END SUBROUTINE cloud_opt

!------------------------------------------------------------------

SUBROUTINE path_lengths(nrad,u,rl,dzl,dl,o3l,vp,pl,eps)

  ! Get the path lengths for the various gases...

  IMPLICIT NONE
  INTEGER :: nrad
  REAL, DIMENSION(nrad) :: rl,dzl,dl,o3l,vp,pl
  REAL, DIMENSION(nrad,3) :: u
  REAL :: rvk0,rvk1,dzl9,rmix,eps
  INTEGER :: k

  u(1,1) = .5 * (rl(2) + rl(1)) * 9.81 * dzl(1)
  u(1,2) = .5 * (dl(2) + dl(1)) * .45575e-3 * 9.81 * dzl(1)
  u(1,3) = o3l(1) * 9.81 * dzl(1)

  rvk0 = rl(1)
  DO k = 2,nrad
     rvk1 = (rl(k) + 1.e-6)
     dzl9 = 9.81 * dzl(k)
     rmix = rvk1 / dl(k)
     vp(k) = pl(k) * rmix / (.622 + rmix)
     u(k,1) = (rvk1 - rvk0) / (log(rvk1 / rvk0) + eps) * dzl9
     u(k,2) = (dl(k) - dl(k-1)) / (log(dl(k) / dl(k-1)) + eps)  &
          * dzl9 * 0.45575e-3
     u(k,3) = 0.5 * dzl9 * (o3l(k) + o3l(k-1))
     rvk0 = rvk1
  END DO

  RETURN
END SUBROUTINE path_lengths

SUBROUTINE radCompUkmo(ngrids,ia,iz,ja,jz,mxp,myp,mzp,ngrid,swFile,lwFile)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: ia,iz,ja,jz,mxp,myp,mzp,ngrid,ngrids
   CHARACTER(LEN=*),INTENT(IN) :: swFile,lwFile

   INTEGER :: ncols
   INTEGER :: kmax
   
   kmax=mzp-1

   ncols=(iz-ia+1)*(jz-ja+1)
   CALL prepRadUKMO_new(ngrids,ngrid,ncols,ia,iz,ja,jz,kMax,swFile,lwFile,mxp,myp)	   
   CALL UKMO_Lw_new(ngrid,kMax,mxp,myp,ncols)
   CALL UKMO_Sw_new(ngrid,kMax,mxp,myp,nCols)

END SUBROUTINE radCompUkmo

SUBROUTINE prepRadUKMO_new(ngrids,ngrid,ncols,ia,iz,ja,jz,kMax,swfile,lwfile,mxp,myp)
  USE ModDateUtils
   USE mem_radiate , ONLY: radiate_g
   USE mem_basic   , ONLY: basic_g
   USE rconstants  , ONLY : cp,cpor,p00,stefan
   USE rconstants  , ONLY: cpi,p00,cpor
   USE micphys     , ONLY: level,pwmas,pwmasi,cparm,emb0,emb1,gnu,dnfac,jhcat, &
                          emb,cx
   USE mem_leaf    , ONLY: leaf_g ! INTENT(IN)
   USE mem_grid    , ONLY : iyear1, imonth1, idate1, itime1, grid_g, npatch
   USE UkmoAdapt   !, ONLY: getOz,sig,mem_UKMO(ngrid)%sigMid,mem_UKMO,sigmid,sig,flip 
   USE ref_sounding, ONLY: pi01dn ! INTENT(IN)
   USE Rad_UKMO    ,ONLY: InitRadUKMO
   USE micphys, ONLY: &
          gnu, level, icloud, irain, ipris, & ! INTENT(IN)
         isnow, iaggr, igraup, ihail         ! INTENT(IN)
   USE mem_micro, ONLY: micro_g
   USE mem_cuparm, ONLY: cuparm_g, nnqparm 
   IMPLICIT NONE
                
   REAL(KIND=r8), PARAMETER :: pi=3.1415926e00_r8
   REAL(KIND=r8), PARAMETER :: solar=1.3533e03_r8 !* 100.0 !test
   REAL(KIND=r8), PARAMETER :: adjust=0.01 !to convert from Pascal to mbar

   INTEGER,INTENT(IN)  :: ngrids,ngrid,ia,iz,ja,jz,kMax,ncols,mxp,myp
   CHARACTER(LEN=*),INTENT(IN) :: swfile
   CHARACTER(LEN=*),INTENT(IN) :: lwfile
     
   REAL(KIND=r8) :: colrad(ncols)
   REAL(KIND=r8) :: press(ncols,kMax+1)

   REAL(KIND=r8) :: tair(ncols,kMax+1)   
   REAL(KIND=r8) :: qe(ncols,kMax+1)    
   REAL(KIND=r8) :: clwp (ncols,kMax)   ! Cloud Liquid Water Path
   REAL(KIND=r8) :: rei  (ncols,kMax)   ! Ice particle Effective Radius (microns)
   REAL(KIND=r8) :: rel  (ncols,kMax)   ! Liquid particle Effective Radius (microns)
   REAL(KIND=r8) :: fice (ncols,kMax)   ! fractional amount of cloud that is ice
   REAL(KIND=r8) :: lmixr(ncols,kMax)   ! ice/water mixing ratio
   REAL(KIND=r8) :: taud (ncols,kMax)
   REAL(KIND=r8) :: date,psurf(ncols),hcpi,prsnz
   INTEGER :: itime(4)
   INTEGER :: i,j,k,m,np,ncol,icount
!   INTEGER :: julday
   REAL(KIND=r8) :: picpi
   REAL(KIND=r8) :: weight

   IF(.not. ukmo_created) THEN
      !PRINT *,'Init Rad_ukmo 1st step with kMax=',kMax;call flush(6)
      CALL Create_mem_ukmo(ngrids)
      ukmo_created=.true.
      CALL Alloc_mem_Rad_UKMO(ngrid,ncols,kMax)
      CALL CreateFlip(kMax,ngrid)!-1)
   END IF
    
   itime(1)=itime1
   itime(2)=imonth1
   itime(3)=idate1
   itime(4)=iyear1
   CALL getco2(itime,mem_UKMO(ngrid)%co2val)
   mem_UKMO(ngrid)%solar=solar
   
   !co2val=0.0_r8
   date=julday(imonth1, idate1, iyear1)
   

   !Pressure at top
   prsnz  = (pi01dn(kMax+1,1)/cp)**cpor*p00
   ncol=0
   DO i=ia,iz
      DO j=ja,jz
         ncol=ncol+1

         !Mapping cols in i,j and vice-versa
         mem_UKMO(ngrid)%MapI(ncol)=i
         mem_UKMO(ngrid)%MapJ(ncol)=j

         !Cos Zenithal and Surface Temperature
         mem_UKMO(ngrid)%coszLocal(ncol)=radiate_g(ngrid)%cosz(i,j)
         mem_UKMO(ngrid)%gtg(ncol)= &
                (radiate_g(ngrid)%rlongup(i,j)/stefan)** 0.25e0_r8

         !Pressure, air temperature and Moisture
         DO k = 1,kMax+1
            picpi = (basic_g(ngrid)%pi0(k,i,j) + basic_g(ngrid)%pp(k,i,j)) * cpi
            press(ncol,k) = (p00 * picpi ** cpor)!/100.0
            tair(ncol,k) = basic_g(ngrid)%theta(k,i,j) * picpi
            qe(ncol,k)=basic_g(ngrid)%rv(k,i,j)
         END DO         

         !surface pressure
         psurf(ncol)=(basic_g(ngrid)%pi0(1,i,j) + basic_g(ngrid)%pp(1,i,j)+ &
                basic_g(ngrid)%pi0(2,i,j) + basic_g(ngrid)%pp(2,i,j))*cpi*0.50
         psurf(ncol)=(p00*psurf(ncol)**cpor)*adjust
         !Flipping the variables
         mem_UKMO(ngrid)%FlipPbot(ncol,kMax)=press(ncol,1)
         DO k=1,kMax
            mem_UKMO(ngrid)%FlipTe(ncol,mem_UKMO(ngrid)%flip(k))=tair(ncol,k)
            mem_UKMO(ngrid)%FlipQe(ncol,mem_UKMO(ngrid)%flip(k))=qe(ncol,k)
            mem_UKMO(ngrid)%FlipPbot(ncol,k)=press(ncol,mem_UKMO(ngrid)%flip(k))
         END DO

         !Delta Pressure 
         mem_UKMO(ngrid)%FlipDP(ncol,1)=mem_UKMO(ngrid)%FlipPbot(ncol,1)-prsnz
         mem_UKMO(ngrid)%FlipPMid(ncol,1)=(mem_UKMO(ngrid)%FlipPbot(ncol,1)+ &
                                           prsnz)*0.50
                                    
         DO k=2,kMax
            mem_UKMO(ngrid)%FlipPMid(ncol,k)= &
                 (mem_UKMO(ngrid)%FlipPbot(ncol,k)+ &
                  mem_UKMO(ngrid)%FlipPbot(ncol,k-1))*0.50
            mem_UKMO(ngrid)%FlipDP(ncol,k)= &
                  mem_UKMO(ngrid)%FlipPbot(ncol,k)- &
                  mem_UKMO(ngrid)%FlipPbot(ncol,k-1)
         END DO
         
         !Land/Ocean mask
         DO np=1,npatch
            IF(leaf_g(ngrid)%leaf_class(i,j,np)==0) THEN
               mem_UKMO(ngrid)%imask(ncol)=0
               cycle
            END IF !Ocean
            IF(leaf_g(ngrid)%leaf_class(i,j,np)==2) THEN
               mem_UKMO(ngrid)%imask(ncol)=13
               cycle
            END IF!land ice         
         END DO
         !Albedo         
         mem_UKMO(ngrid)%AlbVisDiff(ncol)=radiate_g(ngrid)%albedt(i,j)
         mem_UKMO(ngrid)%AlbNirDiff(ncol)=radiate_g(ngrid)%albedt(i,j)
         mem_UKMO(ngrid)%AlbVisBeam(ncol)=radiate_g(ngrid)%albedt(i,j)
         mem_UKMO(ngrid)%AlbNirBeam(ncol)=radiate_g(ngrid)%albedt(i,j)         
         !colatitude - 0 to Pi NPole->SPole
         colrad(ncol)=pi*(1.00e0_r8-(grid_g(ngrid)%glat(i,j)+90.00e0_r8)/180.00e0_r8)
         
        END DO
     END DO

     IF(mem_UKMO(ngrid)%radUkmoNotInitiated) THEN !Doing sigma levels from pressure
                                  ! Initializing from rad files and initializing
                                  ! climatologic ozone

        !Calculating sig and delta sig
         DO i=1,kMax
            mem_UKMO(ngrid)%sig(i)=mem_UKMO(ngrid)%flippbot(15,kMax-i+1)/ &
               mem_UKMO(ngrid)%flippbot(15,kMax)
         END DO
         mem_UKMO(ngrid)%sig(1)=1.
         mem_UKMO(ngrid)%sig(kMax+1)=0.
         DO i=1,kMax
             IF (mem_UKMO(ngrid)%sig(i)<0.1) EXIT
         END DO
         mem_UKMO(ngrid)%nls=kMax-i+1
         DO k=1,kMax
            mem_UKMO(ngrid)%delsig(k)=mem_UKMO(ngrid)%sig(k)-mem_UKMO(ngrid)%sig(k+1)
         END DO
         DO K=1,kMax-1 !-1
             mem_UKMO(ngrid)%sigMid(k)=(mem_UKMO(ngrid)%sig(k)+mem_UKMO(ngrid)%sig(k+1))/2
         END DO

        CALL InitRadUKMO(kMax,mem_UKMO(ngrid)%sig,mem_UKMO(ngrid)%nls,swfile,lwfile)
        CALL InitGetoz(yrl,kMax,mem_UKMO(ngrid)%sig)
         mem_UKMO(ngrid)%radUkmoNotInitiated=.false.
     END IF 

     !Getting O3 from four season climatological values   
     CALL getoz (ncols,Ncols,kMax,mem_UKMO(ngrid)%sigMid,colrad,date,mem_UKMO(ngrid)%flipO3)

     !Adjust pressure values
     mem_UKMO(ngrid)%FlipPbot=mem_UKMO(ngrid)%FlipPbot*Adjust                                   
     mem_UKMO(ngrid)%FlipDP=  mem_UKMO(ngrid)%FlipDP  *Adjust                   
     mem_UKMO(ngrid)%FlipPMid=mem_UKMO(ngrid)%FlipPMid*Adjust
     prsnz=prsnz*adjust
          
     !cloud optical properties CCM3 - Slingo (1989)
     !Getting fice, rei, rel and lmixR
     CALL Cloud_Micro_CCM3(&
             Ncols, kMax , mem_UKMO(ngrid)%sigMid(1:kMax), mem_UKMO(ngrid)%sig(1:kMax) &
             , mem_UKMO(ngrid)%delsig, &
             mem_UKMO(ngrid)%imask   , &
             psurf   , Tair   , Qe    , mem_UKMO(ngrid)%gtg  , &
             mem_UKMO(ngrid)%FlipPBot, prsnz, &
             clwp , lmixr, fice  , rei   , rel   , taud  )
    
     !Flipping optical properties
     !LFR: Note - flip is an array that invert levels
     DO ncol=1,NCols
        DO k=1,kMax
           mem_UKMO(ngrid)%FlipRei  (ncol,mem_UKMO(ngrid)%flip(k))=Rei  (ncol,k)
           mem_UKMO(ngrid)%FlipRel  (ncol,mem_UKMO(ngrid)%flip(k))=Rel  (ncol,k)
           mem_UKMO(ngrid)%FlipFice (ncol,mem_UKMO(ngrid)%flip(k))=Fice (ncol,k)
           !FlipTaud (flip(k))=Taud (k)
           mem_UKMO(ngrid)%FlipLMixR(ncol,mem_UKMO(ngrid)%flip(k))=LMixR(ncol,k)
        END DO
     END DO
     
     !LFR - Just for test - clear sky + other fixed setups
     !mem_UKMO(ngrid)%FlipLMixR=0.0
     !mem_UKMO(ngrid)%flipO3=0.0
     !mem_UKMO(ngrid)%co2val=0.0
     !mem_UKMO(ngrid)%FlipFice=1.0
     !mem_UKMO(ngrid)%FlipRei=10.0
     !mem_UKMO(ngrid)%FlipRel=5.0

        
END SUBROUTINE prepRadUKMO_new

SUBROUTINE UKMO_sw_new(ngrid,kMax,mxp,myp,ncols)
   USE Rad_UKMO    , ONLY: ukmo_swintf
   USE UkmoAdapt, ONLY: mem_UKMO
   USE mem_radiate, ONLY: &
         radiate_g       ! INTENT(INOUT)
   IMPLICIT NONE

   INTEGER,INTENT(IN) :: ngrid,ncols,kMax,mxp,myp
   !REAL,INTENT(INOUT)   :: rshort(mxp,myp)
    
    INTEGER :: i,j,ncol
    INTEGER :: irec=1

    !rshort=0.0
    !Calling short wave radiation
    CALL ukmo_swintf( &
            ! Model Info and flags
            mem_UKMO(ngrid)%nls, ncols, kMax,         &
            ! Solar field
            mem_UKMO(ngrid)%coszLocal, &
            mem_UKMO(ngrid)%solar, &
            ! Atmospheric fields
            mem_UKMO(ngrid)%FlipPBot, &
            mem_UKMO(ngrid)%FlipPMid, &
            mem_UKMO(ngrid)%FlipDP,   &
            mem_UKMO(ngrid)%FlipTe, &
            mem_UKMO(ngrid)%FlipQe, &
            mem_UKMO(ngrid)%FlipO3, &
            mem_UKMO(ngrid)%co2val, &
            mem_UKMO(ngrid)%gtg, &
            ! SURFACE:  albedo
            mem_UKMO(ngrid)%imask       ,&
            mem_UKMO(ngrid)%AlbVisDiff, &
            mem_UKMO(ngrid)%AlbNirDiff, &
            mem_UKMO(ngrid)%AlbVisBeam, &
            mem_UKMO(ngrid)%AlbNirBeam, & 
            mem_UKMO(ngrid)%ToaDown        , &
            mem_UKMO(ngrid)%SfcVisBeam     , &
            mem_UKMO(ngrid)%SfcVisDiff     , &
            mem_UKMO(ngrid)%SfcNirBeam     , &
            mem_UKMO(ngrid)%SfcNirDiff     , &
            mem_UKMO(ngrid)%SfcVisBeamC    , &
            mem_UKMO(ngrid)%SfcVisDiffC    , &
            mem_UKMO(ngrid)%SfcNirBeamC    , &
            mem_UKMO(ngrid)%SfcNirDiffC    , &
            mem_UKMO(ngrid)%ToaNetC        , &
            mem_UKMO(ngrid)%ToaNet         , &
            mem_UKMO(ngrid)%SfcNet         , &
            mem_UKMO(ngrid)%SfcNetC           , &
            mem_UKMO(ngrid)%HeatRateC      , &
            mem_UKMO(ngrid)%HeatRate       , &          
            ! Cloud field
            mem_UKMO(ngrid)%cld, &
            mem_UKMO(ngrid)%clu,  &
            mem_UKMO(ngrid)%FlipFice, &
            mem_UKMO(ngrid)%FlipRei, &
            mem_UKMO(ngrid)%FlipRel, &
            mem_UKMO(ngrid)%FlipLMixR)           

    DO ncol=1,nCols
       i=mem_UKMO(ngrid)%MapI(ncol)
       j=mem_UKMO(ngrid)%MapJ(ncol)
       radiate_g(ngrid)%rshort(i,j)= mem_UKMO(ngrid)%SfcVisBeam(ncol)+ &
                    mem_UKMO(ngrid)%SfcVisDiff(ncol)+ &
                    mem_UKMO(ngrid)%SfcNirBeam(ncol)+ &
                    mem_UKMO(ngrid)%SfcNirDiff(ncol)
    END DO
        
END SUBROUTINE UKMO_sw_new

SUBROUTINE UKMO_lw_new(ngrid,kMax,mxp,myp,ncols)
   USE Rad_UKMO    , ONLY: ukmo_lwintf
   USE UkmoAdapt, ONLY: mem_UKMO
   USE mem_radiate, ONLY: &
         radiate_g       ! INTENT(INOUT)
   IMPLICIT NONE

   INTEGER,INTENT(IN) :: ngrid,kMax,mxp,myp,ncols
!   REAL,INTENT(OUT)   :: rlong(mxp,myp)
!   REAL,INTENT(OUT)   :: rlongUp(mxp,myp)
    
    INTEGER :: i,j,ncol
    INTEGER :: irec=1

    !rlong=0.0
   !LFR  write(*,'(A3,10(2x,A12))') ' k ','sig(k)','flipPBot','flipPMid','flipDp','flipTe','flipQe','flipO3'
   !LFR  do i=1,kMax
      !LFR  write(*,'(I3,10(2x,E12.5))') i,sig(i),mem_UKMO(ngrid)%flippbot(15,i),mem_UKMO(ngrid)%flippmid(15,i),&
           !LFR  mem_UKMO(ngrid)%flipdp(15,i),mem_UKMO(ngrid)%flipte(15,i),mem_UKMO(ngrid)%flipqe(15,i)*1e3, &
           !LFR  mem_UKMO(ngrid)%flipo3(15,i)
   !LFR  enddo
   !LFR  write(*,'(A3,6(2x,A12))') ' k ','cld','clu','FlipFice','FlipRei','FlipRel','FlipMixR'
   !LFR  do i=1,kMax
      !LFR  write(*,'(I3,6(2x,E12.5))') i,mem_UKMO(ngrid)%cld(15,i),mem_UKMO(ngrid)%clu(15,i), &
           !LFR  mem_UKMO(ngrid)%FlipFice(15,i),mem_UKMO(ngrid)%FlipRei(15,i), &
           !LFR  mem_UKMO(ngrid)%FlipRel(15,i),mem_UKMO(ngrid)%FlipLMixR(15,i)
   !LFR  enddo
   !LFR  write(*,*) mem_UKMO(ngrid)%solar, mem_UKMO(ngrid)%co2val, mem_UKMO(ngrid)%gtg(15), mem_UKMO(ngrid)%imask(15), &
            !LFR  mem_UKMO(ngrid)%AlbVisDiff(15), mem_UKMO(ngrid)%AlbNirDiff(15), &
            !LFR  mem_UKMO(ngrid)%AlbVisBeam(15), mem_UKMO(ngrid)%AlbNirBeam(15), &
            !LFR  acos(mem_UKMO(ngrid)%coszLocal(15))*180/3.141592
  
    !Calling short wave radiation
     CALL ukmo_lwintf( &
            ! Model Info and flags
            mem_UKMO(ngrid)%nls,ncols, kMax,         &
            ! Atmospheric fields
            mem_UKMO(ngrid)%FlipPBot, &
            mem_UKMO(ngrid)%FlipPMid, &
            mem_UKMO(ngrid)%FlipDP, &
            mem_UKMO(ngrid)%FlipTe, &
            mem_UKMO(ngrid)%FlipQe, &
            mem_UKMO(ngrid)%FlipO3, &
            mem_UKMO(ngrid)%co2val, &
            mem_UKMO(ngrid)%gtg, &
            ! SURFACE:  albedo
            mem_UKMO(ngrid)%imask       ,&
            !
            mem_UKMO(ngrid)%lw_toa_up_clr , &
            mem_UKMO(ngrid)%lw_toa_up , &
            mem_UKMO(ngrid)%lw_cool_clr    , &
            mem_UKMO(ngrid)%lw_cool    , &
            mem_UKMO(ngrid)%lw_sfc_net_clr,&
            mem_UKMO(ngrid)%lw_sfc_net,&
            mem_UKMO(ngrid)%lw_sfc_down_clr, &
            mem_UKMO(ngrid)%lw_sfc_down, &
            ! Cloud field
            mem_UKMO(ngrid)%cld,&
            mem_UKMO(ngrid)%clu, &
            mem_UKMO(ngrid)%FlipFice, &
            mem_UKMO(ngrid)%FlipRei,&
            mem_UKMO(ngrid)%FlipRel,&
            mem_UKMO(ngrid)%FlipLMixR)           

    DO ncol=1,ncols
       i=mem_UKMO(ngrid)%MapI(ncol)
       j=mem_UKMO(ngrid)%MapJ(ncol)
       radiate_g(ngrid)%rlong(i,j)=  mem_UKMO(ngrid)%lw_sfc_down_clr(ncol)+&
                    mem_UKMO(ngrid)%lw_sfc_down(ncol)                             
    END DO     
        
END SUBROUTINE UKMO_lw_new

