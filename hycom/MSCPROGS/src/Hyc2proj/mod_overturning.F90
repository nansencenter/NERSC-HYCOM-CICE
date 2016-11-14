module mod_overturning

contains


subroutine mosf(u,v,dp,qlon,qlat,idm,jdm,kdm,mostrf,lats,ds,nlats,nds,undef)
use m_spherdist
implicit none

   integer, intent(in) :: idm,jdm,nlats,nds,kdm
   real,    intent(in) :: qlat(0:idm+1,0:jdm+1),qlon(0:idm+1,0:jdm+1)
   real,    intent(in), dimension(idm,jdm,kdm) :: u,v,dp
   real,    intent(out) :: mostrf(nlats,nds)
   real,    intent(in) :: lats(nlats)
   real,    intent(in) :: ds(nds)
   real,    intent(in) :: undef
      
   real depthu(idm,jdm),depthv(idm,jdm) ,scvx(idm,jdm),scuy(idm,jdm), &
        depth(idm,jdm),fumask(idm,jdm), fvmask(idm,jdm)
   real,dimension(idm,jdm,kdm) :: dpu, dpv ,transu,transv
   real, dimension(idm,jdm) :: uint,vint,dptmp, inttransu,inttransv,ftstmask
   real, dimension(idm,jdm,kdm+1) :: sumdpu,sumdpv
   real, dimension(idm,jdm,nds) :: inttransu2,inttransv2,crssecu,crssecv
   logical :: mask(idm,jdm)
   real :: dptmpu, dptmpv
   real :: crssec(nlats,nds)

   integer :: ideep, ilat, i, j, im1, ip1, jm1, jp1,k 
   logical, parameter :: masktest=.false.
   real :: maxdepth, velresidual, trresidual



!
! --- ------------------------------------------------------------------
! --- calculate scale factors and depths in u- and v-points
! --- ------------------------------------------------------------------
!
   !$OMP PARALLEL DO PRIVATE(i,j)
   do j=1,jdm
   do i=1,idm
      depth(i,j) = sum(dp(i,j,:))
   end do
   end do
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO PRIVATE(i,j,ip1,im1,jp1,jm1)
   do j=1,jdm
     jp1=j+1
     jm1=mod(j+jdm-2,jdm)+1
     do i=1,idm
       im1=max(1,i-1)
       ip1=i+1
       scvx(i,j)=spherdist(qlon(ip1,j),qlat(ip1,j), qlon(i  ,j),qlat(i  ,j))
       scuy(i,j)=spherdist(qlon(i,jp1),qlat(i,jp1), qlon(i,j  ),qlat(i,j  ))
       depthu(i,j)=min(depth(i,j),depth(im1,j))
       depthv(i,j)=min(depth(i,j),depth(i,jm1))
       if (depth(i,j).lt.1.) depth(i,j)=0.
       if (depthu(i,j).lt.1.) depthu(i,j)=0.
       if (depthv(i,j).lt.1.) depthv(i,j)=0.
     enddo
   enddo
   !$OMP END PARALLEL DO


!
! --- ------------------------------------------------------------------
! --- calculate dp in u and v points
! --- ------------------------------------------------------------------
!
   sumdpu=0.
   sumdpv=0.
   do k=1,kdm
   !$OMP PARALLEL DO PRIVATE(i,j,im1,jm1)
   do j=1,jdm
     jm1=mod(j+jdm-2,jdm)+1
     do i=1,idm
       im1=max(1,i-1)
       dpu(i,j,k)   = 0.5*(dp(i,j,k) + dp(im1,j,k))
       dpv(i,j,k)   = 0.5*(dp(i,j,k) + dp(i,jm1,k))

       sumdpu(i,j,k+1)=min(sumdpu(i,j,k)+dpu(i,j,k),depthu(i,j))
       sumdpv(i,j,k+1)=min(sumdpv(i,j,k)+dpv(i,j,k),depthv(i,j))

       transu(i,j,k)=dpu(i,j,k)*u(i,j,k)*scuy(i,j)
       transv(i,j,k)=dpv(i,j,k)*v(i,j,k)*scvx(i,j)
     end do
   end do
   !$OMP END PARALLEL DO
   end do
   !print *,maxval(transu),minval(transu)
   !print *,maxval(transv),minval(transv)


   ! Top-down
   inttransu2=0.
   inttransv2=0.
   crssecu=0.
   crssecv=0.
   do ideep=1,nds
   !print *,ideep,nds,ds(ideep)
   do k=1,kdm
   !$OMP PARALLEL DO PRIVATE(i,j,dptmpu,dptmpv)
   do j=1,jdm
   do i=1,idm
      dptmpu= min(sumdpu(i,j,k+1),ds(ideep)) - min(sumdpu(i,j,k),ds(ideep)) 
      dptmpv= min(sumdpv(i,j,k+1),ds(ideep)) - min(sumdpv(i,j,k),ds(ideep)) 
      inttransu2(i,j,ideep) = inttransu2(i,j,ideep)+dptmpu*transu(i,j,k)/(dpu(i,j,k)+1e-4)
      inttransv2(i,j,ideep) = inttransv2(i,j,ideep)+dptmpv*transv(i,j,k)/(dpv(i,j,k)+1e-4)


      ! For computing "cross section area" later
      crssecu(i,j,ideep)= crssecu(i,j,ideep)+dptmpu*scuy(i,j)
      crssecv(i,j,ideep)= crssecv(i,j,ideep)+dptmpv*scvx(i,j)


   end do
   end do
   !$OMP END PARALLEL DO
   end do
   end do
   !print *,maxval(inttransu2),minval(inttransu2)
   !print *,maxval(inttransv2),minval(inttransv2)


!
! --- ------------------------------------------------------------------
! --- calculate  meridional overturning stream function
! --- ------------------------------------------------------------------
!



   ! For every zonal band:
   mostrf=0.
   crssec=0.
   !oldmask=.false.
   !print *,maxval(dpv(:,:,1)),maxval(dpu(:,:,1)),maxval(dp(:,:,1))
   mask=.false.
   do ilat=1,nlats
      !print *,ilat,nlats,lats(ilat)

      ! Find points inside this region (North of lats(i)
      !mask=.false.
      !where (qlat(1:idm,1:jdm)>lats(ilat))
      !   mask=.true.
      !endwhere
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jdm
      do i=1,idm
         mask(i,j)=qlat(i,j)>lats(ilat)
      end do
      end do
      !$OMP END PARALLEL DO


      ! This will only retain transport at the boundary
      ! of the region
      fumask=0.
      fvmask=0.
      ftstmask=0.
      !$OMP PARALLEL DO PRIVATE(i,j,ip1,jp1)
      do j=1,jdm
      do i=1,idm
      if (mask(i,j)) then

         jp1=min(jdm,j+1)
         ip1=min(idm,i+1)

         fumask(i,j)=fumask(i,j)+1.
         fvmask(i,j)=fvmask(i,j)+1.

         fumask(ip1,j)=fumask(ip1,j)-1.
         fvmask(i,jp1)=fvmask(i,jp1)-1.

         ftstmask(i,j)=1.
      end if
      end do
      end do
      !$OMP END PARALLEL DO

      ! Decrease mask to active points
      !mask = (abs(fumask)>1e-4 .or. abs(fvmask)>1e-4)
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jdm
      do i=1,idm
         if (depth(i,j)<1.) then
            fumask(i,j)=0.
            fvmask(i,j)=0.
         end if
         mask(i,j) = abs(fumask(i,j))>1e-4 .or. abs(fvmask(i,j))>1e-4
      end do
      end do
      !$OMP END PARALLEL DO


      !if (ilat==nlats/2.and.masktest) then
      !   call tecfld('fumask',fumask,idm,jdm,depth)
      !   call tecfld('fvmask',fvmask,idm,jdm,depth)
      !   call tecfld('fmask',ftstmask,idm,jdm,depth)
      !endif



       uint=0.
       vint=0.
       !if (ideep==nds.and.ilat==nlats/2) then
       !   call tecfld('uint',inttransu,idm,jdm,depth)
       !   call tecfld('vint',inttransv,idm,jdm,depth)
       !endif



       !$OMP PARALLEL DO PRIVATE(i,j,ideep)
       do ideep=1,nds

          !mostrf(ilat,ideep)= sum(fumask*inttransu2(:,:,ideep),mask=mask) + &
          !                    sum(fvmask*inttransv2(:,:,ideep),mask=mask) 


          ! Integrate over latitude band
          do j=1,jdm
          do i=1,idm
          if(mask(i,j)) then
             mostrf(ilat,ideep) = mostrf(ilat,ideep) + &
                   fumask(i,j)*inttransu2(i,j,ideep) + &
                   fvmask(i,j)*inttransv2(i,j,ideep)

            ! Cross section area
            crssec(ilat,ideep) = crssec(ilat,ideep) +  &
               crssecu(i,j,ideep)*abs(fumask(i,j))  +  &
               crssecv(i,j,ideep)*abs(fvmask(i,j))
          end if
          end do
          end do


          mostrf(ilat,ideep)= mostrf(ilat,ideep)*1e-6
          !if (ideep==1.and.ilat>1) then
          !   mostrf(ilat,ideep)= mostrf(ilat,ideep)+  mostrf(ilat-1,ideep)
          !end if

          !if (ds(ideep)>maxdepth) then
          !    mostrf(ilat,ideep)=undef
          !end if
      end do
      !$OMP END PARALLEL DO

      ! Bottom transport residual
      trresidual  = mostrf(ilat,nds)
      !print *,'residual ',mostrf(ilat,nds)

      ! Which translates into the following barotropic velocity (1e-6 m/s)
      velresidual = trresidual / crssec(ilat,nds)

      ! Make transport "nondivergent" -- Subtract transport residual
      do ideep=1,nds
         mostrf(ilat,ideep) = mostrf(ilat,ideep) - &
            crssec(ilat,ideep)*velresidual
      end do

      ! This should be zero now:
      !print *,'residual fixed ',mostrf(ilat,nds)

      ! Set depths below a maximum depth to undefined
      maxdepth=max( maxval(sumdpu(:,:,kdm+1),mask=abs(fumask)>1e-4), &
                    maxval(sumdpv(:,:,kdm+1),mask=abs(fvmask)>1e-4) )
      do ideep=1,nds
          if (ds(ideep)>maxdepth) then
              mostrf(ilat,ideep)=undef
          end if
      end do

   end do

   !call tecfld('mostrf',mostrf,nlats,nds,mostrf)
   !open(10,file='mostrf.tec',status="unknown")
   !write(10,*)'TITLE = "mostrf"'
   !write(10,*)'VARIABLES = "lat" "depth" "mostrf"'
   !write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',nlats,', J=',nds,', K=1'
   !write(10,'(10e15.5)')   ((lats(ilat)        ,ilat=1,nlats),ideep=1,nds)
   !write(10,'(10e15.5)')   ((ds(ideep)         ,ilat=1,nlats),ideep=1,nds)
   !write(10,'(10e15.5)')((mostrf(ilat,ideep),ilat=1,nlats),ideep=1,nds)
   !close(10)


end subroutine mosf


subroutine mosf_sig0(u,v,sal,tem,dp,qlon,qlat,idm,jdm,kdm,mostrf,lats,sigarray,nlats,nsig,undef)
use m_spherdist
implicit none

   integer, intent(in) :: idm,jdm,nlats,nsig,kdm
   real,    intent(in) :: qlat(0:idm+1,0:jdm+1),qlon(0:idm+1,0:jdm+1)
   real,    intent(in), dimension(idm,jdm,kdm) :: u,v,dp,sal,tem
   real,    intent(out) :: mostrf(nlats,nsig)
   real,    intent(in) :: lats(nlats)
   real,    intent(in) :: sigarray(nsig)
   real, intent(in) :: undef
      
   real depthu(idm,jdm),depthv(idm,jdm) ,scvx(idm,jdm),scuy(idm,jdm), &
        depth(idm,jdm),fumask(idm,jdm), fvmask(idm,jdm)
   real,dimension(idm,jdm,kdm) :: dpu, dpv ,transu,transv
   real, dimension(idm,jdm) :: uint,vint,dptmp, inttransu,inttransv,ftstmask
   real, dimension(idm,jdm,kdm+1) :: sumdpu,sumdpv
   real, dimension(idm,jdm,nsig) :: inttransu2,inttransv2, crssecu,crssecv
   logical :: mask(idm,jdm)
   real :: dptmpu, dptmpv, sigdel
   real :: salu,salv,temu,temv,sigu,sigv
   real :: radian, pi, thref
   integer, dimension(idm,jdm,kdm)  :: usig0_index, vsig0_index
   integer, dimension(idm,jdm,nsig) :: ncnt_usigindex, ncnt_vsigindex
   real :: crssec(nlats,nsig)
   real :: velres, trres

   integer :: ideep, ilat, i, j, im1, ip1, jm1, jp1,k, isig,isigu,isigv,nvalid
   logical, parameter :: masktest=.false.

   include 'stmt_funcs.H'

   pi = acos(0.)*2
   radian=180/pi
   thref=1e-3

   sigdel = (maxval(sigarray)-minval(sigarray))/(nsig-1)
   !print *,'sigdel=',sigdel

!
! --- ------------------------------------------------------------------
! --- calculate scale factors and depths in u- and v-points
! --- ------------------------------------------------------------------
!
   !$OMP PARALLEL DO PRIVATE(i,j)
   do j=1,jdm
   do i=1,idm
      depth(i,j) = sum(dp(i,j,:))
   end do
   end do
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO PRIVATE(i,j,ip1,im1,jp1,jm1)
   do j=1,jdm
     jp1=j+1
     jm1=mod(j+jdm-2,jdm)+1
     do i=1,idm
       im1=max(1,i-1)
       ip1=i+1
       scvx(i,j)=spherdist(qlon(ip1,j),qlat(ip1,j), qlon(i  ,j),qlat(i  ,j))
       scuy(i,j)=spherdist(qlon(i,jp1),qlat(i,jp1), qlon(i,j  ),qlat(i,j  ))
       depthu(i,j)=min(depth(i,j),depth(im1,j))
       depthv(i,j)=min(depth(i,j),depth(i,jm1))
       if (depth(i,j).lt.1.) depth(i,j)=0.
       if (depthu(i,j).lt.1.) depthu(i,j)=0.
       if (depthv(i,j).lt.1.) depthv(i,j)=0.
     enddo
   enddo
   !$OMP END PARALLEL DO


!
! --- ------------------------------------------------------------------
! --- calculate dp in u and v points
! --- ------------------------------------------------------------------
!
   sumdpu=0.
   sumdpv=0.
   ncnt_usigindex=0
   ncnt_vsigindex=0
   do k=1,kdm
   !$OMP PARALLEL DO PRIVATE(i,j,im1,jm1,salu,salv,temu,temv,sigu,sigv,isig)
   do j=1,jdm
     jm1=max(j-1,1)
     do i=1,idm
       im1=max(1,i-1)

       dpu(i,j,k)   = 0.5*(dp(i,j,k) + dp(im1,j,k))
       dpv(i,j,k)   = 0.5*(dp(i,j,k) + dp(i,jm1,k))

       sumdpu(i,j,k+1) = min(depthu(i,j),sumdpu(i,j,k)  + dpu(i,j,k))
       sumdpv(i,j,k+1) = min(depthv(i,j),sumdpv(i,j,k)  + dpv(i,j,k))

       dpu(i,j,k) = sumdpu(i,j,k+1)  - sumdpu(i,j,k)
       dpv(i,j,k) = sumdpv(i,j,k+1)  - sumdpv(i,j,k)

       transu(i,j,k)=dpu(i,j,k)*u(i,j,k)*scuy(i,j)
       transv(i,j,k)=dpv(i,j,k)*v(i,j,k)*scvx(i,j)

       !temu = (dp(i,j,k)*tem(i,j,k)+tem(im1,j,k)*dp(im1,j,k)) / &
       !       (2*max(dpu(i,j,k),1e-4))
       !temv = (dp(i,j,k)*tem(i,j,k)+tem(i,jm1,k)*dp(i,jm1,k)) / &
       !       (2*max(dpv(i,j,k),1e-4))
       !salu = (dp(i,j,k)*sal(i,j,k)+sal(im1,j,k)*dp(im1,j,k)) / &
       !       (2*max(dpu(i,j,k),1e-4))
       !salv = (dp(i,j,k)*sal(i,j,k)+sal(i,jm1,k)*dp(i,jm1,k)) / &
       !       (2*max(dpv(i,j,k),1e-4))

       temu = .5*(tem(i,j,k)+tem(im1,j,k))
       temv = .5*(tem(i,j,k)+tem(i,jm1,k))
       salu = .5*(sal(i,j,k)+sal(im1,j,k))
       salv = .5*(sal(i,j,k)+sal(i,jm1,k))
   
       ! Densities
       sigu = sig(temu,salu)
       sigv = sig(temv,salv)

       ! Put into index array
       usig0_index(i,j,k)= max(1,min(floor((sigu-minval(sigarray))/sigdel)+1,nsig))
       vsig0_index(i,j,k)= max(1,min(floor((sigv-minval(sigarray))/sigdel)+1,nsig))
       !print *,i,j,k,temu,salu,sigu,usig0_index(i,j,k),sal(i,j,k)

       if ( dpu(i,j,k) <1e-4 ) then
          isig=usig0_index(i,j,k)
          ncnt_usigindex(i,j,isig) =  ncnt_usigindex(i,j,isig) +1
       end if
       if ( dpv(i,j,k) <1e-4 ) then
          isig=vsig0_index(i,j,k)
          ncnt_vsigindex(i,j,isig) =  ncnt_vsigindex(i,j,isig) +1
       end if

     end do
   end do
   !$OMP END PARALLEL DO
   end do
   !print *,maxval(transu),minval(transu)
   !print *,maxval(transv),minval(transv)
   !print *,maxval(usig0_index),minval(usig0_index)
   !print *,maxval(vsig0_index),minval(vsig0_index)
   !call tecfld('transu',transu(:,:,5),idm,jdm,depth)
   !call tecfld('usig0_index',float(usig0_index(:,:,5)),idm,jdm,depth)


   ! Top-down
   inttransu2=0.
   inttransv2=0.
   crssecu=0.
   crssecv=0.
   do k=1,kdm
   !$OMP PARALLEL DO PRIVATE(i,j,isigu,isigv)
   do j=1,jdm
   do i=1,idm

      isigu=usig0_index(i,j,k)
      isigv=vsig0_index(i,j,k)
      
      inttransu2(i,j,isigu) = inttransu2(i,j,isigu)+ transu(i,j,k)
      inttransv2(i,j,isigv) = inttransv2(i,j,isigv)+ transv(i,j,k)

      ! Cross section area for layers
      crssecu(i,j,isigu) = crssecu(i,j,isigu) + dpu(i,j,k)*scuy(i,j)
      crssecv(i,j,isigu) = crssecv(i,j,isigv) + dpv(i,j,k)*scvx(i,j)

   end do
   end do
   !$OMP END PARALLEL DO
   end do
   !print *,'inttrans'
   !call tecfld('inttransu_1',inttransu2(:,:,5),idm,jdm,depth)
   !print *,maxval(inttransu2),minval(inttransu2)
   !print *,maxval(inttransv2),minval(inttransv2)

   ! The above gives the transport in each sigma layer. Integrate.
   do isig=2,nsig
   !$OMP PARALLEL DO PRIVATE(i,j)
   do j=1,jdm
   do i=1,idm
      inttransu2(i,j,isig) = inttransu2(i,j,isig)+ inttransu2(i,j,isig-1)
      inttransv2(i,j,isig) = inttransv2(i,j,isig)+ inttransv2(i,j,isig-1)

      crssecu(i,j,isig) = crssecu(i,j,isig) + crssecu(i,j,isig-1)
      crssecv(i,j,isig) = crssecv(i,j,isig) + crssecv(i,j,isig-1)
   end do
   end do
   !$OMP END PARALLEL DO
   end do
   !print *,'inttrans'
   !call tecfld('inttransu_2',inttransu2(:,:,5),idm,jdm,depth)
   !print *,maxval(inttransu2),minval(inttransu2)
   !print *,maxval(inttransv2),minval(inttransv2)


!
! --- ------------------------------------------------------------------
! --- calculate  meridional overturning stream function
! --- ------------------------------------------------------------------
!



   ! For every zonal band:
   mostrf=0.
   crssec=0.
   do ilat=1,nlats
      !print *,ilat,nlats,lats(ilat)

      ! Find points inside this region (North of lats(i)
      mask=.false.
      !where (plat(1:idm,1:jdm)>lats(ilat).and.depth>0.)
      !where (qlat(1:idm,1:jdm)>lats(ilat))
      !   mask=.true.
      !endwhere
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jdm
      do i=1,idm
         mask(i,j)=qlat(i,j)>lats(ilat)
      end do
      end do
      !$OMP END PARALLEL DO


      ! This will only retain transport at the boundary
      ! of the region
      fumask=0.
      fvmask=0.
      ftstmask=0.
      !$OMP PARALLEL DO PRIVATE(i,j,ip1,jp1)
      do j=1,jdm
      do i=1,idm
      if (mask(i,j)) then

         jp1=min(jdm,j+1)
         ip1=min(idm,i+1)

         fumask(i,j)=fumask(i,j)+1.
         fvmask(i,j)=fvmask(i,j)+1.

         fumask(ip1,j)=fumask(ip1,j)-1.
         fvmask(i,jp1)=fvmask(i,jp1)-1.

         ftstmask(i,j)=1.
      end if
      end do
      end do
      !$OMP END PARALLEL DO

      ! Decrease mask to active points
      !mask = (abs(fumask)>1e-4 .or. abs(fvmask)>1e-4) !.and. depth(i,j)>1.
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jdm
      do i=1,idm
         if (depth(i,j)<1.) then
            fumask(i,j)=0.
            fvmask(i,j)=0.
         end if
         mask(i,j) = abs(fumask(i,j))>1e-4 .or. abs(fvmask(i,j))>1e-4
      end do
      end do
      !$OMP END PARALLEL DO



      !if (ilat==nlats/2.and.masktest) then
      !   call tecfld('fumask',fumask,idm,jdm,depth)
      !   call tecfld('fvmask',fvmask,idm,jdm,depth)
      !   call tecfld('fmask',ftstmask,idm,jdm,depth)
      !endif



       !uint=0.
       !vint=0.
       !if (ideep==10) then
       !if (ilat==nlats/2) then
       !   call tecfld('uint',inttransu2(:,:,5),idm,jdm,depth)
       !   call tecfld('vint',inttransv2(:,:,5),idm,jdm,depth)
       !endif

       !$OMP PARALLEL DO PRIVATE(i,j,isig)
       do isig=1,nsig
          ! Integrate over latitude band
          !mostrf(ilat,isig)= sum(fumask*inttransu2(:,:,isig),mask=mask) + &
          !                   sum(fvmask*inttransv2(:,:,isig),mask=mask) 
          do j=1,jdm
          do i=1,idm
          if(mask(i,j)) then
             mostrf(ilat,isig)= mostrf(ilat,isig)+ &
             fumask(i,j)*inttransu2(i,j,isig)    + &
             fvmask(i,j)*inttransv2(i,j,isig)   

             crssec(ilat,isig) = crssec(ilat,isig) + &
                crssecu(i,j,isig)*abs(fumask(i,j)) + &
                crssecv(i,j,isig)*abs(fvmask(i,j))
          end if
          end do
          end do
             
             
          mostrf(ilat,isig)= mostrf(ilat,isig)*1e-6

          !nvalid= sum(abs(nint(fumask))*ncnt_usigindex(:,:,isig)) + &
          !        sum(abs(nint(fvmask))*ncnt_vsigindex(:,:,isig))
          !if (nvalid==0) then
          !   mostrf(ilat,isig)= undef
          !end if
             
      end do
      !$OMP END PARALLEL DO 

      ! Find transport residual
      trres=mostrf(ilat,nsig)

      ! velocity residual
      velres=trres/crssec(ilat,nsig)

      ! Make stream function "nondivergent"
      do isig=1,nsig
         mostrf(ilat,isig)= mostrf(ilat,isig) - velres*crssec(ilat,isig)
      end do

   end do

   !print *,minval(mostrf,mostrf/=undef),maxval(mostrf,mostrf/=undef)

   !call tecfld('mostrf',mostrf,nlats,nds,mostrf)
   !open(10,file='mostrfsig0.tec',status="unknown")
   !write(10,*)'TITLE = "mostrf"'
   !write(10,*)'VARIABLES = "lat" "depth" "mostrf"'
   !write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',nlats,', J=',nsig,', K=1'
   !write(10,'(10e15.5)')   ((lats(ilat)        ,ilat=1,nlats),isig=1,nsig)
   !write(10,'(10e15.5)')   ((sigarray(isig)         ,ilat=1,nlats),isig=1,nsig)
   !write(10,'(10e15.5)')((mostrf(ilat,isig),ilat=1,nlats),isig=1,nsig)
   !close(10)


end subroutine mosf_sig0



subroutine mosf_theta0(u,v,tem,dp,qlon,qlat,idm,jdm,kdm,mostrf,lats,theta0array,nlats,ntheta0,undef)
use m_spherdist
implicit none

   integer, intent(in) :: idm,jdm,kdm,nlats,ntheta0
   real,    intent(in) :: qlat(0:idm+1,0:jdm+1),qlon(0:idm+1,0:jdm+1)
   real,    intent(in), dimension(idm,jdm,kdm) :: u,v,dp,tem
   real,    intent(out) :: mostrf(nlats,ntheta0)
   real,    intent(in) :: lats(nlats)
   real,    intent(in) :: theta0array(ntheta0)
   real, intent(in) :: undef
      
   real,    dimension(idm,jdm) ::  depthu,depthv,scvx,scuy, &
      depth,fumask, fvmask, ftstmask
   real,    dimension(idm,jdm,kdm)   :: dpu, dpv ,transu,transv
   real,    dimension(idm,jdm,kdm+1) :: sumdpu,sumdpv
   real,    dimension(idm,jdm,ntheta0) :: inttransu2,inttransv2, crssecu,crssecv
   integer, dimension(idm,jdm,ntheta0) :: ncnt_uthetaindex, ncnt_vthetaindex
   integer, dimension(idm,jdm,kdm) :: utheta0_index, vtheta0_index
   logical :: mask(idm,jdm)
   real :: thetadel
   real :: temu,temv
   real :: crssec(nlats,ntheta0)
   real :: trres, velres

   integer :: ideep, ilat, i, j, im1, ip1, jm1, jp1,k, itheta,ithetau,ithetav,nvalid
   logical, parameter :: masktest=.false.

   include 'stmt_funcs.H'
   !pi = acos(0.)*2
   !radian=180/pi
   !thref=1e-3

   print *,'mosf_theta entry'
   thetadel = (maxval(theta0array)-minval(theta0array))/(ntheta0-1)
   !print *,'thetadel=',thetadel

!
! --- ------------------------------------------------------------------
! --- calculate scale factors and depths in u- and v-points
! --- ------------------------------------------------------------------
!
   print *,'mosf_theta depth, scuy,scvx'
   !$OMP PARALLEL DO PRIVATE(i,j)
   do j=1,jdm
   do i=1,idm
      depth(i,j) = sum(dp(i,j,:))
   end do
   end do
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO PRIVATE(i,j,ip1,im1,jp1,jm1)
   do j=1,jdm
     jp1=j+1
     jm1=mod(j+jdm-2,jdm)+1
     do i=1,idm
       im1=max(1,i-1)
       ip1=i+1
       scvx(i,j)=spherdist(qlon(ip1,j),qlat(ip1,j), qlon(i  ,j),qlat(i  ,j))
       scuy(i,j)=spherdist(qlon(i,jp1),qlat(i,jp1), qlon(i,j  ),qlat(i,j  ))
       depthu(i,j)=min(depth(i,j),depth(im1,j))
       depthv(i,j)=min(depth(i,j),depth(i,jm1))
       if (depth(i,j).lt.1.) depth(i,j)=0.
       if (depthu(i,j).lt.1.) depthu(i,j)=0.
       if (depthv(i,j).lt.1.) depthv(i,j)=0.
     enddo
   enddo
   !$OMP END PARALLEL DO


!
! --- ------------------------------------------------------------------
! --- calculate dp in u and v points
! --- ------------------------------------------------------------------
!
   sumdpu=0.
   sumdpv=0.
   ncnt_uthetaindex=0
   ncnt_vthetaindex=0
   print *,'dpu,dpv'
   do k=1,kdm
   !$OMP PARALLEL DO PRIVATE(i,j,im1,jm1temu,temv)
   do j=1,jdm
     jm1=max(j-1,1)
     do i=1,idm
       im1=max(1,i-1)

       dpu(i,j,k)   = 0.5*(dp(i,j,k) + dp(im1,j,k))
       dpv(i,j,k)   = 0.5*(dp(i,j,k) + dp(i,jm1,k))

       sumdpu(i,j,k+1) = min(depthu(i,j),sumdpu(i,j,k)  + dpu(i,j,k))
       sumdpv(i,j,k+1) = min(depthv(i,j),sumdpv(i,j,k)  + dpv(i,j,k))

       dpu(i,j,k) = sumdpu(i,j,k+1)  - sumdpu(i,j,k)
       dpv(i,j,k) = sumdpv(i,j,k+1)  - sumdpv(i,j,k)

       transu(i,j,k)=dpu(i,j,k)*u(i,j,k)*scuy(i,j)
       transv(i,j,k)=dpv(i,j,k)*v(i,j,k)*scvx(i,j)

       !temu = (dp(i,j,k)*tem(i,j,k)+tem(im1,j,k)*dp(im1,j,k)) / &
       !       (2*max(dpu(i,j,k),1e-4))
       !temv = (dp(i,j,k)*tem(i,j,k)+tem(i,jm1,k)*dp(i,jm1,k)) / &
       !       (2*max(dpv(i,j,k),1e-4))
       temu = .5*(tem(i,j,k)+tem(im1,j,k))
       temv = .5*(tem(i,j,k)+tem(i,jm1,k))

       ! Put into index array
       utheta0_index(i,j,k)= min(max(1,floor((temu-minval(theta0array))/thetadel)+1),ntheta0)
       vtheta0_index(i,j,k)= min(max(1,floor((temv-minval(theta0array))/thetadel)+1),ntheta0)
       !print *,i,j,k,temu,utheta0_index(i,j,k),sal(i,j,k)

       if ( dpu(i,j,k) <1e-4 ) then
          itheta=utheta0_index(i,j,k)
          ncnt_uthetaindex(i,j,itheta) =  ncnt_uthetaindex(i,j,itheta) +1
       end if
       if ( dpv(i,j,k) <1e-4 ) then
          itheta=vtheta0_index(i,j,k)
          ncnt_vthetaindex(i,j,itheta) =  ncnt_vthetaindex(i,j,itheta) +1
       end if

     end do
   end do
   !$OMP END PARALLEL DO
   end do


   ! Top-down
   print *,' layer transport'
   inttransu2=0.
   crssecu=0.
   crssecv=0.
   inttransv2=0.
   do k=1,kdm
   !$OMP PARALLEL DO PRIVATE(i,j,ithetau,ithetav)
   do j=1,jdm
   do i=1,idm

      ithetau=utheta0_index(i,j,k)
      ithetav=vtheta0_index(i,j,k)
      
      inttransu2(i,j,ithetau) = inttransu2(i,j,ithetau)+ transu(i,j,k)
      inttransv2(i,j,ithetav) = inttransv2(i,j,ithetav)+ transv(i,j,k)

      ! Cross section area for theta layers
      crssecu(i,j,itheta) = crssecu(i,j,itheta) + scuy(i,j)*dpu(i,j,k)
      crssecv(i,j,itheta) = crssecv(i,j,itheta) + scvx(i,j)*dpv(i,j,k)

   end do
   end do
   !$OMP END PARALLEL DO
   end do

   !! The above gives the transport in each layer. Integrate.
   print *,' integrated transport'
   do itheta=ntheta0-1,1,-1
   !$OMP PARALLEL DO PRIVATE(i,j)
   do j=1,jdm
   do i=1,idm
      inttransu2(i,j,itheta) = inttransu2(i,j,itheta)+ inttransu2(i,j,itheta+1)
      inttransv2(i,j,itheta) = inttransv2(i,j,itheta)+ inttransv2(i,j,itheta+1)

      crssecu(i,j,itheta) = crssecu(i,j,itheta) + crssecu(i,j,itheta+1) 
      crssecv(i,j,itheta) = crssecv(i,j,itheta) + crssecv(i,j,itheta+1) 
   end do
   end do
   !$OMP END PARALLEL DO
   end do



!
! --- ------------------------------------------------------------------
! --- calculate  meridional overturning stream function
! --- ------------------------------------------------------------------
!



   ! For every zonal band:
   mostrf=0.
   crssec=0.
   !oldmask=.false.
   !print *,maxval(dpv(:,:,1)),maxval(dpu(:,:,1)),maxval(dp(:,:,1))
   mask=.false.
   do ilat=1,nlats
      !print *,ilat,nlats,lats(ilat)

      ! Find points inside this region (North of lats(i)
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jdm
      do i=1,idm
         mask(i,j)=qlat(i,j)>lats(ilat)
      end do
      end do
      !$OMP END PARALLEL DO


      ! This will only retain transport at the boundary
      ! of the region
      fumask=0.
      fvmask=0.
      ftstmask=0.
      !$OMP PARALLEL DO PRIVATE(i,j,ip1,jp1)
      do j=1,jdm
      do i=1,idm
      if (mask(i,j)) then

         jp1=min(jdm,j+1)
         ip1=min(idm,i+1)

         fumask(i,j)=fumask(i,j)+1.
         fvmask(i,j)=fvmask(i,j)+1.

         fumask(ip1,j)=fumask(ip1,j)-1.
         fvmask(i,jp1)=fvmask(i,jp1)-1.

         ftstmask(i,j)=1.
      end if
      end do
      end do
      !$OMP END PARALLEL DO

      ! Decrease mask to active points
      !mask = (abs(fumask)>1e-4 .or. abs(fvmask)>1e-4)
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jdm
      do i=1,idm
         if (depth(i,j)<1.) then
            fumask(i,j)=0.
            fvmask(i,j)=0.
         end if
         mask(i,j) = abs(fumask(i,j))>1e-4 .or. abs(fvmask(i,j))>1e-4
      end do
      end do
      !$OMP END PARALLEL DO


       !$OMP PARALLEL DO PRIVATE (i,j,itheta)
       do itheta=1,ntheta0

          do j=1,jdm
          do i=1,idm
          if (mask(i,j)) then
             mostrf(ilat,itheta) =  mostrf(ilat,itheta) + &
                fumask(i,j)*inttransu2(i,j,itheta)      + &
                fvmask(i,j)*inttransv2(i,j,itheta)     

             crssec(ilat,itheta) = crssec(ilat,itheta)  + &
                abs(fumask(i,j))*crssecu(i,j,itheta)    + &
                abs(fvmask(i,j))*crssecv(i,j,itheta)
          end if
          end do
          end do

          mostrf(ilat,itheta)= mostrf(ilat,itheta)*1e-6

          !nvalid= sum(abs(nint(fumask))*ncnt_uthetaindex(:,:,itheta)) + &
          !        sum(abs(nint(fvmask))*ncnt_vthetaindex(:,:,itheta))
          !if (nvalid==0) then
          !   mostrf(ilat,itheta)= undef
          !end if
          
      end do
      !$OMP END PARALLEL DO

      print *,ilat,mostrf(ilat,1),crssec(ilat,1)
      ! transport residual at "bottom"
      trres=mostrf(ilat,1)

      ! Velocity residual (1e-6 m/s)
      velres=trres/crssec(ilat,1)

      print *,ilat,trres,velres,crssec(ilat,1)

      ! Make streamfunction "nondivergent"
      do itheta=1,ntheta0
        mostrf(ilat,itheta) = mostrf(ilat,itheta) - velres*crssec(ilat,itheta)
     end do


   end do

   !print *,minval(mostrf),maxval(mostrf)

   !call tecfld('mostrf',mostrf,nlats,nds,mostrf)
   !open(10,file='mostrftheta0.tec',status="unknown")
   !write(10,*)'TITLE = "mostrf"'
   !write(10,*)'VARIABLES = "lat" "depth" "mostrf"'
   !write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',nlats,', J=',ntheta0,', K=1'
   !write(10,'(10e15.5)')   ((lats(ilat)        ,ilat=1,nlats),itheta=1,ntheta0)
   !write(10,'(10e15.5)')   ((theta0array(itheta)         ,ilat=1,nlats),itheta=1,ntheta0)
   !write(10,'(10e15.5)')((mostrf(ilat,itheta),ilat=1,nlats),itheta=1,ntheta0)
   !close(10)


end subroutine mosf_theta0


subroutine mht(u,v,sal,tem,dp,qlon,qlat,idm,jdm,kdm,merht,lats,nlats)
use m_spherdist
implicit none

   integer, intent(in) :: idm,jdm,kdm,nlats
   real,    intent(in) :: qlat(0:idm+1,0:jdm+1),qlon(0:idm+1,0:jdm+1)
   real,    intent(in), dimension(idm,jdm,kdm) :: u,v,dp,tem,sal
   real,    intent(out):: merht(nlats)
   real,    intent(in) :: lats(nlats)
      
   real depthu(idm,jdm),depthv(idm,jdm) ,scvx(idm,jdm),scuy(idm,jdm), &
        depth(idm,jdm),fumask(idm,jdm), fvmask(idm,jdm)
   real,dimension(idm,jdm,kdm) :: dpu, dpv, temu,temv,sigu,sigv
   real, dimension(idm,jdm) :: uint,vint,dptmp
   real, dimension(idm,jdm,kdm+1) :: sumdpu,sumdpv
   logical :: mask(idm,jdm)
   real :: salu,salv

   integer :: ideep, ilat, i, j, im1, ip1, jm1, jp1,k 
   logical, parameter :: masktest=.false.
   real, parameter :: cpsw=3987.
   real :: bartr,barvel,crssec

   real :: pi, radian, thref ! needed in stmt_funcs

   include 'stmt_funcs.H'

   pi = acos(0.)*2
   radian=180./pi
   thref=1e-3

!
! --- ------------------------------------------------------------------
! --- calculate scale factors and depths in u- and v-points
! --- ------------------------------------------------------------------
!
   !$OMP PARALLEL DO PRIVATE(i,j)
   do j=1,jdm
   do i=1,idm
      depth(i,j) = sum(dp(i,j,:))
   end do
   end do
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO PRIVATE(i,j,im1,ip1,jm1,jp1)
   do j=1,jdm
     jp1=j+1
     jm1=mod(j+jdm-2,jdm)+1
     do i=1,idm
       im1=max(1,i-1)
       ip1=i+1
       scvx(i,j)=spherdist(qlon(ip1,j),qlat(ip1,j), qlon(i  ,j),qlat(i  ,j))
       scuy(i,j)=spherdist(qlon(i,jp1),qlat(i,jp1), qlon(i,j  ),qlat(i,j  ))
       depthu(i,j)=min(depth(i,j),depth(im1,j))
       depthv(i,j)=min(depth(i,j),depth(i,jm1))
       if (depth(i,j).lt.1.) depth(i,j)=0.
       if (depthu(i,j).lt.1.) depthu(i,j)=0.
       if (depthv(i,j).lt.1.) depthv(i,j)=0.
     enddo
   enddo
   !$OMP END PARALLEL DO


!
! --- ------------------------------------------------------------------
! --- calculate dp in u and v points
! --- ------------------------------------------------------------------
!
   sumdpu=0.
   sumdpv=0.
   do k=1,kdm
   !$OMP PARALLEL DO PRIVATE(i,j,im1,jm1)
   do j=1,jdm
     jm1=max(1,j-1)
     do i=1,idm
       im1=max(1,i-1)

       dpu(i,j,k)   = 0.5*(dp(i,j,k) + dp(im1,j,k))
       dpv(i,j,k)   = 0.5*(dp(i,j,k) + dp(i,jm1,k))

       sumdpu(i,j,k+1) = min(depthu(i,j),sumdpu(i,j,k)  + dpu(i,j,k))
       sumdpv(i,j,k+1) = min(depthv(i,j),sumdpv(i,j,k)  + dpv(i,j,k))

       dpu(i,j,k) = sumdpu(i,j,k+1)  - sumdpu(i,j,k)
       dpv(i,j,k) = sumdpv(i,j,k+1)  - sumdpv(i,j,k)
       
     end do
   end do
   !$OMP END PARALLEL DO
   end do
      
   do k=1,kdm
      !$OMP PARALLEL DO PRIVATE(i,j,im1,jm1,salu,salv)
      do j=1,jdm
        jm1=max(1,j-1)
        do i=1,idm
          im1=max(1,i-1)
          !temu(i,j,k) = (dp(i,j,k)*tem(i,j,k) + dp(im1,j,k)*tem(im1,j,k))/ &
          !              (2*max(1e-4,dpu(i,j,k)))
          !temv(i,j,k) = (dp(i,j,k)*tem(i,j,k) + dp(i,jm1,k)*tem(i,jm1,k))/ &
          !              (2*max(1e-4,dpv(i,j,k)))
          !salu        = (dp(i,j,k)*sal(i,j,k) + dp(im1,j,k)*sal(im1,j,k))/ &
          !              (2*max(1e-4,dpu(i,j,k)))
          !salv        = (dp(i,j,k)*sal(i,j,k) + dp(i,jm1,k)*sal(i,jm1,k))/ &
          !              (2*max(1e-4,dpv(i,j,k)))
          temu(i,j,k) = .5*(tem(i,j,k) + tem(im1,j,k))
          temv(i,j,k) = .5*(tem(i,j,k) + tem(i,jm1,k))
          salu        = .5*(sal(i,j,k) + sal(im1,j,k))
          salv        = .5*(sal(i,j,k) + sal(i,jm1,k))

          sigu(i,j,k) = sig(temu(i,j,k),salu)+1000.
          sigv(i,j,k) = sig(temv(i,j,k),salv)+1000.
        end do
      end do
      !$OMP END PARALLEL DO
   end do

!
! --- ------------------------------------------------------------------
! --- calculate  meridional overturning stream function
! --- ------------------------------------------------------------------
!



   ! For every zonal band:
   merht=0.
   do ilat=1,nlats

      ! Find points inside this region (North of lats(i)
      !mask=.false.
      !where (qlat(1:idm,1:jdm)>lats(ilat))
      !   mask=.true.
      !endwhere
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jdm
      do i=1,idm
         mask(i,j)=qlat(i,j)>lats(ilat)
      end do
      end do
      !$OMP END PARALLEL DO


      ! This will only retain transport at the boundary
      ! of the region
      fumask=0.
      fvmask=0.
      !$OMP PARALLEL DO PRIVATE(i,j,ip1,jp1)
      do j=1,jdm
      do i=1,idm
      if (mask(i,j)) then

         jp1=min(jdm,j+1)
         ip1=min(idm,i+1)

         fumask(i,j)=fumask(i,j)+1.
         fvmask(i,j)=fvmask(i,j)+1.

         fumask(ip1,j)=fumask(ip1,j)-1.
         fvmask(i,jp1)=fumask(i,jp1)-1.
      end if
      end do
      end do
      !$OMP END PARALLEL DO

      ! Drop points over land
      !where (depth<10.0)
      !   fumask=0.
      !   fvmask=0.
      !endwhere

      !$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jdm
      do i=1,idm
         if (depth(i,j)<1.) then
            fumask(i,j)=0.
            fvmask(i,j)=0.
         end if
         mask(i,j) = abs(fumask(i,j))>1e-4 .or. abs(fvmask(i,j))>1e-4
      end do
      end do
      !$OMP END PARALLEL DO

      !if (ilat==nlats/2.and.masktest) then
      !   call tecfld('fumask',fumask,idm,jdm,depth)
      !   call tecfld('fvmask',fvmask,idm,jdm,depth)
      !endif

      ! We should now have transport masks
      ! for the transport into the region
      ! north of lats(i). Let the streamfunction
      ! be zero at the surface

      ! Integrate over latitude band
      crssec = 0.
      bartr=0.
      barvel=0.

      ! Cross section area and barotropic velocity for transport
      do k=1, kdm
         !$OMP PARALLEL DO PRIVATE(i,j)  &
         !$OMP REDUCTION(+:crssec,bartr)
         do j=1,jdm
         do i=1,idm
         if (mask(i,j)) then

            ! cross section
            crssec= crssec + &
            abs(fumask(i,j))*dpu(i,j,k)*scuy(i,j) + &
            abs(fvmask(i,j))*dpv(i,j,k)*scvx(i,j)

            ! Section transport
            bartr= bartr + &
            fumask(i,j)*u(i,j,k)*dpu(i,j,k)*scuy(i,j) + &
            fvmask(i,j)*v(i,j,k)*dpv(i,j,k)*scvx(i,j)

         end if
         end do
         end do
         !$OMP END PARALLEL DO
      end do

      ! Section velocity component
      barvel=bartr/crssec
      !print *,ilat,barvel,bartr,crssec


      do k=1, kdm
         !$OMP PARALLEL DO PRIVATE(i,j)
         do j=1,jdm
         do i=1,idm
         if (mask(i,j)) then

            ! heat transport
            merht(ilat)= merht(ilat) + &
            fumask(i,j)*(u(i,j,k)-barvel)*temu(i,j,k)*dpu(i,j,k)*sigu(i,j,k)*scuy(i,j) + &
            fvmask(i,j)*(v(i,j,k)-barvel)*temv(i,j,k)*dpv(i,j,k)*sigv(i,j,k)*scvx(i,j)

         end if
         end do
         end do
         !$OMP END PARALLEL DO
      end do

      merht(ilat)=merht(ilat)*cpsw
   end do

end subroutine mht


end module mod_overturning






