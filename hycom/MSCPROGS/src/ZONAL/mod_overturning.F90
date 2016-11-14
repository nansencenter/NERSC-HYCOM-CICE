module mod_overturning

contains


subroutine mosf(u,v,dp,idm,jdm,kdm,mostrf,lats,nlats,deep,ndeep,undef)
use mod_grid
implicit none

   integer, intent(in) :: idm,jdm,nlats,ndeep,kdm
   real,    intent(in), dimension(idm,jdm,kdm) :: u,v,dp
   real,    intent(out) :: mostrf(nlats,ndeep)
   real,    intent(in) :: lats(nlats)
   real,    intent(in) :: deep(ndeep)
   real,    intent(in) :: undef
      
   real,dimension(idm,jdm) :: dpu, dpv, fumask, fvmask, dp2,  &
      ftstmask, sumdpu, sumdpv, inttransu,inttransv,crssecu,crssecv
   logical :: mask(idm,jdm)
   real :: dptmpu, dptmpv
   real :: crssec(nlats,ndeep),strmf(nlats,ndeep)

   integer :: ideep, ilat, i, j, im1, ip1, jm1, jp1,k 
   logical, parameter :: masktest=.false.
   real :: maxdepth, velresidual, trresidual, dz, dlat

   real, external :: spherdist


   !print *,lats
   !print *,deep
   print *,'max dp:',maxval(dp)

 
   ! --- ------------------------------------------------------------------
   ! --- Cycle z-layers
   ! --- ------------------------------------------------------------------
   do ideep=1,ndeep
      inttransu=0.
      inttransv=0.
      crssecu=0.
      crssecv=0.



      ! --- ------------------------------------------------------------------
      ! --- Sum up for this depth interval ( deep(ideeep-1) -> deep(ideep) )
      ! --- ------------------------------------------------------------------
      sumdpu=0.
      sumdpv=0.
      dp2=0.
      do k=1,kdm
         do j=1,jdm
         do i=1,idm
            jm1=max(1,j-1)
            im1=mod(idm+i-2,idm)+1

            ! Layer thickness u/v points
            !dpu(i,j)   = dp(i,j,k) !0.5*(dp(i,j,k) + dp(im1,j,k))
            !dpv(i,j)   = dp(i,j,k) !0.5*(dp(i,j,k) + dp(i,jm1,k))
            dpu(i,j)   = 0.5*(dp(i,j,k) + dp(im1,j,k))
            dpv(i,j)   = 0.5*(dp(i,j,k) + dp(i,jm1,k))

            ! Fraction of layer between ideep-1 and ideep
            if (ideep==1) then
               dptmpu= max( min(sumdpu(i,j)+dpu(i,j),deep(ideep)) - sumdpu(i,j), 0.)
               dptmpv= max( min(sumdpv(i,j)+dpv(i,j),deep(ideep)) - sumdpv(i,j), 0.)
            else
               dptmpu= max( min(sumdpu(i,j)+dpu(i,j),deep(ideep)) - max(sumdpu(i,j),deep(ideep-1)), 0.)
               dptmpv= max( min(sumdpv(i,j)+dpv(i,j),deep(ideep)) - max(sumdpv(i,j),deep(ideep-1)), 0.)
            end if

            ! Test array - should be equal to deep(ideep) - deep(ideep-1)
            dp2(i,j)=dp2(i,j)+dptmpu

            if (ideep==1) then
               if (dp2(i,j)>deep(ideep)) then
                  print '(a,2i5,3f14.2)','Error: ',i,j,dp2(i,j), deep(ideep),deep(ideep-1)
                  print *,'Error: ',sumdpu(i,j),sumdpu(i,j)+dpu(i,j)
                  stop
               end if
            else
               if (dp2(i,j)>deep(ideep)-deep(ideep-1)) then
                  print '(a,2i5,3f14.2)','Error: ',i,j,dp2(i,j), deep(ideep),deep(ideep-1)
                  print *,'Error: ',sumdpu(i,j),sumdpu(i,j)+dpu(i,j)
                  stop
               end if
            end if

            ! u/v transport in z-layers (m^3 s^-1)
            inttransu(i,j) = inttransu(i,j)+dptmpu*u(i,j,k)*scuy(i,j)
            inttransv(i,j) = inttransv(i,j)+dptmpv*v(i,j,k)*scvx(i,j)

            ! For computing "cross section area" later
            crssecu(i,j)= crssecu(i,j)+dptmpu*scuy(i,j)
            crssecv(i,j)= crssecv(i,j)+dptmpv*scvx(i,j)

            ! Augment sum for next k iteration
            sumdpu(i,j)=sumdpu(i,j)+dpu(i,j)
            sumdpv(i,j)=sumdpv(i,j)+dpv(i,j)

         end do
         end do
      end do
      !print *,maxval(dp2),deep(ideep),deep(max(1,ideep-1)),maxval(sumdpu)

      ! --- ------------------------------------------------------------------
      ! --- Set up masks 
      ! --- ------------------------------------------------------------------
      mostrf(:,ideep)=0.
      crssec=0.
      do ilat=1,nlats

         ! Find points inside this region (North of lats(i)
         !$OMP PARALLEL DO PRIVATE(i,j)
         do j=1,jdm
         do i=1,idm
            mask(i,j)=qlat(i,j)>lats(ilat)
         end do
         end do
         !$OMP END PARALLEL DO


         ! This will only retain transport at the boundary of the region
         fumask=0.
         fvmask=0.
         ftstmask=0.
         !$OMP PARALLEL DO PRIVATE(i,j,ip1,jp1)
         do j=1,jdm
         do i=1,idm
         if (mask(i,j)) then
            jp1=min(jdm,j+1)
            if (periodic) then
               ip1=mod(i,idm)+1
            else
               ip1=min(idm,i+1)
            end if
            fumask(i,j)=fumask(i,j)+1.
            fvmask(i,j)=fvmask(i,j)+1.
            fumask(ip1,j)=fumask(ip1,j)-1.
            fvmask(i,jp1)=fvmask(i,jp1)-1.
         end if
         end do
         end do
         !$OMP END PARALLEL DO

         ! Decrease mask to active points
         !$OMP PARALLEL DO PRIVATE(i,j)
         do j=1,jdm
         do i=1,idm
            mask(i,j) = abs(fumask(i,j))>1e-4 .or. abs(fvmask(i,j))>1e-4
         end do
         end do
         !$OMP END PARALLEL DO


          ! --- ------------------------------------------------------------------
          ! --- transport for one depth (ideep) and one latitude band
          ! --- ------------------------------------------------------------------
          ! Integrate over latitude band
          do j=1,jdm
          do i=1,idm
          if(mask(i,j)) then
             mostrf(ilat,ideep) = mostrf(ilat,ideep) + &
                   fumask(i,j)*inttransu(i,j) + &
                   fvmask(i,j)*inttransv(i,j)

            ! Cross section area
            crssec(ilat,ideep) = crssec(ilat,ideep) +  &
               crssecu(i,j)*abs(fumask(i,j))  +  &
               crssecv(i,j)*abs(fvmask(i,j))
          end if
          end do
          end do

          mostrf(ilat,ideep)= mostrf(ilat,ideep)*1e-6

         ! Set depths below a maximum depth to undefined
         maxdepth=max( maxval(sumdpu,mask=abs(fumask)>1e-4), &
                       maxval(sumdpv,mask=abs(fvmask)>1e-4) )
          if (deep(ideep)>maxdepth) mostrf(ilat,ideep)=undef
      end do ! latitude bands ilat
   end do ! depth iterval ideep


   ! --- ------------------------------------------------------------------
   ! --- Convert from transport to meridional overturning streamfunction
   ! --- ------------------------------------------------------------------
   do ideep=2,ndeep
      dz  =deep(ideep)-deep(ideep-1)
      do ilat=1,nlats
        if (mostrf(ilat,ideep-1)/=undef) then
           strmf(ilat,ideep)=strmf(ilat,ideep-1)+mostrf(ilat,ideep-1)
        else
           strmf(ilat,ideep)=undef
        end if
      enddo
   enddo
   mostrf=strmf
end subroutine mosf


subroutine mosf_sig0(u,v,dens,dp,idm,jdm,kdm,mostrf,lats,nlats,sigarray,nsig,undef)
use mod_grid
implicit none

   integer, intent(in) :: idm,jdm,nlats,nsig,kdm
   real,    intent(in), dimension(idm,jdm,kdm) :: u,v,dens,dp
   real,    intent(out) :: mostrf(nlats,nsig)
   real,    intent(in) :: lats(nlats)
   real,    intent(in) :: sigarray(nsig)
   real, intent(in) :: undef
      
   real fumask(idm,jdm), fvmask(idm,jdm)
   logical :: mask(idm,jdm)
   integer :: mostrfcnt(nlats,nsig)
   real :: crssec(nlats,nsig)
   real :: strmf (nlats,nsig)
   real:: dpu, dpv ,transu,transv, densu, densv, dsig
   real :: radian, pi, thref
   integer :: ilat, i, j, im1, ip1, jm1, jp1,k, uindex, vindex, isig
   real, external :: spherdist
   include 'stmt_funcs.H'

   pi = acos(0.)*2
   radian=180/pi
   thref=1e-3

   dsig = (maxval(sigarray)-minval(sigarray))/(nsig-1)


   !print *,'mosf_sig0 needs a tune-up'
   !stop

!
! --- ------------------------------------------------------------------
! --- calculate dp in u and v points
! --- ------------------------------------------------------------------
!
   ! For every zonal band:
   mostrf=0.
   crssec=0.
   do ilat=1,nlats

      ! Find points inside this region (North of lats(i)
      mask=.false.
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jdm
      do i=1,idm
         mask(i,j)=qlat(i,j)>lats(ilat)
      end do
      end do
      !$OMP END PARALLEL DO

      ! This will only retain transport at the boundary of the region
      fumask=0.
      fvmask=0.
      !$OMP PARALLEL DO PRIVATE(i,j,ip1,jp1)
      do j=1,jdm
      do i=1,idm
      if (mask(i,j)) then
         jp1=min(jdm,j+1)
         if (periodic) then
            ip1=mod(i,idm)+1
         else
            ip1=min(idm,i+1)
         end if
         fumask(i,j)=fumask(i,j)+1.
         fvmask(i,j)=fvmask(i,j)+1.
         fumask(ip1,j)=fumask(ip1,j)-1.
         fvmask(i,jp1)=fvmask(i,jp1)-1.
      end if
      end do
      end do
      !$OMP END PARALLEL DO

      ! Decrease mask to active points
      !mask = (abs(fumask)>1e-4 .or. abs(fvmask)>1e-4) !.and. depth(i,j)>1.
      !$OMP PARALLEL DO PRIVATE(i,j)
      do j=1,jdm
      do i=1,idm
         mask(i,j) = abs(fumask(i,j))>1e-4 .or. abs(fvmask(i,j))>1e-4
      end do
      end do
      !$OMP END PARALLEL DO

      do k=1,kdm
         do j=1,jdm
         do i=1,idm
         if (mask(i,j)) then
            jm1=max(j-1,1)
            im1=max(1,i-1)

            dpu = 0.5*(dp(i,j,k) + dp(im1,j,k))
            dpv = 0.5*(dp(i,j,k) + dp(i,jm1,k))

            ! Transport values
            transu=dpu*u(i,j,k)*scuy(i,j)*fumask(i,j)
            transv=dpv*v(i,j,k)*scvx(i,j)*fvmask(i,j)

            densu = .5*(dens(i,j,k)+dens(im1,j,k))
            densv = .5*(dens(i,j,k)+dens(i,jm1,k))
      
            ! Index into vertical dimension
            uindex = max(1,min(floor((densu-minval(sigarray))/dsig)+1,nsig))
            vindex = max(1,min(floor((densv-minval(sigarray))/dsig)+1,nsig))

            ! Update transport array
            if (dpu>.5) then
               mostrf   (ilat,uindex)= mostrf(ilat,uindex) + transu
               mostrfcnt(ilat,uindex)= mostrfcnt(ilat,uindex)+1
            end if
            if (dpv>.5) then
               mostrf   (ilat,vindex)= mostrf(ilat,vindex) + transv
               mostrfcnt(ilat,vindex)= mostrfcnt(ilat,vindex)+1
            end if
         end if
         end do
         end do
      end do ! kdm
   end do ! latitude band ilat

   mostrf=mostrf*1e-6


   ! --- ------------------------------------------------------------------
   ! --- Convert from transport to meridional overturning streamfunction
   ! --- ------------------------------------------------------------------
   do isig=2,nsig
      do ilat=1,nlats
        if (mostrf(ilat,isig-1)/=undef) then
           strmf(ilat,isig)=strmf(ilat,isig-1)+mostrf(ilat,isig-1)
        else
           strmf(ilat,isig)=undef
        end if
      enddo
   enddo
   mostrf=strmf






end subroutine mosf_sig0



subroutine mosf_theta0(u,v,tem,dp,qlon,qlat,idm,jdm,kdm,mostrf,lats,theta0array,nlats,ntheta0,undef)
!use m_spherdist
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
   real, external :: spherdist

   !include 'stmt_funcs.h'
   !pi = acos(0.)*2
   !radian=180/pi
   !thref=1e-3

   print *,'mosf_theta entry'
   thetadel = (maxval(theta0array)-minval(theta0array))/(ntheta0-1)
   !print *,'thetadel=',thetadel

   print *,'mosf_theta needs a tune-up'
   stop
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
!use m_spherdist
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
   real, external :: spherdist

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






