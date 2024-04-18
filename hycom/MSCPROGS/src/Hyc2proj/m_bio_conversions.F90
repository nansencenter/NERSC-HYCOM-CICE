module m_bio_conversions
   real, parameter :: b2=0.04, ny=1.38E-2, N2CHLA=11.0, kd_chl=0.02
   real, parameter :: cnit=14.01,cpho=30.97,csil=28.09, ccar=12.01
   real, parameter :: oxyml=44.6608009,oxygr=32.0,C2NIT=6.625 ! redfield
! _FABM__caglar_
   real, parameter :: kd_chla=0.05737798064012768   ! was 0.4 before Nov23 update, light attenuation coeff. for chlorophyll
   real, parameter :: kd_det=0.16218644594775256   ! was 0.0 before Nov23 update, light attenuation coeff. for detritus
   real, parameter :: kd_dom=0.18206465359221413   ! was 0.0 before Nov23 update, light attenuation coeff. for DOM
   real, parameter :: C2SIL=6.625    ! redfield C:Si mol ratio.
   real, parameter :: C2PHO=106.0    ! redfield C:P mol ratio
   real, parameter :: kd_water=0.041 ! light attenuation coeff. for Arctic sea water 
   real, parameter :: srdet_eco=5.0  ! Sinking rate of detritus
! _FABM__caglar_
   contains

   subroutine chlorophyll_nor(dia,fla,chl_a,idm,jdm,kdm)
!compute chlorophyll: kg m-3
      implicit none
      
      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::dia, fla 
      real, dimension(idm,jdm,kdm)  , intent(out) ::chl_a
      real, dimension(idm,jdm,kdm)  ::boss
      integer :: i,j,k

      chl_a=1.0e-6/N2CHLA*(dia+fla)

   end subroutine chlorophyll_nor

   subroutine chlorophyll_eco(dia,fla,chl_a,idm,jdm,kdm)
!compute chlorophyll: kg m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::dia, fla
      real, dimension(idm,jdm,kdm)  , intent(out) ::chl_a
      real, dimension(idm,jdm,kdm)  ::boss
      integer :: i,j,k

      chl_a=1.0e-6*(dia+fla)

   end subroutine chlorophyll_eco
!------------------------------------------------
   subroutine attenuation_nor(dia,fla,attencoef,idm,jdm,kdm)
!compute the attenuation coefficient: m-1
      implicit none
      
      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::dia, fla 
      real, dimension(idm,jdm,kdm)  , intent(out) ::attencoef
    
      attencoef=b2+ny/N2CHLA*(dia+fla)

   end subroutine attenuation_nor

   subroutine attenuation_eco(dia,fla,attencoef,idm,jdm,kdm)
!compute the attenuation coefficient: m-1
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::dia, fla
      real, dimension(idm,jdm,kdm)  , intent(out) ::attencoef

      attencoef=b2+kd_chl*(dia+fla)

   end subroutine attenuation_eco
!------------------------------------------------
   subroutine nitrate_conv_nor(nit,nitrate,idm,jdm,kdm)
!compute nitrate: mole m-3
      implicit none
      
      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::nit !mmol m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::nitrate

      nitrate=1e-3/cnit*nit

   end subroutine nitrate_conv_nor

   subroutine nitrate_conv_eco(nit,nitrate,idm,jdm,kdm)
!compute nitrate: mole m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::nit !mmol m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::nitrate

      nitrate=1e-3*nit

   end subroutine nitrate_conv_eco
!------------------------------------------------
   subroutine phosphate_conv_nor(pho,phosphate,idm,jdm,kdm)
!compute phosphate: mole m-3
      implicit none
      
      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::pho !mmol m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::phosphate
    
      phosphate=1e-3/cpho*pho

   end subroutine phosphate_conv_nor

   subroutine phosphate_conv_eco(pho,phosphate,idm,jdm,kdm)
!compute phosphate: mole m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::pho !mmol m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::phosphate

      phosphate=1e-3*pho

   end subroutine phosphate_conv_eco
!------------------------------------------------
   subroutine silicate_conv_eco(sil,silicate,idm,jdm,kdm)
!compute silicate: mole m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::sil !mmol m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::silicate

      silicate=1e-3*sil

   end subroutine silicate_conv_eco
!------------------------------------------------
   subroutine pbiomass_nor(dia,fla,biomass,idm,jdm,kdm)
!compute phytoplankton biomass: mole m-3
      implicit none
      
      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::fla, dia !mg P m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::biomass
    
      biomass=1e-3/cnit*(fla+dia)

   end subroutine pbiomass_nor

   subroutine pbiomass_eco(dia,fla,biomass,idm,jdm,kdm)
!compute phytoplankton biomass: mole m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::fla, dia !mg P m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::biomass

      biomass=1e-3/ccar/C2NIT*(fla+dia)

   end subroutine pbiomass_eco
!------------------------------------------------
   subroutine oxygen_conv_nor(oxy,kg_oxy,idm,jdm,kdm)
!compute dissolved oxygen: kg m-3
      implicit none
      
      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::oxy !mg 02 m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::kg_oxy
! original oxygen unit mg/m3
      kg_oxy=1e-6*oxy

   end subroutine oxygen_conv_nor

   subroutine oxygen_conv_eco(oxy,kg_oxy,idm,jdm,kdm)
!compute dissolved oxygen: kg m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::oxy !ml 02 l-1
      real, dimension(idm,jdm,kdm)  , intent(out) ::kg_oxy
! original oxygen unit ml/l
      kg_oxy=1e-6*oxygr*oxyml*oxy

   end subroutine oxygen_conv_eco

!------------------------------------------------
   subroutine primary_production_nor(primprod,pres,gpp_depthint,onem,idm,jdm,kdm)
!compute gross primary production: kg m-2 s-1
      implicit none
      
      integer, intent(in) :: idm,jdm,kdm
      real, intent(in) :: onem
      real, dimension(idm,jdm,kdm)  , intent(in)  ::primprod, pres
      real, dimension(idm,jdm)      , intent(out) ::gpp_depthint

      real, dimension(idm,jdm,kdm)   :: dplayer
      
      integer :: i,j
      real, dimension(idm,jdm,kdm)                ::gpp

! original primprod unit unit: mg N m-3 s-1
      gpp=1e-6*ccar*C2NIT*primprod/cnit !now in kg C m-3 s-1
! calculate layer depth in meters
      do i=1,kdm
       dplayer(:,:,i)=(pres(:,:,i+1)-pres(:,:,i))/onem
      end do

!intergarte over depth
      do i=1,idm
        do j=1,jdm
          gpp_depthint(i,j)=dot_product(gpp(i,j,:),dplayer(i,j,:))
        end do
      end do

   end subroutine primary_production_nor

   subroutine primary_production_eco(primprod,pres,gpp_depthint,onem,idm,jdm,kdm)
!compute gross primary production: kg m-2 s-1
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, intent(in) :: onem
      real, dimension(idm,jdm,kdm)  , intent(in)  ::primprod, pres
      real, dimension(idm,jdm)      , intent(out) ::gpp_depthint

      real, dimension(idm,jdm,kdm)   :: dplayer

      integer :: i,j,k
      real, dimension(idm,jdm,kdm)                ::gpp

! original primprod unit unit: mg C m-3 s-1
      gpp=1e-6*primprod/86400. !now in kg C m-3 s-1
! calculate layer depth in meters
      do i=1,kdm
       dplayer(:,:,i)=(pres(:,:,i+1)-pres(:,:,i))/onem
      end do

       gpp_depthint = 0.0
!intergarte over depth
      do i=1,idm
        do j=1,jdm
!          do k=1,kdm
          gpp_depthint(i,j)=dot_product(gpp(i,j,:),dplayer(i,j,:))
!       gpp_depthint(i,j)=gpp_depthint(i,j)+gpp(i,j,k)*dplayer(i,j,k)
!          end do
        end do
      end do

   end subroutine primary_production_eco

   subroutine primary_production(primprod,pres,gpp_depthint,onem,idm,jdm,kdm)
!compute net primary production: mgC m-2 d-1
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, intent(in) :: onem
      real, dimension(idm,jdm,kdm)  , intent(in)  ::primprod, pres
      real, dimension(idm,jdm)      , intent(out) ::gpp_depthint

      real, dimension(idm,jdm,kdm)   :: dplayer

      integer :: i,j,k
      real, dimension(idm,jdm,kdm)                ::gpp

! original primprod unit unit: mg C m-3 s-1 (gross pp)
      gpp=primprod*86400. !now in mg C m-3 d-1
! calculate layer depth in meters
      do i=1,kdm
       dplayer(:,:,i)=(pres(:,:,i+1)-pres(:,:,i))/onem
      end do

       gpp_depthint = 0.0
!intergarte over depth
      do i=1,idm
        do j=1,jdm
!          do k=1,kdm
          gpp_depthint(i,j)=dot_product(gpp(i,j,:),dplayer(i,j,:))
          gpp_depthint(i,j)=gpp_depthint(i,j) * 0.9 ! assumed 10% respiration
!       gpp_depthint(i,j)=gpp_depthint(i,j)+gpp(i,j,k)*dplayer(i,j,k)
!          end do
        end do
      end do

   end subroutine primary_production


!------------------------------------------------
   subroutine net_primary_production(netpp,dia,fla,pres,npp_euphd,onem,idm,jdm,kdm)
!compute net primary production: g m-2 day-1
      implicit none
      
      integer, intent(in) :: idm,jdm,kdm
      real, intent(in) :: onem
      real, dimension(idm,jdm,kdm)  , intent(in)  ::netpp, pres, fla, dia
      real, dimension(idm,jdm)      , intent(out) ::npp_euphd

      real, dimension(idm,jdm,kdm)   :: dplayer, ldepth, kpar
      
      integer :: i,j,k
      real, dimension(idm,jdm,kdm)   ::npp
      real                           ::partop, parbot, z_euphd 
      logical                        ::reached_bot           

! original primprod unit unit: mg N m-3 s-1
      npp=1e-3*ccar*C2NIT*netpp*86400.0/cnit !now in g C m-3 day-1
! calculate layer depth in meters
      do i=1,kdm
       dplayer(:,:,i)=(pres(:,:,i+1)-pres(:,:,i))/onem
       ldepth(:,:,i)=pres(:,:,i)/onem
       kpar=b2+(dia+fla)*ny/N2CHLA
      end do
      

!intergarte over depth
      npp_euphd=0.0
      do i=1,idm
        do j=1,jdm
          parbot=1.0
          reached_bot=.false.
          k=1
          do while (reached_bot .eqv..false.)
            partop=parbot
            parbot=partop*exp(-dplayer(i,j,k)*kpar(i,j,k))
!if(i>100.and.j>100) then
!            print*,-dplayer(i,j,k), kpar(i,j,k),partop, parbot
!end if
            if (parbot>0.01) then
              npp_euphd(i,j)=npp_euphd(i,j)+npp(i,j,k)*dplayer(i,j,k)
            else
	      z_euphd=-log(0.01/partop)/kpar(i,j,k)
!if(i>100.and.j>100) then
!              print*, z_euphd, ldepth(i,j,k-1:k+1)
!              pause
!end if
              npp_euphd(i,j)=npp_euphd(i,j)+npp(i,j,k)*(z_euphd)
              reached_bot=.true.
            end if
            if (k>=kdm-1) reached_bot=.true.
            k=k+1;
          end do
        end do
      end do

   end subroutine net_primary_production
!------------------------------------------------
   subroutine integrated_chlorophyll(dia,fla,pres,chl_euphd,onem,idm,jdm,kdm)
!compute net primary production: mg m-2 
      implicit none
      
      integer, intent(in) :: idm,jdm,kdm
      real, intent(in) :: onem
      real, dimension(idm,jdm,kdm)  , intent(in)  ::pres, fla, dia 
      real, dimension(idm,jdm)      , intent(out) ::chl_euphd

      real, dimension(idm,jdm,kdm)   :: dplayer, ldepth, kpar
      
      integer :: i,j,k
      real, dimension(idm,jdm,kdm)   ::chl
      real                           ::partop, parbot, z_euphd 
      logical                        ::reached_bot           

! original phytoplankton unit unit: mg N m-3
      chl=(dia+fla)/N2CHLA
! calculate layer depth in meters
      do i=1,kdm
       dplayer(:,:,i)=(pres(:,:,i+1)-pres(:,:,i))/onem
       ldepth(:,:,i)=pres(:,:,i)/onem
       kpar=b2+(dia+fla)*ny/N2CHLA
      end do
      

!intergarte over depth
      chl_euphd=0.0
      do i=1,idm
        do j=1,jdm
          parbot=1.0
          reached_bot=.false.
          k=1
          do while (reached_bot .eqv. .false.)
            partop=parbot
            parbot=partop*exp(-dplayer(i,j,k)*kpar(i,j,k))
            if (parbot>0.37) then
              chl_euphd(i,j)=chl_euphd(i,j)+chl(i,j,k)*dplayer(i,j,k)
            else
	      z_euphd=-log(0.37/partop)/kpar(i,j,k)
              chl_euphd(i,j)=chl_euphd(i,j)+chl(i,j,k)*(z_euphd)
              reached_bot=.true.
            end if
            if (k>=kdm-1) reached_bot=.true.
            k=k+1;
          end do
        end do
      end do

   end subroutine integrated_chlorophyll

   subroutine integrated_chlorophyll_eco(dia,fla,pres,chl_euphd,onem,idm,jdm,kdm)
!compute net primary production: mg m-2 
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, intent(in) :: onem
      real, dimension(idm,jdm,kdm)  , intent(in)  ::pres, fla, dia
      real, dimension(idm,jdm)      , intent(out) ::chl_euphd

      real, dimension(idm,jdm,kdm)   :: dplayer, ldepth, kpar

      integer :: i,j,k
      real, dimension(idm,jdm,kdm)   ::chl
      real                           ::partop, parbot, z_euphd
      logical                        ::reached_bot

! original phytoplankton unit unit: mg Chl m-3
      chl=(dia+fla)
! calculate layer depth in meters
      do i=1,kdm
       dplayer(:,:,i)=(pres(:,:,i+1)-pres(:,:,i))/onem
       ldepth(:,:,i)=pres(:,:,i)/onem
       kpar=b2+(dia+fla)*kd_chl
      end do


!intergarte over depth
      chl_euphd=0.0
      do i=1,idm
        do j=1,jdm
          parbot=1.0
          reached_bot=.false.
          k=1
          do while (reached_bot .eqv. .false.)
            partop=parbot
            parbot=partop*exp(-dplayer(i,j,k)*kpar(i,j,k))
            if (parbot>0.37) then
              chl_euphd(i,j)=chl_euphd(i,j)+chl(i,j,k)*dplayer(i,j,k)
            else
              z_euphd=-log(0.37/partop)/kpar(i,j,k)
              chl_euphd(i,j)=chl_euphd(i,j)+chl(i,j,k)*(z_euphd)
              reached_bot=.true.
            end if
            if (k>=kdm-1) reached_bot=.true.
            k=k+1;
          end do
        end do
      end do

   end subroutine integrated_chlorophyll_eco
!------------------------------------------------

   subroutine det_bottom_flux(det,dsnk,bot_flux,onem,idm,jdm,kdm)
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, intent(in) :: onem
      real, dimension(idm,jdm,kdm)  , intent(in)  :: det,dsnk
      real, dimension(idm,jdm,kdm)      , intent(out) :: bot_flux
      real :: spd(idm,jdm,kdm)
      integer :: i,j,k
! dsnk is the varible detritus sinking speed * detritus concentration
! essentially it is the flux represented as a state variable
! To prevent outlier values, we first calculate the sinking speed (dsnk/det) and
! convert it to 1/d
      spd = (dsnk/det)*86400.
! and set minimum and maximum values, and convert mgC m-2 d-1 --> mmolC m-2 d-1 
      spd = max(spd,0.5)
      spd = min(spd,12.0)
      bot_flux = det * spd / ccar
       

!! compute flux of detritus to the seafloor mgC m-2 d-1 --> mmolC m-2 d-1                                                                                       !      bot_flux=det * srdet_eco / ccar 

   end subroutine det_bottom_flux

!------------------------------------------------

! _FABM__caglar_
      subroutine chlorophyll_fabm(dia,fla,ccl,chl_a,idm,jdm,kdm)
      !compute chlorophyll: mg m-3
        implicit none

        integer, intent(in) :: idm,jdm,kdm
        real, dimension(idm,jdm,kdm)  , intent(in)  ::dia, fla, ccl
        real, dimension(idm,jdm,kdm)  , intent(out) ::chl_a
        real, dimension(idm,jdm,kdm)  ::boss
        integer :: i,j,k

        chl_a=dia+fla+ccl

     end subroutine chlorophyll_fabm

     subroutine chlorophyll(dia,fla,ccl,chl_a,idm,jdm,kdm)
!compute chlorophyll: mg m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::dia, fla, ccl
      real, dimension(idm,jdm,kdm)  , intent(out) ::chl_a
      real, dimension(idm,jdm,kdm)  ::boss
      integer :: i,j,k

      chl_a=(dia+fla+ccl)

     end subroutine chlorophyll


     subroutine nitrate_conv(nit,nitrate,idm,jdm,kdm)
     !compute nitrate: mmole m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::nit !mgC m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::nitrate

      nitrate=nit/ccar/C2NIT

     end subroutine nitrate_conv

     subroutine silicate_conv(sil,silicate,idm,jdm,kdm)
     !compute silicate: mmole m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::sil !mgC m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::silicate

      silicate=sil/ccar/C2SIL

     end subroutine silicate_conv

     subroutine phosphate_conv(pho,phosphate,idm,jdm,kdm)
     !compute phosphate: mmole m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::pho !mgC m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::phosphate

      phosphate=pho/ccar/C2PHO

     end subroutine phosphate_conv

     subroutine pbiomass(dia,fla,ccl,biomass,idm,jdm,kdm)
     !compute phytoplankton biomass: mmoleC m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::fla, dia, ccl !mgC m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::biomass

      biomass=(fla+dia+ccl)/ccar

     end subroutine pbiomass

     subroutine zbiomass(micro,meso,biomass,idm,jdm,kdm)
     !compute zooplankton biomass: mmoleC m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::micro, meso !mgC m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::biomass

      biomass=(micro+meso)/ccar

     end subroutine zbiomass

     subroutine oxygen_conv(oxy,mmol_oxy,idm,jdm,kdm)
!compute dissolved oxygen: mmol m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::oxy !mmol 02 m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::mmol_oxy
! original oxygen unit mmol/m3
      mmol_oxy=oxy

     end subroutine oxygen_conv

     subroutine pp_conv(pp,pp_daily,idm,jdm,kdm)
!compute net primary production (gross PP * 0.9 as net pp): mg m-3 d-1
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::pp ! mg m-3 s-1
      real, dimension(idm,jdm,kdm)  , intent(out) ::pp_daily

      pp_daily=pp*24.*60.*60.

     end subroutine pp_conv

     subroutine attenuation(dia,fla,ccl, det, dom,attencoef,idm,jdm,kdm)
!compute the attenuation coefficient: m-1
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::dia, fla, ccl, det, dom
      real, dimension(idm,jdm,kdm)  , intent(out) ::attencoef

      attencoef=kd_water+kd_chl*(dia+fla+ccl)+(kd_det/1000.0)*det+(kd_dom/1000.0)*dom

     end subroutine attenuation

     subroutine dic_conv(dic,dissic,idm,jdm,kdm)
     !compute dic: mole m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::dic !mmol m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::dissic

      dissic=dic/1000.

     end subroutine dic_conv

     subroutine pco2_conv(spco2_ppm,spco2,idm,jdm,kdm)
     !compute dic: mole m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::spco2_ppm !micromol mol-1 (ppm)
      real, dimension(idm,jdm,kdm)  , intent(out) ::spco2 ! Pa
      ! conversion taken from:
      ! https://acsess.onlinelibrary.wiley.com/doi/pdfdirect/10.2134/asaspecpub53.appendix2

      spco2=spco2_ppm/10.1325

     end subroutine pco2_conv

! _FABM__caglar_
! _vertical_velocity

   subroutine vertical_velocity(u,v,pres,w,scpx,scpy,scux,scvy,plon,plat,depth,onem,idm,jdm,kdm)
!diagnose vertical velocity: m/s
      implicit none
      
      integer, intent(in) :: idm,jdm,kdm
      real, intent(in) :: onem
      real, dimension(idm,jdm,kdm)  , intent(in)    :: u,v
      real, dimension(idm,jdm,kdm+1)  , intent(in)  :: pres
      real, dimension(idm,jdm,kdm+1)  :: lpres
      real, dimension(idm,jdm)  , intent(in)        :: scpx,scpy,scux,scvy,plon,plat,depth
      real, dimension(idm,jdm,kdm)      , intent(out) :: w
      real, dimension(idm,jdm,kdm)                :: dudx,dvdy,dpdx,dpdy,layer_thkn,int_div,thk_adv
      real :: idsum, tasum

      integer :: i,j,k

      lpres=pres/onem;
      ! caulculate derivative:
      do i=2,idm-1
        do j=2,jdm-1
          if (depth(i,j)>0) then
            do k=1,kdm
              dudx(i,j,k)=(u(i+1,j,k)-u(i-1,j,k))/(scux(i-1,j)+scux(i,j));
              dvdy(i,j,k)=(v(i,j+1,k)-v(i,j-1,k))/(scvy(i,j-1)+scvy(i,j));

              dpdx(i,j,k)=(lpres(i+1,j,k+1)-lpres(i-1,j,k+1))/(scpx(i-1,j)+scpx(i,j));
              dpdy(i,j,k)=(lpres(i,j+1,k+1)-lpres(i,j-1,k+1))/(scpy(i,j-1)+scpy(i,j));
            end do
          end if
        end do
      end do
!      print*, dudx(50,50,:)
!      print*, dvdy(50,50,:)
!      print*, dpdx(50,50,:)
!      print*, dpdy(50,50,:)
!  calculate at the bottom of the first layer:
      do i=2,idm-1
        do j=2,jdm-1
          if (depth(i,j)>0) then
            do k=1,kdm
              layer_thkn(i,j,k)=lpres(i,j,k+1)-lpres(i,j,k);
            end do
          end if
        end do
      end do

!  calculate the thickness integrated divergence in each layer
      do i=2,idm-1
        do j=2,jdm-1
          if (depth(i,j)>0) then
            do k=1,kdm
              int_div(i,j,k)=layer_thkn(i,j,k)*(dudx(i,j,k)+dvdy(i,j,k))
            end do
          end if
        end do
      end do

!  calculate thickness advection in each layer:
      thk_adv(:,:,1)=0.0
      do i=2,idm-1
        do j=2,jdm-1
          if (depth(i,j)>0) then
            do k=2,kdm
              thk_adv(i,j,k)=(u(i,j,k)-u(i,j,k-1))*dpdx(i,j,k) + &
                             (v(i,j,k)-v(i,j,k-1))*dpdy(i,j,k);

            end do
          end if
        end do
      end do

!  evaluate the vertical veolocity at the mid-point in each layer
      do i=2,idm-1
        do j=2,jdm-1
          if (depth(i,j)>0) then
            w(i,j,1)=0.5*int_div(i,j,1);  
            idsum=0.0
            tasum=0.0
            do k=2,kdm
              idsum=idsum+int_div(i,j,k)
              tasum=tasum+thk_adv(i,j,k-1)
              w(i,j,k)=idsum-tasum+0.5*int_div(i,j,k)
            end do
          end if
        end do
      end do

   end subroutine vertical_velocity
!   
   subroutine w_velocity(u,v,pres,w,scpx,scpy,scux,scuy,scvx,scvy,plon,plat,depth,onem,idm,jdm,kdm)
!-- from archv2ncdf3      
!diagnose vertical velocity: m/s
      implicit none
      
      integer, intent(in) :: idm,jdm,kdm
      real, intent(in) :: onem
      real, dimension(idm,jdm,kdm)  , intent(in)    :: u,v
      real, dimension(idm,jdm,kdm+1)  , intent(in)  :: pres
      real, dimension(idm,jdm,kdm+1)  :: lpres
      real, dimension(idm,jdm)  , intent(in)        :: scpx,scpy,scux,scuy,scvx,scvy,plon,plat,depth
      real, dimension(idm,jdm,kdm)    , intent(out) :: w
      real, dimension(idm,jdm)                      :: dpdx,dpdy, work
      real             dudxdn,dudxup,dvdydn,dvdyup        
      real, dimension(idm,jdm,kdm)   :: dplayer
      !real :: idsum, tasum

      real, parameter :: flag = 2.0**100
      integer :: i,j,k

      lpres=pres/onem;
! calculate layer depth in meters
      do k=1,kdm
       dplayer(:,:,k)=(lpres(:,:,k+1)-lpres(:,:,k))
      end do
      
!      if     (iowvlin.ne.0) then
      do j= 1,jdm-1
        do i= 1,idm-1
          if (depth(i,j)>0) then
          !if     (ip(i,j).eq.1) then
            dudxdn= &
                 (u(i+1,j  ,1)*scuy(i+1,j  )-u(i,j,1)*scuy(i,j))&
                 /(scpx(i,j)*scpy(i,j))
            dvdydn= &
                 (v(i  ,j+1,1)*scvx(i  ,j+1)-v(i,j,1)*scvx(i,j))&
                 /(scpx(i,j)*scpy(i,j))
            !w(i,j,1)=     dplayer(i,j,1)*(dudxdn+dvdydn) ! layer interface intfwv=1
            w(i,j,1)=0.5*dplayer(i,j,1)*(dudxdn+dvdydn) ! layer center  intfwv=0
          else
            w(i,j,1) = flag
          endif
        enddo
      enddo
      do k= 2,kdm
        do j= 1,jdm-1
          do i= 1,idm-1
            !if     (iu(i,j).eq.1) then
            if (depth(i,j)>0) then
              dpdx(i,j)=&
                     (lpres(i,j,k)*scpy(i,j)-lpres(i-1,j  ,k)*scpy(i-1,j  ))&
                     /(scux(i,j)*scuy(i,j))
            endif
          enddo
        enddo
        do j= 1,jdm-1
          do i= 1,idm-1
            if (depth(i,j)>0) then
            !if     (iv(i,j).eq.1) then
              dpdy(i,j)=&
                     (lpres(i,j,k)*scpx(i,j)-lpres(i  ,j-1,k)*scpx(i  ,j-1))&
                     /(scvx(i,j)*scvy(i,j))
            endif
          enddo
        enddo
        do j=1,jdm
          !if     (iu(2   ,j).eq.1) then
          if (depth(2,j)>0) then
            dpdx(1 ,j)=dpdx(2   ,j)
          endif
          if (depth(idm-1,j)>0) then
          !if     (iu(ii-1,j).eq.1) then
            dpdx(idm,j)=dpdx(idm-1,j)
          endif
        enddo
        do i=1,idm
          if (depth(i,2)>0) then
          !if     (iv(i,2   ).eq.1) then
            dpdy(i,1 )=dpdy(i,2   )
          endif
          if (depth(i,jdm-1)>0) then
          !if     (iv(i,jj-1).eq.1) then
            dpdy(i,jdm)=dpdy(i,jdm-1)
          endif
        enddo
        do j= 1,jdm-1
          do i= 1,idm-1
            if (depth(i,j)>0) then
            !if     (ip(i,j).eq.1) then
              dudxup=&
                   (u(i+1,j  ,k-1)*scuy(i+1,j  )-&
                    u(i  ,j  ,k-1)*scuy(i  ,j  ))&
                   /(scpx(i,j)*scpy(i,j))
              dvdyup=&
                    (v(i  ,j+1,k-1)*scvx(i  ,j+1)-&
                    v(i  ,j  ,k-1)*scvx(i  ,j  ))&
                   /(scpx(i,j)*scpy(i,j))
              dudxdn=&
                   (u(i+1,j  ,k  )*scuy(i+1,j  )-&
                    u(i  ,j  ,k  )*scuy(i   ,j  ))&
                   /(scpx(i,j)*scpy(i,j)) 
              dvdydn=&                    
                   (v(i  ,j+1,k  )*scvx(i   ,j+1)-&
                    v(i  ,j  ,k  )*scvx(i   ,j  ))&
                   /(scpx(i,j)*scpy(i,j))
             !w(i,j,k)=w(i,j,k-1)+0.5*(2.0*dplayer(i,j,k)*(dudxdn+dvdydn)-&  !intfwv=1
             !         (u(i  ,j  ,k)-u(i  ,j  ,k-1))*dpdx(i  ,j  )-&
             !         (u(i+1,j  ,k)-u(i+1,j  ,k-1))*dpdx(i+1,j  )-&
             !         (v(i  ,j  ,k)-v(i  ,j  ,k-1))*dpdy(i  ,j  )-&
             !         (v(i  ,j+1,k)-v(i  ,j+1,k-1))*dpdy(i  ,j+1))
             w(i,j,k)=w(i,j,k-1)+0.5*(dplayer(i,j,k-1)*(dudxup+dvdyup)+& !intfwv=0 
                                      dplayer(i,j,k  )*(dudxdn+dvdydn)-&
                      (u(i  ,j  ,k)-u(i  ,j  ,k-1))*dpdx(i  ,j  )-&
                      (u(i+1,j  ,k)-u(i+1,j  ,k-1))*dpdx(i+1,j  )-&
                      (v(i  ,j  ,k)-v(i  ,j  ,k-1))*dpdy(i  ,j  )-&
                      (v(i  ,j+1,k)-v(i  ,j+1,k-1))*dpdy(i  ,j+1))
            else
              w(i,j,k) = flag
            endif
          enddo
        enddo
        do i= 1,idm
          w(i ,jdm,k) = flag
        enddo
        do j= 1,jdm
          w(idm,j ,k) = flag
        enddo
      enddo !k
!
! --- w is noisy - smooth at least twice.
      do k= 1,kdm
        !call psmoo(w(1,1,k),work)
        !call psmoo(w(1,1,k),work)
        work(:,:)=w(:,:,k)
        w(:,:,k)=psmoo(work(:,:),idm,jdm)
        work(:,:)=w(:,:,k)
        w(:,:,k)=psmoo(work(:,:),idm,jdm)
      enddo
!     !iowvlin 
   end subroutine w_velocity
!---   
   function psmoo(alist,idm,jdm)
!
      integer idm,jdm
      real alist(idm,jdm),blist(idm,jdm), psmoo(idm,jdm)
      real, parameter :: wgt  = 0.25
      real, parameter :: flag = 2.0**100
      integer i,ia,ib,j,ja,jb
!     copy the code that does the smmothing from archv2ncdf3
!     entry psmoo(alist,blist)
!c ---this entry is set up to smooth data carried at -p- points
!c
      do i=1,idm
         do j=1,jdm
            if (alist(i,j).ne.flag) then
               !ja=max(jfp(i,l),j-1)
               !jb=min(jlp(i,l),j+1)
               ja=j-1
               jb=j+1
               if (alist(i,ja).eq.flag) ja=j
               if (alist(i,jb).eq.flag) jb=j
               blist(i,j)=(1.-wgt-wgt)*alist(i,j)+wgt*(alist(i,ja)+alist(i,jb))
            endif
         enddo
      enddo
!
      do j=1,jdm
         do i=1,idm
            if (alist(i,j).ne.flag) then
               !ia=max(ifp(j,l),i-1)
               !ib=min(ilp(j,l),i+1)
               ia=i-1
               ib=i+1
               if (alist(ia,j).eq.flag) ia=i
               if (alist(ib,j).eq.flag) ib=i
               alist(i,j)=(1.-wgt-wgt)*blist(i,j)+wgt*(blist(ia,j)+blist(ib,j))
            endif
         enddo
      enddo
      psmoo(:,:)=alist(:,:)
   end
!  function std(x,ldim)
!  real x(ldim)
!  xmean=sum(x)/ldim
!  std=sqrt(sum((x-xmean)*(x-xmean))/(ldim-1))
!  end
!------------------------------------------------
end module








