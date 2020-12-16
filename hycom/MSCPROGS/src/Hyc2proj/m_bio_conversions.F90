module m_bio_conversions
   real, parameter :: b2=0.04, ny=1.38E-2, N2CHLA=11.0, kd_chl=0.02
   real, parameter :: cnit=14.01,cpho=30.97,csil=28.09, ccar=12.01
   real, parameter :: oxyml=44.6608009,oxygr=32.0,C2NIT=6.625 ! redfield
! _FABM__caglar_
   real, parameter :: kd_fabm=0.04   ! light attenuation coeff. for chlorophyll
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

   subroutine det_bottom_flux(det,bot_flux,onem,idm,jdm,kdm)
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, intent(in) :: onem
      real, dimension(idm,jdm,kdm)  , intent(in)  :: det
      real, dimension(idm,jdm,kdm)      , intent(out) :: bot_flux

      integer :: i,j,k

! compute flux of detritus to the seafloor mgC m-2 d-1 --> mmolC m-2 d-1                                                                                           
      bot_flux=det * srdet_eco / ccar 

   end subroutine det_bottom_flux

!------------------------------------------------

! _FABM__caglar_
      subroutine chlorophyll_fabm(dia,fla,chl_a,idm,jdm,kdm)
      !compute chlorophyll: mg m-3
        implicit none

        integer, intent(in) :: idm,jdm,kdm
        real, dimension(idm,jdm,kdm)  , intent(in)  ::dia, fla
        real, dimension(idm,jdm,kdm)  , intent(out) ::chl_a
        real, dimension(idm,jdm,kdm)  ::boss
        integer :: i,j,k

        chl_a=dia+fla

     end subroutine chlorophyll_fabm

     subroutine chlorophyll(dia,fla,chl_a,idm,jdm,kdm)
!compute chlorophyll: mg m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::dia, fla
      real, dimension(idm,jdm,kdm)  , intent(out) ::chl_a
      real, dimension(idm,jdm,kdm)  ::boss
      integer :: i,j,k

      chl_a=(dia+fla)

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

     subroutine pbiomass(dia,fla,biomass,idm,jdm,kdm)
     !compute phytoplankton biomass: mmoleC m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::fla, dia !mgC m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::biomass

      biomass=(fla+dia)/ccar

     end subroutine pbiomass

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

      pp_daily=pp*24.*60.*60.*0.9

     end subroutine pp_conv

     subroutine attenuation(dia,fla,attencoef,idm,jdm,kdm)
!compute the attenuation coefficient: m-1
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::dia, fla
      real, dimension(idm,jdm,kdm)  , intent(out) ::attencoef

      attencoef=kd_water+kd_fabm*(dia+fla)

     end subroutine attenuation

     subroutine dic_conv(dic,dissic,idm,jdm,kdm)
     !compute dic: mole m-3
      implicit none

      integer, intent(in) :: idm,jdm,kdm
      real, dimension(idm,jdm,kdm)  , intent(in)  ::dic !mmol m-3
      real, dimension(idm,jdm,kdm)  , intent(out) ::dissic

      dissic=dic/1000.

     end subroutine dic_conv


! _FABM__caglar_
end module








