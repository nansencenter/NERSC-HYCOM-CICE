module m_read_field_wrapper

! The main idea here is to clean up the code a bit in eg p_hyc2proj.F90

contains

   subroutine read_field_wrapper_2d(fld,cvar,fbase,idm,jdm,k,hycfile,nrmem,undef,dplayer,silent)
   use mod_read_rstab
   use mod_read_dailyab
   use mod_read_weekab
   use mod_readpak
   implicit none
   integer, intent(in) :: idm,jdm,nrmem,hycfile,k
   character*8, intent(in) :: cvar
   character(len=*), intent(in) :: fbase
   real   , intent(in) :: undef
   real   , intent(out) :: fld(idm,jdm)
   real   , intent(in ), optional :: dplayer(idm,jdm)
   logical, intent(in ), optional :: silent

   logical :: lsilent

   lsilent=.true.
   if (present(silent)) lsilent=silent

   if (hycfile==2) then
      call read_field2d    (trim(fbase),cvar,fld   ,idm,jdm,k,undef)
   elseif (hycfile==3) then
      call read_weekfield2d(trim(fbase),cvar,fld ,idm,jdm,k,undef,lsilent=lsilent)
      if (present(dplayer) .and. trim(cvar)/='pres')  then
         fld=fld/max(dplayer,1e-6)
      end if
   elseif (hycfile==1) then
      ! NB - not present in restart file
      call read_rstfield2d (trim(fbase),cvar,fld  ,idm,jdm,k,undef)
   elseif (hycfile==4) then
      call get_pakfield(trim(fbase),trim(cvar),k,fld,idm,jdm,undef)
   end if
   where(fld/=undef) fld=fld/nrmem
   end subroutine


   subroutine read_field_wrapper_3d(fld,cvar,fbase,idm,jdm,kdm,hycfile,nrmem,undef,pres)
   use mod_read_rstab
   use mod_read_dailyab
   use mod_read_weekab
   use mod_readpak
   implicit none
   integer, intent(in) :: idm,jdm,kdm,nrmem,hycfile
   character*8, intent(in) :: cvar
   character(len=*), intent(in) :: fbase
   real   , intent(in) :: undef
   real   , intent(out) :: fld(idm,jdm,kdm)
   real   , intent(in ) :: pres(idm,jdm,kdm)
   real :: dplayer(idm,jdm)
   integer :: k
   real, parameter :: onem=9806.

   if (hycfile==2) then
      call read_field3d    (trim(fbase),cvar,fld   ,idm,jdm,kdm,undef)
   elseif (hycfile==3) then
      call read_weekfield3d(trim(fbase),cvar,fld   ,idm,jdm,kdm,undef)

      ! Fields are layer-weighted
      if (trim(cvar)/='pres') then
         do k=1,kdm
            dplayer = pres(:,:,k+1)-pres(:,:,k)
            fld(:,:,k)=fld(:,:,k)*onem/max(dplayer,1.)
         end do
      end if

   elseif (hycfile==1) then
      ! NB - not present in restart file
      call read_rstfield3d (trim(fbase),cvar,fld  ,idm,jdm,kdm,undef)
   elseif (hycfile==4) then
      do k=1,kdm
         call get_pakfield(trim(fbase),trim(cvar),k,fld(:,:,k),idm,jdm,undef)
      end do
   end if
   where(fld/=undef) fld=fld/nrmem
   end subroutine

end module m_read_field_wrapper
