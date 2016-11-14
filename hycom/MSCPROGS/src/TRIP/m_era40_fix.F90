module m_era40_fix

! Moved this into seperate subroutine, since this fix is needed for
!synoptic and climatology fields
contains

   subroutine era40_fix(vname,fld,nlon,nlat,flon,flat,dlon,dlat)
      use mod_xc
      implicit none

      character(len=*), intent(in)    :: vname
      real            , intent(in)    :: flon,flat,dlon,dlat
      integer         , intent(in)    :: nlon,nlat
      real            , intent(inout) :: fld(nlon,nlat)

      ! Variable names - ERA40 1.125 by 1.125 
      character(len=*), parameter :: vuwnd='U10M_sfc'
      character(len=*), parameter :: vvwnd='V10M_sfc'
      character(len=*), parameter :: vtair='T2M_sfc'
      character(len=*), parameter :: vdewp='D2M_sfc'
      character(len=*), parameter :: vprec='TP'
      character(len=*), parameter :: vtcc ='TCC_sfc'
      character(len=*), parameter :: vblh ='BLH_sfc'
      character(len=*), parameter :: vlmsk='LSM_sfc'
      character(len=*), parameter :: vpres='MSL_sfc'
      character(len=*), parameter :: vssrd='SSRD_sfc'

      real :: latera40, weight,wrad,wrad2
      integer :: i,j

      if (vname==vprec.or.trim(vname)=='RO') then
         ! ERA40 precipitation is biased in the tropics (See Biasotti et al., J.Clim)
         ! bias is approx 50 % so real precip is ~ 2/3 of precip fields
         ! This puts a latitudinal weighting on precipitation, centered on the
         ! mean ITCZ-position (~5 degrees N). TODO: Adjust to seasonal ITCZ
         ! position
         do j=1,nlat
            latera40=flat+(j-1)*dlat

            ! Weight=1 at 5N, 0 at 30N, 25S
            weight=abs(latera40-5.)
            wrad=20. ! 0-weight distance from 5N
            wrad2=10. ! region around 5N where weight=1
            if (weight<wrad2) then
               weight=1.
            else if (weight<wrad) then
               weight=.5+.5*cos((weight-wrad2)*3.14156/(wrad-wrad2))
            else 
               weight=0.
            end if

            do i=1,nlon
               fld(i,j)=fld(i,j)*(1. - weight*0.33)
               !TEST2 fld(i,j)=fld(i,j)
               !TEST  fld(i,j)=weight*6*3600
            end do
         end do
      else
         if (mnproc==1) print *,'No ERA40 fix for variable '//vname
         call xcstop('(era40_fix)')
         stop '(era40_fix)'
      end if


   END subroutine era40_fix
end module m_era40_fix
