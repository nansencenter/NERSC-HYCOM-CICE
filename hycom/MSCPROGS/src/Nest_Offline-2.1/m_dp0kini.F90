module m_dp0kini
contains 
subroutine dp0kini(dp0k,kdm,blkfilein)
   use m_parse_blkdat
   implicit none
   real, parameter :: onem=9806.
   real, parameter :: qonem=1./onem

   integer, intent(in)  :: kdm
   real   , intent(out) :: dp0k(kdm)
   character(len=*), intent(in) :: blkfilein

   real :: dp00, dp00x, dp00f, ds00, ds00x, ds00f, dp0kf
   real :: ds0kf, thkmin
   real :: dsm,dsms, dpms, dpm, realvar
   real :: ds0k(kdm), dssk(kdm)
   integer :: nhybrd,nsigma,intvar
   integer :: k
   logical :: isopyc, hybrid


   call parse_blkdat('dp00  ','real',dp00 ,intvar,blkfilein)
   call parse_blkdat('dp00x ','real',dp00x,intvar,blkfilein)
   call parse_blkdat('dp00f ','real',dp00f,intvar,blkfilein)
   call parse_blkdat('ds00  ','real',ds00 ,intvar,blkfilein)
   call parse_blkdat('ds00x ','real',ds00x,intvar,blkfilein)
   call parse_blkdat('ds00f ','real',ds00f,intvar,blkfilein)
   call parse_blkdat('thkmin','real',thkmin,intvar,blkfilein)
   call parse_blkdat('nhybrd','integer',realvar,nhybrd,blkfilein)
   call parse_blkdat('nsigma','integer',realvar,nsigma,blkfilein)


   isopyc=nhybrd==0
   if (isopyc) then
      print *,'Can not handle isopycnic yet'
      stop '(dp0kini)'
   end if

   if (nsigma>1) then
      print *,'Can not handle sigma yet '
      stop '(dp0kini)'
   end if

   hybrid = .not. isopyc
   if (hybrid .and. nsigma.le.1) then
     nsigma=1
     dp00 =dp00
     dp00x=dp00x
     dp00f=dp00f
   endif


   ! Consistency checks
   if (dp00f.lt.1.0) then
     write(6,'(/ a /)')  'error - must have dp00f>=1.0'
     stop '(dp0kini)'
   end if
   if (dp00f.eq.1.0 .and. dp00.ne.dp00x) then
     write(6,'(/ a /)') 'error - must have dp00x==dp00 for dp00f==1.0'
         stop '(dp0kini)'
   endif
   if (dp00.gt.dp00x) then
     write(6,'(/ a /)') 'error - dp00x must be at least dp00'
     stop '(dp0kini)'
   endif
   if (ds00.gt.dp00 .or. ds00x.gt.dp00x .or. ds00f.gt.dp00f) then
     write(6,'(/ a /)')  'error - must have ds00,ds00x,ds00f <= dp00,dp00x,dp00f'
     stop '(dp0kini)'
   endif
   if (ds00.le.0.0) then
     write(6,'(/ a /)') 'error - must have ds00>0.0'
     stop '(dp0kini)'
   endif
   if (ds00f.lt.1.0) then
     write(6,'(/ a /)')  'error - must have ds00f>=1.0'
     stop '(dp0kini)'
    endif
   if (ds00f.eq.1.0 .and. ds00.ne.ds00x) then
      write(6,'(/ a /)')  'error - must have ds00x==ds00 for ds00f==1.0'
      stop '(dp0kini)'
   endif
   if (ds00.gt.ds00x) then
     write(6,'(/ a /)') 'error - ds00x must be at least ds00'
     stop '(dp0kini)'
    endif


! --- logorithmic k-dependence of dp0 (deep z's)
      dp00 =onem*dp00
      dp00x=onem*dp00x
      if     (isopyc) then
        dp0k(1)=thkmin*onem
      else
        dp0k(1)=dp00
      endif
      dpm  = dp0k(1)*qonem
      dpms = dpm
      write(6,*)
      write(6,135) 1,dp0k(1)*qonem,dpm,dpms
 135  format('dp0k(',i2,') =',f7.2,' m',  &
                '    thkns =',f7.2,' m',  &
                '    depth =',f8.2,' m')
!
      dp0kf=1.0
      do k=2,kdm
        dp0kf=dp0kf*dp00f
        if     (k.le.nhybrd) then
          dp0k(k)=min(dp00*dp0kf,dp00x)
        else
          dp0k(k)=0.0
        endif
        dpm  = dp0k(k)*qonem
        dpms = dpms + dpm
        write(6,135) k,dp0k(k)*qonem,dpm,dpms
        !write(6,*) 'geopar: dp0kf  = ',dp0kf
        !write(6,*) 'geopar: dp0k   = ',dp0k(k),k
      enddo
!
! --- logorithmic k-dependence of ds0 (shallow z-s)
      ds00 =onem*ds00
      ds00x=onem*ds00x
      if     (isopyc) then
        ds0k(1)=thkmin*onem
      else
        ds0k(1)=ds00
      endif
      dsm  = ds0k(1)*qonem
      dsms = dsm
      write(6,*)
      write(6,130) 1,ds0k(1)*qonem,dsm,dsms
 130  format('ds0k(',i2,') =',f7.2,' m', &
                '    thkns =',f7.2,' m', &
                '    depth =',f8.2,' m')
!
      ds0kf=1.0
      do k=2,nsigma
        ds0kf=ds0kf*ds00f
        ds0k(k)=min(ds00*ds0kf,ds00x)
        dsm  = ds0k(k)*qonem
        dsms = dsms + dsm
        write(6,130) k,ds0k(k)*qonem,dsm,dsms
        !write(6,*) 'geopar: ds0kf  = ',ds0kf
        !write(6,*) 'geopar: ds0k   = ',ds0k(k),k
      enddo
      write(6,*)
!
! --- sigma-depth scale factors
      do k=1,nsigma
        dssk(k)=ds0k(k)/dsms  ! onem * fraction of depths in sigma layer k
      enddo
      do k= nsigma+1,kdm
        ds0k(k)=dp0k(k)
        dssk(k)=0.0           ! these layers are zero in sigma mode
      enddo

end subroutine dp0kini
end module m_dp0kini

