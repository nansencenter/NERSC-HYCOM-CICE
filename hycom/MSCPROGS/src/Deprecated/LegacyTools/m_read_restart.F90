module m_read_restart

use mod_dimensions
type states
   real u(itdm,jtdm,kk)
   real v(itdm,jtdm,kk)
   real d(itdm,jtdm,kk)
   real t(itdm,jtdm,kk)
   real s(itdm,jtdm,kk)
   real ub(itdm,jtdm)
   real vb(itdm,jtdm)
   real pb(itdm,jtdm)
end type states


type (states), save :: mem ! Keeps member read from restart file

contains

function read_restart_mem(filename,imem)
   implicit none
   logical read_restart_mem
   character(len=*),intent(in) :: filename
   integer,         intent(in) :: imem

   character(len=24) filename
   character(len=9), parameter ::  cident='HYCOM_1.1'
   character(len=9), parameter :: cidentA='HYCOM_1.0'
   character(len=9), parameter :: cident0='MIC2HYC00'
   character(len=9), parameter :: cident1='CURVIINT1'
   character(len=9), parameter :: cident2='HYCOM_2.0'
   character(len=9) rident
   logical ex
   integer i,j,k,l


   ! Global arrays
   real, dimension(itdm,jtdm) :: tmp,pmix,pbot,psikk,thkk,dpmixl


   if (mnproc==1) print *,'reading:',filename
   inquire(file=filename,exist=ex)
   if (.not.ex) then
      stop 'restart file for MICOM does not exist'
   endif


   inquire(iolength=j)rident,mem,tpbot,tpsikk,tthkk,tdpmixl
   if (mnproc==1) print *,'iolength=',j
   open(10,file=filename,status='old',form='unformatted',access='direct',recl=j)
   read(10,rec=imem)rident
   if (rident == cident) then
      close(10)
      inquire(iolength=j)rident,mem
      open(10,file=filename,status='old',form='unformatted',access='direct',recl=j)
      if (mnproc==1) print *,'reading HYCOM_1.1 restart file.  record=',imem
      read(10,rec=imem)rident,mem

   elseif (rident == cidentA) then
      if (mnproc==1) print *,'reading HYCOM_1.0 restart file.  record=',imem
      read(10,rec=imem)rident,mem,tmp,tmp,tmp,tmp !,entrn

   elseif (rident == cident0) then
      if (mnproc==1) print *,'reading mic2hyc generated restart file.  record=',imem
      read(10,rec=imem)rident,mem,tmp,tmp,tmp,tmp

   elseif (rident == cident1) then
      if (mnproc==1) print *,'reading curviint generated restart file.  record=',imem
      read(10,rec=imem)rident,mem

   elseif (rident == cident2) then
      if (mnproc==1) print *,'reading HYCOM_2.0 restart file.  record=',imem
      read(10,rec=imem)rident,mem,tpbot,tpsikk,tthkk,tdpmixl

   else
      if (mnproc==1) print *,'ERROR: wrong ident of restart file: ',rident
      if (mnproc==1) print *,'it should be: ',cident
      stop
   endif

   ! Conversion from CGS to SI units when reading old restart files
   if ((rident == cident0).or.(rident == cidentA)) then
      print *,'read_restart: Converting to SI units'
      mem%u=0.01*mem%u
      mem%v=0.01*mem%v
      mem%d=0.1*mem%d
      mem%ub=0.01*mem%ub
      mem%vb=0.01*mem%vb
      mem%pb=0.1*mem%pb
   endif

   if (mnproc==1) print *,'read record ',imem,' from ',filename,' with ident=',rident
   close(10)
   if (mnproc==1) print *,'done'
   read_restart_mem=.true.
end function read_restart_mem




!function read_restart_ice(rungen,rt,imem1,imem2,iceproc)
!   use mod_dimensions
!   use mod_hycomfuncs
!   use mod_year_info
!   use mod_common_ice
!   use mod_xc
!   use mod_ensemble_io
!   use mod_distribute
!   implicit none
!   logical read_restart_ice
!   type(year_info), intent(in)  :: rt
!   character(len=3), intent(in) :: rungen
!   real, dimension(itdm,jtdm) :: tficem,thicem,thsnwm,tticem,ttsrfm,ticeU,ticeV
!   integer, intent(in) :: imem1,imem2,iceproc
!   character(len=27) icename
!   logical ex
!   integer i,j
!   real x(10)
!
!
!   icename=fileice(rungen,rt)
!   print *,'reading:',trim(icename)
!   inquire(file=icename,exist=ex)
!   if (ex) then
!      inquire(iolength=j)tficem,thicem,thsnwm,tticem,ttsrfm  !,ticeU,ticeV
!      open(10,file=icename,status='old',form='unformatted',access='direct',recl=j)
!         ticeU=0.0; ticeV=0.0
!         do i=imem1,imem2
!
!            read(10,rec=i)tficem,thicem,thsnwm,tticem,ttsrfm !,ticeU,ticeV
!
!            ! Distribute ice variables to tiles
!            call distribute_ice(tficem,thicem,thsnwm,tticem,ttsrfm,ticeU,ticeV)
!
!            ! Save in ensemble file
!            if (saveice) call write_ensemble('ICE',i)
!
!         enddo
!      close(10)
!      if (mnproc==1) print *,'done'
!   else 
!      if (mnproc==1) print *,'Restart file for ICE does not exist.'
!      if (mnproc==1) print *,'I continue with fields from iniice'
!      if (saveice) then
!         do i=imem1,imem2
!            call write_ensemble('ICE',i)
!         enddo
!      endif
!   endif
!   read_restart_ice=.true.
!
!end function read_restart_ice
!
!
!
!
!
!function read_restart_eco(rungen,rt,bio,nbc,hbc,imem1,imem2)
!use mod_dimensions
!use mod_ecodim
!use mod_hycomfuncs
!use mod_year_info
!implicit none
!logical read_restart_eco
!type(year_info), intent(in)  :: rt
!character(len=3), intent(in) :: rungen
!real, intent(out) :: bio(nbio,kbc,idm,jdm)
!real, intent(out) :: nbc(kdm,idm,jdm)
!real, intent(out) :: hbc(kbc,idm,jdm)
!integer, intent(in) :: imem1,imem2
!
!character(len=27) econame
!logical ex
!integer i,j
!
!econame=fileeco(rungen,rt)
!print *,'reading:',trim(econame)
!inquire(file=trim(econame),exist=ex)
!if (ex) then
!   inquire(iolength=j)bio,nbc,hbc
!   open(10,file=econame,status='old',form='unformatted',access='direct',recl=j)
!      do i=imem1,imem2
!         read(10,rec=i)bio,nbc,hbc
!         write(53,rec=i)bio,nbc,hbc
!      enddo
!   close(10)
!   print *,'done'
!else 
!   print *,'Restart file for ECO does not exist.'
!   print *,'Continues with fields from biochm_ini'
!#ifdef ECOSYS
!   brestart=.false.
!   vert_flx=0.0  
!   call biochm_ini(brestart,rt)
!   call bclin11(brestart)
!   call biochm_ini2(brestart,rt)
!   do i=imem1,imem2
!      write(53,rec=i)bio,nbc,hbc
!   enddo
!#endif
!endif
!read_restart_eco=.true.
!end function read_restart_eco
!
!function read_restart_ave(rungen,rt,imem1,imem2)
!use mod_dimensions
!use mod_hycomfuncs
!use mod_year_info
!use mod_average
!implicit none
!logical read_restart_ave
!type(year_info), intent(in)  :: rt
!character(len=3), intent(in) :: rungen
!integer, intent(in) :: imem1,imem2
!character(len=27) avename
!logical ex
!integer i,j
!
!print *,'reading:', fileaveres(rungen,rt)
!inquire(file=fileaveres(rungen,rt),exist=ex)
!if (.not.ex) stop 'restart file for AVERAGE does not exist'
!open(10,file=fileaveres(rungen,rt),status='unknown',form='unformatted')
!   read(10)ave_week
!close(10)
!print *,'counter= ',ave_week%counter
!print *,'done'
!read_restart_ave=.true.
!end function read_restart_ave


end module m_read_restart

