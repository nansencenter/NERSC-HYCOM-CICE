program gpdens
! code for computing a density function for the velocity time series
use mod_modtype
use m_density
use m_densitylevels
use m_accplot
use m_mkestat
use m_correlation
use m_exceedence
use m_rotate
use m_time_filter
use m_jot
use m_merge

implicit none

integer, parameter :: maxnr=1000000
type(data_new)  cm(maxnr),gp(maxnr)
type(data_new)  cm2(maxnr),gp2(maxnr)
type(data_new), allocatable :: cmfilt(:),gpfilt(:)
type(data_new), allocatable :: cmtide(:),gptide(:)

integer i,j,m,nr
character(len=100) fnameCM,fnameGP,fname
character(len=140) string
logical ex
integer ncm,ngp,igp,icm,i1,k,nrm
integer deep,istation
real thetacm,thetagp
real cmlon,cmlat,gplon,gplat
integer cmdeep,gpdeep,totsta,stno
real torad,todeg
character(len=6)dum
character(len=50) stname
integer, parameter :: nlev=15   ! <100
real maxA,da



fnameCM=' '; fnameGP=' ';fname=' '
open(10,file='infile.GPCM')
   read(10,'(a)')fnameGP
   read(10,'(a)')fnameCM
   read(10,'(a)')fname
close(10)
print '(4a)','filenames are:',trim(fnameCM), trim(fnameGP), trim(fname)
read(fname(12:15),'(i4)')deep
read(fname(8:10),'(i3)')istation

print *,'Station= ',istation,'  Deep= ',deep

open(10,file='cminfile.in')
   read(10,'(a, i5)')dum,totsta
   do i=1,totsta
      read(10,'(i3, a50)',err=999) stno, stname
      if (istation == stno) exit
   enddo
999 close(10)



! info table
open(10,file='infofile.txt')
do i=1,1000
   read(10,'(a)',end=998)string
   if (string(1:3) == fname(8:10)) exit
enddo
998 close(10)
print '(a)',trim(string)

read(string(6:13) ,'(f8.3)')gplon
read(string(16:23),'(f8.3)')gplat
read(string(26:33),'(f8.3)')cmlon
read(string(36:43),'(f8.3)')cmlat

!read(string(117:120),'(i4)')cmdeep
!read(string(123:126),'(i4)')gpdeep

torad=pi/180.0
todeg=1.0/torad

open(10,file=trim(fnameCM))
   ncm=0
   do m=1,maxnr
      read(10,'(f9.4,2f9.3,4f8.2)',end=101)cm(m)
      cm(m)%dir=atan2(cm(m)%v,cm(m)%u)*todeg
      ncm=ncm+1
   enddo
   print *,'maxnr too small for cm data'
101 close(10)
   if (ncm < 100) then
      call system('touch no_merged_points')
      stop 'empty CM file'
   endif
print '(2a)','files read ok:',trim(fnameCM)



open(10,file=trim(fnameGP))
   ngp=0
   do m=1,maxnr
      read(10,'(f9.4,2f9.3,4f8.2)',end=102)gp(m)  
      gp(m)%dir=atan2(gp(m)%v,gp(m)%u)*todeg
      ngp=ngp+1
   enddo
   print *,'maxnr too small for gp data'
102 close(10)
   if (ngp < 100) then
      call system('touch no_merged_points')
      stop 'empty GP file'
   endif

print '(2a)','files read ok:',trim(fnameGP)

! merging the two data sets
   call merge(gp,cm,gp2,cm2,ngp,ncm,maxnr,nrm)

! generate exceedence diagnostics
   call exceedence(cm2,gp2,nrm,deep)

! generate accumulation plots
   call accplot(cm2,gp2,nrm,deep)

! compute density functions
   maxA=0.0
   call density(cm2,nrm,'CM',trim(fnameCM),thetacm,fname(8:10),fname(12:15),maxA,'U')
   call density(gp2,nrm,'GP',trim(fnameGP),thetagp,fname(8:10),fname(12:15),maxA,'U')
   call densitylevels(fname,maxA,nlev,'U')


! compute Joint occurence tables
   call jot(cm2,nrm,'CM',fname(12:15),stname,cmlon,cmlat)
   call jot(gp2,nrm,'GP',fname(12:15),stname,gplon,gplat)

! compute low-pass filtered timeseries
   allocate(cmfilt(nrm),gpfilt(nrm))
   allocate(cmtide(nrm),gptide(nrm))
   call time_filter(cm2,gp2,cmfilt,gpfilt,cmtide,gptide,nrm,24.0,'AVE',deep)

! compute density functions for filtered velocities
   maxA=0.0
   call density(cmfilt,nrm,'CM',trim(fnameCM)//'_F',thetacm,fname(8:10),fname(12:15),maxA,'F')
   call density(gpfilt,nrm,'GP',trim(fnameGP)//'_F',thetagp,fname(8:10),fname(12:15),maxA,'F')
   call densitylevels(fname,maxA,nlev,'F')

! compute density functions for tidal velocities
   maxA=0.0
   call density(cmtide,nrm,'CM',trim(fnameCM)//'_T',thetacm,fname(8:10),fname(12:15),maxA,'T')
   call density(gptide,nrm,'GP',trim(fnameGP)//'_T',thetagp,fname(8:10),fname(12:15),maxA,'T')
   call densitylevels(fname,maxA,nlev,'T')

! compute MKE statistics
  call mkestat(cm2,gp2,cmfilt,gpfilt,cmtide,gptide,nrm,deep,istation)

! complex correlation
   call correlation(cm2,gp2,cmfilt,gpfilt,cmtide,gptide,nrm,deep)


! print low pass filtered timeseries rotated to principal direction
   call rotate(cmfilt,gpfilt,nrm,thetacm,thetagp,deep)

! call spectral(cm2,gp2)

open(10,file='statloc.dat')
   write(10,'(2f10.2)')cmlon,cmlat
close(10)

open(10,file='statname.dat')
   write(10,'(a)')trim(stname)
close(10)

open(10,file='statinfo.dat')
   write(10,'(i6,"&",3(f10.2,"&"),i6,a,2(I4,a,i2,a))')deep,cm2(1)%day,cm2(nrm)%day,float(nrm-1)/24.0,nrm,&
                                    " & (",int(cmlon),'''',abs(nint((cmlon-int(cmlon))*600.0/10.0)),', ',&
                                           int(cmlat),'''',abs(nint((cmlat-int(cmlat))*600.0/10.0)),') \\\\'
close(10)

open(10,file='infotablecm.txt')
   write(10,'(a)')'\\begin{tabular}[l]'
   write(10,'(3a)')'Station: ',trim(stname),' \\\\'
   write(10,'(a,2(i4,a,i2,a))')'Location: lon:lat $(',int(cmlon),'''',abs(nint((cmlon-int(cmlon))*600.0/10.0)),', ',&
                                                    int(cmlat),'''',abs(nint((cmlat-int(cmlat))*600.0/10.0)),')$ \\\\'
!   write(10,'(a,i5,a)')'Total model water depth: ',cmdeep,' \\'
   write(10,'(a,f8.2,a)')'Start date (day rel. 01.01.95): ',cm2(1)%day,'\\\\'
   write(10,'(a,f8.2,a)')'End date (day rel. 01.01.95): ',cm2(nrm)%day,'\\\\'
   write(10,'(a,i5,a)')'Mooring/model depth: ',deep,'\\\\'
   write(10,'(a,f8.2,a)')'Total averaging time (days): ',float(nrm-1)/24.0,'\\\\'
   write(10,'(a,i8,a)')'Number of merged points: ',nrm,'\\\\'
   write(10,'(a)')'Printed \\today\\\\'
   write(10,'(a)')'\\end{tabular}'
close(10)

open(10,file='infotablegp.txt')
   write(10,'(a)')'\\begin{tabular}[l]'
   write(10,'(3a)')'Station: ',trim(stname),' \\\\'
   write(10,'(a,2(i4,a,i2,a))')'Location: lon:lat $(',int(gplon),'''',abs(nint((gplon-int(gplon))*600.0/10.0)),', ',&
                                                    int(gplat),'''',abs(nint((gplat-int(gplat))*600.0/10.0)),')$ \\\\'
!   write(10,'(a,i5,a)')'Total model water depth: ',gpdeep,' \\\\'
   write(10,'(a,f8.2,a)')'Start date (day rel. 01.01.95): ',gp2(1)%day,'\\\\'
   write(10,'(a,f8.2,a)')'End date (day rel. 01.01.95): ',gp2(nrm)%day,'\\\\'
   write(10,'(a,i5,a)')'Mooring/model depth: ',deep,'\\\\'
   write(10,'(a,f8.2,a)')'Total averaging time (days): ',float(nrm-1)/24.0,'\\\\'
   write(10,'(a,i8,a)')'Number of merged points: ',nrm,'\\\\'
   write(10,'(a)')'Printed \\today\\\\'
   write(10,'(a)')'\\end{tabular}'
close(10)

end program
