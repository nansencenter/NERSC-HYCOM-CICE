module m_jot
contains
subroutine jot(new,nrm,cdata,chdeep,stname,lon,lat)
use mod_data_new
implicit none
integer, intent(in) :: nrm
type(data_new), intent(in) :: new(nrm)
character(len=2), intent(in) :: cdata
character(len=4), intent(in) :: chdeep
real, intent(in) :: lon,lat
character(len=50)stname

integer, parameter :: nx=20
integer, parameter :: ny=16
real, parameter :: da=22.5
real, parameter :: ds=10.0
integer ia(nx,ny)
real     a(nx,ny)
real    as(nx,ny)

real p(2,2),lambda(2),zz(2,2),w(10),t1,t2
integer info
integer i,j,m
real tmp
real time,u,v,speed,dir
real scale,scale1
real varu,varv,covuv,aveu,avev
real x(2,2)
real torad,todeg
character(len=5) cia
character(len=4) ca
integer nri

torad=pi/180.0
todeg=1.0/torad


print '(3a)','JOT for ',cdata,' data'

   ia=0
   do m=1,nrm
      speed=new(m)%speed
      dir=new(m)%dir
      if (dir < 0.0) dir=dir+360.0

      i=int(speed/ds)+1
      j=int(dir/da)+1
      i=min(max(1,i),nx)
      j=min(max(1,j),ny)
      ia(i,j)=ia(i,j)+1
   enddo
101 close(10)

   A=100.0*float(ia)/float(nrm)

   do i=nx,1,-1
      if (sum(ia(i,1:ny)) > 0) then
         nri=i+1
         exit
      endif
   enddo

   cia=' '
   ca=' '
   open(10,file='jot'//cdata//'.tex')
   write(10,'(a)')'\begin{flushleft}'
   write(10,'(a)')'\begin{tabular}{|r||*{16}{p{0.040\textwidth}|}l|}'
   write(10,'(a)')'\hline'
   write(10,'(4a)',advance='no')'\multicolumn{18}{|l|}{\bf ',cdata,'---',trim(stname)
   write(10,'(3a)',advance='no')';   Mooring depth: ',chdeep,' (m)'
   write(10,'(a)',advance='no')';   Location:'
   write(10,'(a,2(i4,a,i2,a))',advance='no')'$(',int(lon),'''',abs(nint((lon-int(lon))*600.0/10.0)),', ',&
                                                 int(lat),'''',abs(nint((lat-int(lat))*600.0/10.0)),')$'
   write(10,'(a,i8)',advance='no')';   Number of merged points: ',nrm
   write(10,'(a)')'}\\'
   write(10,'(a)')'\hline'
   write(10,'(a)',advance='no')'   &'
   do j=1,ny
      write(10,'(f5.1,a)',advance='no')float(j)*da,'&'
   enddo
   write(10,'(a)')'\\'
   write(10,'(a)')'\hline'
   write(10,'(a)')'\hline'
   do i=1,min(nri,nx)
      write(10,'(f5.1,a)',advance='no')float(i)*ds,'&'
      do j=1,ny
         if (ia(i,j) == 0) then
            cia=' '
            ca=' '
         else
            write(cia,'(i5)')ia(i,j)
            write(ca,'(f4.1)')a(i,j)
         endif
         write(10,'(a,a5,a,a4,a)',advance='no')'$',cia,':',ca,'$&' 
      enddo
      write(10,'(a,i6,a,f5.1,a)')'$',sum(ia(i,1:ny)),':',sum(a(i,1:ny)),'$\\'
      write(10,'(a)')'\hline'
   enddo
   write(10,'(a)')'\hline'
   write(10,'(a)',advance='no')'   &'
   do j=1,ny
      write(10,'(i5,a,f4.1,a)',advance='no')sum(ia(1:nx,j)),':',sum(a(1:nx,j)),'&' 
   enddo
   write(10,'(i6,a,f5.1,a)')sum(ia(1:nx,1:ny)),':',sum(a(1:nx,1:ny)),'\\'
   write(10,'(a)')'\hline'
   write(10,'(a)')'\end{tabular}'
   write(10,'(a)')'\end{flushleft}'
   close(10)

end subroutine jot
end module m_jot
