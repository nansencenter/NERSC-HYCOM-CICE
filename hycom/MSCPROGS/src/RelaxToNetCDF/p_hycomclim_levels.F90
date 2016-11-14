! KAL - program reads a hycom-type climatology input file (.d) and 
! KAL - retrieves the number of vertical levels in them. The levels 
! KAL - are dumped to the file hclim_depthlevels
program hycomclim_levels
implicit none

character(len=40) :: ctitle
integer           :: IWI,JWI,KWI
real*4            :: XFIN,YFIN,DXIN,DYIN
real*4,allocatable:: ZLEV(:)

integer :: k


open(41,file='t_m01.d',form='unformatted',status='old')
read(41) ctitle
read(41)  IWI,JWI,KWI,XFIN,YFIN,DXIN,DYIN
if (kwi>0 .and. kwi < 200) then
   allocate(ZLEV(kwi))
   rewind(41)
   read(41) ctitle
   read(41)  IWI,JWI,KWI,XFIN,YFIN,DXIN,DYIN,ZLEV
   close(41)
else
   print *,'depthlevels > 200, must be an error...'
   close(41)
   stop
end if

print *,ctitle
print *,IWI, JWI, KWI, XFIN, YFIN, DXIN, DYIN
print *,ZLEV

open(10,file='hclim_depthlevels',status='replace')
do k=1,kwi
   write(10,*) ZLEV(k)
end do
close(10)
end program
