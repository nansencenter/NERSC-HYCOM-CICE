! skeleton of subroutine - illustrates how to read old latlon.uf file format
! This routine has not yet been compiled...
module m_read_latlonuf
implicit none



subroutine read_latlonuf(qlon,qlat,ulon,ulat,vlon,vlat,plon,plat,idm,jdm,filename,ierr)
   implicit none
   character(len=*), intent(in) :: filename
   real, dimension(idm,jdm), intent(out) ::  &
      qlon, qlat, ulon, ulat, vlon, vlat, plon, plat
   integer,intent(in) :: idm,jdm
   integer,intent(out) :: ierr

   logical :: exuf
   integer :: itst,jtst, ios
   real*8, dimension(0:idm+1,0:jdm+1), intent(out) ::  &
        qlon2, qlat2, ulon2, ulat2, vlon2, vlat2, plon2, plat2



   inquire(file=trim(filename),exist=exuf)
   if (exuf) then
     write(*,*)'Load grid positions from file: latlon.uf'
     open(10,file='./Data/latlon.uf',form='unformatted')
     read(10)itst,jtst
     if (idm/=itst .or. jdm/=jtst) then
        print *,'(read_latlonuf: Grid size mismatch)'
        ierr=1
        return
     end if
     rewind(10)
     read(10,iostat=ios)idm,jdm,qlat2,qlon2,plat2,plon2,ulat2,ulon2,vlat2,vlon2
     close(10)
   else
     print *,'(read_latlonuf: cant find '//trim(filename)//')'
     ierr=2
     return
   end if

   ierr=ios ! NB
   return
end subroutine
end module


