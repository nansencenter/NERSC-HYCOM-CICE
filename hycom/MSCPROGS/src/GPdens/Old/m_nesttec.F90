module m_nesttec
contains
subroutine nesttec(fname,ii,jj,depths,modlon,modlat,factor,sponge,relax,relaxu,relaxv,&
                       dampp,spongu,spongv)
   character(len=*), intent(in) :: fname
   integer, intent(in) :: ii,jj
   real, intent(in) :: depths(ii,jj)
   real, intent(in) :: factor(ii,jj)
   real, intent(in) :: sponge(ii,jj)
   real, intent(in) :: relax(ii,jj)
   real, intent(in) :: relaxu(ii,jj)
   real, intent(in) :: relaxv(ii,jj)
   real, intent(in) :: dampp(ii,jj)
   real, intent(in) :: spongu(ii,jj)
   real, intent(in) :: spongv(ii,jj)
   real, intent(in) :: modlon(ii,jj)
   real, intent(in) :: modlat(ii,jj)

   character*4 tag

   open(10,file='tec'//fname//'.dat',status='unknown')
      write(10,*)'TITLE = "',fname,'"'
      write(10,*)'VARIABLES = "i"  "j" "lon" "lat" "depths" "factor" "sponge" "relax" &
            &"relaxu" "relaxv" "dampp" "spongu" "spongv"'
      write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',ii,', J=',jj,', K=1'


! Positions in Grid index
      write(10,'(30I4)')((i,i=1,ii),j=1,jj)
      write(10,'(30I4)')((j,i=1,ii),j=1,jj)


! Positions in Lon Lat
      write(10,900)((modlon(i,j),i=1,ii),j=1,jj)
      write(10,900)((modlat(i,j),i=1,ii),j=1,jj)
      write(10,900)((depths(i,j),i=1,ii),j=1,jj)
      write(10,900)((factor(i,j),i=1,ii),j=1,jj)
      write(10,900)((sponge(i,j),i=1,ii),j=1,jj)
      write(10,900)((relax(i,j),i=1,ii),j=1,jj)
      write(10,900)((relaxu(i,j),i=1,ii),j=1,jj)
      write(10,900)((relaxv(i,j),i=1,ii),j=1,jj)
      write(10,900)((dampp(i,j),i=1,ii),j=1,jj)
      write(10,900)((spongu(i,j),i=1,ii),j=1,jj)
      write(10,900)((spongv(i,j),i=1,ii),j=1,jj)

   close(10)
900  format(10(1x,e12.5))

end subroutine nesttec
end module m_nesttec
