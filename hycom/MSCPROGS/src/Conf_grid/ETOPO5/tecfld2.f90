subroutine tecfld(fname,ii,jj,modlat,modlon,depths,field)
   character(len=*), intent(in) :: fname
   integer, intent(in) :: ii,jj
   real, intent(in) :: modlat(ii,jj)
   real, intent(in) :: modlon(ii,jj)
   real, intent(in) :: depths(ii,jj)
   real, intent(in) :: field(ii,jj)

   character*4 tag

   open(44,file=fname//'.dat',status='unknown')
      write(44,*)'TITLE = "',fname,'"'
      write(44,*)'VARIABLES = "i"  "j"  "lat" "lon" "depths" "depths2"'
      write(44,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',ii,', J=',jj,', K=1'


! Positions in Grid index
      write(44,'(30I4)')((i,i=1,ii),j=1,jj)
      write(44,'(30I4)')((j,i=1,ii),j=1,jj)


! Positions in Lon Lat
      write(44,900)((modlat(i,j),i=1,ii),j=1,jj)
      write(44,900)((modlon(i,j),i=1,ii),j=1,jj)
      write(44,900)((depths(i,j),i=1,ii),j=1,jj)
      write(44,900)((field(i,j),i=1,ii),j=1,jj)

   close(44)
900  format(10(1x,e12.5))

end subroutine tecfld

