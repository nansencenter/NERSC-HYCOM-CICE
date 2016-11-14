module mod_pakk
contains


subroutine prepak26(rungen,fname,n2d,n3d,varlist,maxfld,idat,rt,sigma,nx,ny,nz)
!  Routine to do preprocessing of data packing
   use mod_year_info
implicit none
   integer, parameter ::maxhdr=1000
   integer :: iw,it,nn,in,n,iy,icc,inc,loop
   integer :: n2d, n3d, maxfld,nz,ny,nx


   character(len=*), intent(in) :: rungen
   character(len=*), intent(in) :: fname
   character(len=*), intent(in) :: varlist(maxfld)
   real, dimension(nz) :: sigma
   type(year_info), intent(in)  :: rt


   integer lines

   character*80 whole(maxhdr)

   integer idat(3)
   INTEGER count

!  Initialise Variables for file header  ???
   whole(:)=" "

!  Create file Header
      !nsday=(86400./baclin)+0.0001

      WRITE(whole(2),'(a,a3)')'Run: ',rungen 
!      WRITE(whole(3),'(a,a4,a,a2,a,a3,a,a2,a,a4)')'y=',rt%cyy,' m=',rt%cmm,' d=',rt%cdd,' h=',rt%chh,' s=',rt%css
      WRITE(whole(3),'(a)')fname
      WRITE(whole(4),103) n2d
      DO count=5,n2d+4
         WRITE (whole(count),'(a5)') varlist(count-4)
      enddo

      WRITE(whole(count),105) n3d
      DO loop=count+1,count+n3d
         WRITE(whole(loop),'(a5)') varlist(loop-count+n2d)
      enddo
      !WRITE(whole(loop),106) IDM,JDM,KDM
      WRITE(whole(loop),106) nx,ny,nz
      WRITE(whole(loop+1),107)

      inc=loop+2
      icc=1
      DO n=1,nz,5
         in=n+4
         IF (in.GT.nz) in=nz

         DO nn=n,in
            WRITE(whole(inc)(icc:icc+10),108) nn,sigma(nn)
            !WRITE(whole(inc)(icc:icc+10),108) nn,1.
            icc=icc+11
         enddo

         inc=inc+1
         icc=1
      enddo

      !WRITE(whole(inc),109) gridn,xpivn
      WRITE(whole(inc),109) 0.0,0.0
      WRITE(whole(inc+1),113)
      WRITE(whole(inc+2),114)
      WRITE(whole(inc+3),115)
      WRITE(whole(inc+4),116)
      DO iw=maxhdr,2,-1
         IF (whole(iw).EQ." ") lines=iw
      enddo

      WRITE(whole(1),'(a3,i4)')'1.1',lines

 103  FORMAT("2D fields: ",I2)
 105  FORMAT("3D fields: ",I2)
 106  FORMAT("Domain: I=",I4," J=",I4," K=",I4)
 107  FORMAT("Layer Sigma-Theta:")
 108  FORMAT(I2,"=",F6.3)
 109  FORMAT("Resolution (gridn)= ",F6.2," Equator is at xpivn=",F6.2) 
 113  FORMAT("Matrix written K(J(I)),all 2D and Layer 1 fields first")
 114  FORMAT("Field separator is '++'.")
 115  FORMAT("Field header format 'Field Name,Layer,x,y,a,b,c,ASCII"," string length'")
 116  FORMAT("in fortran format('++ ',A5,3I2,3I4,I5) and data ","format (40(A2)).")



   open(10,file=trim(fname)//'.hdr')
      do it=1,lines
         iy=index(whole(it),'        ')
         write(10,'(a)') whole(it)(1:iy)
      enddo
   close(10)

end subroutine prepak26




subroutine pakmsk(mask,array,work,ii,jj,util,length)
   implicit  none
! send 'array' to subr. pakk after deleting array elements outside 'mask'
   integer, intent(in)  :: ii,jj, length ! local versions ....
   real   , intent(in)  :: array(ii,jj)
   real   , intent(out) :: work(ii,jj)
   character*2, intent(out) :: util(ii*jj+14)
   integer mask(ii,jj)
   integer :: i,j


   do j=1,jj
   do i=1,ii
      work(i,j)=0.
     if (mask(i,j).ne.0) work(i,j)=array(i,j)
   enddo
   enddo
   call pakk(work,ii,jj,util,length)
end subroutine



subroutine mkfldh(fieldname,ix,iy,layer,ia,ib,ic,length,fieldhead)
   implicit none
   integer, intent(in) :: ix,iy,layer,ia,ib,ic,length
   character(len=*), intent(in)  ::  fieldname
   character(len=*), intent(out) ::  fieldhead
   write(fieldhead,100) fieldname,layer,ix,iy,ia,ib,ic,length 
!KAL --  I6 is too small for huge grids. This approach is backwards compatible
!100 format("++ ",A5,3I2,3I4,I6)
100 format("++ ",A5,3I2,3I4,I10)
end subroutine mkfldh




      subroutine pakk(array,ii,jj,compac,length)
      implicit none
      integer, intent(in) :: ii,jj
      integer :: length,idim
      real    :: base,scal
      integer :: nbits,i,j,l,i1,i2,numb,lngth
!
! --- converts the contents of -array- into an ascii character string which
! --- is stored in character*2 array -compac-. compac(1)...compac(7) contain
! --- the base value, i.e., the minimum value encountered in -array-.
! --- compac(8)...compac(14) contain a scale factor by which the individual
! --- 6-bit integers encoded as ascii character pairs in compac(8),...
! --- compac(length) must be multiplied before the base value is added
! --- during an unpakking operation. base value and scale factor are encoded
! --- in e14.7 format.
!
! --- the printable ascii characters used to encode the integers include
! --- the numbers 0...9, upper- and lower-case letters a...z, a...z, plus
! --- two additional characters '.' and '/' (total of 64).
!
! --- a packing operation fills (ii*jj+14) array elements in -compac- which
! --- must be dimensioned accordingly in the calling program. the total
! --- number of occupied array elements is returned in -length-. in calls to
! --- unpack, -length- is treated as input variable.
!
      real array(ii,jj)
      character*2 char,compac(ii*jj+14),comp2(14)
      character*14 comp14(2)
      equivalence (comp2,comp14)
      data nbits/12/
      base=1.e22
      do 1 i=1,ii
      do 1 j=1,jj
 1    base=min(base,array(i,j))
      scal=0.
      do 2 i=1,ii
      do 2 j=1,jj
 2    scal=max(scal,array(i,j)-base)
      scal=scal/float(2**nbits-1)
      i1=0
      i2=0
      length=14
      do 3 i=1,ii
      do 3 j=1,jj
      if (scal.eq.0.) go to 7
      numb=(array(i,j)-base)/scal+.5
      i1=numb/64
      i2=numb-64*i1
!
! --- map 6-bit numbers onto character set consisting of numbers
! --- 0...9, letters a...z, a...z, and the two characters '.' and '/'.
! --- (if mapping into the character range 32...95 -- which includes the
! --- characters !"#$%&'()*+,-./:;<=>?@[\]^_  -- is deemed safe, delete
! --- the next 6 lines.)
      if (i1.gt.37) i1=i1+6
      if (i1.gt.11) i1=i1+7
      i1=i1+14
      if (i2.gt.37) i2=i2+6
      if (i2.gt.11) i2=i2+7
      i2=i2+14
!
 7    length=length+1
      compac(length)(1:1)=char(i1+32)
      compac(length)(2:2)=char(i2+32)
 100  format (a2)
 3    continue
      write (comp14(1),101) base
      write (comp14(2),101) scal
 101  format (1pe14.7)
      do 8 l=1,14
 8    compac(l)=comp2(l)
!
      return
!
!
      entry unpakk(array,idim,ii,jj,compac,length)
!
      do 9 l=1,14
 9    comp2(l)=compac(l)
      read (comp14(1),101) base
      read (comp14(2),101) scal
      lngth=14
      do 4 i=1,ii
      do 4 j=1,jj
      lngth=lngth+1
      i1=ichar(compac(lngth)(1:1))
      i2=ichar(compac(lngth)(2:2))
!
! --- 6-bit numbers are mapped onto character set consisting of numbers
! --- 0...9, letters a...z, a...z, and the two characters '.' and '/'.
! --- (if mapped into character range 32...95, delete next 6 lines)
      if (i1.gt.96) i1=i1-6
      if (i1.gt.64) i1=i1-7
      i1=i1-14
      if (i2.gt.96) i2=i2-6
      if (i2.gt.64) i2=i2-7
      i2=i2-14
!
 4    array(i,j)=scal*float(64*(i1-32)+(i2-32))+base
      if (lngth.ne.length) stop 'unpack'
      return
      end subroutine pakk

end module mod_pakk
