
c     Simple program which reads the identification of fields from
c     a sequential, unformatted "DNMI felt"-file.

      program se_sekvensiell

      implicit none

      integer*2 ident(20),felt(5)
      integer n
      character*50 navn

c     For reading command line arguments:
      integer iargc,narg

      narg = iargc()
      if (narg .lt. 1) then
         write(*,*)'File name: '
         read(5,'(a50)')navn
      else
         call getarg(1,navn)
      end if

      open(21, file = navn, status = 'old', access = 'sequential',
     >     form = 'unformatted')

      do while (.true.)

         read(21,end=999)ident
         write(*,100)ident
         if (ident(10) .lt. 0 .or. ident(11) .lt. 0)
     >      read(21,end=999)n
         read(21,end=999)felt

         if (ident(10) .lt. 0 .or. ident(11) .lt. 0)
     >      write(*,101)'uncompressed size',
     >                abs(int(ident(10))*int(ident(11))),
     >                'compressed size',n,'algorithm',felt(5)

      end do
 100  format(i4,i5,i3,i5,i5,i4,i5,i2,i5,2i5,3i5,2i6,2i5,i3,i3)
 101  format(a17,1x,i7,3x,a15,1x,i7,3x,a9,1x,i1)
 999  close(21)
      end
