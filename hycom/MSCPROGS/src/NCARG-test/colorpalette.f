CKAL -- Stolen without shame from hycom ALL routines
c
C     subroutine pastel(nbase,ibase)
C     integer nbase,ibase(2)
C
C --- define 2 pastel colors
C
C     parameter (ncol=2)
C     real rgb(3,ncol)
C     data rgb/
C    .         0.8,1.0,0.9,		!  pale blue
C    .         0.8,0.7,0.5 		!  beige
C    .        /
Ccc     .         0.5,0.8,0.7,		!  pale blue
Ccc     .         1.0,0.9,0.7,		!  beige
Ccc     .         0.9,0.8,0.6,		!  beige
C
C     if     (ibase(1).ne.0) then
C       write(*,101) 'existing pastel color table:',' entries',
C    .    ibase(1)+1,' --',ibase(2)
C       return
C     endif
C     ibase(1) = nbase
C     ibase(2) = nbase+ncol
C     nbase    = ibase(2)
C     write (*,101) 'defining pastel color table:',' entries',
C    .   ibase(1)+1,' --',ibase(2)
C     do i=1,ncol
Ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
C       call gscr(1,i+ibase(1),rgb(1,i),rgb(2,i),rgb(3,i))
C     enddo
C     return
C100  format (' color index',i4,'  <==> rgb combination',3f7.3)
C101  format(a32,a,i5,a,i5)
C     end
c
c
      subroutine sstpal(IWKID,nbase,ibase)
      implicit none
      integer IWKID
      integer nbase ! Number of colors
      integer ibase ! Offset in color table 
c
c --- initialize sst color table
c
      integer, parameter :: ncol=64
      integer rgb(3,1:ncol)
      data rgb   /
     . 24, 1,30,  18, 1,30,  14, 1,30,  10, 1,30,   1, 1,30,   1, 1,28,
     .  1, 1,26,   1, 1,24,   1, 1,22,   1, 1,19,   1, 4,15,   1, 7,13,
     .  1,10,13,   1,13,14,   1,12,15,   1,13,19,   1,15,20,   1,17,20,
     .  1,19,21,   1,21,21,   1,23,23,   1,25,25,   1,27,27,   1,28,28,
     .  1,29,29,   1,30,30,   1,30,27,   1,29,25,   1,29,18,   1,28,14,
     .  1,27,15,   1,25,15,   1,23,14,   1,22,10,   1,20,10,   1,18,10,
     .  1,16,10,   1,18, 7,   1,17, 1,   7,19, 1,   7,21, 1,   9,23, 1,
     . 11,25, 1,  14,26, 1,  17,27, 1,  20,28, 1,  24,29, 1,  29,29, 1,
     . 28,28, 1,  28,26, 1,  28,23, 1,  28,20, 1,  28,17, 1,  28,14, 1,
     . 28,11, 1,  28, 8, 1,  28, 5, 1,  27, 1, 1,  25, 1, 1,  22, 1, 1,
     . 19, 1, 1,  16, 1, 1,  13, 1, 1,   1, 1, 1/
      real, parameter ::  scale=.03125
      integer :: i,l0,l1
      real    :: rind,w0,w1
c
      write (*,101) 'defining sst color table:',' entries',
     .   ibase,' --',nbase
c
CKAL  Do a linear interpolation of color table to number of entries
CKAL  specified
      do i=1,nbase

         rind=float(ncol-1)*float(i-1)/float(nbase-1)
         l0=min(floor(rind)+1,ncol-1)
         l1=l0+1
         w0=l1-1-rind
         w1=1.-w0
         print '(i5,f10.2,i4,f10.2,i4,f10.2)',i+ibase,rind,l0,w0,l1,w1
         call gscr(IWKID,i+ibase,
     .            scale*(w0*rgb(1,l0)+w1*rgb(1,l1)),
     .            scale*(w0*rgb(2,l0)+w1*rgb(2,l1)),
     .            scale*(w0*rgb(3,l0)+w1*rgb(3,l1)))
C
      end do
C     do i=1,ncol
Ccc      write (*,100) i+ibase(1),(scale*rgb(l,i),l=1,3)
C       call gscr(IWKID,i+ibase,
C    .            scale*rgb(1,i),scale*rgb(2,i),scale*rgb(3,i))
C     enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
c
c
      subroutine gaudy(IWKID,nbase,ibase)
      implicit none
      integer nbase,ibase,IWKID
 
C --- initialize 'gaudy' color table
 
      integer :: num,ncol
      parameter (num=8,ncol=100)
      real red(num),grn(num),blu(num)
      data red/0.,0.,1.,1.,0.,0.,1.,1./
      data grn/0.,1.,0.,1.,0.,1.,0.,1./
      data blu/0.,0.,0.,0.,1.,1.,1.,1./
C --- 'expo' controls color interpolation. The larger the value, the
C --- more abrupt the transition between the prescribed color values
      real :: expo
      data expo/1.0/
      integer :: i,m
      real :: x,r,g,b,xx
 
      write (*,101) 'defining gaudy color table:',' entries',
     .   ibase,' --',nbase
      !do i=1,ncol
      !  x=1.+float((num-1)*(i-1))/float(ncol-1)
      do i=1,nbase
        x=1.+float((num-1)*(i-1))/float(nbase-1)
        m=min(num-1,int(x))
        x=x-float(m)
        if (x.lt.0.5) then
          xx=.5*(2.*x)**expo
        else
          xx=1.-.5*(2.*(1.-x))**expo
        end if
        r=red(m)*(1.-xx)+red(m+1)*xx
        g=grn(m)*(1.-xx)+grn(m+1)*xx
        b=blu(m)*(1.-xx)+blu(m+1)*xx
ccc        write (*,100) i+ibase(1),r,g,b
        call gscr(1,i+ibase,r,g,b)
      enddo
      return
 100  format (' color index',i4,'  <==> rgb combination',3f7.3)
 101  format(a32,a,i5,a,i5)
      end
C
C
C     subroutine twoton(nbase,ibase)
C     integer nbase,ibase(2)
C
C --- initialize table of 2-tone (blue and red) shades
C
C     parameter (ncol=64)
C     data satur/1./
C
C     if     (ibase(1).ne.0) then
C       write(*,101) 'existing 2-tone color table:',' entries',
C    .    ibase(1)+1,' --',ibase(2)
C       return
C     endif
C     ibase(1) = nbase
C     ibase(2) = nbase+ncol
C     nbase    = ibase(2)
C     write (*,101) 'defining 2-tone color table:',' entries',
C    .   ibase(1)+1,' --',ibase(2)
C     if     (mod(ncol,2).eq.0) then
C       ncolx=ncol-1
C     else
C       ncolx=ncol
C     endif
C     ii=0
C     do i=1,ncolx
C       x=1.-satur*abs(float(i-(ncolx+1)/2))/float((ncolx+1)/2-1)
C       if (i.le.(ncolx+1)/2) then
C         r=x
C         g=x
C         b=1.				!  blue shades
C       else
C         r=1.				!  red  shades
C         g=x
C         b=x
C       endif
C       ii=ii+1
Ccc        write (*,100) ii+ibase(1),r,g,b
C       call gscr(1,ii+ibase(1),r,g,b)
C       if     (ncolx.ne.ncol .and. i.eq.ncol/2) then
C         ii=ii+1
Ccc          write (*,100) ii+ibase(1),r,g,b
C         call gscr(1,ii+ibase(1),r,g,b)
C       endif
C     enddo
C     return
C100  format (' color index',i4,'  <==> rgb combination',3f7.3)
C101  format(a32,a,i5,a,i5)
C     end
C
C
C     subroutine fc100(nbase,ibase)
C     integer nbase,ibase(2)
C
C --- initialize table of 100 false colors.
C
C     parameter (ncol=100)
C     real rgb( 3,1:ncol)
C     real rgb0(3,10)
C     real rgb1(3,10)
C     real rgb2(3,10)
C     real rgb3(3,10)
C     real rgb4(3,10)
C     real rgb5(3,10)
C     real rgb6(3,10)
C     real rgb7(3,10)
C     real rgb8(3,10)
C     real rgb9(3,10)
C     equivalence (rgb(1, 1),rgb0(1,1)),
C    .            (rgb(1,11),rgb1(1,1)),
C    .            (rgb(1,21),rgb2(1,1)),
C    .            (rgb(1,31),rgb3(1,1)),
C    .            (rgb(1,41),rgb4(1,1)),
C    .            (rgb(1,51),rgb5(1,1)),
C    .            (rgb(1,61),rgb6(1,1)),
C    .            (rgb(1,71),rgb7(1,1)),
C    .            (rgb(1,81),rgb8(1,1)),
C    .            (rgb(1,91),rgb9(1,1))
C     data rgb0 /
C    .    1.0000, 1.0000, 1.0000,
C    .    0.9763, 0.9235, 0.9955,
C    .    0.9567, 0.8603, 0.9918,
C    .    0.9371, 0.7970, 0.9880,
C    .    0.9175, 0.7338, 0.9843,
C    .    0.8979, 0.6705, 0.9806,
C    .    0.8783, 0.6073, 0.9769,
C    .    0.8586, 0.5440, 0.9731,
C    .    0.8390, 0.4808, 0.9694,
C    .    0.8194, 0.4175, 0.9656 /
C     data rgb1 /
C    .    0.7998, 0.3543, 0.9619,
C    .    0.7802, 0.2910, 0.9582,
C    .    0.7606, 0.2278, 0.9544,
C    .    0.7410, 0.1645, 0.9507,
C    .    0.7241, 0.0380, 0.9736,
C    .    0.6615, 0.0127, 0.9807,
C    .    0.5988, 0.0000, 0.9860,
C    .    0.5362, 0.0000, 0.9900,
C    .    0.4892, 0.0000, 0.9900,
C    .    0.4422, 0.0000, 0.9900 /
C     data rgb2 /
C    .    0.3800, 0.0000, 0.9900,
C    .    0.3178, 0.0000, 0.9900,
C    .    0.2556, 0.0000, 0.9900,
C    .    0.1934, 0.0000, 0.9900,
C    .    0.1312, 0.0000, 0.9702,
C    .    0.0000, 0.0000, 0.9286,
C    .    0.0030, 0.0287, 0.5463,
C    .    0.0101, 0.1262, 0.2249,
C    .    0.0217, 0.2002, 0.2851,
C    .    0.0333, 0.2741, 0.3453 /
C     data rgb3 /
C    .    0.0537, 0.3375, 0.4587,
C    .    0.0686, 0.3867, 0.5140,
C    .    0.0803, 0.4279, 0.5393,
C    .    0.0983, 0.4944, 0.5544,
C    .    0.1221, 0.5670, 0.5852,
C    .    0.1457, 0.6150, 0.6241,
C    .    0.1692, 0.6629, 0.6629,
C    .    0.1941, 0.7005, 0.7005,
C    .    0.2190, 0.7382, 0.7382,
C    .    0.2439, 0.7758, 0.7758 /
C     data rgb4 /
C    .    0.2892, 0.8330, 0.8330,
C    .    0.3223, 0.8724, 0.8724,
C    .    0.3567, 0.9080, 0.9080,
C    .    0.3911, 0.9435, 0.9435,
C    .    0.4317, 0.9800, 0.9800,
C    .    0.3912, 0.9701, 0.9302,
C    .    0.3478, 0.9537, 0.8740,
C    .    0.2950, 0.9455, 0.7796,
C    .    0.2442, 0.9291, 0.6784,
C    .    0.2083, 0.9067, 0.6100 /
C     data rgb5 /
C    .    0.1755, 0.8390, 0.5624,
C    .    0.1427, 0.7714, 0.5149,
C    .    0.1113, 0.6867, 0.4606,
C    .    0.0661, 0.6170, 0.3500,
C    .    0.0530, 0.5317, 0.2928,
C    .    0.0400, 0.4464, 0.2657,
C    .    0.0282, 0.3873, 0.2226,
C    .    0.0746, 0.4577, 0.1220,
C    .    0.1522, 0.4797, 0.0000,
C    .    0.1886, 0.5369, 0.0000 /
C     data rgb6 /
C    .    0.2021, 0.5941, 0.0000,
C    .    0.2454, 0.6513, 0.0000,
C    .    0.2903, 0.7080, 0.0000,
C    .    0.3362, 0.7643, 0.0000,
C    .    0.3901, 0.7873, 0.0000,
C    .    0.4449, 0.8102, 0.0000,
C    .    0.5006, 0.8330, 0.0000,
C    .    0.5573, 0.8558, 0.0000,
C    .    0.6150, 0.8785, 0.0000,
C    .    0.6700, 0.9012, 0.0000 /
C     data rgb7 /
C    .    0.7334, 0.9238, 0.0000,
C    .    0.8022, 0.9407, 0.0000,
C    .    0.8774, 0.9480, 0.0000,
C    .    0.9500, 0.9500, 0.0000,
C    .    0.9150, 0.8946, 0.0000,
C    .    0.8996, 0.8390, 0.0000,
C    .    0.8996, 0.7797, 0.0000,
C    .    0.8996, 0.7185, 0.0000,
C    .    0.8996, 0.6570, 0.0000,
C    .    0.8996, 0.5959, 0.0000 /
C     data rgb8 /
C    .    0.8996, 0.5350, 0.0000,
C    .    0.8996, 0.4741, 0.0000,
C    .    0.8996, 0.4132, 0.0000,
C    .    0.8996, 0.3522, 0.0000,
C    .    0.8996, 0.2913, 0.0000,
C    .    0.8996, 0.2304, 0.0000,
C    .    0.8996, 0.1695, 0.0000,
C    .    0.8875, 0.1216, 0.0000,
C    .    0.8754, 0.0736, 0.0000,
C    .    0.8392, 0.0000, 0.0000 /
C     data rgb9 /
C    .    0.7820, 0.0000, 0.0000,
C    .    0.7152, 0.0000, 0.0000,
C    .    0.6499, 0.0000, 0.0000,
C    .    0.5846, 0.0000, 0.0000,
C    .    0.5092, 0.0000, 0.0000,
C    .    0.4239, 0.0000, 0.0000,
C    .    0.3386, 0.0000, 0.0000,
C    .    0.2532, 0.0000, 0.0000,
C    .    0.1679, 0.0000, 0.0000,
C    .    0.1000, 0.0000, 0.0000 /
C
C     if     (ibase(1).ne.0) then
C       write(*,101) 'existing fc100 color table:',' entries',
C    .    ibase(1)+1,' --',ibase(2)
C       return
C     endif
C     ibase(1) = nbase
C     ibase(2) = nbase+ncol
C     nbase    = ibase(2)
C     write (*,101) 'defining fc100 color table:',' entries',
C    .   ibase(1)+1,' --',ibase(2)
C     do i=1,ncol
Ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
C       call gscr(1,i+ibase(1),rgb(1,i),rgb(2,i),rgb(3,i))
C     enddo
C     return
C100  format (' color index',i4,'  <==> rgb combination',3f7.3)
C101  format(a32,a,i5,a,i5)
C     end
C
C
C     subroutine ifc100(nbase,ibase)
C     integer nbase,ibase(2)
C
C --- initialize inverted table of 100 false colors.
C
C     parameter (ncol=100)
C     real rgb( 3,1:ncol)
C     real rgb0(3,10)
C     real rgb1(3,10)
C     real rgb2(3,10)
C     real rgb3(3,10)
C     real rgb4(3,10)
C     real rgb5(3,10)
C     real rgb6(3,10)
C     real rgb7(3,10)
C     real rgb8(3,10)
C     real rgb9(3,10)
C     equivalence (rgb(1, 1),rgb0(1,1)),
C    .            (rgb(1,11),rgb1(1,1)),
C    .            (rgb(1,21),rgb2(1,1)),
C    .            (rgb(1,31),rgb3(1,1)),
C    .            (rgb(1,41),rgb4(1,1)),
C    .            (rgb(1,51),rgb5(1,1)),
C    .            (rgb(1,61),rgb6(1,1)),
C    .            (rgb(1,71),rgb7(1,1)),
C    .            (rgb(1,81),rgb8(1,1)),
C    .            (rgb(1,91),rgb9(1,1))
C     data rgb0 /
C    .    1.0000, 1.0000, 1.0000,
C    .    0.9763, 0.9235, 0.9955,
C    .    0.9567, 0.8603, 0.9918,
C    .    0.9371, 0.7970, 0.9880,
C    .    0.9175, 0.7338, 0.9843,
C    .    0.8979, 0.6705, 0.9806,
C    .    0.8783, 0.6073, 0.9769,
C    .    0.8586, 0.5440, 0.9731,
C    .    0.8390, 0.4808, 0.9694,
C    .    0.8194, 0.4175, 0.9656 /
C     data rgb1 /
C    .    0.7998, 0.3543, 0.9619,
C    .    0.7802, 0.2910, 0.9582,
C    .    0.7606, 0.2278, 0.9544,
C    .    0.7410, 0.1645, 0.9507,
C    .    0.7241, 0.0380, 0.9736,
C    .    0.6615, 0.0127, 0.9807,
C    .    0.5988, 0.0000, 0.9860,
C    .    0.5362, 0.0000, 0.9900,
C    .    0.4892, 0.0000, 0.9900,
C    .    0.4422, 0.0000, 0.9900 /
C     data rgb2 /
C    .    0.3800, 0.0000, 0.9900,
C    .    0.3178, 0.0000, 0.9900,
C    .    0.2556, 0.0000, 0.9900,
C    .    0.1934, 0.0000, 0.9900,
C    .    0.1312, 0.0000, 0.9702,
C    .    0.0000, 0.0000, 0.9286,
C    .    0.0030, 0.0287, 0.5463,
C    .    0.0101, 0.1262, 0.2249,
C    .    0.0217, 0.2002, 0.2851,
C    .    0.0333, 0.2741, 0.3453 /
C     data rgb3 /
C    .    0.0537, 0.3375, 0.4587,
C    .    0.0686, 0.3867, 0.5140,
C    .    0.0803, 0.4279, 0.5393,
C    .    0.0983, 0.4944, 0.5544,
C    .    0.1221, 0.5670, 0.5852,
C    .    0.1457, 0.6150, 0.6241,
C    .    0.1692, 0.6629, 0.6629,
C    .    0.1941, 0.7005, 0.7005,
C    .    0.2190, 0.7382, 0.7382,
C    .    0.2439, 0.7758, 0.7758 /
C     data rgb4 /
C    .    0.2892, 0.8330, 0.8330,
C    .    0.3223, 0.8724, 0.8724,
C    .    0.3567, 0.9080, 0.9080,
C    .    0.3911, 0.9435, 0.9435,
C    .    0.4317, 0.9800, 0.9800,
C    .    0.3912, 0.9701, 0.9302,
C    .    0.3478, 0.9537, 0.8740,
C    .    0.2950, 0.9455, 0.7796,
C    .    0.2442, 0.9291, 0.6784,
C    .    0.2083, 0.9067, 0.6100 /
C     data rgb5 /
C    .    0.1755, 0.8390, 0.5624,
C    .    0.1427, 0.7714, 0.5149,
C    .    0.1113, 0.6867, 0.4606,
C    .    0.0661, 0.6170, 0.3500,
C    .    0.0530, 0.5317, 0.2928,
C    .    0.0400, 0.4464, 0.2657,
C    .    0.0282, 0.3873, 0.2226,
C    .    0.0746, 0.4577, 0.1220,
C    .    0.1522, 0.4797, 0.0000,
C    .    0.1886, 0.5369, 0.0000 /
C     data rgb6 /
C    .    0.2021, 0.5941, 0.0000,
C    .    0.2454, 0.6513, 0.0000,
C    .    0.2903, 0.7080, 0.0000,
C    .    0.3362, 0.7643, 0.0000,
C    .    0.3901, 0.7873, 0.0000,
C    .    0.4449, 0.8102, 0.0000,
C    .    0.5006, 0.8330, 0.0000,
C    .    0.5573, 0.8558, 0.0000,
C    .    0.6150, 0.8785, 0.0000,
C    .    0.6700, 0.9012, 0.0000 /
C     data rgb7 /
C    .    0.7334, 0.9238, 0.0000,
C    .    0.8022, 0.9407, 0.0000,
C    .    0.8774, 0.9480, 0.0000,
C    .    0.9500, 0.9500, 0.0000,
C    .    0.9150, 0.8946, 0.0000,
C    .    0.8996, 0.8390, 0.0000,
C    .    0.8996, 0.7797, 0.0000,
C    .    0.8996, 0.7185, 0.0000,
C    .    0.8996, 0.6570, 0.0000,
C    .    0.8996, 0.5959, 0.0000 /
C     data rgb8 /
C    .    0.8996, 0.5350, 0.0000,
C    .    0.8996, 0.4741, 0.0000,
C    .    0.8996, 0.4132, 0.0000,
C    .    0.8996, 0.3522, 0.0000,
C    .    0.8996, 0.2913, 0.0000,
C    .    0.8996, 0.2304, 0.0000,
C    .    0.8996, 0.1695, 0.0000,
C    .    0.8875, 0.1216, 0.0000,
C    .    0.8754, 0.0736, 0.0000,
C    .    0.8392, 0.0000, 0.0000 /
C     data rgb9 /
C    .    0.7820, 0.0000, 0.0000,
C    .    0.7152, 0.0000, 0.0000,
C    .    0.6499, 0.0000, 0.0000,
C    .    0.5846, 0.0000, 0.0000,
C    .    0.5092, 0.0000, 0.0000,
C    .    0.4239, 0.0000, 0.0000,
C    .    0.3386, 0.0000, 0.0000,
C    .    0.2532, 0.0000, 0.0000,
C    .    0.1679, 0.0000, 0.0000,
C    .    0.1000, 0.0000, 0.0000 /
C
C     if     (ibase(1).ne.0) then
C       write(*,101) 'existing ifc100 color table:',' entries',
C    .    ibase(1)+1,' --',ibase(2)
C       return
C     endif
C     ibase(1) = nbase
C     ibase(2) = nbase+ncol
C     nbase    = ibase(2)
C     write (*,101) 'defining ifc100 color table:',' entries',
C    .   ibase(1)+1,' --',ibase(2)
C     do i=1,ncol
C       ii=ncol+1-i
Ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
C       call gscr(1,i+ibase(1),rgb(1,ii),rgb(2,ii),rgb(3,ii))
C     enddo
C     return
C100  format (' color index',i4,'  <==> rgb combination',3f7.3)
C101  format(a32,a,i5,a,i5)
C     end
C
C
C     subroutine fc20(nbase,ibase)
C     integer nbase,ibase(2)
C
C --- initialize table of 20 false colors.
C
C     parameter (ncol=20)
C     real rgb( 3,1:ncol)
C     data rgb   /
C    + .00,.38,.50, .00,.50,.63, .00,.63,.75, .00,.75,.88,
C    + .00,.88,1.0, .00,1.0,1.0, .20,.99,.99, .40,.99,.99,
C    + .60,.99,.99, .80,.99,.99, .99,.99,.00, .99,.88,.00,
C    + .99,.75,.00, .99,.63,.00, .99,.50,.00, .99,.38,.00,
C    + .99,.25,.00, .99,.13,.00, .75,.00,.00, .50,.00,.00 /
C
C     if     (ibase(1).ne.0) then
C       write(*,101) 'existing fc20 color table:',' entries',
C    .    ibase(1)+1,' --',ibase(2)
C       return
C     endif
C     ibase(1) = nbase
C     ibase(2) = nbase+ncol
C     nbase    = ibase(2)
C     write (*,101) 'defining fc20 color table:',' entries',
C    .   ibase(1)+1,' --',ibase(2)
C     do i=1,ncol
Ccc      write (*,100) i+ibase(1),(rgb(l,i),l=1,3)
C       call gscr(1,i+ibase(1),rgb(1,i),rgb(2,i),rgb(3,i))
C     enddo
C     return
C100  format (' color index',i4,'  <==> rgb combination',3f7.3)
C101  format(a32,a,i5,a,i5)
C     end
C
C
