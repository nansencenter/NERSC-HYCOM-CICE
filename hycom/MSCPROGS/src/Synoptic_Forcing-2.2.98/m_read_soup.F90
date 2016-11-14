module m_read_soup

contains

   subroutine read_soup(rt,plon,plat,synuwind,synvwind,synwndspd,syntaux,syntauy) 
      use mod_xc
!      use mod_year_info
      use mod_year_info22
      use m_interpug
      !use m_spherdist
      !use m_tecfld 
      implicit none
     
      
      real,    parameter :: undef=-999.

      ! For the final version
      integer, parameter :: nlon=113,nlat=153
      real,    parameter :: flon=97.0 , flat=-9.00
      real,    parameter :: dlat=0.250, dlon=0.250
      character(len=6),   parameter ::  delta= '3 hr'

      type(year_info) , intent(in) :: rt

      !real, allocatable :: Field(:)
      real :: r4(4),r1,tmp
      character(len=3) :: crec
      integer :: i3(3)
      !integer, allocatable, dimension(:) :: ilons,ilats,point_index
      real, allocatable, dimension(:,:) :: fld1,fld2, &
      windx,windy,lon,lat,weights
      integer :: i,icnt,j
      integer :: ddhhmm
      integer :: yymm
      integer :: ios
      integer :: jul,hh,month,yy,dd,minut
      integer :: whh,lwmonth,lwyy,lwdd,lwminut,lwhh
      integer :: himonth,hiyy,hidd,himinut,hihh
      integer :: totrec,irec,prsrec,ind
      logical :: first, speed, direction,pressure
      integer :: speedrec, dirrec,iol
      character(len=80) :: header
      character(len=52) :: preambl
      integer*8 :: ifirst, ilast
      integer*8 :: time2
      integer :: hirec, lwrec, baserec
  
      real :: d1,d2,d3,d4,d
      real w,d_max     
      real, dimension(1:idm,1:jdm) :: mlon, tmpu, tmpv, tmpw
      real, dimension(idm,1:jdm), intent(in) :: plat,plon
      real, dimension(idm,1:jdm), intent(inout) :: synuwind, synvwind
      real, dimension(idm,1:jdm), intent(inout) :: syntaux, syntauy
      real, dimension(idm,1:jdm), intent(inout) :: synwndspd
      real, dimension(idm,1:jdm) :: tmpuin, tmpvin, tmpwin
      logical :: notfound,ex

      real pi,rad,wfact

      integer :: ip1,im1, jp1,jm1
      integer :: nlat2,nlon2,ios2
      real :: dlon2,dlat2,flat2,flon2
      real :: wndfac, cd_new, w4
      real,parameter :: airdns  =    1.2
      real,parameter :: cd      =    0.0012

      integer :: margin 

      character(len=9) path0
      character(len=30) Sname 
      real, external :: spherdist
      integer, parameter :: itype=0 ! Sets interpolation type bilinear

      mlon = plon ; where (mlon < 0.) mlon=mlon+360.

      pi=acos(-1.)
!      print *,pi
      rad=pi/180.

      !allocate(Field(npoint))
      allocate(windx(nlon,nlat))
      allocate(windy(nlon,nlat))

      allocate(lon(nlon,nlat))
      allocate(lat(nlon,nlat))
      allocate(weights(nlon,nlat))

      allocate(fld1(nlat,nlon))
      allocate(fld2(nlat,nlon))

      speed=.false.
      direction=.false.
      pressure=.false.
 
!      if (winds=='soup') then
         path0(1:9)='./Soup/'
!      endif   1994_Wind.ZIP
! This we don't care about at the moment...
      !inquire(file=trim(path0)//rt%cyy//rt%cmm/'_SOUP_3Hr.Win',exist=ex)
      !if (.not.ex) then WINDX001.dat

      !   if(mnproc==1) print *, 'read_soup:file does not exist:',trim(path0),rt%iyy,rt%imm,'_SOUP_3Hr.Win'
      !   call xcstop ('(read_soup)') 
      !   stop
      !endif
      Sname='_SOUP_3Hr.Win'
      Sname='_SEAFINE_Corr2_3Hr.WIN'

      inquire(exist=ex,file=trim(path0)//rt%cyy//rt%cmm//trim(Sname))
      if (.not.ex) then
         print *,'Cant find '//trim(path0)//rt%cyy//rt%cmm//trim(Sname)
         call flush(lp)
         call xcstop('(read_soup)')
         stop '(read_soup)'
      else
          print *,'Reading '//trim(path0)//rt%cyy//rt%cmm//trim(Sname)
      end if
      open(20,file=trim(path0)//rt%cyy//rt%cmm//trim(Sname),form='formatted')
      Read(20,'(a80)') header
      print *,'Header : ',header
      ind=index(header,'Format')
      read(header(ind+6:),*) ifirst,ilast
      print *,header
      print *,ifirst,ilast

      ! Calculate first record number
      yy=ifirst/1000000
      month=(ifirst-yy*1000000)/10000
      dd=(ifirst-yy*1000000-month*10000)/100
      hh=(ifirst-yy*1000000-month*10000-dd*100)
      lwyy=yy ; lwminut=0 ; lwhh=hh; lwdd=dd; lwmonth=month
      print *,'first=',ifirst,lwyy,lwmonth,lwdd,lwhh


      ! Calculate last record number
      yy=ilast/1000000
      month=(ilast-yy*1000000)/10000
      dd=(ilast-yy*1000000-month*10000)/100
      hh=(ilast-yy*1000000-month*10000-dd*100)
      hiyy=yy ; himinut=0 ; hihh=hh; hidd=dd; himonth=month
      print *,'last=',ilast,hiyy,himonth,hidd,hihh

      ! Data every 30 min
      if (delta=='30 min') then
         jul=datetojulian(lwyy,lwmonth,lwdd,1950,1,1)
         baserec=jul*48+lwhh*2+lwminut/30-1
         lwrec=jul*48+lwhh*2+lwminut/30 - baserec
         jul=datetojulian(hiyy,himonth,hidd,1950,1,1)
         hirec=jul*48+hihh*2+himinut/30 - baserec
      else if (delta=='15 min') then
         jul=datetojulian(lwyy,lwmonth,lwdd,1950,1,1)
         baserec=jul*96+lwhh*4+lwminut/15-1
         lwrec=jul*96+lwhh*4+lwminut/15 - baserec
         jul=datetojulian(hiyy,himonth,hidd,1950,1,1)
         hirec=jul*96+hihh*4+himinut/15 - baserec
      else if (delta=='3 hr') then
         jul=datetojulian(lwyy,lwmonth,lwdd,1950,1,1)
         baserec=jul*8+lwhh/3
         lwrec=jul*8+lwhh/3 - baserec
         jul=datetojulian(hiyy,himonth,hidd,1950,1,1)
         hirec=jul*8+hihh/3 - baserec
      else
         print *,'Unknown delta ...'
         call xcstop(' ')
         stop
      end if

      print *,'baserec,low,hi:',baserec,lwrec,hirec

      ! Create ll grid
      do i=1,nlon
      do j=1,nlat
         lon(i,j)=(i-1)*dlon+flon
         lat(i,j)=(j-1)*dlat+flat
      end do
      end do
     
      ios=0
      notfound=.true.
      do while(ios==0 .and. notfound)
         ! Read Data header
         !iLat= 105iLong= 145DX= 0.125DY= 0.125SWLat=  18.000SWLon= -98.000DT=200209280000
         Read(20,98,iostat=ios) nlat2,nlon2,dlon2,dlat2,flat2,flon2,yy,month,dd,hh,minut
         !write(*,98) nlat2,nlon2,dlon2,dlat2,flat2,flon2,yy,month,dd,hh,minut
         !print *,rt%iyy,rt%imm,rt%idm,rt%ihh

         ! rt%idm is day of wmonth, starting from 0 !
         if (dd==rt%idm+1 .and. hh==rt%ihh .and.yy==rt%iyy .and. rt%imm==month) then 
            write(*,98) nlat2,nlon2,dlon2,dlat2,flat2,flon2,yy,month,dd,hh,minut
            notfound=.false.
         end if


         if (ios==0) then

            ! Consistency check
            if (nlon/=nlon2 .or. nlat/=nlat2 ) then 
               print *,'Grid size mismatch :'
               print *,'file  :',nlon2,nlat2
               print *,'preset:',nlon ,nlat
            else if (dlon/=dlon2 .or. dlat/=dlat2 ) then 
               print *,'Grid spacing mismatch :'
               print *,'file  :',dlon2,dlat2
               print *,'preset:',dlon ,dlat
            else if (flon/=flon2 .or. flat/=flat2 ) then 
               print *,'First Point  mismatch :'
               print *,'file  :',flon2,flat2
               print *,'preset:',flon ,flat
            end if

            read(20,99,iostat=ios) ((fld1(i,j),j=1,nlon),i=1,nlat)
            read(20,99,iostat=ios) ((fld2(i,j),j=1,nlon),i=1,nlat)




            ! Time info
            if (delta=='30 min') then
               jul=datetojulian(yy,month,dd,1950,1,1)
               totrec=jul*48+hh*2+minut/30
            else if (delta=='15 min') then
               jul=datetojulian(yy,month,dd,1950,1,1)
               totrec=jul*96+hh*4+minut/15
            else if (delta=='3 hr') then
               jul=datetojulian(yy,month,dd,1950,1,1)
               totrec=jul*8+hh/3
            end if
            !write(6,'(a8,i5,i3,i3,i3,i3)',advance='no') '    now=',yy,month,dd,hh,minut
            !write(6,*) '   now rec =',baserec,totrec,totrec-baserec



            ! Depending on header info, we create two different files
            irec=totrec-baserec


            windx=transpose(fld1)
            windy=transpose(fld2)
         end if
      end do
      close(20)


      d_max=300*1000 ! 1000 km
      ! Find the smallest distance of the four possible:
      do i=1,nlon
         do j=1,nlat
            d1=spherdist(lon(i,j),lat(i,j),lon(1   ,j),lat(1   ,j))
            d2=spherdist(lon(i,j),lat(i,j),lon(nlon,j),lat(nlon,j))
            d3=spherdist(lon(i,j),lat(i,j),lon(i,nlat),lat(i,nlat))
            d4=spherdist(lon(i,j),lat(i,j),lon(i,   1),lat(i,   1))
            d=min(d1,d2,d3,d4)
            !w=max(0.,min(1.,(d/d_max)**3))
            w=max(0.,min(1.,(d/d_max)))
            weights(i,j)=w
         end do
      end do
           
      tmpu=0.
      tmpv=0.
      tmpw=0.

      call interpug  (windx  ,nlon,nlat,flon,flat,dlon,dlat,tmpu,mlon,plat,itype)   
      call interpug  (windy  ,nlon,nlat,flon,flat,dlon,dlat,tmpv,mlon,plat,itype)   
      call interpug  (weights,nlon,nlat,flon,flat,dlon,dlat,tmpw,mlon,plat,itype)   
      call rotate(tmpu,tmpv,plat,plon,idm,jdm,'l2m') 

      do i=1,idm
         do j=1,jdm
            d1=spherdist(plon(i,j),plat(i,j),plon(1 , j),plat(1 ,j ))
            d2=spherdist(plon(i,j),plat(i,j),plon(idm, j),plat(idm,j ))
            d3=spherdist(plon(i,j),plat(i,j),plon(i ,jdm),plat(i ,jdm))
            d4=spherdist(plon(i,j),plat(i,j),plon(i , 1),plat(i , 1))
            d=min(d1,d2,d3,d4)
            !w=max(0.,min(1.,(d/d_max)**3))
            w=max(0.,min(1.,(d/d_max)))
            tmpw(i,j)=min(w,tmpw(i,j))
         end do
      end do


     
      !print *,irec
      !write(crec,'(i3.3)') irec
      !if (mod(irec,10)==0) then
      !   call tecfld2('INWIND'//crec,nx,ny,mlon,modlat,synuwind,synvwind)
      !   call tecfld2('NWWIND'//crec,nx,ny,mlon,modlat,tmpu*tmpw,tmpv*tmpw)
      !end if

!      !calling xcaput to get tmpu, tmpv and tmpw on the same dimension as synuwind.
!      call xcaput(tmpu,tmpuin,0)
!      call xcaput(tmpv,tmpvin,0)
!      call xcaput(tmpw,tmpwin,0)
!

      synuwind=synuwind*(1.-tmpw)+tmpu*tmpw
      synvwind=synvwind*(1.-tmpw)+tmpv*tmpw

      ! Calculate corresponding drag and wind speed - from read_ecmwf
      !$OMP PARALLEL DO PRIVATE (wndfac,cd_new,w4,i,j) &
#ifdef KARALIGHT
     do j=1,jdm
         do i=1,idm
            im1=i-1
            ip1=i+1
            synwndspd(i,j)=sqrt(synuwind(i,j)*synuwind(i,j)+ &
                                synvwind(i,j)*synvwind(i,j))
            wndfac=MAX(2.5,MIN(32.5,synwndspd(i,j)))
            cd_new = 1.0E-3*(.692 + .0710*wndfac - .000700*wndfac**2)
            wfact=synwndspd(i,j)*airdns*cd_new
            syntaux(i,j)=synuwind(i,j)*wfact
            syntauy(i,j)=synvwind(i,j)*wfact
      end do
      end do
#else
      do j=1,jdm
         jm1=j-1
         jp1=j+1
         do i=1,idm
            im1=i-1
            ip1=i+1

            synwndspd(i,j)=sqrt(synuwind(i,j)*synuwind(i,j)+ &
                                synvwind(i,j)*synvwind(i,j))

            wndfac=(1.+sign(1.,synwndspd(i,j)-11.))*.5
            cd_new=(0.49+0.065*synwndspd(i,j))*1.0e-3*wndfac+cd*(1.-wndfac)

            w4=.25*(synvwind(im1,jp1)+synvwind(i,jp1)+ &
                    synvwind(im1,j  )+synvwind(i,j  ))
            syntaux(i,j)=synuwind(i,j)*sqrt(synuwind(i,j)*synuwind(i,j)+w4*w4)*airdns*cd_new

            w4=.25*(synuwind(i,jm1)+synuwind(ip1,jm1)+ &
                    synuwind(ip1,j)+synuwind(i  ,j  ))
            syntauy(i,j)=synvwind(i,j)*sqrt(synvwind(i,j)*synvwind(i,j)+w4*w4)*airdns*cd_new
         enddo
      enddo
#endif
      !$OMP END PARALLEL DO

      !if (mod(irec,10)==0) then
      !   call tecfld2('OUTWIND'//crec,nx,ny,mlon,modlat,synuwind,synvwind)
      !   call tecfld2('NWWINDWEIGHT'//crec,nx,ny,mlon,modlat,tmpw,tmpw)
      !end if
              
         

      !stop

 98 format("iLat= ",i3,"iLong= ",i3,"DX=",f6.4,"DY=",f6.4,"SWLat=",f8.5,"SWLon=",f8.5, &
	   "DT=",i4,i2,i2,i2,i2)
 99 format(8f10.4)

   end subroutine read_soup
end module m_read_soup 
