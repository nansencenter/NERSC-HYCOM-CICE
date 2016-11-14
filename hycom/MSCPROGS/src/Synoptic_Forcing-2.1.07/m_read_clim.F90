module m_read_clim
contains
subroutine read_clim(nmo,field,filen,diag,clmflag,plon,plat,depths)
! Load data from file 'filen' into array 'field'
! and dumps it in integer format if 'diag'=1.
!KAL -- changed to read one month at a time
   use mod_xc
   use mod_za
   use m_bilin_ecmwf2
   use m_ncvar_read
   use m_era40_fix
   implicit none

   integer, intent(in) :: nmo, diag
   character(len=*), intent(in) :: filen
   character(len=*), intent(in) :: clmflag
   real, intent(out) :: field (idm,jdm)
   real, intent(in), dimension(idm,jdm) :: plon,plat,depths

   integer nxl,nyl,nzl,fact
   integer :: iutil (idm,jdm)
   integer :: ip    (idm,jdm)
   real    :: tmpfld(idm,jdm)
   real    :: gfld  (idm,jdm)
   real grdn,ypn,xpn,ypo,rfact,vmax,vmin,dd,dw,offs,fl
   integer :: i,j,k,iind
   logical :: ex,exa, exb
   real    :: hmin,hmax,hmin2,hmax2
   character*40 :: filena,varid,filenb
   real :: dummy(1,1,1)
   real, dimension(:,:,:), allocatable :: era40data
   real :: mlon (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
   integer :: nlon,nlat
   real    :: dlon,dlat,flon,flat,llon,llat
   integer :: index1, index2

   ip = 1


if (clmflag=='era40') then

   ! Find var id to read - var id is embedded in file name
   index1=index(filen,'climatology_')
   index2=index(filen,'.nc')
   if (index1==0 .or. index2==0) then
      if (mnproc==1) then
         write(lp,*) 'read_clim: Can not get indexes'
         write(lp,*) 'file is '//trim(filen)
         call flush(lp)
      end if
      call xcstop('(read_clim)')
      stop '(read_clim)'
   end if

   varid=filen(index1+12:index2-1)
   !print *,varid

   ! Get lon/lat info:
   call ncvar_read(filen,'Ni' ,dummy(1,1,1)  , 1,1,1,1,1) ; nlon=dummy(1,1,1)
   call ncvar_read(filen,'Di' ,dummy(1,1,1)  , 1,1,1,1,1) ; dlon=dummy(1,1,1)
   call ncvar_read(filen,'Lo1',dummy(1,1,1)  , 1,1,1,1,1) ; flon=dummy(1,1,1)
   call ncvar_read(filen,'Lo2',dummy(1,1,1)  , 1,1,1,1,1) ; llon=dummy(1,1,1)
   call ncvar_read(filen,'Nj' ,dummy(1,1,1)  , 1,1,1,1,1) ; nlat=dummy(1,1,1)
   call ncvar_read(filen,'Dj' ,dummy(1,1,1)  , 1,1,1,1,1) ; dlat=dummy(1,1,1)
   call ncvar_read(filen,'La1',dummy(1,1,1)  , 1,1,1,1,1) ; flat=dummy(1,1,1)
   call ncvar_read(filen,'La2',dummy(1,1,1)  , 1,1,1,1,1) ; llat=dummy(1,1,1)

   allocate(era40data(nlon,nlat,1))

   ! Read data
   call ncvar_read(filen,trim(varid),era40data,nlon,nlat,1,nmo,nmo)

   ! fix for era40 climatology precipitation
   if (trim(varid)=='TP') then
      call era40_fix(varid,era40data,nlon,nlat,flon,flat,dlon,dlat)
   end if


   ! Do bilinear interpolation of data
   mlon = plon ; where (mlon < 0.) mlon=mlon+360.

   !print *,flon,flat
!diag   print *,'read_clim:',trim(varid),nmo,nlon,nlat,minval(era40data),maxval(era40data)
   call bilin_ecmwf2(era40data,nlon,nlat,flon,flat,dlon,dlat, &
                          field,mlon,plat,depths)   

   deallocate(era40data)

!diag   if (trim(varid)=='T2M_sfc') then
!diag      call zaiopf ('tst.a','replace',11)
!diag      call zaiowr(field,ip,.false.,hmin,hmax,11,.true.)
!diag      call zaiocl(11)
!diag   end if
   

elseif (clmflag=='old') then

   vmax=-9999.
   vmin= 9999.

   filena=filen
   filenb=filen
   iind=index(filena,'.forc',.true.)
   filena(iind:iind+4)='.a   '
   filenb(iind:iind+4)='.b   '

   ! Check for presence of original ".forc" files
   inquire(file=filen,exist=ex) 


   if (.not.ex) then
      if (mnproc==1) then
         write(lp,'(a)') 'Can not find file '//filen
         call flush(lp)
      end if
      call xcstop('(read_clim)')
      stop '(read_clim)'
   end if


! KAL - Skipped generation of .ab - files
      
      if (mnproc==1) then
         write(*,'(''I read: '',a,i5)')trim(filen)//' for month ',nmo
         call flush(lp)
      end if
      !close(10)
      open(10,FILE=filen,FORM='formatted',STATUS= 'UNKNOWN')
      read(10,200)nxl,nyl,nzl,fact
      if(idm.NE.nxl.OR.jdm.NE.nyl.OR.nzl.NE.12) then
         CLOSE(10)
         if (mnproc==1) then
            write(lp,*)'Wrong dimension in file ',filen    
            write(lp,*)'File dimension: (nx,ny,nmo)',nxl,nyl,nzl
            write(lp,*)'Model dimension:(nx,ny,nmo)',idm,jdm,12
            call flush(lp)
         end if
         call xcSTOP('read_clim')
         STOP 'read_clim'
      end if

      rfact=1./FLOAT(fact)


      do k=1,12 ! Cycle all months
         read(10,201) iutil

         ! Calculate field
         do i=1,idm
            do j=1,jdm
               tmpfld(i,j)=FLOAT(iutil(i,j))*rfact
               fl=.5+SIGN(.5,tmpfld(i,j)+998.)        !=0 if no dat
               vmin=MIN(vmin,tmpfld(i,j)*fl+(1.-fl)*vmin)
               vmax=MAX(vmax,tmpfld(i,j)*fl+(1.-fl)*vmax)
            enddo
         enddo


         ! if k== nmo, keep it
         if (k==nmo) then
            gfld=tmpfld
            exit ! Exits month loop
         end if


      end do
      close(10)

! --- ---------------------------- Diagnostics -------------
      dd=vmax-vmin
      offs=0.
      if(vmin.LT.0.)offs=-vmin
      if(vmax.LT.1.AND.dd.LT.1.)dw=100.
      if(vmax.GE.1.AND.dd.LT.100.)dw=1.  
      if(dd.GT.100.AND.dd.LT.1000.)dw=.1
      if(dd.GT.1000.AND.dd.LT.10000.)dw=.01
      if(dd.GT.10000.AND.dd.LT.100000.)dw=.001
      if(vmin.GT.200.AND.vmax.LT.400) THEN
        offs=-223.15
        dw=1.
      endif
      !write(*,*)'File: ',filen(:)   ! Already stated
      !write(*,202)vmin,vmax,offs,dw ! Who cares ...

      if(diag.EQ.1)THEN
       do k=1,nzl
        write(*,203)k,filen
        do i=idm,1,-1
         fl=-999.*(.5-SIGN(.5,offs+gfld(i,j)-.1))
         write(*,'(67I2)')i,(INT(-999.*(.5-SIGN(.5,offs+gfld(i,j)-.1))&
                   +(offs+gfld(i,j))*dw) ,j=jdm-3,3,-1)
        enddo
       enddo
      endif

      field=gfld
end if


 200  format(4I8)
 201  format(15I8)
 202  format('Min: ',f9.2,' Max: ',f9.2,' Offs: ',f9.2,' Fact ',f9.2)
 203  format('Month: ',i2,'File: ',a10)
      continue

! KAL -- Based on old "reddat"  -- Changes below
! KAL -- Returns one month. Produces binary file versions of .forc files -- 22062005
! KAL -- Allows for use of  Produces binary file versions of .forc files -- 22062005
end subroutine read_clim
end module m_read_clim
