      module m_har_mdomain
      contains

      SUBROUTINE har_mdomain(jj,kk,nconst,mlon,mlat, depths,
     .                  xmjd,ifreq,tname2)
      use m_csrtide
      use m_admit2
      use netcdf
! --- ----------------------------------------------------------
! --- Utilizing files from the data originator to read data
! --- file and estimate amplitudes and phases in all positions
! --- given by mlon,alta.
! --- ------------------------------------------------------- 
cGE       implicit double precision (a-h,o-z)
       implicit none

       integer, intent(in) :: jj,kk,nconst
       real, intent(in) :: xmjd
       real   , intent(in), dimension(jj,kk) :: mlon, mlat, 
     &    depths

       INTEGER, intent(in) :: ifreq(nconst)

       logical pseudo,radiation,isdata
       !real, dimension(jj,kk) :: mlon,mlat
       real amp, pha,h1, h2, x, y, u, v
       dimension u(3,2), v(3,2)
       dimension x(100),y(100),h1(100),h2(100),amp(100),pha(100)
       real*4 ::locamp(jj,kk,nconst)
       real*4 ::locpha(jj,kk,nconst)
       integer :: idmid, jdmid,ncid,varid
       integer :: i,j, iconst,k
       real :: dlon, dlat, tlp,tds,tide

       CHARACTER*4 tname2(nconst) 
       CHARACTER*9 dname(2*nconst+2),tname

       dname(1)=' "Long  "'
       dname(2)=',"Latt  "'
       DO i=1,nconst
        tname=',"A_xxxx"'
        WRITE(tname(5:8),'(A4)')tname2(i)
        dname(2+i)=tname

        tname=',"P_xxxx"'
        WRITE(tname(5:8),'(A4)')tname2(i)
        dname(2+nconst+i)=tname
!        write(*,*)dname(2+i),dname(2+nconst+i)
       ENDDO

! ---
! --- New - added by knut, dump data in netcdf grid format
! ---
       DO j=1,jj
       DO k=1,kk
       if (depths(j,k)>.1) then
         dlon=mlon(j,k)
         dlat=mlat(j,k)
         call csrtide (xmjd,dlat,dlon,tide,tlp,tds,isdata,
     .              u, v, pseudo, radiation)
         call admit2 (pseudo,radiation,u,v,x,y,h1,h2,amp,pha)

         do iconst=1,nconst
            locamp(j,k,iconst) = amp(ifreq(iconst))
            locpha(j,k,iconst) = pha(ifreq(iconst))
         end do
       else
         do iconst=1,nconst
            locamp(j,k,iconst) = -1d10
            locpha(j,k,iconst) = -1d10
         end do
       end if
       ENDDO
       ENDDO

c --- Netcdf Diagnostics
      print *,'Dumping diagnostics to csrtide.nc'
      call ncerr(NF90_create('csrtide.nc',NF90_CLOBBER,ncid))
      call ncerr(NF90_def_dim(ncid,'idm',jj,idmid))
      call ncerr(NF90_def_dim(ncid,'jdm',kk,jdmid))
      call ncerr(NF90_def_var(ncid,'depths',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_put_att(ncid,varid,'_FillValue',
     &   real(0.,kind=4)))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,depths))

      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'lon',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,mlon))

      call ncerr(NF90_redef(ncid))
      call ncerr(NF90_def_var(ncid,'lat',NF90_FLOAT,
     &      (/idmid,jdmid/),varid))
      call ncerr(NF90_enddef(ncid))
      call ncerr(NF90_put_var(ncid,varid,mlat))

      DO iconst=1,nconst
         call ncerr(NF90_redef(ncid))
        call ncerr(NF90_def_var(ncid,trim(tname2(iconst))//'_amplitude',
     &         NF90_FLOAT,(/idmid,jdmid/),varid))
         call ncerr(NF90_put_att(ncid,varid,'_FillValue',
     &      real(-1e10,kind=4)))
         call ncerr(NF90_enddef(ncid))
         call ncerr(NF90_put_var(ncid,varid,locamp(:,:,iconst)))

         call ncerr(NF90_redef(ncid))
         call ncerr(NF90_def_var(ncid,trim(tname2(iconst))//'_phase',
     &         NF90_FLOAT,(/idmid,jdmid/),varid))
         call ncerr(NF90_put_att(ncid,varid,'_FillValue',
     &      real(-1e10,kind=4)))
         call ncerr(NF90_enddef(ncid))
         call ncerr(NF90_put_var(ncid,varid,locpha(:,:,iconst)))
      end do
      call ncerr(NF90_close(ncid))


! ---
! --- ----------------- Dump data in Tecplot Point format
! ---
       stop

       OPEN(10,FILE='tec_tide.dat',STATUS='UNKNOWN')
       WRITE(10,'(''VARIABLES='',18A9)')(dname(i),i=1,18)
       WRITE(10,'(''ZONE I='',I3,'',J='',I3,'',F=POINT'')')
     .       kk,jj

       DO j=1,jj
        DO k=1,kk
         dlon=mlon(j,k)
         dlat=mlat(j,k)
         call csrtide (xmjd,dlat,dlon,tide,tlp,tds,isdata,
     .              u, v, pseudo, radiation)
         call admit2 (pseudo,radiation,u,v,x,y,h1,h2,amp,pha)

         write(10,'(18f8.2)')
     .             dlon,dlat
     .            ,(amp(ifreq(i)),i=1,nconst)
     .            ,(pha(ifreq(i)),i=1,nconst)
        ENDDO
       ENDDO
       CLOSE(10)
       WRITE(*,*)
       WRITE(*,'(''2D data fields are on file obc_tec.dat'')')

      end subroutine
         subroutine ncerr(errcode)
         use netcdf
         implicit none
         integer, intent(in) :: errcode
         if (errcode/=NF90_NOERR) then
            write(6,'(a)') NF90_STRERROR(errcode)
            stop '(ncerr)'
         end if
         end subroutine


      end module
