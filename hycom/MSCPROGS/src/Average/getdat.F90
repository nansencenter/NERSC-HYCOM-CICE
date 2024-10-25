subroutine getdat(flnm,time,iweight,mntype,iexpt,yrflag,thbase)
 use mod_mean  ! HYCOM mean array interface
 use mod_za    ! HYCOM array I/O interface

 character(len=*) :: flnm
 real*8, dimension(3) :: time
 real :: thbase
 !logical, intent(in) :: icegln,trcout
 integer :: iweight,iexpt,mntype,yrflag

 ! --- read model fields and extract portion of global fields.
 ! --- HYCOM 2.0 array I/O archive or mean/mnsq archive file.

 character(len=80) :: cline
 character(len=79), dimension(5) :: preambl
 character(len=6) :: cvarin
 real :: hminb,hmaxb
 real :: thet
 integer :: i,j,k,iversn,l,layer,sigver
 integer :: ios,ktr,ntr
 logical :: nodens

 real, allocatable :: work(:,:)

 data ni/14/

 allocate( work(idm,jdm) )

 l = len_trim(flnm)
 open (unit=ni,file=flnm(1:l-2)//'.b',form='formatted',status='old',action='read')
 call zaiopf(flnm(1:l-2)//'.a','old', ni)

 read( ni,'(a80/a80/a80/a80)') ctitle
 write(lp,'(a80/a80/a80/a80)') ctitle

 read( ni,*) iversn,cvarin
 write(lp,*) cvarin,' = ',iversn
 if (cvarin.ne.'iversn') then
    write(lp,*)
    write(lp,*) 'error in getdat - input ',cvarin,' but should be iversn'
    write(lp,*)
    call flush(lp)
    stop
 endif

 read( ni,*) iexpt,cvarin
 write(lp,*) cvarin,' = ',iexpt
 if (cvarin.ne.'iexpt ') then
    write(lp,*)
    write(lp,*) 'error in getdat - input ',cvarin,' but should be iexpt '
    write(lp,*)
    call flush(lp)
    stop
 endif

 read( ni,*) yrflag,cvarin
 write(lp,*) cvarin,' = ',yrflag
 if (cvarin.ne.'yrflag') then
    write(lp,*)
    write(lp,*) 'error in getdat - input ',cvarin,' but should be yrflag'
    write(lp,*)
    call flush(lp)
    stop
 endif

 read( ni,*) idmtst,cvarin
 write(lp,*) cvarin,' = ',idmtst
 if (cvarin.ne.'idm   ') then
    write(lp,*)
    write(lp,*) 'error in getdat - input ',cvarin,' but should be idm   '
    write(lp,*)
    call flush(lp)
    stop
 endif

 read( ni,*) jdmtst,cvarin
 write(lp,*) cvarin,' = ',jdmtst
 if (cvarin.ne.'jdm   ') then
    write(lp,*)
    write(lp,*) 'error in getdat - input ',cvarin,' but should be jdm   '
    write(lp,*)
    call flush(lp)
    stop
 endif

 if (idmtst.ne.idm .or. jdmtst.ne.jdm) then
    write(lp,*)
    write(lp,*) 'error in getdat - input idm,jdm',' not consistent with parameters'
    write(lp,*) 'idm,jdm = ',idm,   jdm,   '  (REGION.h)'
    write(lp,*) 'idm,jdm = ',idmtst,jdmtst,'  (input)'
    write(lp,*)
    call flush(lp)
    stop
 endif

 read( ni,'(a)') cline
 write(lp,'(a)') cline(1:len_trim(cline))

 if (cline(24:28).eq.'model') then
    iweight = 1  ! standard archive
    mntype  = 0
 elseif (cline(24:28).eq.' mean') then
    iweight = 0  ! mean archive, see below
    mntype  = 1
 elseif (cline(24:28).eq.' mnsq') then
    iweight = 0  ! mean archive, see below
    mntype  = 2
 else
    write(lp,*)
    write(lp,*) 'error in getdat - wrong archive type.'
    write(lp,*) 'only "model", " mean", or " mnsq" allowed.'
    write(lp,*)
    call flush(lp)
    stop
 endif

 !--------------------------------
 ! 'montg1  ' > montg(:,:)
 !--------------------------------
 read (ni,'(a)',end=6) cline
 write(lp,'(a)') cline(1:len_trim(cline))
 i = index(cline,'=')
 read (cline(i+1:),*) nstep,time(1),layer,thet,hminb,hmaxb
 call getfld(work, ni, hminb,hmaxb, .false.)
 call extrct(work,idm,jdm,iorign,jorign,montg,ii,jj)
 write(lp,'("input  ",a," into ",a)') cline(1:8),'montg   '

 ! --- detect version 2.2 normal archive files
 
 nodens = layer.ne.0
 if (nodens) then
    sigver = layer
    thbase = thet
 endif
 !--------------------------------
 ! 'srfhgt  ' > srfht(:,:)
 !--------------------------------
 read(ni,'(a)',end=6) cline
 write(lp,'(a)') cline(1:len_trim(cline))
 i = index(cline,'=')
 read (cline(i+1:),*)  nstep,time(2),layer,thet,hminb,hmaxb
 call getfld(work, ni, hminb,hmaxb, .false.)
 call extrct(work,idm,jdm,iorign,jorign,srfht,ii,jj)
 write(lp,'("input  ",a," into ",a)') cline(1:8),'srfht   '

 if (lsteric) then
    !--------------------------------
    ! 'steric  ' > steric(:,:)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work,ni,hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,steric,ii,jj)
    write(lp,'("input  ",a," into ",a)') cline(1:8),'steric  '
 endif
 !--------------------------------
 ! 'oneta   ' > oneta(:,:)
 !--------------------------------
 if (loneta) then
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,oneta,ii,jj)
    write(lp,'("input  ",a," into ",a)') cline(1:8),'oneta   '
    !-------------------------------- 
    ! 'mnoneta ' > onetaw(:,:)
    !--------------------------------
    if (mntype.eq.2) then  ! mnsq
       read (ni,'(a)',end=6) cline
       write(lp,'(a)')       cline(1:len_trim(cline))
       i = index(cline,'=')
       read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
       call getfld(work,ni, hminb,hmaxb, .false.)
       call extrct(work,idm,jdm,iorign,jorign,onetaw,ii,jj)
       write(lp,'("input  ",a," into ",a)') cline(1:8),'mnoneta '
       if (cline(1:8).ne.'mnoneta ') then
          write(lp,*)
          write(lp,*) 'error in getdat - must have mnoneta in mnsq archive'
          write(lp,*)
          call flush(lp)
          call xcstop('error - must have mnoneta in mnsq archive')
          stop
       endif
    else
       onetaw(:,:) = oneta(:,:)
    endif !mnsq:else
 else
   oneta(:,:) = 1.0
   onetaw(:,:) = 1.0
 endif
 !--------------------------------
 ! 'surflx  ' > surflx(:,:)
 !--------------------------------
 read (ni,'(a)',end=6) cline
 write(lp,'(a)')       cline(1:len_trim(cline))
 i = index(cline,'=')
 read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
 call getfld(work, ni, hminb,hmaxb, .false.)
 call extrct(work,idm,jdm,iorign,jorign,surflx,ii,jj)
 write(lp,'("input  ",a," into ",a)') cline(1:8),'surflx  '
 !--------------------------------
 ! 'wtrflx  ' > wtrflx(:,:)
 !--------------------------------
 if (lwtrflx) then
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,wtrflx,ii,jj)
    write(lp,'("input  ",a," into ",a)') cline(1:8),'wtrflx  '
 else
    wtrflx(:,:) = 0.0
 endif
 !--------------------------------
 ! 'salflx  ' > salflx(:,:)
 !--------------------------------
 read (ni,'(a)',end=6) cline
 write(lp,'(a)')       cline(1:len_trim(cline))
 i = index(cline,'=')
 read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
 call getfld(work, ni, hminb,hmaxb, .false.)
 call extrct(work,idm,jdm,iorign,jorign,salflx,ii,jj)
 write(lp,'("input  ",a," into ",a)') cline(1:8),'salflx  '
 !--------------------------------
 ! 'bl_dpth ' >  dpbl(:,:)
 !--------------------------------
 read (ni,'(a)',end=6) cline
 write(lp,'(a)')       cline(1:len_trim(cline))
 i = index(cline,'=')
 read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
 call getfld(work, ni, hminb,hmaxb, .false.)
 call extrct(work,idm,jdm,iorign,jorign,dpbl,ii,jj)
 write(lp,'("input  ",a," into ",a)') cline(1:8),'dpbl    '
 !--------------------------------
 ! 'mix_dpth' > dpmixl(:,:)
 !--------------------------------
 read (ni,'(a)',end=6) cline
 write(lp,'(a)')       cline(1:len_trim(cline))
 i = index(cline,'=')
 read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
 call getfld(work, ni, hminb,hmaxb, .false.)
 call extrct(work,idm,jdm,iorign,jorign,dpmixl,ii,jj)
 write(lp,'("input  ",a," into ",a)') cline(1:8),'dpmixl  '

 if (.not. nodens) then
    !--------------------------------
    ! 'tmix    ' > tmix(:,:)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,tmix,ii,jj)
    write(lp,'("input  ",a," into ",a)') cline(1:8),'tmix    '
    !--------------------------------
    ! 'smix    ' > smix(:,:)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,smix,ii,jj)
    write(lp,'("input  ",a," into ",a)') cline(1:8),'smix    '
    !--------------------------------
    ! 'thmix   ' >  thmix(:,:)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    thbase = thet
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,thmix,ii,jj)
    write(lp,'("input  ",a," into ",a)') cline(1:8),'thmix   '
    !--------------------------------
    ! 'umix    ' > umix(:,:)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .true.)
    call extrct(work,idm,jdm,iorign,jorign,umix,ii,jj)
    write(lp,'("input  ",a," into ",a)') cline(1:8),'umix    '
    !--------------------------------
    ! 'vmix    ' > vmix(:,:)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .true.)
    call extrct(work,idm,jdm,iorign,jorign,vmix,ii,jj)
    write(lp,'("input  ",a," into ",a)') cline(1:8),'vmix    '
 endif !.not. nodens

 !--------------------------------
 ! 'kemix   ' > kemix(:,:)
 !--------------------------------
 if (iweight.eq.0) then  ! mean archive
    read (ni,'(a)',end=6) cline
    write(lp,'(a)') cline(1:len_trim(cline)) 
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,kemix,ii,jj)
    write(lp,'("input  ",a," into ",a)') cline(1:8),'kemix   '
 endif

 ! --- is there ice?
 if (icegln) then
    !--------------------------------
    ! 'covice  ' > covice(:,:)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,covice,ii,jj)
    !--------------------------------
    ! 'thcice  ' > thkice(:,:)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)') cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,thkice,ii,jj)
    !--------------------------------
    ! 'temice  ' > temice(:,:)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,temice,ii,jj)
 else
    covice(:,:) = 0.0
    thkice(:,:) = 0.0
    temice(:,:) = 0.0
 endif
 !--------------------------------
 ! 'u_btrop ' > ubaro(:,:)
 !--------------------------------
 read (ni,'(a)',end=6) cline
 write(lp,'(a)')       cline(1:len_trim(cline))
 i = index(cline,'=')
 read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
 call getfld(work, ni, hminb,hmaxb, .true.)
 call extrct(work,idm,jdm,iorign,jorign,ubaro,ii,jj)
 write(lp,'("input  ",a," into ",a)') cline(1:8),'ubaro   '
 !--------------------------------
 ! 'v_btrop ' > vbaro(:,:)
 !--------------------------------
 read (ni,'(a)',end=6) cline
 write(lp,'(a)')       cline(1:len_trim(cline))
 i = index(cline,'=')
 read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
 call getfld(work, ni, hminb,hmaxb, .true.)
 call extrct(work,idm,jdm,iorign,jorign,vbaro,ii,jj)
 write(lp,'("input  ",a," into ",a)') cline(1:8),'vbaro   '
 !--------------------------------
 ! 'kebtrop ' > kebaro(:,:)
 !--------------------------------
 if (iweight.eq.0) then  ! mean archive
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .true.)
    call extrct(work,idm,jdm,iorign,jorign,kebaro,ii,jj)
    write(lp,'("input  ",a," into ",a)') cline(1:8),'kebaro  '
 endif

 do k=1,kk
    !--------------------------------
    ! 'u-vel   ' > u(:,:,k)
    !--------------------------------
    read (ni,'(a)',end=6)   cline
    write(lp,'(a)')         cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .true.)
    call extrct(work,idm,jdm,iorign,jorign,u(1,1,k),ii,jj)
    write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'u       ',k
    if (cline(1:8).ne.'u-vel.  ') then
       write(lp,*)
       write(lp,*) 'error in getdat - layer ',k,' does not exist (kk= ',kk,')'
       write(lp,*)
       call flush(lp)
       stop
    endif
    !--------------------------------
    ! 'v-vel   ' > v(:,:,k)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)') cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .true.)
    call extrct(work,idm,jdm,iorign,jorign,v(1,1,k),ii,jj)
    write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'v       ',k
    !--------------------------------
    ! 'k.e.    ' > ke(:,:,k)
    !--------------------------------
    if (iweight.eq.0) then  ! mean archive
       read (ni,'(a)',end=6) cline
       write(lp,'(a)')       cline(1:len_trim(cline))
       i = index(cline,'=')
       read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
       call getfld(work, ni, hminb,hmaxb, .false.)
       call extrct(work,idm,jdm,iorign,jorign,ke(1,1,k),ii,jj)
       write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'k.e.    ',k
    endif
    !--------------------------------
    ! 'thknss  ' > dp(:,:,k)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)') cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,dp(1,1,k),ii,jj)
    write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'dp      ',k
    !--------------------------------
    ! 'mnthknss' > dw(:,:,k)
    !--------------------------------
    if (mntype.eq.2) then  ! mnsq
       read (ni,'(a)',end=6) cline
       write(lp,'(a)')       cline(1:len_trim(cline))
       i = index(cline,'=')
       read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
       call getfld(work, ni, hminb,hmaxb, .false.)
       call extrct(work,idm,jdm,iorign,jorign,dw(1,1,k),ii,jj)
       write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'dw      ',k
       if (cline(1:8).ne.'mnthknss') then
          write(lp,*)
          write(lp,*) 'error in getdat - must have mnthknss in mnsq archive'
          write(lp,*)
          call flush(lp)
          stop
       endif
    else
       dw(:,:,k) = dp(:,:,k)
    endif
    !--------------------------------
    ! 'temp    ' > temp(:,:,k)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,temp(1,1,k),ii,jj)
    write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'temp    ',k
    !--------------------------------
    ! 'salin   ' > saln(:,:,k)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)') cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,saln(1,1,k),ii,jj)
    write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'saln    ',k
    !--------------------------------
    ! 'density ' > th3d(:,:,k)
    !--------------------------------
    if (.not. nodens) then
       read (ni,'(a)',end=6) cline
       write(lp,'(a)')       cline(1:len_trim(cline))
       i = index(cline,'=')
       read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
       call getfld(work, ni, hminb,hmaxb, .false.)
       call extrct(work,idm,jdm,iorign,jorign,th3d(1,1,k),ii,jj)
       write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'th3d    ',k
    else
       call th3d_p(temp(1,1,k),saln(1,1,k),th3d(1,1,k),ii,jj,sigver,thbase)
       write(lp,'("    ",a8,"calculate ",a,i3)') " ",'th3d    ',k
       if (k.eq.1) then
           tmix(:,:) = temp(:,:,1)
           smix(:,:) = saln(:,:,1)
          thmix(:,:) = th3d(:,:,1)
           umix(:,:) =    u(:,:,1)
           vmix(:,:) =    v(:,:,1)
          write(lp,'("copy   ",a," into ",a)') 'temp.1  ','tmix    '
          write(lp,'("copy   ",a," into ",a)') 'saln.1  ','smix    '
          write(lp,'("copy   ",a," into ",a)') 'th3d.1  ','thmix   '
          write(lp,'("copy   ",a," into ",a)') '   u.1  ','umix    '
          write(lp,'("copy   ",a," into ",a)') '   v.1  ','vmix    '
       endif !k==1
    endif !.not.nodens:else
    !--------------------------------
    ! 'tracer  ' > tracer(:,:,k,ktr)
    !--------------------------------
    if (trcout) then
       do ktr= 1,ntracr
          read (ni,'(a)',end=6) cline
          write(lp,'(a)')       cline(1:len_trim(cline))
          i = index(cline,'=')
          read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
          call getfld(work, ni, hminb,hmaxb, .false.)
          call extrct(work,idm,jdm,iorign,jorign,tracer(1,1,k,ktr),ii,jj)
          write(lp,'("input  ",a," into ",a,2i3)') cline(1:8),'tracer  ',k,ktr
       enddo !ktr
    endif !trcout
    !--------------------------------
    ! 'viscty  ' > visc(:,:,k)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)') cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,visc(1,1,k),ii,jj)
    write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'visc    ',k
    !--------------------------------
    ! 't-diff  ' > tdff(:,:,k)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)') cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,tdff(1,1,k),ii,jj)
    write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'tdff    ',k
    !--------------------------------
    ! 's-diff  ' > sdff(:,:,k)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)') cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,sdff(1,1,k),ii,jj)
    write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'sdff    ',k

    write(lp,'(a,f9.5)') 'finished reading data for layer',thet
    call flush(lp)
    theta(k)=thet

 enddo !k-loop

 if (bgcout) then
    !--------------------------------
    ! > tracer_bgc_3d
    !--------------------------------
    do k=1,kk
       do ktr= 1,ntracr_bgc_3d
          read (ni,'(a)',end=6) cline
          write(lp,'(a)')       cline(1:len_trim(cline))
          i = index(cline,'=')
          read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
          call getfld(work, ni, hminb,hmaxb, .false.)
          call extrct(work,idm,jdm,iorign,jorign,tracer_bgc_3d(1,1,k,ktr),ii,jj)
          write(lp,'("input  ",a," into ",a,2i3)') cline(1:8),'tracer_bgc_3d',k,ktr
       enddo !ktr
    enddo !k-loop
    !--------------------------------
    ! > tracer_bgc_2d
    !--------------------------------
    do ktr= 1,ntracr_bgc_2d
       read (ni,'(a)',end=6) cline
       write(lp,'(a)')       cline(1:len_trim(cline))
       i = index(cline,'=')
       read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
       call getfld(work, ni, hminb,hmaxb, .false.)
       call extrct(work,idm,jdm,iorign,jorign,tracer_bgc_2d(1,1,ktr),ii,jj)
       write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'tracer_bgc_2d'
    enddo !ktr
 endif ! bgcout
    
 if (mntype.eq.1.and.icegln) then
    !--------------------------------
    ! 'si_u    ' > si_u(:,:)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,si_u,ii,jj)
    write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'si_u'
    !--------------------------------
    ! 'si_v    ' > si_v(:,:)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,si_v,ii,jj)
    write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'si_v'
    !--------------------------------
    ! 'surtx   ' > surtx(:,:)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,surtx,ii,jj)
    write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'surtx'
    !--------------------------------
    ! 'surty   ' > surty(:,:)
    !--------------------------------
    read (ni,'(a)',end=6) cline
    write(lp,'(a)')       cline(1:len_trim(cline))
    i = index(cline,'=')
    read (cline(i+1:),*)  nstep,time(3),layer,thet,hminb,hmaxb
    call getfld(work, ni, hminb,hmaxb, .false.)
    call extrct(work,idm,jdm,iorign,jorign,surty,ii,jj)
    write(lp,'("input  ",a," into ",a,i3)') cline(1:8),'surty'
 endif

 if (iweight.eq.0) then
    iweight = nstep  ! mean archive
 endif

 close( unit=ni)
 call zaiocl(ni)

 deallocate( work )

 return

 ! --- unexpected end of file
6 continue
 write (lp,*) '***** unexpected end of archive file *****'
 call flush(lp)
 stop '(e-o-f)'

end subroutine getdat
 
subroutine getdepth(dpthfil)
 use mod_mean  ! HYCOM mean array interface
 use mod_za

 character :: dpthfil*(*)

 ! --- acquire basin depths

 character :: cline*80
 character(len=79), dimension(5) :: preambl
 character(len=6) :: cvarin
 real :: hminb,hmaxb
 integer :: i,j,iversn,l

 real, allocatable :: work(:,:)

 allocate( work(idm,jdm) )

 open(unit=9,file=dpthfil(1:len_trim(dpthfil))//'.b',form='formatted',status='old',action='read')
 read(9, '(a79)') preambl
 write(lp,'(a79)') preambl
 read (9, '(a)')   cline
 write(lp,'(a)')   cline(1:len_trim(cline))
 i = index(cline,'=')
 read (cline(i+1:),*)   hminb,hmaxb
 close(unit=9)

 call zaiopf(dpthfil(1:len_trim(dpthfil))//'.a','old', 9)
 call getfld(work, 9, hminb,hmaxb, .true.)
 call zaiocl(9)
 write(lp,'("read ",a," into ",a)') preambl(1)(1:8),'depths  '
 do j= 1,jj
 do i= 1,ii
            depths(i,j) = work(i,j)
 enddo
 enddo
 do i= 0,ii
          depths(i,0) = 0.0
 enddo
 do j= 0,jj
          depths(0,j) = depths(ii,j)  ! assumed periodic
 enddo

 deallocate( work )

 return

end subroutine getdepth

subroutine getfld(work, iunit, hminb,hmaxb, lzero)
 use mod_za ! HYCOM array I/O interface

 ! --- read a single array

 logical :: lzero
 integer :: iunit
 real, dimension(idm,jdm) :: work
 integer, dimension(idm,jdm) :: mask
 real    :: hminb,hmaxb
 real    :: hmina,hmaxa

 call zaiord(work,mask,.false., hmina,hmaxa, iunit)

 if (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or. abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4 ) then
    write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')      &
            ,'error - .a and .b files not consistent:' &
            ,'.a,.b min = ',hmina,hminb,hmina-hminb    &
            ,'.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb 
    call flush(lp)
    stop
 endif

 if (lzero) then
    do j= 1,jdm
    do i= 1,idm
       if (work(i,j).gt.2.0**99) then
          work(i,j) = 0.0
       endif
    enddo
    enddo
 endif

 return
end subroutine getfld

subroutine th3d_p(temp,saln,th3d,no,mo,sigver,thbase)
 implicit none

 integer no,mo,sigver
 real, dimension(no,mo) :: temp,saln,th3d
 real :: thbase

 ! --- calculate density using the appropriate equation of state.

 if (sigver.eq.1) then
    call th3d_p1(temp,saln,th3d,no,mo,thbase)
 elseif (sigver.eq.2) then
    call th3d_p2(temp,saln,th3d,no,mo,thbase)
 elseif (sigver.eq.3) then
    call th3d_p3(temp,saln,th3d,no,mo,thbase)
 elseif (sigver.eq.4) then
    call th3d_p4(temp,saln,th3d,no,mo,thbase)
 elseif (sigver.eq.5) then
    call th3d_p5(temp,saln,th3d,no,mo,thbase)
 elseif (sigver.eq.6) then
    call th3d_p6(temp,saln,th3d,no,mo,thbase)
 elseif (sigver.eq.7) then
    call th3d_p7(temp,saln,th3d,no,mo,thbase)
 elseif (sigver.eq.8) then
    call th3d_p8(temp,saln,th3d,no,mo,thbase)
 else  !unknown
        th3d(:,:) = 0.0
 endif
 return

end subroutine th3d_p

subroutine th3d_p1(temp,saln,th3d,no,mo,thbase)
 implicit none

 integer :: no,mo
 real, dimension(no,mo) :: temp,saln,th3d
 real :: thbase

 ! --- calculate density using the equation of state.

 integer :: i,j

 include './include/stmt_fns_SIGMA0_7term.h90'

 do j= 1,mo
 do i= 1,no
    th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
 enddo !i
 enddo !j
 return

end subroutine th3d_p1

subroutine th3d_p2(temp,saln,th3d,no,mo,thbase)
 implicit none

 integer :: no,mo
 real, dimension(no,mo) :: temp,saln,th3d
 real :: thbase

 ! --- calculate density using the equation of state.

 integer :: i,j

 include './include/stmt_fns_SIGMA2_7term.h90'

 do j= 1,mo
 do i= 1,no
    th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
 enddo !i
 enddo !j
 return
end subroutine th3d_p2

subroutine th3d_p3(temp,saln,th3d,no,mo,thbase)
 implicit none

 integer :: no,mo
 real, dimension(no,mo) :: temp,saln,th3d
 real    :: thbase

 ! --- calculate density using the equation of state.

 integer :: i,j

 include './include/stmt_fns_SIGMA0_9term.h90'

 do j= 1,mo
 do i= 1,no
    th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
 enddo !i
 enddo !j
 return
end subroutine th3d_p3

subroutine th3d_p4(temp,saln,th3d,no,mo,thbase)
 implicit none

 integer :: no,mo
 real, dimension(no,mo) :: temp,saln,th3d
 real    :: thbase

 ! --- calculate density using the equation of state.

 integer :: i,j

 include './include/stmt_fns_SIGMA2_9term.h90'

 do j= 1,mo
 do i= 1,no
    th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
 enddo !i
 enddo !j
 return
end subroutine th3d_p4

subroutine th3d_p5(temp,saln,th3d,no,mo,thbase)
 implicit none

 integer :: no,mo
 real, dimension(no,mo) :: temp,saln,th3d
 real    :: thbase

 ! --- calculate density using the equation of state.

 integer :: i,j

 include './include/stmt_fns_SIGMA0_17term.h90'

 do j= 1,mo
 do i= 1,no
    th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
 enddo !i
 enddo !j
 return

end subroutine th3d_p5

subroutine th3d_p6(temp,saln,th3d,no,mo,thbase)
 implicit none

 integer :: no,mo
 real, dimension(no,mo) :: temp,saln,th3d
 real    :: thbase

 ! --- calculate density using the equation of state.

 integer :: i,j

 include './include/stmt_fns_SIGMA2_17term.h90'

 do j= 1,mo
 do i= 1,no
    th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
 enddo !i
 enddo !j
 return
end subroutine th3d_p6

subroutine th3d_p7(temp,saln,th3d,no,mo,thbase)
 implicit none

 integer :: no,mo
 real, dimension(no,mo) :: temp,saln,th3d
 real    :: thbase

 ! --- calculate density using the equation of state.

 integer :: i,j

 include './include/stmt_fns_SIGMA0_12term.h90'

 do j= 1,mo
 do i= 1,no
    th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
 enddo !i
 enddo !j

 return
end subroutine th3d_p7

subroutine th3d_p8(temp,saln,th3d,no,mo,thbase)
 implicit none

 integer :: no,mo
 real, dimension(no,mo) :: temp,saln,th3d
 real    :: thbase

 ! --- calculate density using the equation of state.

 integer :: i,j

 include './include/stmt_fns_SIGMA2_12term.h90'

 do j= 1,mo
 do i= 1,no
    th3d(i,j) = sig(r8(temp(i,j)),r8(saln(i,j))) - thbase
 enddo !i
 enddo !j

 return
end subroutine th3d_p8
