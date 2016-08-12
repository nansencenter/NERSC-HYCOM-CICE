      program isuba_field
      use mod_za  ! HYCOM I/O interface
      use mod_zb  ! HYCOM I/O interface for subregion
      implicit none
c
c     create a diferent-grid subregion from a full region hycom file.
c
c     subregion grid is arbitrary, except that any part of this grid
c     that is outside the input grid must be land.
c
c     the existance of file regional.grid.masked is a signal to
c     only consider ocean points on the original grid.  note that
c     regional.grid.masked can be empty, since it is not read in.
c
c     Alan J. Wallcraft,  NRL,  January 2005.
c
      character*80         :: cline,cline_u,cline_out
      character*256        :: flnm_in, flnm_tin,
     &                        flnm_out,flnm_top,flnm_reg
      integer              :: idm_out,jdm_out,khead,kskip,
     &                        iref_out,jref_out,iref_in,jref_in,
     &                        ijgrd
      integer              :: i,ii,ijsafe,ios,ip,itmp,ix,j,jj,jp,jq,l,
     &                        ni,nir,no,nor
      integer              :: k,l0,l1
      integer              :: if_sm,il_sm,jf_sm,jl_sm
      logical              :: laxis(2),laxis2,
     &                        lperiod(2),lcycle,lmask
      real                 :: hmina,hminb,hmaxa,hmaxb,
     &                        deg2rad,dist,distj,dist_max,
     &                        dx,dy,qdx,xp,yp,plat_min,plat_max
      integer, allocatable :: m_sm(:,:),  iv_sm(:,:)
      integer, allocatable :: m_in(:,:),  m_out(:,:),   m_osm(:,:)
      real,    allocatable :: a_in(:,:),  a_out(:,:)
c
      real,    allocatable :: plat_in( :,:),plon_in( :,:)
      real,    allocatable :: plat_out(:,:),plon_out(:,:)
      integer, allocatable :: i_out(:,:),j_out(:,:)
      real,    allocatable :: x_out(:,:),y_out(:,:)
      real,    allocatable :: plat_in_min(:),plat_in_max(:)
c
      real,    parameter   :: hspval=0.5*2.0**100  ! half spval
c
      real*8         ztecnf
      external       ztecnf,ztecng,ztecnp,ztecnb
c
      integer        its
      real*8         acc,err,step
      real*8         x2(2),w(6)
c
      logical        ldebug
      real*8         a,b
      common/zaecnb/ a(0:2,0:2),b(0:2,0:2),ldebug
      save  /zaecnb/
c
      call xcspmd
      call zaiost
      call blkdat(idm_out,jdm_out, khead,kskip,
     &            flnm_reg,flnm_in,flnm_tin,flnm_out,flnm_top)
      call zbiost(idm_out,jdm_out)
c
      allocate(  iv_sm(jdm,2) )
c
      allocate(   m_sm(idm,jdm),     m_osm(idm_out,jdm_out) )
      allocate(   m_in(idm,jdm),     m_out(idm_out,jdm_out) )
      allocate(   a_in(idm,jdm),     a_out(idm_out,jdm_out) )
      allocate( plat_in(idm,jdm), plat_out(idm_out,jdm_out) )
      allocate( plon_in(idm,jdm), plon_out(idm_out,jdm_out) )
c                                                            
      allocate( plat_in_max(jdm), plat_in_min(jdm) )         
c                                                   
      allocate( i_out(idm_out,jdm_out), j_out(idm_out,jdm_out) )
      allocate( x_out(idm_out,jdm_out), y_out(idm_out,jdm_out) )
c
c     get the output p-grid mask from the bathymetry.
c
      if     (flnm_top.eq.'NONE') then
        write(6,'(/a)') 'NO OUTPUT BATHYMETRY'
        a_out(:,:) = 100.0
      else
        l  = len_trim(flnm_top)
        open (unit=13,file=flnm_top(1:l-2)//'.b',form='formatted',
     .        status='old',action='read')
        write(6,'(/a)') 'OUTPUT BATHYMETRY:'
        do i= 1,6
          read(13,'(a)') cline
          write(6,'(a)') cline(1:len_trim(cline))
        enddo
        l = index(cline,'=')
        read (cline(l+1:),*)  hminb,hmaxb
        close(unit=13)
c
        l  = len_trim(flnm_top)
        call zbiopf(flnm_top(1:l-2)//'.a','old', 13)
        call zbiord(a_out,m_out,.false., hmina,hmaxa, 13)
        call zbiocl(13)
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - output bathymetry .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          stop
        endif
      endif !flnm_top
c
      do j= 1,jdm_out
        do i= 1,idm_out
          if     (a_out(i,j).gt.hspval .or.
     &            a_out(i,j).le.0.0        ) then
            m_out(i,j)   = 0
          else
            m_out(i,j)   = 1
          endif
        enddo
      enddo
c
c     get the input p-grid mask from the bathymetry.
c
      if     (flnm_tin.eq.'NONE') then
        write(6,'(/a)') 'NO INPUT BATHYMETRY'
        a_in(:,:) = 100.0
      else
        l  = len_trim(flnm_tin)
        open (unit=13,file=flnm_tin(1:l-2)//'.b',form='formatted',
     .        status='old',action='read')
        write(6,'(/a)') ' INPUT BATHYMETRY:'
        do i= 1,6
          read(13,'(a)') cline
          write(6,'(a)') cline(1:len_trim(cline))
        enddo
        l = index(cline,'=')
        read (cline(l+1:),*)  hminb,hmaxb
        close(unit=13)
c
        l  = len_trim(flnm_tin)
        call zaiopf(flnm_tin(1:l-2)//'.a','old', 13)
        call zaiord(a_in,m_in,.false., hmina,hmaxa, 13)
        call zaiocl(13)
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          write(6,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &      'error - bathymetry .a and .b files not consistent:',
     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          stop
        endif
      endif
c
      do j= 1,jdm
        do i= 1,idm
          if     (a_in(i,j).gt.hspval .or.
     &            a_in(i,j).le.0.0        ) then
            m_in(i,j) = 0
            m_sm(i,j) = 0
          else
            m_in(i,j) = 1
            m_sm(i,j) = 1
          endif
        enddo
      enddo
c
c     should we land-mask the input grid?
c
      inquire(file='regional.grid.masked', exist=lmask)
      if     (lmask) then
        write(6,'(/a/)') ' input domain lon,lat is landmasked'
      else
        write(6,'(/a/)') ' input domain lon,lat is not landmasked'
      endif
c       
c     read the input and output grid locations (no error checking).
c
      nir = 24
      call zaiopf('regional.grid.a','old', nir)
      call zaiord(plon_in,m_in,.false., hmina,hmaxa, nir)
      call zaiord(plat_in,m_in,.false., hmina,hmaxa, nir)
      call zaiocl(nir)
c     
      nor = 25
      call zbiopf(trim(flnm_reg),'old', nor)
      call zbiord(plon_out,m_out,.false., hmina,hmaxa, nor)
      call zbiord(plat_out,m_out,.false., hmina,hmaxa, nor)
      call zbiocl(nor)
c     
c     calculate the output grid location w.r.t. the input grid.
c
      lperiod(1) = maxval(plon_in(:,:))-
     &             minval(plon_in(:,:))  .gt. 350.0
      do j= 1,jdm
        if     (.not.lperiod(1)) then
          exit
        endif
        dy = mod( abs(plon_in(1,j) - plon_in(  3,j)), 360.0 )
        if     (dy.gt.180.0) then
          dy = 360.0 - dy  !abs distance
        endif
        dx = mod( abs(plon_in(1,j) - plon_in(idm,j)), 360.0 )
        if     (dx.gt.180.0) then
          dx = 360.0 - dx  !abs distance
        endif
        lperiod(1) = lperiod(1) .and. dx.lt.dy  !1 and idm closer than 1 and 3
      enddo
      if     (lperiod(1)) then
        write(6,'(a)') ' input domain assumed to be periodic'
      else
        write(6,'(a)') ' input domain assumed to be non-periodic'
      endif
c
      lperiod(2) = maxval(plon_out(:,:))-
     &             minval(plon_out(:,:))   .gt. 350.0
      do j= 1,jdm_out
        if     (.not.lperiod(2)) then
          exit
        endif
        dy = mod( abs(plon_out(1,j) - plon_out(      3,j)), 360.0 )
        if     (dy.gt.180.0) then
          dy = 360.0 - dy  !abs distance
        endif
        dx = mod( abs(plon_out(1,j) - plon_out(idm_out,j)), 360.0 )
        if     (dx.gt.180.0) then
          dx = 360.0 - dx  !abs distance
        endif
        lperiod(2) = lperiod(2) .and. dx.lt.dy  !1 and idm closer than 1 and 3
      enddo
      if     (lperiod(2)) then
        write(6,'(a)') 'output domain assumed to be periodic'
      else
        write(6,'(a)') 'output domain assumed to be non-periodic'
      endif
c       
      laxis(1) = .true.
      do i= 2,idm
        laxis(1) = laxis(1) .and.
     &             maxval(abs(plat_in(1,:)-plat_in(i,:))).le.1.e-2
        if     (.not.laxis(1)) then
          exit
        endif
      enddo
      do j= 2,jdm
        laxis(1) = laxis(1) .and.
     &             maxval(abs(plon_in(:,1)-plon_in(:,j))).le.1.e-2
        if     (.not.laxis(1)) then
          exit
        endif
      enddo
      if     (laxis(1)) then
        write(6,'(a)') ' input domain has 1-d lat/lon axes'
      else
        write(6,'(a)') ' input domain is curvilinear'
      endif
c
      laxis(2) = .true.
      do i= 2,idm_out
        laxis(2) = laxis(2) .and.
     &             maxval(abs(plat_out(1,:)-plat_out(i,:))).le.1.e-2
        if     (.not.laxis(2)) then
          exit
        endif
      enddo
      do j= 2,jdm_out
        laxis(2) = laxis(2) .and.
     &             maxval(abs(plon_out(:,1)-plon_out(:,j))).le.1.e-2
        if     (.not.laxis(2)) then
          exit
        endif
      enddo
      if     (laxis(2)) then
        write(6,'(a)') 'output domain has 1-d lat/lon axes'
      else
        write(6,'(a)') 'output domain is curvilinear'
      endif
      write(6,*)
      call flush(6)
c     
      laxis2 = laxis(1) .and. laxis(2)
c     
      do j= 1,jdm
        plat_in_min(j) = minval(plat_in(:,j))
        plat_in_max(j) = maxval(plat_in(:,j))
      enddo
      plat_min = minval(plat_in_min(:))
      plat_max = maxval(plat_in_max(:))
      dist_max = 0.0
      do j= 1,jdm-1
        do i= 1,idm-1
          dist_max = max( abs(plat_in(i,j)-plat_in(i+1,j)),
     &                    abs(plat_in(i,j)-plat_in(i,j+1)),
     &                    dist_max )
        enddo
      enddo
      dist_max = 2*dist_max  !target must be at least this close in latitude
c
      deg2rad = 4.d0*atan(1.d0)/180.d0  !pi/180
c       
      ldebug = .false.
c     
      do jj= 1,jdm_out
        do ii= 1,idm_out
          if     (ii.eq.1 .or. jj.eq.1 .or.  m_out(ii,jj).eq.1) then
          ldebug = mod(ii,idm_out/10).eq.1 .and.
     &             mod(jj,jdm_out/10).eq.1
          if     (laxis2 .and. ii.ne.1 .and. jj.ne.1) then
c         
c           shortcut for 1-d axes.
c
            i_out(ii,jj) = i_out(ii, 1)
            j_out(ii,jj) = j_out( 1,jj)
            x_out(ii,jj) = x_out(ii, 1)
            y_out(ii,jj) = y_out( 1,jj)
            cycle
          endif
c           
c         find the nearest point by exhaustive search, but improve
c         efficiency by using plat_*_min/_max to exclude far away rows.
c
          xp   = plon_out(ii,jj)
          yp   = plat_out(ii,jj)
          yp   = min(plat_max,max(plat_min,yp))  !in the input latitude range
          qdx  = max(0.001,abs(cos(yp*deg2rad)))
c         
c         start with a nearby point.
c
          if     (ii.eq.1) then
            if     (jj.eq.1) then
              ip = 1
              jp = 1
            else
              ip = i_out(1,jj-1)
              jp = j_out(1,jj-1)
            endif
          endif
          dy =      abs(plat_in(ip,jp) - yp)
          dx = mod( abs(plon_in(ip,jp) - xp), 360.0 )
          if     (dx.gt.180.0) then
            dx = 360.0 - dx
          endif
          dist = qdx*dx+dy
          if     (ldebug) then
            write(6,'(a,4i5,3f9.2)')
     &        'ii,jj,ip,jp,dx,dy,dist = ',
     &         ii,jj,ip,jp,dx,dy,dist
            call flush(6)
          endif
c           
          lcycle = .false.
          do jq= 0,jdm
            if     (jq.eq.0) then
              j = jp  ! search estimated row location first
            else
              j = jq
            endif
            distj = min(dist,dist_max)
            if     (.not. ldebug) then
              if     (yp.lt.plat_in_min(j)-distj .or.
     &                yp.gt.plat_in_max(j)+distj     ) then
                cycle  ! far away row
              endif
            else !debug
              if     (yp.lt.plat_in_min(j)-distj .or.
     &                yp.gt.plat_in_max(j)+distj     ) then
                if     (.not.lcycle) then
                  write(6,'(a,2i5,5x,i5,f9.2)')
     &              'ii,jj,j,dist (cycle-strt)',
     &               ii,jj,j,dist
                  call flush(6)
                elseif (jq.eq.jdm) then
                  write(6,'(a,2i5,5x,i5,f9.2)')
     &              'ii,jj,j,dist (cycle-stop)',
     &               ii,jj,j,dist
                  call flush(6)
                endif
                lcycle = .true.
                cycle  ! far away row
              else
                if     (lcycle) then
                  write(6,'(a,2i5,5x,i5,f9.2)')
     &              'ii,jj,j,dist (cycle-stop)',
     &               ii,jj,j-1,dist
                  call flush(6)
                endif
                lcycle = .false.
              endif
            endif !.not.ldebug;else
            if     (dist.eq.0.0) then
              exit   ! found exact location
            endif
            do i = 1,idm
              if     (lmask .and. m_in(i,j).eq.0) then  !land point
                cycle
              endif
              dy =      abs(plat_in(i,j) - yp)
              dx = mod( abs(plon_in(i,j) - xp), 360.0 )
              if     (dx.gt.180.0) then
                dx = 360.0 - dx
              endif
              if     (qdx*dx+dy.lt.dist) then
                ip   = i
                jp   = j
                dist = qdx*dx+dy
                if     (ldebug) then
                  write(6,'(a,4i5,3f9.2)')
     &              'ii,jj,ip,jp,dx,dy,dist = ',
     &               ii,jj,ip,jp,dx,dy,dist
                  call flush(6)
                endif
              endif
            enddo !i
          enddo !j
c           
c         convert nearest point into bilinear cell and distances.
c
          if     (dist.eq.0.0) then  !exact location
            if     (ip.eq.idm) then
              i_out(ii,jj) = ip-1
              x_out(ii,jj) = 1.0
            else
              i_out(ii,jj) = ip
              x_out(ii,jj) = 0.0
            endif
            if     (jp.eq.jdm) then
              j_out(ii,jj) = jp-1
              y_out(ii,jj) = 1.0
            else
              j_out(ii,jj) = jp
              y_out(ii,jj) = 0.0
            endif
          elseif (dist.gt.dist_max) then !outside grid, use nearest point
            if     (ip.eq.idm) then
              i_out(ii,jj) = ip-1
              x_out(ii,jj) = 1.0
            else
              i_out(ii,jj) = ip
              x_out(ii,jj) = 0.0
            endif
            if     (jp.eq.jdm) then
              j_out(ii,jj) = jp-1
              y_out(ii,jj) = 1.0
            else
              j_out(ii,jj) = jp
              y_out(ii,jj) = 0.0
            endif
          else  !standard case
c       
c           find exact location with napack routine(s).
c           over-kill for rectilinear, but neccessary for curvilinear grids.
c
            if     (ip.eq.1   .and. .not.lperiod(1)) then
              ip = 2
            elseif (ip.eq.idm .and. .not.lperiod(1)) then
              ip = idm-1
            endif
            if     (jp.eq.1) then
              jp = 2
            elseif (jp.eq.jdm) then
              jp = jdm-1
            endif
            do j= 0,2
              do i= 0,2
                ix = ip+i-1
                if     (lperiod(1)) then
                  if     (ix.eq.0) then
                    ix = idm
                  elseif (ix.eq.idm+1) then
                    ix = 1
                  endif
                endif
                b(i,j) =      plat_in(ix,jp+j-1) - yp
                a(i,j) = mod( plon_in(ix,jp+j-1) - xp, 360.0 )
                if     (a(i,j).lt.-180.0) then
                  a(i,j) = 360.0 + a(i,j)
                elseif (a(i,j).gt. 180.0) then
                  a(i,j) = a(i,j) - 360.0
                endif
                a(i,j) = qdx*a(i,j)
              enddo !i
            enddo !j
            if     (b(0,1).eq.b(1,1) .and. b(1,1).eq.b(2,1)) then !rectilinear
              x2(1) = ip - a(1,1)/(a(2,1)-a(1,1))
              x2(2) = jp - b(1,1)/(b(1,2)-b(1,1))
            else  !curvilinear
              step   = 0.0
              x2(1)  = 1.0
              x2(2)  = 1.0
              acc    = 1.e-3
              call cg(x2,err,its,step,acc,10,2,2,
     &                ztecnf,ztecng,ztecnb,ztecnp,w, ldebug)
              if     (its.lt.0) then  !very flat extrema
                x2(1)  = 1.0
                x2(2)  = 1.0
              elseif (min(x2(1),x2(2)).lt.-1.0 .or.
     &                max(x2(1),x2(2)).gt. 3.0     ) then  !very bad cg result
                x2(1)  = 1.0
                x2(2)  = 1.0
              endif
              x2(1) = ip + x2(1)-1.0
              x2(2) = jp + x2(2)-1.0
            endif !rectilinear:curvilinear
            if     (lperiod(1) .and. x2(1).gt.idm) then
              i_out(ii,jj) = idm
              x_out(ii,jj) = max( 0.d0, min( 1.d0, x2(1)-idm ))
            else
              i_out(ii,jj) = max( 1,    min( idm-1, int(x2(1)) ))
              x_out(ii,jj) = max( 0.d0, min( 1.d0, x2(1)-i_out(ii,jj) ))
            endif
            j_out(ii,jj) = max( 1,    min( jdm-1, int(x2(2)) ))
            y_out(ii,jj) = max( 0.d0, min( 1.d0, x2(2)-j_out(ii,jj) ))
          endif !exact point:else
          if     (ldebug) then
            write(6,'(a,4i5,3f9.2)')
     &        'ii,jj,i_,j_,x_,y_,dist = ',
     &         ii,jj,i_out(ii,jj),j_out(ii,jj),
     &               x_out(ii,jj),y_out(ii,jj),dist
            call flush(6)
          endif !ldebug
          endif !i=1,j=1,sea-point
        enddo !ii
      enddo !jj
