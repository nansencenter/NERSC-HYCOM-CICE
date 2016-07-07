      subroutine extrct_p(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
      real, allocatable :: wrot(:,:)
      integer           :: i,j
c
      allocate( wrot(m,n) )
      do j= 1,m
        do i= 1,n
          wrot(m+1-j,i) = work(i,j)
        enddo
      enddo
c
      call extrct_px(wrot,m,n,io,jo,array,no,mo)
c
      deallocate( wrot )
c
      return
      end
      subroutine extrct_pv(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
      real, allocatable :: wrot(:,:)
      integer           :: i,j
c
      allocate( wrot(m,n) )
      do j= 1,m
        do i= 1,n
          wrot(m+1-j,i) = work(i,j)
        enddo
      enddo
c
      call extrct_pvx(wrot,m,n,io,jo,array,no,mo)
c
      deallocate( wrot )
c
      return
      end
      subroutine extrct_q(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
      real, allocatable :: wrot(:,:)
      integer           :: i,j
c
      allocate( wrot(m,n) )
      do j= 1,m
        do i= 1,n
          wrot(m+1-j,i) = work(i,j)
        enddo
      enddo
c
      call extrct_qx(wrot,m,n,io,jo,array,no,mo)
c
      deallocate( wrot )
c
      return
      end
      subroutine extrct_u(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
      real, allocatable :: wrot(:,:)
      integer           :: i,j
c
      allocate( wrot(m,n) )
      do j= 1,m
        do i= 1,n
          wrot(m+1-j,i) = work(i,j)
        enddo
      enddo
c
      call extrct_ux(wrot,m,n,io,jo,array,no,mo)
c
      deallocate( wrot )
c
      return
      end
      subroutine extrct_us(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
      real, allocatable :: wrot(:,:)
      integer           :: i,j
c
      allocate( wrot(m,n) )
      do j= 1,m
        do i= 1,n
          wrot(m+1-j,i) = work(i,j)
        enddo
      enddo
c
      call extrct_usx(wrot,m,n,io,jo,array,no,mo)
c
      deallocate( wrot )
c
      return
      end
      subroutine extrct_v(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
      real, allocatable :: wrot(:,:)
      integer           :: i,j
c
      allocate( wrot(m,n) )
      do j= 1,m
        do i= 1,n
          wrot(m+1-j,i) = work(i,j)
        enddo
      enddo
c
      call extrct_vx(wrot,m,n,io,jo,array,no,mo)
c
      deallocate( wrot )
c
      return
      end
      subroutine extrct_vs(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
      real, allocatable :: wrot(:,:)
      integer           :: i,j
c
      allocate( wrot(m,n) )
      do j= 1,m
        do i= 1,n
          wrot(m+1-j,i) = work(i,j)
        enddo
      enddo
c
      call extrct_vsx(wrot,m,n,io,jo,array,no,mo)
c
      deallocate( wrot )
c
      return
      end

c--------------------------------------------------------------------------
c     ORIGINAL EXTRCT ROUTINES.
c--------------------------------------------------------------------------
      subroutine extrct_px(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the p-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      do j=1,min(mo,m-jo+1)
        jw = jo-1+j
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
c
      do j=m-jo+2,mo  ! arctic patch
        jw = m-1-(j-(m-jo+1))
        do i=1,no
          iw = n-mod(io+i-2,n)
          array(i,j)=work(iw,jw)
        enddo
      enddo
      return
      end
c
      subroutine extrct_pvx(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the p-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c --- for vector fields, use extrct_p for scalar fields.
c --- note the sign change across the bipolar seam.
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      do j=1,min(mo,m-jo+1)
        jw = jo-1+j
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
c
      do j=m-jo+2,mo  ! arctic patch
        jw = m-1-(j-(m-jo+1))
        do i=1,no
          iw = n-mod(io+i-2,n)
          array(i,j)=-work(iw,jw)
        enddo
      enddo
      return
      end
c
      subroutine extrct_qx(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the q-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      do j=1,min(mo,m-jo+1)
        jw = jo-1+j
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
c
      do j=m-jo+2,mo  ! arctic patch
        jw = m-(j-(m-jo+1))
        do i=1,no
          iw = mod(n-mod(io+i-2,n),n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
      return
      end

      subroutine extrct_ux(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the u-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c --- for vector fields, use extrct_us for scalar fields.
c --- note the sign change across the bipolar seam.
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      do j=1,min(mo,m-jo+1)
        jw = jo-1+j
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
c
      do j=m-jo+2,mo  ! arctic patch
        jw = m-1-(j-(m-jo+1))
        do i=1,no
          iw = mod(n-mod(io+i-2,n),n)+1
          array(i,j)=-work(iw,jw)
        enddo
      enddo
      return
      end
c
      subroutine extrct_usx(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the u-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c --- for scalar fields, use extrct_u for vector fields.
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      do j=1,min(mo,m-jo+1)
        jw = jo-1+j
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
c
      do j=m-jo+2,mo  ! arctic patch
        jw = m-1-(j-(m-jo+1))
        do i=1,no
          iw = mod(n-mod(io+i-2,n),n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
      return
      end
c
      subroutine extrct_vx(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the v-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c --- for vector fields, use extrct_vs for scalar fields.
c --- note the sign change across the bipolar seam.
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      do j=1,min(mo,m-jo+1)
        jw = jo-1+j
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
c
      do j=m-jo+2,mo  ! arctic patch
        jw = m-(j-(m-jo+1))
        do i=1,no
          iw = n-mod(io+i-2,n)
          array(i,j)=-work(iw,jw)
        enddo
      enddo
      return
      end
c
      subroutine extrct_vsx(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for the v-grid of global domains with arctic bi-polar patch.
c --- it will also work for closed and near-global domains.
c --- for scalar fields, use extrct_v for vector fields.
c
      integer i,iw,j,jw
c
      if     (m.le.6) then
        call extrct_2d(work,n,m,io,jo,array,no,mo)  ! 2-d domain
        return
      endif
c
      do j=1,min(mo,m-jo+1)
        jw = jo-1+j
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,jw)
        enddo
      enddo
c
      do j=m-jo+2,mo  ! arctic patch
        jw = m-(j-(m-jo+1))
        do i=1,no
          iw = n-mod(io+i-2,n)
          array(i,j)=work(iw,jw)
        enddo
      enddo
      return
      end

      subroutine extrct_2d(work,n,m,io,jo,array,no,mo)
      implicit none
c
      integer n,m,io,jo,no,mo
      real    work(n,m),array(no,mo)
c
c --- array = work(io:io+no-1,jo:jo+mo-1)
c
c --- for 2-d (m<6, infinate f-plane) domains.
c
      integer i,iw,j,jw
c
      do j=1,mo
        do i=1,no
          iw = mod(io+i-2,n)+1
          array(i,j)=work(iw,1)
        enddo
      enddo
      return
      end
c


