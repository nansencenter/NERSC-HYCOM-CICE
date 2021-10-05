










      module m_strmf_eval
      contains

      subroutine strmf_eval(idm,jdm,strmf,ub,vb)
      use mod_grid
      implicit none
c
      integer idm,jdm
      real strmf(idm,jdm)
     .    ,ub(idm,jdm),vb(idm,jdm)
c
      integer i,j,im1,jm1,ip1,jp1
      real, external :: spherdist
c
c --- ------------------------------------------------------------------
c --- integrate the stream function with boundary condition strmf(1,1)=0
c --- ------------------------------------------------------------------
c
      strmf(1,1)=0.
c
      do j=2,jdm
        jm1=j-1
        strmf(1,j)=strmf(1,jm1)
     .            -ub(1,jm1)*depthu(1,jm1)*scuy(1,jm1)
      enddo
c
      do i=2,idm
        im1=i-1
        do j=1,jdm
          strmf(i,j)=strmf(im1,j)
     .              +vb(im1,j)*depthv(im1,j)*scvx(im1,j)
        enddo
      enddo
c
c --- ------------------------------------------------------------------
c --- interpolate the streamfunction to the p-point (also smooths)
c --- ------------------------------------------------------------------
c
      do j=1,jdm
        jp1=mod(j,jdm)+1
        do i=1,idm-1
          ip1=i+1
          strmf(i,j)=.25*(strmf(i,j  )+strmf(ip1,j  )
     .                   +strmf(i,jp1)+strmf(ip1,jp1))
        enddo
      enddo
c
      return
c
      end subroutine

      end module m_strmf_eval
