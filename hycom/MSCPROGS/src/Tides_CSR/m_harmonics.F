      module m_harmonics
      contains
      SUBROUTINE harmonics(jj,kk,nconst,dlon,dlat
     .                           ,xmjd,m,k,ifreq,amppha,deep)
      use m_csrtide
      use m_admit2
! --- ----------------------------------------------------------
! --- Utilizing files from the data originator to read data
! --- file and estimate amplitudes and phases. The result
! --- is put into array amppha
! --- ------------------------------------------------------- 
cGE       implicit double precision (a-h,o-z)
       dimension amppha(2*(jj+kk),2*nconst)
       INTEGER ifreq(nconst)
       real deep
       logical pseudo,radiation,isdata
       dimension mdyhms(6), u(3,2), v(3,2)
       dimension x(100),y(100),h1(100),h2(100),amp(100),pha(100)

       call csrtide (xmjd,dlat,dlon,tide,tlp,tds,isdata,
     .              u, v, pseudo, radiation)
       call admit2 (pseudo,radiation,u,v,x,y,h1,h2,amp,pha)


       DO i=1,nconst
         iamp=2*i-1
         ipha=2*i
         if (deep > 0.0) then
            amppha(m,iamp)=amp(ifreq(i))   !Units mm,dgr 
            amppha(m,ipha)=pha(ifreq(i))
         else
            amppha(m,iamp)=-999.9
            amppha(m,ipha)=-999.9
         endif

       ENDDO
      end subroutine
      end module
