
      SUBROUTINE rmask(iflg,nx,ny)
c -------------------------------------------------------------
c --- Reads the landmask into iflg.
c -------------------------------------------------------------
      DIMENSION iflg(nx,ny)
      INTEGER nx,ny,i,j 
c.ver1      IF(nx.ne.71) STOP 'Test format in rmask' 
      IF(ny.ne.71) STOP 'Test format in rmask' 
      WRITE(*,*)'Loading mask from file flags.dat'
      OPEN(13,FILE='Data/flags.dat',STATUS='unknown')
      READ(13,'(71I1)')((iflg(i,j),j=ny,1,-1),i=nx,1,-1)
c.ver1      READ(13,'(100I1)',END=15)((iflg(i,j),i=1,nx),j=ny,1,-1)
  15  CONTINUE
      CLOSE(13)
      RETURN
      END

