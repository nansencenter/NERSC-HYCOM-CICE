      SUBROUTINE dmpwin(nx,ny,nz,u,v,filen)

      CHARACTER*10 filen
      REAL   u(nx,ny,nz),v(nx,ny,nz)
      INTEGER nx,ny,nz,i,j,k
      OPEN(50,FILE=filen,FORM='UNFORMATTED',STATUS='UNKNOWN',ERR=100)
      WRITE(50)nx,ny,nz
      j=1
      k=1
      WRITE(50)(u(i,j,k),i=1,nx*ny*nz)
      WRITE(50)(v(i,j,k),i=1,nx*ny*nz)
      CLOSE(50)
      GOTO 101
 100  WRITE(*,*)'Failed to open: ',filen 
 101  RETURN 
      END

       
      



