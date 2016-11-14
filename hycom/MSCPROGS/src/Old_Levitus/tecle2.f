      SUBROUTINE tecle2(nx,ny,nz,zl,sl,title,filen,
     +              nx1,nx2,ny1,ny2,nz1,nz2)
c --------------------------------------------------------------
c --- Author: Knud Simonsen, NERSC 14/12 1993
c --------------------------------------------------------------
c --- Dumps the data, or part in 'data' on file in a 
c --- format, which is accepted by 'tecplot'.  
c --- ----------------------------------------------------------
c --- Input:    zl: depths  (1-dim)
c ---           sl: salinity 
c ---     nx,ny,nz: Size of array.
c ---      nx1,nx2: x-indexes which are filed
c ---      ny1,ny2: y-indexes which are filed
c ---      nz1,nz2: z-indexes which are filed.
c ---      title:   Title of the plot
c ---      filen:   Filename
c --- Output: None.
c --------------------------------------------------------------- 
      INTEGER nx,ny,nz,nx1,nx2,ny1,ny2,nz1,nz2,
     +        dx,dy,dz
      REAL zl(nz),sl(nx,ny,nz),dzz,to
      CHARACTER title*10
      character(len=*) filen

      dx=nx2-nx1+1
      dy=ny2-ny1+1
      dz=nz2-nz1+1

      OPEN(10,FILE=filen,STATUS='UNKNOWN')
      WRITE(10,105)                    
      !WRITE(10,101)
      !  WRITE(10,104)dx,dy,dz
      !  WRITE(10,99)(((i,i=nx1,nx2),j=ny1,ny2),k=nz1,nz2)
      !  WRITE(10,99)(((j,i=nx1,nx2),j=ny1,ny2),k=nz1,nz2)
      !  WRITE(10,100)(((zl(k),i=nx1,nx2),j=ny1,ny2),k=nz1,nz2)
      !  WRITE(10,100)(((sl(i,j,k),i=nx1,nx2),j=ny1,ny2),k=nz1,nz2)
      !CLOSE(10)
      WRITE(10,201)
      do k=nz1,nz2
        WRITE(10,204)k,dx,dy
        if (k==nz1) then
           WRITE(10,99)((i,i=nx1,nx2),j=ny1,ny2)
           WRITE(10,99)((j,i=nx1,nx2),j=ny1,ny2)
        else
           WRITE(10,205)
        end if
        WRITE(10,100)((sl(i,j,k),i=nx1,nx2),j=ny1,ny2)
      end do
      CLOSE(10)


  105 FORMAT('TITLE = "TITI"')
  101 FORMAT('VARIABLES = "X","Y","Z","S"')
  102 FORMAT('ZONE T= "Layer ',I2,'", I=',I3,', J=',I3,', F=BLOCK')
  103 FORMAT(2I4,F8.1,2F9.3)
  104 FORMAT('ZONE I=',I3,', J=',I3,
     +                  ', K=',I3,', F=BLOCK')
  201 FORMAT('VARIABLES = "X","Y","S"')
  204 FORMAT('ZONE T="Month ',i2, '", I=',I3,', J=',I3, ', F=BLOCK')
  205 FORMAT('D=(1,2)')
   99 FORMAT(20I4)
  100 FORMAT(10(1x,e10.3))

      RETURN
      END 

