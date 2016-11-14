      SUBROUTINE dmpdat(nx,ny,nz,field,fact,filen)

      !KAL - HMMM CHARACTER*18 filen
      character(len=*) filen
      REAL field(nx,ny,nz),rfact
      INTEGER nx,ny,nz,i,j,k,fact

      rfact=FLOAT(fact)
 
      
      OPEN(10,FILE=filen,FORM='FORMATTED'
     +        ,STATUS='UNKNOWN',ERR=100)
      WRITE(10,'(4I8)')nx,ny,nz,fact
      DO k=1,nz
       WRITE(10,99)((INT(field(i,j,k)*rfact),i=1,nx),j=1,ny)
      ENDDO
      CLOSE(10)
      WRITE(*,*)'Data written to file: ',filen
  99  FORMAT(15I8)
      GOTO 101
 100  WRITE(*,*)'Failed to open ',filen
 101  RETURN 
      END

       
      



