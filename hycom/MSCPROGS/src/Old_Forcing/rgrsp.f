      SUBROUTINE rgrsp(nxe,nye,nze,ypivo,xpivn,ypivn,gridn)
c ----------------------------------------------------------------
c --- Author Knud Simonsen, 10/12 1993, NERSC
c ----------------------------------------------------------------
c --- Reads the file 'gridspec.inp', 
c --- Check if the dimensions in the 
c --- file corresponds to dimensions in the model.
c ----------------------------------------------------------------
c --- Input: nxe,nye,nze: Dim. in the program
c --- Output: ypivo The meridiane, which is considered 
c --- as model 'equator' (ME)
c ---         xpivn x-index for the ME meridian
c ---         ypivn y-index for the true equator
c ---         gridn Gridresolution.
c ----------------------------------------------------------------
      INTEGER nxe,nye,nze,nx,ny,nz,i,j
      REAL ypivo,xpivn,ypivn,gridn
      CHARACTER*70 tull



      OPEN(10,FILE='gridspec.inp')
      DO i=1,4
       READ(10,'(A70)')tull
      ENDDO
      READ(10,100)nx,tull
      READ(10,100)ny,tull
      READ(10,100)nz,tull
      READ(10,200)ypivo,tull
      READ(10,200)xpivn,tull
      READ(10,200)ypivn,tull
      READ(10,200)gridn,tull
      CLOSE(10)
      write(*,*)
      write(*,*)'================================================'
      write(*,*) 'Specifications from "gridspec.inp":'
      write(*,201) nx,ny,nz
      write(*,202)ypivo,xpivn,ypivn,gridn
      write(*,*)'------------------------------------------------'  
      IF(nx.NE.nxe.OR.ny.NE.nye.OR.nz.NE.nze) THEN
       write(*,*) 'Dim. in model'
       write(*,201) nxe,nye,nze   
      write(*,*)'------------------------------------------------'
       STOP 'The grid dim. in program and file are not the same'
      ENDIF
      write(*,*)'================ Lets Go ======================='
 100  FORMAT(I3,A70)
 200  FORMAT(f6.1,A70)    
 201  FORMAT('nx=',I3,' ny=',I3,' nz=',I3)
 202  FORMAT('ypivo=',f6.1,' xpivn=',
     +             f6.1,' ypivn=',f6.1,' gridn=',f4.1)
      RETURN
      END 
    
