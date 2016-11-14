      SUBROUTINE cice(ice,iw,nx,ny,mo,
     +        iflg,mlon,mlat)

      INTEGER iw,nx,ny,mo,iflg(nx,ny)
      REAL ice(nx,ny,mo),mlon(nx,ny),mlat(nx,ny),
     +    a,b,vmin,vmax
      CHARACTER*20 filen

      a=0.     !Const. not used for this dataset
      b=0.
      vmin=0.
      vmax=0.
       
     
      IF(iw.EQ.1) THEN 
       filen= 'Data/ice_5190.dat'    !Read Mean ice conc 1951-90 
       WRITE(*,*)'Mean ice conc. 1951--90 (0-1) (cice.f)'
       CALL lwalsh(ice,iflg,nx,ny,mlon,mlat,filen,a,b,
     +            vmin,vmax)
       write(*,*) 'lwalsh'
                                      !Fillup missing data
                                      !and land mask.
       CALL fillup(ice,nx,ny,iflg,12,0.)

      ELSEIF(iw.EQ.2)THEN
       filen= 'Data/ice_0190.dat'    !Read Mean ice conc 1951-90 
       WRITE(*,*)'Mean ice conc. 1901--90 (0-1) (cice.f)'
       CALL lwalsh(ice,iflg,nx,ny,mlon,mlat,filen,a,b,
     +            vmin,vmax)
                                      !Fillup missing data
                                      !and land mask.
       CALL fillup(ice,nx,ny,iflg,12,0.)
   
      ELSE
       WRITE(*,*)'Ice conc does not exist'
      ENDIF
c
      write(*,*) 'out of cice'
c
      RETURN
      END



