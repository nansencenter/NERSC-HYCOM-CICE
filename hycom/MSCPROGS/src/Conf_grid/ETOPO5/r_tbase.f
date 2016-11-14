      program r_tbase

c     read the TerrainBase data
c
      dimension itopo(4320)

      open(10,file='ntbase')

c     specific a region to be printed
      alonl = 270.0
      alonr = 300.0
      alatt = 60.0
      alatb = 30.0


      write(*,49)alatb,alonl,alatt,alonr
49    format('Sample Data             :'/,
     *       'Column and Row spacing  : 2.0 degrees'/,
     *       'SW corner (lat/lon)     :',2f6.0/,
     *       'NE corner (lat/lon)     :',2f6.0)


      ilonl = nint(alonl)*12 +1
      ilonr = nint(alonr)*12 +1
      ilatb = nint(90.- alatb)*12 + 1
      ilatt = nint(90.- alatt)*12 + 1
      if(ilonr.eq.4321)ilonr=4320
      if(ilatb.eq.2161)ilatb=2160


      do j=1,2160
         read(10,'(20i6)')(itopo(i),i=1,4320)
         if(j.gt.ilatb)go to 60
         if(j.ge.ilatt.and.j.le.ilatb.and.mod(j,24).eq.0)then
          write(*,'(24i6)')(itopo(i),i=ilonl,ilonr,24)
         endif
      enddo

60    print*,'Processing completed: '

      stop 
      
      end


