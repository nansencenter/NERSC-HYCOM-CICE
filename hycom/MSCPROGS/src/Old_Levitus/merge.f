      PROGRAM merge
c -------------------------------------------------------------
c -------------------------------------------------------------
      parameter(nzd=33,nxd=181,nyd=240)

      INTEGER iflg(nxd,nyd)
      REAL rno3_reid(nxd,nyd,nzd),
     .     rno3_nodc(nxd,nyd,nzd),
     .     rno3     (nxd,nyd),
     .     zl(nzd)

      OPEN(10,FILE='esop_no3_nodc.dat',form='formatted')
      OPEN(11,FILE='esop_no3_reid.dat',form='formatted')
c
      DO k=1,nzd                      !Fillup remaining areas
      write(*,*) 'reading layer ',k
       read(10,'(10f9.4)')rno3
      do i=1,nxd
      do j=1,nyd
        rno3_nodc(i,j,k)=rno3(i,j)
      enddo
      enddo
       read(11,'(10f9.4)')rno3
      do i=1,nxd
      do j=1,nyd
        rno3_reid(i,j,k)=rno3(i,j)
      enddo
      enddo
      enddo
      CLOSE(10)
      CLOSE(11)
c
      open(50,file='mask.check')
      read(50,'(2i6)') ii,jj
      if ((ii.ne.nxd).or.(jj.ne.nyd)) stop
      read(50,*)
      read(50,*)
      read(50,*)
      do 76 i=nxd,1,-1
  76  read(50,'(i3,1x,150i1)') i1,(iflg(i,j),j=nyd,1,-1)
      close(50)
c
      do i=1,nxd
      do j=1,nyd
      if (iflg(i,j).eq.2) then
      do k=1,nzd
      rno3_nodc(i,j,k)=rno3_reid(i,j,k)
      enddo
      endif
      enddo
      enddo
c
      close(11)
      OPEN(11,FILE='esop_no3_merge.dat',form='formatted')
c
      DO k=1,nzd                      !Fillup remaining areas
      write(*,*) 'write layer ',k
       DO i=1,nxd
        DO j=1,nyd
         rno3(i,j)=rno3_nodc(i,j,k)
        ENDDO
       ENDDO
       WRITE(11,'(10f9.4)')rno3
      ENDDO
      CLOSE(11)
c
      dums= 0.
      dumt= 10.
      DO i=1,nzd
       IF(i.GE.4)  dumt= 20.
       IF(i.GE.5)  dumt= 25.
       IF(i.GE.9)  dumt= 50.
       IF(i.GE.12) dumt= 100.
       IF(i.GE.24) dumt = 250.
       IF(i.GE.26) dumt=500.
       zl(i) = dums
       dums = dums + dumt
      ENDDO
c
      CALL tecle2(nxd,nyd,nzd,zl,rno3_nodc,
     +              'Initial','merge_tecplot.dat',
     +              1,nxd,1,nyd,1,nzd)
c
      STOP 'Finito'
      END  

