      subroutine geopar
      use mod_xc  ! HYCOM communication interface
      use mod_rest ! HYCOM restart arrays
      use mod_za  ! HYCOM I/O interface
c
c --- set up model parameters related to geography
c
      implicit none
c
      real      hmina,hminb,hmaxa,hmaxb
      integer   i,j
      character preambl(5)*79,cline*80
c
c --- read basin depth array
c
      if     (mnproc.eq.1) then
      write (lp,'(3a)') ' reading bathymetry file from ',
     &                  flnmdep(1:len_trim(flnmdep)),'.[ab]'
      endif
      call xcsync(flush_lp)
      open (unit=9,file=flnmdep(1:len_trim(flnmdep))//'.b',status='old')
      read (     9,'(a79)')  preambl
      read (     9,'(a)')    cline
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      close(unit=9)
      if     (mnproc.eq.1) then
      write (lp,'(/(1x,a))') preambl,cline
      endif
c
      call zaiopf(flnmdep(1:len_trim(flnmdep))//'.a','old', 9)
      call zaiord(depths,ip,.false., hmina,hmaxa, 9)
      call zaiocl(9)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        endif
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j= 1,jdm
        do i= 1,idm
          if     (depths(i,j).gt.0.5*huge) then
            depths(i,j) = 0.0
          endif
        enddo
      enddo
c
c --- determine masks for u,v,p,q points
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,idm
        do i=1,jdm
          ip(i,j)=0
          iu(i,j)=0
          iv(i,j)=0
        enddo
      enddo
c
c --- mass points are defined where water depth is greater than zero
c --- u,v points are located halfway between any 2 adjoining mass points
c
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jdm-1
        do i=1,idm-1
          if (depths(i,j).gt.0.) then
            ip(i,j)=1
          endif
        enddo
        do i=2,idm-1
          if (ip(i-1,j).gt.0.and.ip(i,j).gt.0) then
            iu(i,j)=1
          endif
        enddo
      enddo
!$OMP PARALLEL DO PRIVATE(j,i)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=2,jdm-1
        do i=1,idm-1
          if (ip(i,j-1).gt.0.and.ip(i,j).gt.0) then
            iv(i,j)=1
          endif
        enddo
      enddo
c
      return
      end
