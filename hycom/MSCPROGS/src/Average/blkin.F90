module mod_blkin
  use mod_za    ! HYCOM array I/O interface
  implicit none

  contains

  subroutine blkind(dvar,cvar,cfmt)
  implicit none

  real(kind=8) :: dvar
  character(len=*) :: cfmt
  character(len=6) :: cvar

  ! integer ::lp
  ! common/linepr/lp

  ! read in one d.p. value from stdin

  character(len=6) :: cvarin

  read(*,*) dvar,cvarin
  write(lp,cfmt) cvarin,dvar
  call flush(lp)

  if (cvar.ne.cvarin) then
    write(lp,*) 
    write(lp,*) 'error in blkind - input ',cvarin,' but should be ',cvar
    write(lp,*) 
    call flush(lp)
    stop
  endif
  return
  end subroutine blkind

  subroutine blkinr(rvar,cvar,cfmt)
  implicit none

  real(kind=4) :: rvar
  character(len=*) :: cfmt
  character(len=6) :: cvar

  !integer :: lp
  !common/linepr/lp

  !read in one real value from stdin

  character(len=6) :: cvarin

  read(*,*) rvar,cvarin
  write(lp,cfmt) cvarin,rvar
  call flush(lp)

  if (cvar.ne.cvarin) then
     write(lp,*) 
     write(lp,*) 'error in blkinr - input ',cvarin,' but should be ',cvar
     write(lp,*) 
     call flush(lp)
     stop
  endif
  return
  end subroutine blkinr

  subroutine blkinr2(rvar,nvar,cvar1,cfmt1,cvar2,cfmt2)
  implicit none

  real(kind=4) :: rvar
  integer :: nvar
  character(len=*) :: cfmt1, cfmt2
  character(len=6) :: cvar1,cvar2

  ! integer :: lp
  ! common/linepr/lp

  !read in one real value from stdin,
  !identified as either cvar1 (return nvar=1) or cvar2 (return nvar=2)

  character*6 cvarin

  read(*,*) rvar,cvarin

  if (cvar1.eq.cvarin) then
    nvar = 1
    write(lp,cfmt1) cvarin,rvar
    call flush(lp)
  elseif (cvar2.eq.cvarin) then
    nvar = 2
    write(lp,cfmt2) cvarin,rvar
    call flush(lp)
  else
    write(lp,*) 
    write(lp,*) 'error in blkinr2 - input ',cvarin,' but should be ',cvar1,' or ',cvar2
    write(lp,*) 
    call flush(lp)
    stop
  endif
  return
  end subroutine blkinr2
    
  subroutine blkinr3(rvar,nvar,cvar1,cfmt1,cvar2,cfmt2,cvar3,cfmt3)
  implicit none

  real(kind=4) :: rvar
  integer :: nvar
  character(len=*) :: cfmt1,cfmt2,cfmt3
  character(len=6) :: cvar1,cvar2,cvar3

  ! integer :: lp
  ! common/linepr/lp

  !     read in one of three possible real values from stdin,
  !     identified as either cvar1 (return nvar=1) or
  !                          cvar2 (return nvar=2) or
  !                          cvar3 (return nvar=3)

  call blkinr9(rvar,nvar,cvar1,cfmt1,cvar2,cfmt2,cvar3,cfmt3, &
                  'XXXXXX','("blkinr: ",a6," =",f11.4," ")', &
                  'XXXXXX','("blkinr: ",a6," =",f11.4," ")', &
                  'XXXXXX','("blkinr: ",a6," =",f11.4," ")', &
                  'XXXXXX','("blkinr: ",a6," =",f11.4," ")', &
                  'XXXXXX','("blkinr: ",a6," =",f11.4," ")', &
                  'XXXXXX','("blkinr: ",a6," =",f11.4," ")')
  return
  end subroutine blkinr3

  subroutine blkinr9(rvar,nvar, &
                        cvar1,cfmt1,cvar2,cfmt2,cvar3,cfmt3, &
                        cvar4,cfmt4,cvar5,cfmt5,cvar6,cfmt6, &
                        cvar7,cfmt7,cvar8,cfmt8,cvar9,cfmt9)
  implicit none

  real(kind=4) :: rvar
  integer :: nvar
  character(len=6) :: cvar1,cvar2,cvar3,cvar4,cvar5,cvar6,cvar7,cvar8,cvar9
  character(len=*) :: cfmt1,cfmt2,cfmt3,cfmt4,cfmt5,cfmt6,cfmt7,cfmt8,cfmt9

  ! integer :: lp
  ! common/linepr/lp

  !     read in one of nine possible real values from stdin,
  !     identified as either cvar1 (return nvar=1) or
  !                          cvar2 (return nvar=2) or
  !                          ...
  !                          cvar9 (return nvar=9)

  character*6 cvarin

  read(*,*) rvar,cvarin

  if (cvar1.eq.cvarin) then
     nvar = 1
     write(lp,cfmt1) cvarin,rvar
     call flush(lp)
  elseif (cvar2.eq.cvarin) then
     nvar = 2
     write(lp,cfmt2) cvarin,rvar
     call flush(lp)
  elseif (cvar3.eq.cvarin) then
     nvar = 3
     write(lp,cfmt3) cvarin,rvar
     call flush(lp)
  elseif (cvar4.eq.cvarin) then
     nvar = 4
     write(lp,cfmt4) cvarin,rvar
     call flush(lp)
  elseif (cvar5.eq.cvarin) then
     nvar = 5
     write(lp,cfmt5) cvarin,rvar
     call flush(lp)
  elseif (cvar6.eq.cvarin) then
     nvar = 6
     write(lp,cfmt6) cvarin,rvar
     call flush(lp)
  elseif (cvar7.eq.cvarin) then
    nvar = 7
    write(lp,cfmt7) cvarin,rvar
    call flush(lp)
  elseif (cvar8.eq.cvarin) then
     nvar = 8
     write(lp,cfmt8) cvarin,rvar
     call flush(lp)
  elseif (cvar9.eq.cvarin) then
     nvar = 9
     write(lp,cfmt9) cvarin,rvar
     call flush(lp)
  else
     write(lp,*) 
     write(lp,*) 'error in blkinr9 - input ',cvarin,' but should be:'
     write(lp,*) cvar1,' or ',cvar2,' or ',cvar3,' or'
     write(lp,*) cvar4,' or ',cvar5,' or ',cvar6,' or'
     write(lp,*) cvar7,' or ',cvar8,' or ',cvar9
     write(lp,*) 
     call flush(lp)
     stop
  endif
  return
  end subroutine blkinr9
    
  subroutine blkini(ivar,cvar)
  implicit none

  integer :: ivar
  character(len=6) :: cvar

  ! integer :: lp
  ! common/linepr/lp

  ! read in one integer value from stdin

  character(len=6) :: cvarin

   read(*,*) ivar,cvarin
   write(lp,6000) cvarin,ivar
   call flush(lp)

  if (cvar.ne.cvarin) then
     write(lp,*) 
     write(lp,*) 'error in blkini - input ',cvarin,' but should be ',cvar
     write(lp,*) 
     call flush(lp)
     stop
  endif
  return
  6000 format('blkini: ',a6,' =',i6)
  end subroutine blkini
    
  subroutine blkini2(ivar,nvar,cvar1,cvar2)
  implicit none

  integer :: ivar,nvar
  character(len=6) :: cvar1,cvar2

  ! integer :: lp
  ! common/linepr/lp

  ! read in one integer value from stdin
  ! identified as either cvar1 (return nvar=1) or cvar2 (return nvar=2)

  character(len=6) :: cvarin

  read(*,*) ivar,cvarin
  write(lp,6000) cvarin,ivar
  call flush(lp)

  if (cvarin.eq.cvar1) then
     nvar = 1
  elseif (cvarin.eq.cvar2) then
     nvar = 2
  else
     write(lp,*) 
     write(lp,*) 'error in blkini2 - input ',cvarin,' but should be ',cvar1,' or ',cvar2
     write(lp,*) 
     call flush(lp)
     stop
  endif
  return
6000  format('blkini: ',a6,' =',i6)
  end subroutine blkini2
    
  subroutine blkini3(ivar,nvar,cvar1,cvar2,cvar3)
  implicit none

  integer ::ivar,nvar
  character(len=6) :: cvar1,cvar2,cvar3

  ! integer :: lp
  ! common/linepr/lp

  ! read in one of three possible integer values from stdin
  ! identified as either cvar1 (return nvar=1) or
  !                      cvar2 (return nvar=2) or
  !                      cvar3 (return nvar=3)

  call blkini9(ivar,nvar,cvar1,cvar2,cvar3, &
                            'XXXXXX','XXXXXX','XXXXXX', &
                            'XXXXXX','XXXXXX','XXXXXX')
  return
  end subroutine blkini3
    
  subroutine blkini9(ivar,nvar,cvar1,cvar2,cvar3, &
                                   cvar4,cvar5,cvar6, &
                                   cvar7,cvar8,cvar9)
  implicit none

  integer :: ivar,nvar
  character(len=6) :: cvar1,cvar2,cvar3,cvar4,cvar5,cvar6,cvar7,cvar8,cvar9

  ! integer :: lp
  ! common/linepr/lp

  ! read in one of nine possible integer values from stdin
  ! identified as either cvar1 (return nvar=1) or
  !                      cvar2 (return nvar=2) or
  !                      cvar3 (return nvar=3) or
  !                      ...
  !                      cvar9 (return nvar=9)

  character(len=6) :: cvarin

  read(*,*) ivar,cvarin
  write(lp,6000) cvarin,ivar
  call flush(lp)

  if     (cvarin.eq.cvar1) then
     nvar = 1
  elseif (cvarin.eq.cvar2) then
     nvar = 2
  elseif (cvarin.eq.cvar3) then
     nvar = 3
  elseif (cvar4.eq.cvarin) then
     nvar = 4
  elseif (cvar5.eq.cvarin) then
     nvar = 5
  elseif (cvar6.eq.cvarin) then
     nvar = 6
  elseif (cvar7.eq.cvarin) then
     nvar = 7
  elseif (cvar8.eq.cvarin) then
     nvar = 8
  elseif (cvar9.eq.cvarin) then
     nvar = 9
  else
     write(lp,*) 
     write(lp,*) 'error in blkini9 - input ',cvarin,' but should be:'
     write(lp,*) cvar1,' or ',cvar2,' or ',cvar3,' or'
     write(lp,*) cvar4,' or ',cvar5,' or ',cvar6,' or'
     write(lp,*) cvar7,' or ',cvar8,' or ',cvar9
     write(lp,*) 
     call flush(lp)
     stop
  endif
  return
 6000 format('blkini: ',a6,' =',i6)
  end subroutine blkini9
    
  subroutine blkinl(lvar,cvar)
  implicit none

  logical :: lvar
  character(len=6) :: cvar

  ! integer :: lp
  ! common/linepr/lp

  ! read in one logical value from stdin
  ! due to a SGI bug for logical I/O: read in an integer 0=F,1=T

  character*6 cvarin
  integer     ivar

  read(*,*) ivar,cvarin
  lvar = ivar .ne. 0
  write(lp,6000) cvarin,lvar
  call flush(lp)

  if (cvar.ne.cvarin) then
     write(lp,*) 
     write(lp,*) 'error in blkinr - input ',cvarin,' but should be ',cvar
     write(lp,*) 
     call flush(lp)
     stop
  endif
  return
  6000 format('blkinl: ',a6,' =',l6)
  end subroutine blkinl
    
end module mod_blkin
