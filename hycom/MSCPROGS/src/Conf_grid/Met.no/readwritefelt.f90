PROGRAM readwritefelt

!// Example where the first field in a met.no sequential model file is read
!// and written to a new file.
!// Author: Jon Albretsen, met.no

IMPLICIT NONE

CHARACTER(LEN=99)  :: ifile        ! Input file with model data
CHARACTER(LEN=99)  :: ofile        ! Output file

INTEGER, PARAMETER :: im=720, jm=578 ! Grid dimensions (alternativ, read std. inp. and allocate arrays)
REAL               :: field(im,jm) ! Field
INTEGER            :: mode         ! Choose how to read felt-file
REAL               :: gparam(6)    ! Grid parameters
INTEGER            :: ident(20)    ! Identity string for the sequential data
INTEGER            :: igtype       ! Grid type
INTEGER            :: ierr         ! Error indicator

!ifile = "ifile.seq"
ifile = "CONMAN4km-3.seq_sgi"
ofile = "ofile.seq"

OPEN(UNIT=11,FILE=ifile,ACCESS='SEQUENTIAL',STATUS='UNKNOWN',FORM='UNFORMATTED')
OPEN(UNIT=21,FILE=ofile,ACCESS='SEQUENTIAL',STATUS='UNKNOWN',FORM='UNFORMATTED')

mode = 0   ! Read fields sequentially
! mode = 201 ! Find specified field (see getflt in inout.f90)

CALL getflt(mode,11,0,ident,im*jm,field,igtype,gparam,0,ierr)
IF (ierr /= 0) THEN
  PRINT *,"Error with reading file ",ifile
  STOP
ENDIF

print*, 'Ident is ', ident
print*, ' Topo ', field(1:4,1:4)

CALL putflt(21,ident,igtype,gparam,im*jm,field,1.,ierr)
IF (ierr /= 0) THEN
  PRINT *,"Error with reading file ",ofile
  STOP
ENDIF

CLOSE(11)
CLOSE(21)

END PROGRAM readwritefelt

