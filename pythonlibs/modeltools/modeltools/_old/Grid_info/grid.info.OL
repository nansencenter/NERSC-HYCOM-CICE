  78.0    90.0     ! Position of pole N (lattitude, longitude):
  48.0    00.0     ! Position of pole S (lattitude, longitude):
 -27.0  -02.0 250   ! Longitude interval West lon,   East lon, idim
 -50.0   30.0 325   ! Lattitude interval south limit, north limit, jdim
.true.            ! Generate topography
.true.            ! dump teclatlon.dat
.true.            ! dump micom latlon.dat
.true.            ! mercator grid (true, false)
 -50.0 .false.     ! Lat' of southern boundary
.false.            ! Smoothing, Shapiro filter
 8   2            ! Order of Shapiro filter,  number of passes
