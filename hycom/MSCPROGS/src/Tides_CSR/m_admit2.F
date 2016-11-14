      module m_admit2
      contains
      subroutine admit2 (pseudo,radiation,u,v,x,y,h1,h2,amp,pha)
c
c...to compute admittance and in-phase and quadrature components
c   from orthoweights
c
c...coded by Richard Eanes, UT/CSR, 1993 and 1994
c
c Inputs:
c  pseudo   -- Logical variable, true for CR91 style pseudo-orthoweights
c              Use output from tptide for this
c  radiation - Logical variable, true if the orthoweights are corrected
c              for the radiation potential at T2,S2,R2  (R2 not in this list)
c  u,v       - The ortwhoweights for the diurnal and semidiurnal bands
c
c Outputs:
c  x,y       - The tidal admittances at a list of frequencies spanning
c              the diurnal and semidiurnal bands (see data below)
c  h1,h2     - In-phase and quadrature amplitudes at the same frequencies (cm)
c  amp,pha   - Amplitude and phase at the same list of frequencies. (cm,deg)
c
cGE      implicit double precision (a-h,o-z)
      parameter (nf=30)
      logical first,pseudo,radiation
      character*11 tname
      dimension u(3,2), v(3,2)
      dimension amp(100),pha(100),x(100),y(100),h1(100),h2(100)
      dimension tname(nf),freq(nf),per(nf),hsm(nf),chi(nf)
      dimension i1(0:2),i2(0:2),xp1(nf),xp2(nf),fxy(nf),sp(6,2)
      save

c...the orthotide weight factors (signs with + were wrong in C&R 90 Table A1)
      data ((sp(i,m),i=1,6),m=1,2) /
     . 0.0298,  0.1408, +0.0805,  0.6002, +0.3025,  0.1517,
     . 0.0200,  0.0905, +0.0638,  0.3476, +0.1645,  0.0923/
c
      data (tname(i),freq(i),per(i),hsm(i),chi(i),i=1,nf) /
c              name        f(cycle/d)       period      hs (m)    phase
c long period
     .     'LP   55.565',  0.000147094,    6798.382,    0.0279,   180.0, ! 1
     .     'Sa   56.554',  0.002737779,     365.260,   -0.0049,     0.0, ! 2
     .     'Ssa  57.555',  0.005475819,     182.621,   -0.0310,     0.0, ! 3
     .     'TERa 58.554',  0.008213597,     121.749,   -0.0018,     0.0, ! 4
     .     'Mm   65.455',  0.036291647,      27.555,   -0.0352,     0.0, ! 5
     .     'Mf   75.555',  0.073202203,      13.661,   -0.0666,     0.0, ! 6
     .     'TERm 85.455',  0.109493850,       9.133,   -0.0128,     0.0, ! 7
     .     '93a  93.555',  0.140928587,       7.096,   -0.0020,     0.0, ! 8
c diurnal
     .     '2Q1 125.755',  0.856952384,      -6.904,   -0.0066,   270.0, ! 9
     .     'Q1  135.655',  0.893244048,      -9.213,   -0.0502,   270.0, !10
     .     'O1  145.555',  0.929535695,     -13.841,   -0.2622,   270.0, !11
     .     'M1  155.655',  0.966446233,     -28.297,    0.0206,    90.0, !12
     .     'P1  163.555',  0.997262079,    -221.060,   -0.1220,   270.0, !13
     .     'S1  164.556',  1.000000000,    -559.897,    0.0029,    90.0, !14
     .     'K1  165.555',  1.002737898,    1050.246,    0.3688,    90.0, !15
     .     'PHI1167.555',  1.008213699,     155.580,    0.0053,    90.0, !16
     .     'J1  175.455',  1.039029527,      26.850,    0.0206,    90.0, !17
     .     'OO1 185.555',  1.075940083,      13.485,    0.0113,    90.0, !18
     .     'NU1 195.455',  1.112231730,       9.054,    0.0022,    90.0, !19
c semi-diurnal
     .     '227 227.655',  1.828255527,      -5.704,    0.0047,     0.0, !20
     .     '2N2 235.755',  1.859690264,      -6.950,    0.0160,     0.0, !21
     .     'MU2 237.555',  1.864547173,      -7.193,    0.0193,     0.0, !22
     .     'N2  245.655',  1.895981946,      -9.295,    0.1210,     0.0, !23
     .     'NU2 247.455',  1.900838820,      -9.734,    0.0230,     0.0, !24
     .     'M2  255.555',  1.932273593,     -14.026,    0.6319,     0.0, !25
     .     'L2  265.455',  1.968565204,     -28.566,   -0.0179,   180.0, !26
     .     'T2  272.556',  1.997262163,    -158.475,    0.0172,     0.0, !27
     .     'S2  273.555',  2.000000000,    -279.994,    0.2940,     0.0, !28
     .     'K2  275.555',  2.005475795,     525.123,    0.0800,     0.0, !29
     .     '285 285.455',  2.041767407,      26.181,    0.0045,     0.0/ !30
      data i1 / 1, 9, 20/,  i2 / 8, 19, 30/
      data dt /2.d0/
      data rad_fac /0.97d0/, rad_phase /5.9d0/   ! Ray-Sanchez-Cartrwight radiation factors
      data first /.true./

c...compute constants depending only on frequency
      first = .true.
      if (first) then
         first = .false.
         twopi = 8.d0*atan(1.d0)
         deg = 360.d0/twopi
         do m = 1,2
           mf = (-1)**m
           do i = i1(m),i2(m)
              theta = twopi*freq(i)*dt
              cw = 2.*cos(theta)
              sw = 2.*sin(theta)
              xp1(i) = sp(2,m) - sp(3,m)*cw
              if (pseudo) then
* Cartwright&Ray pseudo-orthoweights
                 xp2(i) = sp(4,m) - sp(5,m)*cw - sp(6,m)*sw
              else
* True orthoweights
                 xp2(i) = sp(4,m) - sp(5,m)*cw + sp(6,m)*sw
              endif
              fxy(i) = mf*abs(hsm(i))*1000.d0   ! scale to return values in mm
c             write (12,200) tname(i),fxy(i)*sp(1,m),fxy(i)*xp1(i),
c    .                                               fxy(i)*xp2(i)
           enddo
         enddo
      endif

      do m = 1,2
         do i = i1(m),i2(m)
            x(i) = sp(1,m)*u(1,m) + xp1(i)*u(2,m) + xp2(i)*u(3,m)
            y(i) = sp(1,m)*v(1,m) + xp1(i)*v(2,m) + xp2(i)*v(3,m)
            h1(i) =  fxy(i)*x(i)
            h2(i) = -fxy(i)*y(i)
            amp(i) = sqrt(h1(i)**2 + h2(i)**2)
            if (amp(i).ne.0.) then
               pha(i) = atan2(h2(i),h1(i))*deg
               if (pha(i).lt.0.) pha(i) = pha(i) + 360.
            else
               pha(i) = 0.
            endif
c...correct T2 and S2 for radiation potential
            if (radiation .and. (i.eq.27 .or. i.eq.28) ) then
               amp(i) = amp(i)*rad_fac
               pha(i) = pha(i) + rad_phase
               h1(i) = amp(i)*cos(pha(i)/deg)
               h2(i) = amp(i)*sin(pha(i)/deg)
            endif
c           write (6,100) i,m,tname(i),hsm(i),freq(i),x(i),y(i),
c    .                    h1(i),h2(i),amp(i),pha(i)
        enddo
      enddo

  100 format (1x,2i2,1x,a11,1x,f8.4,f8.5,1x,4f9.4,1x,2f8.2)
  200 format (a11,3x,3f15.6)

      end subroutine
      end module

