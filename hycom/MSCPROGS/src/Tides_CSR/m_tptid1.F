      module m_tptid1
      contains
      subroutine tptid1( u,v,ts,pseudo,radiation,tide )

*  computes the tide at time ts from the orthoweights u,v.
*  this is a kernel routine, meant to be called only by tptide.

cGE      implicit double precision (a-h,o-z)
      dimension   indx(30,2,4),amp(30,2),freq(30,2),
     *            pha(30,2),ct(30,2),st(30,2),at(3),bt(3),
     *            u00(2),u20(2),u21(2),u40(2),u41(2),v41(2),
     *            p(3,2),q(3,2),u(3,2),v(3,2),
     *            shpn(4),phc(4),dpd(4)
      logical     init,radiation,pseudo
      parameter  (rad_factor=0.97, rad_phase=5.9)
      save
*
      data phc/290.21, 280.12, 274.35, 343.51/,
     *     dpd/13.1763965,0.9856473,0.1114041,0.0529539/,
     *     u00/0.0298,0.0200/, u20/0.1408,0.0905/,
     *     u21/0.0805,0.0638/, u40/0.6002,0.3476/,
     *     u41/0.3025,0.1645/, v41/0.1517,0.0923/
      data tc/40431744.d2/, tslast/-9.99e10/
      data init/.true./
      data indx                                             !(  l  ,m,n)
     + /-3,-3,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0,       !( 1:15,1,1)
     +   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3,       !(16:30,1,1)
     +  -3,-3,-2,-2,-2,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0,       !( 1:15,2,1)
     +   0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3,       !(16:30,2,1)
     +   0, 2, 0, 0, 2, 2, 0, 0, 0, 0, 2, 0, 0, 0, 2,       !( 1:15,1,2)
     +  -3,-2,-2, 0, 0, 0, 1, 2,-2, 0, 0,-2, 0, 0, 0,       !(16:30,1,2)
     +   0, 2, 0, 2, 3,-1, 0, 0, 1, 2, 3,-2,-1, 0, 0,       !( 1:15,2,2)
     +   1,-2, 0, 0, 0,-3,-2,-1, 0, 0, 0, 0,-2, 0, 0,       !(16:30,2,2)
     +   2, 0, 1, 1,-1,-1, 0, 0, 0, 2, 0,-1, 1, 1,-1,       !( 1:15,1,3)
     +   0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0,-2, 0, 0,       !(16:30,1,3)
     +   3, 1, 2, 0, 0, 1, 1, 1, 1,-1,-1, 2, 0, 0, 0,       !( 1:15,2,3)
     +   0, 1,-1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1,       !(16:30,2,3)
     +   0, 0,-1, 0,-1, 0,-2,-1, 0, 0, 0, 0, 0, 1, 0,       !( 1:15,1,4)
     +   0,-1, 0,-1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1,       !(16:30,1,4)
     +   0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0,-1, 0,       !( 1:15,2,4)
     +   0, 0, 0, 0, 1, 0, 0, 0,-1, 0, 1, 2, 0, 0, 1/       !(16:30,2,4)
      data amp                                                !(  l  ,m)
     + /0.00663,0.00802,0.00947,0.05019,0.0018 ,0.00954,      !( 1: 6,1)
     +  0.00152,0.04946,0.26218,0.00171,0.00343,0.00741,      !( 7:12,1)
     +  0.02062,0.00414,0.00394,0.00713,0.00137,0.122  ,      !(13:18,1)
     +  0.0073 ,0.36874,0.05002,0.00293,0.00524,0.00395,      !(19:24,1)
     +  0.02061,0.00409,0.00342,0.00169,0.01128,0.00723,      !(25:30,1)
     +  0.0018 ,0.00467,0.01601,0.01932,0.0013 ,0.00102,      !( 1: 6,2)
     +  0.00451,0.121  ,0.00113,0.02298,0.00106,0.0019 ,      !( 7:12,2)
     +  0.00218,0.02358,0.63194,0.00193,0.00466,0.01786,      !(13:18,2)
     +  0.00447,0.00197,0.01718,0.29402,0.00305,0.00102,      !(19:24,2)
     +  0.07994,0.02382,0.00259,0.00086,0.00447,0.00195/      !(25:30,2)
      data pha                                                !(  l  ,m)
     + / 90., 90., 90., 90., 90., 90.,270., 90., 90.,270.,    !( 1:10,1)
     +  270.,270.,270.,270.,270., 13.,270., 90., 90.,270.,    !(11:20,1)
     +  270.,-13.,270.,270.,270.,270.,270.,270.,270.,270.,    !(21:30,1)
     +    0.,  0.,  0.,  0., 77.,103.,180.,  0., 77.,  0.,    !( 1:10,2)
     +   77.,180.,103.,180.,  0., 77.,180.,180.,  0.,  0.,    !(11:20,2)
     +  283.,  0.,262.,180.,  0.,  0.,  0.,  0.,  0.,  0./    !(21:30,2)
*
      if (init) then
      init = .false.
      rad = atan(1.)/45.
      dpld  = 36.d1 - dpd(1) + dpd(2)
      twopi = atan(1.)*8.
*
      fdpd  = 0.
c     write(6,100)
c 100 format(/'   semidiurnal & diurnal tidal potential terms:'/)
      do 201 m=1,2
         fdpd = fdpd + dpld
         do 101 l=1,30
c  the data formerly read here is now in the data statements above
c           read(9,*) (indx(l,m,n),n=1,4), amp(l,m), pha(l,m)
            theta = fdpd
            do 1 n=1,4
               theta = theta + real(indx(l,m,n))*dpd(n)
    1       continue
            freq(l,m) = theta
            omegat = 2.*theta*rad
            ct(l,m) = 2.*cos(omegat)
            st(l,m) = 2.*sin(omegat)
            if (radiation .and. m.eq.2 .and. 
     *           (l.eq.21.or.l.eq.22.or.l.eq.23)) then
               amp(l,m) = amp(l,m)*rad_factor
               pha(l,m) = pha(l,m) - rad_phase
               write(6,503) l,m,(indx(l,m,n),n=1,4),freq(l,m),amp(l,m),
     *                      pha(l,m)
            else
c              write(6,502) l,m,(indx(l,m,n),n=1,4),freq(l,m),amp(l,m),
c    *                      pha(l,m)
            endif
  101    continue
c        write(6,*) ' '
  201 continue
  502 format(1x,i3,i4,2x,4i2,f12.6,2(f10.5,f10.1,5x))
  503 format(1x,i3,i4,2x,4i2,f12.6,f10.5,f10.1,3x,'incl.rad.')
      endif

      if (ts.ne.tslast) then
         td = (ts - tc)/864.d2
*    compute 4 principal mean longitudes for a given time td
         do 107 n=1,4
            ph = phc(n) + td*dpd(n)
            shpn(n) = mod(ph,360.0)
  107    continue
         fd = td - int(td)
         e  = 36.d1*fd - shpn(1) + shpn(2)
      endif
*
      otide = 0.
      do 300 m=1,2
         if (ts.ne.tslast) then
*
*  compute orthotides p(3),q(3) for given time ts in terms
*  of potential amplitudes at(3), bt(3)
            do 207 k=1,3
               at(k) = 0.
               bt(k) = 0.
  207       continue
            do 407 l=1,30
               ph = m*e + pha(l,m)
               do 307 n=1,4
                  i = indx(l,m,n)
                  if (i.eq.0) go to 307
                  ph = ph + real(i)*shpn(n)
  307          continue
               theta = ph*rad
               a = amp(l,m)*100.
               cr = a*cos(theta)
               sr = a*sin(theta)
               crct = cr*ct(l,m)
               crst = cr*st(l,m)
               srct = sr*ct(l,m)
               srst = sr*st(l,m)
               at(1) = at(1) + cr
               bt(1) = bt(1) - sr
               at(2) = at(2) + crct
               bt(2) = bt(2) - srct
               at(3) = at(3) + crst
               bt(3) = bt(3) - srst
  407       continue
            p(1,m) = u00(m)*at(1)
            q(1,m) = u00(m)*bt(1)
            p(2,m) = u20(m)*at(1) - u21(m)*at(2)
            q(2,m) = u20(m)*bt(1) - u21(m)*bt(2)
            if (pseudo) then
*  cartwright and ray pseudo-orthoweights
               p(3,m) = u40(m)*at(1) - u41(m)*at(2) - v41(m)*at(3)
               q(3,m) = u40(m)*bt(1) - u41(m)*bt(2) - v41(m)*bt(3)
            else
*  correct orthoweight definition
               p(3,m) = u40(m)*at(1) - u41(m)*at(2) + v41(m)*at(3)
               q(3,m) = u40(m)*bt(1) - u41(m)*bt(2) + v41(m)*bt(3)
            endif
         endif
         do 30 k=1,3
            otide = otide + u(k,m)*p(k,m) + v(k,m)*q(k,m)
   30    continue
  300 continue
      tide = otide
      tslast = ts
      end subroutine
      end module
