module m_get_erfc
contains

FUNCTION get_ERFC (X)
!DEC$ ATTRIBUTES NO_ARG_CHECK :: X   
!  erfc = 1.0 - erf (), erf etant la fonction d'erreur
!       definie, si 'erfci' est la fonction inverse, par :
!       0.5 * erfc [ - erfci(p) / 1.414 ] = p
!  Calculs dus a mon collegue J.Labeyrie
REAL, INTENT (IN) :: X
REAL              :: ERFC
! __________________________________________________
!
!  Poids pour l'evaluation quand  0 < abs(dy) < .477
!
    DOUBLE PRECISION, PARAMETER, DIMENSION (5) :: DP0 = (/ &
         113.8641541510502D0,    377.485237685302D0,       &
        3209.377589138469D0,       0.1857777061846032D0,   &
           3.161123743870566D0                          /)
    DOUBLE PRECISION, PARAMETER, DIMENSION (4) :: DQ0 = (/ &
         244.0246379344442D0,   1282.616526077372D0,       &
        2844.236833439171D0,      23.60129095234412D0   /)
!
!              ..........  quand  .477 <= abs(dy) < 4.0
!
    DOUBLE PRECISION, PARAMETER, DIMENSION (9) :: DP1 = (/ &
           8.883149794388376D0,   66.11919063714163D0,     &
         298.6351381974001D0,    881.9522212417691D0,      &
        1712.047612634071D0,    2051.078377826071D0,       &
        1230.339354797997D0,       2.153115354744038D-8,   &
           0.5641884969886701D0                         /)
    DOUBLE PRECISION, PARAMETER, DIMENSION (8) :: DQ1 = (/ &
         117.6939508913125D0,    537.1811018620099D0,      &
        1621.38957456669D0,     3290.79923573346D0,        &
        4362.619090143247D0,    3439.367674143722D0,       &
        1230.339354803749D0,      15.74492611070983D0   /)
!
!              ..........  quand  4.0 <= abs(dy)
!
    DOUBLE PRECISION, PARAMETER, DIMENSION (6) :: DP2 = (/ &
          -3.603448999498044D-01, -1.257817261112292D-01,  &
          -1.608378514874228D-02, -6.587491615298378D-04,  &
          -1.631538713730210D-02, -3.053266349612323D-01/)
    DOUBLE PRECISION, PARAMETER, DIMENSION (5) :: DQ2 = (/ &
           1.87295284992346D0,     5.279051029514284D-01,  &
           6.051834131244132D-02,  2.335204976268692D-03,  &
           2.568520192289822D0                          /)
!
    DOUBLE PRECISION, PARAMETER :: DPI = 0.5641895835477563D0
    DOUBLE PRECISION, PARAMETER :: DBIG  = 9.0D0
    DOUBLE PRECISION, PARAMETER :: DLARG = 6.375D0
    DOUBLE PRECISION, PARAMETER :: DSMAL = 1.0D-10
    DOUBLE PRECISION :: DX, DINV, DRES, DSQ, DNUM, DDEN
!
!  Traitement du signe
!
      IF (X > 0.0) THEN
          ISW = 1
          DX = DBLE (X)
      ELSE
          ISW = -1
          DX = -DBLE (X)
      ENDIF
!
!  Integration par points de Gauss adaptes a l'intervalle
!
      IF (DX < 0.477D0) THEN
!
          IF (DX < DSMAL) THEN
              DRES = DX * DP0(3) / DQ0(3)
          ELSE
              DSQ = DX * DX
              DNUM = DP0(4) * DSQ + DP0(5)
              DDEN = DSQ + DQ0(4)
              DO I = 1, 3
                  DNUM = DNUM * DSQ + DP0(I)
                  DDEN = DDEN * DSQ + DQ0(I)
              ENDDO
              DRES = DX * DNUM / DDEN
          ENDIF
          IF (ISW < 0) DRES = - DRES
          DRES = 1.0D0 - DRES
!
      ELSEIF (DX < 4.0D0) THEN
!
          DNUM = DP1(8) * DX + DP1(9)
          DDEN = DX * DQ1(8)
          DO I = 1, 7
              DNUM = DNUM * DX + DP1(I)
              DDEN = DDEN * DX + DQ1(I)
          ENDDO
          DRES = DNUM / DDEN
          DRES = DRES * DEXP ( -DX*DX)
          IF (ISW < 0) DRES = 2.0D0 - DRES
!
      ELSEIF (ISW > 0 .AND. DX > DBIG) THEN
!
          DRES = 0.0D0
!
      ELSEIF (ISW < 0 .AND. DX >= DLARG) THEN
!
          DRES = 2.0D0
!
      ELSE
!
          DSQ = DX * DX
          DINV = 1.0D0 / DSQ
          DNUM = DP2(5) * DINV + DP2(6)
          DDEN = DINV + DQ2(5)
          DO I = 1, 4
              DNUM = DNUM * DINV + DP2(I)
              DDEN = DDEN * DINV + DQ2(I)
          ENDDO
          DRES = (DPI + DINV * DNUM / DDEN) / DX
          DRES = DRES * DEXP ( -DSQ)
          IF (ISW < 0) DRES = 2.0D0 - DRES
!
      ENDIF
!
      get_ERFC = DRES
      RETURN
END FUNCTION get_ERFC

end module m_get_erfc
