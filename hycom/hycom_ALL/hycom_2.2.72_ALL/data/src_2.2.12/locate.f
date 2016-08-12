      FUNCTION locate(xx,x,IDM)
      IMPLICIT NONE

      REAL, DIMENSION(IDM),INTENT(IN) :: xx
      REAL :: x
      INTEGER :: locate,IDM
      INTEGER :: n,jl,jm,ju
      LOGICAL :: ascnd
      
      n=size(xx)
      ascnd =(xx(n)>=xx(1))
      jl=0
      ju=n+1

      do
      if (ju-jl <=1) exit

      jm=(ju+jl)/2

      if (ascnd .eqv. (x>=xx(jm))) then
        jl=jm
      else
        ju=jm
      end if
      enddo
      if (x==xx(1)) then
        locate=1
      else if (x = =xx(n)) then
        locate=n-1
      else
        locate=jl
        end if

      end FUNCTION locate

      FUNCTION concat(s1, s2)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: s1, s2
      CHARACTER(LEN=LEN_TRIM(s1)+LEN_TRIM(s2)) :: concat ! func name
      concat = TRIM(s1) // TRIM(s2)
      END FUNCTION concat
             
