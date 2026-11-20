
 SUBROUTINE SpecifyPhiAtxy( x, y, PhiAtxy )
 USE Kinds, ONLY : idknd
 IMPLICIT NONE
 REAL(KIND= idknd) :: r,xs,ys,PI
 REAL(KIND= idknd), INTENT(IN) :: x, y
 REAL(KIND= idknd), INTENT(OUT) :: PhiAtxy
    PI = 4.d0*atan(1.d0)
    xs = 0d0; ys = 5d-3
    r = sqrt((x-xs)**2+(y-ys)**2)
    PhiAtxy= (1/(2.d0*PI))*log(abs(r))

 END SUBROUTINE SpecifyPhiAtxy
