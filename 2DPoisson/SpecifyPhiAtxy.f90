
 SUBROUTINE SpecifyPhiAtxy( x, y, PhiAtxy )
 USE Kinds, ONLY : idknd
 IMPLICIT NONE
 REAL(KIND= idknd) :: r,rs,PI
 REAL(KIND= idknd), INTENT(IN) :: x, y
 REAL(KIND= idknd), INTENT(OUT) :: PhiAtxy
    PI = 4.d0*atan(1.d0)
    r = sqrt(x**2+y**2)
    rs = 1e-3
    PhiAtxy= (1/(2.d0*PI))*log(abs(rs-r))

 END SUBROUTINE SpecifyPhiAtxy
