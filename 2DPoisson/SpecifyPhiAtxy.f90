
 SUBROUTINE SpecifyPhiAtxy( x, y, PhiAtxy )
 USE Kinds, ONLY : idknd
 IMPLICIT NONE

 REAL(KIND= idknd), INTENT(IN) :: x, y
 REAL(KIND= idknd), INTENT(OUT) :: PhiAtxy

    PhiAtxy= 3.d0

 END SUBROUTINE SpecifyPhiAtxy
