
 MODULE NaNs_etc
 USE Kinds

 PUBLIC :: IsNan
 INTERFACE IsNan
    MODULE PROCEDURE IsNan_ir4
    MODULE PROCEDURE IsNan_ir8
 END INTERFACE
 
 CONTAINS
 SUBROUTINE IsNan_ir8( x, Flagge )
 REAL(KIND= ir8), INTENT (IN) :: x
 LOGICAL(KIND= ii4), INTENT (OUT) :: Flagge

 INTEGER(KIND= iiknd) :: wbisnan

    Flagge= wbisnan( x ) /= 0

 END SUBROUTINE IsNan_ir8
!*******
 SUBROUTINE IsNan_ir4( x, Flagge )
 REAL(KIND= ir4), INTENT (IN) :: x
 LOGICAL(KIND= ii4), INTENT (OUT) :: Flagge

    CALL IsNan_ir8( REAL(x, KIND= ir8), Flagge)
 END SUBROUTINE IsNan_ir4
!*******
 END MODULE NaNs_etc
