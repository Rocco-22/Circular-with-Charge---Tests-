
!2345678901234567890123456789012345678901234567890123456789012345678901234567890

 SUBROUTINE LSZEIC( xr, yz, iFlag )
 USE Kinds
 USE ModulePlots, ONLY : LineColor

 IMPLICIT NONE

 REAL(KIND= idknd),    INTENT(IN) :: xr, yz
 INTEGER(KIND= iiknd), INTENT(IN) :: iFlag

 REAL(KIND= idknd) :: xA, yA, rcv
 INTEGER :: ifz
 common /clszc/ xA(1:3), yA(1:3), rcv, ifz
 save   /clszc/

 INTEGER(KIND= iiknd) :: il, ixr, iyz
 REAL(KIND= idknd) :: fxr, fyz, x, y

!c      WRITE (UNIT= *,*) 'lszeic:',xr,yz,iFlag
!!    IH= iFlag/10
    il= MOD(iFlag,10_iiknd)
    IF (il > 5) THEN
       RETURN
    END IF

    ixr= INT(xr, KIND= iiknd)
    fxr= xr-ixr
    iyz= INT(YZ, KIND= iiknd)
    fyz= YZ-iyz
    x= (1-fxr)*xA(ixr) + fxr*xA(ixr+1)
    y= (1-fyz)*yA(iyz) + fyz*yA(iyz+1)

    IF ((il /= 2_iiknd) .AND. (il /= 3_iiknd)) THEN
       !
       ! Zeichnen einer Verbindungslinie zum Zielpunkt
       ! Line to the Target.
       !
       WRITE (UNIT= ifz, FMT= '(1P,3(E14.6))') x, y, rcv
    ELSE
       !
       ! Start einer Konturlinie
       ! Start of Contourline.
       !
       ! Positionieren auf den neuen Zielpunkt
       ! Positioning to the next Target.
       !
       WRITE (UNIT= ifz, FMT= 72) LineColor, x, y, rcv
    END IF
 72 FORMAT( /, &
            ' % linecolor= ', I0, /, &
            1P, 3(E14.6) )

 END SUBROUTINE LSZEIC

