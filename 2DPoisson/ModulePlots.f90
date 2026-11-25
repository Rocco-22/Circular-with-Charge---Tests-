
 MODULE ModulePlots
 USE Kinds, ONLY : iiknd, idknd

 PUBLIC
   CHARACTER(LEN= 400), SAVE :: c2dPlotOpts = '-geometry 690x560+10+10 '

   CHARACTER(LEN= 400), SAVE :: PlotOpts = '-geometry 1024x768+410+40 -noclip'
   REAL(KIND= idknd), DIMENSION(1:3), SAVE :: EyePosition= (/ -1., -2.3, 0.5 /)
   REAL(KIND= idknd), SAVE :: rotx= 0, roty= 0, rotz= 0

   CHARACTER(LEN= 3), SAVE :: OnlyPlotfile= 'no '
 CHARACTER(LEN= 3), SAVE :: plotshowtext= 'yes'

   INTEGER(KIND= iiknd), SAVE :: iPlot= 0

   CHARACTER(LEN= 8), SAVE :: PlotForeground= 'black'
   CHARACTER(LEN= 8), SAVE :: PlotBackground= 'white'
   INTEGER(KIND= iiknd), SAVE :: LineColor= 0_iiknd
                 ! 0: Foreground Color, 1: Yellow, 2: Blue, 3: Green,
                 ! 4: Red, 5: Dark Blue, 6: Orange, 7: Pink, 9: Cyan, 10: Brown
! SUBROUTINE getmtv
!          mymtv(:LEN_TRIM(mymtv)), &
!         '-pfg ', PlotForeground(:LEN_TRIM(PlotForeground)), &
!         '-pbg ', Plotbackground(:LEN_TRIM(Plotbackground)), &


 CONTAINS
!********
 SUBROUTINE WritePotentialAsMTVFile( &
     Title, FileName, &
     Phi, xPosition, yPosition)
 USE Kinds
 IMPLICIT NONE

 CHARACTER(LEN= *),                 INTENT(IN) :: Title, FileName
 REAL(KIND= idknd), DIMENSION(:,:), INTENT(IN) :: Phi
 REAL(KIND= idknd), DIMENSION(:),   INTENT(IN) :: xPosition
 REAL(KIND= idknd), DIMENSION(:),   INTENT(IN) :: yPosition
!
! Aequipotential plot
!
 CHARACTER(LEN= *), PARAMETER :: MyName= 'WritePotentialAsMTVFile'
 REAL(KIND= idknd), DIMENSION(1:2,1:2) :: za
 INTEGER(KIND= iiknd), DIMENSION(1:1000) :: iBitMap
 INTEGER :: ix, iy
 REAL(KIND= idknd) :: x, y, xLow, xHigh, yLow, yHigh, pMin, pMax, Phi0
 REAL(KIND= idknd) :: Scale2D, eFactor, Extension, ooe, rDum, est
 INTEGER :: mtv, NumberOfContourValues, ii, Iter
 INTEGER :: ix1, ix2, iy1, iy2

 REAL(KIND= idknd), DIMENSION(1:3) :: xa, ya
 REAL(KIND= idknd), DIMENSION(1:1) :: rcv
 INTEGER :: ifz
 COMMON /clszc/ xa, ya, rcv, ifz
 SAVE   /clszc/
 EXTERNAL lszeic

!     copy the plot information to the common block /CLSZC/

    iBitMap= 0
    ifz= 14
    OPEN(UNIT= ifz, &
         FILE= FileName, &
         FORM= 'FORMATTED', &
       STATUS= 'UNKNOWN' )

    ix1= LBOUND(Phi, DIM= 1)
    ix2= UBOUND(Phi, DIM= 1)
    iy1= LBOUND(Phi, DIM= 2)
    iy2= UBOUND(Phi, DIM= 2)

    xLow = +1.E+30
    xHigh= -1.E+30
    DO ix= ix1, ix2
       x= xPosition(ix)
       xLow = MIN(x, xLow )
       xHigh= MAX(x, xHigh)
    END DO
!! write (*,*) MyName, ': xLow, xHigh:', xLow, xHigh

    yLow = +1.E+30
    yHigh= -1.E+30
    DO iy= iy1, iy2
       y= yPosition(iy)
       yLow = MIN(y, yLow )
       yHigh= MAX(y, yHigh)
    END DO
!! write (*,*) MyName, ': yLow, yHigh:', yLow, yHigh

    pMin= +1.E+30
    pMax= -1.E+30
    DO iy= iy1, iy2
       DO ix= ix1, ix2
          Phi0= Phi(ix,iy)
          pMin= MIN(phi0, pMin)
          pMax= MAX(phi0, pMax)
       END DO
    END DO

    Extension= MAX(ABS(xHigh-xLow), ABS(yHigh-yLow))
    eFactor= 1
    est= 2
    Iter= 0
 11 CONTINUE
    Iter= Iter + 1
    IF (Iter < 30) THEN
!      write (*,*) ' efakt:',efakt
       IF (eFactor < 0.5*Extension) THEN
          eFactor= eFactor * est
          GOTO 11
       ELSE IF (eFactor > 0.5*Extension) THEN
          eFactor= eFactor / est
          GOTO 11
       END IF
    END IF

    WRITE (UNIT= ifz, FMT= 81) Title
 81 FORMAT( &
       '$ DATA= CURVE3D Name=Potential', /, &
       '% leftworld= false', /, &
       '% filltype= 1     #Solid fill for each face', /, &
       '% hiddenline= true', /, &
       '% toplabel= "', A, '"', /, &
       '% xlabel="x-axis"', /, &
       '% ylabel="y-axis"', /, &
       '% zlabel="Potential"', /, &
       '% comment= "no comment"', /, &
       '% boundary= False', /, &
       '% equalscale= True', /, &
       '% axisscale= False', /, &
       '% window.xmin= -0.4 window.ymin= -0.4', /, &
       '% window.xmax= 0.4 window.ymax= 0.4' &
       )
    WRITE (UNIT= ifz, FMT= *) '% zscale= ', eFactor
    WRITE (UNIT= ifz, FMT= 61) xLow, yLow, pMin
    WRITE (UNIT= ifz, FMT= *)
    WRITE (UNIT= ifz, FMT= 61) xHigh, yHigh, pMax
 61 FORMAT(1P,3E20.7)

    NumberOfContourValues= 30
    Scale2D= 0.3
    ooe= Scale2D * eFactor / ( pMax - pMin )
    DO iy= iy1, iy2-1
       yA(1)= yPosition(iy  )
       yA(2)= yPosition(iy+1)
       yA(3)= yPosition(iy+1) + yA(2) - yA(1)
       DO ix= ix1, ix2-1
          xA(1)= xPosition(ix  )
          xA(2)= xPosition(ix+1)
          xA(3)= xPosition(ix+1) + xA(2) - xA(1)

          za(1,1)= Phi(ix  , iy  ) * ooe
          za(1,2)= Phi(ix  , iy+1) * ooe
          za(2,1)= Phi(ix+1, iy  ) * ooe
          za(2,2)= Phi(ix+1, iy+1) * ooe
          DO ii= 1, NumberOfContourValues
             rcv(1)= ooe * ( pMin + (ii-1) * (pMax-pMin) &
                           / FLOAT(NumberOfContourValues-1) )
             rdum= 2*pMax * ooe
             CALL gcontr( za, 2, 2, 2, &
                          rcv, 1, &
                          rDum, iBitMap, lszeic )
          END DO
       END DO
    END DO

    CLOSE (UNIT= ifz, STATUS= 'KEEP')

 END SUBROUTINE WritePotentialAsMTVFile
!******
 END MODULE ModulePlots
