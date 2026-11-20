
 INTEGER FUNCTION NrofThreads()
    NrofThreads= 20                    ! Adjust this.
 END FUNCTION NrofThreads
!------
 MODULE InputData
 USE Kinds, ONLY : idknd
 IMPLICIT NONE

 PUBLIC
    REAL(KIND= idknd), ALLOCATABLE, DIMENSION(:) :: xPosition, yPosition
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Material

    INTEGER :: NrOfCharges
    REAL(KIND= idknd), DIMENSION(1:100) :: xOfCharge, yOfCharge, ChargeOfCharge

    REAL, DIMENSION(0:10) :: Epsr= (/ 1.0, &    !  Epsr(0)
                                      1e20, &   !  Epsr(1)   ! Metal
                                      1e-10, &  !  Epsr(2)   ! Magnetic.
                                      1.0, &    !  Epsr(3)
                                      1e20, &   !  Epsr(4)  ! Metal
                                      1.0, &    !  Epsr(5)
                                      1.0, &    !  Epsr(6)
                                      1.0, &    !  Epsr(7)
                                      1.0, &    !  Epsr(8)
                                      1.0, &    !  Epsr(9)
                                      1.0 /)    !  Epsr(10)
 END MODULE InputData
!------
 SUBROUTINE BuildAndSolve2DPoisson
 USE Kinds, ONLY : idknd
 USE Multigrid, ONLY : &
       PoissonGrid2D, &
       MGAllocateData, &
       MGDeAllocateData, &
       MGData_FirstIx_etc, &
       Multigrid_2D_Type2, &
       MGDebugPlots

 USE InputData, ONLY : &
       xPosition, yPosition, &
       Material, Epsr, &
       NrOfCharges, xOfCharge, yOfCharge, ChargeOfCharge

 USE ModulePlots, ONLY : &
       WritePotentialAsMTVFile

 IMPLICIT NONE

 CHARACTER(LEN= *), PARAMETER :: MyName= 'BuildAndSolve2DPoisson'

 TYPE (PoissonGrid2D) :: Data
 REAL(KIND= idknd), ALLOCATABLE, DIMENSION(:,:) :: Phi
 REAL(KIND= idknd) :: ChargeSquare, LastChargeSquare, ChargeSquareMaximum, rDum
 REAL :: Omega, Omega2
 INTEGER :: i, iter, NIter, ix, iy, ix1, iy1, ix2, iy2
 REAL(KIND= idknd) :: Factor1, Factor2, Factor3, Factor4, Summe, x, y, PhiAtxy
 INTEGER :: ThisNThreads
 INTEGER, EXTERNAL :: NrOfThreads

 integer :: i11

 MGDebugPlots= .false. ! No Debug-Data
       !
       ! Aufstellen der Datastruktur fuer MultiGrid.
       !
       ix1= 0-1
       iy1= 0-1
       ix2= UBOUND( xPosition, DIM= 1) - 1
       iy2= UBOUND( yPosition, DIM= 1) - 1
       CALL MGAllocateData( MyName, Data, ix1, ix2, iy1, iy2 )
!! write (*,*) MyName, ': Passed MGAllocateData.'; call wbflsh( 6 )

       DO ix= LBOUND(Data% dx, DIM= 1)+1, UBOUND(Data% dx, DIM= 1)
          Data% dx(ix)= xPosition(ix+1) - xPosition(ix)
       END DO
       ix= LBOUND(Data% dx, DIM= 1)
       Data% dx(ix)= Data% dx(ix+1)

       DO iy= LBOUND(Data% dy, DIM= 1)+1, UBOUND(Data% dy, DIM= 1)
          Data% dy(iy)= yPosition(iy+1) -yPosition(iy)
       END DO
       iy= LBOUND(Data% dy, DIM= 1)
       Data% dy(iy)= Data% dy(iy+1)
    rDum= 1
       Data% nVariable= 1
       DO iy= LBOUND(Data% Factor, DIM= 3)+1, UBOUND(Data% Factor, DIM= 3)-1
       DO ix= LBOUND(Data% Factor, DIM= 2)+1, UBOUND(Data% Factor, DIM= 2)-1
          IF (dsoEpsA(3,ix,iy) < 1e-5) THEN
             Data% Variable(ix,iy)= .FALSE.
             x= xPosition(ix) - Data% dx(ix) / 2
             y= yPosition(iy) - Data% dy(iy) / 2
             CALL SpecifyPhiAtxy( x, y, PhiAtxy )
             Data% Phi(ix,iy)= PhiAtxy
          ELSE
             Factor1= 0
             Factor2= 0
             Factor3= 0
             Factor4= 0
             IF (dsoEpsA(1,ix,iy) > 1e-4) THEN
                Factor1= rDum / dsoEpsA(1,ix,iy)
                Data% EpsA(1,ix,iy)= Factor1 * Data% dx(ix)
             END IF
             IF (dsoEpsA(2,ix,iy) > 1e-4) THEN
                Factor2= rDum / dsoEpsA(2,ix,iy)
                Data% EpsA(2,ix,iy)= Factor2 * Data% dy(iy)
             END IF
             IF (dsoEpsA(1,ix-1,iy) > 1e-4) THEN
                Factor3= rDum / dsoEpsA(1,ix-1,iy)
                Data% EpsA(1,ix-1,iy)= Factor3 * Data% dx(ix-1)
             END IF
             IF (dsoEpsA(2,ix,iy-1) > 1e-4) THEN
                Factor4= rDum / dsoEpsA(2,ix,iy-1)
                Data% EpsA(2,ix,iy-1)= Factor4 * Data% dy(iy-1)
             END IF
             Summe= Factor1 + Factor2 + Factor3 + Factor4
             IF (Summe /= 0) THEN
                Data% Variable(ix,iy)= .TRUE.
                Data% Factor(1,ix,iy)= Factor1 / Summe
                Data% Factor(2,ix,iy)= Factor2 / Summe
                Data% Factor(3,ix,iy)= Factor3 / Summe
                Data% Factor(4,ix,iy)= Factor4 / Summe
         !!       Data% Factor(5,ix,iy)= 1 / Summe
                Data% ooSumme(ix,iy)= 1 / Summe
                Data% EpsAods(1,ix  ,iy  )= Factor1
                Data% EpsAods(2,ix  ,iy  )= Factor2
                Data% EpsAods(1,ix-1,iy  )= Factor3
                Data% EpsAods(2,ix  ,iy-1)= Factor4
                Data% nVariable= Data% nVariable + 1
             ELSE
                Data% Variable(ix,iy)= .FALSE.
                Data% ooSumme(ix,iy)= 0
             END IF
          END IF
       END DO
       END DO

       !
       ! Clean up somewhat.
       !
       DO iy= LBOUND(Data% Variable, DIM= 2)+1, UBOUND(Data% Variable, DIM= 2)-1
       DO ix= LBOUND(Data% Variable, DIM= 1)+1, UBOUND(Data% Variable, DIM= 1)-1
          IF (.NOT. Data% Variable(ix,iy)) THEN
             IF (.NOT. Data% Variable(ix-1,iy)) THEN
                IF (.NOT. Data% Variable(ix+1,iy)) THEN
                   IF (.NOT. Data% Variable(ix,iy-1)) THEN
                      IF (.NOT. Data% Variable(ix,iy+1)) THEN
                         Data% Phi(ix,iy)= 0
                      END IF
                   END IF
                END IF
             END IF
          END IF
       END DO
       END DO

 if (0 == 1) then ! Debug.
 i11= 11
 OPEN( UNIT= i11, FILE= 'Variable' )
       WRITE (i11,'(A)') ' % toplabel= Variable"'
       WRITE (i11,*) ' $DATA= CONTOUR'
       WRITE (i11,*) ' % nx=', SIZE(Data% Variable, DIM= 1)-2
       WRITE (i11,*) ' % ny=', SIZE(Data% Variable, DIM= 2)-2
       DO iy= LBOUND(Data% Variable, DIM= 2)+1, UBOUND(Data% Variable, DIM= 2)-1, 1
       DO ix= LBOUND(Data% Variable, DIM= 1)+1, UBOUND(Data% Variable, DIM= 1)-1, 1
          IF (Data% Variable(ix,iy)) THEN
             WRITE (i11,*) 1
          ELSE
             WRITE (i11,*) 0
          END IF
       END DO
       END DO
       CLOSE (UNIT= i11, STATUS= 'KEEP')
 end if

!! write (*,*) MyName, ': Vor MGData_FirstIx_etc.'; call wbflsh( 6 )
       CALL MGData_FirstIx_etc( Who= MyName, D= Data )

       !
       ! The Charges.
       !
       Data% Charge= 0
       DO i= 1, NrOfCharges
          CALL FindThePosition( xOfCharge(i), xPosition, ix )
          CALL FindThePosition( yOfCharge(i), yPosition, iy )
          Data% Charge(ix,iy)= ChargeOfCharge(i)
          Data% ChargeTerm(ix,iy)= Data% Charge(ix,iy) * Data% ooSumme(ix,iy)
       END DO

  NIter= 10; Omega= 1.60; Omega2= 1.10 ! ?? Iterations
!! write (*,*) MyName, ': nVariable:', Data% nVariable; call wbflsh(6) ! 1435801
       IF (Data% nVariable > 1e5) THEN
          ! 2011-01-24
          Omega= 1.93; Omega2= 1.60   ! 300, 2e-23
          ThisNThreads= MIN( 8, NrOfThreads() )
       ELSE IF (Data% nVariable > 1e4) THEN
          ThisNThreads= MIN( 4, NrOfThreads() )
       ELSE IF (Data% nVariable > 1e3) THEN
          ThisNThreads= MIN( 3, NrOfThreads() )
       ELSE
          ThisNThreads= MIN( 3, NrOfThreads() )
       END IF
! Mit N Threads: 1 Sekunde
! ThisNThreads= 1 !! Open64, 7 Sekunden.
       LastChargeSquare= HUGE( LastChargeSquare )
       ChargeSquareMaximum= 0
       DO iter= 1, 1001
!! write (*,*) MyName, ': Vor MG.'; call wbflsh( 6 )
          CALL Multigrid_2D_Type2 &
                         ( MyName, Data, &
                           NIter= NIter, Omega= Omega, Omega2= Omega2, &
                           iPlane= 0, Iteration= iter, &
                           ChargeSquare= ChargeSquare, &
                           nThreads= ThisNThreads )
!! write (*,*) MyName, ': Nach MG.'
 if ((mod(iter,10) == 0) .or. (iter<10)) then
    IF (ChargeSquareMaximum /= 0) THEN
       WRITE (UNIT= *, FMT= '(A,I5,1P,2E9.2)') 'i, Residual:', &
         iter, ChargeSquare/ChargeSquareMaximum, ChargeSquare
    END IF
 end if
          ChargeSquareMaximum= MAX( ChargeSquareMaximum, &
                                         ChargeSquare )
          IF (Iter > 5) THEN
             IF (ChargeSquare < 1e-20 * ChargeSquareMaximum * Data% nVariable) THEN
                EXIT
             END IF
          END IF
       END DO

       IF (iter > 1000) THEN
          WRITE (UNIT= *, FMT= '(2A,I5,1P,E12.4)') &
                    MyName, ': i, Residuum:', &
                    iter, ChargeSquare/ChargeSquareMaximum
       END IF

 if (1==1) then
       DO iter= iter+1, iter+10
 write (*,*) MyName, ': Again with large NIter.'
          LastChargeSquare= ChargeSquare
          NIter= 100; Omega= 1.5; Omega2= 1
          CALL Multigrid_2D_Type2 &
                         ( MyName, Data, &
                           NIter= NIter, Omega= Omega, Omega2= Omega2, &
                           iPlane= 0, Iteration= iter, &
                           ChargeSquare= ChargeSquare, &
                           nThreads= ThisNThreads )
 write (*,*) MyName, ': iter:', &
                        iter, ChargeSquare/ChargeSquareMaximum
          IF (ChargeSquare > 0.8 * LastChargeSquare) THEN
             EXIT
          END IF

       END DO
 end if

       ix1= LBOUND( Data% Phi, DIM= 1)+1
       ix2= UBOUND( Data% Phi, DIM= 1)
       iy1= LBOUND( Data% Phi, DIM= 2)+1
       iy2= UBOUND( Data% Phi, DIM= 2)
       ALLOCATE( Phi(ix1:ix2, iy1:iy2) )
       DO iy= iy1, iy2
          DO ix= ix1, ix2
             Phi(ix,iy)= Data% Phi(ix,iy)
          END DO
       END DO
       CALL WritePotentialAsMTVFile( 'Potential', &
                                     'Potential.mtv', &
                                     Phi, xPosition, yPosition )
       CALL WritePhiData( 'PhiData', Phi, xPosition, yPosition )
       DEALLOCATE( Phi )

       CALL MGDeAllocateData( Data )
!----
 CONTAINS
     REAL(KIND= idknd) FUNCTION dsoEpsA(ic, ix, iy)
     INTEGER, INTENT(IN) :: ic, ix, iy
     ! dz is '1'.
        IF (ic == 1) THEN
           dsoEpsA= dx(ix) &
                / (( dy(iy  ) * Epsr(Material(ix  ,iy  )) &
                   + dy(iy-1) * Epsr(Material(ix  ,iy-1)) ) * 0.5 )
        ELSE IF (ic == 2) THEN
           dsoEpsA= dy(iy) &
                / (( dx(ix  ) * Epsr(Material(ix  ,iy  )) &
                   + dx(ix-1) * Epsr(Material(ix-1,iy  )) ) * 0.5 )
        ELSE IF (ic == 3) THEN ! This is essentially just a Flag.
           dsoEpsA= 4 / ( Epsr(Material(ix  ,iy  )) &
                        + Epsr(Material(ix-1,iy  )) &
                        + Epsr(Material(ix-1,iy-1)) &
                        + Epsr(Material(ix  ,iy-1)) )
        END IF
     END FUNCTION dsoEpsA
!-----
     REAL(KIND= idknd) FUNCTION dx(i)
     INTEGER, INTENT(IN) :: i
        dx= xPosition(i-1) - xPosition(i)
     END FUNCTION dx
!-----
     REAL(KIND= idknd) FUNCTION dy(i)
     INTEGER, INTENT(IN) :: i
        dy= yPosition(i-1) - yPosition(i)
     END FUNCTION dy
!-----
     SUBROUTINE FindThePosition( FindThis, InThisArray, iPosition )
     REAL(KIND= idknd), INTENT(IN) :: FindThis
     REAL(KIND= idknd), DIMENSION(:), INTENT(IN) :: InThisArray
     INTEGER, INTENT(OUT) :: iPosition

     REAL(KIND= idknd) :: SmallestDistance, ThisDistance
     INTEGER :: i

        !
        ! Linear Search.
        !
        SmallestDistance= HUGE( SmallestDistance )
        DO i= LBOUND(InThisArray, DIM= 1), UBOUND(InThisArray, DIM= 1)
           ThisDistance= ABS( FindThis - InThisArray(i) )
           IF (ThisDistance <= SmallestDistance) THEN
              iPosition= i
              SmallestDistance= ThisDistance
           END IF
        END DO

     END SUBROUTINE FindThePosition
!-----
     SUBROUTINE WritePhiData( FileName, Phi, xPosition, yPosition )
     CHARACTER(LEN= *),                 INTENT(IN) :: FileName
     REAL(KIND= idknd), DIMENSION(:,:), INTENT(IN) :: Phi
     REAL(KIND= idknd), DIMENSION(:),   INTENT(IN) :: xPosition, yPosition

     CHARACTER(LEN= *), PARAMETER :: MyName= 'WritePhiData'
     INTEGER :: ix, iy, i11, i

       i11= 11
       OPEN(UNIT= i11, &
            FILE= FileName, &
          STATUS= 'UNKNOWN' )
       WRITE (UNIT= i11, FMT= 71) &
               ( LBOUND(Phi, DIM= i), UBOUND(PHi, DIM= i), i= 1, 2 )
 71    FORMAT( 4I9, ' <- nx0, nx, ny0, ny', /, &
               '# Gridplanes nx0..nx, ny0..ny are as follows:' )
 72    FORMAT( 1P, E14.7, 0P, I9, A )
       DO i= LBOUND(Phi, DIM= 1), UBOUND(Phi, DIM= 1)
          WRITE (UNIT= i11, FMT= 72) xPosition(i), i, ' <- x(i), i'
       END DO
       DO i= LBOUND(Phi, DIM= 2), UBOUND(Phi, DIM= 2)
          WRITE (UNIT= i11, FMT= 72) yPosition(i), i, ' <- y(i), i'
       END DO

       WRITE (UNIT= i11, FMT= 73)
    73 FORMAT( '# The Potential Data is written by a loop nest like:', /, &
               '# DO iy= ny0, ny', /, &
               '#    DO ix= nx0, nx', /, &
               '#       WRITE (i11,*) Phi(ix,iy)', /, &
               '#    END DO', /, &
               '# END DO' )

      DO iy= LBOUND(Phi, DIM= 2), UBOUND(Phi, DIM= 2)
         DO ix= LBOUND(Phi, DIM= 1), UBOUND(Phi, DIM= 1)
            WRITE (i11,*) Phi(ix,iy)
         END DO
      END DO

!! write (*,*) MyName, ': before CLOSE.'
       CLOSE( UNIT= i11, STATUS= 'KEEP')

     END SUBROUTINE WritePhiData
 END SUBROUTINE BuildAndSolve2DPoisson
!***********
 SUBROUTINE ReadChargeData
   USE InputData, ONLY : &
         NrOfCharges, xOfCharge, yOfCharge, ChargeOfCharge

   IMPLICIT NONE
   CHARACTER(LEN= *), PARAMETER :: MyName= 'ReadChargeData'

   INTEGER :: i11, i

      i11= 11
      OPEN(UNIT= i11, &
            FILE= 'ChargeData', &
         STATUS= 'OLD' )

      READ (UNIT= i11, FMT= *) NrOfCharges
      DO i= 1, NrOfCharges
         READ (UNIT= i11, FMT= *) xOfCharge(i), yOfCharge(i), ChargeOfCharge(i)
      END DO

   !! write (*,*) MyName, ': before CLOSE.'
      CLOSE( UNIT= i11, STATUS= 'KEEP')

 END SUBROUTINE ReadChargeData
!***********
 SUBROUTINE ReadEpsData
 USE Kinds, ONLY : idknd
 USE InputData, ONLY : &
        Epsr

 IMPLICIT NONE
 CHARACTER(LEN= *), PARAMETER :: MyName= 'ReadEpsData'

 INTEGER :: i11, i, NrOfMaterials, MaterialNumber
 REAL(KIND= idknd) :: ThisEps

    i11= 11
    OPEN(UNIT= i11, &
         FILE= 'EpsData', &
       STATUS= 'OLD' )

    READ (UNIT= i11, FMT= *) NrOfMaterials
    DO i= 1, NrOfMaterials
       READ (UNIT= i11, FMT= *) MaterialNumber, ThisEps
       Epsr(MaterialNumber)= ThisEps
    END DO

!! write (*,*) MyName, ': before CLOSE.'
    CLOSE( UNIT= i11, STATUS= 'KEEP')

 END SUBROUTINE ReadEpsData
!***********
 PROGRAM ReadInputDataAndSolve
 USE InputData, ONLY : &
     xPosition, yPosition, Material

 IMPLICIT NONE

 CHARACTER(LEN= *), PARAMETER :: MyName= 'ReadInputDataAndSolve'
 INTEGER :: nx, ny, nz, ix, iy, iz, ix0, jy, jz
 INTEGER :: i11, i, nPop, jMat
 INTEGER :: MinMat, MaxMat
 REAL, ALLOCATABLE, DIMENSION(:) :: z

 i11= 11
    OPEN(UNIT= i11, &
         FILE= 'Gitterfuellung', &
       STATUS= 'OLD' )

    READ (UNIT= i11, FMT= *) nx, ny, nz
    READ (UNIT= i11, FMT= *)            ! Skip a Line

    ALLOCATE(xPosition(0:nx+1), yPosition(0:ny+1), z(0:nz+1))
    ALLOCATE(Material(-1:nx, -1:ny))
    Material= 0

    ! Read the Positions.
    DO ix= 0, nx+1
       READ (UNIT= i11, FMT= *) xPosition(ix)
    END DO
    DO iy= 0, ny+1
       READ (UNIT= i11, FMT= *) yPosition(iy)
    END DO
    DO iz= 0, nz+1
       READ (UNIT= i11, FMT= *) z(iz)  ! Never used.
    END DO

 if (0==1) then
    ! Skip some Lines.
    DO i= 1, 15
       READ (UNIT= i11, FMT= *)
    END DO
 end if

    ! Read the nPop and jMat, fill the Material-Array at iz= nz/2.
    DO iz= 0, nz/2
       MinMat= 10000
       MaxMat= -1
       DO iy= 0, ny
          ix= 0
          DO WHILE(ix < nx)
             READ (UNIT= i11, FMT= *) nPop, jMat ! , ix0, jy, jz
 if (0==1) then
 write (*,*) MyName, nPop, jMat, ix0, jy, jz
 if (ix0 /= ix) then
    write (*,*) ' Bad: ix, ix0, jy, jz:', ix, ix0, jy, jz
    stop
 end if
 if (jy /= iy) then
    write (*,*) ' Bad: iy, jy:', iy, jy, ix0, jz
    stop
 end if
 if (jz /= iz) then
    write (*,*) ' Bad: iz, jz:', iz, jz, ix0, jy
    stop
 end if
 end if
             Material(ix:ix+nPop-1,iy)= jMat
             ix= ix + nPop
 MinMat= MIN(MinMat, jMat)
 MaxMat= MAX(MaxMat, jMat)
          END DO
       END DO
 write (*,*) ' iz, MinMat, MaxMat:', iz, MinMat, MaxMat
!!       IF (MinMat /= MaxMat) EXIT
    END DO
 write (*,*) ' iz, MinMat, MaxMat:', iz, MinMat, MaxMat

!! write (*,*) MyName, ': before CLOSE.'
    CLOSE (UNIT= i11, STATUS= 'KEEP' )

    CALL WriteMaterialAsMTVFile  ! Not needed.

    CALL ReadEpsData
    CALL ReadChargeData
    CALL BuildAndSolve2DPoisson


 CONTAINS
 SUBROUTINE WriteMaterialAsMTVFile
 INTEGER :: ix, iy, iz, Mat, iOldMat
 INTEGER, PARAMETER :: mtv= 11

    OPEN (UNIT  = mtv, &
          FILE  = 'Material.mtv',  &
          ACCESS= 'SEQUENTIAL',  &
          FORM  = 'FORMATTED',  &
          ACTION= 'WRITE',  &
          STATUS= 'UNKNOWN')
    WRITE (UNIT= mtv, FMT= 72)
 72 FORMAT(  &
       ' $ DATA= CURVE3D', /,  &
       ' % equalscale= True', /,  &
       ' % axisscale= False', /,  &
       ' % hiddenline= True', /,  &
       ' % filltype= 1', /,  &
       ' % toplabel= "Material"', /,  &
       ' % leftworld= False', /,  &
       ' % comment= "No Comment"', / )
 if (1==1) then
    WRITE (UNIT= mtv, FMT= 73) &
         xPosition(1), xPosition(nx), &
         yPosition(1), yPosition(ny), &
         z(nz/2-1), z(nz/2+1)
 73 FORMAT( 1P, &
       ' % xmin= ', E14.6,' xmax= ', E14.6, /,  &
       ' % ymin= ', E14.6,' ymax= ', E14.6, /,  &
       ' % zmin= ', E14.6,' zmax= ', E14.6, / )
 end if

     8010 FORMAT( &
          '$ DATA= CURVE3D Name=Material', /, & 
          '% filltype= 1      # Solid fill for each face', /, &
          '% hiddenline= true # Apply hiddenline/surface sorting', /, &
          '% leftworld= False', /, &
          '% toplabel= "GdfidL, 3D Arrowplot"', /, &
          '% subtitle= "', A, '"', /, &
          '% sidelabel= False', /, &
          '% equalScale= True', /, & 
          '% axisScale= False', /, &
          '% comment= "', A, '"', /, &
          '% xlabel= "x/m", ylabel= "y/m", zlabel= "z/m"', /, &
          '% boundary= False', /, &
          1P, &
          '% eyepos.x=', E12.4, ' eyepos.y=', E12.4, ' eyepos.z=', E12.4, /, &
          '% xmin= ', E14.6, ' xmax= ', E14.6, /, &
          '% ymin= ', E14.6, ' ymax= ', E14.6, /, &
          '% zmin= ', E14.6, ' zmax= ', E14.6, / )

    iz= nz/2
    DO iy= 1, ny
       iOldMat= -1
       DO ix= 1, nx
          Mat= Material(ix,iy)
       !!   IF (Mat /= iOldMat) THEN
          IF (.true.) THEN
             WRITE (UNIT= mtv, FMT= 74) &
                Mat, &
                xPosition(ix  ), yPosition(iy  ), z(iz), &
                xPosition(ix+1), yPosition(iy  ), z(iz), &
                xPosition(ix+1), yPosition(iy+1), z(iz), &
                xPosition(ix  ), yPosition(iy+1), z(iz), &
                xPosition(ix  ), yPosition(iy  ), z(iz)
          END IF
          iOldMat= Mat
       END DO
    END DO
 74 FORMAT( '% fillcolor= ', I0, ' linecolor= 5 filltype= 1', /, &
            1P, 5(3E15.7,/) )

    CLOSE (UNIT  = mtv, &
           STATUS= 'KEEP' )

 END SUBROUTINE WriteMaterialAsMTVFile

 END PROGRAM ReadInputDataAndSolve
! ----
