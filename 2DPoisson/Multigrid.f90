!2345678901234567890123456789012345678901234567890123456789012345678901234567890

!
! Multigrid.fParallel
!

 MODULE Multigrid
 USE Kinds, ONLY : &
        idknd, SparseRKind, &
        dknd=> idknd, aknd=> iaknd, &
        iiknd, ii4, iik1
 USE NaNs_etc, ONLY : &
        IsNaN
 IMPLICIT NONE

 PRIVATE

 PUBLIC :: Multigrid_2D_Type2

 PUBLIC :: Poissongrid2D
 PUBLIC :: MGAllocateData
 PUBLIC :: MGData_FirstIx_etc
 PUBLIC :: MGDeAllocateData

 PUBLIC :: MGDebugPlots

 INTEGER(KIND= KIND(0)), PUBLIC, PARAMETER :: MGUBoundOfFactor= 4

 TYPE Poissongrid2D
    TYPE (Poissongrid2D), POINTER :: Coarse
    REAL(KIND= dknd),     POINTER, DIMENSION(:,:) :: Phi, Charge, ChargeTerm, &
                                                     ooSumme
!! LOGICAL(KIND= ii4),   POINTER, DIMENSION(:,:)   :: Variable
    LOGICAL(KIND= iik1),   POINTER, DIMENSION(:,:)   :: Variable
    REAL(KIND= dknd),     POINTER, DIMENSION(:,:,:) :: Factor, EpsA, EpsAods
    REAL(KIND= dknd),     POINTER, DIMENSION(:)     :: dx, dy
    INTEGER(KIND= iiknd), POINTER, DIMENSION(:)     :: FirstIx, LastIx
    INTEGER(KIND= iiknd), POINTER, DIMENSION(:,:)   :: Number          !<-- LU
    REAL(KIND= dknd),     POINTER, DIMENSION(:)     :: RHS             !<-- LU
    INTEGER(KIND= iiknd) :: FirstIy, LastIy, nVariable
    INTEGER(KIND= iiknd) :: IstParallelerSatz, LUSet
    LOGICAL(KIND= ii4)   :: LUIsFactorised
 END TYPE

 LOGICAL(KIND= ii4), SAVE :: MGDebugPlots= .FALSE._ii4

 CONTAINS
!**********
!**********
 SUBROUTINE MGAllocateData( Who, D, il1, iu1, il2, iu2 )
 CHARACTER(LEN= *), INTENT(IN) :: Who
 TYPE (Poissongrid2D) :: D
 INTEGER(KIND= iiknd), INTENT(IN) :: il1, iu1, il2, iu2

 CHARACTER(LEN= *), PARAMETER :: MyName= 'MGAllocateData'
 INTEGER(KIND= KIND(0)) :: iStat

    ALLOCATE( D% Phi        (il1:iu1,il2:iu2), &
              D% ooSumme    (il1:iu1,il2:iu2), &
              D% Variable   (il1:iu1,il2:iu2), &
              D% Factor (1:MGUBoundOfFactor, &
                             il1:iu1,il2:iu2), &
              D% EpsA   (1:2,il1:iu1,il2:iu2), &
              D% EpsAods(1:2,il1:iu1,il2:iu2), &
              D% ChargeTerm(il1:iu1,il2:iu2), &
              D% Charge     (il1:iu1,il2:iu2), &
              D% dx         (il1:iu1), &
              D% dy                 (il2:iu2), &
       STAT= iStat )
    IF (iStat /= 0) THEN
       write (*,*) Who, ': I Cannot allocate. il1, iu1, il2, iu2:', il1, iu1, il2, iu2
       CALL StopWithError( Who )
    END IF
    NULLIFY( D% Coarse )
    D% ChargeTerm= 0
    D% Charge= 0
    D% ooSumme= 0
    D% Phi= 0
    D% Variable= .FALSE._iik1
    D% Factor= 0
    D% EpsA= 0
    D% EpsAods= 0

 END SUBROUTINE MGAllocateData
!**********
 SUBROUTINE MGData_FirstIx_etc( Who, D )
 CHARACTER(LEN= *), INTENT(IN) :: Who
 TYPE (Poissongrid2D) :: D

 CHARACTER(LEN= *), PARAMETER :: MyName= 'MGData_FirstIx_etc'
 INTEGER(KIND= iiknd) :: ix, iy, iy1, iy2
 INTEGER(KIND= KIND(0)) :: iStat

!! write (*,*) MyName, ': Who:', Who

    D% FirstIy= HUGE(D% FirstIy) / 4
    D% LastIy= -(HUGE(D% LastIy) - 1) / 4
    DO iy= LBOUND(D% Variable, DIM= 2), UBOUND(D% Variable, DIM= 2)
       DO ix= LBOUND(D% Variable, DIM= 1), UBOUND(D% Variable, DIM= 1)
          IF (D% Variable(ix,iy)) THEN
             D% FirstIy= MIN(D% FirstIy, iy)
             D% LastIy= MAX(D% LastIy, iy)  ! MAX() is not needed.
             !
             ! We are only searching whether there is any Variable at this ix.
             !
             EXIT  ! Process next iy
          END IF
       END DO
    END DO

    iy1= LBOUND(D% Variable, DIM= 2)
    iy2= UBOUND(D% Variable, DIM= 2)

    ALLOCATE(D% FirstIx (iy1:iy2), &
             D% LastIx(iy1:iy2), &
             STAT= iStat )
    IF (iStat /= 0) THEN
       write (*,*) MyName, ': I cannot allocate FirstIx etc.', iy1, iy2
       CALL StopWithError( MyName )
    END IF
!! write (*,*) MyName, ': associated(D% Variable):', associated(D% Variable)
!! write (*,*) MyName, ': iy1, iy2:', iy1, iy2
    D% FirstIx= HUGE(1_iiknd)
    D% LastIx= -(HUGE(1_iiknd) - 1)
!! write (*,*) MyName, ': D% FirstIy, D% LastIy:', D% FirstIy, D% LastIy
    DO iy= D% FirstIy, D% LastIy
!! write (*,*) MyName, ': iy:', iy
       DO ix= LBOUND(D% Variable, DIM= 1), UBOUND(D% Variable, DIM= 1)
          IF (D% Variable(ix,iy)) THEN
             D% FirstIx(iy)= ix
             EXIT
          END IF
       END DO
       DO ix= UBOUND(D% Variable, DIM= 1), LBOUND(D% Variable, DIM= 1), -1
          IF (D% Variable(ix,iy)) THEN
             D% LastIx(iy)= ix
             EXIT
          END IF
       END DO
    END DO

    ! 2010-08-13
    DO iy= LBOUND(D% FirstIx, DIM= 1), UBOUND(D% FirstIx, DIM= 1)
       IF (D% FirstIx(iy) > D% LastIx(iy)) THEN
        !!  D% FirstIx(iy)= -(HUGE(D% FirstIx(iy))-2) !! DEC-Nag mag das nicht.
          D% FirstIx (iy)= -(HUGE(1_iiknd) - 2)
          D% LastIx(iy)= D% FirstIx(iy) - 1
       END IF
    END DO

    ! 2015-04-26:
    D% LUIsFactorised= .FALSE._ii4

    NULLIFY(D% RHS) ! 2015-04-27
    NULLIFY(D% Coarse) ! 2015-07-17_2

 END SUBROUTINE MGData_FirstIx_etc
!**********
 RECURSIVE SUBROUTINE MGDeAllocateData( D )
 TYPE (Poissongrid2D) :: D

    IF (ASSOCIATED(D% Coarse)) THEN
       CALL MGDeAllocateData( D% Coarse )
       DEALLOCATE( D% Coarse )
    END IF
    DEALLOCATE( D% Phi )
    DEALLOCATE( D% Charge )
    DEALLOCATE( D% ChargeTerm )
    DEALLOCATE( D% ooSumme )
    DEALLOCATE( D% Variable )
    DEALLOCATE( D% Factor )
    DEALLOCATE( D% EpsA )
    DEALLOCATE( D% EpsAods )
    DEALLOCATE( D% dx )
    DEALLOCATE( D% dy )
    DEALLOCATE( D% FirstIx )
    DEALLOCATE( D% LastIx )

 END SUBROUTINE MGDeAllocateData
!**********
 SUBROUTINE SOR2D( Phi, Variable, Factor, Omega, NIter )
 REAL(KIND= aknd), DIMENSION(0:,0:),   INTENT(INOUT) :: Phi
 LOGICAL(KIND= ii4), DIMENSION(0:,0:), INTENT(IN)    :: Variable
 REAL(KIND= aknd), DIMENSION(:,0:,0:), INTENT(IN)    :: Factor
 REAL(KIND= dknd),                     INTENT(IN)    :: Omega
 INTEGER(KIND= iiknd),                 INTENT(IN)    :: NIter

 CHARACTER(LEN= *), PARAMETER :: MyName= 'SOR2D'
 REAL(KIND= dknd) :: PhiAlt, PhiNeu
 INTEGER(KIND= iiknd) :: jx1, ix1, nx, ny, ix, iy, iter, nThreads
 INTEGER :: NrofThreads

 write (*,*) MyName,': Tach.'

    NThreads= NrofThreads()

    nx= UBOUND(Phi, DIM= 1)-1
    ny= UBOUND(Phi, DIM= 2)-1
    jx1= 1
    DO iter= 1, 2*NIter
!!       ix1= jx1
!$OMP PARALLEL DO NUM_THREADS(NThreads) DEFAULT(SHARED) &
!$OMP PRIVATE(iy,ix1,ix,PhiAlt,PhiNeu)
       DO iy= 1, ny
          IF (MOD(iy,2_iiknd) == 1_iiknd) THEN
             ix1= jx1
          ELSE
             ix1= 3-jx1
          END IF
          DO ix= ix1, nx, 2
             IF (Variable(ix,iy)) THEN
                PhiAlt= Phi(ix,iy)
                PhiNeu= Phi(ix+1,iy  ) * Factor(1,ix,iy)  &
                      + Phi(ix  ,iy+1) * Factor(2,ix,iy)  &
                      + Phi(ix-1,iy  ) * Factor(3,ix,iy)  &
                      + Phi(ix  ,iy-1) * Factor(4,ix,iy)
                Phi(ix,iy)= PhiAlt + Omega * (PhiNeu-PhiAlt)
             END IF
          END DO
!!          ix1= 3-ix1  ! 2=3-1; 1=3-2
       END DO
       jx1= 3-jx1  ! 2=3-1; 1=3-2
    END DO
 END SUBROUTINE SOR2D
!**********
 SUBROUTINE Gauss_Seidel_2D_00( Phi, ChargeTerm, Variable, Factor, NIter, Omega )
 REAL(KIND= aknd), DIMENSION(0:,0:),   INTENT(INOUT) :: Phi
 REAL(KIND= aknd), DIMENSION(0:,0:),   INTENT(IN)    :: ChargeTerm
 LOGICAL(KIND= ii4),          DIMENSION(0:,0:),   INTENT(IN)    :: Variable
 REAL(KIND= aknd), DIMENSION(:,0:,0:), INTENT(IN)    :: Factor
 INTEGER,                              INTENT(IN)    :: NIter
 REAL(KIND= aknd),                     INTENT(IN)    :: Omega

 CHARACTER(LEN= *), PARAMETER :: MyName= 'Gauss_Seidel_2D_00'
 INTEGER :: ix1, jx1, nx, ny, ix, iy, iter
 REAL(KIND= aknd) :: PhiAlt, PhiNeu

 write (*,*) MyName,': Tach.'
    nx= UBOUND(Phi, DIM= 1)-1
    ny= UBOUND(Phi, DIM= 2)-1
    jx1= 1
    DO iter= 1, 2*NIter
       ix1= jx1
       DO iy= 1, ny 
          DO ix= ix1, nx, 2
             IF (Variable(ix,iy)) THEN
                PhiAlt= Phi(ix,iy)
                PhiNeu= Phi(ix+1,iy  ) * Factor(1,ix,iy)  &
                       + Phi(ix  ,iy+1) * Factor(2,ix,iy)  &
                       + Phi(ix-1,iy  ) * Factor(3,ix,iy)  &
                       + Phi(ix  ,iy-1) * Factor(4,ix,iy)  &
                       - ChargeTerm(ix,iy)
                Phi(ix,iy)= PhiAlt + Omega * (PhiNeu-PhiAlt)
             END IF
          END DO
          ix1= 3-ix1  ! 2=3-1; 1=3-2
       END DO
       jx1= 3-jx1  ! 2=3-1; 1=3-2
    END DO

 END SUBROUTINE Gauss_Seidel_2D_00
!**********
 SUBROUTINE Gauss_Seidel_2D( Phi, ChargeTerm, Variable, Factor, NIter, Omega )
 REAL(KIND= dknd), DIMENSION(:,:),   INTENT(INOUT) :: Phi
 REAL(KIND= dknd), DIMENSION(:,:),   INTENT(IN)    :: ChargeTerm
 LOGICAL(KIND= ii4), DIMENSION(:,:),   INTENT(IN)    :: Variable
 REAL(KIND= dknd), DIMENSION(:,:,:), INTENT(IN)    :: Factor
 INTEGER(KIND= iiknd),               INTENT(IN)    :: NIter
 REAL(KIND= aknd),                   INTENT(IN)    :: Omega

 CHARACTER(LEN= *), PARAMETER :: MyName= 'Gauss_Seidel_2D'
 INTEGER(KIND= iiknd) :: ix1, jx1, ix, iy, iter, nThreads
 REAL(KIND= dknd) :: PhiAlt, PhiNeu
 INTEGER :: NrofThreads
  
!! write (*,*) MyName,': lbound(Phi,1):', lbound(Phi,1)  ! 1
 write (*,*) MyName,': Tach.'
    NThreads= NrofThreads()
    jx1= 1
    DO iter= 1, 2*NIter
!$OMP PARALLEL DO NUM_THREADS(NThreads) DEFAULT(SHARED) &
!$OMP PRIVATE(iy,ix1,ix,PhiAlt,PhiNeu)
       DO iy= LBOUND(Phi, DIM= 2)+1, UBOUND(Phi, DIM= 2)-1, 1
 IF (MOD(iy,2_iiknd) == 1_iiknd) THEN
    ix1= jx1
 ELSE
    ix1= 3-jx1
 END IF
          DO ix= LBOUND(Phi, DIM= 1)+ix1, UBOUND(Phi, DIM= 1)-1, 2
             IF (Variable(ix,iy)) THEN
                PhiAlt= Phi(ix,iy)
                PhiNeu= Phi(ix+1,iy  ) * Factor(1,ix,iy)  &
                      + Phi(ix  ,iy+1) * Factor(2,ix,iy)  &
                      + Phi(ix-1,iy  ) * Factor(3,ix,iy)  &
                      + Phi(ix  ,iy-1) * Factor(4,ix,iy)  &
                      - ChargeTerm(ix,iy)
                Phi(ix,iy)= PhiAlt + Omega * (PhiNeu-PhiAlt)
             END IF
          END DO
!!          ix1= 3-ix1  ! 2=3-1; 1=3-2
       END DO
       jx1= 3-jx1  ! 2=3-1; 1=3-2
    END DO

 END SUBROUTINE Gauss_Seidel_2D
!********** 
 RECURSIVE SUBROUTINE Multigrid_2D_Type2 &
    ( Who, Data, NIter, Omega, Omega2, iPlane, Iteration, ChargeSquare, &
      nThreads )

 CHARACTER(LEN= *),    INTENT(IN)  :: Who
 TYPE (Poissongrid2D)              :: Data
 INTEGER(KIND= iiknd), INTENT(IN)  :: NIter
 REAL(KIND= aknd),     INTENT(IN)  :: Omega, Omega2
 INTEGER(KIND= iiknd), INTENT(IN)  :: iPlane, Iteration, nThreads
 REAL(KIND= dknd),     INTENT(OUT) :: ChargeSquare

 CHARACTER(LEN= *), PARAMETER :: MyName= 'Multigrid_2D_Type2'

 REAL(KIND= dknd), POINTER, DIMENSION(:,:) :: Phi, Charge, ChargeTerm, &
                                              ooSumme, &
                                              Phi2, Charge2, ChargeTerm2
 REAL(KIND= dknd), POINTER, DIMENSION(:,:,:) :: Factor, Factor2
 LOGICAL(KIND= iik1), POINTER, DIMENSION(:,:)   :: Variable, Variable2
 REAL(KIND= dknd), POINTER, DIMENSION(:)     :: dx, dy, dx2, dy2
 REAL(KIND= dknd), POINTER, DIMENSION(:,:,:) :: EpsA, EpsA2, EpsAods, EpsAods2

 REAL(KIND= aknd) :: ThisOmega
 REAL(KIND= dknd) :: FluxIntegral, Fxp, Fxm, Fyp, Fym, &
                     Phi0, Error, SumOfCharges, &
                     EpsArea, Length, Summe, &
                     aFactor1, aFactor2, aFactor3, aFactor4, &
                     Phi_Coarse, Blub, b1, b2, b3, b4, xNorm, yNorm
 INTEGER(KIND= iiknd) :: nx, ny, ix, iy, jx, jy, nVariable, iGamma, jIter
 INTEGER(KIND= iiknd) :: nx2, ny2, ix2, iy2, ixA, ixE, iyA, iyE, &
            NPotentials, ix0, iy0, jx2, jy2, jj, jThreads
 LOGICAL(KIND= ii4) :: SomeFlag
 INTEGER(KIND= KIND(0)) :: i, iStat

 !
 ! Wenige Iterationen auf dem feinen Gitter.
 ! Bestimmung der ErrorCharge.
 !
 ! Bestimung der Loesung auf einem groeberen Gitter, indem die ErrorCharge 
 ! des feineren Gitters als rechte Seite verwendet wird.
 !
 ! Interpolation der Loesung des Coarseen Gitters auf das feine Gitter.
 ! Korrektur der feinen Loesung.
 ! u_fein := u_fein + u_Coarse
 !
 ! Wieder einige Iterationen auf dem feinen Gitter.
 !
 !!    NThreads= NrofThreads()

 jThreads= nThreads
!! if (Dieses_Fortransystem(1:LEN('Open64')) /= 'Open64') then !! 84 & 85
 if (1 == 1) then  ! 2013-01-31, ohne Aktivierung:  10x/100, 6.x GFlop (32CPUs)
                   ! 2013-01-31, mit Aktivierung:  09x/100, 6.x GFlop (32CPUs)
    DO jj= nThreads, 1, -1
       IF ((Data% LastIy - Data% FirstIy+1) > 3*jj) THEN ! 2*jj, 3*jj, 4*jj: 09x/100
          ! > 1*jj gibt schlechtere Performance als 2*jj etc.
          EXIT
       END IF
    END DO
    jThreads= MAX( 1_iiknd, jj )
 end if


    Phi=> Data% Phi
    ix0= LBOUND(Phi, DIM= 1)
    iy0= LBOUND(Phi, DIM= 2)    ! iy0= 0
    nx= UBOUND(Phi, DIM= 1)-1  ! 0:nx+1
    ny= UBOUND(Phi, DIM= 2)-1
    ! Bestimmung der Anzahl der Gitterebenen im Coarseen Gitter
    nx2= 1+(nx-ix0)/2
    ny2= 1+(ny-iy0)/2
!! if (iPlane == 0) then
!!    IF (.NOT. ASSOCIATED(Data% Coarse)) THEN
!!       write (*,*) MyName,': ix0, iy0, nx, ny:', ix0, iy0, nx, ny
!!    END IF
!! end if
    ChargeTerm=> Data% ChargeTerm
    Variable   => Data% Variable
    Factor     => Data% Factor

 fGamma : &
 DO iGamma= 1, 2
    IF (iGamma == 1) THEN  ! hierher
       !
       ! Am Ende eines Gamma Durchlaufs wird bereits GaussSeidel ausgefuehrt.
       !

 !!   write (*,*) MyName,': lbound(Phi,1):', lbound(Phi,1)
 !!  CALL Gauss_Seidel_2D(Phi, ChargeTerm, Variable, Factor, NIter, Omega)
       CALL GaussSeidel2Dd( lbound(Phi, 1), lbound(Phi, 2), &
                            ubound(Phi, 1), ubound(Phi, 2), &
                            Data% FirstIx, Data% LastIx, &
                            Data% FirstIy, Data% LastIy, &
                            Phi, ChargeTerm, Variable, Factor, &
                            nIter, Omega, &
                            Data% nVariable, jThreads )

        IF (iPlane /= 0) THEN
        IF ((nx2 < 5) &
       .OR. (ny2 < 5)) THEN
        !!   CALL GaussSeidel2D( SIZE(Phi, DIM= 1)-2, SIZE(Phi, DIM= 2)-2, &
        !!                       Phi, ChargeTerm, Variable, Factor, 2*NIter, Omega )
 if (MGDebugPlots) then
 CALL DebugWritePhi('Nach Gauss-Seidel', 0)
 CALL DebugWriteCharge('Charge', Charge, 0)
 END IF
           RETURN
        END IF
        END IF
     END IF

     dx=>   Data% dx
     dy=>   Data% dy
     EpsA=> Data% EpsA

     IF (.NOT. ASSOCIATED(Data% Coarse)) THEN
!! write (*,*) MyName,': Data% Coarse wird angelegt.', iPlane, nx, ny
        ALLOCATE(Data% Coarse)
        NULLIFY(Data% Coarse% Variable)
        NULLIFY(Data% Coarse% Coarse)
     END IF
     Variable2=> Data% Coarse% Variable
     IF (.NOT. ASSOCIATED(Variable2)) THEN
        ALLOCATE( Variable2  (0:nx2+1, 0:ny2+1)); Variable2= .FALSE._ii4
        ALLOCATE(   Charge2  (0:nx2+1, 0:ny2+1)); Charge2= 0
        ALLOCATE(   EpsA2(1:2,0:nx2+1, 0:ny2+1)); EpsA2= 0
        ALLOCATE(ChargeTerm2(0:nx2+1, 0:ny2+1)); ChargeTerm2= 0
        Data% Coarse% Variable   => Variable2
        Data% Coarse% EpsA       => EpsA2
        Data% Coarse% Charge     => Charge2
        Data% Coarse% ChargeTerm=> ChargeTerm2

        !
        ! Bestimmung der Koeffizienten des Coarseen Gitters
        ! durch Integration der Koeffizienten des feinen Gitters.
        !
        DO iy2= 1, ny2
           iyA= 2*(iy2-1)+1 + iy0  ! 3= 2*(2-1)+1-0 = 2*1+1 = 3
           iyE= MIN(ny, iyA+1_iiknd)
        DO ix2= 1, nx2
           ixA= 2*(ix2-1)+1 + ix0
           ixE= MIN(nx, ixA+1_iiknd)

           !
           ! Bestimmung, welche Punkte im Coarseen Gitter variable Punkte sind.
           ! Wenn einer der benachbarten Punkte des feinen Gitters
           ! festgelegtes Potential hat, hat auch der Coarsee Gitterpunkt
           ! festgelegtes Potential.
           NPotentials= 0
           DO jy= iyA, iyE
              DO jx= ixA, ixE
                 IF (.NOT. Variable(jx,jy)) THEN
                    NPotentials= NPotentials + 1
                 END IF
              END DO
           END DO
           Variable2(ix2,iy2)= NPotentials == 0

           !
           ! Summation der Suszeptibilitaeten (Leitfaehigkeiten).
           ! Der Coarsee Gitterpunkt liegt zwischen den feinen Gitterpunkten.
           !

           EpsArea= 0
           DO jy= iyA, iyE
              EpsArea= EpsArea + EpsA(1,ixE,jy)
           END DO
           EpsA2(1,ix2,iy2)= EpsArea

           EpsArea= 0
           DO jx= ixA, ixE
              EpsArea= EpsArea + EpsA(2,jx,iyE)
           END DO
           EpsA2(2,ix2,iy2)= EpsArea

        END DO
        END DO

     !!   ALLOCATE(  Factor2(4, 0:nx2+1, 0:ny2+1))
        ALLOCATE( Factor2 (1:MGUBoundOfFactor, &
                             0:nx2+1,0:ny2+1), & !! 5 fuer externes GaussSeidel2D
                  EpsAods2(2,0:nx2+1,0:ny2+1), &
                  ooSumme   (0:nx2+1,0:ny2+1), &
                  dx2       (0:nx2+1), &
                  dy2               (0:ny2+1), &
          STAT= iStat )
        IF (iStat /= 0) THEN
           WRITE (*,*) MyName, ': I cannot allocate coarse Data.'
           CALL StopWithError( MyName )
        END IF

        !
        ! Bestimmung der Coarseen Schrittweiten.
        !
        dx2(0    )= 2*dx(ix0+0)
        dx2(nx2+1)= 2*dx( nx+1)
        DO ix2= 1, nx2
           ix= 2*(ix2-1)+1 + ix0
           Length= 0.5*dx(ix)
           IF (ix < nx) THEN
              Length= Length + dx(ix+1)
              IF (ix+1 < nx) THEN
                 Length= Length + 0.5*dx(ix+1)
              END IF;
           END IF;
           dx2(ix2)= Length
        END DO

        dy2(0    )= 2*dy(iy0+0)
        dy2(ny2+1)= 2*dy( ny+1)
        DO iy2= 1, ny2
           iy= 2*(iy2-1)+1 + iy0
           Length= 0.5*dy(iy)
           IF (iy < ny) THEN
              Length= Length + dy(iy+1)
              IF (iy+1 < ny) THEN
                 Length= Length + 0.5*dy(iy+1)
              END IF
           END IF
           dy2(iy2)= Length
        END DO

        Data% Coarse% dx=> dx2
        Data% Coarse% dy=> dy2

        Factor2= 0
        Data% Coarse% Factor=> Factor2
        EpsAods2= 0
        Data% Coarse% EpsAods=> EpsAods2
        ooSumme= 0
        Data% Coarse% ooSumme=> ooSumme
        nVariable= 0
        DO iy= 1, ny2
        DO ix= 1, nx2
           IF (Variable2(ix,iy)) THEN
              aFactor1= EpsA2(1,ix  ,iy  ) / dx2(ix)
              aFactor2= EpsA2(2,ix  ,iy  ) / dy2(iy)
              aFactor3= EpsA2(1,ix-1,iy  ) / dx2(ix-1)
              aFactor4= EpsA2(2,ix  ,iy-1) / dy2(iy-1)
              Summe= aFactor1+aFactor2+aFactor3+aFactor4
              IF (Summe /= 0) THEN
 call IsNaN( 1 / ((Summe + 1) / 2), SomeFlag )
 if (SomeFlag) then
    write (*,*) MyName,': ix, iy, Summe:', ix, iy, Summe
    write (*,'(2A,4E14.7)') MyName,': aFactor(1:4):', aFactor1, aFactor2, aFactor3, aFactor4
    CALL StopWithError( MyName )
 END IF
                 EpsAods2(1,ix  ,iy  )= aFactor1
                 EpsAods2(2,ix  ,iy  )= aFactor2
                 EpsAods2(1,ix-1,iy  )= aFactor3
                 EpsAods2(2,ix  ,iy-1)= aFactor4
! write (*,*) MyName,': iPlane, Summe:', iPlane, Summe
                 Factor2(1,ix,iy)= aFactor1 / Summe
                 Factor2(2,ix,iy)= aFactor2 / Summe
                 Factor2(3,ix,iy)= aFactor3 / Summe
                 Factor2(4,ix,iy)= aFactor4 / Summe
                 ooSumme(ix,iy)= 1 / Summe
                 nVariable= nVariable + 1
              ELSE
                 ooSumme(ix,iy)= 0
                 Variable2(ix,iy)= .FALSE.
              END IF
           END IF
        END DO
        END DO
        Data% Coarse% nVariable= nVariable
        CALL MGData_FirstIx_etc( Who= MyName, D= Data% Coarse )

        ALLOCATE(Phi2(0:nx2+1, 0:ny2+1)); Phi2= 0
        Data% Coarse% Phi=> Phi2

!! write (*,*) MyName,': Ende der Initialisierung.'
     END IF

     EpsAods=> Data% EpsAods
     Charge=> Data% Charge

     Charge2     => Data% Coarse% Charge
     ooSumme     => Data% Coarse% ooSumme
     Variable2   => Data% Coarse% Variable
     ChargeTerm2=> Data% Coarse% ChargeTerm

     ChargeSquare= 0
 !$OMP PARALLEL DO NUM_THREADS(jThreads) DEFAULT(SHARED) &
 !$OMP    PRIVATE(iy2,iyA,iyE,ix2,ixA,ixE,SumOfCharges,jy,jx, &
 !$OMP            Phi0,Fxp,Fyp,Fxm,Fym,FluxIntegral,Error,SomeFlag ) &
 !$OMP    REDUCTION(+:ChargeSquare)
     DO iy2= 1, ny2
        iyA= iy0+2*(iy2-1)+1
        iyE= MIN(ny, iyA+1_iiknd)
     DO ix2= 1, nx2
        ixA= ix0+2*(ix2-1)+1
        ixE= MIN(nx, ixA+1_iiknd)
        !
        ! Summation der Chargeen.
        !
 Charge2(ix2,iy2)= 0
 ChargeTerm2(ix2,iy2)= 0
        IF (Variable2(ix2,iy2)) THEN
           SumOfCharges= 0
           DO jy= iyA, iyE
           DO jx= ixA, ixE
              IF (Variable(jx,jy)) THEN
                 Phi0= Phi(jx,jy)
                 Fxp= ( Phi(jx+1,jy) - Phi0 ) * EpsAods(1,jx  ,jy  )
                 Fyp= ( Phi(jx,jy+1) - Phi0 ) * EpsAods(2,jx,jy)
                 Fxm= ( Phi0 - Phi(jx-1,jy) ) * EpsAods(1,jx-1,jy)
                 Fym= ( Phi0 - Phi(jx,jy-1) ) * EpsAods(2,jx,jy-1)

                 FluxIntegral= Fxp-Fxm + Fyp-Fym

                 Error= FluxIntegral - Charge(jx,jy)
                 ChargeSquare= ChargeSquare + Error**2
 if (0==1) then
    CALL IsNaN( ChargeSquare, SomeFlag )
    IF (SomeFlag) THEN
 !$OMP CRITICAL
       write (*,*) MyName // ': NaN in MultiGrid.'
       write (*,*) MyName,': ChargeSquare:', ChargeSquare, Error
 write (*,*) MyName,': Fxp,Fxm,Fyp,Fym:', Fxp,Fxm,Fyp,Fym
 write (*,*) MyName,': EpsAods(:,,):', EpsAods(:,jx,jy)
 write (*,*) MyName,': EpsAods(:,,):', EpsAods(1,jx-1,jy), EpsAods(2,jx,jy-1)
 write (*,*) MyName,': iPlane:', iPlane
 write (*,*) MyName,': jx,jy:', jx,jy
 write (*,*) MyName,': 5 Phi:', &
           Phi0, Phi(jx+1,jy), Phi(jx-1,jy), Phi(jx,jy+1), Phi(jx,jy-1)
 write (*,*) MyName,': ChargeTerm:', &
                          ChargeTerm(jx,jy), &
                          ChargeTerm(jx-1,jy), &
                          ChargeTerm(jx+1,jy), &
                          ChargeTerm(jx,jy-1), &
                          ChargeTerm(jx,jy+1)
    CALL StopWithError( MyName )
 !$OMP END CRITICAL
    END IF
 end if

! if (iPlane .GT. 1) &
! write (*,*) ' Rot, Error:', Rot, Error

                 SumOfCharges= SumOfCharges + Error
              END IF
           END DO
           END DO
           Charge2(ix2,iy2)=     SumOfCharges
           ChargeTerm2(ix2,iy2)= SumOfCharges * ooSumme(ix2,iy2)
        END IF
     END DO
     END DO
     CALL IsNaN( ChargeSquare, SomeFlag )
     IF (SomeFlag) THEN
        write (*,*) MyName // ': NaN in MultiGrid.'
        CALL StopWithError( MyName )
     END IF

!! write (*,*) MyName,': iPlane, ChargeSquare:', iPlane, ChargeSquare, Data% Coarse% nVariable
     IF (ChargeSquare < 1d-70) THEN
        EXIT fGamma
     END IF

 if (MGDebugPlots) then
 CALL DebugWriteCharge('Charge2', Charge2, 2)
 end if

     nVariable= Data% Coarse% nVariable

     IF (nVariable == 0) THEN
     !! CALL Gauss_Seidel_2D(Phi, ChargeTerm, Variable, Factor, 2*NIter, Omega)
        CALL GaussSeidel2D( SIZE(Phi, DIM= 1)-2, SIZE(Phi, DIM= 2)-2, &
                            Phi, ChargeTerm, Variable, Factor, &
                            INT(2*NIter, KIND= iiknd), Omega )
 if (MGDebugPlots) then
 CALL DebugWritePhi('Nach Gauss-Seidel', 0)
 CALL DebugWriteCharge('Charge', Charge, 0)
 END IF
     ELSE
        Phi2=> Data% Coarse% Phi
        Phi2= 0
!! write (*,*) MyName, ': iPlane, Omega:', iPlane, Omega
        CALL Multigrid_2D_Type2( MyName, Data% Coarse, &
                                 NIter+0_iiknd, &   !! +1 vor 2011-12-16
                                 REAL((1+22*Omega)/23, KIND= aknd), &  ! Hierher
                                 REAL((1+Omega2)/2, KIND= aknd), &     ! Hierher
                                 iPlane+1_iiknd, Iteration, &
                                 ChargeSquare= Blub, &
                                 nThreads= jThreads )

        !
        ! Jetzt Korrektur des feinen Ergebnisses durch Subtraktion des 
        ! groben Ergebnisses vom feinen Ergebnis.
        !
 if (MGDebugPlots) then
 CALL DebugWritePhi('Before Correction', 1)
 end if

!! write (*,*) MyName, ': Before Correction.', jThreads
        Variable2=> Data% Coarse% Variable
     !$OMP PARALLEL DO NUM_THREADS(jThreads) DEFAULT(SHARED) &
     !$OMP    PRIVATE(iy2,iyA,iyE,ix2,ixA,ixE,Phi_Coarse, &
     !$OMP            iy,ix,b1,b2,b3,b4,xNorm,yNorm,jx2,jy2)
        DO iy2= 1, ny2
           iyA= iy0+2*(iy2-1)
           iyE= MIN(ny, iyA+1_iiknd)
           DO ix2= 1, nx2
              IF (Variable2(ix2,iy2)) THEN
                 ixA= ix0+2*(ix2-1)
                 ixE= MIN(nx, ixA+1_iiknd)
                 IF ((iyA /= iyE) .AND. (ixA /= ixE)) THEN
               !!  IF (0==1) THEN
                    DO iy= iyA, iyE   !! 31, 1.6e-4
                       IF (iy==iyA) THEN
                          jy2= iy2-1
                          yNorm= 1.5 / 2.0
                       ELSE
                          jy2= iy2
                          yNorm= 0.5 / 2.0
                       END IF
                    DO ix= ixA, ixE
                       IF (Variable(ix,iy)) THEN
                          IF (ix==ixA) THEN
                             jx2= ix2-1
                             xNorm= 1.5 / 2.0
                          ELSE
                             jx2= ix2
                             xNorm= 0.5 / 2.0
                          END IF
                          b1= Phi2(jx2  ,jy2  )
                          b2= Phi2(jx2+1,jy2  ) - b1
                          b3= Phi2(jx2  ,jy2+1) - b1
                          b4= Phi2(jx2+1,jy2+1) &
                             -Phi2(jx2  ,jy2+1) &
                             -Phi2(jx2+1,jy2  ) + b1
                          Phi(ix,iy)= Phi(ix,iy) &
                              - Omega2 * (b1 + b2*xNorm &
                                        + b3*ynorm + b4*xnorm*ynorm)
                       END IF
                    END DO
                    END DO
                 ELSE
                    Phi_Coarse= Phi2(ix2,iy2)   !! Gleich.
                    DO iy= iyA, iyE
                    DO ix= ixA, ixE
                       IF (Variable(ix,iy)) THEN
                          Phi(ix,iy)= Phi(ix,iy) - Omega2 * Phi_Coarse
                       END IF
                    END DO
                    END DO
                 END IF
              END IF
           END DO
        END DO
!! write (*,*) MyName, ': After Correction.'
 if (MGDebugPlots) then
 CALL DebugWritePhi( 'After Correction', 2 )
 end if

        ! Wieder einige Iterationen auf dem feinen Gitter.
        ! Again some Iterations on the fine Grid.

     !! CALL Gauss_Seidel_2D(Phi, ChargeTerm, Variable, Factor, NIter, Omega)
 jIter= nIter      !! 8 3.50E-38
!! jIter= nIter + 1  !! 7 3.61E-38
 ThisOmega= Omega  !! 7 3.61E-38
!! ThisOmega= (1+Omega)/2  !! 8 4.44E-39
        IF (iGamma == 1) THEN
            jIter= 2 + jIter  !! 8 3.55E-38
        END IF
        CALL GaussSeidel2Dd( lbound(Phi, 1), lbound(Phi, 2), &
                             ubound(Phi, 1), ubound(Phi, 2), &
                             Data% FirstIx, Data% LastIx, &
                             Data% FirstIy, Data% LastIy, &
                             Phi, ChargeTerm, Variable, Factor, &
                             jIter, ThisOmega, &
                             Data% nVariable, jThreads )
     END IF
 END DO fGamma

 CONTAINS
    SUBROUTINE DebugWritePhi( Titel, i1 )
    CHARACTER(LEN= *), INTENT(IN) :: Titel
    INTEGER,           INTENT(IN) :: i1

    CHARACTER(LEN= 400) :: FileName
    INTEGER :: i11

    i11= 23
       WRITE (FileName,'(3(A,I3.3),A)') 'Potential-Plane=', iPlane, &
                                       '-Iteration=', Iteration, &
                                       '-i1=', i1, '.mtv'
       OPEN (UNIT= i11, FILE= FileName)
       WRITE (i11,*) ' $DATA= CONTOUR'
       WRITE (i11,'(3A)') ' % toplabel= "', Titel, '"'
       WRITE (i11,'(3A)') ' % subtitle= "', TRIM(FileName), '"'
       WRITE (i11,*) ' % nx=', SIZE(Phi, DIM= 1)-2
       WRITE (i11,*) ' % ny=', SIZE(Phi, DIM= 2)-2
       DO iy= LBOUND(Phi, DIM= 2)+1, UBOUND(Phi, DIM= 2)-1, 1
       DO ix= LBOUND(Phi, DIM= 1)+1, UBOUND(Phi, DIM= 1)-1, 1
          WRITE (i11,*) Phi(ix,iy)
       END DO
       END DO
       CLOSE (UNIT= i11, STATUS= 'KEEP')

    END SUBROUTINE DebugWritePhi

    !**********************

    SUBROUTINE DebugWriteCharge( Titel, Charge, i1 )
    CHARACTER(LEN= *),                  INTENT(IN) :: Titel
    REAL(KIND= dknd), DIMENSION(0:,0:), INTENT(IN) :: Charge
    INTEGER,                            INTENT(IN) :: i1

    CHARACTER(LEN= 400) :: FileName
    INTEGER :: i11

    i11= 23
       WRITE (FileName,'(3(A,I3.3),A)') 'Charge-Plane=', iPlane, &
                                       '-Iteration=', Iteration, &
                                       '-i1=', i1, '.mtv'
       OPEN (UNIT= i11, FILE= FileName)
       WRITE (i11,*) ' $DATA= CONTOUR'
       WRITE (i11,'(3A)') ' % toplabel= "', Titel, '"'
       WRITE (i11,'(3A)') ' % subtitle= "', TRIM(FileName), '"'
       WRITE (i11,*) ' % nx=', SIZE(Charge, DIM= 1)-2
       WRITE (i11,*) ' % ny=', SIZE(Charge, DIM= 2)-2
       DO iy= LBOUND(Charge, DIM= 2)+1, UBOUND(Charge, DIM= 2)-1, 1
       DO ix= LBOUND(Charge, DIM= 1)+1, UBOUND(Charge, DIM= 1)-1, 1
          WRITE (i11,*) Charge(ix,iy)
       END DO
       END DO
       CLOSE (UNIT= i11, STATUS= 'KEEP')

    END SUBROUTINE DebugWriteCharge

 END SUBROUTINE Multigrid_2D_Type2

!**********
 subroutine StopWithError( Bla )
 character(len= *) :: Bla
    write (*,*) Bla
    Stop
 end subroutine StopWithError
 END MODULE Multigrid
!##########

