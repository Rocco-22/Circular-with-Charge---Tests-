
 SUBROUTINE GaussSeidel2D &
               ( nx, ny, &
                 Phi, ChargeTerm, Variable, Factor, NIter, Omega )
 USE Kinds
 USE Multigrid, ONLY : MGUBoundOfFactor
 IMPLICIT NONE

 INTEGER(KIND= KIND(0)),            INTENT(IN)  :: nx, ny
 REAL(KIND= idknd), DIMENSION(0:nx+1,0:ny+1),   INTENT(INOUT) :: Phi
 REAL(KIND= idknd), DIMENSION(0:nx+1,0:ny+1),   INTENT(IN)    :: ChargeTerm
 LOGICAL(KIND= iik1), DIMENSION(0:nx+1,0:ny+1),   INTENT(IN)    :: Variable
 REAL(KIND= idknd), DIMENSION(1:MGUBoundOfFactor, &
                                0:nx+1,0:ny+1), INTENT(IN)    :: Factor
 INTEGER(KIND= iiknd),                          INTENT(IN)    :: NIter
 REAL(KIND= iaknd),                             INTENT(IN)    :: Omega

 CHARACTER(LEN= *), PARAMETER :: MyName= 'GaussSeidel2D'
 LOGICAL(KIND= iik1), DIMENSION(0-1:ny+1) :: PhiRotBestimmt, PhiSchwarzBestimmt
 INTEGER(KIND= iiknd) :: ix1, jx1, ix, iy, iter, jy, lix1
 LOGICAL(KIND= ii4), PARAMETER :: SingleSweep= .TRUE._ii4
 REAL(KIND= idknd) :: PhiOld, PhiNew
 INTEGER(KIND= iiknd) :: NThreads
 INTEGER :: NrofThreads

    NThreads= 1 * NrofThreads()

    DO iter= 1, 2*NIter
       PhiRotBestimmt= .FALSE.
       PhiRotBestimmt(0)= .TRUE.
       PhiSchwarzBestimmt= PhiRotBestimmt
 !$OMP PARALLEL DO NUM_THREADS(NThreads) &
 !$OMP   PRIVATE(iy, lix1, ix, PhiOld, PhiNew, jy) 
       DO iy= 1, ny, 1
          IF (MOD(iy-1, 2) == 0) THEN
             lix1= 1
          ELSE
             lix1= 2
          END IF
          DO ix= 0 + lix1, nx, 2
             IF (Variable(ix,iy)) THEN
                PhiOld= Phi(ix,iy)
                PhiNew= Phi(ix+1,iy  ) * Factor(1,ix,iy) &
                      + Phi(ix  ,iy+1) * Factor(2,ix,iy) &
                      + Phi(ix-1,iy  ) * Factor(3,ix,iy) &
                      + Phi(ix  ,iy-1) * Factor(4,ix,iy) &
                      - ChargeTerm(ix,iy)
                Phi(ix,iy)= PhiOld + Omega * (PhiNew-PhiOld);
             END IF
          END DO
          PhiRotBestimmt(iy)= .TRUE.
 !$OMP FLUSH(PhiRotBestimmt)
          IF (SingleSweep) THEN
             IF (iy-1-1 >= LBOUND(Phi, DIM= 2)) THEN
                jy= iy-1
                IF (PhiRotBestimmt(jy)) THEN
                   DO ix= 0 + lix1, nx, 2
                      IF (Variable(ix,jy)) THEN
                         PhiOld= Phi(ix,jy)
                         PhiNew= Phi(ix+1,jy  ) * Factor(1,ix,jy) &
                               + Phi(ix  ,jy+1) * Factor(2,ix,jy) &
                               + Phi(ix-1,jy  ) * Factor(3,ix,jy) &
                               + Phi(ix  ,jy-1) * Factor(4,ix,jy) &
                               - ChargeTerm(ix,jy)
                         Phi(ix,jy)= PhiOld + Omega * (PhiNew-PhiOld)
                      END IF
                   END DO
                   PhiSchwarzbestimmt(jy)= .TRUE.
                END IF
             END IF
          END IF
       END DO

       IF (SingleSweep) THEN
          !$OMP PARALLEL DO NUM_THREADS(NThreads) &
          !$OMP   PRIVATE(iy, lix1, ix, PhiOld, PhiNew, jy)
          DO jy= 1, ny, 1
             IF (.NOT. PhiSchwarzbestimmt(jy)) THEN
                IF (MOD(jy-2, 2) == 0) THEN
                   lix1= 1
                ELSE
                   lix1= 2
                END IF
                DO ix= 0 + lix1, nx, 2
                   IF (Variable(ix,jy)) THEN
                      PhiOld= Phi(ix,jy)
                      PhiNew= Phi(ix+1,jy  ) * Factor(1,ix,jy) &
                            + Phi(ix  ,jy+1) * Factor(2,ix,jy) &
                            + Phi(ix-1,jy  ) * Factor(3,ix,jy) &
                            + Phi(ix  ,jy-1) * Factor(4,ix,jy) &
                            - ChargeTerm(ix,jy)
                      Phi(ix,jy)= PhiOld + Omega * (PhiNew-PhiOld)
                   END IF
                END DO
             END IF
          END DO
       END IF

    END DO

 END SUBROUTINE GaussSeidel2D
!***************
 SUBROUTINE GaussSeidel2Dd &
               ( lbx, lby, iubx, iuby, &
                 FirstIx, LastIx, FirstIy, LastIy, &
                 Phi, ChargeTerm, Variable, Factor, NIter, Omega, nVariable, &
                 nThreads )
 USE Kinds
 USE Multigrid, ONLY : MGUBoundOfFactor

 IMPLICIT NONE

 INTEGER(KIND= KIND(0)), INTENT(IN) :: lbx, lby, iubx, iuby
 INTEGER(KIND= iiknd), INTENT(IN), DIMENSION(lby:iuby)  :: FirstIx, LastIx
 INTEGER(KIND= iiknd), INTENT(IN)                       :: FirstIy, LastIy
 REAL(KIND= idknd),    INTENT(INOUT), DIMENSION(lbx:iubx,lby:iuby) :: Phi
 REAL(KIND= idknd),    INTENT(IN), DIMENSION(lbx:iubx,lby:iuby)  :: ChargeTerm
 LOGICAL(KIND= iik1),  INTENT(IN), DIMENSION(lbx:iubx,lby:iuby) :: Variable
 REAL(KIND= idknd),    INTENT(IN), DIMENSION(1:MGUBoundOfFactor, &
                                              lbx:iubx,lby:iuby) :: Factor
 INTEGER(KIND= iiknd), INTENT(IN)    :: NIter, nThreads
 REAL(KIND= iaknd),    INTENT(IN)    :: Omega
 INTEGER(KIND= iiknd), INTENT(IN)    :: nVariable

 CHARACTER(LEN= *), PARAMETER :: MyName= 'GaussSeidel2Dd'
 INTEGER(KIND= iiknd) :: ix1, ix, iy, iter, lix1
 REAL(KIND= idknd) :: PhiOld, PhiNew
 INTEGER(KIND= iiknd) :: nNIter

    nNIter= 2*NIter

 !$OMP PARALLEL NUM_THREADS(nThreads) &
 !$OMP   DEFAULT (SHARED) &
 !$OMP   PRIVATE(iy, lix1, ix1, ix, PhiOld, PhiNew )
    DO iter= 1, nNIter

!! !$OMP BARRIER

 !$OMP DO
       DO iy= FirstIy, LastIy
          IF (MOD(iy-1, 2) == MOD(iter-1,2)) THEN
             lix1= lbx+1  !! 1
          ELSE
             lix1= lbx+2  !! 2
          END IF
          ix1= FirstIx(iy) - ABS( MOD(FirstIx(iy),2_iiknd) - MOD(lix1,2_iiknd) )
          DO WHILE (ix1 < FirstIx(iy))
             ix1= ix1 + 2
          END DO
          DO ix= ix1, LastIx(iy), 2
             IF (Variable(ix,iy)) THEN
                PhiOld= Phi(ix,iy)
                PhiNew= Phi(ix+1,iy  ) * Factor(1,ix,iy) &
                      + Phi(ix  ,iy+1) * Factor(2,ix,iy) &
                      + Phi(ix-1,iy  ) * Factor(3,ix,iy) &
                      + Phi(ix  ,iy-1) * Factor(4,ix,iy) &
                      - ChargeTerm(ix,iy)
                Phi(ix,iy)= PhiOld + Omega * (PhiNew-PhiOld);
             END IF
          END DO
       END DO

    END DO
 !$OMP END PARALLEL

 END SUBROUTINE GaussSeidel2Dd

