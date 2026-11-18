
 PROGRAM ReadPhiData
 USE Kinds, ONLY : idknd
 USE ModulePlots, ONLY : &
       WritePotentialAsMTVFile

 IMPLICIT NONE

 INTEGER :: nx0, nx, ny0, ny, ix, iy, i
 REAL(KIND= idknd), ALLOCATABLE, DIMENSION(:) :: xPosition, yPosition
 REAL(KIND= idknd), ALLOCATABLE, DIMENSION(:,:) :: Phi

    READ (*,*) nx0, nx, ny0, ny
    ALLOCATE( xPosition(nx0:nx), &
              yPosition(ny0:ny), &
              Phi(nx0:nx, ny0:ny) )
    READ (*,*)     ! Skip a Line.
    DO ix= nx0, nx
       READ (*,*) xPosition(ix)
    END DO
    DO iy= ny0, ny
       READ (*,*) yPosition(iy)
    END DO

    ! Skip 6 Lines.
    DO i= 1, 6
       READ (*,*)
    END DO
    DO iy= ny0, ny
       DO ix= nx0, nx
          READ (*,*) Phi(ix,iy)
       END DO
    END DO
    CALL WritePotentialAsMTVFile( 'Potential', &
                                  'Potential2.mtv', &
                                  Phi, xPosition, yPosition )

 END PROGRAM ReadPhiData

