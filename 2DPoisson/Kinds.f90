!
! Ich brauche einige verschiedene INTEGER Kinds, zum Beispiel fuer
! die generischen Unterprogramme 'Register', 'Send Data' und so weiter.
! 'iiknd' sollte etwa 32 bit sein,
!
! iik1, iik2, iik3
!
!
! iGigL : Muss Maximalzahl von Gitterzellen auf einem Rechner fassen.
! iGigG : Muss Maximalzahl von Gitterzellen Global fassen.
!
!23456789012345678901234567890123456789012345678901234567890123456789012
 MODULE Kinds
 USE OMP_LIB

 IMPLICIT NONE
 PUBLIC

!! INCLUDE 'omp_lib.h' ! Fuer OMP_LOCK_KIND

 INTEGER, PARAMETER :: ir4  = SELECTED_REAL_KIND( 6)
 INTEGER, PARAMETER :: ir8  = SELECTED_REAL_KIND(12)
 INTEGER, PARAMETER :: ii4  = SELECTED_INT_KIND(9)
 INTEGER, PARAMETER :: ii8  = SELECTED_INT_KIND(16)
 INTEGER, PARAMETER :: iik1 = SELECTED_INT_KIND(2)
 INTEGER, PARAMETER :: iik2 = SELECTED_INT_KIND(4)
 INTEGER, PARAMETER :: iik3 = SELECTED_INT_KIND(8)
 INTEGER, PARAMETER :: iik4 = SELECTED_INT_KIND(16)
 INTEGER, PARAMETER :: iLen8= SELECTED_INT_KIND(16)

 INTEGER, PARAMETER :: iCellOffset= iik3 ! 2^31 Zellen pro Ebene.
!! INTEGER, PARAMETER :: iwCellOffset= iik3 ! 2^31 Zellen pro Ebene.
!! INTEGER, PARAMETER :: iwCellOffset= iik2 ! 2^15 Zellen pro Ebene.
 INTEGER, PARAMETER :: iw2CellOffset= iik2 ! 2^15 Zellen pro Stange.
 INTEGER, PARAMETER :: i2CellOffset= iik2 ! 2^15 Zellen pro Stange.

!
! iiknd
! iGig2    =>  ixyCol etc
! iGig     =>  NrofCell
! "iGig", as soon as we have memories larger than 8* 2**31
! INTEGER(KIND=iGig) has to be capable of storing numbers
!       larger than 2**31
!
 INTEGER, PARAMETER :: iiDefaultkind  = KIND(0)
!!      INTEGER, PARAMETER :: iiknd  = KIND(0)
 INTEGER, PARAMETER :: iiknd  = ii4
!!      INTEGER, PARAMETER :: iiknd  = iik2 ! EXTCHK


!
!     Falls 64 bit integers existieren,
!     ansonsten die kleinstmoeglichen Integer
!
 INTEGER, PARAMETER :: iPseudoI8 = MAX(ii8, iik1)

!! INTEGER, PARAMETER :: iLen8 = KIND(0) ! ubound(, KIND= iLen8)
!! INTEGER, PARAMETER :: iLen8 = MAX(ii8, KIND(0)) ! ubound(, KIND= iLen8)

!! INTEGER, PARAMETER :: iGig   = ii4  ! 32 Bit
 INTEGER, PARAMETER :: iGig   = ii8  ! 64 Bit
    INTEGER, PARAMETER :: iGigL  = iGig ! Lokale NrofCell
    INTEGER, PARAMETER :: iGigG  = iGig ! Globale NrofCell
 INTEGER, PARAMETER :: iGig2  = MAX(ii4, iiknd) ! ixyCol etc
!! INTEGER, PARAMETER :: iGigAda= MAX(ii4, iiknd) ! Ada NrofCell

 INTEGER, PARAMETER :: ikMEdge  = ii4  ! iMEdge, Index into other Fields.

 INTEGER, PARAMETER :: DefInt= iGig ! Ada-DefInt
!! INTEGER, PARAMETER :: DefInt= SELECTED_INT_KIND(16)  ! Ada-DefInt

 INTEGER, PARAMETER :: iFMTK  = ii8
 INTEGER, PARAMETER :: i64  = MAX(ii8, ii4)
 INTEGER, PARAMETER :: iAssignK  = ii8
 INTEGER, PARAMETER :: il1knd = 1                     !! logical*1

 !! INTEGER, PARAMETER :: WBINT= ii8 ! SuperLU

!!      INTEGER, PARAMETER :: idknd  = KIND(0E0)   !! for CRAY ??
!!      INTEGER, PARAMETER :: idknd  = ir4    !! Default-Kind
 INTEGER, PARAMETER :: idknd  = ir8    !! Default-Kind
 INTEGER, PARAMETER :: idknd_bits= 64

 INTEGER, PARAMETER :: iaknd= ir4  !! for large arrays
 INTEGER, PARAMETER :: iaknd_bits= 32

 INTEGER, PARAMETER :: izknd= iaknd  !! Zellart
!! INTEGER, PARAMETER :: izknd= idknd  !! Zellart
!! INTEGER, PARAMETER :: izknd= ir4

!!      INTEGER, PARAMETER :: iaknd= ir8  !! for large arrays
!!      INTEGER, PARAMETER :: iaknd_bits= 64

!!      INTEGER, PARAMETER :: ibvknd = ir4    !! For the basis-vectors
 INTEGER, PARAMETER :: ibvknd = iaknd    !! For the basis-vectors

 INTEGER, PARAMETER :: isapknd= ir8    !! For v,w of SAP
!!      INTEGER, PARAMETER :: isapknd= idknd  !! For v,w of SAP

 INTEGER, PARAMETER :: SparseRKind = KIND(0D0)
 INTEGER, PARAMETER :: ExternalRKind= ir4

! Mit idknd ( 10000 * dt ) : 175.21 MFlops
! Mit iaknd                : 208.46 MFlops
!! INTEGER, PARAMETER :: iPMLknd= iaknd  ! 6/100, 5.25 GFlop
 INTEGER, PARAMETER :: iPMLknd= idknd  ! ebenfalls.
 INTEGER, PARAMETER :: iIMPknd= idknd  ! 2025-05-18
 INTEGER, SAVE :: dummyintegerofKinds

 TYPE bvA
    REAL(KIND= ibvknd), POINTER, DIMENSION(:) :: A
 END TYPE bvA

 REAL(KIND= idknd), PARAMETER :: &
           UNDEFINED = - 51063._idknd * (2._idknd ** 40), &
           INFINITY  =   51063._idknd * (2._idknd ** 40), &
           AUTO      = UNDEFINED

 END MODULE Kinds

