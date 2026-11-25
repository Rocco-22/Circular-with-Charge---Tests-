
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE fill0(bitmap, n)
      USE Kinds
      IMPLICIT INTEGER(KIND= iiknd) (i-n)
      IMPLICIT REAL(KIND= idknd) (a-h,o-z)
      INTEGER(KIND= iiknd) BITMAP(*), N

!     fill the first n bits of bitmap with zeroes.

      INTEGER(KIND= iiknd), PARAMETER :: NBPW= 30
      loop= n/NBPW
      nblw= MOD(n,NBPW)
      DO i= 1, LOOP, 1
         bitmap(i) = 0
      END DO
      IF (NBLW /= 0) THEN
         bitmap(loop+1)= MOD(bitmap(loop+1), INT(2**(NBPW-NBLW), KIND= iiknd) )
      END IF
      END SUBROUTINE fill0
! *******
      SUBROUTINE mark1(bitmap, n)
      USE Kinds
      IMPLICIT INTEGER(KIND= iiknd) (i-n)
      IMPLICIT REAL(KIND= idknd) (a-h,o-z)
      INTEGER(KIND= iiknd), DIMENSION(*), INTENT(INOUT) :: bitmap
      INTEGER(KIND= iiknd),               INTENT(IN) :: n

!     put a one in the nth bit of bitmap.

      INTEGER(KIND= iiknd), PARAMETER :: NBPW= 30
      nword= (n-1)/NBPW
      nbit= MOD(n-1_iiknd,NBPW)
      i= 2**(nbpw-nbit-1)
      bitmap(nword+1) = bitmap(nword+1) &
       + i*(1-MOD(bitmap(nword+1)/I,2_iiknd))
      END SUBROUTINE mark1
! ********
      FUNCTION iget(bitmap, n)
      USE Kinds
      IMPLICIT INTEGER(KIND= iiknd) (i-n)
      IMPLICIT REAL(KIND= idknd) (a-h,o-z)
      INTEGER(KIND= iiknd), DIMENSION(*), INTENT(IN) :: bitmap
      INTEGER(KIND= iiknd),               INTENT(IN) :: n

!     iget == 0 if the nth bit of bitmap is zero, else iget is one.

      INTEGER(KIND= iiknd), PARAMETER :: NBPW= 30
      nword = (n-1)/NBPW
      nbit = MOD(N-1_iiknd,NBPW)
      iget = MOD(bitmap(nword+1)/2**(NBPW-nbit-1),2)
      END FUNCTION iget
! ******
      SUBROUTINE gcontr(z, nrz, nx, ny, cv, ncv, zmax, bitmap, draw)
      USE Kinds
      IMPLICIT INTEGER(KIND= iiknd) (i-n)
      IMPLICIT REAL(KIND= idknd) (a-h,o-z)
      INTEGER(KIND= iiknd),                 INTENT(IN) :: nrz, nx, ny
      REAL(KIND= idknd), DIMENSION(nrz,ny), INTENT(IN) :: z
      INTEGER(KIND= iiknd),                 INTENT(IN) :: ncv
      REAL(KIND= idknd), DIMENSION(ncv),    INTENT(IN) :: cv
      REAL(KIND= idknd),                    INTENT(IN) :: zmax
      INTEGER(KIND= iiknd), DIMENSION(*),   INTENT(INOUT) :: bitmap
      EXTERNAL Draw

      CHARACTER(LEN= *), PARAMETER :: myname= 'gcontr'
!/*
!     This subroutine draws a contour through equal values of an array.
!
!     *****     Formal arguments     ***********************************
!
!     Z is the array for which contours are to be drawn.  The elements
!     of z are assumed to lie upon the nodes of a topologiCALLy
!     rectangular coordinate system - e.g. cartesian, polar (except
!     the origin), etc.
!
!     NRZ is the number of rows declared for Z in the CALLing program.
!
!     NX is the limit for the first subscript of Z.
!
!     NY is the limit for the second subscript of Z.
!
!     CV are the values of the contours to be drawn.
!
!     NCV is the number of contour values in CV.
!
!     ZMAX is the maximum value of z for consideration.  A value of
!     Z(I,J) greater than ZMAX is a signal that that point and the
!     grid line segments radiating from that point to it's neighbors
!     are to be excluded from contouring.
!
!     BITMAP is a work area large enough to hold 2*NX*NY*NCV bits.  it
!     is accessed by low-level routines, which are described below.
!     Let J be the number of useful bits in each word of BITMAP,
!     as determined by the user machine and implementation of
!     the BITMAP manipulation subprograms described below.  Then
!     the number of words required for the BITMAP is the floor of
!         (2*NX*NY*NCV+J-1)/J.
!
!     DRAW is a user-provided subroutine used to draw contours.
!     The CALLing sequence for DRAW is:
!
!         CALL DRAW (x,y,iflag)
!         let nx = int part of x, fx = fractional part of x.
!         Then x should be interpreted such that increases in nx
!         correspond to increases in the first subscript of Z, and
!         fx is the fractional distance from the abscissa corresponding
!         to nx to the abscissa corresponding to nx+1,
!         and y should be interpreted similarly for the second
!         subscript of Z.
!         The low-order digit of IFLAG will have one of the values:
!             1 - Continue a contour,
!             2 - Start a contour at a boundary,
!             3 - Start a contour not at a boundary,
!             4 - Finish a contour at a boundary,
!             5 - Finish a closed contour (not at a boundary).
!                 Note that requests 1, 4 and 5 are for pen-down
!                 moves, and that requests 2 and 3 are for pen-up
!                 moves.
!             6 - Set x and y to the approximate 'pen' position, using
!                 the notation discussed above.  This CALL may be
!                 ignored, the result being that the 'pen' position
!                 is taken to correspond to z(1,1).
!         iflag/10 is the contour number.
!
!     *****     EXTERNAL SUBPROGRAMS     *******************************
!
!     DRAW is the user-supplied line drawing subprogram described above.
!     DRAW may be sensitive to the host computer and to the plot device.
!     FILL0 is used to fill a BITMAP with zeroes.  CALL FILL0 (BITMAP,N)
!     fills the first n bits of BITMAP with zeroes.
!     MARK1 is used to place a 1 in a specific bit of the BITMAP.
!     CALL MARK1 (BITMAP,N) puts a 1 in the nth bit of the BITMAP.
!     IGET is used to determine the setting of a particular bit in the
!     BITMAP.  I=IGET(BITMAP,N) sets i to zero if the nth bit of the
!     BITMAP is zero, and sets i to one if the nth bit is one.
!     FILL0, MARK1 and IGET are machine sensitive.

!     ******************************************************************
!*/
      INTEGER(KIND= iiknd) :: L1(4), L2(4), IJ(2)

!     L1 AND L2 CONTAIN LIMITS USED DURING THE SPIRAL SEARCH FOR THE
!     BEGINNING OF A CONTOUR.
!     IJ STORES SUBCRIPTS USED DURING THE SPIRAL SEARCH.

      INTEGER(KIND= iiknd) :: I1(2), I2(2), I3(6)
! This save statements are needed for the shit Sun-compiler on "watt"
      save l2,imax,jmax,imin,jmin, ij,i,j,xy,x,y
!/*
!     I1, I2 AND I3 ARE USED FOR SUBSCRIPT COMPUTATIONS DURING THE
!     EXAMINATION OF LINES FROM Z(I,J) TO IT'S NEIGHBORS.
!*/
      DIMENSION XINT(4)

!     XINT IS USED TO MARK INTERSECTIONS OF THE CONTOUR UNDER
!     CONSIDERATION WITH THE EDGES OF THE CELL BEING EXAMINED.

      DIMENSION XY(2)

!     XY IS USED TO COMPUTE COORDINATES FOR THE DRAW SUBROUTINE.

      EQUIVALENCE (L2(1),IMAX), (L2(2),JMAX), (L2(3),IMIN), &
        (L2(4),JMIN)
      EQUIVALENCE (IJ(1),I), (IJ(2),J)
      EQUIVALENCE (XY(1),X), (XY(2),Y)

      DATA L1(3) /-1/, L1(4) /-1/
      DATA I1 /1,0/, I2 /1,-1/, I3 /1,0,0,1,1,0/

      L1(1) = NX
      L1(2) = NY
      DMAX = ZMAX

!     SET THE CURRENT PEN POSITION.  THE DEFAULT POSITION CORRESPONDS
!     TO Z(1,1).

      X = 1.0
      Y = 1.0
      CALL DRAW(X, Y, 6_iiknd)
      ICUR = MAX(1_iiknd,MIN(INT(X, KIND= iiknd),NX))
      JCUR = MAX(1_iiknd,MIN(INT(Y, KIND= iiknd),NY))

!     CLEAR THE BITMAP

      CALL FILL0(BITMAP, INT(2*NX*NY*NCV, KIND= iiknd))

!     SEARCH ALONG A RECTANGULAR SPIRAL PATH FOR A LINE SEGMENT HAVING
!     THE FOLLOWING PROPERTIES:
!          1.  THE END POINTS ARE NOT EXCLUDED,
!          2.  NO MARK HAS BEEN RECORDED FOR THE SEGMENT,
!          3.  THE VALUES OF Z AT THE ENDS OF THE SEGMENT ARE SUCH THAT
!              ONE Z IS LESS THAN THE CURRENT CONTOUR VALUE, AND THE
!              OTHER IS GREATER THAN OR EQUAL TO THE CURRENT CONTOUR
!              VALUE.

!     SEARCH ALL BOUNDARIES FIRST, THEN SEARCH INTERIOR LINE SEGMENTS.
!     NOTE THAT THE INTERIOR LINE SEGMENTS NEAR EXCLUDED POINTS MAY BE
!     BOUNDARIES.

      IBKEY = 0
   10 I = ICUR
      J = JCUR
   20 IMAX = I
      IMIN = -I
      JMAX = J
      JMIN = -J
      IDIR = 0
!     DIRECTION ZERO IS +I, 1 IS +J, 2 IS -I, 3 IS -J.
   30 NXIDIR = IDIR + 1
      K = NXIDIR
      IF (NXIDIR>3) NXIDIR = 0
   40 I = ABS(I)
      J = ABS(J)
      IF (Z(I,J)>DMAX) GO TO 140
      L = 1
!     L=1 MEANS HORIZONTAL LINE, L=2 MEANS VERTICAL LINE.
   50 IF (IJ(L)>=L1(L)) GO TO 130
      II = I + I1(L)
      JJ = J + I1(3-L)
      IF (Z(II,JJ)>DMAX) GO TO 130
!c      ASSIGN 100 TO JUMP
      Jump= 100
!     THE NEXT 15 STATEMENTS (OR SO) DETECT BOUNDARIES.
   60 IX = 1
      IF (IJ(3-L)==1) GO TO 80
      II = I - I1(3-L)
      JJ = J - I1(L)
      IF (Z(II,JJ)>DMAX) GO TO 70
      II = I + I2(L)
      JJ = J + I2(3-L)
      IF (Z(II,JJ)<DMAX) IX = 0
   70 IF (IJ(3-L)>=L1(3-L)) GO TO 90
   80 II = I + I1(3-L)
      JJ = J + I1(L)
      IF (Z(II,JJ)>DMAX) GO TO 90
!c      IF (Z(I+1,J+1)<DMAX) GO TO JUMP, (100, 280)
      IF (Z(I+1,J+1)<DMAX) then
         if (jump == 100) then
            goto 100
         else if (jump == 280) then
            goto 280
         else
            write (*,*) myname, ': bad jump:', jump
            STOP
         end if
      end if
   90 IX = IX + 2
!c      GO TO JUMP, (100, 280)
         if (jump == 100) then
            goto 100
         else if (jump == 280) then
            goto 280
         else
            write (*,*) myname, ': bad jump:', jump
            STOP
         end if

  100 IF (IX==3) GO TO 130
      IF (IX+IBKEY==0) GO TO 130
!     NOW DETERMINE WHETHER THE LINE SEGMENT IS CROSSED BY THE CONTOUR.
      II = I + I1(L)
      JJ = J + I1(3-L)
      Z1 = Z(I,J)
      Z2 = Z(II,JJ)
      DO 120 ICV=1,NCV
        IF (IGET(BITMAP,INT(2*(NX*(NY*(ICV-1)+J-1)+I-1)+L, KIND= iiknd))/=0) GO TO 120
        IF (cv(icv) <= MIN(z1,z2)) GOTO 110
        IF (cv(icv) <= MAX(z1,z2)) GOTO 190
         CALL MARK1(BITMAP, INT(2*(NX*(NY*(ICV-1)+J-1)+I-1)+L, KIND= iiknd))
  110 CONTINUE
  120 CONTINUE
  130 L = L + 1
      IF (L<=2) GO TO 50
  140 L = MOD(IDIR,2_iiknd) + 1
      IJ(L) = ISIGN( int(IJ(L)), int(L1(K)))

!     LINES FROM Z(I,J) TO Z(I+1,J) AND Z(I,J+1) ARE NOT SATISFACTORY.
!     CONTINUE THE SPIRAL.

  150 IF (IJ(L)>=L1(K)) GO TO 170
      IJ(L) = IJ(L) + 1
      IF (IJ(L)>L2(K)) GO TO 160
      GO TO 40
  160 L2(K) = IJ(L)
      IDIR = NXIDIR
      GO TO 30
  170 IF (IDIR==NXIDIR) GO TO 180
      NXIDIR = NXIDIR + 1
      IJ(L) = L1(K)
      K = NXIDIR
      L = 3 - L
      IJ(L) = L2(K)
      IF (NXIDIR>3) NXIDIR = 0
      GO TO 150
  180 IF (IBKEY/=0) RETURN
      IBKEY = 1
      GO TO 10

!     AN ACCEPTABLE LINE SEGMENT HAS BEEN FOUND.
!     FOLLOW THE CONTOUR UNTIL IT EITHER HITS A BOUNDARY OR CLOSES.

  190 IEDGE = L
      CVAL = CV(ICV)
      IF (IX/=1) IEDGE = IEDGE + 2
      IFLAG = 2 + IBKEY
      XINT(IEDGE) = (CVAL-Z1)/(Z2-Z1)
  200 XY(L) = REAL(IJ(L)) + XINT(IEDGE)
      XY(3-L) = REAL(IJ(3-L))
      CALL MARK1(BITMAP, INT(2*(NX*(NY*(ICV-1)+J-1)+I-1)+L, KIND= iiknd))
      CALL DRAW(X, Y, IFLAG+10_iiknd*ICV)
      IF (IFLAG<4) GO TO 210
      ICUR = I
      JCUR = J
      GO TO 20

!     CONTINUE A CONTOUR.  THE EDGES ARE NUMBERED CLOCKWISE WITH
!     THE BOTTOM EDGE BEING EDGE NUMBER ONE.

  210 NI = 1
      IF (IEDGE<3) GO TO 220
      I = I - I3(IEDGE)
      J = J - I3(IEDGE+2)
  220 DO 250 K=1,4
        IF (K==IEDGE) GO TO 250
        II = I + I3(K)
        JJ = J + I3(K+1)
        Z1 = Z(II,JJ)
        II = I + I3(K+1)
        JJ = J + I3(K+2)
        Z2 = Z(II,JJ)
        IF (cval <= MIN(z1,z2)) GOTO 250
        IF (cval > MAX(z1,z2)) GOTO 250
        IF (K==1) GO TO 230
        IF (K/=4) GO TO 240
  230   ZZ = Z1
        Z1 = Z2
        Z2 = ZZ
  240   XINT(K) = (CVAL-Z1)/(Z2-Z1)
        NI = NI + 1
        KS = K
  250 CONTINUE
      IF (NI==2) GO TO 260

!     THE CONTOUR CROSSES ALL FOUR EDGES OF THE CELL BEING EXAMINED.
!     CHOOSE THE LINES TOP-TO-LEFT AND BOTTOM-TO-RIGHT IF THE
!     INTERPOLATION POINT ON THE TOP EDGE IS LESS THAN THE INTERPOLATION
!     POINT ON THE BOTTOM EDGE.  OTHERWISE, CHOOSE THE OTHER PAIR.  THIS
!     METHOD PRODUCES THE SAME RESULTS IF THE AXES ARE REVERSED.  THE
!     CONTOUR MAY CLOSE AT ANY EDGE, BUT MUST NOT CROSS ITSELF INSIDE
!     ANY CELL.

      KS = 5 - IEDGE
      IF (XINT(3)<XINT(1)) GO TO 260
      KS = 3 - IEDGE
      IF (KS<=0) KS = KS + 4

!     DETERMINE WHETHER THE CONTOUR WILL CLOSE OR RUN INTO A BOUNDARY
!     AT EDGE KS OF THE CURRENT CELL.

  260 L = KS
      IFLAG = 1
!c      ASSIGN 280 TO JUMP
      Jump= 280
      IF (KS<3) GO TO 270
      I = I + I3(KS)
      J = J + I3(KS+2)
      L = KS - 2
  270 IF (IGET(BITMAP, INT(2*(NX*(NY*(ICV-1)+J-1)+I-1)+L, KIND= iiknd))==0) GO TO 60
      IFLAG = 5
      GO TO 290
  280 IF (IX/=0) IFLAG = 4
  290 IEDGE = KS + 2
      IF (IEDGE>4) IEDGE = IEDGE - 4
      XINT(IEDGE) = XINT(KS)
      GO TO 200

      END SUBROUTINE gcontr

