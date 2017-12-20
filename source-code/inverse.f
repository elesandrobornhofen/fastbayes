      module inverse
      
      implicit none
      
      contains

************************************************************************
*************************** SUBROUTINE DKMWHF ********************
************************************************************************
	
      SUBROUTINE DKMWHF(A,V,W,det,zero,iflag,nrang,
     &				n,Dim_imap,iopt)
C=======================================================================

C     PURPOSE :  ROUTINE TO OBTAIN A GENERALISED INVERSE OF A
C                HALFSTORED SYMMETRIC MATRIX

C     STRATEGY : USE GAUSSIAN ALGORITHM,SELECTING THE LARGEST
C                DIAGONAL ELEMENTS (ABSOLUTE VALUE) AS PIVOTS;

C                WHEN A PIVOT SMALLER THAN A SPECIFIED VALUE
C                "zero" IS ENCOUNTERED, THE MATRIX IS NOT OF FULL
C                RANK AND THE ROWS & COLUMNS PERTAINING TO THE
C                REMAINING (zero) DIAGONAL ELEMENTS ARE SET TO
C                zero (0.D0); I.E., THE GENERALISED INVERSE OBTAINED
C                IS THE INVERSE OFF THE FULL RANK SUBMATRIX WITH A
C                NO. OF ROWS AND COLUMNS EQUAL TO THE NO. OF
C                DEPENDENCIES 'zeroED' OUT.

C                THIS IS AN ADAPTATION OF C.R. HENDERSON'S ROUTINE
C                "DJNVHF". IT DIFFERS FROM "DKMVHF" IN SUCH, THAT AN
C                EXTRA WORK VECTOR IS USED IN ORDER TO MINIMISE PAGE
C                TURNS; IT REQUIRES ABOUT 90% OF THE CPU TIME USED
C                BY "DKMVHF". IN ADDITION, THE LOG detERMINANT OF
C                THE FULL-RANK (SUB)MATRIX IS CALCULATED.

C     PARAMETERS
C                A     : DOUBLE PRECISION VECTOR OF LENGTH N*(N+1)/2,
C                        CONTAINING THE UPPER TRIANGLE OF THE MATRIX TO
C                        BE INVERTED STORED ROWWISE ON ENTRY AND ITS
C                        INVERSE ON EXIT.
C                V     : DOUBLE PRECISION VECTOR OF LENGTH N,
C                        USED AS WORK SPACE
C                W     : DOUBLE PRECISION VECTOR OF LENGTH N,
C                        CONTAINING THE DIAGONAL ELEMENTS OF THE
C                        INVERSE ON EXIT.
C                det   : DOUBLE PRECISION VARIABLE
C                        RETURNING THE LOG detERMINANT OF THE MATRIX
C                        (THE SIGN OF NEGATIVE PIVOTS IS IGNORED, BUT A
C                        WARNING MESSAGE IS PRINTED)
C                zero  : DOUBLE PRECISION VARIABLE
C                        SPECIFYING THE SMALLEST PIVOT TO BE TREATED AS
C                        NON-zero;
C                        THE APPROPRIATE VALUE DEPENDS ON THE SIZE OF
C                        ELEMENTS OF A AND ALSO ON THE AMOUNT OF PRIOR
C                        ADJUSTMENT FOR ANY KNOWN DEPENDENCIES.
C                        IF DEPENDENCIES ARE EXPECTED, THE APPROPRIATE
C                        ROWS AND COLUMNS SHOULD BE zeroED OUT BEFORE
C                        CALLING THE ROUTINE, THEN A VALUE OF 1.D-10
C                        TO 1.D-12 IS SUFFICIENT,
C                        IF SEVERAL DEPENDENCIES ARE TO BE IDENTIFIED
C                        BY THE ALGORITHM (AND IF ELEMENTS OF A ARE
C                        LARGE), A "zero" AS LARGE AS 1.D-8 OR 1.D-7
C                        MAY BE REQUIRED. THE PROGRAM PRINTS OUT ANY
C                        PIVOTS SMALLER THAN 1.D-5 TO ASSIST IN THE
C                        CHOICE OF A SUITABLE VALUE.
C                iflag : INTEGER VECTOR OF LENGTH N,
C                        USED AS WORK SPACE, CONTAINS ORDER IN WHICH
C                        PIVOTS WERE SELECTED ON EXIT
C                nrang : INTEGER VARIABLE,
C                        RETURNING RANK OF THE MATRIX
C                n     : ORDER OF THE MATRIX TO BE INVERTED.
C                iopt  : INTEGER VARIABLE, TO detERMINE WHETHER WARNING
C                        MESSAGES, ETC. ARE PRINTED ; SET : 0=NO, 1=YES.

C     ROUTINES REQUIRED : NONE

C     COPYRIGHT : KARIN MEYER 1986
C------------------------------------------------------KM--5/86------------
* Matrix wird in Vektor A (Laenge Dim_imap) gespeichert: [1 2 3 4 5 6 ...]

      REAL*8,INTENT(OUT)	:: det
      REAL*8,INTENT(IN)		:: zero
      INTEGER,INTENT(IN)	:: Dim_imap,n,iopt
      REAL*8				:: WW,XX,DMAX,AMAX,BMAX,DIMAX,V(n),W(n)
      INTEGER,INTENT(OUT)	:: nrang
      REAL*8,INTENT(INOUT)	:: A(Dim_imap)
      INTEGER	:: NEG,N1,I,II,IMAX,IMAXM1,J,IMAXP1,IL,IJ,iflag(n)
	
      V=0.0d0
      W=0.0d0
      iflag=0
	
C     ------------------
C     MATRIX IS A SCALAR
C     ------------------
      IF(N.EQ.1) THEN
        XX=A(1)
        IF(DABS(XX).GT.zero) THEN
          A(1)=1.D0/XX
          nrang=1
        ELSE
          A(1)=0.D0
          nrang=0
        END IF
        W(1)=A(1)
        IF(XX.GT.zero)THEN
          det=DLOG(XX)
        ELSE IF(XX.LT.-zero)THEN
          det=DLOG(-XX)
          IF(iopt.EQ.1)PRINT *,'NEGATIVE DIAG. SIGN IGNORED',XX
        ELSE
          det=0.D0
        END IF
        RETURN
      END IF

C     ----------
C     INITIALIZE
C     ----------
      NEG=0
      det=0.D0
      N1=N+1
      DO 1 I=1,N
 1    iflag(I)=0

C     PICK OUT DIAGONAL ELEMENTS
      II=-N
      DO 101 I=1,N
        II=II+N1
        W(I)=A(II)
 101  II=II-I

C     --------------------------
C     GAUSSIAN ELIMINATION STEPS
C     --------------------------
      DO 100 II=1,N

C     FIND DIAG. ELEMENT WITH LARGEST ABSOLUTE VALUE (PIVOT)
        DMAX=0.D0
        AMAX=0.D0
        DO 2 I=1,N
          IF(iflag(I).NE.0)GO TO 2
          BMAX=DABS(W(I))
          IF(BMAX.GT.AMAX)THEN
            DMAX=W(I)
            AMAX=BMAX
            IMAX=I
          END IF
 2      CONTINUE

C     CHECK FOR SINGULARITY
        IF(AMAX.LE.zero) GO TO 11
        IF(iopt.EQ.1.AND.AMAX.LT.1.D-5) PRINT *,'SMALL PIVOT ',II,DMAX
C     SET FLAG
        iflag(IMAX)=II
C     ACCUMULATE LOG detERMINANT
        det=det+DLOG(AMAX)
        IF(DMAX.LT.0.D0)THEN
          NEG=NEG+1
          IF(iopt.EQ.1)PRINT *,'NEGATIVE PIVOT, IGNORE SIGN FOR LOG det'
     *                       ,II,IMAX,DMAX
        END IF
	
        DIMAX=1.D0/DMAX
        IMAXM1=IMAX-1
        IMAXP1=IMAX+1

C     PICK OUT ELEMENTS FOR ROW IMAX
        IL=IMAX-N
        DO 3 I=1,IMAXM1
          IL=IL+N1-I
 3      V(I)=A(IL)
        IL=IL+N1-IMAX
        DO 4 I=IMAXP1,N
          IL=IL+1
 4      V(I)=A(IL)

C     TRANSFORM MATRIX
        IJ=0
        DO 7 I=1,IMAXM1
          WW=V(I)
          IF(WW.EQ.0.D0)THEN
            IJ=IJ+N1-I
          ELSE
            XX=WW*DIMAX
            IJ=IJ+1
            W(I)=W(I)-WW*XX
            DO 71 J=I+1,IMAXM1
              IJ=IJ+1
71          A(IJ)=A(IJ)-XX*V(J)
C        ELEMENT A(I,IMAX)
            IJ=IJ+1
            A(IJ)=XX
            DO 72 J=IMAXP1,N
              IJ=IJ+1
72          A(IJ)=A(IJ)-XX*V(J)
          END IF
 7      CONTINUE

C     ROW IMAX
        IJ=IJ+1
        W(IMAX)=-DIMAX
        DO 73 J=IMAXP1,N
          IJ=IJ+1
 73     A(IJ)=V(J)*DIMAX

        DO 74 I=IMAXP1,N
          WW=V(I)
          IF(WW.EQ.0.D0)THEN
            IJ=IJ+N1-I
          ELSE
            XX=WW*DIMAX
            IJ=IJ+1
            W(I)=W(I)-WW*XX
          DO 75 J=I+1,N
           IJ=IJ+1
 75       A(IJ)=A(IJ)-XX*V(J)
          END IF
 74     CONTINUE

 100  CONTINUE

      WRITE(*,*)'ln(Determinante):',det

 
C     ------------------------------------
C     STORE DIAGONALS BACK & REVERSE SIGNS
C     ------------------------------------
      IJ=0
      DO 9 I=1,N
        IJ=IJ+1
        WW=-W(I)
        W(I)=WW
        A(IJ)=WW
        DO 90 J=I+1,N
          IJ=IJ+1
 90     A(IJ)=-A(IJ)
 9    CONTINUE
      nrang=N
 
 300  IF(iopt.EQ.1.AND.NEG.GT.0)PRINT *,
     *                  'NO. OF NEGATIVE PIVOTS =',NEG
      RETURN

C      ---------------------------------------------------
C      MATRIX NOT OF FULL RANK, RETURN GENERALISED INVERSE
C      ---------------------------------------------------
 11   nrang=II-1
      IJ=0

      DO 14 I=1,N
        IF(iflag(I).EQ.0) THEN
C        ... zero OUT ROW/COLUMN
          W(I)=0.D0
          DO 12 J=I,N
            IJ=IJ+1
 12       A(IJ)=0.D0
        ELSE
          IJ=IJ+1
          WW=-W(I)
          W(I)=WW
          A(IJ)=WW
          DO 13 J=I+1,N
            IJ=IJ+1
            IF(iflag(J).NE.0) THEN
              A(IJ)=-A(IJ)
            ELSE
              A(IJ)=0.D0
            END IF
 13       CONTINUE
        END IF
 14   CONTINUE
 
      IF(iopt.EQ.1)PRINT 15,N,nrang
 15   FORMAT(' GENERALISED INVERSE OF MATRIX WITH ORDER =',I5,
     *                                     '   AND RANK =',I5)

      GO TO 300

      END SUBROUTINE DKMWHF
      
      
      end module inverse
