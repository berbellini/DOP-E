      SUBROUTINE FARB2D (X,LX,Y,LY,Z,NXDIM,CN,ICOL,NC,MODE)
C
      DIMENSION X(LX),Y(LY),Z(NXDIM,LY),CN(NC),ICOL(NC+1)
C
C     FILL AREA WITH BICUBICS FOR 2D CONTOUR PLOTTING
C     -----------------------------------------------
C     FARB-E-2D  VERSION 2.1, 10/1988
C
C     T R I P   ALGORITHM              A. PREUSSER
C
C     AUTHOR: A. PREUSSER
C             FRITZ-HABER-INSTITUT DER MPG
C             FARADAYWEG 4-6
C             D-1000 BERLIN 33
C
C     INPUT PARAMETERS
C     X       ARRAY OF LENGTH LX FOR X-COORDINATES OF
C             A REGULAR GRID
C             IN ASCENDING ORDER.
C
C                     X- AND Y-COORDINATES MUST BE GIVEN
C                             IN CENTIMETERS
C                             ==============
C
C     LX      NUMBER OF GRID LINES X= X(I), I=1,LX
C             PARALLEL TO Y-AXIS.
C     Y       ARRAY OF LENGTH LY FOR Y-COORDINATES
C             IN ASCENDING ORDER.
C     LY      NUMBER OF GRID LINES Y= Y(I), I=1,LY
C             PARALLEL TO X-AXIS.
C     Z       2-DIMENSIONAL ARRAY DIMENSIONED Z(NXDIM,...)
C             DEFINING THE Z-VALUES AT THE GRID POINTS.
C             THE POINT WITH THE COORDINATES X(K), Y(L)
C             RECEIVES THE VALUE Z(K,L), K=1,LX, L=1,LY.
C     NXDIM   FIRST DIMENSION OF ARRAY Z
C     CN      ARRAY OF LENGTH NC FOR THE Z-VALUES OF
C             THE CONTOURS (CONTOUR LEVELS)
C             IN ASCENDING ORDER
C     ICOL    INTEGER ARRAY OF LENGTH NC+1 FOR
C             THE COLOURS TO BE USED FOR THE LINES OR AREAS.
C             VALUES FROM THIS ARRAY ARE PASSED TO
C             THE USER SUPPLIED SUBROUTINE USRPLT.
C             ICOL(I) IS USED FOR THE AREA, WHERE
C                  Z .GT. CN(I-1)        AND
C                  Z .LE. CN(I),
C             FOR I=2,NC.
C             AREAS, WHERE Z.LE.CN(1)
C             ARE FILLED WITH COLOUR ICOL(1),
C             AND AREAS, WHERE Z.GT.ICOL(NC)
C             ARE FILLED WITH COLOUR ICOL(NC+1).
C     NC      NUMBER OF CONTOUR LEVELS, NC.LE.100
C     MODE          0, FILL AREA ONLY
C                   1, LINES ONLY
C                   2, FILL AREA AND LINES
C
C
C     OUTPUT
C     IS PERFORMED BY CALLS TO THE SUBROUTINE    USRPLT
C     TO BE SUPPLIED BY THE USER (AN EXAMPLE FOR USRPLT
C     IS INCLUDED.)
C
C     PARAMETERS OF USRPLT
C                SUBROUTINE USRPLT (X,Y,N,NCOL,MODE)
C                X,Y     REAL ARRAYS OF LENGTH N FOR
C                        THE COORDINATES OF THE POLYGON
C                        TO BE PLOTTED.
C                N       NUMBER OF POINTS OF POLYGON
C                NCOL    COLOUR TO BE USED
C                        FOR THE AREA OR THE LINE.
C                        FOR NCOL, THE PROGRAM PASSES
C                        VALUES OF ICOL AS DESCRIBED ABOVE.
C                MODE    1, LINE DRAWING
C                        0, FILL AREA
C     -------------------------------------------------------------
C
C     THIS MODULE (FARB2D) IS BASED ON SUBROUTINE SFCFIT OF
C          ACM ALGORITHM 474 BY H.AKIMA
C
      DIMENSION ZA(4,2),ZB(5),ZAB(2,3),ZX(2),ZY(2),ZXY(2)
     A, XX(4),YY(4),ZZ(4),ZZX(4),ZZY(4),ZZXY(4)
C
      DATA X4,Z43,X5,Z53,Z63,Z5B1,Z5B2,Z5B3,Z5B4,Z5B5
     A, A5,ZA5B2,ZA5B3,ZA5B4,X6,Z64,Z6B1,Z6B2,Z6B3,Z6B4,Z6B5
     B, ZX,ZY,ZXY,ZAB,ZA,ZB /46*0/
C
C
C    PRELIMINARY PROCESSING
C
      IFA= 0
      LX0= LX
      LXM1= LX0-1
      LXM2= LXM1 - 1
      LY0= LY
      LYM1= LY0 -1
      LYM2= LYM1-1
C
C     ERROR CHECK
C
      IF(LXM2.LT.0) GOTO 400
      IF(LYM2.LT.0) GOTO 410
*
c      DO 20 IX=2,LX0
      DO 20 IX=2,LX
C        IF (X(IX-1)-X(IX)) 20,460,470
        IF (X(IX-1) .lt. X(IX)) then
            go to 20
        ELSE IF (X(IX-1) .eq. X(IX)) then
            go to 460
        ELSE
            go to 470
        ENDIF
   20 CONTINUE
C
c      DO 70 IY=2,LY0
      DO 70 IY=2,LY
C        IF (Y(IY-1)-Y(IY)) 70,490,500
        IF (Y(IY-1) .lt. Y(IY)) then
            go to 70
        ELSE IF (Y(IY-1) .eq. Y(IY)) then
            go to 490
        ELSE
            go to 500
        ENDIF
   70 CONTINUE
C
      DO 80 I=2,NC
C        IF (CN(I-1)-CN(I)) 80,530,540
        IF (CN(I-1) .le. CN(I)) then
            go to 80
        ELSE IF (CN(I-1) .eq. CN(I)) then
            go to 530
        ELSE
            go to 540
        ENDIF
   80 CONTINUE
C
      IF (NC.LE.0) GOTO 560
C
C
C  MAIN DO-LOOPS
C
      DO 390 IY=2,LY0
        IYM2= IY-2
        IYM3= IY-3
        IYML= IY-LY0
        IYML1= IYML+1
        IX6= 0
        DO 380 IX=1,LX0
          IXM1= IX-1
          IXML= IX- LX0
C
C ROUTINES TO PICK UP NECESSARY X,Y, AND Z VALUES TO
C COMPUTE THE ZA,ZB, AND ZAB VALUES, AND TO ESTIMATE
C THEM WHEN NECESSARY
C PRELIMINARY WHEN IX.EQ.1
C
          IF(IXM1.NE.0) GOTO 150
          Y3= Y(IY-1)
          Y4= Y(IY)
          B3= 1./(Y4-Y3)
          IF (IYM2.GT.0) B2= 1./(Y3-Y(IY-2))
          IF (IYM3.GT.0) B1= 1./(Y(IY-2)-Y(IY-3))
          IF (IYML.LT.0) B4= 1./(Y(IY+1)-Y4)
          IF (IYML1.LT.0) B5= 1./(Y(IY+2)-Y(IY+1))
          GOTO 180
C
C  TO SAVE THE OLD VALUES
C
  150     ZA(1,1)= ZA(2,1)
          ZA(1,2)=   ZA(2,2)
          X3= X4
          Z33= Z43
          ZA(2,1)= ZA(3,1)
          ZA(2,2)= ZA(3,2)
          ZAB(1,1)= ZAB(2,1)
          ZAB(1,2)= ZAB(2,2)
          ZAB(1,3)= ZAB(2,3)
  160     X4= X5
          Z43= Z53
          ZB(1)= Z5B1
          ZB(2)= Z5B2
          ZB(3)= Z5B3
          ZB(4)= Z5B4
          ZB(5)= Z5B5
          ZA(3,1)= ZA(4,1)
          ZA(3,2)= ZA(4,2)
          ZAB(2,1)= ZA5B2
          ZAB(2,2)= ZA5B3
          ZAB(2,3)= ZA5B4
  170     X5= X6
          Z53= Z63
          Z54= Z64
          Z5B1= Z6B1
          Z5B2= Z6B2
          Z5B3= Z6B3
          Z5B4= Z6B4
          Z5B5= Z6B5
C TO COMPUTE THE ZA, ZB, AND ZAB VALUES AND
C TO ESTIMATE THE ZB VALUES
C WHEN (IY.LE.3).OR.(IY.GE.LY-1)
C
  180     IX6= IX6 + 1
          IF (IX6.GT.LX0) GOTO 260
          X6= X(IX6)
          Z63= Z(IX6,IY-1)
          Z64= Z(IX6,IY)
          Z6B3= (Z64-Z63)*B3
          IF (LYM2.EQ.0) GOTO 200
          IF (IYM2.EQ.0) GOTO 190
          Z62= Z(IX6,IY-2)
          Z6B2= (Z63-Z62)*B2
          IF (IYML.NE.0) GOTO 190
          Z6B4= Z6B3 + Z6B3 -Z6B2
          GOTO 210
  190     Z65= Z(IX6,IY+1)
          Z6B4= (Z65-Z64)*B4
          IF (IYM2.NE.0) GOTO 210
          Z6B2= Z6B3 + Z6B3 -Z6B4
          GOTO 210
  200     Z6B2= Z6B3
          Z6B4= Z6B3
  210     IF (IYM3.LE.0) GOTO 220
          Z6B1= (Z62-Z(IX6,IY-3))*B1
          GOTO 230
  220     Z6B1= Z6B2 + Z6B2 -Z6B3
  230     IF (IYML1.GE.0) GOTO 240
          Z6B5= (Z(IX6,IY+2) - Z65)*B5
          GOTO 250
  240     Z6B5= Z6B4 + Z6B4 -Z6B3
  250     IF (IX6.EQ.1) GOTO 170
          A5= 1./(X6-X5)
          ZA(4,1)= (Z63-Z53)*A5
          ZA(4,2)= (Z64-Z54)*A5
          ZA5B2= (Z6B2-Z5B2)*A5
          ZA5B3= (Z6B3-Z5B3)*A5
          ZA5B4= (Z6B4-Z5B4)*A5
          IF (IX6.EQ.2) GOTO 160
          GOTO 280
C     TO ESTIMATE THE ZA AND ZAB VALUES
C     WHEN (IX.GE.LX-1).AND.(LX.GT.2)
  260     IF (LXM2.EQ.0) GOTO 270
          ZA(4,1)= ZA(3,1) +ZA(3,1) - ZA(2,1)
          ZA(4,2)= ZA(3,2) +ZA(3,2) - ZA(2,2)
          IF (IXML.EQ.0) GOTO 290
          ZA5B2= ZAB(2,1) + ZAB(2,1) - ZAB(1,1)
          ZA5B3= ZAB(2,2) + ZAB(2,2) - ZAB(1,2)
          ZA5B4= ZAB(2,3) + ZAB(2,3) - ZAB(1,3)
          GOTO 290
C     TO ESTIMATE THE ZA AND ZAB VALUES
C     WHEN (IX.GE.LX-1).AND.(LX.EQ.2)
  270     ZA(4,1)= ZA(3,1)
          ZA(4,2)= ZA(3,2)
          IF (IXML.EQ.0) GOTO 290
          ZA5B2= ZAB(2,1)
          ZA5B3= ZAB(2,2)
          ZA5B4= ZAB(2,3)
C     TO ESTIMATE THE ZA AND ZAB VALUES WHEN IX EQ 1
  280     IF (IXM1.NE.0) GOTO 290
          ZA(2,1)= ZA(3,1) +ZA(3,1) - ZA(4,1)
          ZA(1,1)= ZA(2,1) +ZA(2,1) - ZA(3,1)
          ZA(2,2)= ZA(3,2) +ZA(3,2) - ZA(4,2)
          ZA(1,2)= ZA(2,2) +ZA(2,2) - ZA(3,2)
          ZAB(1,1)= ZAB(2,1) + ZAB(2,1) -ZA5B2
          ZAB(1,2)= ZAB(2,2) + ZAB(2,2) -ZA5B3
          ZAB(1,3)= ZAB(2,3) + ZAB(2,3) -ZA5B4
          GOTO 300
C     NUMERICAL DIFFERENTATION  ---- TO DETERMINE
C     PARTIAL DERIV. ZX,ZY, AND ZXY AS WEIGHTED MEANS OF
C     DIVIDED DIFFERENCES ZA, ZB, AND ZAB, RESPECTIVELY
C
C TO SAVE THE OLD VALUES WHEN IX.NE.1
  290     ZX33= ZX(1)
          ZX34= ZX(2)
          ZY33= ZY(1)
          ZY34= ZY(2)
          ZXY33= ZXY(1)
          ZXY34= ZXY(2)
C
C NEW COMPUTATION
  300     DO 350 JY=1,2
            W2= ABS(ZA(4,JY)-ZA(3,JY))
            W3= ABS(ZA(2,JY)-ZA(1,JY))
            SW= W2 + W3
            IF (SW.EQ.0) GOTO 310
            WX2= W2/SW
            WX3= W3/SW
            GOTO 320
  310       WX2= 0.5
            WX3= 0.5
  320       ZX(JY)= WX2*ZA(2,JY) + WX3*ZA(3,JY)
            W2= ABS(ZB(JY+3)-ZB(JY+2))
            W3= ABS(ZB(JY+1)-ZB(JY))
            SW= W2 + W3
            IF (SW.EQ.0.) GOTO 330
            WY2= W2/SW
            WY3= W3/SW
            GOTO 340
  330       WY2= 0.5
            WY3= 0.5
  340       ZY(JY)= WY2*ZB(JY+1) + WY3*ZB(JY+2)
            ZXY(JY)= WY2*(WX2*ZAB(1,JY) + WX3*ZAB(2,JY))+
     A               WY3*(WX2*ZAB(1,JY+1) + WX3*ZAB(2,JY+1))
  350     CONTINUE
          IF (IXM1.EQ.0) GOTO 380
C
C
C         DEFINITION OF COORDINATES FOR INTERFACE TO FARBRC
          XX(1)= X4
          XX(2)= X3
          XX(3)= X3
          XX(4)= X4
          YY(1)= Y4
          YY(2)= Y4
          YY(3)= Y3
          YY(4)= Y3
          ZZ(1)= Z(IX,IY)
          ZZ(2)= Z(IX-1,IY)
          ZZ(3)= Z33
          ZZ(4)= Z43
          ZZX(1)= ZX(2)
          ZZY(1)= ZY(2)
          ZZXY(1)= ZXY(2)
          ZZX(2)= ZX34
          ZZY(2)= ZY34
          ZZXY(2)= ZXY34
          ZZX(3)= ZX33
          ZZY(3)= ZY33
          ZZXY(3)= ZXY33
          ZZX(4)= ZX(1)
          ZZY(4)= ZY(1)
          ZZXY(4)= ZXY(1)
C
C
          NSIDES= 3
          IF (IX.EQ.2) NSIDES= 4
          IFA= IFA +1
          CALL FARBRC (XX,YY,ZZ,ZZX,ZZY,ZZXY,CN,ICOL,NC,MODE,NSIDES)
  380   CONTINUE
        CALL FRBFCL(ICOL)
  390 CONTINUE
C
C     NORMAL EXIT
      RETURN
C
C     ERROR EXIT
  400 WRITE (*,99999)
      GOTO 600
C
  410 WRITE (*,99998)
      GOTO 600
C
  460 WRITE (*,99993)
      GOTO 480
  470 WRITE (*,99992)
  480 WRITE (*,99991) IX, X(IX)
      GOTO 600
C
  490 WRITE (*,99990)
      GOTO 510
  500 WRITE (*,99989)
  510 WRITE (*,99988)IY,Y(IY)
      GOTO 600
C
  530 WRITE (*,99973)
      GOTO 550
  540 WRITE (*,99972)
  550 WRITE (*,99971) I, CN(I)
      GOTO 600
C
  560 WRITE (*,99970) NC
C
  600 CONTINUE
      RETURN
C     FORMAT STATEMENTS
99999 FORMAT (1X/23H  ***   LX = 1 OR LESS./)
99998 FORMAT (1X/23H  ***   LY = 1 OR LESS./)
99993 FORMAT (1X/27H  ***   IDENTICAL X VALUES./)
99992 FORMAT (1X/33H  ***   X VALUES OUT OF SEQUENCE./)
99991 FORMAT (7H   IX= , I6, 10X, 7HX(IX) =, E12.3)
99990 FORMAT (1X/27H  ***   IDENTICAL Y VALUES./)
99989 FORMAT (1X/33H  ***   Y VALUES OUT OF SEQUENCE./)
99988 FORMAT (7H   IY= , I6, 10X, 7HY(IY) =, E12.3)
99973 FORMAT (1X/28H  ***   IDENTICAL CN VALUES./)
99972 FORMAT (1X/33H  ***  CN VALUES OUT OF SEQUENCE./)
99971 FORMAT (7H    I= , I6, 10X, 7HCN(I) =, E12.3)
99970 FORMAT (1X/18H  ***  NC .LE. 0  /,4H NC=,I10/)
      END
      SUBROUTINE FAR2D (X,LX,Y,LY,Z,NXDIM,CN,ICOL,NC)
C
      DIMENSION X(LX+1),Y(LY+1),Z(NXDIM,LY),CN(NC),ICOL(NC+1)
C
C     FILL AREA OF RECTANGLES FOR A 2D-ARRAY
C     --------------------------------------
C     FARB-E-2D  VERSION 2.1, 10/1988
C
C     AUTHOR: A. PREUSSER
C             FRITZ-HABER-INSTITUT DER MPG
C             FARADAYWEG 4-6
C             D-1000 BERLIN 33
C
C     INPUT PARAMETERS
C     X       ARRAY OF LENGTH LX+1 FOR X-COORDINATES OF
C             A REGULAR GRID
C             IN ASCENDING ORDER.
C     LX      NUMBER OF VALUES IN X-DIRECTION
C     Y       ARRAY OF LENGTH LY+1 FOR Y-COORDINATES
C             IN ASCENDING ORDER.
C     LY      NUMBER OF VALUES IN Y-DIRECTION
C     Z       2-DIMENSIONAL ARRAY DIMENSIONED Z(NXDIM,...)
C             DEFINING THE Z-VALUES FOR THE RECTANGLES DEFINED
C             BY THE GRID LINES X= X(K), Y= Y(L), X= X(K+1),
C             Y= Y(L+1),  K=1,LX, L=1,LY.
C             RECTANGLE K,L RECEIVES VALUE Z(K,L).
C     NXDIM   FIRST DIMENSION OF ARRAY Z
C     CN      ARRAY OF LENGTH NC FOR THE Z-VALUES OF
C             THE LEVELS SEPARATING AREAS OF DIFFERENT
C             COLOURS (IN ASCENDING ORDER).
C     ICOL    INTEGER ARRAY OF LENGTH NC+1 FOR
C             THE COLOURS TO BE USED FOR THE LINES OR AREAS.
C             VALUES FROM THIS ARRAY ARE PASSED TO
C             THE USER SUPPLIED SUBROUTINE USRPLT.
C             ICOL(I) IS USED FOR THE RECTANGLE, WHERE
C                  Z(K,L) .GT. CN(I-1)        AND
C                  Z(K,L) .LE. CN(I),
C             FOR I=2,NC.
C             RECTANGLES, WHERE Z(L,K).LE.CN(1)
C             ARE FILLED WITH COLOUR ICOL(1),
C             AND AREAS, WHERE Z(L,K).GT.ICOL(NC)
C             ARE FILLED WITH COLOUR ICOL(NC+1).
C     NC      NUMBER OF CONTOUR LEVELS
C
C
C     OUTPUT
C     IS PERFORMED BY CALLS TO THE SUBROUTINE    USRPLT
C     TO BE SUPPLIED BY THE USER (AN EXAMPLE FOR USRPLT
C     IS INCLUDED AND WILL BE USED IN CASE THE USER DOES
C     NOT SUPPLY HIS OWN ROUTINE).
C
C     PARAMETERS OF USRPLT
C                SUBROUTINE USRPLT (X,Y,N,NCOL,MODE)
C                X,Y     REAL ARRAYS OF LENGTH N FOR
C                        THE COORDINATES OF THE POLYGON
C                        TO BE PLOTTED.
C                N       NUMBER OF POINTS OF POLYGON
C                NCOL    COLOUR TO BE USED
C                        FOR THE AREA OR THE LINE.
C                        FOR NCOL, THE PROGRAM PASSES
C                        VALUES OF ICOL AS DESCRIBED ABOVE.
C                MODE    1, LINE DRAWING
C                        0, FILL AREA
C     -------------------------------------------------------------
C
      DIMENSION XX(4),YY(4)
      COMMON /FRBCOB/ NFABU,NCOLBU,XFABU(4),YFABU(4)
C     INITIALIZE FILL AREA BUFFER (SET TO CLOSED)
      NFABU= 0
C     INITIALIZE COLOUR OF FILL AREA BUFFER
      NCOLBU= 0
C
C     ERROR CHECKS
C
      IF(LX.LT.1) GOTO 3400
      IF(LY.LT.1) GOTO 3410
*
      DO 20 IX=2,LX+1
C        IF (X(IX-1)-X(IX)) 20,3460,3470
        IF (X(IX-1) .lt. X(IX)) then
            go to 20
        ELSE IF (X(IX-1) .eq. X(IX)) then
            go to 3460
        ELSE
            go to 3470
        ENDIF
   20 CONTINUE
C
      DO 70 IY=2,LY+1
C        IF (Y(IY-1)-Y(IY)) 70,3490,3500
        IF (Y(IY-1) .lt. Y(IY)) then
            go to 70
        ELSE IF (Y(IY-1) .eq. Y(IY)) then
            go to 3490
        ELSE
            go to 3500
        ENDIF
   70 CONTINUE
C
      DO 80 I=2,NC
C        IF (CN(I-1)-CN(I)) 80,3530,3540
        IF (CN(I-1) .lt. CN(I))then
            go to 80
        ELSE IF (CN(I-1) .eq. CN(I))then
            go to 3530
        ELSE
            go to 3540
        ENDIF
   80 CONTINUE
C
      IF (NC.LE.0) GOTO 3560
C
      DO 2000 IY= 1,LY
        YY(3)= Y(IY)
        YY(4)= Y(IY)
        YY(1)= Y(IY+1)
        YY(2)= Y(IY+1)
C
        DO 1000 IX= 1,LX
          XX(2)= X(IX)
          XX(3)= X(IX)
          XX(1)= X(IX+1)
          XX(4)= X(IX+1)
C         DETERMINE COLOUR
          NCOL= 1
          DO 500 I=1,NC
            IF (Z(IX,IY).LE.CN(I)) GOTO 510
            NCOL= NCOL+1
  500     CONTINUE
  510     CONTINUE
C
          IF (NCOL.EQ.NCOLBU) CALL FRBFUP (XX,YY)
          IF (NCOL.NE.NCOLBU) CALL FRBFOP (XX,YY,ICOL,NCOL)
C
 1000   CONTINUE
        CALL FRBFCL(ICOL)
        NCOLBU= 0
 2000 CONTINUE
C
      RETURN
C
C     ERROR EXIT
 3400 WRITE (*,99999)
      GOTO 3600
C
 3410 WRITE (*,99998)
      GOTO 3600
C
 3460 WRITE (*,99993)
      GOTO 3480
 3470 WRITE (*,99992)
 3480 WRITE (*,99991) IX, X(IX)
      GOTO 3600
C
 3490 WRITE (*,99990)
      GOTO 3510
 3500 WRITE (*,99989)
 3510 WRITE (*,99988)IY,Y(IY)
      GOTO 3600
C
 3530 WRITE (*,99973)
      GOTO 3550
 3540 WRITE (*,99972)
 3550 WRITE (*,99971) I, CN(I)
      GOTO 3600
C
 3560 WRITE (*,99970) NC
C
 3600 CONTINUE
      RETURN
C     FORMAT STATEMENTS
99999 FORMAT (1X/23H  ***   LX = 0 OR LESS./)
99998 FORMAT (1X/23H  ***   LY = 0 OR LESS./)
99993 FORMAT (1X/27H  ***   IDENTICAL X VALUES./)
99992 FORMAT (1X/33H  ***   X VALUES OUT OF SEQUENCE./)
99991 FORMAT (7H   IX= , I6, 10X, 7HX(IX) =, E12.3)
99990 FORMAT (1X/27H  ***   IDENTICAL Y VALUES./)
99989 FORMAT (1X/33H  ***   Y VALUES OUT OF SEQUENCE./)
99988 FORMAT (7H   IY= , I6, 10X, 7HY(IY) =, E12.3)
99973 FORMAT (1X/28H  ***   IDENTICAL CN VALUES./)
99972 FORMAT (1X/33H  ***  CN VALUES OUT OF SEQUENCE./)
99971 FORMAT (7H    I= , I6, 10X, 7HCN(I) =, E12.3)
99970 FORMAT (1X/18H  ***  NC .LE. 0  /,4H NC=,I10/)
      END
      SUBROUTINE FARBRC(X,Y,Z,ZX,ZY,ZXY,CN,ICOL,NC,MODE,NSIDES)
C
      DIMENSION X(4),Y(4),Z(4),ZX(4),ZY(4),ZXY(4),CN(NC),ICOL(NC+1)
C
C     F ILL  AR EA  FOR A  B ICUBIC FUNCTION ON A  R E C TANGLE
C     *      **            *                       *   *
C
C     T R I P   ALGORITHM   A.PREUSSER   FARB-E-2D  VERSION 2.1 10/1988
C
C     AUTHOR: A. PREUSSER
C             FRITZ-HABER-INSTITUT DER MPG
C             FARADAYWEG 4-6
C             D-1000 BERLIN 33
C
C
C     THIS SUBROUTINE COMPUTES A BICUBIC FUNCTION FROM THE
C     VALUES X,Y,Z,ZX,ZY,ZXY GIVEN AT THE FOUR VERTICES OF
C     A RECTANGLE, AND PLOTS CONTOURS FOR THE Z-VALUES CN(I),
C     I=1,NC, USING THE COLOURS ICOL(I).
C     AREA FILLING, SET BY PARAMETER MODE, IS AN OPTIONAL FEATURE.
C
C     INPUT PARAMETERS
C     ================
C     X,Y,Z         COORDINATES OF THE VERTICES
C                   IN CENTIMETERS (OR INCHES, IF CMSCAL=2.54).
C                   REAL ARRAYS OF LENGTH (4).
C                   X(I),Y(I),Z(I), I=1,4 DEFINE THE
C                   POSITION OF VERTEX (I).
C                   THE SIDES OF THE RECTANGLE MUST BE PARALLEL
C                   TO THE X- AND Y-AXIS,
C                   AND THE VERTICES MUST BE ORDERED
C                   COUNTER-CLOCKWISE AS IS INDICATED BELOW
C                   (VERTEX 1 IN THE UPPER RIGHT CORNER).
C     ZX,ZY,ZXY     DERIVATIVES OF Z AT THE VERTICES.
C                   REAL ARRAYS OF LENGTH (4).
C     CN            Z- VALUES IN ASCENDING ORDER
C                   FOR THE CONTOUR LEVELS.
C                   REAL ARRAY OF LENGTH (NC)
C     ICOL          INDICES OF THE COLOURS OR PATTERNS
C                   FOR THE AREAS BETWEEN THE CONTOUR LINES.
C                   THE VALUES OF ICOL ARE PASSED TO USRPLT.
C                   ICOL(I) IS USED FOR THE AREA, WHERE
C                        Z .GT. CN(I-1)        AND
C                        Z .LE. CN(I),
C                   FOR I=2,NC.
C                   AREAS, WHERE Z.LE.CN(1)
C                   ARE FILLED WITH COLOUR ICOL(1),
C                   AND AREAS, WHERE Z.GT.ICOL(NC)
C                   ARE FILLED WITH COLOUR ICOL(NC+1).
C     NC            NUMBER OF CONTOUR LEVELS, NC.LE.100
C     MODE          0, FILL AREA ONLY
C                   1, LINES ONLY
C                   2, FILL AREA AND LINES
C     NSIDES        4, COMPUTE ZEROS FOR 4 SIDES. (SHOULD BE USED
C                      AS A DEFAULT)
C                   3, COMPUTE ZEROS FOR 3 SIDES ONLY.
C                      ZEROS FOR SIDE 4 ARE COPIED FROM SIDE 2
C                      (ONLY APPLICABLE, IF THE RECTANGLE OF
C                       THIS CALL IS THE RIGHT HAND NEIGHBOR OF
C                       THE RECTANGLE OF THE PREVIOUS CALL)
C
C     OUTPUT
C     ======
C     IS PERFORMED BY CALLS TO THE SUBROUTINE    USRPLT
C     TO BE SUPPLIED BY THE USER (AN EXAMPLE FOR USRPLT
C     IS INCLUDED).
C
C     PARAMETERS OF USRPLT
C                SUBROUTINE USRPLT (X,Y,N,NCOL,MODE)
C                X,Y     REAL ARRAYS OF LENGTH N FOR
C                        THE COORDINATES OF THE POLYGON
C                        TO BE PLOTTED.
C                N       NUMBER OF POINTS OF POLYGON
C                NCOL    INDEX  DEFINING THE COLOUR FOR
C                        THE AREA OR THE LINE
C                MODE    1, LINE DRAWING
C                        0, FILL AREA
C
C     IF A RECTANGLE RECEIVES ONLY ONE COLOUR,
C     THE AREA IS NOT FILLED AT ONCE.
C     INSTEAD, A 'FILL AREA BUFFER' IS OPENED
C     OR UPDATED, UNTIL A RECTANGLE WITH A
C     DIFFERENT COLOUR IS ENCOUNTERED.
C     THEREFORE, IF THE NEXT CALL TO FARBRC
C     IS NOT FOR A RIGHT-HAND-NEIGHBOR,
C     OR IF IT IS THE LAST CALL, SUBROUTINE
C            FRBFCL
C     MUST BE CALLED BY THE USER
C            CALL FRBFCL(ICOL)  ,
C     IN ORDER TO CLEAR THE FILL AREA BUFFER,
C     AND TO FILL THE AREA OF THE RECTANGLE.
C
C            DENOMINATION OF THE VERTICES AND SIDES OF THE
C                           RECTANGLE
C            Y
C                            SIDE(3)
C  VERTEX(2) * -------------------------------0-------- * VERTEX(1)
C            (                             .            )
C            (                           .              )
C            (                          .               )
C    SIDE(4) (                          . RIDE          ) SIDE(2)
C            (                           .              )
C            (                             .            )
C            (                                .         )
C  VERTEX(3) * ----------------------------------0----- * VERTEX(4)
C                            SIDE(1)                        X
C
C     THE SIDES ARE PARALLEL TO THE CARTESIAN X-Y-SYSTEM.
C
C   -----------------------------------------------------------------
C   END OF USER DOCUMENTATION
C   -----------------------------------------------------------------
C
C           SOME NOMENCLATURE
C
C     STATION      ZERO ON A SIDE
C     RIDE         MOVE FROM ONE STATION TO ANOTHER INSIDE RECT.
C     TRANSFER     MOVE FROM ONE STATION TO THE NEXT ON SIDE
C     TRIP         SEQUENCE OF RIDES AND TRANSFERS
C     ROUND TRIP   SUCCESSFUL TRIP THAT ENDED AT ITS START
C     HORROR TRIP  TRIP THAT DOES NOT FIND AN END
C     JOURNEY      SEQUENCE OF TRIPS STARTING FROM THE SAME
C                  TYPE OF STATIONS (SAME VALUE OF ISTATZ)
C                  AND HAVING THE SAME ORIENTATION.
C
C                  THERE MAY BE THREE JOURNEYS.
C                  THE FIRST TWO ARE COUNTER-CLOCKWISE AND
C                  START AT STATIONS WITH ISTATZ=0 AND =2,
C                  RESPECTIVELY. THE THIRD JOURNEY IS CARRIED
C                  OUT ONLY IN CASE OF NUMERICAL DIFFICULTIES,
C                  WHEN AREAS ARE UNFILLED AFTER THE FIRST TWO.
C                  IT STARTS AT STATIONS WITH ISTATZ=1 OR =0 AND
C                  IS CLOCKWISE.
C
C
C
      COMMON /FRBCOB/ NFABU,NCOLBU,XFABU(4),YFABU(4)
C     NFABU       0, FILL AREA BUFFER CLOSED
C                 1, FILL AREA BUFFER OPEN
C     NCOLBU      COLOUR OF FILL AREA BUFFER
C     XFABU,YFABU X-Y COORDINATES OF FILL AREA BUFFER
C
      COMMON /FRBCOC/ SACMIN,CMSCAL,MAXPOL,NCPMAX,MAXSTA,NCMAXS,NPP,PI
     1,               MAXRID, SIGS(4), NCMAX
C     SACMIN     MINIMAL DISTANCE OF TWO POINTS TO BE STORED
C *** CMSCAL     VARIABLE FOR SWITCHING BETWEEN CM AND INCH
C     MAXPOL     MAXIMUM NUMBER OF POINTS FOR A TRIP
C     NCPMAX     MAXIMUM NUMBER OF POINTS TO BE COMPUTED FOR A RIDE
C     MAXSTA     MAXIMUM NUMBER OF POINTS TO BE STORED FOR A RIDE
C     NCMAXS     MAXIMUM NUMBER OF CONTOURS CROSSING A RECTANGLE SIDE
C     NPP        ACCUMULATED NUMBER OF POINTS FOR A SEQUENCE OF
C                RECTANGLES
C     PI         3.141...
C     MAXRID     MAXIMUM NUMBER OF RIDES FOR A TRIP
C     SIGS       SIGN FOR SIDES (+1 OR -1)
C     NCMAX      MAXIMUM NUMBER OF CONTOUR LEVELS
C
      COMMON /FRBCOF/ KK,KSE,XX4F,YY4F,SIR,COR,CL
C
C     /FRBCOF/ CONTAINS VARIABLES WHICH ARE PASSED TO FUNCTION
C              FRBEVA AS PARAMETERS
C
C     KK          INDEX  OF FUNCTION TO BE EVALUATED BY FRBEVA
C     KSE         ACTUAL SIDE INDEX
C     XX4F,YY4F   COORDINATES FOR POINT P4F (PRELEMINARY POSITION
C                 OF POINT P4)
C     SIR,COR     COSINUS OF DIRECTION NORMAL TO CURVE DIRECTION
C     CL          ACTUAL CONTOUR LEVEL
C
      COMMON /FRBCOP/  P0(4),P1(4),P2(4),P3(4)
     1,     Q0(4),Q1(4),Q2(4)
     2,     R0(4),R1(4)
     3,     P11,P12,P13,P21,P22,P23,P31,P32,P33
C
C      P0,P1,P2,P3    COEFFICIENTS FOR THE POLYNOMIALS
C                     ON THE 4 SIDES .
C                     VARIABLES ON SIDES 1 AND 2 ARE
C                     COUNTER-CLOCKWISE, ON SIDES 3 AND 4
C                     CLOCKWISE.
C      Q0,Q1,Q2       COEFFICIENTS FOR THE DERIVATIVES
C                     OF THE POLYNOMIALS ON THE 4 SIDES.
C      R0,R1          COEFF. FOR THE SECOND DERIVATIVES.
C      P11...P33      COEFF. OF POLYNOMIALS USED TOGETHER
C                     WITH P0(I)...P3(I), I=1 AND I=4, FOR THE
C                     REPRESENTATION INSIDE THE RECTANGLE.
C
      COMMON /FRBCRD/ X0(50,4),Y0(50,4),NCLZR(50,4),SDER(50,4),TZR(50,4)
     A,               NZ(4),SI(4),CO(4),SA(4),SE(4),DX(4),DY(4),SL(4)
     B,               HMIN,SLMAX,KRIDE,NPREC,POSERR,DERNOR,NDIR3
     C,               ZMAX(4),ZMIN(4),ISTATZ(50,4),X3,Y3,ZSOLD
C
C     FRBCRD CONTAINS VARIABLES THAT ARE PASSED TO FRBRID
C            OR THAT ARE RETAINED FOR THE NEXT CALL TO FARBRC (NSIDE=3)
C     X0,Y0       X-Y-COORDINATES OF ZEROS ON SIDES
C     NCLZR       CONTOUR LEVEL FOR ZEROS
C     SDER        DERIVATIVE IN DIRECTION OF SIDES
C                 (SIDE DIRECTION= COUNTER CLOCKWISE)
C     TZR         COORDINATES FOR ZEROS (STATIONS) ON SIDES
C     NZ          NUMBER OF ZEROS ON SIDES
C     SI,CO       COSINUS OF DIRECTION FOR SIDES
C     SA,SE       VALUES OF VARIABLES AT START AND END OF SIDES
C     DX,DY       DIFFERENCES OF X AND Y
C     SL          SIDE LENGTHS
C     HMIN        LENGTH OF SHORTEST SIDE OF RECTANGLE
C     SLMAX       LENGTH OF LONGEST SIDE
C     KRIDE       COUNTS THE CALLS TO FRBRID
C     NPREC       NUMBER OF POINTS FOR THE RECTANGLE
C     POSERR      PERMITTED POSITION ERROR
C     DERNOR      SIGN OF NORMAL DERIVATIVE FOR A RIDE (+1 OR -1)
C     NDIR3       1, FOR COUNTER-CLOCKWISE TRIP
C                 -1, FOR CLOCKWISE TRIP
C     ZMIN,ZMAX    MIN. AND MAX. VALUES OF Z ON SIDES
C     ISTATZ       STATUS OF ZEROS ON SIDES
C                  THE INITIAL STATUS IS = 0
C                  WHEN A ZERO HAS SERVED AS START +1 IS ADDED
C                  WHEN A ZERO HAS SERVED AS END +2 IS ADDED
C     X3,Y3        COORDINATES X(3),Y(3)
C     ZSOLD        ZS OF LAST CALL
C
      DIMENSION CC(100),ZZ(4),T1(4,4),Z1(4,4),TS2(4),IN(4)
     1,         XPOL(300),YPOL(300),XSTACK(100),YSTACK(100)
     2,         XX(2),YY(2)
     3,         JSAR(6),JZAR(6),JSER(6),JZER(6),JSTOPR(6)
     4,         SL1(4),SL12(4)
C
C     CC           SCALED CONTOUR LEVELS
C     ZZ           SCALED Z-VALUES AT VERTICES
C     T1           COORDINATES ON SIDES AT ENDPOINTS OF INTERVALS
C     Z1           Z-VALUES AT ENDPOINTS OF INTERVALS
C     TS2          WORKING ARRAY FOR COMPUTING T1
C     IN           NUMBER OF INTERVALS ON SIDES
C     XPOL,YPOL    X,Y-COORDINATES FOR FILL AREA POLYGON
C     XSTACK,YSTACK X,Y COORDINATES OF STACK
C                  (COORDINATES OF THE FIRST RIDE OF A TRIP)
C     XX,YY        LOCAL COPIES OF X,Y
C     JSAR,JZAR    SIDE AND ZERO FOR START OF RIDES
C     JSER,JZER    SIDE AND ZERO FOR END OF RIDES
C     JSTOPR       STOP MODES OF RIDES
C     SL1          1./SIDE LENGTH
C     SL12         SL1**2
C
C
      SAVE IT
      DATA XSTACK,YSTACK /200*0./
      DATA IT /0/, ZS/0./
C
C
C     INITIALISATION FOR FIRST RECTANGLE
      IT= IT + 1
C
C     SET INSTALLATION PARAMETERS
C
      IF (IT.NE.1) GOTO 10
      CMSCAL= 1.
C *** SET  CMSCAL= 2.54  FOR INCH CALIBRATED PLOTTERS
      SACMIN= 0.02/CMSCAL
      MAXPOL= 300
      NCMAXS= 50
      NCMAX= 100
      MAXSTA= 100
      MAXRID= 6
      NCPMAX= MAXSTA*10
      NPP= 0
      NFABU= 0
      PI= 4.*ATAN(1.)
      SIGS(1)=  1.
      SIGS(2)=  1.
      SIGS(3)= -1.
      SIGS(4)= -1.
   10 CONTINUE
C
C     CHECK NUMBER OF CONTOUR LEVELS
      IF (NC.GT.NCMAX) WRITE (*,8999) NC,NCMAX
 8999 FORMAT ('0***ERROR*** IN FARBRC'/
     1        ' NUMBER OF CONTOURLEVELS NC=',I10/
     2        ' .GT. NCMAX=',I10)
C
      NFAR= 0
C         = NUMBER OF FILL AREA CALLS
      KRIDE= 0
C       = NUMBER OF CALLS TO FRBRID
      NPREC= 0
C          = NUMBER OF CURVE POINTS COMPUTED FOR RECTANGLE
C
C     CHECK CONTOUR VALUES FOR MONOTONY,
C     IF CONTOUR LINE PASSES THROUGH DATA POINT ON
C     SIDE 4, SET NSIDE TO 4 (NO COPY FROM LAST RECTANGLE)
C     SCALE CONTOUR LEVELS
      NSIDE= NSIDES
      IF (CN(1).EQ.Z(2) .OR.
     A    CN(1).EQ.Z(3)) NSIDE= 4
      ZSOLD= ZS
      ZS= (Z(1)+Z(2)+Z(3)+Z(4))/4.
      CC(1)= CN(1) - ZS
      DO 20 KCL=2,NC
        CC(KCL)= CN(KCL) - ZS
        IF (CN(KCL).EQ.Z(2) .OR.
     A      CN(KCL).EQ.Z(3)) NSIDE= 4
        IF (CN(KCL)-CN(KCL-1).GT.0.) GOTO 20
        WRITE (*,9000) KCL,CN(KCL),KCL-1,CN(KCL-1)
   20 CONTINUE
 9000 FORMAT ('0***ERROR*** IN FARBRC'/
     A        '0     CONTOURLEVEL', I5,'=',E15.7/
     B        '0.LE. CONTOURLEVEL', I5,'=',E15.7)
C
C
C     SOME BASIC GEOMETRY FOR THE RECTANGLE
C
      SL(1)= X(4)-X(3)
      SL(2)= Y(1)-Y(4)
      SL(3)= X(1)-X(2)
      SL(4)= Y(2)-Y(3)
      X3= X(3)
      Y3= Y(3)
      DO 50 J=1,4
        ZZ(J)= Z(J) - ZS
        NP1= MOD(J+1,4) + 1
        NP2= MOD(J+2,4) + 1
        DX(J)= X(NP2) - X(NP1)
        DY(J)= Y(NP2) - Y(NP1)
        SL1(J)= 1./SL(J)
        SL12(J)= SL1(J)*SL1(J)
        CO(J)= DX(J)/SL(J)
        SI(J)= DY(J)/SL(J)
        SA(J)= 0.
        IF (J.GT.2) SA(J)= SL(J)
        SE(J)= SL(J)
        IF (J.GT.2) SE(J)= 0.
   50 CONTINUE
      SLMAX= AMAX1(SL(1),SL(2))
      DI2=  - (DY(1)*CO(2) - DX(1)*SI(2))
      DI3= ABS(DY(2)*CO(3) - DX(2)*SI(3))
C
C     CHECK COORDINATES OF VERTICES
      IF (Y(1).NE.Y(2) .OR. X(2).GE.X(1))
     1   WRITE (*,9010) IT
 9010 FORMAT ('0***ERROR***'/
     1        ' Y(1).NE.Y(2)  OR   X(2).GE.X(1)'/
     2        ' IN RECTANGLE NO.', I10,
     3        ' VERTICES MUST BE ORDERED COUNTER-CLOCKWISE'/
     4        ' STARTING IN THE UPPER RIGHT CORNER.'/
     5        ' SIDES MUST BE PARALLEL TO X- AND Y-AXIS')
C
C     CHECK IF VERTICES ARE NUMBERED COUNTER-CLOCKWISE
      IF (DI2.LT.0.) WRITE (*,9020) IT
 9020 FORMAT ('0***ERROR***, VERTICES OF RECTANGLE NO.',I10,
     1        ' NOT IN COUNTER-CLOCKWISE ORDER')
C
      HMIN= AMIN1(DI2,DI3)
C         = SHORTEST SIDE LENGTH
C
C     CHECK HMIN
c      IF (HMIN.LT.0.01/CMSCAL .OR. HMIN.GT.100./CMSCAL)
      if (hmin.lt.0.0001/cmscal .or. hmin.gt.10000./cmscal)
     1   WRITE (*,9030) IT,HMIN
      IF (HMIN.EQ.0.) GOTO 5000
 9030 FORMAT ('0***WARNING***, CHECK RECTANGLE NO.',I10/
     1 ' DIFFERENCE IN X- OR Y- COORDINATES TOO LARGE OR TOO SMALL'/
     2 ' XYDIF=', E20.5/
     3 ' SCALE X AND/OR Y TO CM (OR INCH)')
C
      POSERR= AMIN1(1.E-03/CMSCAL,1.E-03*HMIN*HMIN/SLMAX)
C           = PERMITTED POSITION ERROR
C
C     COPY INFORMATION FOR SIDE 4
      IF (NSIDE.EQ.4) GOTO 80
      JZN= NZ(2)
      NZ(4)= JZN
      ZMIN(4)= ZMIN(2) + ZSOLD - ZS
      ZMAX(4)= ZMAX(2) + ZSOLD - ZS
      IF (JZN.EQ.0) GOTO 80
      JJZ= JZN
      DO 70 JZR=1,JZN
        X0(JJZ,4)= 0.
        Y0(JJZ,4)= Y0(JZR,2)
        NCLZR(JJZ,4)= NCLZR(JZR,2)
        SDER(JJZ,4)= -SDER(JZR,2)
        TZR(JJZ,4)= TZR(JZR,2)
        ISTATZ(JJZ,4)= 0
        JJZ= JJZ - 1
   70 CONTINUE
   80 CONTINUE
C
C     COMPUTE COEFFICIENTS FOR POLYNOMIALS ALONG SIDES
C
          Z3A3= (ZZ(4)- ZZ(3))*SL1(1)
          P0(1)= ZZ(3)
          P1(1)= ZX(3)
          P2(1)= (2.0*(Z3A3-ZX(3))+Z3A3-ZX(4))*SL1(1)
          P3(1)= (-2.*Z3A3+ZX(4)+ZX(3))*SL12(1)
C
          Z3A3= (ZZ(2)- ZZ(3))*SL1(4)
          P0(4)= ZZ(3)
          P1(4)= ZY(3)
          P2(4)= (2.*(Z3A3-ZY(3))+Z3A3-ZY(2))*SL1(4)
          P3(4)= (-2.*Z3A3+ZY(2)+ZY(3))*SL12(4)
C
          Z3A3= (ZZ(1)- ZZ(2))*SL1(3)
          P0(3)= ZZ(2)
          P1(3)= ZX(2)
          P2(3)= (2.*(Z3A3-ZX(2)) + Z3A3 - ZX(1))*SL1(3)
          P3(3)= (-2.*Z3A3+ZX(1)+ZX(2))*SL12(3)
C
          Z3A3= (ZZ(1)- ZZ(4))*SL1(2)
          P0(2)= ZZ(4)
          P1(2)= ZY(4)
          P2(2)= (2.*(Z3A3-ZY(4)) + Z3A3 - ZY(1))*SL1(2)
          P3(2)= (-2.*Z3A3+ZY(1)+ZY(4))*SL12(2)
C
C     DETERMINE POLYNOMIAL COEFF FOR DERIVATIVES ALONG SIDES
            DO  90 J= 1,4
              Q0(J)=    P1(J)
              Q1(J)= 2.*P2(J)
              Q2(J)= 3.*P3(J)
              R0(J)=    Q1(J)
              R1(J)= 2.*Q2(J)
   90       CONTINUE
C
C     SET CONTOUR LEVEL TO BE PASSED TO FRBEVA
      CL= 0.
C
C
C     FIND POINTS ON SIDES OF RECTANGLE,
C     WHERE FIRST DERIVATIVE IS ZERO
C
C     LOOP OVER SIDES
      DO 200 JSA=1,NSIDE
        KSE= JSA
C       SET INITIAL ENDPOINTS OF INTERVALS
        T1(1,KSE)= SA(KSE)
        T1(2,KSE)= SE(KSE)
        I= 2
C       LOOP OVER DERIVATIVES
        DO 150 K=3,4
C
C         SET FUNCTION TO BE EVALUATED BY FRBEVA
          KK= K
C
          TS2(1)= T1(1,KSE)
          II= 2
          TB= T1(1,KSE)
          F2= FRBEVA(TB)
C         LOOP OVER ENDPOINTS OF INTERVALS
          DO 100 J=2,I
            TA= TB
            F1= F2
            TB= T1(J,KSE)
            F2= FRBEVA(TB)
            IF (F1*F2.GT.0.) GOTO 100
            IF (F1.EQ.0. .AND. F2.EQ.0.) GOTO 100
            TS2(II)= FRBZER(TA,TB,F1,F2,POSERR)
            II= II + 1
  100     CONTINUE
          TS2(II)= T1(I,KSE)
          DO 120 J=1,II
  120       T1(J,KSE)= TS2(J)
          I= II
  150   CONTINUE
C       IN(KSE)= NUMBER OF INTERVALS
        I= I-1
        IN(KSE)= I
C       (E.G. IF IN(KSE)=1, THERE IS NO POINT FOR WHICH 1ST DER.=0)
C
C       COMPUTE MAXIMA AND MINIMA FOR EACH SIDE
        NP1= MOD(KSE+1,4) + 1
        NP2= MOD(KSE+2,4) + 1
        ZMAX(KSE)= AMAX1(ZZ(NP1),ZZ(NP2))
        ZMIN(KSE)= AMIN1(ZZ(NP1),ZZ(NP2))
        Z1(1,KSE)= ZZ(NP1)
        Z1(I+1,KSE)= ZZ(NP2)
        IF (I.EQ.1) GOTO 170
        KK= 1
        DO 160 J=2,I
          Z1(J,KSE)= FRBEVA(T1(J,KSE))
          IF (Z1(J,KSE).GT.ZMAX(KSE)) ZMAX(KSE)= Z1(J,KSE)
          IF (Z1(J,KSE).LT.ZMIN(KSE)) ZMIN(KSE)= Z1(J,KSE)
  160   CONTINUE
  170   CONTINUE
  200 CONTINUE
C
C     CHECK, IF RECTANGLE HAS ONE COLOUR,
C     BECAUSE THE MINIMUM IS OVER THE MAX. CONTOUR LEVEL,
C     OR MAXIMUM UNDER MIN. CONTOUR LEVEL
      ZMAXT= AMAX1(ZMAX(1),ZMAX(2),ZMAX(3),ZMAX(4))
      ZMINT= AMIN1(ZMIN(1),ZMIN(2),ZMIN(3),ZMIN(4))
      CN1= CC(1)
      CNN= CC(NC)
      IF (CNN.GE.ZMINT .AND. CN1.LT.ZMAXT) GOTO 500
      NCOL =1
      IF (CNN.LT.ZMINT) NCOL= NC+1
C
C     RECTANGLE HAS ONE COLOUR ONLY
  250 CONTINUE
      IF (MODE.EQ.1) GOTO 400
      IF (NFABU.EQ.1 .AND. NCOL.EQ.NCOLBU) GOTO 300
      CALL FRBFOP(X,Y,ICOL,NCOL)
      GOTO 400
  300 CONTINUE
      CALL FRBFUP (X,Y)
  400 CONTINUE
      NZ(2)= 0
      GOTO 5000
  500 CONTINUE
C
C     FIND MIN. AND MAX. CONTOUR LEVEL FOR RECTANGLE
      DO 600 KCL=1,NC
        IC1= KCL
        IF (CC(KCL).GE.ZMINT) GOTO 610
  600 CONTINUE
  610 CONTINUE
      ICN= NC + 1
      DO 650 KCL=1,NC
        ICN= ICN - 1
        IF (CC(ICN).LE.ZMAXT) GOTO 660
  650 CONTINUE
  660 CONTINUE
C     IC1= FIRST CONTOUR LEVEL
C     ICN= LAST CONTOUR LEVEL
      ICN1= ICN - IC1 + 1
      NCOL= ICN +1
      IF (CC(ICN).EQ.ZMAXT) NCOL= ICN
      IF (ICN1.EQ.0) GOTO 250
C
C     COMPUTE ZEROS ON SIDES FOR ALL CONTOUR LEVELS
C     IN COUNTER-CLOCKWISE ORDER
C
      DO 900 JSA=1,NSIDE
        KSE= JSA
        JN= 0
        IF (CC(ICN).LT.ZMIN(KSE) .OR. CC(IC1).GT.ZMAX(KSE))
     A     GOTO 850
        NI= IN(KSE)
        DO  800 JIN=1,NI
C
          NDIR= SIGN(1.,Z1(JIN+1,KSE)-Z1(JIN,KSE))
          KCL= IC1 - 1
          IF (NDIR.EQ.-1) KCL= ICN + 1
C
          DO 700 KCLL=1,ICN1
            KCL= KCL + NDIR
            CL= CC(KCL)
            F1= Z1(JIN,KSE) - CL
            F2= Z1(JIN+1,KSE) - CL
            F1F2= F1*F2
            IF (F1F2.GT.0.) GOTO  700
            IF (F1F2.LT.0.) GOTO  690
C
C           SPECIAL SITUATIONS
C
            IF (NI.EQ.1 .AND. F1.EQ.0. .AND. F2.EQ.0.) GOTO  670
C           IF () THEN CONTOURLINE = SIDE KSE
C
            IF (F1.EQ.0. .AND.
     1         ABS(T1(JIN,KSE)-T1(1,KSE)).LE.POSERR*3.) GOTO 700
C           IF () THEN LINE PASSES THROUGH A VERTEX AT START OF SIDE
C                 THIS CASE IS HANDLED ON PREVIOUS SIDE
C
            IF (F2.EQ.0. .AND.
     1        ABS(T1(JIN+1,KSE)-T1(NI+1,KSE)).LE.POSERR*3.) GOTO 680
C           IF () THEN LINE PASSES THROUGH VERTEX AT END OF SIDE
            GOTO  690
C
C           CONTOUR LINE = SIDE JSA
  670       CONTINUE
            KVERT= MOD(JSA+1,4) + 1
            XX(1)= X(KVERT)
            YY(1)= Y(KVERT)
            KVERT= MOD(JSA+2,4) + 1
            XX(2)= X(KVERT)
            YY(2)= Y(KVERT)
            IF (MODE.GT.0)
     1      CALL USRPLT(XX,YY,2,KCL,1)
            GOTO 850
  680       CONTINUE
C
C           LINE PASSES THROUGH DATA POINT
C                       (AT END OF SIDE)
C
C           INHIBIT MULTIPLE ZERO AT VERTEX
            IF (JN.EQ.0) GOTO 685
            IF (KCL.NE.NCLZR(JN,KSE)) GOTO 685
            IF (ABS(T1(NI+1,KSE)-TZR(JN,KSE)).LE.POSERR*3.) GOTO 700
  685       CONTINUE
C
C           COMPUTE VALUE ON SIDE JSA
            KK= 1
            EPS= 0.01*SL(KSE)
            TA= T1(NI+1,KSE)-EPS*SIGS(KSE)
            FA= FRBEVA(TA)
C
C           COMPUTE VALUE ON NEXT SIDE
            NSE= MOD(KSE,4) + 1
            KSE= NSE
            EPSN= 0.01*SL(KSE)
            TB= T1(1,NSE)+ EPSN*SIGS(NSE)
            FB= FRBEVA(TB)
C
            KSE= JSA
            IF (FA*FB.GT.0.) GOTO 700
C           IF () THEN CONTOUR LINE IS DEGENERATED TO A POINT
C
C           CONTOUR LINE STARTS FROM VERTEX INTO RECTANGLE
            JN= JN+ 1
            IF (JN.GT.NCMAXS) GOTO 6000
            TZR(JN,KSE)= T1(NI+1,KSE)
            NCLZR(JN,KSE)= KCL
            SDER(JN,KSE)= -FA/EPS
            GOTO 700
C
C           COMPUTE ZERO ON SIDE (STATION)
  690       JN= JN+1
            IF (JN.GT.NCMAXS) GOTO 6000
            KK= 1
            TZR(JN,KSE)= FRBZER(T1(JIN,KSE),T1(JIN+1,KSE),F1,F2,
     A                   POSERR)
C
C           COMPUTE DERIVATIVE AT ZERO
            KK= 4
            SDER(JN,KSE)= FRBEVA(TZR(JN,KSE))*SIGS(KSE)
C           STORE INDEX OF CONTOUR LEVEL
            NCLZR(JN,KSE)= KCL
C
C           CHECK SIGN OF DER., IF SAME LEVEL
            IF (JN.LT.2) GOTO 700
            IF (KCL.NE.NCLZR(JN-1,KSE)) GOTO 700
            IF (ABS(TZR(JN,KSE)-TZR(JN-1,KSE)).GT.POSERR*3.) GOTO 700
            IF (SDER(JN,KSE)*SDER(JN-1,KSE) .LT. 0.) GOTO 700
C
C           IF SIGN IS WRONG OR =0, COMPUTE DER. BY DIFFERENCES
            KK= 1
            EPS= 0.01*SL(KSE)
            TA= TZR(JN-1,KSE) - EPS*SIGS(KSE)
            SDER(JN-1,KSE)= -FRBEVA(TA)/EPS
            TB= TZR(JN,KSE) + EPS*SIGS(KSE)
            SDER(JN,KSE)= FRBEVA(TB)/EPS
C
  700     CONTINUE
  800   CONTINUE
  850   NZ(KSE)= JN
C                = NUMBER OF ZEROS ON SIDE KSE
  900 CONTINUE
C
C     EVERY RIDE SHOULD HAVE START AND END
      IF (NZ(1)+NZ(2)+NZ(3)+NZ(4).LT.2) GOTO 250
C
C
C     CLEAR FILL AREA BUFFER
      CALL FRBFCL (ICOL)
C
C     COMPUTE X0,Y0 FOR EACH ZERO (RELATIVE TO X(3),Y(3)),
C     SET STATUS OF ALL ZEROS TO 0
      DO 1300 JSA=1,NSIDE
        JN= NZ(JSA)
        IF (JN.EQ.0) GOTO 1300
        NP1= MOD(JSA+1,4) + 1
        DO 1280 JZA=1,JN
          ISTATZ(JZA,JSA)= 0
          T= TZR(JZA,JSA)
          X0(JZA,JSA)= X(NP1) -X3 + SA(JSA)*CO(JSA) +T*ABS(CO(JSA))
          Y0(JZA,JSA)= Y(NP1) -Y3 + SA(JSA)*SI(JSA) +T*ABS(SI(JSA))
 1280   CONTINUE
 1300 CONTINUE
C
C      COMPUTE COEFFICIENTS FOR REPRESENTATION INSIDE RECTANGLE
          ZX3B3= (ZX(2)-ZX(3))*SL1(4)
          ZX4B3= (ZX(1)-ZX(4))*SL1(4)
          ZY3A3= (ZY(4)-ZY(3))*SL1(1)
          ZY4A3= (ZY(1)-ZY(2))*SL1(1)
          A= (ZZ(1)-ZZ(4)-ZZ(2)+ZZ(3))*SL1(4)*SL1(1)
     A       - ZX3B3 - ZY3A3 + ZXY(3)
          B= ZX4B3 - ZX3B3 - ZXY(4) + ZXY(3)
          C= ZY4A3 - ZY3A3 - ZXY(2) + ZXY(3)
          D= ZXY(1) - ZXY(4) - ZXY(2) + ZXY(3)
          E= A+A-B-C
          P11= ZXY(3)
          P12= (2.*(ZX3B3-ZXY(3))+ZX3B3-ZXY(2))*SL1(4)
          P13= (-2.*ZX3B3+ZXY(2)+ZXY(3))*SL12(4)
          P21= (2.*(ZY3A3-ZXY(3))+ZY3A3-ZXY(4))*SL1(1)
          P22= (3.*(A+E)+D)*SL1(1)*SL1(4)
          P23= (-3.*E-B-D)*SL1(1)*SL12(4)
          P31= (-2.*ZY3A3+ZXY(4)+ZXY(3))*SL12(1)
          P32= (-3.*E-C-D)*SL1(4)*SL12(1)
          P33= (D+E+E)*SL12(1)*SL12(4)
C
C     INITIALIZE STACK (FOR FIRST RIDE OF A TRIP)
      JSE1ST= 0
      JZE1ST= 0
      JSA1ST= 0
      JZA1ST= 0
      NP1ST=  0
      DERNOS= 0
C
C
C     START 'TRIPS' USING THE ZEROS ON THE SIDES AS 'STATIONS'
C
C     SET PARAMETERS FOR FIRST JOURNEY
C     START AT UNUSED STATIONS (ISTATZ =0)
      ISTART= 0
      NDIR3= 1
C     NDIR3= 1 MEANS COUNTER-CLOCKWISE TRIP
      NDIRV= 2
      NDIRS= 0
C
C     LOOP OVER JOURNEYS
      NJOUR= 3
      IF (MODE.EQ.1) NJOUR=1
      DO 4000 JOURNY= 1,NJOUR
C
C       SET PARAMETERS FOR JOURNY=3
        IF (JOURNY.NE.3) GOTO 1310
        ISTART= 1
        NDIR3= -1
        NDIRV= 1
        NDIRS= 2
 1310   CONTINUE
C
C       LOOP OVER SIDES, JSA= SIDE INDEX
        DO 3000 JSA=1,4
          JZN= NZ(JSA)
          IF (JZN.EQ.0) GOTO 2010
C
C         LOOP OVER ZEROS, JZA= STARTING ZERO
          DO 2000 JZA= 1,JZN
C
C           FOR THIRD JOURNEY, START ALSO AT ISTATZ=0
            IF (JOURNY.EQ.3 .AND. ISTATZ(JZA,JSA).EQ.0) GOTO 1320
C
            IF (ISTATZ(JZA,JSA).NE.ISTART)   GOTO 2000
C           IF ()  THEN THIS STATION WILL NOT SERVE AS START
 1320       CONTINUE
C
C
C           START TRIP FROM STATION  SIDE JSA, ZERO JZA
C
C           IF JOURNY.NE.1 CHECK STACK FIRST
            IF (JOURNY.EQ.1 .OR.
     A          JSA.NE.JSE1ST .OR. JZA.NE.JZE1ST) GOTO 1340
C           COORDINATES FOR NEXT RIDE ARE IN STACK
            JJ= NP1ST
            DO 1330 J=1,NP1ST
              XPOL(J)= XSTACK(JJ)
              YPOL(J)= YSTACK(JJ)
              JJ= JJ-1
 1330       CONTINUE
            NPOL1= NP1ST
            JSA2 = JSA1ST
            JZA2 = JZA1ST
            NSTOP= 0
C           SET NORMAL DERIV. AND COLOR
            DERNOR= -DERNOS
            DERNO1= DERNOR
            NCL1= NCLZR(JZA,JSA)
            NCOL= AMIN1(-DERNOR*NDIR3+1.,1.) + NCL1
            GOTO 1350
C
 1340       CONTINUE
C
C           SET NORMAL DERIVATIVE AND COLOUR FOR FIRST RIDE
            IF (SDER(JZA,JSA).EQ.0.) GOTO 2000
            DERNOR= SIGN(1.,SDER(JZA,JSA))
            DERNO1= DERNOR
            NCL1= NCLZR(JZA,JSA)
            NCOL= AMIN1(-DERNOR*NDIR3+1.,1.) +
     1            NCL1
C           NCOL= COLOUR FOR TRIP
C           THE FIRST LINE OF THE LAST STATEMENT IS 1 OR 0,
C           DEPENDING ON THE SIGN OF SDER
C           (DERIVATIVE IN DIRECTION OF SIDE)
C           AND NDIR3
C           NDIR3= 1, FOR JOURNY=1,2 (COUNTER-CLOCKWISE TRIP)
C           NDIR3= -1, FOR JOURNY=3  (CLOCKWISE TRIP)
C
          CALL FRBRID(JSA,JZA,CC,NC,XPOL,YPOL,MAXPOL,JSA2,JZA2,NPOL1,
     A                NSTOP)
C           FIRST RIDE ENDED ON STATION SIDE JSA2, ZERO JZA2
C
C           CALL FOR LINE DRAWING
            IF (MODE.NE.1 ) GOTO 1350
                CALL USRPLT(XPOL,YPOL,NPOL1,NCL1,1)
C
C           BOOK KEEPING FOR START AND END OF RIDE
 1350       CONTINUE
            JPOL= NPOL1
            JPOLL1= 1
            NP= NPOL1
            IRIDE= 1
            JSER(1)= JSA2
            JZER(1)= JZA2
            JSAR(1)= JSA
            JZAR(1)= JZA
            JSTOPR(1)= NSTOP
            JSIDES= 0
            IF (NSTOP.EQ.2) GOTO 2000
C
C           FIND NEXT STATION FOR CONTINUATION OF TRIP (=TRANSFER)
C
 1360       CONTINUE
            JZA2 = JZA2 + NDIR3
            IF (JZA2.LE.NZ(JSA2).AND.JZA2.GT.0) GOTO 1400
C
C           TAKE NEXT SIDE, ADD VERTEX TO POLYGON
C           SET JZA2 TO FIRST OR LAST ZERO OF THE NEW SIDE
 1370       CONTINUE
            JSIDES= JSIDES + 1
            IF (JSIDES.EQ.5) GOTO 2000
            JPOL= JPOL + 1
            KVERT= MOD(JSA2+NDIRV,4) +1
            XPOL(JPOL)= X(KVERT)
            YPOL(JPOL)= Y(KVERT)
            JSA2= MOD(JSA2+NDIRS,4) + 1
            IF (NZ(JSA2) .EQ. 0) GOTO 1370
            JZA2= 1
            IF (NDIR3.EQ.-1) JZA2= NZ(JSA2)
C
C           CHECK FOR REGULAR END OF TRIP
 1400       CONTINUE
            IF (JSA2.EQ.JSA .AND. JZA2.EQ.JZA) GOTO 1900
C
C           CHECK DIFFERENCE OF CONTOUR LEVELS
C           BETWEEN JSA,JZA AND JSA2,JZA2,
C           CHECK SIGN OF DERIVATIVE OF STATION JZA2,JSA2,
C           SET NEW DERNOR
C
            NCDIF= NCLZR(JZA2,JSA2) - NCL1
            NCDIFA= ABS(NCDIF)
            NC1= -NCDIFA
            IF (NC1.EQ.0) NC1= 1
            DERNOR= DERNO1*NC1
            SDCHEK= SDER(JZA2,JSA2)*DERNOR
            IF (SDCHEK.GE.0. .AND. NCDIFA.LE.1) GOTO 1420
C
C           DO NOT STOP AT VERTEX
            IF (TZR(JZA2,JSA2).LT.POSERR*3. .OR.
     2          TZR(JZA2,JSA2).GT.SL(JSA2)-POSERR*3.) GOTO 1360
C
C           STOP TRIP IN ALL OTHER CASES
            GOTO 2000
C
 1420       CONTINUE
C
C           START NEW RIDE FROM SIDE JSA2, ZERO JZA2
C
C           CHECK STACK FIRST
            IF (JSA2.NE.JSE1ST .OR. JZA2.NE.JZE1ST) GOTO 1500
C           COORDINATES FOR NEXT RIDE ARE IN STACK
C           THE FOLLOWING REPLACES A CALL TO FRBRID
            IF (JPOL+NP1ST.GT.MAXPOL) GOTO 7000
            JJ= NP1ST
            DO 1450 J=1,NP1ST
              XPOL(JPOL+J)= XSTACK(JJ)
              YPOL(JPOL+J)= YSTACK(JJ)
              JJ= JJ-1
 1450       CONTINUE
            JSE= JSA1ST
            JZE= JZA1ST
            NP= NP1ST
            NSTOP= 0
            GOTO 1600
C
 1500       CONTINUE
          CALL FRBRID (JSA2,JZA2,CC,NC,
     A      XPOL(JPOL+1),YPOL(JPOL+1),MAXPOL-JPOL,JSE,JZE,NP,NSTOP)
C           RIDE ENDED AT STATION  SIDE JSE, ZERO JZE
C
C           CALL FOR LINE DRAWING
            IF (MODE.NE.1 .OR. ISTATZ(JZA2,JSA2).NE.0 .OR.
     A          NP.LT.2)                         GOTO 1595
                CALL USRPLT(XPOL(JPOL+1),YPOL(JPOL+1),
     A                      NP,NCLZR(JZA2,JSA2),1)
 1595       CONTINUE
C
C           BOOK KEEPING FOR CONTINUATION RIDE
 1600       CONTINUE
            IF (NSTOP.EQ.2) NP= 0
            JPOLL1= JPOL + 1
            JPOL= JPOL + NP
            IRIDE= IRIDE + 1
            JSAR(IRIDE)= JSA2
            JZAR(IRIDE)= JZA2
            JSER(IRIDE)= JSE
            JZER(IRIDE)= JZE
            JSTOPR(IRIDE)= NSTOP
            JSA2= JSE
            JZA2= JZE
            IF (IRIDE+1.GT.MAXRID) GOTO 2000
C           CONTINUE TRIP
            GOTO 1360
C
C           POLYGON FOR FILL AREA IS COMPLETE
C           (SUCCESSFULL ROUND TRIP)
C
 1900       CONTINUE
C
C           WRITE COORDINATES OF FIRST RIDE TO STACK
            IF (NPOL1.GT.MAXSTA) GOTO 1960
            IF (JSTOPR(1).EQ.1) GOTO 1960
            DO 1950 J=1,NPOL1
              XSTACK(J)= XPOL(J)
              YSTACK(J)= YPOL(J)
 1950       CONTINUE
            JSA1ST= JSA
            JZA1ST= JZA
            JSE1ST= JSER(1)
            JZE1ST= JZER(1)
            DERNOS= DERNO1
            NP1ST= NPOL1
C
C           CALL FOR FILL AREA
 1960       CONTINUE
            IF (JPOL.LT.3) GOTO 2000
            NFAR= NFAR + 1
            IF (MOD(MODE,2).EQ.0)
     A      CALL USRPLT(XPOL,YPOL,JPOL,ICOL(NCOL),0)
C
C           SET FLAGS FOR START AND END OF RIDES
            DO 1990 JR= 1,IRIDE
              JS= JSAR(JR)
              JZ= JZAR(JR)
              IF (MOD(ISTATZ(JZ,JS),2).EQ.0)
     A           ISTATZ(JZ,JS)= ISTATZ(JZ,JS) + 1
              IF (JSTOPR(JR).GT.0) GOTO 1990
              JS= JSER(JR)
              JZ= JZER(JR)
              IF (ISTATZ(JZ,JS).LT.2)
     A          ISTATZ(JZ,JS)= ISTATZ(JZ,JS) + 2
 1990       CONTINUE
C
C         DRAW LINE IF MODE=2
          JS= JSAR(IRIDE)
          JZ= JZAR(IRIDE)
          IF (MODE.EQ.2 .AND. NP.GT.1)
     A       CALL USRPLT(XPOL(JPOLL1),YPOL(JPOLL1),NP,NCLZR(JZ,JS),1)
C
 2000     CONTINUE
C         END OF LOOP OVER ZEROS
 2010     CONTINUE
 3000   CONTINUE
C       END OF LOOP OVER SIDES
C
C       SET ISTART FOR JOURNY=2
        ISTART= 2
C
 4000 CONTINUE
      NPP= NPP + NPREC
C
C     TREAT CASE WHEN WHOLE RECTANGLE HAS TO BE FILLED
C     BECAUSE THERE WAS NO SUCCESSFULL TRIP
      IF (NFAR.NE.0 .OR. MODE.EQ.1) GOTO 4050
      NCOL= ICN + 1
      IF (CC(ICN).EQ.ZMAXT) NCOL= ICN
      GOTO 250
 4050 CONTINUE
C
 5000 CONTINUE
      RETURN
C
C     ERROR EXIT
 6000 CONTINUE
      WRITE (*,6001) NCMAXS
 6001 FORMAT (' ***ERROR IN FARBRC'/
     A        ' MORE THAN', I8, ' CONTOURS CROSSING A',
     B        ' SIDE OF A RECTANGLE.'/
     C        ' INCREASE INSTALLATION PARAMETER NCMAXS')
      RETURN
C
 7000 CONTINUE
      WRITE (*,7001) MAXPOL,IT
 7001 FORMAT (' *** ERROR *** IN FARBRC'/
     A        ' OVERFLOW OF WORKING STORAGE XPOL,YPOL'/
     B        ' MAXPOL= ',I7/' RECT.NO. ',I10)
      RETURN
      END
      SUBROUTINE FRBRID (JSA,JZA,CN,NC,XPOL,YPOL,MAXPOL,JSA2,JZA2,NP,
     A                   NSTOP)
C
C     TRACE CONTOUR FROM SIDE JSA TO SIDE JSA2
C     (RIDE FROM JSA,JZA TO JSA2,JZA2)
C
C     T R I P   ALGORITHM   A.PREUSSER   FARB-E-2D  VERSION 2.1 10/1988
C
C
C     AUTHOR: A. PREUSSER
C             FRITZ-HABER-INSTITUT DER MPG
C             FARADAYWEG 4-6
C             D-1000 BERLIN 33
C
C
C     INPUT PARAMETERS
C     JSA         SIDE INDEX  CONTOUR STARTS FROM
C     JZA         ZERO INDEX  CONTOUR STARTS FROM
C     CN          CONTOUR LEVELS
C     NC          NUMBER OF CONTOUR LEVELS
C     MAXPOL      MAXIMUM NUMBER OF POINTS IN XPOL,YPOL
C
C     OUTPUT PARAMETERS
C     XPOL,YPOL    X-Y-COORDINATES OF THE POINTS OF A RIDE
C     JSA2         SIDE INDEX WHERE CONTOUR ENDS
C     JZA2         ZERO (STATION) INDEX WHERE CONTOUR ENDS
C     NP           NUMBER OF POINTS STORED TO XPOL,YPOL
C     NSTOP        =0, RIDE ENDED AT STATION
C                  =1, RIDE ENDED ON SIDE, NO STATION FOUND
C                  =2, RIDE ENDED INSIDE RECTANGLE
C
C                  IF NSTOP.EQ.1, JSA2,JZA2 INDICATE THE PREVIOUS
C                  STATION ON THE ENDING SIDE.
C                  (NOTE POSITIVE OR NEGATIVE SENSE OF TRIP
C                  INDICATED BY NDIR3). JZA2 MAY BE ZERO.
C
C                  IF NSTOP.EQ.2, JSA2= JSA;    JZA2= JZA  .
C
C
C
      DIMENSION CN(NC),XPOL(MAXPOL),YPOL(MAXPOL)
C
      COMMON /FRBCOC/ SACMIN,CMSCAL,NPMAX,NCPMAX,MAXSTA,NCMAXS,NPP,PI
     1,               MAXRID, SIGS(4), NCMAX
C
      COMMON /FRBCOF/ KK,KSE,XX4F,YY4F,SIR,COR,CL
C
      COMMON /FRBCRD/ X0(50,4),Y0(50,4),NCLZR(50,4),SDER(50,4),TZR(50,4)
     A,               NZ(4),SI(4),CO(4),SA(4),SE(4),DX(4),DY(4),SL(4)
     B,               HMIN,SLMAX,KRIDE,NPREC,POSERR,DERNOR,NDIR3
     C,               ZMAX(4),ZMIN(4),ISTATZ(50,4),X3,Y3,ZSOLD
C
      COMMON /FRBCOR/ RMA,RMAX,DSMAX,DSMIN,FSTEP
     A,               THETAS(4),RACMIN
      SAVE /FRBCOR/
C
      KSE= JSA
      KRIDE= KRIDE +1
      NSAD= 0
      NSTOP= 0
C
C
C     INITIALIZATION, IF FIRST RIDE OF A RECTANGLE
C
      IF (KRIDE.NE.1) GOTO 1400
      RMAX  = AMIN1(0.010/CMSCAL,HMIN*0.010)
C           = DISTANCE NORMAL TO CURVE DIRECTION WITHIN WHICH
C             A ZERO MUST BE FOUND
      DSMAX = HMIN*0.2
C           = MAXIMUM STEP SIZE
      DSMIN = AMIN1(RMAX*0.03,POSERR*8.)
C           = MINIMUM STEP SIZE
      FSTEP = AMIN1(RMAX*8.,DSMAX)
C           = STARTING STEP SIZE
      RACMIN= SACMIN*0.1
C             A POINT IS STORED ONLY IF THE ACCUMULATED R'S (RACC)
C             (CHANGE IN DIRECTION) HAVE REACHED RACMIN
C     SET DIRECTION OF SIDES
      THETAS(1)= 0.
      THETAS(2)= PI*0.5
      THETAS(3)= PI
      THETAS(4)= PI*1.5
 1400 CONTINUE
C
C
C     DEFINE CONTOUR LEVEL
      KCL= NCLZR(JZA,KSE)
      CL= CN(KCL)
C
C     ESTIMATE STARTING DIRECTION
C
C     COMPUTE F1 ON SIDE
      NSE= MOD(KSE,4) + 1
      EPS= 0.01*SL(NSE)
      F1= SDER(JZA,KSE)*EPS
C
C     COMPUTE F2 NORMAL TO SIDE
      KK= 2
      XX4F= X0(JZA,KSE)
      YY4F= Y0(JZA,KSE)
      SIR= CO(KSE)
      COR= - SI(KSE)
      F2= FRBEVA(EPS)
C
C     COMPUTE ANGLE FOR STARTING DIRECTION
      IF (F2.EQ.0.) GOTO 1470
      THETSC= ATAN (-F1/F2)
      IF (THETSC.LT.0.) THETSC= THETSC + PI
      IF (THETSC.EQ.0.) THETSC= (NDIR3-1)*PI*0.5
      GOTO 1480
 1470 CONTINUE
C     IF (F1.EQ.0.) GOTO 1690
      THETSC= PI*0.5
 1480 THETAC= THETAS(KSE) + THETSC
C
C     COMPUTE POINTS
C
C     STORE FIRST POINT
      JP= 1
C     JP= NUMBER OF POINTS STORED FOR THIS CONTOUR LINE
      XPOL(JP)= X0(JZA,KSE)
      YPOL(JP)= Y0(JZA,KSE)
C
C     INITIALIZE TRACING
 1490 CONTINUE
      R= 0.
      DS= FSTEP
      DX12= DS*COS(THETAC)
      DY12= DS*SIN(THETAC)
      XX3= XX4F
      YY3= YY4F
      XX2= XX3 - DX12
      YY2= YY3 - DY12
      XX1= XX2 - DX12
      YY1= YY2 - DY12
      DS23= DS
      DS12= DS
C     DS01= DS
C     POINTS P1,P2,P3,P4 WITH COORDINATES XX1...XX4, YY1...YY4
C     ARE REFERRED TO AS *QUEUE*. DS12...DS34 ARE THE DISTANCES
C     BETWEEN POINTS IN THE QUEUE.
C     XX4F,YY4F ARE PRELEMINARY COORDINATES FOR THE NEXT POINT P4
C     WHICH WILL BE COMPUTED BY THE REGULA FALSI (FRBZER).
C     FOR A DERIVATION OF THE FORMULAS FOR PL0...PL2 SEE
C     PREUSSER,A. COMPUTING AREA FILLING CONTOURS FOR SURFACES
C                 DEFINED BY PIECEWISE POLYNOMIALS.
C                 COMPUTER AIDED GEOMETRIC DESIGN 3,
C                 (1986), P. 267-279
C     THERE IS ALSO AN EXPLANATION FOR THE FOLLOWING PART OF
C     THE ALGORITHM.
      SACC= 0.
C         = ACCUMULATED DISTANCES TO LAST POINT STORED
      RACC= RACMIN
C         = ACCUMULATED R
C     A POINT IS ONLY STORED TO XPOL,YPOL IF SACC.GE.SACMIN
C                                  AND       RACC.GE.RACMIN
      NCP= 0
C         = NUMBER OF POINTS COMPUTED
      NOST= 0
C         = NUMBER OF STEPS NORMAL TO CURVE
      RMA= RMAX
C
C     COMPUTE NEW POINT FOR CONTOUR LINE
C
C     COMPUTE CURVE DIRECTION
 1500 CONTINUE
      DS13= DS23 + DS12
      PL0=  DS23/(DS12*DS13)
      PL1= -DS13/(DS12*DS23)
      PL2= (DS13+DS23)/(DS13*DS23)
      DXDS= PL0*XX1 + PL1*XX2 + PL2*XX3
      DYDS= PL0*YY1 + PL1*YY2 + PL2*YY3
      SQ= SQRT(DXDS*DXDS+DYDS*DYDS)
      DXDS= DXDS/SQ
      DYDS= DYDS/SQ
      COR= -DYDS
      SIR=  DXDS
C
C     SEARCH FOR TWO POINTS WITH OPPOSITE SIGN
      RMA= SIGN (RMAX,RMA)
 1550 CONTINUE
      XX4F= XX3 + DXDS*DS
      YY4F= YY3 + DYDS*DS
      F1= FRBEVA(0.)
      RMA= SIGN(RMA,DERNOR*F1)
      F2= FRBEVA(RMA)
      IF (F1*F2 .LE. 0.) GOTO 1600
C
 1560 CONTINUE
      IF (DS*0.5 .LT. DSMIN) GOTO 1570
C     DIVIDE STEPSIZE IN CURVE DIRECTION BY 2.
      DS= DS*0.5
      GOTO 1550
C
C     DIVIDE STEPSIZE NORMAL TO CURVE BY 2.
 1570 CONTINUE
      NOST= NOST + 1
      RMA= RMA*0.5
      IF (ABS(RMA).LE.POSERR) GOTO 1580
      F2= FRBEVA(RMA)
      IF (F1*F2.GT.0.) GOTO 1570
      GOTO 1600
C
C     SADDLE POINT
C
C     SET NEW DIRECTION
 1580 NSAD= NSAD + 1
      IF (NSAD.GT.1) GOTO 1690
      DXDS1= DXDS
      DYDS1= DYDS
      THETAC= ATAN2(DYDS,DXDS) + PI*0.5
C
C     STORE SADDLE POINT
      XX4F= XX3
      YY4F= YY3
      JPP= JP +1
      IF (JPP.GT.MAXPOL) GOTO 4000
      JP= JPP
      XPOL(JP)= XX3
      YPOL(JP)= YY3
      GOTO 1490
C
C
C     FIND ZERO FOR NEW POINT
C
 1600 CONTINUE
      R= FRBZER(0.,RMA,F1,F2,POSERR)
      NCP= NCP + 1
      IF (NCP.GT.NCPMAX) GOTO 1690
      DS34= SQRT(DS*DS + R*R)
      XX4= XX4F + COR*R
      YY4= YY4F + SIR*R
C
C     CHECK IF POINT IS OUTSIDE THE RECTANGLE
      JSA2=1
      IF (YY4.LT.0.     ) GOTO 1700
      JSA2=2
      IF (XX4.GT.DX(1)       ) GOTO 1700
      JSA2= 3
      IF (YY4.GT.DY(2)       ) GOTO 1700
      JSA2= 4
      IF (XX4.LT.0.     ) GOTO 1700
C
C     POINT IS INSIDE
C
C     STORE POINT TO XPOL,YPOL
      SACC= SACC + DS34
      RACC= RACC + ABS(R)
      IF (SACC.LT.SACMIN .OR. RACC.LT.RACMIN) GOTO 1650
      JPP= JP + 1
      IF (JPP.GT.MAXPOL) GOTO 4000
      JP= JPP
      XPOL(JP)= XX4
      YPOL(JP)= YY4
      SACC= 0.
      RACC= 0.
C
C     UPDATE QUEUE
 1650 CONTINUE
C     DS01= DS12
      DS12= DS23
      DS23= DS34
      XX1= XX2
      YY1= YY2
      XX2= XX3
      YY2= YY3
      XX3= XX4
      YY3= YY4
C
C     SET NEW STEP SIZE
      SOLL1= 2.
      IF (ABS(R).GT.POSERR) SOLL1= ABS(RMA*0.8/R)
      IF (SOLL1.GT.1) DS= AMIN1(DSMAX,DS*SQRT(SOLL1))
C
      GOTO 1500
C
C     TRACING STOPPED
 1690 NSTOP= 1
      JSA2= 1
      DXDS= DXDS1
      DYDS= DYDS1
C
C     POINT IS OUTSIDE
C
C     SEARCH FOR CORRESPONDING ZERO ON SIDES
C     START WITH SIDE JSA2
 1700 CONTINUE
      R30MIN= 99999.
      S30MAX= DSMAX
      JSE= JSA2
      DO 1780 N=1,2
        DO 1770 JESE=1,4
          JJ= NZ(JSE)
          IF (JJ.EQ.0) GOTO 1760
          DO 1750 J=1,JJ
C
C           LEVEL CHECK
            IF (NCLZR(J,JSE).NE.NCLZR(JZA,KSE)) GOTO 1750
C
C           DERIVATIVE CHECK
            IF (SDER(J,JSE)*DERNOR.GT.0.) GOTO 1750
C
C           R-CHECK
            DX30= X0(J,JSE) - XX3
            DY30= Y0(J,JSE) - YY3
            R30= (DY30*DXDS - DX30*DYDS)
            IF (R30*R.LT.0 .AND. N.EQ.1) GOTO 1750
            R30= ABS(R30)
            IF (R30.GE.R30MIN) GOTO 1750
C
C           S-CHECK
            DS30= DX30*DXDS + DY30*DYDS
            IF (DS30.GT.S30MAX .OR. DS30.LT.-DSMIN*2.) GOTO 1750
C
            JZA2= J
            JSA2= JSE
            R30MIN= R30
 1750     CONTINUE
 1760     CONTINUE
          JSE= MOD(JSE,4) + 1
 1770   CONTINUE
        IF (R30MIN.LT.RMAX) GOTO 1800
C       FOR N=2, DO NOT CHECK SIGN OF R30
 1780 CONTINUE
C
C     NO ACCEPTABLE ZERO ON ALL THREE SIDES
C
C     REDUCE STEP SIZE
      IF (NSTOP.NE.1) GOTO 1560
C
      GOTO 1850
C
C     STORE END STATION OF THE RIDE
 1800 CONTINUE
      NSTOP= 0
 1850 CONTINUE
      IF (JP+1.GT.MAXPOL) GOTO 4000
      IF (NSTOP.EQ.1) GOTO 3000
      JP= JP+1
      XPOL(JP)= X0(JZA2,JSA2)
      YPOL(JP)= Y0(JZA2,JSA2)
C
C     ADD COORDINATES OF LOWER LEFT VERTEX
C     NORMAL RETURN
 2000 CONTINUE
      NP= JP
      DO 2050 JP=1,NP
        XPOL(JP)= XPOL(JP) + X3
        YPOL(JP)= YPOL(JP) + Y3
 2050 CONTINUE
      NPREC= NPREC + NP
      RETURN
C
C     ERROR HANDLING
 3000 CONTINUE
      IF (JP.LE.2) GOTO 3900
C
C     CHECK IF LAST POINT IS NEAR BOUNDARY
      JSA2= 1
      TTZR= XX3
      IF (ABS(YY3).LT.SACMIN) GOTO 3300
C
      JSA2= 2
      TTZR= YY3
      IF (ABS(XX3-DX(1)).LT.SACMIN) GOTO 3300
C
      JSA2= 3
      TTZR= XX3
      IF (ABS(YY3-DY(2)).LT.SACMIN) GOTO 3300
C
      JSA2= 4
      TTZR= YY3
      IF (ABS(XX3).LT.SACMIN) GOTO 3300
      GOTO 3900
C
C     TRACING STOPPED NEAR BOUNDARY
C     NSTOP= 1, JSA2 SIDE NUMBER,
C     JZA2 ZERO INDEX  OF LAST OR NEXT ZERO, OR =0
C     DEPENDING ON NDIR3.
C
C     FIND JZA2
 3300 UGR= SA(JSA2)
      JZA2= 0
      JN= NZ(JSA2)
      IF (JN.EQ.0) GOTO 3500
      DO 3400 JJ=1,JN
         JZA2= JZA2 + 1
         OGR= TZR(JZA2,JSA2)
         IF (TTZR.GE.UGR .AND. TTZR.LE.OGR .OR.
     A       TTZR.GE.OGR .AND. TTZR.LE.UGR) GOTO 3450
         UGR= OGR
 3400 CONTINUE
      JZA2= JZA2 + 1
 3450 CONTINUE
      IF (NDIR3.EQ.1) JZA2= JZA2-1
 3500 CONTINUE
      GOTO 2000
C
C     TRACING STOPPED WITHIN RECTANGLE
 3900 CONTINUE
      NSTOP= 2
      JSA2= KSE
      JZA2= JZA
      GOTO 2000
C
 4000 CONTINUE
      WRITE (*,9991)
 9991 FORMAT (1X,'***WARNING*** IN FRBRID'/
     A        ' OVERFLOW OF WORKING STORAGE XPOL,YPOL')
      GOTO 3900
      END
      SUBROUTINE USRPLT (X,Y,N,NCOL,MODE)
      DIMENSION X(N),Y(N)
C
C     ROUTINE FOR CALLS TO PLOT SYSTEM
C
C     MODE= 1, LINE DRAWING
C           0, FILL AREA
C
C     X,Y      COORDINATES FOR POLYGON
C     N        NUMBER OF POINTS FOR POLYGON
C     NCOL     COLOUR FOR POLYGON (INDEX OF CONTOUR LEVEL)
C
      IF (MODE.EQ.0) GOTO 1000
C
C     DRAW LINE
        call newpen(1)
      CALL GPL (N,X,Y)
      RETURN
C
C     FILL AREA
 1000 CONTINUE
C     SET FILL AREA INDEX
c      CALL GSFAI(NCOL)
        call newpen(ncol)
C     FILL AREA
      CALL GFA (N,X,Y)
C
      RETURN
      END
        subroutine gfa(n,x,y)
        real x(n), y(n)
        call shadep(n,x,y)
        return
        end
        subroutine gpl(n,x,y)
        real x(n), y(n)
        call plot(x(1),y(1),3)
        do 1000 i=2,n
            call plot(x(i),y(i),2)
 1000   continue
        call plot(x(1),y(1),3)
        return
        end

      FUNCTION FRBEVA(T)
C
C     T R I P   ALGORITHM   A.PREUSSER   FARB-E-2D  VERSION 2.1 10/1988
C
C     FUNCTION EVALUATION
C
C     AUTHOR      : A. PREUSSER
C                   FRITZ-HABER-INSTITUT
C                   DER MAX-PLANCK-GESELLSCHAFT
C                   FARADAYWEG 4-6
C                   D-1000 BERLIN 33
C
      COMMON /FRBCOF/ KK,KSE,XX4F,YY4F,SIR,COR,CL
C
C     VARIABELS IN /FRBCOF/ ARE USED AS ARGUMENTS
C     FOR AN EXPLANATION SEE SUBROUTINE FARBRC
C
C     KK      NUMBER OF FUNCTION TO BE EVALUATED
C     KK=1    ORIGINAL POLYNOMIAL ALONG SIDE KSE
C     KK=2    BIVARIATE POLYNOMIAL INSIDE RECTANGLE
C     KK=3    2ND DERIVATIVE ALONG SIDE KSE
C     KK=4    1ST DERIVATIVE ALONG SIDE KSE
C
      COMMON /FRBCOP/  P0(4),P1(4),P2(4),P3(4)
     1,     Q0(4),Q1(4),Q2(4)
     2,     R0(4),R1(4)
     3,     P11,P12,P13,P21,P22,P23,P31,P32,P33
C
C
      IF (KK.EQ.2) GOTO 20
      IF (KK.EQ.4) GOTO 40
      IF (KK.EQ.3) GOTO 30
      FRBEVA= P0(KSE)+T*(P1(KSE)+T*(P2(KSE)+T*P3(KSE))) - CL
      RETURN
   20 XX4= XX4F + COR*T
      YY4= YY4F + SIR*T
      S0= P0(1) + YY4*(P1(4)+YY4*(P2(4)+YY4*P3(4)))
      S1= P1(1) + YY4*(P11+YY4*(P12+YY4*P13))
      S2= P2(1) + YY4*(P21+YY4*(P22+YY4*P23))
      S3= P3(1) + YY4*(P31+YY4*(P32+YY4*P33))
      FRBEVA= S0 + XX4*(S1+XX4*(S2+XX4*S3)) - CL
      RETURN
   30 FRBEVA= R0(KSE) + T*R1(KSE)
      RETURN
   40 FRBEVA= Q0(KSE) + T*(Q1(KSE)+T*Q2(KSE))
      RETURN
      END
      SUBROUTINE FRBFOP(X,Y,ICOL,NCOL)
      DIMENSION X(4),Y(4),ICOL(*)
C
C     OPEN FILL AREA BUFFER
C
      COMMON /FRBCOB/ NFABU,NCOLBU,XFABU(4),YFABU(4)
C     IF ALREADY OPEN, CLEAR BUFFER FIRST
      IF (NFABU.EQ.1) CALL FRBFCL (ICOL)
C     FILL BUFFER
      DO 100 J=1,4
        XFABU(J)= X(J)
        YFABU(J)= Y(J)
  100 CONTINUE
C     SET COLOUR OF FILL AREA BUFFER
      NCOLBU= NCOL
C     DECLARE FILL AREA BUFFER OPEN
      NFABU= 1
      RETURN
      END
      SUBROUTINE FRBFUP(X,Y)
      DIMENSION X(4),Y(4)
C
C     UPDATE FILL AREA BUFFER
C
      COMMON /FRBCOB/ NFABU,NCOLBU,XFABU(4),YFABU(4)
      XFABU(4)= X(4)
      YFABU(4)= Y(4)
      XFABU(1)= X(1)
      YFABU(1)= Y(1)
      RETURN
      END
      SUBROUTINE FRBFCL(ICOL)
      DIMENSION ICOL(*)
C
C     CLEAR FILL AREA BUFFER
C
      COMMON /FRBCOB/ NFABU,NCOLBU,XFABU(4),YFABU(4)
C     RETURN, IF FILL AREA BUFFER IS CLOSED
      IF (NFABU.NE.1) GOTO 100
C     DECLARE FILL AREA BUFFER CLOSED
      NFABU= 0
C     CALL FILL AREA ROUTINE
      CALL USRPLT(XFABU,YFABU,4,ICOL(NCOLBU),0)
  100 CONTINUE
      RETURN
      END
      FUNCTION FRBZER (TA,TB,F1,F2,ER)
C
C     T R I P   ALGORITHM   A.PREUSSER   FARB-E-2D  VERSION 2.1 10/1988
C
C     COMPUTE ZERO BETWEEN TA AND TB
C
C     F1= FUNCTION VALUE AT TA
C     F2= FUNCTION VALUE AT TB
C         F1 AND F2 MUST HAVE OPPOSITE SIGN
C         THIS MUST BE CHECKED BEFORE ENTRY
C     ER= PERMITTED ERROR FOR SOLUTION FRBZER
C     NAME OF FUNCTION = FRBEVA
C
C     THE METHOD IS A COMBINATION OF THE REGULA FALSI
C     AND THE MIDPOINT METHOD
C
C     IT IS A MODIFIED VERSION OF THE VIM- (CONTROL DATA
C     USER GROUP) ROUTINE WITH CATALOG IDENTIFICATION
C                C2BKYZERO
C     WRITTEN BY LOREN P. MEISSNER, 1965
C
      A=TA
      B=TB
      FA=F1
      FB=F2
      C=A
      FC=FA
      S=C
      FS=FC
C
   10 CONTINUE
      H=0.5*(B+C)
      IF(ABS(H-B) .LE.ER) GO TO 110
      IF (ABS(FB) .LE. ABS(FC)) GO TO 15
      Y=B
      FY=FB
      G=B
      FG=FB
      S=C
      FS=FC
      GO TO 20
   15 Y=S
      FY=FS
      G=C
      FG=FC
      S=B
      FS=FB
   20 CONTINUE
      IF (FY .NE. FS) GO TO 21
      B=H
      GO TO 29
   21 CONTINUE
      E=(S*FY-Y*FS)/(FY-FS)
      IF (ABS(E-S) .LE.ER) E=S+SIGN(ER,G-S)
      IF ((E-H)*(S-E) .LT. 0.0) GO TO 28
      B=E
      GO TO 29
   28 B=H
C
C *** FUNCTION CALL
   29 FB=FRBEVA(B)
C
      IF (FG*FB .LT. 0.0) GO TO 35
      C=S
      FC=FS
      GO TO 10
   35 CONTINUE
      C=G
      FC=FG
      GO TO 10
C
  110 FRBZER= H
C
      RETURN
      END
      SUBROUTINE GTEXT (X,Y,H,TEXT,W,N)
      CHARACTER*(*) TEXT
C
C     TEXT PLOTTING
C     PARAMETERS ARE IDENTICAL WITH CALCOMP-ROUTINE
C                  SYMBOL
C
C     X,Y        COORDINATES OF LOWER LEFT CORNER OF TEXT
C     H          TEXT HEIGHT
C     TEXT       TEXT
C     W          ANGLE FOR TEXT
C     N          NUMBER OF CHARACTERS IN TEXT
C
        call symbol(x,y,h,text,w,n)
        return
        end
