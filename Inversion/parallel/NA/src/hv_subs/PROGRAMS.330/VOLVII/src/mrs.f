c-----
chttp://gams.nist.gov/serve.cgi/ModuleComponent/1070/Fullsource/ITL/RS.f
c http://gams.nist.gov/serve.cgi/Module/EISPACK/SVD/7956/
c
c     Modified RBHerrmann, Saint Louis University 07 JAN 2003 to
c     use DOUBLE PRECISION INSTEAD OF REAL. THIS ALSO
c     REQUIRED AMIN1 -> DMIN1 AMAX1 -> DMAX1
c-----
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C***BEGIN PROLOGUE  TQL2
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4A5,D4C2A
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Compute eigenvalues and eigenvectors of symmetric
C            tridiagonal matrix.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TQL2,
C     Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     This subroutine finds the eigenvalues and eigenvectors
C     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
C     The eigenvectors of a FULL SYMMETRIC matrix can also
C     be found if  TRED2  has been used to reduce this
C     full matrix to tridiagonal form.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E contains the subdiagonal elements of the input matrix
C          in its last N-1 positions.  E(1) is arbitrary.
C
C        Z contains the transformation matrix produced in the
C          reduction by  TRED2, if performed.  If the eigenvectors
C          of the tridiagonal matrix are desired, Z must contain
C          the identity matrix.
C
C      On Output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct but
C          unordered for indices 1,2,...,IERR-1.
C
C        E has been destroyed.
C
C        Z contains orthonormal eigenvectors of the symmetric
C          tridiagonal (or full) matrix.  If an error exit is made,
C          Z contains the eigenvectors associated with the stored
C          eigenvalues.
C
C        IERR is set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  PYTHAG
C***END PROLOGUE  TQL2
C
        implicit none
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
      DOUBLE PRECISION B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
      DOUBLE PRECISION PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  TQL2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0E0
      B = 0.0E0
      E(N) = 0.0E0
C
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B .LT. H) B = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            IF (B + ABS(E(M)) .EQ. B) GO TO 120
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0E0 * E(L))
         R = PYTHAG(P,1.0D+00)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0E0
         C2 = C
         EL1 = E(L1)
         S = 0.0E0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0E0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0E0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0E0)
            E(I+1) = S * E(I) * R
            S = 1.0E0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
C***BEGIN PROLOGUE  TRED2
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C1B1
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Reduce real symmetric matrix to symmetric tridiagonal
C            matrix using and accumulating orthogonal transformation
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TRED2,
C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     This subroutine reduces a REAL SYMMETRIC matrix to a
C     symmetric tridiagonal matrix using and accumulating
C     orthogonal similarity transformations.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        A contains the real symmetric input matrix.  Only the
C          lower triangle of the matrix need be supplied.
C
C     On Output
C
C        D contains the diagonal elements of the tridiagonal matrix.
C
C        E contains the subdiagonal elements of the tridiagonal
C          matrix in its last N-1 positions.  E(1) is set to zero.
C
C        Z contains the orthogonal transformation matrix
C          produced in the reduction.
C
C        A and Z may coincide.  If distinct, A is unaltered.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  TRED2
C
        implicit none
      INTEGER I,J,K,L,N,II,NM,JP1
      DOUBLE PRECISION A(NM,N),D(N),E(N),Z(NM,N)
      DOUBLE PRECISION F,G,H,HH,SCALE
C
C***FIRST EXECUTABLE STATEMENT  TRED2
      DO 100 I = 1, N
C
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE
C
      IF (N .EQ. 1) GO TO 320
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 2) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I,K))
C
         IF (SCALE .NE. 0.0E0) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
C
  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
C
         F = Z(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0E0
C
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / H
            G = 0.0E0
C     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
C     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE
C
         HH = F / (H + H)
C     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
C
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
C
  290    D(I) = H
  300 CONTINUE
C
  320 D(1) = 0.0E0
      E(1) = 0.0E0
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0E0) GO TO 380
C
         DO 360 J = 1, L
            G = 0.0E0
C
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
C
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
C
  380    D(I) = Z(I,I)
         Z(I,I) = 1.0E0
         IF (L .LT. 1) GO TO 500
C
         DO 400 J = 1, L
            Z(I,J) = 0.0E0
            Z(J,I) = 0.0E0
  400    CONTINUE
C
  500 CONTINUE
C
      RETURN
      END

Cc      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
Cc      implicit none
CcC***BEGIN PROLOGUE  PYTHAG
CcC***REFER TO  EISDOC
CcC
CcC     Finds sqrt(A**2+B**2) without overflow or destructive underflow
CcC***ROUTINES CALLED  (NONE)
CcC***END PROLOGUE  PYTHAG
Cc      DOUBLE PRECISION A,B
CcC
Cc      DOUBLE PRECISION P,Q,R,S,T
CcC***FIRST EXECUTABLE STATEMENT  PYTHAG
Cc      P = DMAX1(ABS(A),ABS(B))
Cc      Q = DMIN1(ABS(A),ABS(B))
Cc      IF (Q .EQ. 0.0D+00) GO TO 20
Cc   10 CONTINUE
Cc         R = (Q/P)**2
Cc         T = 4.0E0 + R
Cc         IF (T .EQ. 4.0E0) GO TO 20
Cc         S = R/T
Cc         P = P + 2.0E0*P*S
Cc         Q = Q*S
Cc      GO TO 10
Cc   20 PYTHAG = P
Cc      RETURN
Cc      END
      double precision function pythag(a,b)
      double precision a,b
c
c   finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end

      subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)
c
      integer i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr
      double precision a(nm,n),w(n),u(nm,n),v(nm,n),rv1(n)
      double precision c,f,g,h,s,x,y,z,tst1,tst2,scale,pythag
      logical matu,matv
c
c   this subroutine is a translation of the algol procedure svd,
c   num. math. 14, 403-420(1970) by golub and reinsch.
c   handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
c
c   this subroutine determines the singular value decomposition
c        t
c   a=usv  of a real m by n rectangular matrix.  householder
c   bidiagonalization and a variant of the qr algorithm are used.
c
c   on input
c
c      nm must be set to the row dimension of two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.  note that nm must be at least
c        as large as the maximum of m and n.
c
c      m is the number of rows of a (and u).
c
c      n is the number of columns of a (and u) and the order of v.
c
c      a contains the rectangular input matrix to be decomposed.
c
c      matu should be set to .true. if the u matrix in the
c        decomposition is desired, and to .false. otherwise.
c
c      matv should be set to .true. if the v matrix in the
c        decomposition is desired, and to .false. otherwise.
c
c   on output
c
c      a is unaltered (unless overwritten by u or v).
c
c      w contains the n (non-negative) singular values of a (the
c        diagonal elements of s).  they are unordered.  if an
c        error exit is made, the singular values should be correct
c        for indices ierr+1,ierr+2,...,n.
c
c      u contains the matrix u (orthogonal column vectors) of the
c        decomposition if matu has been set to .true.  otherwise
c        u is used as a temporary array.  u may coincide with a.
c        if an error exit is made, the columns of u corresponding
c        to indices of correct singular values should be correct.
c
c      v contains the matrix v (orthogonal) of the decomposition if
c        matv has been set to .true.  otherwise v is not referenced.
c        v may also coincide with a if u is not needed.  if an error
c        exit is made, the columns of v corresponding to indices of
c        correct singular values should be correct.
c
c      ierr is set to
c        zero       for normal return,
c        k          if the k-th singular value has not been
c                   determined after 30 iterations.
c
c      rv1 is a temporary storage array.
c
c   calls pythag for  dsqrt(a*a + b*b) .
c
c   questions and comments should be directed to burton s. garbow,
c   mathematics and computer science div, argonne national laboratory
c
c   this version dated august 1983.
c
c   ------------------------------------------------------------------
c
      ierr = 0
c
      do 100 i = 1, m
c
         do 100 j = 1, n
            u(i,j) = a(i,j)
  100 continue
c   .......... householder reduction to bidiagonal form ..........
      g = 0.0d0
      scale = 0.0d0
      x = 0.0d0
c
      do 300 i = 1, n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m) go to 210
c
         do 120 k = i, m
  120    scale = scale + dabs(u(k,i))
c
         if (scale .eq. 0.0d0) go to 210
c
         do 130 k = i, m
            u(k,i) = u(k,i) / scale
            s = s + u(k,i)**2
  130    continue
c
         f = u(i,i)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,i) = f - g
         if (i .eq. n) go to 190
c
         do 150 j = l, n
            s = 0.0d0
c
            do 140 k = i, m
  140       s = s + u(k,i) * u(k,j)
c
            f = s / h
c
            do 150 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  150    continue
c
  190    do 200 k = i, m
  200    u(k,i) = scale * u(k,i)
c
  210    w(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m .or. i .eq. n) go to 290
c
         do 220 k = l, n
  220    scale = scale + dabs(u(i,k))
c
         if (scale .eq. 0.0d0) go to 290
c
         do 230 k = l, n
            u(i,k) = u(i,k) / scale
            s = s + u(i,k)**2
  230    continue
c
         f = u(i,l)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,l) = f - g
c
         do 240 k = l, n
  240    rv1(k) = u(i,k) / h
c
         if (i .eq. m) go to 270
c
         do 260 j = l, m
            s = 0.0d0
c
            do 250 k = l, n
  250       s = s + u(j,k) * u(i,k)
c
            do 260 k = l, n
               u(j,k) = u(j,k) + s * rv1(k)
  260    continue
c
  270    do 280 k = l, n
  280    u(i,k) = scale * u(i,k)
c
  290    x = dmax1(x,dabs(w(i))+dabs(rv1(i)))
  300 continue
c   .......... accumulation of right-hand transformations ..........
      if (.not. matv) go to 410
c   .......... for i=n step -1 until 1 do -- ..........
      do 400 ii = 1, n
         i = n + 1 - ii
         if (i .eq. n) go to 390
         if (g .eq. 0.0d0) go to 360
c
         do 320 j = l, n
c   .......... double division avoids possible underflow ..........
  320    v(j,i) = (u(i,j) / u(i,l)) / g
c
         do 350 j = l, n
            s = 0.0d0
c
            do 340 k = l, n
  340       s = s + u(i,k) * v(k,j)
c
            do 350 k = l, n
               v(k,j) = v(k,j) + s * v(k,i)
  350    continue
c
  360    do 380 j = l, n
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
  380    continue
c
  390    v(i,i) = 1.0d0
         g = rv1(i)
         l = i
  400 continue
c   .......... accumulation of left-hand transformations ..........
  410 if (.not. matu) go to 510
c   ..........for i=min(m,n) step -1 until 1 do -- ..........
      mn = n
      if (m .lt. n) mn = m
c
      do 500 ii = 1, mn
         i = mn + 1 - ii
         l = i + 1
         g = w(i)
         if (i .eq. n) go to 430
c
         do 420 j = l, n
  420    u(i,j) = 0.0d0
c
  430    if (g .eq. 0.0d0) go to 475
         if (i .eq. mn) go to 460
c
         do 450 j = l, n
            s = 0.0d0
c
            do 440 k = l, m
  440       s = s + u(k,i) * u(k,j)
c   .......... double division avoids possible underflow ..........
            f = (s / u(i,i)) / g
c
            do 450 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  450    continue
c
  460    do 470 j = i, m
  470    u(j,i) = u(j,i) / g
c
         go to 490
c
  475    do 480 j = i, m
  480    u(j,i) = 0.0d0
c
  490    u(i,i) = u(i,i) + 1.0d0
  500 continue
c   .......... diagonalization of the bidiagonal form ..........
  510 tst1 = x
c   .......... for k=n step -1 until 1 do -- ..........
      do 700 kk = 1, n
         k1 = n - kk
         k = k1 + 1
         its = 0
c   .......... test for splitting.
c              for l=k step -1 until 1 do -- ..........
  520    do 530 ll = 1, k
            l1 = k - ll
            l = l1 + 1
            tst2 = tst1 + dabs(rv1(l))
            if (tst2 .eq. tst1) go to 565
c   .......... rv1(1) is always zero, so there is no exit
c              through the bottom of the loop ..........
            tst2 = tst1 + dabs(w(l1))
            if (tst2 .eq. tst1) go to 540
  530    continue
c   .......... cancellation of rv1(l) if l greater than 1 ..........
  540    c = 0.0d0
         s = 1.0d0
c
         do 560 i = l, k
            f = s * rv1(i)
            rv1(i) = c * rv1(i)
            tst2 = tst1 + dabs(f)
            if (tst2 .eq. tst1) go to 565
            g = w(i)
            h = pythag(f,g)
            w(i) = h
            c = g / h
            s = -f / h
            if (.not. matu) go to 560
c
            do 550 j = 1, m
               y = u(j,l1)
               z = u(j,i)
               u(j,l1) = y * c + z * s
               u(j,i) = -y * s + z * c
  550       continue
c
  560    continue
c   .......... test for convergence ..........
  565    z = w(k)
         if (l .eq. k) go to 650
c   .......... shift from bottom 2 by 2 minor ..........
         if (its .eq. 30) go to 1000
         its = its + 1
         x = w(l)
         y = w(k1)
         g = rv1(k1)
         h = rv1(k)
         f = 0.5d0 * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
         g = pythag(f,1.0d0)
         f = x - (z / x) * z + (h / x) * (y / (f + dsign(g,f)) - h)
c   .......... next qr transformation ..........
         c = 1.0d0
         s = 1.0d0
c
         do 600 i1 = l, k1
            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = pythag(f,h)
            rv1(i1) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
            if (.not. matv) go to 575
c
            do 570 j = 1, n
               x = v(j,i1)
               z = v(j,i)
               v(j,i1) = x * c + z * s
               v(j,i) = -x * s + z * c
  570       continue
c
  575       z = pythag(f,h)
            w(i1) = z
c   .......... rotation can be arbitrary if z is zero ..........
            if (z .eq. 0.0d0) go to 580
            c = f / z
            s = h / z
  580       f = c * g + s * y
            x = -s * g + c * y
            if (.not. matu) go to 600
c
            do 590 j = 1, m
               y = u(j,i1)
               z = u(j,i)
               u(j,i1) = y * c + z * s
               u(j,i) = -y * s + z * c
  590       continue
c
  600    continue
c
         rv1(l) = 0.0d0
         rv1(k) = f
         w(k) = x
         go to 520
c   .......... convergence ..........
  650    if (z .ge. 0.0d0) go to 700
c   .......... w(k) is made non-negative ..........
         w(k) = -z
         if (.not. matv) go to 700
c
         do 690 j = 1, n
  690    v(j,k) = -v(j,k)
c
  700 continue
c
      go to 1001
c   .......... set error -- no convergence to a
c              singular value after 30 iterations ..........
 1000 ierr = k
 1001 return
      end
