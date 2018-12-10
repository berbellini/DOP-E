c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME IV                                                       c
c                                                                      c
c      PROGRAM: SRFINV                                                 c
c                                                                      c
c      COPYRIGHT 1986, 1991                                            c
c      D. R. Russell, R. B. Herrmann                                   c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c  
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
      program srfinv96
c
c       CHANGES
c       31 OCT 2002 - dimension of NL was changes to ensure 
c           that sequential
c           accumulation works as number of layers is changes
c
c
c     Program SURFINV calculates a singular value decomposition
c     of an arbitrary matrix with n rows and m columns.  The
c     maximum number of columns is 80.  The maximum number of
c     rows is UNLIMITED.  Sequential Householder transformations
c     are used to reduce the matrix to upper triangular form,
c     and then a singular value decomposition is performed on
c     the triangular matrix.  See Lawson and Hanson, 'Solving
c     Least Squares Problems', Prentice-Hall, 1974.
c
c     Column weighting is allowed, and either no or
c     differential smoothing constraints can be used.  'Differential
c     smoothing' implies putting a constraint on the L2 norm of
c     the DIFFERENCE between adjacent model elements, rather than
c     on the L2 norm of the elements themselves, as in no
c     smoothing.
c
c     Program developed by David R. Russell, St. Louis University,
c     Jan. 1984.
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        parameter (LIN=5,LOT=6,LER=0,NL=200,NL2=NL+NL,NROW=NL+NL+NL+1)
c-----
c       LIN - logical unit for FORTRAN read from terminal
c       LOT - logical unit for FORTRAN write to terminal
c       LER - logical unit for FORTRAN error write
c       NL  - number of layers in model
c       NL2 - number of columns in model (1st NL=velocity, 
c             2nd NL=Q inverse)
c       NROW    - maximum number of rows to read in at a time
c             after first set is read in then reads one at a time
c                 For efficient read in NROW >> NL+NL+1
c             of course this trades off against matrix size.
c           The idea is that we initially read in NROWs of data, 
c             where NROW >= NL2
c           Then the H12 routine works on the entire matrix and 
c             stores the result
c           in the first NL2 rows of the A matrix. Subsequent 
c             reads fill in the
c           NL2 +1 -> NROW rows, and augment the first NL2 rows.
c-----
      real a(NROW,NL2+2),x(NL2),b(NROW),h(3*NL2),dd(NL2)
        logical wc(NL2)
        real w(NL2,NL2), winv(NL2,NL2), ta(NL2,NL2)
        real tta(NL2,NL2)
      double precision sum
c-----
c     machine dependent initialization
c-----
      call mchdep()
c-----
      open(1,file='tmpsrfi.09',form='unformatted',access='sequential')
      rewind 1
c-----
c       m   number of unknowns in the matrix should be
c           twice the number of layers
c       nfilt   0   no weighting no smoothing
c           1      weighting no smoothing
c           2       no weighting    smoothing
c           3          weighting    smoothing
c----
      read(1) m,nfilt
      do 5 i=1,m
   5  dd(i)=1.0
      read(1)(dd(i),i=1,m)
      read(1)(wc(i),i=1,m)
C       WRITE(6,*)(wc(i),i=1,m)


c-----
c       form the w and winv matrices for application of the smoothing
c       condition
c-----
        call mkw(W,WC,m,NL2)
        call mkwinv(WINV,WC,m,NL2)
c
c.....set up rows of a and b for first iteration.
c
      n1=NROW
      m2=m/2
      ntot=0
      l=0
      itst=0
      iend=0
      i1=1
      m1=m+1
      n2=n1
      go to 12
  10  l=m+1
      i1=m+2
  12  do 15 j=i1,n1
      read(1,end=18)(a(j,k),k=1,m1)
      l=l+1
      ntot=ntot+1
   15 continue
      go to 20
   18 if(l.lt.i1) go to 30
      iend=1
      n2=l
   20 continue
c
c.....use householder transformations to reduce a to an upper
c     triangular matrix.
c
        if(itst.eq.0)then
            itst=1
            do 23 j=1,m1
                call h12(1,j,j+1,n2,a(1,j),1,t,a(1,j+1),1,n1,m1-j)
  23        continue
        else
            do 25 j=1,m1
                call h12(1,j,m+2,n2,a(1,j),1,t,a(1,j+1),1,n1,m1-j)
  25        continue
        endif
      if(iend.eq.0) go to 10
c
c.....end of sequential iterations.
c.....force upper diagonal to remove slight numerical error
c
  30  continue
      do 3020 i=2,m
      k=i-1
      do 3020 j=1,k
      a(i,j)=0.0
 3020 continue
      do 32 i=1,m
      b(i)=a(i,m1)
  32  continue
      err=a(m1,m1)
      open(2,file='tmpsrfi.10',form='unformatted',access='sequential')
      open(3,file='tmpsrfi.11',form='unformatted',access='sequential')
      rewind 2
      rewind 3
c-----
c       copy a into taa
c-----
        do 4100 j=1,m
            do 4200 i=1,m
                tta(i,j) = a(i,j)
 4200       continue
 4100   continue
c-----
c       multiply by the W inverse matrix
c-----
        call mymatmul(ttA,WINV,tA,m,NL2)
c
c....... a priori column weighting
c
        do 667 j=1,m
            do 666 i=1,m
                ta(i,j)=ta(i,j)*dd(j)
 666        continue
 667    continue
c-----
c       get norm of upper left triangular matrix
c-----
      s1=0.0
      s2=0.0
      do 2000 i=1,m2
      do 2000 j=i,m2
2000  s1=s1+ta(i,j)**2
      s1=sqrt(s1)
c-----
c       get norm of lower right triangular matrix
c-----
      do 2100 i=m2+1,m
      do 2100 j=i,m
2100  s2=s2+ta(i,j)**2
      s2=sqrt(s2)
c-----
c       apply column weighting to control size of eigenvalues
c       We assume that the eigenvalues associated with the
c       first m/2 model parameters can be of vastly different
c       size than those of the second m/2. We do this so
c       that we only require a single damping value to
c       control the solution
c
c       To accomplish this we use the Euclidian or Frobenius
c       norm to adjust the matrix and model space
c-----
c       s1  = norm of the velocity only portion
c       s2  = norm of q only portion
c
c-----
        if(s1.gt.0.0 .and. s2.eq.0.0)then
            srat12 = 1.0
        elseif(s1.eq.0.0 .and. s2.gt.0.0)then
            srat12 = 1.0
        elseif(s1.eq.0.0 .and. s2.eq.0.0)then
            srat12 = 1.0
        else
            srat12 = s1/s2
        endif
C       write(0,*)'surfinv:s1,s2,srat1,srat2',s1,s2,srat12
        do 2200 j=1,m
            if(j.le.m2)then
                srat = 1.0
            else
                srat = srat12
            endif
            dd(j)=dd(j)*srat
            do 2201 i=1,j
                ta(i,j)=ta(i,j)*srat
 2201       continue
 2200   continue
c-----
c       get norm of rescaled matrix
c-----
        s1 = 0.0
        do 2300 i=1,m
            do 2301 j=i,m
                s1 = s1 + ta(i,j)**2
 2301       continue
 2300   continue
        s1 = sqrt(s1)
c-----
c       rescale resultant matrix and solution vector 
c       so that mean norm is
c       10.0. This will mean that a damping of 1.0 will cause the
c       model to move slightly
c-----
        if(s1.eq.0.0)s1 = 1.0
        s1 = s1/10.0
        do 2302 i=1,m
            dd(i) = dd(i)/s1
            do 2303 j=i,m
                a(i,j) = ta(i,j)/s1
 2303       continue
 2302   continue
c
c.......single value decomposition
c
        call svdrs(a,n1,m,m,b,n1,1,h)
c-----
c       The current resolution matrix is
c                     2     2  -1    2  T
c       R(y) = V ( lam + sig  )   lam  V
c
c       since  y    = R  y
c               est       true
c
c       and y = Wx  in terms of x we have
c               -1       2     2  -1    2  T
c       R(x) = W  V ( lam + sig  )   lam  V W
c                                               -1
c       Since solution for x involves the term W  V, that    T
c       will be formed and saved below together with the    U b
c       and the eigenvalues.
c-----
c       To compute the resolving kernel we need to form
c       the matrix
c        T
c       V W W
c          1 2
c       and output the results by row
c-----
c                    T
c       first form  V W  using the fact that W  is a diagonal matrix
c                      1                      1
c-----
        do 770 i=1,m
            do 760 j=1,m
                ta(i,j) = a(j,i)/dd(j)
  760       continue
  770   continue
c-----
c       now form
c        T
c       V W W
c          1 2
c-----
        do 800 i=1,m
            do 790 j=1,m
                sum = 0.0d+00
                do 780 k=1,m
                    sum = sum + ta(i,k)*W(k,j)
 780            continue
                x(j) = sum
 790        continue
            write(3)(x(k),k=1,m)
 800    continue
c-----
c       to obtain the solution vector we need
c            -1  -1        2     2   -1      T
c       x = W   W   V ( lam + sig I )  lam  U b
c            2   1
c       since the weighting matrix was W = W1 W2 where W1 is a diagonal
c       matrix consisting of the 1/dd(i) and W2 is the 
c       general smoothing matrix
c                      -1            T
c       we will pass  W  V, lam and U b to the calling program
c                     
c-----               -1
c       first form  W   V
c                    1
c-----
        do 690 i=1,m
            do 695 j=1,m
                ta(i,j)=dd(i)*a(i,j)
  695       continue
  690   continue
c-----             -1  -1
c       now form  W   W  V
c                   2  1
c-----
c-----
c
c       multiply by W inverse matrix
c       This is now W sup -1 V
c-----
        call mymatmul(WINV,tA,ttA,m,NL2)
      klim=0
      do 970 i=1,m
      if(h(i).lt.1.e-15) go to 980
      klim=klim+1
 970  continue
c
c     ntot: total number of rows
c     klim: number of non-zero eigenvalues
c     m:    number of columns
c     err:  norm of minimum error between observed and computed
c             values:  err = min || b - Ax ||
c     h:    vector of eigenvalues
c     b:    vector of data eigenvector matrix U times residual
c           vector b: 'utb vector'.
c     a:    model eigenvector matrix V.
c
 980  write(2) ntot,klim,m,err
      write(2) (h(i),i=1,m)
      write(2) (b(i),i=1,m)
      do 80 i=1,m
      write(2) (tta(i,j),j=1,m)
  80  continue
      do 85 i=1,3
      close(i,status='keep')
  85  continue
c     stop
      end
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        subroutine svdrs(a,mda,mm,nn,b,mdb,nb,s)
c-----
c     c.l. lawson and r.j. hanson jet propulsion laboratory, 1 mar 74.
c     see 'solving least squares problems', prentice-hall, 1974.
c
c     the array s occupies 3*n cells.
c     a occupies m*n cells
c     b occupies m*nb cells.
c
c     special singular value decomposition subroutine.
c     we have the mxn matrix a and the system a*x=b to solve.
c     either m .ge. n or m .lt. n is permitted.
c     the singular value decomposition
c                  a = u*s*v**(t)
c     is made in such a way that one gets
c     (1) the matrix v in the first n rows and columns of a.
c     (2) the diagonal matrix or ordered singular values in
c         the first n cells of the array s(ip).  ip .ge. 3*n.
c     (3) the matrix product u**(t)*b=g gets placed back in b.
c     (4) the user must complete the solution and do his own
c         singular value analysis.
c-----
        parameter (LIN=5,LOT=6,LER=0)
      dimension a(mda,nn),b(mdb,nb),s(nn,3)
      zero=0.
      one=1.
      n=nn
      if(n.le.0.or.mm.le.0) return
      j=n
  10  continue
      do 20 i=1,mm
C      if(a(i,j)) 50,20,50
      if(a(i,j) .ne. 0.0) go to 50
  20  continue
      if(j.eq.n) go to 40
      do 30 i=1,mm
  30  a(i,j)=a(i,n)
  40  continue
      a(1,n)=j
      n=n-1
  50  continue
      j=j-1
      if(j.ge.1) go to 10
      ns=0
      if(n.eq.0) go to 240
      i=1
      m=mm
  60  if(i.gt.n.or.i.ge.m) go to 150
C      if(a(i,i)) 90,70,90
      if(a(i,i) .ne. 0.0) go to 90
C  70  do 80 j=1,n
      do 80 j=1,n
C      if(a(i,j)) 90,80,90
      if(a(i,j) .ne.0.0) go to 90
  80  continue
      go to 100
  90  i=i+1
      go to 60
 100  if(nb.le.0) go to 115
      do 110 j=1,nb
      t=b(i,j)
      b(i,j)=b(m,j)
 110  b(m,j)=t
 115  do 120 j=1,n
 120  a(i,j)=a(m,j)
      if(m.gt.n) go to 140
      do 130 j=1,n
 130  a(m,j)=zero
 140  continue
      m=m-1
      go to 60
 150  continue
      l=min0(m,n)
c
c     the following loop reduces a to upper bidiagonal and
c     also applies the premultiplying transformations to b.
c
      do 170 j=1,l
      if (j.ge.m) go to 160
      call h12(1,j,j+1,m,a(1,j),1,t,a(1,j+1),1,mda,n-j)
      call h12(2,j,j+1,m,a(1,j),1,t,b,1,mdb,nb)
 160  if (j.ge.n-1) go to 170
      call h12(1,j+1,j+2,n,a(j,1),mda,s(j,3),a(j+1,1),mda,1,m-j)
 170  continue
c
c     copy the bidiagonal matrix into the array s() for qrbd.
c
      if (n.eq.1) go to 190
      do 180 j=2,n
      s(j,1)=a(j,j)
 180  s(j,2)=a(j-1,j)
 190  s(1,1)=a(1,1)
c
      ns=n
      if (m.ge.n) go to 200
      ns=m+1
      s(ns,1)=zero
      s(ns,2)=a(m,m+1)
 200  continue
c
c     construct the explicit nxn by product matrix, w=q1*q2*..*ql*i
c     in the array a().
c
      do 230 k=1,n
      i=n+1-k
      if(i.ge.n-1) go to 210
      call h12(2,i+1,i+2,n,a(i,1),mda,s(i,3),a(1,i+1),1,mda,n-i)
 210  do 220 j=1,n
 220  a(i,j)=zero
 230  a(i,i)=one
c
c     compute the svd of the bidiagonal matrix
c
      call qrbd(ipass,s(1,1),s(1,2),ns,a,mda,n,b,mdb,nb)
c
      go to (240,310), ipass
 240  continue
      if(ns.ge.n) go to 260
      nsp1=ns+1
      do 250 j=nsp1,n
 250  s(j,1)=zero
 260  continue
      if(n.eq.nn) return
      np1=n+1
      do 280 j=np1,nn
      s(j,1)=a(1,j)
      do 270 i=1,n
 270  a(i,j)=zero
 280  continue
      do 300 k=np1,nn
      i=s(k,1)
      s(k,1)=zero
      do 290 j=1,nn
      a(k,j)=a(i,j)
 290  a(i,j)=zero
      a(i,k)=one
 300  continue
      return
 310  write(LOT,320)
      stop
 320  format(' convergence failure in qr bidiagonal svd routine')
      end
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine h12(mode,lpivot,l1,m,u,iue,up,c,ice,icv,ncv)
c
c     construction and/or application of a single
c     Householder transformation..     q = i+u*(u**t)/b
c
c     mode    = 1 or 2   to select algorithm h1 or h2.
c     lpivot is the index of the pivot element.
c     l1,m   if l1 .le. m   the transformation will be constructed to
c            zero elements indexed from l1 through m.   if l1 gt. m
c            the subroutine dose an identity transformation.
c     u(),iuw,up    on entry to h1 u() contains the pivot vector.
c                   iue is the storage increment between elements.
c                                       on exit from h1 u() and up
c                   contain quantities defining the vector u of the
c                   Householder transformation.   on entry to h2 u()
c                   and up should contain quantities previously computed
c                   by h1.  these will not be modified by h2.
c     c()    on entry to h1 or h2 c() contains a matrix which will be
c            regarded as a set of vectors to which the Householder
c            transformation is to be applied.  on exit c() contains the
c            set of transformed vectors.
c     ice    storage increment between elements of vectors in c().
c     icv    storage increment between vectors in c().
c     ncv    number of vectors in c() to be transformed. if ncv .le. 0
c            no operations will be done on c().
c
      dimension u(iue,m),c(*)
      double precision sm,b
      one=1.
c
      if (0.ge.lpivot.or.lpivot.ge.l1.or.l1.gt.m) return
      cl=abs(u(1,lpivot))
      if (mode.eq.2) go to 60
c            ****** construct the transformation. ******
          do 10 j=l1,m
  10      cl=amax1(abs(u(1,j)),cl)
C      if (cl) 130,130,20
      if (cl .le. 0.0) go to 130
C  20  clinv=one/cl
      clinv=one/cl
      sm=(dble(u(1,lpivot))*clinv)**2
          do 30 j=l1,m
  30      sm=sm+(dble(u(1,j))*clinv)**2
c             convert dble. prec. sm to sngl. prec. sm1
      sm1=sm
      cl=cl*sqrt(sm1)
C      if (u(1,lpivot)) 50,50,40
      if (u(1,lpivot) .le. 0.0) go to 50
C  40  cl=-cl
      cl=-cl
  50  up=u(1,lpivot)-cl
      u(1,lpivot)=cl
      go to 70
c     ****** apply the transformation 1+u*(u**t)/b to c. ******
c
C  60  if (cl) 130,130,70
  60  if (cl .le. 0.0) go to 130
  70  if (ncv.le.0) return
      b=dble(up)*u(1,lpivot)
c          b  must be nonpositive here.  if b=0., return.
c
C      if (b) 80,130,130
      if (b .ge. 0.0) go to 130
C  80  b=one/b
      b=one/b
      i2=1-icv+ice*(lpivot-1)
      incr=ice*(l1-lpivot)
          do 120 j=1,ncv
          i2=i2+icv
          i3=i2+incr
          i4=i3
          sm=c(i2)*dble(up)
              do 90 i=l1,m
              sm=sm+c(i3)*dble(u(1,i))
  90          i3=i3+ice
C          if (sm) 100,120,100
          if (sm .eq. 0.0) go to 120
C 100      sm=sm*b
          sm=sm*b
          c(i2)=c(i2)+sm*dble(up)
              do 110 i=l1,m
              c(i4)=c(i4)+sm*dble(u(1,i))
 110          i4=i4+ice
 120      continue
 130      return
          end
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine qrbd(ipass,q,e,nn,v,mdv,nrv,c,mdc,ncc)
c
c     c.l. lawson and r.j. hanson, jet propulsion laboratory, 
c          1973 jun 12.
c     see 'solving least squares problems', prentice-hall, 1974.
c          qr algorithm for singular values of a bidiagonal matrix.
c
      real diff
      logical wntv,havers,fail
      dimension q(nn),e(nn),v(mdv,nn),c(mdc,ncc)
      zero=0.
      one=1.
      two=2.
c
      n=nn
      ipass=1
      if (n.le.0) return
      n10=10*n
      wntv=nrv.gt.0
      havers=ncc.gt.0
      fail=.false.
      nqrs=0
      e(1)=zero
      dnorm=zero
      do 10 j=1,n
  10  dnorm=amax1(abs(q(j))+abs(e(j)),dnorm)
      do 200 kk=1,n
      k=n+1-kk
c
c     test for splitting or rank deficiencies
c     first make test for last diagonal term, q(k), being small.
c
  20  if(k.eq.1) go to 50
C      if(diff(dnorm+q(k),dnorm)) 50,25,50
      if( diff(dnorm+q(k),dnorm) .ne. 0.0) go to 50
c
c     since q(k) is small we will make a special pass to
c     transform e(k) to zero.
c
C  25  cs=zero
      cs=zero
      sn=-one
      do 40 ii=2,k
      i=k+1-ii
      f=-sn*e(i+1)
      e(i+1)=cs*e(i+1)
      call g1(q(i),f,cs,sn,q(i))
c     transformation constructed to zero position (i,k).
c
      if (.not.wntv) go to 40
      do 30 j=1,nrv
  30  call g2(cs,sn,v(j,i),v(j,k))
c
c     accumulate rt. transformations in v.
c
  40  continue
c
c     the matrix is now bidiagonal, and of lower order
c     since e(k) .eq. zero..
c
  50  do 60 ll=1,k
      l=k+1-ll
C      if(diff(dnorm+e(l),dnorm)) 55,100,55
      if(diff(dnorm+e(l),dnorm) .eq.0.0) go to 100
C  55  if(diff(dnorm+q(l-1),dnorm)) 60,70,60
      if( diff(dnorm+q(l-1),dnorm) .eq.0.0 ) go to 70
  60  continue
c
c     this loop can't completesince e(1)=zero.
c
      go to 100
c
c     cancellation of e(l), l.gt.1.
c
  70  cs=zero
      sn=-one
      do 90 i=l,k
      f=-sn*e(i)
      e(i)=cs*e(i)
C      if(diff(dnorm+f,dnorm)) 75,100,75
      if(diff(dnorm+f,dnorm) .eq. 0.0) go to 100
C  75  call g1(q(i),f,cs,sn,q(i))
      call g1(q(i),f,cs,sn,q(i))
      if (.not.havers) go to 90
      do 80 j=1,ncc
  80  call g2(cs,sn,c(i,j),c(l-1,j))
  90  continue
c
c     test for convergence..
c
 100  z=q(k)
      if(l.eq.k) go to 170
c
c     shift from bottom 2 by 2 minor of b**(t)*b.
c
      x=q(l)
      y=q(k-1)
      g=e(k-1)
      h=e(k)
      f=((y-z)*(y+z)+(g-h)*(g+h))/(two*h*y)
      g=sqrt(one+f**2)
      if(f.lt.zero) go to 110
      t=f+g
      go to 120
  110 t=f-g
  120 f=((x-z)*(x+z)+h*(y/t-h))/x
c
c     next qr sweep
c
      cs=one
      sn=one
      lp1=l+1
      do 160 i=lp1,k
      g=e(i)
      y=q(i)
      h=sn*g
      g=cs*g
      call g1(f,h,cs,sn,e(i-1))
      f=x*cs+g*sn
      g=-x*sn+g*cs
      h=y*sn
      y=y*cs
      if (.not.wntv) go to 140
c
c     accumulate rotations (from the right) in 'v'
c
      do 130 j=1,nrv
 130  call g2(cs,sn,v(j,i-1),v(j,i))
 140  call g1(f,h,cs,sn,q(i-1))
      f=cs*g+sn*y
      x=-sn*g+cs*y
      if(.not.havers) go to 160
      do 150 j=1,ncc
 150  call g2(cs,sn,c(i-1,j),c(i,j))
c
c     apply rotations from the left to
c     right hand sides in 'c'..
c
 160  continue
      e(l)=zero
      e(k)=f
      q(k)=x
      nqrs=nqrs+1
      if(nqrs.le.n10) go to 20
c
c     return to 'test for splitting'.
c
      fail=.true.
c
c     cutoff for convergence failure. 'nqrs' will be 2*n usually.
c
 170  if(z.ge.zero) go to 190
      q(k)=-z
      if(.not.wntv) go to 190
      do 180 j=1,nrv
 180  v(j,k)=-v(j,k)
 190  continue
c
c     convergence. q(k) is made nonnegative..
c
 200  continue
      if(n.eq.1) return
      do 210 i=2,n
      if(q(i).gt.q(i-1)) go to 220
 210  continue
      if (fail) ipass=2
      return
c
c     every singular value is in order..
c
 220  do 270 i=2,n
      t=q(i-1)
      k=i-1
      do 230 j=i,n
      if(t.ge.q(j)) go to 230
      t=q(j)
      k=j
 230  continue
      if(k.eq.i-1) go to 270
      q(k)=q(i-1)
      q(i-1)=t
      if(.not.havers) go to 250
      do 240 j=1,ncc
      t=c(i-1,j)
      c(i-1,j)=c(k,j)
 240  c(k,j)=t
 250  if(.not.wntv) go to 270
      do 260 j=1,nrv
      t=v(j,i-1)
      v(j,i-1)=v(j,k)
 260  v(j,k)=t
 270  continue
c
c     end of ordering algorithm
c
      if (fail) ipass=2
      return
      end
c
c- - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine g1(a,b,cos,sin,sig)
c
c     compute orthogonal rotation matrix.
c     compute.. matrix (c,s) so that (c,s)(a)=(sqrt(a**2+b**2))
c                      (-s,c)        (-s,c)(b) (0)
c     compute sig=sqrt(a**2+b**2)
c        sig is computed last to allow for the possibility that
c        sig may be in the same location as a or b.
c
       zero=0.
       one=1.
       if(abs(a).gt.abs(b)) then
           xr=b/a
           yr=sqrt(one+xr**2)
           cos=sign(one/yr,a)
           sin=cos*xr
           sig=abs(a)*yr
       else
           if (b.ne.0.0) then
               xr=a/b
               yr=sqrt(one+xr**2)
               sin=sign(one/yr,b)
               cos=sin*xr
               sig=abs(b)*yr
           else
               sig=zero
                   cos=zero
                   sin=one
           endif
       endif
       return
       end
c
c- - - - - - - - - - - - - - - - - - - - - -
c
      subroutine g2(cos,sin,x,y)
c
c     apply the rotation computed bu g1.f to (x,y).
c
      xr=cos*x+sin*y
      y=-sin*x+cos*y
      x=xr
      return
      end
c
c- - - - - - - - - - - - - - - - - - - - - - -
c
      function diff(x,y)
      diff=x-y
      return
      end

        subroutine mkw(W,WC,m,NL2)
        real W(NL2,NL2)
        logical WC(NL2)
        integer NL2
        integer m
        do 2000 i=1,m
            do 2100 j=1,m
                W(i,j) = 0.0
 2100       continue
            W(i,i) = 1.0
            if(i.lt.m .and. wc(i))then
                W(i,i+1) = -1.0
                
            endif
 2000   continue
        return
        end

        subroutine mkwinv(WI,WC,m,NL2)
c-----
c       form the W inverse matrix
c-----
        real WI(NL2,NL2)
        logical WC(NL2)
        integer m
        integer NL2
c-----
c       create the W^-1 matrix
c-----
c       first fill as differential
c-----
        do 3000 i=1,m
            do 3100 j=1,m
                if(j.lt.i)then
                    WI(i,j) = 0.0
                else
                    WI(i,j) = 1.0
                endif
 3100       continue
 3000   continue
c-----
c       now put in the constraint
c-----
        do 3500 i=1,m-1
            if(.not.wc(i))then
                do 3600 ii=1,i
                    do 3610 jj=i+1,m
                    WI(ii,jj) = 0.0
 3610               continue
 3600           continue
            endif
 3500   continue
        return
        end

        subroutine mymatmul(A,B,C,m,NARR)
c-----
c       C = A B
c-----
        real A(NARR,NARR)
        real B(NARR,NARR)
        real C(NARR,NARR)
        integer NARR, m
        do 1000 i=1,m
            do 2000 j=1,m
                sum = 0.0
                do 3000 k=1,m
                    sum = sum + A(i,k)*B(k,j)
 3000           continue
                C(i,j) = sum
 2000       continue
 1000   continue
        return
        end
