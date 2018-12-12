      program tdisp96
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME IV                                                       c
c                                                                      c
c      PROGRAM: TDISP96                                                c
c                                                                      c
c      COPYRIGHT 1992 ,2010                                            c
c      R. B. Herrmann                                                  c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c       Changes
c       13 SEP 2000 - initialized the fpar and ipar arrays
c           before output
c       14 OCT 2006 - corrected initial conditions in
c           subroutine dltar4 for lfluid case
c       28 DEC 2007 - changed the Earth flattening to now use layer
c           midpoint and the Biswas (1972: PAGEOPH 96, 61-74, 1972)
c           density mapping for P-SV  - note a true comparison
c           requires the ability to handle a fluid core for SH and SV
c       23 AUG 2010 - created
c       20 AMR 2012 - corrected compiler warnings
c       09 SEP 2012 - correct last line of bsort from n=n-m to nn=n-m
c             Thanks to ruoshan at ustc
c-----
c                                                                      c
c                                                                      c
c     The program calculates the dispersion values for any             c
c     layered model, any frequency, and any mode.                      c
c                                                                      c
c     This program will accept one liquid layer at the surface.        c
c     In such case ellipticity of rayleigh wave is that at the         c
c     top of solid array.  Love wave communications ignore             c
c     liquid layer.                                                    c
c                                                                      c
c     The number of periods is not limitted (try 16384), but           c
c     the max mode number = 2000. This number appears at common/phas/  c
c     and common/stor/ , which can be changed by user before compilation
c     Other places in 'surface wave' series which have this restrictionc
c     are in 'srfcomb85', after which the total mode can become 2000,  c
c     then in 'reigen85' and 'leigen85'.                               c
c                                                                      c
c     Pole searching can be done either on C-T or F-K domain.          c
c     Pole searching is improved around C=S- or P-velocities.          c
c                                                                      c
c     Program developed by Robert B Herrmann Saint Louis               c
c     univ. Nov 1971, and revised by C. Y. Wang on Oct 1981.           c
c     Version '85' is developed at Feb 1985.                           c
c                                                                      c
c                                                                      c
c----------------------------------------------------------------------c
c                                                                      c
c     REFERENCE:                                                       c
c                                                                      c
c      Aki & Richards (1980) Chapter 7.                                c
c      Haskell, N.A. (1964) BSSA 54, 377-393.                          c
c      Dunkin, J.W.  (1965) BSSA 55, 335-358.                          c
c      Herrmann,R.B. (1979) BSSA 69, 1-16.                             c
c      Wang, C.Y. & R.B. Herrmann (1980) BSSA 70, 1015-1036.           c
c      Wang, C.Y. (1981) Saint Louis U. Ph.D. Dissertation.            c
c                                                                      c
c----------------------------------------------------------------------c
c       MODIFICATION HISTORY
c
c       24 MAY 2000 - permit all fluid layers
c-----
        parameter (LIN=5, LOT=6, NL=1000, NL2=NL+NL)
        implicit double precision (a-h,o-z)
        real*4 vth(NL2),vtp(NL2)
        real*4 vts(NL2,2)
        equivalence(vth(1),vts(1,1)),(vtp(1),vts(1,2))
        common/pari/ mmax,mode
        common/water/iwat(NL)
        common/pard/ twopi,displ,dispr
        common/vels/ mvts(2),vts

        common/earth/radius
        real radius

        parameter(NPERIOD=2049)
        real fval(NPERIOD)
        integer nfval

        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep
        
        common/timodel/od(NL),oTA(NL),oTC(NL),oTL(NL),oTN(NL),oTF(NL),
     1      oTRho(NL),
     2      oqa(NL),oqb(NL),oetap(NL),oetas(NL), 
     3      ofrefp(NL), ofrefs(NL)
        real od,oTA,oTC,oTN,oTL,oTF,oTRho,oqa,oqb,
     1        oetap,oetas,ofrefp,ofrefs

        real*4 ccmin, ccmax

c-----
c       parameters for model file
c-----
        character title*80
        logical ext
c-----
c       parameters from tdisp96.dat
c-----
        character mname*80
        logical dolove, dorayl
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid
c-----
c       command flags from tdisp96.dat
c-----
c
        real*4 dt
c-----
c       command line variables
c-----
        logical verby
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       initialize
c-----
        radius = 6371.
c-----
c       parse command line flags
c-----
        call gcmdln(verby,ccmin,ccmax)
c-----
c       main program just does the control.
c-----
        twopi=6.283185307179586d+00
c-----
c       get control parameters from tdisp96.dat
c-----
        inquire(file='tdisp96.dat',exist=ext)
        if(.not.ext)then
            call usage('Control file tdisp96.dat'//
     1          ' does not exist')
        endif
        open(1,file='tdisp96.dat',form='formatted',status='unknown',
     1      access='sequential')
        rewind 1
        read(1,*,end=9999,err=9999)dt
C        WRITE(*,*)dt
        read(1,*,end=9999,err=9999)npts,n1,n2
C        WRITE(*,*)npts,n1,n2
        read(1,'(a)',end=9999,err=9999)mname
C        WRITE(*,*)mname
        read(1,*)dolove, dorayl
C        WRITE(*,*)dolove, dorayl
        read(1,*)hs, hr
C        WRITE(*,*)hs,hr
        read(1,*)nmodes
C        WRITE(*,*)nmodes
        read(1,*)faclov,facray
C        WRITE(*,*)faclov,facray
c-----
c       get user specified frequencies, but also so an error check
c-----
        if(n1.lt.0 )then
            nfval = npts
            npts = 0
            do 1010 i=1,nfval
                read(1,*,end=1011,err=1011)xx
                if(xx.ge.0.0)then
                    npts = npts + 1
                    fval(npts) = xx
                endif
 1010       continue
 1011       continue
c-----
c       safety
c-----
            if(npts.gt.1)then
                call bsort(npts,fval,-1)
            endif
            nfval = npts
        else
            nfval = -1
        endif
        close(1)
c-----
c       get the earth model
c-----
        inquire(file=mname,exist=ext)
        if(.not. ext)call usage('Model file does not exist')
        l = lgstr(mname)

        write(LOT,*)'Model name: ',mname(1:l)
        call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.true.)
CCCCC drop out if not flat 1-D constant velocity isotropic
c-----
c       check for water only
c-----
        allfluid = .true.
        do 1200 i=1,mmax
            if(oTN(i).gt.0.0)then
                allfluid = .false.
            endif
 1200   continue
        if(allfluid)then
            dolove = .false.
        endif
c-----
c       set up controls
c-----
        if(dolove)then
            idispl = 1
        else
            idispl = 0
        endif
        if(dorayl)then
            idispr = 1
        else
            idispr = 0
        endif
        mode = nmodes
c-----
c       mmax = number of layer including the half-space.
c       mode = maximum number of modes to be found at a frequency,
c              can be as large as 1500.
c-----
        if(mode.gt.2000)then
            mode=2000
        endif
c-----
c       process dispersion
c-----
        if(faclov.lt.1.0)faclov = 5.0
        if(facray.lt.1.0)facray = 5.0
        displ = faclov
        dispr = facray
c-----
c       get the Love and Rayleigh dispersion as two separate
c       operations.
c       when we are done update the state files. 
c-----
        do 1000 ilvry=1,2
            if(ilvry.eq.1)then
                if(idispl.eq.1)then
                    call disprs(ilvry,dt,npts,iret,mname,
     1                  verby,nfval,fval,ccmin,ccmax)
                endif
            else
                if(idispr.eq.1)then
                    call disprs(ilvry,dt,npts,iret,mname,
     1                  verby,nfval,fval,ccmin,ccmax)
                endif
            endif
 1000   continue
        go to 999
 9999   continue
            write(LOT,*)'Error in tdisp96.dat'
            stop
  999   continue
        end

        subroutine mdsrch(cmin,cmax)
c-----
c       Examine the velocity model. Order the velocity model in terms
c       of increasing velocity. This is done twice, once for S only
c       and then for S and P combined. The objective is to determine
c       regions where a denser search in phase velocity should
c       be made to avoid missing modes
c-----
        parameter (LIN=5, LOT=6, NL=1000, NL2=NL+NL)
        implicit double precision (a-h,o-z)
        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
        real*4 vth(NL2),vtp(NL2)
        real*4 vts(NL2,2)
        equivalence(vth(1),vts(1,1)),(vtp(1),vts(1,2))
        common/pari/ mmax,mode
        common/water/iwat(NL)
        common/vels/ mvts(2),vts
        real*4 cmin, betmn, betmx, cmax
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid
c-----
c       store the layer velocity to be used to initiate 
c       denser pole searching.
c-----
            m=0
            do 310 i=1,mmax
                if(sqrt(TL(i)/Trho(i)).gt.0.000) then
                    m=m+1
                    vth(m)=sqrt(TL(i)/Trho(i))
                endif
  310       continue
            call bsort(m,vth,1)
            m=m+1
            vth(m)=1.0e+30
            mvts(1)=m
            m=0
            do 320 i=1,mmax
                if(sqrt(TA(i)/TRho(i)).gt.0.0) then
                    m=m+1
                    vtp(m)=sqrt(TA(i)/TRho(i))
                endif
                if(sqrt(TL(i)/Trho(i)).gt.0.000) then
                    m=m+1
                    vtp(m)=sqrt(TL(i)/Trho(i))
                endif
  320       continue
            call bsort(m,vtp,1)
            m=m+1
            vtp(m)=1.0e+30
            mvts(2)=m
c-----
c       now establish minimum phase velocity for starting off the
c       solution.
c       This code is from srfdis(IV) and considers whether a water layer
c       is or is not present. If the P velocity in the water layer
c       is less than the S velocity in any other layer, this is the
c       starting velocity.
c       Otherwise, the least shear wave velocity is used. However the
c       starting value is obtained from the Rayleight wave velocity
c       in a halfspace, by an interative solution
c-----
        betmn = 1.0e+20
        betmx = 0.0
        jmn = 1
        jsol = 1
        do 20 i=1,mmax
            if(sqrt(TL(i)/Trho(i)).gt.0.010 .and. 
     1               sqrt(TL(i)/Trho(i)).lt.betmn)then
                betmn = sqrt(TL(i)/Trho(i))
                jmn = i
                jsol = 1
            else if(sqrt(TL(i)/Trho(i)).le.0.01 .and. 
     1              sqrt(TA(i)/Trho(i)).lt.betmn)then
                betmn = sqrt(TA(i)/Trho(i))
                jmn = i
                jsol = 0
            endif
            if(allfluid)then
                if(sqrt(TA(i)/Trho(i)).gt.betmx) 
     1              betmx=sqrt(TA(i)/Trho(i))
            else
                if(sqrt(TL(i)/Trho(i)).gt.betmx) 
     1              betmx=sqrt(TL(i)/Trho(i))
            endif
   20   continue
c-----
c       get starting value for phase velocity, which will 
c       correspond to the VP/VS ratio
c-----
        if(jsol.eq.0)then
c-----
c       water layer
c-----
            cmin = betmn
            cmin=.90*cmin
            cmax = betmx
        else
c-----
c       solid layer solve halfspace period equation
c-----
            cmin = betmn
C            call gtsolh(jmn,cmin)
            cmin=.91*cmin
            cmax = betmx
        endif
        WRITE(6,*)'Searching for ',cmin, ' <= c <= ',cmax
c-----
c       back off a bit to get a starting value at a lower phase velocity
c-----
        return
        end

c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        function dltar(wvno,omega,kk)
c-----
c       control the way to P-SV or SH.
c-----
        implicit double precision (a-h,o-z)
        if(kk.eq.1)then
c-----
c           love wave period equation
c-----
            dltar = dltar1(wvno,omega)
        else if(kk.eq.2)then
c-----
c           rayleigh wave period equation
c-----
            dltar = dltar4(wvno,omega)
C        WRITE(6,*)6.2831853/omega, omega/wvno,dltar
        endif
        end

        function dltar1(wvno,omega)
c-----
c       find SH dispersion values.
c-----
        parameter (NL=1000)
        implicit double precision (a-h,o-z)

        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
        common/pari/ mmax,modew
        common/pard/ twopi,displ,dispr
        common/water/iwat(NL)

        complex*16 esh(2,2), einvsh(2,2)
        double precision hl(2,2)
        logical lshimag
        double precision rsh
        complex*16 e10, e20

c-----
c       Haskell-Thompson love wave formulation from halfspace
c       to surface.
c-----
C        WRITE(6,*)mmax,modew,twopi,displ,dispr

c-----
c       get the eigenfunctions for the halfspace
c-----

		call gtegsh(mmax,wvno, omega, rsh, lshimag)
        call gtesh(esh,einvsh,rsh,wvno,dble(TL(mmax)),lshimag,
     1        .false.,iwat(mmax))

        e1=dreal(einvsh(1,1))
        e2=dreal(einvsh(1,2))
        mmm1 = mmax - 1
        exe = 0.0
        do 600 m=mmm1,1,-1
            if(iwat(m).eq.0)then
               call gtegsh(m,wvno, omega, rsh, lshimag)
               call gtesh(esh,einvsh,rsh,wvno,dble(TL(m)),lshimag,
     1         .false.,iwat(m))
               call varsh(dble(d(m)),rsh,lshimag,
     1            cossh,rsinsh,sinshr,exm)
               call hskl(cossh,rsinsh,sinshr,dble(TL(m)),
     1            iwat(m),hl,exm,exe)
                e10=e1*hl(1,1)+e2*hl(2,1)
                e20=e1*hl(1,2)+e2*hl(2,2)
                xnor=dabs(dreal(e10))
                ynor=dabs(dreal(e20))
                if(ynor.gt.xnor) xnor=ynor
                if(xnor.lt.1.d-40) xnor=1.0d+00
                e1=dreal(e10)/xnor
                e2=dreal(e20)/xnor
            endif
  600   continue
        dltar1=e1
        return
        end

c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        subroutine bsort(nn,x,isign)
c       do bubble sort.
c       isign=+1  increase   =-1 decrease.
c
          real*4 x(nn)
c
        n = nn
        m=0
        do 50 i=2,n
          ii=i-m-1
          do 40 j=ii,1,-1
            x0=x(j+1)-x(j)
              if(isign.le.0) x0=-x0
            if(abs(x0).lt.1.e-7) go to 20
C              if(x0) 10,20,50
            if(x0 .lt. 0.0)then
                go to 10
            else if(x0 .eq. 0.0)then
                go to 20
            else
                go to 50
            endif
10        continue
            x0=x(j)
            x(j)=x(j+1)
            x(j+1)=x0
            go to 40
20        continue
            m=m+1
            do 30 k=j,n-m
              x(k)=x(k+1)
30        continue
            go to 50
40      continue
50    continue
        nn=n-m
        return
        end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        subroutine disprs(ilvry,dt,nn,iret,mname,
     1          verby,nfval,fval,ccmin,ccmax)
c-----
c       search poles in either C-T or F-K domains.
c       To generate synthetic seismograms, F-K domain searching is
c       recommended.
c
c-----
        implicit double precision (a-h,o-z)
        integer NL, NL2
        parameter (LIN=5, LOT=6, NL=1000,NL2=NL+NL)
        integer*4 mmaxx,kmaxx
        integer*4 ifuncc(2)
        real*4 vts(NL2,2)
        real*4 dt, cmin, cmax
        real*4 ccmin,ccmax
        logical verby
        parameter(NPERIOD=2049)
        real fval(NPERIOD)

        parameter(MAXMOD=2000)
        common/phas/ cp(MAXMOD)
        common/pari/ mmax,mode
        common/water/iwat(NL)
        common/pard/ twopi,displ,dispr
        common/stor/ dcc(MAXMOD),c(MAXMOD),nlost,index,nroot1
        common/vels/ mvts(2),vts


        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep
        common/timodel/od(NL),oTA(NL),oTC(NL),oTL(NL),oTN(NL),oTF(NL),
     1      oTRho(NL),
     2      oqa(NL),oqb(NL),oetap(NL),oetas(NL), 
     3      ofrefp(NL), ofrefs(NL)
        real od,oTA,oTC,oTN,oTL,oTF,oTRho,oqa,oqb,
     1        oetap,oetas,ofrefp,ofrefs


        character mname*80
        integer ipar(20)
        real*4 fpar(20)
        character title*80
        integer iunit, idimen,iiso,iflsph,ierr
c-----
c       initialize
c-----
        data ipar/20*0/
        data fpar/20*0.0/
c
   11 format(1x,' improper initial value  no zero found. cmin=',f6.2)
   10 format(' ' ,5x,7f10.4)
   20 format(' ' ,15x,6f10.4)
   30 format(' ' ,16x,'FLAT CRUSTAL MODEL USED  ' )
   40 format(' ' ,7x,'   THICK     TA        TC        TF',
     1    '        TL        TN      DENSITY')
c-----
c       get the earth model
c       This may seem redundant, but if we apply earth flattening to
c       a spherical model, then Love and Rayleigh flattened models
c       are different
c-----
c-----
c       get the earth model
c-----
        call getmod(1,mname,mmmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
        mmax = mmmax
        nsph = iflsph
        ipar(1) = nsph
        fpar(1) = refdep
        fpar(2) = refdep
c-----
c       set fpar(1) = refdep in flat earth model
c   ?   set fpar(2) = refdep in original flat or spherical model
c   ?       fpar(1) may be different because of flattening
c-----
c       for TI we have
c-----
c       set ipar(1)   = 1 if medium is spherical
c       set ipar(2)   = 1 if source is in fluid
c       set ipar(3)   = 1 if receiver is in fluid
c       set ipar(4)   = 1 if eigenfunctions are output with -DER flag
c       set ipar(5)   = 1 if dc/dh are output with -DER flag
c       set ipar(6)   = 1 if dc/dav are output with -DER flag
c       set ipar(7)   = 1 if dc/dbv are output with -DER flag
c       set ipar(8)   = 1 if dc/dr are output with -DER flag
c       set ipar(9)   = 1 if dc/dah are output with -DER flag
c       set ipar(10)  = 1 if dc/dn are output with -DER flag
c       set ipar(11)  = 1 if dc/dbh are output with -DER flag
c-----
        do 339 i=1,mmax
            d(i) = od(i)
            TA(i) = oTA(i)
            TC(i) = oTC(i)
            TF(i) = oTF(i)
            TN(i) = oTN(i)
            TL(i) = oTL(i)
            Trho(i) = oTrho(i)
            if(TN(i).le.1.0e-4*TA(i))then
                iwat(i) = 1
            else
                iwat(i) = 0
            endif
C        WRITE(6,*)d(i),TA(i),TC(i),TF(i),TN(i),TL(i),Trho(i),iwat(i)
  339   continue
c-----
c       d = thickness of layer in kilometers
c       a = compressional wave velocity in km/sec
c       b = transverse wave velocity in km/sec
c       rho = density in gm/cc
c-----

c-----
c       if the model is spherical, then perform a wave type dependent
c       earth flattening
c-----
        if(nsph.gt.0)then
            call sphere(ilvry)
        endif
            write(LOT,30)
            write(LOT,40)
            if(nsph.eq.0)then
                  WRITE(LOT,*)'Velocity model is flat'
            else
                  WRITE(LOT,*)'Velocity model is spherical. ',
     1                'The following is the flattened model'
            endif
            do 341 i=1,mmax-1
                write(LOT,10) d(i),TA(i),TC(i),
     1           TF(i),TL(i),TF(i),Trho(i)
  341       continue
            write(LOT,20) TA(mmax),TC(mmax),TF(mmax),TL(mmax),
     1         TN(mmax),TRho(mmax)
c-----
c       look at the model, determined the range of velocities
c       for the SH or P-SV problems
c-----
        call mdsrch(cmin,cmax)
        if(ccmin.gt.0.0)cmin = ccmin
        if(ccmax.gt.0.0)cmax = ccmax
        write(LOT,*)'cmin,cmax',cmin,cmax
c-----
        if(ilvry.eq.1)then
            write(LOT,*)'Love Computations'
        else if(ilvry.eq.2)then
            write(LOT,*)'Rayleigh Computations'
        endif
        write(LOT,*) ' Number of layers=',mmax,' Maximum modes=',mode
c-----
c       read in control parameters.
c-----
        iret = 1
        write(LOT,*) ' NPTS=',nn,' DT=',dt,' CMIN=',cmin, ' nfval=',
     1     nfval
        df=1./(nn*dt)
        if(nfval.gt.0)then
            kmax = nfval
        else
            kmax=nn/2+1
        endif
        mmaxx=mmax
        kmaxx=kmax
c-----
c       I/O file set up.
c-----
        if(ilvry.eq.1 )then
c-----
c       Love wave
c-----
            open(1,file='tdisp96.lov',status='unknown',
     1          form='unformatted', access='sequential')
            disp = displ
        else
c-----
c       Rayleigh Wave
c-----
            open(1,file='tdisp96.ray',status='unknown',
     1          form='unformatted', access='sequential')
            disp = dispr
        endif
        rewind 1
        call ptsmdt(1,mmaxx,d,ta,tc,tf,tl,tn,Trho,qa,qb,kmaxx,
     1      mname,ipar,fpar)
c-----
c       The main processing section
c-----
        lyr0=0
        do 410 i=1,mvts(ilvry)
                if(cmin.gt.vts(i,ilvry)) lyr0=i
  410   continue
        if(lyr0.eq.0) lyr0=1
        k1=1
        kmode=mode
        c1=cmin
        ndisp=dabs(disp)
        if(ndisp.le.0) ndisp=1
        ifunc=ilvry
c---------------------------------------
c       search poles in F-K domain.
c---------------------------------------
        nlost = 0
        index = 0
        do 600 i=1,MAXMOD
            dcc(i)=0.0
            c(i)=0.0
  600   continue
        nroot1=MAXMOD
        do 2800 k = k1,kmax
c---------------------------------------
c           search poles in F-K domain.
c---------------------------------------
            if(nfval.le.0)then
                t1=dfloat(kmax-k)*df
                if(k.eq.kmax) t1=0.01*df
                t11=dfloat(kmax-k+1)*df
                index=index+1
            else
c-----
c       force a complete search for user provided frequencies
c       do not attempt to follow a curve
c-----
                t1 = fval(k)
                t11 = 0.99*fval(k)
                index = 0
            endif
            omega1=twopi*t1
            omega0=twopi*t11
c-----
            kmode=0
            call poles(ifunc,omega0,omega1,cmin,cmax)
            if(verby)then
                if(k.eq.1)then
                    write(LOT,4443)k,t1,nroot1
                else
                    write(LOT,4444)k,t1,nroot1
                endif
            endif
 4443   format('  Frequency(',i5,')=',1pe11.4,' Modes=',i5)
 4444   format('+ Frequency(',i5,')=',1pe11.4,' Modes=',i5)
            kmode=nroot1
c-----
            if(k.eq.1.and.kmode.eq.0) then
                write(LOT,11) cmin
                iret = 0
            endif
            if(k.eq.1.and.kmode.eq.0) go to 3000
            do 2100 i=1,kmode
                cp(i)=omega1/dcc(i)
 2100       continue
            t0=1./t1
c-----
c       print out possible dc values.
c-----
            if(k.eq.1) then
                lmode=kmode
                if(lmode.le.1) go to 2700
                dcmax=-1.e+30
                dcmin=1.e+30
                dcsum=0.0
                do 2600 i=1,lmode-1
                    dc=cp(i+1)-cp(i)
                    if(dcmax.le.dc) dcmax=dc
                    if(dcmin.ge.dc) dcmin=dc
                    dcsum=dcsum+dc
 2600           continue
                dcsum=dcsum/float(lmode-1)
                write(LOT,*) ' '
        write(LOT,*) '============================================='
                if(ilvry.eq.1)
     1  write(LOT,*) 'dc value for LOVE at the lowest period=',t0
                if(ilvry.eq.2)
     1  write(LOT,*) 'dc value for RAYLEIGH at the lowest period=',t0
        write(LOT,*) ' between c(1)=',cp(1),' and'
        write(LOT,*) '         c(',lmode,')=',  cp(lmode),' :'
        write(LOT,*) ' mean=',dcsum
        write(LOT,*) ' max =',dcmax
        write(LOT,*) ' min =',dcmin
        write(LOT,*) ' '
            endif
 2700           continue
c-----
c       Output.
c-----
            ifuncc(ilvry)=ifunc
C       write(LOT,*)'DISPRS: lmode, kmode, kmodee',
C     1     lmode, kmode, kmodee
            kmodee=kmode
        call ptshed(1,ifuncc(ilvry),kmodee,t0)
C       write(LOT,*) ifuncc(ilvry),kmodee,t0
            if(kmode.gt.0) then
            call ptsval(1,cp,kmode)
C       write(LOT,*) (cp(i),i=1,kmode)
            endif
 2800   continue
        kmodee=-1
        call ptshed(1,kmodee,kmodee,t0)
 3000   continue
        close(1,status='keep')
        return
        end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        subroutine poles(ifunc,omega0,omega,cmin,cmax)
c*****  This is the most important subroutine in this program.
c       search pole on wavenumber-frequency domain.
c       Each pole is followed by jumping method.
c       At first two periods, poles are located very precisely.
c
c       However, there doesn't exist a perfect scheme for such 
c       a complicated problem. The program 'scomb96.f' is 
c       designed to repair this unavoidable fault.
c
c       [eroot]  store the wavenumbers for angular freq [omega0]
c       [cphase] store the phase velocies for [omega0+domega]
c       and they will be replaced by the i'th mode's wavenumber 
c       (at omega) and i'th mode's phase velocity (at omega0)
c-----
        implicit double precision (a-h,o-z)
        integer LIN, LOT, NL, NL2
        parameter (LIN=5, LOT=6,NL=1000,NL2=NL+NL)
        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
        common/pari/ mmax,mode
        common/pard/ twopi,displ,dispr
        common/water/iwat(NL)
        real*4 vts(NL2,2)
        dimension wvm(NL2)
        integer MAXMOD
        parameter(MAXMOD=2000)
        common/stor/ eroot(MAXMOD),cphase(MAXMOD),nlost,index,nroot1
        common/vels/ mvts(2),vts
        real*4 cmin, cmax
        integer*4 nmx,mtest,i,j,k,ii,jj,kk
c-----
c       this routine finds the roots of period equation using
c       regular halving method to initiate the first two
c       frequencies, then followed by jumping method.
c-----
        epi  = 1.d-10
        freq = omega/twopi
        wvmx = omega/dble(cmin)
        wvmn = omega/dble(cmax)
        wvmn1= omega/(dble(cmax)+0.03)
        tp=twopi/omega
c-----
c       nmx evenly divides the real k axis between Kmin and Kmax.
c       To begin we need an initial estimate of the number of modes
c       at the frequency freq. 
c
c       For Love waves in a single layer over a halfspace
c
c       n = omega H sqrt( 1/b(1)**2 - 1/b(2)**2 ) / pi
c-----
        sum = 0.0
        do 30 i=1,mmax
            if(TL(i).eq.0.0)then
                v = sqrt(TA(i)/TRho(i))
                fac = 1.0/(v**2) - 1.0/(cmax**2)
            else
                v = sqrt(TL(i)/TRho(i))
                fac = 1.0/(v**2) - 1.0/(cmax**2)
            endif
            if(fac.lt.0.0)fac = 0.0
            sum = sum + d(i)* sqrt(fac)
   30   continue
        sum = 2.0*freq*sum
        nmx = 200 + 10.0*sum
        ra=10.0
        if(TL(mmax-1).ne.0.0) ra=sqrt(TL(mmax)/TL(mmax-1))
        if(ra.gt.3.0) ra=3.0
        if(ra.gt.1.0) nmx=nmx*ra
    
        disp=displ
        if(ifunc.eq.2) disp=dispr
        ndisp=dabs(disp)
        if(ndisp.eq.0) ndisp=1
        nmx=nmx*ndisp
        rangeK = wvmx-wvmn 
        dk = rangeK/nmx
        do 80 i=1,mvts(ifunc)
            wvm(i)=omega/dble(vts(i,ifunc))
   80   continue
c-----
c       This is the most important routine
c       We know that the phase velocity dispersion is bounded by
c       the maximum shear wave velocity and the
c       minimum of P wave velocity (water layer, shear wave velocity,
c           or surface Rayleigh wave velocity
c
c       The number of modes increases with frequency in the manner
c       of Love waves, e.g.,
c           omega H sqrt ( 1. - B1*B1/B2*B2 )/B1 = N pi
c           (this is just the SH higher mode cutoff for a
c           single layer over a halfspace)
c
c       We also know that the higher modes are very dense at fixed
c       frequency, just above a layer velocity. The array vts
c       is a list of all possible layer velocities in the 
c       phase velocity limits between cmin and bmax.
c
c       We initially estimate the density of roots and bracket
c       zeros in the period equation, with further refinement
c       by interval halving (slow but always works) or linear
c       interpolation. Linear interpolation has problems since the
c       behavior of the period equation is not linear near a root.
c
c       In an extreme case we may have in the region of a root
c                   ^
c                   |
c       DEL         X
c        |
c        |
c        |
c        -----------|---|---|---------
c                       x   x
c
c                   c1  c3  c4
c       linear interpolation between c1 and c4 would not converge fast
c       since the next estimate would be near c4
c
c       linear interpolation between c3 and c1 would also not work.
c       However, if c3 and c4 are not identical, a root estimate
c       may be able to be made by linear interpolation on c3 and c4
c       if the result lies between c1 and c4. If this is the case,
c       we use a modified interval halving.
c
c       The other aspect concerns estimating the roots at the next
c       frequency. Assuming that the phase velocity dispersion
c       always decreases with frequency, 
c-----
        if(nlost.eq.1001) go to 2000
        if(index.gt.2) go to 3000
 2000   continue
c*****  The regular halving method.
c-----
c       At the very first two periods or the places jumping method fails
c       we find the poles using the regular halving method.
c       nmx is chosen for a 40 km crustal model. for shallower
c       thickness a proportionately smaller nmx can be used
c       search for roots of period equation
c-----
        nroot = 0
        c2 = wvmx
        del2 = dltar(c2,omega,ifunc)
        lyr=1
        jj=1
        do 500 i=1,nmx
            jj=jj-1
            if(jj.ne.0) go to 500
            c10 = wvmx-float(i)*dk
            if(i.eq.nmx) c10=wvmn+0.01*dk
            jj = 1
            kk = 1
            if(c10.gt.wvm(lyr)) go to 90
c-----
c           kk and jj represent a denser searching when phase velocity
c           = S wave velocity. Their values can be changed as kk=3*lyr
c           jj=8.
c-----
            if(vts(lyr,ifunc).ne.0.0) 
     1          ra=vts(lyr+1,ifunc)/vts(lyr,ifunc)
            if(ra.gt.3.0) ra=3.01
            if(ra.lt.1.0) ra=1.01
            nra=ra
            kk = 10.0*ra
            kk = kk*ndisp
            lyr = lyr+1
            jj = 4*nra+ndisp
   90       continue
            dk0 = dk/float(kk)
            do 400 j=1,jj
            do 400 k=1,kk
                if(nroot.eq.mode) go to 510
                jk = kk*(j-1)+k
                c1 = c10+dk-dk0*float(jk)
                if(c1.lt.wvmn1) go to 510
                del1 = dltar(c1,omega,ifunc)
                if(dsign(1.0d+00,del1)*
     1              dsign(1.0d+00,del2).ge.0) go to 350
                nroot = nroot + 1
                c4 = c2
                del4 = del2
c-----
c       perform interval halving, result is guaranteed to be within
c       be 2**(-15) abs(c4-c1) of the true zero
c-----
                do 200 ii=1,15
                    c3 = 0.5*(c1+c4)
                    del3 = dltar(c3,omega,ifunc)
                    if(dsign(1.0d+00,del1)
     1                  *dsign(1.0d+00,del3).ge.0)then
                        del1 = del3
                        c1 = c3
                    else 
                        del4 = del3
                        c4 = c3
                    endif
                    if(dabs(c4-c1).lt.epi*c1) go to 250
  200           continue
  250           continue
                c3 = 0.5*(c1+c4)
                if(index.eq.1) cphase(nroot)=omega/c3
                if(nlost.eq.1) cphase(nroot)=omega/c3
                eroot(nroot) = c3
  350           c2 = c1
                del2 = del1
  400       continue
  500   continue
  510   continue
        nlost=nlost+1000
        go to 1250
 3000   continue
c-----
c***    The jumping method.
c       jump(del k)= -(del omega)/c - (del c)*wvno/c
c                         t                tadd
c       if this method fails, the control will return
c       to the regular halving method.
c       This method has been largely improved in order to handle
c       locked mode approximation. Much denser searching is enforced
c       when jumping over c=v(layer).
c-----
        nroot   = 0
        if(nroot1.eq.0) go to 1250
        dk0 = 0.1*dk
        domega  = omega0-omega
        eroot(nroot1+1)=wvmn
        if(nroot1.eq.mode) eroot(nroot1+1)=eroot(nroot1)-5.*dk
        nlost   = 0

        i = 0
 1200   continue
        i = i + 1
        if(i.gt.nroot1)go to 1201
        iNEG = 0

        if (i.ne.1) then
            Cw00 = cphase(i-1)
            Cw01 = omega0/eroot(i)
            Cw02 = omega0/eroot(i+1)
            dCint= Cw01-Cw00
            if ((i+1).le.nroot1) dCint1=Cw02-Cw01
            if ((i+2).le.nroot1) dCint2=omega0/eroot(i+2)-Cw02
        endif

        s   = eroot(i)
        vphase  = omega0/s
        t   = -domega/vphase
        dc  = vphase-cphase(i)
        dci1    = 2.0*omega0/eroot(i+1) - cphase(i+1)
        dci0    = 2.0*omega0/eroot(i) - cphase(i)
        Climit  = 0.0
        if (i.lt.(nroot1-2).and.index.gt.3) then
            Climit = omega/omega0*eroot(i+2)
        endif
        if (i.gt.5 .and. i.lt.(nroot1-2) .and. dc.gt.dCint1) dc=dCint1
        if (i.ge.(nroot1-2)) dc=dc*0.3
        tadd    = -dc*s/vphase

c       write(LOT,*)'i=',i
c       if (i.eq.1) then  
c       write(LOT,*)' '
c       write(LOT,*)'wvmx=',wvmx,' wvmn=',wvmn
c       write(LOT,*)'dc=',dc,' dCint1=',dCint1
c       write(LOT,*)'vphase=',vphase,' cphase=',cphase(i)
c       endif

        dc  = s-eroot(i+1)

        if(i.eq.nroot1.and.nroot1.gt.1) dc=dc0
CRBH    cs = dabs(tadd)/tadd
CRBH    never divide by zero, just check the sign
        if(tadd .gt. 0.0)then
            cs = 1.0
        else if(tadd.lt.0.0)then
            cs = -1.0
        else
            cs = 0.0
        endif
        if(i.eq.1) then
            icare=0
            if(dabs(tadd).gt.dabs(0.05*rangeK)) 
     1          tadd=-dabs(0.05*rangeK)
            if(dabs(t).gt.dabs(0.1*rangeK)) icare=1 
            if(dabs(t).gt.dabs(0.3*rangeK)) icare=2 
            go to 520
        endif
        if(dabs(tadd).lt.dk0) tadd=cs*dk0
        if(dabs(tadd).gt.0.8*dc) tadd=cs*0.5*dc
  520   continue
        dc0=dc
c-----
c       s+t makes the first jump.
c-----
        c2  = s+t
        ncont   = 0
        ifirst  = 0
c       write(LOT,*)' '
c       write(LOT,*)'tp=',tp,'omega=',omega,' tadd=',tadd
c       write(LOT,*)'dk=',dk,' s=',s,' t=',t 
c-----
c       here a backward jump done for the fundamental mode.
c-----
        if(i.eq.1) then
            c3=s-2.*t
            if(c3.gt.wvmx) then
                c3=s-t
                if(c3.gt.wvmx) c3=s-0.5*t
                if(c3.gt.wvmx) c3=s-0.25*t
                if(c3.gt.wvmx) then
                    c3=wvmx*0.9999
                    if (s.gt.wvmx) then
                        if (c2.lt.wvmn) then
                            c4=0.5*(c3+c2)
                        else
                            c4=0.5*(c3+wvmn)
                        endif
                    endif
                endif
            endif
            c4=s
            if (icare.eq.1) then
                c3=c4+0.9*t
                c4=c4+0.95*t
            else if (icare.eq.2) then
                c3=c2+0.01*rangeK
                c4=c2+0.02*rangeK
            endif
            del3=dltar(c3,omega,ifunc)
            del2=dltar(c2,omega,ifunc)
            del4=dltar(c4,omega,ifunc)
            if(dsign(1.0d+00,del2)*
     1          dsign(1.0d+00,del4).le.0) then
                c1=c2
                del1=del2
                c2=c4
                del2=del4
                ihalf=20
                go to 850
            endif
            if(dsign(1.0d+00,del3)*
     1          dsign(1.0d+00,del4).le.0) then
                c1=c4
                del1=del4
                c2=c3
                del2=del3
                ihalf=20
                go to 850
            endif
            if (abs(del2).le.abs(del4) .and. abs(del4).le.abs(del3)
     *          .and. tadd.gt.0.0) tadd=-tadd
        endif
c-----
c       a protection for the first wild jump.
c-----
        nctrl=25
        if(i.ne.1) then
            c3=eroot(i-1)-0.05*dk0
            c20=c2
                c40=s+0.8*t
c       write(LOT,*)'eroot=',eroot(i-1),'  dk0=',dk0
            if(c20.gt.c40) then
                cx=c40
                c40=c20
                c20=cx
            endif
            if(c20.lt.c3.and.c3.le.c40) then
                del2=dltar(c20,omega,ifunc)
                del3=dltar(c3,omega,ifunc)
c       write(LOT,*)'c20=',c20,'  del2=',del2
c       write(LOT,*)'c3 =',c3 ,'  del3=',del3
                if(dsign(1.0d+00,del2)*
     1          dsign(1.0d+00,del3).le.0) then
                    c1=c20
                    del1=del2
                    c2=c3
                    del2=del3
                    ihalf=20
                    go to 850
                endif
            endif
            if(c3.gt.c40) then
                del2=dltar(c20,omega,ifunc)
                del3=dltar(c3,omega,ifunc)
                del4=dltar(c40,omega,ifunc)
c       write(LOT,*)'c20=',c20,'  del2=',del2
c       write(LOT,*)'c40=',c40,'  del4=',del4
c       write(LOT,*)'c3 =',c3 ,'  del3=',del3
                if(dsign(1.0d+00,del3)
     1              *dsign(1.0d+00,del4).le.0) then
                    c1=c40
                    del1=del4
                    c2=c3
                    del2=del3
                    ihalf=20
                    go to 850
                endif
                if(dsign(1.0d+00,del2)
     1              *dsign(1.0d+00,del4).le.0) then
                    c1=c20
                    del1=del2
                    c2=c40
                    del2=del4
                    ihalf=20
                    go to 850
                endif
            endif
            if(c3.lt.c20) then
                c2=c3
            endif

            if (i.le.4 .and. i.ge.(nroot1-2)) goto 528
            fact=1.0
            dCnow=omega/(c2+tadd)-omega/c2
            ddci=dci1-dci0
            if (ddci.le.0 .or. ddci.le.(0.3*dCint1))then 
                Fdense=1
            else if (ddci.le.0 .or. ddci.le.(0.1*dCint1))then 
                Fdense=2
            else
                Fdense=0
            endif
c       write(LOT,*)'Fdense=',Fdense

            if (dCint1 .lt. (0.6*dCint)) then
                fact=3./(3.+(dCint/dCint1)**2)
                if (ifunc.eq.2) fact=0.6*fact
                if (fact.lt.0.25) fact=0.25
                if (Fdense.eq.1) fact=fact*0.5
                if (Fdense.eq.2) fact=fact*0.25
                tadd9=omega/(omega/c2+fact*dCint1)-c2
                if (abs(tadd9).lt.abs(tadd)) tadd=tadd9
                nctrl=50
            else if (dCint1.gt.(1.5*dCint) .and. 
     *          dCint1.gt.(1.5*dCint2)) then
                fact=2.0/(2.0+dCint/dCint1)
                if (Fdense.eq.1) fact=0.2
                if (Fdense.eq.2) fact=0.1
                tadd=tadd*fact
                if (ifunc.eq.2) tadd=tadd*0.5
                nctrl=50
            else if (dCnow.gt.0.5*dCint1) then
                fact=0.5
                        if (ifunc.eq.2) fact=0.4
                if (Fdense.eq.1) fact=0.2
                if (Fdense.eq.2) fact=0.1
                        tadd9=omega/(omega/c2+fact*dCint1)-c2
                        if (abs(tadd9).lt.abs(tadd)) tadd=tadd9
                        nctrl=50
            else if (Fdense.ne.1) then
                if (Fdense.eq.1) tadd=tadd*0.5
                if (Fdense.eq.2) tadd=tadd*0.2
                nctrl=50
            endif
        endif
  528   c22=c2
        cKEP2=c2
        iAGAIN=0
  529   if (iAGAIN.eq.1) then
            tadd=tadd*0.3
            c2=cKEP2
            nctrl=50
c       write(LOT,*)'iAGAIN=',iAGAIN,'  tadd=',tadd 
        else if (iAGAIN.eq.2) then
            tadd=tadd*0.2
            c2=cKEP2
            nctrl=50
c       write(LOT,*)'iAGAIN=',iAGAIN,'  tadd=',tadd 
            if (mtest.gt.20) nctrl=100
            iAGAIN=3
        endif
c-----
c       ntest,nctrl control the number of jumps using tadd.
c       ihalf is the times to do a pole refinement.
c       ncont controls which pole searching mechanism being used.
c-----
        ihalf  = 20
        ncont  = 0
  550   if(c2.le.wvmn1) go to 1250
        del2 = dltar(c2,omega,ifunc)
        ntest  = 0
        mtest  = 0
        tadd1=tadd
c-----
c       mtest,tadd0 control a denser pole searching around c=v(layer).
c       jmp denotes the wild jump (c2=s+t) 
c       just jumping over a c=v(layer).
c-----
  600   continue
        mtest=mtest+1
        tadd0=tadd
        c1=c2+tadd0
        if(c1.le.wvmn1) go to 800
c-----
c       Here add another constrant for missing a small 2_root gap, and
c       reach an unreasonable point
c-----
        if (c1.le.Climit .and. index.gt.3) then
            if (iAGAIN.eq.0) then
                iAGAIN=1
                goto 529
            else if (iAGAIN.eq.1) then
                iAGAIN=2
                goto 529
            endif
        endif
c-----
c       here to catch a reverse dispersion.
c-----
        if(i.ne.1 .and. c1.ge.eroot(i-1)) then
            c2=c22
            tadd=-dabs(tadd)
c       write(LOT,*)'catch a reverse dispersion  tadd=',tadd
            if(c2.eq.c22) go to 550
        endif
        del1    = dltar(c1,omega,ifunc)
c       write(LOT,*)'c1=',c1,' del1=',del1
        if(dsign(1.0d+00,del1)*dsign(1.0d+00,del2).le.0) go to 850
c-----
c       this part do some convergence direction check and try to
c       fix some smaller root_interval
c       the normal convergence is from (-) to (+) or from (+) to (-)
c       if (-) --> (-0) --> (-)   or   (+) --> (+0) --> (+) happen,
c       then we can sure that we miss a very small 2 root_interval gap
c-----
        if (i.ge.nroot1) goto 630 
        if (mtest.eq.1) then
            ckeep1=c1
            dkeep1=del1
            ckeep2=c2
            dkeep2=del2
            jjdir0=0
            if (del1.ge.0.9999999) then
                jjdir0=1
                dire0 =-1.0
            else if (del1.le.-0.9999999) then
                jjdir0=1
                dire0 =+1.0
            endif
            if (abs(del1).lt.abs(del2)) then
                jjdir0=1
                dire0 =del1-del2
            endif
        else if (mtest.ge.2) then
            if((c1+tadd0).le.wvmn1) go to 800
            if (jjdir0.eq.0) then
                if (abs(del1).lt.abs(dkeep1)) then
                    dire0=del1-dkeep1
                    jjdir0=1
                endif
            else
                jjdir1=1
            endif
            dire1=del1-dkeep1

            if ((jjdir0*jjdir1.ne.0).and.(dire1*dire0.lt.0)) then
                iNEG=0
                xLOW=abs(dkeep1)
                ifound=0
                ckeep0=c1
                dkeep0=del1

                do 626 kk=1,10
                    cU1=(ckeep0*2+ckeep1)/3
                    cU2=(ckeep0+ckeep1*2)/3
                    dU1=dltar(cU1,omega,ifunc)
                    dU2=dltar(cU2,omega,ifunc)

                    cD1=(ckeep1*2+ckeep2)/3
                    cD2=(ckeep1+ckeep2*2)/3
                    dD1=dltar(cD1,omega,ifunc)
                    dD2=dltar(cD2,omega,ifunc)

                    direU1=dkeep0-dU1
                    direU2=dU1-dU2
                    direD1=dD1-dD2
                    direD2=dD2-dkeep2

                    if ((dU1*dkeep1).le.0) then
                        iNEG=1
                        cc21=ckeep0
                        cc22=cU1
                        dd21=dkeep0
                        dd22=dU1
                        if((dU2*dkeep1).gt.0) then
                            cc11=cU1
                            cc12=cU2
                            dd11=dU1
                            dd12=dU2
                        endif
                    endif
                    if ((dU2*dkeep1).le.0) then
                        if (iNEG.eq.0) then
                            cc21=cU1
                            cc22=cU2
                            dd21=dU1
                            dd22=dU2
                            cc11=cU2
                            cc12=ckeep1
                            dd11=dU2
                            dd12=dkeep1
                        else
                            cc11=cU2
                            cc12=ckeep1
                            dd11=dU2
                            dd12=dkeep1
                        endif
                        iNEG=iNEG+1
                    endif
                    if ((dD1*dkeep1).le.0) then
                        cc21=ckeep1
                        cc22=cD1
                        dd21=dkeep1
                        dd22=dD1
                        if ((dD2*dkeep1).gt.0) then
                            cc11=cD1
                            cc12=cD2
                            dd11=dD1
                            dd12=dD2
                        endif
                        iNEG=1
                    endif
                    if ((dD2*dkeep1).le.0) then
                        if (iNEG.eq.0) then
                            cc21=cD1
                            cc22=cD2
                            dd21=dD1
                            dd22=dD2
                            cc11=cD2
                            cc12=ckeep2
                            dd11=dD2
                            dd12=dkeep2
                        else
                            cc11=cD2
                            cc12=ckeep2
                            dd11=dD2
                            dd12=dkeep2
                        endif
                        iNEG=iNEG+1
                    endif
                    if (iNEG.ne.0) goto 848
                    if (abs(dU1).lt.xLOW) then
                        ckeep2=ckeep1
                        dkeep2=dkeep1
                        ckeep1=(ckeep0+ckeep1)/2
                        dkeep1=dltar(ckeep1,omega,ifunc)
                    else
                        ckeep0=cU1
                        dkeep0=dU1
                        if (abs(dU2).lt.xLOW) then
                            ckeep2=ckeep1
                            dkeep2=dkeep1
                            ckeep1=cU2
                            dkeep1=dU2
                        else
                            ckeep0=cU2
                            dkeep0=dU2
                        endif
                    endif
                    if (abs(dD2).lt.xLOW) then
                        ckeep0=ckeep1
                        dkeep0=dkeep1
                        ckeep1=(ckeep0+ckeep2)/2
                        dkeep1=dltar(ckeep1,omega,ifunc)
                    else
                        if (abs(dD1).lt.xLOW) then
                            ckeep2=cD2
                            dkeep2=dD2
                        else
                            ckeep2=cD1
                            dkeep2=dD1
                        endif
                    endif
  626           continue
            endif

            ckeep1=c1
            dkeep1=del1
            ckeep2=c2
            dkeep2=del2
        endif
c-----
c
c-----
  630   ntest   = ntest+1
        if(ntest.ge. nctrl) go to 650
        del2    = del1
        c2  = c1
        go to 600
  650   ncont   = ncont+1
        go to (700,720,750),ncont
c-----
c       This is another kind of jumping procedure, which is
c       a remedy when the first jumping method fails.
c-----
  700   continue
c       write(LOT,*)'LABEL 700'
        ihalf  = 20*ndisp
        tadd   = -dc/float(ihalf)
        c2     = c22
        nctrl  = 40*ndisp
        if(i.eq.1) nctrl = 20*ndisp
        go to 550
c-----
c       This is the third kind of jumping procedure.
c-----
  720   continue
c       write(LOT,*)'LABEL 720'
        tadd   = t/float(10*ndisp)
        c2     = s+5.d+00*tadd
        if(i.gt.1.and.c2.gt.eroot(i-1)) c2=eroot(i-1)-0.01*dk0
        nctrl  = 15*ndisp
        ihalf  = 20
        go to 550
c-----
c       If all jumping methods fail, it goes back to the regular
c       halving method.
c-----
  750   continue
c       write(LOT,*)'LABEL 750'
        if(i.eq.mode) go to 1201
        nlost=nlost+1
        wvmn2=wvmn+dk
        if(c1.lt.wvmn2.or.c2.lt.wvmn2) go to 1250
        c3=0.0
        if(eroot(nroot).ne.0.0) c3=omega/eroot(nroot)
        write(LOT,*)
     * 'At perd=',twopi/omega,' mode=',nroot,' c=',c3,' JUMP fail.'
        go to 2000
  800   c1 = wvmn+0.01*dk0
        del1 = dltar(c1,omega,ifunc)
        if(dsign(1.0d+00,del1)*dsign(1.0d+00,del2).le.0) go to 850
        go to 1250
  848   iROOT=1
  849   if (iNEG.ne.0) then
            if (iROOT.eq.1) then
                iROOT=2
                c1=cc11
                c2=cc12
                del1=dd11
                del2=dd12
            else
                iNEG=0
                c1=cc21
                c2=cc22
                del1=dd21
                del2=dd22
                i=i+1
                vphase=omega0/eroot(i)
            endif
        endif
  850   c4 = c2
        del4 = del2
        itype=1
        if (del1.eq.0) then
            c4=c1
            goto 1050
        else if (del4.eq.0) then
            c1=c4
            goto 1050
        endif
c-----
c       refine the roots using halving method.
c       09/28/95 MODIFIED_HALVING_METHOD can reduce 40% computing time
c-----
        aa1 = abs(del1)
        aa4 = abs(del4)
        do 1000 ii=1,ihalf
            if (itype.eq.1) then
                factor=abs(del1)/(abs(del1)+abs(del4))
                if (factor.lt.0.05) factor=0.05
                if (factor.gt.0.95) factor=0.95
                c3 = c1*(1.0-factor)+c4*factor
            else if (itype.eq.0) then
                c3 = 0.5*(c1+c4)
            else if (itype.eq.-1) then
                c3 = c333
            endif

            del3 = dltar(c3,omega,ifunc)
            aa3 = abs(del3)
c-----
c       Now consider the possible cases resulting from this computation
c       Note that sign of DEL4 always = - sign DEL1
c-----
c       Case    DEL1    DEL3    DEL4
c       A   <0  <0  >0
c       B   <0  >=0 >0
c       C   >0  <0  <0
c       D   >0  >=0 <0
c-----

c-----
c       define the default rule
c-----
            itype=0
            if (del1.lt.0.0) then
c-----
c               thus del4 >= 0
c-----
                if (del3.lt.0.0) then
c-----
c       case A
c-----
                    if (aa1.ge.del4.and.aa3.ge.del4) then
                        itype=1
                    else if(aa1.ne.aa3)then
                        c333=c1-aa1*(c1-c3)/(aa1-aa3)
                        if (c333.lt.0.5*(c3+c4) .and.
     *                  (c333.gt.c3)) then
                            itype=-1
                        endif
                    endif
                else if (del3.ge.0.0) then
c-----
c       case B
c-----
                    if (aa1.lt.aa4.and.aa1.lt.aa3) then
                        itype=1
                    else if(del4.ne.del3) then
                        c333=c4-del4*(c4-c3)/(del4-del3)
                        if (c333.gt.0.5*(c1+c3) .and.
     *                  (c333.lt.c3)) then
                            itype=-1
                        endif
                    endif
                endif
            else if (del1.ge.0.0) then
c-----
c           del4 < 0
c-----
                if (del3.lt.0.0) then
c-----
c       case C
c-----
                    if (aa4.ge.del1.and.aa3.ge.del1) then
                        itype=1
                    else if(aa4.ne.aa3)then
                        c333=c4-aa4*(c4-c3)/(aa4-aa3)
                        if (c333.gt.0.5*(c1+c3) .and.
     *                  (c333.lt.c3)) then
                            itype=-1
                        endif
                    endif
                else if (del3.ge.0.0) then
c-----
c       case D
c-----
                    if (aa4.lt.del1.and.del3.ge.aa4) then
                        itype=1
                    else if(del1.ne.del3)then
                        c333=c1+del1*(c3-c1)/(del1-del3)
                        if (c333.lt.0.5*(c3+c4) .and.
     *                  (c333.gt.c3)) then
                            itype=-1
                        endif
                    endif
                endif
            endif
            if(dsign(1.0d+00,del1)*dsign(1.0d+00,del3).ge.0.0) then 
                del1 = del3
                c1 = c3
                aa1 = aa3
            else 
                del4 = del3
                c4 = c3
                aa4 = aa3
            endif 
            if(dabs(c4-c1).lt.epi*c4) go to 1050
 1000   continue
 1050   continue
        jjroot=0
 1060   jjroot=jjroot+1
c-----
c       I like the root position close to small wavenumber c1 
c-----
        if ((c4-c1).lt.5.0D-13) then
            c3=c1
c       write(6,*)'root=',omega/c3
            goto 1070
        endif
        if(abs(del1).lt.abs(del4)) then
            rr = abs(del4)/(abs(del1)+abs(del4))
            if (rr.gt.0.95) rr=0.95
            c3 = c4+rr*(c1-c4)
        else
            c3 = 0.5*(c1+c4)
        endif
        del3= dltar(c3,omega,ifunc)
c       write(6,*)'  c1=',c1,' del1=',del1
c       write(6,*)'  c3=',c3,' del3=',del3
c       write(6,*)'  c4=',c4,' del4=',del4
c       write(6,*)'root=',omega/c3
        if ((del3*del1).lt.0.0) then 
            c4=c3
            del4=del3
            if (jjroot.eq.3) then
                c3=c1
            else
                goto 1060
            endif
        endif
 1070   if(c3.lt.wvmn) go to 1201
c-----
c       examine the result if using the second or third method,
c       or when disp .gt. 4.9.
c-----
        if(disp.gt.4.91) ncont=-1
        if(ncont.ne.0.and.i.ne.1.and.ifirst.ne.1) then
            ifirst=1
            c1=c3
            c2=eroot(nroot)
            if(c1.gt.c2) then
                cx=c1
                c1=c2
                c2=cx
            endif
            ihalf=20
            dc=(c2-c1)/10.0
            c10=c1+0.01*dc
            c2=c2-0.01*dc
            del2=dltar(c2,omega,ifunc)
 1100       continue
            c1=c2-dc
            if(c1.lt.c10) go to 1150
            del1=dltar(c1,omega,ifunc)
            if(dsign(1.0d+00,del1)*dsign(1.0d+00,del2).lt.0) go to 850
            c2=c1
            del2=del1
            go to 1100
 1150       continue
        endif
        nroot = nroot+1
        eroot(nroot)  = c3
        cphase(nroot) = vphase
        if (iNEG.ne.0) goto 849
        go to 1200
 1201   continue
 1250   continue
        if(nroot.gt.nroot1) nroot=nroot1
        nroot1=nroot
        return
        end

        subroutine gtsolh(jmn,c)
c-----
c       using the properties of layer jmn, get the 
c       Rayleigh wave phase velocity
c-----
c-----
c       starting solution
c-----
C        real*4 kappa, k2, gk2
C        c = 0.95*b
C        do 100 i=1,5
C            gamma = b/a
C            kappa = c/b
C            k2 = kappa**2
C            gk2 = (gamma*kappa)**2
C            fac1 = sqrt(1.0 - gk2)
C            fac2 = sqrt(1.0 - k2)
C            fr = (2.0 - k2)**2 - 4.0*fac1*fac2
C            frp = -4.0*(2.0-k2) *kappa
C     1          +4.0*fac2*gamma*gamma*kappa/fac1
C     2          +4.0*fac1*kappa/fac2
C            frp = frp/b
C            c = c - fr/frp
C  100   continue
        return
        end

        subroutine sphere(ifunc)
        integer ifunc
c-----
c       This is subroutine tdosph() in VOLVI/src/tspec96
c       modified to work for the Love Rayleigh problem
c
c       Transform spherical earth to flat earth
c
c       Schwab, F. A., and L. Knopoff (1972). 
c           Fast surface wave and free
c       mode computations, in  
c           Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c           B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  44, 45 , 41 pp 112-113
c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c
c       We will treat all as P-SV for the heck of it
c       This requires more work
c-----
c       mmax    I*4 number of layers
c       TA     R   A 
c       TC     R   C 
c       TF     R   F 
c       TL     R   L 
c       TN     R   N 
c                  note  density not required
c       TD     R   layer thickness
c       v() R   array of velocities
c       h() R   array of layer thicknesses
c       ipsvsh  I       1 - get P time
c                       2 - get SV time
c                       3 - get SH time
c       refdep R   Reference depth for the model specification
c
c       Note we need the constants here.  Since the velocities
c       must increase with depth, e.g., vf = vs (a/r)
c       and that density  varies
c       as rhof = rhos (a/r)^-P, [not the TI surface wave code has not yet
c        been written], then using the model that m = rho beta^2, we have
c
c       TA = rho VA^2,
c       TAf = rhof * VAf^2 = rhos (a/r)^-P VAs^2 (a/r)^2
c           = (a/r)^2-P TAs
c-----
        integer NL
        parameter (NL=1000)
        common/timod/TD(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real*4 TD, TA, TC, TF, TL, TN, TRho
        real*4 qa, qb, etap, etas, frefp, frefs
        common/depref/refdep
        real refdep


        common/pari/ mmax,mode
        integer mmax, mode
        common/water/iwat(NL)
        double precision z0,z1,r0,r1,ar,tmp
C       double precision va,vc,cl,vn

        common/earth/radius
        real radius

        ar=radius
        r0=ar + refdep
        td(mmax)=1.0
        do 10 i=1,mmax
            r1=r0-dble(td(i))
            z0=ar*dlog(ar/r0)
            z1=ar*dlog(ar/r1)
            td(i)=z1-z0
c-----
c        attempt 7 15 2007 - use standard rule but at mid layer depth as per DGH
c-----
            TMP=(ar+ar)/(r0+r1)
         if(ifunc.eq.2)then
c-----
c                SV
c-----
                 rhosph    = trho(i)
                 trho(i)   = rhosph * tmp**(-2.275)

                 ta(i)=ta(i)*tmp**(-0.2750)
                 tc(i)=tc(i)*tmp**(-0.2750)
C                 tf(i)=tf(i)*tmp**(-0.2750)

                 elsph = tl(i)
                 tl(i)  =elsph*tmp**(-0.2750)
                 ensph = tn(i)
                 tn(i)=ensph*tmp**(-0.2750)

            else if(ifunc.eq.1)then
                 rhosph    = trho(i)
C                 vl = sqrt(tl(i)/rhosph)*tmp
C                 vn = sqrt(tn(i)/rhosph)*tmp
                 trho(i) = rhosph * tmp**(-5)
                 elsph = tl(i)
                 tl(i)=elsph*tmp**(-3.0)
                 ensph = tn(i)
                 tn(i)=ensph*tmp**(-3.0)
C              WRITE(6,*)i,tl(i),vl*vl*trho(i),tn(i),vn*vn*trho(i)
            endif
            r0 = r1
   10   continue
        td(mmax)=0.0
        return
        end

        subroutine gcmdln(verby,cmin,cmax)
        logical verby
        character names*40
        real*4 cmin,cmax

        verby = .false.
        cmin = -1.0
        cmax = -1.0
        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 1100
            call mgtarg(i,names)
            if(names(1:2).eq.'-v')then
                verby = .true.
c----
c      fudge to force the dispersion search to start at a specific value
c-----
            else if(names(1:5).eq.'-cmin')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')cmin
            else if(names(1:5).eq.'-cmax')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')cmax
            else if(names(1:2).eq.'-?')then
                call usage(' ')
            else if(names(1:2).eq.'-h')then
                call usage(' ')
            endif
        go to 1000
 1100   continue
            
        return
        end

        subroutine usage(str)
        character str*(*)
        parameter(LER=0)
        write(LER,*)'Usage:',str
        write(LER,*)
     1  'tdisp96 -v -? -h'
        write(LER,*)
     1  '  -v   (default false) Screen output of current computation'
        write(LER,*)
     1  '  -?   this help message'
        write(LER,*)
     1  '  -h   this help message'
        stop
        end
        
c-----
c       notes:
c       For the surface wave the propagator is always positive
c       However the E and E^-1 matrices can be complex
c-----

        subroutine varsh(h,rsh,lshimag,cossh,rsinsh,sinshr,ex)
        implicit none
        double precision h
        double precision rsh
        logical lshimag
        double precision cossh, rsinsh, sinshr
        double precision ex
        
        double precision q, fac, sinq
        q = rsh*h
        ex =  0.0
        if(lshimag)then
            if(rsh.gt.0.0)then
                 cossh = dcos(q)
                 sinq = sin(q)
                 sinshr = sinq/rsh
                 rsinsh = - rsh*sinq
            else
                 cossh  = 1.0d+00
                 sinshr = dble(h)
                 rsinsh = 0.0
            endif
        else
            ex = q
            if(q.lt.16.0d+00)then
                fac = dexp(-2.0d+00*q)
            else
                fac = 0.0d+00
            endif
            cossh = (1.0d+00 + fac) * 0.5d+00
            sinq   = (1.0d+00 - fac) * 0.5d+00
            sinshr = sinq/rsh
            rsinsh = sinq*rsh
        endif
        return
        end
        

        subroutine hskl(cossh,rsinsh,sinshr,TL,iwat,hl,exm,exe)
c-----
c       True cosh( rd b )    = exp(exm) cossh
c       True sinh( rd b )/rb = exp(exm) sinshr
c       True rb sinh( rd b ) = exp(exm) rsinsh
c-----
        implicit none
        double precision cossh, rsinsh, sinshr
        double precision hl(2,2),exm,exe
        double precision TL
        integer iwat

        if(iwat.eq.0)then   
            hl(1,1) = cossh
            hl(2,1) = TL*rsinsh
            hl(1,2) = sinshr/(TL)
            hl(2,2) = cossh
            exe = exe + exm
        else
            hl(1,1) = 1.0
            hl(1,2) = 0.0
            hl(2,1) = 0.0
            hl(2,2) = 1.0
        endif
        return
        end 

        subroutine gtesh(esh,einvsh,rsh,wvno,L,lshimag,ltrueinv,iwat)
        complex*16 esh(2,2), einvsh(2,2)
        double precision wvno, rsh
        double precision L
        logical lshimag,ltrueinv
        integer iwat
        double complex rb
c-----
c       get E and E^-1 for Quasi SH
c       if ltrueinv == .true. then return E^-1
c       else return the normalized avoid divide by zero
c       (2k L rsh) E^-1 to avoid dividing by zero when rsh = 0
c-----
        if(iwat.eq.1)then
c-----
c           for a water layer use an identity matrix. The
c           results have validity only for the extreme case of
c           a model that is fluid-solid-fluid, e.g., ocean,mantle,
c           outer core.  The propagator matrix technique cannot 
c           handle anything more complicated
c-----
             esh(1,1) = dcmplx(1.0d+00, 0.0d+00)
             esh(1,2) = dcmplx(0.0d+00, 0.0d+00)
             esh(2,1) = dcmplx(0.0d+00, 0.0d+00)
             esh(2,2) = dcmplx(1.0d+00, 0.0d+00)
             einvsh(1,1) = dcmplx(1.0d+00, 0.0d+00)
             einvsh(1,2) = dcmplx(0.0d+00, 0.0d+00)
             einvsh(2,1) = dcmplx(0.0d+00, 0.0d+00)
             einvsh(2,2) = dcmplx(1.0d+00, 0.0d+00)
        else
             if(lshimag)then
                 rb =  dcmplx(0.0d+00,rsh)
             else
                 rb = dcmplx(rsh,0.0d+00)
             endif
             esh(1,1) =   wvno
             esh(1,2) =   wvno
             esh(2,1) =   wvno * L * rb
             esh(2,2) = - wvno * L * rb
             einvsh(1,1) =   L * rb
             einvsh(2,1) =   L * rb
             einvsh(1,2) =   1.0
             einvsh(2,2) =  -1.0
             if(ltrueinv)then
                 einvsh(1,1) = einvsh(1,1)/(2.*wvno*L*rb)
                 einvsh(1,2) = einvsh(1,2)/(2.*wvno*L*rb)
                 einvsh(2,1) = einvsh(2,1)/(2.*wvno*L*rb)
                 einvsh(2,2) = einvsh(2,2)/(2.*wvno*L*rb)
             endif
        endif
        return
        end

        subroutine gtegsh(m,wvno,omega,rsh,lshimag)
        implicit none

        integer m
        double precision wvno, omega, rsh
        logical lshimag

        integer NL
        parameter (NL=1000)
        common/timod/H(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real H, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/water/iwat(NL)
        integer iwat


c-----
c       internal variables
c-----
        double precision wvno2, omega2
        wvno2 = wvno*wvno
        omega2 = omega*omega
        if(iwat(m).eq.1)then
c-----
c           check for water layer
c-----
            rsh = 0.0
            lshimag = .false.
        else
            rsh = TN(m)*wvno2/(TL(m)) - TRho(m)*omega2/(TL(m))
            if(rsh .lt. 0.0)then
                rsh = dsqrt(dabs(rsh))
                lshimag = .true.
            else
                rsh = dsqrt(dabs(rsh))
                lshimag = .false.
            endif
        endif
        
        return
        end

        subroutine lmult(d11,d12,d21,d22,hl,iwat,exel,exb,icomp)
        implicit none
c-----
c       multiply SH matrix by a row vector on left
c-----
        complex*16 d11,d12,d21,d22,hl(2,2),e1,e2
        integer iwat
        real*8 exel, exb
        logical icomp
c-----
c       fluid layer do nothing, just return, 
c           equivalent to multiplying by
c       identity matrix
c-----
        if(iwat.eq.0)then
c-----
c       elastic layer
c-----
            e1=d11
            e2=d12
c-----
c           a11 = cosql
c           a12 = yl
c           a21 = zl
c           a22 = cosql
c-----
            d11=e1*hl(1,1) + e2*hl(2,1)
            d12=e1*hl(1,2) + e2*hl(2,2)
            exel = exel + exb
            if(icomp)then
                e1=d21
                e2=d22
                d21=e1*hl(1,1) + e2*hl(2,1)
                d22=e1*hl(1,2) + e2*hl(2,2)
            endif
        endif
        return
        end

        function dltar4(wvno,omga)
c-----
c       find P-SV dispersion values.
c-----
        parameter (NL=1000)
        implicit double precision (a-h,o-z)
        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
        common/pari/ mmax,mode
        integer mmax,mode
        common/water/iwat(NL)
        integer iwat
        real*8 e(5),ee(5),ca(5,5)

        complex*16 gbr(2,5)
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv
        complex*16 p,q
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        real*8 cosp , cossv
        real*8 rsinp, rsinsv
        real*8 sinpr, sinsvr
        real*8 pex, svex
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid

CRBH    integer i,j
c-----
c       set up starting values for bottom halfspace
c-----
        wvno2=wvno*wvno
        om2 = omga*omga
        call evalg(0,mmax,mmax-1,gbr,1,
     1      wvno,omga,om2,wvno2)

c-----
c       CORRECTED 14 OCT 2006
c-----
            e(1) = dreal(gbr(1,1) )
            e(2) = dreal(gbr(1,2) )
            e(3) = dreal(gbr(1,3) )
            e(4) = dreal(gbr(1,4) )
            e(5) = dreal(gbr(1,5) )
C        WRITE(6,'(i4,5e11.3)')mmax,e
c-----
c-----
c       matrix multiplication from bottom layer upward
c-----
        mmm1 = mmax-1
        do 500 m = mmm1,1,-1
           call  gettiegn(Za,Zb,Zc,Zd,Ze,Zf,om2,wvno2,rp, rsv,
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,omga,wvno,
     2      dcmplx(1.0d+00,0.0d+00),dcmplx(1.0d+00,0.0d+00))
                p = rp  * dble(d(m))
                q = rsv * dble(d(m))
                call varsv(p,q,rp, rsv, 
     1              cosp, cossv, rsinp, rsinsv, 
     1              sinpr, sinsvr, pex,svex,iwat(m))
                call dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,
     1              NP,NSV,
     1              x11,x21,x31,x41,x12,x22,x32,x42,
     1              TRho(m),iwat(m),pex+svex,om2)

            nmat = 5
            do 200 i=1,nmat
                cr=0.0d+00
                do 100 j=1,nmat
                    cr=cr+e(j)*ca(j,i)
  100           continue
                ee(i)=cr
  200       continue
            call normc(ee,exa,nmat)
            do 300 i = 1,nmat
                e(i)=ee(i)
  300       continue
  500   continue
        dltar4 = e(1)
        return
        end

        subroutine dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,NP,NSV,
     1      x11,x21,x31,x41,x12,x22,x32,x42,
     1      TRho,iwat,ex,om2)
        implicit none
        real*8 CA(5,5)
        real*8 cosp , cossv 
        real*8 rsinp, rsinsv
        real*8 sinpr, sinsvr
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        real TRho
        integer iwat
        real*8 ex, dfac
        real*8  om2
        complex*16 zrho

        integer i, j
        complex*16 TCA(6,6)
c-----
c       introduce conciseness to reduce operations
c-----
        COMPLEX*16 NPNP
        COMPLEX*16 NPNSV
        COMPLEX*16 NSVNSV
        COMPLEX*16 CPCSV, CPSVR, CPRSV, SPRCSV, SPRSVR
        COMPLEX*16 SPRRSV, RSPCSV,  RSPSVR, RSPRSV
        COMPLEX*16 X11X11, X11X21, X11X31, X11X41
        COMPLEX*16         X21X21, X21X31, X21X41
        COMPLEX*16                 X31X31, X31X41
        COMPLEX*16                         X41X41
        COMPLEX*16 X12X12, X12X22, X12X32, X12X42
        COMPLEX*16         X22X22, X22X32, X22X42
        COMPLEX*16                 X32X32, X32X42
        COMPLEX*16                         X42X42

        COMPLEX*16 FAC01, FAC02, FAC03, FAC04, FAC05
        COMPLEX*16 FAC06, FAC07, FAC08, FAC09, FAC10
        COMPLEX*16 FAC11, FAC12, FAC13, FAC14, FAC15
        COMPLEX*16 FAC16, FAC17, FAC18, FAC19, FAC20
c-----
c        A11     A12     A13    -A13     A15     A16
c        A21     A22     A23    -A23     A25     A15
c        A31     A32     A33    1-A33   -A23    -A13
c       -A31    -A32    1-A33    A33     A23     A13
c        A51     A52    -A32     A32     A22     A12
c        A61     A51    -A31     A31     A21     A11
c-----
c       this will be multipled on the left by the G matrix
c
c       [ G11   G12 G13 -G13    G15 G16 ]
c
c-----
c       or on the right by
c
c       [ H11   H21 H31 -H31    H51 H61  ] ^T
c-----
c       the number of multiplications can be reduced from 
c            36 to 25 if we define a new matrices
c       related to the original matrices by
c-----
c         A11     A12     A13         A15     A16
c         A21     A22     A23         A25     A15
c       2 A31   2 A32    2 A33 -1   -2 A23  -2 A13
c         A51     A52    -A32         A22     A12
c         A61     A51    -A31         A21     A11
c-----
c
c       [ G11   G12  G13    G15 G16  ]
c       [ H11   H21 2 H31   H51 H61  ] ^T
c
c-----
c       this means that some of the original definitions of the 
c           Aij elements must be changed for the
c       definition of the modified 5x5 compount A matrix
c
c       old 6x6                 new 5x5
c       A11 = 1 - A33 + A22     1 - (1/2)(new A33 + 1) + new A2)
c       A53 = -A32              A43 = - (1/2) new A32
c       A63 = -A31              A53 = - (1/2) new A31
c-----
c       To recover the needed elements, we note that the 
c           old G14 = -old G14 = new G13
c-----
c-----
        zrho = dcmplx(dble(TRho),0.0d+00)
        if(ex.gt.35.0d+00)then
            dfac = 0.0d+00
        else
            dfac = dexp(-ex)
        endif
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,5
                do 101 i=1,5
                    ca(i,j) = 0.0d+00
  101           continue
  100       continue
            ca(3,3) = dfac 
            ca(1,1) = cosp 
            ca(5,5) = cosp 
            ca(1,2) = -rsinp/(zrho*om2) 
            ca(2,1) = - zrho*sinpr*om2 
            ca(2,2) = cosp 
            ca(4,4) = cosp 
            ca(4,5) = ca(1,2) 
            ca(5,4) = ca(2,1) 
        else
c-----
c           introduce some repeated expressions to 
c           reduce multiplications
c-----
            NPNP   = dfac/(NP*NP)
            NSVNSV = dfac/(NSV * NSV)
            NPNSV  = 1.0d+00/(NP * NSV)
            CPCSV  = cosp *cossv
            CPSVR  = cosp *sinsvr
            CPRSV  = cosp *rsinsv
            SPRCSV = sinpr*cossv
            SPRSVR = sinpr*sinsvr
            SPRRSV = sinpr*rsinsv
            RSPCSV = rsinp*cossv
            RSPSVR = rsinp*sinsvr
            RSPRSV = rsinp*rsinsv

            X11X11 = x11*x11
            X11X21 = x11*x21    
            X11X31 = x11*x31    
            X11X41 = x11*x41    
            X21X21 = x21*x21
            X21X31 = x21*x31
            X21X41 = x21*x41
            X31X31 = x31*x31
            X31X41 = x31*x41
            X41X41 = x41*x41

            X12X12 = x12*x12
            X12X22 = x12*x22
            X12X32 = x12*x32
            X12X42 = x12*x42
            X22X22 = x22*x22
            X22X32 = x22*x32
            X22X42 = x22*x42
            X32X32 = x32*x32
            X32X42 = x32*x42
            X42X42 = x42*x42

            FAC01= X11X11*X22X32*RSPCSV*(NPNSV)
            FAC02=X11X11*X22X42*RSPRSV*(NPNSV)
            FAC03=X11X11*X32X42*RSPCSV*(NPNSV)
            FAC04=X11X21*X12X32*CPSVR*(NPNSV)
            FAC05=X11X21*X12X42*CPCSV*(NPNSV)
            FAC06=X11X31*X12X22*RSPCSV*(NPNSV)
            FAC07=X11X31*X12X32*RSPSVR*(NPNSV)
            FAC08=X11X31*X12X42*RSPCSV*(NPNSV)
            FAC09=X11X41*X12X22*CPCSV*(NPNSV)
            FAC10=X11X41*X12X32*CPSVR*(NPNSV)
            FAC11=X11X41*X22X32*CPCSV*(NPNSV)
            FAC12=X21X31*X12X12*CPSVR*(NPNSV)
            FAC13=X21X31*X12X42*CPCSV*(NPNSV)
            FAC14=X21X31*X42X42*CPRSV*(NPNSV)
            FAC15=X21X41*X12X12*SPRSVR*(NPNSV)
            FAC16=X21X41*X32X42*SPRCSV*(NPNSV)
            FAC17=X31X41*X12X12*CPSVR*(NPNSV)
            FAC18=X31X41*X22X42*CPRSV*(NPNSV)
            FAC19=X31X41*X32X42*CPCSV*(NPNSV)
            FAC20=X41X41*X22X32*SPRCSV*(NPNSV)        

c-----
c       repeated terms
c       X11X11*X22X32*RSPCSV*(NPNSV)
c       X11X11*X22X42*RSPRSV*(NPNSV)
c       X11X11*X32X42*RSPCSV*(NPNSV)
c       X11X21*X12X32*CPSVR*(NPNSV)
c       X11X21*X12X42*CPCSV*(NPNSV)
c       X11X31*X12X22*RSPCSV*(NPNSV)
c       X11X31*X12X32*RSPSVR*(NPNSV)
c       X11X31*X12X42*RSPCSV*(NPNSV)
c       X11X41*X12X22*CPCSV*(NPNSV)
c       X11X41*X12X32*CPSVR*(NPNSV)
c       X11X41*X22X32*CPCSV*(NPNSV)
c       X21X31*X12X12*CPSVR*(NPNSV)
c       X21X31*X12X42*CPCSV*(NPNSV)
c       X21X31*X42X42*CPRSV*(NPNSV)
c       X21X41*X12X12*SPRSVR*(NPNSV)
c       X21X41*X32X42*SPRCSV*(NPNSV)
c       X31X41*X12X12*CPSVR*(NPNSV)
c       X31X41*X22X42*CPRSV*(NPNSV)
c       X31X41*X32X42*CPCSV*(NPNSV)
c       X41X41*X22X32*SPRCSV*(NPNSV)        
c-----

c-----
c       ALSO NOTE THAT NONE OF THE TCA(?,4) or TCA(4,?) ARE REQUIRED
c       SO DO NOT COMPUTE
c-----
            
c-----
c       elastic layer
c-----
c CA11 12 12 = 11 22 - 12 21
c CA12 12 13 = 11 23 - 13 21
c CA13 12 14 = 11 24 - 14 21 = - 11 13 - 14 21
c CA14 12 23 = 12 23 - 22 13
c CA15 12 24 = 12 24 - 22 14
c CA16 12 34 = 13 24 - 23 14
c CA21 13 12 = 11 32 - 12 31
c CA22 13 13 = 11 33 - 13 31
c CA23 13 14 = 11 34 - 14 31
c CA24 13 23 = 12 33 - 13 32 = 12 22 - 13 32
c CA25 13 24  = 12 34 - 14 32
c CA26 13 34  = 13 34 - 14 33
c CA31 14 12 = 11 42 - 12 41 = - 11 31 - 12 41
c CA32 14 13 = 11 43 - 13 41
c CA33 14 14 = 11 44 - 14 41 = 11 11 - 14 41
c CA34 14 23 = 12 43 - 13 42 = - 12 21 + 13 31
c CA35 14 24 = 12 44 - 14 42 = 12 11 + 14 31
c CA36 14 34 = 13 44 - 14 43 = 13 11 + 14 21  = - CA13
c CA41 23 12 = 21 32 - 22 31 
c CA42 23 13 = 21 33 - 23 31 
c CA43 23 14 = 21 34 - 24 31 = - 21 12 + 13 31
c CA44 23 23 = 22 33 - 23 32 
c CA45 23 24 = 22 34 - 24 32 = - 22 12 + 13 32 
c                            = - 33 12 + 13 32 = - CA24
c CA46 23 34 = 23 34 - 24 33 = - 23 12 + 13 22 = - CA14 
c CA51 24 12 = 21 42 - 22 41 = - 21 31 - 22 41
c CA52 24 13 = 21 43 - 23 41 = - 21 21 - 23 41
c CA53 24 14 = 21 44 - 24 41 = - 11 43 + 13 41 = - CD32
c CA54 24 23 = 22 43 - 23 42 = - 22 21 + 23 31 
c                            = - 33 21 + 23 31 = -CD42
c CA55 24 24 = 22 44 - 24 42 = 33 11 - 13 31 = CD22 
c CA56 24 34 = 23 44 - 24 43 = 23 11 - 13 21 = CD12
c CA61 34 12 = 31 42 - 32 41 = - 31 31 - 32 41
c CA62 34 13 = 31 43 - 33 41 = - 21 31 - 22 41
c CA63 34 14 = 31 44 - 34 41 = - 11 43 + 13 41 = - CD32
c CA64 34 23 = 32 43 - 33 42 = - 22 21 + 23 31 
c                            = - 33 21 + 23 31 = -CD42
c CA65 34 24 = 32 44 - 34 42 = 32 11 - 12 31 = CD21 
c CA66 34 34 = 33 44 - 34 43 = 22 11 - 12 21 = CD11

        
            TCA(1,1) = - X11X21*X31X41*(NPNP ) 
     1            - X12X22*X32X42*(NSVNSV) 
     1            - FAC11
     1            - FAC13
     1            + X11X31*X22X42*RSPRSV*(NPNSV)
     1            + X21X41*X12X32*SPRSVR*(NPNSV)
            TCA(1,2) = - X11X41*X22X22*CPRSV*(NPNSV) 
     1            - X21X21*X12X42*SPRCSV*(NPNSV) 
     1            + X11X21*X22X42*CPRSV*(NPNSV)
     1            + X21X41*X12X22*SPRCSV*(NPNSV)
            TCA(1,3) =   X11X11*X21X41*(NPNP )
     1            + X12X12*X22X42*(NSVNSV)
     1            + FAC09
     1            + FAC05
     1            - FAC15
     1            - FAC02
C           TCA(1,4) = - X11X21*X21X31*(NPNP )
C     1           - X12X22*X22X32*(NSVNSV)
C     1           + X11X31*X22X22*RSPRSV*(NPNSV)
C     1           + X21X21*X12X32*SPRSVR*(NPNSV)
C     1           - X21X31*X12X22*CPCSV*(NPNSV)
C     1           - X11X21*X22X32*CPCSV*(NPNSV)
            TCA(1,5) =  - FAC06
     1             - FAC04
     1             + FAC12
     1             + FAC01
            TCA(1,6) = - X11X11*X21X21*(NPNP )
     1            - X12X12*X22X22*(NSVNSV)
     1            - 2.0*X11X21*X12X22*CPCSV*(NPNSV)
     1            + X21X21*X12X12*SPRSVR  *(NPNSV)
     1            + X11X11*X22X22*RSPRSV  *(NPNSV)
            TCA(2,1) = - X11X41*X32X32*CPSVR*(NPNSV)
     1            - X31X31*X12X42*RSPCSV*(NPNSV)
     1            + X11X31*X32X42*RSPCSV*(NPNSV)
     1            + X31X41*X12X32*CPSVR*(NPNSV)
            TCA(2,2) = - FAC11
     1            - FAC13
     1            + X11X21*X32X42*CPCSV*(NPNSV)
     1            + X31X41*X12X22*CPCSV*(NPNSV)
            TCA(2,3) =   FAC10
     1            + FAC08
     1            - FAC03
     1            - FAC17
C           TCA(2,4) = + X11X31*X22X32*RSPCSV*(NPNSV)
C     1           + X21X31*X12X32*CPSVR*(NPNSV)
C     1           - X11X21*X32X32*CPSVR*(NPNSV)
C     1           - X31X31*X12X22*RSPCSV*(NPNSV)
            TCA(2,5) = - 2.0d+00 * FAC07
     1            + X11X11*X32X32*RSPSVR*(NPNSV)
     1            + X31X31*X12X12*RSPSVR*(NPNSV)
            TCA(2,6) = - FAC04
     1            - FAC06
     1            + FAC01
     1            + FAC12
            TCA(3,1) = - X11X31*X41X41*(NPNP )
     1            - X12X32*X42X42*(NSVNSV)
     1            - X11X41*X32X42*CPCSV*(NPNSV)
     1            - X31X41*X12X42*CPCSV*(NPNSV)
     1            + X11X31*X42X42*RSPRSV*(NPNSV)
     1            + X41X41*X12X32*SPRSVR*(NPNSV)
            TCA(3,2) = - X11X41*X22X42*CPRSV*(NPNSV)
     1            - X21X41*X12X42*SPRCSV*(NPNSV)
     1            + X11X21*X42X42*CPRSV*(NPNSV)
     1            + X41X41*X12X22*SPRCSV*(NPNSV)
            TCA(3,3) =   X11X11*X41X41*(NPNP )
     1            + X12X12*X42X42*(NSVNSV)
     1            + 2.0*X11X41*X12X42*CPCSV*(NPNSV)
     1            - X11X11*X42X42*RSPRSV*(NPNSV)
     1            - X41X41*X12X12*SPRSVR*(NPNSV)
C           TCA(3,4) = - X11X21*X31X41*(NPNP )
C     1           - X32X42*X12X22*(NSVNSV)
C     1           + X11X31*X22X42*RSPRSV*(NPNSV)
C     1           + X21X41*X12X32*SPRSVR*(NPNSV)
C     1           - X11X21*X32X42*CPCSV*(NPNSV)
C     1           - X31X41*X12X22*CPCSV*(NPNSV)
            TCA(3,5) = - FAC08
     1            - FAC10
     1            + FAC03
     1            + FAC17
            TCA(3,6) =  - TCA(1,3)
C           TCA(4,1) =   X21X41*X31X31*(NPNP )
C     1           + X22X42*X32X32*(NSVNSV)
C     1           - X21X41*X32X32*SPRSVR*(NPNSV)
C     1           - X31X31*X22X42*RSPRSV*(NPNSV)
C     1           + X21X31*X32X42*CPCSV  *(NPNSV)
C     1           + X31X41*X22X32*CPCSV  *(NPNSV)
C           TCA(4,2) = - X21X41*X22X32*SPRCSV*(NPNSV)
C     1           - X21X31*X22X42*CPRSV*(NPNSV)
C     1           + X31X41*X22X22*CPRSV*(NPNSV)
C     1           + X21X21*X32X42*SPRCSV*(NPNSV)
C           TCA(4,3) = - X31X41*X11X21*(NPNP )
C     1           - X12X22*X32X42*(NSV*NSV)
C     1           - X31X41*X12X22*CPCSV*(NPNSV)
C     1           - X11X21*X32X42*CPCSV*(NPNSV)
C     1           + X11X31*X22X42*RSPRSV*(NPNSV)
C     1           + X21X41*X12X32*SPRSVR*(NPNSV)
C           TCA(4,4) =   X21X31*X21X31*(NPNP )
C     1           + X22X32*X22X32*(NSVNSV)
C     1           + X21X31*X22X32*CPCSV  *(NPNSV)
C     1           + X21X31*X22X32*CPCSV  *(NPNSV)
C     1           - X21X21*X32X32*SPRSVR*(NPNSV)
C     1           - X31X31*X22X22*RSPRSV*(NPNSV)
C           TCA(4,5) =   TCA(2,3)
C           TCA(4,6) = - TCA(1,4)
            TCA(5,1) = - FAC16
     1            - FAC18
     1            + FAC14
     1            + FAC20
            TCA(5,2) = - 2.0 * X21X41*X22X42*SPRRSV*(NPNSV)
     1            + X21X21*X42X42*SPRRSV*(NPNSV)
     1            + X41X41*X22X22*SPRRSV*(NPNSV)
            TCA(5,3) = - TCA(3,2)
C           TCA(5,4) = - TCA(4,2)
            TCA(5,5) =   TCA(2,2)
            TCA(5,6) =   TCA(1,2)
            TCA(6,1) = - X31X41*X31X41*(NPNP )
     1            - X32X42*X32X42*(NSVNSV)
     1            - 2.0d+00 *FAC19
     1            + X31X31*X42X42*RSPRSV *(NPNSV)
     1            + X41X41*X32X32*SPRSVR *(NPNSV)
            TCA(6,2) = - FAC18
     1            - FAC16
     1            + FAC14
     1            + FAC20
            TCA(6,3) = - TCA(3,1)
C           TCA(6,4) =   TCA(3,1)
            TCA(6,5) =   TCA(2,1)
            TCA(6,6) =   TCA(1,1)
c-----
c       for development only - clean up later
c       define the CA(5,5)
c-----
            CA(1,1) = dreal( TCA(1,1) )
            CA(1,2) = dreal( TCA(1,2) )
            CA(1,3) = dreal( TCA(1,3) )
            CA(1,4) = dreal( TCA(1,5) )
            CA(1,5) = dreal( TCA(1,6) )

            CA(2,1) = dreal( TCA(2,1) )
            CA(2,2) = dreal( TCA(2,2) )
            CA(2,3) = dreal( TCA(2,3) )
            CA(2,4) = dreal( TCA(2,5) )
            CA(2,5) = dreal( TCA(1,5) )

            CA(3,1) =   2.0d+00*dreal( TCA(3,1) )
            CA(3,2) =   2.0d+00*dreal( TCA(3,2) )
c-----
c       beware of normalization
c-----
            CA(3,3) =   2.0d+00*dreal( TCA(3,3)) - 1.0d+00*dfac 
            CA(3,4) =  -2.0d+00*dreal( TCA(2,3) )
            CA(3,5) =  -2.0d+00*dreal( TCA(1,3) )

            CA(4,1) =   dreal( TCA(5,1) )
            CA(4,2) =   dreal( TCA(5,2) )
            CA(4,3) =  -dreal( TCA(3,2) )
            CA(4,4) =   dreal( TCA(2,2) )
            CA(4,5) =   dreal( TCA(1,2) )

            CA(5,1) =   dreal( TCA(6,1) )
            CA(5,2) =   dreal( TCA(5,1) )
            CA(5,3) =  -dreal( TCA(3,1) )
            CA(5,4) =   dreal( TCA(2,1) )
            CA(5,5) =   dreal( TCA(1,1) )
        endif
        return
        end

        subroutine evalg(jbdry,m,m1,gbr,in,
     1      wvno,om,om2,wvno2)
        implicit none
        integer jbdry, m, m1, in
        complex*16 gbr(2,5)
        real*8 wvno,om,wvno2,om2
        integer NL
        parameter(NL=1000)
        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
        common/modlly/mmax
        integer mmax
        common/modspec/allfluid
        logical allfluid
        complex*16 CDSQRT

        complex*16 atna, atnb
        complex*16 cg(6)
        complex*16 g(4,4)
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        integer iwat

        complex*16 e(4,4), einv(4,4)

        integer i,j,k
        complex*16 zsum


C        WRITE(6,*)'allfluid:',allfluid
c-----
c       set up halfspace conditions
c-----
            if(TL(m).eq. 0.0 .or.TN(m).eq.0.0)then
                iwat = 1
            else
                iwat = 0
            endif
c-----
c       HALFSPACE
c-----
            call gettiegn(za,zb,zc,zd,ze,zf,om2,wvno2,rp, rsv, 
     1          x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,om,wvno,
     2          dcmplx(1.0d+00,0.0d+00),dcmplx(1.0d+00,0.0d+00))
            if(iwat.eq.0)then
c-----
c               ELASTIC HALFSPACE
C       WRITE(6,*)'X:'
C       WRITE(6,*)x11,x21,x31,x41
C       WRITE(6,*)x12,x22,x32,x42
C        WRITE(6,*)'RP,RSV:',RP,RSV
C        WRITE(6,*)'NP,NSV,m,om,wvno:',NP,NSV,m,om,wvno
C       WRITE(6,*)'om2,wvno2:',om2,wvno2


                G(1,1) =    x41/(2.0*NP  *rp  )
                G(1,2) = -  x31/(2.0*NP       )
                G(1,3) = -  x21/(2.0*NP  *rp  )
                G(1,4) =    x11/(2.0*NP       )
                G(2,1) =    x42/(2.0*NSV      )
                G(2,2) = -  x32/(2.0*NSV *rsv )
                G(2,3) = -  x22/(2.0*NSV      )
                G(2,4) =    x12/(2.0*NSV *rsv )
                G(3,1) =    x41/(2.0*NP  *rp  )
                G(3,2) =    x31/(2.0*NP       )
                G(3,3) = -  x21/(2.0*NP  *rp  )
                G(3,4) = -  x11/(2.0*NP       )
                G(4,1) = -  x42/(2.0*NSV      )
                G(4,2) = -  x32/(2.0*NSV *rsv )
                G(4,3) =    x22/(2.0*NSV      )
                G(4,4) =    x12/(2.0*NSV *rsv )
          E(1,1)=  X11*rp
          E(2,1)=  X21
          E(3,1)=  X31*rp
          E(4,1)=  X41
          E(1,2)=  X12
          E(2,2)=  X22*rsv
          E(3,2)=  X32
          E(4,2)=  X42*rsv
          E(1,3)=  X11*rp
          E(2,3)= -X21
          E(3,3)=  X31*rp
          E(4,3)= -x41
          E(1,4)= -X12
          E(2,4)=  X22*rsv
          E(3,4)= -x32
          E(4,4)=  X42*rsv
C          do j=1,4
C              do i = 1,4
C                   write(6,*)'E   (',I ,',' , J, ')=',E(i,J)
C             enddo
C          enddo
C          do j=1,4
C              do i = 1,4
C                   write(6,*)'EINV(',I ,',' , J, ')=',G(i,J)
C             enddo
C          enddo
C          do i=1,4
C              do j = 1,4
C                zsum = dcmplx(0.0d+00,0.0d+00)
C                do k=1,4
C                zsum = zsum + E(i,k)*g(k,j)
C                enddo
C                write(6,*)'E INV(',I ,',' , J, ')=',ZSUM
C             enddo
C          enddo

c CG(1) = G 12 12 = 11 22 - 12 21
c CG(2) = G 12 13 = 11 23 - 13 21
c CG(3) = G 12 14 = 11 24 - 14 21
c CG(4) = G 12 23 = 12 23 - 13 22
c CG(5) = G 12 24 = 12 24 - 14 22
c CG(6) = G 12 34 = 13 24 - 14 23
                CG(1) = G(1,1)*G(2,2) - G(1,2)*G(2,1)
                CG(2) = G(1,1)*G(2,3) - G(1,3)*G(2,1)
                CG(3) = G(1,1)*G(2,4) - G(1,4)*G(2,1)
                CG(4) = G(1,2)*G(2,3) - G(1,3)*G(2,2)
                CG(5) = G(1,2)*G(2,4) - G(1,4)*G(2,2)
                CG(6) = G(1,3)*G(2,4) - G(1,4)*G(2,3)
                gbr(in,1) = CG(1)
                gbr(in,2) = CG(2)
                gbr(in,3) = CG(3)
                gbr(in,4) = CG(5)
                gbr(in,5) = CG(6)

            else if(iwat.eq.1)then
c-----
c               FLUID HALFSPACE
c-----
                if(allfluid)then
                    gbr(in,1) = dble(TRho(m))*om2
                    gbr(in,2) = -rp
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                else
                    gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = -dble(TRho(m))*om2
                    gbr(in,5) = rp
                endif
            endif
        return
        end

        subroutine gettiegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, 
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,omg,wvn,
     2      atna,atnb)
        implicit none
        COMPLEX*16 A,B,C,D,E,F
        real*8 wvno2, omega2
        COMPLEX*16 rp, rsv
c-----
c       get the modified eigen vectors x11 and x31 have 
c           lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 have 
c           lambda2 (rsv) factored out
c-----
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
c-----
c       norms
c-----
        COMPLEX*16 NP, NSV
        integer m
        real*8 omg, wvn
        complex*16 atna, atnb

        COMPLEX*16 xka2, xkb2

        integer NL
        parameter (NL=1000)
        common/timod/td(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real*4 td,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
c-----
c       internal variables
c-----
        COMPLEX*16 L2(2)
        COMPLEX*16 bb, cc
        COMPLEX*16 CDSQRT
        COMPLEX*16 SRTBB4AC
        COMPLEX*16 ddef, aabc

        COMPLEX*16 ZFAC
c-----
c       first test to see if a fluid layer - if it is fluid, the
c       eigenfunctions are specially computed and we need only the
c       rp
c-----
        if(TL(m).eq.0.0 .or. TN(m).eq.0.0)then
            rp = cdsqrt(dcmplx(wvno2 -omega2*TRho(m)/(TA(m)), 0.0d+00))
            rsv = dcmplx(0.0d+000, 0.0d+00)
            return
        endif


        a = wvn * TF(m) / (TC(m))
        b = 1.0/(TC(m))
        c = - TRho(m)*omg*omg + wvn*wvn *
     1      (TA(m) -TF(m)*TF(m)/(TC(m)))
        d = - wvn
        e = 1.0/(TL(m))
        f = - TRho(m)*omg*omg

c-----
c       do algebra first to avoid numerical problems
c-----
        ddef = wvn*wvn - TRho(m)*omg*omg/(TL(m))
        aabc = wvn*wvn*TA(m)/TC(m) - TRho(m)*omg*omg/(TC(m))

c-----
c       Do the QUASI P and SV - WE MUST BE CAREFUL HERE CONCERNING
c       BRANCH CUTS OF THE SQUARE ROOT
c-----
c-----
c       The characteristic equation to be solved is
c
c       L^4 - L^2[ 2 ad +ec +fb ] + [ (d^2+ef)(a^2+bc)] = 0
c-----
        bb = 2.0d+00 * a*d + e*c +f*b
        cc = ddef * aabc
c----
c       ensure that the SQRT(bb*bb - 4.0D+00*cc) is in the
c       I and II quadrants
c-----

        SRTBB4AC = CDSQRT(bb*bb - 4.0D+00*cc)
        IF(DIMAG(SRTBB4AC) .lt.0.0D+00)THEN
            SRTBB4AC = - SRTBB4AC
        ENDIF
c-----
c       Determine L^2 with care for roundoff
c-----
        IF(DREAL(BB) .LT.0.0D+00 .AND. DREAL(SRTBB4AC).LT.0.0D+00)THEN
            L2(2) = ( bb - SRTBB4AC) / 2.0d+00
            L2(1) = cc/L2(2)
        ELSE
            L2(1) = ( bb + SRTBB4AC) / 2.0d+00
            L2(2) = cc/L2(1)
        ENDIF
c-----
c       Use the Lambda^2 values to form
c       xka^2 == k^2 - L(1)^2
c       xkb^2 == k^2 - L(2)^2
c       Associate the smallest xka, xkb with the P!
c-----
        xka2 = wvno2 - L2(1)
        xkb2 = wvno2 - L2(2)
        if(cdabs(xkb2) .lt. cdabs(xka2))THEN
                ZFAC = L2(1)
                L2(1) = L2(2)
                L2(2) = ZFAC
        endif
        rp  = CDSQRT(L2(1))
        rsv = CDSQRT(L2(2))

c-----
c       get the norms - note that the true norm will be 
c            2  NP amd 2 L(2) NSV
c       The factorization permits us to use the sin nz/n or n sin nz
c-----
C        NP  = (  L2(1)*(-2*a*b*d + 2*a*a*e + b*c*e - b*b*f)
C     1      + (a*a+b*c)*(2*b*d*d - 2*a*d*e + b*e*f - c*e*e) )
C        NSV = (- L2(2)*(2*b*d*d - 2*a*d*e - c*e*e + b*e*f)
C     1      + (d*d+e*f)*(2*a*b*d - 2*a*a*e + b*b*f - b*c*e) )
c-----
c       get the modified eigen vectors x11 and x31 have 
c           lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 have 
c           lambda2 (rsv) factored out
c-----
        x11 =              (b*d - a*e)
        x21 =  b*L2(1) - e*(b*c + a*a)
        x31 =    L2(1) -   (a*d + c*e)
        x41 = -a*L2(1) + d*(b*c + a*a)

        x12 = -e*L2(2) + b*(d*d + e*f)
        x22 = ( b*d - a*e)
        x32 = d*L2(2) - a*(d*d + e*f)
        x42 = - ( L2(2) -  a*d - b*f)
c-----
c       TEST
c       Force the eigenfunctions to be as given in 7.4.4
c-----
        zfac = rp / x21
        x11  = x11 *zfac
        x21  = rp
        x31  = x31 *zfac
        x41  = x41 *zfac

        zfac = rsv / x12
        x12  = rsv
        x22  = x22 * zfac
        x32  = x32 * zfac
        x42  = x42 * zfac
c-----
c       REDEFINE HERE USING THE adjusted eigenvectors
c       Note that TRUE NP  = 2 * np  * rp
c       Note that TRUE NSV = 2 * nsv * rsv     , where
c       but also note that since the 11 31 22 42 13 33 24 44 have
c       an eigenvalue factored out the normalization using these
c       elements will not use the rp and rsv
c-----  
        np   = x11*x41 - x21*x31
        nsv  = x12*x42 - x22*x32

        return
        end

        subroutine varsv(p,q, rp, rsv, 
     1      cosp, cosq, rsinp, rsinq, 
     1      sinpr, sinqr, pex,svex,iwat)
c-----
c       p = rp  * h
c       q = rsv * h
c       rp  vertical wave number for P
c       rsv vertical wave number for SV
c       cosp=cosh(p)  rsinp =rp *sinh(p)  sinpr = sinh(p)/rp
c       cosq=cosh(q)  rsinsv=rsv*sinh(p)  sinpq = sinh(p)/rsv
c-----
        implicit none
        COMPLEX*16 p, q
        COMPLEX*16 rp, rsv
        real*8 cosp, cosq
        real*8 rsinp, rsinq
        real*8 sinpr, sinqr
        REAL *8 pex,svex
        integer iwat

        REAL*8 pr, pi, qr, qi
        COMPLEX*16 epp, epm, eqp, eqm, erp, erm
        COMPLEX*16 sinp, sinq, sinr

        REAL*8 PFAC, SVFAC, SHFAC
        
        pex  = 0.0d+00
        svex = 0.0d+00
        pr = dreal(p)
        pi = dimag(p)
        qr = dreal(q)
        qi = dimag(q)
        pex   = pr
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            epp = dcmplx(dcos(pi), dsin(pi))/2.0
            epm = dconjg(epp)
            if(pr.lt.15.) then
                pfac=dexp(-2.*pr)
            else
                pfac  = 0.0d+00
            endif
            cosp = epp + pfac*epm
            sinp = epp - pfac*epm
            rsinp = rp *sinp
            sinpr = sinp/rp
            cosq  = 1.0d+00
            rsinq = 0.0d+00
            sinqr = 0.0d+00
        else
c-----
c       elastic layer
c-----
            svex = qr
            epp = dcmplx(dcos(pi), dsin(pi))/2.0
            epm = dconjg(epp)
            eqp = dcmplx(dcos(qi), dsin(qi))/2.0
            eqm = dconjg(eqp)
            if(pr.lt.15.) then
                pfac=dexp(-2.*pr)
            else
                pfac  = 0.0d+00
            endif
            cosp = epp + pfac*epm
            sinp = epp - pfac*epm
            rsinp = rp *sinp
            sinpr = sinp/rp

            if(qr.lt.15.) then
                svfac=dexp(-2.*qr)
            else
                svfac  = 0.0d+00
            endif
            cosq = eqp + svfac*eqm
            sinq = eqp - svfac*eqm
            rsinq = rsv*sinq
            sinqr = sinq/rsv

        endif
        return
        end
        subroutine hska(AA,tcosp,trsinp,tsinpr,tcossv,trsinsv,tsinsvr,
     1      NP,NSV,
     1      X11, X21, X31, X41,X12, X22, X32, X42,
     2      TRho,iwat,ex,om2)
c-----
c       Changes
c
c       01 May 2002  - defined cosp = tcosp/NP to 
c             reduce nmber of complex divisions
c-----
        implicit none
        real*8 AA(4,4)
        real*8 tcosp , tcossv 
        real*8 trsinp, trsinsv
        real*8 tsinpr, tsinsvr
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        real TRho
        integer iwat
        real*8 ex, dfac
        real*8  om2
        complex*16 zrho

c-----
c       introduce shorthand to reduce work
c-----
        COMPLEX*16 cosp , cossv 
        COMPLEX*16 rsinp, rsinsv
        COMPLEX*16 sinpr, sinsvr

        integer i, j
        zrho = dcmplx(dble(TRho),0.0d+00)
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,4
                do 101 i=1,4
                    AA(i,j) = 0.0d+00
  101           continue
  100       continue
            if(ex.gt.35.0d+00)then
                dfac = 0.0d+00
            else
                dfac = dexp(-ex)
            endif
            AA(1,1) = dfac
            AA(4,4) = dfac
            AA(2,2) = tcosp
            AA(3,3) = tcosp
            AA(2,3) = -trsinp/(zrho*om2)
            AA(3,2) = - zrho*om2*tsinpr
        else
c-----
c       elastic layer
c-----
            cosp   = tcosp/NP
            sinpr  = tsinpr/NP
            rsinp  = trsinp/NP
            cossv  = tcossv/NSV
            sinsvr = tsinsvr/NSV
            rsinsv = trsinsv/NSV

            AA(1,1) = dreal(  x11*x41*cosp  + x12*x42*cossv  )
            AA(1,2) = dreal(- x11*x31*rsinp - x12*x32*sinsvr )
            AA(1,3) = dreal(- x11*x21*cosp  - x12*x22*cossv  )
            AA(1,4) = dreal(  x11*x11*rsinp + x12*x12*sinsvr )
            AA(2,1) = dreal(  x21*x41*sinpr + x22*x42*rsinsv )
            AA(2,2) = dreal(- x21*x31*cosp  - x22*x32*cossv  )
            AA(2,3) = dreal(- x21*x21*sinpr - x22*x22*rsinsv )
            AA(3,1) = dreal(  x31*x41*cosp  + x32*x42*cossv  )
            AA(3,2) = dreal(- x31*x31*rsinp - x32*x32*sinsvr )
            AA(4,1) = dreal(  x41*x41*sinpr + x42*x42*rsinsv )
            AA(2,4) = - AA(1,3)
            AA(3,3) =   AA(2,2)
            AA(3,4) = - AA(1,2)
            AA(4,2) = - AA(3,1)
            AA(4,3) = - AA(2,1)
            AA(4,4) =   AA(1,1)
        endif
        return
        end

        subroutine normc(ee,ex,nmat)
c       This routine is an important step to control over- or
c       underflow.
c       The Haskell or Dunkin vectors are normalized before
c       the layer matrix stacking.
c       Note that some precision will be lost during normalization.
c
        implicit double precision (a-h,o-z)
        dimension ee(5)
        ex = 0.0d+00
        t1 = 0.0d+00
        do 10 i = 1,nmat
          if(dabs(ee(i)).gt.t1) t1 = dabs(ee(i))
   10 continue
        if(t1.lt.1.d-40) t1=1.d+00
        do 20 i =1,nmat
          t2=ee(i)
          t2=t2/t1
          ee(i)=t2
   20 continue
c-----
c       store the normalization factor in exponential form.
c-----
        ex=dlog(t1)
        return
        end
