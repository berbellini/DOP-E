        program rftnpr96
c---------------------------------------------------------------------:
c                                                                   c
c    COMPUTER PROGRAMS IN SEISMOLOGY                                c
c    VOLUME IV                                                      c
c                                                                   c
c    PROGRAM: RFTNPR96                                              c
c                                                                   c
c    COPYRIGHT 1986, 1991, 2001                                     c
c    D. R. Russell, R. B. Herrmann                                  c
c    Department of Earth and Atmospheric Sciences                   c
c    Saint Louis University                                         c
c    221 North Grand Boulevard                                      c
c    St. Louis, Missouri 63103                                      c
c    U. S. A.                                                       c
c                                                                   c
c---------------------------------------------------------------------c
c
c     This program checks the input control file 'robs.d' and
c     converts the input data into unformatted binary files
c     to be used in other programs.  The unformatted files
c     are labeled as 'tmpsrfi.xx' where xx is a number from
c   0 to 14.
c
c     Developed by David R. Russell, St. Louis University, Jan. 1984.
c     Restructured Input R. B. Herrmann 8 August 1986
c
c     Restructured August 1986 to simplify input format RBH
c     Modified  September 1986 to include gamma values  RBH
c     Modified  November 2001 to use SURF96 dispersion format and
c     model96 earth model format
c     Modified  January 2002 for clean up 
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     CHANGES
c     23 NOV 2002 - set up default weights of 1.0 for 
c          each receiver function
c-----

        parameter(LER=0,LIN=5,LOT=6)
        parameter(NL=200,NL2=NL+NL)
c-----
c     LIN   - unit for FORTRAN read from terminal
c     LOT   - unit for FORTRAN write to terminal
c     LER   - unit for FORTRAN error output to terminal
c     NL    - number of layers in model
c     NL2   - number of columns in model (first NL2/2 are
c         - velocity parameters, second NL2/2 are Q values)
c-----
        integer nf(13)
        real dd(NL2)
        logical wc(NL2)
        character*80 nmmodl, nmrftn
        logical ext
        common/param/qaqb,itype,dlam,invdep

        integer nf10(NL)
        data nf10/NL*1/
c-----
c     machine dependent initialization
c-----
        call mchdep()
c-----
c     test for the existence of robs.d
c     if it does not exist interactively set it up
c-----
        inquire(file='robs.d',exist=ext)
        if(.not.(ext))    call setrobs()
c-----
        open(1,file='robs.d',status='old',access='sequential')
        rewind 1
c-----
c     nf() is the control array
c     nf(1) = 1 estimated stdev computed from residuals
c           0 no scaling by residuals
c     nf(8) = Model Weighting
c           0 No Weighting
c           1 Read In Layer Velocity Weights
c     nf(9) = Number of Layers in Starting Model (from model file)
c     nf(10)= Input Format (from model file)
c           0 - Inversion a, rho fixed
c           1 - Inversion Poisson Ratio Fixed, Rho computed from Vp
c     nf(11)= Type of Smoothing
c           0 None
c           1 Differential
c     nf(12)= Earth Flattening
c           0 No Earth Flattening - Model is Flat Earth Model
c           1 Earth Flattening - Model is Spherical Earth Model
c           1 frequency
c-----
        read(1,*) (nf(i),i=1,8),nf(11),nf(12)
C       write(LOT,*) (nf(i),i=1,8),nf(11),nf(12)

        read(1,'(a)')nmmodl
C       write(LOT,*)nmmodl
c-----
c     nmmodl= name of file containing model information
c-----
        call getmdl(nmmodl,nf(9),nf(10),iunit,iflsph)
        nf(12) = iflsph
        m=nf(9)
        m2 = m + m
c-----
c     nmrftn= name of file containing list of receiver 
c          function file names
c     nurftn= number of receiver functions to be inverted
c     Open this file, then read each SAC file to insure 
c          that it is a correct SAC
c     file
c-----
        read(1,'(a)')nmrftn
C       write(LOT,*)nmrftn
        call chksac(nmrftn,nurftn)
c-----
c     initialize model weights
c-----
        do 2 i=1,m
            dd(i) = 1.0
            dd(i+m) = dd(i)
            wc(i) = .true.
            wc(i+m) = .true.
            if(i.eq.m)then
                wc(i) = .false.
                wc(i+m) = .false.
            else
                wc(i) = .true.
                wc(i+m) = .true.
            endif
    2   continue
        close (1)
        open(8,file='tmpsrfi.12',form='unformatted',access='sequential')
        rewind 8
        write(8) nf(8),m
        write(8)(dd(i),i=1,m2)
        write(8)(wc(i),i=1,m2)
        nf34=max(nf(3),nf(4))
        nf67=max(nf(6),nf(7))
        close(8,status='keep')
        if(nf(8).eq.0.and.nf(11).eq.0) nfilt=0
        if(nf(8).eq.1.and.nf(11).eq.0) nfilt=1
        if(nf(8).eq.0.and.nf(11).eq.1) nfilt=2
        if(nf(8).eq.1.and.nf(11).eq.1) nfilt=3
        nf34 = max(nf(3),nf(4))
        nf67 = max(nf(6),nf(7))
c-----
c     this value of nup forces a run of dispersion first
c     Also set up default internal values that can be overridden by menu
c-----
        nup=2
        nzer=0
        qaqb = 2.25
        invcsl = 0
        wref = 1.0
        lstinv = -1
        invdep = 1
        itype = 0
        dlam = 1.0
        twnmin = -5.0
        twnmax = 20.0
        iter = 0
        pval = 0.5
        sigv = 0.05
        sigr = 0.05
        sigg = 0.00005
        idtwo = 0
        idum2 = 0
        idum3 = 0
        idum4 = 0
        idum5 = 0
        rdum1 = 0.0
        rdum2 = 0.0
        rdum3 = 0.0
        rdum4 = 0.0
        rdum5 = 0.0
c-----
c     iprog is a binary OR: 2= rftn, 1=surf
c-----
        iprog = 2 + 0
        call pttmp0(iprog,nzer,nf(1),nf(2),nf34,nf(5),nf67,nf10,nfilt,
     1            nup,dlam,qaqb,wref,invcsl,lstinv,
     2      twnmin,twnmax,iter,nurftn,invdep,pval,sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        end

        subroutine pttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,
     1      nfilt,nup,dlam,qaqb,wref,invcsl,lstinv,
     2      twnmin,twnmax,iter,nurftn,invdep,pval,sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        integer NL
        parameter (NL=200)
        integer nf10(NL)
c-----
c     update control file
c-----
        open(1,file='tmpsrfi.00',form='unformatted'
     1            ,access='sequential')
        rewind 1
        write(1) iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup,dlam,qaqb,wref,invcsl,lstinv,
     2      twnmin,twnmax,iter,nurftn,invdep,pval,sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5
        close(1,status='keep')
        return
        end

        subroutine getmdl(nmmodl,mmax,nfmod,iunit,iflsph)
c-----
c     igetmod common
c-----
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)

        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        character nmmodl*(*)
        common/modtit/title
        character*80 title
c-----
c     open the model file. These values will be saved
c     in internal files so that the original model is not
c     modified in any way
c-----
        call getmod(2,nmmodl,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
c----
c     mmax    nf9
c     
c     nfmod   nf10    1 invert for S - P, Rho calculated from S
c             0 invert for S but P, Rho fixed
c-----
        nfmod = 1
        iunit = 0
        LT = LGSTR(TITLE)
        call putmod(2,'tmpsrfi.17',mmax,title(1:lt),iunit,iiso,iflsph,
     1      idimen,icnvel,.false.)
        call putmod(2,'tmpmod96.000',mmax,title(1:lt),iunit,iiso,iflsph,
     1      idimen,icnvel,.false.)
        return
        end

        subroutine setrobs()
c-----
c     set up the control file robs.d interactively
c-----
        parameter (LIN=5, LOT=6, LER=0)
        parameter (NL=200,NL2=NL+NL)
        character*50 nmmodl, nmrftn
        integer nf(12)
        logical ext
        data nf/12*0/
c-----
c     open the file robs.d
c-----
        open(1,file='robs.d',status='unknown',form='formatted',
     1            access='sequential')
        rewind 1
c-----
c     get control for dispersion search
c-----
c-----
c     get second line of robs.d
c-----
        write(LOT,*)
C     1 ' Enter 1 if SW variance based on residuals or'
C       write(LOT,*)
C     1 '       0 if SW variance based on obs std err'
C       read(LIN,*)nf(1)
        nf(1) = 1
C       write(LOT,*)
C     1 ' Model Layer Weighting'
C       write(LOT,*)
C     1 '      0 No weighting'
C       write(LOT,*)
C     1 '      1 Read in weights for layers'
C       read(LIN,*)nf(8)
        nf(8) = 0
C       write(LOT,*)' Enter inversion technique'
C       write(LOT,*)'    0   invert for Vs :Va,rho fixed'
C       write(LOT,*)'    1 : invert for Vs :Poisson fixed, rho from Vp'
C       read(LIN,*)nf(10)
C       write(LOT,*)
C     1 ' Enter Type of Smoothing'
C       write(LOT,*)
C     1 '      0 None'
C       write(LOT,*)
C     1 '      1 Differential'
C       read(LIN,*)nf(11)
        nf(11) = 1
c-----
c     nf(12) on flat/spherical model is now build into the 
c          earth model file
c-----
        nf(12) = 0
        write(1,'(10i5)')(nf(i),i=1,8),nf(11),nf(12)
c-----
c     get model file information
c-----
        write(LOT,*)
     1 ' Enter name of model file'
        read(LIN,'(a)')nmmodl
        write(1,'(a)')nmmodl
c-----
c     see whether the model file exists
c-----
        inquire(file=nmmodl,exist=ext)
c-----
c     if it does not, interactively create it
c-----
        if(.not.(ext))call setmod(nmmodl,mmax)
c-----
c     get receiver function file
c-----
        write(LOT,*)
     1 ' Enter name of receiver function  file list'
        read(LIN,'(a)')nmrftn
        write(1,'(a)')nmrftn
c-----
c     if the file does not exist, then interactively build it
c-----
        inquire(file=nmrftn,exist=ext)
        if(.not.(ext))call setrft(nmrftn)
        return
        end

        subroutine setrft(nmrftn)
c-----
c     interactively set up dispersion file
c-----
        parameter (LIN=5, LOT=6, LER=0)
        character nmrftn*(*)
        character fname*80

        write(LOT,*)
     1 ' Interactively setting up receiver function file list:',nmrftn
c-----
c     open file
c-----
        open(2,file=nmrftn,status='unknown',form='formatted',
     1           access='sequential')
        rewind 2
        write(LOT,*)
     1 'Enter receiver function SAC binary file name, EOF to end'
 1000   continue
        read(LIN,'(a)',end=1001,err=1001)fname
            lf = lgstr(fname)
            write(2,'(a)') fname(1:lf)
        goto 1000
 1001   continue
        close (2)
        return
        end

        subroutine perr(str)
        character str*(*)
        parameter(LER=0,LIN=5,LOT=6)
        dimension nf(10)
        common/param/qaqb,itype,dlam,invdep
        integer NL
        parameter (NL=200)
        integer nf10(NL)
        data nf10/NL*1/
        data nf/10*0/
            nzer=-1
            nf34 = 0
            nfilt = 0
            nup=0
            wref = 1.0
            invcsl = 0
            lstinv = -1
            twnmin = -5.0
            twnmax =  20.0
            iter = 0 
            nurftn = 0 
            invdep = 1
            pval = 0.5
            iprog = 2
        sigv = 0.05
        sigr = 0.05
        sigg = 0.00005
        idtwo = 0
        idum2 = 0
        idum3 = 0
        idum4 = 0
        idum5 = 0
        rdum1 = 0.0
        rdum2 = 0.0
        rdum3 = 0.0
        rdum4 = 0.0
        rdum5 = 0.0
        write(LER,*)str
c-----
c     iprog is a integer  2= rftn, 1=surf 3 = joint
c-----
        iprog = 2 + 0
        call pttmp0(iprog,nzer,nf(1),nf(2),nf34,nf(5),nf67,nf10,
     1      nfilt,nup,
     1      dlam,qaqb,wref,invcsl,lstinv,
     2      twnmin,twnmax,iter,nurftn,invdep,pval,sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        stop
        end


        subroutine chksac(nmrftn,nurftn)
c-----
c     check the SAC binary files to see if they have all the information
c     required for the receiver function computation
c-----
        implicit none
        character nmrftn*80, fname*80
        integer nurftn

        integer LER, LIN, LOT
        parameter(LER=0,LIN=5,LOT=6)
        integer NPTS
        parameter (NPTS=16384)
        real x(NPTS)
        integer n, nerr
        real dt, rayp, gauss, delay, invwgt
        integer lf
        integer lgstr

        integer nzyear, nzjday, nzhour, nzmin, nzmon, nzday
        character*8 kstnm
c-----
c     open control file with list of receiver function names
c-----
        open(2,file=nmrftn,status='unknown',form='formatted',
     1           access='sequential')
        rewind 2
c-----
c     open receiver function control file
c-----
        open(3,file='tmpsrfi.15',access='sequential',form='formatted',
     1      status='unknown')
        rewind 3
        nurftn = 0
 1000   continue
        read(2,'(a)',end=9000)fname
        lf = lgstr(fname)
        call brsac(4,NPTS,fname,x,nerr)
        if(nerr.lt.0)then
            WRITE(LER,*)'Error reading SAC file nerr=',nerr
            lf = lgstr(fname)
            WRITE(LER,*)'SAC FILE:',fname(1:lf)
            if(nerr .eq.-2)then
            WRITE(LER,*)'Number of data points in file exceeds ',NPTS
            endif
            go to 1000
        endif
        call getnhv('NPTS    ',n,nerr)
        call getfhv('DELTA   ',dt,nerr)
        call getfhv('USER0   ',gauss,nerr)
        call getfhv('USER4   ',rayp,nerr)
        call getfhv('B       ',delay,nerr)
        call getnhv('NZYEAR  ',nzyear,nerr)
        call getnhv('NZJDAY  ',nzjday,nerr)
        call getnhv('NZHOUR  ',nzhour,nerr)
        call getnhv('NZMIN   ',nzmin ,nerr)
        call getkhv('KSTNM   ',kstnm ,nerr)
        call mnthdy(nzyear,nzjday,nzmon,nzday)
        if(delay.ne. -12345.)then
            delay = - delay
        endif
        if(n.le.0 .or. dt.le.0.0 .or. gauss.le.0.0 .or. rayp.lt.0.0
     1      .or. delay.lt.0.0)then
            write(LER,*)'Problems with RFTN SAC FILE:',fname(1:lf)
            write(LER,*)'NPTS       :',n
            write(LER,*)'DT         :',dt
            write(LER,*)'Gauss USER0:',gauss
            write(LER,*)'Rayp  USER4:',rayp
            write(LER,*)'Delay B    :',delay
        else
            write(3,'(a)')fname(1:lf)
            write(3,'(a)')kstnm
            write(3,'(i4,1x,i3,1x,i2,1x,i2,1x,i2,1x,i2)')
     1          nzyear, nzjday, nzhour, nzmin, nzmon, nzday
            if(gauss.gt.0.0)then
                invwgt = 1.0/gauss
            else
                invwgt = 1.0
            endif
            write(3,'(i5,f10.3,f10.3,f10.3,f10.3,f10.3)')
     1          n,dt,rayp,gauss,delay,invwgt
            nurftn = nurftn + 1
        endif
        
        go to 1000
 9000   continue
        close (3)
        close (2)
        return
        end

