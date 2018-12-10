        program joint96
c----------------------------------------------------------------------c
c                                                                      c
c        COMPUTER PROGRAMS IN SEISMOLOGY                               c
c        VOLUME IV                                                     c
c                                                                      c
c        PROGRAM: JOINT96                                              c
c                                                                      c
c        COPYRIGHT 1986, 1991, 2002                                    c
c        D. R. Russell, R. B. Herrmann                                 c
c        Department of Earth and Atmospheric Sciences                  c
c        Saint Louis University                                        c
c        221 North Grand Boulevard                                     c
c        St. Louis, Missouri 63103                                     c
c        U. S. A.                                                      c
c                                                                      c
c----------------------------------------------------------------------c
c       CHANGES
c       07 JUL 2002 - add comments, also wc control to permit mixed
c           smoothing
c       17 JUL 2002 - format change for id.eq.48
c       23 NOV 2002 - introduce RFTN weights and individual 
c          layer specific
c           index of how to get VS - option 30 changed, 49 50 introduced
c       09 JAN 2003 - removed extraneous gzap and pzap subroutines
c       28 JAN 2003 - output files preserve SPHERICAL/FLAT EARTH
c       30 JAN 2003 - Katy Raven at Bullard Laboratories 
c               caught bugs in ww
c               also error in do 1221 for id.eq.36
c       06 JAN 2007 - always run srfdrr96 - else Love only data will not work
c       09 FEB 2011 - srfdrr06 now outputs duda for the Rayleigh wave
c        the jsamat is modified accordingly but the duda is not used
c        in these codes but rather in the shallow96 code of VOLX
c
c
c       Interactive inversion of surface wave dispersion
c         data.  Currently accepts both Rayleigh and Love
c         waves with a maximum of 512 periods, 20 phase modes,
c       and 20 group modes, in any combination, for each
c       type of wave.  Programs used in conjunction with SURF (via
c       ssytem calls) are:
c
c
c          JNTPRE96: Checks input data and sets up unformatted binary
c                  files.
c          SRFPRE96: Checks input data and sets up unformatted binary
c                  files.
c          SRFDIS96: Finds theoretical phase velocities v.s. periods for
c                  both Love and Rayleigh modes.
c          SRFDRL96: Finds group velocities of Love modes and calculates
c                  partial derivatives of phase and group velocities
c                  with respect to layer shear velocities.
c          SRFDRR96: Same as DERIVL, except for Rayleigh modes.
c          SRFINV96: Performs a singular value decomposition of
c                   DERIVR and DERIVL, with either no
c                   or differential smoothing constraints.
c          SRFPHV96:  Plots observed and predicted dispersion curves
c          SRFPHR96: Plot model and resolution kernel
c          RFTNDR96: compute receiver function partials
c          RFTNPV96:  Plots observed and predicted receiver functions
c          SRFPHR96: Plots model and resolution kernel
c
c       Program developed by David R. Russell, St. Louis
c       University, Jan. 1984.  Programs disper.f, derivr.f
c         derivl.f used in ssytem calls of surf.f were developed
c       by R. B. Herrmann and C. Y. Wang, St. Louis University.
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        integer LER, LIN, LOT
        parameter(LER=0,LIN=5,LOT=6)
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)
c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c       NL2 - number of columns in model (first NL2/2 are
c           - velocity parameters, second NL2/2 are Q values)
c-----
        integer nf10(NL)
        common/ctrl/numa,d(NL),a(NL),b(NL),r(NL),rat(NL),dd(NL2),x(NL2),
     1      h(NL2),u(NL2),ct(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        character*20 name
        logical ext
        logical lquiet
        integer mnmarg
        integer ssytem
        integer invdep, lstinv
        save invdep, lstinv
c-----
c       invdep  = 0 last inversion was for depth
c           = 1 last inversion was for velocity and Q inverse
c       lstinv  = 2,3,4 depending on the last inversion
c           invdep = 1 for velocity/Q inversion
c                        = 0 for layer thickness inversion
c-----
        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modtit/title
        character*80 title

        real fval
        integer ival, kval
        logical lval
        real pval
c-----
c       machine dependent initialization
c-----
        call mchdep()
        numa=mnmarg()
c-----
c       if the first argument is 39, then clean up tmp files
c       and exit
c-----
        if(numa .gt. 0)then
            call mgtarg(1,name)
            if(name(1:2).eq.'39')then
                kerr=ssytem('rm -f tmpsrfi.*' )
                kerr=ssytem('rm -f tmpmod96.*')
                kerr=ssytem('rm -f tmpmrg*.*')
                call trmnt()
            endif
        endif
c-----
c
c       Check for existence of temporary files.  If not, assume
c         program is starting and run PRESET.
c-----
        inquire(file='tmpsrfi.00',exist=ext)
        if(.not.(ext)) kerr=ssytem('jntpre96')
c-----
c       get control parameters
c-----
        call gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,nup,
     1      dlam,qaqb,wref,invcsl,lstinv,twnmin,twnmax,iter,
     2      nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        nft = nfilt / 2
        if(iprog.ne.3)then
            WRITE(LOT,*)'joint96 requires a joint inversion'
            if(iprog.eq.1)then
            WRITE(LOT,*)'Do surf96 39 to clean up, then joint96'
            else if(iprog.eq.2)then
            WRITE(LOT,*)'Do rftn96 39 to clean up, then joint96'
            else if(iprog.eq.3)then
            WRITE(LOT,*)'Do joint96 39 to clean up, then joint96'
            endif
            STOP
        endif
        if(itot.lt.0)then
            kerr=ssytem('rm tmpsrfi.*')
            call trmnt()
        endif
c-----
c       get earth model
c-----
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
        m = mmax
c-----
c       get current weighting
c-----
        call ddcur(m,m+m)
c-----
c       parse menu input
c-----
        it=0
c-----
c       if numa > 0 assume command line structure being used.
c-----
        numa=mnmarg()
        if(numa .gt. 0)then
            lquiet = .true.
        else
            lquiet = .false.
        endif
   12   i=0
        it=it+1
   14   if(itot.eq.0) go to 18
c-----
c       terminate iteration loop and reset counter
c-----
        if(it.gt.itot)then
            itot = 0
            call pttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,nup,
     1          dlam,qaqb,wref,invcsl,lstinv,
     2          twnmin,twnmax,iter,nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        endif
c-----
c       never exceed the number of entries on the command line
c-----
        if(it.gt.itot) call trmnt()
  18  i=i+1
        if(numa.eq.0) go to 22
c-----
c       check for number of iterations.
c-----
        if(itot.gt.0.and.i.gt.numa) go to 12
        if(itot.eq.0.and.i.gt.numa) call trmnt()
c-----
c       Read in flag (id) of desired command, either from
c       command line or from terminal.
c-----
        call mgtarg(i,name)
        read(name,'(i2)') id
        if(id.ge.0.and.id.le.50) go to 24
        write(LOT,*)'command string not valid - job aborted'
        call trmnt()
  22    if(i.eq.1) call menu()
        write(LOT,*)'ready'
        read(LIN,*) id
c-----
c       execute command based on id.
c-----
   24   continue
        if(id.eq.1)then
c-----
c           Calculate theoretical values and
c               partial derivatives.
c-----
            kerr=ssytem('srfdis96')
            open(2,file='tmpsrfi.08')
            close(2,status='delete')
c-----
c           The order of these next two statements is
c           ESSENTIAL, since srfdrr96 merges with the
c           Love wave partials if they are computed
c-----
            if(nf34.ne.0) kerr=ssytem('srfdrl96')
c-----
c     06 JAN 2007 - always run srfdrr96
c-----
C            if(nf67.ne.0) kerr=ssytem('srfdrr96')
            kerr=ssytem('srfdrr96')
            open(1,file='tmpsrfi.16',status='unknown',
     1          form='formatted',access='sequential')
            rewind 1
            write(1,'(i5,2f20.10,i5)')invdep,twnmin,twnmax,idtwo
            close (1)
            kerr=ssytem('rftndr96')
            nup=1
            call pttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,nup,
     1          dlam,qaqb,wref,invcsl,lstinv,
     2          twnmin,twnmax,iter,nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        else if(id.eq.2.or.id.eq.3.or.id.eq.4)then
c-----
c       safety - do not permit more than one update to model
c-----
            if(nup.gt.1 .and. id.ne.28 .and. id.ne.23 ) then
                if(nup.eq.2)then
                    write(LOT,*)'Model has been updated:'
                    write(LOT,*)
     1          'Run:  (1,2) inversion'
                endif
            else
                lstinv = id
c-----
c               Calculate singular value decomposition.
c-----
                call jramat(id,nfilt,i,wref,invcsl,nf1,
     1              invdep,lstinv,sigr)
                call jsamat(id,nfilt,i,wref,invcsl,nf1,
     1              invdep,lstinv,sigv,sigg)
                call jamat(id,nfilt,i,wref,invcsl,nf1,
     1              invdep,lstinv,pval,sigv,sigr,sigg,
     1              nurftn)
CRFTN96         call ramat(id,nfilt,i,wref,invcsl,nf1,
CRFTN96     1               invdep,lstinv)
                kerr=ssytem('srfinv96')
                if(nup.eq.1) then
                    nup=0
                endif
                call pttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,
     1              nfilt,nup,dlam,qaqb,wref,invcsl,
     2              lstinv,twnmin,twnmax,iter,nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
            endif
        else if(id.eq.6 .or. id.eq.10 .or. id.eq.18 .or. id.eq.19
     1      .or. id.eq.28 .or. id.eq.29 .or. id.eq.23
     2      .or. id.eq.13 .or. id.eq.14 .or. id.eq.24)then
c-----
c           update shear velocities for this iteration, or,
c           list singular value decomposition, or,
c           list shear velocity values, or,
c           list resolving kernels, or,
c           create file with compressional, shear velocities,
c           standard deviations of shear velocities, densities,
c           Poisson's ratio, and resolution information.
c-----
            if(nup.gt.0 .and. id.ne.28 .and. id.ne.23 ) then
                if(nup.eq.2)then
                    write(LOT,*)'Model has been updated:'
                    write(LOT,*)
     1              'Run: (1,2) inversion'
                else if(nup.eq.1)then
                    write(LOT,*)'Model partials computed'
                    write(LOT,*)
     1              'Run: (2) inversion of last update'
                endif
            else
c-----
c       if model is to be updated, it is only done once
c       for any velocity update, and for Q update for
c       coupled inversion
c
c       for uncoupled, we can update q model as many times
c       as desired to make model converge. This is because the
c       necessary partials with respect to layer velocity
c       do not change
c-----
                if( id.eq.6 )then
                    if(lstinv.ne.3 .or. invcsl.gt.1)then
                        nup=2
                    endif
                    iter = iter + 1
                    call modls(6,i,nf1,nf10,dlam,invdep,
     1                  lstinv,iter,rmsv,rmsq)
                    call pttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,
     1              nf10,nfilt,nup,dlam,qaqb,wref,
     2              invcsl,lstinv,twnmin,twnmax,iter,
     3              nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
                    if(invdep.eq.0)then
                        WRITE(LOT,'(a,f9.4,a)')
     1      'RMS change in layer thickness model    :',rmsv, ' km'
                        write(LOT,'(a,i5,a)')
     1                      'ITERATION ',iter,
     1                      ' done: UPDATING H'
                    else
                        if(lstinv.eq.2)then
                        WRITE(LOT,'(a,f9.4,a)')
     1      'RMS change in S-wave velocity model    :',rmsv, ' km/sec'
                        write(LOT,'(a,i5,a)')
     1                      'ITERATION ',iter,
     1                      ' done: UPDATING V'
                        else if(lstinv.eq.3)then
                        WRITE(LOT,'(a,g9.2,a)')
     1      'RMS change in QS inverse      model    :',rmsq, ' '
                        write(LOT,'(a,i5,a)')
     1                      'ITERATION ',iter,
     1                      ' done: UPDATING Q'
                        else if(lstinv.eq.4)then
                        WRITE(LOT,'(a,f9.4,a)')
     1      'RMS change in S-wave velocity model    :',rmsv, ' km/sec'
                        WRITE(LOT,'(a,g9.2,a)')
     1      'RMS change in QS inverse      model    :',rmsq, ' '
                        write(LOT,'(a,i5,a)')
     1                      'ITERATION ',iter,
     1                      ' done: UPDATING V & Q'
                        endif
                    endif
                    write(LOT,1)
    1   format('-------------------------------------',
     1  '---------------------------------------')
                else
                    call modls(id,i,nf1,nf10,dlam,invdep,
     1                  lstinv,iter,rmsv,rmsq)
                endif
            endif
        else if(id.eq.11 .or. id.eq.12 .or. id.eq.16 .or. id.eq.17
     1      .or. id.eq.22 .or. id.eq.27)then
c-----
c           list dispersion values or partial derivatives.
c-----
            call jsamat(id,nfilt,i,wref,invcsl,nf1,invdep,
     1          lstinv,sigv,sigg)
        else if(id.eq.31)then
c-----
c       modify model parameter weighting
c-----
            call ddcur(m,m2)
            if(numa.ne.0)then
                i=i+1
                name = ' '
                call mgtarg(i,name)
                read(name,'(i10)')ival
                i=i+1
                name=' '
                call mgtarg(i,name)
                read(name,'(f10.0)')fval
                if(fval.lt.0.0001)fval=0.0001
            else
                write(LOT,'(i2,a,i2,a,i2,a1,i2,a)')m
     1              ,' layers: 1-',m,
     1              ' for Vs ',m,'-',m2,' for Qbinv'
                write(LOT,*)'Enter i'
                read(LIN,*)ival
                call ddcur(m,m2)
                if(ival.lt.1.or.ival.gt.m2)then
                    write(LOT,*)'ival out of range'
                else
                write(LOT,*)'Current dd(',ival,')=',dd(ival)
                write(LOT,*)'Enter New dd(',ival,')'
                read(LIN,*)fval
                if(fval.lt.0.0001)fval=0.0001
                endif
            endif
            if(ival.gt.0 .and. ival.le.m2)then
                call ddupdt(ival,fval)
            endif
        else if(id.eq.7)then
c-----
c           call partials data 
c           velocity model plotting program
c-----
            if(nup.eq.2  ) then
                write(LOT,*)'model has been updated:'
                write(LOT,*)
     1          'run: (1) for new partials'
            else
                kerr=ssytem('srfvph > /dev/null 2>&1' )
                kerr=ssytem('rftnvp > /dev/null 2>&1' )
            endif
        else if(id.eq.8)then
c-----
c           call gamma Qb inverse plotting program
            if(nup.eq.2  ) then
                write(LOT,*)'model has been updated:'
                write(LOT,*)
     1          'run: (1) for new partials'
            else
                kerr=ssytem('srfgph > /dev/null 2>&1')
            endif
        else if(id.eq.9)then
c-----
c           call model resolution program
c-----
            if(nup.eq.1  ) then
                write(LOT,*)'new partials have been found:'
                write(LOT,*)
     1          'run: (2) inversion '
            else if(nup.eq.2)then
                write(LOT,*)
     1          'run: (1,2) inversion'
            else
                call modls(id,i,nf1,nf10,dlam,invdep,
     1              lstinv,iter,rmsv,rmsq)
                kerr=ssytem('srfrph > /dev/null 2>&1')
            endif
        else if(id.eq.38)then
c-----
c           temporarily end program, leaving all temporary
c           files in place
c-----
            call trmnt()
        else if(id.eq.39)then
c-----
c           permanently end program, removing all tmeporary files
c-----
            kerr=ssytem('rm -f tmpsrfi.*')
            kerr=ssytem('rm -f tmpmod96.*')
            kerr=ssytem('rm -f tmpmrg*.*' )
            call trmnt()
        else if(id.eq.0)then
            call menu()
        else if(id.eq.5)then
            call iquery(
     1      'Invert for Velocity(1)--Depth(0): Currently is',
     2          invdep,lquiet,i)    
            call ibound(0,1,invdep,1)
        else if(id.eq.30)then
c-----
c       modify how VP is computed in each layer
c-----
            if(numa.ne.0)then
c-----
c               command line 30 LAYER VP_compute
c-----
                i=i+1
                name = ' '
                call mgtarg(i,name)
                read(name,'(i10)')ival
                i=i+1
                name=' '
                call mgtarg(i,name)
                read(name,'(i10)')jval
                if(ival.ge.1.or.ival.le.m)then
                    nf10(ival) = jval
                    call ibound(0,1,nf10(ival),1)
                endif
            else
                write(LOT,'(i2,a,i2,a)')m
     1              ,' layers: 1-',m,
     1              ' Vp fixed(0), Vp/Vs fixed(1)'
                write(LOT,*)'Enter i'
                read(LIN,*)ival
                if(ival.lt.1.or.ival.gt.m)then
                    write(LOT,*)'ival out of range'
                else
                    write(LOT,*)'Current nf10(',ival,')='
     1                  ,nf10(ival)
                    write(LOT,*)'Enter New nf10(',ival,')'
                    read(LIN,*)jval
                    nf10(ival) = jval
                    call ibound(0,1,nf10(ival),1)
                endif
            endif
            if(ival.gt.0 .and. ival.le.m)then
                call pttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,
     1              nfilt,nup,dlam,qaqb,wref,invcsl,
     2              lstinv,twnmin,twnmax,iter,nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3              idtwo,idum2,idum3,idum4,idum5,
     4              rdum1,rdum2,rdum3,rdum4,rdum5)
            endif
        else if(id.eq.32)then
            call fquery(
     1      'enter damping factor: Currently is',
     2          dlam, lquiet, i)
            call fbound(0.0, 10000.0, dlam, 1.0)
        else if(id.eq.33)then
            call fquery(
     1      'enter minimum RFTN time window (<=0) : twnmin=',
     2          twnmin,lquiet,i)
            call fbound(-1000.0,0.0,twnmin,-5.0)
        else if(id.eq.34)then
            call fquery(
     1      'enter maximum RFTN time window (> 0) : twnmax=',
     2          twnmax,lquiet,i)
            call fbound(0.0,1000.0,twnmax,+20.0)
        else if(id.eq.35)then
            ival = invcsl
            call iquery(
     1      'non-causal(0),uncpled cau(1),coup caus(2): invcsl=',
     2      invcsl, lquiet, i)
            call ibound(0,2,invcsl,0)
        else if(id.eq.36)then
            call iquery(
     1      'Reset global smoothing none(0),diff(1): Currently is',
     2          nft,lquiet,i)
            call ibound(0,1,nft,1)
            if(nft.eq.0)then
            do 1221 ii=1,m
                call  wcupdt(ii,.false.)
                call  wcupdt(ii+m,.false.)
 1221       continue
            else if(nft.eq.1)then
            do 1222 ii=1,m-1
                call  wcupdt(ii,.true.)
                call  wcupdt(ii+m,.true.)
 1222       continue
            endif
        else if(id.eq.37)then
            call iquery(
     1      'enter number of iterations: itot=',
     2          itot,lquiet,i)
            call ibound(0,1000,itot,0)
        else if(id.eq.40)then
            call fquery(
     1      'Enter minimum velocity dispersion std error',
     2          sigv,lquiet,i)
            call fbound(0.0,1000.0,sigv,1.0)
        else if(id.eq.41)then
            call fquery(
     1      'Enter minimum gamma dispersion std error',
     2          sigg,lquiet,i)
            call fbound(0.0,1000.0,sigg,1.0)
        else if(id.eq.42)then
            call fquery(
     1      'Enter minimum RFTN  std error',
     2          sigr,lquiet,i)
            call fbound(0.0,1000.0,sigr,1.0)
        else if(id.eq.43)then
            call fquery(
     1      'Enter Joint Control (p): RFTN=0 <= p <= 1=SURF',
     2          pval,lquiet,i)
            call fbound(0.0,1.0,pval,0.5)
        else if(id.eq.44)then
            call iquery(
     1      '2x RFTN computation: (0)no (1)yes',
     2          idtwo,lquiet,i)
            call ibound(0,1,idtwo,0)
        else if(id.eq.45)then
            call ddcur(mmax,mmax+mmax)
            call shwwtv(mmax,dd,wc,nf10)
        else if(id.eq.46)then
            call ddcur(mmax,mmax+mmax)
            call shwwtq(mmax,dd,wc)
        else if(id.eq.47)then
            call shwctl()
        else if(id.eq.48)then
c-----
c       Select layer smoothing condition
c       enter S s for smooth or D d for differential
c-----
            call ddcur(m,m2)
            if(numa.ne.0)then
                i=i+1
                name = ' '
                call mgtarg(i,name)
                read(name,'(i10)')ival
                i=i+1
                name=' '
                call mgtarg(i,name)
                read(name,'(i10)')kval
                if(kval.eq.0)then
                    lval = .false.
                else if(kval.eq.1)then
                    lval = .true.
                endif
            else
                write(LOT,'(i2,a,i2,a,i2,a1,i2,a)')m
     1              ,' layers: 1-',m,
     1              ' for Vs ',m+1,'-',m2,' for Qbinv'
                write(LOT,*)'Enter layer/boundary i'
                read(LIN,*)ival
                call ddcur(m,m2)
                if(ival.lt.1.or.ival.gt.m2)then
                    write(LOT,*)'ival out of range'
                else
                write(LOT,*)'Current Smoothing wc(',ival,')=',
     1                  wc(ival)
                write(LOT,*)'Enter New wc(',ival,')',
     1          ':0 no smoothing, 1 smoothing '
                read(LIN,*)kval
                if(kval.eq.0)then
                    lval = .false.
                else if(kval.eq.1)then
                    lval = .true.
                endif
                endif
            endif
            if(ival.gt.0 .and. ival.le.m2)then
                call wcupdt(ival,lval)
            endif
        else if(id.eq.49)then
            call shwrfw()
        else if(id.eq.50)then
            call chgrfw(numa,i)
        endif
            nfilt = mod(nfilt,2) + nft*2
            call pttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,nup,
     1          dlam,qaqb,wref,invcsl,lstinv,
     1          twnmin,twnmax,iter,nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        go to 14
        end

        subroutine trmnt()
c-----
c       subroutine to exit program
c-----
            call exit(0)
        end

        function lnstr(str)
        character*80 str
c-----
c       knowing that a string is at most 80 characters long,
c       determine the minimum length of string by stripping
c       characters off the end. Return 1 at the minimum
c-----
        do 100 i=80,1,-1
            lnstr = i
            if(str(i:i).ne.' ')goto 101
  100   continue
  101   continue
        return
        end

        subroutine gtrofa(r,a)
c-----
c       obtain density (r,gm/cc) as a function of P-vel (a,km/sec)
c-----
        dimension rp(18)
        data rp/1.65,1.85,2.05,2.15,2.23,2.32,2.39,2.50,2.60,2.70,2.85,
     # 2.98,3.14,3.31,3.49,3.66,3.86,4.05/
c-----
c determine density
c-----
        if(a.gt.1.5 .and. a.lt.10.0)then
            xa = (a-1.5)*2.
            na = xa+1
            dr = rp(na+1)-rp(na)
            da=a-na*.5-1.0
            r = rp(na)+dr*da*2.
        else if(a.le.1.5)then
            r = 1.0
        else if(a.ge.10.0)then
            r = 4.10
        else if(a.le.1.5)then
            r = rp(1)
        endif
        return
        end

        subroutine gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup,dlam,qaqb,wref,invcsl,
     2      lstinv,twnmin,twnmax,iter,nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        integer NL
        parameter (NL=200)
        integer nf10(NL)
c-----
c       read control file
c-----
        open(1,file='tmpsrfi.00',form='unformatted',access='sequential')
        rewind 1
        read(1) iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,nup,dlam
     1      ,qaqb,wref,invcsl,lstinv,twnmin,twnmax,iter,
     2      nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5
        close(1,status='keep')
        return
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
c       update control file
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

        subroutine ddupdt(ival,fval)
c-----
c       update layer weighting value in the array dd(NL2)
c       by reading the tmpsrfi.12 data file completely,
c       and writing completely to tmprsrf.14, changing
c       only dd(ival) = fval, and then copying back
c       to tmpsrfi.12 for use by the subroutine amat.f
c-----
        integer NL, NL2
        parameter(NL=200,NL2=NL+NL)
        common/ctrl/numa,d(NL),a(NL),b(NL),r(NL),rat(NL),dd(NL2),x(NL2),
     $      h(NL2),u(NL2),ct(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        open(1,file='tmpsrfi.12',form='unformatted',access='sequential'
     1      ,status='unknown')
        rewind 1
c-----
c       tmpsrfi.01 is used as a temporary storage file
c       to permit updating tmpsrfi.12
c-----
        open(2,file='tmpsrfi.01',form='unformatted',access='sequential'
     1      ,status='unknown')
        rewind 2
        read(1) nd,m
        write(2) nd,m
        m2 = m + m
        read(1)(dd(i),i=1,m2)
        read(1)(wc(i),i=1,m2)
c-----
c       update the dd(ival) value
c-----
        dd(ival)=fval
        dd(ival+m)=fval
        write(2)(dd(i),i=1,m2)
        write(2)(wc(i),i=1,m2)
C 1000  continue
C       read(1,end=1001)ifn,k,md,tp,c,sd
C       write(2)ifn,k,md,tp,c,sd
C       go to 1000
C 1001  continue
        rewind 1
        rewind 2
        read(2)nd,m
        write(1)nd,m
        read(2)(dd(i),i=1,m2)
        write(1)(dd(i),i=1,m2)
        read(2)(wc(i),i=1,m2)
        write(1)(wc(i),i=1,m2)
C 2000  continue
C       read(2,end=2001)ifn,k,md,tp,c,sd
C       write(1)ifn,k,md,tp,c,sd
C       goto 2000
C 2001  continue
        close(2,status='delete')
        close(1)
        return
        end

        subroutine wcupdt(ival,lval)
c-----
c       update layer weighting value in the array dd(NL2)
c       by reading the tmpsrfi.12 data file completely,
c       and writing completely to tmprsrf.14, changing
c       only dd(ival) = fval, and then copying back
c       to tmpsrfi.12 for use by the subroutine amat.f
c-----
        integer ival
        logical lval

        integer NL, NL2
        parameter(NL=200,NL2=NL+NL)
        common/ctrl/numa,d(NL),a(NL),b(NL),r(NL),rat(NL),dd(NL2),x(NL2),
     $      h(NL2),u(NL2),ct(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        open(1,file='tmpsrfi.12',form='unformatted',access='sequential'
     1      ,status='unknown')
        rewind 1
        open(2,file='tmpsrfi.01',form='unformatted',access='sequential'
     1      ,status='unknown')
        rewind 2
        read(1) nd,m
        write(2) nd,m
        m2 = m + m
        read(1)(dd(i),i=1,m2)
        read(1)(wc(i),i=1,m2)
c-----
c       update the dd(ival) value
c-----
        wc(ival)=lval
        wc(ival+m)=lval
        write(2)(dd(i),i=1,m2)
        write(2)(wc(i),i=1,m2)
C 1000  continue
C       read(1,end=1001)ifn,k,md,tp,c,sd
C       write(2)ifn,k,md,tp,c,sd
C       go to 1000
C 1001  continue
        rewind 1
        rewind 2
        read(2)nd,m
        write(1)nd,m
        read(2)(dd(i),i=1,m2)
        write(1)(dd(i),i=1,m2)
        read(2)(wc(i),i=1,m2)
        write(1)(wc(i),i=1,m2)
C 2000  continue
C       read(2,end=2001)ifn,k,md,tp,c,sd
C       write(1)ifn,k,md,tp,c,sd
C       goto 2000
C 2001  continue
        close(2,status='delete')
        close(1)
        return
        end

        subroutine ddcur(m,m2)
c-----
c       get array of current dd values and array limit
c-----
        parameter(NL=200,NL2=NL+NL)
        common/ctrl/numa,d(NL),a(NL),b(NL),r(NL),rat(NL),dd(NL2),x(NL2),
     $      h(NL2),u(NL2),ct(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        open(1,file='tmpsrfi.12',form='unformatted',access='sequential'
     1      ,status='unknown')
CSURF96 open(1,file='tmpsrfi.04',form='unformatted',access='sequential'
CSURF96     1       ,status='unknown')
        rewind 1
        read(1) nd,m
        m2 = m + m
        read(1)(dd(i),i=1,m2)
        read(1)(wc(i),i=1,m2)
        close(1,status='keep')
        return
        end

