      program srfphv96
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME IV                                                       c
c                                                                      c
c      PROGRAM: SRFPHV96                                               c
c                                                                      c
c      COPYRIGHT 1986, 2002                                            c
c      D. R. Russell, R. B. Herrmann                                   c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c-----
c       CHANGES
c       21 MAY 2002 Corrected scaling for depth - previously
c               if layer thickness inversion used, only
c               the second model thicknesses were used
c       10 JUL 2002 Updated to include wc array for layer c
c               inversion control
c       26 APR 2003 Modified plot to place predicted ON TOP of preserved
c       10 FEB 2011 Change to handle the new tmpsrfi.08 format which
c                   has duda partial derivatives
c-----
c
c       This program plots the observed and predicted
c       velocity or gamma dispersion corresponding to the
c       current model
c-----
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
        common/ctrl/ numa,d(NL),a(NL),b(NL),r(NL),rat(NL),
     1      dd(NL2),x(NL2),
     1      h(NL2),u(NL2),c(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        character *20 ibcx,ibcy
c-----
c       common for igetmod
c-----
        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        common/modtit/title
        character title*80
c-----
c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       get the command line 
c-----
        call gcmdln(nid)
c-----
        call gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup,dlam,qaqb,wref,invcsl,
     2      lstinv,twnmin,twnmax,iter,nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
c-----
c       iprog   = inversion type. Logical or of  rftn 2, 
c                 surf 1 - if does
c               not match, terminate this run
c       itot    = total number of iterations [37]
c       nf1     = 1 estimated stdev computed from residuals
c                 0 no scaling by residuals
c       nf2     = TOTAL number of Love Wave Gamma Modes
c                 0 DO NOT PROCESS Love Wave Gamma Data for Q
c       nf34    = TOTAL number of Love Wave modes to process 
c                 (from C and U)
c       nf5     = TOTAL number of Rayleigh Wave Gamma Modes
c                 0 DO NOT PROCESS Rayleigh Wave Gamma Data for Q
c       nf67    = TOTAL number of Rayleigh Wave modes to process 
c                 (from C and U)
c       nf10    = Input Format (from model file)
c                 0 - Inversion a, rho fixed
c                 1 - Inversion Poisson Ratio Fixed, 
c                 Rho computed from Vp
c       nfilt   = smoothing parameter 
c                 0  No model Weight  No smoothing
c                 1  Model Weight     No smoothing
c                 2  No model weight  Smoothing
c                 3  Model weight     Smoothing
c       nup     = state
c               =1 partials computed 
c               =0 surfinv run =1 before
c               =2 on update of invcsl>=1 or 3
c       dlam    = damping [32]
c       qaqb    = 2.25 default  [34]
c       wref    = reference frequency 1.0 is default [33]
c       invcsl  = 0 acausal, 1 uncoupled causel, 2 coupled causal [35]
c       invdep  = 0 last inversion was for depth
c               = 1 last inversion was for velocity and Q inverse
c       lstinv  = 2,3,4,5 depending on the last inversion
c               invdep = 1 for 2,3,4 and 0 for 5
c       twnmin  These give the receiver function window
c       twnmax
c       iter    Current iteration
c       nurftn  Number of receiver functions to be read
c       invdep  0 invert for layer thickness
c           1 invert for velocity
c       pval    0 invert for RFTN, 1 invert for SURF
c-----
        if(mod(iprog,2).ne.1)then
            WRITE(LOT,*)'srfphv96 requires a surface-wave ',
     1       'dispersion inversion'
            STOP
        endif
c-----
c       open graphics
c-----
        if(nid.eq.0)then
            call pinitf('SRFPHV96.PLT')
        else
            call pinitf('SRFPHG96.PLT')
        endif
c-----
c       now plot the velocity/Q inverse model(s) for reference
c       we do this first since the computation of gamma requires
c       the velocity
c-----
        call pltmod(nid)
c-----
c       get earth and q model to set up causal partials
c-----
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
        m = mmax
c-----
c       save current values
c-----
        do 39 i=1,mmax
            if(qb(i).gt.1.0)then
                qbinv(i) = 1.0/qb(i)
            else
                qbinv(i) =     qb(i)
            endif
            if(qa(i).gt.1.0)then
                qainv(i) = 1.0/qa(i)
            else
                qainv(i) =     qa(i)
            endif
            b(i) = vb(i)
            a(i) = va(i)
            d(i) = dl(i)
            r(i) = rho(i)
   39   continue

c-----
        ibcx = 'PERIOD'
        ncx = 6

        if(nid.eq.0)then
            id = 27
            ibcy='VELOCITY (KM/S)'
            ncy = 16
        else if(nid.eq.1)then
            id = 22
            ibcy='GAMMA (1/KM)'
            ncy = 12
        endif

        if(nf34.eq.0 .or. nf67.eq.0)then
            if(nf34.gt.0)ilvry=1
            if(nf67.gt.0)ilvry=2
            x0 = 2.7
            y0 = 1.0
            xaxlen = 6.0
            yaxlen = 6.0
        call vmat(id,wref,invcsl,nf1,x0,y0,xaxlen,yaxlen
     1      ,ilvry,ibcx,ibcy,ncx,ncy,mmax)
        else
c-----
c           plot Love
c-----
            x0 = 2.7
            y0 = 2.0
            xaxlen = 3.0
            yaxlen = 4.0
        call vmat(id,wref,invcsl,nf1,x0,y0,xaxlen,yaxlen
     1      ,1,ibcx,ibcy,ncx,ncy,mmax)
c-----
c           plot Rayleigh
            x0 = 6.5
            y0 = 2.0
            xaxlen = 3.0
            yaxlen = 4.0
        call vmat(id,wref,invcsl,nf1,x0,y0,xaxlen,yaxlen
     1      ,2,ibcx,ibcy,ncx,ncy,mmax)
        endif
        call pend()
        close(1,status='keep')
        end

        subroutine vmat(id,wref,invcsl,nf1,x0,y0,xaxlen,yaxlen
     1      ,ilvry,ibcx,ibcy,ncx,ncy,m)
c-----
c     a general purpose reformatting file
c-----
        parameter(LER=0,LIN=5,LOT=6)
        parameter(NL=200,NL2=NL+NL)
c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - numebr of layers in model
c       NL2 - number of columns in model (first NL/2 are
c           - velocity parameters, second NL/2 are Q values)
c-----
        common/ctrl/ numa,d(NL),a(NL),b(NL),r(NL),rat(NL),
     1      dd(NL2),x(NL2),
     1      h(NL2),u(NL2),ct(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        character ibcx*(*), ibcy*(*)

c-----
c       get inversion control
c-----
        open(1,file='tmpsrfi.04',form='unformatted',access='sequential')
        rewind 1
c-----
c       open dispersion file
c-----
        open(2,file='tmpsrfi.08',form='unformatted',access='sequential')
        rewind 2
c-----
c     read in output of disper, derivl, derivr (unit 1) tmpsrfi.04
c     read in observed data (unit 2) tmpsrfi.08
c
c     Note there is no one-to-one correspondence in the between
c     observed and theoretical. There can be more theoretical than
c     observed because of missing periods at various modes.
c     In addition, observed data can have multiple observations
c     at each period for a given mode
c
c     The mitigating factor is that both are arranged in order
c     of increasing period, as
c     LOVE
c        PHASE
c           MODES
c        GROUP
c           MODES
c     RAYLEIGH
c        PHASE
c           MODES
c        GROUP
c           MODES
c
c-----
        call plot( x0,  y0, -3)
c-----
c       perform this three times, once to get bounds
c       the second to plot theoretical, and the third
c       to plot observed
c-----
        vmin = 1e+20
        vmax = 0.0
        tmin = 1.0e+20
        tmax = 0.0
        nval = 0
        do 100 jrun = 1,3
c-----
c           jrun = 1  - get bounds
c           jrun = 2  - plot observed
c           jrun = 3  - plot predicted
c-----
            if(jrun.eq.3)then
                call newpen(2)
            else 
                call newpen(1)
            endif
            kold = 0
            md1old = 0
c-----
c           read observed dispersion
c-----
            rewind 1
            read(1) nd,m
            m2 = m + m
            read(1)(dd(i),i=1,m2)
            read(1)(wc(i),i=1,m2)
            rewind 2
c-----
c
c           initialize the search
c-----
            tp1 = 0.0
            md1 = 0
            k1 = 0
   10       continue
c-----
c           get observed data
c-----
            read(1,end=40)ifn,k,md,tp,c,sd
            if(c.eq.0.0)go to 10
            if(ifn.ne.ilvry)go to 10
            iwant = ifn*1000 +mod(k,2)*100 + md
            call gttheo(ifn,k,md,tp,itst,k1,md1,tp1,c1,iret,
     1          m,wref,invcsl,invdep,iwant,ihave)
            if(ihave.gt.iwant)then
                rewind 2
                go to 10
            endif
            if(iret.lt.0)then
                rewind 2
                tp1 = 0.0
                md1 = 0
                k1 = 0
                ihave = -1
                go to 10
            endif

c-----
c       if model does not support a higher mode drop from listing
c-----
            if(itst.eq.0)goto 10
            if(ifn .ne.ilvry)goto 10
            if(id.eq.27 .and. k.gt.2 .or. id.eq.22 .and. k.le.2)goto 10
            if(jrun .eq.1)then
c-----
c               get bounds for plot axes
c-----
                if(c.gt.vmax)vmax = c
                if(c1.gt.vmax)vmax = c1
                if(c.lt.vmin)vmin=c
                if(c1.lt.vmin)vmin = c1
                if(tp1.gt.tmax)tmax=tp1
                if(tp1.lt.tmin)tmin=tp1
                nval = nval + 1
            else if(jrun.eq.2)then
c-----
c               plot the observed data
c-----  
C               xx = xaxlen *(tp1 - tmin)/(tmax-tmin)
                xx = xaxlen * alog10(tp1/tmin)/alog10(tmax/tmin)
                yy = yaxlen *(c - vmin)/(vmax-vmin)
                if(k.eq.1)then
                    isym = 1
                else if(k.eq.2)then
                    isym = 2
                else if(k.eq.3) then
                    isym = 1
                endif
                call newpen(md)
                call symbol(xx,yy,0.05,char(isym),0.0,-1)
                call newpen(1)
                if(nf1.eq.0)then
                    dy = sd*yaxlen/(vmax-vmin)
                    dy = dy/2.
                    call plot(xx,yy+dy,3)
                    call plot(xx,yy-dy,2)
                    call plot(xx,yy,3)
                endif
            else if(jrun.eq.3)then
c-----
c               plot the predicted data
c-----
C               xx = xaxlen *(tp1 - tmin)/(tmax-tmin)
                xx = xaxlen * alog10(tp1/tmin)/alog10(tmax/tmin)
                yy = yaxlen *(c1 - vmin)/(vmax-vmin)
                if(k.eq.kold .and. md1.eq.md1old)then
                    call plot(xx,yy,2)
                else
                    call plot(xx,yy,3)
                endif
                kold = k
                md1old = md1
            endif
            goto 10
        write(LOT,*)'in main program - job aborted'
        stop
   40   continue
            if(nval.eq.0)then
                write(LER,*)'NO SUCH VALUES TO PLOT'
                go to 100
            endif
            if(jrun.eq.1)then
                call vmxmn(tmin,5,-1)
                call vmxmn(vmin,5,-1)
                call vmxmn(tmax,5,+1)
                call vmxmn(vmax,5,+1)
                if(tmax/tmin .lt. 10.0)tmax = 10.0*tmin
C               if(tmin.eq.tmax)tmax=tmin + 1.
                if(vmin.eq.vmax)vmax=vmin + 1.
            fstx = tmin
            dvx = (tmax-tmin)/xaxlen
            fsty = vmin
            dvy = (vmax-vmin)/yaxlen
        call dologx(0.0,0.0,xaxlen,tmax,tmin,
     1      0.10,.false.,.false.,.true.,ncx,ibcx)
C           call axis(0.0,0.0,ibcx,-ncx,xaxlen,0.0,fstx,dvx)
            call axis(0.0,0.0,ibcy,+ncy,yaxlen,90.,fsty,dvy)
        if(ilvry.eq.1)then
           call symbol(0.5,yaxlen+0.1,0.10,'LOVE',0.0,4)
        else if(ilvry.eq.2)then
           call symbol(0.5,yaxlen+0.1,0.10,'RAYLEIGH',0.0,8)
        endif
        endif
        call newpen(1)
  100   continue
        close(1,status='keep')
        close(2,status='keep')
        call plot(-x0, -y0, -3)
        return
        end

        subroutine vmxmn(val,inc,jmnmx)
c-----
c       change the number val to have two significant
c       figures evently divisible by inc.
c       jmnmx = negative go for minimum
c       jmnmx = positive go for maximum
c
c       Example for inc = 5, 0.13 -> 0.10 for min or 0.15 for a max
c       for inc = 10  0.13 -> 0.10 for a min or 0.20 for a max
c-----
        if(val.eq.0.0)return
        if(val.gt.0.0)then
            iplmn = +1
            imnmx = jmnmx
        else
            iplmn = -1
            imnmx = - jmnmx
            val = -val
        endif
        if(val.ge.1.0)then
            ival = alog10(val)
        else
            ival = alog10(val) - 1
        endif
c-----
c       convert val to an integer with two significant figures
c-----
        xval = val * 10.0**(1-ival)
        jval = xval
        jval = jval / inc
        if(imnmx.lt.0)then
            jval = jval * inc
        else
            jval = (jval+1)*inc
        endif
        val = iplmn * jval * 10.0**(ival-1)
        return
        end

        subroutine gttheo(ifn,k,md,tp,itst,k1,md1,tp1,c1,iret,
     1  m,wref,invcsl,invdep,iwant,ihave)
        parameter (NL=200,NL2=NL+NL)
        common/ctrl/numa,d(NL),a(NL),b(NL),r(NL),rat(NL),dd(NL2),x(NL2),
     $      h(NL2),u(NL2),ct(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        real*4 dcdb(NL2),dcda(NL2),dgdq(NL2),dudb(NL2),duda(NL2),
     1      dudq(NL2),dcdq(NL2),dgdv(NL2),dcdh(NL2), dudh(NL2)
        save dcdb,dcda,dgdq,dudb,duda,dudq,
     1      dcdq,dgdv,dcdh,dudh
        save cvel,gvel,gam
c-----
c       search through theoretical file to find a match
c       to a particular observation
c
c       OBSERVED
c           ifn : 1 Love 2 Rayl
c           k   : 1 Phase 2 Group
c           md  : Mode 
c           tp  : Period
c       THEORETICAL
c           itst    : Whether dispersion point found
c               : 0 mode does not exist, 1 = Love, 2 = Rayleigh
c           k1  : 1 Phase 2 Group
c           md1 : Mode
c           tp1 : Period
c           c1  : Velocity
c           dd  : Partial Derivatives
c-----
c
c       check for possibility of a double observation
c
c-----
c-----
c       initialize the dd array
c-----
        call zero(dd,1,m+m)
c-----
c       test to see if request is for information already in hand
c       NOTE: the test mod(k,2) == mod(k1,2) arises since
c       we only need dcdb, phase velocity and perehaps dcda to get
c       the anelastic attenutation coefficient gamma
c-----
        iret = 0
        dif = abs( (tp-tp1)/tp )
        jret = -1
        if(ifn.eq.itst .and. md.eq.md1 .and. dif.lt.1.0e-5)then
            if(mod(k,2).eq.mod(k1,2))then
c-----
c       should be .eq. 1 but to permit velocity or gamma
c       dispersion immediately after dispersion computation
c       use .ne. 0
c-----
                if(invdep.ne.0)then
                    if(k.eq.1)then
                        call cpy(dcdb,dd,0,m)
                        call cpy(dcdq,dd,m,m)
                        c1 = cvel
                    else if(k.eq.2)then
                        call cpy(dudb,dd,0,m)
                        call cpy(dudq,dd,m,m)
                        c1 = gvel
                    else if(k.eq.3)then
                        call cpy(dgdv,dd,0,m)
                        call cpy(dgdq,dd,m,m)
                        c1 = gam
                    endif
                else if(invdep.eq.0)then
                    if(k.eq.1)then
                        call zero(dd,1,m+m)
                        call cpy(dcdh,dd,0,m)
                        c1 = cvel
                    else if(k.eq.2)then
                        call zero(dd,1,m+m)
                        call cpy(dudh,dd,0,m)
                        c1 = gvel
                    else if(k.eq.3)then
                        call zero(dd,1,m+m)
                        c1 = gam
                    endif
                endif
                jret=1
            endif
        endif
c-----
c       test to see if information is in hand for a mode
c       not computed
c       E.g., 1 st mode not found, but we are asked for
c       phase and gamma, so read in another observed data
c       value rather than reading another theoretical value
c       and getting the two lists out of synchronization
c-----
        if(itst.eq.0 .and. md.eq.md1 .and. dif.lt.1.0e-5
     1      .and. mod(k,2).eq.mod(k1,2) )then
            jret = 1
        endif
        if(jret.eq.1)return
c-----
c       since we will not reuse previous value
c       read in new theoretical observation
c-----
   10   continue
        iret = 0
        read(2,end=30) itst,k1,md1,tp1
        ihave = 1000*itst + 100*mod(k1,2) + md1
        dif = abs ( (tp-tp1)/tp)
        if(k1.eq.2)then
            read(2,end=30)gvel,(dudb(j),j=1,m)
            read(2,end=30) (dudh(j),j=1,m)
            read(2,end=30) cvel,(dcdb(j),j=1,m)
            read(2,end=30) (dcdh(j),j=1,m)
            if(itst.eq.2)then
                 read(2,end=30)(dcda(j),j=1,m)
                 read(2,end=30)(duda(j),j=1,m)
            endif
        else
            read(2,end=30)cvel,(dcdb(j),j=1,m)
            read(2,end=30) (dcdh(j),j=1,m)
            if(itst.eq.2)read(2,end=30)(dcda(j),j=1,m)
        endif
        if(mod(k,2).ne.mod(k1,2).or.md.ne.md1.or.dif.gt.1.0e-5)then
            go to 10
        endif
c-----
c       define correct parameter
c-----
        call getgam(dcdb,b,qbinv,dcda,a,qainv,m,itst,cvel,tp,
     1      gam,dgdq)
        call zero(dcdq,1,m)
        call zero(dudq,1,m)
        call zero(dgdv,1,m)
        if(invcsl.ge.1)call cslmod(gvel,cvel,gam,wref,tp,m,
     1      qbinv,qainv,dcdb,dudb,itst,dcdq,dudq,dgdv,dcda
     2      ,a,b,k1)
c-----
c       setup causal, uncoupled if invcsl = 1
c-----
        if(invcsl.eq.1 .or. invdep.eq.0)then
            call zero(dudq,1,m)
            call zero(dcdq,1,m)
            call zero(dgdv,1,m)
        endif
        if(invdep.ne.0)then
            if(k1.eq.2)then
                call cpy(dudb,dd,0,m)
                call cpy(dudq,dd,m,m)
                c1 = gvel
            else if(k1.eq.1)then
                if(k.eq.1)then
                    c1 = cvel
                    call cpy(dcdb,dd,0,m)
                    call cpy(dcdq,dd,m,m)
                else if(k.eq.3)then
                    c1 = gam
                    call cpy(dgdv,dd,0,m)
                    call cpy(dgdq,dd,m,m)
                endif
            endif
        else if(invdep.eq.0)then
            if(k1.eq.2)then
                call zero(dd,1,m+m)
                call cpy(dudh,dd,0,m)
                c1 = gvel
            else if(k1.eq.1)then
                if(k.eq.1)then
                    call zero(dd,1,m+m)
                    c1 = cvel
                    call cpy(dcdh,dd,0,m)
                else if(k.eq.3)then
                    call zero(dd,1,m+m)
                    c1 = gam
                endif
            endif
        endif
        if(ihave.lt.iwant)go to 10
        return
   30   continue
            iret = -1
        return
        end

        subroutine cslmod(u,c,gam,wref,tp,m,
     1      qbinv,qainv,dcdb,dudb,itst,dcdq,dudq,dgdv,dcda
     2      ,a,b,k1)
        real*4 dcdb(*),dudb(*),qbinv(*),dcdq(*),dudq(*)
     1      ,dgdv(*),dcda(*),a(*),b(*),qainv(*)
c-----
c       modify phase, group velocities and partials for
c       constant Q causality
c
c       itst 1 = L, 2 = R
c-----
        c0 = c
        u0 = u
        f  = 1./tp
        pi = 3.1415927
        omega = 6.2831853*f
        faclog = alog(f/wref)/pi
        facg = 2.*gam*u0/(pi*omega)
c-----
c       get correct phase velocity
c-----
        c = c0 + 2.*gam*c0*c0*faclog/omega
c-----
c       get correct group velocity
c-----
        cmc0c0 = (c-c0)/c0
        if(k1.eq.2)then
            u0c0 = u0/c0
            u = u0 *( 1. + (2. - u0c0)*cmc0c0 + facg ) 
            uu0 = u/u0
        endif
c-----
c       get correct partials
c-----
        do 100 i=1,m
            if(b(i).eq.0.0)then
                qaiqbi = 1.0e-6
            else 
                if(qbinv(i).gt.0.0)then
                    qaiqbi = qainv(i)/qbinv(i)
                else
                    qaiqbi = 1.0/(0.75*(a(i)/b(i))**2)
                endif
            endif
            dcdp = 0
            if(itst.eq.2)dcdp=dcda(i)
            dcds = dcdb(i)
            dcdb(i) = dcds*(1. + faclog*qbinv(i) )
            if(k1.eq.2)then
            dudv = dudb(i)
            dudb(i) = dudv*(uu0 -u0c0*cmc0c0 + facg )
     1       + dcds*u0c0*(-2.*facg +u0c0*qbinv(i)/pi +
     2          (2.-u0c0)*(faclog*qbinv(i) -cmc0c0) +
     3          u0c0 *cmc0c0 )
            dudq(i) = u0c0*(2.-u0c0)*dcdq(i) +
     1          u0c0*u0c0*(dcds*b(i)+dcdp*a(i)*qaiqbi)/pi
     2          - 2*gam*dcds/c0
            endif
            dcdq(i) = faclog*(dcds*b(i) +dcdp*a(i)*qaiqbi)
            dgdv(i) = 0.5*omega*dcds*qbinv(i)/(c0*c0)
  100   continue
        return
        end

        subroutine getgam(dcdb,b,qbinv,dcda,a,qainv,m,ilvry,cc,tp,
     1          gam,dgdq)
        real*4 dcdb(*),b(*),dcda(*),a(*),dgdq(*),qbinv(*), qainv(*)
c-----
c       determine spatial attenuation at constant frequency
c       and partial derivatives
c-----
        omega = 6.2831853/tp
        factr = 0.5*omega/(cc*cc)
        sum = 0.0
        do 100 i=1,m
            sum = sum + dcdb(i)*b(i)*qbinv(i)
            dgdq(i) = factr * dcdb(i)*b(i)
            if(ilvry.eq.2)then
                if(b(i).eq.0.0)then
                    qaiqbi = 1.0e-6
                else 
                    if(qbinv(i).gt.0.0)then
                        qaiqbi = qainv(i)/qbinv(i)
                    else
                        qaiqbi = 0.75*(a(i)/b(i))**2
                    endif
                endif
                sum = sum + dcda(i)*a(i)*qainv(i)
                dgdq(i) = dgdq(i) + factr*dcda(i)*a(i)*qaiqbi
            endif
  100   continue
        gam = factr * sum
        return
        end

        subroutine zero(x,m,n)
        real*4 x(*)
c-----
c       set x(m) = 0, x(m+1)=0, ..., x(n)=0
c-----
        do 100 i=m,n
            x(i) = 0.0
  100   continue
        return
        end

        subroutine cpy(x,y,m,n)
        real*4 x(*),y(*)
c-----
c       copy contents of array x into array y with offset m
c-----
        do 100 i=1,n
            y(i+m) = x(i)
  100   continue
        return
        end

        subroutine pltmod(nid)
c-----
c       plot the velocity or the Q inverse model
c       
c       nid I*4 0 for velocity 1 for gamma/Q
c-----
c-----
c       common for igetmod
c-----
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)
        common/isomod/d(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        common/modtit/title
        character title*80
c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c-----
        real barr(4), qbarr(4), darr(4)
c-----
c       set the plot
c-----  
        x0 = 0.5
        y0 = 2.0
        call plot(x0,y0,-3)
        xlen = 1.5
        ylen = 4.0
c-----
c       get the velocity model twice, once for the
c       current model and then for the initial model
c       however we use both to get the plot extremes
c-----
c-----
c-----
c       get bounds for the plot
c-----
        bmin = 1.0e+38
        bmax = 0.0
        qbmin = 1.0e+38
        qbmax = 0.0
        zmin = 0.0
        zmax = 0.0
        do 2000 jmd=1,2
c-----
c           get velocity model and q(beta) inverse model
c-----
            if(jmd.eq.1)then
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            else
        call getmod(2,'tmpmod96.000',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            endif
        mlyr = mmax
            depth = 0.0
        do 1000 i=1,mlyr
            if(nid.eq.0)then
                if(vb(i).gt.bmax)bmax = vb(i)
                if(vb(i).lt.bmin)bmin = vb(i)
            else if(nid.eq.1)then
                if(qb(i) .gt.1.0)qb(i) = 1.0/qb(i)
                if(qb(i).gt.qbmax)qbmax = qb(i)
                if(qb(i).lt.qbmin)qbmin = qb(i)
            endif
            if(i.lt.mlyr)depth = depth + d(i)
 1000   continue
        if(depth+d(mlyr-1).gt.zmax)zmax = depth+d(mlyr-1)
 2000   continue
c-----
c       safety test - force a minimum separation in km/sec units
c-----
            bmax = bmax + 0.125
            bmin = bmin - 0.125
        
        darr(1) = 0.0
        darr(2) = zmax
c-----
c       now scale the velocity/Q inverse axis
c-----
        barr(1) = bmin
        barr(2) = bmax
        call dvscale(bmin, bmax)
        barr(3) = bmin
        barr(4) = (bmax - bmin)/xlen
        qbarr(1) = 0.95*qbmin
        qbarr(2) = 1.05*qbmax
        call gscale(qbarr,xlen,2,1)
        qbmin = qbarr(3)
        qbmax = qbarr(3) + xlen * qbarr(4)

        yylen = 10.0
        call gscale(darr,yylen,2,1)
        dmin = 0.0
C       dmax = darr(3) + yylen*darr(4)
        dmax = darr(2)
c-----
c       create bounding box
c-----
        call gbox(0.0,0.0,xlen,ylen)
c-----
c       label the depth axis
c-----
            call labz(dmin,dmax,10,0.0,ylen,ylen,0,dd)
c-----
c       plot the velocity/q inverse axis
c-----
        if(nid.eq.0)then
            call dolinx(0.0,ylen,xlen,bmax,bmin,
     1          0.05,.true.,.true.,.true.,9,'VS (KM/S)')
C           call myaxis(0.0,ylen,'VS (KM/S)', 9,xlen,0.0, 
C     1             barr(3),barr(4))
        else
            call myaxis(0.0,ylen,'QBINV'    , 5,xlen,0.0,
     1          qbarr(3),qbarr(4))
        endif
c-----
c       plot the velocity and Qb models
c-----
c-----
c       now read the model file two times, once for original
c       and then for current
c-----
        do 3000 jmd=1,2
            if(jmd.eq.1)then
        call getmod(2,'tmpmod96.000',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            call newpen(4)
            else
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            call newpen(2)
            endif
        mlyr = mmax
        if(nid .eq. 0)then
            dx = 1.0/barr(4)
            call vplt(0.0,ylen,mlyr,vb,d,bmin,bmax,dx,
     1          dmin,dmax,dd,jmd)
        else if(nid.eq.1)then
            do 2910 i=1,mlyr
                if(qb(i).gt.1.0)then
                    qb(i) = 1.0/qb(i)
                endif
 2910       continue
            dx = 1.0/qbarr(4)
            call vplt(0.0,ylen,mlyr,qb,d,qbmin,qbmax,dx,
     1          dmin,dmax,dd,jmd)
        endif
 3000   continue
c-----
c       annotate
c-----
        call newpen(2)
        call plot(0.0,-0.25,3)
        call plot(0.25*xlen, -0.25, 2)
        call symbol(0.30*xlen, -0.25, 0.10, ' Current',0.0,8)
        call newpen(4)
        call plot(0.0,-0.50,3)
        call plotd(0.25*xlen,-0.50,21,0.05)
        call symbol(0.30*xlen, -0.50, 0.10, ' Initial',0.0,8)
c-----
c       reset the plot
c-----
        call plot(-x0,-y0,-3)
        return
        end

c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: AXIS                                                  c
c                                                                     c
c      COPYRIGHT (C)  1986 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
        subroutine myaxis(xpage,ypage,ititl,nchar,axlen,angle,
     1          firstv,deltav)
        real HTNUM,SYMHT,TICHT
       parameter (HTNUM=0.07,SYMHT=0.10,TICHT=0.07)
c-----
c       xpage,ypage  coordinates of starting point of axis in
c                    inches relative to current origin
c       ititl   axis title
c       nchar   number of characters in title
c               >0 for tic marks, numbers and title on
c                  counterclockwise side
c               <0 for tic marks, numbers and title on
c                  clockwise side
c       axlen   floating point axis length in inches
c       angle   angle (degrees) that x axis makes with
c               horizontal x-direction
c       firstv  scale value at (xpage,ypage)
c       deltav  change in slace between tic marks per unit
c               length
c               <0 value at xpage,ypage is a max, and values
c               decrease along the axis
c-----
c-----
c       we do not worry about the vagaries of storage of characters
c       since the address is passed on down to subroutine symbol
c-----
        character*(*) ititl
        a=1.0
        kn=nchar
        if(kn.lt.0)then
                a= -a
                kn= -kn
        endif
c-----
c       if deltav is too large invoke scientific notation
c-----
        ex = 0.0
        if(deltav.ne.0.0)then
                yex= alog10(abs(deltav))
                if(yex.lt.-2.0)then
                ex = aint(yex+0.01) - 1.0
            else if(yex.lt.-1.0)then
                ex = aint(yex+0.01) 
            endif
                if(yex.ge.2.0) ex = aint(yex + 0.01)
        endif
        xval = firstv*10.0**(-ex)
        xdel = deltav*10.0**(-ex)
        ct = cos(angle*0.01745329)
        st = sin(angle*0.01745329)
        ntic = axlen + 1.0
c-----
c       first put in numbers and title
c       adjust offset for numbers
c-----
        dx = - HTNUM
        dy = 1.5*a*HTNUM - 0.5*HTNUM
c-----
c       find initial position given rotation
c-----
        xn = xpage + dx*ct - dy*st
        yn = ypage + dy*ct + dx*st
        mtic = ntic/2
        do 100 i=1,ntic
            if(i.eq.1)then
                if(xval.lt.0.0)then
                    dll = 2.0*HTNUM
                else
                    dll = 1.0*HTNUM
                endif
                if(abs(xval).ge.10.0)dll = dll + HTNUM
                dxx = dll*ct
                dyy = dll*st
                    call number(xn+dxx,yn+dyy,HTNUM,xval,angle,2)
            else if(i.eq.ntic)then
                if(xval.lt.0.0)then
                    dll = 4.0*HTNUM
                else
                    dll = 3.0*HTNUM
                endif
                if(abs(xval).ge.10.0)dll = dll + HTNUM
                dxx = -dll*ct
                dyy = -dll*st
                    call number(xn+dxx,yn+dyy,HTNUM,xval,angle,2)
            else
                    call number(xn,yn,HTNUM,xval,angle,2)
            endif
                xval=xval+xdel
                xn=xn+ct
                yn=yn+st
c-----
c       halfway down axis put in title
c-----
                if(i.eq.mtic)then
                        z=kn
                        if(ex.ne.0.0)z=z+7.0
                        dx = -0.5*SYMHT*z + axlen*0.5
                        dy = (2.5*a-0.5)*SYMHT
                        xt=xpage +dx*ct-dy*st
                        yt=ypage +dy*ct+dx*st
                        call symbol(xt,yt,SYMHT,ititl,angle,kn)
                        if(ex.ne.0.0)then
                                z=kn+2
                                xt=xt+z*ct*SYMHT
                                yt=yt+z*st*SYMHT
                                call symbol(xt,yt,SYMHT,'*10',angle,3)
                                xt=xt+(3.0*ct-0.8*st)*SYMHT
                                yt=yt+(3.0*st+0.8*ct)*SYMHT
                                call number(xt,yt,0.7*SYMHT,ex,angle,
     1                                  -1)
                        endif
                endif
  100   continue
c-----
c       now put in tic marks
c-----
        call plot(xpage+axlen*ct,ypage+axlen*st,3)
        dx = - TICHT*st*a
        dy =   TICHT*ct*a
        a = ntic -1
        xn = xpage + a*ct
        yn = ypage + a*st
        do 200 i=1,ntic
                call plot(xn,yn,2)
                call plot(xn+dx,yn+dy,2)
                call plot(xn,yn,2)
                xn=xn-ct
                yn=yn-st
  200   continue
        return
        end

        subroutine labz(dmin,dmax,nd,x0,y0,ylen,id,ddd)
c-----
c       label the depth axis - however 
c-----
c       dmin    R*4 - minimum depth
c       dmax    R*4 - maximum depth
c       nd  I*4 - number of segments
c       x0  R*4
c       Y0  R*4 - upper left corner
c       ylen    R*4 - length of Y-axis
c       id  I*4 - 0 single
c                 1 double
c-----
        character ostr*10, fmt*10
        dd = (dmax - dmin) / real (nd)
        if(id.eq.1)then 
            ndo = 2
        else
            ndo = 1
        endif
        if(ndo.eq.1)then
            dz = ylen / real(nd)
        else
            dz = 0.5 * ylen / real(nd)
        endif
        ddd = dz / dd
        if(dmax.gt.10 .and. dmax.lt.100)then
            ndec = 2
            fmt = '(f10.2)'
        else if(dmax.ge.100)then
            ndec = 0
            fmt = '(f10.0)'
        else
            ndec = 3
            fmt = '(f10.3)'
        endif
        do 100 n=1,ndo
            do 110 m=1,nd-1
                ypos = ylen - real(m)*dz - (n-1)*0.5*ylen
                yval = dmin + m*dd
                call plot(x0,ypos,3)
                call plot(x0+0.10,ypos,2)
                write(ostr,fmt)yval
                ht = 0.05   
                call symbol(x0-11.0*ht,ypos-0.5*ht,ht,ostr,0.0,10)
  110       continue
  100   continue
        hht = 0.07
            nchar = 10
            call symbol(x0-7.0*ht,ylen/2.0 - nchar*0.5*hht,
     1          hht,'DEPTH (KM)',90.0,nchar)
        return
        end
                
        subroutine vplt(xpos,ypos,n, x, y,xmin,xmax,dx,
     1          ymin,ymax,dy,iskp)
c-----
c       y   R*4 - array of parameters to be plotted along y-axis
c       x   R*4 - array of parameters to be plotted along x-axis
c       iskp    I*4 - 2 solid line
c                 1 dashed line
c-----
c       plotting velocity model is odd since have to force the steps
c-----
        parameter(NL=200,NL2=NL+NL)
        real*4 x(*), y(*)
        dphl = 0.0
        dphu = y(1)
        do 100 i=2,n
            xx1 = xpos + (x(i-1) - xmin)*dx
            xx2 = xpos + (x(i  ) - xmin)*dx
            yy1 = ypos - (dphl - ymin)*dy
            yy2 = ypos - (dphu - ymin)*dy
            call plot(xx1,yy1,3)
            if(iskp.eq.2)then
                call plot(xx1,yy2,2)
                call plot(xx2,yy2,2)
            else
                if(xx1.ne.xx2 .or. yy1.ne.yy2)then
                call plotd(xx1,yy2,21,0.05)
                call plotd(xx2,yy2,21,0.05)
                endif
            endif
            dphl = dphu
            dphu = dphl + y(i)
  100   continue
        dphu = 1.10*dphl
        yy2 = ypos - (dphu - ymin)*dy
        if(iskp.eq.2)then
            call plot(xx2,yy2,2)
        else
            call plotd(xx2,yy2,21,0.05)
        endif
        return
        end

        subroutine gcmdln(nid)
c-----
c       parse the command line parameters
c-----
c       nid I*4 0 for velocity dispersion and model (default)
c               1 for gamma dispersion and Qb inv model
c-----
        integer*4 nid
        character names*20
        nmarg = mnmarg()
c-----
c       initialize defaults
c-----
        nid = 0
c-----
c       loop through command line arguments
c-----

        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,names)
            if(names(1:2).eq.'-V' .or. names(1:2).eq.'-v')then
                nid = 0
            else if (names(1:2).eq.'-G' .or. names(1:2).eq.'-g')then
                nid = 1
            endif
        go to 1000
 2000   continue
        return
        end

        subroutine gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup,dlam,qaqb,wref,invcsl,
     2      lstinv,twnmin,twnmax,iter,nurftn,invdep,pval,
     3      sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
c-----
c       read control file
c-----
        integer NL
        parameter (NL=200)
        integer nf10(NL)
        open(1,file='tmpsrfi.00',form='unformatted',access='sequential')
        rewind 1
        read(1) iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,nup,dlam
     1      ,qaqb,wref,invcsl,lstinv,twnmin,twnmax,iter,
     2      nurftn,invdep,pval,
     3      sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5
        close(1,status='keep')
        return
        end

        subroutine dvscale(vmin,vmax)
c-----
c       rescale xval so that a plot will look nice
c       this is currently for velocity which can never be negative
c-----
c       vmin    R   - value to be scaled
c       vmax    R   - value to be scaled
c-----
c-----
c       first get size
c-----
        size = max(vmin,vmax)
        if(size.gt.1)then
            lpow = alog10(size)
            tmp = 10.0**lpow
            if(size.lt.tmp)then
                tmp = 0.1 * tmp
            endif
        else
            lpow = alog10(size) + 1
            tmp = 10.0**lpow
            if(size.lt.tmp)then
                tmp = 0.1 * tmp
            endif
        endif
        vmn = vmin
        do 1000 i=0,20
            xx = real(i) * 0.5 *tmp
            if(xx .lt. vmin)then
                vmn = xx 
            endif
 1000   continue
        vmx = vmax
        do 2000 i=20,0,-1
            xx = real(i) * 0.5 *tmp
            if(xx .gt. vmax)then
                vmx = xx 
            endif
 2000   continue
C       write(6,*)vmin,vmax,vmn,vmx,size,lpow,tmp
        if(vmx.gt.vmn)then
            vmin = vmn
            vmax = vmx
        endif
C       write(6,*)vmin,vmax,vmn,vmx,size,lpow,tmp

        return
        end
