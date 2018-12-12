        program tspec96p
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME VI                                                      c
c                                                                     c
c      PROGRAM: TSPEC96P                                              c
c                                                                     c
c      COPYRIGHT 2002                                                 c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63013                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       CHANGES
c       11 SEP 2000 - build in P, SV and SH first arrival times
c       14 JAN 2001 - put in A, C, F, L and N constants into
c           trace header
c       08 JUN 2001 - error in mapping of ksrc - fixed to
c           ksrc(i) = jsrc(lsrc(i))
c       18 JUL 2001 - fixed swpecification for KSRC when source and
c           receiver are in fluids. 
c       09 APR 2002 - another hack at fixing the 
c           JSRC when fluids are involved.
c           commented out lines in gethsp which was not correct, 
c           modified output loop
c       25 APR 2002 - implement TI media with vertical symmetry
c       05 MAY 2002 - corrected putting A C F L and N into header
c       29 OCT 2002 - Fixed evalg and evalh for fluid layers
c       20 JUN 2003 - changed getegn and section evlmat before call hska
c           STILL MUST CORRECT DC CORRECTION SINCE FOR TI MEDIA
c           RP !-> k, RSV !-> k for large k but rather
c       02 FEB 2008 - corrected error in use of atna atnb in getegn/rshof
c       08 SEP 2008 - correctly defined mmax as an integer in subroutine
c               dezero
c-----
        parameter(LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        common/damp/alpha,ieqex
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/modelt/tD(NL),tTA(NL),tTC(NL),tTL(NL),tTN(NL),tTF(NL),
     1      tTRho(NL),
     2      tqa(NL),tqb(NL),tetap(NL),tetas(NL), 
     3      tfrefp(NL), tfrefs(NL)
        real tD, tTA, tTC, tTF, tTL, tTN, tTRho
        real tqa, tqb, tetap, tetas, tfrefp, tfrefs
        common/jout/jsrc(21) , jbdrys, jbdryh
c-----
c       matrix components in layers and boundaries saved
c-----
        complex*16 har(NL,4,4), dar(NL,5,5), hsr(2,5), gbr(2,5), 
     1      hal(NL,2,2), hsl(2,2), gbl(2,2)
        real*8 hex(NL), lex(NL), dex(NL), hexw(NL)
        common/hamat/har
        common/damat/dar
        common/hsrfr/hsr
        common/gbrfr/gbr
        common/hlmat/hal
        common/hsrfl/hsl
        common/gbrfl/gbl
        common/hexex/hex
        common/hexexw/hexw
        common/dexex/dex
        common/lexex/lex 
        common/water/iwater(NL),iwats(2),iwatb(2)
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
        common/lyrctl/lyrins
        logical lyrins
        common/rlimit/rlim
        real*4 rlim
c-----c
        real*4 ffreq(8)
        parameter(NDST=500)
        real*4 r(NDST), tshift(NDST), vred(NDST)
        complex smo(21)
        complex ztmp, zdata
        common/c/pmin,pmax,dp,pcntrl
        common/frlim/fl,fu,df,fwhich
        common/cntrl/ishank,hnkarg,dstcor,dokjar


        real SSA, SSC, SSF, SSL, SSN, SSR
        real RRA, RRC, RRF, RRL, RRN, RRR
c-----
c       ishank  L   - .true. use Hankel function and not Bessel
c                    not that asymptotic tricks are not used
c       hnkarg  R*4 - (default 6.0) For kr > hnkarg use the Hankel
c                   function, otherwise the Bessel
c       isfref  L   - .true. reference frequency uis not the 
c                   default of 1.0 Hz
c       fref    R*4 - reference frequency for Causal Q
c       dstcor  R*4 - 0 first time sample is tshift + r/vred
c                 1 first time sample is tshift + z/vred where
c                   z = abs (source depth - receiver depth)
c                 2 first time sample is tshift + 
c                   sqrt(z*z + r*r) /vred
c       dokjar  L   - .false. use Futterman causal Q
c                 .true.  use Kjartansson Causal Q
c-----
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        logical ext
        character mname*80, title*80

        common/depref/refdep
        common/srcinf/vsa, vsb, vsr

c-----
c       ksrc = temporary array for jsrc for output
c       here jsrc != 0 reflects source information
c       if receiver is not in a fluid DO NOT output pressure field
c-----
        integer ksrc(21)
c-----
c       lsrc maps jsrc to output Green's functions. e.g., if
c       jsrc(8) = radial explosion term, but in final output it
c       occupies position 10, or jsrc(lsrc(10)) = computed
c-----
        integer*4 lsrc(21)
        data lsrc/1,2,3,4,13,5,6,14,7,8,9,10,11,12,15,16,17,18,19,20,21/
        fwhich = -1.0
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln()
c-----
c       set up tolerances
c           rlim is the distance which is effectively zero
c-----
            rlim = 1.0e-07
c-----
c       See if the file hspec96p.dat exists, if it does
c       open it for all control information
c-----
        inquire(file='hspec96p.dat',exist=ext)
        if(.not. ext)then
                write(LER,*)'Control file ',
     1          'hspec96p.dat does not exist'
                go to 9999
        endif
        open(3,file='hspec96p.dat',access='sequential',
     1      form='formatted',status='unknown')
        rewind 3
        call gethsp(3,mname,fl,fu,delt,n,n1,n2,xleng,xfac,
     1      ndist,r,tshift,vred)
        close (3)
c-----
c       process
c-----
        df = 1./(n*delt)
        nyq = n/2 + 1
        nyq2 = 2*nyq
        write(LOT,2)  fl,fu,df,n1,n2,n,
     2      vmin,vamin,vamax,vbmin,vbmax
c-----
c       UNIX output - no carriage control
c-----
    2   format('fl =',f10.5,5x,'fu =',f10.5,5x,'df =',f10.5,/
     1      4x,'n1 =',i4,5x,'n2 =',i4,5x, ' n =',i5/
     2  'vmin  =',f10.5,' vamin =',f10.5,' vamax =',f10.5/
     3  '       ',10x  ,' vbmin =',f10.5,' vbmax =',f10.5)
 2021   format('depths =',f20.6)
 2031   format('depthr =',f14.6)
 2041   format('     r =',f14.6,' tshift=',f10.2,' vred=',f10.2)
    3   format(8f10.5)
    4   format('frequencies for which response computed     ')
    5   format('alpha =',f10.5,5x,'dt =',f10.3)
c-----
c       MSDOS MICROSOFT FORTRAN OUTPUT
c-----
c    2   format(1x,'fl =',f10.5,5x,'fu =',f10.5,5x,'df =',f10.5,/
c     1          1x,4x,'n1 =',i4,5x,'n2 =',i4,5x, ' n =',i5/
c     2  1x,'vmin  =',f10.5,' vamin =',f10.5,' vamax =',f10.5/
c     3  1x,'       ',10x  ,' vbmin =',f10.5,' vbmax =',f10.5)
c 2021   format(1x,'depths =',f20.6)
c 2031   format(1x,'depthr =',f14.6)
c 2041   format(1x,'     r =',f14.6,' tshift=',f10.2,' vred=',f10.2)
c    3   format(1x,8f10.5)
c    4   format(1x,'frequencies for which response computed     ')
c    5   format(1x,'alpha =',f10.5,5x,'dt =',f10.3)
c-----
        write(LOT,*)'SOURCE DEPTH (',mdpths,')'
        do 2020 i=1,mdpths
            write(LOT,2021)depths(i)
 2020   continue
        write(LOT,*)'RECEIVER DEPTH (',mdpthr,')'
        do 2030 i=1,mdpthr
            write(LOT,2031)depthr(i)
 2030   continue
        write(LOT,*)'RECEIVER DISTANCES (',ndist,')'
        do 2040 i=1,ndist
            write(LOT,2041)r(i), tshift(i), vred(i)
 2040   continue
        write(LOT,5)alpha,delt
        if(dokjar)then
        write(LOT,*)'Kjartansson Constant Q operator used'
        endif
        if(dosud)then
            write(LOT,*)'Wavefield separation at Source'
            write(LOT,*)'    Upward P-wave  ',spup
            write(LOT,*)'    Upward S-wave  ',ssup
            write(LOT,*)'    Downward P-wave',spdn
            write(LOT,*)'    Downward S-wave',ssdn
        endif
        if(dorud)then
            write(LOT,*)'Wavefield separation at Receiver'
            write(LOT,*)'    Upward P-wave  ',rpup
            write(LOT,*)'    Upward S-wave  ',rsup
            write(LOT,*)'    Downward P-wave',rpdn
            write(LOT,*)'    Downward S-wave',rsdn
        endif
        write(LOT,4)
c-----
c     open output file for hspec96
c-----
      open(unit=2,status='scratch',form=
     1            'unformatted',access='sequential')
      rewind 2
c-----
c       process the frequencies
c-----
        do 101 i=1,8
  101   ffreq(i)=-1.0
        ndist = 1
        call bufini(1,ierr)
        do 100 ii = n1,n2
            freq=(ii-1)*df
            if(freq.lt.df) freq = 0.01*df
            call excit(freq,xleng,xfac,dk,nk,omega)
            index=mod(ii,8)
            if(index.eq.0)index=8
            ffreq(index)=freq
            if (index.eq.8) then
                    write(LOT,3)ffreq
                    do 102 ij=1,8
                        ffreq(ij)=-1.
  102               continue
            endif
  100   continue
        call buflsh()
c-----
c       output the final spectrum as a function of distance
c-----
        open(unit=4,file='hspec96.grn',status='unknown',
     1      form='unformatted',access='sequential')
        rewind 4
c-----
c       iprog   1 hspec96
c           2 hspec96p
c           3 hwhole96
c-----
        iprog = 2
        write(4)iprog
        write(4) alpha,fl,fu,delt,n,n1,n2,df,nyq2
        write(4)mname
c-----
c       now output the spectrum for each distance
c
c       The order in the temporary file 'hspec91.tmp' is
c           FREQ
c               RAY PARAMETER
c                   SOURCE_DEPTH
c                       RECEIVER_DEPTH
c       This must be rearranged to form
c           RAY PARAMETER
c               SOURCE_DEPTH
c                   RECEIVER_DEPTH
c                       FREQ
c
c-----
c       because first arrival may be non-causal, shift time series
c       by 40*DT seconds
c-----
        tshft = -40.0*delt
c-----
c       TP TSV and TSH have no meaning in ordeinary sense
c-----
        TP = -12345.
        TSV = -12345.
        TSH = -12345.
c-----
c       output by ray parameter
c-----
        do 5000 jd=1,nk
        p  = (jd-1)*dp + pmin
        k = 0
        do 5005 js=1,mdpths
        do 5010 jr=1,mdpthr
            k = k + 1
            rewind 2
            call bufini(0,ierr)
            call srclyr(depths(js),lmax,dph)
c-----
c       define TI constants
c-----
            TP =-12345.
            TSV=-12345.
            TSH=-12345.
            SSR = TRho(lmax)
            SSL = TL(lmax)
            SSN = TN(lmax)
            SSA = TA(lmax)
            SSC = TC(lmax)
            SSF = TF(lmax)
            
            write(4)p,tshft,depths(js)-refdep,
     1          depthr(jr)-refdep,
     2          TP, TSV, TSH, 
     3          SSA, SSC, SSF, SSL, SSN, SSR
            call srclyr(depthr(jr),lmax,dph)
            RRN = TN(lmax)
            RRA = TA(lmax)
c-----
c       Dislocations and forces must act in a solid source region
c       if the receiver is in a fluid, 
c           then permit pressure field output
c-----
            do 5011 i=1,21
                ksrc(i) = jsrc(lsrc(i))
                if(ksrc(i).gt.0)then
                if(i.ge.1.and.i.le.8)then
                    if(SSN .le. 0.0001*SSA)then
                        ksrc(i) = 0
                    endif
                else if(i.ge.11.and.i.le.15)then
                    if(SSN .le. 0.0001*SSA)then
                        ksrc(i) = 0
                    endif
                elseif(i.eq.16)then
                    if(RRN .lt. 0.0001*RRA)then
                        ksrc(i) = 1
                    else
                        ksrc(i) = 0
                    endif
                else if(i.ge.17)then
                    if(RRN.lt.0.0001*RRA.and.SSN.gt.0.0)then
                        ksrc(i) = 1
                    else
                        ksrc(i) = 0
                    endif
                endif
                endif
 5011       continue
            write(4)ksrc
            do 5200 i=n1,n2
                freq=(i-1)*df
                fac = 6.2831853*freq*tshft
                ztmp = cmplx(cos(fac), sin(fac) )
                do 5300 jjd=1,nk
                kk = 0
                do 5301 jjs=1,mdpths
                do 5302 jjr=1,mdpthr
                    kk = kk + 1
                    do 5401 jj=1,21
                        if(jsrc(jj).eq.1)then
                            call bufrd(xr,ierr)
                            call bufrd(xi,ierr)
                            smo(jj)=cmplx(xr,xi)
                        endif
 5401               continue
                    do 5400 jj=1,21
                        if(jsrc(lsrc(jj)).eq.1)then
                        if(jjd.eq.jd .and. k.eq.kk.and.
     1                      ksrc(jj).ne.0)then
                            zdata=ztmp*smo(lsrc(jj))
                            datar= real(zdata)
                            datai=aimag(zdata)
                            write(4)datar,datai
                        endif
                        endif
 5400               continue
 5302           continue
 5301           continue
 5300           continue
 5200       continue
 5010   continue
 5005   continue
 5000   continue
        rr = -1.0
        tt0 = 0.0
        write(4)rr,tt0
        close (4)
        close(2)
 9999   continue
        end

        subroutine gethsp(lun,mname,fl,fu,delt,n,n1,n2,xleng,xfac,
     1      ndist,r,tshift,vred)
c-----
c       read in data file of hspec8 commands
c-----
        parameter(LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        common/damp/alpha,ieqex
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        parameter(NDST=500)
        real*4 r(NDST), tshift(NDST), vred(NDST)
        common/c/pmin,pmax,dp,pcntrl
        character ostr*80
        character mname*80, title*80
        common/lyrctl/lyrins
        logical lyrins
        logical ext

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar

        common/depref/refdep
c-----
c-----
c       UNIX FORTRAN - NO CARRIAGE CONTROL
c-----
   21   format(11i5)
   22   format(9e11.4)
   24   format('XLENG=',e15.7,' XFAC=',e15.7)
   30   format(2e15.7/1x,3i10)
c-----
c       MSDOS MICROSOFT FORTRAN OUTPUT - CARRIAGE CONTROL
c-----
c   21  format(1x,11i5)
c   22  format(1x,9e11.4)
c   24  format(1x,'XLENG=',e15.7,' XFAC=',e15.7)
c   30  format(1x,2e15.7/1x,3i10)
c-----
c       read in distance correction
c-----
        read(lun,*)idcor
        if(idcor.ge.0 .and. idcor.le.2)then
            dstcor = idcor
        else
            dstcor = 0
        endif
c-----
c       read in time domain damping, sampling interval
c-----
        read(lun,*)alphat,delt
c-----
c       read in number of time samples, frequency limits
c-----
        read(lun,*)n,n1,n2  
        alpha = alphat/(n*delt)
        write(LOT,30) alpha,delt,n,n1,n2
        df = 1.0/(n*delt)
        fl = (n1-1)*df
        fu = (n2-1)*df
c-----
c       Specify desired output Green's functions
c-----
        read(lun,*)ieqex
        if(ieqex.gt.6 .or. ieqex.lt.0)ieqex = 2
        if(ieqex .eq. 0)then
            write(LOT,*)'ieqex= ',ieqex,' EARTHQUAKE + EXPLOSION'
        else if(ieqex.eq.1)then
            write(LOT,*)'ieqex= ',ieqex,' POINT FORCES + EXPLOSION'
        else if(ieqex.eq.2)then
            write(LOT,*)'ieqex= ',ieqex,' ALL GREEN'
        else if(ieqex.eq.3)then
            write(LOT,*)'ieqex= ',ieqex,' EXPLOSION ONLY'
        else if(ieqex.eq.4)then
            write(LOT,*)'ieqex= ',ieqex,' EARTHQUAKE ONLY'
        else if(ieqex.eq.5)then
            write(LOT,*)'ieqex= ',ieqex,' POINT FORCES ONLY'
        else if(ieqex.eq.6)then
            write(LOT,*)'ieqex= ',ieqex,' SH ONLY'
        endif
c-----
c       provide names for output Green's functions in order of output
c-----
        if(ieqex.eq.0)then
c-----
c           EARTHQUAKE + EXPLOSION
c-----
            do 1234 i=1,21
                if(i.le.8)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1234       continue
            jsrc(13) = 1
            jsrc(14) = 1
            jsrc(16) = 1
            jsrc(17) = 1
            jsrc(18) = 1
            jsrc(19) = 1
            jsrc(20) = 1
            jsrc(21) = 1
        else if(ieqex.eq.1)then
c-----
c           POINT FORCES + EXPLOSION
c-----
            do 1235 i=1,21
                if(i.ge.7)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1235       continue
            jsrc(13) = 0
            jsrc(14) = 0
            jsrc(17) = 0
            jsrc(18) = 0
            jsrc(19) = 0
        else if(ieqex.eq.2)then
c-----
c           EXPLOSION ONLY
c-----
            do 1236 i=1,21
                jsrc(i) = 1
 1236       continue
        else if(ieqex.eq.3)then
            do 1237 i=1,21
                if(i.eq.7 .or. i.eq.8 .or. i.eq.16)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1237       continue
        else if(ieqex.eq.4)then
c-----
c           EARTHQUAKE ONLY
c-----
            do 1238 i=1,21
                if(i.le.6 .or. i.eq.13 .or. i.eq.14)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1238       continue
        else if(ieqex.eq.5)then
c-----
c           POINT FORCES ONLY
c-----
            do 1239 i=1,21
                if(i.ge.9)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1239       continue
            jsrc(13) = 0
            jsrc(14) = 0
            jsrc(16) = 0
        else if(ieqex.eq.6)then
c-----
c           SH ONLY
c-----
            do 1240 i=1,21
                if(i.eq.13 .or. i.eq.14 .or. i.eq.15)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1240       continue
        endif
c-----
c       jsrc - array giving Green's functions to be evaluated
c           this controls computations, gives far field terms
c           or course true solution for radial may involve transverse
c       jsrc(lsrc) maps into Green,s functions, e.g.,
c           For REP=10, lsrc(10) = 8, and if jsrc(8) = 1
c           P-SV contribution to explosion 
c           radial time history is computed
c
c       ieqex = 0 
c       Earthquake + Explosion
c       1-ZDD   2-RDD   3-ZDS   4-RDS   5-TDS   6-ZSS
c       7-RSS   8-TSS   9-ZEP   10-REP
c
c       ieqex = 1 
c       Point Forces + Explosion
c       11-ZVF  12-RVF  13-ZHF  14-RHF  15-THF  16-PEP
c       9-ZEP   10-REP (others have no meaning)
c
c       ieqex = 2
c       All Green's functions
c       1-ZDD   2-RDD   3-ZDS   4-RDS   5-TDS   6-ZSS
c       7-RSS   8-TSS   9-ZEP   10-REP
c       11-ZVF  12-RVF  13-ZHF  14-RHF  15-THF  16-PEP
c
c       ieqex = 3
c       Explosion Only
c       9-ZEP   10-REP
c
c       ieqex = 4 
c       Earthquake
c       1-ZDD   2-RDD   3-ZDS   4-RDS   5-TDS   6-ZSS
c       7-RSS   8-TSS
c
c       ieqex = 5 
c       Point Forces Only
c       11-ZVF  12-RVF  13-ZHF  14-RHF  15-THF  
c
c       ieqex = 6
c       SH Solution Only
c       5-TDS 8-TSS 15-THF
c
c       If fluid layer for receiver, 16 is forced to be fluid 
c           stress due to explosion
c-----
c       input jbdry = 10*surface + halfspace
c       surface   = 0 - elastic halfspace   = 0 - elastic
c                   1 - free              1 - free 
c                   2 - rigid             2 - rigid
c-----
        read(lun,*)jbdry
        write(LOT,21) jbdry
        if(jbdry.lt.0)jbdry = 0
        ibdrys = jbdry / 10
        if(ibdrys.eq.0)then
            jbdrys = 0
        else if(ibdrys.eq.1)then
            jbdrys = 1
        else if(ibdrys.eq.2)then
            jbdrys = -1
        else
            jbdrys = 0
        endif
        ibdryh = mod(jbdry,10)
        if(ibdryh.eq.0)then
            jbdryh = 0
        else if(ibdryh.eq.1)then
            jbdryh = 1
        else if(ibdryh.eq.2)then
            jbdryh = -1
        else
            jbdryh = 0
        endif
c-----
c       jbdrys  =  surface boundary condition
c           = -1 top surface is rigid
c           =  0 really a halfspace with parameters of top layer    
c           =  1 free surface
c       jbdryh  = halfspace boundary condition
c           = -1  RIGID
c           =  0  ELASTIC
c           = +1  FREE SURFACE
c-----
        if(jbdrys.eq.1)then
            write(LOT,*)' TOP  OF MODEL IS FREE SURFACE  '
        else if(jbdrys.eq.0)then
            write(LOT,*)' TOP  OF MODEL IS HALFSPACE WITH',
     1          ' PROPERTIES OF FIRST LAYER'
        else if(jbdrys.eq.-1)then
            write(LOT,*)' TOP  OF MODEL IS RIGID'
        endif
        if(jbdryh.eq.0)then
            write(LOT,*)' BASE OF MODEL IS HALFSPACE WITH',
     1          ' PROPERTIES OF BOTTOM LAYER'
        else if(jbdryh.eq.-1)then
            write(LOT,*)' BASE OF MODEL IS RIGID'
        else if(jbdryh.eq.1)then
            write(LOT,*)' BASE OF MODEL IS FREE'
        endif


c-----
c       read in the earth model name
c-----
        read(lun,'(a)')mname
        lmnm = lgstr(mname)
        write(LOT,*)mname(1:lmnm)
c-----
c       read in the earth model in model96 format
c-----
        inquire(file=mname,exist=ext)
        if(.not. ext)then
            write(LER,*)'Model file not located'
            stop
        endif
        call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.true.)
c-----
c       check the appropriateness of the model file
c-----
c               - -1 file does not exist
c               - -2 file is not a model file
c                 -3 error in the model file
        if(ierr.eq. -1)then
            call usage('Model file does not exist')
        else if(ierr.eq. -2)then
            call usage('Model file given is not a model96 file')
        else if(ierr.eq. -3)then
            call usage('Error in model file')
        endif
c-----
c       error checking
c-----
        if(idimen.ne.1)then
            call usage('1-D velocity model required')
        endif
        if(icnvel.ne.0)then
            call usage('Constant velocity model required')
        endif
        if(iiso.lt.0 .or. iiso.gt.1)then
            call usage('Isotropic or TI velocity model required')
        endif
        if(iflsph.ne.0)then
            call usage('Flat earth velocity model required')
        endif
c-----
c       do not permit Q < 1 If qa or qb is entered > 1 
c       invert to form q inverse
c-----
        do 3007 i=1,mmax
            if(qa(i) .gt. 1.0)qa(i) = 1.0/qa(i)
            if(qb(i) .gt. 1.0)qb(i) = 1.0/qb(i)
            if(frefp(i) .le. 0.0)frefp(i) = 1.0
            if(frefs(i) .le. 0.0)frefs(i) = 1.0
 3007   continue
        call modcpy(.true.)
        call velbnd()
c-----
c       check model for inconsistencies
c-----
        call chkmod()
c-----
c       read in controls for wavenumber integration
c-----
        read(lun,*) xleng, xfac
        write(LOT,24)xleng, xfac
c-----
c       read in the  source depth
c-----
        read(lun,*)mdpths
        mtmp = NSOURCE
        if(mdpths .gt. NSOURCE)then
            write(ostr,3008)mdpths,mtmp
 3008   format(' NUMBER OF SOURCE DEPTHS',i5,' EXCEEDS DIMENSION ',i5)
            call werror(ostr)
        endif
c-----
c       get maximum source or receiver depth
c-----
        do 2008 i=1,mdpths
            read(lun,*)depths(i)
            depths(i) = depths(i) + refdep
 2008   continue
c-----
c       read in the receiver depth
c-----
        read(lun,*)mdpthr
        mtmp = NRECEIVER
        if(mdpthr .gt. NRECEIVER)then
            write(ostr,3009)mdpthr,mtmp
 3009   format(' NUMBER OF RECEIVER DEPTHS',i5,' EXCEEDS DIMENSION ',i5)
            call werror(ostr)
        endif
        do 2009 i=1,mdpthr
            read(lun,*)depthr(i)
            depthr(i) = depthr(i) + refdep
 2009   continue
c-----
c       check for filling the final depth array
c-----
        mtmp = mdpths * mdpthr
        ntmp = NSR
        if(mtmp .gt. ntmp)then
            write(ostr,3011)mtmp,ntmp
 3011   format(' NUMBER SOURCE-RECEIVER COMB',
     1      i5,' EXCEEDS DIMENSION ',i5)
            call werror(ostr)
        endif
c-----
c       read in the distances
c-----
        read(lun,*)ndist
        mtmp = NDST
        if(mdpthr .gt. NDST)then
            write(ostr,3010)ndist,mtmp
 3010   format(' NUMBER OF DISTANCES',i5,' EXCEEDS DIMENSION ',i5)
            call werror(ostr)
        endif
        do 2010 i=1,ndist
            read(lun,*)rr,tshf,vr
            r(i)=rr
            tshift(i) = tshf
            vred(i) = vr
 2010   continue
c-----
c       get phase velocity limits
c-----
        read(lun,*)pmin,pmax,dp,pcntrl
        write(LOT,*)'[pmin,pmax,dp,pcntrl]=[' ,pmin, ',' ,pmax, ',' ,dp, 
     1      pcntrl, ']'
c-----
c       here pmin = mimimum ray parameter
c       pmax = maximum ray parameter
c       dp = ray parameter increment
c       pcntrl <= 0.0  modified for time series
c           > 0.0 true p-tau response
c-----
        if(pcntrl .le. 0.0)then
            write(LOT,*)' Modifed p-tau, for causal P arrivals'
        else
            write(LOT,*)' True p-tau'
        endif
        if(dp.le.0.0)dp = 1.0
c-----
c       For reasons of efficiency, decide whether to
c       add all layers at once, to the model
c       or to evaluate each layer source-receiver
c       combination separately. 
c
c       Roughly if the number of unique source and receiver depths
c       are mdpths+mdpthr if we insert layers, then we
c       end up with roughly mmax+mdpths+mdpthr layers, and
c       hence layer multiplication of this many matrices
c       for each source-receiver combination. Of course, for
c       equally spaced depth points, some economy arises
c       in avoiding matrix recomputation.
c
c       So if mdpths+mdpthr > 2*mmax we do not make a big model
c       other wise we do
c-----
c       adjust the model so that additional layers are added
c       to permit source and receiver at top of a give layer
c-----
c       lyrins = .true.
c       if(mdpths+mdpthr .gt. 2.0*mmax .and. lyrins.eq. .false.)then
        if(mdpths+mdpthr .gt. 2.0*mmax )then
            lyrins = .false.
            write(LOT,*)' LAYER INSERTION NOT DONE'
        else
            lyrins = .true. 
            write(LOT,*)' LAYER INSERTION DONE'
            do 2108 i=1,mdpths
                call insert(depths(i)+refdep)
 2108       continue
            do 2109 i=1,mdpthr
                call insert(depthr(i)+refdep)
 2109       continue
            call dezero()
c-----
c           check whether neighboring layers are identical
c           to avoid redundant evaluation
c----
            call equlck()
        endif
c-----
c       verify the new model parameters
c-----
        write(LOT,*)'mmax=',mmax
        write(LOT,22)(d(i),TA(i),TC(i),TF(i),TL(i),TN(i),
     1      Trho(i),qa(i),qb(i),i=1,mmax)
C-----
C removed 01 APR 2002 since messes up a mixed fluid/solid modium
Cc-----
Cc      Guarantee that no time wasted if any source in in the water
Cc      since there can only be a center of expansion source
Cc-----
C removed 01 APR 2002 since messes up a mixed fluid/solid modium
C       do 2019 i=1,mdpths
C           call srclay(depths(i), lmaxs(i), dphs)
C           if(b(lmaxs(i)).le.1.0e-04)then
C               do 2091 ii=1,21
C                   if(ii.ne.7 .and. ii.ne.8 
C     1                 .and. ii.ne.16)then
C                       jsrc(ii) = 0
C                   endif
C 2091          continue
C           endif
C 2019  continue
c-----
c       determine position of source and receive in new layers
c-----
        if(lyrins)then
            do 2020 i=1,mdpths
                call srclyr(depths(i), lmaxs(i), dphs)
 2020       continue
            do 2021 i=1,mdpthr
                call srclyr(depthr(i), lmaxr(i), dphs)
 2021       continue
        endif
        return
        end

        subroutine gcmdln()
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        integer*4 mnmarg
        character*50 name
        dokjar = .false.
        dosud = .false.
        spup = .false.
        spdn = .false.
        ssup = .false.
        ssdn = .false.
        dorud = .false.
        rpup = .false.
        rpdn = .false.
        rsup = .false.
        rsdn = .false.
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:3).eq.'-SU')then
                spup = .true.
                ssup = .true.
                ssdn = .false.
                spdn = .false.
                dosud = .true.
            else if(name(1:3).eq.'-SD')then
                spup = .false.
                ssup = .false.
                ssdn = .true.
                spdn = .true.
                dosud = .true.
            else if(name(1:5).eq.'-SPUP')then
                spup = .true.
                dosud = .true.
            else if(name(1:5).eq.'-SSUP')then
                ssup = .true.
                dosud = .true.
            else if(name(1:5).eq.'-SPDN')then
                spdn = .true.
                dosud = .true.
            else if(name(1:5).eq.'-SSDN')then
                ssdn = .true.
                dosud = .true.
            else if(name(1:5).eq.'-RPUP')then
                rpup = .true.
                dorud = .true.
            else if(name(1:5).eq.'-RSUP')then
                rsup = .true.
                dorud = .true.
            else if(name(1:5).eq.'-RPDN')then
                rpdn = .true.
                dorud = .true.
            else if(name(1:5).eq.'-RSDN')then
                rsdn = .true.
                dorud = .true.
            else if(name(1:3).eq.'-RD')then
                rpup = .false.
                rsup = .false.
                rsdn = .true.
                rpdn = .true.
                dorud = .true.
            else if(name(1:3).eq.'-RU')then
                rpup = .true.
                rsup = .true.
                rsdn = .false.
                rpdn = .false.
                dorud = .true.
            else if(name(1:2).eq.'-K')then
                dokjar = .true.
            else if(name(1:2).eq.'-?')then
                call usage(' ')
            else if(name(1:2).eq.'-h')then
                call usage(' ')
            endif
        goto 11
   13   continue
c-----
c       safety check
c-----
        return
        end

        subroutine usage(ostr)
c------
c       write out program syntax
c-----
        character ostr*(*)
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
        if(ostr.ne. ' ' )then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'USAGE: ',
     1  'tspec96p [-K] ',
     2      '[-SU] [-SD] [-SPUP] [-SSUP] [-SPDN] [-SSDN] ',
     3      '[-RU] [-RD] [-RPUP] [-RSUP] [-RPDN] [-RSDN] [-?] [-h]'
        write(LER,*)
     1  '-K      (default Futterman) use Kjartansson Causal Q'
        write(LER,*)
     1  'The following govern wavefield at source. The default is',
     2  ' the entire wavefield'
        write(LER,*)
     1  '-SU      (default whole wavefield) Compute only upgoing ',
     2  '               wavefield from the source'
        write(LER,*)
     1  '-SD      (default whole wavefield) Compute only downgoing ',
     2  '               wavefield from the source'
        write(LER,*)
     1  ' -SPUP  Include upward P at source'
        write(LER,*)
     1  ' -SSUP  Include upward S at source'
        write(LER,*)
     1  ' -SPDN  Include downward P at source'
        write(LER,*)
     1  ' -SSDN  Include downward S at source'
        write(LER,*)
     1  'The following govern wavefield at receiver. The default is',
     2  ' the entire wavefield'
        write(LER,*)
     1  ' -RD    Only downgoing waves at receiver'
        write(LER,*)
     1  ' -RU    Only upgoing waves at receiver'
        write(LER,*)
     1  ' -RPUP  Include upward P at receiver'
        write(LER,*)
     1  ' -RSUP  Include upward S at receiver'
        write(LER,*)
     1  ' -RPDN  Include downward P at receiver'
        write(LER,*)
     1  ' -RSDN  Include downward S at receiver'
        write(LER,*)
     1  '-?                   Display this usage message'
        write(LER,*)
     1  '-h                   Display this usage message'
        stop 
        end

        subroutine aten(om,qa,qb,alpha,atna,atnb,frefp,frefs)
        implicit none
c-----
c       om  C*16    complex frequency
c       
c-----
c       make velocities complex, using Futterman or Kjartansson
c       causality operator
c-----
        complex*16 om, atna, atnb
        real*4 qa, qb, frefp, frefs, alpha

        complex*16 at
        real*8 pi, om1p, om1s, oml, fac, pi2
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        real gama, gamb, rfac
        real*8 CDABS
        complex*16 CDLOG
c-----
c       reference frequency is fref hz
c-----
        om1p=6.2831853*frefp
        om1s=6.2831853*frefs
        pi2 = 1.5707963
        pi=3.1415927d+00
        if(dokjar)then
c-----
c       Kjartansson Constant Q, causal Q operator
c       Kjartansson, E. (1979). 
c           Constant Q-wave propagation and attenuation,
c       J. Geophys. Res. 84, 4737-4748.
c-----
            gama = atan(qa)/pi
            gamb = atan(qb)/pi
            if(gama.le.0.0)then
                atna = dcmplx(1.0d+00,0.0d+00)
            else
                fac = pi2*gama
                rfac = sin(fac)/cos(fac)
                atna = dcmplx(1.0d+00,0.0d+00)/
     1              (( (om/om1p)**dble(-gama) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
            endif
            if(gamb.le.0.0)then
                atnb = dcmplx(1.0d+00,0.0d+00)
            else
                fac = pi2*gamb
                rfac = sin(fac)/cos(fac)
                atnb = dcmplx(1.0d+00,0.0d+00)/
     1          (( (om/om1s)**dble(-gamb) ) *
     2           dcmplx(1.0d+00,-dble(rfac)))
            endif
        else
c-----
c       Futterman Causal Q
c-----
c           low frequency cutoff is 0.01 hz
c-----
            oml=0.062831853d+00
            atna=dcmplx(1.0d+00,0.0d+00)
            atnb=dcmplx(1.0d+00,0.0d+00)
            if(qa.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(om).gt.oml) at=CDLOG(om/om1p)/pi
              if(CDABS(om).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1p*fac))/pi
              endif
              atna=(1.+dble(qa)*at+dcmplx(0.0d+00,dble(qa/2.)))
            endif
            if(qb.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(om).gt.oml) at=CDLOG(om/om1s)/pi
              if(CDABS(om).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1s*fac))/pi
              endif
               atnb=(1.+dble(qb)*at+dcmplx(0.0d+00,dble(qb/2.)))
            endif
        endif
        return
        end
        subroutine bufini(irdwr,ierr)
c-----
c       initialize buffer pointer
c       irdwr = 0 read initialize
c       irdwr = 1 write initialize
c-----
        implicit none
        integer irdwr, ierr
        integer BUFMAX
        parameter(BUFMAX=2000)
        common/buf/iptr,maxb,buffer(BUFMAX)
        integer iptr, maxb
        real buffer
        iptr = 1
        if(irdwr.eq.0)call getbuf(ierr)
        return
        end

        subroutine buflsh()
c-----
c       flush output buffer
c-----
        implicit none
        integer BUFMAX
        parameter(BUFMAX=2000)
        common/buf/iptr,maxb,buffer(BUFMAX)
        integer iptr, maxb
        real buffer
        integer ipt
        integer i
        ipt = iptr -1
        if(ipt.gt.0)write(2)ipt,(buffer(i),i=1,ipt)
        iptr = 1
        return
        end

        subroutine bufwr(x)
c-----
c       fill buffer with floating point variable x,
c       flush buffer as necessary
c-----
        implicit none
        real x
        integer BUFMAX
        parameter(BUFMAX=2000)
        common/buf/iptr,maxb,buffer(BUFMAX)
        integer iptr, maxb
        real buffer
        buffer(iptr) = x
        iptr = iptr + 1
        if(iptr.gt.BUFMAX)call buflsh
        return
        end

        subroutine getbuf(ierr)
c-----
c       read in file contents into buffer, taking care not to
c       read beyond the contents of the file
c-----
        implicit none
        integer ierr
        integer BUFMAX
        parameter(BUFMAX=2000)
        common/buf/iptr,maxb,buffer(BUFMAX)
        integer iptr, maxb
        real buffer
        integer i
c-----
c       ierr = 0 successful read
c            = 1 read error
c            = 2 end of file
c-----
        read(2,err=1000,end=2000)maxb,(buffer(i),i=1,maxb)
        iptr = 1
        ierr = 0
        return
 1000   ierr = 1
        return
 2000   ierr = 2
        return
        end

        subroutine bufrd(x,ierr)
c-----
c       retrieve a value from buffer array, red in new array
c       as necessary
c       iptr is here the next array element to be read
c       it is always >= 1. We do not worry the upper limit
c       since the calling program must worry about this
c       because read always follows a complete write
c-----
        implicit none
        real x
        integer ierr
        integer BUFMAX
        parameter(BUFMAX=2000)
        common/buf/iptr,maxb,buffer(BUFMAX)
        integer iptr, maxb
        real buffer
c-----
c       only yank in new data if actually required
c-----
        if(iptr.gt.maxb)call getbuf(ierr)
        x = buffer(iptr)
        iptr = iptr + 1
        return
        end

        subroutine chkmod()
        implicit none
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        integer jsrc, jbdrys, jbdryh
        common/lwater/lfluid
        logical lfluid

        integer iw, isoldd, isoldu, i
c-----
c       check model for inconsistencies
c-----
c-----
c       Model cannot consist entirely of water layers
c       Also determine first solid layer from top
c-----
        iw = 0  
        isoldd = 0
        do 100 i=1,mmax
            if(TL(i).eq.0.0 .or. TN(i).eq.0.0)then
                iw = iw + 1
            else
                if(isoldd .eq.0)isoldd=i
            endif
  100   continue
        if(iw .eq. mmax)then
            lfluid = .true.
C           call werror('MODEL CONSISTS ONLY OF LIQUID LAYERS')
        else
            lfluid = .false.
        endif
c-----
c       Determine first solid layer from bottom
c-----
        iw = 0  
        isoldu = 0
        do 101 i=mmax,1,-1
            if(TL(i).eq.0.0 .or. TN(i).eq.0.0)then
                iw = iw + 1
            else
                if(isoldu .eq.0)isoldu=i
            endif
  101   continue
c-----
c       Check for interior water layer
c-----
        if(iw.gt.0 .and. .not. lfluid)then
            do 102 i=isoldd,isoldu
                if(TL(i).eq.0.0.or.TN(i).eq.0.0)then
                call werror('MODEL HAS INTERIOR  FLUID LAYERS')
                endif
  102       continue
        endif
c-----
c       If boundary condition is rigid, and the adjacent layer is
c       fluid, terminate 
c-----
C       if(TN(1).le.1.0e-04 .and. jbdrys.eq.-1 .and. lfluid)then
C           call werror('TOP LAYER IS FLUID AND RIGID BOUNDARY')
C       endif
C       if(TN(mmax).le.1.0e-04 .and. jbdryh.eq.-1 .and. lfluid)then
C           call werror('BOTTOM LAYER IS FLUID AND RIGID BOUNDARY')
C       endif
C       if(TL(1).le.1.0e-04 .and. jbdrys.eq.-1 .and. lfluid)then
C           call werror('TOP LAYER IS FLUID AND RIGID BOUNDARY')
C       endif
C       if(TL(mmax).le.1.0e-04 .and. jbdryh.eq.-1 .and. lfluid)then
C           call werror('BOTTOM LAYER IS FLUID AND RIGID BOUNDARY')
C       endif
        return
        end

        subroutine cmult(e,ca,exa,exe)
        implicit none
        common/lwater/lfluid
        logical lfluid
c-----
c       FORM EC where e(1x5) c(5x5)
        complex*16 e(5)
        complex*16 ca(5,5)
        real*8 exa,exe,eval
        real *8 xnorm
        complex*16 c, ee(5)
        integer IUP, i, j

        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 1350 i=1,IUP
            c = dcmplx(0.0d+00,0.0d+00)
            do 1349 j=1,IUP
                c=c+ e(j) * ca(j,i)
 1349       continue
            ee(i)=c
 1350   continue
        exe = exe + exa
        call normc(ee,eval,xnorm)
        do 1351 i=1,IUP
            e(i) = ee(i)*xnorm
 1351   continue
        exe = exe + eval
        return
        end

        subroutine copy5(ca,dar,m,itofrm,dex,exa)
        implicit none
        integer NL
        parameter (NL=200)
        complex*16 ca(5,5)
        complex*16 dar(NL,5,5)
        integer m, itofrm
        real*8 dex(NL)
        real*8 exa

        integer i, j
c-----
c       copy from ca to dar
c-----
        if(itofrm.eq.0)then
            do 100 j=1,5
                do 110 i=1,5
                    dar(m,i,j) = ca(i,j)
  110           continue
                dex(m) = exa
  100       continue
c-----
c       copy from dar to ca
c-----
        else
            do 200 j=1,5
                do 210 i=1,5
                    ca(i,j) = dar(m,i,j)
  210           continue
                exa = dex(m)
  200       continue
        endif
        return
        end

        subroutine copy4(aa,har,m,itofrm,hex,ex)
c-----
c       copy between aa and har arrays for m'th layer
c
c       aa  R   - 4x4 Haskell matrix array
c       har R   - NLx4x4 storage array
c       m   I   - layer index
c       itofrm  I   - 0 copy aa to har, ex to hex
c                 !=0 copy from har to aa, hex to ex
c       hex R   - NL array - storage for exponent
c       ex  R   - exponent value
c-----
        implicit none
        integer NL
        parameter (NL=200)
        complex*16 aa(4,4)
        complex*16 har(NL,4,4)
        integer m, itofrm
        real*8 hex(NL)
        real*8 ex

        integer i, j
c-----
c       copy from aa to har
c-----
        if(itofrm.eq.0)then
            do 100 j=1,4
                do 110 i=1,4
                    har(m,i,j) = aa(i,j)
  110           continue
                hex(m) = ex
  100       continue
c-----
c       copy from har to aa
c-----
        else
            do 200 j=1,4
                do 210 i=1,4
                    aa(i,j) = har(m,i,j)
  210           continue
                ex = hex(m)
  200       continue
        endif
        return
        end

        subroutine copy2(hl,hal,m,itofrm,lex,exb)
        implicit none
        integer NL
        parameter (NL=200)
        complex*16 hl(2,2)
        complex*16 hal(NL,2,2)
        integer m, itofrm
        real*8 lex(NL)
        real*8 exb

        integer i, j
c-----
c       copy from hl to hal
c-----
        if(itofrm.eq.0)then
            do 100 j=1,2
                do 110 i=1,2
                    hal(m,i,j) = hl(i,j)
  110           continue
                lex(m) = exb
  100       continue
c-----
c       copy from hal to hl
c-----
        else
            do 200 j=1,2
                do 210 i=1,2
                    hl(i,j) = hal(m,i,j)
  210           continue
                exb = lex(m)
  200       continue
        endif
        return
        end

        subroutine dezero()
        implicit none
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        integer i
        real dmin
c-----
c       ultimately get rid of zero thickness layers - this
c       will require readjusting the model from top down, and
c       also readjusting the source and receiver indices.
c----
c       Here just guarantee that the halfspace is not of zero thickness
c-----
        dmin = 1.0e+30
        do 100 i=1,mmax-1
            if(d(i) .lt. dmin .and. d(i).gt.0.0)dmin = d(i)
  100   continue
c       if(d(mmax).le.0.0)then
c           d(mmax) = 0.1*dmin
c       endif
        return
        end

        subroutine dmult(da,aa)
c-----
c       propagator up
c       FORM D = DA
c-----
        implicit none
        complex*16 aa(4,4)
        complex*16 sumd,ea(4,4),da(4,4)
        integer i, j, jj
        do 1360 i=1,4
            do 1361 j=1,4
                sumd = dcmplx(0.0d+00,0.0d+00)
                do 1362 jj=1,4
                    sumd=sumd+da(i,jj) * aa(jj,j)
 1362           continue
                ea(i,j)=sumd
 1361       continue
 1360   continue
        do 1363 j=1,4
            do 1364 i=1,4
                da(i,j)=ea(i,j)
 1364       continue
 1363   continue
        return
        end

        subroutine dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,NP,NSV,
     1      x11,x21,x31,x41,x12,x22,x32,x42,
     1      rho,iwat,ex,om2)
        implicit none
        complex*16 CA(5,5)
        COMPLEX*16 cosp , cossv 
        COMPLEX*16 rsinp, rsinsv
        COMPLEX*16 sinpr, sinsvr
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        real rho
        integer iwat
        real*8 ex, dfac
        complex*16  om2
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
c       the number of multiplications can be reduced from 36 to 25 
c           if we define a new matrices
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
        zrho = dcmplx(dble(rho),0.0d+00)
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
                    ca(i,j) = dcmplx(0.0d+00, 0.0d+00)
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
            CA(1,1) = TCA(1,1)
            CA(1,2) = TCA(1,2)
            CA(1,3) = TCA(1,3)
            CA(1,4) = TCA(1,5)
            CA(1,5) = TCA(1,6)

            CA(2,1) = TCA(2,1)
            CA(2,2) = TCA(2,2)
            CA(2,3) = TCA(2,3)
            CA(2,4) = TCA(2,5)
            CA(2,5) = TCA(1,5)

            CA(3,1) =   2.0d+00*TCA(3,1)
            CA(3,2) =   2.0d+00*TCA(3,2)
c-----
c       beware of normalization
c-----
            CA(3,3) =   2.0d+00*TCA(3,3) - 1.0d+00*dfac
            CA(3,4) =  -2.0d+00*TCA(2,3)
            CA(3,5) =  -2.0d+00*TCA(1,3)

            CA(4,1) =   TCA(5,1)
            CA(4,2) =   TCA(5,2)
            CA(4,3) =  -TCA(3,2)
            CA(4,4) =   TCA(2,2)
            CA(4,5) =   TCA(1,2)

            CA(5,1) =   TCA(6,1)
            CA(5,2) =   TCA(5,1)
            CA(5,3) =  -TCA(3,1)
            CA(5,4) =   TCA(2,1)
            CA(5,5) =   TCA(1,1)
        endif
        return
        end

        subroutine equlck()
        implicit none
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald

        integer m
c-----
c       To avoid repeated computation, check to see 
c           if neighboring layers
c       are identical, once for going up and another for going down
c-----
c       First check top down
c-----
        do 100 m=1,mmax
            if(m.eq.1)then
                equald(m) = .false.
            else if(m.gt.1
     1          .and. TA(m).eq.TA(m-1) 
     1          .and. TC(m).eq.TC(m-1) 
     1          .and. TF(m).eq.TF(m-1) 
     1          .and. TL(m).eq.TL(m-1) 
     1          .and. TN(m).eq.TN(m-1) 
     1          .and.  D(m).eq. D(m-1) 
     4          .and. TRho(m).eq.TRho(m-1)
     5          .and. qa(m).eq.qa(m-1)
     6          .and. qb(m).eq.qb(m-1) )then
                equald(m) = .true.
            else
                equald(m) = .false.
            endif
  100   continue
c-----
c       check bottom up
c-----
        do 200 m=1,mmax
            if(m.eq.mmax)then
                equalu(m) = .false.
            else if(m.lt.mmax
     1          .and. TA(m).eq.TA(m+1) 
     1          .and. TC(m).eq.TC(m+1) 
     1          .and. TF(m).eq.TF(m+1) 
     1          .and. TL(m).eq.TL(m+1) 
     1          .and. TN(m).eq.TN(m+1) 
     1          .and.  D(m).eq. D(m+1) 
     4          .and. Trho(m).eq.Trho(m+1)
     5          .and. qa(m).eq.qa(m+1)
     6          .and. qb(m).eq.qb(m+1) )then
                equalu(m) = .true.
            else
                equalu(m) = .false.
            endif
  200   continue
        return
        end

        subroutine evalg(jbdry,m,m1,gbr,gbl,in,
     1      wvno,om,om2,wvno2)
        implicit none
        integer jbdry, m, m1, in
        complex*16 gbr(2,5), gbl(2,2)
        complex*16 wvno,om,wvno2,om2
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
        common/lwater/lfluid
        logical lfluid
        complex*16 CDSQRT

        complex*16 atna, atnb
        complex*16 cg(6)
        complex*16 g(4,4)
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv, rsh
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        integer iwat


c-----
c       set up halfspace conditions
c-----
            if(TL(m).eq. 0.0 .or.TN(m).eq.0.0)then
                iwat = 1
            else
                iwat = 0
            endif
c-----
c       set up halfspace boundary conditions
c
c       jbdry   = -1  RIGID
c           =  0  ELASTIC
c           = +1  FREE SURFACE
c
c-----
        if(jbdry.lt.0)then
c-----
c       RIGID - check properties of layer above
c-----
            if(TL(m) .gt. 0.0 .and. TN(m).gt.0.0)then
c-----
c               ELASTIC ABOVE - RIGID
c-----
                gbr(in,1) = dcmplx(1.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,1) = dcmplx(1.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(0.0d+00,0.0d+00)
            else
c-----
c               FLUID ABOVE - RIGID
c-----
                gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                if(lfluid)then
                    gbr(in,1) = dcmplx(1.0d+00,0.0d+00)
                else
                    gbr(in,4) = dcmplx(1.0d+00,0.0d+00)
                endif
c-----
c               (pseudo SH)
c-----
                gbl(in,1) = dcmplx(1.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(0.0d+00,0.0d+00)
            endif
        else if(jbdry.eq.0)then
c-----
c       HALFSPACE
c-----
            call aten(om,qa(m),qb(m),alpha,atna,atnb,frefp(m),frefs(m))
            call getegn(za,zb,zc,zd,ze,zf,om2,wvno2,rp, rsv, rsh,
     1          x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,om,wvno,
     2          atna,atnb)
            if(iwat.eq.0)then
c-----
c               ELASTIC HALFSPACE
C       WRITE(0,*)'X:'
C       WRITE(0,*)x11,x21,x31,x41
C       WRITE(0,*)x12,x22,x32,x42
C       WRITE(0,*)'RP,RSV,RSH:',RP,RSV,RSH
C       WRITE(0,*)'NP,NSV,m,om,wvno:',NP,NSV,m,om,wvno
C       WRITE(0,*)'om2,wvno2:',om2,wvno2


                G(1,1) =    x41/(2.0*NP  *rp  )
                G(1,2) = -  x31/(2.0*NP       )
                G(1,3) = -  x21/(2.0*NP  *rp  )
                G(1,4) =    x11/(2.0*NP       )
                G(2,1) =    x42/(2.0*NSV      )
                G(2,2) = -  x32/(2.0*NSV *rsv )
                G(2,3) = -  x22/(2.0*NSV      )
                G(2,4) =    x12/(2.0*NSV *rsv )
                G(3,1) = -  x41/(2.0*NP  *rp  )
                G(3,2) = -  x31/(2.0*NP       )
                G(3,3) =    x21/(2.0*NP  *rp  )
                G(3,4) =    x11/(2.0*NP       )
                G(4,1) =    x42/(2.0*NSV      )
                G(4,2) =    x32/(2.0*NSV *rsv )
                G(4,3) = -  x22/(2.0*NSV      )
                G(4,4) = -  x12/(2.0*NSV *rsv )
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

                gbl(in,1) = atnb*atnb*TL(m)*rsh
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
            else if(iwat.eq.1)then
c-----
c               FLUID HALFSPACE
c-----
                if(lfluid)then
                    gbr(in,1) = dble(TRho(m))*om2
                    gbr(in,2) = -rp
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                else
                    gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = -dble(Trho(m))*om2
                    gbr(in,5) = rp
                endif
                gbl(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
            endif
        else if(jbdry.eq.1)then
c-----
c       FREE - check properties of layer above
c-----
            if(TL(m) .gt. 0.0 .and. TN(m).gt.0.0)then
                gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(1.0d+00,0.0d+00)
                gbl(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
                
            else
                gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                if(lfluid)then
                    gbr(in,2) = dcmplx(1.0d+00,0.0d+00)
                else
                    gbr(in,5) = dcmplx(1.0d+00,0.0d+00)
                endif
                gbl(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
            endif
        endif
        return
        end

        subroutine evalh(jbdry,m,m1,hsr,hsl,in,wvno,om,om2,wvno2)
        implicit none
        integer jbdry,m,m1,in
        complex*16 hsr(2,5), hsl(2,2)
        complex*16 wvno,om,wvno2,om2
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
        complex*16 CDSQRT

        complex*16 atna, atnb
        complex*16 HM(4,4)
        complex*16 CH(6)
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv, rsh
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        integer iwat


c-----
c       set up surface conditions
c-----
        call aten(om,qa(m),qb(m),alpha,atna,atnb,frefp(m),frefs(m))
            if(TL(m).eq. 0.0 .or.TN(m).eq.0.0)then
                iwat = 1
            else
                iwat = 0
            endif
        
c-----
c       do top surface condition
c
c           = -1  RIGID
c           =  0  ELASTIC
c           = +1  FREE SURFACE
c-----
        if(iwat.eq.0)then
            if(jbdry.eq.0)then
c-----
c           ELASTIC ABOVE SOLID
c-----
c-----
c       Now give the H matrix
c-----
        call getegn(za,zb,zc,zd,ze,zf,om2,wvno2,rp, rsv, rsh,
     1  x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,om,wvno,
     2  atna,atnb)
                HM(1,1) =   x11*rp
                HM(1,2) =   x12
                HM(1,3) = - x11*rp
                HM(1,4) =   x12
                HM(2,1) =   x21
                HM(2,2) =   x22*rsv
                HM(2,4) =  -x22*rsv
                HM(2,3) =   x21
                HM(3,1) =   x31*rp
                HM(3,2) =   x32
                HM(3,3) = - x31*rp
                HM(3,4) =   x32
                HM(4,1) =   x41
                HM(4,2) =   x42*rsv
                HM(4,3) =   x41
                HM(4,4) =  -x42*rsv
C Compound matrices
c CH(1) = H 12 12 = 11 22 - 12 21
c CH(2) = H 13 12 = 11 32 - 12 31
c CH(3) = H 14 12 = 11 42 - 12 41
c CH(4) = H 23 12 = 21 32 - 31 22
c CH(5) = H 24 12 = 21 42 - 41 22
c CH(6) = H 34 12 = 31 42 - 41 32
                CH(1) = HM(1,1)*HM(2,2) - HM(1,2)*HM(2,1)
                CH(2) = HM(1,1)*HM(3,2) - HM(1,2)*HM(3,1)
                CH(3) = HM(1,1)*HM(4,2) - HM(1,2)*HM(4,1)
                CH(4) = HM(2,1)*HM(3,2) - HM(3,1)*HM(2,2)
                CH(5) = HM(2,1)*HM(4,2) - HM(4,1)*HM(2,2)
                CH(6) = HM(3,1)*HM(4,2) - HM(4,1)*HM(3,2)
                hsr(in,1) =         CH(1)
                hsr(in,2) =         CH(2)
                hsr(in,3) = 2.0d+00*CH(3)
                hsr(in,4) =         CH(5)
                hsr(in,5) =         CH(6)
                hsl(in,1) = dcmplx(1.0d+00,0.0d+00)
                hsl(in,2) = atnb*atnb*TL(m)*rsh
            else if(jbdry.eq.-1)then
c-----
c           RIGID ABOVE SOLID
c-----
                hsr(in,1) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(1.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(0.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(1.0d+00,0.00d+00)
            else if(jbdry.eq.1)then
c-----
c           FREE ABOVE SOLID
c-----
                hsr(in,1) = dcmplx(1.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            endif
        else if(iwat.gt.0)then
            if(jbdry.eq.0)then
c-----
c           HALFSPACE ABOVE FLUID
c-----
                hsr(in,1) = rp
                hsr(in,2) = -dble(TRho(m))*om2
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            else if(jbdry.eq.-1)then
c-----
c           RIGID ABOVE FLUID
c-----
                hsr(in,1) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(1.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            else if(jbdry.eq.1)then
c-----
c           FREE ABOVE FLUID
c-----
                hsr(in,1) = dcmplx(1.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            endif
        endif
        return
        end

        subroutine evlmat(om,wvno,jbdrys,jbdryh,wvno2,om2)
        implicit none
        complex*16 om, wvno, om2, wvno2
        integer jbdrys,jbdryh
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        integer NSOURCE, NRECEIVER, NSR
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        real*4 depths
        integer lmaxs, mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depthr
        integer lmaxr, mdpthr
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
        complex*16 p,q,r
        complex*16  ca(5,5), hl(2,2)
        complex*16 aa(4,4)
        complex*16 zone
c-----
c       matrix components in layers and boundaries saved
c-----
        common/hamat/har
        complex*16 har(NL,4,4)
        common/damat/dar
        complex*16 dar(NL,5,5)
        common/hsrfr/hsr
        complex*16 hsr(2,5)
        common/gbrfr/gbr
        complex*16 gbr(2,5)
        common/hlmat/hal
        complex*16 hal(NL,2,2)
        common/hsrfl/hsl
        complex*16 hsl(2,2)
        common/gbrfl/gbl
        complex*16 gbl(2,2)
        common/hexex/hex
        real*8 hex(NL)
        common/hexexw/hexw
        real*8 hexw(NL)
        common/dexex/dex
        real*8 dex(NL)
        common/lexex/lex 
        real*8 lex(NL)
        common/water/iwater(NL),iwats(2),iwatb(2)
        integer iwater, iwats, iwatb
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
        logical compute
        complex*16 CDSQRT

        complex*16 atna,atnb
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv, rsh
        real*8 pex,svex,shex, cpex
        COMPLEX*16 cosp , cossv , cossh
        COMPLEX*16 rsinp, rsinsv, rsinsh
        COMPLEX*16 sinpr, sinsvr, sinshr
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42


        real*8 dfac
        integer m, iwat
        zone  = dcmplx(1.0d+00,0.0d+00)
c-----
c       evaluate the G matrix 
c           gbr(1,x) is for normal stack
c           gbr(2,x) is for inverted stack
c-----
        call evalg(jbdryh,mmax,mmax-1,gbr,gbl,1,wvno,om,om2,wvno2)
        call evalg(jbdrys,1,   2,     gbr,gbl,2,wvno,om,om2,wvno2)
c-----
c       evaluate the H matrix
c           hsr(1,x) is for normal stack
c           hsr(2,x) is for inverted stack
c-----
        call evalh(jbdrys,1,   2,     hsr,hsl,1,wvno,om,om2,wvno2)
        call evalh(jbdryh,mmax,mmax-1,hsr,hsl,2,wvno,om,om2,wvno2)
c-----

c-----
c       matrix multiplication from bottom layer upward
c-----
        do 1340 m = 1,mmax,1
c-----
c       first check to see if computations already done
c-----
            if(equald(m))then
                compute = .false.
            else
                compute = .true.
            endif
            if(compute)then
            if(TL(m).eq. 0.0 .or.TN(m).eq.0.0)then
                iwat = 1
            else
                iwat = 0
            endif
                call aten(om,qa(m),qb(m),alpha,atna,atnb,
     1          frefp(m),frefs(m))
                call getegn(za,zb,zc,zd,ze,zf,om2,wvno2,rp, rsv, rsh,
     1          x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,om,wvno,
     2          atna,atnb)
c-----
c       now evaluate the various hyperbolic functions
c-----
                p = rp  * dble(d(m))
                q = rsv * dble(d(m))
                r = rsh * dble(d(m))
                call var(p,q,r, rp, rsv, rsh,
     1              cosp, cossv, cossh, rsinp, rsinsv, rsinsh,
     1              sinpr, sinsvr, sinshr,pex,svex,shex,iwat)
                call dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,
     1              NP,NSV,
     1              x11,x21,x31,x41,x12,x22,x32,x42,
     1              Trho(m),iwat,pex+svex,om2)
c-----
c       For isotropic, Re p > Re sv do adjsut only the
c           SV part of the Haskell matrix
c       For TI this is not always true, so carefully adjust
c-----
            if(pex .gt. svex)then
c-----
c               PEX > SVEX, factor out PEX
c-----
                if((pex-svex).gt. 40.0d+00)then
                    dfac = 0.0d+00
                else
                    dfac = dexp(-(pex-svex))
                endif
                cpex = pex
                call hska(AA,cosp,rsinp,sinpr,
     1              dfac*cossv,dfac*rsinsv,dfac*sinsvr,NP,NSV,
     1              X11, X21, X31, X41,X12, X22, X32, X42, 
     2              Trho(m),iwat,pex,om2)
        else
c-----
c               SVEX > PEX, factor out SVEX
c-----
                if((svex-pex).gt. 40.0d+00)then
                    dfac = 0.0d+00
                else
                    dfac = dexp(-(svex-pex))
                endif
                cpex = svex
                call hska(AA,dfac*cosp,dfac*rsinp,dfac*sinpr,
     1              cossv,rsinsv,sinsvr,NP,NSV,
     1              X11, X21, X31, X41,X12, X22, X32, X42, 
     2              Trho(m),iwat,pex,om2)
        endif
                call hskl(cossh,rsinsh,sinshr,TL(m),
     1              iwat,hl,atnb)
            endif
            iwater(m) = iwat
            call copy5(ca,dar,m,0,dex,pex+svex)
            call copy4(aa,har,m,0,hex,cpex)
            call copy2(hl,hal,m,0,lex,shex)
 1340   continue
        return
        end

        subroutine excit(freq,xleng,xfac,dk,nk,omega)
c-----
c     sample response for all wavenumbers at a given frequency
c     using Bouchon equal wavenumber sampling = dk
c     with offset of 0.218dk
c-----
        implicit none

        real freq, xleng, xfac, dk, omega
        integer nk

        integer LER, LIN, LOT
        parameter(LER=0,LIN=5,LOT=6)
        integer NL
        parameter(NL=200)
        integer NSOURCE, NRECEIVER, NSR
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        real vmin,vamin,vamax,vbmin,vbmax
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        integer jsrc, jbdrys, jbdryh
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
c-----
c       set up common blocks for wavenumber sampling at
c       suitable depths. This is necessary since the matrix
c       evaluation is done here for all source-receiver pairs
c       The source-receiver distance is important for the
c       wavenumber sampling at low frequencies
c-----
        common/kint1/gasymp
            logical gasymp(NSR)
        common/kint2/mkup
            integer mkup(NSR)
        common/kint3/wave
            real*4 wave(NSR,2)
        common/kint4/aa,bb,cc
            complex aa(NSR,21),bb(NSR,21),cc(NSR,21)

        complex*16 wvn,om, wvn2, om2
        complex gg(21)
        complex zeye
c-----
c       matrix components in layers and boundaries saved
c-----
        common/hamat/har
        complex*16 har(NL,4,4)
        common/damat/dar
        complex*16 dar(NL,5,5)
        common/hsrfr/hsr
        complex*16 hsr(2,5)
        common/gbrfr/gbr
        complex*16 gbr(2,5)
        common/hlmat/hal
        complex*16 hal(NL,2,2)
        common/hsrfl/hsl
        complex*16 hsl(2,2)
        common/gbrfl/gbl
        complex*16 gbl(2,2)
        common/hexex/hex
        real*8 hex(NL)
        common/hexexw/hexw
        real*8 hexw(NL)
        common/dexex/dex
        real*8 dex(NL)
        common/lexex/lex 
        real*8 lex(NL)
        common/water/iwater(NL),iwats(2),iwatb(2)
        integer iwater, iwats, iwatb
        common/lyrctl/lyrins
        logical lyrins
        common/c/pmin,pmax,dp,pcntrl
        real pmin,pmax,dp,pcntrl


        integer ii, js, jr, k, j
        real p
        real dphr, dphs


        zeye = cmplx(0.0,1.0)
        omega=6.2831853*freq
        om=dcmplx(dble(omega),-dble(alpha))
        om2 = om * om
c-----
c       output by ray parameter
c-----
        nk = (pmax - pmin)/dp
        if(nk.lt.0)nk = 0
        nk = nk + 1
        do 3998 ii=1,nk
            p  = (ii-1)*dp + pmin
c-----
c       safety
c-----
            if(p.eq.0.0)p = 0.01*dp
            wvn=dcmplx(dble(p),0.0d+00)*om
            wvn2 = wvn*wvn
c-----
c       evaluate matrices first
c-----
            if(lyrins)then
                call evlmat(om,wvn,jbdrys,jbdryh,wvn2,om2)
            endif
c-----
c       now evaluate for a specific source, receiver position
c-----
            k = 0
            do 4000 js=1,mdpths
                do 4010 jr=1,mdpthr
                    k = k + 1
                    if(.not.lyrins)then
c-----
c                   evaluate matrices first
c                   for currently defined layering
c-----
                    call modcpy(.false.)
                    call insert(depths(js))
                    call insert(depthr(jr))
                    call srclyr(depths(js), 
     1                  lmaxs(js), dphs)
                    call srclyr(depthr(jr), 
     1                  lmaxr(jr), dphr)
                    call dezero()
                    call equlck()
                call evlmat(om,wvn,jbdrys,jbdryh,wvn2,om2)
                    endif
                    call rshof(gg,om,wvn,lmaxs(js),
     1                  lmaxr(jr),wvn2,om2)
c-----
c       To make radial look pulse like for small ray parameters
c       Also take time derivative of the point force and pressure fields
c       The multiplication by 'i' accomplishes this for 
c           responses with a linear term in wavenumber
c       Note that these are the integrands and not the 
c           output Green s functions. Thus
c       the TDS contribution (output=5) is integrand (13)
c-----
                    if(pcntrl .lt.0.0)then
                        GG(2)  = - GG(2) * zeye
                        GG(3)  =   GG(3) * zeye
                        GG(5)  = - GG(5)
                        GG(6)  =   GG(6) * zeye
                        GG(14) =   GG(14) * zeye
                        GG(8)  = - GG(8) * zeye
                        GG(9)  =   GG(9) * zeye * om
                        GG(10) = - GG(10) * zeye * om * zeye
                        GG(11) =   GG(11) * zeye * om * zeye
                        GG(12) =   GG(12) * zeye * om 
                        GG(15) =   GG(15) * zeye * om 
                    endif
c-----
c       output
c-----
                    do 3999 j=1,21
                        if(jsrc(j).eq.1)then
                            call bufwr(real(gg(j)))
                            call bufwr(aimag(gg(j)))
                        endif
 3999               continue
 4010           continue
 4000       continue
 3998   continue
        return
        end

        subroutine getegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,omg,wvn,
     2      atna,atnb)
        implicit none
        COMPLEX*16 A,B,C,D,E,F
        COMPLEX*16 wvno2, omega2
        COMPLEX*16 rsh, rp, rsv
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
        COMPLEX*16 omg, wvn
        complex*16 atna, atnb

        COMPLEX*16 xka2, xkb2

        integer NL
        parameter (NL=200)
        common/timodel/H(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real H, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
c-----
c       internal variables
c-----
        COMPLEX*16 L2(2)
        COMPLEX*16 bb, cc
        COMPLEX*16 CDSQRT 
        COMPLEX*16 SRTBB4AC
        COMPLEX*16 ddef, aabc

        COMPLEX*16 ZFAC

        real mu
        complex*16 CTF
        complex*16 atna2, atnb2

c-----
c       first test to see if a fluid layer - if it is fluid, the
c       eigenfunctions are specially computed and we need only the
c       rp
c-----
        atna2= atna*atna
        atnb2= atnb*atnb
        if(TL(m).eq.0.0 .or. TN(m).eq.0.0)then
            rp = cdsqrt(wvno2 -omega2*TRho(m)/(TA(m)*atna2))
            rsv = dcmplx(0.0d+000, 0.0d+00)
            rsh = dcmplx(0.0d+000, 0.0d+00)
            return
        endif
        
c-----
c       Do the SH
c-----
        rsh = CDSQRT(TN(m)*atnb2*wvno2/(TL(m)*atnb2) 
     1      - Trho(m)*omega2/(TL(m)*atnb2))
c-----
c       I use the artifice of a mu instead of the direct TC
c       since I wish to associate the attenuation with the
c       P and S directly
c-----
C        mu=(TC(m)+TA(m)+6.0*TL(m)
C     1      +5.0*TN(m)-2.0*TF(m))/15.0
C        CTF=0.5d+00*
C     1    ((TC(m)+TA(m))*atna*atna
C     2      +(6.0*TL(m)+5.0*TN(m)-15.0*mu)*atnb*atnb)

        CTF = TF(m)*(atna2 - 2.*(TL(m)/TA(m))*atnb2)/
     1         (1. - 2.*(TL(m)/TA(m)))

        a = wvn * CTF / (TC(m)*atna2)
        b = 1.0/(TC(m)*atna2)
        c = - TRho(m)*omg*omg + wvn*wvn *
     1      (TA(m)*atna2 -CTF*CTF/(TC(m)*atna2))
        d = - wvn
        e = 1.0/(TL(m)*atnb2)
        f = - TRho(m)*omg*omg

c-----
c       do algebra first to avoid numerical problems
c-----
        ddef = wvn*wvn - TRho(m)*omg*omg/(TL(m)*atnb2)
        aabc = wvn*wvn*TA(m)/TC(m) - TRho(m)*omg*omg/(TC(m)*atna2)

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
c       get the norms - note that the true norm will be 2  NP and 
c           2 L(2) NSV
c       The factorization permits us to use the sin nz/n or n sin nz
c-----
        NP  = (  L2(1)*(-2*a*b*d + 2*a*a*e + b*c*e - b*b*f)
     1      + (a*a+b*c)*(2*b*d*d - 2*a*d*e + b*e*f - c*e*e) )
        NSV = (- L2(2)*(2*b*d*d - 2*a*d*e - c*e*e + b*e*f)
     1      + (d*d+e*f)*(2*a*b*d - 2*a*a*e + b*b*f - b*c*e) )
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
c       Force the eigenfunctions to be as given in 5.4.4
c-----
        zfac = rp / x21
        x11  = x11 *zfac
        x21  = x21 *zfac
        x31  = x31 *zfac
        x41  = x41 *zfac

        zfac = rsv / x12
        x12  = rsv
        x22  = x22 * zfac
        x32  = x32 * zfac
        x42  = x42 * zfac
        
        np   = x11*x41 - x21*x31
        nsv  = x12*x42 - x22*x32

        return
        end
        subroutine hska(AA,tcosp,trsinp,tsinpr,tcossv,trsinsv,tsinsvr,
     1      NP,NSV,
     1      X11, X21, X31, X41,X12, X22, X32, X42, 
     2      rho,iwat,ex,om2)
c-----
c       Changes
c
c       01 May 2002  - defined cosp = tcosp/NP to reduce 
c           number  of complex divisions
c-----
        implicit none
        complex*16 AA(4,4)
        COMPLEX*16 tcosp , tcossv 
        COMPLEX*16 trsinp, trsinsv
        COMPLEX*16 tsinpr, tsinsvr
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        real rho
        integer iwat
        real*8 ex, dfac
        complex*16  om2
        complex*16 zrho

c-----
c       introduce shorthand to reduce work
c-----
        COMPLEX*16 cosp , cossv 
        COMPLEX*16 rsinp, rsinsv
        COMPLEX*16 sinpr, sinsvr

        integer i, j
        zrho = dcmplx(dble(rho),0.0d+00)
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,4
                do 101 i=1,4
                    AA(i,j) = dcmplx(0.0d+00,0.0d+00)
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

            AA(1,1) =   x11*x41*cosp  + x12*x42*cossv 
            AA(1,2) = - x11*x31*rsinp - x12*x32*sinsvr
            AA(1,3) = - x11*x21*cosp  - x12*x22*cossv 
            AA(1,4) =   x11*x11*rsinp + x12*x12*sinsvr
            AA(2,1) =   x21*x41*sinpr + x22*x42*rsinsv
            AA(2,2) = - x21*x31*cosp  - x22*x32*cossv 
            AA(2,3) = - x21*x21*sinpr - x22*x22*rsinsv
            AA(3,1) =   x31*x41*cosp  + x32*x42*cossv 
            AA(3,2) = - x31*x31*rsinp - x32*x32*sinsvr
            AA(4,1) =   x41*x41*sinpr + x42*x42*rsinsv
            AA(2,4) = - AA(1,3)
            AA(3,3) =   AA(2,2)
            AA(3,4) = - AA(1,2)
            AA(4,2) = - AA(3,1)
            AA(4,3) = - AA(2,1)
            AA(4,4) =   AA(1,1)
        endif
        return
        end

        subroutine hskl(cossh,rsinsh,sinshr,TL,iwat,hl,atnb)
        implicit none
        complex*16 cossh, rsinsh, sinshr
        real TL
        integer iwat
        complex*16 hl(2,2)
        complex*16 atnb
        if(iwat.eq.0)then   
            hl(1,1) = cossh
            hl(2,1) = TL*rsinsh*atnb*atnb
            hl(1,2) = sinshr/(TL*atnb*atnb)
            hl(2,2) = cossh
        else
            hl(1,1) = dcmplx(1.0d+00,0.0d+00)
            hl(1,2) = dcmplx(0.0d+00,0.0d+00)
            hl(2,1) = dcmplx(0.0d+00,0.0d+00)
            hl(2,2) = dcmplx(1.0d+00,0.0d+00)
        endif
        return
        end 

        subroutine hsupdn(aa,saa,wvno2,om,om2,wvno,m,
     1      lmaxr,lmaxs,isr)
c-----
c       aa  C*16    P-SV matrix
c       saa C*16    SH matrix
c       wvno2   C*16    wavenumber squared
c       om  C*16    complex angular frequency
c       wvno    C*16    wavenumber 
c       m   I   layer index for matrix computation
c               (perform wavefield separation in this layer)
c       lmaxr   I   index of receiver layer
c       lmaxs   I   index of source layer
c               (used for definition of up/down)
c       isr I   0 - compete for receiver layer
c               1 - compute for source
c-----
        implicit none
        complex*16 aa(4,4), saa(2,2), wvno2, om, wvno,om2
        integer m, lmaxr, lmaxs, isr
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud
        complex*16 zero, zone 
        logical pup, pdn, sup, sdn, doud
        complex*16 CDSQRT

        complex*16 atna, atnb
        complex*16 atna2, atnb2
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv, rsh
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42

        integer iwat, i, j
c-----
c       define necessary constants
c-----
        zero = dcmplx(0.0d+00,0.0d+00)
        zone = dcmplx(1.0d+00,0.0d+00)
c-----
c       if receiver is beneath source, it is necessary to redefine
c       local concept of up/down
c-----
        if(isr.eq.0)then
            if(lmaxs .gt. lmaxr)then
                pup = rpup
                pdn = rpdn
                sup = rsup
                sdn = rsdn
            else
                pup = rpdn
                pdn = rpup
                sup = rsdn
                sdn = rsup
            endif
            doud = dorud
        else
            if(lmaxs .gt. lmaxr)then
                pup = spup
                pdn = spdn
                sup = ssup
                sdn = ssdn
            else
                pup = spdn
                pdn = spup
                sup = ssdn
                sdn = ssup
            endif
            doud = dosud
        endif
        call aten(om,qa(m),qb(m),alpha,atna,atnb,frefp(m),frefs(m))
        atna2 = atna*atna
        atnb2 = atnb*atnb
            if(TL(m).eq. 0.0 .or.TN(m).eq.0.0)then
                iwat = 1
            else
                iwat = 0
            endif
        call getegn(za,zb,zc,zd,ze,zf,om2,wvno2,rp, rsv, rsh,
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,om,wvno,
     2      atna,atnb)
        do 100 j=1,4
            do 101 i = 1,4
                aa(i,j) = zero
  101       continue
  100   continue
        saa(1,1) = zero
        saa(1,2) = zero
        saa(2,1) = zero
        saa(2,2) = zero
c-----
c       now test wave field conditions
c-----
        if(pup .and. pdn .and. sup .and. sdn)then
            aa(1,1)  = zone
            aa(2,2)  = zone
            aa(3,3)  = zone
            aa(4,4)  = zone
            saa(1,1) = zone
            saa(2,2) = zone
        else
c-----
c           water layer
c-----
            if(iwat.ne.0)then
c-----
c               coefficients of exp(nua*h)
c-----
                if(pup)then
                    aa(2,2) = aa(2,2) + 0.5d+00*zone
                    aa(2,3) = aa(2,3) - 
     1                  rp/(2.0d+00*Trho(m)*om2)
                    aa(3,2) = aa(3,2) - 
     1                  Trho(m)*om2/(2.0d+00*rp)
                    aa(3,3) = aa(3,3) + 0.5d+00*zone
                endif
c-----
c               coefficients of exp(-nua*h)
c-----
                if(pdn)then
                    aa(2,2) = aa(2,2) + 0.5d+00*zone
                    aa(2,3) = aa(2,3) + 
     1                  rp/(2.0d+00*Trho(m)*om2)
                    aa(3,2) = aa(3,2) + 
     1                  Trho(m)*om2/(2.0d+00*rp)
                    aa(3,3) = aa(3,3) + 0.5d+00*zone
                endif

c-----
c           elastic layer
c-----
            else
c-----
c               coefficients of exp(nua*h)
c-----
                if(pup)then
                    AA(1,1) = AA(1,1) + x11*x41/NP 
                    AA(1,2) = AA(1,2) - x11*x31*rp/NP 
                    AA(1,3) = AA(1,3) - x11*x21/NP 
                    AA(1,4) = AA(1,4) + x11*x11*rp/NP 
                    AA(2,1) = AA(2,1) + x21*x41/(NP*rp )
                    AA(2,2) = AA(2,2) - x21*x31/NP 
                    AA(2,3) = AA(2,3) - x21*x21/(NP*rp )
                    AA(3,1) = AA(3,1) + x31*x41/NP 
                    AA(3,2) = AA(3,2) - x31*x31*rp/NP 
                    AA(4,1) = AA(4,1) + x41*x41/(NP*rp )

                endif
c-----
c               coefficients of exp(-nua*h)
c-----
                if(pdn)then
                    AA(1,1) = AA(1,1) + x11*x41/NP 
                    AA(1,2) = AA(1,2) + x11*x31*rp/NP 
                    AA(1,3) = AA(1,3) - x11*x21/NP 
                    AA(1,4) = AA(1,4) - x11*x11*rp/NP 
                    AA(2,1) = AA(2,1) - x21*x41/(NP*rp )
                    AA(2,2) = AA(2,2) - x21*x31/NP 
                    AA(2,3) = AA(2,3) + x21*x21/(NP*rp )
                    AA(3,1) = AA(3,1) + x31*x41/NP 
                    AA(3,2) = AA(3,2) + x31*x31*rp/NP 
                    AA(4,1) = AA(4,1) - x41*x41/(NP*rp )
                endif
c-----
c               coefficients of exp(nub*h)
c-----
                if(sup)then
                    AA(1,1) = AA(1,1) + x12*x42/NSV
                    AA(1,2) = AA(1,2) - x12*x32/(NSV*rsv)
                    AA(1,3) = AA(1,3) - x12*x22/NSV
                    AA(1,4) = AA(1,4) + x12*x12/(NSV*rsv)
                    AA(2,1) = AA(2,1) + x22*x42*rsv/NSV
                    AA(2,2) = AA(2,2) - x22*x32/NSV
                    AA(2,3) = AA(2,3) - x22*x22*rsv/NSV
                    AA(3,1) = AA(3,1) + x32*x42/NSV
                    AA(3,2) = AA(3,2) - x32*x32/(NSV*rsv)
                    AA(4,1) = AA(4,1) + x42*x42*rsv/NSV

                    saa(1,1) = saa(1,1) + zone
                    saa(1,2) = saa(1,2) 
     1                  + zone/(TL(m)*rsh*atnb2)
                    saa(2,1) = saa(2,1) + TL(m)*rsh*atnb2
                    saa(2,2) = saa(2,2) + zone
                endif
c-----
c               coefficients of exp(-nub*h)
c-----
                if(sdn)then
                    AA(1,1) = AA(1,1) + x12*x42/NSV
                    AA(1,2) = AA(1,2) + x12*x32/(NSV*rsv)
                    AA(1,3) = AA(1,3) - x12*x22/NSV
                    AA(1,4) = AA(1,4) - x12*x12/(NSV*rsv)
                    AA(2,1) = AA(2,1) - x22*x42*rsv/NSV
                    AA(2,2) = AA(2,2) - x22*x32/NSV
                    AA(2,3) = AA(2,3) + x22*x22*rsv/NSV
                    AA(3,1) = AA(3,1) + x32*x42/NSV
                    AA(3,2) = AA(3,2) + x32*x32/(NSV*rsv)
                    AA(4,1) = AA(4,1) - x42*x42*rsv/NSV

                    saa(1,1) = saa(1,1) + zone
                    saa(1,2) = saa(1,2) 
     1                  - zone/(TL(m)*rsh*atnb2)
                    saa(2,1) = saa(2,1) - TL(m)*rsh*atnb2
                    saa(2,2) = saa(2,2) + zone
                endif
                aa(2,4) = -aa(1,3)
                aa(3,3) =  aa(2,2)
                aa(3,4) = -aa(1,2)
                aa(4,2) = -aa(3,1)
                aa(4,3) = -aa(2,1)
                aa(4,4) =  aa(1,1)
c-----
c       special case clean up
c-----
                if(pup.and.sup .or. pdn.and.sdn)then
                    aa(1,1) = zone
                    aa(2,2) = zone
                    aa(3,3) = zone
                    aa(4,4) = zone
                    aa(1,3) = zero
                    aa(2,4) = zero
                    aa(3,1) = zero
                    aa(4,2) = zero
                endif
                do 200 j=1,4
                    do 201 i = 1,4
                        aa(i,j) = aa(i,j)/dble(2.0)
  201               continue
  200           continue
                saa(1,1) = saa(1,1) / dble(2.0)
                saa(1,2) = saa(1,2) / dble(2.0)
                saa(2,1) = saa(2,1) / dble(2.0)
                saa(2,2) = saa(2,2) / dble(2.0)
            endif
        endif
        return
        end

        subroutine insert(dph)
        implicit none
        real dph
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax

        integer m
        real dep, dp, dphh, hsave
        integer ls
c-----
c       Insert a depth point into the model by splitting a layer
c       so that the point appears at the top boundary of the layer
c       dph = depth of interest
c-----
c       determine the layer in which the depth dph lies.
c       if necessary, adjust  layer thickness at the base
c-----
c       Here determine layer corresponding to specific depth dph
c       If the bottom layer is not thick enough, extend it
c
c       dep - depth to bottom of layer
c       dphh    - height of specific depth above bottom of the layer
c-----
        if(dph.le.0)then
            d(1) = d(1) - dph
            return
        else if(dph.ge.0)then
            dep = 0.0 
            dp = 0.0 
            dphh = -1.0
            do 100 m = 1,mmax 
                dp = dp + d(m) 
                dphh = dp - dph 
                if(m.eq.mmax)then
                    if(d(mmax).le.0.0 .or. dphh.lt.0.0)then
                        d(mmax) = (dph - dp)
                    endif
                endif
                dep = dep + d(m) 
                dphh = dep - dph 
                ls = m 
                if(dphh.ge.0.0) go to 101 
  100       continue 
  101       continue 
        endif
c-----
c       In the current model, the depth point is in the ls layer
c       with a distance dphh to the bottom of the layer
c
c       Do not create unnecessary layers, e.g., 
c           at surface and internally
c       However do put in a zero thickness layer 
c           at the base if necessary
c-----
        if(dphh .eq. 0.0 .and. ls.ne.mmax)then
            return
        else
c-----
c           adjust layering
c-----
             do 102 m = mmax,ls,-1
                d(m+1) = d(m)
                TA(m+1) = TA(m)
                TC(m+1) = TC(m)
                TF(m+1) = TF(m)
                TL(m+1) = TL(m)
                TN(m+1) = TN(m)
                TRho(m+1) = TRho(m)
                qa(m+1) = qa(m)
                qb(m+1) = qb(m)
                etap(m+1) = etap(m)
                etas(m+1) = etas(m)
                frefp(m+1) = frefp(m)
                frefs(m+1) = frefs(m)
  102       continue
            hsave = d(ls)
            d(ls) = hsave - dphh
            d(ls+1) = dphh
            ls = ls + 1
            mmax = mmax + 1
            if(d(mmax).lt.0.0)d(mmax)=0.0
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

        subroutine modcpy(totmp) 
        implicit none
        logical totmp
c-----
c       copy model to temporary array
c-----
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/modelt/tD(NL),tTA(NL),tTC(NL),tTL(NL),tTN(NL),tTF(NL),
     1      tTRho(NL),
     2      tqa(NL),tqb(NL),tetap(NL),tetas(NL), 
     3      tfrefp(NL), tfrefs(NL)
        real tD, tTA, tTC, tTF, tTL, tTN, tTRho
        real tqa, tqb, tetap, tetas, tfrefp, tfrefs
        integer i, mmaxt
c-----
c       copy to temporary array
c-----
        if(totmp)then
            do 20 i = 1,mmax 
                td(i) = d(i)
                tTA(i) = TA(i)
                tTC(i) = TC(i)
                tTF(i) = TF(i)
                tTL(i) = TL(i)
                tTN(i) = TN(i)
                tTRho(i) = TRho(i)
                tqa(i) = qa(i)
                tqb(i) = qb(i)
                tetap(i) = etap(i)
                tetas(i) = etas(i)
                tfrefp(i) = frefp(i)
                tfrefs(i) = frefs(i)
   20       continue 
            mmaxt = mmax
        else
            do 30 i = 1,mmaxt 
                d(i) = td(i)
                TA(i) = tTA(i)
                TC(i) = tTC(i)
                TF(i) = tTF(i)
                TL(i) = tTL(i)
                TN(i) = tTN(i)
                TRho(i)  = tTRho(i)
                qa(i)    = tqa(i)
                qb(i)    = tqb(i)
                etap(i)  = tetap(i)
                etas(i)  = tetas(i)
                frefp(i) = tfrefp(i)
                frefs(i) = tfrefs(i)
   30       continue 
            mmax = mmaxt
        endif
        return 
        end 

        subroutine normc(e,ex,xnorm)
        implicit none

        complex*16 e(*)
        real*8 ex, xnorm

        common/lwater/lfluid
        logical lfluid
        real*8 test,testt,x,y,fac
        real*8 DREAL
        integer i, IUP
        test = 0.0D+00
        testt = 0.0D+00
        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 2 i = 1,IUP
            if(dabs(dreal(e(i))).gt.testt) testt =dabs(dreal(e(i)))
            if(dabs(dimag(e(i))).gt.testt) testt =dabs(dimag(e(i)))
    2   continue
        if(testt.lt.1.0e-30)testt=1.0
        do 1 i =1,IUP
            x=dreal(e(i))/testt
            y=dimag(e(i))/testt
            fac = x*x + y*y
            if(test.lt.fac) test = fac
    1   continue
        test = testt*dsqrt(test)
        if(test.lt.1.0d-30) test=1.0
        xnorm = 1./test
        ex =-dlog(xnorm)
        return
        end

        subroutine rcmult(y,c,exa,exe)
        implicit none
        complex*16 y(5,5)
        complex*16 c(5,5)
        real*8 exa, exe
        common/lwater/lfluid
        logical lfluid
c-----
c       FORM YC where y(5x5) c(5x5) RETURN Y
c-----
        real*8 eval
        real*8 xnorm
        complex*16 ztmp, ee(5,5)
        integer IUP, i, j, k
        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 1350 i=1,IUP
            do 1351 j=1,IUP
                ztmp = dcmplx(0.0d+00,0.0d+00)
                do 1349 k=1,IUP
                    ztmp=ztmp+ y(i,k) * c(k,j)
 1349           continue
                ee(i,j)=ztmp
 1351       continue
 1350   continue
        exe = exe + exa
        call rnormc(ee,eval,xnorm)
        do 1353 j=1,IUP
            do 1352 i=1,IUP
                y(i,j) = ee(i,j)*xnorm
 1352       continue
 1353   continue
        exe = exe + eval
        return
        end

        subroutine rnormc(e,ex,xnorm)
        implicit none

        common/lwater/lfluid
        logical lfluid
        real*8 ex
        real*8 DREAL
        real *8 test,testt,x,y,fac,xnorm
        complex*16 e(5,5)

        integer i, j, IUP

        test = 0.0D+00
        testt = 0.0D+00
        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 3 j=1,IUP
            do 2 i = 1,IUP
            if(dabs(dreal(e(i,j))).gt.testt) testt =dabs(dreal(e(i,j)))
            if(dabs(dimag(e(i,j))).gt.testt) testt =dabs(dimag(e(i,j)))
    2       continue
    3   continue
        if(testt.lt.1.0e-30)testt=1.0
        do 4 j=1,IUP
            do 1 i =1,IUP
                x=dreal(e(i,j))/testt
                y=dimag(e(i,j))/testt
                fac = x*x + y*y
                if(test.lt.fac) test = fac
    1       continue
    4   continue
        test = testt*dsqrt(test)
        if(test.lt.1.0d-30) test=1.0
        xnorm = 1./test
        ex =-dlog(xnorm)
        return
        end

        subroutine rshof(gg,om,wvno, lmaxs, lmaxr, wvno2, om2) 
        implicit none
        complex gg(21) 
        complex*16 om, wvno, om2, wvno2
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        real vmin,vamin,vamax,vbmin,vbmax
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
c-----
c       gus - surface displacements or potentials or top of layer
c-----
        complex*16 gus(21)
        complex*16 cd(5),da(4,4),fr,y(4,4)
        complex*16 fourpo
        complex*16 d11,d12,fl 
        complex*16 s21,s32,s14,s34,s32e,s34e 
        complex*16 s24,s33
        complex*16 wv4pi
        complex*16 zero
        real *8 fact,exe,exl,exel,exll,elj
        real *8 exwu
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud
        complex*16 sdd(4), sds(4), sss(4), sep(4), svf(4), shf(4)
        complex*16 zone
        common/jout/jsrc(21), jbdrys, jbdryh
        integer jsrc, jbdrys, jbdryh
        complex*16 haa(4,4), saa(2,2)
        real*8 fourpi
        real*8 CDABS
        complex*16 CDSQRT
        real*8 DREAL

        complex*16 atna, atnb
        complex*16 atna2, atnb2
        real mu
        complex*16 CTF

        integer i, j, k, iwats, lmaxr, lmaxs, jj, iwatr



c-----
c       Initialization
c-----
        fourpi=12.5663706d+00
        fourpo=12.5663706d+00*om*om
        wv4pi = 2.0d+00 * wvno / fourpi
        zero  = dcmplx(0.0d+00,0.0d+00)
        zone  = dcmplx(1.0d+00,0.0d+00)
c-----
c       do not evaluate for wvno = 0.0
c-----
        if(CDABS(wvno).eq.0.0d+00) then
            do 102 i=1,21
                gg(i) = cmplx(0.0,0.0)
  102       continue
        else
c-----
c       process for this wavenumber and frequency
c-----
            do 101 i = 1,21
                gus(i) = zero
  101       continue
            call aten(om,qa(lmaxs),qb(lmaxs),alpha,atna,atnb,
     1          frefp(lmaxs),frefs(lmaxs))
C            mu=(TC(lmaxs)+TA(lmaxs)+6.0*TL(lmaxs)
C     1          +5.0*TN(lmaxs)-2.0*TF(lmaxs))/15.0
C            CTF=0.5d+00*
C     1      (  (TC(lmaxs)+TA(lmaxs))*atna*atna
C     2        +(6.0*TL(lmaxs)+5.0*TN(lmaxs)-15.0*mu)*atnb*atnb)
            atna2 = atna * atna
            atnb2 = atnb * atnb
        CTF = TF(lmaxs)*(atna2 - 2.*(TL(lmaxs)/TA(lmaxs))*atnb2)/
     1         (1. - 2.*(TL(lmaxs)/TA(lmaxs)))

            if(TL(lmaxs).eq. 0.0 .or.TN(lmaxs).eq.0.0)then
                iwats = 1
            else
                iwats = 0
            endif
            if(TL(lmaxr).eq. 0.0 .or.TN(lmaxr).eq.0.0)then
                iwatr = 1
            else
                iwatr = 0
            endif
            if(dosud)then
                call hsupdn(haa,saa,wvno2,om,om2,wvno,lmaxs,
     1              lmaxr,lmaxs,1)
            endif
c-----
c           this follows Wang and Herrmann (1980), 
c               "A numerical study of P-, SV-, and SH-wave 
c               generation in a plane layerd medium,
c               Bull. Seism. Soc. Am. 70, 1015-1036.
c           For p-SV we use equation (8). Note below in "do 60" 
c               that the Z component (odd index) is reversed 
c               from (8) so that +z is up
c
c       using the trick of reducing the 6x6 compound matrices to 5x5, 
c           the
c       following correspondence between what is given her and that in 
c           Wang
c       and Herrmann is in effect
c-----
c       Compound Matrix (6x6)   Here (5x5)
c
c       X|12/12 ==      cd(1)
c       X|12/13 ==      cd(2)
c       X|12/14 ==      cd(3)
c       X|12/23 ==      -cd(3)
c       X|12/24 ==      cd(4)
c       X|12/34 ==      cd(5)
c       X|12/21 = -X|12/12
c       X|12/ij = -X|12/ji
c
c       where
c
c        |12
c       X|  = X|12/ij
c        |ij
c-----
            call scoef(cd,da,fr,om,exe,exl,exwu,wvno,
     1          fl,d11,d12,exel,exll,lmaxs,lmaxr,wvno2,om2) 
c-----
c       Form X|12/ij x Zj2 in Equation 8 of Wang and Herrmann
c----
c KLUDGE to CHANGE ORDER AND ALSO GET UR CORRECT FOR WATER LAYER
            do 50 k=1,2
c-----
c       k =1 UR     k=2 UZ
c       jj=1 = UZ, jj=2 = UR
c-----
                if(k.eq.1)then
                    jj = 2
                else if(k.eq.2)then
                    jj = 1
                endif
                if(iwatr.eq.1 .and. jj.eq.2)then
c-----
c       compute Tz and then get Ur from UZ
c-----
c-----
                    j = 3
                else
                    j = k
                endif
                y(1,jj)= cd(1)*da(2,j) + cd(2)*da(3,j) 
     1              + cd(3)*da(4,j)
                y(2,jj)= -cd(1)*da(1,j) - cd(3)*da(3,j) 
     1              + cd(4)*da(4,j)
                y(3,jj)=-cd(2)*da(1,j) + cd(3)*da(2,j) 
     1              + cd(5)*da(4,j)
                y(4,jj)= -cd(3)*da(1,j) - cd(4)*da(2,j) 
     1              - cd(5)*da(3,j)
   50       continue
            if(iwatr.eq.1)then
                y(1,2) = - wvno*y(1,2)/(TRho(lmaxr)*om*om)
                y(2,2) = - wvno*y(2,2)/(TRho(lmaxr)*om*om)
                y(3,2) = - wvno*y(3,2)/(TRho(lmaxr)*om*om)
                y(4,2) = - wvno*y(4,2)/(TRho(lmaxr)*om*om)
            endif
c-----
c           evaluate different Green's functions
c           apply source terms
c----- 
c-----
c           START OF P-SV
c-----
c           First compute the DELTA displacement-stress source terms
c           for inverted model, the UZ, TR elements change, These will
c           be only those required for dipoles and forces
c           
c           Stress-displacement discontinuities for Green's functions
c           Green   dUr dUz dTz dTr
c           DD      s32     s34
c           DS  s21
c           SS              s14
c           EX      s32e        s34e
c           VF          s33
c           HF              s24
c-----
            if(iwats.eq.1)then
                s14  = zero
                s21  = zero
                s24  = zero
                s32  = zero
                s33  = zero
                s34  = zero
                s34e = zero
            else
                s14  = -wv4pi
                s21  =  2.0d+00/(fourpi*TL(lmaxs)*atnb*atnb)
                s24  = -2.0d+00/fourpi
                s32  =  4.0d+00/(fourpi*TC(lmaxs)*atna*atna)
                s33  = -2.0d+00/fourpi
                s34  = -wv4pi*(1.0d+00+2.0d+00*CTF/
     1              (TC(lmaxs)*atna*atna))
                s34e =  wv4pi*(1.0d+00-CTF/(TC(lmaxs)*atna*atna))
            endif
            s32e=2.0d+00/(fourpi*TC(lmaxs)*atna*atna)
c-----
c           receiver beneath the source
c-----
            if(lmaxr .gt. lmaxs)then
                s14  = - s14
                s24  = - s24
                s32  = - s32
                s32e = - s32e
                s34  = - s34
                s34e = - s34e
            endif
c-----
c           For complete wavefield computation do not
c           waste cycles computing W matrix elements
c-----
            if(.not. dosud)then
                do 61 j=1,2
c       DD
                    gus(j   )=s32 *y(2,j)+ s34*y(4,j)
c       DS
                    gus(j+ 2)=s21 *y(1,j)             
c       SS
                    gus(j+ 4)=             s14*y(4,j)
c       EX
                    gus(j+ 6)=s32e*y(2,j)+s34e*y(4,j)
c       VF
                    gus(j+ 8)=s33 *y(3,j)
c       HF
                    gus(j+10)=s24 *y(4,j)
   61           continue
            ELSE
                sdd(1) = zero
                sdd(2) = s32
                sdd(3) = zero
                sdd(4) = s34
                sds(1) = s21
                sds(2) = zero
                sds(3) = zero
                sds(4) = zero
                sss(1) = zero
                sss(2) = zero
                sss(3) = zero
                sss(4) = s14
                sep(1) = zero
                sep(2) = s32e
                sep(3) = zero
                sep(4) = s34e
                svf(1) = zero
                svf(2) = zero
                svf(3) = s33
                svf(4) = zero
                shf(1) = zero
                shf(2) = zero
                shf(3) = zero
                shf(4) = s24
c-----
c               Change the source vector for up/down going wavefields
c-----

                call svupdn(sdd,haa)
                call svupdn(sds,haa)
                call svupdn(sss,haa)
                call svupdn(sep,haa)
                call svupdn(svf,haa)
                call svupdn(shf,haa)
c-----
c           compute the response
c-----
                do 64 j=1,2
                    gus(j   ) = sdd(1)*y(1,j) + sdd(2)*y(2,j) 
     1                  + sdd(3)*y(3,j) + sdd(4)*y(4,j)
                    gus(j+ 2) = sds(1)*y(1,j) + sds(2)*y(2,j) 
     1                  + sds(3)*y(3,j) + sds(4)*y(4,j)
                    gus(j+ 4) = sss(1)*y(1,j) + sss(2)*y(2,j) 
     1                  + sss(3)*y(3,j) + sss(4)*y(4,j)
                    gus(j+ 6) = sep(1)*y(1,j) + sep(2)*y(2,j) 
     1                  + sep(3)*y(3,j) + sep(4)*y(4,j)
                    gus(j+ 8) = svf(1)*y(1,j) + svf(2)*y(2,j) 
     1                  + svf(3)*y(3,j) + svf(4)*y(4,j)
                    gus(j+10) = shf(1)*y(1,j) + shf(2)*y(2,j) 
     1                  + shf(3)*y(3,j) + shf(4)*y(4,j)
   64           continue
            ENDIF
c-----
c       invert the vertical
c-----
            do 62 j=1,12,1
                gus(j) = -gus(j)
   62       continue
c-----
c           if receiver beneath the source unflip radial
c-----
            if(lmaxr .gt. lmaxs)then
                do 63 j=2,12,2
                    gus(j) = - gus(j)
   63           continue
            endif
c-----
c           If the receiver is in the water and the source 
c           is an explosion generate pressure time history at receiver
c           Also in water layer radial time history 
c           is generated differently
c-----
            if(iwatr.eq.1)then
                gus(16) =  dble(TRho(lmaxr))*gus(8) /wvno
                gus(17) =  dble(TRho(lmaxr))*gus(2) /wvno
                gus(18) =  dble(TRho(lmaxr))*gus(4) /wvno
                gus(19) =  dble(TRho(lmaxr))*gus(6) /wvno
                gus(20) =  dble(TRho(lmaxr))*gus(10) /wvno
                gus(21) =  dble(TRho(lmaxr))*gus(12) /wvno
            endif
c-----
c           END OF P-SV
c-----
c           START OF SH
c-----
            if(iwats.eq.0 .and. iwatr.eq.0)then
                sds(1) = 2.0d+00/(fourpi*TL(lmaxs)*atnb*atnb)
                sds(2) = zero
                sss(1) = zero
                sss(2) = -2.0d+00*wvno/12.5663706d+00
                shf(1) = zero
                shf(2) = -2.0d+00/12.5663706d+00
                if(.not. dosud )then
                    gus(13) = ( d11*sds(1)           )
                    gus(14) = (            d12*sss(2))
                    gus(15) = (            d12*shf(2))
                ELSE
                    call shupdn(sds,saa)
                    call shupdn(sss,saa)
                    call shupdn(shf,saa)
                    gus(13) = ( d11*sds(1) + d12*sds(2))
                    gus(14) = ( d11*sss(1) + d12*sss(2))
                    gus(15) = ( d11*shf(1) + d12*shf(2))
                ENDIF
                if(lmaxr .gt. lmaxs)then
                    gus(13) = - gus(13)
                endif
            endif
c-----
c           END OF SH
c-----
c-----
c           do final scaling for exponential
c-----
c-----
c           SV
c-----
c           fix for radial derived from vertical for fluid
c-----
            do 71 k=1,2
                elj = -exe + exl 
                fact = 0.0D+00
                if(elj.gt.-55.) fact=dexp(elj)
                do 72 i=0,10,2
                    j = i + k
                    gg(j) = ( gus(j) * fact/fr)
c-----
c           flip UZ to make vertical positive up
c-----
                    if(k.eq.1)then
                        gg(j) = -gg(j)
                    endif
   72           continue
c----
c           do pressure field
c----
                if(k.eq.1)then
                    gg(16) = - (gus(16)*fact/fr)
                    gg(17) = - (gus(17)*fact/fr)
                    gg(18) = - (gus(18)*fact/fr)
                    gg(19) = - (gus(19)*fact/fr)
                    gg(20) = - (gus(20)*fact/fr)
                    gg(21) = - (gus(21)*fact/fr)
                endif
   71       continue
c-----
c           SH
c-----
            elj=-exel+exll
            if(iwats.eq.0)then
                if(elj.gt.-55.) then
                    fact = dexp(elj)
                    gg(13)=(gus(13)*fact)/(fl)
                    gg(14)=(gus(14)*fact)/(fl)
                    gg(15)=(gus(15)*fact)/(fl)
                else
                    gg(13) = cmplx(0.0,0.0)
                    gg(14) = cmplx(0.0,0.0)
                    gg(15) = cmplx(0.0,0.0)
                endif
            else
                gg(13) = cmplx(0.0,0.0)
                gg(14) = cmplx(0.0,0.0)
                gg(15) = cmplx(0.0,0.0)
            endif
        endif
        return
        end

        subroutine scoef(cd,da,fr,om,exe,exl,exwu,wvno,
     1      fl,d11,d12,exel,exll,llmaxs,llmaxr,wvno2, om2)
c-----
c       llmaxs  I*4 - source layer index in original model
c       llmaxr  I*4 - receiver layer index in original model
c-----
        implicit none
        complex*16 cd(5), da(4,4), fr, om,wvno, om2, wvno2
        real*8 exe, exl, exwu
        complex*16 fl, d11, d12
        real*8 exel, exll
        integer llmaxs, llmaxr
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/damp/alpha,ieqex
        real alpha  
        integer ieqex
        complex*16  ca(5,5)
        complex*16  e(5)
        complex*16 e1,e2,e21, e22
        real *8 ex,exa,exb
        real*8 dzero
        complex*16 zdum
        complex*16 aa(4,4)
        complex*16 cy(5,5)
        complex*16 zero, zone
        complex*16 y11, y12, y21, y22, sd11, sd21
        complex*16 hl(2,2)
c-----
c       matrix components in layers and boundaries saved
c-----
        common/hamat/har
        complex*16 har(NL,4,4)
        common/damat/dar
        complex*16 dar(NL,5,5)
        common/hsrfr/hsr
        complex*16 hsr(2,5)
        common/gbrfr/gbr
        complex*16 gbr(2,5)
        common/hlmat/hal
        complex*16 hal(NL,2,2)
        common/hsrfl/hsl
        complex*16 hsl(2,2)
        common/gbrfl/gbl
        complex*16 gbl(2,2)
        common/hexex/hex
        real*8 hex(NL)
        common/hexexw/hexw
        real*8 hexw(NL)
        common/dexex/dex
        real*8 dex(NL)
        common/lexex/lex 
        real*8 lex(NL)
        common/water/iwater(NL),iwats(2),iwatb(2)
        integer iwater, iwats, iwatb
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
        logical retrieve
c-----
c       check for decomposition at wavefield at receiver
c-----
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud
        complex*16 haa(4,4), saa(2,2)

        integer i, j, mm, m
        integer iwat, lmaxr, lmaxs, in

c-----
c       this routine computes the layer response. 
c           To simplify the mathematics
c       of the case of receiver above or beneath the source, the
c       layer is internally flipped.
c
c       llmaxs  I*4 - source layer index in original model
c       llmaxr  I*4 - receiver layer index in original model
c       mm  I*4 - pointer to layer in original model
c
c       lmaxs   I*4 - source layer index in current model
c       lmaxr   I*4 - receiver layer index in current model
c       mm  I*4 - pointer to layer in current model
c
c       in  I*4 - 1 use current model (source beneath receiver)
c               - 2 use inverted model (receiver beneath source)
c-----
c       initialize matrices
c-----
        zero = dcmplx(0.0d+00,0.0d+00)
        zone  = dcmplx(1.0d+00,0.0d+00)
        dzero = 0.0d+00
        exe=0.0
        exl=0.0
        exwu = 0.0
        do 2 j = 1,4
            do 3 i = 1,4
                da(i,j)=zero
    3       continue
            da(j,j) = zone
    2   continue
        do 12 j=1,5
            do 13 i=1,5
                cy(i,j) = zero
   13       continue
            cy(j,j) = zone
   12   continue
        y11 = zone
        y12 = zero
        y21 = zero
        y22 = zone
        exel = 0.0
        exll = 0.0
c-----
c     set up halfspace conditions
c-----
        if(llmaxs .ge. llmaxr)then
            in = 1
        else
            in = 2
        endif
        do 100 i=1,5
            e(i) = gbr(in,i)
  100   continue
        e1 = gbl(in,1)
        e2 = gbl(in,2)
        do 11 i=1,5
            cd(i)=e(i)
   11   continue
        d11=e1
        d12=e2
c-----
c       set up limits on the layer stacking
c-----
        if(llmaxs .ge. llmaxr)then
            lmaxs = llmaxs
            lmaxr = llmaxr
        else
            lmaxs = mmax - llmaxs + 2
            lmaxr = mmax - llmaxr + 2
        endif
c-----
c       matrix multiplication from bottom layer upward
c-----
        do 1340 mm = mmax,1,-1
            if(llmaxs .ge. llmaxr)then
                m = mm
                if(equalu(m))then
                    retrieve = .false.
                else
                    retrieve = .true.
                endif
            else
                m = mmax + 1 - mm
                if(equald(m))then
                    retrieve = .false.
                else
                    retrieve = .true.
                endif
            endif
            iwat = iwater(m)
            if(retrieve)then
                call copy5(ca,dar,m,1,dex,exa)
                call copy2(hl,hal,m,1,lex,exb)
                call copy4(aa,har,m,1,hex,ex)
            endif
            call cmult(e,ca,exa,exe)
            call lmult(e1,e2,e21,e22,hl,iwat,exel,exb,.false.)
            if(mm.lt.lmaxr)then
                call rcmult(cy,ca,exa,exl)
                call lmult(y11,y12,y21,y22,hl,iwat,
     1              exll,exb,.true.)
            else if(mm.ge.lmaxr .and. mm.lt.lmaxs) then
                call dmult(da,aa)
                    exl = exl + ex
c-----
c       save values at top of source layer
c-----
            else if(mm.eq.lmaxs) then
                    do 1352 i=1,5
                    cd(i)=e(i)
 1352           continue
                exl=exe
                exll = exel
                d11=e1
                d12=e2
            endif
            if(mm.eq.1)then
                do 200 i=1,5
                    ca(i,1) = hsr(in,i)
  200           continue
                sd11 = hsl(in,1)
                sd21 = hsl(in,2)
                zdum = e1
                e1 = zdum*sd11 + e2*sd21
c               e2 = zdum*sd11 - e2*sd21
                zdum = y11
                y11 = zdum*sd11 + y12*sd21
                zdum = y21
                y21 = zdum*sd11 + y22*sd21
                zdum = dcmplx(0.0,0.0)
                do 1402 i=1,5
                    zdum = zdum + e(i)*ca(i,1)
 1402           continue
                e(1) = zdum
                call rcmult(cy,ca,dzero,exl)
            endif
 1340 continue
c-----
c       get final matrices
c-----
c-SH
        fl=e1
c-P-SV
c       form x(l,m)y(ij|12)
c-----take care of x(i,j) y(1j|12) and replace the da
            aa(1,1) =   zero
            aa(2,1) =   cy(1,1)
            aa(3,1) =   cy(2,1)
            aa(4,1) =   cy(3,1) /(2.0d+00 )
c change sign 0430 1200
            aa(1,2) = - aa(2,1)
            aa(2,2) =   zero
            aa(3,2) = - aa(4,1)
            aa(4,2) =   cy(4,1)
            aa(1,3) = - aa(3,1)
            aa(2,3) = - aa(3,2)
            aa(3,3) =   zero
            aa(4,3) =   cy(5,1)
            aa(1,4) = - aa(4,1)
            aa(2,4) = - aa(4,2)
            aa(3,4) = - aa(4,3)
            aa(4,4) =   zero
            call dmult(da,aa)
        fr=e(1)
c-----
c       if decomposion of wavefield at receiver modify the
c       returned matrices
c-----
        if(dorud)then
            
            call hsupdn(haa,saa,wvno2,om,om2,wvno,llmaxr, 
     1          llmaxr,llmaxs,0)
c-----
c           SH
c           B = SAA B
c-----
            d11 = d11*(saa(1,1)*y11 + saa(1,2)*y21)
            d12 = d12*(saa(1,1)*y11 + saa(1,2)*y21)
c-----
c           SV
c           BT = ST CD X HAAT
c           Form X HAAT
c-----
            call trans4(haa)
            call dmult(da,haa)
        else
            d11 = y11*d11
            d12 = y11*d12
        endif
        return
        end

        subroutine shupdn(s,saa)
        implicit none
        complex*16 s(2), saa(2,2)
        complex*16 tmp
            tmp  = saa(1,1)*s(1) + saa(1,2)*s(2)
            s(2) = saa(2,1)*s(1) + saa(2,2)*s(2)
            s(1) = tmp
        return
        end

        subroutine srclay(depth,lmax,dph)
        implicit none
        real depth, dph
        integer lmax
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/lyrctl/lyrins
        logical lyrins
        if(.not.lyrins)then
            call modcpy(.false.)
            call insert(depth)
        endif
        call srclyr(depth,lmax,dph)
        return
        end

        subroutine srclyr(depth,lmax,dph)
        implicit none
        real depth, dph
        integer lmax
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        integer mmax
        common/modlly/mmax

        integer m
        real dep
c-----
c       Find source/receiver boundary. It is assumed that
c       it will lie upon a boundary
c
c       lmax = source layer 
c            = 0 is the free surface 
c       depth = source depth 
c       dph = height of  source above lmax + 1 interface 
c       lmax = 0 is the free surface 
c-----
        if(depth.le.0.0)then
            lmax = 1
            dph = 0.0
        else
            dep = 0.0 
            do 100 m = 2,mmax
                dep = dep + d(m-1) 
                dph = dep - depth 
                lmax = m 
                if(abs(dph).lt. 0.0001*d(m-1) .or.
     1              abs(dph).lt.1.0e-6)go to 101
  100       continue 
  101   continue 
        endif
        return 
        end 

        subroutine svupdn(s,haa)
        implicit none
        complex*16 s(4), haa(4,4)
        complex*16 t(4), tmp
        integer i, j
        do 100 i=1,4
            tmp = dcmplx(0.0d+00,0.0d+00)
            do 101 j=1,4
                tmp = tmp + haa(i,j)*s(j)
  101       continue
            t(i) = tmp
  100   continue
        do 200 i=1,4
            s(i) = t(i)
  200   continue
        return
        end
        
        subroutine trans4(a)
c-----
c       from a transpose 
c-----
        implicit none
        complex*16 a(4,4), zdum
        integer i, j
        do 100 i=1,4
            do 101 j=i,4
                zdum = a(i,j)
                a(i,j) = a(j,i)
                a(j,i) = zdum
  101       continue
  100   continue
        return
        end

        subroutine var(p,q,r, rp, rsv, rsh,
     1      cosp, cosq, cosr, rsinp, rsinq, rsinr,
     1      sinpr, sinqr, sinrr,pex,svex,shex,iwat)
c-----
c       p = rp  * h
c       q = rsv * h
c       r = rsh * h
c       rp  vertical wave number for P
c       rsv vertical wave number for SV
c       rsh vertical wave number for SH
c       cosp=cosh(p)  rsinp =rp *sinh(p)  sinpr = sinh(p)/rp
c       cosq=cosh(q)  rsinsv=rsv*sinh(p)  sinpq = sinh(p)/rsv
c       cosr=cosh(r)  rsinsh=rsh*sinh(p)  sinpr = sinh(p)/rsh
c-----
        implicit none
        COMPLEX*16 p, q, r
        COMPLEX*16 rp, rsv, rsh
        COMPLEX*16 cosp, cosq, cosr
        COMPLEX*16 rsinp, rsinq, rsinr
        COMPLEX*16 sinpr, sinqr, sinrr
        REAL *8 pex,svex,shex
        integer iwat

        REAL*8 pr, pi, qr, qi, rr, ri
        COMPLEX*16 epp, epm, eqp, eqm, erp, erm
        COMPLEX*16 sinp, sinq, sinr

        REAL*8 PFAC, SVFAC, SHFAC
        
        pex  = 0.0d+00
        svex = 0.0d+00
        shex = 0.0d+00
        pr = dreal(p)
        pi = dimag(p)
        qr = dreal(q)
        qi = dimag(q)
        rr = dreal(r)
        ri = dimag(r)
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
            cosr  = 1.0d+00
            rsinr = 0.0d+00
            sinrr = 0.0d+00
            shex  = 0.0d+00
        else
c-----
c       elastic layer
c-----
            svex = qr
            shex = rr
            epp = dcmplx(dcos(pi), dsin(pi))/2.0
            epm = dconjg(epp)
            eqp = dcmplx(dcos(qi), dsin(qi))/2.0
            eqm = dconjg(eqp)
            erp = dcmplx(dcos(ri), dsin(ri))/2.0
            erm = dconjg(erp)
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

            if(rr.lt.15.) then
                shfac=dexp(-2.*rr)
            else
                shfac  = 0.0d+00
            endif
            cosr = erp + shfac*erm
            sinr = erp - shfac*erm
            rsinr = rsh*sinr
            sinrr = sinr/rsh
        endif
        return
        end

        subroutine velbnd() 
        implicit none
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
c-----
c       get bounds on earth model 
c-----
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        real vmin,vamin,vamax,vbmin,vbmax

        real va, vc, vf, vl, vn
        integer i
c-----
c       initialize bound search
c-----
        vamin = 1.0e+38
        vbmin = 1.0e+38
        vmin  = 1.0e+38
        vamax = 0.0
        vbmax = 0.0
        write(LOT,2) 
    2   format(' Working model'/
     1      ' ',7x,'d',8x,'va',8x,'vc',8x,'vf',8x,'vl',8x,'vn'
     1          ,8x,'rho',6x,'1/qa',6x,'1/qb')
    3   format(' ',7f10.3,2f10.6) 
        do 20 i = 1,mmax 
            va = sqrt(TA(i)/TRho(i))
            vc = sqrt(TC(i)/TRho(i))
            vf = sqrt(TF(i)/TRho(i))
            vl = sqrt(TL(i)/TRho(i))
            vn = sqrt(TN(i)/TRho(i))
            if(va.gt.vamax)vamax=va
            if(vc.gt.vamax)vamax=vc
            if(vl.gt.vbmax)vbmax=vl
            if(vn.gt.vbmax)vbmax=vn
            if(va.lt.vamin)vamin=va
            if(vc.lt.vamin)vamin=vc
            if(vl.lt.vbmin .and. vl.gt.0.0)vbmin=vl
            if(vn.lt.vbmin .and. vn.gt.0.0)vbmin=vn
            if(vl.gt.0.1 .or. vn.gt.0.1)then
                if(vl.lt.vmin)vmin=vl
                if(vn.lt.vmin)vmin=vn
            else
                if(vc.lt.vmin)vmin=vc
                if(va.lt.vmin)vmin=va
            endif
            if(i.lt.mmax)then
            write(LOT,3)d(i),va,vc,vf,vl,vn,TRho(i),qa(i),qb(i)
            else
            write(LOT,5)     va,vc,vf,vl,vn,TRho(i),qa(i),qb(i)
            endif
   20   continue 
    5   format(' ',10x,6f10.3,2f10.6/' ') 
c-----
c     obtain extreme velocity limits
c-----
      return 
      end 

        subroutine werror(ostr)
c-----
c       output error message and terminate program
c-----
        implicit none
        integer LER, LIN, LOT
        parameter(LER=0, LIN=5, LOT=6)
        character ostr*(*)
        write(LER,*)'PROGRAM TERMINATION'
        write(LER,*)ostr
        stop
        end

