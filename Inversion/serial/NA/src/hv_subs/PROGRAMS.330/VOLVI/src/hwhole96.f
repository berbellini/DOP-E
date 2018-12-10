      program hwhole96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME VI                                                      c
c                                                                     c
c      PROGRAM: HWHOLE96                                              c
c                                                                     c
c      COPYRIGHT 1985                                                 c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       CHANGES
c       11 SEP 2000 - build in P, SV and SH first arrival times
c       24 OCT 2000 - build in new formalism based on RBH notes
c       14 JAN 2001 - put in A, C, F, L and N constants into
c           trace header
c       27 AUG 2003 - put in far-field approximation
c-----
c-----
c     THIS PROGRAM COMPUTES THE SPECTRA OF THE HALFSPACE GREEN's
c     FUNCTION FOR DOUBLE COUPLE and EXPLOSION GREEN's FUNCTION
c     BY FOLLOWING HASKELL, N. A. (1963). Radiation pattern of
c         Rayleigh waves from a fault of arbitrary dip and 
c         direction of motion in a homogeneous medium, Bull.
c         Seism. Soc. Am. 53, 619-642.
c-----
        parameter(LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/damp/alpha,ieqex
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/modelt/dt(NL),at(NL),bt(NL),rhot(NL),mmaxt,
     1      qat(NL),qbt(NL)
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
        parameter(NDST=500)
        real*4 r(NDST), tshift(NDST), vred(NDST)
        complex*16 cresp(16)
        complex ztmp, zdata
        common/c/cmax,c1,c2,cmin
        common/frlim/fl,fu,df,fwhich
c-----
c       ishank  L   - .true. use Hankel function and not Bessel
c                    not that asymptotic tricks are not used
c       hnkarg  R*4 - (default 6.0) For kr > hnkarg use the Hankel
c                   function, otherwise the Bessel
c       dstcor  R*4 - 0 first time sample is tshift + r/vred
c                 1 first time sample is tshift + z/vred where
c                   z = abs (source depth - receiver depth)
c                 2 first time sample is tshift + 
c                   sqrt(z*z + r*r) /vred
c       dokjar  L   - .false. use Futterman causal Q
c                 .true.  use Kjartansson Causal Q
c-----
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar

        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        logical ext
        character mname*80, title*80

        real VSA, VSB, VSR
        REAL*4 SA, SC, SF, SL, SN, SR
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
c       See if the file hspec96.dat exists, if it does
c       open it for all control information
c-----
        inquire(file='hspec96.dat',exist=ext)
        if(.not. ext)then
                write(LER,*)'Control file ',
     1          'hspec96.dat does not exist'
                go to 9999
        endif
        open(3,file='hspec96.dat',access='sequential',
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
        write(LOT,2)  fl,fu,df,n1,n2,n
c-----
c       UNIX output - no carriage control
c-----
    2   format('fl =',f10.5,5x,'fu =',f10.5,5x,'df =',f10.5,/
     1      4x,'n1 =',i4,5x,'n2 =',i4,5x, ' n =',i5)
 2021   format('depths =',f20.6)
 2031   format('depthr =',f14.6)
 2041   format('     r =',f14.6,' tshift=',f10.2,' vred=',f10.2)
    4   format('frequencies for which response computed     ')
    5   format('alpha =',f10.5,5x,'dt =',f10.3)
c-----
c       MSDOS MICROSOFT FORTRAN OUTPUT
c-----
c    2   format(1x,'fl =',f10.5,5x,'fu =',f10.5,5x,'df =',f10.5,/
c     1          1x,4x,'n1 =',i4,5x,'n2 =',i4,5x, ' n =',i5)
c 2021   format(1x,'depths =',f20.6)
c 2031   format(1x,'depthr =',f14.6)
c 2041   format(1x,'     r =',f14.6,' tshift=',f10.2,' vred=',f10.2)
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
c       output the final spectrum as a function of distance
c-----
        open(unit=4,file='hspec96.grn',status='unknown',
     1      form='unformatted',access='sequential')
        rewind 4
c-----
c       oprog   1 hspec96
c           2 hspec96p
c           3 hwhole96
c-----
        iprog = 3
        write(4)iprog
        write(4) alpha,fl,fu,delt,n,n1,n2,df,nyq2
        write(4)mname
c-----
c       $now out the spectrum for each distacne
c
c           DIST
c               SOURCE_DEPTH
c                   RECEIVER_DEPTH
c                       FREQ
c
c-----
        do 5000 jd=1,ndist
        do 5005 js=1,mdpths
        do 5010 jr=1,mdpthr
c-----
c       new formalism has Z increasing downward
c-----
            depth = - (depths(js) - depthr(jr))
            call setup(depth,r(jd) )
            call gett0(t0,r(jd),depths(js),
     1          depthr(jr),tshift(jd),vred(jd),dstcor)
            DIST = SQRT(R(JD)**2 + ( DEPTHS(JS) - DEPTHR(JR))**2)
            TP  = DIST/A(1)
            if(B(1).gt.0.0)then
                TSV = DIST/B(1)
                TSH = DIST/B(1)
            else
                TSV = -12345.0
                TSH = -12345.0
            endif
            VSR = rho(1)
            VSA = a(1)
            VSB = b(1)
            SR = VSR
            SL = VSR*VSB*VSB
            SN = VSR*VSB*VSB
            SA = VSR*VSA*VSA
            SC = VSR*VSA*VSA
            SF = SA - 2.0*SN
            write(4)r(jd),t0,depths(js),
     1              depthr(jr),
     2              TP,TSV,TSH, 
     3              SA, SC, SF, SL, SN, SR
c-----
c       if the receiver is in a fluid, 
c       then permit pressure field output
c-----
            do 5011 i=1,21
                if(B(1) .gt. 0.0)then
                    if(i.lt.16)then
                        ksrc(i) = jsrc(lsrc(i))
                    else
                        ksrc(i) = 0
                    endif
                else
                    if(i.eq.9 .or. i.eq. 10 .or. i.eq.16)then
                        ksrc(i) = jsrc(lsrc(i))
                    else
                        ksrc(i) = 0
                    endif
                endif
 5011       continue
            write(4)ksrc
            do 5020 ii=n1,n2
                freq=(ii-1)*df
                if(freq.lt.df) freq = 0.01*df
                call wvint(r(jd),cresp,alpha,freq,a(1),b(1),
     1              rho(1),qa(1),qb(1),ieqex,depth,
     1              frefp(1),frefs(1))
                fac = 6.2831853*freq*t0
                ztmp = cmplx(cos(fac), sin(fac) )
                do 5200 jj=1,16
                    if(jsrc(lsrc(jj)).eq.1)then
                        zdata = ztmp * cresp(jj)
                        datar= real(zdata)
                        datai=aimag(zdata)
                        write(4)datar,datai
                    endif
 5200           continue
 5020       continue
 5010   continue
 5005   continue
 5000   continue
        rr = -1.0
        tt0 = 0.0
        write(4)rr,tt0
        close (4)
 9999   continue
        end

        subroutine gethsp(lun,mname,fl,fu,delt,n,n1,n2,xleng,xfac,
     1      ndist,r,tshift,vred)
c-----
c       read in data file of hspec8 commands
c-----
        implicit none
        integer lun
        character mname*80
        real fl, fu, delt, xleng, xfac
        integer n, n1, n2, ndist
        integer LER, LIN, LOT
        parameter(LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        integer NSOURCE, NRECEIVER, NSR 
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        integer jsrc, jbdrys, jbdryh
        integer NDST
        parameter(NDST=500)
        real*4 r(NDST), tshift(NDST), vred(NDST)
        common/c/cmax,c1,c2,cmin
        real cmax,c1,c2,cmin  
        character ostr*80
        character title*80
        common/lyrctl/lyrins
        logical lyrins
        logical ext
        integer ibdrys, ibdryh, jbdry, idcor, lmnm
        integer lgstr
        integer iunit, iiso, iflsph, idimen,icnvel,ierr
        real rr,tshf,vr
        integer i, ntmp, mtmp 
        real alphat, df

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar

c-----
c       UNIX FORTRAN - NO CARRIAGE CONTROL
c-----
   21   format(11i5)
   22   format(6e11.4)
   24   format('XLENG=',e15.7,' XFAC=',e15.7)
   30   format(2e15.7/1x,3i10)
 4012       format('WAVENUMBER FILTERING bounded in phase velocites'/
     1          '[cmax,c1,c2,cmin]=','[', f10.3, ',', 
     1          f10.3, ',', f10.3, ',', f10.3,']' /
     3          '( -1.0 means 0.0 for cmin and infinity for cmax)')
 4013   format('WAVENUMBER FILTERING NOT DONE')
c-----
c       MSDOS MICROSOFT FORTRAN OUTPUT - CARRIAGE CONTROL
c-----
c   21  format(1x,11i5)
c   22  format(1x,6e11.4)
c   24  format(1x,'XLENG=',e15.7,' XFAC=',e15.7)
c   30  format(1x,2e15.7/1x,3i10)
c 4012      format(1x,'WAVENUMBER FILTERING bounded in phase velocites'/
c     1         1x,'[cmax,c1,c2,cmin]=',
c     2         1x,'[', f10.3, ',', f10.3, ',', f10.3, ',', f10.3,']' /
c     3         1x,'( -1.0 means 0.0 for cmin and infinity for cmax)')
c 4013  format(1x,'WAVENUMBER FILTERING NOT DONE')
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
c       provide names for output Green s functions in order of output
c-----
        do 1233 i=16,21
            jsrc(i) = 0
 1233   continue
        if(ieqex.eq.0)then
            do 1234 i=1,16
                if(i.le.8)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1234       continue
            jsrc(13) = 1
            jsrc(14) = 1
            jsrc(16) = 1
        else if(ieqex.eq.1)then
            do 1235 i=1,16
                if(i.ge.7)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1235       continue
            jsrc(13) = 0
            jsrc(14) = 0
        else if(ieqex.eq.2)then
            do 1236 i=1,16
                jsrc(i) = 1
 1236       continue
        else if(ieqex.eq.3)then
            do 1237 i=1,16
                if(i.eq.7 .or. i.eq.8 .or. i.eq.16)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1237       continue
        else if(ieqex.eq.4)then
            do 1238 i=1,16
                if(i.le.6 .or. i.eq.13 .or. i.eq.14)then
                    jsrc(i) = 1
                else
                    jsrc(i) = 0
                endif
 1238       continue
        else if(ieqex.eq.5)then
            do 1239 i=1,16
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
            do 1240 i=1,16
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
        if(iiso.ne.0)then
            call usage('Isotropic velocity model required')
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
 2009   continue
c-----
c       check for filling the final depth array
c-----
        mtmp = mdpths * mdpthr
        ntmp = NSR
        if(mtmp .gt. ntmp)then
            write(ostr,3011)mtmp,ntmp
 3011   format(' NUMBER SOURCE-RECEIVER COMB',
     1          i5,' EXCEEDS DIMENSION ',i5)
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
c       read in phase velocity limits for wavenumber filtering
c       Wavenumber filtering will consist of following
c
c       |*ZERO*|-COSINE TAPER-|*ALL PASS*|-COSINE TAPER-|*ZERO
c       |      |              |          |              |
c            omega          omega      omega          omega
c k =   0    -----          -----      -----          -----   infinity
c            cmax            c1         c2             cmin
c-----
c       If c2 or cmin <= 0, then upper wavenumber limit is infinite
c       If c1 or cmax <= 0, then lower wavenumber limit is zero
c-----
        read(lun,*,end=4010,err=4010)cmax,c1,c2,cmin
            if(c1.le.0.0)cmax = -1.0
            if(c2.le.0.0)cmin = -1.0
            write(LOT,4012)cmax,c1,c2,cmin
        goto 4011
 4010   continue
            cmax = -1.0
            c1 = -1.0
            c2 = -1.0
            cmin = -1.0
        write(LOT,4013)
 4011   continue
c-----
c       verify the new model parameters
c-----
        write(LOT,*)'mmax=',mmax
        write(LOT,22)(d(i),a(i),b(i),rho(i),qa(i),qb(i),i=1,mmax)
c-----
c       Guarantee that no time wasted if any source is in the water
c-----
        if(b(1).le.1.0e-04)then
            do 2091 i=1,16
                jsrc(i) = 0
 2091       continue
            jsrc(7) = 1
            jsrc(8) = 1
            jsrc(16) = 1
        else
            jsrc(16) = 0
        endif
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
     
        subroutine setup(z,r)
        implicit real*8 (a-h,o-z)
        real*4 z,r
        common/zrinit/fac(34),r1
        r1=dsqrt(dble(r)**2 + dble(z)**2)
c-----
c       if distance == 0 , force small answers
c-----
        if(r1.le.0.0001d+00)r1=0.0001d+00
        r2 = r1**2
        r3 = r1*r2
        r4 = r2*r2
        r5 = r3*r2
        r6=r3*r3
        r7 = r3*r4
        zr3 = z/r3
        rr3 = r/r3
        z3r4 = z*z*z/r4
        zr4  = z/r4
        rzzr4 = r*z*z/r4
        rrzr4 = r*r*z/r4
        rr4 = r/r4
        r3r4 = r*r*r/r4
        r3r5 = r3r4/r1
        z3r5 = z3r4/r1
        zr5 = zr4/r1
        rr5 = rr4/r1
        rzzr5 = rzzr4/r1
        rrzr5 = rrzr4/r1
        z3r6 = z3r4/r2
        r3r6 = r3r4/r2
        rrzr6 = rrzr4/r2
        rzzr6 = rzzr4/r2
        z3r7 = z3r4/r3
        r3r7 = r3r4/r3
        rzzr7 = rzzr4/r3
        rrzr7 = rrzr4/r3
c-----
c       Fz
c-----
        fac(1) = zr3
        fac(2) = z/r2
c-----
c       Fr
c-----
        fac(3) = rr3
        fac(4) = r/r2
c-----
c       Frrr
c-----
        fac(5) = (-9.*rr5 + 15.*r3r7)
        fac(6) = (-9.*rr4 + 15.*r3r6)
        fac(7) = (-3.*rr3 + 6.*r3r5)
        fac(8) = ( r3r4)
c-----
c       Fzzz
c-----
        fac(9) = (-9.*zr5 + 15.*z3r7)
        fac(10)= (-9.*zr4 + 15.*z3r6)
        fac(11)= (-3.*zr3 + 6.*z3r5)
        fac(12)= ( z3r4)
c-----
c       Frrz
c-----
        fac(13)= (-3.*zr5 + 15.*rrzr7)
        fac(14)= (-3.*zr4 + 15.*rrzr6)
        fac(15)= (6.*rrzr5 - zr3)
        fac(16)= (rrzr4 )
c-----
c       Frzz
c-----
        fac(17)= (-3.*rr5 + 15.*rzzr7)
        fac(18)= (-3.*rr4 + 15.*rzzr6)
        fac(19)= (6.*rzzr5 - rr3)
        fac(20)= (rzzr4)
c-----
c       Frz
c-----
        rzr3 = r*zr3
        rzr4 = rzr3/r1
        rzr5 = rzr3/r2
        fac(21) = -3.*rzr5
        fac(22) = -3.*rzr4
        fac(23) = -rzr3
c-----
c       Fzz
c-----
        z2r5 = z * zr5
        z2r4 = z * zr4
        fac(24) = (3.*z2r5 - 1./r3)
        fac(25) = (3.*z2r4 - 1./r2)
        fac(26) = z * zr3
c-----
c       Frr
c-----
        r2r5 = r * rr5
        r2r4 = r * rr4
        fac(27) = (3.*r2r5 - 1./r3)
        fac(28) = (3.*r2r4 - 1./r2)
        fac(29) = r * rr3
c-----
c       Fror
c-----
        fac(30) = 1.0 / r3
        fac(31) = 1.0 / r2
c-----
c       Frzor
c-----
        fac(32) = -3.*zr5
        fac(33) = -3.*zr4
        fac(34) = -zr3
        return
        end


        subroutine wvint(r,cresp,alpha,freq,a,b,rho,qa,qb,ieqex,h,
     1      frefp,frefs)
        implicit none
        real r, h, alpha, freq, frefp, frefs
        integer ieqex
        
        complex*16 omega,fz(2),fr(2),frrr(2),frrz(2),frzz(2),fzzz(2)
        complex*16 frz(2),fzz(2),frr(2),f(2),fror(2),frzor(2)
        complex*16 cresp(16)
        complex*16 pfac(2),xka,xkb,atna,atnb
        complex*16 kv,eye
        complex*16 pi4ro
        complex*16 zero, zone
        real*8 pi4r
        real*4 a,b,rho,qa,qb


        common/zrinit/fac(34),r1
        real*8 fac,r1

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        common/solu/ipsvsh
        integer ipsvsh

        real*8 dr,di

        complex*16 ibr, sbr, dsbdhr
        complex*16 npfac(2)
        complex*16 nfsh, ffsh

        common/cdhoop/docagn
        logical docagn

        common/approx/dofar
        logical dofar

        real sgn
        integer i, iwat
c-----
c       initialize
c-----
        eye = dcmplx(0.0d+00,1.0d+00)
        zero = dcmplx(0.0d+00,0.0d+00)
        zone = dcmplx(1.0d+00,0.0d+00)
        do 100 i=1,16
            cresp(i) = zero
  100   continue

        omega = dcmplx(dble(6.2831853*freq),-dble(alpha))
        pi4r = dble(4.*3.1415927*rho)
        pi4ro = - pi4r*omega*omega

        call aten(omega,qa,qb,xka,xkb,alpha,a,b,atna,atnb,iwat,
     1      frefp,frefs)
        pfac(1) = -eye * xka * dcmplx(r1,0.0d+00)
        pfac(2) = -eye * xkb * dcmplx(r1,0.0d+00)
        dr = dreal(pfac(1))
        di = dimag(pfac(1))
        pfac(1) = dexp(dr)*dcmplx(dcos(di),dsin(di))
        dr = dreal(pfac(2))
        di = dimag(pfac(2))
        pfac(2) = dexp(dr)*dcmplx(dcos(di),dsin(di))

        npfac(1) = -eye * xka * dcmplx(dabs(dble(h)),0.0d+00)
        npfac(2) = -eye * xkb * dcmplx(dabs(dble(h)),0.0d+00)
        dr = dreal(npfac(1))
        di = dimag(npfac(1))
        npfac(1) = dexp(dr)*dcmplx(dcos(di),dsin(di))
        dr = dreal(npfac(2))
        di = dimag(npfac(2))
        npfac(2) = dexp(dr)*dcmplx(dcos(di),dsin(di))
c-----
c       For comparison to Cagniard-de Hoop, note that the non-causal
c       arrival never has an effect, so it could be 
c           dropped here for computation
c
c       However for comparison to omega-k integration we need it
c       and also need it for r = 0
c----- 
        if(docagn)npfac(2) = zero
        if(h.gt.0.0)then
            sgn = +1.0
        else if(h.lt.0.0)then
            sgn = -1.0
        else
            sgn = 0.0
        endif
        if(r.eq.0.0)then
            sbr = npfac(2)/(2.0*h)
            ibr = zero
            if(b.gt.0.0)then
            dsbdhr = sgn * npfac(2)*(-eye*xkb/2.0d+00)
            else
            dsbdhr = 0.0
            endif
        else
            if(b.gt.0.0)then
                sbr = (npfac(2) - pfac(2))/(eye*xkb*r*r)
                dsbdhr = sgn*(sgn*( npfac(2)
     1              -dabs(dble(h))
     2              *pfac(2)/dabs(r1)))/(r*r)
            else
                dsbdhr = 0.0
                sbr = 0.0
            endif
            ibr = (- pfac(2)/r1 
     1          +(2.0d+00/dabs(dble(r)))*sbr)/dble(r)
        endif
c-----
c       get partials for P potential
c-----
        kv = eye * xka
        call deriv(pfac(1),kv,fr(1),fz(1),frrr(1),frrz(1),
     1      frzz(1),fzzz(1),frz(1),fzz(1),frr(1),f(1),
     2      fror(1), frzor(1) )
c-----
c       get partials for S potential
c-----
        kv = eye * xkb
        call deriv(pfac(2),kv,fr(2),fz(2),frrr(2),frrz(2),
     1      frzz(2),fzzz(2),frz(2),fzz(2),frr(2),f(2),
     2      fror(2), frzor(2) )
            if(b.gt.0.0)then
c-----
c       P contribution
c-----
                if(ipsvsh.eq.0 .or. ipsvsh.eq.1)then
c-----
c       ZDD
c-----
                    cresp(1) =  cresp(1) - 
     1                  (3.0d+00*fzzz(1) + 
     1                  xka*xka*fz(1) 
     2                   )/pi4ro
c-----
c       RDD
c-----
                    cresp(2) = cresp(2) - 
     1                  (3.0d+00*frzz(1) 
     2                  + xka*xka*fr(1))/pi4ro
c-----
c       ZDS
c-----
                    cresp(3) = cresp(3) + 
     1                   (-2.0d+00*frzz(1)
     1                  )/pi4ro
c-----
c       RDS
c-----
                    cresp(4) = cresp(4) +
     1                  (-2.0d+00*(frrz(1)) 
     1                   )/pi4ro
c-----
c       TDS
c-----
                    cresp(5) = cresp(5) +
     1                  ( 2.0d+00*(frzor(1)) 
     1                  )/pi4ro
                    if(r.gt.0.0)then
c-----
c       ZSS
c-----
                    cresp(6) = cresp(6) - 
     1                  (2.0d+00*frrz(1)+fzzz(1)
     1                  +xka*xka*fz(1)
     1                  )/pi4ro
                    endif
c-----
c       RSS
c-----
                    cresp(7) = cresp(7) -
     1                  (2.0d+00*frrr(1)+frzz(1)+xka*xka*fr(1)
     2                  )/pi4ro
c-----
c       TSS
c-----
                    cresp(8) = cresp(8) +
     1                  ( -2.0d+00*frrr(1)-2.0d+00*frzz(1)
     1                  -2.0d+00*xka*xka*fr(1)
     3                  )/pi4ro
c-----
c       ZEX
c-----
                    cresp(9) = cresp(9) + 
     1                  (- fz(1))/
     2                  (4.0*3.1415927*rho*a*a*atna*atna)
c-----
c       REX
c-----
                    cresp(10) = cresp(10) + 
     1                   (- fr(1))/
     2                  (4.0*3.1415927*rho*a*a*atna*atna)
c-----
c       ZVF
c-----
                    cresp(11) = cresp(11) + 
     1                  (fzz(1) )/(pi4ro)
c-----
c       RVF
c-----
                    cresp(12) = cresp(12) + 
     1                    (frz(1) )/pi4ro
c-----
c       ZHF
c-----
                    cresp(13) = cresp(13) + 
     1                   (frz(1) )/pi4ro
c-----
c       RHF
c-----
                    cresp(14) = cresp(14) +
     1                  (  frr(1) )/pi4ro
c-----
c       THF
c-----
                    cresp(15) = cresp(15) - fror(1)/pi4ro


                endif
c-----
c       SV contribution
c-----
                if(ipsvsh.eq.0 .or. ipsvsh.eq.2 .or. ipsvsh.eq.4 )then
                    cresp(1) = cresp(1) - 
     1                  (
     1                  - 3.0d+00*fzzz(2)
     2                  - 3.0d+00*xkb*xkb*fz(2) )/pi4ro
                    cresp(2) = cresp(2) - 
     1                  (
     1                  - 3.0d+00*frzz(2) 
     2                  )/pi4ro
                    cresp(3) = cresp(3) + 
     1                   (
     1                  2.0d+00*frzz(2) 
     2                  + xkb*xkb*fr(2))/pi4ro
                    ffsh =   xkb*xkb*fz(2)/pi4ro
                    nfsh =   xkb*xkb*dsbdhr /(pi4ro)
                    cresp(4) = cresp(4) -
     1                   (2.0d+00*(-frrz(2)) 
     1                  - xkb*xkb*fz(2) )/pi4ro
     1                  -nfsh 
                    cresp(5) = cresp(5) +
     1                  ( 2.0d+00*(-frzor(2)) 
     1                  -xkb*xkb*fz(2))/pi4ro
     1                  -nfsh -ffsh
                    if(r.gt.0.0)then
                    cresp(6) = cresp(6) - 
     1                  (
     1                  -2.0d+00*frrz(2)-fzzz(2)
     1                  -xkb*xkb*fz(2))/pi4ro
                    endif
                    ffsh =  xkb*xkb*fr(2)/pi4ro
                    nfsh =  xkb*xkb*ibr*(2.0d+00)/pi4ro
                    cresp(7) = cresp(7) -
     1                  (
     1                  -2.0d+00*frrr(2)-frzz(2)
     2                  -2.0d+00*xkb*xkb*fr(2))/pi4ro
     1                  -nfsh
                    cresp(8) = cresp(8) +
     1                   ( 
     2                  + 2.0d+00*frrr(2)+2.0d+00*frzz(2)
     3                  +xkb*xkb*fr(2))/pi4ro
     1                  -nfsh - ffsh
                    cresp(11) = cresp(11) + 
     1                  (- fzz(2) - xkb*xkb*f(2))/(pi4ro)
                    cresp(12) = cresp(12) + 
     1                   ( - frz(2))/pi4ro
                    cresp(13) = cresp(13) + 
     1                    (- frz(2))/pi4ro
                    ffsh = - xkb*xkb*f(2)/pi4ro
                    nfsh =   xkb*xkb*sbr/(pi4ro)
                    cresp(14) = cresp(14) +
     1                  ( - frr(2) - xkb*xkb*f(2))/pi4ro
     1                  -nfsh
                    cresp(15) = cresp(15) -
     1                  ( -fror(2) - xkb*xkb*f(2))/pi4ro
     1                  -nfsh - ffsh 
                endif
c-----
c       SH contribution
c-----
                if(ipsvsh.eq.0 .or. ipsvsh.eq.3 .or. ipsvsh.eq.4 )then
                    ffsh =   xkb*xkb*fz(2)/pi4ro
                    nfsh =   xkb*xkb*dsbdhr /(pi4ro)
                    cresp(4) = cresp(4) +
     1                  nfsh 
                    cresp(5) = cresp(5) +
     1                  nfsh + ffsh
                    ffsh =  xkb*xkb*fr(2)/pi4ro
                    nfsh =  xkb*xkb*ibr*(2.0d+00)/pi4ro
                    cresp(7) = cresp(7) +
     1                  nfsh 
                    cresp(8) = cresp(8) +
     1                  nfsh + ffsh
                    ffsh = - xkb*xkb*f(2)/pi4ro
                    nfsh =   xkb*xkb*sbr/(pi4ro)
                    cresp(14) = cresp(14) 
     1                  +nfsh
                    cresp(15) = cresp(15) 
     1                  +nfsh +ffsh
                endif
                cresp(16) = dcmplx(0.0d+00,0.0d+00)
            else
                if(ipsvsh.eq.0 .or. ipsvsh.eq.1)then
                cresp(1) = dcmplx(0.0d+00,0.0d+00)
                cresp(2) = dcmplx(0.0d+00,0.0d+00)
                cresp(3) = dcmplx(0.0d+00,0.0d+00)
                cresp(4) = dcmplx(0.0d+00,0.0d+00)
                cresp(5) = dcmplx(0.0d+00,0.0d+00)
                cresp(6) = dcmplx(0.0d+00,0.0d+00)
                cresp(7) = dcmplx(0.0d+00,0.0d+00)
                cresp(8) = dcmplx(0.0d+00,0.0d+00)
                cresp(9) = - fz(1)/(4.0*3.1415927*rho*a*a*atna*atna)
                cresp(10) = - fr(1)/(4.0*3.1415927*rho*a*a*atna*atna)
                cresp(11) = dcmplx(0.0d+00,0.0d+00)
                cresp(12) = dcmplx(0.0d+00,0.0d+00)
                cresp(13) = dcmplx(0.0d+00,0.0d+00)
                cresp(14) = dcmplx(0.0d+00,0.0d+00)
                cresp(15) = dcmplx(0.0d+00,0.0d+00)
                cresp(16) = + f(1)/(4.*3.1415927*a*a*atna*atna)
                endif
            endif
c-----
c       the responses given above were based on a 
c           cylindrical coordinate system
c       with Z positive down. Convert all the Z coordinates 
c           to positive up
c       The tangential and radial displacements do not change
c-----
        cresp(1)  = - cresp(1)
        cresp(3)  = - cresp(3)
        cresp(6)  = - cresp(6)
        cresp(9)  = - cresp(9)
        cresp(11) = - cresp(11)
        cresp(13) = - cresp(13)
        return
        end

        subroutine deriv(pfac,kv,fr,fz,frrr,frrz,frzz,fzzz,frz,
     1  fzz,frr,f,fror,frzor)
        common/zrinit/fac(34),r1

        common/approx/dofar
        logical dofar

        real*8 fac,r1
        complex*16 pfac,kv,fr,fz,frrr,frrz,frzz,fzzz,frz,fzz,frr,f
        complex*16 fror, frzor

        if(dofar)then
            f  = pfac/r1
            fz = -pfac* ( kv * fac(2))
            fr = -pfac* ( kv * fac(4))
            fror = cmplx(0.0d+00,0.0d+00)
            frz = -pfac*(  kv*kv*fac(23))
            frzor = cmplx(0.0d+00,0.0d+00)
            fzz =  pfac*(  kv*kv*fac(26))
            frr =  pfac*(  kv*kv*fac(29))
            frrr = -pfac*( 
     1           kv*kv*kv*fac(8) )
            fzzz = -pfac*( 
     1           kv*kv*kv*fac(12) )
            frrz = -pfac*( 
     1           kv*kv*kv*fac(16) )
            frzz = -pfac*( 
     1           kv*kv*kv*fac(20) )
        else
            f  = pfac/r1
            fz = -pfac* (fac(1) + kv * fac(2))
            fr = -pfac* (fac(3) + kv * fac(4))
            fror = -pfac * ( fac(30) + kv * fac(31))
            frz = -pfac*( fac(21) + kv*fac(22) + kv*kv*fac(23))
            frzor = -pfac*( fac(32) + kv*fac(33) + kv*kv*fac(34))
            fzz =  pfac*( fac(24) + kv*fac(25) + kv*kv*fac(26))
            frr =  pfac*( fac(27) + kv*fac(28) + kv*kv*fac(29))
            frrr = -pfac*( fac(5) + kv*fac(6) + kv*kv*fac(7)
     1          + kv*kv*kv*fac(8) )
            fzzz = -pfac*( fac(9) + kv*fac(10) + kv*kv*fac(11)
     1          + kv*kv*kv*fac(12) )
            frrz = -pfac*( fac(13) + kv*fac(14) + kv*kv*fac(15)
     1          + kv*kv*kv*fac(16) )
            frzz = -pfac*( fac(17) + kv*fac(18) + kv*kv*fac(19)
     1          + kv*kv*kv*fac(20) )
        endif
        return
        end

        subroutine aten(omega,qa,qb,xka,xkb,alpha,a,b,atna,atnb,iwat,
     2      frefp,frefs)
c-----
c       make velocities complex, using Futterman causality operator
c-----
        real*4 qa,qb,alpha,a,b
        complex*16 omega,at,atna,atnb,xka,xkb
        real*8 pi, om1p, om1s, oml, fac, pi2
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
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
                atna = cmplx(1.0,0.0)
            else
                fac = pi2*gama
                rfac = sin(fac)/cos(fac)
                atna = dcmplx(1.0d+00,0.0d+00)/
     1              (( (omega/om1p)**dble(-gama) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
            endif
            if(b.gt.1.0e-04*a)then
                if(gamb.le.0.0)then
                    atnb = cmplx(1.0,0.0)
                else
                    fac = pi2*gamb
                    rfac = sin(fac)/cos(fac)
                    atnb = dcmplx(1.0d+00,0.0d+00)/
     1              (( (omega/om1s)**dble(-gamb) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
                endif
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
              if(CDABS(omega).gt.oml) at=CDLOG(omega/om1p)/pi
              if(CDABS(omega).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1p*fac))/pi
              endif
              atna=(1.+dble(qa)*at+dcmplx(0.0d+00,dble(qa/2.)))
            endif
            if(qb.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(omega).gt.oml) at=CDLOG(omega/om1s)/pi
              if(CDABS(omega).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1s*fac))/pi
              endif
               atnb=(1.+dble(qb)*at+dcmplx(0.0d+00,dble(qb/2.)))
            endif
        endif
        xka=omega/(dble(a)*atna)
        if(b.le.1.0e-04*a)then
            iwat = 1
            xkb = dcmplx(0.0d+00,0.0d+00)
        else
            iwat = 0
            xkb=omega/(dble(b)*atnb)
        endif
        return
        end

        subroutine gcmdln()
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
c       ishank  L   - .true. use Hankel function and not Bessel
c                    not that asymptotic tricks are not used
c       hnkarg  R*4 - (default 6.0) For kr > hnkarg use the Hankel
c                   function, otherwise the Bessel
c       dokjar  L   - .false. use Futterman causal Q
c                 .true.  use Kjartansson Causal Q
c-----
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        common/cdhoop/docagn
        logical docagn

        common/approx/dofar
        logical dofar

        common/solu/ipsvsh

        integer*4 mnmarg
        character*50 name

        ishank = .false.
        hnkarg = 6.0
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
        docagn = .false.
        dofar = .false.

        ipsvsh = 0
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-K')then
                dokjar = .true.
            else if(name(1:2).eq.'-C')then
                docagn = .true.
            else if(name(1:2).eq.'-F')then
                dofar = .true.
            else if(name(1:2).eq.'-P')then
                ipsvsh = 1
            else if(name(1:3).eq.'-SV')then
                ipsvsh = 2
            else if(name(1:3).eq.'-SH')then
                ipsvsh = 3
            else if(name(1:2).eq.'-S' .and. name(1:3).ne.'-SV'
     1          .and. name(1:3).ne.'-SH')then
                ipsvsh = 4
C           else if(name(1:3).eq.'-SU')then
C               spup = .true.
C               ssup = .true.
C               ssdn = .false.
C               spdn = .false.
C               dosud = .true.
C           else if(name(1:3).eq.'-SD')then
C               spup = .false.
C               ssup = .false.
C               ssdn = .true.
C               spdn = .true.
C               dosud = .true.
C           else if(name(1:5).eq.'-SPUP')then
C               spup = .true.
C               dosud = .true.
C           else if(name(1:5).eq.'-SSUP')then
C               ssup = .true.
C               dosud = .true.
C           else if(name(1:5).eq.'-SPDN')then
C               spdn = .true.
C               dosud = .true.
C           else if(name(1:5).eq.'-SSDN')then
C               ssdn = .true.
C               dosud = .true.
C           else if(name(1:5).eq.'-RPUP')then
C               rpup = .true.
C               dorud = .true.
C           else if(name(1:5).eq.'-RSUP')then
C               rsup = .true.
C               dorud = .true.
C           else if(name(1:5).eq.'-RPDN')then
C               rpdn = .true.
C               dorud = .true.
C           else if(name(1:5).eq.'-RSDN')then
C               rsdn = .true.
C               dorud = .true.
C           else if(name(1:3).eq.'-RD')then
C               rpup = .false.
C               rsup = .false.
C               rsdn = .true.
C               rpdn = .true.
C               dorud = .true.
C           else if(name(1:3).eq.'-RU')then
C               rpup = .true.
C               rsup = .true.
C               rsdn = .false.
C               rpdn = .false.
C               dorud = .true.
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
        if(hnkarg.lt.3.0)hnkarg=3.0
        return
        end

        subroutine usage(str)
c------
c       write out program syntax
c-----
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
        character str*(*)
        lstr = lgstr (str)
        write(LER,*)str(1:lstr)
        write(LER,*)'USAGE: ',
     1  'hwhole96  [-K] [-?] [-h] '
C     2     '[-SU] [-SD] [-SPUP] [-SSUP] [-SPDN] [-SSDN] ',
C     3     '[-RU] [-RD] [-RPUP] [-RSUP] [-RPDN] [-RSDN] '
        write(LER,*)
     1  '-K      (default Futterman) use Kjartansson Causal Q'
C       write(LER,*)
C     1 'The following govern wavefield at source. The default is',
C     2 ' the entire wavefield'
C       write(LER,*)
C     1 '-SU      (default whole wavefield) Compute only upgoing ',
C     2 '               wavefield from the source'
C       write(LER,*)
C     1 '-SD      (default whole wavefield) Compute only downgoing ',
C     2 '               wavefield from the source'
C       write(LER,*)
C     1 ' -SPUP  Include upward P at source'
C       write(LER,*)
C     1 ' -SSUP  Include upward S at source'
C       write(LER,*)
C     1 ' -SPDN  Include downward P at source'
C       write(LER,*)
C     1 ' -SSDN  Include downward S at source'
C       write(LER,*)
C     1 'The following govern wavefield at receiver. The default is',
C     2 ' the entire wavefield'
C       write(LER,*)
C     1 ' -RD    Only downgoing waves at receiver'
C       write(LER,*)
C     1 ' -RU    Only upgoing waves at receiver'
C       write(LER,*)
C     1 ' -RPUP  Include upward P at receiver'
C       write(LER,*)
C     1 ' -RSUP  Include upward S at receiver'
C       write(LER,*)
C     1 ' -RPDN  Include downward P at receiver'
C       write(LER,*)
C     1 ' -RSDN  Include downward S at receiver'
        write(LER,*)
     1  '-?                   Display this usage message'
        write(LER,*)
     1  '-h                   Display this usage message'
        stop 
        end

        subroutine gett0(t0,r,depths,depthr,tshift,vred,dstcor)
        real*4 t0, r, depths, depthr, tshift, vred
        integer*4 dstcor
c-----
c       compute time of first sample of the time series
c-----
            if(dstcor.eq.0)then
                rr = r
            else if(dstcor.eq.1)then
                rr = abs(depthr - depths)
            else if(dstcor.eq.2)then
                rr = sqrt(r*r + (depthr-depths)*(depthr-depths))
            endif
            if(vred.eq.0.0)then
                t0 = tshift 
            else
                t0 = tshift + rr/vred
            endif
        return
        end


