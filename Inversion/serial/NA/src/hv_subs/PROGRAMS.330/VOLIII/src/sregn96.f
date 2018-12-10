        program sregn96
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: SREGN96                                                c
c                                                                      c
c      COPYRIGHT 1996, 2010                                            c
c      R. B. Herrmann                                                  c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c
c
c     This program calculates the group velocity and partial
c     derivatives of Love waves for any plane multi-layered
c     model.  The propagator-matrix, instead of numerical-
c     integration method is used, in which the Haskell rather
c     than Harkrider formalisms are concerned.
c
c     Developed by C. Y. Wang and R. B. Herrmann, St. Louis
c     University, Oct. 10, 1981.  Modified for use in surface
c     wave inversion, with addition of spherical earth flattening
c     transformation and numerical calculation of group velocity
c     partial derivatives by David R. Russell, St. Louis
c     University, Jan. 1984.
c
c     Rewrite of theory to agree more with hspec96 wavenumber
c     integration code
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c Revision history:
c       07 AUG 2002 - make string lengths 120 characters from 80 
c       13 OCT 2006 - verbose output of energy integrals
c       26 SEP 2008 - fixed undefined LOT in subroutine up
c       14 JUN 2009 - reformulate in general terms as per book
c       01 AUG 2010 - change code to agree with book  for
c                     migration to TI. Also clean up code using
c                     implicit none
c       20 SEP 2012 - cleaned up some compiler warnings
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit none
c-----
c       common blocks
c-----
        integer NL
        parameter(NL=1000)
c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c-----
        common/isomod/di(NL),ai(NL),bi(NL),rhoi(NL),
     1      qai(NL),qbi(NL),etapi(NL),etasi(NL), 
     2      frefpi(NL), frefsi(NL)
        real*4 di, ai, bi, rhoi, qai, qbi, etapi, etasi, frefpi, frefsi

        common/model/  d(NL),a(NL),b(NL),rho(NL),qa(NL), qb(NL)
        real*4 d, a, b, rho, qa, qb
        common/depref/refdep
        real*4 refdep
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs, xmu, xlam
        integer mmax

        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0, dcda, dcdb, dcdr, dcdh

        real*4 sdcda(NL), sdcdb(NL), sdcdh(NL), sdcdr(NL) 
        real*4 spur(NL), sptr(NL), spuz(NL), sptz(NL)

        common/wateri/iwat(NL)
        integer iwat
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid

        common/sumi/   sumi0,sumi1,sumi2,sumi3,flagr,are,ugr
        real*8 sumi0, sumi1, sumi2, sumi3, flagr, are, ugr

        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        real*4 vtp,dtp,rtp

c-----
c       function prototype
c-----
        integer lgstr

c-----
c       local variables
c-----
        integer LER, LIN, LOT
        parameter(LER=0,LIN=5,LOT=6)
        integer ipar(10)
        integer nipar(10)
        integer i,j,k,l
        integer iunit,iiso,iflsph,idimen,icnvel,ierr
        integer MAXMOD
        parameter(MAXMOD=2000)
        integer lss, lsso, lrr, lrro, n1, n2, npts, ifirst
        integer ifunc, mode, nsph, nper, nmodes, mmaxot

        logical dotmp
        logical ext
        logical dogam
        logical dderiv
        logical verbose
        logical nwlyrs, nwlyrr
        logical dolove, dorayl

        character mname*120
        character hsfile*120, hrfile*120, title*120 
        character*12 fname(2)

        real*4 fpar(10)
        real*4 s1, dephs, dephr, faclov, facray, hr, hs, dt
        real*4 deplw, depthr, depths, depup
        real*4 ohs, ohr
        real*4 rare, rtr, rtz, rur0, rur, ruz
        real*4 sare, sd2ur, sd2uz, sduz, sur, sur0, suz, sdur
        real*4 sum, sumgr, sumgv, sumkr
        real*4 wvnsrc, wvnrec, xl2m
        real*4 twopi

        real*8 durdz, duzdz, d2urdz, d2uzdz
        real*8 rurdz, ruzdz
        real*8 c, omega, wvno, gammar, csph, usph
        real*8 cp(MAXMOD), t


c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line information
c-----  
        call gcmdln(hsfile,hrfile,hs,hr,dotmp,dogam,dderiv,
     1      nipar,verbose)
c-----
c       wired in file names     - output of this program 
c                           is always in the
c               binary file sregn96.egn or sregn96.der
c                               - input from sdisp96 is always
c               binary file sdisp96.ray
c-----
        if(dotmp)then
            fname(1) = 'tsdisp96.ray'
        else
            fname(1) = 'sdisp96.ray'
        endif
        if(dderiv)then
            fname(2) = 'sregn96.der'
        else
            fname(2) = 'sregn96.egn'
        endif
c-----
c       get control parameters from sdisp96.dat
c-----
        inquire(file='sdisp96.dat',exist=ext)
        if(.not.ext)then
            call usage('Control file sdisp96.dat'//
     1          ' does not exist')
        endif
        open(1,file='sdisp96.dat',form='formatted',status='unknown',
     1      access='sequential')
        rewind 1
        read(1,*,end=9999,err=9999)dt
        read(1,*,end=9999,err=9999)npts,n1,n2
        read(1,'(a)',end=9999,err=9999)mname
        read(1,*)dolove, dorayl
        read(1,*)ohs, ohr
        read(1,*)nmodes
        read(1,*)faclov,facray
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
c-----
c       check for water only
c-----
        allfluid = .true.
        do 1200 i=1,mmax
            if(bi(i).gt.0.0)then
                allfluid = .false.
            endif
 1200   continue
c-----
c       get the Q information into the program
c       since this is not carried through in the sdisp96 output
c-----
        do 1234 i=1,mmax
C            qai(i) = qa(i)
C            qbi(i) = qb(i)
            if(.not.dogam)then
                qai(i) = 0.0
                qbi(i) = 0.0
            endif
            if(qai(i).gt.1.0)qai(i)=1.0/qai(i)
            if(qbi(i).gt.1.0)qbi(i)=1.0/qbi(i)
            zqai(i) = qai(i)
            zqbi(i) = qbi(i)
            if(frefpi(i).le.0.0)frefpi(i) = 1.0
            if(frefsi(i).le.0.0)frefsi(i) = 1.0
            zfrefp(i) = frefpi(i)
            zfrefs(i) = frefsi(i)
            zetap(i) = etapi(i)
            zetas(i) = etasi(i)
 1234   continue
        nsph = iflsph
c-----
c       get source depth, and sphericity control
c-----
        if(hs.lt.-1.0E+20)then
            depths = ohs
        else
            depths = hs
        endif
        if(hr.lt.-1.0E+20)then
            depthr = ohr
        else
            depthr = hr
        endif
        depthr = depthr + refdep
        depths = depths + refdep
c-----
c       if there is a spherical model, map the source 
c            and receiver depth 
c       in the spherical model to the equivalent depth in the flat model
c-----
        if(nsph.gt.0)then
                dephs = 6371.0 *alog( 6371.0/(6371.0 - (depths-refdep)))
                dephr = 6371.0 *alog( 6371.0/(6371.0 - (depthr-refdep)))
        else
                dephs = depths
                dephr = depthr
        endif
c----
c       see if the file sdisp96.ray exists
c-----
        inquire(file=fname(1),exist=ext)
        if(.not.ext)then
            call usage('Dispersion file: '//fname(1)//
     1          ' does not exist')
        endif
c-----
c       open output file of sdisp96 and of sregn96
c-----
        open(1,file=fname(1),form='unformatted',status='unknown',
     1          access='sequential')
        open(2,file=fname(2),form='unformatted',status='unknown',
     1          access='sequential')
        rewind 1
        rewind 2
c-----
c       obtain the earth model: note if the original was spherical,
c       this will be the transformed model
c-----
        call gtsmdl(1,mmax,d,a,b,rho,qa,qb,nper,
     2              mname,ipar,fpar)
c-----
c       define modulus of rigidity, also get current transformed model
c       parameters
c       set fpar(1) = refdep
c       set ipar(1) = 1 if medium is spherical
c       set ipar(2) = 1 if source is in fluid
c       set ipar(3) = 1 if receiver is in fluid
c       set ipar(4) = 1 if eigenfunctions are output with -DER flag
c       set ipar(5) = 1 if dc/dh are output with -DER flag
c       set ipar(6) = 1 if dc/da are output with -DER flag
c       set ipar(7) = 1 if dc/db are output with -DER flag
c       set ipar(8) = 1 if dc/dr are output with -DER flag
c-----
        deplw = 0.0
        depup = 0.0
        ipar(2) = 0
        ipar(3) = 0
        do 185 i=1,mmax
            zd(i) = d(i)
            za(i) = a(i)
            zb(i) = b(i)
            zrho(i) = rho(i)
                xmu(i)=dble(rho(i)*b(i)*b(i))
            xlam(i)=dble(rho(i)*(a(i)*a(i)-2.*b(i)*b(i)))
            depup = deplw + d(i)
            if(b(i) .lt. 0.0001*a(i))then
                iwat(i) = 1
                if(depths.ge.deplw .and. depths.lt.depup)then
                    ipar(2) = 1
                endif
                if(depthr.ge.deplw .and. depthr.lt.depup)then
                    ipar(3) = 1
                endif
            else
                iwat(i) = 0
            endif
            deplw = depup
  185   continue
        do 186 i=4,8
            ipar(i) = nipar(i)
  186   continue
        call putmdl(2,mmax,d,a,b,rho,qa,qb,nper,depths-refdep,
     2      depthr-refdep,mname,ipar,fpar)
c-----
c               If a spherical model is used, reconstruct the thickness
c               of the original model to get the mapping from
c               spherical to flat. We need this for correct
c               spherical partial derivatives.
c
c               THIS MUST BE DONE BEFORE SOURCE LAYER INSERTION
c-----
        if(nsph.gt.0)then
                call bldsph()
        endif
c-----
c       split a layer at the source depth
c       the integer ls gives the interface at which the eigenfunctions
c       are to be defined
c-----
        call insert(dephs,nwlyrs,lsso)
        call insert(dephr,nwlyrr,lrro)
        call srclyr(dephs,lss)
        call srclyr(dephr,lrr)
c-----DEBUG
c       output the new model with source layer
c-----
C       write(6,*)'lss,lrr:',lss,lrr
C        do 9902 i=1,mmax
C                write(6,*)iwat(i),zd(i),za(i),zb(i),zrho(i),xmu(i)
C 9902   continue
c-----
c       EnD DEBUG
c-----
        twopi=2.*3.141592654
        ifirst = 0
  400   continue
        call gtshed(1,ifunc,mode,t,ierr)
        if(ierr.ne.0)go to 700
        s1=t
c-----
c       DEBUG OUTPUT
c-----
c        write(6,*) ifunc,mode,s1
c-----
c       END DEBUG OUTPUT
c-----
        call puthed(2,ifunc,mode,s1)
        omega=twopi/t
        if(ifunc.lt.0) go to 700
        if(mode.le.0) go to 400
        read(1) (cp(k),k=1,mode)
        do 600 k=1,mode
                c=cp(k)
c-----
c       main part.
c-----
            wvno=omega/c
            call svfunc(omega,wvno)
            call energy(omega,wvno,mmax)
c-----
c       the gamma routine will use the spherical model, but the
c       frequency dependence and Q of the original model
c-----
C        WRITE(6,*)'c(noq)=',c
            if(dogam)then
                call gammap(omega,wvno,gammar)
                c = omega/wvno
C        WRITE(6,*)'c(q)=',c
            else
                gammar = 0.0d+00
            endif
c------
c     also check for possible conversion errors in IEEE
c     conversion from double precision to single precision
c-----
            if(dabs(uu0(1)).lt.1.0d-36)uu0(1)=0.0d+00
            if(dabs(uu0(3)).lt.1.0d-36)uu0(3)=0.0d+00
            if(dabs(c).lt.1.0d-36)c=0.0d+00
            if(dabs(ugr).lt.1.0d-36)ugr=0.0d+00
            if(dabs(are).lt.1.0d-36)are=0.0d+00
            if(dabs(flagr).lt.1.0d-36)flagr=0.0d+00
            if(dabs(gammar).lt.1.0d-36)gammar=0.0d+00

c-----
c       output necessary eigenfunction values for
c       source excitation
c-----
            mmaxot = mmax
            if(nsph.gt.0)then
c-----
c           sphericity correction for partial derivatives 
c           of the original model
c-----
C        WRITE(6,*)'c(flat)=',c,' u(flat)=',ugr
                call sprayl(omega,c,mmaxot,csph,usph,ugr)
                wvno = omega / csph
C        WRITE(6,*)'c(sph )=',csph,' u(sph )=',usph

            endif

c-----
c      get the derivatives of the eitenfunctions required
c      for source excitation from the definition of stress. For
c      completeness get the second derivative from the first
c      derivatives and the equation of motion for the medium
c-----
            xl2m = xlam(lss) + 2.0d+00 * xmu(lss)
            duzdz = ( tz(lss) + wvno*xlam(lss)*ur(lss))/xl2m
            if(iwat(lss).eq.1)then
                durdz = wvno*uz(lss)
                d2urdz = wvno*duzdz
            else
                durdz = -wvno*uz(lss) + tr(lss)/xmu(lss)
                d2urdz = (-wvno*duzdz*(xlam(lss)+xmu(lss)) 
     1              - ur(lss)*(zrho(lss)*omega*omega -
     2              wvno*wvno*xl2m))/xmu(lss)
            endif
            d2uzdz = ( - uz(lss) *( zrho(lss)*omega*omega - 
     1          xmu(lss)*wvno*wvno) + 
     2          wvno*durdz*(xlam(lss)+xmu(lss)))/ xl2m
            if(dabs(duzdz).lt.1.0d-36)duzdz=0.0d+00
            if(dabs(durdz).lt.1.0d-36)durdz=0.0d+00
            if(dabs(d2uzdz).lt.1.0d-36)d2uzdz=0.0d+00
            if(dabs(d2urdz).lt.1.0d-36)d2urdz=0.0d+00

            xl2m = xlam(lrr) + 2.0d+00 * xmu(lrr)
            ruzdz = ( tz(lrr) + wvno*xlam(lrr)*ur(lrr))/xl2m
            if(iwat(lrr).eq.1)then
                rurdz = wvno*uz(lrr)
            else
                rurdz = -wvno*uz(lrr) + tr(lrr)/xmu(lrr)
            endif
            if(dabs(ruzdz).lt.1.0d-36)ruzdz=0.0d+00
            if(dabs(rurdz).lt.1.0d-36)rurdz=0.0d+00

            sur = sngl(ur(lss))
            sdur = sngl(durdz)
            sd2ur = sngl(d2urdz)
            suz = sngl(uz(lss))
            sduz = sngl(duzdz)
            sd2uz = sngl(d2uzdz)
            sare = sngl(are)
            wvnsrc = sngl(wvno)
            sur0 = sngl(uu0(1))

            if( verbose ) then
                 WRITE(LOT,2)t,c,ugr,gammar,
     1                sumi0,sumi1,sumi2,sumi3,
     2                flagr,are,sur0
    2    format(' T=',e15.7,'  C=',e15.7,'  U=',e15.7,'  G=',e15.7/
     1          'I0=',e15.7,' I1=',e15.7,' I2=',e15.7,' I3=',e15.7/
     2          ' L=',e15.7,' AR=',e15.7,'  E=',f15.7)
            endif

            rur = sngl(ur(lrr))
            rtr = sngl(tr(lrr))
            ruz = sngl(uz(lrr))
            rtz = sngl(tz(lrr))
            rare = sngl(are)
            wvnrec = sngl(wvno)
            rur0 = sngl(uu0(1))

            sumkr = 0.0
            sumgr = 0.0
            sumgv = 0.0
            if(nsph.gt.0)then
                wvno = omega/ csph
                ugr = usph
            endif
        if(dderiv)then
c-----
c               if a layer was inserted, get the partial
c               derivative for the original layer.
c           The sequence is cannot be changed.
c           Originally the model grows because the source layer
c           is added and then the receiver layer
c           This we must first strip the receiver and lastly the
c           source
c----
c       initialize
c-----
C           WRITE(6,*)'nwlyrr,lrro:',nwlyrr,lrro
C           WRITE(6,*)'nwlyrs,lsso:',nwlyrs,lsso
            if(nwlyrr)then
                call collap(lrro+1,mmaxot)
            endif
            if(nwlyrs)then
                call collap(lsso+1,mmaxot)
            endif
            call chksiz(ur,spur,mmaxot)
            call chksiz(tr,sptr,mmaxot)
            call chksiz(uz,spuz,mmaxot)
            call chksiz(tz,sptz,mmaxot)
            call chksiz(dcdh,sdcdh,mmaxot)
            call chksiz(dcdb,sdcdb,mmaxot)
            call chksiz(dcda,sdcda,mmaxot)
            call chksiz(dcdr,sdcdr,mmaxot)
c-----
c           up to this point the dcdh are changes to phase velocity if
c           if the layer boundary changes. Here we change this to mean
c           the dc/dh for a change in layer thickness
c
c           A layer becomes thicker if the base increases and the top
c           decreases its position. The dcdh to this point indicates 
c           the effect of moving a boundary down. Now we convert to
c           the effect of changing a layer thickness.
c-----
            do 505 i=1,mmaxot-1
 
                sum = 0.0
                do 506 j=i+1,mmaxot
                    sum = sum + sdcdh(j)
  506           continue
                sdcdh(i) = sum
  505       continue
            sdcdh(mmaxot) = 0.0

            call putder(2,6,sngl(wvno),sngl(ugr), 
     1          sngl(gammar), 
     1          sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0,
     2          rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3          sumkr,sumgr,sumgv,mmaxot,
     4          sdcdh,sdcda,sdcdb,sdcdr,
     5          spur,sptr,spuz,sptz,ipar)
        else
            call putegn(2,2,1,sngl(wvno),sngl(ugr),
     1          sngl(gammar),
     1          sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2          rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3          sumkr,sumgr,sumgv)
        endif
c-----
c       DEBUG OUTPUT
c-----
c           write(6,*) wvno,c,ugr,are
c           write(6,*) uu(ls),dut
cc-----
c       END DEBUG OUTPUT
c-----
  600   continue
        go to 400
  700   continue
c-----
c       close input file from sdisp96 and output file of this program
c-----
        do 900 i=1,2
                close(i,status='keep')
  900   continue

 9999   continue
        end

        subroutine insert(dph,newlyr,ls)
        implicit none
c-----
c       procdure arguments
c-----
        real*4 dph
        logical newlyr
        integer ls
c-----
c       common blocks
c-----
        integer NL
        parameter(NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs, xmu, xlam
        integer mmax
        integer iwat
        common/wateri/iwat(NL)
c-----
c       local variables
c------
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        real*8 dep, dp, dphh, hsave
        integer m, m1

c-----
c       Insert a depth point into the model by splitting a layer
c       so that the point appears at the top boundary of the layer
c       dph = depth of interest
c       newlyr  - L .true. layer added
c               .false. no layer added to get source eigenfunction
c-----
c-----
c       determine the layer in which the depth dph lies.
c       if necessary, adjust  layer thickness at the base
c-----
c       Here determine layer corresponding to specific depth dph
c       If the bottom layer is not thick enough, extend it
c
c       dep     - depth to bottom of layer
c       dphh    - height of specific depth above bottom of the layer
c-----
        dep = 0.0
        dp = 0.0
        dphh=-1.0
        ls = 1
        do 100 m =1, mmax
                dp = dp + zd(m)
                dphh = dp - dph
                if(m.eq.mmax)then
                        if(zd(mmax).le.0.0d+00 .or. dphh.lt.0.0)then
                                zd(mmax) = (dph - dp)
                        endif
                endif
                dep = dep + zd(m)
                dphh = dep - dph
                ls = m
                if(dphh.ge.0.0) go to 101
  100   continue
  101   continue
c-----
c       In the current model, the depth point is in the ls layer
c       with a distance dphh to the bottom of the layer
c
c       Do not create unnecessary layers, e.g., 
c            at surface and internally
c       However do put in a zero thickness layer 
c            at the base if necessary
c-----
        if(dph .eq. 0.0)then
            newlyr = .false.
                return
        else if(dphh .eq. 0.0 .and. ls.ne.mmax)then
            ls = ls + 1
            newlyr = .false.
                return
        else
            newlyr = .true.
c-----
c               adjust layering
c-----
                 do 102 m = mmax,ls,-1
                       m1=m+1
                        zd(m1) = zd(m)
                        za(m1) = za(m)
                        zb(m1) = zb(m)
                        zrho(m1) = zrho(m)
                        zqai(m1) = zqai(m)
                        zqbi(m1) = zqbi(m)
                        zfrefp(m1) = zfrefp(m)
                        zfrefs(m1) = zfrefs(m)
                        zetap(m1) = zetap(m)
                        zetas(m1) = zetas(m)
                        xmu(m1)=xmu(m)
                        xlam(m1)=xlam(m)
                iwat(m1) = iwat(m)
  102           continue
                hsave=zd(ls)
                zd(ls) = hsave - dphh
                zd(ls+1) = dphh
                mmax = mmax + 1
        endif
        return
        end

        subroutine srclyr(depth,lmax)
        parameter (LER=0, LIN=5, LOT=6)
        parameter(NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs, xmu, xlam
        common/wateri/iwat(NL)
c-----
c       Find source/receiver boundary. It is assumed that
c       it will lie upon a boundary
c
c       lmax = source layer 
c       depth = source depth 
c-----
        dep = 0.0
        do 100 i=1,mmax
            if(abs(depth - dep).le.0.001*zd(i))then
                lmax = i
                return
            endif
            dep = dep + zd(i)
  100   continue
        return
        end 

        subroutine gcmdln(hsfile,hrfile,hs,hr,dotmp,dogam,dderiv,
     1      nipar,verbose)
c-----
c       parse the command line arguments
c-----
c       hsfile  C*120   - name of source depth file
c       hrfile  C*120   - name of receiver depth file
c       hs  R*4 source depth (single one specified
c       hr  R*4 receiver depth (single one specified
c       dotmp   L   - .true. use file tsdisp96.lov
c       dogam   L   - .true. incorporate Q
c       dderiv  L   - .true. output depth dependent values
c       nipar   I*4 - array o integer controls
c           set nipar(4) = 1 if eigenfunctions are output with -DER flag
c           set nipar(5) = 1 if dc/dh are output with -DER flag
c           set nipar(6) = 1 if dc/da are output with -DER flag
c           set nipar(7) = 1 if dc/db are output with -DER flag
c           set nipar(8) = 1 if dc/dr are output with -DER flag
c       verbose L   - .true. output information on energy integrals
c-----
        implicit none
c-----
c       procedure arguments
c-----
        character hrfile*120, hsfile*120
        real hs, hr
        logical dotmp, dogam, dderiv
        integer nipar(10)
        logical verbose

c-----
c       procedure prototypes
c-----
        integer mnmarg
c-----
c       local arguments
c-----
        character name*40
        integer i, nmarg

        hrfile = ' '
        hsfile = ' '
        hs = -1.0E+21
        hr = -1.0E+21
        dotmp = .false.
        dogam = .true.
        dderiv = .false.
        nipar(4) = 0
        nipar(5) = 0
        nipar(6) = 0
        nipar(7) = 0
        nipar(8) = 0
        verbose = .false.
        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
                call mgtarg(i,name)
                if(name(1:3).eq.'-HR')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')hr
                else if(name(1:3).eq.'-HS')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')hs
C               else if(name(1:4).eq.'-FHR')then
C                   i = i + 1
C                   call mgtarg(i,hrfile)
C               else if(name(1:4).eq.'-FHS')then
C                   i = i + 1
C                   call mgtarg(i,hsfile)
                else if(name(1:2) .eq. '-T')then
                    dotmp = .true.
                else if(name(1:4) .eq. '-NOQ')then
                    dogam = .false.
                else if(name(1:4) .eq. '-DER')then
                    nipar(4) = 1
                    nipar(5) = 1
                    nipar(6) = 1
                    nipar(7) = 1
                    nipar(8) = 1
                    dderiv = .true.
                else if(name(1:3).eq.'-DE')then
                    dderiv = .true.
                    nipar(4) = 1
                else if(name(1:3).eq.'-DH')then
                    dderiv = .true.
                    nipar(5) = 1
                else if(name(1:3).eq.'-DB')then
                    dderiv = .true.
                    nipar(7) = 1
                else if(name(1:3).eq.'-DR')then
                    dderiv = .true.
                    nipar(8) = 1
                else if(name(1:2).eq.'-V')then
                    verbose = .true.
                else if(name(1:2) .eq. '-?')then
                    call usage(' ')
                else if(name(1:2) .eq. '-h')then
                    call usage(' ')
                endif
        go to 1000
 2000   continue
        return
        end

        subroutine usage(ostr)
        implicit none
c-----
c       procedure arguments
c-----
        character ostr*(*)
c-----
c       local variables
c-----
        integer LER
        parameter (LER=0)
        integer lostr
c-----
c       function prototype
c-----
        integer lgstr

        if(ostr .ne. ' ')then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'Usage: sregn96 ',
C    1      ' -FHR recdep -FHS srcdep -HS hs -HR hr ' ,
     1      ' -HS hs -HR hr ' ,
     2       ' [-NOQ] [-T] [-DER -DE -DH -DB -DR ] [-?] [-h]'
C       write(LER,*)
C     1 '-FHS srcdep (overrides -HS )  Name of source depth  file'
C       write(LER,*)
C     1 '-FHR recdep (overrides -HR )  Name of receiver depth  file'
        write(LER,*)
     1  '-HS hs      (default 0.0 )  Source depth '
        write(LER,*)
     1  '-HR hr      (default 0.0 )  Receiver depth'
        write(LER,*)
     1  '-NOQ        (default Q used) perfectly elastic'
        write(LER,*)
     1  '-T          (default false) use tdisp96.lov not disp96.lov'
        write(LER,*)
     1  '-DER        (default false) output depth dependent values'
        write(LER,*)
     1  '-DE         (default false) output eigenfunctions(depth)'
        write(LER,*)
     1  '-DH         (default false) output DC/DH(depth)'
        write(LER,*)
     1  '-DB         (default false) output DC/DB(depth)'
        write(LER,*)
     1  '-DR         (default false) output DC/DR(depth)'
        write(LER,*)
     1  '-V          (default false) list energy integrals'
        write(LER,*)
     1  '-?          (default none )  this help message '
        write(LER,*)
     1  '-h          (default none )  this help message '
        stop
        end

        subroutine collap(ls,mmaxot)
        implicit none
c-----
c       routine arguments
c-----
        integer ls, mmaxot
c-----
c       common blocks
c-----
        integer NL
        parameter(NL=1000)
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0, dcda, dcdb, dcdr, dcdh
c-----
c       local arguments
c-----
        integer i

        do 501 i = ls-1,mmaxot
            if(i .eq. ls -1)then
                dcda(i) = dcda(i) + dcda(i+1)
                dcdb(i) = dcdb(i) + dcdb(i+1)
                dcdr(i) = dcdr(i) + dcdr(i+1)
            endif
            if(i.gt.ls)then
                dcda(i-1) = dcda(i)
                dcdb(i-1) = dcdb(i)
                dcdh(i-1) = dcdh(i)
                dcdr(i-1) = dcdr(i)
                ur(i-1) = ur(i)
                tr(i-1) = tr(i)
                uz(i-1) = uz(i)
                tz(i-1) = tz(i)
            endif
  501   continue
        mmaxot = mmaxot - 1
        return
        end

        subroutine chksiz(dp,sp,mmaxot)
c-----
c       correctly convert double precision to single precision
c-----
        real*8 dp(mmaxot)
        real*4 sp(mmaxot)
            do 610 i=1,mmaxot
                if(dabs(dp(i)).lt.1.0d-36)then
                    sp(i) = 0.0
                else
                    sp(i) = sngl(dp(i))
                endif
  610       continue
        return
        end

        subroutine sprayl(om,c,mmax,csph,usph,ugr)
c-----
c       Transform spherical earth to flat earth
c       and relate the corresponding flat earth dispersion to spherical
c
c       Schwab, F. A., and L. Knopoff (1972). 
c          Fast surface wave and free
c       mode computations, in  
c          Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c          B. A. Bolt (ed),
c       Academic Press, New York
c
c       Rayleigh Wave Equations 111, 114 p 144
c
c       Partial with respect to parameter uses the relation
c       For phase velocity, for example,
c
c       dcs/dps = dcs/dpf * dpf/dps, c is phase velocity, p is
c       parameter, and f is flat model
c
c       om      R*8     angular frequency
c       c       R*8     phase velocity
c       mmax    I*4     number of layers 
c-----
        implicit none
c-----
c       procedure arguments
c-----
        real*8 om, c, csph,usph,ugr
        integer mmax
c-----
c       common blocks
c-----
        integer NL
        parameter(NL=1000)
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0
        real*8 dcda, dcdb, dcdr, dcdh
        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        real*4 vtp,dtp,rtp
c-----
c       local arguments
c-----
        real*8 ar, tm
        integer i

        ar = 6370.0d0
        tm=sqrt(1.+(c/(2.*ar*om))**2)
            do 20 i=1,mmax
                dcda(i)=dcda(i)*  vtp(i)/(tm**3)
                dcdb(i)=dcdb(i)*  vtp(i)/(tm**3)
                dcdh(i)=dcdh(i)*  dtp(i)/(tm**3)
                dcdr(i)=dcdr(i)*  rtp(i)/(tm**3)
   20       continue
c       write(6,*)'c flat=',c,' csph=',c/tm
        usph = ugr*tm
        csph = c/tm
        return
        end

        subroutine bldsph()
c-----
c       Transform spherical earth to flat earth
c
c       Schwab, F. A., and L. Knopoff (1972). 
c          Fast surface wave and free
c       mode computations, in  
c          Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c          B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  44, 45 , 41 pp 112-113
c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c     Revised 28 DEC 2007 to use mid-point, assume linear variation in
c     slowness instead of using average velocity for the layer
c     Use the Biswas (1972:PAGEOPH 96, 61-74, 1972) density mapping
c
c       This is for Rayleigh Waves
c
c-----
        implicit none
        integer NL
        parameter (NL=1000)
c-----
c       common blocks
c-----
        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        real*4 vtp,dtp,rtp
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
        integer mmax
c-----
c       local arguments
c-----
        real*8 ar, dr, r0, r1, z0, z1
        real*8 tmp
        integer i
c-----
c       vtp is the factor used to convert spherical 
c          velocity to flat velocity
c       rtp is the factor used to convert spherical 
c          density  to flat density 
c       dtp is the factor used to convert spherical 
c          boundary to flat boundary
c-----
c-----
c       duplicate computations of srwvds(IV)
c-----
        ar=6370.0d0
        dr=0.0d0
        r0=ar
        z0 = 0.0d+00
        zd(mmax)=1.0
        do 10 i=1,mmax
            r1 = r0 * dexp(-zd(i)/ar)
            if(i.lt.mmax)then
                z1 = z0 + zd(i)
            else
                z1 = z0 + 1.0d+00
            endif
            TMP=(ar+ar)/(r0+r1)
            vtp(i) = TMP
            rtp(i) = TMP**(-2.275)
            dtp(i) = sngl(ar/r0)
            r0 = r1
            z0 = z1
   10   continue
C        write(6,*)'vtp:',(vtp(i),i=1,mmax)
C        write(6,*)'rtp:',(rtp(i),i=1,mmax)
C        write(6,*)'dtp:',(dtp(i),i=1,mmax)
c-----
c       at this point the model information is no longer used and
c       will be overwritten
c-----
        return
        end

        function ffunc(nub,dm)
        implicit none
        complex*16 ffunc
        complex*16 nub
        real*8 dm
        complex*16 exqq
        complex*16 argcd
c-----
c       get the f function
c-----
        if(cdabs(nub).lt. 1.0d-08)then
            ffunc = dm
        else
            argcd = nub*dm
            if(dreal(argcd).lt.40.0)then
                  exqq = cdexp(-2.0d+00*argcd)
            else
                  exqq=0.0d+00
            endif
            ffunc = (1.d+00-exqq)/(2.d+00*nub)
        endif
        return
        end

        function gfunc(nub,dm)
        implicit none
        complex*16 gfunc
        complex*16 nub
        real*8 dm
        complex*16 argcd
        argcd = nub*dm
        if(dreal(argcd).lt.75)then
             gfunc = cdexp(-argcd)*dm
        else
             gfunc = dcmplx(0.0d+00, 0.0d+00)
        endif
        return
        end

        function h1func(nua,nub,dm)
        implicit none
        complex*16 h1func
        complex*16 nua, nub
        real*8 dm
        complex*16 argcd
        complex*16 exqq
        if(cdabs(nub+nua).lt. 1.0d-08)then
            h1func = dm
        else
            argcd = (nua+nub)*dm
            if(dreal(argcd).lt.40.0)then
                  exqq = cdexp(-argcd)
            else
                  exqq=0.0d+00
            endif
            h1func = (1.d+00-exqq)/(nub+nua)
        endif
        return
        end

        function h2func(nua,nub,dm)
        implicit none
        complex*16 h2func
        complex*16 nua, nub
        real*8 dm
        complex*16 argcd
        complex*16 exqq, exqp
        if(cdabs(nub-nua).lt. 1.0d-08)then
c-----
c           this should never occur for surface waves
c-----
            h2func =  dm
        else
            argcd = nua*dm
            if(dreal(argcd).lt.40.0)then
                  exqp = cdexp(-argcd)
            else
                  exqp=0.0d+00
            endif
            argcd = nub*dm
            if(dreal(argcd).lt.40.0)then
                  exqq = cdexp(-argcd)
            else
                  exqq=0.0d+00
            endif
            h2func = (exqq - exqp)/(nua-nub)
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
        dimension ee(nmat)
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

c-----
c      this section of code is specific to the isotropic problem
c      the TI code will use the same subroutine calls. The
c      conversion of isotropic to TI involves replacement of
c      subroutines and/or pieces of code, especially the common
c      block describing the model parameters
c-----

      subroutine svfunc(omega,wvno)
c
c     This combines the Haskell vector from sub down and
c     Dunkin vector from sub up to form the eigenfunctions.
c
      implicit none
c-----
c     command arguments
c-----
        real*8 omega, wvno
c-----
c     common blocks
c-----
        integer NL
        parameter (NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1    zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2    xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
        integer mmax
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     1                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0
        real*8 dcda, dcdb, dcdr, dcdh
        common/dunk/   cd(NL,5),exe(NL),exa(NL)
        real*8 cd
        real*8 exe, exa
        common/wateri/iwat(NL)
        integer iwat
        common/hask/   vv(NL,4)
        real*8 vv
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid

c-----
c       unknown common check for clean up
c-----
        common/hwat/wh(NL,2),hex(NL)
        real*8 wh, hex
c-----
c       internal variables
c-----
        real*8 fr
        integer i
        real*8 uu1, uu2, uu3, uu4, ext, f1213, fact
        real*8 cd1, cd2, cd3, cd4, cd5, cd6
        real*8 tz1, tz2, tz3, tz4
        integer jwat
c-----
c       get compound matrix from bottom to top
c-----
        call up(omega,wvno,fr)
        call down(omega,wvno)

c-----
c       get propagator from top to bottom 
c-----
        f1213 = -cd(1,2)
        ur(1) = cd(1,3)/cd(1,2)
        uz(1) = 1.0d+00
        tz(1) = 0.0d+00
        tr(1) = 0.0d+00
        uu0(1) = ur(1)
        uu0(2) = 1.0
        uu0(3) = fr
        uu0(4) = fr
c------
c      following
c                12       -1
c       Bk =    X|       Z
c                5-l,5-k  4l
c             -----------------
c                  12
c              - R |
c                  13
c-----
c       Compound Matrix (6x6)   Here (5x5)
c
c       X|12/12 ==      cd(1)
c       X|12/13 ==      cd(2)
c       X|12/14 ==      cd(3)
c       X|12/23 ==     -cd(3)
c       X|12/24 ==      cd(4)
c       X|12/34 ==      cd(5)
c       X|12/21 =  -X|12/12
c       X|12/ij =  -X|12/ji
c
c       where
c
c        |12
c       X|  = X|12/ij == X(1,i)*X(2,j) - X(1,j)*X(2,i)
c        |ij

c   and using the shorthand
c       True [11 12 13 14 15 16 ] = Reduced [ 11 12 13 -13 14 15 ]
c     we have
c-----
c   In long hand with full notation and using the complex matrix short hand, e.g.,
C   that the denominator is -R12 and recalling R11 == 0
c
c
c    |    |         |                               | |      |
c    |    |         |           12      12      12  | |   -1 |
c    | B1 |         |   0      X|     -X|      X|   | |  Z   |
c    |    |         |           34      24      23  | |   41 |
c    |    |         |                               | |      |
c    |    |         |   12              12      12  | |   -1 |
c    | B2 |         | -X|       0      X|     -X|   | |  Z   |
c    |    |     1   |   34              14      13  | |   42 |
c    |    | = _____ |                               | |      |
c    |    |      12 |   12      12              12  | |   -1 |
c    | B3 |   -R|   |  X|     -X|       0      X|   | |  Z   |
c    |    |      13 |   24      14              12  | |   43 |
c    |    |         |                               | |      |
c    |    |         |   12      12      12          | |   -1 |
c    | B4 |         | -X|      X|     -X|       0   | |  Z   |
c    |    |         |   23      13      12          | |   44 |
c
c
c-----
c    Here the vv are the first column of the  Haskell down, and the
c    cd are the first row if the X up. At the surface R=X 
c------
        ext = 0.0
        do 200 i=2,mmax
         cd1=  cd(i,1)
         cd2=  cd(i,2)
         cd3=  cd(i,3)
         cd4= -cd(i,3)
         cd5=  cd(i,4)
         cd6=  cd(i,5)
         tz1 = -vv(i,4)
         tz2 = -vv(i,3)
         tz3 =  vv(i,2)
         tz4 =  vv(i,1)


            uu1 =
     1                  tz2*cd6 - tz3*cd5 + tz4*cd4
            uu2 =
     1        -tz1*cd6          + tz3*cd3 - tz4*cd2
            uu3 =
     1         tz1*cd5 - tz2*cd3          + tz4*cd1
            uu4 =
     1        -tz1*cd4 + tz2*cd2 - tz3*cd1

            ext=exa(i)+exe(i) -exe(1)
            if(ext.gt.-80.0 .and. ext.lt.80.0 ) then
                fact=dexp(ext)

                ur(i)=uu1*fact/f1213
                uz(i)=uu2*fact/f1213
                tz(i)=uu3*fact/f1213
                tr(i)=uu4*fact/f1213
            else
                ur(i) = 0.0
                uz(i) = 0.0
                tz(i) = 0.0
                tr(i) = 0.0
            endif
  200   continue
c-----
c       correction for fluid layers on top if not
c       all fluid
c-----
        if(.not. allfluid)then
                jwat = 0
                do 300 i=1,mmax
                    if(iwat(i).gt.0)then
                            jwat = i
                    else
                            go to 301
                    endif
  300   continue
  301   continue
                if(jwat.gt.0)then
                        do i=1,jwat
                            ur(i) = 0.0
                            tr(i) = 0.0
                        enddo
                endif
        endif

CRBHc-----
CRBHc       Continue KLUDGE if top layers are water
CRBHc       The problem here was that even though the 
CRBHc          phase velocity is determined,
CRBHc       the computed eigenfunctions do not match perfectly,
CRBHc       and we do not match the free surface requirement that Tz=0 at
CRBHc       top of fluid. So RBH used the fluid propagator from top
CRBHc       down to the solid-fluid interface, the Haskell from
CRBHc       base up to the same interface, and then adjust the amplitudes
CRBHc       to match. The purpose is to get a better estimate of the
CRBHc       eigenfunction in the fluid.
CRBHc
CRBHc       Ultimately use Kennett
CRBHc-----
CRBHc       first get the largest exponent
CRBHc-----
CRBH        if(iwat(1).eq.1)then
CRBHc-----
CRBHc           OK first layer is fluid
CRBHc-----
CRBH            hmax = 0.0
CRBH            jwat=0
CRBH            do 300 i=1,mmax-1
CRBH                if(iwat(i).eq.1)then
CRBH                    jwat = i
CRBH                    hmax = hex(i+1)
CRBH                endif
CRBH  300       continue
CRBHc-----
CRBHc           now get correction factor
CRBHc-----
CRBH            corrd = sqrt(uz(jwat+1)**2 + tz(jwat+1)**2 )
CRBH            corrw = sqrt(wh(jwat+1,1)**2 + wh(jwat+1,2)**2 )
CRBH            if(corrd.eq.0.0)then
CRBHc-----
CRBHc       give precedence to surface values
CRBHc-----
CRBH                fac = 0.0
CRBH            else
CRBH                fac = corrw/corrd
CRBH            endif
CRBH            fac = fac * exp(hmax - hmax)
CRBH            fac = fac * 2.0
CRBHc-----
CRBHc           eventually normalize
CRBHc-----
CRBH        IF(.not. ALLFLUID)then
CRBH            do 400 i=1,mmax
CRBH                if(i.le.jwat)then
CRBH                    efac = exp(hex(i) - hmax)
CRBH                    efac = efac * 2.0
CRBH                    ur(i) = 0.0
CRBH                    uz(i) = wh(i,1)*efac
CRBH                    tz(i) = wh(i,2)*efac
CRBH                    tr(i) = 0.0
CRBH                else
CRBH                    ur(i) = ur(i) * fac
CRBH                    uz(i) = uz(i) * fac
CRBH                    tr(i) = tr(i) * fac
CRBH                    tz(i) = tz(i) * fac
CRBH                endif
CRBH  400       continue
CRBH        ENDIF
CRBH        endif
        
        return
        end

        subroutine up(omega,wvno,fr)
c-----
c       This finds the values of the Dunkin vectors at
c       each layer boundaries from bottom layer upward.
c
c       Note that these are the reduced vectors, e.g.,
c       The true compound (G An-1 ... A1 H) is
c       True [11 12 13 14 15 16 ] = Reduced [ 11 12 13 -13 14 15 ]
c-----
        implicit none
c-----
c       command line arguments
c-----
        real*8 omega, wvno, fr
c-----
c       common blocks
c-----
        integer NL
        parameter (NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
        integer mmax
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0
        real*8 dcda, dcdb, dcdr, dcdh
        common/dunk/   cd(NL,5),exe(NL),exa(NL)
        real*8 cd
        real*8 exe, exa
        common/wateri/iwat(NL)
        integer iwat
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid
c-----
c       unknown common check for clean up
c-----
        common/save/   dd(5,5),aa(4,4),ex1,ex2
        real*8 dd, aa, ex1, ex2
        common/aamatx/ ww(NL),xx(NL),yy(NL),zz(NL),
     1      cospp(NL),cosqq(NL)
        real*8 ww, xx, yy, zz, cospp, cosqq
        common/engerw/ wra,wd,wba
        real*8 wra,wd,wba
c-----
c       internal arguments
c-----

c-----
c       The save labeled common is used to pass the
c       Dunkin 5x5 and Haskell 4x4 matrices to the main program
c       exp(ex1) is the scaling for the Dunkin-Thrower compound matrix
c       exp(ex2) is the scaling for the Haskell matrices
c-----

        complex*16 rp, rsv, p, q
        real*8 wvno2, om2
        complex*16 gbr(2,5)
        real*8 cr, ee(5), exn
        real*8 ca(5,5)
        real*8 xka, xkb, pex, svex

        real*8 cosp,rsinp,sinpr,cossv,rsinsv,sinsvr
        real*8 exsum

        integer i,j,m,mmm1,nmat
c-----
c       initialize base vector
c-----
c-----
c       set up starting values for bottom halfspace
c-----
        wvno2=wvno*wvno
        om2 = omega*omega
        call evalg(0,mmax,mmax-1,gbr,1,
     1      wvno,omega,om2,wvno2)
c-----
c       note for surface waves, the ra, rb are real for the halfspace as long
c       as the phase velocity is less than the respective halfspace velocity
c-----

            cd(mmax,1)= dreal(gbr(1,1))
            cd(mmax,2)= dreal(gbr(1,2))
            cd(mmax,3)= dreal(gbr(1,3))
            cd(mmax,4)= dreal(gbr(1,4))
            cd(mmax,5)= dreal(gbr(1,5))
            exe(mmax) = 0.0d+00
c------
c       matrix multiplication from bottom layer upward
c------
        mmm1=mmax-1
        exsum = 0.0
        do 500 m = mmm1,1,-1
                xka = omega/za(m)
                if(zb(m).gt.0.0)then
                  xkb = omega/zb(m)
                else
                  xkb = 0.0
                endif
                rp =CDSQRT(dcmplx(wvno2-xka*xka,0.0d+00))
                rsv=CDSQRT(dcmplx(wvno2-xkb*xkb,0.0d+00))
                p = rp  * zd(m)
                q = rsv * zd(m)
                call varsv(p,q,rp, rsv,
     1              cosp, cossv, rsinp, rsinsv,
     1              sinpr, sinsvr, pex,svex,iwat(m))
                call dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,
     1              sngl(zrho(m)),sngl(zb(m)),iwat(m),pex,pex+svex,
     2              wvno,wvno2,om2)

            nmat = 5
            do 200 i=1,nmat
                cr=0.0d+00
                do 100 j=1,nmat
                    cr=cr+cd(m+1,j)*ca(j,i)
  100           continue
                ee(i)=cr
  200       continue
            exn= 0.0d+00
            call normc(ee,exn,nmat)
            exsum = exsum + pex + svex + exn
            exe(m) = exsum
            do 300 i = 1,nmat
                cd(m,i)=ee(i)
  300       continue
  500   continue
c-----
c       define period equation
c-----
        fr=cd(1,1)
        return
        end

        subroutine dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,
     1              rho,b,iwat,ex,exa,wvno,wvno2,om2)
        implicit none
        real*8 ca(5,5), wvno, wvno2, om2
        real*4 rho,b
        real*8 cosp,rsinp,sinpr,cossv,rsinsv,sinsvr
        integer iwat

        real rho2
        real*8 gam
        real*8 ex, exa
        real*8 a0
        real*8 cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz

        real*8 gam2,gamm1,gamm2,a0c,xz2,wy2,temp
        real*8 dfac
        real*8 cqww2, cqxw2, g1wy2, gxz2, g2wy2, g2xz2
        real*8 gg1, a0cgg1
        integer i, j
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
c       if we define a new matrices
c       related to the original matrices by
c-----
c         A11     A12     A13         A15     A16
c         A21     A22     A23         A25     A15
c        2 A31   2 A32   2 A33 -1   -2 A23  -2 A13
c         A51     A52    -A32         A22     A12
c         A61     A51    -A31         A21     A11
c-----
c
c       [ G11   G12  G13    G15 G16  ]
c       [ H11   H21 2 H31   H51 H61  ] ^T
c
c-----
c       this means that some of the original definitions of the 
c       Aij elements must be changed for the
c       definition of the modified 5x5 compount A matrix
c
c       old 6x6                 new 5x5
c       A11 = 1 - A33 + A22     1 - (1/2)(new A33 + 1) + new A2)
c       A53 = -A32              A43 = - (1/2) new A32
c       A63 = -A31              A53 = - (1/2) new A31
c-----
c       To recover the needed elements, 
c          we note that the old G14 = -old G14 = new G13
c-----
c-----
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,5
                do 101 i=1,5
                    ca(i,j) = 0.0d+00
  101           continue
  100       continue
            if(ex.gt.35.0d+00)then
                dfac = 0.0d+00
            else
                dfac = dexp(-ex)
            endif
            ca(3,3) = dfac
            ca(1,1) = cosp
            ca(5,5) = cosp
            ca(1,2) = - rsinp/(rho*om2)
            ca(2,1) = - rho*sinpr*om2
            ca(2,2) = cosp
            ca(4,4) = cosp
            ca(4,5) = ca(1,2)
            ca(5,4) = ca(2,1)
        else
            if ( exa.lt. 60.0)then
                a0=dexp(-exa)
            else
                a0 = 0.0d+00
            endif
            cpcq = cosp * cossv
            cpy  = cosp*sinsvr
            cpz  = cosp*rsinsv
            cqw  = cossv*sinpr
            cqx  = cossv*rsinp
            xy   = rsinp*sinsvr
            xz   = rsinp*rsinsv
            wy   = sinpr*sinsvr  
            wz   = sinpr*rsinsv
c-----
c       elastic layer
c-----
            rho2= rho*rho
            gam = 2.0*b*b*wvno2/om2
            gam2  = gam*gam
            gamm1 = gam-1.
            gamm2 = gamm1*gamm1
            cqww2 = cqw * wvno2
            cqxw2 = cqx / wvno2
            gg1 = gam*gamm1
            a0c  = dcmplx(2.0d+00,0.0d+00)*
     1          (dcmplx(a0,0.0d+00)-cpcq)
            xz2  = xz/wvno2
            gxz2 = gam*xz2
            g2xz2 = gam2 * xz2
            a0cgg1 = a0c*(gam+gamm1)
            wy2  = wy*wvno2
            g2wy2 = gamm2 * wy2
            g1wy2 = gamm1 * wy2

c-----
c       OK by symmetry
c----
            temp = a0c*gg1 + g2xz2 + g2wy2
            ca(3,3) = a0 + temp + temp
            ca(1,1) = cpcq-temp
            ca(1,2) = (-cqx + wvno2*cpy)/(rho*om2)
            temp = dcmplx(0.5d+00,0.0d+00)*a0cgg1 + gxz2 + g1wy2
            ca(1,3) = wvno*temp/(rho*om2)

            ca(1,4) = (-cqww2+cpz)/(rho*om2)
            temp = wvno2*(a0c + wy2) + xz
            ca(1,5) = -temp/(rho2*om2*om2)

            ca(2,1) = (-gamm2*cqw + gam2*cpz/wvno2)*rho*om2
            ca(2,2) = cpcq
            ca(2,3) = (gamm1*cqww2 - gam*cpz)/wvno
            ca(2,4) = -wz
            ca(2,5)=ca(1,4)


            temp =dcmplx(0.5d+00,0.0d+00)*a0cgg1*gg1 
     1          + gam2*gxz2 + gamm2*g1wy2
            ca(3,1) = -dcmplx(2.0d+00,0.0d+00)*temp*rho*om2/wvno
            ca(3,2) = -wvno*(gam*cqxw2 - gamm1*cpy)*
     1          dcmplx(2.0d+00,0.0d+00)

            ca(3,4)=-2.0d+00*ca(2,3)
            ca(3,5)=-2.0d+00*ca(1,3)

            ca(4,1) = (-gam2*cqxw2 + gamm2*cpy)*rho*om2
            ca(4,2) = -xy
            ca(4,3)= -ca(3,2)/2.0d+00
            ca(4,4)=ca(2,2)
            ca(4,5)=ca(1,2)

            temp = gamm2*(a0c*gam2 + g2wy2) + gam2*g2xz2
            ca(5,1) = -rho2*om2*om2*temp/wvno2
            ca(5,2)=ca(4,1)
            ca(5,3)=-ca(3,1)/2.0d+00
            ca(5,4)=ca(2,1)
            ca(5,5)=ca(1,1)
        endif
        return
        end

        subroutine evalg(jbdry,m,m1,gbr,inp,
     1      wvno,om,om2,wvno2)
c-----
c       this is from hspec96, but not everything is required
c       tot he layered halfspace prroblem
c-----
        implicit none
c-----
c       procedure arguments
c-----
        integer jbdry, m, m1, inp
        complex*16 gbr(2,5)
        real*8 wvno, wvno2, om, om2
c-----
c       common blocks
c-----
        integer NL
        parameter(NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1    zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2    xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas,
     1      zfrefp, zfrefs,
     1      xmu, xlam
        integer mmax
        common/modspec/allfluid
        logical allfluid
        common/emat/e,einv,ra,rb
        complex*16 ra,rb
        complex*16 e(4,4), einv(4,4)
        common/wateri/iwat(NL)
        integer iwat

c-----
c       internal parameters
c-----
        real*8 xka,xkb,gam,gamm1
        complex*16 CDSQRT


        integer i,j
        
c-----
c       set up halfspace conditions
c-----
        xka = om/za(m)
        if(zb(m).gt.0.0)then
          xkb = om/zb(m)
        else
          xkb = 0.0
        endif
        ra=CDSQRT(dcmplx(wvno2-xka*xka,0.0d+00))
        rb=CDSQRT(dcmplx(wvno2-xkb*xkb,0.0d+00))
        gam = dble(zb(m))*wvno/om
        gam = 2.0d+00 * (gam * gam)
        gamm1 = gam - dcmplx(1.0d+00,0.0d+00)
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
            if(zb(m) .gt. 0.0)then
c-----
c               ELASTIC ABOVE - RIGID
c-----
                gbr(inp,1) = dcmplx(1.0d+00,0.0d+00)
                gbr(inp,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,5) = dcmplx(0.0d+00,0.0d+00)
            else
c-----
c               FLUID ABOVE - RIGID
c-----
                gbr(inp,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,5) = dcmplx(0.0d+00,0.0d+00)
                if(allfluid)then
                    gbr(inp,1) = dcmplx(1.0d+00,0.0d+00)
                else
                    gbr(inp,4) = dcmplx(1.0d+00,0.0d+00)
                endif
c-----
c               (pseudo SH)
c-----
            endif
        else if(jbdry.eq.0)then
c-----
c       HALFSPACE
c-----
            if(iwat(m).eq.0)then
c-----
c               ELASTIC HALFSPACE
c-----
c       multiply G of Herrmann 2001 by - rho^2 om^4 k^2 ra rb
c       should have no effect since it is in both numerator and
c       denominator -- however will not give the correct potential
c       coefficient -- so rethink?
c-----
            E(1,1)= wvno
            E(1,2)= rb
            E(1,3)= wvno
            E(1,4)= -rb

            E(2,1)= ra
            E(2,2)= wvno
            E(2,3)= -ra
            E(2,4)= wvno

            E(3,1)= zrho(m)*om2*gamm1
            E(3,2)= zrho(m)*om2*gam*rb/wvno
            E(3,3)= zrho(m)*om2*gamm1
            E(3,4)= -zrho(m)*om2*gam*rb/wvno

            E(4,1)= zrho(m)*om2*gam*ra/wvno
            E(4,2)= zrho(m)*om2*gamm1
            E(4,3)= -zrho(m)*om2*gam*ra/wvno
            E(4,4)= zrho(m)*om2*gamm1

            EINV(1,1)= 0.5*gam/wvno
            EINV(1,2)= -0.5*gamm1/ra
            EINV(1,3)= -0.5/(zrho(m)*om2)
            EINV(1,4)= 0.5*wvno/(zrho(m)*om2*ra)

            EINV(2,1)= -0.5*gamm1/rb
            EINV(2,2)= 0.5*gam/wvno
            EINV(2,3)= 0.5*wvno/(zrho(m)*om2*rb)
            EINV(2,4)= -0.5/(zrho(m)*om2)

            EINV(3,1)= 0.5*gam/wvno
            EINV(3,2)=  0.5*gamm1/ra
            EINV(3,3)= -0.5/(zrho(m)*om2)
            EINV(3,4)= -0.5*wvno/(zrho(m)*om2*ra)

            EINV(4,1)= 0.5*gamm1/rb
            EINV(4,2)= 0.5*gam/wvno
            EINV(4,3)= -0.5*wvno/(zrho(m)*om2*rb)
            EINV(4,4)= -0.5/(zrho(m)*om2)

                gbr(inp,1)=dble(zrho(m)*zrho(m))*om2*om2*
     1              (-gam*gam*ra*rb+wvno2*gamm1*gamm1)
                gbr(inp,2)=-dble(zrho(m))*(wvno2*ra)*om2
                gbr(inp,3)=-dble(zrho(m))*(-gam*ra*rb+wvno2*gamm1)
     1              *om2*wvno
                gbr(inp,4)=dble(zrho(m))*(wvno2*rb)*om2
                gbr(inp,5)=wvno2*(wvno2-ra*rb)
             gbr(inp,1)=0.25*gbr(inp,1)/
     1               (-zrho(m)*zrho(m)*om2*om2*wvno2*ra*rb)
             gbr(inp,2)=0.25*gbr(inp,2)/
     1               (-zrho(m)*zrho(m)*om2*om2*wvno2*ra*rb)
             gbr(inp,3)=0.25*gbr(inp,3)/
     1               (-zrho(m)*zrho(m)*om2*om2*wvno2*ra*rb)
             gbr(inp,4)=0.25*gbr(inp,4)/
     1               (-zrho(m)*zrho(m)*om2*om2*wvno2*ra*rb)
             gbr(inp,5)=0.25*gbr(inp,5)/
     1               (-zrho(m)*zrho(m)*om2*om2*wvno2*ra*rb)

            else if(iwat(m).eq.1)then
c-----
c               FLUID HALFSPACE
c-----
                if(allfluid)then
                    gbr(inp,1) = dble(0.5) / ra
                    gbr(inp,2) = dcmplx(0.5d+00,0.0d+00)/
     1                  (-dble(zrho(m))*om2)
                    gbr(inp,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(inp,4) = dcmplx(0.0d+00,0.0d+00)
                    gbr(inp,5) = dcmplx(0.0d+00,0.0d+00)
                else
                    gbr(inp,1) = dcmplx(0.0d+00,0.0d+00)
                    gbr(inp,2) = dcmplx(0.0d+00,0.0d+00)
                    gbr(inp,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(inp,4) = dble(0.5*zrho(m)*om2) / ra
                    gbr(inp,5) = dcmplx(-0.5d+00,0.0d+00)
                endif
c-----
c           for safety null the matrices and then fill with a 2x2
c-----
            do i=1,4
               do j=1,4
                  E(i,j) = dcmplx(0.0d+00, 0.0d+00)
                  EINV(i,j) = dcmplx(0.0d+00, 0.0d+00)
               enddo
            enddo
            E(1,1)=  ra
            E(1,2)= -ra
            E(2,1)= -zrho(m)*om2
            E(2,2)= -zrho(m)*om2

            EINV(1,1)= 0.5/ra
            EINV(1,2)= -0.5/(zrho(m)*om2)
            EINV(2,1)= -0.5/ra
            EINV(2,2)= -0.5/(zrho(m)*om2)
            endif
        else if(jbdry.eq.1)then
c-----
c       FREE - check properties of layer above
c-----
            if(zb(m) .gt. 0.0)then
                gbr(inp,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,5) = dcmplx(1.0d+00,0.0d+00)
                
            else
                gbr(inp,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(inp,5) = dcmplx(0.0d+00,0.0d+00)
                if(allfluid)then
                    gbr(inp,2) = dcmplx(1.0d+00,0.0d+00)
                else
                    gbr(inp,5) = dcmplx(1.0d+00,0.0d+00)
                endif
            endif
        endif
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
c
c       for pure real or imaginary rp, and rsv, the
c       functions returned are pure real
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
        COMPLEX*16 epp, epm, eqp, eqm
        COMPLEX*16 sinp, sinq

        REAL*8 PFAC, SVFAC
        
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
            if(pr.lt.30.) then
                pfac=dexp(-2.*pr)
            else
                pfac  = 0.0d+00
            endif
            cosp = dreal(epp + pfac*epm)
            sinp = epp - pfac*epm
            rsinp = dreal(rp *sinp)
            sinpr = dreal(sinp/rp)
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
            if(pr.lt.30.) then
                pfac=dexp(-2.*pr)
            else
                pfac  = 0.0d+00
            endif
            cosp = dreal(epp + pfac*epm)
            sinp = epp - pfac*epm
            rsinp = dreal(rp *sinp)
            sinpr = dreal(sinp/rp)

            if(qr.lt.30.) then
                svfac=dexp(-2.*qr)
            else
                svfac  = 0.0d+00
            endif
            cosq = dreal(eqp + svfac*eqm)
            sinq = eqp - svfac*eqm
            rsinq = dreal(rsv*sinq)
            sinqr = dreal(sinq/rsv)

        endif
        return
        end

        subroutine hska(AA,cosp,rsinp,sinpr,tcossv,trsinsv,tsinsvr,
     2      rho,b,iwat,pex,svex,wvno,wvno2,om2)
c-----
c       Changes
c
c-----
        implicit none
c-----
c       subroutine variables
c-----
        real*8 AA(4,4)
        real*8 cosp , tcossv 
        real*8 rsinp, trsinsv
        real*8 sinpr, tsinsvr
        real rho,b
        integer iwat
        real*8 pex, svex
        real*8  om2, wvno, wvno2

c-----
c       local variables
c-----
        real*8 gam, gamm1
        real*8 cossv 
        real*8 rsinsv
        real*8 sinsvr
        real*8 dfac

        integer i, j

        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,4
                do 101 i=1,4
                    AA(i,j) = 0.0d+00
  101           continue
  100       continue
            if(pex.gt.35.0d+00)then
                dfac = 0.0d+00
            else
                dfac = dexp(-pex)
            endif
            AA(1,1) = dfac
            AA(4,4) = dfac
            AA(2,2) = cosp
            AA(3,3) = cosp
            AA(2,3) = -rsinp/(rho*om2)
            AA(3,2) = - rho*om2*sinpr
        else
c-----
c       elastic layer
c-----
c-----
c           adjust the definitions of the SV trig functions
c           Note that the exp(pex) and exp(svex) have
c           been factored from the these functions, so that
c           cosp     = exp(pex) ( 1 + exp(-2 pex)/2 = exp(pex) cosp
c               True                                               local
c           cossv     = exp(svex) ( 1 + exp(-2 svex)/2 = exp(svex) cossv
c               True                                                    local
c           This separation is fine for the combound matrices which have
c           terms such as cosp cosq, however for the Haskell propagator we need
c           cosp - cossv, or since pex >= svex
c
c           (cosp - cossv)    = exp(pex) [ cosp      - exp(svex-pex) cossv    ]
c                       true                  local                       local
c           It is for this reason that we must modify the tcossv trsinsv rsinsvr in
c           the command line 
c-----
            if( (pex-svex) .gt. 70.0)then
                dfac = 0.0d+00
            else
                dfac = dexp(svex-pex)
            endif
            cossv = dfac * tcossv
            rsinsv = dfac * trsinsv
            sinsvr = dfac * tsinsvr
                
            gam = 2.0*b*b*wvno2/om2
            gamm1 = gam -1.0
C            AA(1,1) =  gam*cosp - gamm1*cossv
            AA(1,1) =  cossv + gam*(cosp - cossv)
            AA(1,2) =  -wvno*gamm1*sinpr + gam*rsinsv/wvno
            AA(1,3) =  -wvno*(cosp-cossv)/(rho*om2)
            AA(1,4) =   (wvno2*sinpr - rsinsv)/(rho*om2)
            AA(2,1) =   gam*rsinp/wvno - wvno*gamm1*sinsvr
            AA(2,2) =   cosp - gam*(cosp- cossv)
            AA(2,3) =   ( -rsinp + wvno2*sinsvr)/(rho*om2)
            AA(2,4) = - AA(1,3)
            AA(3,1) =   rho*om2*gam*gamm1*(cosp-cossv)/wvno
            AA(3,2) = rho*om2*(-gamm1*gamm1*sinpr+gam*gam*rsinsv/wvno2)
            AA(3,3) =   AA(2,2)
            AA(3,4) = - AA(1,2)
            AA(4,1) =   rho*om2*(gam*gam*rsinp/wvno2-gamm1*gamm1*sinsvr)
            AA(4,2) = - AA(3,1)
            AA(4,3) = - AA(2,1)
            AA(4,4) =   AA(1,1)
        endif
        return
        end

        subroutine down(omega,wvno)
c-----
c       This finds the values of the Haskell vectors at
c       each layer boundaries from top layer downward.
c-----
        implicit none
c-----
c     command arguments
c-----
        real*8 omega, wvno
c-----
c     common blocks
c-----
        integer NL
        parameter (NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1    zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2    xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas,
     1      zfrefp, zfrefs,
     1      xmu, xlam
        integer mmax
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     1                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0
        real*8 dcda, dcdb, dcdr, dcdh
        common/dunk/   cd(NL,5),exe(NL),exa(NL)
        real*8 cd
        real*8 exe, exa
        common/wateri/iwat(NL)
        integer iwat
        common/hask/   vv(NL,4)
        real*8 vv
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid
c-----
c      internal variables
c-----
       real*8 wvno2, om2
       integer m
       real*8 aa(4,4)
       integer i,j
       real*8 xka, xkb, pex, svex
       real*8 ex2
       complex*16 rp, rsv, p, q
       real*8 cosp,rsinp,sinpr,cossv,rsinsv,sinsvr
       real*8 cc, aa0(4)
       real*8 exsum

       om2 = omega*omega
       wvno2 = wvno*wvno
c-----
c      propagate down to get the first column of the
c      Am .. A2 A1
c-----
c-----
c      initialize the top surface for the
c      first column of the Haskell propagator
c-----
      do i=1,4
         if(i.eq.1)then
             vv(1,i) = 1.0d+00
         else
             vv(1,i) = 0.0d+00
         endif
      enddo
      exa(1) = 0.0
      
c------
c       matrix multiplication from top layer downward
c------
       exsum  = 0.0
       do 500 m= 1, mmax -1
                xka = omega/za(m)
                if(zb(m).gt.0.0)then
                  xkb = omega/zb(m)
                else
                  xkb = 0.0
                endif
                rp =CDSQRT(dcmplx(wvno2-xka*xka,0.0d+00))
                rsv=CDSQRT(dcmplx(wvno2-xkb*xkb,0.0d+00))
                p = rp  * zd(m)
                q = rsv * zd(m)
                call varsv(p,q,rp, rsv,
     1              cosp, cossv, rsinp, rsinsv,
     1              sinpr, sinsvr, pex,svex,iwat(m))

                call hska(AA,cosp,rsinp,sinpr,
     1                cossv,rsinsv,sinsvr,
     2                sngl(zrho(m)),sngl(zb(m)),iwat(m),
     3                pex,svex,wvno,wvno2,om2)
            do 300 i=1,4
                cc=0.0d+00
                do 200 j=1,4
                    cc=cc+aa(i,j)*vv(m,j)
  200           continue
                aa0(i)=cc
  300       continue
            ex2 = 0.0
            call normc(aa0,ex2,4)
            exsum = exsum + pex + ex2
            exa(m+1)= exsum
            
            do 400 i=1,4
                vv(m+1,i)=aa0(i)
  400       continue
  500   continue
c-----
c       vv is  column 1 of the Haskell propagator at the top of layer m
c-----
        return
        end


        subroutine gammap(omega,wvno,gammar)
c-----
c       This routine finds the attenuation gamma value.
c
c-----
        implicit none
c-----
c       routine arguments
c-----
        real*8 omega, wvno, gammar
c-----
c       common blocks
c-----
        integer NL
        parameter (NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
        integer mmax
        common/wateri/iwat(NL)
        integer iwat
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     1                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0
        real*8 dcda, dcdb, dcdr, dcdh
c-----
c       internal variables
c-----
        real*8 dc, pi, x, omgref, c
        integer i

        gammar=0.0
        dc = 0.0
        pi = 3.141592653589493d+00
        do 100 i=1,mmax
            if(iwat(i).eq.  0)then
                x=dcdb(i)*zb(i)*zqbi(i)
                gammar = gammar + x
                omgref=2.0*pi*zfrefs(i)
                dc = dc + dlog(omega/omgref)*x/pi
            endif
            x=dcda(i)*za(i)*zqai(i)
            gammar = gammar + x
            omgref=2.0*pi*zfrefp(i)
            dc = dc + dlog(omega/omgref)*x/pi
  100   continue
        c=omega/wvno
        gammar=0.5*wvno*gammar/c
        c=c+dc
        wvno=omega/c
        return
        end

        subroutine energy(om,wvno,mmax)
c-----
c       determine the energy integrals and also the
c       phase velocity partials
c-----
        implicit none
c-----
c       subroutine arguments
c-----
        real*8 om, wvno
        integer mmax
c-----
c       common blocks
c-----
        integer NL
        parameter (NL=1000)
        common/sumi/   sumi0,sumi1,sumi2,sumi3,flagr,are,ugr
        real*8 sumi0, sumi1, sumi2, sumi3, flagr, are, ugr
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     1                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0
        real*8 dcda, dcdb, dcdr, dcdh
        common/wateri/iwat(NL)
        integer iwat
        
c-----
c       external function call prototypes
c-----
        real*8 intijr
c-----
c       internal variables
c-----
        real*8 om2, wvno2
        real*8 TA, TC, TF, TL, TN
        real*8 INT11, INT13, INT22, INT44, INT33, INT24
        real*8 URUR, UZUZ, DURDUR, DUZDUZ, URDUZ, UZDUR
        real*8 c
        real*8 fac
        real*8 facah, facav, facbh, facbv,  facr
c       real*8 face

        integer m
c-----
c       coefficients of the ODE
c-----
        real*8 a12, a14, a21, a23
        real *8 ah, av, bh, bv, eta, rho

        integer TYPELYR
c-----
c       TYPELYR = -1  layer represents upper halfspace
c                  0  layer is true intenal layer
c                 +1  layer represents lower halfspace
c-----

c----
c       initialize
c----

        sumi0 = 0.0
        sumi1 = 0.0
        sumi2 = 0.0
        sumi3 = 0.0
        c = om/wvno
        om2 = om*om
        wvno2 = wvno*wvno

        do m=1,mmax
           call getmat(m,wvno,om,a12, a14, a21, a23,
     1      ah,av,bh,bv,eta,rho,
     2      TA,TC,TF,TL,TN,iwat(m))
c------
c       get the integrals for layers over a halfspace
c       this is here is we ever adopt the code to the coal 
c       seam problem
c-----
           if(m.eq.mmax)then
                   typelyr = 1
           else
                   typelyr = 0
           endif
           INT11 = intijr(1,1,m,typelyr,om,om2,wvno,wvno2)
           INT13 = intijr(1,3,m,typelyr,om,om2,wvno,wvno2)
           INT22 = intijr(2,2,m,typelyr,om,om2,wvno,wvno2)
           INT24 = intijr(2,4,m,typelyr,om,om2,wvno,wvno2)
           INT33 = intijr(3,3,m,typelyr,om,om2,wvno,wvno2)
           INT44 = intijr(4,4,m,typelyr,om,om2,wvno,wvno2)


           if(iwat(m).eq.1)then
c-----
c       fluid - not these are for the 4x4 formulation
c-----
             URUR   = INT22*(wvno/(rho*om2))**2
             UZUZ   = INT11
             URDUZ  =  - (wvno/(rho*om2))*a12*INT22
             DUZDUZ = a12*a12*INT22
             sumi0  = sumi0 + rho*(URUR + UZUZ)
             sumi1  = sumi1 + TA*URUR
             sumi2  = sumi2 - TF*URDUZ
             sumi3  = sumi3 + TC*DUZDUZ
             facah = rho*ah*(URUR -2.*eta*URDUZ/wvno)
             facav = rho*av*DUZDUZ / wvno2
             dcda(m) = facah + facav
             facr = -0.5*c*c*(URUR + UZUZ)
             dcdr(m) = 0.5*(av*facav + ah*facah )/rho + facr
c-----
c            define the ur in the fliud from the Tz - this will be at
c            the top of the layer
c-----
             ur(m) = - wvno*tz(m)/(rho*om2)

           else
c-----
c       solid
c-----
             URUR   = INT11
             UZUZ   = INT22
             DURDUR = a12*a12*INT22 + 2.*a12*a14*INT24 + a14*a14*INT44
             DUZDUZ = a21*a21*INT11 + 2.*a21*a23*INT13 + a23*a23*INT33
             URDUZ  = a21*INT11 + a23*INT13
             UZDUR  = a12*INT22 + a14*INT24
             sumi0  = sumi0 + rho*(URUR + UZUZ)
             sumi1  = sumi1 + TL*UZUZ + TA*URUR
             sumi2  = sumi2 + TL*UZDUR - TF*URDUZ
             sumi3  = sumi3 + TL*DURDUR + TC*DUZDUZ

c-----
c            partial derivatives of phase velocity with
c            respect to medium parameters. Note that these are
c            later divided by (U sumi0) when these are finalized
c
c            note that the code distainguishes between
c            ah and av, bh and bv and uses eta. These are
c            actually TI parameters, and for isotropic media
c            ah = av, bh = bv aned eta = 1. The code is written this
c            way for easier conversion to a TI case
c-----
             facah = rho*ah*(URUR -2.*eta*URDUZ/wvno)
             facav = rho*av*DUZDUZ / wvno2
             facbh = 0.0
             facbv = rho*bv*(UZUZ + 2.*UZDUR/wvno + DURDUR/wvno2 +
     1           4.*eta*URDUZ/wvno )
             dcda(m) = facah + facav
             dcdb(m) = facbv + facbh
C       WRITE(6,*)'ah,av,bh,bv,eta:',ah,av,bh,bv,eta
C       WRITE(6,*)'facah:',facah
C       WRITE(6,*)'facav:',facav
C       WRITE(6,*)'facbh:',facbh
C       WRITE(6,*)'facbv:',facbv
C       WRITE(6,*)'m,dcda(m),dcdb(m):',m,dcda(m),dcdb(m)

c            face = - TF*URDUZ/(wvno*eta)
c            dcdah(m) = facah
c            dcdav(m) = facav
c            dcdbh(m) = facbh
c            dcdbv(m) = facbv
c            dcdeta(m) = faceta


             facr = -0.5*c*c*(URUR + UZUZ)
c-----
c      this is correct for TI
             dcdr(m) = 0.5*(av*facav + ah*facah + bv*facbv)/rho + facr

c-----
c       partial with layer thickness needs the value at the layer, e.g.,
c       dcdh(1) is top surface
c-----
           endif
        enddo
c-----
c       determine final parameters
c-----

        flagr=om2*sumi0-wvno2*sumi1-2.d+00*wvno*sumi2-sumi3
        ugr=(wvno*sumi1+sumi2)/(om*sumi0)
        are=wvno/(2.d+00*om*ugr*sumi0)
        fac = are*c/wvno2
c-----
c       use the final factors to apply the 1/(U I0) factor
c----
        do m=1,mmax
           dcda(m) = dcda(m) /(ugr*sumi0)
           dcdb(m) = dcdb(m) /(ugr*sumi0)
           dcdr(m) = dcdr(m) /(ugr*sumi0)
        enddo
c-----
c       determine the dcdh
c-----
        call getdcdh(om2,wvno,wvno2,fac)
        return
        end

        subroutine getmat(m,wvno,om,a12, a14, a21, a23,
     1         ah,av,bh,bv,eta,rho,
     2         TA,TC,TF,TL,TN,iwat)
        
c-----
c       get matrix elements of the ODE
c       We do not require all eight non-zero elements
c       also return the model characterization parameters
c-----
        implicit none
c-----
c       procedure arguments
c-----
        real*8 wvno,om,a12, a14, a21, a23
        real*8 ah,av,bh,bv,eta,rho
        real*8 TA,TC,TF,TL,TN
        integer m, iwat
c-----
c       common blocks
c-----

        integer NL
        parameter (NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
        integer mmax

        if(iwat.eq.1)then
c-----
c            non-gravitating fluid - using 2x2 ODE
c-----
             ah= za(m)
             av= za(m)
             bh= 0.0
             bv= 0.0
             rho=zrho(m)
             eta = 1.0
     
             TL = 0.0
             TN = 0.0
             TC = zrho(m)*za(m)*za(m)
             TA = zrho(m)*za(m)*za(m)
             TF = TA - 2.*TN
     
             a12 = - ( wvno*wvno - om*om/(ah*ah))/(rho*om*om)
        else
c-----
c            elastic - using 4x4 ODE
c-----
             ah= za(m)
             av= za(m)
             bh= zb(m)
             bv= zb(m)
             rho=zrho(m)
             eta = 1.0
     
     
             TL = zrho(m)*zb(m)*zb(m)
             TN = zrho(m)*zb(m)*zb(m)
             TC = zrho(m)*za(m)*za(m)
             TA = zrho(m)*za(m)*za(m)
             TF = TA - 2.*TN
     
             a12 = -wvno
             a14 = 1.0/TL
             a21 = wvno * TF/TC
             a23 = 1.0/TC
        endif

        return
        end

        function intijr(i,j,m,typelyr,om,om2,wvno,wvno2)
        implicit none
c-----
c       procedure arguments
c-----
        real*8 intijr
        integer i,j,m
        real*8 om, om2, wvno,wvno2
        integer TYPELYR
c-----
c       TYPELYR = -1  layer represents upper halfspace
c                  0  layer is true intenal layer
c                 +1  layer represents lower halfspace
c-----
c NOTE DO NOT PERMIT RA RB to be ZERO add a small number TAKEN CARE OF
c IN FUNC CALLS
c beware of E matrix for fluid - internally I use a 4x4
c but hte book uses a 2x2 in the 8-10 problem
c need typelyr, wvno, om, om2, wvno2, m, mmax, do not need medium
c-----
c       Potential coefficients
c       kmpu, kmsu are the upward coefficients at bottom of layer
c       km1pd, km1sd are the downward coefficients at top of layer
c-----
        integer NL
        parameter (NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
        integer mmax
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     1                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0
        real*8 dcda, dcdb, dcdr, dcdh
        common/wateri/iwat(NL)
        integer iwat
        common/emat/e,einv,ra,rb
        complex*16 e(4,4), einv(4,4)
        complex*16 ra, rb

c-----
c       external function definitions
c-----
        complex*16 ffunc, gfunc, h1func, h2func
c-----
c       internal variables
c-----
        complex*16 cintijr, FA, GA, FB, GB, H1, H2
        complex*16 gbr(2,5)
        complex*16 kmpu, kmsu, km1pd, km1sd


c-----
c       call evalg to get the E and EINV matrices
c       evalg knows about water - we ignore everything else
c-----
        call evalg(0,m,m-1,gbr,1,
     1      wvno,om,om2,wvno2)
c-----
c       for an elastic solid
c
c                       T        PU  SvU  PD  SvD T
c       [Ur, Uz, Tr, Tz]  = E [ K  ,K   ,K  ,K   ]
c and
c          PU  SvU  PD  SvD T   -1                 T
c       [ K  ,K   ,K  ,K   ] = E   [Ur, Uz, Tr, Tz]
c
c       The text used the notation K m-1 to represent the potentials
c       at the top of a layer and K m those at the bottom of the layer
c       Because of Fortran indexing UZ(1) is the Z displacement at the
c       top of layer 1 and UZ(2) is the value at the bottom.  This is
c       reason for the slight incondistency in presentation below
c-----       

        if(iwat(m).eq.1)then
c-----
c       fluid
c-----
            if(typelyr .lt. 0)then
              kmpu = einv(1,1)*uz(m+1)
     1             + einv(1,2)*tz(m+1) 
              cintijr = e(i,1)*e(j,1)*kmpu*kmpu/(2.0*ra)
            else if(typelyr .eq. 0)then
              km1pd = einv(2,1)*uz(m)
     1              + einv(2,2)*tz(m) 
              kmpu = einv(1,1)*uz(m+1)
     1             + einv(1,2)*tz(m+1)  
              FA=ffunc(ra,zd(m))
              GA=gfunc(ra,zd(m))
              cintijr = e(i,1)*e(j,1)*kmpu*kmpu*FA
     1          +(e(i,1)*e(j,2)+e(i,2)*e(j,1))*kmpu*km1pd*GA       
     2          + e(i,2)*e(j,2)*km1pd*km1pd*FA
            else
              km1pd = einv(2,1)*uz(m)
     2              + einv(2,2)*tz(m)
              cintijr = e(i,2)*e(j,2)*km1pd*km1pd/(2.0*ra)
            endif
        else
c-----
c       solid
c-----
            if(typelyr .lt. 0)then
c-----
c             potential coefficients for upper halfspace
c             based on the displacement, stress at top
c             Basically m = 1
c-----
              kmpu = einv(1,1)*ur(m) + einv(1,2)*uz(m)
     2             + einv(1,3)*tz(m) + einv(1,4)*tr(m)
              kmsu = einv(2,1)*ur(m) + einv(2,2)*uz(m)
     2             + einv(2,3)*tz(m) + einv(2,4)*tr(m)
              cintijr = e(i,1)*e(j,1)*kmpu*kmpu/(2.0*ra) 
     1         + (e(i,1)*e(j,2)+e(i,2)*e(j,1))*kmpu*kmsu/(ra+rb)
     2         + e(i,2)*e(j,2)*kmsu*kmsu/(2.0*rb)
            else if(typelyr .eq. 0)then
c-----
c             downward potentials coefficients at top of layer m
c-----
              km1pd = einv(3,1)*ur(m) + einv(3,2)*uz(m)
     1              + einv(3,3)*tz(m) + einv(3,4)*tr(m)
              km1sd = einv(4,1)*ur(m) + einv(4,2)*uz(m)
     1              + einv(4,3)*tz(m) + einv(4,4)*tr(m)
c-----
c             upward potentials coefficients at bottom of layer m
c-----
              kmpu = einv(1,1)*ur(m+1) + einv(1,2)*uz(m+1)
     1             + einv(1,3)*tz(m+1) + einv(1,4)*tr(m+1)
              kmsu = einv(2,1)*ur(m+1) + einv(2,2)*uz(m+1)
     1             + einv(2,3)*tz(m+1) + einv(2,4)*tr(m+1)
              FA=ffunc(ra,zd(m))
              GA=gfunc(ra,zd(m))
              FB=ffunc(rb,zd(m))
              GB=gfunc(rb,zd(m))
              H1=h1func(ra,rb,zd(m))
              H2=h2func(ra,rb,zd(m))
              cintijr = e(i,1)*e(j,1)*kmpu*kmpu*FA 
     1          + e(i,3)*e(j,3)*km1pd*km1pd*FA
     1          + e(i,2)*e(j,2)*kmsu*kmsu*FB 
     1          + e(i,4)*e(j,4)*km1sd*km1sd*FB
     4          + H1*((e(i,1)*e(j,2)+e(i,2)*e(j,1))*kmpu*kmsu +
     5             (e(i,3)*e(j,4)+e(i,4)*e(j,3))*km1pd*km1sd)
     6          + H2*((e(i,1)*e(j,4)+e(i,4)*e(j,1))*kmpu*km1sd +
     7             (e(i,2)*e(j,3)+e(i,3)*e(j,2))*km1pd*kmsu)
     8          + GA*(e(i,1)*e(j,3)+e(i,3)*e(j,1))*kmpu*km1pd
     9          + GB*(e(i,2)*e(j,4)+e(i,4)*e(j,2))*kmsu*km1sd
            else
c-----
c             downward potential coefficients for lower halfspace
c             based on the displacement, stress at bottom
c             Basically m = mmax
c-----
              km1pd = einv(3,1)*ur(m)+einv(3,2)*uz(m)
     2            +einv(3,3)*tz(m)+einv(3,4)*tr(m)
              km1sd = einv(4,1)*ur(m)+einv(4,2)*uz(m)
     2            +einv(4,3)*tz(m)+einv(4,4)*tr(m)
              cintijr = e(i,3)*e(j,3)*km1pd*km1pd/(2.0*ra) 
     1         +(e(i,3)*e(j,4)+e(i,4)*e(j,3))*km1pd*km1sd/(ra+rb)
     2         +e(i,4)*e(j,4)*km1sd*km1sd/(2.0*rb)
            endif
        endif

        intijr = dreal(cintijr)
        return
        end

        subroutine getdcdh(om2,wvno,wvno2,fac)
        implicit none
c-----
c       procedure parameters
c-----
        real*8 om2,wvno,wvno2,fac
c-----
c       common blocks
c-----
        integer NL
        parameter (NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
        integer mmax
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     1                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0
        real*8 dcda, dcdb, dcdr, dcdh
        common/wateri/iwat(NL)
        integer iwat
c-----
c       internal variables
c-----
        real*8 tur, tuz, ttz, ttr
        real*8 dfac, gfac1, gfac2, gfac3, gfac4, gfac5, gfac6
c  DEVELOPMENT
        real*8 drho, dmu, dlm
        real*8 duzdzp, dlur2, xl2mp, xl2mm
        real*8 dl2mu, duzdzm, drur2
        real*8 URB, DURDZM, DURDZP


        integer m
c-----

        do  m=1,mmax
            
            tuz = uz(m)
            ttz = tz(m)
            ttr = tr(m)
            if(iwat(m).eq.1)then
                tur = -wvno*ttz/(zrho(m)*om2)
            else
                tur = ur(m)
            endif
c-----
c       this assumes that the top is a halfspace
c-----
            if(m.eq.1)then
                drho = zrho(1) - 0.0
                dmu  = xmu(1) - 0.0
                dlm = xlam(1) - 0.0
                dl2mu = dlm + dmu + dmu
                xl2mp = xlam(m) + xmu(m) + xmu(m)
                duzdzp = (ttz + wvno*xlam(m)*tur)/xl2mp
                if(iwat(m) .eq.1)then
                    durdzp = wvno*tuz
                else
                    durdzp = (ttr/xmu(m)) - wvno*tuz
                endif
                drur2 = tur*tur*drho
                dlur2 = tur*tur*dl2mu

                gfac1 =  om2*drho*tuz**2
                gfac2 =  om2*drur2
                gfac3 = -wvno2*dmu*tuz**2
                gfac4 = -wvno2*dlur2
                gfac5 =  (xl2mp*duzdzp**2)
                gfac6 =  (xmu(m)*durdzp**2 )
            else
                drho = zrho(m) - zrho(m-1)
                dmu = xmu(m) - xmu(m-1)
                dlm = xlam(m) - xlam(m-1)
                dl2mu = dlm + dmu + dmu
                xl2mp = xlam(m)   + xmu(m)   + xmu(m)
                xl2mm = xlam(m-1) + xmu(m-1) + xmu(m-1)
                duzdzp = (ttz + wvno*xlam(m)*tur)/xl2mp
                if(xmu(m).eq.0.0)then
                    durdzp = wvno*tuz
                else
                    durdzp = (ttr/xmu(m)) - wvno*tuz
                endif
                if(xmu(m-1).eq.0.0 )then
                    durdzm = wvno*tuz
                else
                    durdzm = (ttr/xmu(m-1)) - wvno*tuz
                endif
c-----
c       attempt to fix for water layer, since Ur is not continuous
c       across fluid - solid interface or fluid-fluid interface
c-----
                if(iwat(m-1).eq.1 .and. iwat(m).eq.0 )then
                    URB = -wvno*tz(m)/(zrho(m-1)*om2)
                    drur2 = tur*tur*zrho(m)-URB*URB*zrho(m-1)
                    dlur2 = tur*tur*xl2mp -URB*URB*xl2mm
                    duzdzm = (ttz + wvno*xlam(m-1)*URB)/xlam(m-1)
                else if(iwat(m-1).eq.1 .and. iwat(m).eq.1 )then
                    URB = -wvno*tz(m)/(zrho(m-1)*om2)
                    drur2 = tur*tur*zrho(m)- URB*URB*zrho(m-1)
                    dlur2 = tur*tur*xl2mp  - URB*URB*xl2mm
                    duzdzm = (ttz + wvno*xlam(m-1)*URB)/xl2mm
                else
                    drur2 = tur*tur*drho
                    dlur2 = tur*tur*dl2mu
                    duzdzm = (ttz + wvno*xlam(m-1)*tur)/xl2mm
                endif


                gfac1 =  om2*drho*tuz**2
                gfac2 =  om2*drur2
                gfac3 = -wvno2*dmu*tuz**2
                gfac4 = -wvno2*dlur2
                gfac5 =  (xl2mp*duzdzp**2-xl2mm*duzdzm**2)
                gfac6 =  (xmu(m)*durdzp**2 - xmu(m-1)*durdzm**2)
            endif
            dfac = fac * (
     1          gfac1 + gfac2 + gfac3 + gfac4
     2          + gfac5 + gfac6 )
            if(dabs(dfac).lt.1.0d-38)then
                dfac = 0.0d+00
            endif
            dcdh(m) = dfac
        enddo
        return
        end
