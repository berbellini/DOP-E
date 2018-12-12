        program sregn96
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME IV                                                       c
c                                                                      c
c      PROGRAM: SREGN96                                                c
c                                                                      c
c      COPYRIGHT 1996                                                  c
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
c
c     This program calculates the group velocity and partial
c     derivatives of Rayleigh waves for any plane multi-layered model.
c     The propagator-matrix, instead of numerical-integration
c     method is used, in which the Haskell rather than
c     Harkrider formalisms are concered.
c
c     Developed by C. Y. Wang and R. B. Herrmann, St. Louis University,
c     Oct. 10, 1981.  Modified for use in surface wave inversion, with
c     addition of spherical earth flattening transformation and
c     numerical calculation of group velocity partial derivatives by
c     David R. Russell, St. Louis University, Jan. 1984.
c
c     Layer thickness partial derivatives added following the formalism 
c     of
c     Keilis-Borok et al "Surface waves in laterally heterogeneous media
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c Revision history:
c       04 JUN 2001 - iwat.eq.0.0 -> iwat.eq.0 in subroutine gammap
c       07 AUG 2002 - make string lengths 120 characters from 80 
c       13 OCT 2006 - verbose output of energy integrals
c           introduced proper surface normalization for all fluid case
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      parameter(NL=1000)
      parameter(LER=0,LIN=5,LOT=6)
c-----
c     LIN   - unit for FORTRAN read from terminal
c     LOT   - unit for FORTRAN write to terminal
c     LER   - unit for FORTRAN error output to terminal
c     NL    - number of columns in model 
c-----
        common/isomod/di(NL),ai(NL),bi(NL),rhoi(NL),
     1      qai(NL),qbi(NL),etapi(NL),etasi(NL), 
     2      frefpi(NL), frefsi(NL)
        common/model/  d(NL),a(NL),b(NL),rho(NL),qa(NL), qb(NL)
        real*4 d, a, b, rho, qa, qb
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam


        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0, dcda, dcdb, dcdr, dcdh

        real*4 sdcda(NL), sdcdb(NL), sdcdh(NL), sdcdr(NL) 
        real*4 spur(NL), sptr(NL), spuz(NL), sptz(NL)

        common/wateri/iwat(NL)


        common/sumi/   sumi0,sumi1,sumi2,sumi3,flagr,are,ugr
        real*8 sumi0, sumi1, sumi2, sumi3, flagr, are, ugr
        common/depref/refdep

        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        real*4 vtp,dtp,rtp
        real*4 s1
        real*8 durdz, duzdz, d2urdz, d2uzdz
        real*8 rurdz, ruzdz
        parameter(MAXMOD=2000)
        real*8 cp(MAXMOD), t
        real*8 c, omega, wvno, gamma, csph, usph
        logical nwlyrs, nwlyrr
        character*12 fname(2)
c-----
c       wired in file names     - output of this program is 
c          always in the
c               binary file sregn96.egn or sregn96.der
c                               - input from sdisp96 is always
c               binary file sdisp96.ray
c-----
        logical dolove, dorayl
        character hsfile*120, hrfile*120, title*120
        logical dotmp
        logical ext
        logical dogam
        logical dderiv

        character mname*120
        integer ipar(10)
        real*4 fpar(10)

        integer nipar(10)

        logical verbose
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid

c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line information
c-----  
        call gcmdln(hsfile,hrfile,hs,hr,dotmp,dogam,dderiv,
     1      nipar,verbose)
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
c-----
        do 1234 i=1,mmax
C           qai(i) = qa(i)
C           qbi(i) = qb(i)
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
c       if there is a spherical model, map the source and 
c          receiver depth 
c       in the spherical model to the equivalent depth in the flat model
c-----
        if(nsph.gt.0)then
                dephs = 6370.0 * alog  ( 6370.0/(6370.0 - depths))
                dephr = 6370.0 * alog  ( 6370.0/(6370.0 - depthr))
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
c           If a spherical model is used, reconstruct the thickness
c           of the original model to get the mapping from
c           spherical to flat. We need this for correct
c           spherical partial derivatives. 
c
c           THIS MUST BE DONE BEFORE SOURCE LAYER INSERTION
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
c-----
c       EnD DEBUG
c-----
        twopi=2.*3.141592654
        ifirst = 0
  400   continue
        call gtshed(1,ifunc,mode,t,ierr)
        if(ierr.ne.0)go to 700
        s1=t
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
            call energy(omega,wvno)
c-----
c       the gamma routine will use the spherical model, but the
c       frequency dependence and Q of the original model
c-----
            if(dogam)then
                call gammap(omega,wvno,gamma)
            else
                gamma = 0.0
            endif
c-----
c-----
c       output the derivatives.
c-----
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
            if(dabs(gamma).lt.1.0d-36)gamma=0.0d+00
            do 500 i=1,mmax
c-----
c       correct for differences between Haskell and Levshin and Yanson
c-----
                ur(i) = wvno * ur(i)
C               if(iwat(i).eq.1)then
C                   ur(i) = - zrho(i)*tz(i)/wvno
C               endif
                tz(i) = omega*omega * tz(i)
                tr(i) = wvno*omega*omega * tr(i)
  500       continue
c-----
c       everything is in Levshin and Yanson format
c-----
            if(nsph.gt.0)then
c-----
c           sphericity correction for partial derivatives 
c           of the original model
c-----
                call sprayl(omega,dble(c),mmaxot,csph,usph,ugr)
            endif

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
                 WRITE(LOT,2)t,c,ugr,gamma,
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
                wvno = omega/csph
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
            mmaxot = mmax
C       write(6,*)'0:',6.2831853/omega,c,( dcdh(j),j=1,mmaxot)
            if(nwlyrr)then
                call collap(lrro+1,mmaxot)
            endif
C       write(6,*)'1:',6.2831853/omega,c,( dcdh(j),j=1,mmaxot)
            if(nwlyrs)then
                call collap(lsso+1,mmaxot)
            endif
C       write(6,*)'2:',6.2831853/omega,c,( dcdh(j),j=1,mmaxot)
            call chksiz(ur,spur,mmaxot)
            call chksiz(tr,sptr,mmaxot)
            call chksiz(uz,spuz,mmaxot)
            call chksiz(tz,sptz,mmaxot)
            call chksiz(dcdh,sdcdh,mmaxot)
            call chksiz(dcda,sdcda,mmaxot)
            call chksiz(dcdb,sdcdb,mmaxot)
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
C       write(6,*)'3:',6.2831853/omega,c,(sdcdh(j),j=1,mmaxot)

            call putder(2,6,sngl(wvno),sngl(ugr), 
     1          sngl(gamma), 
     1          sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0,
     2          rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3          sumkr,sumgr,sumgv,mmaxot,
     4          sdcdh,sdcda,sdcdb,sdcdr,
     5          spur,sptr,spuz,sptz,ipar)
        else
            call putegn(2,2,1,sngl(wvno),sngl(ugr),
     1          sngl(gamma),
     1          sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2          rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3          sumkr,sumgr,sumgv)
        endif
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
c       stop
        end

      subroutine svfunc(omega,wvno)
c
c     This combines the Haskell vector from sub down and
c     Dunkin vector from sub up to form the eigenfunctions.
c
      implicit double precision (a-h,o-z)
      integer*4 NL
      parameter (NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
      common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
      common/dunk/   uu(NL,5),exe(NL),exa(NL),exh(NL)
      common/hask/   vv(NL,4)
        common/hwat/wh(NL,2),hex(NL)
        common/wateri/iwat(NL)
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid
c
        call up(omega,wvno,fr)
        call down(omega,wvno)
        wvno2=wvno*wvno
        omega2=omega*omega
c------
c     f5 is the fifth component is compound R| 12 ij.
c     uu0(1) is the ellipticity.
c     uu0(2) is Uz at the surface. It is used for normalization.
c     uu0(3)=tz is actually the period equation.
c     uu0(4)=tr should be zero.
c------
        f5=uu(1,4)
        uu0(1)=wvno*uu(1,3)/f5
        uu0(2)=1.0d+00
        uu0(3)=fr
        uu0(4)=fr
        ur(1)=uu(1,3)/f5
        uz(1) = 1.0d+00
        tz(1) = 0.0d+00
        tr(1) = 0.0d+00
c------
c     uu is Z and vv is X in Wang s dissertation eq. III-1-6 (p.76).
c     i.e., Z a Haskell vector from top layer down and X a Dunkin
c     compound matrix from bottom layer up.
c------
        do 200 i=2,mmax
            i1=i-1
            uu1 =
     *        vv(i1,2)*uu(i,1)+vv(i1,3)*uu(i,2)      +vv(i1,4)*uu(i,3)
            uu2 =
     *       -vv(i1,1)*uu(i,1)-vv(i1,3)*uu(i,3)*wvno2+vv(i1,4)*uu(i,4)
            uu3 =
     *       -vv(i1,1)*uu(i,2)+vv(i1,2)*uu(i,3)*wvno2+vv(i1,4)*uu(i,5)
            uu4 =
     *       -vv(i1,1)*uu(i,3)-vv(i1,2)*uu(i,4)      -vv(i1,3)*uu(i,5)
            ext=0.0d+00
            do 100 k=1,i1
                ext=ext+exa(k)-exe(k)
  100       continue
            if(ext.gt.-80.0 .and. ext.lt.80.0 ) then
                fact=dexp(ext)

                ur(i)=uu1*fact/f5
                uz(i)=uu2*fact/f5
                tz(i)=uu3*fact/f5
                tr(i)=uu4*fact/f5
            else
                ur(i) = 0.0
                uz(i) = 0.0
                tz(i) = 0.0
                tr(i) = 0.0
            endif
  200   continue
c-----
c       Continue KLUDGE if top layers are water
c       The problem here was that even though the 
c          phase velocity is determined,
c       the computed eigenfunctions do not match perfectly,
c       and we do not match the free surface requirement that Tz=0 at
c       top of fluid. So RBH used the fluid propagator from top
c       down to the solid-fluid interface, the Haskell from
c       base up to the same interface, and then adjust the amplitudes
c       to match. The purpose is to get a better estimate of the
c       eigenfunction in the fluid.
c
c       Ultimately use Kennett
c-----
c       first get the largest exponent
c-----
        if(iwat(1).eq.1)then
c-----
c           OK first layer is fluid
c-----
            hmax = 0.0
            jwat=0
            do 300 i=1,mmax-1
                if(iwat(i).eq.1)then
                    jwat = i
                    hmax = hex(i+1)
                endif
  300       continue
c-----
c           now get correction factor
c-----
            corrd = sqrt(uz(jwat+1)**2 + tz(jwat+1)**2 )
            corrw = sqrt(wh(jwat+1,1)**2 + wh(jwat+1,2)**2 )
            if(corrd.eq.0.0)then
c-----
c       give precedence to surface values
c-----
                fac = 0.0
            else
                fac = corrw/corrd
            endif
            fac = fac * exp(hmax - hmax)
            fac = fac * 2.0
c-----
c           eventually normalize
c-----
        IF(.not. ALLFLUID)then
            do 400 i=1,mmax
                if(i.le.jwat)then
                    efac = exp(hex(i) - hmax)
                    efac = efac * 2.0
                    ur(i) = 0.0
                    uz(i) = wh(i,1)*efac
                    tz(i) = wh(i,2)*efac
                    tr(i) = 0.0
                else
                    ur(i) = ur(i) * fac
                    uz(i) = uz(i) * fac
                    tr(i) = tr(i) * fac
                    tz(i) = tz(i) * fac
                endif
  400       continue
        ENDIF
        endif
        
        return
        end

        subroutine up(omega,wvno,fr)
c-----
c       This finds the values of the Dunkin vectors at
c       each layer boundaries from bottom layer upward.
c-----
        implicit double precision (a-h,o-z)
        integer*4 NL
        parameter (NL=1000)
        real*8 ee0(5)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        common/dunk/   uu(NL,5),exe(NL),exa(NL),exh(NL)
c-----
c       The save labeled common is used to pass the
c       Dunkin 5x5 and Haskell 4x4 matrices to the main program
c       e(ex1) is the scaling for the Dunkin
c       e(ex2) is the scaling for the Haskell matrices
c-----
        common/save/   dd(5,5),aa(4,4),ex1,ex2
        common/aamatx/ ww(NL),xx(NL),yy(NL),zz(NL),
     1      cospp(NL),cosqq(NL)
        common/engerw/ wra,wd,wba
        common/wateri/iwat(NL)
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid
c-----
c       initialize base vector
c-----
        wvno2=wvno*wvno
        xka=omega/za(mmax)
        wvnop=wvno+xka
        wvnom=dabs(wvno-xka)
        ra=dsqrt(wvnop*wvnom)
        rrho=zrho(mmax)
c-----
c       reversed from Herrmann and sdisp96 because of 
c          way matrices are ordered
c-----
        if(allfluid)then
            uu(mmax,1)=0.0
            uu(mmax,2)=0.0
            uu(mmax,3)=0.0
            uu(mmax,4)= -ra/0.5d+00
            uu(mmax,5)= -rrho/0.5d+00
        else
            xkb=omega/zb(mmax)
            wvnop=wvno+xkb
            wvnom=dabs(wvno-xkb)
            rb=dsqrt(wvnop*wvnom)
            t = zb(mmax)/omega
            gammk = 2.d+00*t*t
            gam = gammk*wvno2
            gamm1 = gam - 1.d+00
c------
c       since symmetry of compund matrix (component 3 and 4)
c       only 5 terms needed.
c------
            uu(mmax,1)=wvno2-ra*rb
            uu(mmax,2)=-rrho*rb
            uu(mmax,3)=rrho*(gamm1-gammk*ra*rb)
            uu(mmax,4)=rrho*ra
            uu(mmax,5)=rrho*rrho*(gamm1*gamm1-gam*gammk*ra*rb)
        endif
c------
c       matrix multiplication from bottom layer upward
c------
        mmx1=mmax-1
        do 400 k=mmx1,1,-1
            k1=k+1
            dpth=zd(k)
            rrho=zrho(k)
            xka = omega/za(k)
            t = zb(k)/omega
            gammk = 2.d+00*t*t
            gam = gammk*wvno2
            wvnop=wvno+xka
            wvnom=dabs(wvno-xka)
            ra=dsqrt(wvnop*wvnom)
            p=ra*dpth
            if(iwat(k).eq.1)then
                nmat = 5
                call varw(k,p,ra,wvno,xka,dpth,w,x,cosp,ex)
                ex2 = ex 
                ex1 = ex
                call dnkaw(rrho,w,x,cosp)
                exh(k) = ex
            else
                nmat = 5
                xkb = omega/zb(k)
                wvnop=wvno+xkb
                wvnom=dabs(wvno-xkb)
                rb=dsqrt(wvnop*wvnom)
                q=rb*dpth
                call var(k,p,q,ra,rb,wvno,xka,xkb,dpth)
                call dnka(wvno2,gam,gammk,rrho)
            endif
            do 200 i=1,nmat
                cc=0.0d+00
                do 100 j=1,nmat
                    cc=cc+dd(i,j)*uu(k1,j)
  100           continue
                ee0(i)=cc
  200       continue
c------
c       normalization of ee is important to prevent over- or
c            under-flow.
c------
            rab = 0.0
            call normc(ee0,rab,nmat)
c------
c       exe stores the exponential terms and carefully controls
c       the precision problem for P-SV Haskell matrix.
c------
            exe(k)=ex1+rab
            exa(k)=ex2
            do 300 i=1,nmat
                uu(k,i)=ee0(i)
  300       continue
  400   continue
c-----
c       define period equation
c-----
        fr=uu(1,5)
        return
        end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        subroutine down(omega,wvno)
c-----
c       This finds the values of the Haskell vectors at
c       each layer boundaries from top layer downward.
c-----
        implicit double precision (a-h,o-z)
        integer*4 NL
        parameter (NL=1000)
        dimension aa0(5)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
        common/dunk/   uu(NL,5),exe(NL),exa(NL),exh(NL)
        common/hask/   vv(NL,4)
        common/save/   dd(5,5),aa(4,4),ex1,ex2
        common/aamatx/ ww(NL),xx(NL),yy(NL),zz(NL),
     *                 cospp(NL),cosqq(NL)
        common/wateri/iwat(NL)
        common/hwat/wh(NL,2),hex(NL)
        wvno2=wvno*wvno
        do 100 j=1,4
            vv(1,j)=0.0d+00
  100   continue
        vv(1,4)=1.0d+00
        aa0(5)=0.0d+00
        mmx1=mmax-1
c-----
c       Kludge around fluids
c       wh(1,1) = free surface displacement
c       wh(1,2) = free surface vertical stress
c-----
        wh(1,1) = 1.0
        wh(1,2) = 0.0
        hex(1) = 0.0
        do 500 k=1,mmx1
            k1=k-1
            if(k.eq.1) k1=1
            t=zb(k)/omega
            gammk=2.d+00*t*t
            gam=gammk*wvno2
            w=ww(k)
            x=xx(k)
            y=yy(k)
            z=zz(k)
            cosp=cospp(k)
            cosq=cosqq(k)
            rrho=zrho(k)
            ex = exh(k)
            if(iwat(k).eq.1)then
                call hskaw(w,x,y,z,cosp,cosq,
     1              wvno2,gam,gammk,rrho,ex)
                wh(k+1,1) = cosp * wh(k,1) - x * wh(k,2)/rrho
                wh(k+1,2) = - rrho*w*wh(k,1) + cosp*wh(k,2)
                hex(k+1) = hex(k) + ex
            else
                call hska (w,x,y,z,cosp,cosq,wvno2,gam,gammk,rrho)
            endif
            do 300 j=1,4
                cc=0.0d+00
                do 200 i=1,4
                    cc=cc+vv(k1,i)*aa(i,j)
  200           continue
            aa0(j)=cc
  300       continue
            call normc(aa0,ex2,4)
            exa(k)=exa(k)+ex2
            
            do 400 i=1,4
                vv(k,i)=aa0(i)
  400       continue
  500   continue
        return
        end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        subroutine hska (w,x,y,z,cosp,cosq,wvno2,gam,gammk,rho)
c-----
c       Haskell matrix
c-----
        implicit double precision (a-h,o-z)
        common/save/ dd(5,5),aa(4,4),ex1,ex2
        gamm1 = gam-1.d+00
        temp = x-wvno2*y
        aa(2,3) = temp/rho
        temp = temp*gammk
        aa(4,3) = temp+y
        aa(4,1) = -(gamm1*y+gam*aa(4,3))*rho
        aa(2,1) = -wvno2*aa(4,3)
        temp = cosq-cosp
        aa(1,3) = temp/rho
        aa(2,4) = -wvno2*aa(1,3)
        aa(4,2) = rho*gammk*gamm1*temp
        aa(3,1) = -wvno2*aa(4,2)
        temp = temp*gam
        aa(1,1) = cosq-temp
        aa(2,2) = cosp+temp
        aa(3,3) = aa(2,2)
        aa(4,4) = aa(1,1)
        temp = z-wvno2*w
        aa(1,4) = temp/rho
        aa(1,2) = -w-gammk*temp
        aa(3,4) = -wvno2*aa(1,2)
        aa(3,2) = rho*(gam*aa(1,2)-gamm1*w)
        return
        end

        subroutine hskaw(w,x,y,z,cosp,cosq,wvno2,gam,gammk,rho,ex)
c-----
c       Haskell matrix for water layer
c-----
        implicit double precision (a-h,o-z)
        common/save/ dd(5,5),aa(4,4),ex1,ex2
        if(ex.lt.80.0)then
            fac = dexp(-ex)
        else
            fac = 0.0
        endif
        aa(1,1) = fac
        aa(1,2) = 0.0d+00
        aa(1,3) = 0.0d+00
        aa(1,4) = 0.0d+00
        aa(2,1) = 0.0d+00
        aa(2,2) = cosp
        aa(2,3) =   x / rho
        aa(2,4) = 0.0d+00
        aa(3,1) = 0.0d+00
        aa(3,2) =   rho * w
        aa(3,3) = cosp
        aa(3,4) = 0.0d+00
        aa(4,1) = 0.0d+00
        aa(4,2) = 0.0d+00
        aa(4,3) = 0.0d+00
        aa(4,4) = fac
        return
        end

        subroutine dnka(wvno2,gam,gammk,rho)
c-----
c       Dunkin matrix.
c-----
        implicit double precision (a-h,o-z)
        common/save/ dd(5,5),aa(4,4),ex1,ex2
        common/ovrflw/ a0,cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        gamm1 = gam-1.d+00
        twgm1=gam+gamm1
        gmgmk=gam*gammk
        gmgm1=gam*gamm1
        gm1sq=gamm1*gamm1
        rho2=rho*rho
c------
c       a0 takes care of the factor being exponentially controlled.
c------
        a0pq=a0-cpcq
        dd(1,1)=cpcq-2.d+00*gmgm1*a0pq-gmgmk*xz-wvno2*gm1sq*wy
        dd(2,1)=(gm1sq*cqw-gmgmk*cpz)*rho
        dd(3,1)=-(gammk*gamm1*twgm1*a0pq+gam*gammk*gammk*xz+
     *          gamm1*gm1sq*wy)*rho
        dd(4,1)=(gmgmk*cqx-gm1sq*cpy)*rho
        dd(5,1)=-(2.d+00*gmgmk*gm1sq*a0pq+gmgmk*gmgmk*xz+
     *          gm1sq*gm1sq*wy)*rho2
        dd(1,2)=(cqx-wvno2*cpy)/rho
        dd(2,2)=cpcq
        dd(3,2)=gammk*cqx-gamm1*cpy
        dd(4,2)=-xy
        dd(5,2)=dd(4,1)
        dd(1,4)=(wvno2*cqw-cpz)/rho
        dd(2,4)=-wz
        dd(3,4)=gamm1*cqw-gammk*cpz
        dd(4,4)=dd(2,2)
        dd(5,4)=dd(2,1)
        dd(1,5)=-(2.d+00*wvno2*a0pq+xz+wvno2*wvno2*wy)/rho2
        dd(2,5)=dd(1,4)
        dd(3,5)=-(twgm1*a0pq+gammk*xz+wvno2*gamm1*wy)/rho
        dd(4,5)=dd(1,2)
        dd(5,5)=dd(1,1)
        t=-2.d+00*wvno2
        dd(1,3)=t*dd(3,5)
        dd(2,3)=t*dd(3,4)
        dd(3,3)=a0+2.d+00*(cpcq-dd(1,1))
        dd(4,3)=t*dd(3,2)
        dd(5,3)=t*dd(3,1)
        return
        end

        subroutine dnkaw(rho1,w,x,cosp)
        implicit double precision (a-h,o-z)
        common/save/ dd(5,5),aa(4,4),ex1,ex2
            do 100 j=1,5
                do 101 i=1,5
                    dd(i,j) = 0.0d+00
  101           continue
                dd(j,j) = 1.0d+00
  100       continue
        if(ex1.lt.90.0d+00)then
            fac = dexp(-ex1)
        else if(ex1.eq.0.0d+00)then
            fac = 1.0d+00
            fac = 0.0
        endif
            dd(1,1) = cosp
            dd(1,2) =   x/ rho1
            dd(2,1) =   rho1*w
            dd(2,2) = cosp
            dd(3,3) = fac
            dd(4,4) = cosp
            dd(4,5) = dd(1,2)
            dd(5,4) = dd(2,1)
            dd(5,5) = cosp
        return
        end

      subroutine var(m,p,q,ra,rb,wvno,xka,xkb,dpth)
c-----
c     find variables cosP, cosQ, sinP, sinQ, etc.
c     as well as cross products required for compound matrix
c-----
c     To handle the hyperbolic functions correctly for large
c     arguments, we use an extended precision procedure,
c     keeping in mind that the maximum precision in double
c     precision is on the order of 16 decimal places.
c
c     So  cosp = 0.5 ( exp(+p) + exp(-p))
c              = exp(p) * 0.5 * ( 1.0 + exp(-2p) )
c     becomes
c         cosp = 0.5 * (1.0 + exp(-2p) ) with an exponent p
c     In performing matrix multiplication, we multiply the modified
c     cosp terms and add the exponents. At the last step
c     when it is necessary to obtain a true amplitude,
c     we then form exp(p). For normalized amplitudes at any depth,
c     we carry an exponent for the numerator and the denominator, and
c     scale the resulting ratio by exp(NUMexp - DENexp)
c
c     The propagator matrices have three basic terms
c
c     HSKA        cosp  cosq
c     DUNKIN      cosp*cosq     1.0
c
c     When the extended floating point is used, we use the
c     largest exponent for each, which is  the following:
c
c     Let pex = p exponent > 0 for evanescent waves = 0 otherwise
c     Let sex = s exponent > 0 for evanescent waves = 0 otherwise
c     Let exa = pex + sex
c
c     Then the modified matrix elements are as follow:
c
c     Haskell:  cosp -> 0.5 ( 1 + exp(-2p) ) exponent = pex
c               cosq -> 0.5 ( 1 + exp(-2q) ) * exp(q-p)
c                                            exponent = pex
c            (this is because we are normalizing all elements in the
c             Haskell matrix )
c    Compound:
c              cosp * cosq -> normalized cosp * cosq exponent 
c                     = pex + qex
c               1.0  ->    exp(-exa)
c-----
      implicit double precision (a-h,o-z)
      integer*4 NL
      parameter (NL=1000)
      common/save/ dd(5,5),aa(4,4),exa,ex
      common/aamatx/w0(NL),x0(NL),y0(NL),z0(NL),
     1              cosp0(NL),cosq0(NL)
      common/ovrflw/   a0,cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
      exa=0.0d+00
      ex =0.0d+00
      a0=1.0d+00
c-----
c     examine P-wave eigenfunctions
c        checking whether c> vp c=vp or c < vp
c-----
      pex = 0.0d+00
      sex = 0.0d+00
      if(wvno.lt.xka)then
             sinp = dsin(p)
             w=sinp/ra
             x=-ra*sinp
             cosp=dcos(p)
      elseif(wvno.eq.xka)then
             cosp = 1.0d+00
             w = dpth
             x = 0.0d+00
      elseif(wvno.gt.xka)then
             pex = p
             fac = 0.0d+00
             if(p.lt.18.0d+00)fac = dexp(-2.0d+00*p)
             cosp = ( 1.0d+00 + fac) * 0.5d+00
             sinp = ( 1.0d+00 - fac) * 0.5d+00
             w=sinp/ra
             x=ra*sinp
      endif
c-----
c     examine S-wave eigenfunctions
c        checking whether c > vs, c = vs, c < vs
c-----
      if(wvno.lt.xkb)then
             sinq=dsin(q)
             y=sinq/rb
             z=-rb*sinq
             cosq=dcos(q)
      elseif(wvno.eq.xkb)then
             cosq=1.0d+00
             y=dpth
             z=0.0d+00
      elseif(wvno.gt.xkb)then
             sex = q
             fac = 0.0d+00
             if(q.lt.18.0d+00)fac = dexp(-2.0d+0*q)
             cosq = ( 1.0d+00 + fac ) * 0.5d+00
             sinq = ( 1.0d+00 - fac ) * 0.5d+00
             y = sinq/rb
             z = rb*sinq
      endif
c-----
c     form eigenfunction products for use with compound matrices
c-----
      exa = pex + sex
      ex = pex
      a0=0.0d+00
      if(exa.lt.80.0d+00) a0=dexp(-exa)
      cpcq=cosp*cosq
      cpy=cosp*y
      cpz=cosp*z
      cqw=cosq*w
      cqx=cosq*x
      xy=x*y
      xz=x*z
      wy=w*y
      wz=w*z
      qmp = sex - pex
      fac = 0.0d+00
      if(qmp.gt.-40.0d+00)fac = dexp(qmp)
      cosq = cosq*fac
      y=fac*y
      z=fac*z
      w0(m)=w
      x0(m)=x
      y0(m)=y
      z0(m)=z
      cosp0(m)=cosp
      cosq0(m)=cosq
      return
      end

        subroutine varw(m,p,ra,wvno,xka,dpth,w,x,cosp,exa)
c-----
c       find variables cosP, sinP, etc.
c       as well as cross products required for compound matrix
c-----
c       To handle the hyperbolic functions correctly for large
c       arguments, we use an extended precision procedure,
c       keeping in mind that the maximum precision in double
c       precision is on the order of 16 decimal places.
c
c       So  cosp = 0.5 ( exp(+p) + exp(-p))
c                = exp(p) * 0.5 * ( 1.0 + exp(-2p) )
c       becomes
c           cosp = 0.5 * (1.0 + exp(-2p) ) with an exponent p
c       In performing matrix multiplication, we multiply the modified
c       cosp terms and add the exponents. At the last step
c       when it is necessary to obtain a true amplitude,
c       we then form exp(p). For normalized amplitudes at any depth,
c       we carry an exponent for the numerator and the denominator, and
c       scale the resulting ratio by exp(NUMexp - DENexp)
c
c       The propagator matrices have three basic terms
c
c       HSKA        cosp  cosq
c       DUNKIN      cosp*cosq     1.0
c
c       When the extended floating point is used, we use the
c       largest exponent for each, which is  the following:
c
c       Let pex = p exponent > 0 for evanescent waves = 0 otherwise
c       Let exa = pex + sex
c
c       Then the modified matrix elements are as follow:
c
c       Haskell:  cosp -> 0.5 ( 1 + exp(-2p) ) exponent = pex
c                                              exponent = pex
c              (this is because we are normalizing all elements in the
c               Haskell matrix )
c    Compound:
c                cosp * cosq -> normalized cosp * cosq exponent 
c                     = pex + qex
c                 1.0  ->    exp(-exa)
c-----
        implicit double precision (a-h,o-z)
        parameter(NL=1000)
        common/aamatx/w0(NL),x0(NL),y0(NL),z0(NL),
     1      cosp0(NL),cosq0(NL)
        exa=0.0d+00
        a0=1.0d+00
c-----
c       examine P-wave eigenfunctions
c          checking whether c> vp c=vp or c < vp
c-----
        pex = 0.0d+00
        if(wvno.lt.xka)then
               sinp = dsin(p)
               w=sinp/ra
               x=-ra*sinp
               cosp=dcos(p)
        else if(wvno.eq.xka)then
               cosp = 1.0d+00
               w = dpth
               x = 0.0d+00
        else if(wvno.gt.xka)then
               pex = p
               fac = 0.0d+00
               if(p.lt.18.00d+00)fac = dexp(-2.0d+00*p)
               cosp = ( 1.0d+00 + fac) * 0.5d+00
               sinp = ( 1.0d+00 - fac) * 0.5d+00
               w=sinp/ra
               x=ra*sinp
        endif
c-----
c       form eigenfunction products for use with compound matrices
c-----
        exa = pex 
        ex = pex
        w0(m)=w
        x0(m)=x
        cosp0(m)=cosp
        return
        end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        subroutine normc(ee,ex,nmat)
c
c     This is an important step to control over- or
c     underflow.
c     The Haskell or Dunkin vectors are normalized before
c     the layer matrix stacking.
c     Note that some precision might be lost during normalization.
c
        implicit double precision (a-h,o-z)
        dimension ee(5)
        t1 = 0.0d+00
        do 10 i = 1,nmat
            if(dabs(ee(i)).gt.t1) t1 = dabs(ee(i))
   10   continue
        if(t1.lt.1.d-30) t1=1.d+00
        do 20 i =1,nmat
            ee(i)=ee(i)/t1
   20   continue
        ex=dlog(t1)
        return
        end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        subroutine energy(omega,wvno)
c------
c     This takes the energy integral by an analytic
c     way using the eigenfuntions found above.
c     The amplitude factor (sum of energy integrals) thus can be
c     accurately determined without taking any sort of derivatives.
c     The first version of this subroutine is provided by Dr. Harkrider.
c------
        implicit double precision (a-h,o-z)
        integer*4 NL
        parameter (NL=1000)
        complex*16 t(4,4),tt(4,4),ff(6),pp(4)
        complex*16 nua,nub,c1,c2
        common/coef/   t,tt,ff,pp
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
        common/eigfun/ uu(NL,4),uu0(4),dcda(NL),dcdb(NL),
     *                 dcdr(NL),dcdh(NL)
        common/sumi/   sumi0,sumi1,sumi2,sumi3,flagr,are,ugr
        common/wateri/iwat(NL)
c-----
c       may have added a layer to put the source into the model
c       we do not want this to affect the concept of partial
c       derivatives
c-----
        omega2=omega*omega
        wvomg2=wvno*omega2
        wvno2=wvno*wvno
        sumi0=0.0d+00
        sumi1=0.0d+00
        sumi2=0.0d+00
        sumi3=0.0d+00
        do 300 k=1,mmax
            k1=k+1
            kk=k
            if(k.eq.mmax) kk=NL+1
            dpth=zd(k)
            daa=za(k)
            drho=zrho(k)
            dlam=xlam(k)
            dbb=zb(k)
            dmu=xmu(k)
            dlamu=dlam+2.d+00*dmu
            xka=omega/daa
            wvnop=wvno+xka
            wvnom=dabs(wvno-xka)
            ra=dsqrt(wvnop*wvnom)
            if(ra.lt.1.d-15) ra=1.d-15
            if(wvno.lt.xka) then
                nua=dcmplx(0.0d+00,ra)
            else
                nua=dcmplx(ra,0.0d+00)
            endif
            if(iwat(k).eq.1)then
                uzup = uu(k,2)
                tzup = uu(k,3)
                urup = -uu(k,3)/zrho(k)
                uzdn = uu(k1,2)
                tzdn = uu(k1,3)
                uu(k,1) = - tzup/drho
                urdn = -uu(k1,3)/zrho(k)
      if(k.eq.mmax)WRITE(6,*)'uzup,tzup,urup,uzdn,tzdn,urdn:',
     1      uzup,tzup,urup,uzdn,tzdn,urdn
                call winteg(urur,uzuz,urduz,duzduz,
     1              dpth,nua,drho,
     1              omega2,wvomg2,wvno2,
     2              uzup,tzup,urup,uzdn,tzdn,urdn,kk)
                urur = urur * wvno2 
                duzduz = duzduz 
                urduz = urduz * wvno 
C       WRITE(6,*)'urur,uzuz,urduz,duzduz:',urur,uzuz,urduz,duzduz
                sumi0=sumi0+drho*(urur+uzuz)
                sumi1=sumi1+dlam*urur
                sumi2=sumi2-dlam*urduz
                sumi3=sumi3+dlam*duzduz
                dldl=-wvno2*urur+2.d+00*wvno*urduz-duzduz
                dldr=omega2*(urur+uzuz)
                dcda(k)=2.d+00*drho*daa*omega*dldl/wvno2
                dcdb(k)=0.0d+00
                dcdr(k)=dldr+dlam*dldl/drho
                dcdr(k)=dcdr(k)*omega/wvno2
            else
                xkb=omega/dbb
                gammk=dbb/omega
                gammk=2.d+00*gammk*gammk
                gam=gammk*wvno2
                wvnop=wvno+xkb
                wvnom=dabs(wvno-xkb)
                rb=dsqrt(wvnop*wvnom)
c----
c subterfuge to avoid division by zero
c the matrices should be changed
c-----
                if(rb.lt.1.d-15) rb=1.d-15
                if(wvno.lt.xkb) then
                    nub=dcmplx(0.0d+00,rb)
                else
                    nub=dcmplx(rb,0.0d+00)
                endif
                call tminus(nua,nub,wvno2,gam,gammk,drho)
                call tplus(nua,nub,wvno2,gam,gammk,drho)
                call intval(kk,nua,nub,dpth)
                do 200 i=1,2
                    i2=i+2
                    c1=0.0d+00
                    c2=0.0d+00
                    do 100 j=1,4
                        if(kk.ne.(NL+1)) then
                        c1=c1+tt(i,j)*uu(k1,j)
                        endif
                        c2=c2+tt(i2,j)*uu(k,j)
  100               continue
                    pp(i)=c1
                    pp(i2)=c2
  200           continue
                urur=ssum(kk,1,1)*wvno2
                urtz=ssum(kk,1,3)*wvomg2
                uzuz=ssum(kk,2,2)
                uztr=ssum(kk,2,4)*wvomg2
                tztz=ssum(kk,3,3)*omega2*omega2
                trtr=ssum(kk,4,4)*wvomg2*wvomg2
                urduz=(wvno*dlam*urur+urtz)/dlamu
                uzdur=-wvno*uzuz+uztr/dmu
                durdur=wvno2*uzuz -
     1              2.d+00*wvno*uztr/dmu+trtr/(dmu*dmu)
                duzduz=(wvno2*dlam*dlam*urur
     1              +2.d+00*wvno*dlam*urtz+tztz)/(dlamu*dlamu)
                sumi0=sumi0+drho*(uzuz+urur)
                sumi1=sumi1+dlamu*urur+dmu*uzuz
                sumi2=sumi2+dmu*uzdur-dlam*urduz
                sumi3=sumi3+dlamu*duzduz+dmu*durdur
                dldl=-wvno2*urur+2.d+00*wvno*urduz-duzduz
                dldm=2.d+00*(wvno2*urur+wvno*uzdur+duzduz)+
     1              wvno2*uzuz+durdur
                dldm=-dldm
                dldr=omega2*(urur+uzuz)
                dcda(k)=2.d+00*drho*daa*omega*dldl/wvno2
                dcdb(k)=2.d+00*drho*dbb*omega*(dldm-2.*dldl)/wvno2
                dcdr(k)=dldr+dlam*dldl/drho+dmu*dldm/drho
                dcdr(k)=dcdr(k)*omega/wvno2
            endif
  300   continue
        dldk=-2.d+00*(wvno*sumi1+sumi2)
        do 400 k=1,mmax
            dcda(k)=dcda(k)/dldk
            dcdb(k)=dcdb(k)/dldk
            dcdr(k)=dcdr(k)/dldk
  400   continue
c------
c     flagr is the Lagrangian, which should be zero in perfect case.
c     ugr is the group velocity.
c     are is the amplitude factor.
c
c       Note all partials can be derived from the Lagrangian
c       
c------
        flagr=omega2*sumi0-wvno2*sumi1-2.d+00*wvno*sumi2-sumi3
        ugr=(wvno*sumi1+sumi2)/(omega*sumi0)
        are=wvno/(2.d+00*omega*ugr*sumi0)
        c = omega/wvno
        fac = are*c/wvno2
C       WRITE(6,*)'ARE,C,UGR,SUMI0,I1,I2,I3,LR:',
c          are,c,ugr,sumi0,sumi1,sumi2,sumi3,flagr
c-----
c       to convert to Levshin and Yanson form we must
c
c       Ur (LJ) = k Ur(Haskell)
c       Uz (LJ) = Uz (Haskell)
c       Tz (LJ) = omega**2 * Tz (Haskell)
c       Tr (LJ) = k omega**2 * Tr (Haskell)
c
c-----
c-----
c       now a lot of effort to get partial with respect to 
c          layer thickness
c-----
        do 500 k=1,mmax
            uz = uu(k,2)
            tz = omega2*uu(k,3)
            tr = wvno*omega2*uu(k,4)
            if(iwat(k).eq.1)then
                ur = -wvno*uu(k,3)/zrho(k)
            else
                ur = wvno * uu(k,1)
            endif
C       WRITE(6,*)'k,wvno,omega,uz,tz,tr,ur:',k,wvno,omega,uz,tz,tr,ur
c-----
c       this assumes that the top is a halfspace
c-----
            if(k.eq.1)then
                drho = zrho(1) - 0.0
                dmu  = xmu(1) - 0.0
                dlm = xlam(1) - 0.0
                dl2mu = dlm + dmu + dmu
                xl2mp = xlam(k) + xmu(k) + xmu(k)
                duzdzp = (tz + wvno*xlam(k)*ur)/xl2mp
                if(xmu(k) .eq.0.0)then
                    durdzp = wvno*uz
                else
                    durdzp = (tr/xmu(k)) - wvno*uz
                endif
                drur2 = ur*ur*drho
                dlur2 = ur*ur*dl2mu

                gfac1 =  omega2*drho*uz**2
                gfac2 =  omega2*drur2
                gfac3 = -wvno2*dmu*uz**2
                gfac4 = -wvno2*dlur2
                gfac5 =  (xl2mp*duzdzp**2)
                gfac6 =  (xmu(k)*durdzp**2 )
C       WRITE(6,*)'drho,dmu,dl2mu,drur2,dlur2,duzdzp,durdzp:',
c          drho,dmu,dl2mu,drur2,dlur2,duzdzp,durdzp
            else
                drho = zrho(k) - zrho(k-1)
                dmu = xmu(k) - xmu(k-1)
                dlm = xlam(k) - xlam(k-1)
                dl2mu = dlm + dmu + dmu
                xl2mp = xlam(k)   + xmu(k)   + xmu(k)
                xl2mm = xlam(k-1) + xmu(k-1) + xmu(k-1)
                duzdzp = (tz + wvno*xlam(k)*ur)/xl2mp
                if(xmu(k).eq.0.0)then
                    durdzp = wvno*uz
                else
                    durdzp = (tr/xmu(k)) - wvno*uz
                endif
                if(xmu(k-1).eq.0.0 )then
                    durdzm = wvno*uz
                else
                    durdzm = (tr/xmu(k-1)) - wvno*uz
                endif
c-----
c       attempt to fix for water layer, since Ur is not continuous
c       across fluid - solid interface or fluid-fluid interface
c-----
                if(iwat(k-1).eq.1 .and. iwat(k).eq.0 )then
                    URB = -wvno*uu(k,3)/zrho(k-1)
                    drur2 = ur*ur*zrho(k)-URB*URB*zrho(k-1)
                    dlur2 = ur*ur*xl2mp -URB*URB*xl2mm
                    duzdzm = (tz + wvno*xlam(k-1)*URB)/xlam(k-1)
                else if(iwat(k-1).eq.1 .and. iwat(k).eq.1 )then
                    URB = -wvno*uu(k,3)/zrho(k-1)
                    drur2 = ur*ur*zrho(k)- URB*URB*zrho(k-1)
                    dlur2 = ur*ur*xl2mp  - URB*URB*xl2mm
                    duzdzm = (tz + wvno*xlam(k-1)*URB)/xl2mm
                else
                    drur2 = ur*ur*drho
                    dlur2 = ur*ur*dl2mu
                    duzdzm = (tz + wvno*xlam(k-1)*ur)/xl2mm
                endif

C       WRITE(6,*)'drho,dmu,dl2mu,drur2,dlur2,duzdzp,durdzp:',
c          drho,dmu,dl2mu,drur2,dlur2,duzdzp,durdzp
                gfac1 =  omega2*drho*uz**2
                gfac2 =  omega2*drur2
                gfac3 = -wvno2*dmu*uz**2
                gfac4 = -wvno2*dlur2
                gfac5 =  (xl2mp*duzdzp**2-xl2mm*duzdzm**2)
                gfac6 =  (xmu(k)*durdzp**2 - xmu(k-1)*durdzm**2)
            endif
C       WRITE(6,*)'fac,gfac1,gfac2,gfac3,gfac4,gfac5,gfac6:',
c          fac,gfac1,gfac2,gfac3,gfac4,gfac5,gfac6
            dfac = fac * (
     1          gfac1 + gfac2 + gfac3 + gfac4
     2          + gfac5 + gfac6 )
            if(dabs(dfac).lt.1.0d-38)then
                dfac = 0.0d+00
            endif
            dcdh(k) = sngl(dfac)
  500   continue
        return
        end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      function ssum(kk,i,j)
c
c     The analytic forms of the solution of integral:
c     Integral U*U dz = T-matrix * eigenfnction * integral-coefs
c
      double precision ssum
      integer*4 NL
      parameter (NL=1000)
      complex*16 t(4,4),tt(4,4),ff(6),pp(4)
      complex*16 sum0,sum1,sum2,sum3,sum4,sum5,sum6
      common/coef/ t,tt,ff,pp
      if(kk.eq.(NL+1)) go to 100
      sum1=t(i,1)*t(j,1)*pp(1)*pp(1)+t(i,3)*t(j,3)*pp(3)*pp(3)
      sum1=sum1*ff(1)
      sum2=t(i,2)*t(j,2)*pp(2)*pp(2)+t(i,4)*t(j,4)*pp(4)*pp(4)
      sum2=sum2*ff(2)
      sum3=(t(i,1)*t(j,2)+t(i,2)*t(j,1))*pp(1)*pp(2)
     *      +(t(i,3)*t(j,4)+t(i,4)*t(j,3))*pp(3)*pp(4)
      sum3=sum3*ff(3)
      sum4=(t(i,1)*t(j,4)+t(i,4)*t(j,1))*pp(1)*pp(4)
     *      +(t(i,3)*t(j,2)+t(i,2)*t(j,3))*pp(2)*pp(3)
      sum4=sum4*ff(4)
      sum5=(t(i,1)*t(j,3)+t(i,3)*t(j,1))*pp(1)*pp(3)
      sum5=sum5*ff(5)
      sum6=(t(i,2)*t(j,4)+t(i,4)*t(j,2))*pp(2)*pp(4)
      sum6=sum6*ff(6)
      sum0=sum1+sum2+sum3+sum4+sum5+sum6
      ssum=(sum0)
      if(dabs(ssum).lt.1.0d-30)sum=0.0d+00
      return
100   continue
      sum1=pp(3)*t(i,3)*t(j,3)*pp(3)*ff(1)
      sum2=pp(4)*t(i,4)*t(j,4)*pp(4)*ff(2)
      sum3=pp(3)*(t(i,3)*t(j,4)+t(i,4)*t(j,3))*pp(4)
      sum3=sum3*ff(3)
      sum0=sum1+sum2+sum3
      ssum=(sum0)
      if(dabs(ssum).lt.1.0d-30)sum=0.0d+00
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine intval(k,nua,nub,dpth)
c
c     This finds the coefficients needed for integrals.
c     see Wang's dissertation eq. III-1-11 (p.80).
c
      implicit double precision (a-h,o-z)
      integer*4 NL
      parameter (NL=1000)
      complex*16 t(4,4),tt(4,4),ff(6),pp(4)
      complex*16 nua,nub,p,q,pq,expp,exqq,pdum
      common/coef/ t,tt,ff,pp
      if(k.eq.(NL+1)) then
            ff(1)=0.5d+00/nua
            ff(2)=0.5d+00/nub
            ff(3)=1.0d+00/(nua+nub)
      else
            p=nua*dpth
            q=nub*dpth
            pq=(nua+nub)*dpth
            call ifpq(p,p+p,expp)
            ff(1)=(1.0d+00-expp)/(2.0d+00*nua)
            call ifpq(q,q+q,exqq)
            ff(2)=(1.0d+00-exqq)/(2.0d+00*nub)
            pdum = pq/2.0d+00
            call ifpq(pdum,pq,expp)
            ff(3)=(1.0d+00-expp)/(nua+nub)
            pdum = p/2.0d+00
            call ifpq(pdum,p,expp)
            pdum = q/2.0d+00
            call ifpq(pdum,q,exqq)
            ff(4)=(exqq-expp)/(nua-nub)
            ff(5)=dpth*expp
            ff(6)=dpth*exqq
      endif
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine ifpq(p,pq,expq)
c
c     finds exp(pq)
c
      implicit double precision (a-h,o-z)
        COMPLEX*16 CDEXP
      complex*16 p,pq,expq
        double precision p0
      p0=dreal(p)
      if(dabs(p0).lt.40.0) go to 100
      expq=0.0d+00
      go to 200
100   continue
      expq=cdexp(-pq)
200   continue
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine tplus(nua,nub,wvno2,gam,gammk,rho)
c
c     E matrix.  B(motion-stress)=E*K(potential-constant)
c
      implicit double precision (a-h,o-z)
      complex*16 t(4,4),tt(4,4),ff(6),pp(4)
      complex*16 nua,nub
      common/coef/ t,tt,ff,pp
      t(1,1)=-1.d+00/rho
      t(1,2)=nub/rho
      t(1,3)=t(1,1)
      t(1,4)=-t(1,2)
      t(2,1)=-nua/rho
      t(2,2)=wvno2/rho
      t(2,3)=-t(2,1)
      t(2,4)=t(2,2)
      t(3,1)=1.d+00-gam
      t(3,2)=gam*nub
      t(3,3)=t(3,1)
      t(3,4)=-t(3,2)
      t(4,1)=-gammk*nua
      t(4,2)=-t(3,1)
      t(4,3)=-t(4,1)
      t(4,4)=t(4,2)
      do 100 i=1,4
      do 100 j=1,4
      t(i,j)=0.5d+00*t(i,j)
100   continue
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine tminus(nua,nub,wvno2,gam,gammk,rho)
c
c     E-inverse matrix
c
      implicit double precision (a-h,o-z)
      complex*16 t(4,4),tt(4,4),ff(6),pp(4)
      complex*16 nua,nub
      common/coef/ t,tt,ff,pp
      gamm1=gam-1.0d+00
      tt(1,1)=-rho*gam
      tt(1,2)=rho*gamm1/nua
      tt(1,3)=1.0d+00
      tt(1,4)=-wvno2/nua
      tt(2,1)=-rho*gamm1/nub
      tt(2,2)=rho*gammk
      tt(2,3)=1.d+00/nub
      tt(2,4)=-1.0d+00
      tt(3,1)=tt(1,1)
      tt(3,2)=-tt(1,2)
      tt(3,3)=1.0d+00
      tt(3,4)=-tt(1,4)
      tt(4,1)=-tt(2,1)
      tt(4,2)=tt(2,2)
      tt(4,3)=-tt(2,3)
      tt(4,4)=-1.0d+00
      return
      end

        subroutine insert(dph,newlyr,ls)
        implicit double precision (a-h,o-z)
        parameter(NL=1000)
        parameter (LER=0, LIN=5, LOT=6)
        logical newlyr
        real*4 dph

        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
        common/wateri/iwat(NL)

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
c       Do not create unnecessary layers, e.g., at 
c          surface and internally
c       However do put in a zero thickness layer at the 
c          base if necessary
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
     1      zfrefp, zfrefs,
     1      xmu, xlam
        common/wateri/iwat(NL)
c-----
c       Find source/receiver boundary. It is assumed that
c       it will lie upon a boundary
c
c       lmax = source layer 
c       depth = source depth 
c       dph = height of  source above lmax + 1 interface 
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

        subroutine gammap(omega,wvno,gamma)
c-----
c       This routine finds the attenuation gamma value.
c
c-----
        implicit double precision (a-h,o-z)
        parameter (NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
        common/wateri/iwat(NL)
        common/eigfun/ uu(NL,4),uu0(4),dcda(NL),dcdb(NL),
     *                 dcdr(NL),dcdh(NL)
        gamma=0.0
        dc = 0.0
        pi = 3.141592653589493d+00
        do 100 i=1,mmax
            if(iwat(i).eq.  0)then
                x=dcdb(i)*zb(i)*zqbi(i)
                gamma = gamma + x
                omgref=2.0*pi*zfrefs(i)
                dc = dc + dlog(omega/omgref)*x/pi
            endif
            x=dcda(i)*za(i)*zqai(i)
            gamma = gamma + x
            omgref=2.0*pi*zfrefp(i)
            dc = dc + dlog(omega/omgref)*x/pi
  100   continue
        c=omega/wvno
        gamma=0.5*wvno*gamma/c
        c=c+dc
        wvno=omega/c
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

        parameter(NL=1000)
        implicit double precision (a-h,o-z)
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*4 vtp,dtp,rtp
        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        ar = 6370.0d0
        tm=sqrt(1.+(c/(2.*ar*om))**2)
            do 20 i=1,mmax
                dcda(i)=dcda(i)*  vtp(i)/(tm**3)
                dcdb(i)=dcdb(i)*  vtp(i)/(tm**3)
                dcdh(i)=dcdh(i)*  dtp(i)/(tm**3)
                dcdr(i)=dcdr(i)*  rtp(i)/(tm**3)
   20       continue
c       write(6,*)'c flat=',c,' csph=',c/tm
        usph = ugr/tm
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
c       This is for Rayleigh Waves
c
c-----
        parameter (NL=1000)
        implicit double precision (a-h,o-z)
        real*4 vtp,dtp,rtp
        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), xlam(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs,
     1      xmu, xlam
        real*8 ar, dr, r0, r1, z0, z1
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
            vtp(i) = sngl((ar/r1 - ar/r0)*(ar/(z1-z0)))
            rtp(i) = (r0*(r0/ar)-r1*(r1/ar))/(zd(i)*2.0d0)
            dtp(i) = sngl(ar/r0)
            r0 = r1
            z0 = z1
   10   continue
c       write(6,*)'vtp:',(vtp(i),i=1,mmax)
c       write(6,*)'rtp:',(rtp(i),i=1,mmax)
c       write(6,*)'dtp:',(dtp(i),i=1,mmax)
c-----
c       at this point the model information is no longer used and
c       will be overwritten
c-----
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
c       dotmp   L   - .true. use file tsdisp96.ray
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
        character hrfile*120, hsfile*120
        logical dotmp, dogam, dderiv
        integer nipar(10)
        logical verbose


        character name*40
        integer mnmarg

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
                else if(name(1:3).eq.'-DA')then
                    dderiv = .true.
                    nipar(6) = 1
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
        integer LER
        parameter (LER=0)
        character ostr*(*)
        if(ostr .ne. ' ')then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'Usage: sregn96 ',
C     1     ' -FHR recdep -FHS srcdep -HS hs -HR hr ' ,
     1      ' -HS hs -HR hr ' ,
     2       ' [-NOQ] [-T] [-DER -DE -DH -DA -DB -DR ] [-?] [-h]'
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
     1  '-T          (default false) use tdisp96.ray not disp96.ray'
        write(LER,*)
     1  '-DER        (default false) output depth dependent values'
        write(LER,*)
     1  '-DE         (default false) output eigenfunctions(depth)'
        write(LER,*)
     1  '-DH         (default false) output DC/DH(depth)'
        write(LER,*)
     1  '-DA         (default false) output DC/DA(depth)'
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

        subroutine winteg(urur,uzuz,urduz,duzduz,
     1              dpth,nua,drho,
     1              omega2,wvomg2,wvno2,
     2              uzup,tzup,urup,uzdn,tzdn,urdn,kk)
        implicit double precision (a-h,o-z)
        complex*16 nua

        complex*16 dp,dpp,dppex,expp
        complex*16 f1, f2, f3
        REAL*8 DREAL
        COMPLEX*16 CDEXP
c-----
c       to convert to Levshin and Yanson form we must
c
c       Ur (LJ) = k Ur(Haskell)
c       Uz (LJ) = Uz (Haskell)
c       Tz (LJ) = omega**2 * Tz (Haskell)
c       Tr (LJ) = k omega**2 * Tr (Haskell)
c
c-----
c       Get integrals for energy terms for fluid
c       
c       urur    R*8 = INT Ur Ur dz
c       uzuz    R*8 = INT Uz Uz dz
c       urduz   R*8 = INT Ur dUz/dz dz
c       duzduz  R*8 = INT dUz/dz dUz/dz dz
c       dpth    R*8 layer thickness
c       nua C*16    complex vertical wavenumber
c       drho    R*8 layer density
c       omega2  R*8 omega squared
c       wvomg2  R*8 k omega omega
c       wvno2   R*8 k*k
c       uzup    R*8 Uz at top of layer
c       urup    R*8 Ur at top of layer, note this is Ur(Haskell)
c       tzup    R*8 Tz at top of layer, note this is Tz(Haskell)
c       uzdn    R*8 Uz at bottom of layer
c       urdn    R*8 Ur at bottom of layer, note this is Ur(Haskell)
c       tzdn    R*8 Tz at bottom of layer, note this is Tz(Haskell)
c       kk      I   if this = NL+1 then we have only the halfspace
c                   contribution
c-----
c       The trick here is to note that (Harkrider, personal 
c          communication)
c       
c       | Uz |   | (Nu/Rho) exp(-Nu z)     -(Nu/Rho) exp(Nu z) | | D' |
c       |    | = |                                             | |    |
c       | Tz | = |      exp(-Nu z)                exp(Nu z)    | | D" |
c
c       D' is downgoing, D" is upgoing
c-----
c
c       Ur =  (-1 / Rho) Tz
c       
c       D' =  (1/2)(Rho/Nu)[ Uz + (Nu/Rho)Tz ]
c                                             top
c
c       D" =  (1/2)(Rho/Nu)[-Uz + (Nu/Rho)Tz ] exp(-Nu H)
c                                                        bottom
c
c       D" exp(Nu H) =  (1/2)(Rho/Nu)[-Uz + (Nu/Rho)Tz ] 
c                                                       bottom
c       for halfspace, there is just exponential decay and then
c       INT  Uz  Uz dz = (Uz^2)    (1/2 Nu)
c                              top
c
c
c-----
      integer*4 NL
      parameter (NL=1000)
        f1=nua*dpth
        if(dreal(f1).lt.40.0) then
            expp=cdexp(-2.d+00*f1)
        else
            expp=dcmplx(0.0d+00,0.0d+00)
        endif
        f3=(1.d+00-expp)/(2.d+00*nua)
        if(dreal(f1).lt.75.0) then
            expp=cdexp(-f1)
        else
            expp=0.0d+00
        endif
        f2=expp
        
        if(kk .eq. (NL+1))then
c-----
c              halfspace
c-----
               uzuz   = uzup*uzup/(2.*nua)
               duzduz = nua*nua*uzup*uzup/(2.*nua)
        else
c-----
c              layer
c-----

               dp = (0.5d+00 * drho/nua)*
     1             (uzup + (nua/drho)*tzup)
               dppex = (0.5d+00 * drho/nua)*
     1             (-uzdn + (nua/drho)*tzdn)
               dpp = (0.5d+00 * drho/nua)*
     1             (-uzdn + (nua/drho)*tzdn)*f2
               uzuz = (nua/drho)*(nua/drho) *
     1             ( dp*dp*f3 - 2.0d+00*dp*dpp*dpth + 
     2          (dppex*dppex - dpp*dpp)/(2.0d+00*nua))
               urur = (-1.0d+00/drho)*(-1.0d+00/drho)*
     1             ( dp*dp*f3 + 2.0d+00*dp*dpp*dpth + 
     2          (dppex*dppex - dpp*dpp)/(2.0d+00*nua))
               duzduz = (-nua*nua/drho)*(-nua*nua/drho)*
     1             ( dp*dp*f3 + 2.0d+00*dp*dpp*dpth + 
     2          (dppex*dppex - dpp*dpp)/(2.0d+00*nua))
               urduz = (-1.0d+00/drho)*(-nua*nua/drho)*
     1             ( dp*dp*f3 + 2.0d+00*dp*dpp*dpth + 
     2          (dppex*dppex - dpp*dpp)/(2.0d+00*nua))
        endif
        return
        end

        subroutine collap(ls,mmaxot)
        parameter(NL=1000)
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0, dcda, dcdb, dcdr, dcdh
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
