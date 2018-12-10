        program slegn96
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: SLEGN96                                                c
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
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c Revision history:
c       07 AUG 2002 - make string lengths 120 characters from 80 
c       13 OCT 2006 - verbose output of energy integrals
c       26 SEP 2008 - fixed undefined LOT in subroutine up
c       14 JUN 2009 - reformulate in general terms as per book
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter(LER=0,LIN=5,LOT=6)
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
        common/model/  d(NL),a(NL),b(NL),rho(NL),qa(NL), qb(NL)
        real*4 d, a, b, rho, qa, qb
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs, xmu

        common/eigfun/ uu(NL),tt(NL),dcdb(NL),dcdh(NL),dcdr(NL),uu0(4)
        real*8 uu, tt, dcdb, dcdh, dcdr, uu0
        real*4 sdcda(NL), sdcdb(NL), sdcdh(NL), sdcdr(NL) 
        real*4 spur(NL), sptr(NL), spuz(NL), sptz(NL)
        common/wateri/iwat(NL)
        common/sumi/   sumi0,sumi1,sumi2,flagr,ale,ugr
        real*8 sumi0, sumi1, sumi2, flagr, ale, ugr
        common/depref/refdep

        real*4 vtp,dtp,rtp
        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        real*4 s1
        parameter(MAXMOD=2000)
        real*8 cp(MAXMOD), t
        real*8 c, omega, wvno, gamma, csph, usph
        real*8 Eut, Edut, Eut0, Ett0, Ed2ut
        logical nwlyrs, nwlyrr
        character*12 fname(2)
c-----
c       wired in file names     - output of this program 
c                           is always in the
c               binary file slegn96.egn or slegn96.der
c                               - input from sdisp96 is always
c               binary file sdisp96.lov
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
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line information
c-----  
        call gcmdln(hsfile,hrfile,hs,hr,dotmp,dogam,dderiv,
     1      nipar,verbose)
        if(dotmp)then
            fname(1) = 'tsdisp96.lov'
        else
            fname(1) = 'sdisp96.lov'
        endif
        if(dderiv)then
            fname(2) = 'slegn96.der'
        else
            fname(2) = 'slegn96.egn'
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
c       see if the file sdisp96.lov exists
c-----
        inquire(file=fname(1),exist=ext)
        if(.not.ext)then
            call usage('Dispersion file: '//fname(1)//
     1          ' does not exist')
        endif
c-----
c       open output file of sdisp96 and of slegn96
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
                call bldsph(mmax)
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
                call shfunc(omega,wvno)
            call energy(omega,wvno,Eut,Edut,Ed2ut,Eut0,Ett0,
     1          lss,lrr)
c-----
c       the gamma routine will use the spherical model, but the
c       frequency dependence and Q of the original model
c-----
C        WRITE(6,*)'c(noq)=',c
            if(dogam)then
                call gammaq(omega,wvno,gamma)
                c = omega/wvno
C        WRITE(6,*)'c(q)=',c
            else
                gamma = 0.0d+00
            endif
c-----
c       output necessary eigenfunction values for
c       source excitation
c-----
            mmaxot = mmax
                if(nsph.gt.0)then
c-----
c               sphericity correction for partial derivatives
c               of the original model
c-----
C        WRITE(6,*)'c(flat)=',c,' u(flat)=',ugr
                 call splove(omega,c,mmaxot,csph,usph,ugr)
                 wvno = omega/csph
C        WRITE(6,*)'c(sph )=',csph,' u(sph )=',usph
                endif
c-----
c           check for underflow
c-----
            if(dabs(Eut).lt. 1.0d-36)Eut = 0.0d+00 
            if(dabs(Edut).lt. 1.0d-36)Edut = 0.0d+00 
            if(dabs(Ed2ut).lt. 1.0d-36)Ed2ut = 0.0d+00 
            if(dabs(Eut0).lt. 1.0d-36)Eut0 = 0.0d+00 
            if(dabs(Ett0).lt. 1.0d-36)Ett0 = 0.0d+00 
            if(dabs(ugr).lt.1.0d-36)ugr=0.0d+00
            if(dabs(sumi0).lt.1.0d-36)sumi0=0.0d+00
            if(dabs(sumi1).lt.1.0d-36)sumi1=0.0d+00
            if(dabs(sumi2).lt.1.0d-36)sumi2=0.0d+00
            if(dabs(ale).lt.1.0d-36)ale=0.0d+00
            if(dabs(flagr).lt.1.0d-36)flagr=0.0d+00
            if(dabs(gamma).lt.1.0d-36)gamma=0.0d+00

c-----
c      get the derivatives of the eitenfunctions required
c      for source excitation from the definition of stress. For
c      completeness get the second derivative from the first
c      derivatives and the equation of motion for the medium
c-----
            sut = sngl(Eut)
            sdut = sngl(Edut)
            sd2ut = sngl(Ed2ut)
            suz = 0.0
            sduz = 0.0
            sd2uz = 0.0
            sale = sngl(ale)
            wvnsrc = sngl(wvno)
            sur0 = 0.0

            if( verbose ) then
                 WRITE(LOT,2)t,c,ugr,gamma,
     1                sumi0,sumi1,sumi2,
     2                flagr,ale
    2    format(' T=',e15.7,'  C=',e15.7,'  U=',e15.7,' G=',e15.7/
     1          'I0=',e15.7,' I1=',e15.7,' I2=',e15.7/
     2          ' L=',e15.7,' AL=',e15.7)
C                 WRITE(LOT,'(3e26.17/3e26.17/3e26.17)')t,c,ugr,gamma,
C     1                sumi0,sumi1,sumi2,
C     2                flagr,ale
            endif

            rut = sngl(Eut0)
            rtt = sngl(Ett0)
            ruz = 0.0
            rtz = 0.0
            rale = sngl(ale)
            wvnrec = sngl(wvno)
            rur0 = 0.0

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
C            do i = 1,mmax
C            WRITE(6,'(i5,6f10.6)')
C    1            I,zd(i),uu(i),tt(i),dcdh(i),dcdb(i),dcdr(i)
C            enddo
C            WRITE(6,*)'nwlyrr,lrro:',nwlyrr,lrro
C            WRITE(6,*)'nwlyrs,lsso:',nwlyrs,lsso
            if(nwlyrr)then
                call collap(lrro+1,mmaxot)
            endif
            if(nwlyrs)then
                call collap(lsso+1,mmaxot)
            endif
            call chksiz(uu,spur,mmaxot)
            call chksiz(tt,sptr,mmaxot)
            call chksiz(dcdh,sdcdh,mmaxot)
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
C           WRITE(6,*)'mmaxot:',mmaxot
C           do i = 1,mmaxot
C            WRITE(6,'(i5,5f10.6)')
C    1            I,uu(i),tt(i),sdcdh(i),sdcdb(i),sdcdr(i)
C           enddo

            call putder(2,5,sngl(wvno),sngl(ugr), 
     1          sngl(gamma), 
     1          sut,sdut,sd2ut,suz,sduz,sd2uz,sale,wvnsrc,sur0,
     2          rut,rtt,ruz,rtz,rale,wvnrec,rur0,
     3          sumkr,sumgr,sumgv,mmaxot,
     4          sdcdh,sdcda,sdcdb,sdcdr,
     5          spur,sptr,spuz,sptz,ipar)
        else
            call putegn(2,1,1,sngl(wvno),sngl(ugr),
     1          sngl(gamma),
     1          sut,sdut,suz,sduz,sale,wvnsrc,sur0,
     2          rut,rtt,ruz,rtz,rale,wvnrec,rur0,
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

        subroutine shfunc(omega,wvno)
c-----
c       This routine evaluates the eigenfunctions by calling sub up.
c-----
        implicit double precision (a-h,o-z)
        parameter(NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs, xmu
        common/wateri/iwat(NL)
        common/eigfun/ uu(NL),tt(NL),dcdb(NL),dcdh(NL),dcdr(NL),uu0(4)
        common/save/   exl(NL)
c-----
c       get the eigenfunctions in an extended floating point form
c       by propagating from the bottom upward. Note since we move in this
c       direction the propagator matrix A(d) is evaluated as A(-d)
c-----
        call up(omega,wvno,fl)
        uu0(1)=1.0
c-----
c       uu0(2)=stress0 is actually the value of period equation.
c       uu0(3) is used to print out the period equation value before
c       the root is refined.
c-----
        uu0(2)=fl
        uu0(3)=0.0
        uu0(4)=0.0
c-----
c       convert to actual displacement from the extended floating point
c       working top to down, where the amplitude should be small
c       Also normalize so that the V(0) = 1
c-----
        ext=0.0
        umax = uu(1)
        tt(1) = 0.0
        do 100 k=2,mmax
            if(iwat(k).eq.0)then
                ext=ext+exl(k-1)
                fact=0.0
                if(ext.lt.80.0) fact=1./dexp(ext)
                uu(k)=uu(k)*fact
                tt(k)=tt(k)*fact
            else
                uu(k) = 0.0
                tt(k) = 0.0
            endif
            if(abs(uu(k)).gt.abs(umax))then
                umax = uu(k)
            endif
  100   continue
        if(uu(1).ne.0.0)then
                umax = uu(1)
        endif
        if(abs(umax).gt.0.0)then
            do 200 k=1,mmax
                if(iwat(k).eq.0)then
                    uu(k) = uu(k) / umax
                    tt(k) = tt(k) / umax
                endif
  200       continue
        endif
c-----
c       force boundary condition at the free surface
c       free surface is stress free
c-----
        return
        end

        subroutine up(omega,wvno,fl)
c-----
c       This routine calculates the elements of Haskell matrix,
c       and finds the eigenfunctions by analytic solution.
c-----
c       History:
c       
c       16 JUN 2009 - generalized the routine so that  the
c         nature of the propagator matrices is buried in the
c         subroutine hskl
c         This will permit simpler extension to TI
c-----
        implicit double precision (a-h,o-z)
        integer LER, LIN, LOT
        parameter(LER=0,LIN=5,LOT=6)
        parameter(NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs, xmu
        common/wateri/iwat(NL)
        common/eigfun/ uu(NL),tt(NL),dcdb(NL),dcdh(NL),dcdr(NL),uu0(4)
        common/save/   exl(NL)
        real*8 a11, a12, a21, a22

        real*8 rb, omega, wvno, xkb, dpth, exl
        real*8 hl(2,2)
c
c-----
c       apply the boundary conditions to the halfspace
c-----
c-----
c     kludge for fluid core
c-----
        if(zb(mmax).gt.0.01)then
                dpth = 0.0
                call varl(mmax,rb,omega,wvno,xkb,dpth,eexl)
               if(wvno.lt.xkb)then
                    write(LOT,*) ' imaginary nub derivl'
                    write(LOT,*)'omega,wvno,b(mmax)',omega,wvno,zb(mmax)
                 endif
               uu(mmax)=1.0
               tt(mmax)=-xmu(mmax)*rb
        else
               uu(mmax)=1.0
               tt(mmax)=0.0
        endif
        exl(mmax) = 0.0
        mmx1=mmax-1
        do 500 k=mmx1,1,-1
            if(iwat(k).eq.0)then
                dpth=zd(k)
                call varl(k,rb,omega,wvno,xkb,dpth,eexl)
                call hskl(hl,iwat(k))
                k1=k+1
c-----
c               We actually use A^-1 since we go from the bottom to top
c-----
                a11 = hl(1,1)
                a22 = hl(2,2)
                a12 = - hl(1,2)
                a21 = - hl(2,1)
C            WRITE(6,*)'eexl:',eexl
C            WRITE(6,*)'hl:',hl

   
                amp0 = a11 * uu(k1) + a12*(tt(k1))
                str0 = a21 * uu(k1) + a22*(tt(k1))
c-----
c               now normalize and save the normalization
c-----
                rr=dabs(amp0)
                ss=dabs(str0)
                if(ss.gt.rr) rr=ss
                if(rr.lt.1.d-30) rr=1.d+00
                exl(k)=dlog(rr)+eexl
                uu(k)=amp0/rr
                tt(k)=str0/rr
                ttlast = tt(k)
            endif
  500   continue
        fl=ttlast
        return
        end

        subroutine energy(omega,wvno,Eut,Edut,Ed2ut,Eut0,Ett0,
     1      lss,lrr)
c-----
c       This routine calculates the values of integrals I0, I1,
c       and I2 using analytic solutions. It is found
c       that such a formulation is more efficient and practical.
c
c       This is based on a suggestion by David Harkrider in 1980.
c       Given the eigenfunctions, which are properly defined,
c       define the potential coefficients at the bottom and
c       top of each layer. If the V(z) = A exp (nu z) + B exp (-nu z)
c       We want the B at the top of the layer for z=0, and the A at the
c       bottom of the layer, from the relation K=[Kup Kdn]^T
c       = [A B]^T = E^-1 U
c-----
c       History:
c       
c       16 JUN 2009 - generalized the routine so that  the
c         nature of the propagator matrices is buried in the
c         subroutine emat which provides E and E^-1
c         This will permit simpler extension to TI
c-----
        implicit double precision (a-h,o-z)
        integer LER, LIN, LOT
        complex*16 nub,xnub,exqq,top,bot,f1,f2,f3
        parameter(NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs, xmu
        common/wateri/iwat(NL)
        common/eigfun/ uu(NL),tt(NL),dcdb(NL),dcdh(NL),dcdr(NL),uu0(4)
        common/sumi/   xi0,xi1,xi2,flagr,ale,ugr
        COMPLEX*16 CDEXP
        complex*16 kmup, km1dn
        complex*16 e11,e12,e21,e22
        complex*16 einvl(2,2), el(2,2)
        complex*16 f, g
        real*8 TN, TL
        real*8 DCDBH, DCDBV, DCDRSH
        real*8 VSHH, VSHV

        c=omega/wvno
        omega2=omega*omega
        wvno2=wvno*wvno
        sumi0=0.0d+00
        sumi1=0.0d+00
        sumi2=0.0d+00
c-----
c RBH beware when nub = 0
c-----
        do 300 k=1,mmax
            if(iwat(k).eq.0)then
            TN = zrho(k)*zb(k)*zb(k)
            TL = zrho(k)*zb(k)*zb(k)
            VSHH = zb(k)
            VSHV = zb(k)
            k1=k+1
            drho=zrho(k)
            dpth=zd(k)
            call varl(k,rb,omega,wvno,xkb,dpth,eexl)
            dmu=xmu(k)
            if(rb .lt. 1.0d-10)rb = 1.0d-10
c-----
c      for safety do not permit rb = 0
c-----
            if(k.eq.mmax) then
                upup  =(0.5d+00/rb)*uu(mmax)*uu(mmax)
                dupdup=(0.5d+00*rb)*uu(mmax)*uu(mmax)
            else
c-----
c               use Einv to get the potential coefficients
c               from the U and T
c-----
                nub=dcmplx(rb,0.0d+00)
                if(wvno.lt.xkb) nub=dcmplx(0.0d+00,rb)
C                WRITE(6,*)'nub:',nub
                xnub=dmu*nub
                call emat(xnub,wvno,el,einvl)
                km1dn = einvl(2,1)*uu(k)  + einvl(2,2)*tt(k)
                kmup  = einvl(1,1)*uu(k1) + einvl(1,2)*tt(k1)

c-----
c               get the f function
c-----
                f3=nub*dpth
                exqq=0.0d+00
                if(dreal(f3).lt.40.0) exqq=cdexp(-2.d+00*f3)
                f=(1.d+00-exqq)/(2.d+00*nub)
c-----
c               get the g function
c-----
                exqq=0.0d+00
                if(dreal(f3).lt.75.0) exqq=cdexp(-f3)
                g=dpth*exqq

                f1 = f*(el(1,1)*el(1,1)*kmup*kmup 
     1              + el(1,2)*el(1,2)*km1dn*km1dn)
                f2 = g*(el(1,1)*el(1,2)+el(1,1)*el(1,2))*kmup*km1dn
c-----
c               cast to a real  upup
c-----
                upup = f1 + f2

c-----
c               cast to a real dupdup
c-----
                dupdup = nub*nub*(f1 - f2)

c-----
            endif
C            WrItE(*,*)k,upup,dupdup
            sumi0=sumi0+drho*upup
            sumi1=sumi1+TN*upup
            sumi2=sumi2+TL*dupdup
            dcr=-0.5d+00*c*c*c*upup
            dcb=0.5d+00*c*(upup+dupdup/wvno2)
            DCDBH = c*drho*VSHH*upup
            DCDBV = c*drho*VSHV*dupdup/wvno2
            dcdb(k)=DCDBH + DCDBV
            DCDRSH= 0.5*c*(-c*c*upup + VSHH*VSHH*upup +
     1       VSHV*VSHV*dupdup/wvno2)
            dcdr(k)=DCDRSH
            else
                  dcdb(k)=0.0
                  dcdr(k)=0.0
            endif
  300   continue
        do 400 k=1,mmax
            if(iwat(k).eq.0)then
                dcdb(k)=dcdb(k)/sumi1
                dcdr(k)=dcdr(k)/sumi1
            else
                dcdb(k)=0.0
                dcdr(k)=0.0
            endif
  400   continue
        flagr=omega2*sumi0-wvno2*sumi1-sumi2
        ugr=sumi1/(c*sumi0)
        ale=0.5d+00/sumi1
        xi0=sumi0
        xi1=sumi1
        xi2=sumi2
c-----
c       define partial with respect to layer thickness
c-----
c       fac = 0.5d+00*c**3/(omega2*sumi1)
        fac = ale*c/wvno2
        c2 = c*c
        llflag = 0
        do 500 k=1,mmax
        if(iwat(k).eq.0)then
            if(llflag.eq.0)then
                drho = zrho(k)
                dmu  = xmu(k)
                dvdz = 0.0
            else 
                drho = zrho(k) - zrho(k-1)
                dmu  = xmu(k) - xmu(k-1)
                dvdz = tt(k)*tt(k)*(1.0/xmu(k) - 1.0/xmu(k-1))
            endif
            dfac = fac * ( uu(k)*uu(k)*
     1          (omega2*drho - wvno2*dmu) + dvdz)
            if(dabs(dfac).lt.1.0d-38)then
                dcdh(k) = 0.0
            else
                dcdh(k) = sngl(dfac)
            endif
            llflag = llflag + 1
        else
            dcdh(k) = 0.0
        endif
  500   continue
c-----
c       compute the eigenfuntions and depth derivatives
c       at the source depth
c-----
        if(iwat(lss).eq.0)then
            Eut = uu(lss)
            Edut = tt(lss)/xmu(lss)
            Ed2ut = ( - (omega/zb(lss))**2 + wvno2)*Eut
        else
            Eut = 0.0d+00
            Edut = 0.0d+00
            Ed2ut = 0.0d+00
        endif
        if(iwat(lrr).eq.0)then
            Eut0 = uu(lrr)
            Ett0 = tt(lrr)
        else
            Eut0 = 0.0d+00
            Ett0 = 0.0d+00
        endif
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
     2      xmu(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs, xmu
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
     2      xmu(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs, xmu
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

        subroutine gammaq(omega,wvno,gamma)
c-----
c       This routine finds the attenuation gamma value.
c
c-----
        implicit double precision (a-h,o-z)
        parameter (NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs, xmu
        common/wateri/iwat(NL)
        common/eigfun/ ut(NL),tt(NL),dcdb(NL),dcdh(NL),dcdr(NL),uu0(4)
        gamma=0.0
        dc = 0.0
        pi = 3.141592653589493d+00
        do 100 i=1,mmax
            if(iwat(i).eq.0)then
                x=dcdb(i)*zb(i)*zqbi(i)
                gamma = gamma + x
                omgref=2.0*pi*zfrefs(i)
                dc = dc + dlog(omega/omgref)*x/pi
            endif
  100   continue
        c=omega/wvno
        gamma=0.5*wvno*gamma/c
        c=c+dc
        wvno=omega/c
        return
        end

        subroutine splove(om,c,mmax,csph,usph,ugr)
c-----
c       Transform spherical earth to flat earth
c       and relate the corresponding flat earth dispersion to spherical
c
c       Schwab, F. A., and L. Knopoff (1972). Fast surface wave 
c            and free
c       mode computations, in  
c            Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c            B. A. Bolt (ed),
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
        real*8 om, tm
        real*8 c, usph, csph, ugr
        parameter(NL=1000)
        common/eigfun/ ut(NL),tt(NL),dcdb(NL),dcdh(NL),dcdr(NL),uu0(4)
        real*8 ut, tt, dcdb, dcdh, dcdr, uu0
        real*4 vtp,dtp,rtp
        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        a = 6371.0d0
        tm=sqrt(1.+(3.0*c/(2.*a*om))**2)
            do 20 i=1,mmax
                dcdb(i)=dcdb(i)*  vtp(i)/(tm**3)
                dcdh(i)=dcdh(i)*  dtp(i)/(tm**3)
                dcdr(i)=dcdr(i)*  rtp(i)/(tm**3)
   20       continue
        csph = c / tm
        usph = ugr * tm
c       write(6,*)'c flat=',c,' csph=',c/tm
        return
        end

        subroutine bldsph(mmax)
c-----
c       Transform spherical earth to flat earth
c
c       Schwab, F. A., and L. Knopoff (1972). 
c            Fast surface wave and free
c       mode computations, 
c            in  Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c            B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  44, 45 , 41 pp 112-113
c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c     Revised 28 DEC 2007 to use mid-point, assume linear variation in
c     slowness instead of using average velocity for the layer
c     Use the Biswas (1972:PAGEOPH 96, 61-74, 1972) density mapping
c
c       This is for Love Waves
c
c-----
        parameter (NL=1000)
        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        real*4 vtp,dtp,rtp
        common/model/  d(NL),a(NL),b(NL),rho(NL),qa(NL), qb(NL)
        real*4 d, a, b, rho, qa, qb
        real*8 ar, dr, r0, r1, z0, z1
c-----
c       vtp is the factor used to convert 
c            spherical velocity to flat velocity
c       rtp is the factor used to convert 
c            spherical density  to flat density 
c       dtp is the factor used to convert 
c            spherical boundary to flat boundary
c-----
c-----
c       duplicate computations of srwvds(IV)
c-----
        ar=6371.0d0
        dr=0.0d0
        r0=ar
        z0 = 0.0d+00
        d(mmax)=1.0
        do 10 i=1,mmax
            r1 = r0 * dexp(dble(-d(i))/ar)
            if(i.lt.mmax)then
                z1 = z0 + dble(d(i))
            else
                z1 = z0 + 1.0d+00
            endif
            TMP=(ar+ar)/(r0+r1)
             vtp(i) = TMP
             rtp(i) = TMP**(-5)
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
                    nipar(6) = 0
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
        integer LER
        parameter (LER=0)
        character ostr*(*)
        if(ostr .ne. ' ')then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'Usage: slegn96 ',
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
        parameter(NL=1000)
        common/eigfun/ uu(NL),tt(NL),dcdb(NL),dcdh(NL),dcdr(NL),uu0(4)
        real*8 uu, tt, dcdb, dcdh, dcdr, uu0
        do 501 i = ls-1,mmaxot
            if(i .eq. ls -1)then
                dcdb(i) = dcdb(i) + dcdb(i+1)
                dcdr(i) = dcdr(i) + dcdr(i+1)
            endif
            if(i.gt.ls)then
                dcdb(i-1) = dcdb(i)
                dcdh(i-1) = dcdh(i)
                dcdr(i-1) = dcdr(i)
                uu(i-1) = uu(i)
                tt(i-1) = tt(i)
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

        subroutine varl(m,rb,omega,wvno,xkb,dpth,eexl)
c-----
c     find variables cosQ,  sinQ, etc.
c-----
c     To handle the hyperbolic functions correctly for large
c     arguments, we use an extended precision procedure,
c     keeping in mind that the maximum precision in double
c     precision is on the order of 16 decimal places.
c
c     So  cosq = 0.5 ( exp(+q) + exp(-q))
c              = exp(q) * 0.5 * ( 1.0 + exp(-2q) )
c     becomes
c         cosq = 0.5 * (1.0 + exp(-2q) ) with an exponent q
c     In performing matrix multiplication, we multiply the modified
c     cosp terms and add the exponents. At the last step
c     when it is necessary to obtain a true amplitude,
c     we then form exp(p). For normalized amplitudes at any depth,
c     we carry an exponent for the numerator and the denominator, and
c     scale the resulting ratio by exp(NUMexp - DENexp)
c
c     The propagator matrices have three basic terms
c
c     HSKA        cosq  sinq
c
c     When the extended floating point is used, we use the
c     largest exponent for each, which is  the following:
c
c     Let sex = s exponent > 0 for evanescent waves = 0 otherwise
c
c     Then the modified matrix elements are as follow:
c
c     Haskell:  cosq -> 0.5 ( 1 + exp(-2q) ) exponent = pex
c-----
      implicit double precision (a-h,o-z)
      integer m
      real*8 rb, xkb, omega, wvno, dpth, eexl
      integer*4 NL
      parameter (NL=1000)
        common/dmodl/  zd(NL),za(NL),zb(NL),zrho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL),
     2      xmu(NL), mmax
        real*8 zd, za, zb, zrho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs, xmu
        common/lvar/cosq,sinq,yl,zl,mu
        real*8 cosq,sinq,yl,zl,mu

      real*8 wvnop, wvnom     , fac

c-----
c     define the  horizontal wavenumber for the S-wave
c-----
      xkb = omega/zb(m)
c-----
c     define the vertical wavenumber for the given wvno
c-----
      wvnop = wvno + xkb
      wvnom = dabs(wvno - xkb)
      fac = wvnop*wvnom
C      WrItE(*,*)'fac:',m,fac
      rb=dsqrt(wvnop*wvnom)
      q = rb*dpth
      mu = zrho(m)*zb(m)*zb(m)
c-----
c     examine S-wave eigenfunctions
c        checking whether c > vs, c = vs, c < vs
c-----
      eexl = 0.0
      if(wvno.lt.xkb)then
             sinq=dsin(q)
             yl=sinq/rb
             zl=-rb*sinq
             cosq=dcos(q)
      else if(wvno.eq.xkb)then
             cosq=1.0d+00
             yl=dpth
             zl=0.0d+00
      else if(wvno.gt.xkb)then
             eexl = q
             fac = 0.0d+00
             if(q.lt.18.0d+00)fac = dexp(-2.0d+0*q)
             cosq = ( 1.0d+00 + fac ) * 0.5d+00
             sinq = ( 1.0d+00 - fac ) * 0.5d+00
             yl = sinq/rb
             zl = rb*sinq
      endif

      return
      end

      subroutine hskl(hl,iwat)
      implicit none
c-----
c     procedure arguments
c-----
      real*8 hl(2,2)
      integer iwat
c-----
c     common blocks
c-----
      common/lvar/cosq,sinq,yl,zl,mu
      real*8 cosq,sinq,yl,zl,mu

      if(iwat.eq.0)then
            hl(1,1) = cosq
            hl(1,2) = yl/mu
            hl(2,1) = zl*mu
            hl(2,2) = cosq
      else
            hl(1,1) = 1.0d+00
            hl(1,2) = 0.0d+00
            hl(2,1) = 0.0d+00
            hl(2,2) = 1.0d+00
      endif
      return
      end

      subroutine emat(xnub,wvno,el,einvl)
      implicit none
c-----
c     compute the E and E^-1 matrices
c-----
      real*8 wvno
      complex*16 xnub, el(2,2), einvl(2,2)
        
          einvl(1,1) =  0.5/wvno
          einvl(1,2) =  0.5/(wvno*xnub)
          einvl(2,1) =  0.5/wvno
          einvl(2,2) = -0.5/(wvno*xnub)
          el(1,1) =  wvno
          el(1,2) =  wvno
          el(2,1) =  wvno*xnub
          el(2,2) = -wvno*xnub
      return
      end
