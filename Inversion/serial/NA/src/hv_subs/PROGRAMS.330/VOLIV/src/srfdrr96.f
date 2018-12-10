c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME IV                                                       c
c                                                                      c
c      PROGRAM: SRFDRR                                                 c
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
c
      program srfdrr96
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
c     Layer thickness partial derivatives added following 
c         the formalism of
c-----
c     Changes:
c     09 FEB 2011 - compute duda
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      parameter(NL=200,NL2=NL+NL)
      parameter(LER=0,LIN=5,LOT=6)
c-----
c     LIN   - unit for FORTRAN read from terminal
c     LOT   - unit for FORTRAN write to terminal
c     LER   - unit for FORTRAN error output to terminal
c     NL    - number of columns in model (first NL2/2 are
c           - velocity parameters, second NL2/2 are Q values)
c-----
        implicit double precision (a-h,o-z)
        real*4 d,a,b,rho,qa1,qb1
        common/model/  d(NL),a(NL),b(NL),rho(NL),qa1(NL),
     *                 qb1(NL),xmu(NL),xlam(NL),mmax,ll
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL2),dcdb(NL2),dcdr(NL)
        common/water/ urb,dphw(NL),water0,dphw0,kw
        common/sumi/   sumi0,sumi1,sumi2,sumi3,flagr,are,ugr
        real*4 btp,dtp
        common/sphere/btp(NL),dtp(NL),dudb(NL2),duda(NL2)
        real*8 dcdbsv(NL2),dcdasv(NL2)
        real*4 h,t,t1,c,cn
        real*4 st0
        character*10 fname(2)
        data fname/'tmpsrfi.07','tmpsrfi.08'/
c-----
c     machine dependent initialization
c-----
        call mchdep()
c-----
        ip=1
        ig=2
        open(1,file=fname(1),form='unformatted',access='sequential')
        rewind 1
        open(2,file=fname(2),form='unformatted',access='sequential')
        rewind 2
c-----
c       obtain the earth model:
c-----
        rewind 1
        read(1) mmax,nsph
        mmax2 = mmax + mmax
        read(1)(btp(i),i=1,mmax)
        read(1)(dtp(i),i=1,mmax)
        do 185 i=1,mmax
            read(1) d(i),a(i),b(i),rho(i)
            xmu(i)=rho(i)*b(i)*b(i)
            xlam(i)=rho(i)*(a(i)*a(i)-2.*b(i)*b(i))
  185   continue
c-----
c       We will merge the output of the Love and Rayleigh partials
c       onto the file tmpsrfi.08 by first rewinding an 
c           existing tmpsrfi.08,
c       then writing the tmpsrfi.05 from srfdrl96 onto it, 
c           and then adding the
c       Rayleigh values
c-----
        call getlove(2,mmax)
        if(b(1).le.0.0)then
            kw = 1
            ll = kw + 1
            dphw(1) = 0.0
            dphw(2) = d(1)
            dphw0   = d(1)
        else
            kw = 0
            ll = kw + 1
        endif
        read(1) nper,igr,h
c-----
c       nper    number of frequencies (periods)
c       igr 0 phase velocity only
c           1 group velocity only
c           2 phase and group velocity data
c       h   increment for period to get partial with respect to period
c-----
        if(igr.ge.2) then
            open(4,file='tmpsrfi.09',form='unformatted',
     1          access='sequential')
            rewind 4
        endif
c-----
c     read in the dispersion values.
c-----
  400   continue
        read(1,end=700) itst,mode,t,t1,c,cn
        if(itst.eq.0)go to 400
c-----
c       itst    0 - higher mode or wave type does not exist 
c           but make dummy
c               entry
c       mode    surface wave mode
c       t   period
c       t1  slightly different period for partial
c       c   phase velocity at t
c       cn  phase velocity at tn
c-----
        t0=t
        if(itst.gt.0) then
            if(igr.gt.0) t0=t*(1.+h)
c-----
c       main part.
c-----
            twopi=2.*3.141592654
            om=twopi/t0
            omega=twopi/t
            wvno=omega/c
            call svfunc(omega,wvno)
            call energy(omega,wvno)
            if(igr.gt.0)then
                cp=c
                omp=omega
                ugp=ugr
                omega=twopi/t1
                omn=omega
                c=cn
                wvno=omega/c
c-----
c       save previous results
c-----
                do 420 i=1,mmax2
                    dcdasv(i)=dcda(i)
                    dcdbsv(i)=dcdb(i)
  420           continue
                call svfunc(omega,wvno)
                call energy(omega,wvno)
                ugr=(ugr+ugp)/2.
                c=(cp+cn)/2.
                do 430 i=1,mmax2
                    dcdn=dcdb(i)
                    dcdb(i)=(dcdbsv(i)+dcdn)/2.
                    uc1=ugr/c
                    ta=uc1*(2.-uc1)*dcdb(i)
                    tb=uc1*uc1*((dcdbsv(i)-dcdn)/(2.*h))
                    dudb(i)=ta+tb
c-----
c     added 09 FEB 2011
c-----
                    dcdn=dcda(i)
                    dcda(i)=(dcda(i)+dcdasv(i))/2.
                    ta=uc1*(2.-uc1)*dcda(i)
                    tb=uc1*uc1*((dcdasv(i)-dcdn)/(2.*h))
                    duda(i)=ta+tb
c-----
c     end new
c-----
  430           continue
            endif
c-----
c       sphericity correction
c-----
            if(nsph.gt.0)then
                if(igr.eq.0)then
                    call sprayl(om,c,ugr,mmax2,0)
                else if(igr.eq.1)then
                    call sprayl(om,c,ugr,mmax2,1)
                else if(igr.eq.2)then
                    call sprayl(om,c,ugr,mmax2,0)
                    call sprayl(om,c,ugr,mmax2,1)
                endif
            endif
        endif
c-----
c       output the derivatives.
c       the first mmax elements are partial with respect to velocity
c       the next mmax are the partial with respect to moving the
c       boundary, e.g., takes two to change layer thickness
c-----
        m=mode
        if(igr.eq.0 .or. igr.eq.2)then
                st0 = sngl(t0)
                write(2) itst,ip,mode,st0
            if(itst.ne.0) then
                write(2) c,(sngl(dcdb(i)),i=1,mmax)
                write(2) (sngl(dcdb(i)),i=mmax+1,mmax2)
                write(2) (sngl(dcda(i)),i=1,mmax)
            
            endif
            if(igr.eq.2)then
                write(4) itst,ig,mode,t0
                if(itst.ne.0) then
                    write(4) ugr,(dudb(i),i=1,mmax)
                    write(4) (dudb(i),i=mmax+1,mmax2)
                    write(4) c,(dcdb(i),i=1,mmax)
                    write(4) (dcdb(i),i=mmax+1,mmax2)
                    write(4) (dcda(i),i=1,mmax)
                    write(4) (duda(i),i=1,mmax)
                endif
            endif
        else if(igr.eq.1)then
                st0 = sngl(t0)
                write(2) itst,ig,mode,st0
            if(itst.ne.0) then
                write(2) sngl(ugr),(sngl(dudb(i)),i=1,mmax)
                write(2) (sngl(dudb(i)),i=mmax+1,mmax2)
                write(2) c,(sngl(dcdb(i)),i=1,mmax)
                write(2) (sngl(dcdb(i)),i=mmax+1,mmax2)
                write(2) (sngl(dcda(i)),i=1,mmax)
                write(2) (sngl(duda(i)),i=1,mmax)
            endif
        endif
        go to 400
  700   continue
c-----
c       end of data read and processing
c       do final clean up, outputing group velocity stuff
c       after all the phase velocity stuff
c-----
        if(igr.ge.2) then
            rewind 4
            do 800 i=1,9000
                read(4,end=900) itst,ig,mode,t0
                st0 = sngl(t0)
                write(2) itst,ig,mode,st0
                if(itst.ne.0) then
                    read(4) ugr,(dudb(j),j=1,mmax)
                    read(4) (dudb(j),j=mmax+1,mmax2)
                    read(4) c,(dcdb(j),j=1,mmax)
                    read(4) (dcdb(j),j=mmax+1,mmax2)
                    read(4) (dcda(j),j=1,mmax)
                    read(4) (duda(j),j=1,mmax)

                    write(2) sngl(ugr),(sngl(dudb(j)),j=1,mmax)
                    write(2) (sngl(dudb(j)),j=mmax+1,mmax2)
                    write(2) c,(sngl(dcdb(j)),j=1,mmax)
                    write(2) (sngl(dcdb(j)),j=mmax+1,mmax2)
                    write(2) (sngl(dcda(j)),j=1,mmax)
                    write(2) (sngl(duda(j)),j=1,mmax)

                endif
  800       continue
  900       close(4,status='delete')
        endif
c-----
c       enter and indicator that this is the end of the data set
c-----
C       itst = 0
C       mode = -1
C       t0 = 0.0
C       write(2) itst,ig,mode,t0
        do 950 i=1,2
            close(i,status='keep')
  950   continue
c       stop
        end

        subroutine sprayl(om,c,u,mmax2,iflag)
c-----
c       Transform spherical earth to flat earth
c       and relate the corresponding flat earth dispersion to spherical
c
c       Schwab, F. A., and L. Knopoff (1972). Fast surface wave and free
c       mode computations, in  Methods in Computational Physics, 
c           Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c           B. A. Bolt (ed),
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
c       u       R*8     group velocity
c       mmax2    I*4     number of layers* 2, first mmax/2 values are
c                       partial with respect to velocity, second are
c                       partial with respect to layer thickness
c       iflag   I*4     0 - phase velocity
c                       1 - group velocity
c-----

        parameter(NL=200,NL2=NL+NL)
        implicit double precision (a-h,o-z)
        real c
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *      dcda(NL2),dcdb(NL2),dcdr(NL)
        real*4 btp,dtp
        common/sphere/btp(NL),dtp(NL),dudb(NL2),duda(NL2)
        a=6370.0
        mmax = mmax2/2
        rval = a
        tm=sqrt(1.+(c/(2.*a*om))**2)
        if(iflag.eq.1) then
            do 10 i=1,mmax
                tmp=dudb(i)*tm+u*c*dcdb(i)/(tm*(2.*a*om)**2)
                dudb(i)=btp(i)*tmp
                tmp=dudb(i+mmax)*tm+u*c*dcdb(i+mmax)
     1              /(tm*(2.*a*om)**2)
                dudb(i+mmax) = (a/rval)*tmp
                rval = rval - dtp(i)
c-----
c     new 09 FEB 2011
c-----
                tmp=duda(i)*tm+u*c*dcda(i)/(tm*(2.*a*om)**2)
                duda(i)=btp(i)*tmp
                tmp=duda(i+mmax)*tm+u*c*dcda(i+mmax)
     1              /(tm*(2.*a*om)**2)
                duda(i+mmax) = (a/rval)*tmp
c-----
c    end new
c-----
   10       continue
            u=u*tm
        else
            do 20 i=1,mmax
                dcdb(i)=dcdb(i)*btp(i)/(tm**3)
                dcdb(i+mmax)=dcdb(i+mmax)*(a/rval)/(tm**3)
                rval = rval - dtp(i)
c-----
c     new 09 FEB 2011
c-----
                dcda(i)=dcda(i)*btp(i)/(tm**3)
                dcda(i+mmax)=dcda(i+mmax)*(a/rval)/(tm**3)
c-----
c    end new
c-----
   20       continue
            c=c/tm
        endif
        return
        end
c
      subroutine svfunc(omega,wvno)
c
c     This combines the Haskell vector from sub down and
c     Dunkin vector from sub up to form the eigenfunctions.
c
      implicit double precision (a-h,o-z)
      integer*4 NL
      parameter (NL=200,NL2=NL+NL)
      real*4 d,a,b,rho,qa1,qb1
      common/model/  d(NL),a(NL),b(NL),rho(NL),qa1(NL),
     *                 qb1(NL),xmu(NL),xlam(NL),mmax,ll
      common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL2),dcdb(NL2),dcdr(NL)
      common/water/ urb,dphw(NL),water0,dphw0,kw
      common/dunk/   uu(NL,5),exe(NL),exa(NL)
      common/hask/   vv(NL,4)
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
      f5=uu(ll,4)
      uu0(1)=wvno*uu(ll,3)/f5
      uu0(2)=1.0d+00
      uu0(3)=fr
      uu0(4)=fr
      ur(ll)=uu(ll,3)/f5
      uz(ll)=1.0d+00
      tz(ll)=-water0
      tr(ll)=0.0d+00
c------
c     uu is Z and vv is X in Wang's dissertation eq. III-1-6 (p.76).
c     i.e., Z a Haskell vector from top layer down and X a Dunkin
c     compound matrix from bottom layer up.
c------
      do 200 i=ll+1,mmax
      i1=i-1
      uu1 =
     *  vv(i1,2)*uu(i,1)+vv(i1,3)*uu(i,2)      +vv(i1,4)*uu(i,3)
      uu2 =
     * -vv(i1,1)*uu(i,1)-vv(i1,3)*uu(i,3)*wvno2+vv(i1,4)*uu(i,4)
      uu3 =
     * -vv(i1,1)*uu(i,2)+vv(i1,2)*uu(i,3)*wvno2+vv(i1,4)*uu(i,5)
      uu4 =
     * -vv(i1,1)*uu(i,3)-vv(i1,2)*uu(i,4)      -vv(i1,3)*uu(i,5)
      ext=0.0d+00
      do 100 k=ll,i1
      ext=ext+exa(k)-exe(k)
100   continue
      fact=0.0d+00
      if(ext.gt.-80.0 .and. ext.lt.80.0 ) fact=dexp(ext)
      ur(i)=uu1*fact/f5
      uz(i)=uu2*fact/f5
      tz(i)=uu3*fact/f5
      tr(i)=uu4*fact/f5
200   continue
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine up(omega,wvno,fr)
c
c     This finds the values of the Dunkin vectors at
c       each layer boundaries from bottom layer upward.
c
      implicit double precision (a-h,o-z)
      integer*4 NL
      parameter (NL=200,NL2=NL+NL)
      dimension ee0(5)
      real*4 d,a,b,rho,qa1,qb1
      common/model/  d(NL),a(NL),b(NL),rho(NL),qa1(NL),
     *                 qb1(NL),xmu(NL),xlam(NL),mmax,ll
      common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL2),dcdb(NL2),dcdr(NL)
      common/dunk/   uu(NL,5),exe(NL),exa(NL)
      common/save/   dd(5,5),aa(4,4),ex1,ex2
      common/aamatx/ ww(NL),xx(NL),yy(NL),zz(NL),
     *                 cospp(NL),cosqq(NL)
      common/water/ urb,dphw(NL),water0,dphw0,kw
        common/engerw/ wra,wd,wba
      wvno2=wvno*wvno
      xka=omega/dble(a(mmax))
      xkb=omega/dble(b(mmax))
      wvnop=wvno+xka
      wvnom=dabs(wvno-xka)
      ra=dsqrt(wvnop*wvnom)
      wvnop=wvno+xkb
      wvnom=dabs(wvno-xkb)
      rb=dsqrt(wvnop*wvnom)
      t = dble(b(mmax))/omega
      gammk = 2.d+00*t*t
      gam = gammk*wvno2
      gamm1 = gam - 1.d+00
      rrho=dble(rho(mmax))
c------
c     since symmetry of compund matrix (component 3 and 4)
c     only 5 terms needed.
c------
      uu(mmax,1)=wvno2-ra*rb
      uu(mmax,2)=-rrho*rb
      uu(mmax,3)=rrho*(gamm1-gammk*ra*rb)
      uu(mmax,4)=rrho*ra
      uu(mmax,5)=rrho*rrho*(gamm1*gamm1-gam*gammk*ra*rb)
c------
c     matrix multiplication from bottom layer upward
c------
      mmx1=mmax-1
      do 400 k=mmx1,ll,-1
      k1=k+1
      dpth=dble(d(k))
      rrho=dble(rho(k))
      xka = omega/dble(a(k))
      xkb = omega/dble(b(k))
      t = dble(b(k))/omega
      gammk = 2.d+00*t*t
      gam = gammk*wvno2
      wvnop=wvno+xka
      wvnom=dabs(wvno-xka)
      ra=dsqrt(wvnop*wvnom)
      p=ra*dpth
      wvnop=wvno+xkb
      wvnom=dabs(wvno-xkb)
      rb=dsqrt(wvnop*wvnom)
      q=rb*dpth
      call var(k,p,q,ra,rb,wvno,xka,xkb,dpth)
      call dnka(wvno2,gam,gammk,rrho)
      exa(k)=ex2
      do 200 i=1,5
      cc=0.0d+00
      do 100 j=1,5
      cc=cc+dd(i,j)*uu(k1,j)
  100 continue
      ee0(i)=cc
  200 continue
c------
c     normalization of ee is important to prevent over- or
c     under-flow.
c------
      call normc(ee0,rab)
c------
c     exe stores the exponential terms and carefully controls
c     the precision problem for P-SV Haskell matrix.
c------
      exe(k)=ex1+rab
      do 300 i=1,5
      uu(k,i)=ee0(i)
  300 continue
  400 continue
c------
c     consider the water layer.
c------
      water0=0.0d+00
      if(ll.eq.1) go to 500
      xka=omega/dble(a(1))
      wvnop=wvno+xka
      wvnom=dabs(wvno-xka)
      ra0=dsqrt(wvnop*wvnom)
      wra=ra0
      p=ra0*dphw0
      znul=1.0d-05
      call var(NL,p,znul,ra0,znul,wvno,xka,znul,dphw0)
c------
c     water0 describes the surface water layer effect.
c------
      rrho=dble(rho(1))
      water0=rrho*ww(NL)/cospp(NL)
      ex3=ex2
      cospp3=cospp(NL)
      if(dabs(cospp3).lt.1.d-30) cospp3=dsign(1.d+00,cospp3)*1.d-30
      ww3=ww(NL)
      tmp=0.0d+00
      if(ex2.gt.-80.0 .and. ex2.lt.80.0) tmp=1.d+00/dexp(ex2)
      uz(1)=tmp/cospp3
        ur(1) = 0.0
        tz(1) = 0.0
c------
c     prepare for subroutine wenerg which takes the energy
c     integral over water layer.
c------
      q=2.d+00*ex3
      tmp=0.0d+00
      if(q.gt.-80.0.and.q.lt.80.0) tmp=1./dexp(q)
      wd=tmp*dphw0/(2.d+00*cospp3*cospp3)
      wba=ww3/(2.d+00*cospp3)
      if(wvno.lt.xka) wra=-wra
      if(kw.le.1) go to 501
      do 450 i=1,kw
        dpth=dphw(i)
        p=ra0*dpth
        call var(NL,p,znul,ra0,znul,wvno,xka,znul,dpth)
        q=ex2-ex3
        tmp=0.0d+00
        if(q.gt.-80.0.and.q.lt.80.0) tmp=dexp(q)
        ur(i)=wvno*ww(NL)/cospp3
        uz(i)=cospp(NL)/cospp3
        tz(i)=-rrho*omega*omega*ww(NL)/cospp3
        ur(i)=tmp*ur(i)
        uz(i)=tmp*uz(i)
        tz(i)=tmp*tz(i)
  450       continue
  501   continue
c-----
c       get radial displacement at base of water
c-----
        dpth=dphw(2)
        p=ra0*dpth
        call var(NL,p,znul,ra0,znul,wvno,xka,znul,dpth)
        q=ex2-ex3
        tmp=0.0d+00
        if(q.gt.-80.0d+00 .and.q.lt.80.0d+00 ) tmp=dexp(q)
        urb=tmp*wvno*ww(NL)/cospp3
c-----
c       define period equation
c-----
  500 continue
      fr=uu(ll,5)+water0*uu(ll,4)
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine down(omega,wvno)
c
c     This finds the values of the Haskell vectors at
c       each layer boundaries from top layer downward.
c
      implicit double precision (a-h,o-z)
      integer*4 NL
      parameter (NL=200,NL2=NL+NL)
      dimension aa0(5)
      real*4 d,a,b,rho,qa1,qb1
      common/model/  d(NL),a(NL),b(NL),rho(NL),qa1(NL),
     *                 qb1(NL),xmu(NL),xlam(NL),mmax,ll
      common/dunk/   uu(NL,5),exe(NL),exa(NL)
      common/hask/   vv(NL,4)
      common/save/   dd(5,5),aa(4,4),ex1,ex2
      common/aamatx/ ww(NL),xx(NL),yy(NL),zz(NL),
     *                 cospp(NL),cosqq(NL)
      wvno2=wvno*wvno
      do 100 j=1,4
      vv(ll,j)=0.0d+00
100   continue
      vv(ll,4)=1.0d+00
      aa0(5)=0.0d+00
        mmx1=mmax-1
      do 500 k=ll,mmx1
      k1=k-1
      if(k.eq.ll) k1=ll
      t=b(k)/omega
      gammk=2.d+00*t*t
      gam=gammk*wvno2
      w=ww(k)
      x=xx(k)
      y=yy(k)
      z=zz(k)
      cosp=cospp(k)
      cosq=cosqq(k)
      rrho=dble(rho(k))
      call hska(w,x,y,z,cosp,cosq,wvno2,gam,gammk,rrho)
      do 300 j=1,4
      cc=0.0d+00
      do 200 i=1,4
      cc=cc+vv(k1,i)*aa(i,j)
200   continue
      aa0(j)=cc
300   continue
      call normc(aa0,ex2)
      exa(k)=exa(k)+ex2
      do 400 i=1,4
      vv(k,i)=aa0(i)
400   continue
500   continue
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine hska(w,x,y,z,cosp,cosq,wvno2,gam,gammk,rho)
c
c     Haskell matrix
c
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
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine dnka(wvno2,gam,gammk,rho)
c
c     Dunkin matrix.
c
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
c     a0 takes care of the factor being exponentially controlled.
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
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
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
c           = pex + qex
c               1.0  ->    exp(-exa)
c-----
      implicit double precision (a-h,o-z)
      integer*4 NL
      parameter (NL=200,NL2=NL+NL)
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
             if(p.lt.16)fac = dexp(-2.0d+00*p)
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
             if(q.lt.16)fac = dexp(-2.0d+0*q)
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
      if(exa.lt.60.0d+00) a0=dexp(-exa)
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
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine normc(ee,ex)
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
      do 10 i = 1,5
      if(dabs(ee(i)).gt.t1) t1 = dabs(ee(i))
   10 continue
      if(t1.lt.1.d-30) t1=1.d+00
      do 20 i =1,5
      ee(i)=ee(i)/t1
   20 continue
      ex=dlog(t1)
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine energy(omega,wvno)
c
c     This takes the energy integral by an analytic
c     way using the eigenfuntions found above.
c     The amplitude factor (sum of energy integrals) thus can be
c     accurately determined without taking any sort of derivatives.
c     The first version of this subroutine is provided by Dr. Harkrider.
c
      implicit double precision (a-h,o-z)
      integer*4 NL
      parameter (NL=200,NL2=NL+NL)
      complex*16 t(4,4),tt(4,4),ff(6),pp(4)
      complex*16 nua,nub,c1,c2
      common/coef/   t,tt,ff,pp
      real*4 d,a,b,rho,qa1,qb1
      common/model/  d(NL),a(NL),b(NL),rho(NL),qa1(NL),
     *                 qb1(NL),xmu(NL),xlam(NL),mmax,ll
      common/eigfun/ uu(NL,4),uu0(4),dcda(NL2),dcdb(NL2),
     *                 dcdr(NL)
      common/sumi/   sumi0,sumi1,sumi2,sumi3,flagr,are,ugr
      common/water/ urb,dphw(NL),water0,dphw0,kw
        omega2=omega*omega
        wvomg2=wvno*omega2
        wvno2=wvno*wvno
        sumi0=0.0d+00
        sumi1=0.0d+00
        sumi2=0.0d+00
        sumi3=0.0d+00
        do 300 k=ll,mmax
            k1=k+1
            kk=k
            if(k.eq.mmax) kk=NL+1
            daa=dble(a(k))
            dbb=dble(b(k))
            drho=dble(rho(k))
            dlam=xlam(k)
            dmu=xmu(k)
            dlamu=dlam+2.d+00*dmu
            dpth=dble(d(k))
            xka=omega/daa
            xkb=omega/dbb
            gammk=dbb/omega
            gammk=2.d+00*gammk*gammk
            gam=gammk*wvno2
            wvnop=wvno+xka
            wvnom=dabs(wvno-xka)
            ra=dsqrt(wvnop*wvnom)
            wvnop=wvno+xkb
            wvnom=dabs(wvno-xkb)
            rb=dsqrt(wvnop*wvnom)
c----
c subterfuge to avoid division by zero
c the matrices should be changed
c-----
            if(ra.lt.1.d-15) ra=1.d-15
            if(rb.lt.1.d-15) rb=1.d-15
            nua=dcmplx(ra,0.0d+00)
            nub=dcmplx(rb,0.0d+00)
            if(wvno.lt.xka) nua=dcmplx(0.0d+00,ra)
            if(wvno.lt.xkb) nub=dcmplx(0.0d+00,rb)
            call tminus(nua,nub,wvno2,gam,gammk,drho)
            call tplus(nua,nub,wvno2,gam,gammk,drho)
            call intval(kk,nua,nub,dpth)
            do 200 i=1,2
                i2=i+2
                c1=0.0d+00
                c2=0.0d+00
                do 100 j=1,4
                    if(kk.ne.(NL+1))c1=c1+tt(i,j)*uu(k1,j)
                    c2=c2+tt(i2,j)*uu(k,j)
100             continue
                pp(i)=c1
                pp(i2)=c2
200         continue
            urur=dsum(kk,1,1)*wvno2
            urtz=dsum(kk,1,3)*wvomg2
            uzuz=dsum(kk,2,2)
            uztr=dsum(kk,2,4)*wvomg2
            tztz=dsum(kk,3,3)*omega2*omega2
            trtr=dsum(kk,4,4)*wvomg2*wvomg2
            urduz=(wvno*dlam*urur+urtz)/dlamu
            uzdur=-wvno*uzuz+uztr/dmu
            durdur=wvno2*uzuz-2.d+00*wvno*uztr/dmu+trtr/(dmu*dmu)
            duzduz=(wvno2*dlam*dlam*urur+2.d+00*wvno*dlam*urtz+
     *          tztz)/(dlamu*dlamu)
            sumi0=sumi0+drho*(uzuz+urur)
            sumi1=sumi1+dlamu*urur+dmu*uzuz
            sumi2=sumi2+dmu*uzdur-dlam*urduz
            sumi3=sumi3+dlamu*duzduz+dmu*durdur
            dldl=-wvno2*urur+2.d+00*wvno*urduz-duzduz
            dldm=2.d+00*(wvno2*urur+wvno*uzdur+duzduz)+
     *          wvno2*uzuz+durdur
            dldm=-dldm
            dldr=omega2*(urur+uzuz)
            dcda(k)=2.d+00*drho*daa*omega*dldl/wvno2
            dcdb(k)=2.d+00*drho*dbb*omega*(dldm-2.*dldl)/wvno2
            dcdr(k)=dldr+dlam*dldl/drho+dmu*dldm/drho
            dcdr(k)=dcdr(k)*omega/wvno2
300     continue
c-----
c       get water layer energy
c-----
        if(ll.ne.1) call wenerg(omega,wvno)
c-----
c       normalize according to total energy
c-----
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
c------
        flagr=omega2*sumi0-wvno2*sumi1-2.d+00*wvno*sumi2-sumi3
        ugr=(wvno*sumi1+sumi2)/(omega*sumi0)
        are=wvno/(2.d+00*omega*ugr*sumi0)
c-----
c       to convert to Levshin and Yanson form we must
c
c       Ur (LJ) = k Ur(Haskell)
c       Uz (LJ) = Uz (Haskell)
c       Tz (LJ) = omega**2 * Tz (Haskell)
c       Tr (LJ) = k omega**2 * Tr (Haskell)
c
c       Except for the water layer, which is OK
c-----
        c = omega/wvno
c       face = 0.5*c**3/(omega*(omega*sumi1 + c*sumi2))
        fac = are*c/wvno2
c-----
c       now get partials with respect to boundary change
c-----
        do 500 k=1,mmax
            if(k.ge.ll)then
                uz = uu(k,2)
                ur = wvno * uu(k,1)
                tz = omega2*uu(k,3)
                tr = wvno*omega2*uu(k,4)
            else
c----- Water is not Haskell, but LJ
                uz = uu(k,2)
                ur = uu(k,1)
                tz = uu(k,3)
                tr = uu(k,4)
            endif
            if(k.eq.1)then
                drho = rho(1) - 0.0
                dmu  = xmu(1) - 0.0
                dlm = xlam(1) - 0.0
                dlmu = dlm + dmu + dmu
                xl2mp = xlam(k) + xmu(k) + xmu(k)
                duzdzp = (tz + wvno*xlam(k)*ur)/xl2mp
                if(xmu(k) .eq.0.0)then
                    durdzp = wvno*uz
                else
                    durdzp = (tr/xmu(k)) - wvno*uz
                endif
                drur2 = ur*ur*drho
                dlur2 = ur*ur*dlmu

                gfac1 = omega2*drho*uz**2
                gfac2 = omega2*drur2
                gfac3 = -wvno2*dmu*uz**2
                gfac4 = -wvno2*dlur2
                gfac5 = (xl2mp*duzdzp**2)
                gfac6 = (xmu(k)*durdzp**2 )
            else
                drho = rho(k) - rho(k-1)
                dmu = xmu(k) - xmu(k-1)
                dlm = xlam(k) - xlam(k-1)
                dlmu = dlm + dmu + dmu
                xl2mp = xlam(k) + xmu(k) + xmu(k)
                xl2mm = xlam(k-1) + xmu(k-1) + xmu(k-1)
                duzdzp = (tz + wvno*xlam(k)*ur)/xl2mp
                if(xmu(k).eq.0.0)then
                    durdzp = wvno*uz
                else
                    durdzp = (tr/xmu(k)) - wvno*uz
                endif
                if(xmu(k-1).eq.0.0)then
                    durdzm = wvno*uz
                else
                    durdzm = (tr/xmu(k-1)) - wvno*uz
                endif
c-----
c       attempt to fix for water layer
c-----
                if(k.eq.ll)then
                    drur2 = ur*ur*rho(k)-urb*urb*rho(k-1)
                    dlur2 = ur*ur*xl2mp -urb*urb*xl2mm
                    duzdzm = (tz + wvno*xlam(k-1)*urb)/xl2mm
                    
                else
                    drur2 = ur*ur*drho
                    dlur2 = ur*ur*dlmu
                    duzdzm = (tz + wvno*xlam(k-1)*ur)/xl2mm
                endif

                gfac1 = omega2*drho*uz**2
                gfac2 = omega2*drur2
                gfac3 = -wvno2*dmu*uz**2
                gfac4 = -wvno2*dlur2
                gfac5 = (xl2mp*duzdzp**2-xl2mm*duzdzm**2)
                gfac6 = (xmu(k)*durdzp**2 - xmu(k-1)*durdzm**2)
            endif
            dfac = fac * (
     1          gfac1 + gfac2 + gfac3 + gfac4
     2          + gfac5 + gfac6 )
            if(dabs(dfac).lt.1.0d-38)then
                dfac = 0.0d+00
            endif
            dcdb(k+mmax) = sngl(dfac)
  500   continue
c-----
c       ultimate kludge
c       As with all of the propagator matrices implemented here, 
c           assume that
c       the fluid layers lie on top of the elastic layers, e.g., the
c       fluid and solids are not interbedded.  
c           Then get the correct dcdh for
c       the fluid/solid interface, we force the total SUM dc/dh = 0
c-----
        if(ll.ne.1)then
        smdcdh = dcdb(mmax+1)
        jbdy = 0
        do 600 k=2,mmax
            if(b(k-1).eq.0.0 .and. b(k).ne.0.0)then
                jbdy = k
            else
                smdcdh = smdcdh + dcdb(mmax+k)
            endif
  600   continue
        if(jbdy.ne.0)then
            dcdb(mmax+jbdy) = -smdcdh
        endif
        endif
c-----
c       now convert from partials with respect to changes in layer
c       boundary to changes in layer thickness
c-----
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
            do 505 i=1,mmax-1
                sm = 0.0
                do 506 j=i+1,mmax
                    sm = sm + dcdb(mmax+j)
  506           continue
                dcdb(mmax+i) = sm
  505       continue
            dcdb(mmax+mmax) = 0.0
        return
        end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      function dsum(kk,i,j)
c
c     The analytic forms of the solution of integral:
c     Integral U*U dz = T-matrix * eigenfnction * integral-coefs
c
      double precision dsum
      integer*4 NL
      parameter (NL=200)
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
      dsum=(sum0)
      if(dabs(dsum).lt.1.0d-30)dsum=0.0d+00
      return
100   continue
      sum1=pp(3)*t(i,3)*t(j,3)*pp(3)*ff(1)
      sum2=pp(4)*t(i,4)*t(j,4)*pp(4)*ff(2)
      sum3=pp(3)*(t(i,3)*t(j,4)+t(i,4)*t(j,3))*pp(4)
      sum3=sum3*ff(3)
      sum0=sum1+sum2+sum3
      dsum=(sum0)
      if(dabs(dsum).lt.1.0d-30)dsum=0.0d+00
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine intval(k,nua,nub,dpth)
c
c     This finds the coeficients needed for integrals.
c     see Wang's dissertation eq. III-1-11 (p.80).
c
      implicit double precision (a-h,o-z)
      integer*4 NL
      parameter (NL=200,NL2=NL+NL)
      complex*16 t(4,4),tt(4,4),ff(6),pp(4)
      complex*16 nua,nub,p,q,pq,expp,exqq,pdum
      common/coef/ t,tt,ff,pp
      if(k.eq.(NL+1)) go to 100
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
      return
100   continue
      ff(1)=0.5d+00/nua
      ff(2)=0.5d+00/nub
      ff(3)=1.0d+00/(nua+nub)
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
      complex*16 p,pq,expq
      p0=dreal(p)
      if(dabs(p0).lt.40.0) go to 100
      expq=0.0d+00
      go to 200
100   continue
      expq=CDEXP(-pq)
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
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine wenerg(omega,wvno)
c
c     calculate energy trapped in the top water layer.
c
      implicit double precision (a-h,o-z)
      integer*4 NL
      parameter (NL=200,NL2=NL+NL)
      real*4 d,a,b,rho,qa1,qb1
      common/model/  d(NL),a(NL),b(NL),rho(NL),qa1(NL),
     *                 qb1(NL),xmu(NL),xlam(NL),mmax,ll
      common/sumi/   sumi0,sumi1,sumi2,sumi3,flagr,are,ugr
      common/eigfun/ uu(NL,4),uu0(4),dcda(NL2),dcdb(NL2),
     *                 dcdr(NL)
      common/water/ urb,dphw(NL),water0,dphw0,kw
        common/engerw/ wra,wd,wba
c
      omega2=omega*omega
      wvno2=wvno*wvno
      wra2=wra*dabs(wra)
      urur  =(wba-wd)/wra2
      urur  =urur*wvno2
      uzuz  =wba+wd
      urduz =wba-wd
      urduz =urduz*wvno
      duzduz=(wba-wd)*wra2
      daa=dble(a(1))
      dr=dble(rho(1))
      dlam=xlam(1)
      sumi0=sumi0+dr*(urur+uzuz)
      sumi1=sumi1+dlam*urur
      sumi2=sumi2-dlam*urduz
      sumi3=sumi3+dlam*duzduz
      dldl=-wvno2*urur+2.d+00*wvno*urduz-duzduz
      dldr=omega2*(urur+uzuz)
      do 100 i=1,kw
      dd=dphw(i+1)-dphw(i)
      rat=dd/dphw0
      dcda(i)=rat*2.d+00*dr*daa*omega*dldl/wvno2
      dcdr(i)=rat*(dldr+dlam*dldl/dr)
      dcdr(i)=dcdr(i)*omega/wvno2
      dcdb(i)=0.0d+00
100   continue
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        subroutine getlove(lun,m)
        implicit none
        integer NL, NL2
        parameter(NL=200,NL2=NL+NL)
        integer lun,m
        integer itst,k1,md1
        real*4 tp1
        integer j
        real gvel,cvel, dudb(NL),dudh(NL),dcdb(NL),dcdh(NL)
        real dcda(NL)
        
        open(3,file='tmpsrfi.05',form='unformatted',access='sequential')
        rewind 3
c-----
c     general purpose routine to read to the end of a sequential
c     file in order to position the file pointer at the end
c     so that the next write statement 'appends' to the file
c     this is necessary since the place at which the file is
c     opened is not necessarily standardized
c-----
      rewind lun
 1000 continue
        read(3,end=9999) itst,k1,md1,tp1        
        write(lun) itst,k1,md1,tp1        
        if(k1.eq.2)then
            read(3,end=9999)gvel,(dudb(j),j=1,m)
            read(3,end=9999) (dudh(j),j=1,m)
            read(3,end=9999) cvel,(dcdb(j),j=1,m)
            read(3,end=9999) (dcdh(j),j=1,m)
            if(itst.eq.4)then
                read(3,end=9999)(dcda(j),j=1,m)
            endif
            write(lun)gvel,(dudb(j),j=1,m)
            write(lun) (dudh(j),j=1,m)
            write(lun) cvel,(dcdb(j),j=1,m)
            write(lun) (dcdh(j),j=1,m)
            if(itst.eq.4)then
                write(lun)(dcda(j),j=1,m)
            endif
        else
            read(3,end=9999)cvel,(dcdb(j),j=1,m)
            read(3,end=9999) (dcdh(j),j=1,m)
            if(itst.eq.2)then
                read(3,end=9999)(dcda(j),j=1,m)
            endif
            write(lun)cvel,(dcdb(j),j=1,m)
            write(lun) (dcdh(j),j=1,m)
            if(itst.eq.2)then
                write(lun)(dcda(j),j=1,m)
            endif
        endif
        go to 1000
 9999   continue
      close (4)
      return
      end
