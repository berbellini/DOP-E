c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME IV                                                       c
c                                                                      c
c      PROGRAM: SRFDRL                                                 c
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
      program srfdrl96
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
c     Changes
c
c     28 DEC 2007 - permit bottom layer to be a fluid
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter(LER=0,LIN=5,LOT=6)
        parameter(NL=200,NL2=NL+NL)
c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c       NL2 - number of columns in model (first NL2/2 are
c           - velocity parameters, second NL2/2 are Q values)
c-----
        common/model/  d(NL),a(NL),b(NL),rho(NL),qa1(NL),
     *                 qb1(NL),xmu(NL),xlam(NL),mmax,ll
        common/eigfun/ ut(NL),tt(NL),dcdb(NL2),dcdr(NL),uu0(4),
     *                 dudb(NL2),btp(NL),dtp(NL)
        common/water/ dphw(NL),water0,dphw0,kw
        common/sumi/   sumi0,sumi1,sumi2,flagr,ale,ugr
        character*10 fname(2)
        data fname/'tmpsrfi.06','tmpsrfi.05'/
c-----
c     machine dependent initialization
c-----
        call mchdep()
c-----
        ip=1
        ig=2
        open(1,file=fname(1),form='unformatted',access='sequential')
        open(2,file=fname(2),form='unformatted',access='sequential')
        rewind 1
        rewind 2
c-----
c       obtain the earth model:
c-----
        read(1) mmax,nsph
        mmax2 = mmax + mmax
        read(1)(btp(i),i=1,mmax)
        read(1)(dtp(i),i=1,mmax)
        do 185 i=1,mmax
            read(1) d(i),a(i),b(i),rho(i)
            xmu(i)=sngl(dble(rho(i))*dble(b(i))*dble(b(i)))
c           xlam(i)=rho(i)*(a(i)*a(i)-2.*b(i)*b(i))
  185   continue
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
c               but make dummy
c               entry
c       mode    surface wave mode 1=FUND 2=1st if 0 mode 
c               does not exist here
c       t   period
c       t1  slightly different period for partial
c       c   phase velocity at t
c       cn  phase velocity at tn
c-----
        t0=t
        if(itst.gt.0) then
            if(igr.gt.0) t0=t*(1.+h)
c-----
c           main part.
c-----
            twopi=2.*3.141592654
            om=twopi/t0
            omega=twopi/t
            wvno=omega/c
            call shfunc(omega,wvno)
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
c           save previous results
c-----
                do 420 i=1,mmax2
                    dudb(i)=dcdb(i)
  420           continue
                call shfunc(omega,wvno)
                call energy(omega,wvno)
                ugr=(ugr+ugp)/2.
                c=(cp+cn)/2.
                do 430 i=1,mmax2
                    dcdn=dcdb(i)
                    dcdb(i)=(dudb(i)+dcdn)/2.
                    uc1=ugr/c
                    ta=uc1*(2.-uc1)*dcdb(i)
                    tb=uc1*uc1*((dudb(i)-dcdn)/(2.*h))
                    dudb(i)=ta+tb
  430           continue
            endif
c-----
c           sphericity correction
c-----
            if(nsph.gt.0)then
                if(igr.eq.0)then
                    call splove(om,c,ugr,mmax2,0)
                else if(igr.eq.1)then
                    call splove(om,c,ugr,mmax2,1)
                else if(igr.eq.2)then
                    call splove(om,c,ugr,mmax2,0)
                    call splove(om,c,ugr,mmax2,1)
                endif
            endif
        else
            mode = 0
            t0 = 0
        endif
c-----
c       output the derivatives.
c       the first mmax elements are partial with respect to velocity
c       the next mmax are the partial with respect to moving the
c       boundary, e.g., takes two to change layer thickness
c-----
        m=mode
        if(igr.eq.0 .or. igr.eq.2) then
            write(2) itst,ip,mode,t0
            if(itst.ne.0) then
                write(2) c,(dcdb(i),i=1,mmax)
                write(2)   (dcdb(i),i=mmax+1,mmax2)
            endif
            if(igr.eq.2)then
                write(4) itst,ig,mode,t0
                if(itst.ne.0) then
                    write(4) ugr,(dudb(i),i=1,mmax)
                    write(4) (dudb(i),i=mmax+1,mmax2)
                    write(4) c,(dcdb(i),i=1,mmax)
                    write(4) (dcdb(i),i=mmax+1,mmax2)
                endif
            endif
        else if(igr.eq.1) then
            write(2) itst,ig,mode,t0
            if(itst.ne.0) then
                write(2) ugr,(dudb(i),i=1,mmax)
                write(2)   (dudb(i),i=mmax+1,mmax2)
                write(2) c,(dcdb(i),i=1,mmax)
                write(2)   (dcdb(i),i=mmax+1,mmax2)
            endif
        endif
        go to 400
  700   continue
c-----
c       end of data read and processing
c       do final clean up, outputting group velocity stuff
c       after all the phase velocity stuff
c-----
        if(igr.ge.2) then
            rewind 4
            do 800 i=1,9000
                read(4,end=900) itst,ig,mode,t0
                write(2) itst,ig,mode,t0
                if(itst.ne.0) then
                    read(4) ugr,(dudb(j),j=1,mmax)
                    read(4)(dudb(j),j=mmax+1,mmax2)
                    read(4)c, (dcdb(j),j=1,mmax)
                    read(4)(dcdb(j),j=mmax+1,mmax2)

                    write(2) ugr,(dudb(j),j=1,mmax)
                    write(2)(dudb(j),j=mmax+1,mmax2)
                    write(2)c, (dcdb(j),j=1,mmax)
                    write(2)(dcdb(j),j=mmax+1,mmax2)

                endif
  800       continue
  900       close(4,status='delete')
        endif
c-----
c       output an indicator that this is the end of the Love Wave data
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

        subroutine splove(om,c,u,mmax2,iflag)
c-----
c       Transform spherical earth to flat earth
c       and relate the corresponding flat earth dispersion to spherical
c
c       Schwab, F. A., and L. Knopoff (1972). Fast surface wave and free
c       mode computations, in  Methods in Computational Physics, 
c               Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c               B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  49, 52  pp 114
c
c       Partial with respect to parameter uses the relation
c       For phase velocity, for example,
c
c       dcs/dps = dcs/dpf * dpf/dps, c is phase velocity, p is
c       parameter, and f is flat model
c
c       om  R*4 angular frequency
c       c   R*4 phase velocity
c       u   R*4 group velocity
c       mmax2   I*4 number of layers* 2, first mmax/2 values are
c               partial with respect to velocity, second are
c               partial with respect to layer thickness
c       iflag   I*4 0 - phase velocity
c               1 - group velocity
c-----

        parameter(NL=200,NL2=NL+NL)
        common/eigfun/ ut(NL),tt(NL),dcdb(NL2),dcdr(NL),uu0(4),
     1      dudb(NL2),btp(NL),dtp(NL)
        a=6370.0
        mmax = mmax2/2
        rval = a
        tm=sqrt(1.+(3.*c/(2.*a*om))**2)
        if(iflag.eq.1) then
            do 10 i=1,mmax
                tmp=dudb(i)*tm+u*c*dcdb(i)*(3./(2.*a*om))**2/tm
                dudb(i)=btp(i)*tmp
                tmp=dudb(i+mmax)*tm+u*c*dcdb(i+mmax)
     1              *(3./(2.*a*om))**2/tm
                dudb(i+mmax) = (a/rval)*tmp
                        rval = rval - dtp(i)
   10       continue
            u=u*tm
        else
            do 20 i=1,mmax
                dcdb(i)=dcdb(i)*btp(i)/(tm**3)
                dcdb(i+mmax)=dcdb(i+mmax)*(a/rval)/(tm**3)
                        rval = rval - dtp(i)
   20       continue
            c=c/tm
        endif
        end

      subroutine shfunc(omega,wvno)
c-----
c     This routine evaluates the eigenfunctions by calling sub
c       up.
c-----
      parameter(NL=200,NL2=NL+NL)
      real*8 exl(NL),ext,fact
      common/model/  d(NL),a(NL),b(NL),rho(NL),qa1(NL),
     *                 qb1(NL),xmu(NL),xlam(NL),mmax,ll
      common/eigfun/ uu(NL,2),dcdb(NL2),dcdr(NL),uu0(4),
     *                 dudb(NL2),btp(NL),dtp(NL)
      common/save/   exl
      call up(omega,wvno,fl)
c-----
c       uu0(2)=stress0 is actually the value of period equation.
c       uu0(3) is used to print out the period euation value before
c       the root is refined.
c-----
        uu0(1)=1.0
        uu0(2)=fl
        uu0(3)=0.0
        uu0(4)=0.0
        ext=0.0
        do 100 k=ll+1,mmax
            ext=ext+exl(k-1)
            fact=0.0
            if(ext.lt.85.0) fact=1./dexp(ext)
            uu(k,1)=uu(k,1)*fact/uu(ll,1)
            uu(k,2)=uu(k,2)*fact/uu(ll,1)
  100   continue
        uu(ll,1)=1.0
        uu(ll,2)=0.0
        return
        end
c
c  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine up(omega,wvno,fl)
c     This routine calculates the elements of Haskell matrix,
c     and finds the eigenfunctions by analytic solution.
c
      parameter(LER=0,LIN=5,LOT=6)
      parameter(NL=200,NL2=NL+NL)
      real*8 exl(NL),qq,rr,ss,exqm,exqp,sinq,cosq
      common/model/  d(NL),a(NL),b(NL),rho(NL),qa1(NL),
     *                 qb1(NL),xmu(NL),xlam(NL),mmax,ll
      common/eigfun/ uu(NL,2),dcdb(NL2),dcdr(NL),uu0(4),
     *                 dudb(NL2),btp(NL),dtp(NL)
      common/save/   exl
      wvno2=wvno*wvno
      xkb=omega/b(mmax)
c-----
c     kludge for fluid core
c-----
      if(b(mmax).gt.0.01)then
             rb=sqrt(abs(wvno2-xkb*xkb))
             if(wvno.lt.xkb)then
                   write(LOT,*) ' imaginary nub derivl'
                   write(LOT,*)'omega,wvno,b(mmax)',omega,wvno,b(mmax)
               endif
             uu(mmax,1)=1.0
             uu(mmax,2)=-xmu(mmax)*rb
      else
             uu(mmax,1)=1.0
             uu(mmax,2)=0.0
      endif
      mmx1=mmax-1
      do 500 k=mmx1,ll,-1
      k1=k+1
      dpth=d(k)
      xkb=omega/b(k)
      rb=abs(wvno2-xkb*xkb)
      rr=dble(rb)
      rr=dsqrt(rr)
      ss=dble(dpth)
      qq=rr*ss
C      if(wvno-xkb) 100,200,300
        if( wvno .lt. xkb)then
            go to 100
        else if(wvno .eq. xkb)then
            go to 200
        else
            go to 300
        endif
100   sinq=dsin(qq)
      cosq=dcos(qq)
      y=sinq/rr
      z=-rr*sinq
      qq=0.0
      go  to 400
200   qq=0.0
      cosq=1.0d+0
      y=dpth
      z=0.0
      go to 400
300   if(qq.gt.40.0) go to 350
      exqp=1.
      exqm=1./dexp(qq+qq)
      sinq=(exqp-exqm)*0.5
      cosq=(exqp+exqm)*0.5
      y=sinq/rr
      z=rr*sinq
      go to 400
350   continue
      y=0.5/rr
      z=0.5*rr
      cosq=0.5
400   continue
      amp0=cosq*uu(k1,1)-y*uu(k1,2)/xmu(k)
      str0=cosq*uu(k1,2)-z*xmu(k)*uu(k1,1)
      rr=abs(amp0)
      ss=abs(str0)
      if(ss.gt.rr) rr=ss
      if(rr.lt.1.d-30) rr=1.d+00
      exl(k)=dlog(rr)+qq
      uu(k,1)=amp0/rr
      uu(k,2)=str0/rr
500   continue
      fl=uu(ll,2)
      return
      end
c
c  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine energy(omega,wvno)
c     This routine calculates the values of integrals I0, I1,
c     and I2 using analytic solutions. It is found
c     that such a formulation is more efficient and practical.
c
      parameter(NL=200,NL2=NL+NL)
      real*8 wvno0,omega0,c,sumi0,sumi1,sumi2
      real*8 xkb,rb,dbb,drho,dpth,dmu,wvno2,omega2
      real*8 upup,dupdup,dcb,dcr
      complex*16 nub,xnub,exqq,top,bot,f1,f2,f3,zarg
      common/model/  d(NL),a(NL),b(NL),rho(NL),qa1(NL),
     *                 qb1(NL),xmu(NL),xlam(NL),mmax,ll
      common/eigfun/ uu(NL,2),dcdb(NL2),dcdr(NL),uu0(4),
     *                 dudb(NL2),btp(NL),dtp(NL)
      common/sumi/   xi0,xi1,xi2,flagr,ale,ugr
        real*8 c2, fac, dvdz, dfac
      wvno0=dble(wvno)
      omega0=dble(omega)
      c=omega0/wvno0
      omega2=omega0*omega0
      wvno2=wvno0*wvno0
      sumi0=0.0d+00
      sumi1=0.0d+00
      sumi2=0.0d+00
      do 300 k=ll,mmax
            k1=k+1
            dbb=dble(b(k))
            drho=dble(rho(k))
            dpth=dble(d(k))
            dmu=dble(xmu(k))
            xkb=omega0/dbb
            rb=dsqrt(dabs(wvno2-xkb*xkb))
            if(k.eq.mmax) then
                  upup  =(0.5/rb)*uu(mmax,1)*uu(mmax,1)
                  dupdup=0.5*rb*uu(mmax,1)*uu(mmax,1)
            else
                  if(wvno0.lt.xkb)then
                        nub=dcmplx(0.0d+00,rb)
                  else
                        nub=dcmplx(rb,0.0d+00)
                  endif
                  xnub=dmu*nub
                  top=uu(k,1)-uu(k,2)/xnub
                  bot=uu(k1,1)+uu(k1,2)/xnub
                  f3=nub*dpth
                  if(dreal(f3).lt.40.0d+00) then
                        zarg = -2.0d+00*f3
                        exqq=dexp(dreal(zarg))*
     1                     dcmplx(dcos(dimag(zarg)),dsin(dimag(zarg)))
c                       exqq=cdexp(-2.*f3)
                  else
                        exqq=dcmplx(0.0d+00,0.0d+00)
                  endif
                  f1=(1.-exqq)/(2.*nub)
                  if(dreal(f3).lt.80.0d+00)then
                        zarg = -f3
                        exqq=dexp(dreal(zarg))*
     1                     dcmplx(dcos(dimag(zarg)),dsin(dimag(zarg)))
c                       exqq=cdexp(-f3)
                  else
                        exqq=dcmplx(0.0d+00,0.0d+00)
                  endif
                  f1=0.25*f1*(top*top+bot*bot)
                  f2=0.5 *dpth*exqq*top*bot
                  upup=dreal(f1+f2)
                  f3=xnub*xnub*(f1-f2)
                  dupdup=dreal(f3)/(dmu*dmu)
            endif
            sumi0=sumi0+drho*upup
            sumi1=sumi1+dmu*upup
            sumi2=sumi2+dmu*dupdup
            dcr=-0.5*c*c*c*upup
            dcb=0.5*c*(upup+dupdup/wvno2)
            dcdb(k)=2.*drho*dbb*dcb
            dcdr(k)=dcr+dbb*dbb*dcb
  300   continue
c-----
c       now that the energy integral I1 is defined get final partial
c-----
        do 400 k=ll,mmax
            dcdb(k)=dcdb(k)/sumi1
            dcdr(k)=dcdr(k)/sumi1
  400   continue
            
c-----
c       get lagrangian, group velocity, energy integral
c-----
            flagr=omega2*sumi0-wvno2*sumi1-sumi2
            ugr=sumi1/(c*sumi0)
            ale=0.5/sumi1
            xi0=sumi0
            xi1=sumi1
            xi2=sumi2
c-----
c       define partial with respect to layer thickness
c-----
c       fac = 0.5d+00*c**3/(omega2*sumi1)
        fac = ale*c/wvno2
        c2 = c*c
c-----
c       for initial layer
c-----
        if(ll.ne.1)then
            dcdb(1) = 0.0
            dcdb(mmax+1) = 0.0
        endif
        do 500 k=ll,mmax
            if(k.eq.ll)then
                drho = rho(k)
                dmu  = xmu(k)
            else
                drho = rho(k) - rho(k-1)
                dmu  = xmu(k) - xmu(k-1)
            endif
            if(k.eq.ll)then
                dvdz = 0.0
            else
                dvdz = uu(k,2)*uu(k,2)*(1.0/xmu(k) - 1.0/xmu(k-1))
            endif
            dfac = fac * ( uu(k,1)*uu(k,1)*
     1          (omega2*drho - wvno2*dmu) + dvdz)
            if(dabs(dfac).lt.1.0d-38)then
                dcdb(k+mmax) = 0.0
            else
                dcdb(k+mmax) = dble(dfac)
            endif
  500   continue
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
                sum = 0.0
                do 506 j=i+1,mmax
                    sum = sum + dcdb(mmax+j)
  506           continue
                dcdb(mmax+i) = sum
  505       continue
            dcdb(mmax+mmax) = 0.0
        return
        end
