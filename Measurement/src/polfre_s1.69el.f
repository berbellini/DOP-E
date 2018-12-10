      program polfre_s1
c
c-----------------------------------------------------------
c This version: polfre_s1.69el:
c =============================
c This is the same as polfre_s1.68el.f but output zdeg and wdeg
c optionally.
c Below the changes of polfre_s1.68el:
c This is the same as polfre_s1.66el.f with some addition:
c (1) The output file contains now 6 more columns for the normalized
c semi-major and semi-minor vector length.
c (2) Set DOP to zero whenever the particle motion is in a plane
c which  deviates from the vertical plane by more than "zdeg" deg.
c (3) Set DOP to zero whenever the semi major vector is not
c in the horizontal plane or parallel to the vertical axis. 
c A tolerance/deviation of "wdeg" deg is permitted.
c "zdeg" and "wdeg" are the variables used in the code.
c (4) Option to output 6 more columns for the normalized
c semi-major and semi-minor vector.
c (5) Option to output H/V ratio.
c-----------------------------------------------------------
c
c WORKING PG for linear/elliptical particle motion.
c
c This program computes the degree of polarization (dop)
c for three component sensors employing sliding data windows
c as function of frequency.
c
c The employed dop strategy has been introduced in 
c Schimmel & Gallart (2003) using analytic signal theory in 
c the t-domain. Later, in Schimmel & Gallart (2004), the dop 
c strategy has been extended to the tf-domain, but 
c using now a different strategy to compute the polarization 
c attributes. Main reason for changing the way to determine 
c polarization attributes has been to get experiences with 
c other approaches and to show that the dop strategy is 
c independent on the polarization approach.
c
c Thus, this program uses an eigen approach (eigen values and 
c eigen vectors) to construct the dop following the
c time-domain definition. A freq-dependent Gaussian window 
c can be used. This version is up-graded using the S-Transform 
c which is very similar to my 1st approach. 
c
c The polarization approach has been used for filtering 
c Schimmel and Gallart (2003,2004,2005) and signal characterization
c Stutzmann et al. (2009), Schimmel et al. (2011), Obrebski et al. (2012),
c Sergeant et al. (2013), among others.
c
c REFERENCES for methodology and examples:
c Schimmel M., and J. Gallart, The use of instantaneous 
c    polarization attributes for seismic signal detection 
c    and image enhancement , Geophys.J.Int.,, 155, 653-668, 
c    doi:10.1046/j.1365-246X.2003.02077.x, 2003.
c Schimmel & Gallart, Degree of polarization filter for
c    frequency-dependent signal enhancement through noise
c    suppression, Bull.Seism.Soc.Am., 94, 1016-1035, 
c    doi: 10.1785/0120030178, 2004.
c Schimmel & Gallart, The inverse S Transform in filters 
c    with time-frequency localization , IEEE Transactions on 
c    Signal Processing, 53 (11), 4417 - 4422, 
c    doi:10.1109/TSP.2005.857065, 2005.
c Stutzmann, E., Schimmel, M., Patau, G., Maggi, A., Global 
c    climate imprint on seismic noise , Geochem. Geophys. Geosyst., 
c    10, Q11004, doi:10.1029/2009GC002619, 2009
c Schimmel, M., Stutzmann, E., Ardhuin, F., Gallart, J., 
c    Polarized Earth's Ambient Microseismic Noise , Geochem. Geophys. 
c    Geosyst., doi:10.1029/2011GC003661, 2011.
c Obrebski, M.J., Ardhuin, F., Stutzmann, E., Schimmel, M., 
c    How moderate sea states can generate loud seismic noise in 
c    the deep ocean, Geophys. Res. Lett., 39, L11601, 
c    doi: 10.1029/2012GL051896, 2012.
c Sergeant A., Stutzmann E., Maggi A., Schimmel M., Ardhuin F., 
c    Obrebski M., Frequency-dependent noise sources in the North 
c    Atlantic Ocean, Geochem. Geophys. Geosyst., 14, 
c    doi:10.1002/2013GC004905, 2013.
c
c (PDFs can be downloaded from my home-page.)
c---------------------------------------------------------------
c This a WORKING (RESEARCH) VERSION rather than a final program.
c I constantly adapt it. ...
c Simplifications are required for the straight final code.
c The program contains options which are not useful any more
c or which were coded for a very specific purpose. These will be
c removed in a future version of this program.
c This research program is not for distribution.
c Copyright & Author: Martin Schimmel (schimmel@ictja.csic.es)
c
c HOW to start:
c 1. compile together with sac-library.
c 2. just execute polfre_s1.66el without any argument to
c    obtain the list of parameters and their usage.
c 3. see my example and play with the scripts.
c
c BUGS & COMMENTS: Please, report remaining bugs! Also, do
c not hesitate to send me your comments on interesting,
c modifications or results.
c----------------------------------------------------------------
      include 'include_pol.h'
c use: nftf=nstf/2 + 1, nstf=2**n, maxns<nstf
c     parameter (nstf=1024,nftf=513 )
      parameter (nstderr=0,ioutdat=10)

c deg of polarization 
      real degpol(maxtr,maxns,maxfr)
      real sum_maj(maxns),sum_pla(maxns)
      real mocyc
c state vector
      real r1(maxns,maxfr,3),r2(maxns,maxfr,3),r3(maxns,maxfr,3)
      double precision  r1state(3),r2state(3)
      real r1amp(maxns,maxfr),r2amp(maxns,maxfr)
c spectrum
      complex cspec, cmean, c1spec, c2spec, c3spec
      complex ctf1(nstf,nftf),ctf2(nstf,nftf),ctf3(nstf,nftf)
c gaussian window
      real rycle,nocyc
c spectral matrix:
      complex c11,c12,c13,c21,c22,c23,c31,c32,c33
c complex variables
      complex c1,c2,c3
c eigen values & vector
      double precision dsr(3,3),dsi(3,3)
      double precision devectr(3,3),devecti(3,3),devalues(3)
c dummies for eigen routines:
      double precision dfv1(3),dfv2(3),fm1(2,3)
c semi major (real part of state vector):
      real rsig1(maxns),rsig2(maxns),rsig3(maxns)
c planarity:
      real plan1(maxns),plan2(maxns),plan3(maxns)
c linearity:
      real linear(maxns,maxfr),lin(maxns)
c dummy:
      real rdum1(maxns)
c for azimuth output:
      real azim(maxns,maxfr)
      real sig1(maxns),sig2(maxns),sig3(maxns)
c file names:
      character*80 name1, name2, name3
c info file:
       character par*14 
       character kstnm*8
c      will be changed soon
       character ahour*2
c logical switches:
      logical lasc,lsac,lbin,lmed,lbug,lave
      logical lfre,lwlenf
      logical llin,llim,lell,llh,lelm
c i/o parameters:
      real pow,f1,f2,dt
      integer ns,nt,nf,nsamp,ntrac,nfreq
      integer wlen,wlenf
c input data:
      real trace1(maxns,maxtr),trace2(maxns,maxtr)
      real trace3(maxns,maxtr)
c function
      complex dftc
c Version 1.68:
      real zdeg,wdeg,zrad,wrad
      logical lhv,lsemi,lzwout

      twopi=6.28318530717958
      pi=twopi/2.
      rad2deg=57.29577951
      lbug=.false.

      narg=iargc()
      if (narg.lt.3) then
        write(*,*)
        write(*,*)'Dimensions:'
        write(*,*)'MAX: # samples, traces, frequencies:'
        write(*,*)maxns,maxtr,maxfr
        write(*,*)'--------------------------------------------'
        write(*,*)
        call usepol
      endif
cccccccccccccc
c Zero arrays:
cccccccccccccc
      write(*,*)'  ... warming up ... '
      write(*,*)'  '
      do ns=1,maxns
         sum_maj(ns)=0.
         sum_pla(ns)=0.
         rdum1(ns)=0.
         rsig1(ns)=0.
         rsig2(ns)=0.
         rsig3(ns)=0.
         plan1(ns)=0.
         plan2(ns)=0.
         plan3(ns)=0.
         lin(ns)=0.
         do nf=1,maxfr
           r1amp(ns,nf)=0.
           r2amp(ns,nf)=0.
           linear(ns,nf)=0.
           do i=1,3
             r1(ns,nf,i)=0.
             r2(ns,nf,i)=0.
             r3(ns,nf,i)=0.
           enddo
         enddo
      enddo
      do i=1,3
        dfv1(i)=0.D0
        dfv2(i)=0.D0
        fm1(1,i)=0.D0
        fm1(2,i)=0.D0
        devalues(i)=0.D0
        do j=1,3
          dsr(i,j)=0.D0
          dsi(i,j)=0.D0
          devectr(i,j)=0.D0
          devecti(i,j)=0.D0
        enddo
      enddo
      do nt=1,maxtr
         do ns=1,maxns
           trace1(ns,nt)=0.
           trace2(ns,nt)=0.
           trace3(ns,nt)=0.
           do nf=1,maxfr
              degpol(nt,ns,nf)=0.
           enddo
         enddo
      enddo

ccccccccccccccccccccccccccccc
c Handle the parameter input:
ccccccccccccccccccccccccccccc
      call getpar(name1,name2,name3,f1,f2,dt,nfreq,ntrac,nsamp,wlen,
     &   wlenf,rycle,nocyc,mocyc,pow,lasc,lbin,lsac,lmed,lave,nflen,
     &   lfre,fre,lwlenf,llin,llim,lelm,lell,dopm,llh,rlll,zdeg,wdeg,
     &   lhv,lsemi,lzwout)
c
cccccccccccc
c READ DATA:
cccccccccccc
      if (lsac) then
        call sac_read(name1,trace1,1,nsampa,dt1,bega)
        call getnhv('nzhour',nhour1,nerr)
        call getnhv('nzmin',nmin1,nerr)
        call getnhv('nzsec',nsec1,nerr)
        call getnhv('nzmsec',nmsec1,nerr)
        call sac_read(name2,trace2,1,nsampb,dt2,begb)
        call getnhv('nzhour',nhour2,nerr)
        call getnhv('nzmin',nmin2,nerr)
        call getnhv('nzsec',nsec2,nerr)
        call getnhv('nzmsec',nmsec2,nerr)
        call sac_read(name3,trace3,1,nsampc,dt3,begc)
        call getnhv('nzhour',nhour3,nerr)
        call getnhv('nzmin',nmin3,nerr)
        call getnhv('nzsec',nsec3,nerr)
        call getnhv('nzmsec',nmsec3,nerr)

c check start time:
        if (nhour1.ne.nhour2.and.nhour1.ne.nhour3) then
           stop 'different start time (hour)'
        else if (nmin1.ne.nmin2.and.nmin1.ne.nmin3) then
           stop 'different start time (min)'
        else if (nsec1.ne.nsec2.and.nsec1.ne.nsec3) then
           stop 'different start time (sec)'
        else if (dt1.lt.1) then
           if (nmsec1.ne.nmsec2.and.nmsec1.ne.nmsec3) then
             stop 'different start time (msec)'
           endif
        endif

        if (dt.le.0.) dt=dt1
        if (abs(dt-dt1).gt.0.0001) then
          write(*,*)'dt,dt1:',dt,dt1
          stop' sac files with different sample intervals?'
        endif
        if (abs(dt1-dt2).gt.0.0001) then
          write(*,*)'dt,dt1:',dt1,dt2
          stop' sac files with different sample intervals?'
        endif
        if (abs(dt1-dt3).gt.0.0001) then
          write(*,*)'dt1,dt3:',dt1,dt3
          stop' sac files with different sample intervals?'
        endif
        if (bega.ne.begb.or.bega.ne.begc) then
          write(*,*)' Data with different beg!',bega,begb,begc
          if ((abs(bega-begb).ge.(dt1/10.)).or.
     &          (abs(bega-begc).ge.(dt1/10))) then
          write(*,*)'beg difference gt dt1 div 10!',bega,begb,begc
          stop
          endif
        endif
        beg=bega
        if (nsamp.gt.0) then
         nsamp=min(nsamp,nsampa,nsampb,nsampc)
        else
         nsamp=min(nsampa,nsampb,nsampc)
        endif
        call getkhv('kstnm',kstnm,nerr)
        ll=leng(kstnm)
        call getnhv('nzyear',iyear,nerr)
        call getnhv('nzjday',jjj,nerr)
        if (jjj.lt.10) then
         write(par,107)jjj
 107     format("00",i1)
        else if (jjj.lt.100) then
         write(par,108)jjj
 108     format("0",i2)
        else
         write(par,109)jjj
 109     format(i3)
        endif

c include hour (will be changed soon)
        if (nhour1.lt.10)  then
         write(ahour,'(a1,i1)') '0',nhour1
        else
         write(ahour,'(i2)') nhour1
        endif
        par(1:ll+7)=par(1:3)//"."//ahour//"."//kstnm(1:ll)

c       par(1:ll+4)=par(1:3)//"."//kstnm(1:ll)
ccc     open(21,file='info_data.tmp')
ccc     rewind(21)
        if (nsamp.lt.50) then
          nsamp=1
ccc       write(21,*)par,nsamp,dt
          stop 'hard coded stop: nsamp<50'
        endif
c       write(21,110)iyear,'.',par(1:ll+4),nsamp,dt
c include hour (will be done different)
ccc     write(21,110)iyear,'.',par(1:ll+7),nsamp,dt
c110    format(i4,a1,A,1x,i5,1x,f10.7)

ccc     close(21)
      endif
      if (lasc) stop ' Option not working yet!'
      if (lbin) then
c Handle bin files:
c------------------
        call bin_read(name1,trace1,ntrac,0,1,nsamp)
        call bin_read(name2,trace2,ntrac,0,1,nsamp)
        call bin_read(name3,trace3,ntrac,0,1,nsamp)
      endif
c bin dop output:
      if (lsac) lbin=.true.
c
ccccccccccccccccc
c Set parameters: 
ccccccccccccccccc
c
      if (rycle.gt.0.01) then
c get power for fft,freqs, ... :
c-------------------------------
        npow=1
        nsmp=2
 5      if (nsamp.le.nsmp) goto 10
        npow=npow+1
        nsmp=nsmp*2
        goto 5
 10     nsmp2=nsmp/2
        if (nsmp2.gt.nftf) stop 'change dims of ctf'
        if (nsmp.gt.nstf) stop 'change dims of ctf'
 
        fnyq=0.5/dt
        if (f2.gt.fnyq) f2=fnyq

        df=1./(float(nsmp)*dt)
        f1=nint(f1/df)*df
        f2=nint(f2/df)*df
        dfstep=(f2-f1)/float(nfreq-1)
        if (dfstep.lt.df) dfstep=df
        ndfstep=nint(dfstep/df)
        dfstep=ndfstep*df
        nfreq=int((f2-f1)/dfstep) + 2
        f2=(nfreq-1)*dfstep+f1
        if (nfreq.gt.maxfr) then
          nfreq=maxfr
          f2=(nfreq-1)*dfstep+f1
        endif
        nflen2=2*nflen+1
      else
        nsmp=nsamp
        nsmp2=nsamp/2
        df=1./(float(nsamp)*dt)
        f1=nint(f1/df)*df
        f2=nint(f2/df)*df
        dfstep=(f2-f1)/float(nfreq-1)
        if (dfstep.lt.df) dfstep=df
        ndfstep=nint(dfstep/df)
        dfstep=ndfstep*df
        nfreq=int((f2-f1)/dfstep) + 2
        f2=(nfreq-1)*dfstep+f1
        if (nfreq.gt.maxfr) then
          nfreq=maxfr
          f2=(nfreq-1)*dfstep+f1
        endif
        nflen2=2*nflen+1
        truncat=0.1
      endif

      if (f2.gt.fnyq) stop ' f2 > fnyq '

c Others:
      zrad=zdeg/rad2deg
      wrad=wdeg/rad2deg
      wrad2=twopi/4.-wrad

c Parameter check:
      write(*,*)'  '
      write(*,*)'DATA PROPERTIES: '
      write(*,*)'  # samples, 3-C traces            : ',nsamp,ntrac
      write(*,*)'  time increment dt                : ',dt
      write(*,*)'  natural frequency increment df   : ',df
      write(*,*)'PROCESSING:'
      if (lwlenf) then
      write(*,*)'  use freq-dependent dop window!'
      write(*,*)'  window at f2, power              : ',wlenf,pow
      else
      write(*,*)'  window,power                     : ',wlen,pow
      endif
      write(*,*)'  frequencies f1,f2                : ',f1,f2
      write(*,*)'  # frequencies                    : ',nfreq
      write(*,*)'  # frequencies for spec average   : ',nflen2
      write(*,*)'  Set dop=0 for angle(motion plane, Z) > zdeg'
      write(*,*)'  and/or angle(semi major, H) > wdeg or'
      write(*,*)'  angle(semi major, Z) > wdeg '
      write(*,*)'  zdeg, wdeg                       : ',zdeg,wdeg 
      if (lave) then
      write(*,*)'  (average spectral matrices)'
      else
      write(*,*)'  (average frequency spectrum)'  
      endif
      write(*,*)'  increment in frequency range     : ',dfstep
      write(*,*)'  ( ',ndfstep,' * ',df,' = ',dfstep,' )'
      write(*,*)'  median                           : ',lmed 
      if (rycle.gt.0.01) then
      write(*,*)'  cycles within 2sigma window      : ',rycle
      write(*,*)'  # samples in 2sigma window for f1: ',rycle/f1/dt+1
      write(*,*)'  # samples in 2sigma window for f2: ',rycle/f2/dt+1
      else
      write(*,*)'  2*std of Gaussian window (time)  : ',nocyc
      write(*,*)'  # samples in 2sigma window       : ',nint(nocyc/dt)+1
      write(*,*)'  truncate Gaussian window at      : ',truncat 
cc    write(*,*)'  (Grrrr, still old fashion!)'
      endif
      if (lfre) write(*,*)'  write fre.asc with dop at fre: ',fre
      write(*,*)'  '

cccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (dopm.lt.0.) then
       open(131,file="azi_spec.asc")
       rewind(131)
       write(131,202)beg,(nsamp-1)*dt+beg,dt,f1,f1+(nfreq-1)*dfstep,
     & dfstep
       open(136,file="azi_plan.asc")
       rewind(136)
       write(136,202)beg,(nsamp-1)*dt+beg,dt,f1,f1+(nfreq-1)*dfstep,
     & dfstep
       open(137,file="azi_majo.asc")
       rewind(137)
       write(137,202)beg,(nsamp-1)*dt+beg,dt,f1,f1+(nfreq-1)*dfstep,
     & dfstep
      endif

 202  format(3(F10.2,1x),3(F14.10,1x))


c====================================
c   S T A R T  P R O C E S S I N G :
c====================================
ccc   n10=ntrac*0.1
ccc   n20=ntrac*0.2
ccc   n40=ntrac*0.4
ccc   n60=ntrac*0.6
ccc   n80=ntrac*0.8
ccc   n90=ntrac*0.9
ccc   if (n10.eq.0) n10=1
      write(*,*)'      .... start working hard!'

ccccccccccccccccccc
c Loop over traces:
ccccccccccccccccccc
      do nt=1,ntrac
c Provide user with progress in %:
ccc     if (nt.eq.n10) write(*,*)'  10% done '
ccc     if (nt.eq.n20) write(*,*)'  20% done '
ccc     if (nt.eq.n40) write(*,*)'  40% done '
ccc     if (nt.eq.n60) write(*,*)'  60% done '
ccc     if (nt.eq.n80) write(*,*)'  80% done '
ccc     if (nt.eq.n90) write(*,*)'  90% done '

        do ns=1,nsamp
          sig1(ns)=trace1(ns,nt)
          sig2(ns)=trace2(ns,nt)
          sig3(ns)=trace3(ns,nt)
        enddo

c Get localized spectra:
c------------------------
        if (rycle.gt.0.01) then
           call st_4ward(sig1,nsamp,rycle,dt,ctf1,nsmp,npow)
           call st_4ward(sig2,nsamp,rycle,dt,ctf2,nsmp,npow)
           call st_4ward(sig3,nsamp,rycle,dt,ctf3,nsmp,npow)
        endif

        stdwidth=nocyc/2.
ccccccccccccccccc
c Loop over time:
ccccccccccccccccc
        do ns=1,nsamp

ccccccccccccccccccccccccccccccc
c Loop over center frequencies:
ccccccccccccccccccccccccccccccc
          bbeg=beg+(ns-1)*dt
          do nf=1,nfreq
           fcent=f1+(nf-1)*dfstep 
           ifcent=nint(fcent/df)+1
           iff=ifcent+1+nflen
           if (nt.eq.1.and.ns.eq.1) then
             write(*,*)' Center frequencies:',fcent
           endif

cccccccccccccccccccccccccccccccccccccccc
c determine azimuth of major axis
cccccccccccccccccccccccccccccccccccccccc
       if (dopm.lt.0.) then
           tttt=twopi*fcent*bbeg
       zzz3=real(ctf3(ns,ifcent)*cexp(cmplx(0.,-tttt)))
       zzz2=real(ctf2(ns,ifcent)*cexp(cmplx(0.,-tttt)))
           azi=rad2deg*atan2(zzz3,zzz2)
           if (azi.gt.180.) azi=azi-180.
           if (azi.lt.0.) azi=azi+180.
           write(131,*)bbeg,fcent,azi
       endif

cccccccccccccccccccccccccc
c Average spectral matrix
cccccccccccccccccccccccccc
           if (lave) then
             c11=0.
             c22=0.
             c33=0.
             c12=0.
             c13=0.
             c23=0.
             mflen2=0
             do nnf=1,nflen2
               nff=iff-nnf
               if (nff.lt.1.or.nff.gt.nsmp2) goto 11
               mflen2=mflen2+1
               if (rycle.gt.0.01) then
                 c1spec=ctf1(ns,nff)
                 c2spec=ctf2(ns,nff)
                 c3spec=ctf3(ns,nff)
               else
                 freq=fcent-(nnf-1-nflen)*df
cccc             stdwidth=nocyc/2.
cccc             nstd=nint(stdwidth/dt)
                 do icomp=1,3
                   if (icomp.eq.1) call window2(sig1,rdum1,nsamp,
     &                nlenw,stdwidth,dt,ns,truncat)
                   if (icomp.eq.2) call window2(sig2,rdum1,nsamp,
     &                nlenw,stdwidth,dt,ns,truncat)
                   if (icomp.eq.3) call window2(sig3,rdum1,nsamp,
     &                nlenw,stdwidth,dt,ns,truncat)

                   rnorm=nlenw
                   if (nlenw.lt.mocyc) rnorm=mocyc
                    
                    cspec=dftc(rdum1,nlenw,freq,dt)
                    if (icomp.eq.1) c1spec=cspec/rnorm
                    if (icomp.eq.2) c2spec=cspec/rnorm
                    if (icomp.eq.3) c3spec=cspec/rnorm
                 enddo
               endif

c spectral matrix:
c-----------------
               c11=c11+c1spec*conjg(c1spec)
               c22=c22+c2spec*conjg(c2spec)
               c33=c33+c3spec*conjg(c3spec)
               c12=c12+c1spec*conjg(c2spec)
               c13=c13+c1spec*conjg(c3spec)
               c23=c23+c2spec*conjg(c3spec)
 11            continue
             enddo
c
c The averaged spectrum:
c-----------------------
             c11=c11/float(mflen2)
             c22=c22/float(mflen2)
             c33=c33/float(mflen2)
             c12=c12/float(mflen2)
             c13=c13/float(mflen2)
             c23=c23/float(mflen2)
             c21=conjg(c12)
             c31=conjg(c13)
             c32=conjg(c23)
           else
cccccccccccccccccccccccccc
c Average over frequency :
cccccccccccccccccccccccccc
             if (rycle.gt.0.01) then
               c1spec=cmplx(0.,0.)
               c2spec=cmplx(0.,0.)
               c3spec=cmplx(0.,0.)
               mflen2=0
               do nnf=1,nflen2
                 nff=iff-nnf
                 if (nff.ge.1.and.nff.le.nsmp2) then
                   mflen2=mflen2+1
                   c1spec=c1spec+ctf1(ns,nff)
                   c2spec=c2spec+ctf2(ns,nff)
                   c3spec=c3spec+ctf3(ns,nff)
                 endif 
               enddo
               if (mflen2.ne.0) then
                 c1spec=c1spec/float(mflen2)
                 c2spec=c2spec/float(mflen2)
                 c3spec=c3spec/float(mflen2)
               endif
             else
               do icomp=1,3
                 cmean=cmplx(0.,0.)
                 do nnf=1,nflen2
c still to check for frequencies aliasing!!!
                   freq=fcent-(nnf-1-nflen)*df
cccc               stdwidth=nocyc/2.
cccc               nstd=nint(stdwidth/dt)
                   if (icomp.eq.1) call window2(sig1,rdum1,nsamp,
     &                nlenw,stdwidth,dt,ns,truncat)
                   if (icomp.eq.2) call window2(sig2,rdum1,nsamp,
     &                nlenw,stdwidth,dt,ns,truncat)
                   if (icomp.eq.3) call window2(sig3,rdum1,nsamp,
     &                nlenw,stdwidth,dt,ns,truncat)

                   rnorm=nlenw
                   if (nlenw.lt.mocyc) rnorm=mocyc
                    
                   cspec=dftc(rdum1,nlenw,freq,dt)
                   cmean=cmean+cspec/rnorm
                 enddo
                 if (icomp.eq.1) c1spec=cmean/float(nflen2)
                 if (icomp.eq.2) c2spec=cmean/float(nflen2)
                 if (icomp.eq.3) c3spec=cmean/float(nflen2)
               enddo
             endif
  
c spectral matrix:
c-----------------
             c11=c1spec*conjg(c1spec)
             c22=c2spec*conjg(c2spec)
             c33=c3spec*conjg(c3spec)
             c12=c1spec*conjg(c2spec)
             c13=c1spec*conjg(c3spec)
             c23=c2spec*conjg(c3spec)
             c21=conjg(c12)
             c31=conjg(c13)
             c32=conjg(c23)
           endif

c  determine eigen values and vectors explicitely:
               dsr(1,1)=dble( real(c11) )
               dsi(1,1)=0.0D0
               dsr(1,2)=dble( real(c12) )
               dsi(1,2)=dble( aimag(c12) )
               dsr(1,3)=dble( real(c13) )
               dsi(1,3)=dble( aimag(c13) )
               dsr(2,1)=dble( real(c21) )
               dsi(2,1)=dble( aimag(c21) )
               dsr(2,2)=dble( real(c22) )
               dsi(2,2)=0.0D0
               dsr(2,3)=dble( real(c23) )
               dsi(2,3)=dble( aimag(c23) )
               dsr(3,1)=dble( real(c31) )
               dsi(3,1)=dble( aimag(c31) )
               dsr(3,2)=dble( real(c32) )
               dsi(3,2)=dble( aimag(c32) )
               dsr(3,3)=dble( real(c33) )
               dsi(3,3)=0.0D0
    
               matz=1
               nm=3
               call ch(nm,nm,dsr,dsi,devalues,matz,devectr,
     &                 devecti,dfv1,dfv2,fm1,ierr)
               acr=sngl(devectr(1,3))
               aci=sngl(devecti(1,3))
               c1=cmplx(acr,aci)
               acr=sngl(devectr(2,3))
               aci=sngl(devecti(2,3))
               c2=cmplx(acr,aci)
               acr=sngl(devectr(3,3))
               aci=sngl(devecti(3,3))
               c3=cmplx(acr,aci)

               call state_rot(c1,c2,c3,r1state,r2state,r1ampl,
     &             r2ampl)
  
c State vector r1+ir2 and its amplitudes:
            r1amp(ns,nf)=r1ampl
            r2amp(ns,nf)=r2ampl
            r1(ns,nf,1)=sngl(r1state(1))
            r1(ns,nf,2)=sngl(r1state(2))
            r1(ns,nf,3)=sngl(r1state(3))
            r3(ns,nf,1)=sngl(r2state(1))
            r3(ns,nf,2)=sngl(r2state(2))
            r3(ns,nf,3)=sngl(r2state(3))
            call crosspn(r1state,r2state,dfv1)
c new state vector: ==> r1 + i r2 = r1 + i pla
            r2(ns,nf,1)=sngl(dfv1(1))
            r2(ns,nf,2)=sngl(dfv1(2))
            r2(ns,nf,3)=sngl(dfv1(3))
            linear(ns,nf)=1.-r2ampl/r1ampl

            continue
          enddo

        enddo

        if (lbug) write(*,*)' Get deg of pol:'
        awlenf=float(wlenf)*f2
        do nf=1,nfreq
          fcent=f1+(nf-1)*dfstep
          if (lwlenf) wlen=nint(awlenf/fcent)
          do ns=1,nsamp
             rsig1(ns)=r1(ns,nf,1)
             rsig2(ns)=r1(ns,nf,2)
             rsig3(ns)=r1(ns,nf,3)
             plan1(ns)=r2(ns,nf,1)
             plan2(ns)=r2(ns,nf,2)
             plan3(ns)=r2(ns,nf,3)
             lin(ns)=linear(ns,nf)

             if (dopm.lt.0.) then
c azi following planarity vector:
              azi=rad2deg*atan2(plan3(ns),plan2(ns)) - 90.0
              if (azi.gt.180.) azi=azi-360.
              if (azi.lt.-180.) azi=azi+360.
c 180 deg ambiguity:
              begdt=beg+(ns-1)*dt
              write(136,*)begdt,fcent,azi
c azi following major vector:
              azi=rad2deg*atan2(rsig3(ns),rsig2(ns))
              if (azi.gt.180.) azi=azi-180.
              if (azi.lt.0.) azi=azi+180.
              write(137,*)begdt,fcent,azi
             endif

             if (lell) then
               azi=rad2deg*atan2(plan3(ns),plan2(ns))-90.
               if (azi.gt.180.) azi=azi-360.
               if (azi.lt.-180.) azi=azi+360.
               azim(ns,nf)=azi
             elseif (llh) then
               azi=rad2deg*atan2(rsig3(ns),rsig2(ns))-90.
               if (azi.gt.180.) azi=azi-180.
               if (azi.lt.0.) azi=azi+180.
               azim(ns,nf)=azi
             elseif (dopm.gt.0) then
c Still provide an azimuth although this does make
c sence only for the elliptically polarized signals.
               azi=rad2deg*atan2(plan3(ns),plan2(ns))-90.
               if (azi.gt.180.) azi=azi-360.
               if (azi.lt.-180.) azi=azi+360.
               azim(ns,nf)=azi
             endif
c check what happens to the others!

          enddo
           
          if (lbug) write(*,*)' smooth linearity:'
          if (lmed) then
            call movmedian(wlen,lin,rdum1,nsamp)
          else
            call movmean(wlen,lin,rdum1,nsamp)
          endif
          if (lbug) write(*,*)' before vector projection'
cccc      if (lmed) then
cccc       call mov_proj_med(rsig1,rsig2,rsig3,sum_maj,nsamp,
cccc &         wlen,pow)
cccc       call mov_proj_med(plan1,plan2,plan3,sum_pla,nsamp,
cccc &         wlen,pow)
cccc      else
           call mov_proj_mean(rsig1,rsig2,rsig3,sum_maj,
     &         nsamp,wlen,pow)
           call mov_proj_mean(plan1,plan2,plan3,sum_pla,nsamp,
     &         wlen,pow)
cccc      endif
          if (lbug) write(*,*)' after vector projection'

          do ns=1,nsamp
            rdu=rdum1(ns)
            if (llin) then
               rrr=rdu
            else if (llim) then
               rrr=sum_maj(ns)
            else if (lelm) then
               rrr=sum_pla(ns)
            else if (lell) then
c            if (rdu.lt.rlll) then
             if (lin(ns).lt.rlll) then
c weight by deviation of vertical plane
c vertical component of the p (p is normalized!):
               aaaa=plan1(ns)
               aaa1=asin(aaaa)
               aaaa=cos(aaa1)
               rrr=sum_pla(ns)*aaaa**8
             else
               rrr=0.
             endif
            else if (llh) then
c            if (rdu.gt.rlll) then
             if (lin(ns).gt.rlll) then
c vertical component of the a (a is normalized!):
               aaaa=rsig1(ns)
               aaa2=asin(aaaa)
               aaaa=cos(aaa2)
c rdu is not squared since aaaa more important !
               rrr=sum_maj(ns)*rdu*aaaa**8
             else
               rrr=0.
             endif
            else
               if (rdu.lt..75) then
                  rrr=sum_pla(ns)
               else
                  rrr=sum_maj(ns)
               endif
            endif
c Here I include the 1.68 DOP filter options:
c  angle with vertical plane: 
            if (zdeg.lt.90) then
              aaaa=plan1(ns)
              aaa1=asin(aaaa)
              if (abs(aaa1).gt.zrad) rrr=0.
            endif
c  angle of semi major with Z and H:
            if (wdeg.lt.90) then
              aaaa=rsig1(ns)
              aaa2=asin(aaaa)
              aaa2=abs(aaa2)
              if (aaa2.gt.wrad.and.aaa2.lt.wrad2) rrr=0.
            endif 

c to avoid numerical problems:
            if (rrr.ge.0..and.rrr.le.1.) then
               degpol(nt,ns,nf)=rrr
            else
               degpol(nt,ns,nf)=0.
            endif

          enddo
        enddo
  
        if (lbug) write(*,*)' trace loop',nt
      enddo


c========================================================
c      W   R   I   T   E       O   U   T   P   U   T    : 
c========================================================
c
      if (lfre) then
         nf=nint((fre-f1)/dfstep) + 1
         fcent=f1+(nf-1)*dfstep 
         open(87,file='fre.asc')
         rewind(87)
         nt=1
         if (lwlenf) wlen=awlenf/fcent
c        write(87,*)fcent,wlen,mocyc,pow,nflen,dt,nt
         do ns = 1,nsamp
           write(87,*)(ns-1)*dt+beg,degpol(nt,ns,nf)
         enddo
         close(87)
      endif
          
cccccccccccccccc
c Binary output:
cccccccccccccccc
      if (lbin.and.dopm.lt.0.) then
        open(unit=ioutdat,form='UNFORMATTED',
     &      file='degpolmat.bin',err=825)
        go to 830
825     write(nstderr,'(a)') 'Can''t open file : degpolmat.bin'
        stop
830     rewind ioutdat
 
        write(nstderr,'(/6x,a)')'Output: degpolmat.bin '
        write(nstderr,'(/6x,a,1x,i4)')
     :       'Number of records =',ntrac
        write(nstderr,'(/6x,a,1x,i4)')
     :       'Number of samples/trace =',nsamp
        write(nstderr,'(/6x,a,1x,i4)')
     :       'Number of frequencies   =',nfreq
 
c write one header line with important parameters.
        if (lwlenf) then
          write(ioutdat)ntrac,nsamp,nfreq,wlenf,pow,f1,f2,
     &     rycle,nocyc,lmed,dfstep
        else
          write(ioutdat)ntrac,nsamp,nfreq,wlen,pow,f1,f2,
     &     rycle,nocyc,lmed,dfstep
        endif
        do nt = 1,ntrac
           do ns = 1,nsamp
             write(ioutdat)(degpol(nt,ns,nf),nf=1,nfreq)
           enddo
        enddo
        close(unit=ioutdat)
      endif
ccccccccccccccccccccccccccccccccc
c ASCII output for test purposes:
ccccccccccccccccccccccccccccccccc
      if (lbug) then
        open(ioutdat,file='degpolmat.asc')
        if (lwlenf) then
          write(ioutdat,*)ntrac,nsamp,nfreq,wlenf,pow,f1,f2,
     &     rycle,nocyc,lmed,dfstep
        else
          write(ioutdat,*)ntrac,nsamp,nfreq,wlen,pow,f1,f2,
     &     rycle,nocyc,lmed,dfstep
        endif
        do nt = 1,ntrac
           do ns = 1,nsamp
             write(ioutdat,*)(degpol(nt,ns,nf),nf=1,nfreq)
           enddo
        enddo
        close(unit=ioutdat)                             
      endif

cccccccccccccccccccccccccc
c Output linear dop:
cccccccccccccccccccccccccc
      write(*,*)
      if (dopm.ge.0.) then
        write(*,*)'OUTPUT: azi_dopm.asc'
        open(19,file='azi_dopm.asc')
      else
        if (llin) then
          open(19,file='degpolmat_lin.asc')
        else if (llim) then
          open(19,file='degpolmat_lim.asc')
        else if (lelm) then
          open(19,file='degpolmat_elm.asc')
        else if (lell) then
          open(19,file='degpolmat_ell.asc')
        else
          open(19,file='degpolmat_dop.asc')
        endif
      endif
      write(*,*)

      if (dopm.ge.0.) then
 191   format(I4,1x,I3,1x,3(I2,1x),I3)
 19    format(F8.3,1x,F14.10,1x,F8.5,1x,F14.5,1x,F6.4,6(1x,F12.7))
 29    format(F8.3,1x,F14.10,1x,F8.5,1x,F14.5,1x,F6.4,6(1x,F12.7),
     & 2(1x,F14.6))
 39    format(F8.3,1x,F14.10,1x,F8.5,1x,F14.5,1x,F6.4,2(1x,F14.6))
 49    format(F8.3,1x,F14.10,1x,F8.5,1x,F14.5,1x,F6.4)
 81    format(F8.3,1x,F14.10,1x,F8.5,1x,F14.5,1x,F6.4,2(1x,F9.4))
 82    format(F8.3,1x,F14.10,1x,F8.5,1x,F14.5,1x,F6.4,2(1x,F14.6),
     & 2(1x,F9.4))
       write(19,191)iyear,jjj,nhour1,nmin1,nsec1,nmsec1
       do nf=1,nfreq
         fff=f1+(nf-1)*dfstep
         do ns = 1,nsamp
          if (degpol(1,ns,nf).gt.dopm) then
            if (lsemi.and..not.lhv) then
            write(19,19)azim(ns,nf),fff,degpol(1,ns,nf),beg+(ns-1)*dt,
     & linear(ns,nf),r1(ns,nf,1),r1(ns,nf,2),r1(ns,nf,3),
     & r3(ns,nf,1),r3(ns,nf,2),r3(ns,nf,3)
            else if (lsemi.and.lhv) then
             ra1=r1(ns,nf,1)
             ra2=r1(ns,nf,2)
             ra3=r1(ns,nf,3)
             rb1=r3(ns,nf,1)
             rb2=r3(ns,nf,2)
             rb3=r3(ns,nf,3)
             if (abs(ra1).gt.0.7071067) then
               rbh=sqrt(rb2*rb2+rb3*rb3)
               rba=1.0-linear(ns,nf)
               rhv=rbh*rba/abs(ra1)
             else
               rah=sqrt(ra2*ra2+ra3*ra3)
               rba=1.0/(1.0-linear(ns,nf))
               rhv=rah*rba/abs(rb1)
             endif
            write(19,29)azim(ns,nf),fff,degpol(1,ns,nf),beg+(ns-1)*dt,
     & linear(ns,nf),r1(ns,nf,1),r1(ns,nf,2),r1(ns,nf,3),
     & r3(ns,nf,1),r3(ns,nf,2),r3(ns,nf,3),rhv,rba
            else if (.not.lsemi.and..not.lzwout.and.lhv) then
             ra1=r1(ns,nf,1)
             ra2=r1(ns,nf,2)
             ra3=r1(ns,nf,3)
             rb1=r3(ns,nf,1)
             rb2=r3(ns,nf,2)
             rb3=r3(ns,nf,3)
             if (abs(ra1).gt.0.7071067) then
               rbh=sqrt(rb2*rb2+rb3*rb3)
               rba=1.0-linear(ns,nf)
               rhv=rbh*rba/abs(ra1)
             else
               rah=sqrt(ra2*ra2+ra3*ra3)
               rba=1.0/(1.0-linear(ns,nf))
               rhv=rah*rba/abs(rb1)
             endif
            write(19,39)azim(ns,nf),fff,degpol(1,ns,nf),beg+(ns-1)*dt,
     & linear(ns,nf),rhv,rba
            else if (.not.lsemi.and..not.lzwout.and..not.lhv) then
            write(19,49)azim(ns,nf),fff,degpol(1,ns,nf),beg+(ns-1)*dt,
     & linear(ns,nf)
            else if (lzwout.and..not.lhv) then
             ra1=r1(ns,nf,1)
             rp1=r2(ns,nf,1)
             aaa1=asin(rp1)
             aaa2=asin(ra1)
             aaa1=abs(aaa1)*rad2deg
             aaa2=abs(aaa2)*rad2deg
            write(19,81)azim(ns,nf),fff,degpol(1,ns,nf),beg+(ns-1)*dt,
     & linear(ns,nf),aaa1,aaa2
            else if (lzwout.and.lhv) then
             ra1=r1(ns,nf,1)
             rp1=r2(ns,nf,1)
             aaa1=asin(rp1)
             aaa2=asin(ra1)
             aaa1=abs(aaa1)*rad2deg
             aaa2=abs(aaa2)*rad2deg
             ra1=r1(ns,nf,1)
             ra2=r1(ns,nf,2)
             ra3=r1(ns,nf,3)
             rb1=r3(ns,nf,1)
             rb2=r3(ns,nf,2)
             rb3=r3(ns,nf,3)
             if (abs(ra1).gt.0.7071067) then
               rbh=sqrt(rb2*rb2+rb3*rb3)
               rba=1.0-linear(ns,nf)
               rhv=rbh*rba/abs(ra1)
             else
               rah=sqrt(ra2*ra2+ra3*ra3)
               rba=1.0/(1.0-linear(ns,nf))
               rhv=rah*rba/abs(rb1)
             endif
            write(19,82)azim(ns,nf),fff,degpol(1,ns,nf),beg+(ns-1)*dt,
     & linear(ns,nf),rhv,rba,aaa1,aaa2
            endif
          endif
         enddo
        enddo 
      else
c write "my header" for gmt scripts:
        write(19,202)beg,(nsamp-1)*dt+beg,dt,f1,f1+(nfreq-1)*dfstep,
     &   dfstep
c
        do nf=1,nfreq
         fff=f1+(nf-1)*dfstep
         do ns = 1,nsamp
          write(19,*)beg+(ns-1)*dt,fff,degpol(1,ns,nf)
         enddo
        enddo
      endif
      close(19)
c
ccccccccccccccccccccccc
c Write info to screen:
ccccccccccccccccccccccc
c  check carefully !
c Parameter check:
      write(*,*)'check carefully ! :'
      write(*,*)'  '
      write(*,*)'DATA PROPERTIES: '
      write(*,*)'  # samples, 3-C traces            : ',nsamp,ntrac
      write(*,*)'  time increment dt                : ',dt
      write(*,*)'  natural frequency increment df   : ',df
      write(*,*)'PROCESSING:'
      if (lwlenf) then
      write(*,*)'  use freq-dependent dop window!'
      write(*,*)'  freq-dep window, power           : ',wlenf,pow
      else
      write(*,*)'  window, power                    : ',wlen,pow
      endif
      write(*,*)'  frequencies f1,f2                : ',f1,f2
      write(*,*)'  # frequencies                    : ',nfreq
      write(*,*)'  # of samples for spec average    : ',nflen2
      write(*,*)'  increment in frequency range     : ',dfstep
      write(*,*)'  ( ',ndfstep,' * ',df,' = ',dfstep,' )'
      write(*,*)'  median                           : ',lmed
      if (lave) then
      write(*,*)'  The spec. matrix is averaged !!! '
      write(*,*)'  # of freq components for average : ',nflen2
      else
      write(*,*)'  The DFT spectra are averaged !!!'
      write(*,*)'  # of freq components for average : ',nflen2
      endif
      if (nocyc.le.0) then
      write(*,*)'  cycles, gauss truncation         : ',rycle,truncat
      write(*,*)'  # samples in 2sigma window for f1: ',rycle/f1/dt+1
      write(*,*)'  # samples in 2sigma window for f2: ',rycle/f2/dt+1
      else
      write(*,*)'  2*std of Gaussian window (time)  : ',nocyc
      write(*,*)'  # samples in 2sigma window       : ',nint(nocyc/dt)+1
      write(*,*)'  truncate Gaussian window at      : ',truncat
      endif
      write(*,*)'  '
      write(*,*)'check carefully the used parameters!'
      write(*,*)
      write(*,*)"---------------- polfre_s1.6 finished ------"//
     & "--------------"
      write(*,*)

      end
c----------------------------------------------------------------------
c                     F U N C T I O N S
c----------------------------------------------------------------------
      real function median(y,num)
c  Author: M.Schimmel (schimmel@ictja.csic.es)
cccccccccccccccccccccccccccccccccccccccc

      parameter (memax=9000)
      real y(1),x(memax)
      integer num
      logical li1,li2,leven

      if (num.gt.memax) stop ' Change dims in function median.'

c array copy:
      do i=1,num
         x(i)=y(i)
      enddo

      rmed=0.
      meven=mod(num,2)
      leven=.false.
      if (meven.eq.0) leven=.true.

      mid=nint(num/2.)
 1    n1=1
      n2=num

 5    rtest=x(mid)
 7    li1=.false.
      do i1=n1,mid
        r1=x(i1)
        if (r1.gt.rtest) then
           r1exchge=r1
           i1exchge=i1
           li1=.true.
           goto 11
        endif
      enddo
 11   continue

      li2=.false.
      do i2=n2,mid,-1
        r2=x(i2)
        if (r2.lt.rtest) then
          r2exchge=r2 
          i2exchge=i2 
          li2=.true.
          goto 12
        endif
      enddo
 12   continue

c exchange values and continue:
      if (li1.and.li2) then
        x(i1exchge)=r2exchge
        x(i2exchge)=r1exchge
        n1=i1exchge+1
        n2=i2exchge-1
        goto 7
      endif

c only i2 reached center:
      if (li1.and..not.li2) then       
        x(i1exchge)=x(mid)
        x(mid)=r1exchge
        n1=i1exchge+1
        n2=num
        goto 5
      endif

c only i1 reached center:
      if (li2.and..not.li1) then       
        x(i2exchge)=x(mid)
        x(mid)=r2exchge
        n1=1
        n2=i2exchge-1
        goto 5
      endif

c i1 & i2 reached center
      if (.not.li1.and..not.li2) rmed=rtest

      if (leven) then
        if (meven.eq.0) then
          mid=mid+1
          meven=1
          rmed1=rmed 
          goto 1
        endif
        rmed=(rmed+rmed1)/2.0
      endif
     
      median=rmed
 
      return
      end
c__________________________________________________________________
        complex FUNCTION dftc(y,nlen,freq,dt)
c
c Discrete Fourier Transformation: (use natural frequencies only!)
c
c  y        is array with time series.
c  nlen     is the length of data vector y.
c  freq     is desired frequency.
c  dt       is time interval between samples in y.
c
c  Author: M.Schimmel (schimmel@ictja.csic.es)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        include 'include_pol.h'
c maxns
        real y(maxns),dt,freq
        double precision twpi,omeg,omega,dfreq,ddt,dy
        double complex cdum,camp,cphs
 
        twpi=8.D0*datan(1.0D0)
        ddt=dt
        dfreq=freq
 
        omeg=-twpi*dfreq*ddt
        cdum=cmplx(0.D0,0.D0)
        do n=1,nlen
           dy=y(n)
           omega=omeg*float((n-1))
           camp=dcmplx(dy,0.)
           cphs=dcmplx(0.,omega)
           cdum=cdum+camp*cdexp(cphs)
        enddo
        dftc=cdum
 
       end 
c------------------------------------------------------------------------
      integer FUNCTION leng(string)
c     leng is the position of the first non blank character
c     in the character variable char
      character string*(*)
      integer l
      l=len(string)
         do i=1,l
           if (string(i:i).eq.' ') go to 10
         enddo
         i=0
  10      leng=i-1
        return
        end

c-------------------------------------------------------------------
c                     S U B R O U T I N E S
c-------------------------------------------------------------------
      SUBROUTINE crosspn(x,y,z)
c
c z = (x cross y)/|x cross y|
c Author: M. Schimmel (schimmel@ictja.csic.es)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision x(3),y(3),z(3),dd,z1,z2,z3
 
      x1=x(1)
      x2=x(2)
      x3=x(3)
      y1=y(1)
      y2=y(2)
      y3=y(3)
      z1=x2*y3-x3*y2
      z2=x3*y1-x1*y3
      z3=x1*y2-x2*y1
      dd=dsqrt(z1*z1+z2*z2+z3*z3)
      if (dd.gt.0.0000001) then
        z(1)=z1/dd
        z(2)=z2/dd
        z(3)=z3/dd
      else
        z(1)=0.0000001
        z(2)=0.0000001
        z(3)=0.0000001
      endif
 
      return
      end
c-------------------------------------------------------------------
      SUBROUTINE mov_proj_med(rs1,rs2,rs3,rout,nsmp,wlen,pow)
c
c  Get med vector in the moving window (wlen samples).
c  Sum the projections of the vector rs? onto the unit med vector.
c  vector rs? is expected to be normalized.
c
c  Input: rs1,rs2,rs3,nsmp,wlen
c  Output: rout
c  Author: M. Schimmel (schimmel@ictja.csic.es)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'include_pol.h'
c  nax=1024
      real rs1(nax),rs2(nax),rs3(nax),rout(nax)
      real d1(nax),d2(nax),d3(nax),pow
      double precision dd,r1,r2,r3,ddum
      double precision rr1,rr2,rr3
      integer nsmp,wlen,nv2
 
      if (nsmp.gt.nax) stop' Change dimensions in mov_proj_med!'
      do i=1,nsmp
         rout(i)=0.
         d1(i)=0.
         d2(i)=0.
         d3(i)=0.
      enddo
      call movmedian_vec_pi(wlen,rs1,rs2,rs3,d1,d2,d3,nsmp)
      nv2=int((wlen-1)/2)
      do i=1,nsmp
         r1=dble(d1(i))
         r2=dble(d2(i))
         r3=dble(d3(i))
         dd=dsqrt(r1*r1+r2*r2+r3*r3)
         if (dd.lt.1.e-7) then
           write(*,*)'Problem with vector norm.'
           dd=1.e-7
         endif
         r1=r1/dd
         r2=r2/dd
         r3=r3/dd
         i1=i-nv2
         i2=i+nv2
         if (i1.lt.1) i1=1
         if (i2.gt.nsmp) i2=nsmp
         dd=0.D0
         do ic=i1,i2
           rr1=r1*dble(rs1(ic))
           rr2=r2*dble(rs2(ic))
           rr3=r3*dble(rs3(ic))
           dd=dd+abs(rr1+rr2+rr3)**pow
         enddo
         ddum=dd/dble(i2-i1+1)
         dd=ddum**pow
         rout(i)=sngl(dd)
      enddo
 
      return
      end
c-----------------------------------------------------------------
      SUBROUTINE mov_proj_mean(rs1,rs2,rs3,rout,nsmp,wlen,pow)
c
c  Get mean vector in the moving window (wlen samples).
c  Sum the projections of the vector rs? onto the unit mean vector.
c  vector rs? is expected to be normalized.
c
c  Input: rs1,rs2,rs3,nsmp,wlen
c  Output: rout
c  Author: M. Schimmel (schimmel@ictja.csic.es)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      include 'include_pol.h'
c nax
      real rs1(nax),rs2(nax),rs3(nax),rout(nax)
      real d1(nax),d2(nax),d3(nax),pow
      double precision dd,r1,r2,r3,ddum
      double precision rr1,rr2,rr3
      integer nsmp,wlen,nv2
 
      if (nsmp.gt.nax) stop ' Change dimensions in mov_proj_mean!'
      do i=1,nsmp
         rout(i)=0.
         d1(i)=0.
         d2(i)=0.
         d3(i)=0.
      enddo
      call movmean_vec_pi(wlen,rs1,rs2,rs3,d1,d2,d3,nsmp)
      nv2=int((wlen-1)/2)
      do i=1,nsmp
         r1=dble(d1(i))
         r2=dble(d2(i))
         r3=dble(d3(i))
         dd=dsqrt(r1*r1+r2*r2+r3*r3)
         r1=r1/dd
         r2=r2/dd
         r3=r3/dd
         i1=i-nv2
         i2=i+nv2
         if (i1.lt.1) i1=1
         if (i2.gt.nsmp) i2=nsmp
         dd=0.D0
         do ic=i1,i2
           rr1=r1*dble(rs1(ic))
           rr2=r2*dble(rs2(ic))
           rr3=r3*dble(rs3(ic))
           dd=dd+abs(rr1+rr2+rr3)**pow
         enddo
         i123=i2-i1+1
         ddum=dd/dble(i123)
         dd=ddum**pow
         rout(i)=sngl(dd)
      enddo
 
      return
      end
c--------------------------------------------------------- 
      SUBROUTINE movmean_vec_pi(nm,x1,x2,x3,y1,y2,y3,nx)
c
c  Consider that the sign of the vector x(x1,x2,x3) is not known.
c  Therefore switch vectors within window whenever their angle is
c  larger than 90 deg with the other vectors in the window.
c
c  INPUT:
c     nm is the window size for the averaging.
c     The window is centered at each sample.
c     x(nx) contains the amplitudes to be averaged.
c  OUTPUT:
c     y(nx) mean values.
c
c   Author: M. Schimmel (schimmel@ictja.csic.es)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'include_pol.h'
c mx=nax
      integer nm,nx
      real x1(nax),x2(nax),x3(nax),y1(nax),y2(nax),y3(nax)
 
      if (nx.gt.nax) stop ' Change dimensions in movmean_vec_pi!'
 
      na=int((nm-1)/2)
      do ix=1,nx
        i1=ix-na
        i2=ix+na
        if (i1.lt.1) i1=1
        if (i2.gt.nx) i2=nx
        in=i2-i1+1
        call mean_pi(in,x1(i1),x2(i1),x3(i1),rm1,rm2,rm3)
        y1(ix)=rm1
        y2(ix)=rm2
        y3(ix)=rm3
      enddo
 
      return
      end
c__________________________________________________________________
      SUBROUTINE movmedian_vec_pi(nm,x1,x2,x3,y1,y2,y3,nx)
c
c  Consider that the sign of the vector x(x1,x2,x3) is not known.
c  Therefore switch vectors within window whenever their angle es
c  larger than 90 deg with the other vectors in the window.
c
c  INPUT:
c     nm is the window size for the averaging.
c     The window is centered at each sample.
c     x(nx) contains the amplitudes to be averaged.
c  OUTPUT:
c     y(nx) median values.
c
c  Author: M. Schimmel (schimmel@ictja.csic.es)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'include_pol.h'
      integer nm,nx
      real x1(nax),x2(nax),x3(nax),y1(nax),y2(nax),y3(nax)
 
      if (nx.gt.nax) stop ' Change dimensions in movmedian_vec_pi!'
 
      na=(nm-1)/2
      do ix=1,nx
        i1=ix-na
        i2=ix+na
        if (i1.lt.1) i1=1
        if (i2.gt.nx) i2=nx
        in=i2-i1+1
        call median_pi(in,x1(i1),x2(i1),x3(i1),rm1,rm2,rm3)
        y1(ix)=rm1
        y2(ix)=rm2
        y3(ix)=rm3
      enddo
 
      return
      end 
c_________________________________________________________________
      SUBROUTINE median_pi(nx,x1,x2,x3,rm1,rm2,rm3)
c
c  find the median vector when the vectors are known with
c  phase pi uncertainty.
c
c  input: nx samples, x1,x2,x3 data vector
c  output: median vector components rm1,rm2,rm3
c
c  Author: M. Schimmel (schimmel@ictja.csic.es)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'include_pol.h'
      real x1(nax),x2(nax),x3(nax)
      real z1(50),z2(50),z3(50)
      real rm1,rm2,rm3,median
      logical lsign
 
      if (nx.gt.50) stop ' Change dims in median_pi.'
      do mm=1,nx
        z1(mm)=x1(mm)
        z2(mm)=x2(mm)
        z3(mm)=x3(mm)
      enddo
      r1=x1(1)
      r2=x2(1)
      r3=x3(1)
      ncount=0
 10   continue
      lsign=.true.
      ncount=ncount+1
      do mm=1,nx
        xx1=z1(mm)
        xx2=z2(mm)
        xx3=z3(mm)
        rr=r1*xx1+r2*xx2+r3*xx3
        if (rr.lt.0.0) then
           lsign=.false.
           z1(mm)=-xx1
           z2(mm)=-xx2
           z3(mm)=-xx3
        endif
      enddo
      r1=median(z1,nx)
      r2=median(z2,nx)
      r3=median(z3,nx)
      if (ncount.ge.4.or.lsign) goto 30
      goto 10
 30   continue
      rm1=r1
      rm2=r2
      rm3=r3
 
      return
      end
c-------------------------------------------------------  
      SUBROUTINE mean_pi(nx,x1,x2,x3,rm1,rm2,rm3)
c
c  find the mean vector when the vectors are known with
c  phase pi uncertainty.
c
c  input: nx samples, x1,x2,x3 data vector
c  output: mean vector components rm1,rm2,rm3
c
c  Author: M. Schimmel (schimmel@ictja.csic.es)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'include_pol.h'
      real x1(nax),x2(nax),x3(nax)
      real z1(nax),z2(nax),z3(nax)
      real rm1,rm2,rm3
      logical lsign
 
      if (nx.gt.nax) stop ' Change dims in mean_pi.'
      do mm=1,nx
        z1(mm)=x1(mm)
        z2(mm)=x2(mm)
        z3(mm)=x3(mm)
      enddo
      r1=x1(1)
      r2=x2(1)
      r3=x3(1)
      ncount=0
 10   continue
      lsign=.true.
      ncount=ncount+1
      rmean1=0.0
      rmean2=0.0
      rmean3=0.0
      do mm=1,nx
        xx1=z1(mm)
        xx2=z2(mm)
        xx3=z3(mm)
        rr=r1*xx1+r2*xx2+r3*xx3
        if (rr.lt.0.0) then
           lsign=.false.
           xx1=-xx1
           xx2=-xx2
           xx3=-xx3
        endif
        rmean1=rmean1+xx1
        rmean2=rmean2+xx2
        rmean3=rmean3+xx3
        z1(mm)=xx1
        z2(mm)=xx2
        z3(mm)=xx3
      enddo
      r1=rmean1/float(nx)
      r2=rmean2/float(nx)
      r3=rmean3/float(nx)
      if (ncount.ge.4.or.lsign) goto 30
      goto 10
 30   continue
      rm1=r1
      rm2=r2 
      rm3=r3
 
      return
      end
cc-------------------------------------------------------------------
       SUBROUTINE state_rot(c1,c2,c3,r1,r2,r1amp,r2amp)
c
c  Perform the phase rotation of the complex state vector to
c  obtain the semi-major y semi-minor axis if the spectral matrix is
c  a pure state.
c  This is in analogy to the computation of the  
c  instantaneous phase for a complex valued 3-comp vector
c  using Morozov & Smithson maximization.
c  Thus, we follow Morozov & Smithson here:
c
c  r1(3)+ i r2(3) is the complex normalized state vector.
c  r1amp and r2amp are the semi-major and semi minor.
c  Author: M. Schimmel (schimmel@ictja.csic.es)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       complex c1,c2,c3
       double complex cd1,cd2,cd3,ca,cb
       double precision pp,pi2,a,b,a1,a2,a3,b1,b2,b3
       double precision r1(3),r2(3)
       real rr,r1amp,r2amp

       pi2=2.D0*datan(1.D0)
       reg=0.0000001
       cd1=c1
       cd2=c2
       cd3=c3
       ca=dcmplx(0.5D0,0.D0)*(cd1*cd1+cd2*cd2+cd3*cd3)
       cb=cd1+cd2+cd3
       cb=cb*cb*dcmplx(0.5D0,0.D0)
       ca=ca+cmplx(reg,0.)*cb
       pp=0.5D0*datan2(dimag(ca),dreal(ca))
       ca=cdexp(dcmplx(0.D0,-pp))
       a1=dreal(ca*cd1)
       a2=dreal(ca*cd2)
       a3=dreal(ca*cd3)
       b1=dimag(ca*cd1)
       b2=dimag(ca*cd2)
       b3=dimag(ca*cd3)

       a=dsqrt(a1*a1+a2*a2+a3*a3)
       b=dsqrt(b1*b1+b2*b2+b3*b3)
       if (b.gt.a) stop 'since minor > major'
       rr=sngl(a1*b1+a2*b2+a3*b3)
       if (abs(rr).gt.0.01) then
         write(*,*)' Problem with state vector.'
         write(*,*)' real and imag parts are not perpendicular!'
         write(*,*)a1,a2,a3,' * ',b1,b2,b3,'  = ',rr
         stop
       endif

       r1(1)=a1/a
       r1(2)=a2/a
       r1(3)=a3/a
       r2(1)=b1/b
       r2(2)=b2/b
       r2(3)=b3/b
       r1amp=sngl(a)
       r2amp=sngl(b)

       end
c-------------------------------------------------------------------------
       SUBROUTINE window2(x,y,nx,ny,std,dt,nc,truncat)
c
c  x   : trace to be filtered
c  y   : windowed trace = x(i) * gauss(nc)
c  std : std of Gaussian (time unit)
c  dt  : sample interval of x.
c  truncat: truncate time series at truncat.
c
c  note y(1) corresponds a x(nc-r1)*truncat
c  Author: M. Schimmel (schimmel@ictja.csic.es)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       include 'include_pol.h' 
       real x(maxns),std, y(maxns)

       rstd=std/dt
       r1=sqrt(log(truncat)*(-2.)) * rstd
       n1=nc-nint(r1)
       n2=nc+nint(r1)
       if (n1.lt.1) n1=1
       if (n2.gt.nx) n2=nx
       ny=n2-n1+1
       rstd=rstd**2
 
       do ii=n1,n2
          rrr=float(ii-nc)
          arg=-0.5*rrr*rrr/rstd
          kk=ii-n1+1
          y(kk)=x(ii)*exp(arg)
       enddo
 
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE getpar(name1,name2,name3,f1,f2,dt,nfreq,ntrac,
     & nsamp,wlen,wlenf,rycle,nocyc,mocyc,pow,lasc,lbin,lsac,lmed,
     & lave,nflen,lfre,fre,lwlenf,llin,llim,lelm,lell,dopm,llh,rlll,
     & zdeg,wdeg,lhv,lsemi,lzwout)
c  Author: M. Schimmel
c===============================================================
c
      character*80 name1, name2, name3
      character par*68, par1*68
      logical lbin,lasc,lsac,lmed,lave,lfre,lwlenf
      logical lhv,lsemi,lzwout
      logical llin,llim,lelm,lell,llh
      integer ntrac, nsamp, nfreq, nflen,wlen,wlenf
      real pow,f1,f2,dt,rycle,nocyc,mocyc,rlll,zdeg,wdeg

cccccccccccccccccccccc
c get the input files:
cccccccccccccccccccccc
c number of arguments: narg
      narg=iargc()

      if(narg.lt.3) call usepol 

      call getarg(1,name1)
      call getarg(2,name2)
      call getarg(3,name3)

ccccccccccccccccccccccccccccccccccccccccc
c Now get the remaining input parameters:
ccccccccccccccccccccccccccccccccccccccccc
c defaults et al.: 
      dt=0.
      f1=0.
      f2=0.
      lasc=.false.
      lbin=.false.
      lsac=.true.
      lmed=.false.
      lave=.false.
      lfre=.false.
      llin=.false.
      llim=.false.
      lelm=.false.
      lell=.true.
      lzwout=.false.
      rlll=1.0
      llh=.false.
      lwlenf=.false.
      lhv=.false.
      lsemi=.false.
      dopm=0.7
      wlen=7
      wlenf=0
      pow=0.
      rycle=3.
      nocyc=0.
      mocyc=0.
      nflen=2
      nfreq=50
      nsamp=0
      ntrac=1
      zdeg=15.
      wdeg=15.

      do iarg=4,narg
        call getarg(iarg,par)
        l=leng(par)
        if (par(1:3).eq.'f1=') then
           par1=par(4:l)
           read(par1,*)f1
        else if (par(1:3).eq.'f2=') then
           par1=par(4:l)
           read(par1,*)f2
        else if (par(1:3).eq.'dt=') then
           par1=par(4:l)
           read(par1,*)dt
        else if (par(1:5).eq.'wlen=') then
           par1=par(6:l)
           read(par1,*)wlen
           lwlenf=.false.
        else if (par(1:6).eq.'wlenf=') then
           par1=par(7:l)
           read(par1,*)wlenf
           lwlenf=.true.
        else if (par(1:4).eq.'ntr=') then
           par1=par(5:l)
           read(par1,*)ntrac
        else if (par(1:4).eq.'nsp=') then
           par1=par(5:l)
           read(par1,*)nsamp
        else if (par(1:4).eq.'nfr=') then
           par1=par(5:l)
           read(par1,*)nfreq
        else if (par(1:4).eq.'fre=') then
           par1=par(5:l)
           read(par1,*)fre
           lfre=.true.
        else if (par(1:5).eq.'wdeg=') then
           par1=par(6:l)
           read(par1,*)wdeg
           wdeg=abs(wdeg)
        else if (par(1:5).eq.'zdeg=') then
           par1=par(6:l)
           read(par1,*)zdeg
           zdeg=abs(zdeg)
        else if (par(1:5).eq.'zwout') then
           lzwout=.true.
           lsemi=.false.
        else if (par(1:5).eq.'dopm=') then
           par1=par(6:l)
           read(par1,*)dopm
           if (dopm.lt.0.) stop 'dopm should be >= 0'
        else if (par(1:4).eq.'pow=') then
           par1=par(5:l)
           read(par1,*)pow
        else if (par(1:6).eq.'cycle=') then
           par1=par(7:l)
           read(par1,*)rycle
        else if (par(1:6).eq.'nocyc=') then
           par1=par(7:l)
           read(par1,*)nocyc
           rycle=0.
        else if (par(1:6).eq.'mocyc=') then
           par1=par(7:l)
           read(par1,*)mocyc
           rycle=0.
        else if (par(1:6).eq.'nflen=') then
           par1=par(7:l)
           read(par1,*)nflen
        else if (par(1:3).eq.'ave') then
           lave=.true.
        else if (par(1:3).eq.'med') then
           lmed=.true.
        else if (par(1:3).eq.'bin') then
           lbin=.true.
        else if (par(1:3).eq.'asc') then
           lasc=.true.
        else if (par(1:3).eq.'sac') then
           lsac=.true.
        else if (par(1:2).eq.'hv') then
           lhv=.true.
        else if (par(1:4).eq.'semi') then
           lsemi=.true.
           lzwout=.false.
        else if (par(1:4).eq.'ell=') then
             par1=par(5:l)
             read(par1,*)rlll
        else
           call usepol
        endif
      enddo

      if (.not.lsac.and.dt.le.0.) stop ' Need sample interval dt= '
      if (mocyc.ne.0.0) nocyc=(mocyc-1.)*dt

      return
      end
c========================================================================
      SUBROUTINE usepol
c  Author: M. Schimmel (schimmel@ictja.csic.es)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      write(*,*)
      write(*,*)'USAGE:  polfre in1 in2 in3  [parameters]'
      write(*,*)
      write(*,*)' in?    :  bin/sac files.'
      write(*,*)'             (use 1=Z, 2=N, 3=E)'
      write(*,*)' '
      write(*,*)' parameters: wlen= wlenf= pow= dt= nsp= ntr= nfr='
      write(*,*)'             nflen= f1= f2= nocyc= cycle= med '
      write(*,*)'             ave fre='
      write(*,*)' '
      write(*,*)'  f1,f2 :  Min & max frequency to be analysed.'
      write(*,*)'  nfr   :  # of frequency samples.'
      write(*,*)'           (pg alters nfr to adjust freq range.)' 
      write(*,*)'  wlen  :  Length of moving DOP window (in samples).'
      write(*,*)'  wlenf :  Freq dependent DOP window: '//
     & 'wlen=wlenf*f2/f '
      write(*,*)'           (wlen,wlenf are integers!)'
      write(*,*)'  bin   :  Use bin i/o (rather than default: sac)'
      write(*,*)'           i: without header, o: 1 line header.'
      write(*,*)'  dt    :  Data time increment (not needed for SAC).'
      write(*,*)'  ntr   :  # of traces in bin files.'
      write(*,*)'  nsp   :  # of samples per trace in bin files.'
      write(*,*)'           Not needed for SAC data whenever the entire'
      write(*,*)'           trace is analyzed. But, can be used to '//
     & 'restrict'
      write(*,*)'           # of samples in SAC data.'
      write(*,*)'  cycle :  # of periods within 2*std width of moving'
      write(*,*)'           Gauss window used to move to the tf-domain.'
      write(*,*)'  nocyc :  2*std width of freq independent Gauss '//
     & 'window.'
      write(*,*)'           A time value (positive & real) is required'
      write(*,*)'           (use mocyc to provide numb. of samples).'
      write(*,*)'  med   :  Use median rather than mean.'//
     & ' This applies to'
      write(*,*)'           smooth the linearity'//
     & ' (eq. 3 in Schimmel &'
      write(*,*)'           Gallart, 2004) and to determine'//
     & ' the vector'
      write(*,*)'           for the vector projections.' //
     & ' Using the mean is'
      write(*,*)'           faster and the differences are'//
     & ' (usually)'
      write(*,*)'           very small. NOTE: median vector '//
     & 'projections'
      write(*,*)'           have been replaced by mean due to a bug.'
      write(*,*)'  ave   :  Average spectral matrix rather than'
      write(*,*)'           spectra (default).'
      write(*,*)'  nflen :  2*nflen+1 frequency components used'
      write(*,*)'           to average (mean) spectrum.'
      write(*,*)'  fre=  :  Write ascii file with DOP for nearest '//
     & 'frequency'
      write(*,*)'           (time, DOP).'
      write(*,*)'  ell=x :  With this option the output is restricted'//
     & ' to' 
      write(*,*)'           linearity < x. On default, the program '//
     & 'outputs' 
      write(*,*)'           the elliptical DOP for all linearity '//
     & 'values.' 
      write(*,*)'  dopm= :  Output is further limited by ell. '//
     & 'dop > dopm'
      write(*,*)'           (dop=0.7 is default).' 
      write(*,*)'  zdeg= :  Set DOP to zero whenever the plane of '//
     & 'particle motion '
      write(*,*)'           deviates from the vertical by more than '//
     & 'zdeg deg.'
      write(*,*)'           Default: zdeg=15. '
      write(*,*)'  wdeg= :  Set DOP to zero whenever the semi-major'//
     & ' vector '
      write(*,*)'           deviates from the horizontal plane or '//
     & 'vertical ax'
      write(*,*)'           by more than wdeg deg. Default: wdeg='//
     & '15.'
      write(*,*)'           zdeg=90, wdeg=90 or more will '// 
     & 'desactivate these'
      write(*,*)'           options.'
      write(*,*)'  zwout :  Output zdeg & wdeg values. Default: no'//
     & ' output.'
      write(*,*)'           semi & zwout do not work together!!!' 
      write(*,*)'  semi  :  Output 6 more columns: semi-major &'//
     & ' semi-minor vector.'
      write(*,*)'           Default: no output.'
      write(*,*)'  hv    :  Output H/V ratio: either |bh|/|a1| or '//
     & ' |ah|/|b1|.'
      write(*,*)'           ah=sqrt(a2*a2+a3*a3) bh=sqrt(b2*b2+b3*b3).'
      write(*,*)'           Further, |a|/|b| or |b|/|a| is output '//
     & 'depending whether'
      write(*,*)'           the semi-major a is parallel to V or H.'  
      write(*,*)'           Use with wdeg. Default: no output.'
      write(*,*)' '
      write(*,*)'OUTPUT file: azi_dopm.asc'
      write(*,*)' '
      write(*,*)'    baz [deg], freq [Hz], dop, time [s], linearity,'//
     &  ' a1, a2, a3, b1, b2, b3, hv, ab, zdeg, wdeg'
      write(*,*)' '
      write(*,*)'    baz : Retrograde particle motion is assumed.'
      write(*,*)'    time: Time in seconds with respect to 1st '//
     & 'sample' 
      write(*,*)'          (time=(ns-1)dt+beg, beg is begin time '//
     & '1st '
      write(*,*)'          sample, ns is sample index, dt is sample '
      write(*,*)'          interval).'
      write(*,*)'    dop : Ell. DOP for motion in a vertical plain.' 
      write(*,*)'          I.e., downweight DOP for motion which is '//
     & 'not in a' 
      write(*,*)'          vertical plain.'
      write(*,*)'    a,b : Normalized semi-major and semi-minor '//
     & 'vector' 
      write(*,*)'          as function of time and frequency. a1,b1 '//
     & 'vertical components.'
      write(*,*)'    hv  : Is H/V ratio as defined above.'
      write(*,*)'    ab  : Is a/b or b/a ratio as defined above.'
      write(*,*)'    zdeg: measured angel in deg as defined above.'
      write(*,*)'    wdeg: measured angel in deg as defined above.'
      write(*,*)
      write(*,*)'DEFAULTS:' 
      write(*,*)' '    
      write(*,*)'  wlen  : 7'
      write(*,*)'  med   : not set.'
      write(*,*)'  ave   : not set.'
      write(*,*)'  nflen : 2.'
      write(*,*)'  nfr   : 50.'
      write(*,*)'  cycle : 3.'
      write(*,*)'  pow   : 4.'
      write(*,*)'  dopm  : 0.7'
      write(*,*)'  zdeg  : 15'
      write(*,*)'  wdeg  : 15'
      write(*,*)' '
      write(*,*)'FINAL REMARKS: This code is a research algorithm '//
     & '(working version)'
      write(*,*)'and has been distributed in the hope to '//
     & 'be useful for your' 
      write(*,*)'research. Please, note I can not take any '//
     & 'warranties, i.e., use'
      write(*,*)'this program on your own risk. If you use '//
     & 'results from this'
      write(*,*)'program in publications then, please, cite '//
     & 'Schimmel & Gallart'
      write(*,*)'(2004) and Schimmel et al. (2011). Thanks! '
      write(*,*)
      write(*,*)'(1) Main theoretical development:'
      write(*,*)' Schimmel & Gallart, Degree of polarization '// 
     & 'filter for frequency-'
      write(*,*)' dependent signal enhancement through noise '//
     & 'suppression, Bull.'
      write(*,*)' Seism.Soc.Am., 94, 1016-1035, doi: '//
     & '10.1785/0120030178, 2004.'
      write(*,*)'(2) Main adaptation/application to study noise:'
      write(*,*)' Schimmel, M., Stutzmann, E., Ardhuin, F., '//
     & 'Gallart, J., Polarized' 
      write(*,*)' Earths Ambient Microseismic Noise, Geochem. '//
     & 'Geophys. Geosyst., '
      write(*,*)' doi:10.1029/2011GC003661, 2011.'
      write(*,*)' '
      write(*,*)'BUGS & COMMENTS: Please, report '//
     & 'bugs! Do not hesitate to'
      write(*,*)'contact me for any questions,'//
     & ' improvements, collaborations,'
      write(*,*)'nice results, ...'

      write(*,*)' '
      write(*,*)'LAST MODIFICATION: 22.09.2015 '
      write(*,*)' '
      write(*,*)'AUTHOR:  M. Schimmel (schimmel@ictja.csic.es) '
      write(*,*)' '
      stop
      return
      end
c---------------------------------------------------------------
       SUBROUTINE sac_read(traces,trace,ntr,ns,dt1,beg)
c
c  Read sac data from file traces into vector.
c  Author: M. Schimmel (schimmel@ictja.csic.es)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'include_pol.h'
      parameter (nstderr=0,indat=10)
      character*80 traces
      integer ntr,ns
      real trace(maxns,maxtr),sig(maxns)

      call rsac1(traces,sig,ns,beg,dt1,maxns,nerr)
      if (nerr.eq.-803) then
        write(*,*)' READ only first maxns samples!',maxns
      else if (nerr.ne.0) then
         write(*,*)' '
         write(*,*)' Problem with SAC trace.'
         write(*,*)' '
         stop
      endif

      if (ns.gt.maxns) ns=maxns
      if (ntr.gt.maxtr) stop 'BIN_READ: decrease # of traces.'

      k=1
      do j=1,ns
        trace(j,k)=sig(j)
      enddo

      return
      end
c__________________________________________________________________
       SUBROUTINE bin_read(traces,trace,ntr,nskip,n1,ns)
c
c  Read time section data from file traces into array trace.
c  n1,ns - index 1st sample and # of samples to read per trace.
c  ntr   - number of traces to be read.
c  nskip - number of records to skip at the beginning of
c          the data set.
c  Author: M. Schimmel (schimmel@ictja.csic.es)
c  ------------------------------------------------------
      include 'include_pol.h'
      parameter (nstderr=0,indat=10)
      character*80 traces
      integer ntr,nskip,n1,ns
      real trace(maxns,maxtr)
 
      n0=n1-1
      if (ns.gt.maxns) stop 'BIN_READ: decrease # of samples.'
      if (ntr.gt.maxtr) stop 'BIN_READ: decrease # of traces.'
 
      open(unit=indat,form='UNFORMATTED',file=traces,err=825)
      goto 830
825   write(nstderr,'(a)') 'Can''t open file : ',traces
      stop
830   rewind indat
 
      write(nstderr,'(/6x,a,a)') 'Time section is file ',traces
      write(nstderr,'(/6x,a,a,1x,i4)')
     :     'Number of records to skip at the ',
     :     'beginning of the data set =', nskip
      write(nstderr,'(/6x,a,1x,i4)')
     :     'Number of traces to read =',ntr
      write(nstderr,'(/6x,a,a,1x,i4)')
     :     'Number of samples to skip at the ',
     :     'beginning of each trace=', n0
      write(nstderr,'(/6x,a,1x,i4)')
     :     'Number of samples/trace =',ns
 
      do k = 1, nskip
         read(indat)skipit
      enddo
      do k = 1,ntr
         read(indat,err=999)(skipit,j=1,n0),(trace(j,k),j=1,ns)
      enddo
      close(unit=indat)
      write(nstderr,'(/6x,a,1x,i4,1x,a,1x,i4)')
     :     'Read traces',nskip+1,' to ',nskip+ntr
 
      return
 999  stop ' Problems in reading bin file.'
      end 
c-----------------------------------------------------------------
      SUBROUTINE movmean(nm,x,y,nx)
c
c  INPUT:
c     nm is the window size for the averaging.
c     The window is centered at each sample.
c     x(nx) contains the amplitudes to be averaged.
c  OUTPUT:
c     y(nx) mean values.
c
c  Author: M. Schimmel (schimmel@ictja.csic.es)
c  ------------------------------------------------------
      include 'include_pol.h'
      integer nm,nx
      real x(nax),y(nax),amean
 
      if (nx.gt.nax) stop ' Change dimensions in movmean!'
 
      na=(nm-1)/2
      do ix=1,nx
        amean=0.
        i1=ix-na
        i2=ix+na
        if (i1.lt.1) i1=1
        if (i2.gt.nx) i2=nx
        do iy=i1,i2
           amean=amean + x(iy)
        enddo
        in=i2-i1+1
        y(ix)=amean/float(in)
      enddo
 
      return
      end
c------------------------------------------------------------------------ 
      SUBROUTINE movmedian(nm,x,y,nx)
c
c  INPUT:
c     nm is the window size for the averaging.
c     The window is centered at each sample.
c     x(nx) contains the amplitudes to be averaged.
c  OUTPUT:
c     y(nx) median values.
c
c  Author: M. Schimmel (schimmel@ictja.csic.es)
c  ------------------------------------------------------
      include 'include_pol.h'
      integer nm,nx
      real x(nax),y(nax),median
      real dum(nax)
 
      if (nx.gt.nax) stop ' Change dimensions in movmedian!'
 
      do ikl=1,nx
        dum(ikl)=x(ikl)
      enddo
      na=(nm-1)/2
      do ix=1,nx
        i1=ix-na
        i2=ix+na
        if (i1.lt.1) i1=1
        if (i2.gt.nx) i2=nx
        in=i2-i1+1
        y(ix)=median(dum(i1),in)
      enddo
 
      return
      end
c_______________________________________________________________________
      SUBROUTINE st_4ward(sig,nsamp,rycle,dt,ctf,nsmp,npow)
c
c INPUT: sig,nsamp,rycle,dt,nsmp,npow
c * sig(nsamp) is real time series with sample interval dt.
c * rycle is the number of periods that equals the 2 std of
c   the Gaussian window.
c * nsmp=2**npow  (==> for FFT)
c * dt is sample interval (time).
c OUTPUT: ctf
c * ctf(nsamp,nsmp2+1) is the complex t-f representation of
c   sig(nsamp).
c   nsamp in ctf is the center sample of 'moving Gauss-window'.
c
c NOTES:
c Please, send bugs, comments, improvements, other routines
c to schimmel@ija.csic.es .
c
c Author: M. Schimmel (schimmel@ictja.csic.es)
c last change 22.08.06
c=================================================================
c nftf=nstf/2 + 1, nstf=2**n, nstf>maxns
      include 'include_pol.h'
      complex csig(nstf),ctf(nstf,nftf),cdum(nstf)
      real sig(maxns),pi,fact,rycle
      integer npow,nsmp,nsamp

      pi=4.0*atan(1.0)
      fact=rycle/2.
      nsmp2=nsmp/2

      do ns=1,nsamp
        csig(ns)=cmplx(sig(ns),0.0)
      enddo
c to go the safe way:
      do ns=nsamp+1,nstf
        csig(ns)=cmplx(0.0,0.0)
      enddo
c
c Perform FFT:
c-------------
      rn=1.
      call clogc(npow,csig,rn,dt)

c=======================
c Start frequency loop:
c=======================
      pp=2.0*pi*pi*fact*fact

c for n=0 (zero frequency nf=1):
c (zero mean t-series assumed!!!!)
c---------------------------------
      do ns=1,nsamp
        ctf(ns,1)=cmplx(0.0,0.0)
      enddo

c and now remaining frequencies:
c--------------------------------
      do nf=2,nsmp2+1
        do mf=1,nsmp
          cdum(mf)=cmplx(0.0,0.0)
        enddo
        do mf=1,nsmp
          argp=-pp*float((mf-1)**2)/float((nf-1)**2)
          argn=-pp*float((nsmp-mf+1)**2)/float((nf-1)**2)
          if (argp.gt.-25.or.argn.gt.-25) then
            mn=mf+nf-1
            if (mn.gt.nsmp) mn=mn-nsmp
             rr=exp(argp)+exp(argn)
             cdum(mf)=csig(mn)*rr
          endif
        enddo

        rn=-1.
        call clogc(npow,cdum,rn,dt)
        do ns=1,nsamp
          ctf(ns,nf)=cdum(ns)
        enddo
      enddo

      return
      end
c---------------------------------------------------------------- 
      SUBROUTINE clogc(n,x,sigm,dt)
 
c--- performs fft on signals with length 2**n and sampling interval
c--- of dt seconds (if in the time domain; notice that dt*df=1/2**n).
c--- the signal is stored in x. it may be complex.
c--- the spectrum is returned in x. it is almost always complex.
c--- a time-to-frequency transform is done with sign=+1. (conform
c--- the convention adopted in aki and richards - the alternative
c--- convention may be obtained by taking complex conjugates after
c--- the call to clogc).
c--- the normalization factor 1./twopi occurs in the frequency-to
c--- time transform (again aki&richards).
c--- normalization is such that physical dimensions are respected.
c--- thus, if the time signal is dimensioned in meters, the
c--- resulting spectral density in x is in meters/hz. for example,
c--- if the time signal is the unit sinc function of width dt, centered
c--- at t=0, the spectral density is dt for all values of the frequency.
c
c--- array locations: if x contains the spectrum, it has the spectrum
c--- for positive frequencies in the first 2**n/2+1 elements, such that
c--- x(1) is at 0 hz, x(2) at df hertz, and x(2**n/2+1) at the nyquist,
c--- where df=1./(2**n*dt) and the nyquist is 1./(2*dt) hz.
c--- the second half of x contains the spectrum for negative frequencies
c--- such that x(2**n) is at -df, x(2**n-1) at -2*df hz etcetera.
c--- if x contains the time signal, x(1) is at time 0, x(2)
c--- at time dt etc.
c
      dimension x(1),m(25)
      complex x,wk,hold,q
      if(sigm.ge.0.) then
        sign=1.
      else
        sign=-1.
      endif
      lx=2**n
      do 1 i=1,n
    1 m(i)=2**(n-i)
      do 4 l=1,n
      nblock=2**(l-1)
      lblock=lx/nblock
      lbhalf=lblock/2
      k=0
      do 4 iblock=1,nblock
      fk=k
      flx=lx
      v=sign*6.283185308*fk/flx
      wk=cmplx(cos(v),sin(v))
      istart=lblock*(iblock-1)
      do 2 i=1,lbhalf
      j=istart+i
      jh=j+lbhalf
      q=x(jh)*wk
      x(jh)=x(j)-q
      x(j)=x(j)+q 
    2 continue
      do 3 i=2,n
      ii=i
      if(k.lt.m(i)) go to 4
    3 k=k-m(i)
    4 k=k+m(ii)
      k=0
      do 7 j=1,lx
      if(k.lt.j) go to 5
      hold=x(j)
      x(j)=x(k+1)
      x(k+1)=hold
    5 do 6 i=1,n
      ii=i
      if(k.lt.m(i)) go to 7
    6 k=k-m(i)
    7 k=k+m(ii)
      if(sign.gt.0.) go to 9
      flx=flx*dt
      do 8 i=1,lx
    8 x(i)=x(i)/flx
      return
    9 do 10 i=1,lx
   10 x(i)=x(i)*dt
      return
      end 
c---------------------------------------------------------------------
c    S U B R O U T I N E S    F R O M    S E I S P A C K
c_____________________________________________________________________
c---------------------------------------------------------------------
c
c---------------------------------------------------------------------
c  double precision eispack routines:
c---------------------------------------------------------------------

      subroutine rsd(nm,n,a,w,matz,z,fv1,fv2,ierr)
c
      integer n,nm,ierr,matz
      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1d(nm,n,a,w,fv1,fv2)
c  tqlrat encounters catastrophic underflow on the Vax
c     call  tqlrat(n,w,fv2,ierr)
      call  tql1d(n,w,fv1,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2d(nm,n,a,w,fv1,z)
      call  tql2d(nm,n,w,fv1,z,ierr)
   50 return
      end

c------------------------------------------------------------------------
      subroutine tred1d(nm,n,a,d,e,e2)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),e2(n)
      double precision f,g,h,scale
c
c     this subroutine is a translation of the algol procedure tred1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix
c     to a symmetric tridiagonal matrix using
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction in its strict lower
c          triangle.  the full upper triangle of a is unaltered.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
c
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0d0
  125    continue
c
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 300
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         e2(i) = scale * scale * h
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         h = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
c
  280    continue
c
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
c
  300 continue
c
      return
      end
c---------------------------------------------------------------------------

      subroutine tred2d(nm,n,a,d,e,z)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end
c------------------------------------------------------------------------

      subroutine tql1d(n,d,e,ierr)
c
      integer i,j,l,m,n,ii,l1,l2,mml,ierr
      double precision d(n),e(n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql1,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the ql method.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        e has been destroyed.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 290 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  210    p = d(l) + f
c     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
c----------------------------------------------------------------

      subroutine tql2d(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc  eigensystem of a complex hermitian matrix. ccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
c
      integer i,j,n,nm,ierr,matz
      double precision ar(nm,n),ai(nm,n),w(n),zr(nm,n),zi(nm,n),
     x       fv1(n),fv2(n),fm1(2,n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex hermitian matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a=(ar,ai).
c
c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex hermitian matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        zr  and  zi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1, fv2, and  fm1  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  htridi(nm,n,ar,ai,w,fv1,fv2,fm1)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
c
         do 30 j = 1, n
            zr(j,i) = 0.0d0
   30    continue
c
         zr(i,i) = 1.0d0
   40 continue
c
      call  tql2(nm,n,w,fv1,zr,ierr)
      if (ierr .ne. 0) go to 50
      call  htribk(nm,n,ar,ai,fm1,n,zr,zi)
   50 return
      end
c=======================================================================
      subroutine htridi(nm,n,ar,ai,d,e,e2,tau)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
      double precision f,g,h,fi,gi,hh,si,scale,pythag
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure tred1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a complex hermitian matrix
c     to a real symmetric tridiagonal matrix using
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex hermitian input matrix.
c          only the lower triangle of the matrix need be supplied.
c
c     on output
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction in their full lower
c          triangles.  their strict upper triangles and the
c          diagonal of ar are unaltered.
c
c        d contains the diagonal elements of the the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c        tau contains further information about the transformations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      tau(1,n) = 1.0d0
      tau(2,n) = 0.0d0
c
      do 100 i = 1, n
  100 d(i) = ar(i,i)
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(ar(i,k)) + dabs(ai(i,k))
c
         if (scale .ne. 0.0d0) go to 140
         tau(1,l) = 1.0d0
         tau(2,l) = 0.0d0
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 290
c
  140    do 150 k = 1, l
            ar(i,k) = ar(i,k) / scale
            ai(i,k) = ai(i,k) / scale
            h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
  150    continue
c
         e2(i) = scale * scale * h
         g = dsqrt(h)
         e(i) = scale * g
         f = pythag(ar(i,l),ai(i,l))
c     .......... form next diagonal element of matrix t ..........
         if (f .eq. 0.0d0) go to 160
         tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
         si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
         h = h + f * g
         g = 1.0d0 + g / f
         ar(i,l) = g * ar(i,l)
         ai(i,l) = g * ai(i,l)
         if (l .eq. 1) go to 270
         go to 170
  160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         ar(i,l) = g
  170    f = 0.0d0
c
         do 240 j = 1, l
            g = 0.0d0
            gi = 0.0d0
c     .......... form element of a*u ..........
            do 180 k = 1, j
               g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
               gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
  180       continue
c
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
               gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
  200       continue
c     .......... form element of p ..........
  220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
  240    continue
c
         hh = f / (h + h)
c     .......... form reduced a ..........
         do 260 j = 1, l
            f = ar(i,j)
            g = e(j) - hh * f
            e(j) = g
            fi = -ai(i,j)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
c
            do 260 k = 1, j
               ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k)
     x                           + fi * tau(2,k) + gi * ai(i,k)
               ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k)
     x                           - fi * e(k) - gi * ar(i,k)
  260    continue
c
  270    do 280 k = 1, l
            ar(i,k) = scale * ar(i,k)
            ai(i,k) = scale * ai(i,k)
  280    continue
c
         tau(2,l) = -si
  290    hh = d(i)
         d(i) = ar(i,i)
         ar(i,i) = hh
         ai(i,i) = scale * dsqrt(h)
  300 continue
c
      return
      end
c=====================================================================
      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
c====================================================================
**** for old version, "send otqlrat from eispack"
** From dana!moler Tue, 1 Sep 87 10:15:40 PDT
** New TQLRAT
      SUBROUTINE TQLRAT(N,D,E2,IERR)
C
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      DOUBLE PRECISION D(N),E2(N)
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG
C
C     This subroutine is a translation of the Algol procedure tqlrat,
C     Algorithm 464, Comm. ACM 16, 689(1973) by Reinsch.
C
C     This subroutine finds the eigenvalues of a symmetric
C     tridiagonal matrix by the rational QL method.
C
C     On input
C
C        N is the order of the matrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E2 contains the squares of the subdiagonal elements of the
C          input matrix in its last N-1 positions.  E2(1) is arbitrary.
C
C      On output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct and
C          ordered for indices 1,2,...IERR-1, but may not be
C          the smallest eigenvalues.
C
C        E2 has been destroyed.
C
C        IERR is set to
C          zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG for  DSQRT(A*A + B*B) .
C
C     Questions and comments should be directed to Burton S. Garbow,
C     Mathematics and Computer Science Div, Argonne National Laboratory
C
C     This version dated August 1987.
C     Modified by C. Moler to fix underflow/overflow difficulties,
C     especially on the VAX and other machines where epslon(1.0d0)**2
C     nearly underflows.  See the loop involving statement 102 and
C     the two statements just before statement 200.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0D0
      T = 0.0D0
      E2(N) = 0.0D0
C
      DO 290 L = 1, N
         J = 0
         H = DABS(D(L)) + DSQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = EPSLON(T)
         C = B * B
         if (c .ne. 0.0d0) go to 105
C        Spliting tolerance underflowed.  Look for larger value.
         do 102 i = l, n
            h = dabs(d(i)) + dsqrt(e2(i))
            if (h .gt. t) t = h
  102    continue
         b = epslon(t)
         c = b * b
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         S = DSQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * S)
         R = PYTHAG(P,1.0D0)
         D(L) = S / (P + DSIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0D0) G = B
         H = G
         S = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
C           Avoid division by zero on next pass
            if (g .eq. 0.0d0) g = epslon(d(i))
            h = g * (p / r)
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0D0) GO TO 210
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0D0) GO TO 130
  210    P = D(L) + F
C     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
c=======================================================================
      subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)
c
      integer i,j,k,l,m,n,nm
      double precision ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
      double precision h,s,si
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure trbak1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a complex hermitian
c     matrix by back transforming those of the corresponding
c     real symmetric tridiagonal matrix determined by  htridi.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction by  htridi  in their
c          full lower triangles except for the diagonal of ar.
c
c        tau contains further information about the transformations.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     note that the last component of each returned vector
c     is real and that vector euclidean norms are preserved.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
c     .......... transform the eigenvectors of the real symmetric
c                tridiagonal matrix to those of the hermitian
c                tridiagonal matrix. ..........
      do 50 k = 1, n
c
         do 50 j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
   50 continue
c
      if (n .eq. 1) go to 200
c     .......... recover and apply the householder matrices ..........
      do 140 i = 2, n
         l = i - 1
         h = ai(i,i)
         if (h .eq. 0.0d0) go to 140
c
         do 130 j = 1, m
            s = 0.0d0
            si = 0.0d0
c
            do 110 k = 1, l
               s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
               si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
  110       continue
c     .......... double divisions avoid possible underflow ..........
            s = (s / h) / h
            si = (si / h) / h
c
            do 120 k = 1, l
               zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
               zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
  120       continue
c
  130    continue
c
  140 continue
c
  200 return
      end
c=========================================================================
      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to 
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying 
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end
c---------------------------------------------------------------
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
c----------------------------------------------------------------
