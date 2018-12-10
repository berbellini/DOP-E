c
c   12 SEP 2000 - prototype the TP, TSV and TSH times
c       Note these are NO OPS - must search and test all rays
c   17 OCT 2002 - Added description of dfile format to usage routine
c   09 JAN 2005 - subroutine chkdep, hm1 not defined changed to dm1
c       baker@usgs.gov
c   12 DEC 2007 - set header values evlat, evlon, stlat, stlon to -12345
c       for compatibility of resulting SAC traces
c
c NOTES
C       2. Look at dimensions of 20 which used to be maximum number
C           of layers
C       3. Look at dimensions of 40 which may be related to the 20?
c
c  genray96 - asymptotic dislocation response using generalized ray 
c          theory.  modified bessel functions are represented by 
c          their asymptotic series.  explicit near field terms are 
c          also included.  all wave fields are represented(p,sv,sh). 
c          output consists of specified displacements for the 4 
c          (vertical and radial) and 2 (tangential) fundamental 
c          green's functions. 
c
c          the output gives the arrival time of the earliest
c          arrival, but the traces begin 5*dt seconds prior
c          to this
c
c  computational parameters 
c        dt - final time spacing for response 
c        npts - number of points in time series
c        dist - epicentral distance
c        t0 and vred define fime of first sample
c           if vred .gt.0 then  tfirst = t0 + dist/vred
c           else   tfirst = t0
c
c        nasym=number of asymptotic expansion terms computed 
c              (=1, decoupled far-field assumptions) 
c  model parameters 
c        mmax - number of layers
c        vp(i),vs(i),rho(i),th(i) - p, s velocity, density, and 
c            layer thickness( km/sec, gm/cc, km) 
c        nh(j) - ray path ; for each j specify layer 
c            number each segment is in. 
c        nd - number of ray segments in ray description. 
c        nm(j) - mode description for ray path. p=5, sv=3, sh=4. 
c        ndeg - degeneracy of each ray (see hron(1971),BSSA p 765-779) 
c            in the case when the source and receiver are in the same 
c            layer, set ndeg lt 0(zero) for upgoing rays from 
c            the source.(except for direct ray) 
c        hs - source depth. 
c        hr - receiver depth 
c        nsurf =1, receiver on surface 
c              =0, buried receiver. 
c        ttmax - maximum record length relative to start time arrival. 
c                rays with travel times greater than (tfirst+ttmax) 
c                will be disregarded.  
c        ddtt -  parameter for spacing in complex p vs t contour
c                to save computer time the p(t) is not computer
c                at intervals of dt seconds. As a matter of fact
c                the program as received had a very crude search
c                This is a problem near the Rayleigh wave pole
c                where poor sampling yields poor waveform
c                The original code lines of the form
c                      tn=tn+dt*5.
c                and
c                      tn=tn+dt*(5+l*l)
c                are have the 5 replaced by ddtt, e.g.
c                      tn=tn+dt*ddtt
c        tele -  The original program had a quadratic expansion
c                of the time spaceing in contor by the l*l
c                mentioned above.  This is modified
c                so that the code not reads
c                      tn=tn+dt*(ddtt + tele*l)
c                tele=0.0 and ddtt=0.0 yields the default
c                values of tele=0.0 and ddtt=5.0
c                Trial and error will specify the best values
c
c
c    REFERENCE:
c
c    Helmberger, D. V. (1968). The crust-mantle transition in the
c          Bering Sea, Bull. Seism. Soc. Am. 58, 179-214.
c
c
c
c    PROGRAM ORGANIZATION
c
c    For the purpose of those who might need to compile and
c    load this program using overlays, the following is a
c    list of all subroutines and functions used in the program
c    in the order of appearance and call
c
c
c   genray96
c       mchdep
c       gcmdln
c           mnmarg
c           mgtarg
c           usage
c       getmod
c       lgstr
c       chkdep
c       set
c       goamd
c       npow2
c       trav
c       pnot
c           cagcon
c               cr
c           dtdp
c               cr
c       ttime
c           cagcon
c               cr
c       gwamd
c       contor
c           compt
c               cagcon
c                   cr
c               dtdp
c                   cr
c           dcpg
c               put
c               get
c           get
c           put
c           realt
c               cagcon
c                   cr
c       model
c           dcpg
c               put
c               get
c           dtdp
c               cr
c           gramd
c           gwamd
c           integ
c               fcon
c                   dcpg
c                       put
c                       get
c                   four
c           interp
c               get
c               put
c               area2
c           prduct
c               gencof
c                   tranm
c                       cr
c                   refft
c                       cr
c           radpat
c               cr
c           vrtang
c               area3
c                   get
c               bessel
c                   bessl
c               dcpg
c                   put
c                   get
c               put
c               dstdps
c                   cr
c               recev
c                   cr
c       convlv
c           fcon
c               four
c           gramd
c           wrhd96
c           wrtr96
c       setcl
c       gcamd
c
c   Open files
c   4   used to read in ray descriptions, then closed
c   4   used to read in distance file
c   1   used to read in model and then closed
c   2   used for pseudo virtual storage is last lines are
c           uncommented
c   8   output file genray96.grn
c   9   output file genray96.tim
c
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        integer NLAY
        parameter (NLAY=200)
        common/isomod/th(NLAY),vp(NLAY),vs(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep

        common/modlly/mmax
        common/intcon/ddtt,tele
c
c   tfirst  R*8 - time of first arrival
c   r   R*8 - distance
c   p0,t0   R*8 - ray parameter and time of direct arrival
c   p1,t1   R*8 - ray parameter and time of fastest head wave
c   p2,t2   R*8 - ray parameter and time of second fastest head wave
c
        real*8 tfirst,r,p0,t0,p1,t1,p2,t2,tmax
        common/rays/nh(100),nm(100),ndeg,nd 
        common/resp/a(1024),tfust
        common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
        logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
        common/virt/idbl1,idbl2,idbl3,idbl4,idbl5,idbl6,idbl7
        common/vrt1/jdbl1,jdbl2,jdbl3,jdbl4,jdbl5,jdbl6,jdbl7
        integer*4 set
        integer jsrc(21)
        integer*4 dogrn
        character*80 dfile, mname, title
        logical ext
        integer*4 mmax, iunit, iiso, iflsph
        integer*4 ierr
        integer*4 lgstr
        logical jdosh, jdosv, jdop
        logical srcflu, recflu
        common/coff/it(100),nup1(100) 
        common/travel/alp(NLAY),als(NLAY),ndsv,nup 
c
c   control for selecting up/down going rays at the receiver and source
c
        integer*4 isupdn
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        logical ldoray
        logical ltime
c
c   call machine dependent initialization
c
        call mchdep()
c
c   parse command line arguments
c
        call gcmdln(dogrn,dfile,ddtt,tele,nasym,isupdn,ltime)
c
c   dogrn   I*4 - 0 do ALL (default)
c             1 do only EQ + Explosion
c             2 do only Explosion plus point force
c   isupdn  I*4 - selection of rays from source
c               1 = down, 0 = both, -1 = up
c   ltime   L   - .true. only output arrival times 
c               else synthetic
c
c
c   open genray96.ray
c
        inquire(file='genray96.ray',exist=ext)
        if(.not. ext)then
            write(LER,*)'Ray description file ray.out not located'
            go to 9000
        endif
        open(4,file='genray96.ray',access='sequential',form='formatted',
     1      status='unknown')
        rewind 4
        read(4,'(a)',end=9999,err=9999)mname
        read(4,*,end=9999,err=9999)hs, hr
        lmnm = lgstr(mname)
c
c   get the earth model
c
        inquire(file=mname,exist=ext)
        if(.not. ext)then
            write(LER,*)'Model file not located'
            go to 9000
        endif
        write(LOT,*)'Model name: ',mname(1:lmnm)
        call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.true.)
c
c   check the appropriateness of the model file
c
        if(ierr.lt.0)go to 9999
c
c   error checking
c
        if(idimen.ne.1)then
            write(LER,*)'1-D velocity model required'
            go to 9000
        endif
        if(icnvel.ne.0)then
            write(LER,*)'Constant velocity model required'
            go to 9000
        endif
        if(iiso.ne.0)then
            write(LER,*)'Isotropic velocity model required'
            go to 9000
        endif
        if(iflsph.ne.0)then
            write(LER,*)'Flat earth velocity model required'
            go to 9000
        endif
c
c   IF either the source or the receiver are in the model halfspace,
c   augment the model by adding one more layer
c
        call chkdep(hs,hr,srcflu,recflu)
c
c   UNTIL transmission and reflection coefficients are correct
c   for solid/fluid or fluid/fluid use artifice of making
c   shear velocity very small
c   10 23 99
c  
        do 1222 i=1,mmax
            if(vs(i).eq.0.0)vs(i) = 0.00001
 1222   continue
        j = lgstr(mname)
        l = lgstr(title)
        write(LOT,*)'Model         :',mname(1:j)
        write(LOT,*)'Title         :',title(1:l)
        write(LOT,*)'Depth source  :',hs
        write(LOT,*)'Depth receiver:',hr
        if(hr .le. 0.0)then
            nsurf = 1
        else
            nsurf = 0
        endif
        write(LOT,*)'nsurf   :',nsurf
        if(dosud)then
            write(LOT,*)'isupdn    :',isupdn
            write(LOT,*)'spup      :',spup
            write(LOT,*)'ssup      :',ssup
            write(LOT,*)'spdn      :',spdn
            write(LOT,*)'ssdn      :',ssdn
        endif
c
c   open pseudo virtual large arrays
c   set opens unit1 and creates junkfil for virtual storage
c   idbl1 is for gcb, idbl2 denotes dr, etc
c   gcb = store 2400 complex*16 values
c
        idbl1 = set(4800)
        jdbl1 = 2400
c
c    dr store 2400 real*8 values
c
        idbl2 = set(2400)
        jdbl2 = 2400
c
c    pc store 2400 complex*16 values
c
        idbl3 = set(4800)
        jdbl3 = 2400
c
c    tc store 2400 real*8 values
c
        idbl4 = set(2400)
        jdbl4 = 2400
c
c   set up array for temporary storage of fft values
c
        idbl5 = set(4096)
c
c   open data files 
c
        call goamd()
c
c   open the output file for the ray descriptions
c
        open(8,file='genray96.grn',form='formatted',status='unknown') 
        rewind 8
c
c   open the output file for travel times
c
        open(9,file='genray96.tim',form='formatted',status='unknown',
     1      access='sequential')
        rewind 9
        write(9,'(a)')mname(1:lmnm)
        write(9,'(2f15.5)')hs,hr
        l=1 
c
        write(LOT,201) prnta,prntb,prntc,prntd 
  201   format(' ','prnta=',l5,2x,'prntb=',l5,2x,'prntc=',
     1      l5,2x,'prntd=',l5) 
        write(LOT,202) prnte,prntf,prntg,prnth 
  202   format(' ','prnte=',l5,2x,'prntf=',l5,2x,'prntg=',
     1      l5,2x,'prnth=',l5) 
        write(LOT,205) nasym 
  205   format(' ','nasym=',i5) 
c
c 
c
c
c   set up the jsrc array
c
        do 1000 i=1,16
            jsrc(i) = 1
 1000   continue
c
c   enforce output by source definition
c
        if(dogrn.eq.1)then
            jsrc(11) = 0
            jsrc(12) = 0
            jsrc(13) = 0
            jsrc(14) = 0
            jsrc(15) = 0
        else if(dogrn.eq.2)then
            jsrc( 1) = 0
            jsrc( 2) = 0
            jsrc( 3) = 0
            jsrc( 4) = 0
            jsrc( 5) = 0
            jsrc( 6) = 0
            jsrc( 7) = 0
            jsrc( 8) = 0
        endif
        jsrc(16) = 0
        jsrc(17) = 0
        jsrc(18) = 0
        jsrc(19) = 0
        jsrc(20) = 0
        jsrc(21) = 0
            
        if(ddtt.le.0.0)ddtt=5.0
        if(tele.le.0.0)tele=0.0
c
c   read in the ray description
c   The number of rays, jm, is decided by the descriptions
c   read
c
c   read in the ray definitions once, and store them on unit 3
c   This removes a dimension constraint on the number of rays
c
            open(3,status='scratch',access='sequential',
     1          form='unformatted')
            rewind 3
            jm = 0
            jdosh = .false.
            jdosv = .false.
            jdop  = .false.
 2000   continue
            read(4,103,end=97,err=97) nn, ndeg
            read(4,104,end=97,err=97) (nh(j),j=1,nn)
            read(4,104,end=97,err=97) (nm(j),j=1,nn)
            jm = jm + 1
            write(3) nn, ndeg
            write(3) (nh(j),j=1,nn) 
            write(3) (nm(j),j=1,nn) 
            if(nm(1).eq.3 )jdosv = .true.
            if(nm(1).eq.4 )jdosh = .true.
            if(nm(1).eq.5 )jdop  = .true.
        go to 2000
   97   continue
c
c   jsodh jdosv and jdop will control final output, e.g., the
c   jsrc (16) flags. Only turnoff previously set flags
c
        if(recflu)jsrc(16) = 1
        if(jdosh .and. (.not. jdop) .and. (.not. jdosv) )then
            jsrc(1) = 0
            jsrc(2) = 0
            jsrc(3) = 0
            jsrc(4) = 0
            jsrc(6) = 0
            jsrc(7) = 0
            jsrc(9) = 0
            jsrc(10) = 0
            jsrc(11) = 0
            jsrc(12) = 0
            jsrc(13) = 0
            jsrc(14) = 0
            jsrc(16) = 0
        else if(jdosv .and.  jdop .and. (.not. jdosh) )then
            jsrc(5) = 0
            jsrc(8) = 0
            jsrc(15) = 0
        endif
        write(LOT,*)'dogrn:',dogrn
        write(LOT,*)'jsrc :',jsrc 
        write(LOT,*)'jdop     :',jdop
        write(LOT,*)'jdosv    :',jdosv
        write(LOT,*)'jdosh    :',jdosh
        close (4)
        
  103   format(2i5)
  104   format(5x,14i5) 
c
  108   format(' ',/) 
  109   format(' ','xxxxxx station number',i4,2x,'range=',f10.3,
     1      ' tfirst=',f10.3,' tmax=',f10.3)
  110   format(' ','ray number',i4) 
  106   format(' ',2i5) 
  107   format(' ',5x,14i5) 
  212   format(' ','ray',i4,2x,
     1      'has been disregarded because of travel time ') 
  214   format(' ','p0=',d15.8,2x,'t0=',d15.8,2x,'p1=',
     1      d15.8,2x,'t1=',d15.8 ) 
c
c   begin processing
c
        open(4,file=dfile,access='sequential',form='formatted',
     1      status='unknown')
        rewind 4
        ns = 0
   99   continue
            rewind 3
            read(4,*,err=9999,end=9999)dist,dt,npts,t0,vred
            dist = abs(dist)
c
c   dist    R*4 - desired distance from the source
c           cannot be zero
c   dt  R*4 - desired sampling interval
c   npts    R*4 - desired number of points in the time series
c   t0,vred R*4 - control for the time of the first sample
c           if vred > 0, then the
c               tfirst = t0 + dist/vred
c           else
c               tfirst = t0
c
            if(vred.gt.0.0)then
                tfrst = t0 + dist/vred
            else
                tfrst = t0
            endif
c
c   check for power of 2
c
            call npow2(npts)
            np = npts
            ttmax=(np-1)*dt
            write(LOT,211) ttmax ,ddtt,tele, dt
  211       format(' ','ttmax=',f10.3,' ddtt=',f10.3,' tele=',f10.3,
     1          ' dt=',f10.3)

            write(LOT,108) 
            ns = ns + 1
            r= dist
c
c       add a little overhead, but read all ray descriptions, 
c       compute time of first arrival, use this to cull rays
c
            tfirst=1.0d+038
            do 998 iray=1,jm
                read(3) nn, ndeg
                nd=nn
                read(3) (nh(j),j=1,nn)
                read(3) (nm(j),j=1,nn)
                if(iray.eq.1)then
                    vps = vp(nh(1))
                    vss = vs(nh(1))
                    rhs = rho(nh(1))
                endif
                call trav(hs,hr,.false.,isupdn,ldoray)
                if(.not. ldoray)then
                    write(LOT,*)'Ray:',iray,
     1              ' rejected because of type and direction'
                    go to 998
                endif
                call pnot(p0,t0,r,.false.,dt)
                call ttime(p0,p1,t1,p2,t2,r,.false.)
                if(t1.lt.tfirst) tfirst=t1
c-
c   output travel time information
c
                write(9,'(i5)')iray
                write(9,103) nn, ndeg
                write(9,104) (nh(j),j=1,nn)
                write(9,104) (nm(j),j=1,nn)
                write(9,'(5e15.7)')r,p0,t0,p1,t1
  998       continue
            if(ltime)go to 99
            tfust = tfirst
            if(tfust.gt.tfrst)then
                tfust = tfrst
            endif
            tmax=tfust+ttmax
            do 1313 i=1,1024
                a(i)=0.0
 1313       continue
            do 113 k=1,6 
                do 114 j=1,4 
                    jrec= (j-1)*6 + k
                    call gwamd(a,jrec)
  114           continue
  113       continue
            ncl1=0
            write(LOT,109) ns, dist,tfust,tmax
            rewind 3
            do 98 iray=1,jm 
                write(LOT,110) iray 
                read(3) nn, ndeg
                nd=nn
                read(3) (nh(j),j=1,nn) 
                read(3) (nm(j),j=1,nn) 
c
c           trav determines ray information in relation 
c           to earth model 
c
                call trav(hs,hr,prnta,isupdn,ldoray) 
                if(.not. ldoray)go to 98
c
c           finding p0 
c
                call pnot(p0,t0,r,prntb,dt) 
c
c           determine first arrival time for ray 
c
                call ttime(p0,p1,t1,p2,t2,r,prntc) 
c
c           p0 is ray parameter of direct wave
c           t0 is travel time of this wave
c           p1 is ray parameter of possible refracted arrival
c           t1 is the travel time of this arrival
c           This is important since the Cagniard contour leaves
c           the real p-axis at p=p0, but the coutour should begin
c           at the least of p0,p1 in order to get the refracted
c           arrival contributions
c           cull rays with respect to travel time 
c
                write(LOT,106)nn, ndeg
                write(LOT,107) (nh(j),j=1,nn)
                write(LOT,107) (nm(j),j=1,nn)
                write(LOT,214) p0,t0,p1,t1
                ttot=t1-tfirst
                if(ttot.le.ttmax)then
c
c               ray is successful
c
                    ncl1=ncl1+1
c 
c               finding contour path on the 
c               complex p plane 
c 
                    call contor(p0,t0,p1,
     1                      t1,r,dt,tmax) 
c
c               computing ray response 
c               without 1.0/sqrt(t) dependence. 
c 
                    lis=nh(1) 
                    mrex=nm(1) 
                    call model(dt,np,nsurf,r,nasym,mrex)
                else
                    write(LOT,212) iray 
                endif
   98       continue 
c
c       convolving 1.0/sqrt(t) dependence and 
c       scaling for final step resonse 
c
            call convlv(np,dt,jsrc,r,hr,hs,tfust,mname,
     1          srcflu,refdep,recflu,vps,vss,rhs)
            write(LOT,112) ncl1 
        go to 99
 9999   continue
        write(LOT,108) 
        write(LOT,111) 
  111   format(' ','all cases completed') 
  112   format(' ','number of successful rays=',i4) 
        close (8) 
        close (9) 
        call setcl()
        call gcamd()
        close (3)
 9000   continue
        close (4) 
        end

        subroutine area2(v1,v,v2,dt,a) 
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        real*8 v1,v2,a,tb,ta,s,v 
        common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
        logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
c
c   using the theoretical area a determined from area3, this 
c   subroutine adjusts the amplitude at t0 such that the area 
c   between t1 and t2 agrees with a. this assumes a trapezoidal 
c   technique of evaluating area (such as in four for w=0) 
c
c   Area = a = 0.5*v1*d1 + 0.5*v2*dt + v*dt 
c
        s=a 
        tb=-dt 
        ta=dt 
        v=(2.*s + v1*tb - v2*ta)/(ta-tb) 
        if(prntf) then
            write(LOT,100) a,tb,ta,v1,v2,v 
        endif
  100   format(' ','from area2: a=',d15.8,2x,'tb=',d15.8,2x,
     1      'ta=',d15.8,2x,'v1=',d15.8,2x,'v2=',
     2      d15.8,2x,'v=',d15.8) 
        return 
        end 

        subroutine area3(g,a,nt0)
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        double complex g
        real*8 a1,a2,a,dt,b
        real*8 get
        REAL*8 DREAL
        common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth
        logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth
        common/virt/idbl1,idbl2,idbl3,idbl4,idbl5,idbl6,idbl7
c
c   this subroutine finds the theoretical area of 
c   response between adjacent
c   time points around tO. a1 refers to the area 
c   before tO and a2 to that
c   after tO
c
        a =0.0d+00
        a1=0.0d+00
        a2=0.0d+00
        if(prntg) write(LOT,*) ' area3 =  nt0 ', nt0
        if(nt0.ne.1) then
            dt=get(idbl4,nt0)-get(idbl4,nt0-1)
            if(prntg) write(LOT,*) ' 1 dt = ',dt,'  g  ', g
            if(prntg) write(LOT,*) get(idbl4,nt0) , get(idbl4,nt0-1)
            b=dimag(g)
            a1=2.0d+00*dsqrt(dt)*b
c
c   now get contribution from other points not near
c   singularity
c
            do 11 i=1,9
            a1=a1+(get(idbl2,nt0-11+i)+get(idbl2,nt0-10+i))*
     1            (get(idbl4,nt0-10+i)-get(idbl4,nt0-11+i))*0.5d+00
   11       continue
        endif
        dt=get(idbl4,nt0+1)-get(idbl4,nt0)
        if(prntg) write(LOT,*) ' 2 dt = ',dt,'  g  ', g
        b=dreal(g)
        a2=2.0*dsqrt(dt)*b
c
c   now get contribution from other points not near singularity
c
        do 21 i=1,9
        a2=a2+(get(idbl2,nt0+i)+get(idbl2,nt0+i+1))*
     1      (get(idbl4,nt0+i+1)-get(idbl4,nt0+i))*0.5
   21   continue
        a=a1+a2
c
c   finally use a final correction because singularity correction assume
c   that there is no initial level between nt0-1 and nt0+1
c
c   here assumes that t(nt0)-t(nt0-1) = t(nt0+1)-t(nt0)
c
c   a=a+get(idbl2,nt0+1)*dt
c   if(nt0.gt.1) then
c       a=a+get(idbl2,nt0-1)*dt
c   endif
c
        if(prntf) write(LOT,100) a,a1,a2
  100   format(' ','from area3, a=',d15.8,5x,'a1=',d15.8,5x,'a2=',d15.8)
        return
        end

        subroutine bessl(z,nas,nord,pbsl) 
        double complex pbsl,z
c
c   this routine computes each term of the asymptotic expansion 
c   for modified bessel functions and derivatives of bessel 
c   functions(less factors of 1/s) 
c
c   Abromowitz and Stegun
c
c   Kn(z) = (pi/2z)**(1/2) exp(-z) 
c   { 1 + (u-1)/(8z) + (u-1)(u-9)/2!(8z)**2
c           + (u-1)(u-9)(u-25)/3!(8z)**3 + ...
c   u = 4n**2
c
c
c   r   R*4 - distance
c   p   C*8 - complex ray parameter (z=pr)
c   nas I*4 - number of the term, e.g., in above (u-1)/(8z) is
c               nas=2
c   nord    I*4 - order of Bessel function, e.g., n in equation
c   nterm   I*4 - 1 evaluate Kn(z)
c             2 evaluate dKn(z)/z
c   pbsl    C*8 - compex value of evaluated function
c
            if(nas.le.0)then
                pbsl=dcmplx(0.0d+00,0.0d+00)
            else
                pbsl=dcmplx(1.0d+00,0.0d+00)
            endif
            if(nas.le.1) then
                return 
            endif
            hmu=4.0*nord*nord 
            do 10 i=2,nas 
                k1=i-1 
                k=2*k1-1 
                pbsl=pbsl*(hmu-k**2)/(8.*z*k1) 
   10       continue
        return
        end

        subroutine bessel(r,nas,nord,nterm,p,pbsl) 
        double complex pbsl,p,pbsl1 
        real*8 r 
        double complex z
        
        z = p*r
c   pbsl = dcmplx(1.0d+00,0.0d+00)
c   if(nas.eq.1)return
        if(nterm.eq.1)then
c
c       evaluate K(nord)(pr)
c
            call bessl(z,nas,nord,pbsl) 
        else if(nterm.eq.2)then
c
c       evaluate dK(nord)(pr)/d(pr)
c
            call bessl(z,nas,nord-1,pbsl) 
            call bessl(z,nas-1,nord,pbsl1) 
            pbsl =  pbsl + (nord)*pbsl1/z
        else if(nterm.eq.3)then
c
c       evaluate (nord) K(nord)(pr)/(pr)
c
            call bessl(z,nas-1,nord,pbsl) 
            pbsl =         (nord)*pbsl/z
        endif
        return
        end

        function cagcon(p,r) 
        double complex cagcon,cr,p,a,ea,eb 
c
c   computes complex time as a function of complex 
c   ray parameter. 
c
        integer NLAY
        parameter (NLAY=200)
        common/isomod/th(NLAY),vp(NLAY),vs(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/travel/alp(NLAY),als(NLAY),nd,nup 
        real*8 r 
        a=dcmplx(0.0d+00,0.0d+00) 
        do 10 i=1,nd 
            ea=dcmplx(0.0d+00,0.0d+00)
            eb=dcmplx(0.0d+00,0.0d+00)
            if(alp(i).gt.0.0d+00) ea=cr(p,dble(vp(i)) )
            if(als(i).gt.0.0d+00) eb=cr(p,dble(vs(i)) )
            a=a + (ea*alp(i) +eb*als(i))*th(i) 
   10   continue
        cagcon=p*r + a 
        return 
        end 


      subroutine compt(tc,pc,p0,dt,r)
c
c    subroutine to solve for complex pc corresponding to a particular
c    value of tc
c
      real*4 dt
      real*8 tc,r
      double complex pc,cagcon,p,tp,p0
     $             ,dtdp,dtdpp
c
c    use Newton iteration
c
      p=p0+dcmplx(0.0d+00,0.0001d+00)
  100 continue
            tp=cagcon(p,r)-tc
            if(cdabs(tp).lt.0.01d+00*dble(dt)) go to 1000
            dtdpp=dtdp(p,r)
            p=p-tp/dtdpp
      go to 100
 1000 continue
      pc=p
      return
      end

      subroutine contor(p0,t0,p1,t1,r,dt,tmax)
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
      real*4 dt
      real*8 p0,t0,p1,t1,ddt,tc,tmax,eps,r,get,tsvp,t0mdt
     1      ,ps,pe
      common/intcon/ddtt,tele
      common/ans/ncp,nt0 
      common/virt/idbl1,idbl2,idbl3,idbl4,idbl5,idbl6,idbl7
      common/vrt1/jdbl1,jdbl2,jdbl3,jdbl4,jdbl5,jdbl6,jdbl7
      double complex pc,p 
      common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
      logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
c
c    this subroutine finds the path of the cagniard contour on the 
c    complex p plane. ddtt sets the courseness of 
c   the search when the contour
c    leaves real p-axis. if one is interested in 
c   only refraction not in the
c    surface waves set ddtt larger. This reduces 
c   the number of ray segments to
c    be computed and thus increases computer time. 
c   Langston used ddtt=5 for pg,
c    Pn etc.
c
        if(prntd) write(LOT,100) 
  100   format(' ','from contor') 
c
c    default value of nt0 if t0 is not reached. 
c
      nt0= -1 
      ncownt=0 
      maur=0 
      eps=1.0d-07
      ddt=dble(dt)
      k=1
      if(dabs(p0-p1).le.eps) go to 50
c
c    take special care in the situation that t0-t1 < dt 
c   in this case there
c    is no special sampling
c
      if((t0-t1).ge.ddt) go to 51
      imin=1
      ps=0.0d+00
      tsvp=t0-dt
      go to 151
   51 continue
      imin=2
c
c    find p values on the real axis and at this point we 
c   know the following
c
c    t1       t0
c    p1       p0
c
c    we want points closely spaced about p0 and p1 and 
c   so we sample points more
c    densely about p1 and p0 and interpolate between them.
c
      pc=dcmplx(p1,0.0d00)
      tc=t1
      call dcpg(pc,k,1,idbl3)
      call put(idbl4,k,tc)
      if(prntd) write(LOT,105) pc,tc,k
      k=k+1
c
c    sample densely about t1 but never go beyond t0-dt
c
      t0mdt=t0-dt
      tc=t1
   52 continue
            tc=tc+0.5d+00*dble(dt)
            if(tc.ge.t0mdt) go to 53
            call realt(tc,pc,p1,p0,0.5*dt,r)
            call dcpg(pc,k,1,idbl3)
            call put(idbl4,k,tc)
            if(prntd) write(LOT,105) pc,tc,k
            if(k.ge.jdbl3) go to 80
            k=k+1
      go to 52
   53 continue
c
c    get ray parameter for t0-dt
c
      tc=t0-dt
      call realt(tc,pc,p1,p0,0.5*dt,r)
      tsvp=tc
      call dcpg(pc,k,1,idbl3)
      call put(idbl4,k,tc)
      if(prntd) write(LOT,105) pc,tc,k
      if(k.ge.jdbl3) go to 80
      k=k+1
c
c    get 10 samples from t0-dt to t0-0.1*dt
c
      ps=p1
  151 continue
      pe=p0
      do 56 i=imin,10
            tc=tsvp+0.1*(i-1)*dt
            call realt(tc,pc,ps,pe,0.1*dt,r)
            call dcpg(pc,k,1,idbl3)
            call put(idbl4,k,tc)
            if(prntd) write(LOT,105) pc,tc,k
            if(k.ge.jdbl3) go to 80
            k=k+1
            ps=pc
   56 continue
c
c    save time and ray parameter of direct arrival
c
   50 continue
      pc=dcmplx(p0,0.0d+00)
      call dcpg(pc,k,1,idbl3)
      tc=t0
      call put(idbl4,k,tc)
      if(prntd) write(LOT,105) pc,tc,k
      nt0=k
      if(k.ge.jdbl3) go to 80
      k=k+1
c
c    get 10 samples from t0+dt to t0
c
      do 57 i=1,10
            p=pc
            tc=t0+i*0.1*dt
            call compt(tc,pc,p,0.1*dt,r)
            call dcpg(pc,k,1,idbl3)
            call put(idbl4,k,tc)
            if(prntd) write(LOT,105) pc,tc,k
            if(k.ge.jdbl3) go to 80
            k=k+1
   57 continue
c
c   get samples from rest of contour in imaginary p
c   Also adaptively modify tele so that array dimension is not
c   exceeded
c
c   The desired upper limit is tmax! The 11 continue ... go to 11
c   loop starts at the current tc, and increments time by
c   time(l) = tc0 + (l+1)*dt*ddtt + dt*tele*(0 + 1 + 2 ... + l)
c              = tc0 + (l+1)*dt*ddtt + dt*tele*l*(l+1)/2
c
        nleft = jdbl3 - (k + 1)
c
c   solve for lmax in
c      tmax    = tc0 + (lmax+1)*dt*ddtt + dt*tele*lmax*(lmax+1)/2
c
c   First see if the quadratic term is required assuming that tele = 0
c   If there are, then use the current value of tele. Else solve
c   for a better value of tele, called ttele
c
        xleft = (tmax - tc)/(dt * ddtt) - 1
        if(xleft.lt.0.0)xleft = 0.0
        mleft = xleft
        if(mleft .lt. nleft)then
            ttele = tele
        else
            ttele = (tmax - tc - (nleft+1)*dt*ddtt )/
     1          (0.5*dt*nleft*(nleft+1.0))
        endif
      l=0
   11 continue
           tc=tc+dt*(ddtt+ttele*l)
           p=pc
           call compt(tc,pc,p,dt,r)
           call dcpg(pc,k,1,idbl3)
           call put(idbl4,k,tc)
           if(prntd) write(LOT,105) pc,tc,k
  105   format(' ','pc=',2d15.8,5x,'tc=',d15.8,5x,'k=',i4)
           if(tc.gt.tmax) go to 80
           if(k.ge.jdbl3) go to 80
           k=k+1
           l=l+1
           go to 11
   80 continue
      ncp=k
      write(LOT,103) ncp,get(idbl4,ncp),nt0, ttele
  103 format(' ','number of contour points=',i5,' tmax=',d10.3,
     1           ' nt0=',i5,' ttele=',e10.3 )
      return
      end

        subroutine convlv(npts,dt,ksrc,r,hr,hs,tfust,mname,
     1      srcflu,refdep,recflu,vps,vss,rhs) 
c
c   convolves in 1/sqrt(t) dependence and multiplies final step 
c   function response by sqrt(2/r)/pi.  also prints, plots, 
c   and punches out final response. 
c
c   npts    I*4 - numer of time samples
c   dt  R*4 - sample interval in sec
c   ksrc    I*4 - array to indicate if Green's are computed
c   r   R*4 - epicentral distance
c   hr  R*4 - receiver depth
c   hs  R*4 - source depth
c   tfust   R*4 - time of first sample
c   mname   C*80    - name of earth model file
c   srcflu  L   - .false. source in solid
c   recflu  L   - .false. receiver in solid
c             .true.  source in fluid
c   vps R*4 - P-velocity at source
c   vss R*4 - S-velocity at source
c   rhs R*4 - density at source
c
c   variables in common blocks
c
c   iftype  I*4 File type
c           1 - single trace
c           3 - three component
c           16 - Green's function
c
c   iobsyn  I*4 1 - observed
c           2 - synthetic
c
c   itmfrq  I*4 1 - time series
c           2 - Fourier spectra (not implemented)
c
c   iunit   I*4 1   - counts 
c           2   - cm
c           3   - cm/sec
c           4   - cm/sec/sec
c           5   - m
c           6   - m/sec
c           7   - m/sec/sec
c           8   - microns
c           9   - microns/sec
c           10  - microns/sec/sec
c   junit   I*4
c           11  - Pa  (nt/m^2)
c           12  - MPa  (mega Pascal)
c
c   cfilt   C*80    comment on filtering operations
c
c   keyear  I*4 event year
c   kemon   I*4 event mon
c   keday   I*4 event day
c   kehour  I*4 event hour
c   kemin   I*4 event minute
c   esec    R*4 event second
c   evlat   R*4 event latitude
c   evlon   R*4 event longitude
c   evdep   R*4 event depth
c   
c   stname  C*8 station name
c           NOTE: AH(6), SEED(5), CSS(6), SAC(8)
c   
c   stlat   R*4 station latitude
c   stlon   R*4 station longitude
c   stelev  R*4 station elevation
c
c
c   distkm  R*4 epicentral distance in km
c   distdg  R*4 epicentral distance in degrees (111.195 km/degree)
c   evstaz  R*4 event -> station azimuth, degrees east of north
c   stevaz  R*4 station -> event azimuth, degrees east of north
c   
c   cpulse  C*80    pulse description
c   
c   ccomnt  C*80    comment
c
c   jsrc    I*4 Array of indices to indicate if traces are 
c               present
c           iftype =  1  jsrc(1) indicates if trace is 
c               present
c           iftype =  3  jsrc(i),i=1,5 indicates if trace 
c               is present
c               i = Z, 2=N, 3=E, 4=R, 5=T
c               (we could carry all five, but the real 
c               purpose  is to permit no more than 
c               3 traces, but 2 may all that are 
c               present in a real data set)
c           
c
c

        common/ihdr96/ ftype, iftype, iobsyn, itmfrq, iunit, idecon,
     1      keyear, kemon, keday, kehour, kemin,
     2      ksyear, ksmon, ksday, kshour, ksmin,
     3      jsrc, junit
        integer*4 ftype, iftype, iobsyn, itmfrq, iunit, idecon
        integer*4 keyear, kemon, keday, kehour, kemin
        integer*4 ksyear, ksmon, ksday, kshour, ksmin
        integer*4 jsrc(21), junit

        common/rhdr96/ esec, distkm, distdg, evstaz, stevaz,
     1      evlat, evlon, evdep, stlat, stlon, stelev,
     2      TP, TSV, TSH,
     3      SA, SC, SF, SN, SL, SR
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev
        REAL*4 TP, TSV, TSH
        REAL*4 SA, SC, SF, SL, SN, SR

        common/chdr96/ stname, cfilt, cpulse, ccomnt
        character*8 stname*8, cfilt*80, cpulse*80, ccomnt*80    

        character stcomp*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        real*4 ssec

        character mname*(*)
        logical srcflu, recflu
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        parameter (NSAMP=1024)
        common/resp/a(NSAMP),tfirst
        common/fftcon/data(2,2048)
        real *4 a,tfirst
        real*8 r
        dimension ksi(16),kki(16)
        common/scrat/ts2(NSAMP)
        real*4 ts2
        integer ksrc(16)
        character*8 ost(16)
        data ksi/3,3,2,2,2,1,1,1,4,4,5,5,6,6,6,4/ 
        data kki/1,2,1,2,3,1,2,3,1,2,1,2,1,2,3,4/ 
        data ost/'ZDD     ', 'RDD     ', 'ZDS     ', 'RDS     ',
     1       'TDS     ', 'ZSS     ', 'RSS     ', 'TSS     ',
     2       'ZEX     ', 'REX     ', 'ZVF     ', 'RVF     ',
     3       'ZHF     ', 'RHF     ', 'THF     ', 'PEX     '/
c
c   begin
c
        con=dsqrt(2./r)/3.14159 
        np2=2*npts 
c
        do 100 i=1,21
            if(ksrc(i).le.16)then
            if(ksrc(i).gt.0)then
                if(kki(i).eq.1)then
                    jsrc(i) = 1
                else if(kki(i).eq.2)then
                    jsrc(i) = 4
                else if(kki(i).eq.3)then
                    jsrc(i) = 5
                endif
            else
                jsrc(i) = 0
            endif
            else
                jsrc(i) = 0
            endif
  100   continue
        if(srcflu)then
            do 101 i=1,8
                jsrc(i) = 0
  101       continue
            do 102 i=11,15
                jsrc(i) = 0
  102       continue
C       10 23 99 EVENTUALLY redo jsrc from top level - everything in
c   main is ignored - also no formula for pressure in fluid
C           jsrc(16) = 1
        endif
        if(recflu)then
                jsrc(16) = 1
        endif
c
        write(LOT,104) r,dt,tfirst 
c
c   constructing 1/sqrt(t).  note that the first two points 
c   represents the singularity.  it is constructed in the same 
c   manner as the response at p0. 
c
        do 220 i=1,1024 
            ts2(i)=0. 
  220   continue
        ts2(2)=3./(sqrt(dt)*2.0) 
        do 221 i=3,npts 
            ts2(i)=1.0/sqrt(dt*(i-2)) 
  221   continue
c
c   OUTPUT TO THE FILE genray96.grn
c
c   First set header variables
c
        iftype = 16
        iobsyn = 2
        itmfrq = 1
        iunit = 1
        junit = 11
        cfilt = 'None'
        keyear = 0
        kemon = 0
        keday = 0
        kehour = 0
        kemin = 0
        esec = 0.0
        evlat = -12345
        evlon = -12345
        evdep = hs-refdep

        stname = 'GRN16'
        stlat  = -12345
        stlon = -12345
        stelev = hr-refdep
        distkm = r
        distdg = r/111.195
        evstaz = 0.0
        stevaz = 180.0
        cpulse = 'genray96'
        lmnm = lgstr(mname)
        ccomnt = mname(1:lmnm)

c
c   TP TSV and TSH have no meaning in ordeinary sense
c
        TP  = -12345.
        TSV = -12345.
        TSH = -12345.
        SA = rhs*vps*vps
        SC = rhs*vps*vps
        SL = rhs*vss*vss
        SN = rhs*vss*vss
        SF = rhs*(vps*vps - 2.0*vss*vss)
        SR = rhs
 
c
c   Output header
c
        call wrhd96(08,nerr)
c
c   Output the traces
c
        do 99 ii=1,16 
            if(jsrc(ii).gt.0)then
c
c   kk = 1, 2, 3 -> Z R T Green's functions
c   ks = mapping from Langston Green's function to Herrmann
c
            ks=ksi(ii) 
            kk=kki(ii) 
            do 222 i=1,2048 
                data(1,i)=0.0
                data(2,i)=0.0
  222       continue
            izero = 0
c
c   note a must be dimensioned 1024 since open(2 in genray requires this
c
            if(izero.ne.1)then
                jrec= (kk-1)*6 + ks
                call gramd(a,jrec)
                do 19 i=1,npts
                    data(1,i)=a(i)
c               data(2,i)=0.0
                    data(2,i)=ts2(i)
   19           continue
                call fcon(np2,dt)
            endif
c
c   here again a double length fft is used to avoid wrap around
c   noise
c
            do 25 i=1,npts 
                if(izero.eq.0)then
                    a(i)=data(1,i)*con
                else if(izero.eq.1)then
                    a(i) = 0.0
                endif
   25       continue
            if(ksrc(ii).eq.1) then
                if(kk.eq.1) then
                    write(LOT,2101) ost(ii) 
                    cmpinc = -90.0
                    cmpaz  =   0.0
                else if(kk.eq.2) then
                    write(LOT,2102) ost(ii) 
                    cmpinc = 0.0
                    cmpaz  =   0.0
                else if(kk.eq.3) then
                    write(LOT,2103) ost(ii) 
                    cmpinc = 0.0
                    cmpaz  =  90.0
                else if(kk.eq.4) then
                    write(LOT,2104) ost(ii) 
                    cmpinc = 0.0
                    cmpaz  =  90.0
                endif
                stcomp = ost(ii)
                cmpdt = dt
                ksyear = 0
                ksmon = 0
                ksday = 0
                kshour = 0
                ksmin = 0
                ssec = tfust
                call wrtr96(8,stcomp,cmpinc,cmpaz,
     1              cmpdt, npts, ksyear, ksmon, 
     2              ksday, kshour, ksmin, ssec, 
     3              a,nerr,NSAMP)

            endif
            endif
   99   continue 
  104   format(' ','range=',f10.3,5x,'dt=',e10.3,5x,
     1      'first arrival=',e10.3)
 2101   format(' ','final response - vertical   component for source ',
     1          a8) 
 2102   format(' ','final response - radial     component for source ',
     1          a8) 
 2103   format(' ','final response - tangential component for source ',
     1          a8) 
 2104   format(' ','final response - pressure   component for source ',
     1          a8) 
c
c   FORTRAN WRITE TO OPEN FILE
c   NO CARRIAGE CONTROL FOR UNIX OR MSDOS
c
        return
        end 

      function cr(p,v) 
      double complex cr,p 
      real*8 rsq,pr,pi,a,b,d,phi,e,f,v,t1 
           t1=1.0d-08 
           rsq=1.0/(v*v) 
           pr=p 
           pi=dimag(p) + 0.00001
           a=rsq -pr*pr + pi*pi 
           b=-2.*pi*pr 
           d=dsqrt(dsqrt(a*a + b*b)) 
           if(dabs(pi).ge.t1)then
                phi=datan2(b,a) 
           else
                phi=0. 
                if(a.lt.0.) phi=3.141592654 
           endif
           e=dcos(phi/2.) 
           f=dsin(phi/2.) 
           if(f.gt.t1 .or. e.le.0.0)then
                e = -e
                f = -f
           endif
           a=d*e 
           b=d*f 
           cr=dcmplx(a,b) 
      return 
      end 

      function dstdps(p0) 
      double complex p,cr 
        integer NLAY
        parameter (NLAY=200)
        common/isomod/th(NLAY),vp(NLAY),vs(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/travel/alp(NLAY),als(NLAY),nd,nnup 
      real*8 dstdps,p0,a,ea,eb,b,c
        REAL*8 DREAL
c
c       calculates second derivative of t with respect to 
c       p for first motion approximation. 
c
           a= 0.0d+00
           p=p0 
           do 10 i=1,nd 
                b=0.0d+00 
                c=0.0d+00 
                ea=0.0d+00 
                eb=0.0d+00 
                if(alp(i).ne.0.0d+00) then
                     ea=dreal(cr(p,dble(vp(i))) )
                     b=-th(i)*alp(i)/(ea*ea*ea*vp(i)*vp(i)) 
                endif
                if(als(i).ne.0.0d+00) then
                     eb=dreal(cr(p,dble(vs(i)))) 
                     c=-th(i)*als(i)/(eb*eb*eb*vs(i)*vs(i)) 
                endif
                a=a+b+c 
   10      continue 
           dstdps=a 
      return 
      end 

      function dtdp(p,r) 
      double complex dtdp,p,a,b,c,cr,ea,eb 
        integer NLAY
        parameter (NLAY=200)
        common/isomod/th(NLAY),vp(NLAY),vs(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/travel/alp(NLAY),als(NLAY),nd,nup 
        real*8 r
           a=dcmplx(0.0d+00,0.0d+00)
           do 10 i=1,nd 
                b=dcmplx(0.0d+00,0.0d+00)
                c=dcmplx(0.0d+00,0.0d+00)
                ea=dcmplx(0.0d+00,0.0d+00)
                eb=dcmplx(0.0d+00,0.0d+00)
                if(alp(i).ne.0.0d+00)then
                     ea=cr(p,dble(vp(i)) )
                     b=dcmplx(th(i)*alp(i))/ea 
                endif
                if(als(i).ne.0.0d+00) then
                     eb=cr(p,dble(vs(i)) )
                     c=dcmplx(th(i)*als(i))/eb 
                endif
                a=a+b+c 
   10      continue 
           dtdp=r-p*a 
      return 
      end 

      subroutine fcon(n,dt)
      common/fftcon/data(2,2048)
      common/virt/idbl1,idbl2,idbl3,idbl4,idbl5,idbl6,idbl7
      complex sum
      complex acmp,bcmp
      real*8 valr,vali
      double complex dcv
c
c    this subroutine convolves two time series
c    x and y where data(1,i)=x(i) and
c    data(2,i)=y(i)
c    time-frequency domain
c
c    save y - time series
c
c
c   use trick of time domain damping to reduce periodicity problems
c
        alpha = 5.0/(n*dt)
        fac = 1.0
        dfac = exp(-alpha*dt)
        n21 = n/2 + 1
        do 100 i=1,n
            data(1,i) = data(1,i)*fac
            data(2,i) = data(2,i)*fac
            valr= 0.0d+00
            vali= data(2,i)
            dcv = dcmplx(valr,vali)
            call dcpg(dcv,i,1,idbl5)
            data(2,i)=0.0
            fac = fac * dfac
  100   continue
        call four(data,n,-1,dt,df)
c
c   save complex data values
c   retrieve real y values - set up for fft
c
        do 105 i=1,n
            valr = data(1,i)
            vali = data(2,i)
            call dcpg(dcv,i,0,idbl5)
            data(2,i)=0.0
            data(1,i)= dimag(dcv)
            dcv = dcmplx(valr,vali)
            call dcpg(dcv,i,1,idbl5)
  105   continue
        call four(data,n,-1,dt,df)
c
c   convolve
c
        n21=n/2 + 1
        do 200 i=1,n21
            n2i=n+2-i
            acmp = cmplx(data(1,i),data(2,i))
            call dcpg(dcv,i,0,idbl5)
            bcmp = dcv
            sum = acmp*bcmp
            data(1,i)=real(sum)
            data(2,i)=aimag(sum)
            if(i.gt.1)then
                data(1,n2i)=real(sum)
                data(2,n2i)=-aimag(sum)
            endif
  200   continue
        data(2,n21)=0.0
        call four(data,n,+1,dt,df)
c
c   unattenuate the time series
c
        fac = 1.0
        dfac = 1.0/dfac
        do 300 i=1,n
            data(1,i) = data(1,i)*fac
            fac = fac * dfac
  300   continue
      return
      end


      subroutine four(data,nn,isign,dt,df)
c
c    the cooley-tookey fast fourier transform in usasi basic fortran
c    transform(j) = sum(data(i)*w**((i-1)(j-1)), where i and j run
c    from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  data is a one-
c    dimensional complex array (i.e., the real and imaginary parts of
c    data are located immediately adjacent in storage, such as fortran
c    places them) whose length nn is a power of two.  isign
c    is +1 or -1, giving the sign of the transform.  transform values
c    are returned in array data, replacing the input data.  the time is
c    proportional to n*log2(n), rather than the usual n**2
c    rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c    b is the number of bits in the floating point fraction.
c
c    the program computes df from dt, dt from df and checks to see
c    if they are consistent. In addition, the transforms are multiplied
c    by dt or df to make the results dimensionally correct
c
      real*4 data(*)
            n = 2 * nn
            if(dt.eq.0.0) dt = 1./(nn*df)
            if(df.eq.0.0) df = 1./(nn*dt)
            if(dt.ne.(nn*df)) df = 1./(nn*dt)
            j = 1
            do 5 i=1,n,2
                if(i.lt.j)then
                    tempr = data(j)
                    tempi = data(j+1)
                    data(j) = data(i)
                    data(j+1)=data(i+1)
                    data(i) = tempr
                    data(i+1) = tempi
                endif
                m = n/2
    3           if(j.le.m) goto 4
                    j = j-m
                    m = m/2
                    if(m.ge.2)goto 3
    4           continue
                j=j+m
    5       continue
            mmax = 2
    6       if(mmax.ge.n) goto 10
                istep= 2 *mmax
                theta = 6.283185307/float(isign*mmax)
                sinth=sin(theta/2.)
                wstpr=-2.*sinth*sinth
                wstpi=sin(theta)
                wr=1.0
                wi=0.0
                do 9 m=1,mmax,2
                        do 8 i=m,n,istep
                            j=i+mmax
                            tempr=wr*data(j)-wi*data(j+1)
                            tempi=wr*data(j+1)+wi*data(j)
                            data(j)=data(i)-tempr
                            data(j+1)=data(i+1)-tempi
                            data(i)=data(i)+tempr
                            data(i+1) = data(i+1)+tempi
    8                   continue
                        tempr = wr
                        wr = wr*wstpr-wi*wstpi + wr
                        wi = wi*wstpr+tempr*wstpi + wi
    9           continue
                mmax = istep
            go to 6
   10       continue
            if(isign.gt.0)then
c
c   frequency to time domain
c
                do 1001 iiii = 1,n
                    data(iiii) = data(iiii) * df
 1001           continue
            else
c
c   time to frequency domain
c
                do 1003 iiii = 1,n
                    data(iiii) = data(iiii) * dt
 1003           continue
            endif
        return
        end
      
      function gencof(mr,kr,ml,kl,nup,it,p) 
c
c   mr  I*4 - wave type of incident ray
c   ml  I*4 - layer of incident ray
c   kr  I*4 - wave type of converted ray
c   kl  I*4 - layer of converted ray
c   nup I*4 -  1 upward propagating ray (e.g., toward surface)
c             -1 downward propagating ray
c   it  I*4 -  0 transmission
c              1 reflection
c              2 direct ray
c
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        integer NLAY
        parameter (NLAY=200)
        common/isomod/th(NLAY),vp(NLAY),vs(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modlly/mmax
      double complex gencof,p 
      common/rmode/lv 
      common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
      logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
      double complex tpp,tps,tsp,tss,rpp,rps,rsp,rss 
C      real*8 th,c,s,d 
      m=ml 
      if(it.ne.1) then
c
c    transmission 
c
           call tranm(p,dble(vp(ml)),dble(vs(ml)),dble(rho(ml)),
     1      dble(vp(kl)),dble(vs(kl)),dble(rho(kl)),
     1          tpp,tps,tsp,tss) 
           k=kl 
           if(nup.ne.-1) then
                tps=-tps 
                tsp=-tsp 
           endif
           if(kr.ne.mr) then
                if(kr.ne.5)then
                     gencof=tps 
                else
                     gencof=tsp 
                endif
           else
                if(kr.ne.5)then
                     gencof=tss 
                else
                     gencof=tpp 
                endif
            endif
      else
c
c    reflection 
c
           k=kl-1 
           if(nup.eq.-1) k=kl+1 
           if(k.eq.0)then
           call refft(p,dble(vp(ml)),dble(vs(ml)),dble(rho(ml)),
     1      1.0d-06,1.0d-06,1.0d-06,
     2          rpp,rps,rsp,rss) 
           else
           call refft(p,dble(vp(ml)),dble(vs(ml)),dble(rho(ml)),
     1      dble(vp(k)),dble(vs(k)),dble(rho(k)),
     2          rpp,rps,rsp,rss) 
           endif
           if(k.le.ml) then
                rps=-rps 
                rsp=-rsp 
           endif
           if(kr.ne.mr)then
                 if(kr.ne.5)then
                      gencof=rps
                 else
                      gencof=rsp
                 endif
           else
                 if(kr.ne.5)then
                      gencof=rss
                 else
                      gencof=rpp
                 endif
           endif
      endif
      if(.not.prntf) return 
        if(k.ne.0.0)then
            write(LOT,100) vp(m),vs(m),rho(m),vp(k),vs(k),rho(k) 
        else
            write(LOT,100) vp(m),vs(m),rho(m),1.0d-6,1.0d-6,1.0d-6
        endif
      write(LOT,101) tpp,tps,tsp,tss 
      write(LOT,102) rpp,rps,rsp,rss 
100     format(' ','cb=',d15.8,2x,'sb=',d15.8,2x,'db=',d15.8,
     1      2x,'ca=',d15.8,2x,'sa=',d15.8,2x,'da=',d15.8) 
101     format(' ','tpp=',2d10.3,2x,'tps=',2d10.3,2x,'tsp=',
     1      2d10.3,2x,'tss= ',2d10.3) 
102     format(' ','rpp=',2d10.3,2x,'rps=',2d10.3,2x,'rsp=',
     1      2d10.3,2x,'rss= ',2d10.3) 
      return 
      end 

        subroutine gcmdln(dogrn,dfile,ddtt,tele,nasym,isupdn,ltime)
c
c   dogrn   I*4 - 0 do ALL (default)
c             1 do only EQ + Explosion
c             2 do only Explosion plus point force
c   dfile   C*80    - name of file giving distance, number points
c   ddtt    R*4
c   tele    R*4
c   nasym   I*4 - number of asymptotic terms in modified
c               Bessel function 
c   prtn?   L   = Output control for print statements
c   isupdn  I*4 - selection of rays from source
c               1 = down, 0 = both, -1 = up
c   ltime   L   - .true. only output arrival times,
c               else synthetic also
c
        common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
        logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 

        integer*4 isupdn
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        character*80 dfile
        integer*4 dogrn, nasym
        real*4 ddtt, tele
        logical ltime

        character*80 name
        integer mnmarg
        nmarg = mnmarg()
        dogrn = 0
        dfile = ' '
        ddtt = 1.0
        tele= 0.0
        nasym = 1
        prnta = .false.
        prntb = .false.
        prntc = .false.
        prntd = .false.
        prnte = .false.
        prntf = .false.
        prntg = .false.
        prnth = .false.
        isupdn = 0
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
        ltime = .false.

        i = 0
   11   i = i + 1
        if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:4).eq.'-ALL')then
                dogrn = 0
            else if(name(1:5).eq.'-EQEX')then
                dogrn = 1
            else if(name(1:4).eq.'-EXF')then
                dogrn = 2
            else if(name(1:2).eq.'-d')then
                i = i + 1
                call mgtarg(i,dfile)
            else if(name(1:2).eq.'-n')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,i10)')nasym
                if(nasym.gt.10)nasym=10
            else if(name(1:3).eq.'-va')then
                prnta = .true.
            else if(name(1:3).eq.'-vb')then
                prntb = .true.
            else if(name(1:3).eq.'-vc')then
                prntc = .true.
            else if(name(1:3).eq.'-vd')then
                prntd = .true.
            else if(name(1:3).eq.'-ve')then
                prnte = .true.
            else if(name(1:3).eq.'-vf')then
                prntf = .true.
            else if(name(1:3).eq.'-vg')then
                prntg = .true.
            else if(name(1:3).eq.'-vh')then
                prnth = .true.
            else if(name(1:3).eq.'-SU')then
                isupdn = -1
                spup = .true.
                ssup = .true.
                ssdn = .false.
                spdn = .false.
                dosud = .true.
            else if(name(1:3).eq.'-SD')then
                isupdn =  1
                spup = .false.
                ssup = .false.
                ssdn = .true.
                spdn = .true.
                dosud = .true.
            else if(name(1:5).eq.'-SPUP')then
                spup = .true.
                dosud = .true.
                isupdn = -1
            else if(name(1:5).eq.'-SSUP')then
                isupdn = -1
                ssup = .true.
                dosud = .true.
            else if(name(1:5).eq.'-SPDN')then
                isupdn = 1
                spdn = .true.
                dosud = .true.
            else if(name(1:5).eq.'-SSDN')then
                isupdn = 1
                ssdn = .true.
                dosud = .true.
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
            else if(name(1:5).eq.'-TIME')then
                ltime = .true.
            else if(name(1:2).eq.'-?')then
                call usage()
            else if(name(1:2).eq.'-h')then
                call usage()
            endif
        goto 11
   13   continue
        return
        end

        subroutine usage()
c-
c   write out program syntax
c
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
        write(LOT,*)'USAGE: ',
     1  'genray96 [-ALL] [-EQEX] [-EXF] [-d dfile] [-n nasym]  ',
     2      '[-va] [-vb] [-vc] [-vd] [-ve] [-vf] [ -vg] [-vh]',
     3      '[-SU] [-SD] [-SPUP] [-SSUP] [-SPDN] [-SSDN] ',
C     4     '[-RU] [-RD] [-RPUP] [-RSUP] [-RPDN] [-RSDN] 
     5      '[-TIME] [-?] [-h]'
        write(LOT,*)
     1  '-ALL   (default true)    Output all Green functions'
        write(LOT,*)
     1  '-EQEX  (default false)   Output only earthquake/explosions',
     2       ' Green functions'
        write(LOT,*)
     1  '-EXF   (default false)   Output only explosions/point',
     2  ' force', ' Green functions'
        write(LOT,*)
     1  '-d dfile (required)      Name of distance file '
        write(LER,*)
     1  '   dfile contains one of more lines with following entries'
        write(LER,*)
     1  '       DIST(km) DT(sec) NPTS T0(sec) VRED(km/s)',
     2  '           first time point is T0 + DIST/VRED',
     3  '           VRED=0 means infinite velocity though'
        write(LOT,*)
     1  '-n nasym (default 1)     Number of asymptotic terms'
        write(LOT,*)
     1  'The following govern wavefield at source. The default is',
     2  ' the entire wavefield'
        write(LOT,*)
     1  ' -SU    Only upward waves from source'
        write(LOT,*)
     1  ' -SD    Only downward waves from source'
        write(LOT,*)
     1  ' -SPUP  Include upward P at source'
        write(LOT,*)
     1  ' -SSUP  Include upward S at source'
        write(LOT,*)
     1  ' -SPDN  Include downward P at source'
        write(LOT,*)
     1  ' -SSDN  Include downward S at source'
        write(LOT,*)
     1  ' -TIME  Only output travel times, no synthetic'
        write(LOT,*)
     1  '[-va .. -vh ] Flags to output internal results (def:false)'
        write(LOT,*)
     1  '-?                   Display this usage message'
        write(LOT,*)
     1  '-h                   Display this usage message'
        stop 
        end

        subroutine integ(ntgr,dt,npts,y1) 
        common/fftcon/data(2,2048)
        real*4 y1
        dimension y1(1024)
c
c   this routine integrates each given term in the asymptotic 
c   series for each ray.  uses double of points.  This
c   corresponds to the ,multiplication by 1/s**ntgr 
c        in the Laplace domain
c
        np=2*npts 
c
c   calculate integration operator in time domain 
c
        if(ntgr.le.0)return
        do 9 i=1,2048 
            data(1,i)=0.0
            data(2,i)=0.0
    9   continue
        if(ntgr.le.1) then
            do 6 i=1,npts 
                data(1,i)=y1(i)
                data(2,i)=1.0
    6       continue
        else
            fac = 1
            do 11 i=1,ntgr-1
                fac = fac*i
   11       continue
            do 10 i=1,npts 
                data(1,i)=y1(i)
                data(2,i)= (dt*(i-1))**(ntgr-1)/fac
   10       continue
        endif
c
c   convolve 
c
        call fcon(np,dt)
        do 12 i=1,npts
            y1(i)=data(1,i)
   12   continue
        return 
        end 

      subroutine interp(n,np,dt,nt0,a,y,t1) 
c
c    n     - number of contour points
c    np    - numbr of points in final array
c    dt    - final sampling interval
c    nt0   - point on contour corresponding to first motion
c    a     - first motion area
c    y     - timeseries
c    t1    - time of first sample on timeseries
c
c   this subroutine interpolates the final response before 
c   the convolution with 1/sqrt(t).  it also approximates 
c   the singularity at t0 using the results of area3. 
c
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        common/virt/idbl1,idbl2,idbl3,idbl4,idbl5,idbl6,idbl7
        real*4 y(1024)
        real*8 get
        real*8 t
        real*8 a,slp,b,x1,x2,tm,t1
        real*8 y1,y2,y3
        common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
        logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
c
c   Given area about the direct arrival, e.g., about nt0,
c   adjust the direct arrival amplitude to guarantee the
c   proper amplitude, using area2()
c
c   This can only be done when there is a direct arrival, e.g.,
c   nt0 >= 1
c
        if(nt0.gt.0)then
            if(nt0.gt.10)then
                y1 = get(idbl2,nt0-10)
            else
                y1 = 0.0
            endif
            y3 = get(idbl2,nt0+10)
            ddt = get(idbl4,nt0+10) - get(idbl4,nt0)
            call area2(y1,y2,y3,ddt,a)
            call put(idbl2,nt0,y2)
        endif
c----
c   initialize final array
c----
        do 21 i=1,1024
            y(i) = 0.0
   21   continue
c
c   set time of first sample as close as possible to
c   t1 but adjust slightly so that the sample at the first motion is 
c   interpolated
c
c   As per D. Helmberger, if the first motion point is not interpolated,
c   this gives correct first motion amplitude
c   in the sense that this time series must be convolved with
c   1/sqrt(t), however, the signal can be up to 0.5 dt in
c   error in arrival time. This may be crucial when several arrivals
c   come together almost simultaneously. The ultimate would be
c   to compute this time series so that the first motion
c   time fits exactly at a sample, convolve with 1/sqrt(t), and
c   then sample at even spacings. Right now we get a discrepancly
c   between this and full-wavenumber
c   So for better comparison with full-wavenumber, we do not
c   accept an error of up to 0.5 dt.
c
        t = get(idbl4,n)
        mfin = (t-t1)/dt + 1
        if(mfin.gt.np)mfin=np
        do 30 i=1,mfin
            tm = dt*(i-1) + t1
            if(tm.ge.get(idbl4,1))goto 24
            goto 30
   24       continue
            do 25 j=2,n
            m=j
            if(tm.le.get(idbl4,j)) go to 26
   25       continue
   26       continue
                  x1 = get(idbl2,m-1)
                  x2 = get(idbl2,m  )
                  slp=(x2-x1)/(get(idbl4,m)-get(idbl4,m-1))
                  b = -slp*get(idbl4,m) + x2
                  y(i) = slp*tm + b
   30       continue
            if(prnte)then
                  write(LOT,101)
                  write(LOT,102)dt,np
                  write(LOT,103) (y(i),i=1,np)
            endif
  101 format(' ','interpolated trace including response at t0'//)
  102 format(' ','dt=',f10.3,5x,'np=',i4)
  103 format(' ',8e15.8)
      return
      end

      subroutine model(dt,npts,nsurf,r,nasym,mrex) 
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
      common/ans/ncp,nt0 
      common/virt/idbl1,idbl2,idbl3,idbl4,idbl5,idbl6,idbl7
      double complex p,gc,psq,dpdt,dtdp,gcb 
      common/resp/a(1024),tfirst
      real*8 r,av,tfust 
      common/rmode/love 
      common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
      logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
      common/scrat/y(1024)
      real*4 y
      tfust=tfirst
c
c
c    prduct finds the product of all 
c        reflection and transmission coefficients
c    less receiver functions.  ks=4 is for isotropic source.
c
c
100     format(' ','ndisp=',i4,2x,'ks=',i4) 
101     format(' ','nas=',i4) 

      do 95 ks=1,6 
           if(ks.ge.3 .and. ks.le.5 .and.love.eq.2)go to 95 
           if((ks.eq.4).and.(mrex.le.4)) go to 95 
           do 90 n=1,ncp 
                call dcpg(p,n,0,idbl3)
                call prduct(p,gc) 
c
c    calculate basic term 
c
                gcb=dcmplx(0.0d+00,0.0d+00)
                call dcpg(gcb,n,1,idbl1)
                call radpat(ks,p,psq) 
                dpdt=dcmplx(1.0d+00,0.0d+00)
                if(n.ne.nt0) then
                     dpdt=dcmplx(1.0d+00,0.0d+00)/dtdp(p,r) 
                endif
                gcb=gc*psq*dpdt 
                call dcpg(gcb,n,1,idbl1)
90         continue 
c
c    calculate each wanted displacement(
c   ndisp   = 1, vert
c       = 2, rad
c       = 3, tang
c       = 4, pressure
c   ks  = 1, Strike Slip
c       = 2, Vertical Dip Slip
c       = 3, 45 degree Dip Slip
c       = 4, Explosion
c       = 5, Vertical Force
c       = 6, Horizontal Force
c   love    = 1, P-SV generated at source
c       = 2, SH generated at source
c
           do 99 ndisp=1,4 
                if((ks.eq.4.or.ks.eq.5).and.ndisp.eq.3)go to 99
            if(ndisp.eq.4 .and. ks.ne.4)go to 99
            if(ndisp.eq.1 .and. love.eq.2)go to 99
                if(prntf) write(LOT,100) ndisp,ks 
c
c    per asymptotic term 
c
                jrec= (ndisp-1)*6 + ks
                call gramd(a,jrec)
                do 98 nas=1,nasym 
                     if(prntf) write(LOT,101) nas 
                     call vrtang(ndisp,nsurf,ks,nas,r,av,ntgr,n98) 
                     if(prntg)write(LOT,*)' return vrtang , n98=',n98
                     if(n98.ne.0) then
c
c    interpolate ray 
c
                     call interp(ncp,npts,dt,nt0,av,y,tfust)
                     if(prntg)write(LOT,*)' return interp,npts,dt,nt0',
     1                   npts,dt,nt0
c
c    integrate if neccessary 
c
                     if(ntgr.ne.0) then
                          call integ(ntgr,dt,npts,y) 
                          if(prntg)then
                              write(LOT,*)'return integ ntgr,dt,npts',
     1                             ntgr,dt,npts
                          endif
                     endif
c
c    add response to displacement
c
                     if(prntg)then
                          write(LOT,*)'ks,ndist,jrec,love,mrex',
     1                         ks,ndisp,jrec,love,mrex
                          write(LOT,*)' npts ',npts
                     endif
                     do 26 i=1,npts
                          a(i)=a(i) + y(i)
   26                continue
                     endif
   98           continue 
                call gwamd(a,jrec)
   99      continue 
   95 continue 
      return 
      end 


        subroutine npow2(npts)
c
c   Given npts, determine the N=2**m such that N >= npts
c   return the new ntps
c
        integer*4 nsamp, npts
        nsamp = npts
        npts = 1
 1000   continue
        npts = 2*npts
        if(npts.lt.nsamp)goto 1000
        return
        end

        subroutine pnot(p0,t0,r,prntb,dt) 
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        common/rays/nh(100),nm(100),ndeg,nd 
        common/rmode/love 
        integer NLAY
        parameter (NLAY=200)
        common/isomod/th(NLAY),vp(NLAY),vs(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        real*8 v,pn,pp,a,r,p0,t0 
        common/travel/alp(NLAY),als(NLAY),ndp,nnup 
        logical prntb
        double complex t,cagcon,dtdp,p , za
        REAL*8 DREAL
c
c   Get travel time for arrivals on the real p axis
c   
c

  100   format(' ','from pnot v=',f10.3) 
  101   format(' ','p0=',d15.8,10x,'t0=',d15.8) 
  102   format(' ','k=',i4,10x,'p=',2d15.8) 

        v=0.0d+00 
        do 10 i=1,ndp 
            if(alp(i).gt.0.0d+00) v=dmax1(v,dble(vp(i)))
            if(als(i).gt.0.0d+00) v=dmax1(v,dble(vs(i)))
   10   continue 
        if(prntb) write(LOT,100) v 
        eps=1.0d-05 
        p=dcmplx(1.0d+00/v - eps,0.0d+00)
        za=dtdp(p,r) 
        a = dreal(za)
        if(a .ge. 0.0d+00)then
            p0=dreal(p) 
            t=cagcon(p,r) 
            t0=dreal(t) 
            if(prntb) write(LOT,101) p0,t0 
            return 
        else
            k=0 
            pn=dreal(p) 
            pp=0.0d+00 
   13       continue
                k=k+1 
                p=dcmplx((pn+pp)/2.0d+00 ,0.0d+00)
                if(prntb) write(LOT,102) k,p 
                za=dtdp(p,r)
                a = dreal(za)
                if((dabs(a).le.0.0005d+00*dble(dt)).or.
     1              (k.ge.60)) then
                    p0=dreal(p) 
                    t=cagcon(p,r) 
                    t0=dreal(t) 
                    if(prntb) write(LOT,101) p0,t0 
                    return 
                endif
                if(a.gt.0.0d+00) then
                    pp = dreal(p)
                else
                    pn = dreal(p)
                endif 
            go to 13 
        endif
        end 

        subroutine prduct(p,gc) 
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        double complex p,gencof,gc 
        common/coff/it(100),nup1(100) 
        common/rays/nh(100),nm(100),ndeg,nd 
        common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
        logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
c
c   computes product of generalized reflection and 
c   transmission coefficients, then multiples in 
c   appropriate receiver functions. 
c

  100   format(' ','gc=',2d15.8) 

        n=nd 
        if(n.eq.1) then
            gc = dcmplx(1.0d+00,0.0d+00)
        else
            n1=n-1 
            gc = dcmplx(1.0d+00,0.0d+00)
            do 10 i=1,n1 
                nup=nup1(i) 
                ik=it(i) 
                ml=nh(i) 
                mr=nm(i) 
                kl=nh(i+1) 
                kr=nm(i+1) 
                gc=gc*gencof(mr,kr,ml,kl,nup,ik,p) 
   10       continue 
            if(prntf) write(LOT,100) gc 
        endif
        c=ndeg
        gc=gc*dcmplx(dble(abs(c)) ,0.0d+00)
        return 
        end 

        subroutine radpat(ks,p,psq) 
c
c   this subroutine computes the coefficients for the 
c   p, sv and sh source potentials for an arbitrarily oriented 
c   point dislocation( after langston and helmberger,1975) 
c   and an isotropic source. 
c   ks = 1, vertical strike-slip term 
c      = 2, vertical dip-slip term 
c      = 3, 45 degree dip-slip term 
c      = 4, isotropic source. 
c      = 5, vertical point force source. 
c      = 6, horizontal point force source. 
c
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        double complex p,psq,cr,ea,eb,sr 
        double complex ptemp
        common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
        logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
        common/rays/nh(100),nm(100),ndeg,nd 
        common/coff/it(100),nup1(100) 
        integer NLAY
        parameter (NLAY=200)
        common/isomod/th(NLAY),vp(NLAY),vs(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modlly/mmax
        real*8 fact,eps,vsi2
C       real*8 th,vp,vs,rho 
        COMPLEX*16 CDSQRT
  100   format(' ','from radpat, ks=',i4,5x,'lis=',i4,5x,'mr=',
     1      i4,5x,'eps=',f10.3/
     2      ' ','p=',2d15.8,5x,'psq=',2d15.8,5x,'sr=',2d15.8) 

c
c   lis = layer index of source used to get source wave velocity
c   mr  = code of ray leaving the source: 3(SV), 4(SH), 5(P)
c   eps = code indicating ray going upward or downward from source
c
        lis=nh(1) 
        mr=nm(1) 
        eps=-nup1(1) 
        if(mr.eq.3) then
c
c       sv wave sources 
c
            eb=cr(p,dble(vs(lis)) )
            ptemp = p
            psq=CDSQRT(ptemp)/eb 
            if(ks.eq.1)then
                sr=dcmplx(eps,0.0d+00)*p*eb 
            else if(ks.eq.2)then
                sr=eb*eb - p*p 
            else if(ks.eq.3)then
                sr=dcmplx(-3.0*eps,0.0d+00)*p*eb 
            else if(ks.eq.4)then
                sr=dcmplx(0.0d+00,0.0d+00)
            else if(ks.eq.5)then
                sr= -p 
            else if(ks.eq.6)then
                sr=dcmplx(eps,0.0d+00)*eb 
            else
                sr=dcmplx(0.0d+00,0.0d+00)
            endif
        else if(mr.eq.4) then
c
c       sh wave sources 
c
            eb=cr(p,dble(vs(lis)) )
            ptemp = p
            psq=CDSQRT(ptemp)/eb 
            vsi2 = 1.0d+00/(vs(lis)*vs(lis))
            if(ks.eq.1)then
                sr=dcmplx(-vsi2,0.0d+00)
            else if(ks.eq.2)then
                sr=dcmplx(-eps*vsi2,0.0d+00)*eb/p
            else if(ks.eq.6)then
                sr=dcmplx(-vsi2,0.0d+00)/p
            else
                sr=dcmplx(0.0d+00,0.0d+00)
            endif
            else
c
c       p wave sources 
c
            ea=cr(p,dble(vp(lis)) )
            ptemp = p
            psq=CDSQRT(ptemp)/ea 
            if(ks.eq.1)then
                sr= p*p 
            else if(ks.eq.2)then
                sr=dcmplx(2.0d+00*eps,0.0d+00)*p*ea 
            else if(ks.eq.3)then
                sr= -p*p + dcmplx(2.0d+00,0.0d+00)*ea*ea 
            else if(ks.eq.4)then
                sr=dcmplx(1.0d+00,0.0d+00)/dble(vp(lis)*vp(lis))
            else if(ks.eq.5)then
                sr=dcmplx(eps,0.0d+00)*ea
            else if(ks.eq.6)then
                sr= p
            else
                sr=dcmplx(0.0d+00,0.0d+00)
            endif
        endif
        psq=psq*sr 
c
c   put in normalization for moment of 1.0e+20 dyne-cm
c
        fact=12.5663706*rho(lis)
        psq=psq/dcmplx(fact,0.0d+00)
        if(prnte) write(LOT,100) ks,lis,mr,eps,p,psq,sr 
        return 
        end 

        subroutine realt(tc,pc,ps,pe,dt,r)
c
c   interval halving routine to fine pc(tc) given the fact that
c   tc-cagcon(pc,r) has a zero crossing between ps and pe
c
c
        real*4 dt
        real*8 tc,ps,pe,pm,pp,tm,tp,t,r
        double complex pc,cagcon,p
        REAL*8 DREAL

        pm=ps
        pc=dcmplx(pm,0.0d+00)
        tm=dreal(cagcon(pc,r))-tc
        pp=pe
        pc=dcmplx(pp,0.0d+00)
        tp=dreal(cagcon(pc,r))-tc
   10   continue
            p=dcmplx((pp+pm)/2.0d+00,0.0d+00)
            t=dreal(cagcon(p,r))-tc
            if(dsign(1.0d+00,t).eq.dsign(1.0d+00,tm)) then
                tm=t
                pm=p
            else
                tp=t
                pp=p
            endif
        if(dabs(t).gt.dble(0.01*dt)) go to 10
        pc=p
        return
        end

        subroutine recev(mr,kr,nup,p,rcv,ntype,nsurf) 
c
c   complex receiver functions 
c   
c   mr  I*4 - Incident wave type
c             3 SV
c             4 SH
c             5 P
c   kr  I*4 - receiver layer
c   nup I*4  - ray direction (?)
c             1 Up
c             -1 Down (-ve of epsilon convention)
c   p   C*8 - double complex ray parameter
c   rcv C*8 - double complex receiver function
c   ntype   I*4 - Desired motion component
c             1 Vertical
c             2 Radial
c             3 Tangential
c             4 Pressure
c   nsurf   I*4 - 1 Evaluate receiver functions for free surface
c             0 Evaluate receiver function within a layer
c
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        double complex p,rcv,ea,eb,rp,rpb,sr,cr 
        integer NLAY
        parameter (NLAY=200)
        common/isomod/th(NLAY),vp(NLAY),vs(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modlly/mmax
C       real*8 th,vp,vs,rho 
        common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
        logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
        rcv = dcmplx(0.0d+00, 0.0d+00)
        if(nsurf.ne.0) then
c
c   receiver at surface
c
            if(mr.eq.4) then
c
c           sh waves 
c
                if(ntype.eq.2) then 
                    rcv=dcmplx(2.0d+00,0.0d+00)*p
                else if(ntype.eq.3) then
                    rcv=dcmplx(2.0d+00,0.0d+00)*p 
                endif
            else
c
c           p and sv wave cases
c
                ea=cr(p,dble(vp(1)) )
                eb=cr(p,dble(vs(1)) )
                rp=(eb*eb-p*p)**2+dcmplx(4.0d+00,0.0d+00)
     1              *p*p*ea*eb 
                rpb=rp*dcmplx(dble(vs(1)*vs(1)),0.0d+00) 
                if(mr.eq.3) then
c
c               sv waves incident
c
                    if(ntype.eq.1)then
                        sr=dcmplx(4.0d+00,0.0d+00)
     1                      *ea*eb*p 
                    else if(ntype.eq.2) then
                        sr=dcmplx(-2.0d+00,0.0d+00)
     1                      *eb*(eb*eb - p*p) 
                    else if(ntype.eq.3)then
                        sr=dcmplx(-2.0d+00,0.0d+00)
     1                      *eb*(eb*eb - p*p) 
                    endif
                    rcv=sr/rpb 
                else if(mr.eq.5)then
c
c               p waves 
c
                    if(ntype.eq.1)then
                         sr=dcmplx(2.0d+00,0.0d+00)
     1                      *ea*(eb*eb - p*p) 
                        rcv=sr/rpb 
                    else if(ntype.eq.2) then
                         sr=dcmplx(4.0d+00,0.0d+00)
     1                      *ea*eb*p 
                        rcv=sr/rpb 
                    else if(ntype.eq.3) then
                         sr=dcmplx(4.0d+00,0.0d+00)
     1                      *ea*eb *p
                        rcv=sr/rpb 
                    else if(ntype.eq.4 .and. 
     1                  vs(kr).le.0.00001 ) then
                        rcv= dcmplx(0.0d+00,0.0d+00)
                    endif
                endif
            endif
        else 
c
c   buried receiver 
c
            if(mr.eq.3) then
c
c           sv waves incident
c
                if(ntype.eq.1) then
                    rcv=p 
                else if(ntype.eq.2) then
                    rcv=cr(p,dble(vs(kr)))*
     1                  dcmplx(dble(-nup),0.0d+00) 
                else if(ntype.eq.3) then
                    rcv=cr(p,dble(vs(kr)))*
     1                  dcmplx(dble(-nup),0.0d+00) 
                endif
             else if(mr.eq.4) then
c
c           sh waves incident
c
                if(ntype.eq.2) then
                     rcv=p
                else if(ntype.eq.3) then
                     rcv=p 
                endif
             else if(mr.eq.5) then
c
c           p waves 
c
                if(ntype.eq.1) then
                     rcv=cr(p,dble(vp(kr)))*nup 
                else if(ntype.eq.2) then
                     rcv= p 
                else if(ntype.eq.3) then
                     rcv= p 
                else if(ntype.eq.4 .and. vs(kr).le.0.00001 ) then
                     rcv= -dble(rho(kr))
                endif
            endif
        endif
        if(prnth) write(LOT,100) rcv 
  100   format(' ','receiver function=',2d15.8) 
        return 
        end 

        subroutine refft(p,v3,s1,d1,v4,s2,d2,rpp,rps,rsp,rss) 
c
c   reflection coefficients from medium 1 (v3, s1, d1) to
c       medium 2 (v4, s2, d2)
c
        common/rmode/love 
        real*8 k1,k2,k3,k4,v1,v2,v3,v4,s1,s2,d1,d2 ,b1,b2
        double complex e1,e2,e1p,e2p,c1,c2,c3,c4,c5,c6,ap,bp,cr,bt 
        double complex rpp,aps,bps,rps,a,b,rss,asp,bsp,rsp,p 
c
c   love    1   P-SV at source
c       2   SH at source
c
        if(v3.eq.v4 .and. s1.eq.s2 .and. d1.eq.d2)then
            rpp = dcmplx(0.0d+00, 0.0d+00)
            rsp = dcmplx(0.0d+00, 0.0d+00)
            rps = dcmplx(0.0d+00, 0.0d+00)
            rss = dcmplx(0.0d+00, 0.0d+00)
        else
            e1p=cr(p,s1) 
            e2p=cr(p,s2) 
            if(love.ne.2)then
                v1=v3 
                v2=v4 
                k4=s2*s2*d2/(s1*s1*d1) 
                b1=0.5d+00/(1.0d+00-k4) 
                b2=0.5d+00*k4/(k4-1.0d+00) 
                k1=b1/(s1*s1) 
                k2=b2/(s2*s2 )
                k3=k1+k2 
                e1=cr(p,v1) 
                e2=cr(p,v2) 
                c1=(p*p)*(dcmplx(k3,0.0d+00)-p*p)**2 
                c2=p*p*e1*e1p*e2p 
                c3=(e1*e1p)*(dcmplx(k2,0.0d+00)-p*p)**2 
                c4=e2p*(dcmplx(k1,0.0d+00)-p*p)**2 
                c5=dcmplx(k1*k2,0.0d+00)*e1*e2p 
                c6=dcmplx(k1*k2,0.0d+00)*e1p 
                ap=c1+c3-c5 
                bp=c2+c4-c6 
                a=-c1+c3-c5 
                b=-c2+c4-c6 
                bt=ap+bp*e2 
                rpp=(a-b*e2)/bt 
                aps=dcmplx(2.0d+00,0.0d+00)*p*e1
     1              *(dcmplx(k2,0.0d+00)-p*p)
     2              *(dcmplx(k3,0.0d+00)-p*p) 
                bps=dcmplx(2.0d+00,0.0d+00)*p*e1
     1              *(dcmplx(k1,0.0d+00)-p*p)*e2p 
                rps=(aps-bps*e2)/bt 
                a=-c1 +c3 +c5 
                b=-c2 +c4 +c6 
                rss=(a-b*e2)/bt 
                asp=dcmplx(2.0d+00,0.0d+00)*p*e1p
     1              *(dcmplx(k2,0.0d+00)-p*p)
     2              *(dcmplx(k3,0.0d+00)-p*p) 
                bsp=dcmplx(2.0d+00,0.0d+00)*p*e1p
     1              *(dcmplx(k1,0.0d+00)-p*p)*e2p 
                rsp=-(asp-bsp*e2)/bt 
            else
                r1=d1*s1*s1 
                r2=d2*s2*s2 
                rss=(r1*e1p-r2*e2p) / (r1*e1p+r2*e2p) 
                rpp=dcmplx(0.0d+00,0.0d+00)
                rps=dcmplx(0.0d+00,0.0d+00)
                rsp=dcmplx(0.0d+00,0.0d+00)
            endif
        endif
        return 
        end 

        subroutine tranm(p,v3,s1,d1,v4,s2,d2,tpp,tps,tsp,tss) 
c
c   transmission coefficients from medium 1 (v3, s1, d1) to
c       medium 2 (v4, s2, d2)
c
        common/rmode/love 
        real*8 k1,k2,k3,k4,v1,v2,v3,v4,s1,s2,d1,d2,rho1,rho2 
        double complex cr,e1,e2,e1p,e2p,c1,c2,c3,c4,c5,c6,ap,bp,bs,t 
        double complex tpp,tps,ts,tsp,tss,p 
c
c   love    1   P-SV at source
c       2   SH at source
c
        if(v3.eq.v4 .and. s1.eq.s2 .and. d1.eq.d2)then
            tpp = dcmplx(1.0d+00, 0.0d+00)
            tsp = dcmplx(0.0d+00, 0.0d+00)
            tps = dcmplx(0.0d+00, 0.0d+00)
            tss = dcmplx(1.0d+00, 0.0d+00)
        else
            e1p=cr(p,s1) 
            e2p=cr(p,s2) 
            if(love.ne.2) then
                v1=v3 
                v2=v4 
                rho1=d1 
                rho2=d2 
                k4=rho2*s2*s2/(rho1*s1*s1) 
                b1=0.5/(1.0d+00-k4) 
                b2=0.5*k4/(k4-1.0d+00) 
                k1=b1/(s1*s1) 
                k2=b2/(s2*s2)
                k3=k1+k2 
                e1=cr(p,v1) 
                e2=cr(p,v2) 
                c1=(p*p)*(dcmplx(k3,0.0d+00)-p*p)**2 
                c2=p*p*e1*e1p*e2p 
                c3=(e1*e1p)*(dcmplx(k2,0.0d+00)-p*p)**2 
                c4=e2p*(dcmplx(k1,0.0d+00)-p*p)**2 
                c5=dcmplx(k1*k2,0.0d+00)*e1*e2p 
                c6=dcmplx(k1*k2,0.0d+00)*e1p 
                ap=c1+c3-c5 
                bp=c2+c4-c6 
                bs=ap+e2*bp 
                t=dcmplx(2.0*k1,0.0d+00)*e1
     1              *(e2p*(dcmplx(k1,0.0d+00)-p*p)
     2              -e1p*(dcmplx(k2,0.0d+00)-p*p)) 
                tpp=t/bs 
                t=dcmplx(2.0*k1,0.0d+00)*p*e1
     1              * (e1p*e2-(dcmplx(k3,0.0d+00)-p*p)) 
                tps=t/bs 
                ts=dcmplx(2.0*k1,0.0d+00)*p*e1p
     1              *((dcmplx(k3,0.0d+00)-p*p)-e1*e2p) 
                tsp=ts/bs 
                ts=dcmplx(-2.0*k1,0.0d+00)*e1p
     1              *(e1*(dcmplx(k2,0.0d+00)-p*p)
     2              -e2*(dcmplx(k1,0.0d+00)-p*p)) 
                tss=ts/bs 
            else
                r1=d1*s1*s1 
                r2=d2*s2*s2 
                tss=(dcmplx(2.0d+00*r1,0.0d+00)*e1p)/ 
     1              (r1*e1p+r2*e2p) 
                tpp=dcmplx(0.0d+00,0.0d+00)
                tps=dcmplx(0.0d+00,0.0d+00)
                tsp=dcmplx(0.0d+00,0.0d+00)
            endif
        endif
        return 
        end 

        subroutine trav(hs,hr,prnta,isupdn,ldoray) 
c
c   hs  R*4 - source depth
c   hr  R*4 - receiver depth
c   prnta   L   - .true. print out details of this routine
c   isupdn  I*4 - selection of rays from source
c   ldoray  L   - .true. use this ray, else exclude it
c               This selection is used when only
c               part of the wavefield at source and receiver
c               are desired
c
        logical prnta
        integer*4 isupdn
        logical ldoray
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        common/rays/nh(100),nm(100),ndeg,nd 
        integer NLAY
        parameter (NLAY=200)
        common/isomod/th(NLAY),vp(NLAY),vs(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modlly/mmax
        common/travel/alp(NLAY),als(NLAY),ndeep,nup 
        common/coff/it(100),nup1(100) 
        common/rmode/love 
c
c   control for selecting up/down going rays at the receiver and source
c
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        ldoray = .true.
c
c   specify whether the problem os SH or P-SV
c
        if(nm(1).eq.4) then
            love=2 
        else
            love=1 
        endif
        n=nd
        do 10 i=1,NLAY 
            alp(i)=0. 
            als(i)=0. 
   10   continue
c
c   set up counters of how many times P or S are encountered in a layer
        do 11 i=1,n 
            if(nm(i).eq.5) alp(nh(i))=alp(nh(i))+1 
            if((nm(i).eq.3).or.(nm(i).eq.4)) als(nh(i))
     1          =als(nh(i)) + 1 
   11    continue 
c
c   determining ray direction from source.  
c
c   nup = 1, up going; 
c       =-1, down going.   
c   lis = source layer 
c   lir = receiver  layer 
c   if ndeg is lt 0, then the ray is up going from source. 
c   This is to remove the ambiguity when source and receiver 
c   are in the same layer. 
c
        lis=nh(1) 
        lir=nh(n) 
        nl=1 
c
c   determine how many times the ray goes through the source layer
c
        do 12 i=1,n 
            if(nh(i).eq.lis) nl=nl+1 
   12   continue
c
c   
c
        nup=(-1)**nl 
c
c   if receiver layer beneath the source, force ray to go down
c
        if(lir.gt.lis) nup=-nup 
c
c   if source and receiver in the same layer, use the ndeg flag
c
        if(lir.eq.lis) then
            if(ndeg.lt.0) then
                nup=+1 
            else
                nup = -1
            endif
        endif
c
c   for a direct ray, and receiver beneath the source
c       nup = -1
c   else 
c       nup =  1
c
        if( n.eq.1 ) then
            if(hr.ge.hs) then
                nup=-1 
            else
                nup= 1
            endif
        endif
        if(prnta) write(LOT,100) lis,nup 
  100   format(' ','lis=',i5,5x,'nup=',i5) 
c
c   finding interaction types at each interface and direction 
c   of each ray segment. 
c   nup1    =  1, upgoing ray 
c       = -1, downgoing ray 
c   it  =  0, transmission; 
c       =  1, reflection; 
c       =  2, direct ray. 
c
        n1=n-1 
        nup1(1)=nup 
        if(n.ne.1)then
            do 19 i=1,n1 
                k=nh(i) 
                m=nh(i+1) 
                if(m.eq.k) it(i)=1 
                if(m.ne.k) it(i)=0 
                if((nup1(i).eq.+1).and.(it(i).eq.1)) nup1(i+1)=-1 
                if((nup1(i).eq.-1).and.(it(i).eq.1)) nup1(i+1)=+1 
                if((nup1(i).eq.+1).and.(it(i).eq.0)) nup1(i+1)=+1 
                if((nup1(i).eq.-1).and.(it(i).eq.0)) nup1(i+1)=-1 
   19       continue 
        else if(n.eq.1)then
            it(1) = 2
        endif
        if(prnta) then
            write(LOT,105) (it(i),i=1,n) 
            write(LOT,106) (nup1(i),i=1,n) 
        endif
  105   format(' ','it',30i2/(3x,30i2)) 
  106   format(' ','nup1',30i2/(5x,30i2)) 
c
c   find receiver position 
c
        lir1=lir-1 
        thtot=0. 
        do 22 i=1,lir1 
            thtot=th(i) + thtot 
   22   continue
        hrl=hr - thtot 
        if(prnta) write(LOT,107) hrl 
  107   format(' ','hrl=',f10.3) 
c
c      adjust multiplier for receiver layer to include depth 
c
        a1=hrl/th(lir) 
        a2=(th(lir)-hrl)/th(lir) 
        nupa=nup1(n) 
        if(nm(n).eq.5) then
            if(nupa.eq.1) alp(lir)=alp(lir) - a1 
            if(nupa.eq.-1) alp(lir)=alp(lir) - a2 
        else if((nm(n).eq.3).or.(nm(n).eq.4)) then
            if(nupa.eq.1) als(lir)=als(lir) - a1 
            if(nupa.eq.-1) als(lir)=als(lir) - a2 
        endif
c
c   find source position 
c
        lis1=lis-1 
        thtot=0. 
        do 13 i=1,lis1 
            thtot=th(i)+thtot 
   13   continue
        hsl=hs-thtot 
        if(prnta) write(LOT,101) hsl 
  101   format(' ','hsl=',f10.3) 
c
c         adjust multiplier for source layer to include depth. 
c
        a1=hsl/th(lis) 
        a2=(th(lis)-hsl)/th(lis) 
        if(nm(1).eq.5) then
            if(nup.eq.1) alp(lis)=alp(lis)-a2 
            if(nup.eq.-1) alp(lis)=alp(lis)-a1 
        else if((nm(1).eq.3).or.(nm(1).eq.4)) then
            if(nup.eq.1) als(lis)=als(lis)-a2 
            if(nup.eq.-1) als(lis)=als(lis)-a1 
        endif
c
c   find deepest layer that the ray penetrates 
c
        ndeep=0 
        do 17 i=1,n 
            ndeep=max0(ndeep,nh(i)) 
   17   continue
        if(prnta) then
            write(LOT,102) ndeep 
            write(LOT,103) (alp(i),i=1,ndeep) 
            write(LOT,104) (als(i),i=1,ndeep) 
        endif
  102   format(' ','ndeep=',i4) 
  103   format(' ','alp',5x,10f10.3) 
  104   format(' ','als',5x,10f10.3) 
c
c   Important. For the very special case of a single ray between
c   source and receiver, which are at the same depth, the
c   alp and als will be 0.0. This will cause trouble later.
c   Since alp and als are multipliers fo layer thickness, and thus
c   are dimensionless, augment the values by a small epsilon.
c   Otherwise we get divide by zero somewhere:
c
        if(als(nh(1)).eq.0.0 .and. nm(1).ne.5)then
            als(nh(1)) = 1.0d-05
        endif
        if(alp(nh(1)).eq.0.0 .and. nm(1).eq.5)then
            alp(nh(1)) = 1.0d-05
        endif
c
c   we now know everything about this ray. Now if the
c   upgoing/downgoing ray selection is made, only pass those
c   arrays satisfying these requirements
c
        
        ldoray = .true.
        if(dosud)then
c
c       now check for the wavetype
c
            if(isupdn.eq.1 )then
                if( (ssdn.and.nm(1).ne.5) .or.
     1              (spdn.and.nm(1).eq.5)) then
                    ldoray = .true.
                else
                    ldoray = .false.
                endif
                if(nup .eq.  1)ldoray = .false.
            else if(isupdn.eq. -1)then
                if( (ssup.and.nm(1).ne.5) .or.
     1              (spup.and.nm(1).eq.5)) then
                    ldoray = .true.
                else
                    ldoray = .false.
                endif
                if(nup .eq. -1)ldoray = .false.
            endif
        endif
        return 
        end 

        subroutine ttime(p0,p1,t1,p2,t2,r,prntc) 
c
c   p0  R*8 - ray parameter of direct path
c   p1  R*8 - ray parameter of head wave
c   t1  R*8 - travel time of head wave
c   r   R*8 - distance
c   prntc   L   - .true. output debug information
c
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        integer NLAY
        parameter (NLAY=200)
        common/isomod/th(NLAY),vp(NLAY),vs(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modlly/mmax
        common/coff/it(100),nup1(100) 
        common/rays/nh(100),nm(100),ndeg,nd 
        real*8 p0,p1,t1,p2,t2,r,va,vb 
        double complex p,t,cagcon 
        logical prntc
        REAL*8 DREAL
c
c   For each ray segment,
c   nup1         1 ray is going up in the layer
c           -1 ray is going down in the layer
c   it       0 interface is transmission
c            1 reflection of ray
c            2 direct ray
c

  100   format(' ','p=',2d15.8,5x,'t=',2d15.8) 
  101   format(' ','p1=',d15.8,10x,'t1=',d15.8) 

c
c   determine the smallest possible value of ray parameter on the
c   real p-axis. This means considering the possible branch points
c
        n=nd
        p1=p0 
        do 99 i=1,n 
            nup=nup1(i) 
            vb=vs(nh(i)) 
            va=vb 
            if(nm(1).ne.4) va=vp(nh(i)) 
            p1=dmin1(p1,1.0/va,1.0/vb) 
            if(i.ne.n .and. it(i).ne.0)then
                if(nup.ne.1) then
                    vb=vs(nh(i) + 1) 
                    va=vb 
                    if(nm(1).ne.4) va=vp(nh(i) + 1) 
                else
                    if(nh(i) .ne.1)then
                    vb=vs(nh(i) - 1) 
                    va=vb 
                    if(nm(1).ne.4) va=vp(nh(i) - 1) 
                    endif
                endif
                p1=dmin1(p1,1.0/va,1.0/vb) 
            endif
   99   continue 
        p=dcmplx(p1,0.0d+00)
        t=cagcon(p,r) 
        t1=dreal(t)
        if(prntc) write(LOT,100) p,t 
        if(prntc) write(LOT,101) p1,t1 
        return 
        end 

        subroutine vrtang(ndisp,nsurf,ks,nas,r,av,ntgr,n98) 
c
c   dt  R*4 - sampling interval
c   ndisp   I*4 - 1 vertical
c             2 radial
c             3 tangential
c             4 pressure
c   nsurf   I*4 - 1 Evaluate receiver functions for free surface
c             0 Evaluate receiver function within a layer
c   ks  I*4 - Source type
c             1 Vertical strike slip
c             2 Vertical dip slip
c             3 45 degree dip slip
c             4 Explosion
c             5 Vertical Force
c             6 Horizontal Force
c   nas I*4 - number of asymptotic terms
c   r   R*8 - epicentral distance
c   av
c   ntgr    I*4 - number of times to integrate
c   n98 I*4 - 0 do not compute contribution
c
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        common/ans/ncp,nt0 
        common/rays/nh(100),nm(100),ndeg,nd 
        common/virt/idbl1,idbl2,idbl3,idbl4,idbl5,idbl6,idbl7
        common/lprint/prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
        logical prnta,prntb,prntc,prntd,prnte,prntf,prntg,prnth 
        double complex gcb,pc,p,gc,rcv,pbsl 
        real*8 dr,dsq,dstdps,r ,av,p0
        common/coff/it(100),nup1(100) 
        REAL*8 DREAL

  100   format(' ','vrtang--ntype=',i4,2x,'nterm=',i4,2x,'nord=',i4,2x,
     1      'gam=',f10.3) 
  102   format(' ','av=',d10.3,2x,'ntgr=',i4) 
  103   format(' ','from vrtang---ntype=',i4,2x,
     1      'gc=',2d12.5,' dsq=',d12.5) 
  105   format(' ','*****error from vrtang, ks=4 but ndisp=3 *****') 

        n98=0 
        ntgr=nas - 1 
        nss=nas 
        m=nd
        ms=nm(1) 
        nup=nup1(m) 
        ml=nm(m) 
        kl=nh(m) 
c
c   nterm   I*4 - parameter to control the appearance of the
c             spr term in the Hankel functions.
c           
c             The integral is of the form:
c             = 1   Kn(spr)dp
c             = 2   { Kn-1(spr) + n Kn(spr)/spr dp
c             = 3   {             n Kn(spr)/spr dp
c

        if(ks.eq.4)then
c
c       isotropic source displacements 
c
            if(ndisp.eq.1)then
c
c           vertical response 
c
                ntype=1 
                nterm=1 
                nord=0 
            else if(ndisp.eq.2)then
c
c           radial response 
c
                ntype=2 
                nterm=1 
                nord=1 
            else if(ndisp.eq.3)then
c
c           error check 
c
                write(LOT,105) 
                n98=1 
                return 
            else if(ndisp.eq.4)then
c
c   pressure field due to explosion
c
                ntype=4 
                nterm=1 
                nord=0 
            endif
        else if(ks.ne.4)then
            if(ndisp.eq.1)then
c
c           vertical response 
c
                ntype=1 
                nterm=1 
                if(ks.le.3)then
                    nord=3-ks 
                else if(ks.eq.5)then
                    nord = 0
                else if(ks.eq.6)then
                    nord = 1
                endif
            else if(ndisp.eq.2)then
c
c           radial response 
c
                ntype=2 
                if(ms.eq.4)then
c
c               SH contribution to radial - near field term
c

                    if(ks.le.3)then
                        nord=3-ks 
                    else if(ks.eq.5)then
                        nord = 0
                    else if(ks.eq.6)then
                        nord = 1
                    endif
                    if(nord.eq.0 .or. nas.eq.1) return 
                    nss=nas 
                    nterm=3 
                else
c
c               P-SV contribution to radial-near/far field term
c
                    if(ks.le.3)then
                        nord=3-ks 
                    else if(ks.eq.5)then
                        nord=0
                    else if(ks.eq.6)then
                        nord=1
                    endif
                    if(ks.eq.2 .or. ks.eq.1 .or. ks.eq.6)then
                        nterm = 2 
                    else
                        nterm = 1
                    endif
                endif
            else if(ndisp.eq.3)then
c
c           tangential response
c
                ntype=3 
                if(ms.ne.4)then
c
c               P-SV contribution to tangential-near field term
c
                    if(ks.le.3)then
                        nord=3-ks 
                    else if(ks.eq.5)then
                        nord = 0
                    else if(ks.eq.6)then
                        nord = 1
                    endif
                    if(nord.eq.0 .or. nas.eq.1) return 
                    nss=nas
                    nterm=3 
                else
c
c               SH contribution to tangential - near/far field term
c
                    if(ks.le.3)then
                        nord=3-ks 
                    else if(ks.eq.5)then
                        nord=0
                    else if(ks.eq.6)then
                        nord=1
                    endif
                    if(ks.eq.2 .or. ks.eq.1 .or. ks.eq.6)then
                        nterm = 2 
                    else
                        nterm = 1
                    endif
                    if(nord.eq.0) return 
                endif
            endif
        endif
        if(prntg) write(LOT,100) ntype,nterm,nord
c
c   calculate response 
c
        if(prntg)write(LOT,*)'ncp = ',ncp
        do 99 n=1,ncp 
            call dcpg(gcb,n,0,idbl1)
            gc=gcb
            call dcpg(pc,n,0,idbl3)
            p=pc 
            call recev(ml,kl,nup,p,rcv,ntype,nsurf) 
            call bessel(r,nss,nord,nterm,p,pbsl) 
            gc=gc*rcv*pbsl
            dr=dimag(gc)
            call put(idbl2,n,dr)
   99   continue
        if(prntg)write(LOT,*)' vrtang--after 99 continue'
        if(prntg)write(LOT,*)' nt0 ',nt0
c
c   determine first motion approximation at p0.  find dstdps 
c   at p0 and assume dt separation on either side of t0. 
c
        if(nt0.ge.0) then
            call dcpg(pc,nt0,0,idbl3)
            p0=dreal(pc)
            dsq=dstdps(p0) 
            call dcpg(gcb,nt0,0,idbl1)
            gc=gcb
            call dcpg(pc,nt0,0,idbl3)
            p=pc
            call recev(ml,kl,nup,p,rcv,ntype,nsurf)
            call bessel(r,nss,nord,nterm,p,pbsl)
            gc=gc*rcv*pbsl
            gc=gc/dcmplx(dsqrt(-2.0d+00*dsq),0.0d+00)
            call area3(gc,av,nt0) 
            if(prntg)write(LOT,103) ntype,gc,dsq
            if(prntg) write(LOT,102) av,ntgr 
        endif
c
c   take imaginary part for response 
c
        n98=1 
        if(prntg)write(LOT,*)' returning from vrtang'
        return 
        end 


        subroutine dcpg(dc,i,j,int)
c
c      double put get of complex numbers
c   dc  C*8 - complex numbr
c   i   I*4 - index of array
c   j   I*4 -  0 get value
c              1 put
c   int I*4 -  reference to offset from set
c
        REAL*8 DREAL
        double complex dc
        real*8 xr,xi
        real*8 get
        k=2*i -1
        if(j.eq.0)then
            xr = get(int,k)
            xi = get(int,k+1)
            dc = dcmplx(xr,xi)
        else
            xr = dreal(dc)
            xi = dimag(dc)
            call put(int,k,xr)
            call put(int,k+1,xi)
        endif
        return
        end

c---------------------------------------------------------------------c
c             THE FOLLOWING ROUTINES ARE MACHINE DEPENDENT
c
c             THERE ARE TWO VERSIONS
c
c             1) FOR USE WITH A 32 BIT MINICOMPUTER WITH
c                VIRTUAL ADDRESSING
c
c             2) FOR USE WITH A 16 BIT MINICOMPUTER. DISK
c                IS USED AS AUXILIARY MEMORY FOR STORAGE OF
c                LARGE ARRAYS. A BUFFERED I/O SCHEME IS USED
c                FOR EFFICIENCY. DIRECT ACCESS FORTRAN I/O
c                IS REQUIRED, WHICH MAY BE DIFFICULT ON SOME
c                LARGE MAINFRAMES
c
c---------------------------------------------------------------------c



c---------------------------------------------------------------------c
c
c
c             VIRTUAL MACHINE VERSION
c                TESTED ON RIDGE32, HP9000, MASSCOMP
c
c---------------------------------------------------------------------c
        function get(i,j)
        common/BIG/bigsav(5,5000)
        real*8 get
        real*8 bigsav
            get = bigsav(i,j)
        return
        end

        subroutine put(i,j,x)
        common/BIG/bigsav(5,5000)
        real*8 x
        real*8 bigsav
            bigsav(i,j)=x
        return
        end

        function set(size)
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        integer*4 set,size
        common/BDATA/iv
            iv = iv +1
            set = iv
            write(LOT,*)'set: set=',set
        return
        end

        block data
        common/BIG/bigsav(5,5000)
        real*8 bigsav
        common/BDATA/iv
        data bigsav/25000*0.0d00/
        data iv/0/
        end

        subroutine gramd(a,jrec)
        dimension a(1024)
        common/BIGDUM/b
        real*4 b(1024,24)
        do 100 j=1,1024
            a(j)=b(j,jrec)
  100   continue
        return
        end

        subroutine gwamd(a,jrec)
        dimension a(1024)
        common/BIGDUM/b
        real*4 b(1024,24)
        do 100 j=1,1024
            b(j,jrec)=a(j)
  100   continue
        return
        end

        subroutine goamd()
        return
        end
     
        subroutine gcamd()
        return
        end
     
        subroutine setcl()
        return
        end

c----------------------------------------------------------------------c
c
c           16 BIT VERSIONS
c              TESTED ON RIDGE32, PDP11/70
c          This uses disk as virtual memory
c
c----------------------------------------------------------------------c
c
c       function get(i,j)
c       common/dv/iopen,iv,ivfrst(5),ivlast(5),ipop,irec,
c    1          buffy(128,5),buffi(128),lrec(5)
c       real*8 buffy,buffi
c       integer*4 i,j
c       real*8 get
c       jrec = (j+ipop-1)/ipop
c       krec= jrec-1 + ivfrst(i)
c       k= j -ipop*(jrec-1)
c       if(irec.ne.krec.and.lrec(i).ne.krec)then
c               do 100 ll=1,128
c  100           buffi(ll)=buffy(ll,i)
c               write(1,rec=lrec(i))buffi
c               read(1,rec=krec)buffi
c               do 101 ll=1,128
c  101           buffy(ll,i)=buffi(ll)
c       endif
c       get = buffy(k,i)
c       irec=krec
c       lrec(i)=krec
c       return
c       end
c
c
c       subroutine put(i,j,x)
cc put x(j) into array i
c       common/dv/iopen,iv,ivfrst(5),ivlast(5),ipop,irec,
c    1          buffy(128,5),buffi(128),lrec(5)
c       real*8 buffy,buffi
c       integer*4 i,j
c       real*8 x
c       jrec = (j+ipop-1)/ipop
c       krec= jrec-1 + ivfrst(i)
c       k= j -ipop*(jrec-1)
c       if(irec.ne.krec.and.lrec(i).ne.krec)then
c               do 100 ll=1,128
c  100           buffi(ll)=buffy(ll,i)
c               write(1,rec=lrec(i))buffi
c               read(1,rec=krec)buffi
c               do 101 ll=1,128
c  101           buffy(ll,i)=buffi(ll)
c       endif
c       buffy(k,i)=x
c       irec=krec
c       lrec(i)=krec
c       return
c       end
c
c       function set(size)
c       common/dv/iopen,iv,ivfrst(5),ivlast(5),ipop,irec,
c    1          buffy(128,5),buffi(128),lrec(5)
c       real*8 buffy,buffi
c       integer*4 set,size
c       if(iopen.le.0)then
c               open(1,file='junkfil',access='direct',recl=1024,
c    1                  form='unformatted',status='new')
c               iopen=1
c       endif
c       iv = iv + 1
c       if(iv.eq.1)then
c               ivfrst(iv) = 1
c       else
c               ivfrst(iv) = ivlast(iv-1) + 1
c       endif
c       itemp = size + ipop -1
c       nrecs = itemp/ipop -1
cc       ivlast(iv)=ivfrst(iv)+mrecs -1
c       ivlast(iv)=ivfrst(iv)+nrecs
cc we have to initialize the array with something
c       do 100 irec=ivfrst(iv),ivlast(iv)
c       write(1,rec=irec)buffi
c     lrec(iv)=irec
c  100   continue
c       set = iv
c       return
c       end
c
c       block data
c       common/dv/iopen,iv,ivfrst(5),ivlast(5),ipop,irec,
c    1          buffy(128,5),buffi(128),lrec(5)
c       real*8 buffy,buffi
c       data iopen/0/,iv/0/,ivfrst/0,0,0,0,0/
c       data ivlast/0,0,0,0,0/
c       data ipop/128/
c       data buffy/640*0.0/
c       data buffi/128*0.0/
c       data lrec/5*0/
c       data irec/0/
cc--this is set up to use disk as virtual storage
cc  of double precision arrays
cc  Here we chose 128 real*8 buffering.
cc  to change this change dimensions in block data.f, put.f get.f
cc  ipop and the initializations here and the
cc  recl= in open which is presently 8*128
cc  irec is a pointer to the present disk record
cc  lrec(j) is a pointer to the disk record which
cc  provided data for the buffer buffy(k,j)
cc  buffi is required only to do disk io
cc  
c       end
c
c     subroutine gramd(a,jrec)
c     dimension a(1024)
c           read(2,rec=jrec)a
c     return
c     end
c
c     subroutine gwamd(a,jrec)
c     dimension a(1024)
c           write(2,rec=jrec)a
c     return
c     end
c
c     subroutine goamd()
c           open(2,recl=4096,form='unformatted',
c    1                  access='direct',status='scratch')
c     return
c     end
c
c     subroutine gcamd()
c           close(2,status='delete')
c     return
c     end
c
c     subroutine setcl()
c           close(1,status='delete')
c     return
c     end

        subroutine chkdep(hs,hr,srcflu, recflu)
        integer NLAY
        parameter (NLAY=200)
        common/isomod/th(NLAY),vp(NLAY),vs(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modlly/mmax
        logical srcflu, recflu
c
c   Get the maximum depth of source or receiver
c   hs  R*4 - source depth in model
c   hr  R*4 - receiver depth in model
c   srcflu  L   - .false. source in solid
c             .true.  source in fluid
c   recflu  L   - .false. receiver in solid
c             .true.  receiver in fluid
c   
c
        hmax = amax1(hs,hr)
c
c   get the thickness of the layer stack
c
        dep = 0.0
        dm1 = -100.0
        srcflu = .false.
        recflu = .false.
        do 100 i=1,mmax-1
            dm1 = dep
            dep = dep + th(i)
            if(hs.ge.dm1 .and. hs.le.dep .and.
     1              vs(i).eq.0)then
                srcflu = .true.
            endif
            if(hr.ge.dm1 .and. hr.le.dep .and.
     1              vs(i).eq.0)then
                recflu = .true.
            endif
c
c   NOT BE CARFEUL IF WE EVER GO TO UPPER ATMOSPHERE
c
  100   continue
        if(hmax .gt. dep)then
            th(mmax) = hmax + 1.0
            vp(mmax+1) = vp(mmax)
            vs(mmax+1) = vs(mmax)
            rho(mmax+1) = rho(mmax)
            qa(mmax+1) = qa(mmax)
            qb(mmax+1) = qb(mmax)
            etap(mmax+1) = etap(mmax)
            etas(mmax+1) = etas(mmax)
            frefp(mmax+1) = frefp(mmax)
            frefs(mmax+1) = frefs(mmax)
            mmax = mmax + 1
        endif
        return
        end

