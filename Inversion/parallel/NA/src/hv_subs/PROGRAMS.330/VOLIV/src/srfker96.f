      program srfker96
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME IV                                                       c
c                                                                      c
c      PROGRAM: SRFKER96                                               c
c                                                                      c
c      COPYRIGHT 2012                                                  c
c      R. B. Herrmann                                                  c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c-----
c       CHANGES
c         26 JUL 2012 - added to distribution
c-----
        integer LOT
        integer NL,NLAY, NL2
        parameter(LOT=6)
        parameter (NL=200,NLAY=200,NL2=NL+NL)
c-----
c       NL  - number of layers in model
c       NL2 - number of columns in model (first NL2/2 are
c           - velocity parameters, second NL2/2 are Q values)
c-----
        integer nf10(NL)
        common/ctrl/numa,d(NL),a(NL),b(NL),r(NL),rat(NL),dd(NL2),x(NL2),
     $      h(NL2),u(NL2),ct(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        real*4 dcdb(NL2),dcda(NL2),dgdq(NL2),dudb(NL2),duda(NL2),
     1      dudqa(NL2),dudqb(NL2),dcdqa(NL2),dcdqb(NL2),
     1      dgda(NL2),dgdb(NL2),dcdh(NL2), dudh(NL2), dgdh(NL2),
     1      dgdqa(NL2),dgdqb(NL2)
        real cvel,gvel,gam

        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY),
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modtit/title
        character*80 title

        character ilorr*1, ipug*1


        call gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup,dlam,qaqb,wref,invcsl,
     2      lstinv,twnmin,twnmax,iter,nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
c-----
c       get the current model
c-----
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
        m = mmax
        do i=1,m
              if(qa(i).gt.1.0)qa(i) = 1.0/qa(i)
              if(qb(i).gt.1.0)qb(i) = 1.0/qb(i)
        enddo
c-----
c       get all the dispersion values
c-----
        open(2,file='tmpsrfi.08',form='unformatted',access='sequential')
        rewind 2
  400   continue
        read(2,end=30) itst,k1,md1,tp1
        if (itst.eq.1)then
               ilorr = 'L'
        else if (itst.eq.2)then
               ilorr = 'R'
        endif
        if(k1.eq.1)then
               ipug='C'
        else if(k1.eq.2)then
               ipug='U'
        endif
C        WRITE(6,*)ilorr,ipug,md1,tp1
c-----
c       zero everything
c-----
        call zero(dudh,1,m)
        call zero(dcdh,1,m)
        call zero(duda,1,m)
        call zero(dcdb,1,m)
        call zero(duda,1,m)
        call zero(dcdb,1,m)
        call zero(dcdqa,1,m)
        call zero(dcdqb,1,m)
        call zero(dudqa,1,m)
        call zero(dudqb,1,m)
        call zero(dgda,1,m)
        call zero(dgdb,1,m)
        call zero(dgdqa,1,m)
        call zero(dgdqb,1,m)
        call zero(dgdh,1,m)
        gvel = 0.0
        cvel = 0.0

c-----
c       THEORETICAL
c           itst    : Whether dispersion point found
c               : 0 mode does not exist, 1 = Love, 2 = Rayleigh
c           k1  : 1 Phase 2 Group
c           md1 : Mode
c           tp1 : Period
c-----
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
            if(itst.eq.2)then
                  read(2,end=30)(dcda(j),j=1,m)
             endif
        endif
c-----
c       only tabulate if the group velocity has been computed
c-----
        if(k1.ne.2)go to 400
        call elas(m,dl,cvel,dcda,dcdb,dcdh,
     1     gvel,duda,dudb,dudh,itst,md1,tp1)
c----
c   check this I meay need the transformed model velocities
c   here output the causal parameters
c-----
        call getgam(dcdb,vb,qb,dcda,va,qa,m,itst,cvel,tp1,
     1      gam)
c-----
c     modify to get the dcda duda for causal?
c     how about dcdh dudh causal?
c     After this call the partial derivatives are over written to the
c     causal values
c-----
        call cslmod(gvel,cvel,gam,wref,tp1,m,
     1      qb   ,qa   ,dcdb,dcda,dudb,duda,itst,va,vb,2,
     1      dcdh,dudh,dgdh,
     1      dcdqb,dcdqa,dudqb,dudqa,dgdb,dgda,dgdqb,dgdqa)

        

             call anelas(m,dl,cvel,dcda,dcdb,dcdh,
     1          gvel,duda,dudb,dudh,dgdh,itst,md1,tp1,
     2          gam,dcdqb,dcdqa,dudqb,dudqa,dgdb,dgda,dgdqb,dgdqa)

        go to 400
   30   continue
        write(LOT,9)
    9   format('___________________________________________',
     1    '________________________________________________')
        close(2)
        end

        subroutine elas(m,dl,cvel,dcda,dcdb,dcdh,
     1     gvel,duda,dudb,dudh,itst,md1,tp1)
c-----
c           itst    : Whether dispersion point found
c               : 0 mode does not exist, 1 = Love, 2 = Rayleigh
c           k1  : 1 Phase 2 Group
c           md1 : Mode
c           tp1 : Period
c-----
        integer LOT
        parameter(LOT=6)
        integer NL,NLAY, NL2
        parameter (NL=200,NLAY=200,NL2=NL+NL)
        real*4 dcdb(NL2),dcda(NL2),dgdq(NL2),dudb(NL2),duda(NL2),
     1      dudq(NL2),dcdq(NL2),dgdv(NL2),dcdh(NL2), dudh(NL2)
        real cvel,gvel,gam
        real dl(NL)
c-----
c       internally md1 = 1 corresdponds to the fundamental mode
c-----
        WRITE(LOT, 9)
        if(itst.eq.1)then
        WRITE(LOT,10)tp1,md1 -1,cvel,gvel
        WRITE(LOT,11)
        do i=1,m
        WRITE(LOT,12),i,dl(i),dcdb(i),dudb(i),dcdh(i),dudh(i)
        enddo
        else
        WRITE(LOT,20)tp1,md1 -1,cvel,gvel
        WRITE(LOT,21)
        do i=1,m
        WRITE(LOT,22),i,dl(i),dcda(i),dcdb(i),duda(i),
     1       dudb(i),dcdh(i),dudh(i)
        enddo
        endif
    9   format('___________________________________________',
     1    '________________________________________________')
   10   format('Elastic  Love wave:     Period=',f10.3,' Mode =',i5,
     1    ' C=',f10.3,' U=',f10.3)
   11   format('LAYER    THICK     dc/db      dU/db      dc/dh',
     1    '      dU/dh')
   12   format(i5,f10.3,4(1pe11.3))
   20   format('Elastic  Rayleigh wave: Period=',f10.3,' Mode =',i5,
     1    ' C=',f10.3,' U=',f10.3)
   21   format('LAYER    THICK     dc/da      dc/db      dU/da',
     1    '      dU/db      dC/dh      dU/dh')
   22   format(i5,f10.3,6(1pe11.3))
         
        return
        end

        subroutine anelas(m,dl,cvel,dcda,dcdb,dcdh,
     1          gvel,duda,dudb,dudh,dgdh,ilvry,md1,tp1,
     2          gam,dcdqb,dcdqa,dudqb,dudqa,dgdb,dgda,dgdqb,dgdqa)
        integer LOT
        parameter(LOT=6)
        integer NL,NLAY, NL2
        parameter (NL=200,NLAY=200,NL2=NL+NL)
        real*4 dcdb(NL2),dcda(NL2),dgdq(NL2),dudb(NL2),duda(NL2),
     1      dudqb(NL2),dcdqb(NL2),dgdb(NL2),dgdqb(NL2),
     1      dudqa(NL2),dcdqa(NL2),dgda(NL2),dgdqa(NL2),
     1      dcdh(NL2), dudh(NL2), dgdh(NL2)
        real cvel,gvel,gam
        real dl(NL)
c-----
c       internally md1 = 1 corresdponds to the fundamental mode
c-----
        WRITE(LOT, 9)
        if(ilvry.eq.1)then
        WRITE(LOT,10)tp1,md1 -1,cvel,gvel,gam
        WRITE(LOT,11)
        do i=1,m
        WRITE(LOT,12),i,dl(i),dcdb(i),dudb(i),dcdh(i),dudh(i),
     1      dcdqb(i),dudqb(i),dgdqb(i)
        enddo
        else
        WRITE(LOT,20)tp1,md1 -1,cvel,gvel,gam
        WRITE(LOT,21)
        do i=1,m
        WRITE(LOT,22),i,dl(i),dcda(i),dcdb(i),duda(i),
     1      dudb(i),dcdh(i),dudh(i),
     1      dcdqa(i),dcdqb(i),dudqa(i),dudqb(i),
     1      dgdqa(i),dgdqb(i)
        enddo
        endif
    9   format('___________________________________________',
     1    '________________________________________________')
   10   format('Anelastic Love wave:     Period=',f10.3,' Mode =',i5,
     1    ' C=',f10.3,' U=',f10.3, ' GAMMA=',1pe11.3)
   11   format('LAYER    THICK     dc/db      dU/db      dc/dh',
     1    '      dU/dh     dc/dQbi    dU/dQbi',
     1    '    dg/dQbi')
   12   format(i5,f10.3,7(1pe11.3))
   20   format('Anelastic Rayleigh wave: Period=',f10.3,' Mode =',i5,
     1    ' C=',f10.3,' U=',f10.3 ' GAMMA=',1pe11.3)
   21   format('LAYER    THICK     dc/da      dc/db      dU/da',
     1    '      dU/db      dc/dh      dU/dh',
     2    '     dc/dQai    dc/dQbi    dU/dQai',
     2    '    dU/dQbi',
     2    '    dg/dQai    dg/dQbi')
   22   format(i5,f10.3,12(1pe11.3))
        return
        end

        subroutine gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup,dlam,qaqb,wref,invcsl,
     2      lstinv,twnmin,twnmax,iter,nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        integer NL
        parameter (NL=200)
        integer nf10(NL)
c-----
c    read control file
c-----
        open(1,file='tmpsrfi.00',form='unformatted',access='sequential')
        rewind 1
        read(1) iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,nup,dlam
     1      ,qaqb,wref,invcsl,lstinv,twnmin,twnmax,iter,
     2      nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5
        close(1,status='keep')
        return
        end

        subroutine getmod(rlun,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,listmd)
c-----
c       HISTORY
c
c       09 08 2000  gave ierr an initial default value for g77
c       01 13 2001  put in close(lun) if file is not model file
c       03 MAY 2002     Modify to permit read from standard input
c       06 JUL 2005 moved inquire to permit use of STDIN
c
c-----
c       General purpose model input
c       This model specification is designed to be as 
c           general as possible
c
c       Input lines
c       Line 01: MODEL
c       Line 02: Model Name
c       Line 03: ISOTROPIC or ANISOTROPIC or 
c           TRANSVERSELY ANISOTROPIC
c       Line 04: Model Units, First character 
c           is length (k for kilometer
c           second is mass (g for gm/cc), third is time (s for time)
c       Line 05: FLAT EARTH or SPHERICAL EARTH
c       Line 06: 1-D, 2-D or 3-D
c       Line 07: CONSTANT VELOCITY
c       Line 08: open for future use
c       Line 09: open for future use
c       Line 10: open for future use
c       Line 11: open for future use
c       Lines 12-end:   These are specific to the model
c           For ISOTROPIC the entries are
c           Layer Thickness, P-velocity, S-velocity, Density, Qp, Qs,
c           Eta-P, Eta S (Eta is frequency dependence), 
c           FreqRefP, FreqRefP
c-----
cMODEL
cTEST MODEL.01
cISOTROPIC
cKGS
cFLAT EARTH
c1-D
cCONSTANT VELOCITY
cLINE08
cLINE09
cLINE10
cLINE11
c H  VP  VS   RHO   QP  QS   ETAP   ETAS REFP  REFS
c1.0    5.0 3.0 2.5 0.0 0.0 0.0 0.0 1.0 1.0
c2.0    5.1 3.1 2.6 0.0 0.0 0.0 0.0 1.0 1.0
c7.0    6.0 3.5 2.8 0.0 0.0 0.0 0.0 1.0 1.0
c10.0   6.5 3.8 2.9 0.0 0.0 0.0 0.0 1.0 1.0
c20.0   7.0 4.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0
c40.0   8.0 4.7 3.3 0.0 0.0 0.0 0.0 1.0 1.0
c-----
c-----
c       rlun    I*4 - logical unit for reading model file. This
c                 unit is released after the use of this routine
c       mname   C*(*)   - model name - if this is stdin or 
c           STDIN just read
c                 from standard input
c       mmax    I*4 - number of layers in the model, last layer is
c                    halfspace
c       title   C*(*)   - title of the model file
c       iunit   I*4 - 0 Kilometer, Gram, Sec
c       iiso    I*4 - 0 isotropic 
c                 1 transversely anisotropic 
c                 2 general anisotropic 
c       iflsph  I*4 - 0 flat earth model
c                 1 spherical earth model
c       idimen  I*4 - 1 1-D
c               - 2 2-D
c               - 3 3-D
c       icnvel  I*4 - 0 constant velocity
c                 1 variable velocity
c       ierr    I*4 - 0 model file correctly read in
c               - -1 file does not exist
c               - -2 file is not a model file
c                 -3 error in the model file
c       listmd  L   - .true. list the model
c------

        implicit none
        character mname*(*), title*(*)
        integer rlun
        integer*4 mmax, iunit, iiso, iflsph, idimen, icnvel
        integer*4 ierr
        character string*80
        logical listmd
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

        logical ext
        character ftype*80
        integer lun, j, i, irefdp

c-----
c       test to see if the file exists
c-----
        ierr = 0
c-----
c       test for input
c-----
        if(MNAME(1:5).eq.'stdin' .or. mname(1:5).eq.'STDIN')then
c-----
c           do not open anything, use standard output
c-----
            lun = LIN
        else
            lun = rlun
            inquire(file=mname,exist=ext)
            if(.not.ext)then
                ierr = -1
                write(LER,*)'Model file does not exist'
                return
            endif
c-----
c           open the file
c-----
            open(lun,file=mname,status='old',form='formatted',
     1          access='sequential')
            rewind lun
        endif
c-----
c       verify the file type
c-----
c-----
c       LINE 01
c-----
        read(lun,'(a)')ftype
        if(ftype(1:5).ne.'model' .and. ftype(1:5).ne.'MODEL')then
            ierr = -2
            write(LER,*)'Model file is not in model format'
            close(lun)
            return
        endif
c-----
c       LINE 02
c-----
        read(lun,'(a)')title
c-----
c       LINE 03
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'ISO' .or. string(1:3).eq.'iso')then
            iiso = 0
        else if(string(1:3).eq.'TRA' .or. string(1:3).eq.'tra')then
            iiso = 1
        else if(string(1:3).eq.'ANI' .or. string(1:3).eq.'ani')then
            iiso = 2
        endif
c-----
c       LINE 04
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'KGS' .or. string(1:3).eq.'kgs')then
            iunit = 0
        endif
c-----
c       LINE 05
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'FLA' .or. string(1:3).eq.'fla')then
            iflsph = 0
        else if(string(1:3).eq.'SPH' .or. string(1:3).eq.'sph')then
            iflsph = 1
        endif
c-----
c       LINE 06
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'1-d' .or. string(1:3).eq.'1-D')then
            idimen = 1
        else if(string(1:3).eq.'2-d' .or. string(1:3).eq.'2-D')then
            idimen = 2
        else if(string(1:3).eq.'3-d' .or. string(1:3).eq.'3-D')then
            idimen = 3
        endif
c-----
c       LINE 07
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'CON' .or. string(1:3).eq.'con')then
            icnvel = 0
        else if(string(1:3).eq.'VAR' .or. string(1:3).eq.'var')then
            icnvel = 1
        endif
c-----
c       get lines 8 through 11
c-----
        do 900 i=8,11
            read(lun,'(a)')string
  900   continue
c-----
c       get model specifically for 1-D flat isotropic
c-----
c-----
c       get comment line
c-----
        read(lun,'(a)')string
        mmax = 0
        refdep = 0.0
        irefdp = 0
        if(iiso.eq.0)then
 1000       continue
            j = mmax +1
                read(lun,*,err=9000,end=9000)d(j),a(j),b(j),
     1              rho(j),qa(j),qb(j),etap(j),etas(j),
     2              frefp(j),frefs(j)
                if(d(j).lt.0.0)then
                    d(j) = -d(j)
                    refdep = refdep + d(j)
                    irefdp = j
                endif
            mmax = j
            go to 1000
 9000       continue
        endif
    1   format(' LAYER             H      P-VEL     S-VEL   DENSITY  ')
    2   format(' ',i5,5x,4f10.3)
    3   format(' ','-SURFACE ','- - - - - ','- - - - - ',
     1      '- - - - - ','- - - - - -')
        if(mmax.gt.0)then
            if(listmd)then
            ierr = 0
            write(LOT,1)
            do 2000 i=1,mmax
                write(LOT,2)
     1              i,d(i),a(i),b(i),rho(i)
                if(i.eq.irefdp)write(LOT,3)
 2000       continue
            endif
        else 
            ierr = -3
            write(LER,*)'Error in model file'
        endif
        if(lun.ne.LIN)close (lun)
        return
        end

        subroutine cslmod(u,c,gam,wref,tp,m,
     1      qbinv,qainv,dcdb,dcda,dudb,duda,ilvry,a,b,k1,
     1      dcdh,dudh,dgdh,
     1     dcdqb,dcdqa,dudqb,dudqa,dgdb,dgda,dgdqb,dgdqa)
        integer NL,NLAY, NL2
        parameter (NL=200,NLAY=200,NL2=NL+NL)
        real u,c,gam,wref,tp
        integer m,k1,ilvry
        real*4 qbinv(NL2),qainv(NL2),a(NL2),b(NL2)
        real*4 dcdh(NL2), dudh(NL2), dgdh(NL2)
        real*4 dcdb(NL2),dcda(NL2),dudb(NL2),duda(NL2)
        real*4 dcdqb(NL2),dcdqa(NL2),dudqb(NL2),dudqa(NL2),dgdb(NL2),
     1         dgda(NL2),dgdqb(NL2),dgdqa(NL2)


c-----
c       modify phase, group velocities and partials for
c       constant Q causality
c
c       ilvry 1 = L, 2 = R
c           k1  : 1 Phase 2 Group
c           tp : Period
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
        faca = 0.0
        fars = 0.0
        
c-----
c       get correct partials
c-----
        do 100 i=1,m
c-----
c          save values that will be overwritten
c          dcdh is not changed
c          dudh is not changed
c-----
            dcdp = dcda(i)
            dcds = dcdb(i)
            dcdb(i) = dcds*(1. + faclog*qbinv(i) )
            dcda(i) = dcdp*(1. + faclog*qainv(i) )
            dcdqa(i) = faclog*dcdp*a(i)
            dcdqb(i) = faclog*dcds*b(i)
            dgdqa(i) = 0.5*omega*dcdp*a(i)/(c0*c0)
            dgdqb(i) = 0.5*omega*dcds*b(i)/(c0*c0)
            if(k1.eq.2)then
                 duds = dudb(i)
                 dudp = duda(i)
                 dudb(i) = duds*(uu0 -u0c0*cmc0c0 + facg )
     1            + dcds*u0c0*(-2.*facg +u0c0*qbinv(i)/pi +
     2               (2.-u0c0)*(faclog*qbinv(i) -cmc0c0) +
     3               u0c0 *cmc0c0 )
                 duda(i) = dudp*(uu0 -u0c0*cmc0c0 + facg )
     1            + dcdp*u0c0*(-2.*facg +u0c0*qainv(i)/pi +
     2               (2.-u0c0)*(faclog*qainv(i) -cmc0c0) +
     3               u0c0 *cmc0c0 )
                 dudqa(i) = 
     1                   u0c0*(2.-u0c0)*dcdqa(i) 
     1                   +
     1                u0c0*u0c0*dcdp*a(i)/pi
                 dudqb(i) =               
     1                   u0c0*(2.-u0c0)*dcdqb(i) 
     1                  +
     1                u0c0*u0c0*dcds*b(i)/pi
            endif
  100   continue
        return
        end

        subroutine getgam(dcdb,b,qbinv,dcda,a,qainv,m,ilvry,cc,tp,
     1          gam)
        real*4 dcdb(*),b(*),dcda(*),a(*),qbinv(*), qainv(*)
c-----
c       determine spatial attenuation at constant frequency
c       and partial derivatives
c-----
        omega = 6.2831853/tp
        factr = 0.5*omega/(cc*cc)
        sumg = 0.0
        do 100 i=1,m
            sumg = sumg + dcdb(i)*b(i)*qbinv(i)
            if(ilvry.eq.2)then
                sumg = sumg + dcda(i)*a(i)*qainv(i)
            endif
  100   continue
        gam = factr * sumg
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

