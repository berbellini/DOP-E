        program mkmod96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: MKMOD96                                               c
c                                                                     c
c      COPYRIGHT 2001 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       This program interactively creates model96 models files
c-----
c       Changes
c       07 MAR 2002 - program did not write the density
c       04 APR 2002 - arguments in putmod were incorrect
c
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        
        character*80  mname
        character*80  title

        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/depref/refdep

c-----
c       enter name of mode file
c-----
        write(LOT,*)'Write creating model96 file for isotropic ',
     1      'constant velocity layers, 1-D model'
        write(LOT,*)'Enter name of the earth model file'
        read(LIN,'(a)')mname
        lm = lgstr(mname)
        write(LOT,*)' Model file is :', mname(1:lm)
c-----
c       Get model comment
c-----
        write(LOT,*)'Enter model comment'
        read(LIN,'(a)')title
        lt = lgstr(title)
        write(LOT,*)' Comment is    :', title(1:lt)
c-----
        lun = 1
c-----
c       Isotropic model
c-----
        iiso = 0
c-----
c       Units are KGS
c-----
        iunit = 0
c-----
c       FLAT or SPHERICAL model
c-----
 1000   continue
        write(LOT,*)'Enter 0 for flat earth model'
        write(LOT,*)'      1 for spherical earth model'
        read(LIN,*)iflsph
        if(iflsph.lt.0 .or. iflsph.gt.1)go to 1000
        write(LOT,*)' Model flat/sph:', iflsph
c-----
c       Model dimension
c-----
        idimen = 1
c-----
c       constant velocity layers
c-----
        icnvel = 0
c-----
c       read in the model
c-----
        write(LOT,*)'Enter Velocity Model, EOF (CTRL D) to end:'
        write(LOT,*)' H     VP     VS      RHO    QP QS ETAP ETAS',
     1      ' FREFP FREFS'
        WRITE(LOT,*)'(km) (km/s) (km/s) (gm/cm^3) -- -- ---- ----',
     1      '  (Hz)  (Hz)'
        mmax = 0
 2000   continue
        read(LIN,*,end=2001,err=2000)hh,vp,vs,r,qp,qs,ep,es,fp,fs
        write(LOT,*)hh,vp,vs,r,qp,qs,ep,es,fp,fs
        mmax = mmax + 1
        d(mmax) = hh
        a(mmax) = vp
        b(mmax) = vs
        rho(mmax) = r
        qa(mmax) = qp
        qb(mmax) = qs
        etap(mmax) = ep
        etas(mmax) = es
        frefp(mmax) = fp
        frefs(mmax) = fs
        go to 2000
 2001   continue
c-----
c       set the reference depth to zero - c
c                     it is the user responsibility to
c       put in the negative layer thicknesses c
c                     to indicate above sea level
c-----
        refdep = 0
        WRITE(LOT,*)'Creating the model file:',mname(1:lm)
c-----
c       create the model file
c-----
        call putmod(lun,mname(1:lm),mmax,title(1:lt),
     1      iunit,iiso,iflsph,idimen,icnvel,.true.)
        end
        

