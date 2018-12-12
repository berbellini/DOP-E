        subroutine modls(id,in,nf1,nf10,dlam,invdep,lstinv,iter,
     1      rmsv, rmsq)
c-----
c       Changes
c       22 FEB 2002 - correct model comment
c       21 MAY 2002 - correct model comment for command line
c       08 JUL 2002 - introduce wc control
c       09 JAN 2003 - Katy Raven at Bullard Laborotories caught
c           problem with assigning density water
c           raven@esc.cam.ac.uk
c-----

c-----
        integer LER, LIN, LOT
        parameter(LER=0,LIN=5,LOT=6)
c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c       NL2 - number of columns in model (first NL2/2 are
c           - velocity parameters, second NL2/2 are Q values)
c-----
c     Calculate shear velocity model based on singular value
c     decomposition, standard deviations, and resolving kernels.
c-----
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)
        integer id, in, nf1, invdep, lstinv,iter
        integer nf10(NL)
        real dlam

        common/ctrl/numa,dd(NL),aa(NL),bb(NL),rr(NL),rt(NL),
     1      dw(NL2),x(NL2),
     2      h(NL2),u(NL2),c(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     3      wc(NL2)
        logical wc
        real d(NL),a(NL),b(NL),r(NL),rat(NL)
        real*4 cendep
        character*20 name
        character*50 fmt, fmt1
c-----
c       common for iputmod
c-----
        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        common/modtit/title
        character title*80
        character ntitle*80

        character ostr*12
        save itype
c-----
c       temporary variables
c-----
        real*4 oa(NL), ob(NL), od(NL), oqbinv(NL)


        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
        mlyr = MMAX
c-----
c       save current values
c-----
        do 39 i=1,mmax
            if(qb(i).gt.1.0)then
                qbinv(i) = 1.0/qb(i)
            else
                qbinv(i) =     qb(i)
            endif
            if(qa(i).gt.1.0)then
                qainv(i) = 1.0/qa(i)
            else
                qainv(i) =     qa(i)
            endif
            b(i) = vb(i)
            a(i) = va(i)
            d(i) = dl(i)
            oa(i) = a(i)
            ob(i) = b(i)
            r(i) = rho(i)
            oqbinv(i) = qbinv(i)
            od(i) = d(i)
   39   continue
c-----
c       get singular value decomposition
c-----
        call sgvsol(ntot,klim,m,err,dlam,mlyr,
     1      id,var,rmsv,rmsq,nf1,d)
        if(invdep.eq.1)then
c-----
c       use modified shear velocities to contruct compressional
c       velocity and density profile if desired
c
c       x() is the solution of Ax = y
c       x( 1   ...   NL ) is change is S velocity or layer thickness
c       x(NL+1 ... NL+NL) is change in QBINV
c-----
c
c-----
c           test for liquid surface layer - do not permit
c           a liquid to become a solid
c-----
            if(b(1) .eq.0.0)then
                x(1) = 0.0
                x(mlyr+1) = 0.0
            endif
c-----
c           apply the corrections
c-----
        do 40 i=1,m
            if(i.le.mlyr)then
                b(i)=b(i)+x(i)
                if(nf10(i).eq.0)then
c-----
c       P velocity is FIXED
c       P velocity is fixed if not low shear velocity sediment
c-----
c-----
c       force density  change if VS < 1.5/km/sec
c-----
                    if(b(i).le.1.5.and.b(i).gt.0.0)then
c-----
c       From Paul Mayne, density is
c       rho (gm/cc) = 0.8 log Vs(m/sec) for particulate material to
c                   = 0.8 log Vs(m/sec) - 0.3
c-----  
                        r(i)=0.8*alog10(b(i))+2.3
                    endif
                else
c-----
c     update compressional velocity and density
c     (assuming Poisson s ratio fixed)
c       note density changed according to P and also the S velocity
c-----
                    if(b(i).eq.0.0)then
                        r(i) = rho(i)
                    else
                    if(b(i).le.1.5)then
c-----
c       From Paul Mayne, density is
c       rho (gm/cc) = 0.8 log Vs(m/sec) for particulate material to
c                   = 0.8 log Vs(m/sec) - 0.3
c-----  
                        r(i)=0.8*alog10(b(i))+2.3
                    else 
                        call gtrofa(r(i),a(i))
                    endif
                    if(vb(i).gt.0.0)then
                        oldrat = va(i)/vb(i)
                        a(i) = oldrat * b(i)
                    endif
                    endif
                endif
c-----
c                   Q inversion
c-----
            else if(i.gt.mlyr)then
                if(qbinv(i-mlyr) .gt.0.0)then
                    oldrat = qainv(i-mlyr)/
     1                  qbinv(i-mlyr)
                else
                    oldrat = 1.0/2.25
                endif
                qbinv(i-mlyr) = qbinv(i-mlyr) + x(i)
                qainv(i-mlyr) = oldrat*
     1                  qbinv(i-mlyr)
            endif
   40   continue
c-----
c       update layer thicknesses, but do not permit any to be
c       negative of zero. Also if the thickness wants to
c       do to zero, then decrease the weight by a factor of
c       5
c-----
        else if(invdep.eq.0)then
            do 41 i=1,mlyr
                if(i.lt.mlyr)then
                    d(i) = d(i) + x(i)
                    if(d(i) .le.0.0)then
c-----
c       just decrease the layer thickness
c       slightly adjust the damping
c-----
                        call ddcur(m,m2)
                        dval = 0.95*dw(i)
                        if(dval.lt.0.01)then
                            dval = 0.01
                        endif
                        call ddupdt(i,dval)
                        d(i) = 0.75 * dl(i)
                    endif
                else
                    x(mlyr) = 0.0
                endif
   41       continue
        endif
c-----
c       parse according to value of ID
c-----
        if(id.eq.10)then
            write(LOT,*)'singular values:'
            write(LOT,*) (h(i),i=1,m)
        else if(id.eq.6)then
            open(2,file='tmpsrfi.02',access='sequential',
     1          status='unknown',form='formatted')
            rewind 2
            if(lstinv.eq.2 .or. lstinv.eq.5)then
                mlw = 1
                mup = mlyr
            else if(lstinv.eq.3)then
                mlw = mlyr +1
                mup = mlyr + mlyr
            else if(lstinv.eq.4)then
                mlw = 1
                mup = mlyr + mlyr
            endif
            write(2,'(7i5,2e16.8)')
     1          iunit,mlyr,itype,invdep,mlw,
     2          mup,lstinv,dlam
            do 146 ii=mlw, mup
                if(ii.gt.mlyr)then
                    i = ii - mlyr
                else
                    i = ii
                endif
                ddd = d(i)
                if(i.eq.mlyr)ddd = 0.0
                if(invdep.eq.1)then
                    write(2,'(i5,4e16.8/5x,5e16.8)')i,
     1                  cendep(i,d,mlyr),ddd,a(i),b(i),
     2                  qbinv(i),oa(i),ob(i),oqbinv(i)
                else
                    write(2,'(i5,4e16.8/5x,2e16.8)')i,
     1                  cendep(i,d,mlyr),ddd,a(i),b(i),
     2                  qbinv(i),od(i)
                endif
                write(2,'(5e16.8)')(v(ii,j),j=mlw,mup)
  146       continue
            close (2)
            do 147 i=1,mlyr
                ddd = d(i)
                if(i.eq.mlyr)ddd = 0.0
                dl(i) = ddd
                va(i) = a(i)
                vb(i) = b(i)
                rho(i) = r(i)
            QA(i) = qainv(i)
            QB(I) = qbinv(i)
            frefp(i) = 1.0
            frefs(i) = 1.0
            etap(i) = 0.0
            etas(i) = 0.0
  147       continue
c-----
c       lun I*4 - logical unit for writing model file. This
c                 unit is released after the use of this routine
c       nmmodl  C*(*)   - model name
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
c------
            iunit = 0
            iiso = 0
            mmax = mlyr
C           iflsph = 0
            idimen = 1
            icnvel = 0
            REFDEP = 0.0
            title = ' '
            lt = lgstr(title)
        call putmod(2,'tmpsrfi.17',mmax,title(1:lt),iunit,iiso,iflsph,
     1      idimen,icnvel,.false.)
        write(ostr,1)iter
    1   format('tmpmod96.',i3.3)
        call putmod(2,ostr,mmax,title(1:lt),iunit,iiso,iflsph,
     1      idimen,icnvel,.false.)
        else if(id.eq.13 .or. id.eq.18)then
                if(iunit.eq.0)then
                        fmt = '(3f10.4,2e10.3,f10.4)'
                else if(iunit.eq.1)then
                        fmt = '(3f10.1,2e10.3,f10.1)'
                else if(iunit.eq.2)then
                        fmt = '(3f10.1,2e10.3,f10.1)'
                endif
            if(lstinv.eq.2)then
                write(LOT,*)' INVERSION FOR S-VEL'
            else if(lstinv.eq.3)then
                write(LOT,*)' INVERSION FOR QB(INV)'
            else if(lstinv.eq.4)then
                write(LOT,*)' INVERSION FOR S-VEL and QB(INV)'
            else if(lstinv.eq.5)then
                write(LOT,*)' INVERSION FOR LAYER THICKNESS'
            endif
            write(LOT,*) ' Estimated data standard dev.:',var
            if((lstinv.eq.2.or.lstinv.eq.4).and.id.eq.18)then
                write(LOT,*)'     DEPTH', 
     1              ' THICKNESS',
     2              '     S-VEL',
     3              ' SIG DELVL',
     4              ' RESL in H',
     5              ' DEL (VEL)'
            write(LOT,*) ' RMS Velocity model perturbation   :',rmsv
            else if((lstinv.eq.3.or.lstinv.eq.4).and.id.eq.13)then
                write(LOT,*)'     DEPTH', 
     1              ' THICKNESS',
     2              '   QB(INV)',
     3              '       SIG DELQB',
     4              '     RESL in H',
     5              '     DEL QBINV'
            write(LOT,*) ' RMS Q-inverse model perturbation  :',rmsq
            else if(lstinv.eq.5)then
                write(LOT,*)'     DEPTH', 
     1              ' THICKNESS',
     2              '     S-VEL',
     3              ' SIG DEL H',
     4              '   DEL (H)'
            write(LOT,*) ' RMS layer thickness perturbation  :',rmsv
            endif
            do 185 i=1,mlyr
                if(i.eq.mlyr)then
                    ddd = 0.0
                else
                    ddd = d(i)
                endif
                if(id.eq.18)then
                    if(lstinv.ne.5)then
                    write(LOT,fmt)
     1              cendep(i,d,mlyr),ddd,b(i),c(i),
     2                  u(i),x(i)
                    else
                    write(LOT,fmt)
     1              cendep(i,d,mlyr),ddd,b(i),c(i),
     2                  x(i)
                    endif
                else if(id.eq.13)then
                    write(LOT,'(2f10.5,1p4e14.5)')
     1              cendep(i,d,mlyr),ddd,
     2              qbinv(i),c(mlyr+i),u(i),x(mlyr+i)
                        endif
  185       continue
        else if(id.eq.23 .or. id.eq.28)then
            if(numa.ne.0)then
                in = in + 1
                call getarg(in,name)
                write(ntitle,11)iter
   11           format('Model after ',i5,' iterations')
            else
                if(id.eq.28)then
                    write(LOT,*)'Enter File Name',
     1                  ' for Velocity Model'
                else if(id.eq.23)then
                    write(LOT,*)'Enter File Name',
     1                  ' for Qbinv Model'
                endif
                read(LIN,'(a)') name
                write(LOT,*)'Enter Model Title/comment'
                read(LIN,'(a)')ntitle
            endif
            open(4,file=name,access='sequential',
     1          status='unknown',form='formatted')
            rewind 4
            if(iunit.eq.0)then
                fmt = '(3f10.4,f10.4,e10.3,f10.4,e10.3)'
                fmt1= '(f10.4,2e11.4)'
            else if(iunit.eq.1)then
                fmt = '(3f10.1,f10.4,e10.3,f10.4,e10.3)'
                fmt1= '(f10.1,2e11.4)'
            else if(iunit.eq.2)then
                fmt = '(3f10.1,f10.4,e10.3,f10.4,e10.3)'
                fmt1= '(f10.1,2e11.4)'
            endif
c-----
c       output last updated model
c-----
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            mlyr = mmax
c-----
c       lun I*4 - logical unit for writing model file. This
c                 unit is released after the use of this routine
c       nmmodl  C*(*)   - model name
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
c------
        REFDEP = 0.0
        lt = lgstr(ntitle)
        call putmod(4,name,mmax,ntitle(1:lt),iunit,iiso,iflsph,
     1      idimen,icnvel,.true.)
c-----
c       output resolution matrix by row
c-----
        else if(id.eq.29 .or. id.eq.24)then
            if(numa.ne.0) then
                in=in+1
                call getarg(in,name)
            else
                if(id.eq.29)then
                    write(LOT,*)'Enter File Name for'
     1                  ,' Velocity Resolution by row'
                else if(id.eq.24)then
                    write(LOT,*)'Enter File Name for'
     1                  ,'Qbinv resolution by row'
                endif
                read(LIN,'(a)') name
            endif
            open(4,file=name,access='sequential',status='unknown',
     1          form='formatted')
            rewind 4
            if(id.eq.29)then
                itype = 0
            else if(id.eq.24)then
                itype = 1
            endif
            write(4,'(3i5,f10.3,8x,e12.5,i5)')
     1          iunit,mlyr,itype,dlam,invdep
c-----
c       output last updated model
c-----
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            mlyr = mmax
            dtop = 0.0
            do 136 i=1,mlyr
                ddd=dl(i)
                if(i.eq.m)ddd=0.0
                if(invdep.eq.0)then
                    dval = dtop
                    dtop = dtop + ddd
                else if(invdep.eq.1)then
                    dval = cendep(i,d,mlyr)
                endif
                if(id.eq.29)then
                    write(4,*)i,dval,ddd,a(i),
     1                  b(i),r(i),rat(i)
                    write(4,'(8f10.3)')(v(i,j),j=1,mlyr)
                else if(id.eq.24)then
                    write(4,*)i,dval,ddd,
     1                  qbinv(i)
                    write(4,'(8f10.3)')(v(i+mlyr,j+mlyr),
     1                  j=1,mlyr)
                endif
  136       continue
            close(4,status='keep')
c-----
c       output resolution matrix by row
c-----
        else if(id.eq.14 .or. id.eq.19)then
            write(LOT,*) mlyr, 'layers'
            write(LOT,*)' '
            if(id.eq.19)then
                write(LOT,*)'layer, dpth, thick, vp,'
     1              ,' Vs, resolving kernel by row'
            else if(id.eq.14)then
                write(LOT,*)'layer, dpth, thick,'
     1              ,' Qbinv, Resolving kernel by row'
            endif
            write(LOT,*) ' '
            do 135 i=1,mlyr
                ddd = dl(i)
                if(i.eq.mlyr)ddd = 0.0
                if(id.eq.19)then
                    write(LOT,*)i,cendep(i,d,mlyr),ddd,a(i),b(i)
                    write(LOT,*)' '
                    write(LOT,'(8f10.3)')(v(i,j),j=1,mlyr)
                else if(id.eq.14)then
                    write(LOT,*)i,cendep(i,d,mlyr),ddd,qbinv(i)
                    write(LOT,*)' '
                    write(LOT,'(8f10.3)')(v(i+mlyr,j+mlyr),
     1                  j=1,mlyr)
                endif
  135       continue
        else if(id.eq.9)then
            open(2,file='tmpsrfi.02',access='sequential',
     1          status='unknown',form='formatted')
            rewind 2
            if(lstinv.eq.2 .or. lstinv.eq.5)then
                mlw = 1
                mup = mlyr
            else if(lstinv.eq.3)then
                mlw = mlyr +1
                mup = mlyr + mlyr
            else if(lstinv.eq.4)then
                mlw = 1
                mup = mlyr + mlyr
            endif
            write(2,'(7i5,2e16.8)')
     1          iunit,mlyr,itype,invdep,mlw,
     2          mup,lstinv,dlam
            do 145 ii=mlw, mup
                if(ii.gt.mlyr)then
                    i = ii - mlyr
                else
                    i = ii
                endif
                ddd = d(i)
                if(i.eq.mlyr)ddd = 0.0
                if(invdep.eq.1)then
                    write(2,'(i5,4e16.8/5x,3e16.8)')i,
     1                  cendep(i,d,mlyr),ddd,a(i),b(i),
     2                  qbinv(i),ob(i),oqbinv(i)
                else
                    write(2,'(i5,4e16.8/5x,2e16.8)')i,
     1                  cendep(i,d,mlyr),ddd,a(i),b(i),
     2                  qbinv(i),od(i)
                endif
                write(2,'(5e16.8)')(v(ii,j),j=mlw,mup)
  145       continue
            close (2)
        endif
        return
        end

        function cendep(i,d,mlyr)
c-----
c       determine the depth to the center of layer i
c-----
        parameter (NL=200,NL2=NL+NL)
        real*4 d(NL)
c-----
c       i   - layer index
c       d   - thickness of i'th layer
c           - halfspace thickness means nothing
c       mlyr    - number of layers including halfspace
c-----
        cendep = 0.0
        if(i.lt.mlyr)then
            imx = i
        else
            imx = mlyr - 1
        endif
        do 128 j=1,imx
            cendep = cendep + d(j)
  128   continue
        if(i.lt.mlyr)then
            cendep = cendep - d(i)/2
        else
            cendep = cendep + d(mlyr-1)/2
        endif
        return
        end

        subroutine sgvsol(ntot,klim,m,err,dlam,
     1      mlyr,id,var,rmsv,rmsq,nf1,d)
c-----
c       general subroutine to return singular value decomposition 
c       as well as solution vector
c-----
c       ntot
c       klim
c       m   - number of unknowns
c       err -
c       h   - vector
c       u   - vector
c       v   - array
c       x   - solution vector
c       dlam    - damping factor
c       b   - layer velocity array
c       qbinv   - layer q-beta inverse array
c       mlyr    - number of layers in model
c       id  - menu command
c       c   - vector standard error estimate
c       var - variance of fit to observed data
c       rmsv    - rms change in velocity/thickness model fit
c       rmsq    - rms change in qinverse model fit
c       dd  - dummy array for temporary storage
c       mlyr    - number of layers in model (m = mlyr+mlyr)
c       d   - array of layer thicknesses
c       nf1 - = 1 use model variance for standard error estimates
c             = 0 do not use model variance
c-----
c-----
        integer NL, NL2
        parameter (NL=200,NL2=NL+NL)
        common/ctrl/numa,dd(NL),a(NL),b(NL),r(NL),rat(NL),
     1      d2(NL2),x(NL2),
     $      h(NL2),u(NL2),c(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        real d(NL)
        real*8 sum,sum1v,sum1q
        real*4 cendep
        real dtmp(NL2)
c-----
        open(4,file='tmpsrfi.10',form='unformatted',
     1      access='sequential')
        rewind 4
c
c     ntot: total number of rows
c     klim: number of non-zero eigenvalues
c     m:    number of columns
c     err:  norm of minimum error between observed and computed
c             values:  err = min || b - Ax ||
c     h:    vector of eigenvalues
c     u:    vector of data eigenvector matrix U times residual
c           vector b: 'utb vector'.
c     v:    model eigenvector matrix V.
c
        read(4) ntot,klim,m,err
        read(4) (h(i), i=1,m)
        read(4) (u(i), i=1,m)
        do 5 i=1,m
            read(4) (v(i,j),j=1,m)
c           write(LOT,*)i, (v(i,j),j=1,m)
    5   continue
        close(4,status='keep')
        if(id.eq.10)return
c-----
c       calculate solution
c            -1    2   2  -1    T
c       x = W  V( L + s I )  L U b
c-----
        sum1v = 0.0d+00
        sum1q = 0.0d+00
        do 35 i=1,m
            sum = 0.0d+00
            do 30 j=1,klim
                sum=sum+dble(u(j)*v(i,j)*h(j)/
     1              (h(j)*h(j) + dlam) )
   30       continue
            x(i) = sngl(sum)
            if(i.le.mlyr)then
                sum1v = sum1v + dble(x(i)**2)
            else
                sum1q = sum1q + dble(x(i)**2)
            endif
   35   continue
        rmsv = sngl(sum1v)/float(mlyr)
        rmsv = sqrt(rmsv)
        rmsq = sngl(sum1q)/float(mlyr)
        rmsq = sqrt(rmsq)
c-----
c     calculate model standard deviations
c     c: vector of standard deviations
c-----
        sum=dble(err*err)
        do 60 i=klim+1,m
   60   sum=sum+u(i)**2
        sn=float(ntot)
        sk=float(klim)
        var=sngl(sum)
        if(dlam.ne.0.0) then
            do 65 i=1,klim
   65           sum=sum+dble((u(i)*
     1              dlam/(h(i)*h(i)+dlam) )**2)
            var=sngl(sum)
            sum=0.0d0
            do 70 i=1,klim
   70       sum=sum+dble(h(i)*h(i)/(h(i)**2 + dlam))
            sk=sngl(sum)
        endif
        if(sn.gt.sk) var=var/(sn-sk)
        if(sn.le.sk) var=1.
        var=sqrt(abs(var))
        do 85 i=1,m
            sum=0.0d0
            do 80 j=1,klim
   80           sum=sum+dble((v(i,j)*h(j)/
     1              (h(j)*h(j)+dlam) )**2)
            c(i)=sngl(sum)
            c(i)=sqrt(c(i))
            if(nf1.eq.1) c(i)=var*c(i)
   85   continue
c-----
c       compute resolution matrix
c-----
        open(4,file='tmpsrfi.10',form='unformatted',
     1      access='sequential')
        rewind 4
        open(2,file='tmpsrfi.11',form='unformatted',access='sequential')
        rewind 2
        open(3,file='tmpsrfi.14',form='unformatted',access='sequential')
        rewind 3
c-----
c     calculate resolving kernel matrix
c-----
        read(4) ntot,klim,m,err
        read(4) (h(i),i=1,m)
        read(4) (u(i),i=1,m)
        do 100 i=1,m
  100   read(2) (v(i,j),j=1,m)
        do 110 i=1,klim
            do 105 j=1,m
  105       v(i,j)=v(i,j)*h(i)/(h(i)+dlam/h(i))
  110   continue
        do 125 i=1,m
            read(4) (u(jj),jj=1,m)
            do 120 j=1,m
                sum=0.0d0
                do 115 k=1,klim
  115               sum=sum+u(k)*v(k,j)
  120       dtmp(j)=sngl(sum)
c-----
c       dd is a temporary store of V(L+2+sig I)^-1 L^2
c-----
            write(3) (dtmp(jj),jj=1,m)
  125   continue
        rewind 3
        do 130 i=1,m
        read(3) (v(i,j),j=1,m)
  130   continue
        close(4,status='keep')
        close(2,status='keep')
        close(3,status='delete')
        if(id.ne.13 .and. id.ne.18)return
c-----
c       estimate resolution of a layer by taking second moment of
c       resolution matrix and layer center point depth
c-----
        if(dlam.eq.0.0)then
            do 146 i=1,mlyr
                h(i)=cendep(i,d,mlyr)
                u(i)=d(i)
  146       continue
            u(m) = 0.0
        else
            do 170 i=1,mlyr
                sum=0.0d0
                do 150 j=1,mlyr
                    if(id.eq.18)then
                        k=i
                        l=j
                    else if(id.eq.13)then
                        k=i+mlyr
                        l=j+mlyr
                    endif
                    sum=sum+dble(v(k,l)**2)
  150           continue
                stot=sngl(sum)
                if(stot.eq.0.0)then
                    h(i)=cendep(i,d,mlyr)
                    u(i)=0.0
                else
                    sum=0.0d0
                    do 160 j=1,mlyr
                        if(id.eq.18)then
                            k=i
                            l=j
                        else if(id.eq.13)then
                            k=i+mlyr
                            l=j+mlyr
                        endif
                        sum=sum+dble(v(k,l)**2*
     1                      cendep(j,d,mlyr))
  160               continue
                    h(i)=sngl(sum)/stot
                    sum=0.0d0
                    do 165 j=1,mlyr
                        if(id.eq.18)then
                            k=i
                            l=j
                        else if(id.eq.13)then
                            k=i+mlyr
                            l=j+mlyr
                        endif
                        sum=sum+dble((v(k,l)*
     1                  (cendep(j,d,mlyr)-h(i)))**2)
  165               continue
                    u(i)=sngl(sum)/stot
                    u(i)=sqrt(u(i))*2.
                    if(u(i).lt.d(i)) u(i)=d(i)
                endif
  170       continue
        endif
        return
        end
