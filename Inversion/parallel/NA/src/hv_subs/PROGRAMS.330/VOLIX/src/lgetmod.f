        subroutine getmod(lun,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,listmd,zmin,zmax,xmin,xmax)
c-----
c     General purpose model input
c     This model specification is designed to be as general as possible
c
c     Input lines
c     Line 01:    MODEL
c     Line 02:    Model Name
c     Line 03:    ISOTROPIC or ANISOTROPIC or TRANSVERSELY ANISOTROPIC
c     Line 04:    Model Units, First character is length 
c                 (k for kilometer), second is mass (g for gm/cc), 
c                 third is time (s for time)
c     Line 05:    FLAT EARTH or SPHERICAL EARTH
c     Line 06:    1-D, 2-D or 3-D
c     Line 07:    CONSTANT VELOCITY
c     Line 08: open for future use
c     Line 09: open for future use
c     Line 10: open for future use
c     Line 11: open for future use
c     Lines 12-end:   These are specific to the model
c-----
c     lun I*4 - logical unit for reading model file. This
c               unit is released after the use of this routine
c     mname   C*(*)   - model name
c     mmax    I*4 - number of layers in the model, last layer is
c                  halfspace
c     title   C*(*)   - title of the model file
c     iunit   I*4 - 0 Kilometer, Gram, Sec
c     iiso    I*4 - 0 isotropic 
c               1 transversely anisotropic 
c               2 general anisotropic 
c     iflsph  I*4 - 0 flat earth model
c               1 spherical earth model
c     idimen  I*4 - 1 1-D
c             - 2 2-D
c             - 3 3-D
c     icnvel  I*4 - 0 constant velocity
c               1 variable velocity
c     ierr    I*4 - 0 model file correctly read in
c             - -1 file does not exist
c             - -2 file is not a model file
c               -3 error in the model file
c     listmd  L   - .true. list the model
c------

        character mname*(*), title*(*)
        integer*4 mmax, iunit, iiso, iflsph, idimen
        integer*4 ierr
        character string*80
        character line11*80
        logical listmd
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c       LER I*4 - logical unit for standard error
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
c-----
c       common for laterally varying model
c-----
        integer NLAYER, NNODE, NDIST
        parameter (NLAYER=20, NNODE=20, NDIST=100)
        common/latbnd/xc, zc, nc, ii
        real*4 xc(NLAYER,NNODE), zc(NLAYER,NNODE)
        integer*4 nc(NLAYER), ii(NLAYER,NNODE)

        common/latval/vp,vs,rho,qp,etap,qs,etas
        real*4 vp(2,NLAYER), vs(2,NLAYER), rho(2,NLAYER), 
     1      qp(2,NLAYER), etap(2,NLAYER), 
     2      qs(2,NLAYER), etas(2,NLAYER)
c-----
c       xc(i,j), yc(i,j) j=1,nc(i)
c           coordinate pairs of layer boundary
c       ii(i,j) a parameter for interpolation used by Cerveny Program
c           set to -1 to avoid spline interpolation
c       vp(i)   P veloctiy in Layer
c       dvp(i)  dP/dz in layer   P-veloctiy = vp(i) + dvp(i)*z
c       vs(i)   S velocity in layer
c       dvs(i)  dS/dP in layer
c       rho(i)  density velocity in layer
c       qp(i)   Quality factor in layer
c       qs(i)   Quality factor in layer
c       etap(i)     frequency dependence in layer
c       etas(i)     frequency dependence in layer
c-----

        common/depref/refdep

        logical ext
        character ftype*80

c-----
c       test to see if the file exists
c-----
        ierr = 0
        inquire(file=mname,exist=ext)
        if(.not.ext)then
            ierr = -1
            write(LER,*)'Model file does not exist'
            return
        endif
c-----
c       open the file
c-----
        open(lun,file=mname,status='old',form='formatted',
     1      access='sequential')
        rewind lun
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
c       get lines 8 through 10
c-----
        do 900 i=8,10
            read(lun,'(a)')string
  900   continue
        read(lun,'(a)')line11
c-----
c       error checking
c-----
        if(idimen.lt.0 .or. idimen.gt.2)then
            write(LER,*)'Only 1-D or 2-D model permitted'
            ierr = -1
            return
        endif
        if(iflsph.ne.0)then
            write(LER,*)'Spherical Earth not permitted'
            ierr = -1
            return
        endif
        if(iiso.ne.0)then
            write(LER,*)'Only isotropic model permitted'
            ierr = -1
            return
        endif
c-----
c       here is the complicated stuff
c-----
        if(idimen.eq.1)then
            if(icnvel.eq.0)then
                call gt1dcv(lun,mmax,listmd,xmin,xmax,zmin,zmax)
            else if(icnvel.eq.1)then
                call gt1dvv(lun,mmax,listmd,xmin,xmax,zmin,zmax)
            endif
        else if(idimen.eq.2)then
            if(line11 .eq. 'CERVENY' .or. line11 .eq. 'cerveny')then
            if(icnvel.eq.0)then
                call gt2dcv(lun,mmax,listmd,xmin,xmax,zmin,zmax)
            else if(icnvel.eq.1)then
                call gt2dvv(lun,mmax,listmd,xmin,xmax,zmin,zmax)
            endif
            else
               write(LER,*)'2-D model missing CERVENY in Line 11'
                ierr = -1
            endif
        endif
        close (lun)
        return
        end
