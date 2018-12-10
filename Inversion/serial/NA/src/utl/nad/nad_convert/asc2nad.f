c----------------------------------------------------------------------
c
c	asc2nad - Program to convert ascii file to NAD file.
c
c	Calling sequence:
c
c			 asc2nad file.in file.out
c
c	where,
c	     `file.in'  = input asc file (e.g. `NAV22.mdl')
c	     `file.out' = output direct access `*.nad' file (e.g. `NAV22.nad')
c	
c----------------------------------------------------------------------
c
 	Program asc2nad
c
	parameter	(nemax=11000,ndmax=24,maxlayers=6)
	parameter	(nhmax=24+12*maxlayers+10000)
	
	real*4		models(ndmax*nemax)
	real*4		data(nemax)
c	real*4		vbest(ndmax)
c	real*4		qalpha(maxlayers)
c	real*4		qbeta(maxlayers)
	real*4		ranges(2,ndmax)
	real*4		scales(ndmax+1)
c						length of string has
c						to be at least nhmax
	character	header(nhmax)

      	character*256 	fnme_in
      	character*256 	fnme
      	character*256 	fnme_out
      	character*256 	answer

        write(*,100)

	m = iargc()
	if(m.ne.2)then
	   write(*,*)' Incorrect number of arguments'
	   write(*,*)
	   write(*,*)' Syntax:'
	   write(*,*)'        asc2nad filein fileout'
	   write(*,*)'        (e.g. NAV22.mdl NAV22.nad)'
	   write(*,*)
           stop
        end if
c
c						read in file names 
        call getarg(1,answer)
        read(answer,fmt='(a20)')fnme_in
        call getarg(2,answer)
        read(answer,fmt='(a20)')fnme_out

c						read in ensemble of models
c						from asc file 
	write(*,*)' Reading models...'
        lu_in = 15

	nlayers = 6
	call read_asc
     &       (lu_in,fnme_in,ne,nd,ns1,ns,nr,nh,ndmax,nemax,
     &        nhu,ranges,scales,models,data)

	na_head = 0
c						write out ensemble of models
c						as a direct access file 
        lu_out = 10

	   fnme = 'junk'
           call NA_header
     &          (lu_in,fnme,header,nhmax,nh,nd,
     &           ranges,scales,ns1,ns,nr,nhu)

	write(*,*)
	write(*,*)' Writing direct access file...'
	write(*,*)
c						write direct access nad file

 	call write_nad
     &       (lu_out,fnme_out,nd,ne,nh,nhu,header,1,models,data)

 200    format(/"  Note that the internal NA portion of",
     &          " the header file"/
     &      19x," will not be written to nad header."/)

 100    format(/"  Program asc2nad - converts",
     &          " ascii receiver function model file "/
     &      19x," to direct access nad file."/)

	stop
	end
c
c ---------------------------------------------------------------------
c
c       write_nad - write a direct access file in 
c		    multi-record NAD format
c
c       Input:
c             lu                : logical unit of file
c             fnme              : filename
c             nhmax             : maximum size of array header
c             ndmax             : maximum size of array data
c             nemax             : maximum size of array models
c	      iform		: =0 then single record format
c				  =1 then multi-record format (for large nd)
c
c       Output:
c             nh                : length in bytes of file header
c             nhu               : length in bytes of user portion of header
c             nd                : dimension of parameter space
c             ne                : number of models in ensemble
c             header            : header character string of length nh (char)
c             data(nd)          : array of data values for each model (real*4)
c             models(nd,ne)     : array of model values  (real*4)
c
c       Comments:
c                The direct access NAD file format:
c
c                VARIABLE       TYPE            SIZE IN BYTES
c                nd             int             4
c                ne             int             4
c                nh             int             4
c                nhu            int             4
c                header         character       nh
c                models         real*4          4*nd*ne
c                data           real*4          4*nd
c                tail           real*4          4
c
c		 In single record mode a direct access file of length
c                [4x(4+nd*ne+ne+1) + nh] bytes is produced.
c
c		 In multi record mode a direct access file of length
c                [(ne+1)*(max(20+nh,4(nd+1))] bytes is produced.
c
c       	Calls are made to subroutine read_da.
c
c                This routine assumes that direct access
c                files are opened with the record length specified
c                in bytes. This is the default for most machines
c                but not Dec machines. (Often a compiler option is
c                available on the DEC/compaq to use bytes rather than 
c                4-byte words.
c
c                                       M. Sambridge, RSES, November 2001
c
c ---------------------------------------------------------------------
c
      Subroutine write_nad
     &           (lu,fnme,nd,ne,nh,nhu,header,iform,models,data)
c
      real*4            models(nd,ne)
      real*4            data(ne)
      character         header(nh)
      real*4            tail
      character*256 	fnme
      logical		warn

      warn = .true.
      warn = .false.

c						write new 
c						multi-record format
      if(iform.eq.1)then

c						calculate length of header
	 len1 = 4*5+nh
	 len2 = 4*(nd+1)
         mul  = 1 + (len1-1)/len2
c        write(*,*)mul
         lenh = mul*len2
         num = ne + mul
         is1 = num*len2 
         is2 = 4*(5+nd*ne+ne)+nh

         write(*,*)' Number of models                         :',ne
         write(*,*)' Number of dimensions                     :',nd
         write(*,*)' Original header length in bytes          :',len1
         write(*,*)' Final header length in bytes             :',lenh
         write(*,*)' Direct access file record length         :',len2
         write(*,*)' Number of records                        :',num 
         write(*,*)' Size of nad file in multi-record format  :',is1
         write(*,*)' Size of nad file in single record format :',is2
         
c							write header
         open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=lenh)

c						write out header
c						for multi-record format 

         write(lu,rec=1)-mul,nd,ne,nh,nhu,header

	 close(lu)
c							write models
         open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=len2)

         do i=1,ne
            call wnad(lu,mul+i,nd,models(1,i),data(i))
         end do
         close(lu)


      else
c						write original 
c						single record format
c					
c						set total length of nad file
         len = 4+nd*ne+ne+1
         len = 4*len+nh

         write(*,*)' size of nad file = ',len
c						open direct access nad file

         open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=len)

         tail = -999.0
         write(lu,rec=1)nd,ne,nh,nhu,header,models,data,tail
 
         close(lu)

      end if

      return
      end

      Subroutine wnad(lu,i,nd,models,data)

      real*4            models(nd)
      real*4            data

      write(lu,rec=i)models,data

      return
      end

c ----------------------------------------------------------------------------
c
c       NA_header - writes NA-specific information to NAD header.
c
c		    This routine adds various NA-header info to
c		    the header written by the user.
c
c						M. Sambridge, June 1999
c
c ----------------------------------------------------------------------------
c
        Subroutine NA_header
     &             (lu,fnme,header,nh_max,nh,nd,
     &              range,scales,n1,n2,n3,nhu)

	real*4		  range(2,nd)
	real*4		  scales(nd+1)
        character*(*)     header
        character*256 	  fnme

c                                               calculate total header length
	rlen = 3*nd + 4
        len = 4*rlen+nh
	nhu = nh
	nh_na   = 4*rlen
	nh_tot  = len
c	write(*,*)' nh_user = ',nhu
c	write(*,*)' nh_na   = ',nh_na  
c	write(*,*)' nh_tot  = ',nh_tot 

        if(nh_tot.gt.nh_max)then
           write(*,*)
           write(*,*)' Error - header array too small'
           write(*,*)
           write(*,*)'         current size = ',nh_max
           write(*,*)'        required size = ',nh_tot
           write(*,*)
           write(*,*)' Remedy - adjust nh_max in parameter',
     &                        ' file and recompile'
           write(*,*)
           stop
        end if

c						write out header information
	call write_header
     &       (lu,fnme,len,nd,nh,range,scales,n1,n2,n3,header)

	nh = nh_tot
c						read header information
c						into character string
	call read_header
     &       (lu,fnme,nh,len,header)


c	write(50,*)header(nh_na+1:nh_tot)

	return
	end
c
c ----------------------------------------------------------------------------
c
c       read_header - converts header information into a character
c		      string by writing it to a direct access file
c		      and then reading it back as a character string
c
c       Calls no other routines.
c
c						M. Sambridge, June 1999
c
c ----------------------------------------------------------------------------
c
c                                               open direct access 
c						temporary file
	Subroutine read_header
     &  	   (lu,fnme,nh,len,header)

        character         header(nh)
        character*256 	  fnme

        open(lu,file=fnme,status='unknown',
     &       form='unformatted',access='direct',recl=len)

        read(lu,rec=1)header

        close(lu)

	return
	end
c
c ----------------------------------------------------------------------------
c
c       write_header - converts header information into a character
c		       string by writing it to a direct access file
c		       and then reading it back as a character string
c
c       Calls no other routines.
c
c						M. Sambridge, June 1999
c
c ----------------------------------------------------------------------------
c
c                                               open direct access 
c						temporary file
	Subroutine write_header
     &  	   (lu,fnme,len,nd,nh,
     &              range,scales,n1,n2,n3,header)

	real*4		  range(2,nd)
	real*4		  scales(nd+1)
        character         header(nh)
        character*256 	  fnme

        open(lu,file=fnme,status='unknown',
     &       form='unformatted',access='direct',recl=len)

        write(lu,rec=1)n1,n2,n3,range,scales,header

        close(lu)

	return
	end

c
c ---------------------------------------------------------------------
c
c       write_rfi - write ensemble of models produced by
c                   narfi and garfi programs.
c
c ---------------------------------------------------------------------
c
	Subroutine read_asc
     &       (lu,fnme2,ne,nd,ns1,ns,nr,nh,ndmax,nemax,
     &        nhu,ranges,scales,models,data)

        real*4          models(nd,ne)
        real*4          data(nd)
        real*4          ranges(2,nd)
        real*4          scales(nd+1)

        character*256   fnme2

        open(lu,file=fnme2,status='old')

        read(lu,*)
        read(lu,*)
        read(lu,*)
        read(lu,*)nd
	if(nd.gt.ndmax)then
	   write(*,*)
	   write(*,*)' Program asc2nad'
	   write(*,*)
	   write(*,*)' Maximum number of variables is too small'
	   write(*,*)' Number in input file = ',nd
	   write(*,*)' Maximum possible     = ',ndmax
	   write(*,*)
	   write(*,*)' Remedy: increase size of array parameter'
	   write(*,*)'         ndmax and recompile'
	   write(*,*)
	   stop
	end if
        read(lu,*)ne
	if(ne.gt.nemax)then
	   write(*,*)
	   write(*,*)' Program asc2nad'
	   write(*,*)
	   write(*,*)' Maximum number of models is too small'
	   write(*,*)' Number in input file = ',ne
	   write(*,*)' Maximum possible     = ',nemax
	   write(*,*)
	   write(*,*)' Remedy: increase size of array parameter'
	   write(*,*)'         nemax and recompile'
	   write(*,*)
	   stop
	end if
        read(lu,*)nh
c       read(lu,*)nhu
        read(lu,*)ns1
        read(lu,*)ns
        read(lu,*)nr
        read(lu,*)scales(1)
        read(lu,*)
        read(lu,*)
        do k=1,nd
           read(lu,*)i,ranges(1,k),ranges(2,k),scales(k+1)
        end do
        read(lu,*)
        read(lu,*)
        read(lu,*)
        do i=1,ne
           read(lu,*)k,(models(k,i),k=1,nd),data(i)
           read(lu,*)
        end do
	nhu = 0

        close(lu)

        return
        end
