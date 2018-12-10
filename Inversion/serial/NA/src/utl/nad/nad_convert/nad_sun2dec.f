c----------------------------------------------------------------------
c
c	nad_sun2dec - Program to read a Neighbourhood Algorithm Direct 
c	              access file written on a big endian machine 
c		      (e.g. a sun workstation) and replace it with the same
c		      file which can be read on a little endian machine 
c		      (e.g. a dec/compaq alpha). 
c
c		       This program should be compiled and run on a little
c		       endian machine like a dec alpha, because it uses
c		       the convert option in the open statement.
c
c	Comments:
c		       This program assumes a particular format for the
c		       reserved portion of the NAD header and does nothing to
c		       the user portion of the NAD header. The user
c		       portion will not be copied properly.
c
c		       This program can be used to convert the demo NAD files
c		       provided with the package from sun big endian to
c		       dec little endian format.
c
c	Calling sequence:
c
c			 nad_sun2dec file1 file2out 
c
c	where,
c	     `file1'   = input direct access `*.nad' file (from sun)
c	     `file2'   = output direct access `*.nad' file (for a dec)
c	
c----------------------------------------------------------------------
c
 	Program nad_sun2dec
c
	parameter	(nemax=101000,ndmax=24,nhmax=10000)
	
	real*4		models(ndmax*nemax)
	real*4		data(nemax)
	real*4		ranges(2,ndmax)
	real*4		scales(ndmax+1)
	character	header(nhmax)

      	character*256 	fnme
      	character*256 	fnmeo
      	character*256 	answer

        write(*,200)

        m = iargc()
        if(m.ne.2)then
           write(*,*)' Incorrect number of arguments'
           write(*,*)
           write(*,*)' Syntax:'
           write(*,*)'        nad_sun2dec filein fileout'
           write(*,*)
           stop
        end if

c						read in filenames
        call getarg(1,answer)
        read(answer,fmt='(a20)')fnme
        call getarg(2,answer)
        read(answer,fmt='(a20)')fnmeo
c						read in ensemble of models
c						from direct access nad file 
        lu_in = 15
	write(*,*)' Reading big endian direct access file...'

	call read_nad_be
     &               (lu_in,fnme,nh,nhu,iform,nd,ne,
     &                nhmax,ndmax,nemax,
     &                header,data,models)

	write(*,*)' Finished reading direct access file'
c
c						write out ensemble of models
c						in rfi ascii format
	nha = nh - nhu
	write(*,*)
	write(*,*)' Summary of nad file read in'
	write(*,*)
	write(*,*)' Number of variables                  = ',nd
	write(*,*)' Number of models                     = ',ne
	write(*,*)' Total length of header               = ',nh
	write(*,*)' Length of reserved header            = ',nha
	write(*,*)' Length of user header                = ',nhu
	a = data(1)
	avg = 0.
	mopt = 1
	do i=1,ne
	   if(data(i).lt.a)then
              a = data(i)
	      mopt = i
	   end if
	   avg = avg + data(i)
	end do
	avg = avg/real(ne)
	write(*,*)
	write(*,*)' Minimum misfit                       = ',a
	write(*,*)' Average misfit over all models       = ',avg
	write(*,*)' Index of lowest misfit model         = ',mopt
c
c						read nad header

	if(nhu.eq.nh)then
	   write(*,*)' No NA header in NAD file'
	else
	   lu = 20
 	   call read_NAheader
     &          (lu_in,fnme,nd,nh,nhu,iform,
     &           ranges,scales,ns1,ns,nr)

	   write(*,*)' Number of samples in first iteration = ',ns1
	   write(*,*)' Number of samples per iteration      = ',ns
	   write(*,*)' Number of cells per iteration        = ',nr
	   write(*,*)' ranges '
	   do k=1,nd
	      write(*,*)k,ranges(1,k),ranges(2,k),scales(k+1)
	   end do

	   write(*,*)
           lu = 10
           call NA_header
     &          (lu,fnmeo,header,nhmax,nh,nd,
     &           ranges,scales,ns1,ns,nr,nhu)

	end if
c
c                                               write direct access nad file
        write(*,*)
        write(*,*)' Writing big endian direct access file...'
        write(*,*)
        lu_out = 16
        call write_nad_le
     &       (lu_out,fnmeo,nd,ne,nh,nhu,header,1,models,data)

 100    format(1x,72("-"))
 200    format(/" Program nad_sun2dec - reads direct access nad file "/
     &        18x,"written on a sun and replaces it with one "/
     &        18x," for a big endian machine."/)

	stop
	end
c
c ---------------------------------------------------------------------
c
c       read_nad_be - read a direct access file in NAD format. 
c
c	Input:
c	      lu	      : logical unit of file
c	      fnme	      : filename 
c	      nhmax	      : maximum size of array header
c	      ndmax	      : maximum size of array data
c	      nemax	      : maximum size of array models
c
c	Output:
c	      nh	      : total length in bytes of file header 
c	      nhu	      : length in bytes of user written part of header 
c	      iform	      : format for nad file (1=multi-record,0=single)
c	      mul	      : number of records in the header (with padding)
c	      nd	      : dimension of parameter space
c	      ne	      : number of models in ensemble
c	      header          : header character string of length nh (char)
c	      data(nd)	      : array of data values for each model (real*4)
c	      models(nd,ne)   : array of model values  (real*4)
c
c	Comments:
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
c                File must contain a single record of length
c                [4x(4+nd*ne+ne+1) + nh] bytes
c
c		Calls are made to subroutine read_da. 
c
c                This routine assumes that direct access
c                files are opened with the record length specified
c                in bytes. This is the default for most machines
c                but not Dec machines. (Often a compiler option is
c                available on the DEC/compaq to use bytes rather than 
c                4-byte words.
c
c					M. Sambridge, RSES, Nov. 2001
c
c ---------------------------------------------------------------------
c
      	Subroutine read_nad_be
     &                     (lu,fnme,nh,nhu,iform,nd,ne,
     &                      nhmax,ndmax,nemax,
     &                      header,data,models)
c
      	real*4            models(*)
      	real*4            data(*)
      	character*(*)     header

      	character*256 	  fnme

	mul = 0
	iform = 0
c						read in size of header
        len = 4
        open(lu,file=fnme,status='unknown',
     &       form='unformatted',access='direct',recl=len,
     &       convert='big_endian')
        read(lu,rec=1)mul
        close(lu)
c						read in number of models,
c						number of dimensions and 
c						size of header
	if(mul.gt.0)then
c						we are in single record format
          len = 16
          open(lu,file=fnme,status='unknown',
     &         form='unformatted',access='direct',recl=len,
     &         convert='big_endian')
          read(lu,rec=1)nd,ne,nh,nhu
          mul = 0
          close(lu)
	else
c						we are in multi record format
          len = 20
          open(lu,file=fnme,status='unknown',
     &         form='unformatted',access='direct',recl=len,
     &         convert='big_endian')
          read(lu,rec=1)mul,nd,ne,nh,nhu
          mul = -mul
          close(lu)
	  iform = 1
	end if

c						check array sizes
	if(nd.gt.ndmax)then
          write(*,*)
          write(*,*)
     &    'Error - maxmimum number of dimensions is too small'
          write(*,*)
     &    '        current maximum number of dimensions =',ndmax
          write(*,*)
     &    '        number of dimensions in input file   =',nd
          write(*,*)'        '
          write(*,*)'Remedy - increase size of parameter ndmax'
          write(*,*)'         to at least this value and recompile'
          write(*,*)'        '
          stop
        end if

	if(ne.gt.nemax)then
          write(*,*)
          write(*,*)
     &    'Error - maxmimum number of models is too small'
          write(*,*)
     &    '        current maximum number of models =',nemax
          write(*,*)
     &    '        number of models in input file   =',ne
          write(*,*)'        '
          write(*,*)'Remedy - increase size of parameter nemax'
          write(*,*)'         to at least this value and recompile'
          write(*,*)'        '
          stop
        end if

	if(nh.gt.nhmax)then
          write(*,*)
          write(*,*)
     &    'Error - maxmimum size of header is too small'
          write(*,*)
     &    '        current maximum size of header =',nhmax
          write(*,*)
     &    '        size of header in input file   =',nh
          write(*,*)'        '
          write(*,*)'Remedy - increase size of parameter nhmax'
          write(*,*)'         to at least this value and recompile'
          write(*,*)'        '
          stop
        end if

 	call read_da(lu,fnme,nd,ne,nh,mul,header,models,data)

      	return
      	end

c
c ---------------------------------------------------------------------
c
c       read_da - read a direct access file containing
c                 ne models with dimension nd, ne data 
c                 values and a header of size nh.
c
c	Input:
c	      lu		: logical unit of file
c	      fnme		: filename 
c	      nh		: length in bytes of file header (minus padding)
c	      mul		: number of records in the header (with padding)
c	      nd		: dimension of parameter space
c	      ne		: number of models in ensemble
c
c	Output:
c	      header		: header character string of length nh (char)
c	      data(nd)		: array of data values for each model (real*4)
c	      models(nd,ne)	: array of model values  (real*4)
c
c	Comments:
c		 The direct access NAD file format:
c
c		 VARIABLE	TYPE		SIZE IN BYTES
c		 nd		int     	4
c		 ne		int     	4
c		 nh		int     	4
c		 header		character     	nh
c		 models		real*4     	4*nd*ne
c		 data		real*4     	4*nd
c		 tail		real*4     	4
c	
c	         File must contain a single record of length 
c		 [4x(3+nd*ne+ne+1) + nh] bytes
c
c
c					M. Sambridge, RSES, Nov. 2001
c
c ---------------------------------------------------------------------
c
      	Subroutine read_da(lu,fnme,nd,ne,nh,mul,header,models,data)
c
      	real*4            models(nd,ne)
      	real*4            data(ne)
      	character         header(nh)
c     	real*4            tail

      	character*256 	  fnme
c					
c                                              calculate length of header
 	iform = 0
 	if(mul.ne.0)iform = 1

 	lenh  = 4*5+nh
        len   = 4*(nd+1)

 	if(iform.eq.1)then

           open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=lenh,
     &          convert='big_endian')

       	   read(lu,rec=1)idum1,idum2,idum3,idum4,idum5,header
	
 	   close(lu)

 	   write(*,*)' record length               ',len
 	   write(*,*)' header length               ',lenh
 	   write(*,*)' Number of records in header ',mul
c
           open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=len,
     &          convert='big_endian')

           do i=1,ne
              call rnad(lu,i+mul,nd,models(1,i),data(i))
           end do

           close(lu)

 	else
c					       read in file as 
c					       a single record
	   len = 4+nd*ne+ne
           len = 4*len+nh

           open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=len,
     &          convert='big_endian')

           read(lu,rec=1)i,j,k,kk,header,models,data

           close(lu)

 	end if
c
      	return
      	end

      	Subroutine rnad(lu,i,nd,models,data)

      	real*4            models(nd)
      	real*4            data

      	read(lu,rec=i)models,data

      	return
      	end
c
c ---------------------------------------------------------------------
c
c       read_NAheader - extracts the NA program information placed
c                       in the header block of a NAD file
c
c ---------------------------------------------------------------------
c
        Subroutine read_NAheader
     &             (lu,fnme,nd,nh,nhu,iform,
     &              ranges,scales,ns1,ns,nr)

        real*4          ranges(2,nd)
        real*4          scales(nd+1)
        character*256   fnme

        if(nhu.ge.nh)return

        if(iform.eq.0)then
           len = 7+nd*3 + 1
           len = 4*len
           open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=len,
     &          convert='big_endian')
           read(lu,rec=1)i,j,k,kk,ns1,ns,nr,ranges,scales
           close(lu)
        else
           len = 4*5+nh
           open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=len,
     &          convert='big_endian')
           read(lu,rec=1)idum1,idum2,idum3,idum4,idum5,
     &                   ns1,ns,nr,ranges,scales
           close(lu)
        end if

        return
        end
c
c ----------------------------------------------------------------------------
c
c       read_header - converts header information into a character
c                     string by writing it to a direct access file
c                     and then reading it back as a character string
c
c       Calls no other routines.
c
c                                               M. Sambridge, June 1999
c
c ----------------------------------------------------------------------------
c
c                                               open direct access
c                                               temporary file
        Subroutine read_header
     &             (lu,fnme,nh,len,header)

        character         header(nh)
        character*256     fnme

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
c                      string by writing it to a direct access file
c                      and then reading it back as a character string
c
c       Calls no other routines.
c
c                                               M. Sambridge, June 1999
c
c ----------------------------------------------------------------------------
c
c                                               open direct access
c                                               temporary file
        Subroutine write_header
     &             (lu,fnme,len,nd,nh,
     &              range,scales,n1,n2,n3,header)

        real*4            range(2,nd)
        real*4            scales(nd+1)
        character         header(nh)
        character*256     fnme

        open(lu,file=fnme,status='unknown',
     &       form='unformatted',access='direct',recl=len)

        write(lu,rec=1)n1,n2,n3,range,scales,header

        close(lu)

        return
        end
c
c ----------------------------------------------------------------------------
c
c       NA_header - writes NA-specific information to NAD header.
c
c                   This routine adds various NA-header info to
c                   the header written by the user.
c
c       Calls no other routines.
c
c                                               M. Sambridge, June 1999
c
c ----------------------------------------------------------------------------
c
        Subroutine NA_header
     &             (lu,fnme,header,nh_max,nh,nd,
     &              range,scales,n1,n2,n3,nhu)

        real*4            range(2,nd)
        real*4            scales(nd+1)
        character*(*)     header
        character*256     fnme


c                                               calculate total header length
        rlen = 3*nd + 4
        len = 4*rlen+nh
        nhu = nh
        nh_na   = 4*rlen
        nh_tot  = len
c       write(*,*)' nh_user = ',nhu
c       write(*,*)' nh_na   = ',nh_na
c       write(*,*)' nh_tot  = ',nh_tot

        if(nh_tot.gt.nh_max)then
           write(*,*)
           write(*,*)' Error - header array too small'
           write(*,*)
           write(*,*)'         current size = ',nh_max
           write(*,*)'        required size = ',nh_tot
           write(*,*)
           write(*,*)' Remedy - adjust nh_max in parameter'
           write(*,*)'          to at least this value and recompile'
           write(*,*)
           stop
        end if

c                                               write out header information
        call write_header
     &       (lu,fnme,len,nd,nh,range,scales,n1,n2,n3,header)

        nh = nh_tot
c                                               read header information
c                                               into character string
        call read_header
     &       (lu,fnme,nh,len,header)


c       write(50,*)header(nh_na+1:nh_tot)

        return
        end
c
c ---------------------------------------------------------------------
c
c       write_nad_le - write a direct access file in 
c		       multi-record NAD format for little endian machine
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
      Subroutine write_nad_le
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
     &          form='unformatted',access='direct',recl=lenh,
     &          convert='little_endian')

c						write out header
c						for multi-record format 

         write(lu,rec=1)-mul,nd,ne,nh,nhu,header

	 close(lu)
c							write models
         open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=len2,
     &          convert='little_endian')

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
     &          form='unformatted',access='direct',recl=len,
     &          convert='little_endian')

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
