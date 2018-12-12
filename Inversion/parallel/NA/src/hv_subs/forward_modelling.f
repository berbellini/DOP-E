c------------------------------------------------------------------------
c
c	Subroutine forward_modelling - calculates the array "predicted_
c				       data" from the model "rmodel"
c
c	Note_1: Calls are made to : theo
c
c       Note_2: rmodel(1:nlayer): thickness
c               rmodel(nlayer+1:2*nlayer): velocity at upper interface
c               rmodel(2*nlayer+1:3*nlayer): velocity at lower interface
c               rmodel(3*nlayer+1:4*nlayer): Vp/Vs ratio
c
c-------------------------------------------------------------------------
c
	subroutine forward_modelling(
     &		rmodel, moddim, ndata, nwave, fs, 
     &          predicted_data, vvs, thick_i, vp_i,  vs_i,  rho_i, 
     &		n_layers_input, modelcode)
c
        include 'hv_param.inc'


	real*4		rmodel(maxmoddim)

	real*4          predicted_data(maxdata,maxwave),
     &			thick, vs, vp, rho, vvs(maxlayer),
     &			tot_thickness, tot_thick_i,
     &			thick_i(200), vs_i(200), vp_i(200), rho_i(200),
     &			delta_thick, thick0, thick1, thick2,
     &			thick_a(20)

c
	integer		istatus, ndata(maxwave), n_layers_input
	character*20	modelcode
	character*100	scratch_path
	character*120	directory
	character*500	command
	real		bedrock_depth, ice_density, ice_vp, water_vs, water_vp, water_density

	nlayer=moddim/2


	scratch_path = "../../models/"
	directory = trim(scratch_path) // adjustl(modelcode)
	call system("mkdir " // directory)
	

c	prepare model for herrmann code
	qa=0.0
	qb=0.0
        etap=0.0
	etas=0.0
	frefp=1.0
	frefs=1.0


c	*********************************
c	Set apriori informations
	bedrock_depth = 3.51
	! Rémy, F. & Tabacco, I. E., 2000. Bedrock features and ice flow near 
	! the EPICA Ice Core Site (Dome C, Antarctica), Geophys. Res. Lett., , 27, 405–408.
	ice_density = 0.92     	
	ice_vp = 3.94  
	! Wittlinger, G. & Farra, V., 2015. Evidence of unfrozen liquids 
	! and seismic anisotropy at the base of the polar ice sheets, 
	! Polar Sci., 9(1), 66–79.

	! We assume that there is a liquid water layer at the base of the ice:
	water_vs = 0.0
	water_vp = 1.5
	water_density = 1.0

c	*********************************



	open(2,file=trim(directory)//"/model.d", form='formatted'
     &		, access='sequential', status='replace')
	write(2,fmt="(a)")'MODEL.01'
	write(2,fmt="(a)")'MODELLO DEMO'
	write(2,fmt="(a)")'ISOTROPIC'
	write(2,fmt="(a)")'KGS'
	write(2,fmt="(a)")'SPHERICAL EARTH'
	write(2,fmt="(a)")"1-D"
	write(2,fmt="(a)")'CONSTANT VELOCITY'
	write(2,fmt="(a)")'LINE08'
	write(2,fmt="(a)")'LINE09'
	write(2,fmt="(a)")'LINE10'
	write(2,fmt="(a)")'LINE11'
	write(2,fmt="(a)")'H	VP	VS	RHO	QP	QS	ETAP	ETAS	FREFP	FREFS'
	tot_thickness = 0.0

	thick_a(1) = rmodel(1)
	thick_a(2) = bedrock_depth - (thick_a(1) + rmodel(3)) 
	thick_a(3) = rmodel(3)
 	thick_a(4) = rmodel(4)
	thick_a(5) = rmodel(5)
	thick_a(6) = rmodel(6)

	do i=1,nlayer
	  	thick  = rmodel(i)
		vs=rmodel(nlayer+i)

		! We use by default Brocher relations to obtain Vp and rho given Vs. These relationships are valid
		! for rocks only. For ice and water we define Vp and Rho from a priori information.

		call brocher_Vp_Vs(vs, vp)
		call brocher_rho_Vs(vs, rho)

		! the first two layers are made of ice, then we use the empirical values for the ice
		if (i .lt. 3) then
			rho = ice_density
			vp = ice_vp 
		end if

		! The third layer is the water layer. We keep the thickness of this layer free to vary
		if (i .eq. 3) then
			rho = water_density
			vs = water_vs
			vp = water_vp
		end if

                write(2,100) thick_a(i) , ' ', vp,' ',vs,' ',rho,' ',qa,' ',qb,' ',etap,' ',
     &				etas,' ',frefp,' ',frefs
		vvs(i) = vs
		tot_thickness = tot_thickness + thick_a(i)
	end do

c	Complete the model with the lithosphere from LITHO1.0
c	If the inverte dcrust is thinner than CRUST1.0, then it completes the 
c	model with CRUST1.0 
	tot_thick_i = 0.0
	do j=1,n_layers_input
		tot_thick_i = tot_thick_i + thick_i(j)
		if (tot_thick_i .ge. tot_thickness) then
			delta_thick = tot_thick_i-tot_thickness
			if (delta_thick .gt. 0.0) then   ! add a layer to match depth the moho with LITHO1.0
				write(2,100) , delta_thick, ' ', vp_i(j),' ',vs_i(j),' ',
     &                                 rho_i(j),' ',qa,' ',qb,' ',etap,' ',
     &				       etas,' ',frefp,' ',frefs
			end if
			goto 101
		end if
	end do
	
  101	do k=j+1,n_layers_input    ! Adding LITHO1.0
		write(2,100) 	thick_i(k), ' ', vp_i(k),' ',vs_i(k),' ',
     &				rho_i(k),' ',qa,' ',qb,' ',
     &				etap,' ',etas,' ',frefp,' ',frefs
	end do
		
	close(2)


	! Computing the theoretical ellipticity from the proposed model calling Herrmann's code
	command = "python ../src/hv_subs/forward_problem.py " // directory // 
     &		"   ../../periods.txt    "  
	call system(command)

	! Reading the ellipticity curve and saving in memory
	open(16, file=trim(directory)//"/predicted_ellipticity", status="old", IOSTAT=istatus)
	do i=1,300
		read(16, *, IOSTAT = istatus) T, predicted_data(i, 1)
		if (istatus .lt. 0) goto 92
	end do
	
  92	continue
	close(16)
	
	
  100	format (f6.2,a,f5.3,a,f5.3,a,f5.3,a,f3.1,a,f3.1,a,f3.1,a,f3.1,a,f3.1,a,f3.1)

	return
	end





c--------------------- subroutines --------------------------
	subroutine brocher_rho_Vs(Vs, rho)
		call brocher_Vp_Vs(Vs, Vp)
		call brocher_rho_Vp(Vp, rho)
	return
	end

	subroutine brocher_rho_Vp(Vp, rho)
		rho = 1.6612*Vp - 0.4721*(Vp**2) + 0.0671*(Vp**3) - 0.0043*(Vp**4) + 0.000106*(Vp**5)
	return
	end

	subroutine brocher_Vp_Vs(Vs, Vp)
		Vp = 0.9409 + 2.0947*Vs - 0.8206*(Vs**2) + 0.2683*(Vs**3) - 0.0251*(Vs**4)
	return 
	end


























