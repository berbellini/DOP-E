C  program to demonstrate the use of gread 
C	to overlay graphic files 
C	void gread(char *fname, float X0, float Y0, float Xlow, float Ylow, 
C	float Xhigh, float Yhigh, int PageNum);

	call pinitf('GREAD.PLT')


c-----
c	place PLTTST.2 image on the entire printable screen 
c-----
	call gread('PLTTST.2', 0.0, 0.0, 0.0, 0.0, 10.0, 8.0,1,1.0,1.0)
c-----
c	 place PLTTST.1 image in the region (0,0)->(5,5) 
c-----
	call gread('PLTTST.1', 0.0, 0.0, 0.0, 0.0, 5.0, 5.0,1,2.0,2.0)

c-----
c	fill a central region with blue 
c-----
	call newpen(4)
	call shader(2.65,1.65,4.85,5.85,0,0,0.01,0.01)
	call newpen(0)
c-----
c	 erase a portion of that blue region, leaving a blue border 
c-----
	call shader(2.75,1.75,4.75,5.75,0,0,0.01,0.01)
c-----
c	 reset pen for safety 
c-----
	call newpen(1)
c-----
c	take the (0,0)->(4,4) region of page 1 of PLTTST
c		and offset the lower left corner to 2.75,1.75, 
c		and scale x by 1.5 and y by 0.5
c		this will place everything in the region
c		(2.75,1.75)->(4.75,5.75) 
c-----
	call gread('PLTTST.PLT', 2.75,1.75,0.0,0.0,4.0,4.0,1,0.5,1.0)

c-----
c	fill a central region with red 
c-----
	call newpen(2)
	call shader(3.95,2.95,6.05,5.05,0,0,0.01,0.01)
c-----
c	 set most of the interior to the background color 
c-----
	call newpen(0)
	call shader(4.0,3.0,6.0,5.0,0,0,0.01,0.01)
c-----
c	 reset pen for safety 
c-----
	call newpen(1)
c-----
c	take the (0,1)->(4,5) region of page 1 of PLTTST
c		and offset the lower left corner of original plot to (4,2) 
c		and scale x by 0.5 and y by 0.5
c		this will place everything in the region (4,3)->(6,5) 
c-----
	call gread('PLTTST.PLT', 4.0, 3.0, 0.0, 1.0, 4.0,5.0,2,0.5,0.5)

	call pend()
	end

	
