#include "calplot.h"

/* program to demonstrate the use of gread 
	to overlay graphic files 
	void gread(char *fname, float X0, float Y0, float Xlow, float Ylow, 
		float Xhigh, float Yhigh, int PageNum);
*/

main()
{

	pinitf("GREAD.PLT");

		/* place PLTTST.2 image on the entire printable screen */
	gread("PLTTST.2", 0.0, 0.0, 0.0, 0.0, 10.0, 8.0, 1, 1.0, 1.0);
		/* place PLTTST.1 image in the region (0,0)->(5,5) */
	gread("PLTTST.1", 0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 1, 2.0, 2.0);

		/* fill a central region with blue */
	newpen(4);
	shader(2.65,1.65,4.85,5.85,0,0,0.01,0.01);
		/* erase a portion of that blue region, leaving a
			blue border */
	newpen(0);
	shader(2.75,1.75,4.75,5.75,0,0,0.01,0.01);
		/* reset pen for safety */
	newpen(1);
		/* take the (0,0)->(4,4) region of page 1 of PLTTST
			and offset the lower left corner to 2.75,1.75, 
			and sclae x by 0.5 and y by 1.0
			this will place everything in the region
			(2.75,1.75)->(4.75,5.75) */
	gread("PLTTST.PLT", 2.75, 1.75, 0.0, 0.0, 4.0, 4.0, 1, 0.5, 1.0);

	
		/* fill a central region with red */
	newpen(2);
	shader(3.95,2.95,6.05,5.05,0,0,0.01,0.01);
		/* set most of the interior to the background color */
	newpen(0);
	shader(4.0,3.0,6.0,5.0,0,0,0.01,0.01);
		/* reset pen for safety */
	newpen(1);
		/* take the (0,1)->(4,5) region of page 1 of PLTTST.PLT
			and offset the lower left corner of original
                        plot to (4,2) , the scale by 0.5 in x and y
			this will place everything in the region
			(4,3)->(6,5) */
	gread("PLTTST.PLT", 4.0, 3.0, 0.0, 1.0, 4.0, 5.0, 2, 0.5, 0.5);

	pend();
}
