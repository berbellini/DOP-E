#include <stdio.h>
#include <stdlib.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

main()
{
	visual();
}

visual()
{
	Display *display;
	Colormap colormap;
	Window window;
	int screen;

	XVisualInfo vTemplate;
	XVisualInfo *visualList;
	int visualsMatched;

	int x ;
	int y ;
	unsigned int width , height , border_width ;
	unsigned long mask;
	int depth;
	XSetWindowAttributes attributes;
	x = 0; y = 0 ;
	width = 100; height = 100; border_width = 5;

/* Vol I page 210 - 211 */


	display = XOpenDisplay(NULL);
	screen = DefaultScreen(display);
	depth = DisplayPlanes(display,screen);
	vTemplate.screen = screen;
	vTemplate.depth = depth;
	fprintf(stderr,"depth = %d\n",depth);
	visualList = XGetVisualInfo(display, 
		VisualScreenMask |
		VisualDepthMask  
			, &vTemplate, &visualsMatched);
	fprintf(stderr,"screen	%d\n",visualList[0].screen);
	fprintf(stderr,"depth 	%d\n",visualList[0].depth );
	fprintf(stderr,"class 	%d\n",visualList[0].class );
	fprintf(stderr,"     where \n");
	fprintf(stderr,"     StaticGray    0\n");
	fprintf(stderr,"     GrayScale     1\n");
	fprintf(stderr,"     StaticColor   2\n");
	fprintf(stderr,"     PseudoColor   3\n");
	fprintf(stderr,"     TrueColor     4\n");
	fprintf(stderr,"     DirectColor   5\n");

	fprintf(stderr,"colormap_size 	%d\n",visualList[0].colormap_size );
	fprintf(stderr,"bits_per_rgb 	%d\n",visualList[0].bits_per_rgb );
	if(visualsMatched == 0 )
		fatalError("No matching visuals\n");

	colormap = XCreateColormap (display, RootWindow(display,screen),
		visualList[0].visual, AllocNone);
	window = XCreateWindow (display, RootWindow(display,screen),
		x, y, width, height, border_width, vTemplate.depth,
		InputOutput, visualList[0].visual, mask, &attributes);
	XSetWindowColormap(display,window, colormap);

	XFree(visualList);
}

fatalError(char *str)
{
	fprintf(stderr,"%s",str);
	exit(0);
}
