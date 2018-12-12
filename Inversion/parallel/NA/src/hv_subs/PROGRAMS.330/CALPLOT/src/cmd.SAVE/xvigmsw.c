/* interface program to connect XviG calls to
	GRAPHAPP CALLS

	Robert B. Herrmann
	Computer Programs in Seismology
	09 SEP 2001
*/


#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<ctype.h>


/* CALPLOT parameters */
#define MaXCOL 35	/* note I need this define put it
				must be identical to the one in
				plotxvig.c
			*/
#ifndef INT
#ifdef MSDOS
#define INT long
#else
#define INT int
#endif
#endif

#define ON      1
#define OFF     0

INT dc_ClipRight, dc_ClipTop, dc_ClipBottom, dc_ClipLeft;
INT dc_ClipRight_sv, dc_ClipTop_sv, dc_ClipBottom_sv, dc_ClipLeft_sv;
extern int Ldc_top, Ldc_right, Ldc_bottom, Ldc_left; /* for internal clipping */
extern INT dc_left     ;    /* device plotting limits */
extern INT dc_right    ;
extern INT dc_bottom   ;
extern INT dc_top      ;




/* XviG function definitions */

long XviG_CloseCursor(void);
 int XviG_CloseWindow(char *name);
void XviG_DrawLine(int x1, int y1, int x2, int y2);
void XviG_DrawPoint(int x, int y);
void XviG_Exit(void);
void XviG_FillPolygon(int *coords, int npoints);
void XviG_FillRectangle(int x1, int y1, int x2, int y2);
void XviG_Flush(void);
int XviG_GetChar(void);
int XviG_GetCursor(int type, int *x_pos, int *y_pos); 
int XviG_Init(char *classname, int color_array[][3], int nr_of_colors);
int XviG_OpenWindow(char *name, int x, int y, unsigned int *width, unsigned int *height);
void XviG_SendMessage(int type, int i1, int i2, int i3, int i4);
void XviG_SetColor(int nr);
void XviG_SetCursor(long cursor);
void XviG_SetGC(int gc);
void XviG_WindowPosition(int *x, int *y);
void XviG_WindowSize(unsigned int *width, unsigned int *height);

char Classname[100];

/* XviG defines */
#define XviG_BUTTON1  -1
#define XviG_BUTTON2  -2
#define XviG_BUTTON3  -3

#define XviG_CURSOR_ARROW	 0L
#define XviG_CURSOR_XORARROW	-1L
#define XviG_CURSOR_XHAIR 	-2L
#define XviG_CURSOR_PLUS  	-3L
#define XviG_CURSOR_BOX  	-4L
#define XviG_CURSOR_RUBBER  	-5L
#define XviG_CURSOR_OFF  	-6L
#define XviG_CURSOR_HYPERBOLA  	-7L

/* GRAPHAPP interface */
/* GRAPHAPP */

#include <graphapp.h>
static rgb  logval[MaXCOL];     /* index into MS API colormap */
bitmap u;	/* bitmap for backup store */
window w, wtemp;
int MYinput (char *title, char *text, int l);
static point *Pt;
static point pt;
void mouse_click( drawing w, int buttons, point p);
void mouse_move( drawing w, int buttons, point p);
void win_redraw(drawing w, rect p);
void  win_getkey(drawing w, int key);
void  win_actiongetkey(drawing w, int key);
static int mouse_button;
static int key_char;
static int actionkey_char;
static int contloop = ON;
void Mainloop(void);
static int border, title, dcClipLeft, dcClipBottom, dcClipRight, dcClipTop;
static int width, height;

static int page_dirty = 0;

static void draw_cursor(int curstyp, int xc, int yc, 
		int xl, int yl, int xh, int yh);
static int inclipregion(int xpos, int ypos);


#define True 1
#define False 0
static int curstype = 0;
static int old_curstype = 0;
static int cur_x = -1;	/* coordinates of last call to di_cross */
static int cur_y = -1;
static int mouse_x, mouse_y;
static int user_curstype;
static int xhair_color = 0;
static int cursor_type;
static int  cursor_on, cursor_drawn;
static  int prev_cursor_x, prev_cursor_y;
static unsigned int xhair_width, xhair_height;
static int do_event_loop = True;

/*------------------------------------------------------------------------------
-- Some general macro definitions
------------------------------------------------------------------------------*/

#define ABS(n)    ((n) < 0 ? -(n) : (n))
#define MaX(a,b)  ((a) > (b) ? (a) : (b))
#define MiN(a,b)  ((a) < (b) ? (a) : (b))






/* implementation */
/*	Notes:
	For the event mechanism we will always need a loop
	the functions will do nothing but pass values?
	The Mainloop will have a flag top determine when to
	stop reacting -- Mainloop must always be alive??

	perhaps contloop should have a set of different values?
*/

long XviG_CloseCursor(void)
{
}

int XviG_CloseWindow(char *name)
{
	return(0);
}

void XviG_DrawLine(int xx1, int yy1, int xx2, int yy2)
{
	drawto(u);
	drawline(pt(xx1,yy1), pt(xx2,yy2));
	page_dirty = 1;
}

void XviG_DrawPoint(int x, int y)
{
	drawto(u);
	drawpoint(pt(x,y));
	page_dirty = 1;
}

void XviG_Exit(void)
{
	exitapp();
}

void XviG_FillPolygon(int *coords, int npoints)
{
	int i, j;
	if((Pt = (point  *)calloc(npoints, sizeof(point))) == NULL)
		return;
	for(i=0, j=0; i < npoints; i++, j+=2){
	Pt[i] = pt( coords[j], coords[j+1]) ;
	}
	drawto(u);
	fillpolygon(Pt, npoints);
	free(Pt);
	page_dirty = 1;

}


void XviG_FillRectangle(int x1, int y1, int x2, int y2)
{     

	setdrawmode(S);
	drawto(u);
	fillrect(rect(x1, y1, x2-x1, y2-y1));
	page_dirty = 1;

}

void XviG_Flush(void)
{
	if(page_dirty)win_redraw(w, rect(0,0,width,height)); /* refresh */
	page_dirty = 0;
}

int XviG_GetChar(void)
{
	int action;
	contloop = ON;
	if(page_dirty)win_redraw(w, rect(0,0,width,height)); /* refresh */
	Mainloop();
	/* map key board to ASCII note we shoudl permit a CTRL L and CTRL M
		to get 12 for NL and 13 for CR
		in graphapp we would have to use setkeyaction to note the
		flag if set ?? */
	
	action = getkeystate();
	page_dirty = 0;
	if(action == CtrlKey){
		if( key_char == 'L' || key_char == 'l')
			return(12);
		if( key_char == 'M' || key_char == 'm')
			return(13);
	}
	if(key_char == 10)
		return(12);
	else
		return(key_char);

}

int XviG_GetCursor(int type, int *x_pos, int *y_pos) 
{
	contloop =  ON;
	/* update the current window */
	if(page_dirty)win_redraw(w, rect(0,0,width,height)); /* refresh */
	/* arguments are needed 
			but not actually used */
	/* initialize key */
	key_char = 0 ;
	setdrawmode(DxorS);
	setdrawmode(notD);
	XviG_SetColor(xhair_color);
	drawto(w);
	Mainloop();
	*x_pos = mouse_x;
	*y_pos = mouse_y;
	cur_x = mouse_x;
	cur_y = mouse_y;
	page_dirty = 0;
	cursor_on = False;
	setdrawmode(S);
	drawto(u);
	if(isprint(key_char))
		return(key_char);
	else
		return(mouse_button);
}

int XviG_Init(char *classname, int color_array[][3], int nr_of_colors)
{
	int i;
	int depth;

	/* define the color values */
	for(i=0; i < nr_of_colors ; i++)
		logval[i] = rgb(color_array[i][0], 
				color_array[i][1],
				color_array[i][2]);

	/* start graphapp by creating a window */

	wtemp = newwindow(classname,
		rect(dc_left,dc_top -dc_top,dc_right-dc_left,dc_top-dc_bottom),
		StandardWindow);
	strcpy(Classname, classname);
	depth = getdepth(wtemp);
	hide(wtemp);
	
	/* RBH note that the return = 0 if cannot create the window 
		check this */
	return(depth);
}

int XviG_OpenWindow(char *name, int x, int y, unsigned int *Width, unsigned int *Height)
{
	int depth;
	/* used to redefine the window height */
/* graphapp code does not do a resize 
	so kludge this by opening the window */
/*
	
	resize(w, rect(0, 0, *width, *height));
	resize(u, rect(0, 0, *width, *height));
*/
	/* start graphapp by creating a window */

	width  = *Width;
	height = *Height;

	w = newwindow(Classname,
		rect(0, 0, width, height), 
		StandardWindow);
	depth = getdepth(w);
	drawto(w);

	u = newbitmap(width, height, depth);

	/* define the event mechanisms */
	setmousedown(w, mouse_click);
	setmousemove(w, mouse_move);
	setredraw(w, win_redraw);
	setkeydown(w, win_getkey);

	setdrawmode(S);
	setbackground(w,logval[0]) ;
	show(w);
	setdrawmode(S);
	drawto(u);	/* ? */
}

void XviG_SendMessage(int type, int i1, int i2, int i3, int i4)
{
/* this is used to define the region for the user cursor implementation
as in
XviG_SendMessage( 1, border, title, 0, 0);
XviG_SendMessage( 2, dc_ClipLeft, dc_ClipBottom, dc_ClipRight, dc_ClipTop);
*/
	if(type == 1){
		border = i1;
		title  = i2;
	}  else if(type == 2){
		dcClipLeft   = i1; 
		dcClipBottom = i2; 
		dcClipRight  = i3; 
		dcClipTop    = i4;
	}
}


void XviG_SetColor(int nr)
{
	setcolor(logval[nr] );
}

void XviG_SetCursor(long cursor)
{
	int i;
	/* these are all placeholders now */
	if(cursor == XviG_CURSOR_ARROW){
		user_curstype = 0;
		/* define hw cursor */
		old_curstype = user_curstype;
		setcursor(ArrowCursor);
	} else if(cursor == XviG_CURSOR_XORARROW){
		user_curstype = 1;
		/* set Xhair color, define hw cursor blank */
		old_curstype = user_curstype;
	} else if(cursor == XviG_CURSOR_XHAIR){
		user_curstype = 2;
		/* set Xhair color, define hw cursor blank */
		old_curstype = user_curstype;
	} else if(cursor == XviG_CURSOR_PLUS){
		user_curstype = 3;
		/* set Xhair color, define hw cursor blank */
		old_curstype = user_curstype;
	} else if(cursor == XviG_CURSOR_BOX){
		user_curstype = 4;
		/* set Xhair color, define hw cursor blank */
		old_curstype = user_curstype;
	} else if(cursor == XviG_CURSOR_RUBBER){
		user_curstype = 5;
		/* set Xhair color, define hw cursor blank */
		old_curstype = user_curstype;
	} else if(cursor == XviG_CURSOR_OFF){
		user_curstype = old_curstype;
		cur_x = -1;
		cur_y = -1;
	}
	xhair_color = 0;
	if(cursor)cursor_on = True;
}

void XviG_SetGC(int gc)
{
	if(gc == 0)
		setdrawmode(S);
	else if(gc == 1)
		setdrawmode(DxorS);
}

void XviG_WindowPosition(int *x, int *y)
{
	*x = 0;
	*y = 0;
}

void XviG_WindowSize(unsigned int *width, unsigned int *height)
{
/*
	*width = dc_right -dc_left;
	*height = dc_top -dc_bottom;
*/
}

/* interface cursor stuff */
struct curs {
	int yv;
	int xlv;
	int xhv;
} ;
struct curs curs_arrow[] = {
	{  0, 0, 0 },
	{  1, 0, 1 },
	{  2, 0, 2 },
	{  3, 0, 3 },
	{  4, 0, 4 },
	{  5, 0, 5 },
	{  6, 0, 6 },
	{  7, 0, 7 },
	{  8, 0, 4 },
	{  9, 0, 1 },
	{ 10, 0, 0 },
	{  9, 3, 4 },
	{ 10, 4, 5 },
	{ 11, 4, 5 },
	{ 12, 5, 5 }
};

static void draw_cursor(int curstyp, int xc, int yc, 
		int xl, int yl, int xh, int yh)
{
/* (xc,yc) is the center point
	The xl,yl,xh,yh are screen coordinates of clip region
	curstyp	0	arrow
		1	arrow
		2	cross hair
		3	plus
		4	rubber box  mode
		5	rubber band mode
		6	rubber mode off 
		7	User pixmap
*/
	int nl, i;
	int xp,xm,yp,ym;
	nl = sizeof(curs_arrow)/sizeof(struct curs);
	if(curstyp == 0){
		setcursor(ArrowCursor);
/*
              XDefineCursor(display, window, arrow_cursor);
*/
	} else if (curstyp == 1) {
		for(i=0;i<nl;i++)
                	drawline(
				pt(xc+curs_arrow[i].xlv,
				yc+curs_arrow[i].yv),
				pt(xc+curs_arrow[i].xhv,
				yc+curs_arrow[i].yv));

	} else if (curstyp == 2) {
                drawline(
			pt(xc, yl), pt(xc, yh));
                drawline(
			pt(xl, yc), pt(xh, yc));
	} else if (curstyp == 3) {
		xp = MiN((xc+10),xh);
		xm = MaX((xc-10),xl);
		yp = MiN((yc-10),yh);
		ym = MaX((yc+10),yl);
                drawline(
			pt(xc, ym), pt(xc, yp));
                drawline(
			pt(xm, yc), pt(xp, yc));
	} else if (curstyp == 4) {
                drawline(
			pt(cur_x, cur_y), pt(xc, cur_y));
                drawline(
			pt(xc, cur_y), pt(xc, yc));
                drawline(
			pt(xc, yc), pt(cur_x, yc));
                drawline(
			pt(cur_x, yc), pt(cur_x, cur_y));
	} else if (curstyp == 5) {
                drawline(
			pt(cur_x, cur_y), pt(xc, yc));
	}
	cursor_type = 1;
	
}

static int inclipregion(int xpos, int ypos)
{
	/* the clip region is in user coordinates */
	xpos = xpos - border;
	ypos = Ldc_top - border - ypos;
	if(ypos < dc_bottom)
		return 0;
	else {
		if(xpos <= dc_ClipRight && xpos >= dc_ClipLeft
			&& ypos <=dc_ClipTop && ypos >= dc_ClipBottom)
				return 1;
		else
				return 0;
	}
}

/* --------------------------GRAPHAPP CODE---------------------------------------*/

/* event mechanism for GRAPHAPP */

void mouse_click( drawing w, int buttons, point pp)
{
	contloop = OFF;
/* MSW   LEFT = 1
	CENTER (if permitted) = 2
	RIGHT = 4
	LEFT + RIGHT = 5
	CALPLOT assumes 1 3 2
                    for L C R
check this 09 09 2001 with xvig source code
*/
	if(buttons == 1)
		mouse_button = 1;
	else if(buttons == 2)
		mouse_button = 3;
	else if(buttons == 4)
		mouse_button = 2;
	else if(buttons == 5)
		mouse_button = 3;
	if (cursor_on) {
		if(cursor_drawn){
			drawto(w);
			setdrawmode(DxorS);
			setdrawmode(notD);
			XviG_SetColor(xhair_color);
			if (cursor_drawn) {
			draw_cursor( curstype,  prev_cursor_x,  prev_cursor_y, 
				border+dc_ClipLeft, Ldc_top-border-dc_ClipTop, 
				border+dc_ClipRight, Ldc_top-border-dc_ClipBottom);
			}
		}
		cursor_drawn = False;
	}
	cursor_on = False;
	/* do a bounds check */
	prev_cursor_x = pp.x;
	prev_cursor_y = pp.y;
	mouse_x = pp.x;
	mouse_y = pp.y;
	if(prev_cursor_x < dc_left+border)
		prev_cursor_x = dc_left+border;
	if(prev_cursor_x > dc_right+border)
		prev_cursor_x = dc_right+border;
	if(prev_cursor_y > Ldc_top - border )
		prev_cursor_y = Ldc_top - border ;
	if(prev_cursor_y < border )
		prev_cursor_y = border;
}

void mouse_move( drawing w, int buttons, point pp)
{
	mouse_x = pp.x;
	mouse_y = pp.y;

	if (cursor_on) {
		drawto(w);
		setdrawmode(DxorS);
		setdrawmode(notD);
		XviG_SetColor(xhair_color);
		if (cursor_drawn) {
		draw_cursor( curstype,  prev_cursor_x,  prev_cursor_y, 
			border+dc_ClipLeft, Ldc_top-border-dc_ClipTop, 
			border+dc_ClipRight, Ldc_top-border-dc_ClipBottom);
		}
		prev_cursor_x = pp.x;
		prev_cursor_y = pp.y;
		if(prev_cursor_x < dc_left+border)
			prev_cursor_x = dc_left+border;
		if(prev_cursor_x > dc_right+border)
			prev_cursor_x = dc_right+border;
		if(prev_cursor_y > Ldc_top - border )
			prev_cursor_y = Ldc_top - border ;
		if(prev_cursor_y < border )
			prev_cursor_y = border;
		if(inclipregion(prev_cursor_x,prev_cursor_y)){
			setcursor(BlankCursor);
/*
			XDefineCursor(display, window, empty_cursor);
*/
			curstype = user_curstype;
		} else {
			curstype = 0;
		}
		draw_cursor( curstype,  prev_cursor_x,  prev_cursor_y, 
			border+dc_ClipLeft, Ldc_top-border-dc_ClipTop, 
			border+dc_ClipRight, Ldc_top-border-dc_ClipBottom);
		cursor_drawn = True;
}
}

void Mainloop(void)
{
	while (doevent() && contloop)
		continue;
}

void win_redraw(drawing w, rect p)
{
	bitblt(w,u,pt(0,0), rect(0,0,width,height),S);
	show(w);
	drawto(u);
	
}

void win_getkey(drawing w, int key){
	contloop = OFF;
/* RBH note check mappings into ASCII , e.g., 
	DEL BS to BS
*/
	key_char = key;
	if (cursor_drawn) {
		draw_cursor( curstype,  prev_cursor_x,  prev_cursor_y, 
			border+dc_ClipLeft, Ldc_top-border-dc_ClipTop, 
			border+dc_ClipRight, Ldc_top-border-dc_ClipBottom);
		cursor_drawn = False;
	}
}

void win_actiongetkey(drawing w, int key){
	actionkey_char = key;
}

