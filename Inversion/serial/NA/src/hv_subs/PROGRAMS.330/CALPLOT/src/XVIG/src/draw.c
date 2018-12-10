/* File>>> draw.c
--
-- %M% -- version %I% (IMEC)            last updated: %E%
--
-- Copyright (c) 1993
-- IMEC vzw
-- Kapeldreef 75
-- B-3001 LEUVEN
-- BELGIUM
--
-- Author   : A. Demaree
--
-- Date     : February 1, 1993
--
-- Function :
--
-- Comment  :
--
-- Review   :
--
	revision history

 *	23 JUL 2003 Following Email from Meijian An, 
 *	Department of Geophysics, Institute of Astronomics, 
 *	Geophysical and Atmospheric Sciences, *	University of Sao Paulo
 *	Sao Paulo, Brazil,  carefully free memory
*/


/*------------------------------------------------------------------------------
-- Global include files
------------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <X11/Xlib.h>

/*------------------------------------------------------------------------------
-- Local include files
------------------------------------------------------------------------------*/

#include "xviglocal.h"

/*------------------------------------------------------------------------------
-- Static variable declarations
------------------------------------------------------------------------------*/

static XPoint *points = (XPoint *) NULL;
static int nr_of_points = 0;

/*------------------------------------------------------------------------------
-- Local function declarations
------------------------------------------------------------------------------*/
void XviG_SetGC(int gc)
{
	if(gc == ON){
		XSetFunction(XviG_display, XviG_gc, GXxor);
		XSetFunction(XviG_display, XviG_gc_fill, GXxor);
		 RBH_xoring_state = ON;
	} else if (gc == OFF){
		XSetFunction(XviG_display, XviG_gc, GXcopy);
		XSetFunction(XviG_display, XviG_gc_fill, GXcopy);
		 RBH_xoring_state = OFF;
	}
}


/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_DrawPoint(int x,
                    int y)
{
  XDrawPoint(XviG_display, XviG_window, XviG_gc, x, y);
  XDrawPoint(XviG_display, XviG_pixmap, XviG_gc, x, y);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_DrawLine(int x1,
                   int y1,
                   int x2,
                   int y2)
{
/*
fprintf(stderr,"Line: (%d,%d) -> (%d,%d)\n",x1,y1,x2,y2);
*/
  XDrawLine(XviG_display, XviG_window, XviG_gc, x1, y1, x2, y2);
  XDrawLine(XviG_display, XviG_pixmap, XviG_gc, x1, y1, x2, y2);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_DrawPolyLine(int *coords,
                       int npoints)
{
  int i, j;

  if (!points)
  {
    points = (XPoint *) malloc(npoints * sizeof(XPoint));
    nr_of_points = npoints;
  }
  else
  if (npoints > nr_of_points)
  {
    points = (XPoint *) realloc((void *) points, npoints * sizeof(XPoint));
    nr_of_points = npoints;
  }

  for (i = 0, j = 0; j < npoints; j++)
  {
    points[j].x = coords[i++];
    points[j].y = coords[i++];
  }

  XDrawLines(XviG_display, XviG_window, XviG_gc,
             points, npoints, CoordModeOrigin);
  XDrawLines(XviG_display, XviG_pixmap, XviG_gc,
             points, npoints, CoordModeOrigin);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_DrawRectangle(int x1,
                        int y1,
                        int x2,
                        int y2)
{
  XDrawRectangle(XviG_display, XviG_window, XviG_gc,
                 Min(x1, x2), Min(y1, y2), Abs(x2-x1), Abs(y2-y1));
  XDrawRectangle(XviG_display, XviG_pixmap, XviG_gc,
                 Min(x1, x2), Min(y1, y2), Abs(x2-x1), Abs(y2-y1));
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_DrawPolygon(int *coords,
                      int npoints)
{
  int i, j, npts_new;

  npts_new = ((coords[0] != coords[(npoints-1)*2]) ||
              (coords[1] != coords[(npoints*2)-1])) ? npoints + 1 : npoints;

  if (!points)
  {
    points = (XPoint *) malloc(npts_new * sizeof(XPoint));
    nr_of_points = npts_new;
  }
  else
  if (npts_new > nr_of_points)
  {
    points = (XPoint *) realloc((void *) points, npts_new * sizeof(XPoint));
    nr_of_points = npts_new;
  }

  for (i = 0, j = 0; j < npoints; j++)
  {
    points[j].x = coords[i++];
    points[j].y = coords[i++];
  }

  if (npoints != npts_new)
  {
    points[j].x = coords[0];
    points[j].y = coords[1];
  }

  XDrawLines(XviG_display, XviG_window, XviG_gc,
             points, npts_new, CoordModeOrigin);
  XDrawLines(XviG_display, XviG_pixmap, XviG_gc,
             points, npts_new, CoordModeOrigin);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_DrawArc(int x,
                  int y,
                  unsigned int radius1,
                  unsigned int radius2,
                  int angle1,
                  int angle2)
{
  XDrawArc(XviG_display, XviG_window, XviG_gc,
           x-radius1, y-radius2, 2*radius1, 2*radius2, angle1*64, angle2*64);
  XDrawArc(XviG_display, XviG_pixmap, XviG_gc,
           x-radius1, y-radius2, 2*radius1, 2*radius2, angle1*64, angle2*64);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_FillRectangle(int x1,
                        int y1,
                        int x2,
                        int y2)
{
/*
fprintf(stderr,"FILL: (%d,%d) -> (%d,%d)\n",x1,y1,x2,y2);
*/
  XFillRectangle(XviG_display, XviG_window, XviG_gc_fill,
                 Min(x1, x2), Min(y1, y2), Abs(x2-x1), Abs(y2-y1));
  XFillRectangle(XviG_display, XviG_pixmap, XviG_gc_fill,
                 Min(x1, x2), Min(y1, y2), Abs(x2-x1), Abs(y2-y1));
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_FillPolygon(int *coords,
                      int npoints)
{
  int i, j;

  if (!points)
  {
    points = (XPoint *) malloc(npoints * sizeof(XPoint));
    nr_of_points = npoints;
  }
  else
  if (npoints > nr_of_points)
  {
    points = (XPoint *) realloc((void *) points, npoints * sizeof(XPoint));
    nr_of_points = npoints;
  }

  for (i = 0, j = 0; j < npoints; j++)
  {
    points[j].x = coords[i++];
    points[j].y = coords[i++];
  }

  XFillPolygon(XviG_display, XviG_window, XviG_gc_fill,
               points, npoints, Complex, CoordModeOrigin);
  XFillPolygon(XviG_display, XviG_pixmap, XviG_gc_fill,
               points, npoints, Complex, CoordModeOrigin);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_FillArc(int x,
                  int y,
                  unsigned int radius1,
                  unsigned int radius2,
                  int angle1,
                  int angle2)
{
  XFillArc(XviG_display, XviG_window, XviG_gc_fill,
           x-radius1, y-radius2, 2*radius1, 2*radius2, angle1*64, angle2*64);
  XFillArc(XviG_display, XviG_pixmap, XviG_gc_fill,
           x-radius1, y-radius2, 2*radius1, 2*radius2, angle1*64, angle2*64);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_Flush(void)
{
  XFlush(XviG_display);
}
