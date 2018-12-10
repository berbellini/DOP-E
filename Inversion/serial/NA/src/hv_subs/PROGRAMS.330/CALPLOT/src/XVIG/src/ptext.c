/* File>>> ptext.c
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
-- Date     : June 23, 1993
--
-- Function : Routine to draw text by means of polylines.
--
-- Comment  :
--
-- Review   :
--
*/


/*------------------------------------------------------------------------------
-- Global include files
------------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <X11/Xlib.h>

/*------------------------------------------------------------------------------
-- Local include files
------------------------------------------------------------------------------*/

#include "xviglocal.h"
#include "charcodes.h"

/*------------------------------------------------------------------------------
-- Some macro definitions
------------------------------------------------------------------------------*/

#define PI  ((double)3.14159265358979323846)
#define ROUND(x)  ((int)floor((x) + 0.5))

/*------------------------------------------------------------------------------
-- Static variable declarations
------------------------------------------------------------------------------*/

static double t_text[3][2];
static XPoint points[16];

/*------------------------------------------------------------------------------
-- Local function declarations
------------------------------------------------------------------------------*/

static void get_charmatrix(int x_offset, int y_offset,
                           unsigned int width, unsigned int height,
                           int rotation);
static void transchar(int x, int y,
                      int *xt, int *yt);
static void drawchar(int c, int offset);

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_PolyText(char *contents,
                   int x,
                   int y,
                   unsigned int width,
                   unsigned int height,
                   int rotation)
{
  int len, newx, newy, offset, i;

  if ((len = (int) strlen(contents)) == 0)
    return;

  rotation %= 360;
  if (rotation < 0)
    rotation += 360;

  get_charmatrix(x, y, width, height, rotation);

  if ((rotation > 90) && (rotation <= 270))
  {
    transchar(7*len, 4, &newx, &newy);
    t_text[0][0] = newx * 1.0;
    t_text[0][1] = newy * 1.0;
    t_text[1][0] = -t_text[1][0];
    t_text[1][1] = -t_text[1][1];
    t_text[2][0] = -t_text[2][0];
    t_text[2][1] = -t_text[2][1];
  }

  XSetLineAttributes(XviG_display, XviG_gc,
                     XviG_save_linewidth, XviG_save_linestyle,
                     CapRound, JoinRound);

  for (i = 0, offset = 0; i < len; i++, offset += 7)
    drawchar(contents[i], offset);

  XSetLineAttributes(XviG_display, XviG_gc,
                     XviG_save_linewidth, XviG_save_linestyle,
                     CapButt, JoinMiter);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

static void get_charmatrix(int x_offset, int y_offset,
                           unsigned int width, unsigned int height,
                           int rotation)
{
  double rcos,rsin;

  if (rotation == 0)
  {
    rcos = 1.0;
    rsin = 0.0;
  }
  else
  if (rotation == 90)
  {
    rcos = 0.0;
    rsin = 1.0;
  }
  else
  if (rotation == 180)
  {
    rcos = -1.0;
    rsin = 0.0;
  }
  else
  if (rotation == 270)
  {
    rcos = 0.0;
    rsin = -1.0;
  }
  else
  {
    rcos = cos((rotation*PI)/180.0);
    rsin = sin((rotation*PI)/180.0);
  }

  t_text[0][0] = x_offset * 1.0;
  t_text[0][1] = y_offset * 1.0;
  t_text[1][0] =  rcos * (width/7.0);
  t_text[1][1] =  rsin * (width/7.0);
  t_text[2][0] = -rsin * (height/8.0);
  t_text[2][1] =  rcos * (height/8.0);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

static void transchar(int x, int y,
                      int *xt, int *yt)
{
  *xt = ROUND(t_text[0][0] + t_text[1][0]*x + t_text[2][0]*y);
  *yt = ROUND(t_text[0][1] - t_text[1][1]*x - t_text[2][1]*y);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

static void drawchar(int c, int offset)
{
  int npolys, npoints, charpos, polypos, pointpos, xt, yt;

  charpos = charpointer[((c < 32) || (c > 126)) ? 0 : c - 32];

  for (npolys = 0, polypos = charpos + 1;
       npolys < charcode[charpos];
       npolys++, polypos = pointpos)
  {
    for (npoints = 0, pointpos = polypos + 1;
         npoints < charcode[polypos];
         npoints++, pointpos += 2)
    {
      transchar(charcode[pointpos]+offset, charcode[pointpos+1], &xt, &yt);

      points[npoints].x = xt;
      points[npoints].y = yt;
    }

    if (charcode[polypos] == 2)
    {
      XDrawLine(XviG_display, XviG_window, XviG_gc,
                (int) points[0].x, (int) points[0].y,
                (int) points[1].x, (int) points[1].y);
      XDrawLine(XviG_display, XviG_pixmap, XviG_gc,
                (int) points[0].x, (int) points[0].y,
                (int) points[1].x, (int) points[1].y);
    }
    else
    {
      XDrawLines(XviG_display, XviG_window, XviG_gc,
                 points, charcode[polypos], CoordModeOrigin);
      XDrawLines(XviG_display, XviG_pixmap, XviG_gc,
                 points, charcode[polypos], CoordModeOrigin);
    }
  }
}
