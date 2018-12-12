/* File>>> pattern.c
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
#include "pattern.h"
#include "dither.h"

/*------------------------------------------------------------------------------
-- Static variable declarations
------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------
-- Local function declarations
------------------------------------------------------------------------------*/

static void Create_Fillpattern(int nr,
                               char *bits,
                               unsigned int width,
                               unsigned int height);

static void Create_Filldither(int nr,
                               char *bits,
                               unsigned int width,
                               unsigned int height);

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_SetLineStyle(int nr,
                       unsigned int width)
{
  int local_nr;

  if ((local_nr = nr % XviG_NR_OF_LINESTYLES) == 0)
    XSetLineAttributes(XviG_display, XviG_gc,
                       width, LineSolid, CapButt, JoinMiter);
  else
  {
    local_nr = Abs(local_nr);
    XSetLineAttributes(XviG_display, XviG_gc,
                       width, LineOnOffDash, CapButt, JoinMiter);
    XSetDashes(XviG_display, XviG_gc,
               0, linestyle[local_nr], linestyle_length[local_nr]);
  }

  /*
  -- Save the linewidth and linestyle to be used in the XviG_PolyText call
  */

  XviG_save_linewidth = width;
  XviG_save_linestyle = (local_nr == 0) ? LineSolid : LineOnOffDash;
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_SetFillStyle(int nr)
{
  int local_nr;

  if ((local_nr = nr % XviG_NR_OF_FILLPATTERNS) == 0)
    XSetFillStyle(XviG_display, XviG_gc_fill, FillSolid);
  else
  {
    local_nr = Abs(local_nr);
    XSetFillStyle(XviG_display, XviG_gc_fill, FillOpaqueStippled);
    XSetStipple(XviG_display, XviG_gc_fill, fillpattern[local_nr]);
  }
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_CreateFillpatterns(void)
{
  Create_Fillpattern( 0, (char *) fill00_bits, fill00_width, fill00_height);
  Create_Fillpattern( 1, (char *) fill01_bits, fill01_width, fill01_height);
  Create_Fillpattern( 2, (char *) fill02_bits, fill02_width, fill02_height);
  Create_Fillpattern( 3, (char *) fill03_bits, fill03_width, fill03_height);
  Create_Fillpattern( 4, (char *) fill04_bits, fill04_width, fill04_height);
  Create_Fillpattern( 5, (char *) fill05_bits, fill05_width, fill05_height);
  Create_Fillpattern( 6, (char *) fill06_bits, fill06_width, fill06_height);
  Create_Fillpattern( 7, (char *) fill07_bits, fill07_width, fill07_height);
  Create_Fillpattern( 8, (char *) fill08_bits, fill08_width, fill08_height);
  Create_Fillpattern( 9, (char *) fill09_bits, fill09_width, fill09_height);
  Create_Fillpattern(10, (char *) fill10_bits, fill10_width, fill10_height);
  Create_Fillpattern(11, (char *) fill11_bits, fill11_width, fill11_height);
  Create_Fillpattern(12, (char *) fill12_bits, fill12_width, fill12_height);
  Create_Fillpattern(13, (char *) fill13_bits, fill13_width, fill13_height);
  Create_Fillpattern(14, (char *) fill14_bits, fill14_width, fill14_height);
  Create_Fillpattern(15, (char *) fill15_bits, fill15_width, fill15_height);
  Create_Fillpattern(16, (char *) fill16_bits, fill16_width, fill16_height);
  Create_Fillpattern(17, (char *) fill17_bits, fill17_width, fill17_height);
  Create_Fillpattern(18, (char *) fill18_bits, fill18_width, fill18_height);
  Create_Fillpattern(19, (char *) fill19_bits, fill19_width, fill19_height);
  Create_Fillpattern(20, (char *) fill20_bits, fill20_width, fill20_height);
  Create_Fillpattern(21, (char *) fill21_bits, fill21_width, fill21_height);
  Create_Fillpattern(22, (char *) fill22_bits, fill22_width, fill22_height);
  Create_Fillpattern(23, (char *) fill23_bits, fill23_width, fill23_height);
  Create_Fillpattern(24, (char *) fill24_bits, fill24_width, fill24_height);
  Create_Fillpattern(25, (char *) fill25_bits, fill25_width, fill25_height);
  Create_Fillpattern(26, (char *) fill26_bits, fill26_width, fill26_height);
  Create_Fillpattern(27, (char *) fill27_bits, fill27_width, fill27_height);
  Create_Fillpattern(28, (char *) fill28_bits, fill28_width, fill28_height);
  Create_Fillpattern(29, (char *) fill29_bits, fill29_width, fill29_height);
  Create_Fillpattern(30, (char *) fill30_bits, fill30_width, fill30_height);
  Create_Fillpattern(31, (char *) fill31_bits, fill31_width, fill31_height);
  Create_Fillpattern(32, (char *) fill32_bits, fill32_width, fill32_height);
  Create_Fillpattern(33, (char *) fill33_bits, fill33_width, fill33_height);
  Create_Fillpattern(34, (char *) fill34_bits, fill34_width, fill34_height);
  Create_Fillpattern(35, (char *) fill35_bits, fill35_width, fill35_height);
  Create_Fillpattern(36, (char *) fill36_bits, fill36_width, fill36_height);
  Create_Fillpattern(37, (char *) fill37_bits, fill37_width, fill37_height);
  Create_Fillpattern(38, (char *) fill38_bits, fill38_width, fill38_height);
  Create_Fillpattern(39, (char *) fill39_bits, fill39_width, fill39_height);
  Create_Fillpattern(40, (char *) fill40_bits, fill40_width, fill40_height);
  Create_Fillpattern(41, (char *) fill41_bits, fill41_width, fill41_height);

  Create_Filldither( 0, (char *)patt0400_bits, patt0400_width, patt0400_height);
  Create_Filldither( 1, (char *)patt0401_bits, patt0401_width, patt0401_height);
  Create_Filldither( 2, (char *)patt0402_bits, patt0402_width, patt0402_height);
  Create_Filldither( 3, (char *)patt0403_bits, patt0403_width, patt0403_height);
  Create_Filldither( 4, (char *)patt0404_bits, patt0404_width, patt0404_height);
  Create_Filldither( 5, (char *)patt1600_bits, patt1600_width, patt1600_height);
  Create_Filldither( 6, (char *)patt1601_bits, patt1601_width, patt1601_height);
  Create_Filldither( 7, (char *)patt1602_bits, patt1602_width, patt1602_height);
  Create_Filldither( 8, (char *)patt1603_bits, patt1603_width, patt1603_height);
  Create_Filldither( 9, (char *)patt1604_bits, patt1604_width, patt1604_height);
  Create_Filldither(10, (char *)patt1605_bits, patt1605_width, patt1605_height);
  Create_Filldither(11, (char *)patt1606_bits, patt1606_width, patt1606_height);
  Create_Filldither(12, (char *)patt1607_bits, patt1607_width, patt1607_height);
  Create_Filldither(13, (char *)patt1608_bits, patt1608_width, patt1608_height);
  Create_Filldither(14, (char *)patt1609_bits, patt1609_width, patt1609_height);
  Create_Filldither(15, (char *)patt1610_bits, patt1610_width, patt1610_height);
  Create_Filldither(16, (char *)patt1611_bits, patt1611_width, patt1611_height);
  Create_Filldither(17, (char *)patt1612_bits, patt1612_width, patt1612_height);
  Create_Filldither(18, (char *)patt1613_bits, patt1613_width, patt1613_height);
  Create_Filldither(19, (char *)patt1614_bits, patt1614_width, patt1614_height);
  Create_Filldither(20, (char *)patt1615_bits, patt1615_width, patt1615_height);
  Create_Filldither(21, (char *)patt1616_bits, patt1616_width, patt1616_height);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

static void Create_Fillpattern(int nr,
                               char *bits,
                               unsigned int width,
                               unsigned int height)
{
  if ((fillpattern[nr] = XCreateBitmapFromData(XviG_display, XviG_dummy_window,
                                               bits, width, height)) == None)
    printf("WARNING : Could not create fillpattern %d.\n", nr);
}

static void Create_Filldither(int nr,
                               char *bits,
                               unsigned int width,
                               unsigned int height)
{
  if ((filldither[nr] = XCreateBitmapFromData(XviG_display, XviG_dummy_window,
                                               bits, width, height)) == None)
    printf("WARNING : Could not create filldither %d.\n", nr);
}

void RBH_SetDither(int nr)
{
  int local_nr;

  if ((local_nr = (nr % 35)) == 0)
    XSetFillStyle(XviG_display, XviG_gc_fill, FillSolid);
  else
  {
    local_nr = Abs(nr);
    XSetFillStyle(XviG_display, XviG_gc_fill, FillOpaqueStippled);
    XSetStipple(XviG_display, XviG_gc_fill, filldither[local_nr]);
  }
}

/*------------------------------------------------------------------------------
--
--

------------------------------------------------------------------------------*/

void XviG_CleanupFillpatterns(void)
{
  int i;

  for (i = 0; i < XviG_NR_OF_FILLPATTERNS; i++)
    if (fillpattern[i] != None)
      XFreePixmap(XviG_display, fillpattern[i]);

/* RBH 1/7/98 */
  for (i = 0; i < RBH_NR_OF_DITHERS; i++)
    if (filldither[i] != None)
      XFreePixmap(XviG_display, filldither[i]);
}
