/* File>>> ftext.c
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
-- Date     : September 27, 1993
--
-- Function : Routines to draw text using fonts.
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

/*------------------------------------------------------------------------------
-- Local type definitions
------------------------------------------------------------------------------*/

typedef struct FontItem aFontItem, *FontItem;

struct FontItem {
  char *name;
  XFontStruct *font_struct;
  int status;               /* (-1 = error) (0 = not yet loaded) (1 = ok) */
};

/*------------------------------------------------------------------------------
-- Static variable declarations
------------------------------------------------------------------------------*/

static FontItem fontList = (FontItem) NULL;
static int nr_of_fonts = 0;
static int current_font = 0;
static char *def_fontlist[5] = {
              "-adobe-courier-medium-r-normal--10-100-75-75-m-60-iso8859-1",
              "-adobe-courier-medium-r-normal--12-120-75-75-m-70-iso8859-1",
              "-adobe-courier-medium-r-normal--14-140-75-75-m-90-iso8859-1",
              "-adobe-courier-medium-r-normal--18-180-75-75-m-110-iso8859-1",
              "-adobe-courier-medium-r-normal--24-240-75-75-m-150-iso8859-1" };

/*------------------------------------------------------------------------------
-- Local function declarations
------------------------------------------------------------------------------*/

static void get_fontnames(void);

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

int XviG_SetFont(int nr)
{
  if (!fontList)
    get_fontnames();

  if (nr_of_fonts == 0)
  {
    fprintf(stderr, "ERROR : No fontnames defined to select.\n");
    current_font = 0;
    return 0;
  }

  if ((nr < 1) || (nr > nr_of_fonts))
  {
    printf("WARNING : Invalid font number %d, using font number %d.\n",
           nr, nr_of_fonts);
    current_font = nr_of_fonts;
  }
  else
    current_font = nr;

  nr = current_font - 1;

  if (fontList[nr].status == 1)
  {
    XSetFont(XviG_display, XviG_gc, fontList[nr].font_struct->fid);
    return 1;
  }

  if (fontList[nr].status == 0)
  {
    if (fontList[nr].font_struct =
                           XLoadQueryFont(XviG_display, fontList[nr].name))
    {
      fontList[nr].status = 1;
      XSetFont(XviG_display, XviG_gc, fontList[nr].font_struct->fid);
      return 1;
    }

    fprintf(stderr, "ERROR : Cannot load font '%s'.\n", fontList[nr].name);
    fontList[nr].status = -1;
    current_font = 0;
    return 0;
  }

  if (fontList[nr].status == -1)
  {
    fprintf(stderr, "ERROR : Cannot select font '%s'.\n", fontList[nr].name);
    current_font = 0;
    return 0;
  }

  current_font = 0;
  return 0;
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

static void get_fontnames(void)
{
  int i;

  if ((i = XviG_GetRMMaxFonts()) == 0)
  {
    fontList = (FontItem) malloc(5 * sizeof(aFontItem));
    fontList[0].status = fontList[1].status = fontList[2].status
                       = fontList[3].status = fontList[4].status
                       = 0;
    fontList[0].name = def_fontlist[0];
    fontList[1].name = def_fontlist[1];
    fontList[2].name = def_fontlist[2];
    fontList[3].name = def_fontlist[3];
    fontList[4].name = def_fontlist[4];
    nr_of_fonts = 5;
  }
  else
  {
    fontList = (FontItem) malloc(i * sizeof(aFontItem));
    nr_of_fonts = i;
    for (i = 0; i < nr_of_fonts; i++)
    {
      fontList[i].status = 0;
      if (!(fontList[i].name = XviG_GetRMFontName(i+1)))
      {
        printf("WARNING : Maxfonts = %d, but only %d fontnames defined.\n",
               nr_of_fonts, i);
        nr_of_fonts = i;
      }
    }
  }
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_FontText(char *contents,
                   int x,
                   int y)
{
  if (current_font == 0)
  {
    fprintf(stderr, "ERROR : No font selected.\n");
    return;
  }

  XDrawString(XviG_display, XviG_window, XviG_gc,
              x, y, contents, (int) strlen(contents));
  XDrawString(XviG_display, XviG_pixmap, XviG_gc,
              x, y, contents, (int) strlen(contents));
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_FontTextSize(char *contents,
                       int *x_offset,
                       int *y_offset,
                       unsigned int *width,
                       unsigned int *height)
{
  int dir_rtn, ascent_rtn, descent_rtn;
  XCharStruct overall_rtn;

  if (current_font == 0)
  {
    fprintf(stderr, "ERROR : No font selected.\n");
    *x_offset = *y_offset = 0;
    *width = *height = 0;
    return;
  }

  XTextExtents(fontList[current_font-1].font_struct,
               contents, (int) strlen(contents),
               &dir_rtn, &ascent_rtn, &descent_rtn, &overall_rtn);

  *x_offset = overall_rtn.lbearing;
  *y_offset = overall_rtn.descent;
  *width = overall_rtn.rbearing - overall_rtn.lbearing;
  *height = overall_rtn.ascent + overall_rtn.descent;
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_CleanupFonts(void)
{
  int i;

  if (fontList)
  {
    for (i = 0; i < nr_of_fonts; i++)
    {
      if (fontList[i].status == 1)
        XFreeFont(XviG_display, fontList[i].font_struct);
    }

    free(fontList);

    fontList = (FontItem) NULL;
    nr_of_fonts = 0;
    current_font = 0;
  }
}
