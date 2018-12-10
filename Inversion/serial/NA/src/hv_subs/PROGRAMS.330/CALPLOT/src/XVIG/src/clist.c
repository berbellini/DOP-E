/* File>>> clist.c
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

/*------------------------------------------------------------------------------
-- Static variable declarations
------------------------------------------------------------------------------*/

static CurItem startCursor = (CurItem) NULL,
               stopCursor = (CurItem) NULL;

/*------------------------------------------------------------------------------
-- Local function declarations
------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

CurItem New_CurItem(Pixmap pixmap,
                    unsigned int width,
                    unsigned int height,
                    int hot_x,
                    int hot_y)
{
  CurItem curcur;

  curcur = (CurItem) malloc(sizeof(aCurItem));

  curcur->pixmap = pixmap;
  curcur->width = width;
  curcur->height = height;
  curcur->hot_x = hot_x;
  curcur->hot_y = hot_y;
  curcur->next = (CurItem) NULL;

  if (!startCursor)
  {
    curcur->prev = (CurItem) NULL;
    startCursor = stopCursor = curcur;
  }
  else
  {
    curcur->prev = stopCursor;
    stopCursor->next = curcur;
    stopCursor = curcur;
  }

  return curcur;
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void Delete_CurItem(CurItem curcur)
{
  if (curcur == startCursor)
  {
    if (curcur == stopCursor)
      startCursor = (CurItem) NULL;
    else
    {
      curcur->next->prev = (CurItem) NULL;
      startCursor = curcur->next;
    }
  }
  else
  if (curcur == stopCursor)
  {
    curcur->prev->next = (CurItem) NULL;
    stopCursor = curcur->prev;
  }
  else
  {
    curcur->prev->next = curcur->next;
    curcur->next->prev = curcur->prev;
  }

  /*
  -- Free it
  */

  free(curcur);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_CleanupCursors(void)
{
  CurItem loopcur, nextcur;

  for (loopcur = startCursor; loopcur; loopcur = nextcur)
  {
    nextcur = loopcur->next;
    XviG_DeleteCursor((long) loopcur);
  }
}
