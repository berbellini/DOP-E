/* File>>> rmdb.c
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
#include <sys/file.h>
#include <X11/Xlib.h>
#include <X11/Xresource.h>

#ifdef SOLARIS
#include <X11/Xw32defs.h>
#endif
#include <unistd.h>

/*------------------------------------------------------------------------------
-- Local include files
------------------------------------------------------------------------------*/

#include "xviglocal.h"

/*------------------------------------------------------------------------------
-- Static variable declarations
------------------------------------------------------------------------------*/

static XrmDatabase rm_database = (XrmDatabase) NULL;

/*------------------------------------------------------------------------------
-- Local function declarations
------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_GetRMDatabase(void)
{
  char db_pathname[256];
  char *homedir;

  if (!(homedir = getenv("HOME")))
    return;

  strcpy(db_pathname, homedir);
  strcat(db_pathname, "/.Xdefaults");

  if (access(db_pathname, R_OK) == -1)
    return;

  rm_database = XrmGetFileDatabase(db_pathname);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_DeleteRMDatabase(void)
{
  if (rm_database)
  {
    XrmDestroyDatabase(rm_database);
    rm_database = (XrmDatabase) NULL;
  }
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

char *XviG_GetRMGeometry(char *name)
{
  char str_name[512];
  char *str_type_rtn;
  XrmValue geometry;

  if (!rm_database)
    return (char *) NULL;

  strcpy(str_name, XviG_class);
  strcat(str_name, ".calxvig.");
  strcat(str_name, name);
  strcat(str_name, ".geometry");

  if (XrmGetResource(rm_database, str_name, XviG_class,
                     &str_type_rtn, &geometry))
    return (char *) geometry.addr;

  return (char *) NULL;
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

int XviG_GetRMMaxFonts(void)
{
  char str_name[256];
  char *str_type_rtn, *endptr;
  XrmValue maxfonts;
  int nr;

  if (!rm_database)
    return 0;

  strcpy(str_name, XviG_class);
  strcat(str_name, ".calxvig.maxfonts");

  if (XrmGetResource(rm_database, str_name, XviG_class,
                     &str_type_rtn, &maxfonts))
  {
    nr = (int) strtol((char *) maxfonts.addr, &endptr, 10);

    if ((*endptr != '\0') || (nr < 0))
      return 0;
    else
      return nr;
  }

  return 0;
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

char *XviG_GetRMFontName(int nr)
{
  char str_name[256], str_nr[8];
  char *str_type_rtn;
  XrmValue fontname;

  sprintf(str_nr, "%d", nr);

  strcpy(str_name, XviG_class);
  strcat(str_name, ".calxvig.font");
  strcat(str_name, str_nr);

  if (XrmGetResource(rm_database, str_name, XviG_class,
                     &str_type_rtn, &fontname))
    return (char *) fontname.addr;

  return (char *) NULL;
}
