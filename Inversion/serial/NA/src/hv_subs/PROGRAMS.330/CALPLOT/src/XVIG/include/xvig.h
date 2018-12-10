/* File>>> xvig.h
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
-- Date     : October 1, 1993
--
-- Function : Header file for XviG version 1.1
--
-- Comment  :
--
-- Review   :
--
*/


#ifndef __XVIG_H
#define __XVIG_H


#define ON	1
#define OFF	0
/*
-- The default colors, if not set by the user program
*/

#define XviG_NR_OF_COLORS        16

#define XviG_COLOR_BLACK          0
#define XviG_COLOR_WHITE          1
#define XviG_COLOR_RED            2
#define XviG_COLOR_GREEN          3
#define XviG_COLOR_BLUE           4
#define XviG_COLOR_CYAN           5
#define XviG_COLOR_MAGENTA        6
#define XviG_COLOR_YELLOW         7
#define XviG_COLOR_ORANGE         8
#define XviG_COLOR_GREEN_YELLOW   9
#define XviG_COLOR_GREEN_CYAN    10
#define XviG_COLOR_BLUE_CYAN     11
#define XviG_COLOR_BLUE_MAGENTA  12
#define XviG_COLOR_RED_MAGENTA   13
#define XviG_COLOR_DARK_GRAY     14
#define XviG_COLOR_LIGHT_GRAY    15


/*
-- The linestyles and fillpatterns
*/

#define XviG_NR_OF_LINESTYLES     8
#define XviG_NR_OF_FILLPATTERNS  42


/*
-- The cursor types
*/

#define XviG_CURSOR_ARROW	 0L
#define XviG_CURSOR_XORARROW	-1L
#define XviG_CURSOR_XHAIR 	-2L
#define XviG_CURSOR_PLUS  	-3L
#define XviG_CURSOR_BOX  	-4L
#define XviG_CURSOR_RUBBER  	-5L
#define XviG_CURSOR_OFF  	-6L
#define XviG_CURSOR_HYPERBOLA  	-8L


/*
-- The cursor input enabling
*/

#define XviG_KEY         0
#define XviG_BUTTON      1
#define XviG_KEY_BUTTON  2


/*
-- The mouse button return values
*/

#define XviG_BUTTON1  -1
#define XviG_BUTTON2  -2
#define XviG_BUTTON3  -3
#define XviG_BUTTON4  -4
#define XviG_BUTTON5  -5


/*
-- The functions
*/

extern int XviG_Init(char *classname,
                     int color_array[][3],
                     int nr_of_colors);

extern void XviG_Exit(void);

extern int XviG_OpenWindow(char *name,
                           int x,
                           int y,
                           unsigned int *width,
                           unsigned int *height);

extern int XviG_CloseWindow(char *name);

extern int XviG_SelectWindow(char *name);

extern void XviG_WindowSize(unsigned int *width,
                            unsigned int *height);

extern void XviG_WindowPosition(int *x,
                                int *y);

extern void XviG_ClearWindow(void);

extern void XviG_SetColor(int nr);

extern void XviG_SetLineStyle(int nr,
                              unsigned int width);

extern void XviG_SetFillStyle(int nr);

extern int XviG_SetFont(int nr);

extern void XviG_DrawPoint(int x,
                           int y);

extern void XviG_DrawLine(int x1,
                          int y1,
                          int x2,
                          int y2);

extern void XviG_DrawPolyLine(int *coords,
                              int npoints);

extern void XviG_DrawRectangle(int x1,
                               int y1,
                               int x2,
                               int y2);

extern void XviG_DrawPolygon(int *coords,
                             int npoints);

extern void XviG_DrawArc(int x,
                         int y,
                         unsigned int radius1,
                         unsigned int radius2,
                         int angle1,
                         int angle2);

extern void XviG_FillRectangle(int x1,
                               int y1,
                               int x2,
                               int y2);

extern void XviG_FillPolygon(int *coords,
                             int npoints);

extern void XviG_FillArc(int x,
                         int y,
                         unsigned int radius1,
                         unsigned int radius2,
                         int angle1,
                         int angle2);

extern void XviG_PolyText(char *contents,
                          int x,
                          int y,
                          unsigned int width,
                          unsigned int height,
                          int rotation);

extern void XviG_FontText(char *contents,
                          int x,
                          int y);

extern void XviG_FontTextSize(char *contents,
                              int *x_offset,
                              int *y_offset,
                              unsigned int *width,
                              unsigned int *height);

extern void XviG_Flush(void);

extern void XviG_SetSenseKbd(int sense_char);

extern int XviG_SenseKbd(void);

extern void XviG_OpenCursor(unsigned int width,
                            unsigned int height,
                            int hot_x,
                            int hot_y);

extern long XviG_CloseCursor(void);

extern void XviG_DeleteCursor(long cursor);

extern void XviG_SetCursor(long cursor);

extern int XviG_GetCursor(int type,
                          int *x_pos,
                          int *y_pos);

/* rbh extension */
extern void XviG_SendMessage(int type, int i1, int i2, int i3, int i4);
extern int XviG_GetChar(void);
#define RBH_NR_OF_DITHERS  22
extern void RBH_SetDither(int nr);
void XviG_SetGC(int gc);

#endif  /* __XVIG_H */
