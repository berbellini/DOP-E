# $Id: makefile.in,v 1.62 2008/07/20 18:30:36 tom Exp $
# template makefile for DIALOG
#
SHELL		= /bin/sh

srcdir		= .


prefix		= /usr/local
exec_prefix	= ${prefix}

top_builddir	= .

x		= 
o		= .o
a		= .a

DESTDIR		=

bindir		= ${exec_prefix}/bin
includedir	= ${prefix}/include
libdir		= ${exec_prefix}/lib
mandir		= ${prefix}/man

BINDIR		= $(DESTDIR)$(bindir)
INCLUDEDIR	= $(DESTDIR)$(includedir)
LIBDIR		= $(DESTDIR)$(libdir)
MANDIR		= $(DESTDIR)$(mandir)

MAN1DIR		= $(MANDIR)/man1
MAN3DIR		= $(MANDIR)/man3

# see po/makefile
localedir       = $(prefix)/share/locale

CFLAGS		= -g -O2
CPPFLAGS	=  -D_GNU_SOURCE -DHAVE_CONFIG_H -I. -I$(srcdir) -DLOCALEDIR=\"$(localedir)\"
EXTRA_CFLAGS	= 
CC		= gcc
CPP		= gcc -E
AR		= ar
LDFLAGS		= 
LIBS		= -lncurses  -lm -L/lib # 
RANLIB		= ranlib

RM		= rm -f
LINT		= lint

LIBTOOL		=  
LIBTOOL_OPTS	= 
LIBTOOL_CLEAN	= 
LIBTOOL_COMPILE	= 
LIBTOOL_CREATE	= ${AR} -cr
LIBTOOL_LINK	= ${CC}
LIBTOOL_INSTALL	= 
LIBTOOL_UNINSTALL = 

INSTALL		= /usr/bin/install -c
INSTALL_PROGRAM	= $(LIBTOOL_INSTALL) ${INSTALL}
INSTALL_SCRIPT	= ${INSTALL}
INSTALL_DATA	= ${INSTALL} -m 644

LINK		= $(LIBTOOL_LINK)

PACKAGE		= `echo dialog |sed s,x,x,`
#
# Standard .c to .o compile line.
#
.SUFFIXES: .c .i $o
.c.i :
# compiling
	$(CPP) $(CPPFLAGS) -c $< >$@
.c$o :
# compiling
	$(LIBTOOL_COMPILE) $(CC) $(CFLAGS) $(EXTRA_CFLAGS) $(CPPFLAGS) -c $<

EXTRAOBJS =  trace$o rc$o calendar$o fselect$o timebox$o formbox$o guage$o pause$o progressbox$o tailbox$o mixedform$o mixedgauge$o
OBJECTS = $(EXTRAOBJS) \
	arrows$o \
	buttons$o \
	checklist$o \
	columns$o \
	dlg_keys$o \
	editbox$o \
	inputbox$o \
	inputstr$o \
	menubox$o \
	mouse$o \
	mousewget$o \
	msgbox$o \
	textbox$o \
	ui_getc$o \
	util$o \
	version$o \
	yesno$o

LIB_OBJECT = ${OBJECTS}

SRCS = $(OBJECTS:$o=.c)
HDRS = \
	dialog.h \
	dlg_colors.h \
	dlg_config.h \
	dlg_keys.h

LIB	= libdialog$a
PROG	= dialog$x
ALL	= $(LIB) $(PROG)

BIN_DIRS	= $(BINDIR) $(MAN1DIR)
LIB_DIRS	= $(LIBDIR) $(INCLUDEDIR)

PO_DIR	= $(srcdir)/po

all	:: $(ALL)

#install	:: install-lib
#uninstall	:: uninstall-lib

# :
#	@echo "Building GNU gettext library..."
#	@cd intl && $(MAKE) 

dialog$o \
$(OBJECTS) : $(srcdir)/dialog.h $(srcdir)/dlg_keys.h dlg_config.h VERSION

$(LIB) : $(LIB_OBJECT)
	$(LIBTOOL_CREATE) $(LIB) $(LIB_OBJECT)
	$(RANLIB) $@

dialog$x : $(LIB) dialog$o # 
	$(LINK) -o $@ dialog$o -L. -ldialog $(LDFLAGS) $(LIBS)

clean	\
distclean \
install \
install-full \
install-strip \
uninstall ::
	@echo making $@
#	@-test -f intl/makefile && cd intl && $(MAKE)  $@
#	@-test -f po/makefile && cd po && $(MAKE)  $@

clean	::
	@- $(RM) -r autom4te.cache
	- $(RM) *.*cov *.da *.bb *.bbg
	- $(LIBTOOL_CLEAN) $(RM) *.i *$o $(ALL)
	- $(RM) headers.sed
	- $(RM) core *~ tags TAGS

distclean :: clean
	$(RM) intl/libintl.h intl/po2tbl.sed
	$(RM) makefile dlg_config.h config.cache config.log config.status
	$(RM) samples/install/makefile
	$(RM) headers-sh
	$(RM) dialog-config

install :: install-bin install-man
	@echo "made $@"

install-full :: install-lib install-bin install-man
	@echo "made $@"

install-strip :: $(PROG) $(BINDIR)
	$(INSTALL_PROGRAM) -s $(PROG) $(BINDIR)

uninstall :: uninstall-bin uninstall-man
	@echo "made $@"

uninstall-full :: uninstall-bin uninstall-man uninstall-lib
	@echo "made $@"

install-bin : $(PROG) $(BINDIR) 
	$(INSTALL_PROGRAM) $(PROG) $(BINDIR)/$(PACKAGE)$x

uninstall-bin :
	$(RM) $(BINDIR)/$(PACKAGE)$x

install-man :: $(MAN1DIR)
	$(INSTALL_DATA) $(srcdir)/dialog.1 $(MAN1DIR)/$(PACKAGE).1

uninstall-man ::
	$(RM) $(MAN1DIR)/$(PACKAGE).1

# most users do not want/need the library, so the install rules are distinct.
install-lib :: $(BINDIR) dialog-config
	@ echo "** installing config script in $(BINDIR)"
	@ $(INSTALL_SCRIPT) dialog-config $(BINDIR)

install-lib :: $(LIB_DIRS) $(LIB) headers.sed
	@ echo "** installing library in $(LIBDIR)"
	@ $(LIBTOOL_INSTALL) $(INSTALL_DATA) $(LIB) $(LIBDIR)
	@ echo "** installing headers in $(INCLUDEDIR)"
	@ $(SHELL) $(srcdir)/headers-sh $(INSTALL_DATA) $(INCLUDEDIR) $(srcdir) dialog.h
	@ $(SHELL) $(srcdir)/headers-sh $(INSTALL_DATA) $(INCLUDEDIR) $(srcdir) dlg_colors.h
	@ $(SHELL) $(srcdir)/headers-sh $(INSTALL_DATA) $(INCLUDEDIR) $(srcdir) dlg_keys.h
	@ $(SHELL) $(srcdir)/headers-sh $(INSTALL_DATA) $(INCLUDEDIR) .         dlg_config.h

uninstall-lib :: $(BINDIR)
	$(RM) $(BINDIR)/dialog-config

uninstall-lib :: $(LIB_DIRS)
	- $(LIBTOOL_UNINSTALL) $(RM) $(LIBDIR)/$(LIB)
	$(RM) $(INCLUDEDIR)/dialog.h
	$(RM) $(INCLUDEDIR)/dlg_colors.h
	$(RM) $(INCLUDEDIR)/dlg_keys.h
	$(RM) $(INCLUDEDIR)/dlg_config.h

install-lib :: $(MAN3DIR)
	$(INSTALL_DATA) $(srcdir)/dialog.3 $(MAN3DIR)/$(PACKAGE).3

uninstall-lib ::
	$(RM) $(MAN3DIR)/$(PACKAGE).3

headers.sed : $(srcdir)/headers-sh
	$(SHELL) $(srcdir)/headers-sh $(INCLUDEDIR) $(srcdir)

$(MAN1DIR) \
$(MAN3DIR) \
$(BINDIR) \
$(INCLUDEDIR) \
$(LIBDIR) : ; $(srcdir)/mkdirs.sh $@

#tags :
#	ctags $(SRCS) $(HDRS)

#TAGS :
#	etags $(SRCS) $(HDRS)

update-po:
	rsync -Lrtvz  translationproject.org::tp/latest/dialog/  $(PO_DIR)
	test -f $(PO_DIR)/makefile && cd $(PO_DIR) && $(MAKE) $@

lint:
	$(LINT) $(CPPFLAGS) *.c
