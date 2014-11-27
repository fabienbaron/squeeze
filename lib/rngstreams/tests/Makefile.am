## Process this file with automake to produce Makefile.in

SUBDIRS = .

INCLUDES = -I$(top_srcdir)/src

LDADD = $(top_builddir)/src/librngstreams.la

check_PROGRAMS = \
 test2RngStream \
 testRngStream

TESTS = \
 test2RngStream \
 testRngStream

#EXTRA_DIST = \
# test2RngStream.res

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = Makefile.in

