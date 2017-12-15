#===============================================================================
# Makefile.am : for automake
# Copyright (C) 2006-2017 Jens Kleinjung
#===============================================================================

bin_PROGRAMS =  minset

AM_CPPFLAGS = $(INTI_CFLAGS)
AM_CFLAGS = -Wall

minset_SOURCES = \
alphabet.c alphabet.h cath.c ga.c ga.h gapar.h getseqs.c getseqs.h \
homstrad.c huffman.c huffman.h lz.c lz.h minset.c minset.h \
minsetpar.h parse_args.c parse_args.h suffix_tree.c suffix_tree.h

minset_LDADD = $(INTI_LIBS)

EXTRA_DIST = doxygen.cfg

CLEANFILES = $(TARGET) *.o 
DISTCLEANFILES = libtool config.cache config.log

