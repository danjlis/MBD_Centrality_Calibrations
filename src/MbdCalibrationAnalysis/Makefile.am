AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libmbdana.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib



libmbdana_la_LIBADD = \
  -lfun4all \
  -lcalotrigger \
  -lcalotrigger_io \
  -lcentrality_io \
  -leventplaneinfo_io \
  -lcalo_io \
  -lepd_io \
  -lmbd_io \
  -lzdcinfo_io \
  -lffarawobjects \
  -lglobalvertex

pkginclude_HEADERS = \
  MbdAna.h

libmbdana_la_SOURCES = \
  MbdAna.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_mbdana

testexternals_mbdana_SOURCES = testexternals.cc
testexternals_mbdana_LDADD = libmbdana.la


testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
