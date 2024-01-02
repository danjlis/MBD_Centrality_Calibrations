AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libinttmbd.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib


libinttmbd_la_LIBADD = \
  -lfun4all \
  -lfun4allraw \
  -lEvent \
  -lmbd_io \
  -lffarawobjects \
  -lphool \
  -lSubsysReco

pkginclude_HEADERS = \
  InttMbd.h

libinttmbd_la_SOURCES = \
  InttMbd.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_inttmbd

testexternals_inttmbd_SOURCES = testexternals.cc
testexternals_inttmbd_LDADD = libinttmbd.la


testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
