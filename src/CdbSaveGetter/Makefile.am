AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libcdbsavegetter.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib

libdbsavegetter_la_LIBADD = \
  -lfun4all \
  -phool

pkginclude_HEADERS = \
  CdbSaveGetter.h

libcdbsavegetter_la_SOURCES = \
  CdbSaveGetter.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_cdbsavegetter

testexternals_cdbsavegetter_SOURCES = testexternals.cc
testexternals_cdbsavegetter_LDADD = libcdbsavegetter.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
