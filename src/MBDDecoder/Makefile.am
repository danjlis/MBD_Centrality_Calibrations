AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libmbddecoder.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib


libmbddecoder_la_LIBADD = \
  -lcdbobjects \
  -lCLHEP \
  -lffamodules \
  -lgsl \
  -lgslcblas \
  -lglobalvertex_io \
  -lsph_onnx \
  -lphparameter \
  -lphool \
  -lSubsysReco \
  -lcalo_reco

pkginclude_HEADERS = \
  MBDdecoder.h

libmbddecoder_la_SOURCES = \
  MBDdecoder.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_mbddecoder

testexternals_mbddecoder_SOURCES = testexternals.cc
testexternals_mbddecoder_LDADD = libmbddecoder.la


testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
