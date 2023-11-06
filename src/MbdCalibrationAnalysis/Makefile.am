AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include  \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = \
   libmbdcalibrationanalysis.la \
   libmbdana.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib


libmbdcalibrationanalysis_la_LIBADD = \
  -lfun4all \
  -lcalo_io

libmbdana_la_LIBADD = \
  -lfun4all \
  -lcalo_io \
  -lmbd_io \
  -lglobalvertex

pkginclude_HEADERS = \
  MbdCalibrationAnalysis.h \
  MbdAna.h

libmbdcalibrationanalysis_la_SOURCES = \
  MbdCalibrationAnalysis.cc

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
  testexternals_mbdcalibrationanalysis \
  testexternals_mbdana

testexternals_mbdcalibrationanalysis_SOURCES = testexternals.cc
testexternals_mbdcalibrationanalysis_LDADD = libmbdcalibrationanalysis.la

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
