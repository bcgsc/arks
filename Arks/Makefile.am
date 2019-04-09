bin_PROGRAMS = arks

arks_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

arks_CPPFLAGS = -I$(top_srcdir)/Arks \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)

arks_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a -lz

arks_LDFLAGS = $(OPENMP_CXXFLAGS)

arks_SOURCES = \
	Arks.h \
	Arks.cpp \
	DistanceEst.h \
	kseq.h	
