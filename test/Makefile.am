LIBNCPOL2SDPA = $(top_builddir)/src/libncpol2sdpa-1.0.la
AM_CPPFLAGS = -I$(top_builddir)/src
bin_PROGRAMS = exampleNcPol benchmarkCase
exampleNcPol_SOURCES = exampleNcPol.cpp
exampleNcPol_LDADD = $(LIBNCPOL2SDPA)
benchmarkCase_SOURCES = benchmarkCase.cpp
benchmarkCase_LDADD = $(LIBNCPOL2SDPA)
