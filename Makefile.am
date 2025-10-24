# This file is part of rseg
#
# Copyright (C) 2025 Andrew D. Smith
#
# Authors: Andrew D. Smith
#
# This is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This software is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.

ACLOCAL_AMFLAGS = -I m4

SUBDIRS := src/smithlab_cpp
install installdirs: SUBDIRS := $(filter-out src/smithlab_cpp, $(SUBDIRS))
AM_CPPFLAGS = -I $(top_srcdir)/src/common -I $(top_srcdir)/src/smithlab_cpp -I $(top_srcdir)/src/bamxx

AM_CXXFLAGS = -Wall -Wextra -Wpedantic -Wno-unknown-attributes

EXTRA_DIST = \
             README.md \
             docs/rseg_manual.tex

noinst_LIBRARIES = librseg.a
librseg_a_SOURCES = \
	src/common/BAMFile.cpp \
	src/common/Distro.cpp \
	src/common/EvaluateBoundaries.cpp \
	src/common/LoadReadsByRegion.cpp \
	src/common/ModelParams.cpp \
	src/common/numerical_utils.cpp \
	src/common/ReadCounts.cpp \
	src/common/rseg_utils.cpp \
	src/common/SelectBinSize.cpp \
	src/common/Smoothing.cpp \
	src/common/SortGenomicRegion.cpp \
	src/common/SplitDistro.cpp \
	src/common/ThreeStateScaleSplitHMM.cpp \
	src/common/ThreeStateScaleSplitResolveMixture.cpp \
	src/common/TwoStateScaleHMM.cpp \
	src/common/TwoStateScaleResolveMixture.cpp \
	src/common/TwoStateScaleSplitHMM.cpp \
	src/common/TwoStateScaleSplitResolveMixture.cpp

# librseg_a_SOURCES += \
# 	src/bamxx/bamxx.hpp \
#         src/common/rseg_gaussinv.hpp \
# 	src/common/bam_record_utils.hpp \
# 	src/common/counts_header.hpp \
# 	src/common/xcounts_utils.hpp \
#         src/common/BetaBin.hpp \
#         src/common/EmissionDistribution.hpp \
#         src/common/Epiread.hpp \
#         src/common/EpireadStats.hpp \
#         src/common/LevelsCounter.hpp \
#         src/common/MSite.hpp \
#         src/common/Smoothing.hpp \
#         src/common/ThreeStateHMM.hpp \
#         src/common/TwoStateHMM.hpp \
#         src/common/TwoStateHMM_PMD.hpp \
#         src/common/bsutils.hpp \
#         src/common/numerical_utils.hpp \
#         src/common/rseg_utils.hpp \
#         src/common/dnmt_error.hpp

# ADS: additional radmeth sources to help isolate the parts using GSL for
# optimization
librseg_a_SOURCES += \
        src/radmeth/radmeth_optimize_series.hpp \
        src/radmeth/radmeth_optimize_gamma.hpp \
        src/radmeth/radmeth_optimize_params.hpp \
        src/radmeth/radmeth_model.hpp \
        src/radmeth/radmeth_lgamma.hpp \
        src/radmeth/radmeth_utils.hpp \
        src/radmeth/radmeth_design.hpp

LDADD = librseg.a src/smithlab_cpp/libsmithlab_cpp.a

bin_PROGRAMS = rseg

rseg_SOURCES = src/rseg/rseg.cpp

## ADS: these are the files output by the test scripts. They can be
## identified by doing:
##
## grep "^outfile" test_scripts/*.test
##
# CLEANFILES = \
#     tests/tests.testout
