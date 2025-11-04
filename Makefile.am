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

AM_CXXFLAGS = -Wall -Wextra -Wpedantic

EXTRA_DIST = \
             README.md \
             docs/rseg_manual.tex

noinst_LIBRARIES = librseg.a
librseg_a_SOURCES = \
	src/common/Distro.cpp \
	src/common/EvaluateBoundaries.cpp \
	src/common/LoadReadsByRegion.cpp \
	src/common/ModelParams.cpp \
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

librseg_a_SOURCES += \
	src/common/bam_record_utils.hpp \
	src/common/Distro.hpp \
	src/common/EvaluateBoundaries.hpp \
	src/common/LoadReadsByRegion.hpp \
	src/common/ModelParams.hpp \
	src/common/log_sum_log.hpp \
	src/common/ReadCounts.hpp \
	src/common/rseg_utils.hpp \
	src/common/SelectBinSize.hpp \
	src/common/Smoothing.hpp \
	src/common/SortGenomicRegion.hpp \
	src/common/SplitDistro.hpp \
	src/common/ThreeStateScaleSplitHMM.hpp \
	src/common/ThreeStateScaleSplitResolveMixture.hpp \
	src/common/TwoStateScaleHMM.hpp \
	src/common/TwoStateScaleResolveMixture.hpp \
	src/common/TwoStateScaleSplitHMM.hpp \
	src/common/TwoStateScaleSplitResolveMixture.hpp

LDADD = librseg.a src/smithlab_cpp/libsmithlab_cpp.a

bin_PROGRAMS = rseg rseg-diff deadzones

rseg_SOURCES = src/rseg/rseg.cpp
rseg_diff_SOURCES = src/rseg/rseg-diff.cpp
deadzones_SOURCES = src/rseg/deadzones.cpp
