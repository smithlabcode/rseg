#    Makefile from common code directory of RSEG source tree
#
#    Copyright (C) 2011 University of Southern California and
#                       Andrew D. Smith
#
#    Authors: Qiang Song and Andrew D. Smith
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

ifndef RSEG_ROOT
$(error Must define RSEG_ROOT variable)
endif

SMITHLAB_CPP = ../smithlib
COMMON_LIB = ../common

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
PROGS =  rseg rseg-diff deadzones join

INCLUDEDIRS = $(SMITHLAB_CPP) $(COMMON_LIB)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS += -lgsl -lgslcblas 

CXX = g++
CXXFLAGS = -Wall -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g -lefence -L$(HOME)/lib

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

rseg: $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o RNG.o) \
	$(addprefix $(COMMON_LIB)/, 	numerical_utils.o rseg_utils.o \
	Distro.o SplitDistro.o Smoothing.o \
	SortGenomicRegion.o LoadReadsByRegion.o ReadCounts.o  SelectBinSize.o \
	TwoStateScaleHMM.o TwoStateResolveMixture.o \
	EvaluateBoundaries.o)		

rseg-diff: $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o RNG.o) \
	$(addprefix $(COMMON_LIB)/, numerical_utils.o rseg_utils.o \
	Distro.o SplitDistro.o Smoothing.o \
	SortGenomicRegion.o LoadReadsByRegion.o ReadCounts.o SelectBinSize.o \
	TwoStateScaleSplitHMM.o TwoStateScaleSplitResolveMixture.o \
	ThreeStateScaleSplitHMM.o ThreeStateScaleSplitResolveMixture.o \
	EvaluateBoundaries.o)		

deadzones: $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

install: $(PROGS)
	@mkdir -p $(RSEG_ROOT)/bin
	@install -m 755 $(PROGS) $(RSEG_ROOT)/bin

clean:
	@-rm -f $(PROGS) *.o *~

.PHONY: clean

