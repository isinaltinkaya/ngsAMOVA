CXX ?= g++

#  -g     add debugging info to the executable 
#  -Wall  turn on most compiler warnings
CXXFLAGS  = -g -Wall
LIBS = -lz -lhts

# HTSSRC := $(CURDIR)/htslib
# HTSSRC := /maps/projects/lundbeck/scratch/pfs488/AMOVA/vcfToGlf/htslib


#if htslib source is defined
ifdef HTSSRC

#if hts source is set to systemwide
ifeq ($(HTSSRC),systemwide)
$(info HTSSRC set to systemwide; assuming systemwide installation)
LIBS += -lhts

else

#if hts source path is given
# Adjust $(HTSSRC) to point to your top-level htslib directory
$(info HTSSRC defined: $(HTSSRC))
CXXFLAGS += -I"$(realpath $(HTSSRC))"
LIBHTS := $(HTSSRC)/libhts.a
LIBS := $(LIBHTS) $(LIBS)

endif

#if htssrc not defined
else

$(info HTSSRC not defined; using htslib submodule)
$(info Use `make HTSSRC=/path/to/htslib` to build using a local htslib installation)
$(info Use `make HTSSRC=systemwide` to build using the systemwide htslib installation)


HTSSRC := $(CURDIR)/htslib
CXXFLAGS += -I$(HTSSRC)
LIBHTS := $(HTSSRC)/libhts.a
LIBS := $(LIBHTS) $(LIBS)

all: .activate_module

endif

.PHONY: .activate_module test

.activate_module:
	git submodule update --init --recursive
	$(MAKE) -C $(HTSSRC)




# CPPFLAGS += -I$(HTSSRC)



PROGRAM = cmpgt

CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)


-include $(OBJ:.o=.d)

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d



all: $(PROGRAM) $(OBJ)

$(PROGRAM): $(OBJ)
	$(CXX) $(CPPFLAGS) $(LIBS) -o $(PROGRAM) *.o 

clean:
	$(RM) *.o *.d $(PROGRAM)

test: 
	./vcfgl -in test/t1.vcf -out test/t1_testvcfgl -err 0.01 -seed 42 -depth 1 -isSim 1;
	bash -c "diff test/t1_vcf_vcfgl.vcf test/t1_testvcfgl.vcf" ;
