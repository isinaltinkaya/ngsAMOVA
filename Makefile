####################################################################################################
# ngsAMOVA Makefile


CXX ?= g++

#  -g     add debugging info to the executable 
#  -Wall  turn on compiler warnings
#CXXFLAGS  := -g -Wall
LIBS = -lz -lm -lbz2 -llzma -lcurl -lpthread



####################################################################################################
## [developer mode]
# - to compile in developer mode, use `make dev`
# - to compile in release mode, use `make`

DEV_MODE=0
ifeq ($(filter dev,$(MAKECMDGOALS)),dev)
	DEV_MODE=1
	CXXFLAGS  := -g -Wall
endif

DEV_VAL=$(shell grep "#define DEV [0,1]" dev.h | cut -d " " -f 3)

# - if DEV_MODE is already set to the desired value, do nothing; otherwise, change it in dev.h
ifeq ($(DEV_VAL),$(DEV_MODE))
$(info [INFO] No change in DEV_MODE)
else
$(info [INFO] Changing DEV_MODE)
$(shell sed -i "s/#define DEV [0-1]/#define DEV $(DEV_MODE)/g" dev.h)
endif

dev: all



####################################################################################################
## [cryptolib check]
# - check if cryptolib is available to link
CRYPTO_TRY=$(shell echo 'int main(){}'|$(CXX) -x c++ - -lcrypto 2>/dev/null -o /dev/null; echo $$?)
ifeq "$(CRYPTO_TRY)" "0"
$(info Crypto library is available to link; adding -lcrypto to LIBS)
LIBS += -lcrypto
else
$(info Crypto library is not available to link; will not use -lcrypto)
endif




####################################################################################################
## [htslib source check]
# - to define htslib source, use `make HTSSRC=/path/to/htslib`
# - to use the systemwide htslib installation, use `make HTSSRC=systemwide`
# - to use the submodule, just use `make`

#if htslib source is defined
ifdef HTSSRC

#if hts source is set to systemwide
ifeq ($(HTSSRC),systemwide)
$(info HTSSRC set to systemwide; assuming systemwide installation)
LIBHTS := -lhts

else

#if hts source path is given
# Adjust $(HTSSRC) to point to your top-level htslib directory
$(info HTSSRC defined: $(HTSSRC))
CXXFLAGS += -I"$(realpath $(HTSSRC))"
LIBHTS := $(realpath $(HTSSRC))/libhts.a

endif

#if htssrc not defined
else

$(info HTSSRC not defined; using htslib submodule)
$(info Use `make HTSSRC=/path/to/htslib` to build using a local htslib installation)
$(info Use `make HTSSRC=systemwide` to build using the systemwide htslib installation)

HTSSRC := $(realpath $(CURDIR)/htslib)
CXXFLAGS += -I"$(HTSSRC)"
LIBHTS := $(HTSSRC)/libhts.a


all: .activate_module

endif

.PHONY: .activate_module test clean 

.activate_module:
	git submodule update --init --recursive
	$(MAKE) -C $(HTSSRC)

	
PROGRAM = ngsAMOVA
all: $(PROGRAM)

CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)

-include $(OBJ:.o=.d)

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d


$(PROGRAM): $(OBJ)
	$(CXX) -o $(PROGRAM) *.o $(LIBHTS) $(LIBS) 


####################################################################################################
## [clean]
# - to clean up the directory, use `make clean`
clean:
	$(RM) *.o *.d $(PROGRAM)


####################################################################################################
## [unit tests]
# - to run unit tests, use `make test`
# - tests are defined in test/Makefile
test:
	$(MAKE) -C test
####################################################################################################

