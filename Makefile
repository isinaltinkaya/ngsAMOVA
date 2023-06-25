####################################################################################################
# ngsAMOVA Makefile
####################################################################################################


####################################################################################################
#
# `make` 		- compile in release mode
#
# `make dev` 	- compile in developer mode
# 				-O0 (no optimization)
# 				-g (debugging info)
# 				-Wall (turn on compiler warnings)
# 				-save-temps (save intermediate files; useful for seeing macro expansions)
# 				-v (verbose)
#
# `make clean` 	- clean up the director
#
####################################################################################################


####################################################################################################
# -- BEGIN BLOCK --
# only run the block if not clean mode
ifneq (clean,$(filter clean,$(MAKECMDGOALS)))
#


# Compiler specs
CXX ?= g++

# Libraries to link
LIBS := -lz -lm -lbz2 -llzma -lcurl -lpthread
$(info [INFO]    Make received LIBS="$(LIBS)")


# avoid using assert statements for safety
# define NDEBUG to disable assert statements
CXXFLAGS := -DNDEBUG
$(info [INFO]    Make received CXXFLAGS="$(CXXFLAGS)")

OPTIM_FLAGS := -O3


DEVH_VAL=$(shell grep "define DEV [0-2]" dev.h | cut -d " " -f 3)


ifeq (dev,$(filter dev,$(MAKECMDGOALS)))
##

$(info )
$(info [INFO]    Compiling in developer mode)
$(info _________________________________________________________________________________)

DEV_MODE=1
OPTIM_FLAGS := -O0
CXXFLAGS += -g -Wall $(OPTIM_FLAGS) -save-temps -v
$(info )
$(info [INFO]    Make received CXXFLAGS="$(CXXFLAGS)")
$(info )

else
##


DEV_MODE=0
$(info _________________________________________________________________________________)
$(info )
$(info [INFO]    Compiling in release mode)
$(info _________________________________________________________________________________)
$(info )

OPTIM_FLAGS := -O3
CXXFLAGS += $(OPTIM_FLAGS)

endif
##

ifneq ($(DEVH_VAL),$(DEV_MODE))
##

$(info [INFO]    Changing DEV in dev.h from $(DEVH_VAL) to $(DEV_MODE))
$(shell sed -i "1s/define DEV $(DEVH_VAL)/define DEV $(DEV_MODE)/g" dev.h)

endif
##

#
endif
# -- END BLOCK --
####################################################################################################

dev: all


####################################################################################################
## [cryptolib check]
# - check if cryptolib is available to link
CRYPTO_TRY=$(shell echo 'int main(){}'|$(CXX) -x c++ - -lcrypto 2>/dev/null -o /dev/null; echo $$?)

ifeq "$(CRYPTO_TRY)" "0"
$(info [INFO]    Crypto library is available to link; adding -lcrypto to LIBS)
LIBS += -lcrypto

else

$(info [INFO]    Crypto library is not available to link; will not use -lcrypto)
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

$(info [INFO]    HTSSRC set to systemwide; assuming systemwide installation)
LIBHTS := -lhts

else

#if hts source path is given
# Adjust $(HTSSRC) to point to your top-level htslib directory
$(info [INFO]    HTSSRC defined: $(HTSSRC))
CPPFLAGS += -I"$(realpath $(HTSSRC))"
LIBHTS := $(realpath $(HTSSRC))/libhts.a

endif

#if htssrc not defined
else

$(info [INFO]    HTSSRC not defined; using htslib submodule)
$(info [INFO]       `make HTSSRC=/path/to/htslib` to build using a local htslib installation)
$(info [INFO]       `make HTSSRC=systemwide` to build using the systemwide htslib installation)

HTSSRC := $(realpath $(CURDIR)/htslib)
CPPFLAGS += -I"$(HTSSRC)"
LIBHTS := $(HTSSRC)/libhts.a


all: .activate_module

endif


.PHONY: .activate_module
.activate_module:
	$(info [INFO]    Activating htslib submodule)
	git submodule update --init --recursive
	$(MAKE) -C $(HTSSRC)

	

####################################################################################################



PROGRAM = ngsAMOVA
all: $(PROGRAM)

CXXSRC := $(wildcard *.cpp)

# Preprocessed C++ files
PREP := $(CXXSRC:.cpp=.ii)

# Assembly source files
ASM := $(CXXSRC:.cpp=.s)

# Object files
OBJ := $(CXXSRC:.cpp=.o)

# Dependency files
DEP := $(OBJ:.o=.d)


LIBS += $(LIBHTS)

$(PROGRAM): $(OBJ)
	$(CXX) $(OBJ) -o $(PROGRAM) $(LIBS) 


-include $(DEP)
FLAGS=$(CPPFLAGS) $(CXXFLAGS)
%.o: %.cpp
	$(CXX) -c  $(FLAGS) $*.cpp
	$(CXX) -MM $(FLAGS) $*.cpp >$*.d




####################################################################################################
## [clean]
# - clean up the directory, use `make clean`
.PHONY: clean
clean:
	$(RM) $(OBJ) $(DEP) $(PREP) $(ASM) $(PROGRAM)

####################################################################################################

## [test]
# - to run unit tests, use `make test`
# - tests are defined in test/Makefile
.PHONY: test
test:
	$(MAKE) -C test

####################################################################################################

