###################################################################################################
# ngsAMOVA Makefile
####################################################################################################

default: all

PROGRAM_NAME = ngsAMOVA
PROGRAM_BIN = bin/$(PROGRAM_NAME)

PRINTF_BOLD := $(shell command -v tput > /dev/null && tput bold || echo "")
PRINTF_GREEN := $(shell command -v tput > /dev/null && tput setaf 2 || echo "")
PRINTF_RED := $(shell command -v tput > /dev/null && tput setaf 1 || echo "")
PRINTF_YELLOW := $(shell command -v tput > /dev/null && tput setaf 3 || echo "")
PRINTF_NORMAL := $(shell command -v tput > /dev/null && tput sgr0 || echo "")

# Compiler specs
CXX ?= g++


# No-compile targets
NO_COMPILE = help clean test 

VAL_ADD_LAPACKLIB = -llapack

VAL_ADD_CRYPTOLIB = -lcrypto


####################################################################################################
# [BLOCK START]
# ony run if make will try compiling, i.e. not NO_COMPILE

DEVH = build/dev.h

ifeq ($(filter $(NO_COMPILE),$(MAKECMDGOALS))$(findstring test,$(MAKECMDGOALS)),) #1

$(info ________________________________________________________________________________)
$(info )
$(info [INFO]    Checking for library availability)
$(info )

## ----------------- [LAPACK LIBRARY] ----------------- ##
## -> [lapack availability check]
LAPACK_TRY = $(shell echo 'int main(){}'|$(CXX) -x c++ - -llapack 2>/dev/null -o /dev/null; echo $$?)

ifeq "$(LAPACK_TRY)" "0" #1_0

$(info $(PRINTF_GREEN)[INFO]    -> LAPACK library is available to link$(PRINTF_NORMAL))
THIS_LAPACKLIB = $(VAL_ADD_LAPACKLIB)

# if LAPACK_TRY != 0
else  #1_0

THIS_LAPACKLIB = 

$(info $(PRINTF_YELLOW)[INFO]	  -> LAPACK library is not available to link; will use EIGEN instead$(PRINTF_NORMAL))
endif #1_0


## ----------------- [CRYPTO LIBRARY] ----------------- ##
## -> [cryptolib availability check]
CRYPTO_TRY = $(shell echo 'int main(){}'|$(CXX) -x c++ - -lcrypto 2>/dev/null -o /dev/null; echo $$?)

ifeq "$(CRYPTO_TRY)" "0" #1_1

$(info $(PRINTF_GREEN)[INFO]    -> Crypto library is available to link$(PRINTF_NORMAL))
THIS_CRYPTOLIB = $(VAL_ADD_CRYPTOLIB)

# if CRYPTO_TRY != 0
else  #1_1

THIS_CRYPTOLIB = 

$(info $(PRINTF_YELLOW)[INFO]    -> Crypto library is not available to link; will not use crypto functions$(PRINTF_NORMAL))
endif #1_1

## ----------------- [HTSLIB SOURCE] ------------------ ##
## -> [htslib source check]

#if htslib source is defined
ifdef HTSSRC #1_2

#if hts source is set to systemwide
ifeq ($(HTSSRC),systemwide) #1_2_1

$(info [INFO]    -> HTSSRC set to systemwide; assuming systemwide installation)
HTS_TRY = $(shell echo 'int main(){}'|$(CXX) -x c++ - -lhts 2>/dev/null -o /dev/null; echo $$?)

ifeq "$(HTS_TRY)" "0" #1_2_1_0

$(info $(PRINTF_GREEN)[INFO]    -> Systemwide installation of HTSlib is available to link$(PRINTF_NORMAL))

THIS_LIBHTS := -lhts

else #1_2_1_0

$(info $(PRINTF_RED)[ERROR]   -> Systemwide installation of HTSlib is not available to link$(PRINTF_NORMAL))
$(info $(PRINTF_RED)[ERROR]   -> Please provide the path to the HTSlib source directory or ensure that HTSlib is installed systemwide$(PRINTF_NORMAL))
$(info $(PRINTF_RED)[ERROR]   -> Exiting$(PRINTF_NORMAL))
$(error $(PRINTF_RED)[ERROR]   -> Systemwide installation of HTSlib is not available to link. Please make sure you set the HTSSRC variable correctly. See 'make help' for more information$(PRINTF_NORMAL))

endif #1_2_1_0

else #1_2_1

#if hts source path is given
# Adjust $(HTSSRC) to point to your top-level htslib directory
$(info [INFO]    -> HTSSRC is defined as $(HTSSRC))
CPPFLAGS := -I$(realpath $(HTSSRC))
THIS_LIBHTS := $(realpath $(HTSSRC))/libhts.a

endif #1_2_1

#if htssrc not defined
else #1_2

$(info [INFO]    -> HTSSRC is not defined; using htslib submodule)

#check if htslib submodule exists
ifeq (,$(wildcard $(CURDIR)/htslib)) #1_2_0

$(info $(PRINTF_RED)[ERROR]   -> HTSlib submodule not found$(PRINTF_NORMAL))
$(info $(PRINTF_RED)[ERROR]   -> Exiting$(PRINTF_NORMAL))
$(error $(PRINTF_RED)[ERROR]   -> HTSlib submodule not found. Please make sure you set the HTSSRC variable correctly. See 'make help' for more information$(PRINTF_NORMAL))

else #1_2_0

HTSSRC := $(realpath $(CURDIR)/htslib)
CPPFLAGS := -I$(HTSSRC)
THIS_LIBHTS := $(HTSSRC)/libhts.a

endif #1_2_0

all: .activate_module

endif #1_2

DEV_FLAGS = -g -Wall -Wextra -Wpedantic -Werror
OPTIM_OFF = -O0
OPTIM_ON = -O3

PREV_BUILD_MODE := $(shell grep -oP 'define DEV \K\d' $(DEVH))


ifeq (dev,$(filter dev,$(MAKECMDGOALS))) #1_3

THIS_BUILD_MODE := 1

ifeq (1,$(PREV_BUILD_MODE)) #1_3_1

# PREV_BUILD_MODE is 1
BUILD_MODE_CHANGED := 0

else #1_3_1

# PREV_BUILD_MODE is 0
BUILD_MODE_CHANGED := 1

endif #1_3_1
else #1_3

THIS_BUILD_MODE := 0

ifeq (1,$(PREV_BUILD_MODE)) #1_3_2
# PREV_BUILD_MODE is 1

BUILD_MODE_CHANGED := 1

else #1_3_2
# PREV_BUILD_MODE is 0

BUILD_MODE_CHANGED := 0

endif #1_3_2
endif #1_3




ifeq (1,$(BUILD_MODE_CHANGED)) #1_4
$(info [INFO]    -> Build mode has been changed)
$(info [INFO]    -> Updating $(DEVH))
$(shell sed -i 's/define DEV $(PREV_BUILD_MODE)/define DEV $(THIS_BUILD_MODE)/' $(DEVH))
endif #1_4



ifeq (dev,$(filter dev,$(MAKECMDGOALS))) #1_5

$(info )
$(info ________________________________________________________________________________)
$(info )
$(info [INFO]    Requested compilation in developer/debug mode)
$(info )


OPTIM_FLAGS := $(OPTIM_OFF)
THIS_MODE_FLAGS := $(DEV_FLAGS) $(OPTIM_FLAGS)


else #1_5

$(info )
$(info ________________________________________________________________________________)
$(info )
$(info [INFO]    Requested compilation in release mode)
$(info )

OPTIM_FLAGS := $(OPTIM_ON)
THIS_MODE_FLAGS := $(OPTIM_FLAGS)

endif #1_5


CXXFLAGS ?=
$(info [INFO]    -> CXXFLAGS was "$(CXXFLAGS)")
CXXFLAGS += $(THIS_MODE_FLAGS)
$(info [INFO]    -> Updated CXXFLAGS to "$(CXXFLAGS)")


$(info [INFO]    -> LIBS was "$(LIBS)")
LIBS ?=
LIBS += $(THIS_LIBHTS) $(THIS_CRYPTOLIB) $(THIS_LAPACKLIB) -lz -lm -lbz2 -llzma -lcurl -lpthread
$(info [INFO]    -> Updated LIBS to "$(LIBS)")

$(info )
$(info ________________________________________________________________________________)
$(info )

endif #1

# [BLOCK END]
####################################################################################################


####################################################################################################


VERSION = v0.6

ifneq ($(wildcard .git),)
$(info [INFO]    -> Building from git repository, using git to determine version)
VERSION := $(VERSION)-$(shell git describe --always 2>/dev/null || echo "unknown")
$(info [INFO]    -> Version is determined to be $(VERSION))
endif

VERSIONH = build/version.h

.PHONY: .check_version

.check_version:
	@echo "#define NGSAMOVA_VERSION \"$(VERSION)\"" > $(VERSIONH).tmp
	@if [ ! -f $(VERSIONH) ] || ! cmp -s $(VERSIONH).tmp $(VERSIONH); then \
		mv $(VERSIONH).tmp $(VERSIONH); \
		echo "[INFO]    Updated $(VERSIONH) with version $(VERSION)"; \
	else \
		rm -f $(VERSIONH).tmp; \
	fi



dev: all

all: .check_version $(PROGRAM_BIN)

CXXSRC := $(wildcard src/*.cpp)

# Preprocessed C++ files
PREP := $(CXXSRC:.cpp=.ii)

# Assembly source files
ASM := $(CXXSRC:.cpp=.s)

# Object files
OBJ := $(CXXSRC:.cpp=.o)

# Dependency files
DEP := $(OBJ:.o=.d)

-include $(DEP)

ifeq ($(filter $(NO_COMPILE),$(MAKECMDGOALS))$(findstring test,$(MAKECMDGOALS)),) #1
CXXFLAGS ?=
$(info [INFO]    -> CXXFLAGS was "$(CXXFLAGS)")
ifneq "$(LAPACK_TRY)" "0" 
CXXFLAGS += -DEIGEN_EXTERN_INSTANTIATIONS -DEIGEN_NO_DEBUG -fno-exceptions
endif
$(info [INFO]    -> Updated CXXFLAGS to "$(CXXFLAGS)")

FLAGS = -Ibuild $(CXXFLAGS) $(CPPFLAGS) -Iinclude 

endif #1


# Build info
BUILDH = build/build.h

$(BUILDH): 
	$(info [INFO]    Writing build info to $(BUILDH))
	$(shell echo '#define NGSAMOVA_MAKE_CXX ("$(CXX)")' > $(BUILDH))
	$(shell echo '#define NGSAMOVA_MAKE_LIBS ("$(LIBS)")' >> $(BUILDH))
	$(shell echo '#define NGSAMOVA_MAKE_FLAGS ("$(FLAGS)")' >> $(BUILDH))
	$(shell echo '#define NGSAMOVA_MAKE_HTSSRC ("$(HTSSRC)")' >> $(BUILDH))
	$(shell echo '#define NGSAMOVA_MAKE_CXXFLAGS ("$(CXXFLAGS)")' >> $(BUILDH))
	$(shell echo '#define NGSAMOVA_MAKE_CPPFLAGS ("$(CPPFLAGS)")' >> $(BUILDH))
	$(shell if [ "$(LAPACK_TRY)" = "0" ]; then echo '#define NGSAMOVA_USE_LAPACK 1'; else echo '#define NGSAMOVA_USE_LAPACK 0'; fi >> $(BUILDH))
	$(shell if [ "$(CRYPTO_TRY)" = "0" ]; then echo '#define NGSAMOVA_USE_CRYPTO 1'; else echo '#define NGSAMOVA_USE_CRYPTO 0'; fi >> $(BUILDH))


$(PROGRAM_BIN): $(OBJ) $(DEVH) $(VERSIONH) $(BUILDH)
	@echo "________________________________________________________________________________"
	@echo ""
	@echo "$(PRINTF_YELLOW)[INFO]    -> Finishing up$(PRINTF_NORMAL)"
	mkdir -p bin/
	$(CXX) -o $@ $(OBJ) $(LIBS) 
	@echo "________________________________________________________________________________"
	@echo ""
	@printf "$(PRINTF_BOLD)$(PRINTF_GREEN)"
	@echo "[FINISHED]    $(PROGRAM_NAME) is now ready to use!"
	@echo "              Full path to the program: $(CURDIR)/$(PROGRAM_BIN)"
	@echo ""
	@echo "              To get started, run:"
	@echo "              $(CURDIR)/$(PROGRAM_BIN) -h"
	@echo "              or:"
	@echo "              ./$(PROGRAM_BIN) -h"
	@echo ""
	@printf "$(PRINTF_NORMAL)"

%.o: %.cpp $(DEVH) $(VERSIONH) $(BUILDH) 
	@echo ""
	@echo "________________________________________________________________________________"
	@echo ""
	@echo "$(PRINTF_YELLOW)[INFO]    -> Compiling and generating dependencies for $*.cpp$(PRINTF_NORMAL)"
	$(CXX) -MMD -MP -c $(FLAGS) $< -o $@

####################################################################################################
## [install]
# - install the program

PREFIX ?= /usr/local
EXEC_PREFIX    ?= $(PREFIX)
BINDIR         ?= $(EXEC_PREFIX)/bin

INSTALL         = install
INSTALL_DIR     = $(INSTALL) -dm0755
INSTALL_PROGRAM = $(INSTALL) -m0755

install: $(PROGRAM_BIN)
	@echo "________________________________________________________________________________"
	@echo ""
	@echo "[INFO]    Installing $(PROGRAM_NAME) to $(BINDIR)"
	@echo ""
	@echo "________________________________________________________________________________"
	$(INSTALL_DIR) $(DESTDIR)$(BINDIR)
	$(INSTALL_PROGRAM) $(PROGRAM_BIN) $(DESTDIR)$(BINDIR)/$(PROGRAM_NAME)
	@echo ""
	@echo "[INFO]    $(PROGRAM_NAME) has been installed to $(BINDIR)"

####################################################################################################
## [uninstall]
# - uninstall the program

uninstall:
	@echo "________________________________________________________________________________"
	@echo ""
	@echo "[INFO]    Uninstalling $(PROGRAM_NAME) from $(BINDIR)"
	@echo ""
	@echo "________________________________________________________________________________"
	$(RM) $(DESTDIR)$(BINDIR)/$(PROGRAM_NAME)
	@echo ""
	@echo "[INFO]    $(PROGRAM_NAME) has been uninstalled from $(BINDIR)"

####################################################################################################
## [clean]
# - clean up the directory

.PHONY: clean
clean:
	@echo "________________________________________________________________________________"
	@echo ""
	@echo "[INFO]    Cleaning up the directory"
	@echo ""
	@echo "________________________________________________________________________________"
	$(RM) $(OBJ) $(DEP) $(PREP) $(ASM) $(VERSIONH) $(BUILDH) $(PROGRAM_BIN) vgcore.* tests/testwd/*
	@echo ""
	
####################################################################################################
## [test]
# 

.PHONY: test test% test-%

test:
	@echo "________________________________________________________________________________"
	@echo ""
	@echo "[INFO]    Running all unit tests"
	@echo ""
	@echo "________________________________________________________________________________"
	@echo ""
	/bin/bash tests/runTests.sh -e $(PROGRAM_BIN) -t regular -d tests/data -r tests/reference -w tests/testwd 

test%:
	@echo "________________________________________________________________________________"
	@echo ""
	@echo "[INFO]    Running unit test for test $*"
	@echo ""
	@echo "________________________________________________________________________________"
	@echo ""
	/bin/bash tests/runTests.sh -e $(PROGRAM_BIN) -t regular -d tests/data -r tests/reference -w tests/testwd  -x $*

test-%:
	@echo "________________________________________________________________________________"
	@echo ""
	@echo "[INFO]    Running all unit tests with test type $*"
	@echo "" 
	@echo "________________________________________________________________________________"
	/bin/bash tests/runTests.sh -e $(PROGRAM_BIN) -t $* -d tests/data -r tests/reference -w tests/testwd 

####################################################################################################
## [.activate_module]

.PHONY: .activate_module

.activate_module:
	@echo "________________________________________________________________________________"
	@echo ""
	@echo "[INFO]    Activating HTSlib submodule"
	@echo ""
	@if ! command -v git &> /dev/null; then \
	echo "$(PRINTF_RED)[ERROR]   Git is not installed! Cannot update the submodules. Please make sure you have git installed and try again.$(PRINTF_NORMAL)"; \
	exit 1; \
	fi
	git submodule update --init --recursive
	$(MAKE) -C $(HTSSRC)
	@echo "[INFO]	-> HTSlib submodule is now activated"
	@echo ""
	@echo "________________________________________________________________________________"

####################################################################################################
## [help]

.PHONY: help
help:
	@echo ""
	@echo "----------------------------------------"
	@echo " Program: $(PROGRAM_NAME)"
	@echo " Version: $(VERSION)"
	@echo " License: GNU GPLv3.0"
	@echo "----------------------------------------"
	@echo ""
	@echo " Usage:"
	@echo "   make [target] [FLAG=value...]"
	@echo ""
	@echo " Targets:"
	@echo "   help ------ Print this help message"
	@echo "   dev ------- Compile in developer/debug mode (activates flags: -g -Wall -O0)"
	@echo "   clean ----- Clean up the directory"
	@echo "   test ------ Run unit tests"
	@echo "   test% ----- Run unit test for test % (e.g. test1)"
	@echo "   test-% ---- Run all unit tests with test type % (e.g. test-vg for valgrind test)"
	@echo "   install --- Install the program"
	@echo "   uninstall - Uninstall the program"
	@echo ""
	@echo " Flags:"
	@echo "   HTSSRC : Specifies the source of HTSlib"
	@echo "     (empty) ---------- Use the HTSlib submodule [default]"
	@echo "     systemwide ------- Use the systemwide HTSlib installation"
	@echo "     /path/to/htslib -- Use the HTSlib installation at /path/to/htslib"
	@echo ""
	@echo " Examples:"
	@echo "   make"
	@echo "     Compile in release mode using HTSlib submodule"
	@echo "   make HTSSRC=/path/to/htslib"
	@echo "     Compile in release mode using /path/to/htslib"
	@echo "   make dev HTSSRC=systemwide"
	@echo "     Compile in developer mode using the systemwide HTSlib installation"
	@echo "   make install PREFIX=$(HOME)/.local"
	@echo "     Install the program to $(HOME)/.local/bin"
	@echo "   sudo make install                      "
	@echo "     Install the program to /usr/local/bin"
	@echo "   sudo make uninstall"
	@echo "     Uninstall the program from /usr/local/bin"
	@echo ""
	@echo " Note: If no values are provided for HTSSRC, CXX, CXXFLAGS, or LIBS, defaults will be used."
	@echo " Your default values are:"
	@echo "   - HTSSRC: '$(HTSSRC)'"
	@echo "   - CXX: '$(CXX)'"
	@echo "   - CXXFLAGS: '$(CXXFLAGS)'"
	@echo "   - LIBS: '$(LIBS)'"
	@echo ""

####################################################################################################
