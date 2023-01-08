CXX ?= g++

#  -g     add debugging info to the executable 
#  -Wall  turn on most compiler warnings
CXXFLAGS  := -g -Wall
LIBS = -lz -lm -lbz2 -llzma -lcurl -lpthread


CRYPTO_TRY=$(shell echo 'int main(){}'|$(CXX) -x c++ - -lcrypto 2>/dev/null -o /dev/null; echo $$?)
ifeq "$(CRYPTO_TRY)" "0"
$(info Crypto library is available to link; adding -lcrypto to LIBS)
LIBS += -lcrypto
else
$(info Crypto library is not available to link; will not use -lcrypto)
endif

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

.PHONY: .activate_module test

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

clean:
	$(RM) *.o *.d $(PROGRAM)

test:
	rm -rvf test/testwd;
	mkdir -pv test/testwd;
	./ngsAMOVA -doEM 1 -doAMOVA 3 -doTest 0 -in test/test_s9_d1_1K.vcf -isSim 1 -minInd 2 -printMatrix 0  -doDist 1 -maxIter 100 -nThreads 0 -tole 1e-10  --hascolnames 1 -m test/metadata_with_header.tsv -out test/testwd/test_s9_d1_1K
	bash -c "diff <(cut -d, -f4-12 test/testwd/test_s9_d1_1K.sfs.csv)  test/ref/test_s9_d1_1K.sfs.csv";

# ./ngsAMOVA -in test/test_s9_d1_1K.vcf -m test/metadata_with_header.tsv -doDist 1 -doAMOVA 3 -doEM 1 -isSim 1 -out test/testwd/test_s9_d1_1K.vcf_testput --formula "Individual~Population";
