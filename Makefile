# CC = gcc
CC = g++

# compiler flags:
#  -g     add debugging info to the executable 
#  -Wall  turn on most compiler warnings
# CPPFLAGS  = -g -Wall -I /maps/projects/lundbeck/scratch/pfs488/AMOVA/vcfToGlf/htslib/htslib/ -lz
CPPFLAGS  = -g -Wall 
LIBS = -lz -lhts

# HTSSRC := $(CURDIR)/htslib
HTSSRC := /maps/projects/lundbeck/scratch/pfs488/AMOVA/vcfToGlf/htslib


CPPFLAGS += -I$(HTSSRC)
LIBHTS := $(HTSSRC)/libhts.a
LIBS := $(LIBHTS) $(LIBS)



TARGET = vcfReader

all: $(TARGET)

$(TARGET): $(TARGET).cpp
		  $(CC) $(CPPFLAGS) $(LIBS) -o $(TARGET) $(TARGET).cpp

clean:
		  $(RM) $(TARGET)
