###################################################################
# COPYRIGHT: EUMETSAT
#
# PRODUCED BY:
# Norwegian Meteorological Institute (DNMI)
# Research and Development Department
# P.O.BOX 43 - Blindern, N-0313 OSLO, NORWAY
#        
# This SW was developed by DNMI and DMI within the context of the
# Co-operation Agreement for the development of a pilot SAF on
# Ocean and Sea Ice.
###################################################################

###################################################################
#
# TYPE: 
# Makefile
# 
# PURPOSE:
# To create icedrift_solve_simplex
#
# AUTHOR:
# Thomas Lavergne, met.no/FoU, 22.09.2008
#
# MODIFIED:
# ESN, DMI, 02.10.2010: Removed DMI specific part
# $Id: icedrift_solve_simplex.make 12096 2017-07-01 18:57:35Z thomasl $
#
###################################################################

#include ../../../../../Makefile.in
include /home/emilyjd/osisaf-hl-sw/Makefile.in
include ../../Makefile.in


# It is nice to include yourself and common dir
SELFINC = ./  

# Include common ice libs and incs
COMMON       = $(OSIROOT)/OSI_HL_Ice/common
OSINETCDFSRC = $(COMMON)/libosinetcdf/src
OSINETCDFINC = $(COMMON)/libosinetcdf/include
OSINETCDFLIB = $(COMMON)/libosinetcdf/lib
ICECOMMONLIB = $(COMMON)/libicecommon/lib
ICECOMMONINC = $(COMMON)/libicecommon/include
SIMPLEXLIB   = src_simplex


# Loader FLAGS required, specify as needed
LDFLAGS = \
	$(GLOB_FLG_LD) \
	-L$(LRSID_COMMON_LIB) -llrsid \
	-L$(ICECOMMONLIB) -licecommon \
	-L$(OSINETCDFLIB) -losinetcdf \
	-L$(NETCDFLIB)    -lnetcdf \
	-L$(FMUTILLIB)    -lfmutil \
	-L$(PROJLIB)      -lproj \
	-L$(SIMPLEXLIB)   -lsimplex \
	-L$(RSPRODLIB)    -lrsprod \
	-L$(LISTLIB)      -llist \
        -lm


# Turn on optimization or debugging if required.
OPT =

ifeq ($(OSI_SITE),DMI)
OPT = 
endif


# Set CFLAGS as needed.
CFLAGS = $(GLOB_FLG_CC) $(OPT) \
	-I$(SELFINC) \
	-I$(LRSID_COMMON_INC) \
	-I$(COMMONDEFSINC) \
	-I$(FMUTILINC) \
	-I$(USEPROJINC) \
	-I$(ICECOMMONINC) \
	-I$(NETCDFINC) \
	-I$(OSINETCDFINC) \
	-I$(PROJINC) \
	-I$(SIMPLEXLIB) \
	-I$(LISTINC) \
	-I$(RSPRODINC)

# List the object files. Important!!! Object files not sourcefiles at
# present as the string substitution is not tested on IRIX yet...
LOCOBJS = icedrift_solve_simplex.o
EXTOBJS = $(LRSID_TRACKING)/src/icedrift_prepost.o $(LRSID_TRACKING)/src/icedrift_model.o $(LRSID_TRACKING)/src/icedrift_io.o \
          $(LRSID_TRACKING)/src/icedrift_solve_common.o $(LRSID_TRACKING)/src/icedrift_solve_filter.o \
		  $(LRSID_TRACKING)/src/func_summerhold.o
OBJS    = $(LOCOBJS) $(EXTOBJS)  

# Specify the name of the executable file
# This will be installed properly if make install is executed.
RUNFILE = icedrift_solve_simplex


# Specify name of dependency files (e.g. header files)
DEPS = $(LRSID_COMMON_INC)/liblrsid.a


$(RUNFILE): $(OBJS) 
	$(LINKER) $(OBJS) -o $(RUNFILE) $(LDFLAGS)


# Specify requirements for the object generation.
$(OBJS): $(DEPS)


# Cleaning and installing commands
clean:
	-rm -f $(LOCOBJS)

distclean: clean
ifdef RUNFILE
	-rm -f $(RUNFILE)
endif

install:
ifdef RUNFILE
	install -d ../bin
	install $(RUNFILE) ../bin/
endif

