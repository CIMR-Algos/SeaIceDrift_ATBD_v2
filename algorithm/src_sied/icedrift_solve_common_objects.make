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
# To create some 'common' object files for the solve softwares:
#    - icedrift_prepost.o 
#    - icedrift_model.o
#    - icedrift_io.o
#
# AUTHOR:
#    Thomas Lavergne, met.no/FoU, 22.09.2008
#
# MODIFIED:
# CVS:
# $Id: icedrift_solve_common_objects.make 5642 2012-10-12 15:37:24Z thomasl $
#
###################################################################


#include ../../../../../Makefile.in
include /home/emilyjd/osisaf-hl-sw/Makefile.in
include ../../Makefile.in

COMMON       = $(OSIROOT)/OSI_HL_Ice/common
ICECOMMONINC = $(COMMON)/libicecommon/include
OSINETCDFINC = $(COMMON)/libosinetcdf/include


# It is nice to include yourself and common dir
SELFINC = ./  


# Loader FLAGS required, specify as needed


# Turn on optimization or debugging if required.
OPT =


# Set CFLAGS as needed.
CFLAGS =  $(GLOB_FLG_CC) $(OPT) \
	-I$(SELFINC)   \
	-I$(LRSID_COMMON_INC) \
	-I$(COMMONDEFSINC) \
	-I$(FMUTILINC) \
	-I$(NETCDFINC) \
	-I$(LISTINC) \
	-I$(RSPRODINC) \
	-I$(ICECOMMONINC) \
	-I$(OSINETCDFINC) \
	-I$(USEPROJINC) \
	-I$(PROJINC)   

# List the object files. Important!!! Object files not sourcefiles at
# present as the string substitution is not tested on IRIX yet...
OBJS =  icedrift_prepost.o icedrift_model.o icedrift_io.o icedrift_solve_common.o icedrift_solve_filter.o \
        icedrift_uncertainties.o

# Specify the name of the executable file
# This will be installed properly if make install is executed.
RUNFILE = 


# Specify name of dependency files (e.g. header files)
DEPS = 


all : $(OBJS)


# Specify requirements for the object generation.
$(OBJS): $(DEPS)


# Cleaning and installing commands
clean:
	-rm -f $(OBJS)

distclean: clean
ifdef RUNFILE
	-rm -f $(RUNFILE)
endif

install:
ifdef RUNFILE
    install -d ../bin
	install $(RUNFILE) ../bin/
endif

