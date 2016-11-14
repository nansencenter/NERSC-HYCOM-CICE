VPATH = .:RCS:TMP:../../../XTRA_MODULES_GETPIVOTS

.SUFFIXES:
.SUFFIXES: .o .F90 .F .H   

#path of MSCPROGS relative to this directory;
MSCP=../../../../../..
#path of MSCPROGS relative to ./TMP directory;
MSCP2=$(MSCP)/..
#name of executable;
TARGET=amsre-piv-hyc2.2.12

# include ../Make.Inc/make.inc
include $(MSCP)/src/Make.Inc/make.inc

# These are the hycom utility library and the conformal mapping library;
# NB: Path is relative to TMP dir
# LIBS:= -L../../../lib/ -lhycnersc $(LIBS)
LIBS:= -L$(MSCP2)/lib/ -lhycnersc -lconfmap $(LIBS)

# This is the hycom utility includes (pre-compiled modules - ie *.mod files)
# NB: Path is relative to TMP dir
# INCLUDE:= -I../ -I../../../include/ $(INCLUDE)
INCLUDE:= -I../ -I$(MSCP2)/include/ $(INCLUDE)


# AIX (xlf) peculiarities. NB: Only necessary for Fortran files
ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGS:=$(subst -,$(subs),$(CPPFLAGS))
endif

# Rules for running cpp and updating files in TMP directory
.F90.o:
	cd ./TMP ; $(CF90) -c $(CPPFLAGS) $(FFLAGS) $(F90FLG) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) -c $(CPPFLAGS) $(FFLAGS) $(F77FLG) $(INCLUDE) -o $*.o ../$<




all : $(TARGET) 

########################################################################
OBJECTS= \
	 m_ncvar_read.o \
	 m_pivotp_micom.o \
	 p_getpivots.o

$(TARGET): $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS) 
#############################################################################

all : $(TARGET) 

# We dont trust mosf yet
install : all 
	mkdir -p $(MSCP)/bin
	cp $(TARGET) $(MSCP)/bin
	# mkdir -p ../../bin
	# cp $(TARGET) ../../bin


clean:
	rm TMP/*.o TMP/*.mod $(TARGET) 

# p_ww3.o : m_pivotp_micom.o m_ncvar_read.o
# m_pivotp_micom.o :  m_spherdist.o 
