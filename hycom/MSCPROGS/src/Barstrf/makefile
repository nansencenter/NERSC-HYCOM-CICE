VPATH = .:RCS:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .F 


include ../Make.Inc/make.inc


# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc $(LIBS)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../include/ $(INCLUDE)


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


TARGET1 = hyc2bastrf
all: $(TARGET1)

#############################################################################
OBJECTS1=  \
		  m_handle_err.o  \
		  m_strmf_eval.o  \
		  p_hyc2bastrf.o  

$(TARGET1):  $(OBJECTS1) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET1) $(OBJECTS1) $(LIBS) 
#############################################################################

install : all
	mkdir -p ../../bin
	cp $(TARGET1) ../../bin


clean:
	rm ./TMP/*.o TMP/*.mod TMP/*.f TMP/*.f90 $(TARGET1)
