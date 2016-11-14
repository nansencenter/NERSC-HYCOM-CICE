VPATH = .:Build/

.SUFFIXES:
.SUFFIXES: .o .F90 .F 


include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc -lconfmap $(LIBS)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../include/ $(INCLUDE)



# AIX (xlf) peculiarities. NB: Only necessary for Fortran files
ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGS:=$(subst -,$(subs),$(CPPFLAGS))
endif




# Rules for running cpp and updating files 
.F.o:
	cd Build/ ; $(CF90)  $(CPPFLAGS) $(INCLUDEDIR) $(FFLAGS) $(F77FLG) $(INCLUDE) -c ../$<

.F90.o:
	cd Build/ ; $(CF90)  $(CPPFLAGS) $(INCLUDEDIR) $(FFLAGS) $(F90FLG) $(INCLUDE) -c ../$<


TARGET=rivers
TARGET2=biorivers


#################################################################################
OBJECTS=\
		  m_handle_err.o  \
		  mod_rivers.o  \
		  p_rivers.o

$(TARGET): $(OBJECTS) 
	cd Build ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS)
#################################################################################
OBJECTS2=\
		  m_handle_err.o  \
		  mod_biorivers.o  \
		  p_biorivers.o

$(TARGET2): $(OBJECTS2) 
	cd Build ; $(LD) $(LINKFLAGS) -o ../$(TARGET2) $(OBJECTS2) $(LIBS)
#################################################################################

all: $(TARGET) $(TARGET2) 

install : all
	mkdir -p ../../bin_setup/
	cp $(TARGET) $(TARGET2) ../../bin_setup/

clean:
	rm Build/*.o Build/*.mod Build/*.f Build/*.f90  $(TARGET)

