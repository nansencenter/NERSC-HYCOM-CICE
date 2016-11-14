VPATH = ..:Build

.SUFFIXES:
.SUFFIXES: .o .F90 .F .f


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

.F.o:
	cd Build/ ;$(CF90) $(INCLUDE) $(CPPFLAGS) $(FFLAGS) $(F77FLG) -c  ../$<
.F90.o:
	cd Build/ ;$(CF90) $(INCLUDE) $(CPPFLAGS) $(FFLAGS) $(F90FLG) -c  ../$<
.f.o:
	cd Build/ ;$(CF90) -c $(FFLAGS) $(F77FLG)  ../$<





TARGET=nestport
TARGET2=nestini
TARGET3=checknest
TARGET4=nestpivots
all : $(TARGET) $(TARGET2) $(TARGET3) $(TARGET4)


#################################################################################
OBJECTS = \
		  p_nestport.o 

$(TARGET): $(OBJECTS) 
	cd Build/ ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS)
#################################################################################
OBJECTS2 = \
		  m_parse_blkdat.o   \
		  p_nestini.o 
$(TARGET2): $(OBJECTS2) 
	cd Build/ ; $(LD) $(LINKFLAGS) -o ../$(TARGET2) $(OBJECTS2) $(LIBS)
##################################################################################
OBJECTS3 = \
		  m_parse_blkdat.o   \
		  m_nearestpoint.o   \
		  mod_nest_io.o   \
		  p_check_nest.o  \
		  dayfor.o

$(TARGET3): $(OBJECTS3) 
	cd Build/ ; $(LD) $(LINKFLAGS) -o ../$(TARGET3) $(OBJECTS3) $(LIBS)
##################################################################################
OBJECTS4 = \
		  m_nearestpoint.o   \
		  mod_nest_io.o   \
		  m_parse_blkdat.o   \
		  p_nestpivots.o \
		  dayfor.o

$(TARGET4): $(OBJECTS4) 
	cd Build/ ; $(LD) $(LINKFLAGS) -o ../$(TARGET4) $(OBJECTS4) $(LIBS)
##################################################################################

clean:
	rm Build/*.o Build/*.mod Build/*.f Build/*.f90  $(TARGET) $(TARGET2) $(TARGET3)  $(TARGET4)


install: all
	mkdir -p ../../bin_setup/ 
	cp $(TARGET) $(TARGET2) $(TARGET3) $(TARGET4) ../../bin_setup/
