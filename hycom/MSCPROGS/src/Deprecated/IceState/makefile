VPATH = .:RCS:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h  

#ARCH=x86_64
#ARCH=powerpc
ARCH=$(shell uname -p) #works with gmake
include make.$(ARCH)


# Rules for running cpp and updating files in TMP directory
.H.h:
	rm -f ./TMP/$*.h
	cat  $*.H  > ./TMP/$*.h


.F90.o:
	rm -f ./TMP/$*.f90
	cat  $*.F90 | $(CPP) $(CPPFLAGS)  > ./TMP/$*.f90
	cd ./TMP ; $(CF90) -c $(FFLAGS) $(F90FLG) $(INCLUDE) -o $*.o $*.f90  

.F.o:
	rm -f ./TMP/$*.f
	cat  $*.F | $(CPP) $(CPPFLAGS)  > ./TMP/$*.f
	cd ./TMP ; $(CF77) -c -col72 $(FFLAGS) $(F77FLG) $(INCLUDE) -o $*.o $*.f  


########################################################################
# TARGET1 : fully fledged monster...  Generates section plots from 
#           section and infiles. The spaghetti code should be split up ...
#
#
# TARGET2 : Finds section intersections, and dumps these to a datafile 
#           which can be read by other programs.

TARGET = icestate2ice


#############################################################################
OBJECTS=  mod_xc.o mod_za.o mod_year_info.o p_icestate2ice.o wtime.o


$(TARGET): $(INC2) $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET) $(OBJECTS) $(LIBS) 
#############################################################################


clean:
	cd ./TMP ; rm *.f  *.o *.f90 *.h *.mod

