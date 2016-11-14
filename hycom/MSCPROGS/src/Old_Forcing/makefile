VPATH= .:TMP

.SUFFIXES:
.SUFFIXES: .o .f90 .f



include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc  $(LIBS)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../include/ $(INCLUDE)


.f.o:
	cd TMP; $(CF77)  -c $(F77FLG)  $(FFLAGS) $(INCLUDE) -o $*.o ../$<

OBJS =	\
   cairt.o ccloud.o chumid.o csst.o cuv.o cwind.o cwstde.o\
	dmpdat.o fillup.o inclcc.o intpol.o lcoads.o lecmw2.o\
	lecmwf.o lecmwv.o lherov.o newforc.o newold.o posks2.o\
	rHeRo.o recmw2.o recmwv.o rgrsp.o rmask.o  \
	rotate.o tecall.o spherdist.o dmpwin.o rcoads.o \
	recmwf.o lamip.o mixairt.o lcoads3.o lecmw3.o rstat.o \
	fillall.o tec1.o tecij.o mercator.o \
	filldepth.o shap2d.o dmp2dtec.o crain.o lrain.o rrain.o \
	cice.o rwalsh.o lwalsh.o intpol2.o


TARGET=old_newf

$(TARGET) :	$(OBJS) 
	cd TMP ; $(CF77)  -o ../$(TARGET) $(LINKFLAGS) $(OBJS) $(LIBS)

all : $(TARGET)

install : all
	mkdir -p ../../bin_setup/
	cp $(TARGET) ../../bin_setup/

clean :
	rm TMP/*.o $(TARGET)

