FFILES =  topo.f 
OFILES = 
FFLAGS =  -O2 -freeform  -64  -mips4 -r8 -r10000
LDFLAGS =  -64  -mips4 -r8 -r10000
LIBS = -lcomplib.sgimath

.f.o:
	f90 -c  $(FFLAGS) -o $*.o $*.f

OBJS = $(FFILES:.f=.o)

TARGET = topo

$(TARGET): $(OBJS)
	f90 -o $(TARGET) $(LDFLAGS)  $(OBJS) $(LIBS) 
