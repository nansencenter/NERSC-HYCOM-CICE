

NERSCLIB=\
Nersclib \

SUBDIRS=\
Average \
Barstrf \
Conf_grid/Code/ \
ConfmapRoutines \
Curviint \
DProfile \
DateTools \
Ensstat \
ExtractNC2D \
ExtractNC3D \
FindLayer \
GridToLL \
Hyc2proj \
IceDrift \
Idealized_Grid \
MkEnsemble \
Relax \
RelaxToNetCDF \
River_Forcing \
SSHFromState \
Section \
TRIP \
ZONAL  \
Model_input-2.2.12 \
Model_input-2.2.37


.PHONY: subdirs $(SUBDIRS)


lib:
		$(MAKE) -C Nersclib/ install 

     
subdirs: $(SUBDIRS) 
     
$(SUBDIRS): lib
		 $(MAKE) -C $@ all


all: lib subdirs 

clean:
		 for dir in $(SUBDIRS) $(NERSCLIB) ; do \
			$(MAKE) -C $$dir clean;\
		 done
		 rm -f ../bin/* ../bin_setup/* ../include/*.mod ../lib/*.a

install: all lib
		 for dir in $(SUBDIRS); do \
			$(MAKE) -C $$dir install;\
		 done

