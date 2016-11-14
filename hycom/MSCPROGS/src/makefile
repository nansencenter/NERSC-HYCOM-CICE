

NERSCLIB=\
Nersclib \

SUBDIRS=\
Average \
Barstrf \
Conf_grid/Code/ \
ConfmapRoutines \
Copymem \
Curviint \
DProfile \
DateTools \
Ensstat \
ExtractNC2D \
ExtractNC3D \
FindLayer \
GP \
GPdens \
GridToLL \
Hyc2proj \
IceDrift \
Idealized_Grid \
MkEnsemble \
Multiproc_New \
NORSEXClim \
Nest_Offline-2.1 \
Nestbat \
Nesting-2.2 \
ObsCompare \
Old_Forcing \
Old_Levitus \
Relax \
RelaxToNetCDF \
River_Forcing \
SSHFromState \
Section \
Synoptic_Forcing-2.1.07 \
Synoptic_Forcing-2.2.12 \
Synoptic_Forcing-2.2.37 \
Synoptic_Forcing-2.2.98 \
TRIP \
Tides_CSR \
ZONAL  \
Model_input-2.2.12 \
Model_input-2.2.37
#Tides_FES/nersc_src/ \


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

