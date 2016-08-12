#
# Makefile for archive processing
#
# --- usage: make hycomproc ARCH=sun >& Make_hycomproc
#
# --- Tunable parameters in ../../config/$(ARCH)_setup
#
# --- assumes netcdf module and library are in this directory
#

.SUFFIXES: 
.SUFFIXES: .c .F .f .o

.F:
	@echo "Must have an explicit rule for" $*
	@exit 1
.f:
	@echo "Must have an explicit rule for" $*
	@exit 1
.c:
	@echo "Must have an explicit rule for" $*
	@exit 1

include ../../config/$(ARCH)_setup

#
# ---------------------------------------------------------------------
# Standard part of the Makefile starts here
# ---------------------------------------------------------------------
#

all:
	/bin/csh Make_all.com

clean:
	/bin/rm -f *.o *.a *.mod M*log *2archv *2d *2data *2t *3z *sf *sfl *sfz *sig2a *restart *hv *ncombc *_inc *gen *ncdf

MODS     =	mod_plot.o mod_xc.o mod_za.o               wtime.o
MODB     =	mod_plot.o mod_xc.o mod_za.o mod_zb.o      wtime.o
MODZ     =	mod_plot.o mod_xc.o mod_za.o mod_zb.o      wtime.o mod_ppsw.o
MODR     =	mod_plot.o mod_xc.o mod_za.o mod_restart.o wtime.o

DATAX_OBJ = 	bigrid.o           blkin.o    \
		buoflx.o                      \
		extrct.o  extrot.o fordate.o  \
                getdat.o  getdtm.o getdepth.o \
		          indxi.o  indxj.o    \
                pakk.o             psmoo.o    \
		          zebra.o  zh.o       \
                dum_gks.o

DATAH_OBJ = 	$(DATAX_OBJ) horout.o
DATAR_OBJ = 	$(DATAX_OBJ) putdat.o
DATAZ_OBJ = 	$(DATAH_OBJ) layer2z.o mixlay.o

NCDFH_OBJ = 	$(DATAX_OBJ) horout_nc.o
NCDFZ_OBJ = 	$(NCDFH_OBJ) layer2z.o mixlay.o

ARCHM_OBJ = 	bigrid.o blkin.o  \
		extrct.o extrot.o \
                getdat.o getdtm.o getdepth.f putdat.o \
		indxi.o  indxj.o  \
		zebra.o  zh.o     dum_gks.o \
		pakk.o   
# 		pakk_micom.o   

ARCHV_OBJ = 	bigrid.o blkin.o  \
		extrct.o extrot.o fordate.o \
                getdat.o getdtm.o getdepth.o putdat.o \
		indxi.o  indxj.o  pakk.o     psmoo.o  \
		zebra.o   zh.o    dum_gks.o 

ARCHT_OBJ = 	          zh.o    dum_gks.o 

ARCHV_SRC = 	hycomarchv.f micomarchv.f \
                conv_archv.f trim_archv.f mrgl_archv.f cooperhaines_archv.f \
                hybgen_archv.f     \
                hybgen_archv_old.f \
                ncoda_archv.f     \
                ncoda_archv_inc.f \
                remap_archv.f remapi_archv.f \
                field2data.f \
                archt2archv.f \
                archv2data2d.f archv2data2t.f archv2data3z.f \
                archv2datasf.f archv2datasfl.f archv2datasfz.f \
                archv2restart.f putdat.f

#
# --- executables.
#

archv2data2d:       $(MODB) archv2data2d.o  $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o archv2data2d  archv2data2d.o  $(MODB) $(DATAH_OBJ)

archm2data2d:       $(MODB) archm2data2d.o  $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o archm2data2d  archm2data2d.o  $(MODB) $(DATAH_OBJ)

archv2data2t:       $(MODB) archv2data2t.o  $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o archv2data2t  archv2data2t.o  $(MODB) $(DATAH_OBJ)

archm2data2t:       $(MODB) archm2data2t.o  $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o archm2data2t  archm2data2t.o  $(MODB) $(DATAH_OBJ)

archv2data3z:       $(MODZ) archv2data3z.o  $(DATAZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2data3z  archv2data3z.o  $(MODZ) $(DATAZ_OBJ)

archm2data3z:       $(MODZ) archm2data3z.o  $(DATAZ_OBJ)
	$(LD) $(LDFLAGS) -o archm2data3z  archm2data3z.o  $(MODZ) $(DATAZ_OBJ)

archv2datasf:       $(MODB) archv2datasf.o  $(DATAZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2datasf  archv2datasf.o  $(MODB) $(DATAZ_OBJ)

archv2datasfl:      $(MODB) archv2datasfl.o $(DATAZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2datasfl archv2datasfl.o $(MODB) $(DATAZ_OBJ)

archv2datasfz:      $(MODB) archv2datasfz.o $(DATAZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2datasfz archv2datasfz.o $(MODB) $(DATAZ_OBJ)

field2data:         $(MODB) field2data.o    $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o field2data    field2data.o    $(MODB) $(DATAH_OBJ)

field2ncdf:         $(MODB) field2data.o    $(NCDFH_OBJ)
	$(LD) $(LDFLAGS) -o field2ncdf    field2data.o    $(MODB) $(NCDFH_OBJ) libnetcdf.a $(EXTRANCDF)

archv2ncdf2d:       $(MODB) archv2data2d.o  $(NCDFH_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdf2d  archv2data2d.o  $(MODB) $(NCDFH_OBJ) libnetcdf.a $(EXTRANCDF)

archm2ncdf2d:       $(MODB) archm2data2d.o  $(NCDFH_OBJ)
	$(LD) $(LDFLAGS) -o archm2ncdf2d  archm2data2d.o  $(MODB) $(NCDFH_OBJ) libnetcdf.a $(EXTRANCDF)

archv2ncdf2t:       $(MODB) archv2data2t.o  $(NCDFH_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdf2t  archv2data2t.o  $(MODB) $(NCDFH_OBJ) libnetcdf.a $(EXTRANCDF)

archv2ncdf3z:       $(MODZ) archv2data3z.o  $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdf3z  archv2data3z.o  $(MODZ) $(NCDFZ_OBJ) libnetcdf.a $(EXTRANCDF)

archm2ncdf3z:       $(MODZ) archm2data3z.o  $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archm2ncdf3z  archm2data3z.o  $(MODZ) $(NCDFZ_OBJ) libnetcdf.a $(EXTRANCDF)

archv2ncdfsf:       $(MODB) archv2datasf.o  $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdfsf  archv2datasf.o  $(MODB) $(NCDFZ_OBJ) libnetcdf.a $(EXTRANCDF)

archv2ncdfsfl:      $(MODB) archv2datasfl.o $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdfsfl archv2datasfl.o $(MODB) $(NCDFZ_OBJ) libnetcdf.a $(EXTRANCDF)

archv2ncdfsfz:     $(MODB) archv2datasfz.o $(NCDFZ_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncdfsfz archv2datasfz.o $(MODB) $(NCDFZ_OBJ) libnetcdf.a $(EXTRANCDF)

archv2ncombc:       $(MODB) archv2ncombc.o  $(DATAH_OBJ)
	$(LD) $(LDFLAGS) -o archv2ncombc  archv2ncombc.o  $(MODB) $(DATAH_OBJ)

archv2restart:      $(MODR) archv2restart.o $(DATAX_OBJ)
	$(LD) $(LDFLAGS) -o archv2restart archv2restart.o $(MODR) $(DATAX_OBJ)

archm2restart:      $(MODR) archm2restart.o $(DATAX_OBJ)
	$(LD) $(LDFLAGS) -o archm2restart archm2restart.o $(MODR) $(DATAX_OBJ)

cooperhaines_archv: $(MODS) cooperhaines_archv.o  $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o cooperhaines_archv cooperhaines_archv.o $(MODS) $(ARCHV_OBJ)

archt2archv:        $(MODS) archt2archv.o  $(ARCHT_OBJ)
	$(LD) $(LDFLAGS) -o archt2archv archt2archv.o $(MODS) $(ARCHT_OBJ)

hycomarchv:         $(MODS) hycomarchv.o   $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o hycomarchv  hycomarchv.o $(MODS) $(ARCHV_OBJ)

micomarchv:         $(MODS) micomarchv.o   $(ARCHM_OBJ)
	$(LD) $(LDFLAGS) -o micomarchv  micomarchv.o $(MODS) $(ARCHM_OBJ)

conv_archv:         $(MODS) conv_archv.o   $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o conv_archv  conv_archv.o $(MODS) $(ARCHV_OBJ)

hybgen_archv:       $(MODS) hybgen_archv.o     $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o hybgen_archv      hybgen_archv.o     $(MODS) $(ARCHV_OBJ)

hybgen_archv_old:   $(MODS) hybgen_archv_old.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o hybgen_archv_old  hybgen_archv_old.o $(MODS) $(ARCHV_OBJ)

mrgl_archv:         $(MODS) mrgl_archv.o   $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o mrgl_archv  mrgl_archv.o $(MODS) $(ARCHV_OBJ)

trim_archv:         $(MODS) trim_archv.o   $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o trim_archv  trim_archv.o $(MODS) $(ARCHV_OBJ)

ncoda_archv:        $(MODS) ncoda_archv.o      $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o ncoda_archv     ncoda_archv.o     $(MODS) $(ARCHV_OBJ)

ncoda_archv_inc:    $(MODS) ncoda_archv_inc.o  $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o ncoda_archv_inc ncoda_archv_inc.o $(MODS) $(ARCHV_OBJ)

remap_archv:        $(MODS) remap_archv.o  $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o remap_archv remap_archv.o $(MODS) $(ARCHV_OBJ)

remapi_archv:       $(MODS) remapi_archv.o  $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o remapi_archv remapi_archv.o $(MODS) $(ARCHV_OBJ)

remap_ncombc:       $(MODS) remap_ncombc.o $(ARCHV_OBJ)
	$(LD) $(LDFLAGS) -o remap_ncombc remap_ncombc.o $(MODS) $(ARCHV_OBJ)

restart2archv:      $(MODR) restart2archv.o $(DATAR_OBJ)
	$(LD) $(LDFLAGS) -o restart2archv restart2archv.o $(MODR) $(DATAR_OBJ)

restart2archv_ssh:  $(MODR) restart2archv_ssh.o   $(DATAR_OBJ)
	$(LD) $(LDFLAGS)  -o restart2archv_ssh   restart2archv_ssh.o   $(MODR) $(DATAR_OBJ)

restart2archv_sig2a: $(MODR) restart2archv_sig2a.o $(DATAR_OBJ)
	$(LD) $(LDFLAGS)  -o restart2archv_sig2a restart2archv_sig2a.o $(MODR) $(DATAR_OBJ)

#
# --- explicit dependencies for each .[fF] file.
#

archv2data2d.o:  archv2data2d.f  mod_plot.o mod_za.o
archm2data2d.o:  archm2data2d.f  mod_plot.o mod_za.o
archv2data3z.o:  archv2data3z.f  mod_plot.o mod_za.o
archm2data3z.o:  archm2data3z.f  mod_plot.o mod_za.o
archv2datasf.o:  archv2datasf.f  mod_plot.o mod_za.o
archv2datasfl.o: archv2datasfl.f mod_plot.o mod_za.o
archv2datasfz.o: archv2datasfz.f mod_plot.o mod_za.o
archv2ncombc.o:  archv2ncombc.f  mod_plot.o mod_za.o
archv2restart.o: archv2restart.f mod_plot.o mod_za.o mod_restart.o
archm2restart.o: archm2restart.f mod_plot.o mod_za.o mod_restart.o

hycomarchv.o:    hycomarchv.f    mod_plot.o mod_za.o
micomarchv.o:    micomarchv.f    mod_plot.o mod_za.o
mrgl_archv.o:    mrgl_archv.f    mod_plot.o mod_za.o
trim_archv.o:    trim_archv.f    mod_plot.o mod_za.o
remap_archv.o:   remap_archv.f   mod_plot.o mod_za.o

hybgen_archv.o:     hybgen_archv.f     mod_plot.o mod_za.o
hybgen_archv_old.o: hybgen_archv_old.f mod_plot.o mod_za.o

ncoda_archv.o:      ncoda_archv.f      mod_plot.o mod_za.o
ncoda_archv_inc.o:  ncoda_archv_inc.f  mod_plot.o mod_za.o

restart2archv.o: restart2archv.f mod_plot.o mod_za.o mod_restart.o

bigrid.o:   bigrid.f      mod_plot.o
blkin.o:    blkin.f
buoflx.o:   buoflx.f
dum_gks.o:  dum_gks.f
extrct.o:   extrct.f
extrot.o:   extrot.f
fordate.o:  fordate.f
getdat.o:   getdat.f      mod_plot.o mod_za.o
getdepth.o: getdepth.f    mod_plot.o mod_za.o
getdtm.o:   getdtm.f      mod_plot.o mod_za.o
horout.o:   horout.f      mod_plot.o mod_xc.o mod_zb.o
horout_nc.o: horout_nc.f  mod_plot.o mod_xc.o mod_zb.o
indxi.o:    indxi.f       mod_plot.o
indxj.o:    indxj.f       mod_plot.o
layer2z.o:  layer2z.f
mixlay.o:   mixlay.f
pakk.o:     pakk.f
prtmsk.o:   prtmsk.f
psmoo.o:    psmoo.f       mod_plot.o
putdat.o:   putdat.f      mod_plot.o mod_za.o
wtime.o:    wtime.F
zebra.o:    zebra.f

mod_plot.o: mod_plot.F
mod_ppsw.o: mod_ppsw.F
mod_xc.o:   mod_xc.F   
mod_za.o:   mod_za.F   mod_xc.o
mod_zb.o:   mod_zb.F

zh.o:       zh.F zh_sun.f zh_t3e.f zh_ymp.f
