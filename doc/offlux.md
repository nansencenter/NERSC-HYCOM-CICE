[TOC]

# Prerequisites
	You need to have your efflux [ab] files in "../relax/XXX/offlux.[ab]" in which XXX   
	denotes the experiment number for example 010. 

	in "blkdat.input" file, you have two flags need to be properly set: (1) nrdflg; 	(2) flxoff.

	You should produce the net shortwave radiation and net long wave radiation fields as input available at "../force/synoptic/XXX", 
	here I assumed that our current directory is experiment folder, i.e. expt_X.XX. These extra fields are generated
	if you fire the following command on terminal:
        ../bin/atmo_synoptic.sh erai+all 2007-01-02T00:00:00 2007-12-31T00:00:00

# Values of nrdflg

	0: no net radiation fluxes are involved.
	1: using Bignami 1995 parameterisation or net longwave radiation and the rest from ERA-I. With this value
	   of nrdflg, the net longwave radiation is calculated inside HYCOM. Note that here, we use downwelling shortwave radiation as input and 
	   albedo effect is accounted inside HYCOM-CICE.
        2: 
        3: Old NERSC parameterisation for longewave radiation is perfomed. My trick is to divide the net heat flux into two terms: 
	   SST dependent and not-dependent terms. The SST-independent term is provided as input field, and the SST dependent is handelled using the model SST
           using Budyko (1974) formulation. We use downwelling shortwave radiation from ERA-I.
        4: All radiation forcing are provided from ERA-I.
        5: We use net shortwave radiation and downwelling longwave radiation from ERA-I.
        6: We use net shortwave radiation from ERA-I and downwelling longwavee radiation using Bignami (1995) formulation.
        7: The same as 6, except we use time-variant representation for the ofllux. This flag need to be used when "flxoff" 
	   is set to 1. 

