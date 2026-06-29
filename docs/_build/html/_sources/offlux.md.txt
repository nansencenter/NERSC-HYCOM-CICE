[TOC]

# Prerequisites
	You need to have your efflux [ab] files in "../relax/XXX/offlux.[ab]" in which XXX   
	denotes the experiment number for example 010. 

	in "blkdat.input" file, you have two flags need to be properly set: (1) nrdflg; 	(2) flxoff.

	You should produce the net shortwave radiation and net long wave radiation fields as input available at "../force/synoptic/XXX", 
	here I assumed that our current directory is experiment folder, i.e. expt_X.XX. These extra fields are generated
	if you fire the following command on terminal:
        ../bin/atmo_synoptic.sh erai+all 2007-01-02T00:00:00 2007-12-31T00:00:00
