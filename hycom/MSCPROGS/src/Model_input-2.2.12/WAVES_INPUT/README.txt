README.txt
Author: Timothy Williams
Date:   20130614, 11:01:49 CEST

#########################################################################
Pivot point code for nearest-neighbour interpolation
- adapted from Jon's p_wamnsea.F90 code
- independent of model, and easy to change for new data source
- grid boundary of data source allowed to be included
  (can be important for low-resolution data)
#########################################################################

#########################################################################
In 'WAVES/PIVOTS' folder:

F90 code and makefiles to get pivot points for
- WW3 (WAVEWATCH III)
- WAMNSEA (WAM North Sea)

In 'RUN_HERE' directory:
'get_pivots.sh' - runs executable on all models
'test_pivots.m' - compare interpolated variables to original ones
   *change 'mod_no' to check different model
   *change 'itest','jtest' to give diagnostics at particular point
#########################################################################

#########################################################################
In 'ICECONS/PIVOTS' folder:

F90 code and makefiles to get pivot points for
- SMOS (Ice thickness)
- AMSR-E (Ice concentration)

In 'RUN_HERE' directory:
'get_pivots.sh' - runs executable on all models
'test_pivots.m' - compare interpolated variables to original ones
   *change 'mod_no' to check different model
   *change 'itest','jtest' to give diagnostics at particular point

#########################################################################

#########################################################################
In 'ICECONS/DAILYAVE_WAVES_INPUT' folder:

Script 'dailyave_waves_input.sh' to link and rename daily average files
so ficem and hicem can be used as inputs in lieu of
SMOS/AMSR-E

#########################################################################
