!undef PARTICLE   # Include code for particle tracers
!undef ICE        # Save and read ice restart files
!undef HIBLER     # Include ice dynamics from Hardner
!undef ENSEMBLE   # Include code for ensemble stuff
!undef RIVER      # Include River fluxes (edit the ./Data/river.dat file)
!undef MIXDIAG    # computes and accumulates vertical transports
!undef NEST_OUTER # Save nesting conditions to be used in local models
!undef NEST_INNER # read nesting conditions saved by a global model
!undef OPEN_MP    # Parallelization using OPEN_MP directives.

!undef CLIMPRT    # Dont add just print current field and stop
!undef DEBUG_FRC
!undef DEBUG_RELAX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define OPEN_MP
#undef OWI_WINDS
#undef NEST_OUTER
#undef NEST_INNER
#undef TIDES
#undef HIBLER
#undef ICE
