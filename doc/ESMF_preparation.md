#Introduction
This short note explains some of the steps needed to set up a model for ESMF.
Examples will be taken from the HYCOM and CICE code used in NERSC-HYCOM-CICE.
The documentation is not exhaustive, for that there is the ESMF User and
Reference manuals. However, we hope that the document will be of practical use
for doing the ESMF-ization of an existing code. 

The Necessary steps are as follows, and will be described in greater detail
below:

* Structure your code to fit with the ESMF framework. This means separating the
code into subroutines for initializing, running and finalizing the model.

* Initialize the ESMF Framework

* Register the run, initialize and finalize methods with the ESMF framework.

* Running the ESMF init routines

* * Get numper of cpus, and mpi communicator object from from ESMF system

* * Set up ESMF data structures. This includes import/export fields, the
description of how the model grid is split among CPUs for MPI runs. The
coordinates and masks of the model grid. Also additional tools, such as clocks,
are described.

* Running the ESMF run routines. Setting up exchanges of fields.

* Running the ESMF coupler object


# Structure code to fit in the ESMF framework

In order to  set up your model code to fit with the ESMF framework, you need to split
it between a Init phase, a run phase and a finalize stage. Each of these stages
can have multiple phases, making it possible to have more fine-grained control
of the code execution. The stages do the following

# # Initialize  stage
 
The initialize stage will take care of initialization of the model, such as
setting up model grids, set up a initial state, time variables,  etc etc. All the stuff that
needs to be done before actually running the model.  In addition, the initialize
stage usually takes care of setting up many of the structures needed by ESMF
itself. This involves setting up the description of the grid in a way that
ESMF can understand (ESMF_Grid), setting up and Export State, and an Import
State. The two latter will be needed to exhange fields between different models.

In HYCOM, this stage is done by the subroutine HYCOM_Init inside mod_hycom.F , while in CICE,
it is done by the subroutine ice_init_esmf in esmf/drivers/cice_comp_esmf.F90

# # Run  stage

This step takes care of the actual time integration of the model. In addition it
will  take care of importing variables from other models from the "Import
State", and putting variables to be exported to other models into an "Export
State". Note that this routine should only run the model forward in time for a
limited number of time steps

In HYCOM, this is done by the subroutine HYCOM_Run in mod_hycom.F. While in  in CICE,
it is done by the subroutine ice_run_esmf in esmf/drivers/cice_comp_esmf.F90

# # Finalize stage

This routine is responsible for making sure the model exits cleanly, as well as
making sure that ESMF also is closed down in a clean way.

In HYCOM, this is done by the subroutine HYCOM_Final in mod_hycom.F. In  in CICE,
it is done by the subroutine ice_final_esmf in esmf/drivers/cice_comp_esmf.F90


# Initialize ESMF Framework

This is one of the first steps needed to do in the code. If it is called before
doing MPI_Init() in the code, then ESMF will set up MPI, and create a MPI
communicator object. In the main program  (hycom_cice.), ESMF_Initialize is
called like this 


         call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN,
        &                      logkindflag=ESMF_LOGKIND_MULTI,
        &                                  vm=worldVM,
        &                                  rc=rc)

It is possible to inquire the virtual machine object (here, worldVM), to get
information at runtime. The following gets the number of MPI tasks (petcount), and
the rank of this MPI task (localPet)

         call ESMF_VMGet(worldVM, petCount=petCount, localPET=localPet,
        &                rc=rc)

Once the ESMF f

# Register Start, Run and Finalize stages in ESMF

This need 

All calls of a models Initialize, Run and Finalize stages will be done through
the ESMF framework, so ESMF needs to know what routines to call. This is done
with the setservices call. The following shows the setservices call in HYCOM 


         call ESMF_GridCompSetEntryPoint(
        &     gridComp,
        &     ESMF_METHOD_INITIALIZE,
        &     HYCOM_Init,
        &     phase=1,
        &     rc=rc)
   cKAL  NEW
         call ESMF_GridCompSetEntryPoint(
        &     gridComp,
        &     ESMF_METHOD_INITIALIZE,
        &     HYCOM_Init_2,
        &     phase=2,
        &     rc=rc)
   cKAL  NEW
         call ESMF_GridCompSetEntryPoint(
        &     gridComp,
        &     ESMF_METHOD_RUN,
        &     HYCOM_Run,
        &     phase=1,
        &     rc=rc)
         call ESMF_GridCompSetEntryPoint(
        &     gridComp,
        &     ESMF_METHOD_FINALIZE,
        &     HYCOM_Final,
        &     phase=1,
        &     rc=rc)

and CICE:


    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
      ice_init_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
      ice_run_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
      ice_final_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)



# ESMF ts
