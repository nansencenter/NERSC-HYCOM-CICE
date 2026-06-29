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

* Create gridded components (models), coupler components, and import/export states for these.

* Register the run, initialize and finalize methods with the ESMF framework.

* Running the ESMF init routines

  * Get numper of cpus, and mpi communicator object from from ESMF system

  * Set up ESMF data structures. This includes import/export fields, the
description of how the model grid is split among CPUs for MPI runs. The
coordinates and masks of the model grid. Also additional tools, such as clocks,
are described.

* Running the ESMF run routines. Setting up exchanges of fields.

* Running the ESMF coupler object


# Structure code to fit in the ESMF framework

In order to  set up your model code to fit with the ESMF framework, you need to split
it between a Init phase, a run phase and a finalize routine. Each of these routines
can have multiple phases, making it possible to have more fine-grained control
of the code execution. The routines do the following

## Initialize  routine
 
The initialize routine will take care of initialization of the model, such as
setting up model grids, set up a initial state, time variables,  etc etc. All the stuff that
needs to be done before actually running the model.  In addition, the initialize
routine usually takes care of setting up many of the structures needed by ESMF
itself. This involves setting up the description of the grid in a way that
ESMF can understand (ESMF_Grid), setting up and Export State, and an Import
State. The two latter will be needed to exhange fields between different models.

In HYCOM, this routine is done by the subroutine HYCOM_Init inside mod_hycom.F , while in CICE,
it is done by the subroutine ice_init_esmf in esmf/drivers/cice_comp_esmf.F90

## Run  routine

This step takes care of the actual time integration of the model. In addition it
will  take care of importing variables from other models from the "Import
State", and putting variables to be exported to other models into an "Export
State". Note that this routine should only run the model forward in time for a
limited number of time steps

In HYCOM, this is done by the subroutine HYCOM_Run in mod_hycom.F. While in  in CICE,
it is done by the subroutine ice_run_esmf in esmf/drivers/cice_comp_esmf.F90

## Finalize routine

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

Once initialized, it is possible to inquire the virtual machine object returned by ESMF_Initialize (here, worldVM), to get
information at runtime. 

# Create gridded components (models), coupler components, and import/export states for these.

## Create gridded component

The gridded components are created by the main program in hycom_cice.F. To begin with, the gridded components are empty, but they are populated with more information as the gridded components call their initialize routines (to be registered in the next section). In hycom_cice.F there are two gridcomponents, one for the ocean (HYCOM) and one for the sea-ice (CICE). The gridded component for the sea ice model is created like this in hycom_cice.F:


    c --- Create the SEAICE gridded component
          iceGridComp = ESMF_GridCompCreate(
         &  name='SEAICE Component', rc=rc)


## Create coupler components

The coupler component takes care of the communication between the gridded components describing HYCOM and CICE. It is created in hycom_cice.F. As with the gridded components, the coupler component is initially empty, but will be filled up with information when initialized at a later stage. We will not say a lot about the coupler components in this document, since it is relatively simple in our case. Since the gridded components CICE and HYCOM use the same model grid, the coupled component is mainly responsible for transferring fields from one models import state to another models export state and vice versa. 

If the model grids were different, some interpolation would be necessary, and it would be natural to do this inside the coupler component. The same applies to field processing that depends on more than one gridded component for instance.

The coupler component is created like this in hycom_cice.F

     c --- Create the OICPL coupler component
           o2iCplComp = ESMF_CplCompCreate(
          &     name="OICPL Coupler Component", rc=rc)


## Create import and export states

After the creation of the gridded components, the import and export states of the gridded components are created. The import and export states provide the gridded components with the necessary input and output fields it needs to communicate with the other gridcomponents. The communication between the gridded component is always done through the coupler object.

The following shows the creation of the CICE import and export states in hycom_cice.F:

     c --- Create empty SEAICE import/export states
           iceImpState = ESMF_StateCreate(Name="SEAICE Import",
          &      stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
     c
           iceExpState = ESMF_StateCreate(Name="SEAICE Export",
          &      stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)

Note that these states are provided as input and output arguments to the initialize, run and finalize routines of the gridded model. For the coupler component, the states are added to the coupler component so that it is available for routine inside the coupled component object.

# Register Start, Run and Finalize routines in ESMF


All calls of a models Initialize, Run and Finalize routines will be done through
the ESMF framework, using the ESMF_GridCompInitialize, 
ESMF_GridCompRun and ESMF_GridCompFinalize routines; so ESMF needs to know what routines do what.
This is done
with the setservices call which tells ESMF what routines to run when a model components ESMF_GridCompInitialize, 
ESMF_GridCompRun and ESMF_GridCompFinalize routines are called. The following shows the setservices call in HYCOM 


         call ESMF_GridCompSetEntryPoint(
        &     gridComp,
        &     ESMF_METHOD_INITIALIZE,
        &     HYCOM_Init,
        &     phase=1,
        &     rc=rc)
         call ESMF_GridCompSetEntryPoint(
        &     gridComp,
        &     ESMF_METHOD_INITIALIZE,
        &     HYCOM_Init_2,
        &     phase=2,
        &     rc=rc)
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

In short; When the hycom model components initialize routine
(signalled by argument ESMF_METHOD_INITIALIZE) is called, the
HYCOM_Init routine will be called, when the hycom model components run routine
(signaled by argument ESMF_METHOD_RUN) is called, the HYCOM_Run routine will be
called, etc. etc ..
These routines are located in mod_hycom.F.

For CICE, the corresponding setup is as follows

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
      ice_init_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
      ice_run_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
      ice_final_esmf, phase=1, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

these routines are located in drivers/esmf/cice_comp_esmf.F90.


# Running the ESMF init routines for the gridded components

These routines usually take care of setting up the ESMF gridded component attributes. Things such as the grid description (longitude, latitude and masks) are set her, as well as a description of which MPI takes care of which part of the domain. The ESMF init routines also set up the import and export states.

In addition, the initialization of the model itself is usually done inside this routine. Below is a short description of the typical ESMF-related tasks done in this routine. This comes in addition to whatever initialization is needed for the model itself.

## Get numper of cpus, and mpi communicator object from from ESMF system

One of the first steps when initializing the model is to know how many CPUs are available for it. This is done through the ESMF framework. The framework provides the gridded components with a virtual machine(vm) object that has important information for the model.

Once inside the gridded components initialize routine, the vm object can be retrieved like this (see HYCOM_Init):

    c --- Get VM from gridComp
          call ESMF_GridCompGet(gridComp, vm=vm, rc=rc)

It is possible to inquire the virtual machine object , to get
information at runtime. The following gets the number of MPI tasks (petcount), 
the rank of this MPI task (localPet), as well as the MPI communicator object for this gridded component.

         call ESMF_VMGet(vm, petCount=petCount, localPET=localPet,
        &                rc=rc, mpiCommunicator=mpiCommunicator)

## Set up the distributed grid



### Set up a DELayout

A DELayout is a mapping 
