      program hycom
c
c --- Non-ESMF driver for stand-alone HYCOM ocean model
c
      use mod_hycom, only : end_of_run,
     &                      HYCOM_Init,
     &                      HYCOM_Run,
     &                      HYCOM_Final
c
      implicit none
c
c --- Initialize HYCOM.
      call HYCOM_Init

c --- Run HYCOM.
      do !until end of run
        call HYCOM_Run
        if     (end_of_run) then
          exit
        endif
      enddo
c
c --- Finalize HYCOM.
      call HYCOM_Final
c
      stop
      end program hycom
