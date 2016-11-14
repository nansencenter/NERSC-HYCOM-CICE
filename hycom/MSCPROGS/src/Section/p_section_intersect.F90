program section_intersect
   use mod_xc
   use mod_za
   use mod_sections
   use mod_grid
   implicit none

   ! Initialize IO for .ab files
   CALL XCSPMD  ! -- Requires "regional.grid.b" to be present
   CALL ZAIOST

   ! Get grid 
   call get_grid()

   ! Read section specification from sections.in
   ! period flag set in get_grid
   call read_sections_in()  

   ! This routine sets up the p-cells active along a section. NB - great circles
   call section_nodepoints(plon,plat,idm,jdm,periodic)

   ! Join subsections
   call sections_join()

   ! Save section for later use
   call save_section_nodes()

end program section_intersect
