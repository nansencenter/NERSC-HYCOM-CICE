!=======================================================================

! Spatial parameter conditions
!
! authors: Jiping Xie, NERSC
!
! 2025: uniformed values in global netcdf files
      module ice_domain_para

      use ice_kinds_mod
      use ice_fileunits, only: nu_diag
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size, only: nx_global, ny_global, max_blocks
      use ice_constants, only: c0, c1

      private
      public :: init_cice_para,para_file,cice_para, &
                Pcice1,Pcice2,Pcice3,Pcice4,Pcice5, &
                Pcice6,Pcice7,Pcice8,Pcice9,Pcice10, &
                Pcice11,Pcice12,Pcice13,Pcice14

      public :: Gcice1,Gcice2,Gcice3,Gcice4,Gcice5,  &
                Gcice6,Gcice7,Gcice8,Gcice9,Gcice10, &
                Gcice11,Gcice12,Gcice13,Gcice14

      !! try to keep the consistence with setting in cice_in
      real (kind=dbl_kind),parameter,public ::    &
         P_ice1=330., P_ice2=917., P_ice3=0.95, P_ice4=300., &
         P_ice5=0.0005, P_ice6=0.00536, P_ice7=27500.,       &
         P_ice8=750., P_ice9=0.05, P_ice10=1.2,             &
         P_ice11=0.05, P_ice12=3.0, P_ice13=0.03, P_ice14=4.
  
      save
      
      character(len=char_len_long) :: para_file
      logical (log_kind)           :: cice_para

      real (kind=dbl_kind),  &
         dimension (nx_block,ny_block,max_blocks):: &
         Pcice1  , & ! value of rhos
         Pcice2  , & ! value of rhoi
         Pcice3  , & ! value of emissivity
         Pcice4  , & ! value of floediam
         Pcice5  , & ! value of iceruf
         Pcice6  , & ! value of dragio 
         Pcice7  , & ! value of Pstar
         Pcice8  , & ! value of rsnw_mlt
         Pcice9  , & ! value of hi_ssl
         Pcice10 , & ! value of R_snw
         Pcice11 , & ! value of astar
         Pcice12 , & ! value of mu_rdg 
         Pcice13 , & ! value of hs1 
         Pcice14     ! ice_ref_salt 


      real (kind=dbl_kind), dimension (nx_global, ny_global):: &
         Gcice1  , & ! value of rhos
         Gcice2  , & ! value of rhoi
         Gcice3  , & ! value of emissivity
         Gcice4  , & ! value of floediam
         Gcice5  , & ! value of iceruf
         Gcice6  , & ! value of dragio 
         Gcice7  , & ! value of Pstar
         Gcice8  , & ! value of rsnw_mlt
         Gcice9  , & ! value of hi_ssl
         Gcice10 , & ! value of R_snw
         Gcice11 , & ! value of astar
         Gcice12 , & ! value of mu_rdg 
         Gcice13 , & ! value of hs1 
         Gcice14     ! ice_ref_salt 

      contains

      subroutine init_cice_para
      use ice_communicate,    only: my_task, master_task
      use ice_constants,      only: c0, c1, puny,field_loc_center,   &
                                    field_type_scalar
      use ice_read_write,     only: ice_read_nc, ice_read_global_nc, &
                                    ice_open_nc, ice_close_nc,       &
                                    ice_var_exists
      use ice_broadcast,      only: broadcast_array
      use ice_blocks,         only: nx_block, ny_block
      use ice_domain_size,    only: max_blocks
      use ice_domain,         only: nblocks, distrb_info 
      use ice_gather_scatter, only: scatter_global

      implicit none
        
      integer (kind=int_kind) :: fid, it
      character (char_len)    :: fldname     ! field name in netCDF file

      real (kind=dbl_kind),   dimension(:,:),allocatable :: work1

      Gcice1=P_ice1;   Gcice2=P_ice2;   Gcice3=P_ice3; 
      Gcice4=P_ice4;   Gcice5=P_ice5;   Gcice6=P_ice6; 
      Gcice7=P_ice7;   Gcice8=P_ice8;   Gcice9=P_ice9; 
      Gcice10=P_ice10; Gcice11=P_ice11; Gcice12=P_ice12; 
      Gcice13=P_ice13; Gcice14=P_ice14;

      Pcice1=P_ice1;   Pcice2=P_ice2;   Pcice3=P_ice3; 
      Pcice4=P_ice4;   Pcice5=P_ice5;   Pcice6=P_ice6; 
      Pcice7=P_ice7;   Pcice8=P_ice8;   Pcice9=P_ice9; 
      Pcice10=P_ice10; Pcice11=P_ice11; Pcice12=P_ice12; 
      Pcice13=P_ice13; Pcice14=P_ice14;

      if (my_task == master_task) then
         call ice_open_nc(para_file, fid)
         do it =1,14
            allocate(work1(nx_global, ny_global))
            select case(it)
            case(1)
               fldname='rhos'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice1(:,:)=work1
               endif        
            case(2)
               fldname='rhoi'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice2(:,:)=work1
               endif
            case(3)
               fldname='emissi'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice3(:,:)=work1
               endif
            case(4)
               fldname='floediam'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice4(:,:)=work1
               endif
            case(5)
               fldname='iceruf'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice5(:,:)=work1
               endif
            case(6)
               fldname='dragio'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice6(:,:)=work1
               endif
            case(7)
               fldname='Pstar'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice7(:,:)=work1
               endif
            case(8)
               fldname='rsnw_mlt'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice8(:,:)=work1
               endif
            case(9)
               fldname='hi_ssl'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice9(:,:)=work1
               endif
            case(10)
               fldname='R_snw'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice10(:,:)=work1
               endif
            case(11)
               fldname='astar'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice11(:,:)=work1
               endif
            case(12)
               fldname='mu_rdg'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice12(:,:)=work1
               endif
            case(13)
               fldname='hs1'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice13(:,:)=work1
               endif
            case(14)
               fldname='ice_ref_salt'
               if (ice_var_exists(fid,fldname)) then
                  call ice_read_global_nc(fid,1,fldname,work1,.true.)
                  Gcice14(:,:)=work1
               endif
            end select
            deallocate(work1)
         enddo
         call ice_close_nc(fid)
      endif

      end subroutine


      end module ice_domain_para
