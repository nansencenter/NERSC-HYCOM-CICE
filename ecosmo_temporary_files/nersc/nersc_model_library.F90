module nersc_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: nersc_model_factory

contains

   subroutine create(self,name,model)

      use fabm_nersc_ecosmo
!      use ecosmo_bg
      ! Add new BB models here

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('ecosmo');       allocate(type_nersc_ecosmo::model)
!         case ('ecosmo_bg');       allocate(type_nersc_ecosmo_bg::model)
         ! Add new BB models here
      end select

   end subroutine



end module
