module mod_hycom_fabm
#ifdef FABM
   use fabm, only: type_model
   public
   type (type_model), save :: fabm_model
#endif
end module