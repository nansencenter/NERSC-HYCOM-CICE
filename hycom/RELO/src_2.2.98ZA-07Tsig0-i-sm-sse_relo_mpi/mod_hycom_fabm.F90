module mod_hycom_fabm
#ifdef FABM
   use fabm, only: type_model
   use fabm_config, only: fabm_create_model_from_yaml_file
   public
   type (type_model), save :: fabm_model
#endif
end module