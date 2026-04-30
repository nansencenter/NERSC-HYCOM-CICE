from seaclim_expts import SeaclimExpts
"""
A class to build URLs from CMIP6 file descriptors
"""
class CmipUrlBuilder():
    class DomainId():
        def __init__(self, short_name, domain):
            self.short_name = short_name
            self.domain = domain
        def __str__(self):
            return f"\"{self.domain}\", alias \"{self.short_name}\""
        def __repr__(self):
            return str(self)
        def __call__(self):
            return {self.short_name: self.domain}

    @staticmethod
    def available_domains():
        return []
    @staticmethod
    def get_domain(d):
        return None

    def set_domain(self, d):
        pass

    def is_secure(self):
        return True
    
    def __call__(self, descriptors):
        ps = self._path_segments(descriptors)
        ff = self._file_fragments(descriptors)
        return "/".join(ps) + "/" + "_".join(ff)
    
"""
A class to build paths to files for MPI-ESM1 data.
"""
class MPIBuilder(CmipUrlBuilder):
    de = "esgf3.dkrz.de"
    uk = "esgf.ceda.ac.uk"
    domains = [
        CmipUrlBuilder.DomainId("de", de),
        CmipUrlBuilder.DomainId("uk", uk),
        ]

    security = {
        de: True,
        uk: True,
        }

    service = {
        de: "thredds/fileServer",
        uk: "thredds/fileServer",
        }
        
    cmip_root = {
        de: "cmip6",
        uk: "esg_cmip6/CMIP6",
        }

    default_version = "v20190710"

    # Institution varies by experiment
    institution = {
        SeaclimExpts.historical: "MPI-M",
        SeaclimExpts.ssp126: "DKRZ",
        SeaclimExpts.ssp245: "DKRZ",
        SeaclimExpts.ssp370: "DKRZ",
        }

    variant = "r1i1p1f1"
    grid_label = "gn"
    source = "MPI-ESM1-2-HR"

    @staticmethod
    def available_domains():
        return MPIBuilder.domains
    @staticmethod
    def get_domain(d):
        for did in MPIBuilder.domains:
            if d == did.short_name or d == did.domain:
                return did.domain

    def __init__(self, domain = "uk"):
        self.set_domain(domain)

    def set_domain(self, d):
        for did in MPIBuilder.domains:
            if d == did.short_name or d == did.domain:
                self.domain = did.domain
                self.security = MPIBuilder.security[self.domain]
                self.service = MPIBuilder.service[self.domain]
                self.cmip_root = MPIBuilder.cmip_root[self.domain]
        return self

    def is_secure(self):
        return self.security

    def _version(self, descriptors):
        if (descriptors.experiment_id == SeaclimExpts.historical and
            descriptors.variable_id in ("tas", "psl")):
                return "v20190815"
        else:
            return MPIBuilder.default_version

    def _path_segments(self, descriptors):
         return [
            self.service,
            self.cmip_root,
            descriptors.activity_id,
            MPIBuilder.institution[descriptors.experiment_id],
            MPIBuilder.source,
            descriptors.experiment_id,
            MPIBuilder.variant,
            descriptors.table_id,
            descriptors.variable_id,
            MPIBuilder.grid_label,
            self._version(descriptors),
            ]
    def _file_fragments(self, descriptors):
        return [
            descriptors.variable_id,
            descriptors.table_id,
            MPIBuilder.source,
            descriptors.experiment_id,
            MPIBuilder.variant,
            MPIBuilder.grid_label,
            descriptors.dates, # Start and end dates
            ]

"""
A class to build paths to files for NorESM data.
"""
class NorEsmBuilder(CmipUrlBuilder):
    de = "esgf3.dkrz.de"
    uk = "esgf.ceda.ac.uk"
    domains = [
        CmipUrlBuilder.DomainId("de", de),
        CmipUrlBuilder.DomainId("uk", uk),
        ]

    security = {
        de: True,
        uk: True,
        }

    service = {
        de: "thredds/fileServer",
        uk: "thredds/fileServer",
        }
        
    cmip_root = {
        de: "cmip6",
        uk: "esg_cmip6/CMIP6",
        }

    scenario_version = "v20191108"
    version = {
        SeaclimExpts.historical: "v20191108",
        SeaclimExpts.ssp126: scenario_version,
        SeaclimExpts.ssp245: scenario_version,
        SeaclimExpts.ssp370: scenario_version,
        }

    # Institution varies by experiment
    institution = "NCC"
    variant = "r1i1p1f1"
    grid_label = "gr" #"gn"
    source = "NorESM2-MM"

    @staticmethod
    def available_domains():
        return NorEsmBuilder.domains
    @staticmethod
    def get_domain(d):
        for did in NorEsmBuilder.domains:
            if d == did.short_name or d == did.domain:
                return did.domain

    def __init__(self, domain = "uk"):
        self.set_domain(domain)

    def set_domain(self, d):
        for did in NorEsmBuilder.domains:
            if d == did.short_name or d == did.domain:
                self.domain = did.domain
                self.security = NorEsmBuilder.security[self.domain]
                self.service = NorEsmBuilder.service[self.domain]
                self.cmip_root = NorEsmBuilder.cmip_root[self.domain]
        return self

    def is_secure(self):
        return self.security

    def _path_segments(self, descriptors):
         return [
            self.service,
            self.cmip_root,
            descriptors.activity_id,
            NorEsmBuilder.institution,
            NorEsmBuilder.source,
            descriptors.experiment_id,
            NorEsmBuilder.variant,
            descriptors.table_id,
            descriptors.variable_id,
            descriptors.grid_id if hasattr(descriptors, "grid_id") else NorEsmBuilder.grid_label,
            self.version[descriptors.experiment_id],
            ]
    def _file_fragments(self, descriptors):
        return [
            descriptors.variable_id,
            descriptors.table_id,
            NorEsmBuilder.source,
            descriptors.experiment_id,
            NorEsmBuilder.variant,
            descriptors.grid_id if hasattr(descriptors, "grid_id") else NorEsmBuilder.grid_label,
            descriptors.dates, # Start and end dates
            ]

"""
A class to build paths to files for NorESM data.
"""
class CNRMBuilder(CmipUrlBuilder):
    de = "esgf3.dkrz.de"
    uk = "esgf.ceda.ac.uk"
    fr = "esg1.umr-cnrm.fr"
    domains = [
        CmipUrlBuilder.DomainId("de", de),
        CmipUrlBuilder.DomainId("uk", uk),
        CmipUrlBuilder.DomainId("fr", fr),
        ]

    security = {
        de: False,
        uk: True,
        fr: False,
        }

    service = {
        de: "thredds/fileServer",
        uk: "thredds/fileServer",
        fr: "thredds/fileServer",
        }
        
    cmip_root = {
        de: "cmip6",
        uk: "esg_cmip6/CMIP6",
        fr: "CMIP6_CNRM",
        }
    scenario_version = "v20190328"
    version = {
        SeaclimExpts.historical: "v20181206",
        SeaclimExpts.ssp126: scenario_version,
        SeaclimExpts.ssp245: scenario_version,
        SeaclimExpts.ssp370: "v20191021",
        }

    # Institution varies by experiment
    institution = "CNRM-CERFACS"
    variant = "r1i1p1f2"
    grid_label = lambda table: "gr" if table != "Omon" else "gn"
    source = "CNRM-ESM2-1"

    @staticmethod
    def available_domains():
        return CNRMBuilder.domains
    @staticmethod
    def get_domain(d):
        for did in CNRMBuilder.domains:
            if d == did.short_name or d == did.domain:
                return did.domain

    def __init__(self, domain = "uk"):
        self.set_domain(domain)

    def set_domain(self, d):
        for did in CNRMBuilder.domains:
            if d == did.short_name or d == did.domain:
                self.domain = did.domain
                self.security = CNRMBuilder.security[self.domain]
                self.service = CNRMBuilder.service[self.domain]
                self.cmip_root = CNRMBuilder.cmip_root[self.domain]
        return self

    def is_secure(self):
        return self.security

    def _path_segments(self, descriptors):
         return [
            self.service,
            self.cmip_root,
            descriptors.activity_id,
            CNRMBuilder.institution,
            CNRMBuilder.source,
            descriptors.experiment_id,
            CNRMBuilder.variant,
            descriptors.table_id,
            descriptors.variable_id,
            CNRMBuilder.grid_label(descriptors.table_id),
            self.version[descriptors.experiment_id],
            ]
    def _file_fragments(self, descriptors):
        return [
            descriptors.variable_id,
            descriptors.table_id,
            CNRMBuilder.source,
            descriptors.experiment_id,
            CNRMBuilder.variant,
            CNRMBuilder.grid_label(descriptors.table_id),
            descriptors.dates, # Start and end dates
            ]
