"""
Fetches forcing data from a CMIP6 data mirror.
"""

import copy

from cmip_url_builders import *
from seaclim_expts import SeaclimExpts

"""
A class to hold the various descriptors of CMIP6 data.
"""
class CMIP6Descriptors():
    historical = SeaclimExpts.historical
    ssp126 = SeaclimExpts.ssp126
    ssp245 = SeaclimExpts.ssp245
    ssp370 = SeaclimExpts.ssp370
    
    experiments = {
        historical,
        ssp126,
        ssp245,
        ssp370,
        }

    activity = {
        historical: "CMIP",
        ssp126: "ScenarioMIP",
        ssp245: "ScenarioMIP",
        ssp370: "ScenarioMIP",
        }


    eid = "experiment_id"
    tid = "table_id"
    vid = "variable_id"
    did = "dates"
    aid = "activity"

    def __init__(self, arg = {}):
        self.experiment_id = None
        self.table_id = None
        self.variable_id = None

        arg_map = {
            CMIP6Descriptors.eid: self.set_experiment_id,
            CMIP6Descriptors.tid: self.set_table_id,
            CMIP6Descriptors.vid: self.set_variable_id,
            CMIP6Descriptors.did: self.set_dates
            }
        for i in arg_map:
            if i in arg:
                arg_map[i](arg[i])

    def set_experiment_id(self, eid):
        if not eid in CMIP6Descriptors.experiments:
            raise KeyError(f"Experiment id \"{eid}\" not allowed. Permitted experiment IDs are {CMIP6Descriptors.experiments}")
        self.experiment_id = eid
        self.activity_id = CMIP6Descriptors.activity[eid]
        return self
    def set_variable_id(self, vid):
        self.variable_id = vid
        return self
    def set_table_id(self, tid):
        self.table_id = tid
        return self
    def set_dates(self, start_all, end = None):
        if end is None:
            return self.set_dates(start_all[0], start_all[1])
        self.dates = start_all + "-" + end
        return self
    def copy(self):
        return copy.deepcopy(self)
    def __str__(self):
        string = "{"
        string += f"activity_id: {self.activity_id}, experiment_id: {self.experiment_id}, variable_id: {self.variable_id} ,table_id: {self.table_id}"

        if hasattr(self, "dates"): string += f", dates: {self.dates}"

        return string + "}"

"""
A base class that fetches data from a given server.
"""
class SeaclimURLs():
    def __init__(self, path_builder):
        self.host = path_builder.domain
        self.path_builder = path_builder
        self.scheme = "https" if path_builder.is_secure() else "http"
    def path(self, descriptors):
        return self.path_builder(descriptors)
    def url(self, arg):
        if isinstance(arg, CMIP6Descriptors):
            path = self.path(arg)
        else:
            path = arg
        return self.scheme + "://" + self.host + "/" + path + ".nc"
