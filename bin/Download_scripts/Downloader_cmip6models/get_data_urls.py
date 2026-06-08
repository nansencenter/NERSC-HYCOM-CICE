from cmip_url_builders import MPIBuilder, NorEsmBuilder, CNRMBuilder
from seaclim_urls import CMIP6Descriptors, SeaclimURLs
from seaclim_expts import SeaclimExpts

import datetime

# Model names
MPI = "MPI"
NorESM = "NorESM"
CNRM = "CNRM"

# Experiment names
historical = SeaclimExpts.historical
ssp126 = SeaclimExpts.ssp126
ssp245 = SeaclimExpts.ssp245
ssp370 = SeaclimExpts.ssp370

# Table names
threehr = "3hr"
sixhr = "6hrPlevPt"
e3hr = "E3hr"
omon = "Omon"
oyr = "Oyr"

"""
Generates the data URLS for retrieving SEACLIM data.
"""

"""
Defines tables as a function of variable (only).
"""
def table_MPI(var):
    sixhr_vars = {"tas", "huss", "uas", "vas", "psl"}
    threehr_vars = {"rsds", "rlds", "pr", "prsn"}
    omon_vars = {"uo", "vo", "thetao", "so", "ssh", "zos", "no3", "po4", "si", "o2", "talk", "dissic"}

    if var in sixhr_vars:
        return sixhr
    elif var in threehr_vars:
        return threehr
    elif var in omon_vars:
        return omon
    else:
        return "???"

def table_NorESM(var):
    monthly_vars = {"uo", "vo", "thetao", "so", "zos", "o2", "no3", "po4", "si", "talk", "dissic"}
    annual_vars = {}
    if var in monthly_vars:
        return omon
    elif var in annual_vars:
        return oyr
    else:
        return "???"

def table_CNRM(var):
    sixhr_vars = {"psl"}
    threehr_vars = {"tas", "huss", "rsds", "rlds", "pr", "prsn"}
    e3hr_vars = {"uas", "vas"}
    omon_vars = {"uo", "vo", "thetao", "so", "ssh", "zos", "no3", "po4", "si", "o2", "talk", "dissic"}

    if var in sixhr_vars:
        return sixhr
    elif var in threehr_vars:
        return threehr
    elif var in e3hr_vars:
        return e3hr
    elif var in omon_vars:
        return omon
    else:
        return "???"

def times_MPI(experiment, table):
    canonical_years = {
        historical: range(1960, 2015, 5),
        ssp126: range(2015, 2100, 5), 
        ssp245: range(2015, 2100, 5), 
        ssp370: range(2015, 2100, 5),
    }
    times = []
    # Offset from the base for start dates by table
    start_delta = {
        threehr: datetime.timedelta(hours = 1, minutes = 30),
        sixhr: datetime.timedelta(hours = 6),
        omon: datetime.timedelta(),
    }
    end_delta = {
        threehr: datetime.timedelta(hours = -1, minutes = -30),
        sixhr: datetime.timedelta(),
        omon: datetime.timedelta(hours = -24),
    }
    fmt = "%Y%m%d%H%M" if (table != omon) else "%Y%m"
    for y in canonical_years[experiment]:
        # The base date is the turn of the year
        base_st = datetime.datetime(y, 1, 1, 0, 0)
        base_en = datetime.datetime(y+5, 1, 1, 0, 0)
        start = base_st + start_delta[table]
        end = base_en + end_delta[table]
        times.append((start.strftime(fmt), end.strftime(fmt)))

    return times

def times_NorESM(experiment, table):
    length = 10
    initial_yr = 1960
    break_yr = 2015
    final_yr = 2101

    scenario_years = range(break_yr // length * length + 1, final_yr, length)
    canonical_years = {
        historical: range(initial_yr // length * length, break_yr, length),
        ssp126: scenario_years,
        ssp245: scenario_years,
        ssp370: scenario_years,
    }
    times = []
    # Offset from the base for start dates by table
    end_delta = datetime.timedelta(hours = -24)

    min_yr = initial_yr if experiment == historical else break_yr
    max_yr = break_yr if experiment == historical else final_yr

    fmt = "%Y%m"
    for y in canonical_years[experiment]:
        # The base date is the turn of the year
        base_st = datetime.datetime(max(y, min_yr), 1, 1, 0, 0)
        base_en = datetime.datetime(min(y+length, max_yr), 1, 1, 0, 0)
        start = base_st
        end = base_en + end_delta
        times.append((start.strftime(fmt), end.strftime(fmt)))

    return times

def times_CNRM(descriptors):
    experiment = descriptors.experiment_id
    table = descriptors.table_id
    var = descriptors.variable_id
    full_dur = 2**20
    file_durations = {
        full_dur: {"zos"},
        50: {'psl', "so", "no3", "po4", "si", "o2", "talk", "dissic"},
        25: {'tas', 'rsds', "uo", "vo", "thetao"},
        20: {'huss', 'rlds'},
        10: {'pr', 'prsn', 'uas', 'vas'},
        }
    dur = None
    for file_duration in file_durations:
        if var in file_durations[file_duration]: dur = file_duration
    if dur == None:
        raise KeyError("File duration not found for variable" + str(var))

    # For historical, the start year is always 1850 and the end year is 2015.
    # For scenarios, the start year is 2015 and the end year is 2100.
    # Find the first starting year or last ending year

    origin_yr = 1850
    break_yr = 2015
    final_yr = 2100
    start_yr = ((1960 - origin_yr) // dur) * dur + origin_yr if experiment == historical else break_yr
    end_yr = (break_yr // dur) * dur if experiment == historical else (final_yr // dur) * dur

    # Treat the special case for single file variables (dur = 1000)
    if dur == full_dur:
        (start_years, end_years) = ([origin_yr], [break_yr]) if experiment == historical else ([break_yr], [final_yr+1])
    else:
        start_years = list(range(start_yr, end_yr, dur))
        end_years = list(range(start_yr+dur, end_yr+dur, dur))

        if experiment == historical:
            start_years.append(end_years[-1] if len(end_years) else origin_yr)
            end_years.append(break_yr)
        else:
            if end_years[-1] > (final_yr + 1): end_years[-1] = final_yr + 1

    times = []
    # Offset from the base for start dates by table
    start_delta = {
        threehr: datetime.timedelta(hours = 3),
        sixhr: datetime.timedelta(hours = 6),
        e3hr: datetime.timedelta(hours = 1, minutes = 30),
        omon: datetime.timedelta(),
    }
    end_delta = {
        threehr: datetime.timedelta(),
        sixhr: datetime.timedelta(),
        e3hr: datetime.timedelta(hours = -1, minutes = -30),
        omon: datetime.timedelta(hours = -24),
    }
    time_shift_fields = {"rsds", "rlds", "pr", "prsn"}
    time_key = table if not var in time_shift_fields else e3hr
    fmt = "%Y%m%d%H%M" if (table != omon) else "%Y%m"
    for i, st in enumerate(start_years):
        # The base date is the turn of the year
        base_st = datetime.datetime(st, 1, 1, 0, 0)
        base_en = datetime.datetime(end_years[i], 1, 1, 0, 0)
        start = base_st + start_delta[time_key]
        end = base_en + end_delta[time_key]
        times.append((start.strftime(fmt), end.strftime(fmt)))

    return times

def grid_NorESM(var):
    gr_vars = {"no3", "po4", "si", "o2", "talk", "dissic"}
    return "gr" if var in gr_vars else "gr"  #else "gn"

"""
Returns several copies of the given descriptor with the times set to the set of files matching this descriptor for the given model.
"""
def set_times(model, descriptors):
    if model == MPI:
        times = times_MPI(descriptors.experiment_id, descriptors.table_id)
    elif model == NorESM or model == NorESM+"uk":
        times = times_NorESM(descriptors.experiment_id, descriptors.table_id)
    elif model == CNRM:
        # CNRM file times depend on experiment, table and variable, so just
        # pass the descriptor object
        times = times_CNRM(descriptors)
    else:
        raise NameError("get_data_urls::start_times: Unkown model " + str(model))

    timed_descriptors = []
    for t in times:
        timed_descriptors.append(descriptors.copy().set_dates(t))

    return timed_descriptors

"""
Extends a list and returns that list
"""
def rextend(l_list, r_list):
    l_list.extend(r_list)
    return l_list

"""
Does a Cartesian product over the set of models, experiments and variables with the most obvious values. Provides a base for more conditional changes.
"""
def base_urls():
    # Objects defining the servers for each model
    servers = {
        MPI: MPIBuilder("de"),
        NorESM: NorEsmBuilder("de"),
        NorESM+"uk": NorEsmBuilder("uk"),
        CNRM: CNRMBuilder("fr"),
        }
    table_fn = {
        MPI: table_MPI,
        NorESM: table_NorESM,
        NorESM+"uk": table_NorESM,
        CNRM: table_CNRM,
        }
    grid_fn = {
        NorESM: grid_NorESM,
        NorESM+"uk": grid_NorESM,
        }
    experiments = {
        MPI: SeaclimExpts.all_expts,
        NorESM: [SeaclimExpts.ssp126, SeaclimExpts.ssp245, SeaclimExpts.ssp370],
        NorESM+"uk": [SeaclimExpts.historical],
        CNRM: SeaclimExpts.all_expts
        }
    atm_var = ["tas", "huss", "uas", "vas", "rsds", "rlds", "pr", "prsn", "psl"]
    ocn_var = [
        "uo", "vo", "thetao", "so", "zos",
        "no3", "po4", "si", "o2", "talk", "dissic",
        ]
    all_var = atm_var + ocn_var
    variables = {
        MPI: all_var,
        NorESM: ocn_var,
        NorESM+"uk": ocn_var,
        CNRM: all_var,
        }

    model_name = {
        MPI: MPI,
        NorESM: NorESM,
        NorESM+"uk": NorESM,
        CNRM: CNRM,
        }
    urls = {}
    # Loop over…
    # Model
    for model in servers:
        surl = SeaclimURLs(servers[model])

        if not model_name[model] in urls:
            urls.update({model_name[model]: {}})
        m_map = urls[model_name[model]]
        # Experiment
        for expt in experiments[model]:
            m_map.update({expt: {}})
            e_map = m_map[expt]
            # Variable
            for var in variables[model]:
                e_map.update({var: []})
                v_list = e_map[var]
                desc = CMIP6Descriptors({
                    CMIP6Descriptors.eid: expt,
                    CMIP6Descriptors.vid: var,
                    CMIP6Descriptors.tid: table_fn[model](var),
                    })
                if model in grid_fn:
                    desc.grid_id = grid_fn[model](var)
                # Files/times:
                for t_desc in set_times(model, desc):
                    v_list.append(surl.url(t_desc))

    return urls

if __name__ == "__main__":
    urls = base_urls()
    f = open("urls.txt", "wt")
    # Unpack the URLS and write them serially
    for model in urls:
        for expt in urls[model]:
            for var in urls[model][expt]:
                for url in urls[model][expt][var]:
                    f.write(url+"\n")

    f.close()
