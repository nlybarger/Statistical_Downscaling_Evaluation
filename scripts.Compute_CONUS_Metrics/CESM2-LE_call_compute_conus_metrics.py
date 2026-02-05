import funcs_DS_eval as fdse
import funcs_read_ds_cmip5 as frdc5

def compute_conus_metric_maps_cmip5_ds():
    scenarios = ['hist', 'rcp45', 'rcp85']
    iscen = 0 # CHANGE FOR PARALLELIZATION
    scen = scenarios[iscen]
    imod = 0 # CHANGE FOR PARALLELIZATION
    imeth = 0 # CHANGE FOR PARALLELIZATION
    mod,meth,_,_,_ = frdc5.get_ds_metavars(FULL_LIST=False, imod=imod, imeth=imeth)
    if scen == 'hist':
        TimeStart = '1981-01-01'
        TimeEnd = '2016-12-31'
    else:
        itim = 0
        tmpstrt = ['2076-01-01', '2056-01-01', '2036-01-01']
        tmpend = ['2099-12-31', '2079-12-31', '2059-12-31']
        TimeStart = tmpstrt[itim]
        TimeEnd = tmpend[itim]
    timslic = slice(TimeStart, TimeEnd)
    fdse.Compute_CONUS_Metric_Maps(timslic, mod=mod, meth=meth, scenario=scen, OVERWRITE=True)


def compute_conus_metric_maps_obs():
    obsdsetlist = ['CONUS404', 'PRISM', 'Livneh', 'NLDAS', 'GMET', 'nClimGrid', 'gridMET']
    obsind = 0 # CHANGE FOR PARALLELIZATION
    obsdset = obsdsetlist[obsind]
    TimeStart = '1981-01-01'
    TimeEnd = '2016-12-31'
    timslic = slice(TimeStart, TimeEnd)
    fdse.Compute_CONUS_Metric_Maps(timslic, obsdset=obsdset, OVERWRITE=True)


def compute_conus_metric_maps_gardlens():
    scen = 'hist'
    TimeStart = '1981-01-01'
    TimeEnd = '2016-12-31'
    timslic = slice(TimeStart, TimeEnd)

    ensind = 0 # CHANGE FOR PARALLELIZATION

    enslist = frdc5.CESM2_LE_ensemble_list(GARDLENS=True)
    ens = enslist[ensind]
    fdse.Compute_CONUS_Metric_Maps(timslic, ens=ens, GARDLENS=True, OVERWRITE=True)


def compute_conus_metric_maps_cesm2le_raw():
    scen = 'hist'
    TimeStart = '1981-01-01'
    TimeEnd = '2016-12-31'
    timslic = slice(TimeStart, TimeEnd)
    # ensind = 0 # CHANGE FOR PARALLELIZATION
    enslist = frdc5.CESM2_LE_ensemble_list(GARDLENS=False)
    # ens = enslist[ensind]
    for ens in enslist:
        fdse.Compute_CONUS_Metric_Maps(timslic, ens=ens, GARDLENS=False, OVERWRITE=True)


if __name__ == "__main__":
    print("Starting function calls...")
    # compute_conus_metric_maps_cmip5_ds()
    # compute_conus_metric_maps_obs()
    # compute_conus_metric_maps_gardlens()
    compute_conus_metric_maps_cesm2le_raw()
