import funcs_DS_eval as fdse
import funcs_read_ds_cmip5 as frdc5
import xarray as xr

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
        tmpstrt = ['2064-01-01', '2044-01-01', '2024-01-01']
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


def compute_metrics_from_maps():
    regions = ['Desert Southwest', 'Great Lakes', 'Gulf Coast', 'Mid Atlantic', 
            'Mountain West', 'North Atlantic', 'Northern Plains', 
            'Pacific Northwest', 'Pacific Southwest']
    # nmod = 6
    # nmeth = 7
    imod = 0 # CHANGE FOR PARALLELIZATION
    imeth = 0 # CHANGE FOR PARALLELIZATION
    wtpr_metric_maps = {}
    diri = '/glade/derecho/scratch/nlybarger/WT_centroids_era5/wtpr_metric_maps'
    for region in regions:
        regname = region.replace(' ', '')
        wtpr_metric_maps[region] = xr.open_dataset(f'{diri}/{regname}_wtpr_DS_metric_maps_1981-2016.nc')
        wtpr_metric_maps[region] = wtpr_metric_maps[region].rename({'method': 'methods'})
    # for imod in range(nmod):
        # for imeth in range(nmeth):
    for region in regions:
        fdse.Compute_Regional_DS_Metrics(region, wtpr_metric_maps[region], imod=imod, imeth=imeth,  OVERWRITE=True)


def compute_metrics_from_maps_GARDLENS():
    regions = ['Desert Southwest', 'Great Lakes', 'Gulf Coast', 'Mid Atlantic', 
            'Mountain West', 'North Atlantic', 'Northern Plains', 
            'Pacific Northwest', 'Pacific Southwest']
    regind = 0 # CHANGE FOR PARALLELIZATION
    region = regions[regind]
    regname = region.replace(' ', '')
    wtpr_metric_maps = {}
    diri = '/glade/derecho/scratch/nlybarger/WT_centroids_era5/wtpr_metric_maps/GARDLENS/'
    wtpr_metric_maps = xr.open_dataset(f'{diri}{regname}_wtpr_metric_maps_GARDLENS_CESM2-LE.nc')
    fdse.Compute_Regional_DS_Metrics(region, wtpr_metric_maps, GARDLENS=True)



if __name__ == "__main__":
    print("Starting function calls...")
    # compute_conus_metric_maps_cmip5_ds()
    # compute_conus_metric_maps_obs()
    # compute_conus_metric_maps_gardlens()
    # compute_conus_metric_maps_cesm2le_raw()
    # compute_metrics_from_maps()
    # compute_metrics_from_maps_GARDLENS()
