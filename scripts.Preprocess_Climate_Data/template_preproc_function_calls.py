import funcs_read_obs as fro
import funcs_preproc_clim_data as fpcd
import funcs_read_ds_cmip5 as frdc5


def call_preproc_obs():
    yr0 = 1980
    obsdsetlist = ['CONUS404', 'PRISM', 'Livneh', 'NLDAS', 'GMET', 'nClimGrid', 'gridMET']
    obsind = 0 # CHANGE FOR PARALLELIZATION
    fro.regrid_obs_icargrid(yr0, obsdsetlist[obsind])


def call_concat_obs():
    # obs = 'ERA5'
    obsdsetlist = ['CONUS404', 'PRISM', 'Livneh', 'NLDAS', 'GMET', 'nClimGrid', 'gridMET']
    obsind = 0 # CHANGE FOR PARALLELIZATION
    yr0 = 1980
    fro.combine_daily_obs(obsdsetlist[obsind], year=yr0)


def call_preproc_cmip5():
    models = ['ACCESS1-3','CanESM2','CCSM4','MIROC5','MRI-CGCM3','NorESM1-M']
    imod = 0  # CHANGE FOR PARALLELIZATION
    model = models[imod]
    fpcd.preproc_cmip_forWT(model)


def call_preproc_CESM2_LE():
    ensind = 0  # CHANGE FOR PARALLELIZATION
    fpcd.prep_cesm2le_forWT(ensind=ensind, expind=0)


def call_combine_cesm2le_forWT():
    timslic = slice('1981-01-01', '2016-12-31')
    enslist = frdc5.CESM2_LE_ensemble_list()
    ensind = 0  # CHANGE FOR PARALLELIZATION
    ens = enslist[ensind]
    fpcd.combine_cesm2le_forWT(ens, timslic=timslic, scenario='hist')


def call_era5_download_daily_tp():
    year = 1980  # CHANGE FOR PARALLELIZATION
    fpcd.download_era5_daily_tp(year)


def call_preproc_era5():
    year = 1980  # CHANGE FOR PARALLELIZATION
    month = 1  # CHANGE FOR PARALLELIZATION
    fpcd.make_era5_daily(year, month)


if __name__ == "__main__":
    print("Starting function calls...")
    # call_preproc_obs()
    # call_preproc_cmip5()
    # call_preproc_CESM2_LE()
    # call_combine_cesm2le_forWT()
    # call_era5_download_daily_tp()
    # call_preproc_era5()
    # call_concat_obs()
