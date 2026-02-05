import funcs_weather_typing as fwt
import funcs_DS_eval as fdse


"""
This is a template script to call the various functions in funcs_weather_typing.py. Each function call is wrapped in 
its own function for ease of parallelization. Modify the regind, imod, obsind variables as needed for parallelization.
This and the associated PBS script are duplicated and modified in the Generate*.bash scripts to run on Casper.
"""
def call_Optimal_Cluster_Test():
    regind = 7 # CHANGE FOR PARALLELIZATION
    TimeStart = '1981-01-01'
    TimeEnd = '2016-12-31'
    fwt.Optimal_Cluster_Test(regind, numWTs=6, TimeStart=TimeStart, TimeEnd=TimeEnd)


def call_era5_centroids():
    regions, bukods = fdse.make_buco_region_map()
    TimeStart = '1981-01-01'
    TimeEnd = '2016-12-31'
    nWTs = 6
    _, VariableList = fwt.get_regional_WT_varlist()
    lonmin, lonmax, latmin, latmax = fdse.get_region_edges(regions, bukods)
    for reg in regions:
        latrng = slice(latmin[reg], latmax[reg])
        lonrng = slice(lonmin[reg], lonmax[reg])
        fwt.Compute_WT_centroids_ERA5(TimeStart, TimeEnd, VariableList[reg], lonrng, latrng, reg, numberWTs=nWTs, OVERWRITE=True)


def call_wtpr_cmip5():
    regind = 7 # CHANGE FOR PARALLELIZATION
    imod = 0 # CHANGE FOR PARALLELIZATION
    fwt.compute_wtpr_cmip5(regind, imod, TimeStart='1981-01-01', TimeEnd='2016-12-31', numWTs=6,OVERWRITE=False)


def call_wtpr_obs():
    regind = 7 # CHANGE FOR PARALLELIZATION
    obsind = 0 # CHANGE FOR PARALLELIZATION
    fwt.compute_wtpr_obs(regind, obsind, TimeStart='1981-01-01', TimeEnd='2016-12-31', numWTs=6,OVERWRITE=True)

def call_wtpr_GARDLENS():
    regind = 7 # CHANGE FOR PARALLELIZATION
    fwt.compute_wtpr_GARDLENS(regind, TimeStart='1981-01-01', TimeEnd='2016-12-31', numWTs=6,OVERWRITE=True)


def call_compute_rmse_map():
    regind = 7 # CHANGE FOR PARALLELIZATION
    TimeStart = '1981-01-01'
    TimeEnd = '2016-12-31'
    fwt.compute_rmse_map(regind, TimeStart=TimeStart, TimeEnd=TimeEnd, numWTs=6)


def call_compute_rmse_map_GARDLENS():
    regind = 7 # CHANGE FOR PARALLELIZATION
    fwt.compute_rmse_map_GARDLENS(regind, numWTs=6)


if __name__ == "__main__":
    print("Starting function calls...")
    # call_Optimal_Cluster_Test()
    # call_era5_centroids()
    # call_wtpr_cmip5()
    # call_wtpr_obs()
    # call_wtpr_GARDLENS()
    call_compute_rmse_map()
    # call_compute_rmse_map_GARDLENS()
