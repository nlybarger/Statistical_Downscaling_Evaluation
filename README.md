Suite of Python scripts and Jupyter notebooks for the evaluation and ranking of statistically 
downscaled CMIP5 datasets based on an extensive metric suite, as well as plotting figures 
corresponding to the description in Lybarger et al., (2026).

Each scripts folder performs one aspect of the analysis. Within each of these, there are 
generic template files that call functions stored within scripts.Functions and default 
PBS scripts for job submission on the NCAR HPC. There are also bash Generator scripts that
create multiple copies of these templates for crude parallelization across HPC nodes.


scripts.Preprocess_Climate_Data:
  These scripts have hard-coded paths to observational datasets, downscaled CMIP5 datasets,
  raw CMIP5 data, CESM2-LE data, and GARD-LENS data. These would be very difficult to
  customize for an HPC environment other than NCAR Glade. Essentially, these scripts take
  the raw data you'd get from online repositories and processed them to look the same.
  For downscaled data and high resolution observational data, this includes removing leap days, 
  regridding to a common 1/8 degree grid, and trimming to the common temporal era. For CMIP5
  data, CESM2-LE data, and ERA5, this includes removing leap days, regridding to a common 1
  degree grid (for WT applications), and trimming to the common temporal era.

scripts.Compute_Weather_Types
  These scripts use that preprocessed data to perform the weather typing analysis described
  in the paper. The optimal cluster test is performed to determine the meteorological
  variables that maximize cross-WT variance in WT-average precipitation using ERA5 (1 degree) data.
  Then, centroids are computed for each region using the ERA5 data and the optimal variable
  set of that region using the SciPy k-means algorithm. Raw CMIP5 model data days are then
  categorized into the ERA5 centroids using their corresponding meteorological variables.
  Downscaled CMIP5 data are then sorted into those weather types using the timestamps of the
  raw model data, and the same is done for high resolution observations using ERA5 data.
  The same process is carried out for CESM2-LE and the corresponding GARD-LENS data.
  Then, metrics comparing the day-to-day correspondence between the days within each WT
  and the climatological WT precipitaton are computed.

scripts.Compute_CONUS_metrics
  Uses the preprocessed downscaled CMIP5, GARD-LENS, and high resolution observations to compute
  the Temperature and Precipitation metrics. It also then compares the metric maps computed for
  each downscaled dataset to those computed for observations to get error metrics, that are then
  further analyzed, normalized, and plotted in scripts.Plotting.

scripts.Plotting
  Performs the final analyses and Figure plotting for metric normalization, uncertainty
  partitioning, method-model ranking, and climate change signal comparison.
