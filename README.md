# NASA International Space Apps Challenge - InSAR Change Detectives!

Team LTI Software and Engineering:

 - Alexandre Beaulieu
 - Mathieu Bergeron
 - Gabriel Gibeau Sanchez
 - Olivier Lavergne
 - Maxime Rondeau

## Notebook location
Notebook location: 
./notebooks/Report.ipynb
## Overview
Here are the processing steps.

### Data source
#### Satellite products
We use L1 data from the SENTINEL-1A Interferometric Wide Swath Level 1 S Product. Data are downloaded via the [python scripts](data/download_2015NepalEarthquake.py) provided by the Alaska Satellite Facility's [bulk download service][asf].

#### Weather data
We use weather data from the [Copernicus Climate Change Service (C3S) Climate Data Store][atmopspheric-data]. There is a utility scrip to download the atmospheric data corresponding to three geological events of interest; namely the Ridgecrest CA USA earthquakes of July 4th and 5th 2019 the Kilauea HI USA volcanic eruption of May 4th 2018 and the Kumamoto KYU Japan earthquake of April 16th 2016. The files have also been pushed to the repository for your convinience.

##### ERA5 data download script

Once the steps (found [here](https://confluence.ecmwf.int/display/CKB/How+to+install+and+use+CDS+API+on+Windows) *Note: The .cdspapirc file need to be created manually for Windows user) to create a login for ERA5's API have been completed atmospheric data can be easily downloaded the query(ies) string(s) in the `ERA5APICalls.py` script. Example of queries can be obtained from https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=form.

### Processing steps
To generate an interferogram Level 1 data is preprocessed using the following steps similar the example found at [this repository][sentinel-1-snappy].
 1.  Open source files
 2.  Select chosen subswath
 3.  Apply Geocoding
 4.  Burst co-registration using orbit and (digital elevation map) DEM
 5.  Estimate constant range and azimuth offsets for the whole image
 6.  Compute interferogram
 7.  Debursts the TOPSAR product
 8.  Compute and subtract topographic phase
 9.  Multilook: Average spectra over a set number of bands
 10. Apply terrain corrections to remove shadows and order terrain-induces artefacts


[asf]:               https://asf.alaska.edu/how-to/data-tools/data-tools/#bulk_download
[sentinel-1-snappy]: https://github.com/crisjosil/InSAR_Snappy
[atmospheric-corrections]: https://doi.org/10.1016/j.jappgeo.2009.03.010
[atmopspheric-data]: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=overview \Mu??oz Sabater J. (2021) was downloaded from the Copernicus Climate Change Service (C3S) Climate Data Store\

