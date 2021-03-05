# GeoHydroHub Examples

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cheginit/geohydrohub-examples/HEAD)
[![CI](https://github.com/cheginit/geohydrohub-examples/actions/workflows/test.yml/badge.svg)](https://github.com/cheginit/geohydrohub-examples/actions/workflows/test.yml)

[GeoHydroHub](https://pygeohydro.readthedocs.io) is a software stack consists of six
Python libraries and is designed to aid in watershed analysis through web services.
Currently, they only includes hydrology and climatology data
within the US. Some of the major capabilities of this software stack are as follows:

* Easy access to many web services for subsetting data and returning the requests as masked
  xarrays or GeoDataFrames.
* Splitting large requests into smaller chunks under-the-hood since web services usually limit
  the number of items per request. So the only bottleneck for subsetting the data
  is the local machine memory.
* Navigating and subsetting NHDPlus database (both meduim- and high-resolution) using web services.
* Cleaning up the vector NHDPlus data, fixing some common issues, and computing vector-based
  accumulation through a river network.
* A URL inventory for some of the popular (and tested) web services.
* Some utilities for manipulating the data and visualization.
