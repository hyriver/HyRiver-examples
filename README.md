[![Logo](https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/hyriver_logo_text.png)](https://github.com/cheginit/HyRiver)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cheginit/hyriver-examples/HEAD)
[![CI](https://github.com/cheginit/hyriver-examples/actions/workflows/test.yml/badge.svg)](https://github.com/cheginit/hyriver-examples/actions/workflows/test.yml)

# Examples Notebooks

[HyRiver](https://hyriver.readthedocs.io) is a software stack consisting of six
Python libraries that are designed to aid in watershed analysis through web services.
Currently, this project only includes hydrology and climatology data
within the US. Some of the major capabilities of HyRiver are as follows:

* Easy access to many web services for subsetting data on server-side and returning the requests
  as masked Datasets or GeoDataFrames.
* Splitting large requests into smaller chunks, under-the-hood, since web services often limit
  the number of features per request. So the only bottleneck for subsetting the data
  is your local machine memory.
* Navigating and subsetting NHDPlus database (both medium- and high-resolution) using web services.
* Cleaning up the vector NHDPlus data, fixing some common issues, and computing vector-based
  accumulation through a river network.
* A URL inventory for some of the popular (and tested) web services.
* Some utilities for manipulating the obtained data and their visualization.
