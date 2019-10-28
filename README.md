EXE_analysis
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/EXE_analysis.png)](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/EXE_analysis)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/REPLACE_WITH_APPVEYOR_LINK/branch/master?svg=true)](https://ci.appveyor.com/project/REPLACE_WITH_OWNER_ACCOUNT/EXE_analysis/branch/master)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/EXE_analysis/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/EXE_analysis/branch/master)

## Description 
`EXE_analysis` is a Python package of data analysis tools for expanded ensemble (EXE) simulations.

## Installation
All the Python scripts in this package are written in Python 3. Currently the package can be installed by following the commandes below:
```
git clone https://github.com/wehs7661/EXE_analysis.git
cd EXE_analysis
pip install -e .
```

## Examples and usage
### `EXE_histogram.py`
- To check the inputs, ouputs and the description of the code, run `EXE_histogram -h`
- To perform data analysis on `complex_0.log` using `EXE_histogram.py` and print the results to `results_histo.txt`, run `EXE_histogram -l complex_0.log`. Then two figures named with defaults will be generated, including `WL_t_complex_0.png` and `Final_hist_complex_0.png` as shown below.
- To perform data analysis on both `complex_0.log` and `complex_1.log` using `EXE_histogram.py` and print the results to `results_histo_multi.txt`, run `EXE_histogram -l complex_*.log` (or `EXE_histogram -l complex_0.png complex_1.png` without using wildcards). 4 figures will be generated.
- Example figures as the output of data analysis performed on `complex_0.png`:
![alt-text-1](WL_t_complex_0.png "title-1") ![alt-text-2](Final_hist_complex_0.png "title-2")

### `EXE_stat_plot.py`

## Copyright

Copyright (c) 2019, Wei-Tse Hsu


## Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
