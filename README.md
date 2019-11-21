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
### 1. `EXE_histogram.py`: Data analysis of Wang-Landau incrementor and weighting factors
- To check the inputs, ouputs and the description of the code, run `EXE_histogram -h`.
- The program could be run without specifying any parameters (specifically, just `EXE_histogram`). Under this situation, the program will analyze all the log files (if any) in the current folder with the length of the simulation for weights average calculation as the default (20 ns).

- To perform data analysis on `solvent_0.log` using `EXE_histogram.py`, run `EXE_histogram -l solvent_0.log`. If the weights were equilibrated in the simulation, then three figures named with defaults will be generated, including `WL_t_solvent_0.png`, `Equil_hist_55ns_solvent_0.png`, and `Final_hist_114ns_solvent_0.png` as shown below, given that the last histogram was generated at about 114 ns and the weights were equilibrated at about 55 ns. (The log files are not provided here since the files are too large.) If the weights had not equilibrated in the simulation, then only the final histogram will be generated.

- Wildcards are available. Therefore, say that there are only two log files in the current folder, including `solvent_0.log` and `solvent_1.log`, to perform data analysis on both files, we can run either `EXE_histogram -l solvent_*.log` or `EXE_histogram -l *.log` (or `EXE_histogram` or `EXE_histogram -l solvent_0.log solvent_1.log` without using wildcards). If the weights were equilibrated in both simulations, then 6 figures will be generated.

- Sample STDOUT message as the output of data analysis performed on both `solvent_0.log` and `solvent_1.log`:
```
The log file(s) to be analyzed: solvent_0.log, and solvent_1.log.
Length of the simulation for the weights average calculation: 20 ns.

Data analysis of the file solvent_0.log:
========================================
The uncertainty of the free energy difference is 0.174 kT.

Or at the simulation temperature (298.15 K), the uncertainty is 0.103 kcal/mol

The weights were equilibrated at 54.937 ns

The final weights are:
  0.00000 14.56500 28.46294 41.67213 54.11278 65.79066 76.92339 87.23114 96.83672 105.77297 113.89335 121.35474 128.16594 134.20094 139.62177 144.44218 148.59055 152.19318 155.31149 157.89021 160.03456 161.06497 161.96962 162.79330 163.44638 163.69086 163.86891 163.89583 163.73541 163.55287 163.23662 162.82605 162.19519 161.42480 160.44183 159.40266 158.37007 157.29657 156.56107 156.00758

The average weights of the last 20 ns of the simulation before the weights are equilibrated (from 34.938 to 54.938 ns) are:
  0.0 14.64648 28.57863 41.79821 54.28418 66.01599 77.11028 87.46918 97.07739 106.00707 114.12933 121.61515 128.38682 134.43235 139.8606 144.69019 148.83927 152.46849 155.52978 158.09579 160.25379 161.2892 162.18532 162.97901 163.68603 163.94997 164.0765 164.11113 164.01516 163.84075 163.53011 163.11137 162.47901 161.7277 160.71876 159.73908 158.73955 157.66436 156.94407 156.40345 

Data analysis of the file solvent_1.log:
========================================
The uncertainty of the free energy difference is 0.246 kT.

Or at the simulation temperature (298.15 K), the uncertainty is 0.146 kcal/mol

The weights were equilibrated at 50.864 ns

The final weights are:
  0.00000 14.60170 28.43523 41.58579 54.08459 65.87537 76.98775 87.34401 97.00613 105.89536 114.04492 121.53658 128.34734 134.41112 139.80858 144.59074 148.78474 152.36356 155.42416 158.01062 160.16965 161.19495 162.09843 162.89775 163.52808 163.75729 163.82927 163.90205 163.76106 163.58467 163.28349 162.79807 162.18027 161.25948 160.28998 159.15817 158.14124 157.05890 156.27008 155.67805

The average weights of the last 20 ns of the simulation before the weights are equilibrated (from 30.866 to 50.866 ns) are:
  0.0 14.64105 28.54121 41.76535 54.25945 66.08005 77.21117 87.55368 97.19715 106.09082 114.29638 121.77106 128.55119 134.60321 140.06286 144.82714 149.00974 152.61477 155.64057 158.26884 160.41002 161.40118 162.33332 163.13015 163.71199 163.96416 164.08451 164.09988 163.98061 163.81325 163.52252 163.01236 162.4054 161.54296 160.49669 159.40129 158.34479 157.25842 156.47471 155.92581 

2 file(s) analyzed.
Total time elapsed (including plotting): 20.138840436935425 seconds.
```

- Sample figures as the output of data analysis performed on `solvent_0.log`:
<img src="EXE_analysis/examples/WL_t_solvent_0.png" width="400"/> <br/>
<img src="EXE_analysis/examples/Equil_hist_55ns_solvent_0.png" width="400"/> <img src="EXE_analysis/examples/Final_hist_114ns_solvent_0.png" width="400"/>

### 2. `EXE_state_plot.py`: Data analysis of the exploration of the state
- To check the inputs, ouputs and the description of the code, run `EXE_state_plot -h`.


### 3. Alchemical free energy calculations
- With `dhdl` files, we are able to calculate the free energy difference between two thermodyanmic states. This calculation is now temporarily relying on the package [alchemical-analysis](https://github.com/MobleyLab/alchemical-analysis), but will move to [alchemlyb](https://github.com/alchemistry/alchemlyb) in the future.
- Using `alchemical-analysis`, to calculate the energy difference given one or more `dhdl` files with prefix as `complex`, run `alchemical_analysis -p complex_1 -t 298.15 -u kcal -w -g -x`. 
- To check the meaning of each useful flag used in `alchemical-analysis`, run `alchemical-analysis -h`.


## Copyright

Copyright (c) 2019, Wei-Tse Hsu


## Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
