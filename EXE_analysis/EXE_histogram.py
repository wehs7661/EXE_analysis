#!/usr/bin/env python
"""
This is a Python script for analyzing the log file generated from the weights equilibration simulation. 
The script performs the following analysis:
1. Generate a plot of Wang-Landau incrementor as a function of time
2. Output the histogram at the last time frame of the expanded ensemble simulation
3. Print out equilibrated Wang-Landau weights which can directly be pasted to the .mdp file if needed. 
4. Estimate the uncertainty in free energy difference between the first and the last state from the final histogram.
"""

import sys
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import natsort


def initialize():
    """
    An initializing function
    """

    parser = argparse.ArgumentParser(
        description='This code analyzes the log file generated from weights equilibration simulations.')
    parser.add_argument('-f',
                        '--log',
                        nargs='+',
                        help='The filename(s) of the log file. Wildcards available.')
    parser.add_argument('-k',
                        '--keyword',
                        nargs='+',
                        help='The keyword used in the filename of the figure produced.')
    parser.add_argument('-t',
                        '--temp',
                        nargs='+',
                        help='The temperature of the simulation. This is for estimating the uncertainty of free energy difference.')

    args_parse = parser.parse_args()

    return args_parse


def get_equilibrated_info(logfile):
    """
    This function analyzes the log file and performs the following tasks:
    1. Output the data needed for generating a plot of Wang-Landau incrementor as a function of time
    2. Print out equilibrated Wang-Landau weights which can directly be pasted to the .mdp file. 
       For example, the equilibrated weights of a solvent simulation can often serve as a good initial
       guess of init-lambda-weights in the binding complex simulation.

    Parameters
    ----------

    logfile : str
        The filename of the log file

    Returns
    -------

    time : np.array
        An array showing the time frame that the WL incrementor is changed
    wl_incrementor : np.array
        The array of Wang-Landau incrementors
    n_states : int
        The number of alchemical states
    final_weights : str
        The final equilibrated lambda weights
    equil_time : float
        The time that the weights are equilibrated

    Example
    -------

    >>> get_equilibrated_weights('solvent_0.log')
    (array([23951.4, 24438.6, 25023. , 25581.4, 26125. , 27201.6, 27908.2,
       28462.6, 29076.6, 30122.8, 31729.2, 32799.8, 35157.6, 39770.4,
       41585. , 43209.6, 48910.6, 50988.2, 52211.8, 53420. , 54937.2]), array([0.1      , 0.08     , 0.064    , 0.0512   , 0.04096  , 0.032768 ,
       0.026214 , 0.020972 , 0.016777 , 0.013422 , 0.010737 , 0.0085899,
       0.0068719, 0.0054976, 0.004398 , 0.0035184, 0.0028148, 0.0022518,
       0.0018014, 0.0014412, 0.0011529]), 40, ' 0.00000 14.56500 28.46294 41.67213 54.11278 65.79066 76.92339 87.23114 96.83672 105.77297 113.89335 121.35474 128.16594 134.20094 139.62177 144.44218 148.59055 152.19318 155.31149 157.89021 160.03456 161.06497 161.96962 162.79330 163.44638 163.69086 163.86891 163.89583 163.73541 163.55287 163.23662 162.82605 162.19519 161.42480 160.44183 159.40266 158.37007 157.29657 156.56107 156.00758\n')
    """

    f = open(logfile, 'r')
    lines = f.readlines()
    f.close()

    wl_incrementor = []
    step = []
    line_n = 0
    counter = 0    # a counter for examining if weights are equilibrated

    for l in lines:
        line_n += 1
        if ' dt ' in l:
            time_step = float(l.split('=')[1])

        if 'Wang-Landau incrementor is:' in l and not wl_incrementor:
            # not wl_incrementor returns True if the list is empty
            wl_incrementor.append(float(l.split(':')[1]))  # initial weight

        if 'weights are now: ' in l:
            # this line only appears before the incrementor is about to change
            weights_line = l.split(':')
            step.append(int(weights_line[0].split()[1]))

            # the info of the incrementor typically appear in line line_n + 7 but it depends
            search_lines = lines[line_n + 1: line_n + 7]
            for l_search in search_lines:
                if 'Wang-Landau incrementor is:' in l_search:
                    wl_incrementor.append(float(l_search.split(':')[1]))

        if 'Weights have equilibrated' in l:
            counter += 1
            equil_step = l.split(':')[0].split()[1]    # the step that the weights are equilibrated
            # the info of the incrementor typically appear in line line_n - 3 but it depends
            search_lines2 = lines[line_n - 8: line_n]
            for l_search in search_lines2:
                if 'weights are now: ' in l_search:
                    final_weights = l_search.split(':')[2]
                    # no need to change to float since this is only for easy pasting to the .mdp file
                    float_weights = [float(i) for i in l_search.split(':')[2].split()]
            break

    time = np.array(step) * time_step
    equil_time = float(equil_step) * time_step / 1000   # units: ns

    if counter == 0:
        print('The weights have not equilibrated.')
        print('The Wang-Landau incrementor at the last time frame (%5.3f ns) is %s.' %
              (time[-1] / 1000, str(wl_incrementor[-1])))
        print('Check the log file for more information.')
        sys.exit()

    n_states = len(float_weights)
    wl_incrementor = np.array(wl_incrementor)

    return time, wl_incrementor, n_states, final_weights, equil_time 


def get_final_histogram(n_states, logfile, temp):
    """
    This function analyzes the log file and performs the following tasks:
    1. Output the counts of each lambda state at the last time frame (for plotting histogram)
    2. Estimate the uncertainty of free energy difference from the final histogram

    Paraneters
    ----------

    n_states : int
        Number of lambda states
    logfile : str
        The filename of the log file

    Returns
    -------

    counts : np.array
        The counts of each lambda state

    Example
    -------
    >>> get_final_histogram(40, 'solvent_0.log')
    [8678. 8437. 8680. 9007. 8606. 7642. 8269. 7878. 7689. 7906. 7451. 7416.
 7939. 7470. 7540. 7858. 7664. 7423. 7527. 7322. 7325. 7538. 7173. 7034.
 6943. 6910. 6935. 6805. 6463. 6371. 6249. 6425. 6353. 6618. 6789. 6810.
 6426. 6408. 6675. 6271.]
    """

    f = open(logfile, 'r')
    lines = f.readlines()
    f.close()
    lines.reverse()     # from this point, lines has been reverse

    line_n = 0
    counts = np.zeros(n_states)
    for l in lines:
        line_n += 1
        if 'MC-lambda information' in l:
            for i in range(n_states):
                # start from lines[line_n - 3]
                counts[i] = float(lines[line_n - 3 - i].split()[5])
            break

    kb = 1.38064852E-23                               # Boltzmann constant
    Na = 6.0221409E23                                 # Avogadro's number
    error = np.abs(np.log(counts[0] / counts[-1]))    # dimensionless error

    if temp is None:
        print('The uncertainty of the free energy difference is %5.3f kT.' % error)
        temp = 298.15    # default
        error *= (kb * Na * temp / 1000) * 0.23900573613
        print('Or at 298.15K, the uncertainty is %5.3f kcal/mol' % error)
    else:
        error *= (kb * Na * float(temp) / 1000) * \
            0.23900573613       # unit: kcal/mol
        print('The uncertainty of the free energy difference is %5.3f kcal/mol.' % error)

    return counts


if __name__ == '__main__':
    args = initialize()

    rc('font', **{
        'family': 'sans-serif',
        'sans-serif': ['DejaVu Sans'],
        'size': 10
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='serif')

    if isinstance(args.log, str):                      # the case of only one input
        args.log = list(args.log)
    if isinstance(args.log, str) and '*' in args.log:  # the case of using wildcard
        args.log = natsort.natsorted(glob.glob(args.log), reverse=False)
    if isinstance(args.keyword, str):
        args.keyword = list(args.keyword)
    if isinstance(args.temp, str):
        args.temp = list(args.temp)

    # Check if the default keyword should be used
    if args.keyword is None:
        args.keyword = [logname.split('.')[0] for logname in args.log]

    for i in range(len(args.log)):
        if len(args.log) > 1:
            print('Output data of the log file: %s' % args.log[i])

        [time, wl_incrementor, n_states, weights, equil_time] = get_equilibrated_info(args.log[i])
        time = time / 1000  # convert from ps to ns

        if args.temp is None:
            counts = get_final_histogram(n_states, args.log[i], args.temp)
        else:
            counts = get_final_histogram(n_states, args.log[i], args.temp[i])

        print('The weights are equilibrated at %5.3f ns' % equil_time)
        print('The final weights are: ', weights)

        # Plot WL incrementor as a function of time
        plt.figure()  # ready to plot!
        plt.step(time, wl_incrementor)
        plt.xlabel('Time (ns)')
        plt.ylabel('Wang-Landau incrementor ($ k_{B} T$)')
        plt.title('Wang-Landau incrementor as a function of time')
        plt.grid()
        plt.savefig('WL_t_%s.png' % args.keyword[i], dpi=600)
        plt.show()

        # Plot the final histogram
        plt.figure()
        plt.bar(np.arange(1, 41), height=counts)
        plt.xlabel('States')
        plt.ylabel('Counts')
        plt.title('The final histogram of the simulation')
        plt.grid()
        plt.savefig('Final_hist_%s.png' % args.keyword[i], dpi=600)
        plt.show()
