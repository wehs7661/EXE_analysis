import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

def initialize():
    parser = argparse.ArgumentParser(
        description='This code analyzes the dhdl file generated from expanded \
                    ensemble simulations and plot the state as a function of \
                    time. With the log file given, the time that the weights \
                    are still equilibrating will be truncated.')
    parser.add_argument('-i',
                        '--dhdl',
                        help='The file name of the dhdl file.')
    parser.add_argument('-l',
                        '--log',
                        help='The file name of the log file, which is \
                             required only if the user want to truncate \
                             the stage in which the weights are still \
                             equilibrating.')
    parser.add_argument('-f',
                        '--freq',
                        help='The frequency (every certain steps) to extract \
                              the data to plot. Default: 1 (step).')
    parser.add_argument('-k',
                        '--keyword',
                        help='The keyword used in the filename of the figure \
                             produced.')
    parser.add_argument('-s',
                        '--shift',
                        default=True,
                        action='store_false',
                        help='Whether to shift the first time frame back to \
                             zero if the dhdl file is truncated. Will make no \
                             differences if no log file is provided.')

    args_parse = parser.parse_args()
    
    return args_parse


def equilibration_time(logfile):
    """
    This function analyzes the log file find the time interval for dhdl 
    truncation. Specifically, the interval starts from the time frame that the
    weights are equilibrated and ends at the last time frame of the log file.

    Parameters
    ----------
    logfile : str
        The filename of the log file

    Returns
    -------
    equil_last : float
        The time that the weights are equilibrated
    """

    f = open(logfile, 'r')
    lines = f.readlines()
    f.close()

    line_n = 0
    equil = False    # to examine if weights are equilibrated

    # First find the time frame the the weights are equilibrated
    for l in lines:
        line_n += 1
        if ' dt ' in l:
            time_step = float(l.split('=')[1])

        if 'Weights have equilibrated' in l:
            equil = True
            equil_last_step = l.split(':')[0].split()[1]    
            # the step that the weights are equilibrated
            break

    equil_last = float(equil_last_step) * time_step   # units: ps

    if equil is False:
        print('The weights have not equilibrated.')
        print('Check the log file for more information.')
        sys.exit()

    return equil_last


def make_trunacated_dhdl(dhdl, equil_last):
    """
    This function reads in a dhdl file and truncate the part that the weights are still equilibrating
    """
    prefix = dhdl.split('.xvg')[0]
    out_name = prefix + '_truncated.xvg'
    outfile = open(out_name, 'w')
    with open(dhdl) as infile:
        for line in infile:
            if line[0] == '#' or line[0] == '@':
                outfile.write(line)
            else:
                time = float(line.split()[0])
                if time > equil_last:
                    outfile.write(line)
    

def data_extraction(dhdl, freq):
    """Reads a dhdl file and returns the states and times.

    Parameters
    -----------
    dhdl : str
        The name of the dhdl file to read.
    freq : float
        The 
    
    Returns
    ----------
    states : np.array
        The lambda states.
    times : np.array
        The times. 
    """

    f = open(dhdl,'r')
    lines = f.readlines()
    f.close()

    if freq is None:
        freq = 1

    # Figure out the number of state from the metatexts
    N_states = 0
    n = 0     # line number of the metatext
    for l in lines:  
        n += 1   # so l is lines[n - 1]
        if l[0] == '@' and ' to ' in l:
            N_states += 1
            if lines[n][0] != '@':
                break
            
    # Extract the state-time data
    state, time = [], []    # units of time: ps
    i = 0     # line number (excluding metatexts)
    for l in lines:
        if l[0] != '#' and l[0] != '@':
            i += 1
            if i % freq == 0:
                time.append(float(l.split()[0]))
                state.append(int(float(l.split()[1])))
    state = np.array(state)
    time = np.array(time) / 1000      # units: ns

    return time, state, N_states


def main():
    args = initialize()

    if args.keyword is None:
        args.keyword = args.dhdl.split('_dhdl')[0].split('/')[-1] 

    dhdl_to_plot = args.dhdl    # dhdl file to plot (will be replaced with the truncated file if the log file is provided)
    png_name = 'state_plot_%s_untruncated.png' % args.keyword  # will also be replaced given a log file
    
    if args.log is None:
        print('Nothing to truncate: the log file is not given!')
        title = 'Exploration of states as a function of time (untruncated simulation)'
    else:
        equil_last = equilibration_time(args.log)
        print('The weights are equilibrated at %5.3f ns' % (equil_last / 1000))
        print('Accordingly, the first %5.3f ns of the simulation will be truncated.' % (equil_last / 1000))
        make_trunacated_dhdl(args.dhdl, equil_last)
        dhdl_to_plot = args.dhdl.split('.xvg')[0] + '_truncated.xvg'
        png_name = 'state_plot_%s_truncated.png' % args.keyword

    time, state, N_states = data_extraction(dhdl_to_plot, args.freq)
    print('Note: the length of untruncated simulation is %6.3f ns' % time[-1])

    if args.shift and args.log is not None:
        time -= time[0]
        title = 'Exploration of states as a function of time with equilibrated weights'
    
    if not args.shift and args.log is not None:
        title = 'Exploration of states as a function of time'

    rc('font', **{
        'family': 'sans-serif',
        'sans-serif': ['DejaVu Sans'],
        'size': 10
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='serif')
    
    plt.plot(time, state)
    plt.title(title)
    plt.xlabel('Time (ns)')
    plt.ylabel('State')
    plt.ylim([0, N_states])
    plt.grid(True)
    plt.savefig(png_name)
    plt.show()