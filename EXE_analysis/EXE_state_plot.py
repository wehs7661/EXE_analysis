import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

def initialize():
    parser = argparse.ArgumentParser(
        description='This code analyzes the dhdl file generated from expanded ensemble simulations and plot the state as a function of time.')
    parser.add_argument('-i',
                        '--dhdl',
                        help='The file name of the dhdl file')
    parser.add_argument('-f',
                        '--freq',
                        help='The frequency (every certain steps) to extract the data to plot')
    parser.add_argument('-k',
                        '--keyword',
                        help='The keyword used in the filename of the figure produced.')

    args_parse = parser.parse_args()
    
    return args_parse

def make_trunacated_dhdl(dhdl):
    """
    This 


    """
    return 

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

    return time, state 

def main():
    args = initialize()

    if args.keyword is None:
        args.keyword = args.dhdl.split('_dhdl')[0].split('/')[-1] 

    time, state = data_extraction(args.dhdl, args.freq)

    rc('font', **{
        'family': 'sans-serif',
        'sans-serif': ['DejaVu Sans'],
        'size': 10
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='serif')
    
    plt.plot(time, state)
    plt.title('Exploration of states as a function of time')
    plt.xlabel('Time (ns)')
    plt.ylabel('State')
    plt.grid(True)
    plt.savefig('state_plot_%s.png' % args.keyword)
    plt.show()
