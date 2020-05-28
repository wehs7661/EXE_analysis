import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def initialize():
    parser = argparse.ArgumentParser(
        description="This code analyze the log file from a expanded ensemble \
                    simulation and plot <dH/dL> as a function of lambda (L).")
    parser.add_argument('-l',
                        '--log',
                        type=str,
                        help='The file name of the log file.')

    arg_parse = parser.parse_args()

    return arg_parse

class dHdL:
    def __init__(self, logfile):
        """
        Obtains required parameters for the plotting

        Parameters
        ----------
        logfile : str
            The filename of the log file
        """
        f = open(logfile, 'r')
        lines = f.readlines()
        f.close()

        self.N_states = None
        self.coul, self.vdW = [], []

        for l in lines:
            if 'n-lambdas' in l:
                self.N_states = int(l.split('=')[1])

            if 'Started mdrun' in l:
                break

    def plot_dhdl(self, logfile):
        f = open(logfile, 'r')
        lines = f.readlines()
        f.close()
        lines.reverse()

        final_found = False
        line_n = 0
        N, vdw, coul, dG = np.zeros(self.N_states), np.zeros(self.N_states), np.zeros(self.N_states), np.zeros(self.N_states) 

        for l in lines:  # direcion: upward
            line_n += 1
            if 'dG(in kT)' in l:   # should find this line first
                print('yes')
                coul_idx = l.split().index('CoulL')
                vdw_idx = l.split().index('VdwL')
                dG_idx = l.split().index('dG(in') - 1

            if 'MC-lambda information' in l:
                final_found = True
                for i in range(self.N_states):
                    # start from lines[line_n - 3]
                    # 'MC-lambda information' is lines[line_n - 1]
                    N[i] = int(lines[line_n - 3 - i].split()[0])
                    vdw[i] = float(lines[line_n - 3 - i].split()[vdw_idx])
                    coul[i] = float(lines[line_n - 3 - i].split()[coul_idx])
                    dG[i] = float(lines[line_n - 3 - i].split()[dG_idx])

            if '  Step  ' in l and final_found is True:
                break

        lambda_all = coul + vdw    # x-axis
        plt.figure()
        plt.plot(lambda_all, dG)
        plt.show()

def main():

    rc('font', **{
    'family': 'sans-serif',
    'sans-serif': ['DejaVu Sans'],
    'size': 10
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='serif')

    args = initialize()

    if args.log is None:
        for file in os.listdir('.'):
            if file.endswith('.log'):
                args.log = file 
        if args.log is None:
            print('No log files provided or found! Please check if the directory is correct or specify the file name of the log file.')

    dHdL(args.log).plot_dhdl(args.log)