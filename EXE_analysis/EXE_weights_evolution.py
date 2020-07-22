
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import rc 
from .EXE_histogram import LogInfo

def initialize():
    parser = argparse.ArgumentParser(
        description='This code plot the time evolution of the weights in \
                    EXE or MetaD-EXE.')
    parser.add_argument('-l',
                        '--log',
                        help='The name of the log file.')

    args_parse = parser.parse_args()

    if args_parse.log is None:
        for file in os.listdir('.'):
            if file.endswith('.log'):
                args_parse.log = file

    return args_parse

class WeightsEvolution(LogInfo):
    def __init__(self, logfile):
        LogInfo.__init__(self, logfile)
        self.prefix = logfile.split('.')[0]

    def get_WL_weights(self, logfile):
        f = open(logfile, 'r')
        lines = f.readlines()
        f.close()

        dG = []
        for i in range(self.N_states - 1):
            dG.append([])

        n_lines = 0
        n = 0
        for l in lines:
            n_lines += 1

            if 'dG(in kT)' in l: # lines[n_lines - 1]
                n += 1
                for i in range(self.N_states- 1):
                    if '<<' not in lines[n_lines + i]:
                        dG[i].append(float(lines[n_lines + i].split()[-1]))
                    else:
                        dG[i].append(float(lines[n_lines + i].split()[-2]))
        
        # compute the free energy difference of each state w.r.t the first state (index: 0)
        dG_0 = []
        for i in range(self.N_states - 1):
            if i == 0:
                dG_0.append(np.array(dG[0]))
            else:
                dG_0.append(np.array(dG[i]) + dG_0[-1])

        return dG_0

    def plot_dG(self, dG_0):
        time = (np.arange(len(dG_0[0])) + 1) * self.nstlog * self.dt / 1000  # ns
        plt.figure()
        for i in range(self.N_states - 1):
            plt.plot(time, dG_0[i], label=f'State {i + 2}')
        plt.xlabel('Time (ns)')
        plt.ylabel('WL weight difference w.r.t the 1st state (kT)')
        plt.title('Wang-Landau weight difference as a function of time')
        plt.legend(ncol=2)
        plt.grid()
        plt.savefig(f'{self.prefix}_weights_evolution.png', dpi=600)
        plt.show()
        
def main():
    args = initialize()

    rc('font', **{
        'family': 'sans-serif',
        'sans-serif': ['DejaVu Sans'],
        'size': 10
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='serif')

    WE = WeightsEvolution(args.log)
    dG_0 = WE.get_WL_weights(args.log)
    WE.plot_dG(dG_0)




