import argparse
import os
import sys
import numpy as np
from sklearn import linear_model

def initialize():
    parser = argparse.ArgumentParser(
        description="This code analyze a log file from a fixed-weight EXE \
                    simulation when the weights were not quite right, or \
                    more specifically, when the system gets stuck and could \
                    not sample all the states. The weights are iteratively \
                    adjusted by following g_k' = g_k + ln(N_(k-1)/N_k) if \
                    N_k > 0. For the region of scaling the electrostatic \
                    interactions, the weights are estimated such that ddG \
                    is approximately constant." )
    parser.add_argument('-l',
                        '--log',
                        type=str,
                        help='The file name of the log file.')
    
    args_parse = parser.parse_args()

    return args_parse

class WeightsIteration:
    def __init__(self, logfile):
        """
        Obtains required parameters for the adjustment.

        Parameters
        ----------
        logfile : str
            The filename of the log file
        """
        f = open(logfile, 'r')
        lines = f.readlines()
        f.close()

        self.fixed, self.N_states = None, None

        for l in lines:
            if 'n-lambdas' in l:
                self.N_states = int(l.split('=')[1])

            if 'lmc-stats' in l and self.fixed is None:
                if l.split('=')[1].split()[0] == 'no':
                    self.fixed = True
                else:
                    self.fixed = False
                    print('Wrong simulation type! This method is only applicable to adjusting the weights in a fixed-weight simulation.')
                    sys.exit()
            
            if 'Started mdrun' in l:
                break

    def adjust_weights(self, logfile):
        f = open(logfile, 'r')
        lines = f.readlines()
        f.close()
        lines.reverse()

        final_found = False
        line_n = 0
        N, count = np.zeros(self.N_states), np.zeros(self.N_states)
        G, dG = np.zeros(self.N_states), np.zeros(self.N_states)

        # First read in the data needed
        for l in lines:  # direction: upward
            line_n += 1
            if 'dG(in kT)' in l:   # should find this line first 
                coul_idx = l.split().index('CoulL')
                count_idx = l.split().index('Count')
                G_idx = l.split().index('G(in')
                dG_idx = l.split().index('dG(in') - 1

            if 'MC-lambda information' in l:
                final_found = True
                for i in range(self.N_states):
                    # start from lines[line_n - 3]
                    # 'MC-lambda information' is lines[line_n - 1]
                    N[i] = int(lines[line_n - 3 - i].split()[0])
                    count[i] = float(lines[line_n - 3 - i].split()[count_idx])
                    G[i] = float(lines[line_n - 3 - i].split()[G_idx])
                    dG[i] = float(lines[line_n - 3 - i].split()[dG_idx])

                    # find the state at which the coulombic force is totally off
                    if float(lines[line_n - 3 - i].split()[coul_idx]) == 1 and float(lines[line_n - 2 - i].split()[coul_idx]) < 1:
                        state_coul_off = int(lines[line_n - 3 - i].split()[0])

            if '  Step  ' in l and final_found is True:
                break

        # Start the weights iteration
        G_adjstd = np.zeros(self.N_states)
        n_zeros = 0     # number of 0 weights
        n_adjstd = 0    # number of states that have been adjusted (including 0 weights)
        model = linear_model.LinearRegression()
        for i in range(state_coul_off):
            if count[i] == 0:
                break
            if G[i] == 0:
                n_zeros += 1
                n_adjstd += 1
                G_adjstd[i] == 0
            else:
                n_adjstd += 1
                G_adjstd[i] = G[i] + np.log(count[i - 1] / count[i])

        # estimate the weights of the state with 0 count
        # In the region of scaling coulombic force, dG is approximately linear
        dG_input = dG[n_zeros - 1:n_adjstd]    # data for linear regression
        N_input = N[n_zeros - 1:n_adjstd]
        results = model.fit(N_input.reshape(-1,1), dG_input)
        n_prediction = state_coul_off - (len(dG_input) + n_zeros)

        for i in range(n_prediction):
            G_adjstd[n_adjstd + i] = G_adjstd[n_adjstd - 1 + i] + model.predict([[N_input[-1] + 1 + i]])

        # The remaining weights will be the same as the original ones
        for i in range(state_coul_off - 1, self.N_states):
            G_adjstd[i] = G[i]

        R2 = model.score(N_input.reshape(-1,1), dG_input)

        return G_adjstd, R2

def main():
    args = initialize()

    if args.log is None:
        for file in os.listdir('.'):
            if file.endswith('.log'):
                args.log = file
        if args.log is None:
            print('No log files provided or found! Please check if the directory is correct or specify the name of the log files.')

    WI = WeightsIteration(args.log)
    G_new, R2 = WI.adjust_weights(args.log)

    G_new_str = ''
    for i in range(len(G_new)):
        G_new_str += '%6.5f ' % G_new[i]

    result_str = '\nData analysis of the file %s:' % args.log
    print(result_str)
    print('=' * len(result_str))
    print('As a result of the weights iteration, the new weights for the next iteration are:')
    print(G_new_str)
    print('\nFor the region of scaling electrostatic interactions, if the count of a state is 0, we applied linear regression on dG to estimate the new weights.')
    print('\nThe R-squared value of the linear regression is: %s' % R2)
        