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
    parser.add_argument('-n',
                        '--ignore',
                        type=int,
                        help='The number of the first few states keeping \
                             the same weight. Default: 0.')
    
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

        self.fixed, self.N_states, self.final_t = None, None, None
        self.furthest, self.R2, self.RMSD = None, None, None
        self.n_input, self.n_prediction = None, None

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

    def adjust_weights(self, n, logfile):
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
                self.final_t = float(lines[line_n - 2].split()[1])  # in ps
                break

        # Find the furthest visited state
        for i in range(len(count) - 1):
            if count[i] != 0 and count[i + 1] == 0:
                self.furthest = int(N[i])
            if self.furthest is None:
                self.furthest = self.N_states

        # Calculate RMSD in count
        # Note that RMSD here is proportional to the time, so we have to normalize RMSD with t
        count_avg = sum(count)/len(count)   # theoretical/ideal count for each state  --> ref
        self.RMSD = np.sqrt(sum((np.array(count) - count_avg) ** 2) / len(count)) / self.final_t

        # Start the weights iteration
        G_adjstd = np.zeros(self.N_states)
        n_zeros = 0     # number of 0 weights
        n_coul_adjstd = 0    # number of states that have been adjusted (including 0 weights) in the "coulobic region"
        # the maximum of n_coul_adjustd = state_coul_off

        # All the states with non-zero counts should be adjusted
        # no matter which region the states are
        for i in range(self.N_states):
            if G[i] == 0:
                n_zeros += 1
                n_coul_adjstd += 1
                G_adjstd[i] == 0

            if G[i] !=0 and count[i] != 0:
                if i + 1 <= state_coul_off:
                    n_coul_adjstd += 1
                if count[i - 1] == 0:
                    G_adjstd[i] = G[i]
                else:
                    G_adjstd[i] = G[i] + np.log(count[i - 1] / count[i])
                    #G_adjstd[i] = G[i] + np.log(count[0] / count[i])
            
        # Now deal with the states with 0 count
        # 1. The weight will be predicted by linear regression if the state is for turning off coulubic forces
        # 2. The weight wlll remain the same as the weight of input mdp file if the state is for turning off vdw forces
        
        # Case 1: In the region of scaling coulombic force, dG is approximately linear
        
        dG_input = dG[n_zeros - 1:n_coul_adjstd]    # data for linear regression
        N_input = N[n_zeros - 1:n_coul_adjstd]
        self.n_input = len(N_input)    # number of input data points for linear regression
        
        # Note: n_zeros + self.n_input + self.n_prediction = stat_coul_off
        self.n_prediction = state_coul_off - (self.n_input + n_zeros)
        
        if self.n_prediction > 0:
            model = linear_model.LinearRegression()
            model.fit(N_input.reshape(-1,1), dG_input)
            self.R2 = model.score(N_input.reshape(-1,1), dG_input)
            for i in range(self.n_prediction):
                G_adjstd[n_coul_adjstd + i] = G_adjstd[n_coul_adjstd - 1 + i] + model.predict([[N_input[-1] + 1 + i]])

        # Case 2: The remaining weights will be the same as the original ones
        for i in range(state_coul_off, self.N_states):
            if count[i] == 0:
                G_adjstd[i] = G[i]

        # if n > 0, the weight should remain the same, so modify back if needed
        if n > 0:
            for i in range(n):
                G_adjstd[i] = G[i]

        return G_adjstd

def main():
    args = initialize()

    if args.ignore is None:
        args.ignore = 0

    if args.log is None:
        for file in os.listdir('.'):
            if file.endswith('.log'):
                args.log = file
        if args.log is None:
            print('No log files provided or found! Please check if the directory is correct or specify the name of the log files.')

    WI = WeightsIteration(args.log)
    G_new = WI.adjust_weights(args.ignore, args.log)

    G_new_str = ''
    for i in range(len(G_new)):
        G_new_str += '%6.5f ' % G_new[i]

    result_str = '\nData analysis of the file %s:' % args.log
    print(result_str)
    print('=' * (len(result_str) - 1))  # len(result_str) includes \n 
    print('Simulation length: %s ps (furthest state sampled: %s)' % (WI.final_t, WI.furthest))
    print('RSMD in count (reference: counts averaged over all the states): %s\n' % WI.RMSD)
    print('As a result of the weights iteration, the new weights for the next iteration are:')
    print(G_new_str)
    print('\nFor the region of scaling electrostatic interactions, if the count of a state is 0, we applied linear regression on dG to estimate the new weight of that state.')
    
    if WI.n_prediction <= 0:
        print('Since all the states for turning off electrostatic were sampled, no linear regression is needed.')
    else:
        print('\nNumber of data points for linear regression: %s' % WI.n_input)
        print('The R-squared value of the linear regression is: %s' % WI.R2)
        