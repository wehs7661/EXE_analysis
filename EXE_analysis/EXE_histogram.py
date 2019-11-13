"""
This is a Python script for analyzing the log file generated from the expanded
ensemble simulation. The script performs the following analysis:
1. Generate a plot of Wang-Landau incrementor as a function of time.
2. Show and save the histogram at the last time frame of the simulation.
3. Show and save the histogram at the time that the weights are equilibrated.
4. Print out equilibrated Wang-Landau weights which can directly be pasted to
    the .mdp file if needed.
5. Print out the averaged Wang-Landau weights for a user-definied portion of
    the last iteration.
6. Estimate the uncertainty in free energy difference between the first and
    the last state from the final histogram.
"""

import sys
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import natsort
import time as timer
# time will be a variable name used below


def initialize():
    """
    An argument parser as an initializing function
    """

    parser = argparse.ArgumentParser(
        description='This code analyzes the log file generated from expanded \
                    ensemble simulations.')
    parser.add_argument('-l',
                        '--log',
                        nargs='+',
                        help='The filename(s) of the log file. Wildcards \
                            available.')
    parser.add_argument('-f',
                        '--frac',
                        help='The fraction which is used to calcuate the \
                            average weights. A fraction of 0.25 (default) \
                            means average the weights of the last 25 percent \
                            of the steps before weights equilibration.')
    parser.add_argument('-k',
                        '--keyword',
                        nargs='+',
                        help='The keyword used in the filename of the figure \
                            produced.')

    args_parse = parser.parse_args()

    return args_parse


class LogInfo:
    def __init__(self, logfile):
        """
        Gets the needed parameters from the log file and set up relevant
        attributes to run the analysis

        Parameters
        ----------
        logfile : str
            The filename of the log file
        """
        f = open(logfile, 'r')
        lines = f.readlines()
        f.close()

        self.dt, self.cutoff, self.wl_scale, self.N_states, self.fixed, \
            self.init_wl, self.temp, self.start = None, None, None, None, \
            None, None, None, None
        line_n = 0

        for l in lines:
            line_n += 1
            if 'dt  ' in l and self.dt is None:
                self.dt = float(l.split('=')[1])

            if 'weight-equil-wl-delta' in l and self.cutoff is None:
                self.cutoff = float(l.split('=')[1])

            if 'wl-scale' in l and self.wl_scale is None:
                self.wl_scale = float(l.split('=')[1])

            if 'n-lambdas' in l and self.N_states is None:
                self.N_states = int(l.split('=')[1])

            if 'lmc-stats' in l and self.fixed is None:
                if l.split('=')[1] == 'no':
                    self.fixed = True
                else:
                    self.fixed = False

            if 'init-wl-delta' in l and self.init_wl is None:
                self.init_wl = float(l.split('=')[1])

            if 'ref-t' in l and self.temp is None:
                self.temp = l.split(':')[1]

            if 'Started mdrun' in l:
                self.start = line_n
                # the line number that the simulation got started
                break


class EXEAnalysis(LogInfo):
    """
    A class inheriting from the class LogInfo
    """

    def __init__(self, logfile):
        """
        Sets up the properties of the istance of EXEAnalysis
        """
        self.equil = False  # for examining if weights are equilibrated
        self.equil_time = None
        LogInfo.__init__(self, logfile)
        # So any instance of EXEAnalysis will have all the properties that a
        # LogInfo instance has

    def get_equil_info(self, logfile, frac):
        """
        This function analyzes the log file and performs the following tasks:
        1. Output the data needed for plotting Wang-Landau incrementor as a
            function of time
        2. Print out equilibrated Wang-Landau weights which can directly be
            pasted to the .mdp file.
        3. Estimate the uncertainty of free energy difference from the
            equilibrted histogram

        Parameters
        ----------
        logfile : str
            The filename of the log file
        frac : float
            The fraction of the time before the weights were equilibrated.
            Used for calculating the average weights.

        Returns
        -------
        time : np.array
            An array containing the time frames when WL incrementor changed
        wl_incrementor : np.array
            The array of Wang-Landau incrementors
        final_weights : str
            The final equilibrated lambda weights
        weights_avg :
            The average weights of the last portion of steps before the
            weights are equilibrated
        equil_time : float
            The time that the weights are equilibrated
        """

        # ======================= 1. Initialization ==========================
        f = open(logfile, 'r')
        lines = f.readlines()
        f.close()

        wl_incrementor = []
        step = []
        line_n = self.start   # start counting from self.start
        weights = []          # a list of weights at a certain step
        weights_all = []      # a list of lists of weights at different time
        # steps

        # ============= 2. Extract wl-incrementors and weights ===============
        for l in lines[self.start:]:  # skip the metadata
            line_n += 1
            if 'Wang-Landau incrementor is' in l and not wl_incrementor:
                # not wl_incrementor returns True if the list if
                wl_incrementor.append(self.init_wl)  # initial wl-incrementor

            if 'weights are now: ' in l:
                # this line only appears before the Wang-Landau incrementor is
                # about to change
                weights_line = l.split(':')
                step.append(int(weights_line[0].split()[1]))

                search_lines = lines[line_n + 1: line_n + 15]
                # the info of the incrementor typically appears in
                # line line_n + 15 but it depends

                for l_search in search_lines:
                    if 'Wang-Landau incrementor is:' in l_search:
                        wl_incrementor.append(float(l_search.split(':')[1]))

            # Create the list to calculate the average weights (first find the
            # point of the last update before the weights are equilibrated)
            if 'Wang-Landau incrementor is' in l and len(wl_incrementor) > 0 \
                    and (wl_incrementor[-1] > self.cutoff) \
                    and (wl_incrementor[-1] * self.wl_scale <= self.cutoff):
                for i in range(self.N_states):
                    weights.append(float(lines[line_n + i + 1].split()[6]))
                weights_all.append(weights)
                weights = []

            if 'Weights have equilibrated' in l:
                self.equil = True

                equil_step = l.split(':')[0].split()[1]
                # the step that the weights are equilibrated
                # the info of the incrementor typically appear in
                # line line_n - 8 but it depends

                search_weights = lines[line_n - 8: line_n]
                # lines for searching weights

                for l_search in search_weights:
                    if 'weights are now: ' in l_search:
                        final_weights = l_search.split(':')[2]
                        # This would include '\n' Note: there is no need to
                        # change to float since this is only for easy pasting
                        # to the .mdp file

                search_counts = lines[line_n - (30 + self.N_states): line_n]
                # 30 is the approximate number of lines of metadata between
                # the counts data and the position at which the weights are
                # found equilibrated

                equil_counts = []

                search_n = line_n - (30 + self.N_states)
                # line number of the lines for searching

                for l_search in search_counts:
                    search_n += 1
                    if 'MC-lambda information' in l_search:
                        for i in range(self.N_states):
                            # start from lines[search_n + 2]
                            equil_counts.append(float(lines[search_n + 2 + i].
                                                split()[5]))
                break

        wl_incrementor = np.array(wl_incrementor)
        time = np.array(step) * self.dt
        self.equil_time = float(equil_step) * self.dt / 1000

        # ========== 3. Exit if the weights have not equilibrated ============
        if self.equil is False:
            N_update = int(np.ceil(np.log(self.cutoff / wl_incrementor[-1]) /
                                   np.log(self.wl_scale)))
            print('The weights have not equilibrated.')
            print('The last time frame that the Wang-Landau increTotal time \
                requiedmentor was updated (%5.3f ns) is %s.' %
                  (time[-1] / 1000, str(wl_incrementor[-1])))
            print('The Wang-Landau scale and the cutoff of Wang-Landau \
                incrementor are %s and %s, respectively.' %
                  (str(self.wl_scale), str(self.wl_cutoff)))
            print('Therefore, it requires %s more updates in Wang-Landau \
                incrementor for the weights to be equilibrated.' % N_update)
            print('Check the log file for more information.')
            sys.exit()

        # =================== 4. Uncertainty estimation ======================
        kb = 1.38064852E-23                               # Boltzmann constant
        Na = 6.0221409E23                                 # Avogadro's number
        err_kt = np.abs(np.log(equil_counts[0] / equil_counts[-1]))
        err_kcal = err_kt * (kb * Na * float(self.temp) / 1000) * 0.23900573613
        print('The uncertainty of the free energy difference is %5.3f kT.\n'
              % err_kt)
        print('Or at the simulation temperature (%s K), the uncertainty is \
              %5.3f kcal/mol\n' % (str(float(self.temp)), err_kcal))

        # ================ 5. Average weights calculation ====================
        n_points = np.floor(frac * len(weights_all))
        # number of lists to be averaged

        weights_all = np.array(weights_all[-int(n_points):])
        # only average a last certain portion

        list_sum = 0
        for i in range(len(weights_all)):
            list_sum += weights_all[i]
        weights_avg_float = list_sum/len(weights_all)
        weights_avg_float = [round(x, 5) for x in weights_avg_float]

        weights_avg = ''    # make weights_ave as a string to be easily copied
        for i in range(len(weights_avg_float)):
            weights_avg += (' ' + str(weights_avg_float[i]))

        return time, wl_incrementor, final_weights, weights_avg, equil_counts

    def get_final_counts(self, logfile):
        """
        This function analyzes the log file and output the count of each
        lambda state at the last time frame (for plotting histogram).

        Paraneters
        ----------
        logfile : str
            The filename of the log file

        Returns
        -------
        counts : np.array
            The counts of each lambda state
        """
        f = open(logfile, 'r')
        lines = f.readlines()
        f.close()
        lines.reverse()     # from this point, lines has been reverse

        line_n = 0
        final_counts = np.zeros(self.N_states)
        for l in lines:
            line_n += 1
            if 'MC-lambda information' in l:
                for i in range(self.N_states):
                    # start from lines[line_n - 3]
                    final_counts[i] = float(lines[line_n - 3 - i].split()[5])
                break

        return final_counts


def main():
    time_needed = []
    s0 = timer.time()
    args = initialize()

    rc('font', **{
        'family': 'sans-serif',
        'sans-serif': ['DejaVu Sans'],
        'size': 10
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='serif')

    if isinstance(args.log, str):                # the case of only one input
        args.log = list(args.log)
    if isinstance(args.log, str) and '*' in args.log:  # to enables wildcards
        args.log = natsort.natsorted(glob.glob(args.log), reverse=False)
    if isinstance(args.keyword, str):
        args.keyword = list(args.keyword)

    # Check if the default should be used
    if args.keyword is None:
        args.keyword = [logname.split('.')[0] for logname in args.log]
    if args.frac is None:
        args.frac = 0.25
    elif float(args.frac) <= 0 or float(args.frac) >= 1:
        print('Error: The fraction for average weights calculation should be \
            between 0 and 1!')
        sys.exit()

    print('\n')
    for i in range(len(args.log)):
        result_str = 'Data analysis of the file %s:' % args.log[i]
        print(result_str)
        print('=' * len(result_str))

        log_info = LogInfo(args.log[i])
        EXE = EXEAnalysis(args.log[i])

        if log_info.fixed is False:
            [time, wl_incrementor, weights_f, weights_a, equil_counts] = \
                EXE.get_equil_info(args.log[i], args.frac)
            # If the weights are not equilibrated, the code temrinates here.

            time = time / 1000  # convert from ps to ns

            # Print the results!
            print('The weights were equilibrated at %5.3f ns\n' %
                  EXE.equil_time)
            print('The average weights of the last %s percent' %
                  str(args.frac * 100), 'of steps right before the weights \
                  are equilibrated are: ', weights_a, '\n')
            print('The final weights are: ', weights_f)
            e0 = timer.time()
            time_needed.append(e0 - s0)

            # Plot WL incrementor as a function of time
            s1 = timer.time()
            plt.figure()  # ready to plot!
            plt.step(time, wl_incrementor)
            plt.xlabel('Time (ns)')
            plt.ylabel('Wang-Landau incrementor ($ k_{B} T$)')
            plt.title('Wang-Landau incrementor as a function of time')
            plt.grid()
            plt.savefig('WL_t_%s.png' % args.keyword[i], dpi=600)
            e1 = timer.time()
            time_needed.append(e1 - s1)
            plt.show()

            # Plot the equilbrated histogram
            s2 = timer.time()
            plt.figure()
            plt.bar(np.arange(1, log_info.N_states + 1), height=equil_counts)
            plt.xlabel('States')
            plt.ylabel('Counts')
            plt.title('The equilibrated histogram of the simulation')
            plt.grid()
            plt.savefig('Equil_hist_%s.png' % args.keyword[i], dpi=600)
            e2 = timer.time()
            time_needed.append(e2 - s2)
            plt.show()

        s3 = timer.time()
        # Extract the final counts for plotting the final histogram, no matter
        # the weights are fixed during the simulation
        final_counts = EXE.get_final_counts(args.log[i])

        # Plot the final histogram
        plt.figure()
        plt.bar(np.arange(1, log_info.N_states + 1), height=final_counts)
        plt.xlabel('States')
        plt.ylabel('Counts')
        plt.title('The final histogram of the simulation')
        plt.grid()
        plt.savefig('Final_hist_%s.png' % args.keyword[i], dpi=600)
        e3 = timer.time()
        time_needed.append(e3 - s3)
        plt.show()

        print('Elasped time (including plotting): %s seconds.'
              % sum(time_needed))
