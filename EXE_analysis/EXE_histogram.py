"""
This is a Python script for analyzing the log file generated from the expanded
ensemble simulation. The script performs the following analysis:
1. Generate a plot of Wang-Landau incrementor as a function of time.
2. Show and save the histogram at the last time frame of the simulation.
3. Show and save the histogram at the time that the weights are equilibrated.
4. Print out equilibrated Wang-Landau weights which can directly be pasted to
    the .mdp file if needed.
5. Print out the averaged Wang-Landau weights over a user-definied length of 
    simulation.
6. Estimate the uncertainty in free energy difference between the first and
    the last state from the final histogram.
"""

import os
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
                        type=str,
                        nargs='+',
                        help='The filename(s) of the log file. Wildcards \
                            available.')
    parser.add_argument('-k',
                        '--keyword',
                        type=str,
                        nargs='+',
                        help='The keyword used in the filename of the figure \
                            produced.')
    parser.add_argument('-a',
                        '--avg_len',
                        type=float,
                        default=20,
                        help='The length of the simulation that the calculation \
                            of average weights are based on. -a 20 means that the \
                            weights of last 20 ns before the weights are eauilibrated\
                            will be averaged. Default: 20 ns.')

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
            self.init_wl, self.temp, self.start, self.nstlog = None, None, \
            None, None, None, None, None, None, None
        line_n = 0

        for l in lines:
            line_n += 1
            if 'dt  ' in l and self.dt is None:
                self.dt = float(l.split('=')[1])

            if 'nstlog' in l and self.nstlog is None:
                self.nstlog = float(l.split('=')[1])

            if 'weight-equil-wl-delta' in l and self.cutoff is None:
                self.cutoff = float(l.split('=')[1])

            if 'wl-scale' in l and self.wl_scale is None:
                self.wl_scale = float(l.split('=')[1])

            if 'n-lambdas' in l and self.N_states is None:
                self.N_states = int(l.split('=')[1])

            if 'lmc-stats' in l and self.fixed is None:
                if l.split('=')[1].split()[0] == 'no':
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
        self.avg_start = None  # the starting piont of the weights average calculation
        self.avg_end = None   # the endpoint time frame of weights average calculation
        LogInfo.__init__(self, logfile)
        # So any instance of EXEAnalysis will have all the properties that a
        # LogInfo instance has

    def get_equil_info(self, logfile):
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

        Returns
        -------
        time : np.array
            An array containing the time frames when WL incrementor changed
        wl_incrementor : np.array
            The array of Wang-Landau incrementors
        final_weights : str
            The final equilibrated lambda weights
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

            # Collect revelant information about weights equilibration
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

                wl_incrementor = np.array(wl_incrementor)
                time = np.array(step) * self.dt  # time array for plotting
                self.equil_time = float(equil_step) * self.dt / 1000  # units: ns

                # Note: avg_end will not be taken into account in the average calculation
                avg_endstep = int(equil_step) - (int(equil_step) % self.nstlog) + self.nstlog
                self.avg_end = avg_endstep * self.dt   # units: ps
                
                break

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
        print('Or at the simulation temperature (%s K), the uncertainty is %5.3f kcal/mol\n' % (str(float(self.temp)), err_kcal))

        return time, wl_incrementor, final_weights, equil_counts

    def get_avg_weights(self, logfile, avg_len):
        """
        This function performs the weights average calculation.

        Parameters
        ----------
        logfile : str
            The filename of the log file

        Returns
        -------
        weights_avg : np.array
            The average weights over an user-defined length of the simulation.
        """
        self.avg_start = self.avg_end - avg_len * 1000     # unit: ps
        if self.avg_start <= 0:
            print('The starting point of the weights average calculation is less than 0!')
            print('Please specified a shorter length of simulation for the',
                  'weights average calculation.')
            sys.exit()

        f = open(logfile, 'r')
        lines = f.readlines()
        f.close()

        search_start = False
        line_n = self.start
        weights = []    # to store weights at each time frame to be averaged
        weights_all = []    # a list of list of weights at different time frames

        # collect the data of weights to be averaged
        for l in lines[self.start:]:    # skip the metadata
            line_n += 1
            if str(self.avg_start) in l: 
                search_start = True

            if 'Wang-Landau incrementor is:' in l and search_start is True:
                for i in range(self.N_states):
                    weights.append(float(lines[line_n + i + 1].split()[6]))
                weights_all.append(weights)
                weights = []

            if str(self.avg_end) in l:
                break
        
        # Average the weights
        list_sum = 0
        weights_all = np.array(weights_all)
        for i in range(len(weights_all)):
            list_sum += weights_all[i]
        weights_avg_float = list_sum/len(weights_all)
        weights_avg_float = [round(x, 5) for x in weights_avg_float]

        weights_avg = ''   # make weights_avg as a string to be easily copied
        for i in range(len(weights_avg_float)):
            weights_avg += (' ' + str(weights_avg_float[i]))

        return weights_avg

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
        final_time : float
            The last time frame of the simulation outputting the histogram information.
        final_counts : np.array
            The final counts of each lambda state.
        """
        f = open(logfile, 'r')
        lines = f.readlines()
        f.close()
        lines.reverse()     # from this point, lines has been reverse

        final_found = False
        line_n = 0
        final_counts = np.zeros(self.N_states)
        for l in lines:
            #  print(l)   # this will print from the bottom
            line_n += 1
            if 'MC-lambda information' in l:  # should find this line first
                final_found = True
                for i in range(self.N_states):
                    # start from lines[line_n - 3]
                    # 'MC-lambda information' is lines[line_n - 1]
                    final_counts[i] = float(lines[line_n - 3 - i].split()[5])

            if '  Step  ' in l and final_found is True:
                # '    Step      Time    ' is lines[line_n - 1]
                final_time = float(lines[line_n - 2].split()[1])  # in ps
                break

        return final_time, final_counts


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

    if args.log is None:
        args.log = []
        for file in os.listdir('.'):
            if file.endswith('.log'):
                args.log.append(file)
        if not args.log:
            print('No log files found! Please check if the directory is correct or specify the name of the log files.')
        else:
            args.log = natsort.natsorted(args.log)

    if isinstance(args.log, str):                # the case of only one input
        args.log = list(args.log)
    if isinstance(args.log, str) and '*' in args.log:  # to enables wildcards
        args.log = natsort.natsorted(glob.glob(args.log), reverse=False)
    if isinstance(args.keyword, str):
        args.keyword = list(args.keyword)

    # Check if the default should be used
    if args.keyword is None:
        args.keyword = [logname.split('.')[0] for logname in args.log]

    log_str = ''
    for i in range(len(args.log)):
        if i == 0:
            log_str += args.log[i]         
        elif i != len(args.log) - 1:
            log_str += ', '
            log_str += args.log[i]
        else:
            log_str += ', and '
            log_str += args.log[i]

    print('\n')
    print('The log files to be analyzed: %s.' % log_str)
    print('Length of the simulation for the weights average calculation: %s ns.\n' % args.avg_len)
    e0 = timer.time()
    time_needed.append(e0 - s0)

    for i in range(len(args.log)):
        s1 = timer.time()
        result_str = 'Data analysis of the file %s:' % args.log[i]
        print(result_str)
        print('=' * len(result_str))

        log_info = LogInfo(args.log[i])
        EXE = EXEAnalysis(args.log[i])
        e1 = timer.time()
        time_needed.append(e1 - s1)

        if log_info.fixed is False:
            s2 = timer.time()
            [time, wl_incrementor, weights_f, equil_counts] = EXE.get_equil_info(args.log[i])
            # If the weights are not equilibrated, the code temrinates here.

            time = time / 1000  # convert from ps to ns

            # Print the results!
            print('The weights were equilibrated at %5.3f ns\n' %
                  EXE.equil_time)
            print('The final weights are:\n', weights_f)

            e2 = timer.time()
            time_needed.append(e2 - s2)

            # Plot WL incrementor as a function of time
            s3 = timer.time()
            plt.figure()  # ready to plot!
            plt.step(time, wl_incrementor)
            plt.xlabel('Time (ns)')
            plt.ylabel('Wang-Landau incrementor ($ k_{B} T$)')
            plt.title('Wang-Landau incrementor as a function of time')
            plt.grid()
            plt.savefig('WL_t_%s.png' % args.keyword[i], dpi=600)
            e3 = timer.time()
            time_needed.append(e3 - s3)
            plt.show()

            # Plot the equilbrated histogram
            s4 = timer.time()
            etime_title = str(round(EXE.equil_time, 1))     # first decimal point
            etime_png = str(int(round(EXE.equil_time, 0)))  # round to integrer
            # the keyword corresponding to the equilibrated time in the filename of png
            plt.figure()
            plt.bar(np.arange(1, log_info.N_states + 1), height=equil_counts)
            plt.xlabel('States')
            plt.ylabel('Counts')
            plt.title('The equilibrated histogram of the simulation (%s ns)' % etime_title)
            if max(equil_counts) >= 10000:
                plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            plt.grid()
            plt.savefig('Equil_hist_%sns_%s.png' % (etime_png, args.keyword[i]), dpi=600)
            e4 = timer.time()
            time_needed.append(e4 - s4)
            plt.show()

        s5 = timer.time()
        if log_info.fixed is True:
            print('This is a fixed-weight expanded ensemble simulation.')
            print('Accordingly, only the final histogram will be output and saved.')
        # Extract the final counts for plotting the final histogram, no matter
        # the weights are fixed during the simulation
        final_time, final_counts = EXE.get_final_counts(args.log[i])
        final_time /= 1000    # from ps to ns
        ftime_title = str(round(final_time, 1))  # 1st decimal point
        ftime_png = str(int(round(final_time, 0)))    # round to integer
        # the keyword corresponding to the final time in the filename of png

        # Plot the final histogram
        plt.figure()
        plt.bar(np.arange(1, log_info.N_states + 1), height=final_counts)
        plt.xlabel('States')
        plt.ylabel('Counts')
        plt.title('The final histogram of the simulation (at %s ns)' % ftime_title)
        if max(final_counts) >= 10000:
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.grid()
        plt.savefig('Final_hist_%sns_%s.png' % (ftime_png, args.keyword[i]), dpi=600)
        e5 = timer.time()
        time_needed.append(e5 - s5)
        plt.show()
        
        s6 = timer.time()
        # Average weights calculation
        weights_a = EXE.get_avg_weights(args.log[i], args.avg_len)
        print('The average weights of the last %s ns' % str(args.avg_len), 
        'of the simulation before the weights are equilibrated (from %s to %s ns) are:\n' 
        %(str(EXE.avg_start / 1000), str(EXE.avg_end / 1000) ), weights_a, '\n')
        e6 = timer.time()
        time_needed.append(e6 - s6)
        
    print('%s files analyzed.' % len(args.log))
    print('Total time elapsed (including plotting): %s seconds.\n'
              % sum(time_needed))
