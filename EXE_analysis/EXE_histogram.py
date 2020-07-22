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
                            will be averaged. Default: 20 ns. If multiple files are\
                            given, the same value applies to all.')
    parser.add_argument('-m',
                        '--mdp',
                        type=str,
                        help='The .mdp file as the basis of the newly generated .mdp\
                            file for the fixed-weight simulation. Note that a new .mdp\
                            file will be generated only if there were only one log\
                            file provided and the weights had been equilibrated.')
    parser.add_argument('-w',
                        '--weights',
                        type=str,
                        choices=['equilibrated', 'adjusted', 'average'],
                        default='adjusted',
                        help='Which kind of weights to be used in the new .mdp file.')

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
            self.init_wl, self.temp, self.start, self.nstlog, self.wl_ratio = None, None, \
            None, None, None, None, None, None, None, None
        self.init_w = []
        line_n = 0

        for l in lines:
            line_n += 1
            if 'dt  ' in l and self.dt is None:
                self.dt = float(l.split('=')[1])

            if 'nstlog' in l and self.nstlog is None:
                self.nstlog = float(l.split('=')[1])

            if 'weight-equil-wl-delta' in l and self.cutoff is None:
                self.cutoff = float(l.split('=')[1])

            if 'wl-ratio' in l and self.wl_ratio is None:
                self.wl_ratio = float(l.split('=')[1])

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

            if 'init-lambda-weights[' in l:
                self.init_w.append(float(l.split('=')[1]))

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
        self.max_Nratio = None   # To examine the flatness of the histogram
        self.min_Nratio = None   # To examine the flatness of the histogram 
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
            if 'Wang-Landau incrementor is' in l and not wl_incrementor:
                # not wl_incrementor returns True if the list if
                wl_incrementor.append(self.init_wl)  # initial wl-incrementor

            if 'weights are now: ' in l:
                n += 1
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
                            if lines[search_n + 2 + i].split()[-1] == '<<':
                                equil_counts.append(float(lines[search_n + 2 + i].split()[-4]))
                            else:
                                equil_counts.append(float(lines[search_n + 2 + i].split()[-3]))
                
                self.equil_time = float(equil_step) * self.dt / 1000  # units: ns

                # Note: avg_end will not be taken into account in the average calculation
                avg_endstep = int(equil_step) - (int(equil_step) % self.nstlog) + self.nstlog
                self.avg_end = avg_endstep * self.dt   # units: ps

                break

        time = np.array(step) * self.dt  # time array for plotting
        wl_incrementor = np.array(wl_incrementor)
        # ========== 3. Exit if the weights have not equilibrated ============
        if self.fixed is False and self.equil is False:
            if len(wl_incrementor) != 0 and self.wl_scale < 0.999:  # avoid speical cases that wl_scale > 0.999
                N_updated = len(wl_incrementor) -1 
                N_update = int(np.ceil(np.log(self.cutoff / wl_incrementor[-1]) /
                                np.log(self.wl_scale)))  # number of updates required
            _, final_counts = self.get_final_counts(logfile)  # to get self.final_w
            diff_w = np.array(self.final_w) - np.array(self.init_w)
            diff_w = [round(x, 2) for x in diff_w]
            print('The weights have not equilibrated.')
            print('Initial weights: %s' % (' '.join([str(i) for i in self.init_w])))
            print('Final weights:   %s' % (' '.join([str(i) for i in self.final_w])))
            print('The difference between the initial weights and final weights are:')
            print(' '.join(str(i) for i in diff_w))
            #print('\nThe Wan-Landau incrementor has been updated for %s times.' % N_updated)
            #print('The last time frame that the Wang-Landau incrementor was updated (%5.3f ns) is %s.' %
                  #(time[-1] / 1000, str(wl_incrementor[-1])))
            #print('The Wang-Landau scale and the cutoff of Wang-Landau incrementor are %s and %s, respectively.' %
                  #(str(self.wl_scale), str(self.cutoff)))
            #print('Therefore, it requires %s more updates in Wang-Landau incrementor for the weights to be equilibrated.' % N_update)
            print('Check the log file for more information.')


            
            #weights_list = [float(self.final_w.split()[i]) for i in range(len(self.final_w.split()))]
            weights_list = self.final_w
            weights_adjstd = np.zeros(len(weights_list))
            wght_adjstd_str = ''
            for i in range(len(weights_list)):
                if i == 0:
                    weights_adjstd[i] = 0
                    wght_adjstd_str += '0.00000 ' 
                else:
                    weights_adjstd[i] = weights_list[i] + np.log(final_counts[i - 1] / final_counts[i])
                    wght_adjstd_str += '%6.5f ' % weights_adjstd[i] 
            RMSD = np.sqrt((1/len(weights_list) * sum((np.array(weights_list) - weights_adjstd) ** 2)))
            print(wght_adjstd_str)
            


            return time, wl_incrementor
            #sys.exit()

        # ================== 4. Additional information =======================
        if self.equil is True:
            avg_counts = sum(equil_counts) / len(equil_counts)
            self.max_Nratio = max(equil_counts) / avg_counts
            self.min_Nratio = min(equil_counts) / avg_counts

        # ==================== 5. Weights adjustment =========================
        # Formula: g'_k = g_k + ln(count_(k - 1) / count_k)
        if self.equil is True:    
            weights_list = [float(final_weights.split()[i]) for i in range(len(final_weights.split()))]
            weights_adjstd = np.zeros(len(weights_list))
            wght_adjstd_str = ''
            for i in range(len(weights_list)):
                if i == 0:
                    weights_adjstd[i] = 0
                    wght_adjstd_str += '0.00000 ' 
                else:
                    weights_adjstd[i] = weights_list[i] + np.log(equil_counts[i - 1] / equil_counts[i])
                    wght_adjstd_str += '%6.5f ' % weights_adjstd[i] 
            RMSD = np.sqrt((1/len(weights_list) * sum((np.array(weights_list) - weights_adjstd) ** 2)))
            print(wght_adjstd_str)
        # =================== 6. Uncertainty estimation ======================
        kb = 1.38064852E-23                               # Boltzmann constant
        Na = 6.0221409E23                                 # Avogadro's number
        if self.equil is True:
            err_kt = np.abs(np.log(equil_counts[0] / equil_counts[-1]))
            err_kcal = err_kt * (kb * Na * float(self.temp) / 1000) * 0.23900573613
            print('The uncertainty of the free energy difference is %5.3f kT.\n'
                % err_kt)
            print('Or at the simulation temperature (%s K), the uncertainty is %5.3f kcal/mol\n' %
                (str(float(self.temp)), err_kcal))

        return time, wl_incrementor, final_weights, equil_counts, RMSD, wght_adjstd_str

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
        self.avg_end=5000
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
            search_start=True
            #if 'Wang-Landau incrementor is:' in l and search_start is True:
            if 'MC-lambda information' in l and search_start is True:
                for i in range(self.N_states):
                    if lines[line_n + i + 1].split()[-1] == '<<':
                        weights.append(float(lines[line_n + i + 1].split()[-3]))
                    else:
                        weights.append(float(lines[line_n + i + 1].split()[-2]))
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
        print(weights_avg)
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
        self.final_w = []
        for l in lines:
            #  print(l)   # this will print from the bottom
            line_n += 1
            if 'MC-lambda information' in l:  # should find this line first
                final_found = True
                if self.fixed is True:
                    data_line = line_n - 3
                else:
                    data_line = line_n - 4
                
                for i in range(self.N_states):
                    if lines[data_line - i].split()[-1] == '<<':
                        self.final_w.append(float(lines[data_line - i].split()[-3]))
                        final_counts[i] = float(lines[data_line - i].split()[-4])
                    else:
                        self.final_w.append(float(lines[data_line - i].split()[-2]))
                        final_counts[i] = float(lines[data_line - i].split()[-3])
                
                """
                for i in range(self.N_states):
                    if lines[data_line - i + 1].split()[-1] == '<<':
                        self.final_w.append(float(lines[data_line - i + 1].split()[-3]))
                        final_counts[i] = float(lines[data_line - i + 1].split()[-4])
                    else:
                        self.final_w.append(float(lines[data_line - i + 1].split()[-2]))
                        final_counts[i] = float(lines[data_line - i + 1].split()[-3])
                """ 
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

    gen_mdp = False # whether to generate a new mdp for the fixed-weight simulation

    if args.log is None:
        args.log = []
        for file in os.listdir('.'):
            if file.endswith('.log'):
                args.log.append(file)
        if not args.log:
            print('No log files found! Please check if the directory is correct or specify the name of the log files.')
        else:
            args.log = natsort.natsorted(args.log)
        #if len(args.log) == 1:
        #    gen_mdp = True

    if isinstance(args.log, str) and '*' in args.log:  # to enables wildcards
        args.log = natsort.natsorted(glob.glob(args.log), reverse=False)

    # Check if the default keyword should be used
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
    print('The log file(s) to be analyzed: %s.' % log_str)
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

        #output = EXE.get_equil_info(args.log[i])
        if log_info.fixed is False and EXE.equil is True:
            s2 = timer.time()
            time, wl_incrementor, weights_f, equil_counts, RMSD, wght_adjstd_str = EXE.get_equil_info(args.log[i])
            # If the weights are not equilibrated, the code temrinates here.

            time = time / 1000  # convert from ps to ns

            # Print the results!
            print('The weights were equilibrated at %5.3f ns\n' %
                  EXE.equil_time)
            print('The Wang-Landau ratio was set as %s, which means that all N_ratio should be larger than %s and smaller than %6.5f.\n'
                % (log_info.wl_ratio, log_info.wl_ratio, 1 / float(log_info.wl_ratio)))
            print('At the time that the weights were equilibrated, the largest and smallest N_ratio are %6.5f and %6.5f, respectively.\n'
                % (EXE.max_Nratio, EXE.min_Nratio))
            print('The final/equilibrated weights are:\n', weights_f)
            print('The adjusted weights based on the equilibrated histogram are (RMSD: %6.5f kT):\n %s\n' % (RMSD, wght_adjstd_str))

            e2 = timer.time()
            time_needed.append(e2 - s2)

            # Plot WL incrementor as a function of time
            s3 = timer.time()
            plt.figure()  # ready to plot!
            plt.step(time, wl_incrementor)
            plt.xlabel('Time (ns)')
            plt.ylabel('Wang-Landau incrementor ($ k_{B} T$)')
            plt.minorticks_on()
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
            plt.minorticks_on()
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
            print('Accordingly, only the final histogram will be output and saved.\n')
        elif EXE.equil is False:
            EXE.get_equil_info(args.log[i])

        if log_info.fixed is False and EXE.equil is False or log_info.fixed is True:
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
            plt.minorticks_on()
            plt.title('The final histogram of the simulation (at %s ns)' % ftime_title)
            if max(final_counts) >= 10000:
                plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            plt.grid()
            plt.savefig('Final_hist_%sns_%s.png' % (ftime_png, args.keyword[i]), dpi=600)
            e5 = timer.time()
            time_needed.append(e5 - s5)
            plt.show()

        
        # Average weights calculation
        #if log_info.fixed is False and EXE.equil is True:
        if True:
            s6 = timer.time()
            weights_a = EXE.get_avg_weights(args.log[i], args.avg_len)
            print('The average weights over the last %s ns' % str(args.avg_len),
                'of the simulation before the weights are equilibrated (from %s to %s ns) are:\n'
                % (str(EXE.avg_start / 1000), str(EXE.avg_end / 1000)), weights_a, '\n')
            e6 = timer.time()
            time_needed.append(e6 - s6)
        
        # Generate the new mdp file for the fixed-weight simulation
        # Note that this section of code is unreachable if the weights had not bee equilibrated
        if args.mdp is None and gen_mdp is True:
            for file in os.listdir('.'):
                if file.endswith('.mdp') and file != 'mdout.mdp':
                    args.mdp = file
            if args.mdp is None:
                gen_mdp = False
                print('No mdp file provided or found! Therefore, no mdp file for the fixed-weight simulation will be generated.\n')

        if gen_mdp is True:
            new_mdp = args.mdp.split('.mdp')[0] + '_fixed.mdp'
            print('Generating a new mdp file for the fixed-weight simulation from %s ...\n' % args.mdp)
            f = open(args.mdp, 'r')
            mdp_lines = f.readlines()
            f.close()
            
            # make sure to start from an empty file
            if os.path.isfile(new_mdp):
                os.system("rm %s" % new_mdp)
            else:
                os.system("touch %s" % new_mdp)
        
            for l in mdp_lines:
                with open(new_mdp, "a+") as file:
                    # turn off lmc-stats
                    if 'lmc-stats' in l:
                        space = l.split('=')[0].split('lmc-stats')[1]
                        l = 'lmc-stats' + space + '= no\n'

                    # parameters to be commented out
                    if 'lmc-weights-equil' in l \
                        or 'weight-equil-wl-delta' in l \
                        or 'init-wl-delta' in l \
                        or 'wl-scale' in l \
                        or 'wl-ratio' in l:
                        l = ';' + l

                    # initial guess of the weights
                    if 'init-lambda-weights' in l:
                        space = l.split('=')[0].split('init-lambda-weights')[1]
                        l = 'init-lambda-weights' + space + '= '
                        if args.weights == 'equilibrated':
                            l += weights_f + '\n'
                        if args.weights == 'adjusted':
                            l += wght_adjstd_str + '\n'
                        if args.weights == 'average':
                            l += weights_a + '\n'
                    file.write(l)
            print('Type of weights to be used: %s\n' % args.weights)
            print('A new mdp file for the fixed-weight simulation (%s) is generated!\n' % new_mdp)

    print('%s file(s) analyzed.' % len(args.log))
    print('Total time elapsed (including plotting): %s seconds.\n'
          % sum(time_needed))
