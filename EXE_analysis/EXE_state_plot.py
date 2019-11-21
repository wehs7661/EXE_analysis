import os
import sys
import glob
import time as timer
import argparse
import natsort
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
                        type=str,
                        nargs='+',
                        help='The file name of the dhdl file.')
    parser.add_argument('-l',
                        '--log',
                        type=str,
                        nargs='+',
                        help='The file name of the log file.')
    parser.add_argument('-f',
                        '--freq',
                        type=int,
                        default=1,
                        help='The frequency (every certain steps) to extract \
                              the data to plot. Default: 1 (step). If multiple\
                              files are given, the same value applies to all.')
    parser.add_argument('-k',
                        '--keyword',
                        type=str,
                        nargs='+',
                        help='The keyword used in the filename of the figure \
                             produced.')
    parser.add_argument('-s',
                        '--shift',
                        default=True,
                        action='store_false',
                        help='Whether to shift the first time frame back to \
                             zero if the dhdl file is truncated (-s means shift \
                             the data). Will make no differences if no log file \
                             is provided. If multiple files are given, the same \
                             value of this parameter applies to all.')

    args_parse = parser.parse_args()

    return args_parse


class StateTimeAnalysis:
    """
    A class for state-time analysis of expanded ensemble simulations
    """

    def __init__(self):
        self.N_states = None
        self.fixed = None
        self.equil = False
        self.equil_time = None
        self.dt = None
        self.init_wl = None

    def get_equil_info(self, logfile):
        f = open(logfile, 'r')
        lines = f.readlines()
        f.close()

        step = []
        wl_incrementor = []

        # First find the time frame that the weights were equilibrated
        line_n = 0
        for l in lines:
            line_n += 1
            if ' dt ' in l and self.dt is None:
                self.dt = float(l.split('=')[1])

            if 'n-lambdas' in l and self.N_states is None:
                self.N_states = int(l.split('=')[1])

            if 'lmc-stats' in l and self.fixed is None:
                if l.split('=')[1].split()[0] == 'no':
                    self.fixed = True
                    break
                else:
                    self.fixed = False

            if 'init-wl-delta' in l and self.init_wl is None:
                self.init_wl = float(l.split('=')[1])

            if 'Wang-Landau incrementor is' in l and not wl_incrementor:
                wl_incrementor.append(self.init_wl)   # initial wl-incrementor

            if 'weights are now: ' in l:
                weights_line = l.split(':')
                step.append(int(weights_line[0].split()[1]))
                search_lines = lines[line_n + 1: line_n + 15]
                for l_search in search_lines:
                    if 'Wang-Landau incrementor is:' in l_search:
                        wl_incrementor.append(float(l_search.split(':')[1]))

            if 'Weights have equilibrated' in l:
                self.equil = True
                equil_last_step = l.split(':')[0].split()[1]
                # the step that the weights are equilibrated
                self.equil_time = float(equil_last_step) * self.dt   # units: ps
                break

        wl_incrementor = np.array(wl_incrementor)
        time = np.array(step) * self.dt

        return time, wl_incrementor

    def truncate_dhdl(self, dhdl):
        """
        This function reads in a dhdl file, truncate the part that the weights are still equilibrating and save as a new file
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
                    if time > self.equil_time:
                        outfile.write(line)

    def state_time_data(self, dhdl, freq):
        """Reads a dhdl file to returns the data required for plotting a state-time plot.
        Also, estimate how long it was required for the system to visit all the intermediate
        staes, either before or after the weights were equilibrated.

        Parameters
        -----------
        dhdl : str
            The name of the dhdl file to read.
        freq : float
            The frequency to extract the state-time data

        Returns
        ----------
        state : np.array
            The lambda states.
        time : np.array
            The times. 
        """
        f = open(dhdl, 'r')
        lines = f.readlines()
        f.close()

        # If the number of state is unknown (b.c. that equil_time wasn't used)
        # it will be figured out from the metatexts of the dhdl file
        if self.N_states is None:
            n = 0   # line number of the metatexts
            self.N_states = 0
            for l in lines:
                n += 1   # so l is lines[n - 1]
                if l[0] == '@' and ' to ' in l:
                    self.N_states += 1
                    if lines[n][0] != '@':
                        break

        # Collect the state-time data and estimate the "visit-all time"
        visited = [False for i in range(self.N_states)]
        visit_start = 0         # starting point
        visit_time = []         # a list of visit-all times
        state, time = [], []    # units of time :ps
        i = 0                   # line number (excluding metatexts)
        for l in lines:
            if l[0] != '#' and l[0] != '@':
                i += 1
                n = int(float(l.split()[1]))  # the state just visited
                visited[n] = True
                if all(visited) is True:
                    # all returns True if all of the items are True
                    visited = [False for i in range(self.N_states)]
                    visit_time.append(float(l.split()[0]) - visit_start)
                    visit_start = float(l.split()[0])

                if i % freq == 0:
                    time.append(float(l.split()[0]))
                    state.append(int(float(l.split()[1])))
        state = np.array(state)
        time = np.array(time) / 1000      # units: ns

        return time, state, visit_time


def main():
    time_needed = []
    s0 = timer.time()
    args = initialize()

    if args.dhdl is None:
        args.dhdl = []
        for file in os.listdir('.'):
            if file.endswith('_dhdl.xvg'):
                args.dhdl.append(file)
        if not args.dhdl:
            print('No dhdl files provided or found! Please check if the directory is correct or specify the name of the dhdl file.')
            sys.exit()
        else:
            args.dhdl = natsort.natsorted(args.dhdl)

    if args.log is None:
        args.log = []
        for file in os.listdir('.'):
            if file.endswith('.log'):
                args.log.append(file)
        if not args.log:
            print('No log files provided or found! Pleaes check if the directory is correct or specify the name of the log file.')
            sys.exit()
        else:
            args.log = natsort.natsorted(args.log)

    if isinstance(args.dhdl, str) and '*' in args.dhdl:   # to enable wildcards
        args.dhdl = natsort.natsorted(glob.glob(args.dhdl), reverse=False)
    if isinstance(args.log, str) and '*' in args.log:     # to enable wildcards
        args.log = natsort.natsorted(glob.glob(args.log), reverse=False)

    if len(args.dhdl) != len(args.log):
        print('Different number of log files and dhdl files were provided.')
        print('Please check that each dhdl file has its corresponding log file.')
        sys.exit()

    # Check if the default keyword should be used
    if args.keyword is None:
        args.keyword = [dhdl.split('_dhdl.xvg')[0] for dhdl in args.dhdl]

    rc('font', **{
        'family': 'sans-serif',
        'sans-serif': ['DejaVu Sans'],
        'size': 10
    })
    # Set the font used for MathJax - more on this later
    rc('mathtext', **{'default': 'regular'})
    plt.rc('font', family='serif')

    log_str, dhdl_str = '', ''
    for i in range(len(args.log)):
        if i == 0:
            dhdl_str += args.dhdl[i]
            log_str += args.log[i]
        elif i != len(args.log) - 1:
            dhdl_str += ', '
            dhdl_str += args.dhdl[i]
            log_str += ', '
            log_str += args.log[i]
        else:
            dhdl_str += ', and '
            dhdl_str += args.dhdl[i]
            log_str += ', and '
            log_str += args.log[i]

    print('\n')
    print('The dhdl file(s) to be analyzed: %s.' % dhdl_str)
    if len(args.log) > 0:
        print('The log file(s) to be analyzed: %s.' % log_str)
    
    e0 = timer.time()
    time_needed.append(e0 - s0)

    for i in range(len(args.dhdl)):
        s1 = timer.time()
        result_str = '\nData analysis of the file %s and %s:' % (args.dhdl[i], args.log[i])
        print(result_str)
        print('=' * len(result_str))

        STA = StateTimeAnalysis()
        wl_time, wl_incrementor = STA.get_equil_info(args.log[i])

        # Three possible cases to deal with
        # Case 1: Fixed-weight simulation (result: one state-time plot)
        # Case 2: Equilibrating-weight simulation
        # Case 2-1: the weights were equilibrated at some point (truncation required, result: two state-time plots)
        # Case 2-2: the weights had not been equilibrated (result: one state-time plot)
        # All three cases output a figure of the whole (untruncated) simulation.
        # Therefore, we first deal with the whole (untruncated) simulation.

        dhdl_file = args.dhdl[i]    # dhdl file to be analyzed
        png_name = 'state_plot_%s_whole.png' % args.keyword[i]   # untruncated
        time, state, visit = STA.state_time_data(dhdl_file, args.freq)

        # print out some info and decide the title of the plot
        if STA.fixed is True:
            title = 'Exploration of states as a function of time with fixed weights.'
            print('Note: This is a fixed-weight simulation.\n')
        else:
            print('Note: This is non-fixed-weight simulation.\n')
            title = 'Exploration of states as a function of time'
            if STA.equil is True:
                print('The weights were equilibrated at %5.3f ns (shown as the dash line in the figure).\n' %
                      (STA.equil_time / 1000))
            else:
                print('The weights have not been equilibrate.\n')
        print('The length of the whole simulation: %6.3f ns.\n' % time[-1])

        # Having all the data, we can plot the figures now.
        # Start with the figure based on the whole (untruncated) simulation
        plt.plot(time, state)
        if STA.equil is True and STA.fixed is False:
            plt.axvline(x=STA.equil_time / 1000, color='k', linestyle='--', zorder=10, linewidth=1.2)
        plt.title(title)
        plt.xlabel('Time (ns)')
        plt.ylabel('State')
        plt.minorticks_on()
        plt.ylim([0, STA.N_states])
        plt.grid(True)
        plt.savefig(png_name)
        e1 = timer.time()
        time_needed.append(e1 - s1)
        plt.show()

        # Analysis of the time required to sample all the states
        s2 = timer.time()
        if not visit:
            print('The system could not sample all the states during %s ns of simulation.\n' % time[-1])
        else:
            max_visit = max(visit)
            max_visit_idx = visit.index(max_visit)
            # Note: this is for the part before the weights were equilibrated
            if STA.fixed is False and STA.equil is True:  # truncation required
                print('Before the weights were equilibrated, the longest time required for the '
                    'system to sample all the intermediate states is %s ns, which is from %s to %s ns.\n'
                    % (str(max_visit / 1000), str(sum(visit[:max_visit_idx]) / 1000), str((sum(visit[:max_visit_idx]) + max_visit) / 1000)))
            else:   # either fixed simulation or equilibrating simulation with non-equilibrated weights
                print('In this simulation, the longest time required for the '
                    'system to sample all the intermediate states is %s ns, which is from %s to %s ns.\n'
                    % (str(max_visit / 1000), str(sum(visit[:max_visit_idx]) / 1000), str((sum(visit[:max_visit_idx]) + max_visit) / 1000)))

            # plot the result!
            simulation_time = np.array([sum(visit[:i + 1]) for i in range(len(visit))]) / 1000   # units: ns
            plt.plot(simulation_time, np.array(visit) / 1000)
            plt.text(0, max(visit) / 1000 * 0.98, '(Number of times that all the states were sampled: %s)' %
                    len(visit), fontsize=9)
            plt.xlabel('Simulation time (ns)')
            plt.ylabel('Time required to sample all the states (ns)')
            plt.minorticks_on()
            plt.grid(True)
            plt.savefig('visit_time_%s.png' % args.keyword[i])
            e2 = timer.time()
            time_needed.append(e2 -s2)
            plt.show()

        # If in the equilibrating weight simulation, the weights were equilibrated
        # plot the 'visit-time' as a function of wl-incrementor
        # First use the simulation time to find the corresponding wl-incrementor
        if STA.fixed is False and STA.equil is True:
            s3 = timer.time()
            visit_wl = np.zeros(len(visit))
            for j in range(len(visit)):
                for k in range(len(wl_time)):
                    if simulation_time[j] * 1000 > wl_time[k] and simulation_time[j] * 1000 < wl_time[k + 1]:
                        visit_wl[j] = wl_incrementor[k]
            plt.plot(visit_wl, np.array(visit) / 1000)
            plt.scatter(visit_wl, np.array(visit) / 1000, color='k')
            plt.xscale('log')
            plt.xlabel('Wang-Landau incrementor ($ k_{B} T$)')
            plt.ylabel('Time required to sample all the states (ns)')
            plt.minorticks_on
            plt.grid(True)
            plt.savefig('visit_time_wl_%s.png' % args.keyword[i])
            e3 = timer.time()
            time_needed.append(e3 -s3)
            plt.show()

        if STA.fixed is False and STA.equil is True:   # then truncation is required
            s4 = timer.time()
            STA.truncate_dhdl(args.dhdl[i])
            dhdl_trnctd = dhdl_trnctd = args.dhdl[i].split('.xvg')[0] + '_truncated.xvg'
            png_trnctd = 'state_plot_%s_truncated.png' % args.keyword[i]
            time_t, state_t, visit_t = STA.state_time_data(dhdl_trnctd, args.freq)
            title_t = 'Exploration of states as a function of time with fixed weights.'
            os.remove(dhdl_trnctd)

            # figure based on the truncated data
            if args.shift:
                time_t -= time_t[0]

            plt.plot(time_t, state_t)
            plt.title(title_t)
            plt.xlabel('Time (ns)')
            plt.ylabel('State')
            plt.minorticks_on()
            plt.ylim([0, STA.N_states])
            plt.grid(True)
            plt.savefig(png_trnctd)
            e4 =timer.time()
            time_needed.append(e4 -s4)
            plt.show()            

            # Analysis of the time required to sample all the states
            # Note: this is for the part after the weights were equilibrated
            if not visit_t:    # length = 0
                s5 = timer.time()
                fix_length = round((time[-1] - STA.equil_time / 1000), 3)
                print('After the weights were equilibrated, the system could not sample all '
                      'the states in the rest of the simulation, which was from %s to %s ns (length: %s ns).\n'
                      % (str(round(STA.equil_time / 1000, 3)), str(time[-1]), str(fix_length)))
                e5 = timer.time()
                time_needed.append(e5 -s5)
            else:
                s6 = timer.time()
                max_visit_t = max(visit_t)
                max_visit_idx_t = visit_t.index(max_visit_t)
                print('After the weigths were equilibrated, the longest time required for the '
                      'system to sample all the intermediate states is %s ns, which is from %s to %s ns.'
                      % (str(max_visit_t / 1000), str(sum(visit_t[:max_visit_idx_t]) / 1000), str((sum(visit_t[:max_visit_idx_t]) + max_visit_t) / 1000)))

                # plot the result!
                simulation_time_t = np.array([sum(visit_t[:i + 1]) for i in range(len(visit_t))]) / 1000   # units: ns
                plt.plot(simulation_time_t, np.array(visit_t) / 1000)
                plt.text(0, max(visit_t) / 1000 * 0.98,
                         '(Number of times that all the states were sampled: %s)' % len(visit_t), fontsize=9)
                plt.xlabel('Simulation time (ns)')
                plt.ylabel('Time required to sample all the states (ns)')
                plt.minorticks_on()
                plt.grid(True)
                plt.savefig('')
                e6 = timer.time('visit_time_%s.png' % args.keyword[i])
                time_needed.append(e6- s6)
                plt.show()

    print('%s file(s) (%s simulation(s)) analyzed.' % (len(args.log) * 2, len(args.log)))
    print('Total time elapsed (including plotting): %s seconds.\n'
          % sum(time_needed))
