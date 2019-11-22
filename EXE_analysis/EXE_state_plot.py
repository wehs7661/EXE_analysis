import os
import sys
import glob
import argparse
import natsort
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from .EXE_histogram import LogInfo
from .EXE_histogram import EXEAnalysis


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

    args_parse = parser.parse_args()

    return args_parse


class StateTimeAnalysis(EXEAnalysis):
    """
    A class for state-time analysis of expanded ensemble simulations. When instantiating
    this class, one parameter (logfile) is required.
    """

    def __init__(self, logfile):
        self.sample_all = False
        EXEAnalysis.__init__(self, logfile)

    def get_state_time(self, dhdl, freq):
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
                    # All the states had been sampled at least once
                    # Note: all returns True if all of the items are True
                    self.sample_all = True
                    visited = [False for i in range(self.N_states)]
                    visit_time.append(float(l.split()[0]) - visit_start)
                    visit_start = float(l.split()[0])

                if i % freq == 0:
                    time.append(float(l.split()[0]))
                    state.append(int(float(l.split()[1])))
        state = np.array(state)
        time = np.array(time) / 1000      # units: ns

        return time, state, visit_time

    def plot_data(self, x, y, type, title, png_name, trnctd=False):
        # Common tasks/settings
        plt.figure()
        if type == 'state-time':
            x -= x[0]   # always shift to 0 (only has influence on truncated data)
        plt.plot(x, y)
        plt.title(title)
        plt.minorticks_on
        plt.grid(True)

        if type == 'state-time':
            if self.equil is True and self.fixed is False and trnctd is False:
                plt.axvline(x=self.equil_time, color='k', linestyle='--', zorder=10, linewidth=1.2)
            plt.xlabel('Time (ns)')
            plt.ylabel('State')
            plt.ylim([0, self.N_states])
            plt.savefig(png_name)
            plt.show()
        elif type == 'visit-time':
            plt.text(0, max(y) * 0.98, '(Number of times that all the states were sampled: %s)' %
                     len(y), fontsize=9)
            plt.xlabel('Simulation time (ns)')
            plt.ylabel('Time required to sample all the states (ns)')
            plt.savefig(png_name)
            plt.show()

        elif type == 'visit-wl':
            plt.scatter(x, y, color='k')
            plt.xscale('log')
            plt.xlabel('Wang-Landau incrementor ($ k_{B} T$)')
            plt.ylabel('Time required to sample all the states (ns)')
            plt.savefig(png_name)
            plt.show()

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
                    if time > self.equil_time * 1000:
                        outfile.write(line)


def main():
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

    for i in range(len(args.dhdl)):
        result_str = '\nData analysis of the file %s and %s:' % (args.dhdl[i], args.log[i])
        print(result_str)
        print('=' * len(result_str))

        STA = StateTimeAnalysis(args.log[i])
        # Note that STA inherits from EXE, which inherits from LogInfo
        # By initiating STA, we can use all the attributes/functions in EXE and LogInfo

        if STA.fixed is True:     # Case 1: fixed-weight simulation (EXE.get_equil_info not needed)
            # print out relevant info
            print('This is a fixed-weight simulation.\n')

            # state-time plot
            title = 'Exploration of states as a function of time with fixed weights'
            png_name = 'state_plot_%s_whole.png' % args.keyword[i]
            time, state, visit = STA.get_state_time(args.dhdl[i], args.freq)
            print('The length of the whole simulation: %6.3f ns.\n' % time[-1])
            STA.plot_data(time, state, 'state-time', title, png_name)

            # visit-time plot
            if not visit:
                print('The system could not sample all the states during %s ns of simulation.\n' % time[-1])
            else:
                # Figure 1: visit-time v.s. simulation time
                title = None
                png_name = 'visit_time_%s.png' % args.keyword[i]
                simulation_time = np.array([sum(visit[:i + 1]) for i in range(len(visit))]) / 1000   # units: ns
                STA.plot_data(simulation_time, visit, 'visit-time', title, png_name)
                # No plot of visit-time v.s. wl-incrementor for fixed weight simulation

                max_visit = max(visit)
                max_visit_idx = visit.index(max_visit)
                print('In this simulation, the longest time required for the '
                      'system to sample all the intermediate states is %s ns, which is from %s to %s ns.\n'
                      % (str(max_visit / 1000), str(sum(visit[:max_visit_idx]) / 1000), str((sum(visit[:max_visit_idx]) + max_visit) / 1000)))

        else:
            # Then STA.get_equil_info is needed
            STA = StateTimeAnalysis(args.log[i])   # inherit attributes from EXE
            wl_time, wl_incrementor, _, _ = STA.get_equil_info(args.log[i])   # update attributes like equil_time
            # Since STA inherits from EXE, instead of using STA.get_equil_info, we can use STA.get_equil_info

            # print out relevant info
            print('This is a non-fixed-weight simulation.\n')   # EXE.get_equil_info needed
            title_w = 'Exploration of states as a function of time'  # for the "whole simulation"
            title_t = 'Exploration of states as a function of time with fixed weights'   # for post-equilibration
            png_name_w = 'state_plot_%s_whole.png' % args.keyword[i]
            png_name_t = 'state_plot_%s_truncated.png' % args.keyword[i]

            if STA.equil is True:   # Case 2-1: equilibrating weight simulation with wegiths equilibrated at some point
                print('The weights were equilibrated at %5.3f ns (shown as the dash line in the figure).\n' % (STA.equil_time))

                # First deal with the whole simulation
                time, state, visit = STA.get_state_time(args.dhdl[i], args.freq)
                STA.plot_data(time, state, 'state-time', title_w, png_name_w)

                # Then, deal with the truncated file (post-equilibration)
                STA.truncate_dhdl(args.dhdl[i])
                dhdl_trnctd = args.dhdl[i].split('.xvg')[0] + '_truncated.xvg'
                time_t, state_t, visit_t = STA.get_state_time(dhdl_trnctd, args.freq)
                STA.plot_data(time_t, state_t, 'state-time', title_t, png_name_t, True)

                # visit-time plot (both the whole simulation and the post-equilibration)
                data = [visit, visit_t]
                prefix = ['Before', 'After']
                suffix = ['', '_truncated']

                str1 = 'The system could not sample all the states during %s ns of simulation.\n' % time[-1]
                fix_length = round((time[-1] - STA.equil_time), 3)
                str2 = 'After the weights were equilibrated, the system could not sample all the states in the rest of the simulation, which was from %s to %s ns (length: %s ns).\n' % (
                    str(round(STA.equil_time, 3)), str(time[-1]), str(fix_length))
                msg_not_visit_all = [str1, str2]

                for m in range(2):
                    if not data[m]:
                        print(msg_not_visit_all[m])
                    else:
                        max_visit = max(data[m])
                        max_visit_idx = data[m].index(max_visit)
                        msg_visit_all = '%s the weights were equilibrated, the longest time required for the system to sample all the intermediate states is %s ns, which is from %s to %s ns.\n' % (
                            prefix[i], str(max_visit / 1000), str(sum(data[m][:max_visit_idx]) / 1000), str((sum(data[m][:max_visit_idx]) + max_visit) / 1000))

                        print(msg_visit_all)
                        # Figure 1: visit time v.s. simulation time
                        title = None
                        png_name = 'visit_time_%s%s.png' % (args.keyword[i], suffix[m])
                        simulation_time = np.array([sum(data[m][:i + 1])
                                                    for i in range(len(data[m]))]) / 1000   # units: ns
                        STA.plot_data(simulation_time, np.array(data[m]) / 1000, 'visit-time', title, png_name)

                        # Figure 2: visit time v.s. wl-incrementor
                        visit_wl = np.zeros(len(data[m]))
                        png_name_wl = 'visit_time_wl_%s%s.png' % (args.keyword[i], suffix[m])
                        for j in range(len(data[i])):
                            for k in range(len(wl_time)):
                                if simulation_time[j] * 1000 > wl_time[k] and simulation_time[j] * 1000 < wl_time[k + 1]:
                                    visit_wl[j] = wl_incrementor[k]
                        STA.plot_data(visit_wl, np.array(data[i]) / 1000, 'visit-wl', title, png_name_wl)

            else:
                print('The weights have not been equilibrated.\n')

                # state-time plot
                title = 'Exploration of states as a function of time'
                png_name = 'state_plot_%s_whole.png' % args.keyword[i]
                time, state, visit = STA.get_state_time(args.dhdl[i], args.freq)
                print('The length of the whole simulation: %6.3f ns.\n' % time[-1])
                STA.plot_data(time, state, 'state-time', title, png_name)

                # visit-time plot
                if not visit:
                    print('The system could not sample all the states during %s ns of simulation.\n' % time[-1])
                else:
                    # Figure 1: visit-time v.s. simulation time
                    title = None
                    png_name = 'visit_time_%s.png' % args.keyword[i]
                    simulation_time = np.array([sum(visit[:i + 1]) for i in range(len(visit))]) / 1000   # units: ns
                    STA.plot_data(simulation_time, visit, 'visit-time', title, png_name)
                    # No plot of visit-time v.s. wl-incrementor for fixed weight simulation

                    max_visit = max(visit)
                    max_visit_idx = visit.index(max_visit)
                    print('In this simulation, the longest time required for the '
                          'system to sample all the intermediate states is %s ns, which is from %s to %s ns.\n'
                          % (str(max_visit / 1000), str(sum(visit[:max_visit_idx]) / 1000), str((sum(visit[:max_visit_idx]) + max_visit) / 1000)))

                    # Figure 2: visit time v.s. wl-incrementor
                    visit_wl = np.zeros(len(visit))
                    png_name_wl = 'visit_time_wl_%s.png' % (args.keyword[i])
                    for j in range(len(visit)):
                        for k in range(len(wl_time)):
                            if simulation_time[j] * 1000 > wl_time[k] and simulation_time[j] * 1000 < wl_time[k + 1]:
                                visit_wl[j] = wl_incrementor[k]
                    STA.plot_data(visit_wl, np.array(data[i]) / 1000, 'visit-wl', title, png_name_wl)

        print('%s file(s) (%s simulation(s)) analyzed.' % (len(args.log) * 2, len(args.log)))

