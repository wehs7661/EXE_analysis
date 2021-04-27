import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from .EXE_weights_evolution import WeightsEvolution

def initialize():
    parser = argparse.ArgumentParser(
        description='This code performs a block analysis on the free energy as a \
                    function of time obtained from parsing the log file and output \
                    the uncertainty as a function of block size.')
    parser.add_argument('-l',
                        '--log',
                        type=str,
                        help='The name of the log file.')
    parser.add_argument('-dt',
                        '--dt',
                        type=float,
                        default=2,
                        help='The spacing (units: ps) between the adjacent data points of the free energy difference time series.')
    parser.add_argument('-t',
                        '--truncation',
                        type=float,
                        default=0,
                        help='The fraction of the free energy time series to be truncated before block analysis.')
    parser.add_argument('-b',
                        '--block_info',
                        nargs='+',
                        required=True,
                        help='The minimum, maximum and spacing of the block size for plotting.')
    
    args_parse = parser.parse_args()
    
    if args_parse.log is None:
        for f in os.listdir('.'):
            if f.endswith('.log'):
                args_parse.log = f

    if len(args_parse.block_info) != 3:
        print('The block size information is in the wrong format!')
        sys.exit()

    return args_parse

def block_analysis(data, b_size):
    """
    Parameters
    ----------
    data : (array-like)
        Time series to be analyzed with block analysis.
    b_size : (int)
        Number of data points included in one block
    """
    n_blocks = int(np.ceil(len(data) / b_size))  # number of blocks
    blocks, block_data = [], []
    for i in range(n_blocks):
        blocks.append(data[b_size * i: b_size * (i + 1)])
        block_data.append(np.mean(blocks[-1])) 
    
    # use block_data to calculate sample mean and sample std
    block_avg = np.mean(block_data)
    block_std = np.std(block_data, ddof=1)

    return block_avg, block_std

def logger(*args, **kwargs):
    print(*args, **kwargs)
    with open("EXE_df_block_analysis.txt", "a") as f:
        print(file=f, *args, **kwargs)

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
    df = dG_0[-1]   # free energy difference (between the first and the last state)
    logger(f'Simulation length (before truncation): {len(df) * args.dt / 1000} ns ({len(df)} time frames)')
    
    df = df[int(args.truncation * len(df)):]  # truncate the data
    logger(f'{args.truncation * 100}% of the data ({int(args.truncation * len(df))} time frames, {int(args.truncation * len(df)) * args.dt / 1000}ns) was truncated before the block analysis')
    
    b_min = float(args.block_info[0])
    b_max = float(args.block_info[1])
    b_spacing = float(args.block_info[2])
    block_pts = list(np.arange(b_min, b_max, b_spacing))
    block_pts.append(b_max)

    avg, std = [], []  # avg and std as a function of time
    for i in block_pts:
        block_result = block_analysis(df, b_size=int(i))
        logger(f'Block size: {i} ({i * args.dt} ps), Free energy difference: {block_result[0]: .4f} +/- {block_result[1]: .4f} kT')
        avg.append(block_result[0])
        std.append(block_result[1])

    # Plot the free energy difference and the uncertainty as a function of block size
    plt.figure()
    plt.plot(np.array(block_pts) * args.dt, avg)
    plt.xlabel('Block size (ps)')
    plt.ylabel('Free energy difference ($ k_{B}T $)')
    plt.title('Free energy difference as a function of block size')
    plt.grid()
    plt.savefig('df_block_analysis.png', dpi=600)

    plt.figure()
    plt.plot(np.array(block_pts) * args.dt, std)
    plt.xlabel('Block size (ps)')
    plt.ylabel('Uncertainty ($ k_{B}T $)')
    plt.title('Uncertainty as a function of block size')
    plt.grid()
    plt.savefig('df_err_block_analysis.png', dpi=600)

