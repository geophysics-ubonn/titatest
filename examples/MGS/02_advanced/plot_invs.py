#!/usr/bin/env python
"""
TODO: Documentation
END DOCUMENTATION
"""
import os
import fnmatch
import subprocess
import crlab_py.tomodir as TD
import numpy as np
td = TD.tomodir()
# CRTomo.config is a dictionary with additional functions
import crlab_py.CRTomo as CRTomo
import glob
from multiprocessing import Pool, TimeoutError


def plot_td(directory):
    summary_dir = 'overview'
    outfile = summary_dir + '/' + os.path.basename(directory) + '.png'
    if os.path.isfile(outfile):
        return

    if not os.path.isdir(summary_dir):
        os.makedirs(summary_dir)
    if not os.path.isdir(directory):
        return
    print(directory)

    # plot
    pwd = os.getcwd()
    os.chdir(directory)
    cmd = ''.join((
        'plot_td.py ',
        ' --cminmag 20 --cmaxmag 50 ',
        ' --cminpha -30 --cmaxpha 0 ',
        ' --cminim -4 --cmaxim -3 --cminre -1.5 --cmaxre -1.2',
        ' --cmaglin ',
    ))
    subprocess.call(cmd, shell=True)
    subprocess.call('td_create_overview.py', shell=True)
    os.chdir(pwd)

    ## copy overview ##
    overview_file = directory + os.sep + 'overview.png'
    if(os.path.isfile(overview_file)):
        os.system('cp {0} {1}'.format(overview_file, outfile))


def create_overview(indir, outdir):
    pool = Pool(processes=4)

    directories = sorted(glob.glob('{0}/*'.format(indir)))
    pool.map(plot_td, directories)
    # plots = map(plot_td, directories)
    # for p in plots:
    #     print(p)
    # for directory in directories:
    #     plot_td(directory)


if __name__ == '__main__':
    # create_overview('SIMS', 'output_SIMS')
    create_overview('SIM2', 'output_SIMS2')
