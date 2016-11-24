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


def create_overview(prefix):
    summary_dir = 'overview'

    if(not os.path.isdir(summary_dir)):
        os.makedirs(summary_dir)

    for directory in sorted(glob.glob('SIMS/*')):
        if not os.path.isdir(directory):
            continue
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
            outfile = summary_dir + '/' + os.path.basename(directory) + '.png'
            os.system('cp {0} {1}'.format(overview_file, outfile))

        # ## copy 10 Hz ##
        # tdir = directory + os.sep + 'invmod/05_10.000000'
        # cmplx, fpi = td.get_result_files(tdir)

        # for obj, name in ((cmplx.cre, 'cre'),
        #                   (cmplx.cim, 'cim')):
        #     in_mag_file = obj + '.png'
        #     out_mag_file = single_dir + os.sep + name
        #     out_mag_file += directory + '_10Hz.png'
        #     os.system('cp "{0}" "{1}"'.format(in_mag_file, out_mag_file))


if __name__ == '__main__':
    create_overview('parstudy_')
