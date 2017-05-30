#!/usr/bin/python
"""list the compiler flags and the run time for each of the test cases
"""
import os

with open('settings.txt', 'r') as fid:
    for nr, line in enumerate(fid.readlines()):
        tomodir = 'output_tomodirs/{0:03}_tomodir'.format(nr)
        run_file = tomodir + '/inv/run.ctr'
        if os.path.isfile(run_file):
            run_lines = open(run_file, 'r').readlines()
            time_line = run_lines[-2].strip()
            if time_line.startswith('CPU'):
                run_time = time_line.replace('/', ' ')
            else:
                run_time = ' ' * 35

        print('{0} {1}'.format(run_time, line.strip()))
