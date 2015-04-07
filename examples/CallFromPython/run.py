#!/usr/bin/python
"""Some modules such as pandas change the cpu affinity of the python process,
limiting it to only one cpu. This prevents external programs, as e.g. CRTomo,
to be called using subprocess.call and facility multiple cpus.

This script shows the problem and how to reset the cpu affinity

Note that this problem does not always come up, but just in case...
"""
import subprocess
import os
# forgot what this todes
os.environ['OPENBLAS_MAIN_FREE'] = '1'


def pa():
    for line in open('/proc/self/status'):
        if 'Cpu' in line:
            print(line)
pa()
# this import will reduce the affinity
import pandas
pandas
# change the cpu affinity back to use all CPUs
# os.system("taskset -p 0xff %d" % os.getpid())
pa()
exit()
p = subprocess.Popen('CRTomo_dev_knuet',
                     shell=True)  # , stdout=subprocess.STDOUT)
p.wait()
