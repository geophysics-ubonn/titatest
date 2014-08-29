#!/usr/bin/env python
# -*- coding: utf-8 -*-

from crlab_py.CRTomo import *

new_opts = {}
new_opts['mswitch2'] = (4, 5, 6)
new_opts['lambda'] = (0.1, 0.5, 0.9)

crt_config = config()
init_list = [crt_config, ]
permutated_options = permutate_options(init_list, new_opts)

# write to files
for nr, i in enumerate(permutated_options):
    i.write_to_file('{0:02}_crtomo.cfg'.format(nr))
