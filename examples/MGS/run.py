#!/usr/bin/env python
# *-* coding: utf-8 *-*
import os
import numpy as np
import subprocess
import shutil
import crlab_py.CRTomo as CRTomo
import crlab_py.tomodir as TD


def create_sim(outdir, mgs_type, mgs_beta):
    td = TD.tomodir()
    td.create_directory_structure(outdir)

    cfg = CRTomo.config()

    cfg['cells_z'] = '-1'
    cfg['fpi_inv'] = 'F'

    cfg['mag_rel'] = '2'
    cfg['mag_abs'] = '1e-3'

    cfg['pha_abs'] = '0'

    cfg['d2_5'] = '1'
    cfg['fic_sink'] = 'F'

    cfg['mswitch2'] = '{0}'.format(mgs_type)
    cfg['lambda'] = '{0}'.format(mgs_beta)

    cfg.write_to_file(outdir + os.sep + 'exe/crtomo.cfg')
    print(cfg)

    shutil.copy('GRID/config.dat', outdir + os.sep + 'config/config.dat')
    shutil.copy('GRID/elem.dat', outdir + os.sep + 'grid/elem.dat')
    shutil.copy('GRID/elec.dat', outdir + os.sep + 'grid/elec.dat')
    shutil.copy('GRID/rho.dat', outdir + os.sep + 'rho/rho.dat')

    pwd = os.getcwd()
    os.chdir(outdir + os.sep + 'exe')
    subprocess.call('list_config_files.py crmod.cfg crmod.cfg', shell=True)
    os.chdir(pwd)

    noisefile = outdir + os.sep + 'exe/crt.noisemod'
    print(noisefile)
    with open(noisefile, 'w') as fid:
        # seed
        fid.write('1\n')
        # a = 1%
        fid.write('1\n')
        fid.write('1e-4\n')
        fid.write('0\n')
        fid.write('0\n')
        fid.write('0\n')
        fid.write('0\n')


if __name__ == '__main__':
    for mgs_type in range(5, 10):
        for mgs_beta in np.logspace(-5, -1, 5):
            directory = 'SIMS/mgs_type_{0}_beta_{1}'.format(mgs_type, mgs_beta)
            if not os.path.isdir(directory):
                os.makedirs(directory)
            create_sim(directory, mgs_type, mgs_beta)
