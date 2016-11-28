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

    # cfg['mag_rel'] = '1'
    # cfg['mag_abs'] = '5e-4'
    # cfg['pha_abs'] = '0'
    cfg['mag_rel'] = '0.5'
    cfg['mag_abs'] = '5e-4'
    cfg['pha_abs'] = '0'

    cfg['d2_5'] = '1'
    cfg['fic_sink'] = 'F'

    cfg['mswitch2'] = '{0}'.format(mgs_type)
    cfg['lambda'] = '{0}'.format(mgs_beta)

    cfg.write_to_file(outdir + os.sep + 'exe/crtomo.cfg')
    print(cfg)

    shutil.copy('GRID2/config.dat', outdir + os.sep + 'config/config.dat')
    shutil.copy('GRID2/elem.dat', outdir + os.sep + 'grid/elem.dat')
    shutil.copy('GRID2/elec.dat', outdir + os.sep + 'grid/elec.dat')
    shutil.copy('GRID2/rho4.dat', outdir + os.sep + 'rho/rho.dat')

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
        # fid.write('1\n')
        fid.write('0.4\n')
        fid.write('1e-4\n')
        fid.write('0\n')
        fid.write('0\n')
        fid.write('0\n')
        fid.write('0\n')


def prepare_simulation(simdir):

    mgs_betas = np.hstack((
        np.logspace(-5, -1, 5),
        0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.48,
        0.5,
        0.51, 0.55,
        0.7, 0.71, 0.72, 0.73, 0.74, 0.75,
        0.8, 0.85, 0.9, 1.0,
    ))
    for mgs_type in (7, ):
        for mgs_beta in mgs_betas:
            directory = '{0}/mgs_type_{1}_beta_{2:.8f}'.format(
                simdir,
                mgs_type,
                mgs_beta
            )
            if not os.path.isdir(directory):
                os.makedirs(directory)
                create_sim(directory, mgs_type, mgs_beta)

    # we need ine reference unversion using the regular triangular reg.
    directory = simdir + '/smoothness_regularization'
    if not os.path.isdir(directory):
        os.makedirs(directory)
        create_sim(directory, 1, 0)

if __name__ == '__main__':
    prepare_simulation('SIM2')
