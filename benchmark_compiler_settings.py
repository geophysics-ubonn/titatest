#!/usr/bin/env python
# *-* coding: utf-8 *-*
# create CRTomo binaries with different compiler settings
import os
import itertools
import hashlib
import re
import subprocess
import shutil


# these flags must always be used
compiler_flags_required = [
    '-fopenmp',
]

# these are permutated
compiler_flags_optional = [
    ('', '-O2', '-O3'),
    ('', '-march=native', ),
    ('', '-ftree-vectorize '),
    # check if this does change our results
    ('', '-ffast-math'),
    ('', '-funroll-loops'),
    ('', '-floop-nest-optimize'),
]


def create_permutate_binaries():
    # flags_permutated = itertools.permutations(compiler_flags_optional)
    flags_permutated = itertools.product(*compiler_flags_optional)
    flags_strings = []
    with open('settings.txt', 'w') as fid:
        for nr, flags in enumerate(flags_permutated):
            flag_str = ' '.join(compiler_flags_required) + ' ' + re.sub(
                ' +', ' ', ' '.join(flags).strip()
            )
            flags_strings.append(flag_str)
            fid.write(
                '{0} {1} {2}\n'.format(
                    nr, hashlib.md5(flag_str).hexdigest(), flag_str
                )
            )
    return flags_strings


def compile_flags(output_dir, flags):
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    computer_name = subprocess.check_output('uname -n', shell=True).strip()
    git_branch = subprocess.check_output(
        "git branch", shell=True
    ).strip()[2:]
    binary_name = 'src/CRTomo_' + git_branch + '_' + computer_name
    print('binary_name', binary_name)

    for nr, flag in enumerate(flags):
        output_file = output_dir + '/{0:03}_CRTomo'.format(nr)
        if os.path.isfile(output_file):
            print('skipping {0}'.format(output_file))
            continue

        if os.path.isfile(binary_name):
            os.unlink(binary_name)
        subprocess.call('make clean', shell=True)
        subprocess.call('./clean_autotools.sh', shell=True)
        os.environ['FCFLAGS'] = flag
        os.environ['FFFLAGS'] = flag
        subprocess.call('./autogen.sh')
        subprocess.call('make')
        if os.path.isfile(binary_name):
            shutil.move(
                binary_name,
                output_file,
            )


def main():
    flags = create_permutate_binaries()
    compile_flags('output_benchmark_binaries', flags)


if __name__ == '__main__':
    main()
