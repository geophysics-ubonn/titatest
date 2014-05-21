import copy


class config(dict):
    """
    Write CRTomo configuration files (crtomo.cfg).

    This class is essentially a dict of CRTomo configurations with a few extra
    functions that know how to write a proper crtomo.cfg file.
    """
    def __init__(self, *arg, **kw):
        super(config, self).__init__(*arg, **kw)
        self.set_defaults()

    def set_defaults(self):
        """
        Fill the dictionary with all defaults
        """
        self['mswitch'] = 1
        self['elem'] = '../grid/elem.dat'
        self['elec'] = '../grid/elec.dat'
        self['volt'] = '../mod/volt.dat'
        self['inv_dir'] = '../inv'
        self['diff_inv'] = 'F ! difference inversion?'
        self['iseed_var'] = 'iseed variance'
        self['cells_x'] = '0    ! # cells in x-direction'
        self['cells_z'] = '0    ! # cells in z-direction'
        self['ani_x'] = '1.000  ! smoothing parameter in x-direction'
        self['ani_z'] = '1.000  ! smoothing parameter in z-direction'
        self['max_it'] = '20    ! max. nr of iterations'
        self['dc_inv'] = 'F     ! DC inversion?'
        self['robust_inv'] = 'F     ! robust inversion?'
        self['fpi_inv'] = 'T     ! final phase improvement?'
        self['mag_rel'] = '5'
        self['mag_abs'] = '1e-4'
        self['pha_a1'] = 0
        self['pha_b'] = 0
        self['pha_rel'] = 0
        self['pha_abs'] = 0.5
        self['hom_bg'] = 'T'
        self['hom_mag'] = '10.00'
        self['hom_pha'] = '0.00'
        self['another_ds'] = 'F'
        self['2-5d'] = '0'
        self['fic_sink'] = 'T'
        self['fic_sink_node'] = '6467'
        self['boundaries'] = 'F'
        self['boundaries_file'] = 'boundary.dat'
        self['mswitch2'] = '1'
        self['lambda'] = 'lambda'

    def write_to_file(self, filename):
        """
        Write the configuration to a file. Use the correct order of values.
        """
        fid = open(filename, 'w')

        for key in ('mswitch', 'elem', 'elec', 'volt', 'inv_dir', 'diff_inv',
                    -1, -1, -1, 'iseed_var', 'cells_x', 'cells_z',
                    'ani_x', 'ani_z', 'max_it', 'dc_inv', 'robust_inv',
                    'fpi_inv', 'mag_rel', 'mag_abs', 'pha_a1', 'pha_b',
                    'pha_rel', 'pha_abs', 'hom_bg', 'hom_mag', 'hom_pha',
                    'another_ds', '2-5d', 'fic_sink', 'fic_sink_node',
                    'boundaries', 'boundaries_file', 'mswitch2', 'lambda'):
            if(key == -1):
                fid.write('\n')
            else:
                fid.write('{0}\n'.format(self[key]))

        fid.close()


def permutate_options(configs, new_opts):
    """
    Permutate configuration options

    Parameters
    ----------
    configs : list of initial CRTomo.config() instances
    new_opts : dictionary of CRTomo options to permutate. The keys of this dict
               correspond to keys in CRTomo.config, and the item holds the
               content to be permutated.

    Returns
    -------
    configs : A list of CRTomo.config instances with the permutated options set


    Example
    -------
    To permutate various regularization types and correponding lambdas, use

    new_opts = {}
    new_opts['mswitch2'] = (4, 5, 6)
    new_opts['lambda'] = (0.1, 0.5, 0.9)

    crt_config = CRTomo.config()
    init_list = [crt_config, ]
    permutated_options = permutate_options(init_list, new_opts)

    # write to files
    for nr, i in enumerate(permutated_options):
        i.write_to_file('{0:02}_crtomo.cfg'.format(nr))
    """
    # parameter to loop through
    key = new_opts.keys()[0]
    items = new_opts[key]
    del(new_opts[key])

    # loop through all configs in 'configs' and permutate with options in
    # new_dict
    template_configs = range(0, len(configs))
    for nr in template_configs:
        for item in items:
            tmp_config = copy.deepcopy(configs[nr])
            tmp_config[key] = item
            configs.append(tmp_config)

    # delete the templates
    for nr in reversed(template_configs):
        del(configs[nr])

    # if new_opts is empty...
    if((new_opts and True) or False):
        return permutate_options(configs, new_opts)
    else:
        return configs
