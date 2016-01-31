#!/usr/bin/python

nr_cells_x = 44
depth_of_interface = 4
start_cell = nr_cells_x * (depth_of_interface - 1)
strength = 0

with open('decouplings.dat', 'w') as fid:
    fid.write('{0}\n'.format(nr_cells_x))
    for i in range(0, nr_cells_x):
        cell_id = start_cell + i
        cell2_id = start_cell + i + nr_cells_x
        cell_id += 1
        cell2_id += 1
        fid.write('{0} {1} {2}\n'.format(cell_id, cell2_id, strength))
