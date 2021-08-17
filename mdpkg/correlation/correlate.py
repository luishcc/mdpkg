import numpy as np
import sys
import os
sys.path.insert(0, os.path.expanduser('~')+'/md-projects/lampy')
from grid import Grid
from readLammps import DumpReader
from math import floor

from scipy.fft import fft, fftfreq, fftshift, rfft, rfftfreq

import matplotlib.pyplot as plt


# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.size": 14,
#     "font.sans-serif": ["Helvetica"]})


rd = DumpReader('thread.lammpstrj')
rd.map_snapshot_in_file()


for time in rd.timesteps.keys():
    print(time/100 - 100)
    rd.read_snapshot(time)
    grd = Grid(rd.snapshots[time], size = float(sys.argv[1]))

    rrange = int(sys.argv[2])
    num = floor(grd.num_z/2)

    plt.figure(1)
    for r in range(rrange):

        a = grd.compute_density_correlation(r)
        plt.plot(np.linspace(0, grd.length_z/2, num), a, label=f'R={r}')

        if r == 8:
            # f = fftshift(fft(a))
            # freq = fftshift(fftfreq(len(a)))
            f = rfft(a) / len(a)
            freq = rfftfreq(len(a))

            plt.figure(2)
            plt.plot(freq[:], f.real[:], 'k-', marker='o')
            # plt.plot([min(freq[1:20]), max(freq[1:20])], [0,0], 'r-.')
            # plt.ylim(-0.6,3)
            plt.plot([0, 0.15], [0, 0], 'b--')
            plt.xlim(-0.01,0.15)
            plt.savefig(f'gif2/{time}.png', format='png')
            plt.close(2)

        if float('Nan') in a:
            continue

    plt.xlabel(r'$\delta z$')
    plt.ylabel(r'$G(r,\delta z)$')
    plt.ylim(-0.05, 1.1)
    plt.plot([0, grd.length_z/2], [0, 0], 'k--')
    plt.legend(loc='right')
    plt.savefig(f'gif/{time}.png', format='png')
    plt.close(1)

    del rd.snapshots[time]
