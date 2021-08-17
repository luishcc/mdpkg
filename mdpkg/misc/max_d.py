import numpy as np

from readLammps import DumpReader


def max(list):
    m = 0
    i = 0
    total = len(list)
    for atom in list:
        print(i, total)
        x, y, z = atom.x
        dist = np.sqrt(x**2 + y**2)
        if dist > m:
            m = dist*1
        i+=1
    return m


data = DumpReader('dump.atom')
print(max(data.atoms))
