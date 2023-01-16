import numpy as np
from math import floor


class Gridz:

    ''' Grid on only z axis of liquid thread simulations.
    Atoms are connected over z periodic boundary '''

    def __init__(self, snap, size=1.5, R=6):
        self.cell = {}
        self.size = size

        self.lx = snap.box.get_length_x()
        self.ly = snap.box.get_length_y()
        self.length_z = snap.box.get_length_z()
        self.num_z = round(self.length_z / self.size)
        self.size_z = self.length_z / (self.num_z)

        for atom in snap.atoms.values():
            # if atom.position[0]**2 + atom.position[1]**2 <= R**2:
            idz = self.get_idz(atom.position)
            try:
                self.cell[idz].add_atom(atom)
                self.cell[idz].set_density()
            except:
                self.cell[idz] = Cellz(idz, self.lx, self.ly, self.size_z)
                self.cell[idz].compute_volume(R)
                self.cell[idz].add_atom(atom)
            # else:
            #     continue

        self.ncells = len(self.cell)

    def get_idz(self, pos):
        return floor((pos[2]*.9999) / self.size_z) # .9999 because of edge coord from LAMMPS is possible

    def compute_density_correlation(self):

        def correlate(dz):
            sumsq = 0
            corr = 0
            for z in range(self.num_z - dz):
                try:
                    sumsq += self.cell[z].get_density()**2
                    corr += self.cell[z].get_density() \
                          * self.cell[z+dz].get_density()
                except KeyError:
                    continue
            for i in range(z+1, self.num_z):
                try:
                    sumsq += self.cell[i].get_density()**2
                except KeyError:
                    continue
            nn = (self.num_z-dz)
            nn2 = (self.num_z)
            try:
                corr /= (sumsq / nn2) * nn
            except ZeroDivisionError:
                return float('NaN')
            return corr

        num = floor(self.num_z/2)
        correlation = np.zeros(num)
        dz_list = np.zeros(num)
        for dz in range(num):
            correlation[dz] = correlate(dz)
        return correlation

    def set_forces(self):
        for id, cell in self.cell.items():
            cell.get_force()

    def set_velocities(self):
        for id, cell in self.cell.items():
            cell.get_velocity()


class Cellz:

    def __init__(self, idz, sx, sy, sz):
        self.atoms = []
        self.id = idz
        self.size = [sx, sy, sz]
        # self.density = None
        self.force = None
        self.velocity = None
        self.force_cylindrical = None
        self.velocity_cylindrical = None
        self.volume = None
        self._property = {}

    def add_atom(self, atom):
        self.atoms.append(atom)

    def get_attr(self, attribute):
        return self._property[attribute]

    def set_size(self, sx, sy, sz):
        self.size = [sx, sy, sz]

    def compute_volume(self, R):
        self.volume = 4 * self.size[0] * self.size[1] * self.size[2]
        # self.volume = np.pi * R**2 * self.size[2]
        # self.volume = 4 * R**2 * self.size[2]

    def set_density(self):
        self.density = len(self.atoms) / self.volume
        self._property['density'] = self.density
        return len(self.atoms) / self.volume

    def get_density(self):
        return len(self.atoms) / self.volume

    def get_force(self):
        f = [0, 0, 0]
        for atom in self.atoms:
            f = [a+b for a, b in zip(f, atom.force)]
        self.force = [i/len(self.atoms) for i in f]
        self._property['fx'] = self.force[0]
        self._property['fy'] = self.force[1]
        self._property['fz'] = self.force[2]
        return self.force


    def get_velocity(self):
        v = [0, 0, 0]
        for atom in self.atoms:
            v = [a+b for a, b in zip(v, atom.velocity)]
        self.velocity = [i/len(self.atoms) for i in v]
        self._property['vx'] = self.velocity[0]
        self._property['vy'] = self.velocity[1]
        self._property['vz'] = self.velocity[2]
        return self.velocity
