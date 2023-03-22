import numpy as np
from math import floor


class GridSurf:

    ''' Grid on only z axis of liquid thread surface.
    Atoms are connected over z periodic boundary '''

    def __init__(self, snap, shape, thickness=1.5):
        self.cell = {}

        lx = snap.box.get_length_x()
        ly = snap.box.get_length_y()
        lz = snap.box.get_length_z()
        self.num_z = len(shape)
        self.size_z = lz / (self.num_z)

        for atom in snap.atoms.values():
            z = snap.box.zlo + atom.position[2]*lz
            if z > lz:
                z = lz*0.9999
            idz = self.get_idz(z)
            # print(self.num_z, idz, z, lz)
            Rup = shape[idz] + 0.1*thickness
            Rlow = shape[idz] - 0.9*thickness
            if Rlow < 0:
                Rlow = 0
            x = snap.box.xlo + atom.position[0]*lx
            y = snap.box.ylo + atom.position[1]*ly
            # print(Rlow, x**2 + y**2, Rup)
            if Rlow**2 <= x**2 + y**2 <= Rup**2:
                try:
                    self.cell[idz].add_atom(atom)
                    self.cell[idz].set_density()
                except:
                    self.cell[idz] = CellSurf(idz, thickness, self.size_z)
                    self.cell[idz].compute_volume(Rlow, Rup)
                    self.cell[idz].compute_area(Rup)
                    self.cell[idz].add_atom(atom)
            else:
                continue

        self.ncells = self.num_z

    def get_idz(self, z):
        return floor((z*1) / self.size_z) # .9999 because of edge coord from LAMMPS is possible

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


class CellSurf:

    def __init__(self, idz, sr, sz):
        self.atoms = []
        self.id = idz
        self.size = [sr, sz]
        self.types = {1:0, 2:0, 3:0}
        # self.density = None
        self.force = None
        self.velocity = None
        self.force_cylindrical = None
        self.velocity_cylindrical = None
        self.volume = None
        self._property = {}

    def add_type(self, type):
        self.types[type] += 1

    def add_atom(self, atom):
        self.atoms.append(atom)
        self.add_type(int(atom.type))

    def get_attr(self, attribute):
        return self._property[attribute]

    def set_size(self, sr, sz):
        self.size = [sr, sz]

    def compute_volume(self, Rlow, Rup):
        self.volume = np.pi * (Rup**2 - Rlow**2) * self.size[1]

    def compute_area(self, Rup):
        self.area = 2*np.pi * Rup * self.size[1]

    def get_density_type(self, type):
        return self.types[type] / self.volume

    def get_area_density_type(self, type):
        return self.types[type] / self.area

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
