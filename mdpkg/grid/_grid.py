import numpy as np
from math import floor


class Grid:

    ''' Grid in cylindrical coordinates for liquid thread simulations.
    Atoms are connected over z periodic boundary '''

    def __init__(self, snap, size=2):
        self.cell = {}
        self.size = size

        self.length_z = snap.box.get_length_z()
        self.num_z = round(self.length_z / self.size)
        self.size_z = self.length_z / (self.num_z)


        for atom in snap.atoms.values():
            idr, idp = self.get_idpolar(atom.position)
            idz = self.get_idz(atom.position)
            try:
                self.cell[(idr, idp,  idz)].add_atom(atom)
            except:
                self.cell[(idr, idp,  idz)] = Cell(idr, idp,  idz)
                self.cell[(idr, idp,  idz)].set_nangle(self.get_numphi(idr))
                self.cell[(idr, idp,  idz)].compute_volume(self.size)
                self.cell[(idr, idp,  idz)].add_atom(atom)

        self.ncells = len(self.cell)


    def get_idpolar(self, pos):

        ''' Calculate the IDs of radius coordinate and angular coordinate
        for a cell volume approximately close to cartesian, i.e. different
        angular division for different radius'''

        r = np.sqrt(pos[0]**2 + pos[1]**2)
        idr = self.get_idr(r)
        N = self.get_numphi(idr)
        bin_theta = 2*np.pi / N
        angle = np.arctan2(pos[1], pos[0])+np.pi
        idp = floor(angle / bin_theta)
        return idr, idp


    def get_idr(self,r):
        return floor(r / self.size)

    def get_idz(self, pos):
        return floor(pos[2] / self.size_z)

    def get_numphi(self, idr):
        return round(np.pi*(2*idr+1))

    def compute_density_correlation(self, r):

        idr = self.get_idr(r)

        def correlate(dz):
            num_phi = self.get_numphi(idr)
            sumsq = 0
            corr = 0
            for z in range(self.num_z - dz):
                for phi in range(num_phi):
                    try:
                        sumsq += self.cell[(idr, phi, z)].get_density()**2
                        corr += self.cell[(idr, phi, z)].get_density() \
                              * self.cell[(idr, phi, z+dz)].get_density()
                    except KeyError:
                        continue
            for i in range(z+1, self.num_z):
                for phi in range(num_phi):
                    try:
                        sumsq += self.cell[(idr, phi, i)].get_density()**2
                    except KeyError:
                        continue
            nn = (self.num_z-dz) * num_phi
            nn2 = (self.num_z) * num_phi
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


class Cell:

    def __init__(self, idr, idp, idz):
        self.atoms = []
        self.id = (idr, idp, idz)
        self.volume = None
        self.nangle = None
        self.angle = None
        self.density = None

    def add_atom(self, atom):
        self.atoms.append(atom)

    def set_nangle(self, nangle):
        self.nangle = nangle
        self.angle = 2*np.pi / self.nangle

    def compute_volume(self, size):
        self.volume = self.angle * size**3 * (2*self.id[0] + 1) * 0.5

    def get_density(self):
        return len(self.atoms) / self.volume


if __name__=='__main__':

    import sys
    import os
    sys.path.insert(0, os.path.expanduser('~')+'/md-projects/lampy')
    from readLammps import DumpReader, Reader

    rd = Reader(sys.argv[1])
    rd.map_snapshot_in_file()
    rd.read_snapshot(2500)
    data = rd.snapshots[2500]

    # data = DumpReader(sys.argv[1])


    grd = Grid(data, size = float(sys.argv[2])*0.75)

    idr = []
    idz = []
    d = []
    for key, cell in grd.cell.items():
        idr.append(cell.id[0])
        idz.append(cell.id[2])
        d.append(cell.get_density()/cell.nangle)

    from scipy.sparse import coo_matrix
    coo = coo_matrix((d, (idr, idz)))

    coo = coo.todense().transpose()

    import matplotlib.pyplot as plt


    #plt.figure(1)
    #plt.tricontourf(idr,idz,d)
    #plt.colorbar()
    #plt.show()
    #exit()

    fig, ax = plt.subplots(1,1)
    im = ax.imshow(coo, extent=[0, 1, 0, 1], aspect=10)
    fig.colorbar(im)
    ax.set_xlabel('Radius')
    ax.set_ylabel('Length')
    ax.yaxis.set_major_locator(plt.NullLocator()) # remove y axis ticks
    ax.xaxis.set_major_locator(plt.NullLocator()) # remove x axis ticks
    # fig.set_dpi(460)
    plt.show()
