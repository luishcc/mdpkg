#   Module for reading LAMMPS Dump and Data output files
#


#----------------------------------------------------------------
class Atom:

    ''' Atom object: used to describe describe position and other properties
of a single atom/particle '''

    def __init__(self, id, p, v=None, type=None):
        self.id = id
        self.position = p
        self.velocity = v
        self.type = type
        self.properties = {}

    def set_property(self, name, vale):
        self.properties[name] = value

#----------------------------------------------------------------
class Box:

    ''' Box object: Contains x y z limits of simulation box.
TO DO: Add boundary condition type '''

    def __init__(self, list):
        self.xlo = list[0]
        self.xhi = list[1]
        self.ylo = list[2]
        self.yhi = list[3]
        try:
            self.zlo = list[4]
            self.zhi = list[5]
        except:
            self.zlo = None
            self.zhi = None

    def get_length_x(self):
        return self.xhi - self.xlo

    def get_length_y(self):
        return self.yhi - self.ylo

    def get_length_z(self):
        return self.zhi - self.zlo

#----------------------------------------------------------------
class Snapshot:

    ''' Snapshot object: Contains the collection of atoms/particles in a single
timestep. Objects are instantiated by a Reader object '''

    def __init__(self, time):
        self.time = time
        self.atoms = {}

    def set_box(self, list):
        self.box = Box(list)

#----------------------------------------------------------------
class DumpReader:

    ''' Main class to read Dump files. It reads single snapshots at a time and
works by first mapping where it snapshot begins in the file and then
the .read_snapshot( method can be called on a specific timestep '''

    def __init__(self, file_name):
        self.file_name = file_name
        self.timesteps = {}
        self.snapshots = {}

    def skip_lines(self, f, n_skip):
        for _ in range(n_skip):
            f.readline()

    def read_snapshot_at(self, time):
        ss = Snapshot(time)
        with open(self.file_name, 'r') as file:
            self.skip_lines(file, self.timesteps[time])
            self.read_snapshot(file)

    def map_snapshot_in_file(self):
        with open(self.file_name, 'r') as file:
            reading = True
            idl = 0
            while reading:
                try:
                    file.readline()
                    self.timesteps[int(file.readline())] = idl
                    file.readline()
                    natom = int(file.readline())
                    self.skip_lines(file, natom + 5)
                    idl += natom + 9
                except:
                    reading = False
                    continue

    def read_sequential(self):
        file = open(self.file_name, 'r')
        self.read_snapshot(file)
        self._file = file

    def read_next(self):
        self.read_snapshot(self._file)

    def close_read(self):
        self._file.close()

    def delete_snapshot(self, time):
        del self.snapshots[time]

    def skip_next(self, n=1):
        for i in range(n):
            self.skip_lines(self._file, 3)
            num = int(self._file.readline())
            self.skip_lines(self._file, num + 5)

    def read_snapshot(self, file):
        file.readline()
        time = int(file.readline())
        ss = Snapshot(time)
        file.readline()
        ss.natoms = int(file.readline())
        file.readline()
        xx = file.readline().split()
        yy = file.readline().split()
        zz = file.readline().split()
        xyz_lim = [None] * 6
        for i in range(2):
            xyz_lim[i] = float(xx[i])
            xyz_lim[i+2] = float(yy[i])
            xyz_lim[i+4] = float(zz[i])
        ss.set_box(xyz_lim)
        ss.dump_attributes = file.readline().split()[1:]
        for n in range(ss.natoms):
            line = file.readline().split()
            id = int(line[0])
            t = line[1]
            p = [float(line[2]), float(line[3]), float(line[4])]
            ss.atoms[id] = (Atom(id, p, type=t))
        self.snapshots[time] = ss



#----------------------------------------------------------------
class DataReader:

    def __init__(self, file):
        self.file_name = str(file)
        self.atoms = []
        #self.box = None
        #self.timestep = None
        self.parse()

    def parse(self):
        list = []

        file = open(self.file_name, 'r')

        linenumber = 0
        id = 0
        reading_atoms = False
        for line in file:

            if line.find('ITEM: ATOMS') >= 0:
                reading_atoms = True
                continue

            if reading_atoms:
                d = {}
                l = line.split()
                d['position'] = [float(l[2]), float(l[3]), float(l[4])]
                d['tag'] = l[1]
                d['velocity'] = None
                d['properties'] = None
                self.atoms.append(Atom(d, id))
                id += 1


    def parse_xyz_style(self):
        print('Implementation Incomplete / Not working')
        return None

    def delete_atom(self, id):
        pass


#----------------------------------------------------------------
if __name__=='__main__':

    # a = DumpReader('dump.atom')
    a = DumpReader('dump.many')
    a.map_snapshot_in_file()
