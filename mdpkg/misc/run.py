from mpi4py import MPI
from lammps import lammps

lmp = lammps()
lmp.file('in.thread')

me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

MPI.Finalize()
