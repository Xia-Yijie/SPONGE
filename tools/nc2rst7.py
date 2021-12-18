import argparse
from netCDF4 import Dataset



parser = argparse.ArgumentParser(description="transform a rst7 file from .nc to .rst7")
parser.add_argument("-rst7",required=True, help="the name of the .rst7 restart file.")
parser.add_argument("-nc", required=True, help="the name of the .nc restart file.")
args = parser.parse_args()


data = Dataset(args.nc)

crds = data.variables["coordinates"][:]
crds = list(crds.reshape(-1,3))
vels = data.variables["velocities"][:]
vels = list(vels.reshape(-1,3))

cell_lengths = data.variables["cell_lengths"][:]
cell_angles = data.variables["cell_angles"][:]
title = data.title
time = data.variables["time"][0]
atom_numbers = data.dimensions["atom"].size

towrite="%s\n%d %f\n"%(title, atom_numbers, time)
for i,ci in enumerate(crds):
    towrite += " ".join(map(lambda x: "%12.7f"%x,ci))
    if i%2 == 1:
        towrite += "\n"
if i%2 == 0:
    towrite += "\n"
    
for i,ci in enumerate(vels):
    towrite += " ".join(map(lambda x: "%12.7f"%x,ci))
    if i%2 == 1:
        towrite += "\n"
if i%2 == 0:
    towrite += "\n"

towrite+= " ".join(map(lambda x: "%12.7f"%x,cell_lengths))
towrite+= " ".join(map(lambda x: "%12.7f"%x,cell_angles))

fw=open(args.rst7,"w")
fw.write(towrite)
fw.close()