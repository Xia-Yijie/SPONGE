#!/usr/bin/env python
import argparse
from struct import unpack
from netCDF4 import Dataset
import numpy as np
import os
from itertools import islice

parser = argparse.ArgumentParser(description="convert a traj file from .dat to .nc")
parser.add_argument("-n",required=True,type=int, help="the number of atoms in the traj file",metavar="atom_numbers")
parser.add_argument("-frame",type=int, help="the number of frames you want convert. If not set, a value will be guessed from the file size.",metavar="frame_numbers")
parser.add_argument("-x",required=True, help="the name of the .dat traj file.",metavar="dat_traj_file")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-box", help="the name of the box traj file",metavar="box_traj_file")
group.add_argument("-xyz", type=float, nargs=3, help="the three constant box length",metavar=("X","Y","Z"))
parser.add_argument("-nc", required=True, help="the name of the .nc traj file.",metavar="nc_traj_file")
parser.add_argument("-dt", type=float, default=1.0, help="the time difference between two frames in the unit of ps. 1.0 for default.")
parser.add_argument("--title",default="DEFAULT", help="the title of the traj.")

args = parser.parse_args()

data=Dataset(args.nc,"w",format="NETCDF3_64BIT_OFFSET")
data.createDimension("frame",0)
data.createDimension("spatial",3)
data.createDimension("atom",args.n)
data.createDimension("cell_spatial",3)
data.createDimension("label",5)
data.createDimension("cell_angular",3)

data.title=args.title
data.application="SPONGE"
data.program="SPONGE"
data.programVersion="0.8"
data.Conventions="AMBER"
data.ConventionVersion="1.0"

data.createVariable("time","f",("frame"))
data.variables["time"].unit="picosecond"
data.createVariable("spatial","S1",("spatial"))
data.createVariable("coordinates","f",("frame","atom","spatial"))
data.variables["coordinates"].unit="angstrom"
data.createVariable("cell_spatial","S1",("cell_spatial"))
data.createVariable("cell_angular","S1",("cell_angular","label"))
data.createVariable("cell_lengths","f8",("frame","cell_spatial"))
data.variables["cell_lengths"].unit="angstrom"
data.createVariable("cell_angles","f8",("frame","cell_angular"))
data.variables["cell_angles"].unit="degree"

data.variables["cell_angular"][:]=[[b'a',b'l',b'p',b'h',b'a'],
                                   [b'b',b'e',b't',b'a',b' '],
                                   [b'g',b'a',b'm',b'm',b'a']]
data.variables["cell_spatial"][:]=[b'a',b'b',b'c']
data.variables["spatial"][:]=[b'x',b'y',b'z']


if not args.frame:
    size = os.path.getsize(args.x)
    args.frame = size // 12 // args.n
    print("frame = %d"%args.frame)

if args.box:
    with open(args.box) as fr:
        box = np.loadtxt(islice(fr,args.frame), dtype = np.float32)
    box = box.reshape((-1,6))
    data.variables["cell_lengths"][:] = box[:,:3]
    data.variables["cell_angles"][:] = box[:,3:]
    box = None
else:
    data.variables["cell_lengths"][:] = np.array(args.xyz) 
    data.variables["cell_angles"][:]=[90,90,90]
    
data.variables["time"][:]= np.arange(args.frame, dtype = np.float32) * args.dt

data.variables["coordinates"][:] = np.fromfile(args.x, dtype=np.float32)[:args.frame*args.n*3].reshape((args.frame, args.n, 3))
data.close()

