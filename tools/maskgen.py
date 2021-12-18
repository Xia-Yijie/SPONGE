import argparse
import os

parser = argparse.ArgumentParser(description="Use VMD to generate a file to record the atom indexes of the corresponding mask.\n")
parser.add_argument("-p",required=True, help="the parm file")
parser.add_argument("-c", help="the restart file")
parser.add_argument("-o",required=True, help="the output file")
parser.add_argument("--vmd", metavar = "vmd", default="vmd", help="the command to start vmd")
args = parser.parse_args()

s = input("Please Enter Your Selection Mask:\n")

p = args.p.split(os.path.sep)
p = "/".join(p)

c=""
if args.c:
  c = args.c.split(os.path.sep)
  c = "/".join(c)

temp_write = """set f [open "{0}" "w"]
mol load parm7 {1}
mol addfile {2} type rst7
atomselect top "{3}"
puts $f [atomselect0 list]
close $f
quit
""".format(args.o, p, c, s)

temp = open("maskgen_temp_tcl_file","w")
temp.write(temp_write)
temp.close()

os.system("{0} -dispdev none -e maskgen_temp_tcl_file".format(args.vmd))
os.remove("maskgen_temp_tcl_file")
