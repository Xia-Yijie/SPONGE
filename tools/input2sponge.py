import argparse
import os
import numpy as np
import parmed as pmd
from parmed.constants import TINY
    
parser = argparse.ArgumentParser(description="convert the input files from AMBER or GROMACS to SPONGE")

parser.add_argument("-f", required=True, help="the original program: AMBER or GROMACS", metavar="file_from", choices = ["AMBER","GROMACS"])
parser.add_argument("-p", required=True, help="the original topology file. .parm7 file for AMBER and .top file for GROMACS", metavar="topology_file")
parser.add_argument("-c", required=True, help="the original coordinate file. .rst7 file for AMBER and .gro file for GROMACS", metavar="coordinate_file")
parser.add_argument("-folder", required=True, help="the folder for the output files. If not exist, it will be made.", metavar="folder")
parser.add_argument("-prefix", required=True, help="the prifix for the output files", metavar="prefix")

args = parser.parse_args()

inputs = pmd.load_file(args.p, xyz = args.c) 
inputs2 = pmd.load_file(args.c)

folder = os.path.abspath(args.folder)
if not os.path.exists(folder):
    os.mkdir(folder)

atom_numbers = len(inputs.atoms)

coordinate_file_name = os.path.join(folder,"%s_coordinate.txt"%args.prefix)

coordinate_file = open(coordinate_file_name,"w")
towrite = "%d\n"%(atom_numbers)
if args.f == "AMBER":
    coordinates = inputs.coordinates
elif args.f == "GROMACS":
    coordinates = inputs.coordinates

for coordinate in coordinates:
    towrite += "%f %f %f\n"%(coordinate[0],coordinate[1],coordinate[2])
towrite += "%f %f %f %f %f %f\n"%(inputs2.box[0],inputs2.box[1],inputs2.box[2],inputs2.box[3],inputs2.box[4],inputs2.box[5])
coordinate_file.write(towrite)
coordinate_file.close()
print("c =",coordinate_file_name)
del coordinates

write_velocity = True
velocities = None
if args.f == "AMBER":
    if not inputs2.hasvels:
        write_velocity = False
    else:
        velocities = inputs2.velocities[0]
elif args.f == "GROMACS":
    velocities = inputs2.velocities

if write_velocity and (type(velocities) != type(None)):
    velocity_file_name = os.path.join(folder,"%s_velocity.txt"%args.prefix)
    velocity_file = open(velocity_file_name,"w")
    towrite = "%d\n"%(atom_numbers)
    for velocity in velocities:
        velocity /= 20.455
        towrite += "%f %f %f\n"%(velocity[0], velocity[1], velocity[2])
    velocity_file.write(towrite)
    velocity_file.close()
    print("v0 =",velocity_file_name)
del write_velocity
del velocities

residue_numbers = len(inputs.residues)
residue_file_name = os.path.join(folder,"%s_residue.txt"%args.prefix)
residue_file = open(residue_file_name,"w")
towrite = "%d %d\n"%(atom_numbers, residue_numbers)
for residue in inputs.residues:
    towrite += "%d\n"%(len(residue))
residue_file.write(towrite)
print("residue_in_file =",residue_file_name)
residue_file.close()

mass_file_name = os.path.join(folder,"%s_mass.txt"%args.prefix)
mass_file = open(mass_file_name,"w")
towrite = "%d\n"%(atom_numbers)
for atom in inputs.atoms:
    towrite += "%f\n"%(atom.mass)
mass_file.write(towrite)
print("mass_in_file =",mass_file_name)
mass_file.close()

charge_file_name = os.path.join(folder,"%s_charge.txt"%args.prefix)
charge_file = open(charge_file_name,"w")
towrite = "%d\n"%(atom_numbers)
for atom in inputs.atoms:
    towrite += "%f\n"%(atom.charge * 18.2223)
charge_file.write(towrite)
print("charge_in_file =",charge_file_name)
charge_file.close()


#把GROMACS中virtual atom和constraint的假键给删掉
#注意这会把键删去，但是不会把bond_partner给删去
if args.f == "GROMACS":
    i = 0
    while i < len(inputs.bonds):
        bond = inputs.bonds[i]
        if bond.type.k < TINY:
            inputs.bonds.remove(bond)
            i -= 1
        i += 1


#处理AMBER中的虚拟位点信息
#AMBER中只有水支持虚拟位点
#AMBER中的constraint在mdin中控制，而非拓扑文件中
elif args.f == "AMBER":
    virtual_atoms = []
    for residue in inputs.residues:
        if residue.name in ("WAT", "HOH"):
            if len(residue.atoms) == 4:
                virtual_bond = residue.atoms[3].bonds[0]
                if residue.atoms[1].bonds[0].atom1 == residue.atoms[0] or residue.atoms[1].bonds[0].atom2 == residue.atoms[0]:
                    ho_bond = residue.atoms[1].bonds[0].type.req
                    hh_bond = residue.atoms[1].bonds[1].type.req
                else:
                    ho_bond = residue.atoms[1].bonds[1].type.req
                    hh_bond = residue.atoms[1].bonds[0].type.req
                factor = virtual_bond.type.req / (4 * ho_bond * ho_bond - hh_bond * hh_bond) ** 0.5
                
                
                virtual_atoms.append([2, residue.atoms[3].idx, residue.atoms[0].idx, 
                    residue.atoms[1].idx, residue.atoms[2].idx, factor ,factor])
                inputs.bonds.remove(virtual_bond)
                del virtual_bond
                del ho_bond
                del hh_bond
    if virtual_atoms:
        virtual_atom_file_name = os.path.join(folder,"%s_vatom.txt"%args.prefix)
        virtual_atom_file = open(virtual_atom_file_name,"w")
        towrite = ""
        for virtual_atom in virtual_atoms:
            temp = "{} "*len(virtual_atom)
            temp = temp[:-1] + "\n"
            towrite += temp.format(*virtual_atom)

        virtual_atom_file.write(towrite)
        virtual_atom_file.close()
        print("virtual_atom_in_file = %s"%(virtual_atom_file_name))
    del virtual_atoms
        
    
bond_numbers = len(inputs.bonds)
if bond_numbers:
    bond_file_name = os.path.join(folder,"%s_bond.txt"%args.prefix)
    bond_file = open(bond_file_name,"w")
    towrite = "%d\n"%(bond_numbers)
    for bond in inputs.bonds:
        towrite += "%d %d %f %f\n"%(bond.atom1.idx, bond.atom2.idx, bond.type.k, bond.type.req)
    bond_file.write(towrite)
    bond_file.close()
    print("bond_in_file =",bond_file_name)

angle_numbers = len(inputs.angles)
if angle_numbers:
    angle_file_name = os.path.join(folder,"%s_angle.txt"%args.prefix)
    angle_file = open(angle_file_name,"w")
    towrite = "%d\n"%(angle_numbers)
    for angle in inputs.angles:
        towrite += "%d %d %d %f %f\n"%(angle.atom1.idx, angle.atom2.idx, angle.atom3.idx, angle.type.k, angle.type.utheteq.value_in_unit(pmd.unit.radian))
    angle_file.write(towrite)
    angle_file.close()
    print("angle_in_file =",angle_file_name)

dihedral_numbers = len(inputs.dihedrals)
if dihedral_numbers:
    dihedral_numbers = 0
    towrite = ""
    for dihedral in inputs.dihedrals:
        if dihedral.type.phi_k != 0:
            dihedral_numbers += 1
            towrite += "%d %d %d %d %d %f %f\n"%(dihedral.atom1.idx, dihedral.atom2.idx, dihedral.atom3.idx, dihedral.atom4.idx,
                dihedral.type.per, dihedral.type.phi_k, dihedral.type.uphase.value_in_unit(pmd.unit.radian))
    if dihedral_numbers:
        dihedral_file_name = os.path.join(folder,"%s_dihedral.txt"%args.prefix)
        dihedral_file = open(dihedral_file_name,"w")
        towrite = "%d\n"%(dihedral_numbers) + towrite
        dihedral_file.write(towrite)
        dihedral_file.close()
        print("dihedral_in_file =",dihedral_file_name)

    #may need fix for gromacs
    nb14_numbers = 0
    towrite = ""
    for dihedral in inputs.dihedrals:
        if dihedral.type.scnb != 0 and dihedral.type.scee != 0 and not dihedral.ignore_end:
            nb14_numbers += 1
            towrite += "%d %d %f %f\n"%(dihedral.atom1.idx, dihedral.atom4.idx,
                1.0 / dihedral.type.scnb, 1.0 / dihedral.type.scee)
    if nb14_numbers:
        towrite = "%d\n"%(nb14_numbers) + towrite
        nb14_file_name = os.path.join(folder,"%s_nb14.txt"%args.prefix)
        nb14_file = open(nb14_file_name,"w")
        nb14_file.write(towrite)
        nb14_file.close()
        print("nb14_in_file =",nb14_file_name)



if args.f == "AMBER":
    LJ_depth = inputs.LJ_depth
    LJ_radius = inputs.LJ_radius
    atom_type_numbers = len(LJ_depth)
    towrite = "%d %d\n\n"%(atom_numbers, atom_type_numbers)
    LJ_idx = []
    for atom in inputs.atoms:
        LJ_idx.append(atom.nb_idx - 1)
elif args.f == "GROMACS":
    atom_type_numbers = 0
    atom_types = []
    temp_factor = 2 ** (-5/6)
    LJ_depth = []
    LJ_radius = []
    LJ_idx = []
    for atom in inputs.atoms:
        atom_type = atom.atom_type
        found = atom_type_numbers
        for i in range(atom_type_numbers):
            if atom_type == atom_types[i]:
                found = i
                break
        if found == atom_type_numbers:
            atom_type_numbers += 1
            atom_types.append(atom_type)
        LJ_depth.append(atom_type.epsilon)
        LJ_radius.append(temp_factor * atom_type.sigma)
        LJ_idx.append(found)
    del atom_types
    del temp_factor
    del found
    towrite = "%d %d\n\n"%(atom_numbers, atom_type_numbers)    
        
    
for i in range(atom_type_numbers):
    temp = []
    for j in range(i+1):
        temp.append("%16.8e"%(np.sqrt(LJ_depth[i] * LJ_depth[j]) * ((LJ_radius[i] + LJ_radius[j])) ** 12))
    towrite += " ".join(temp)
    towrite += "\n"
towrite += "\n"
for i in range(atom_type_numbers):
    temp = []
    for j in range(i+1):
        temp.append("%16.8e"%(np.sqrt(LJ_depth[i] * LJ_depth[j]) * 2 * ((LJ_radius[i] + LJ_radius[j])) ** 6))
    towrite += " ".join(temp)
    towrite += "\n"
towrite += "\n"
for i in range(atom_numbers):
    towrite += "%d\n"%(LJ_idx[i])

lj_file_name = os.path.join(folder,"%s_lj.txt"%args.prefix)
lj_file = open(lj_file_name,"w")
lj_file.write(towrite)
lj_file.close()
print("lj_in_file =",lj_file_name)
del atom_type_numbers
del LJ_depth
del LJ_radius
del LJ_idx

towrite = ""
excluded_atom_numbers = 0
for atom in inputs.atoms:
    temp = 0
    exclusions_temp = atom.bond_partners + atom.angle_partners + atom.dihedral_partners + atom.tortor_partners + atom.exclusion_partners
    exclusions = []
    for atom_e in exclusions_temp:
        if (atom_e > atom):
            exclusions.append("%d"%atom_e.idx)
            temp += 1
    exclusions.sort(key=lambda x:int(x))
    towrite += "%d %s\n"%(temp, " ".join(exclusions))
    excluded_atom_numbers += temp
    del exclusions_temp
    del exclusions
    del temp
towrite = "%d %d\n"%(atom_numbers, excluded_atom_numbers) + towrite
exclude_file_name = os.path.join(folder,"%s_exclude.txt"%args.prefix)
exclude_file = open(exclude_file_name,"w")
exclude_file.write(towrite)
exclude_file.close()
print("exclude_in_file =",exclude_file_name)



#其他复杂情况的处理
#GROMACS的处理
if args.f == "GROMACS":
    from contextlib import closing
    from parmed.gromacs._gromacsfile import GromacsFile
    from parmed import gromacs as gmx
    molecules = {}
    system = []
    with closing(GromacsFile(args.p, includes=[gmx.GROMACS_TOPDIR], defines=None)) as f:
        current_section = None
        current_molecule = None
        constraints = []
        virtual_atoms = []
        for line in f:
            line = line.strip()
            if not line: 
                continue
            if line[0] == '[':       
                current_section = line[1:-1].strip()
                if current_molecule and current_section not in ('moleculetype', 'system', 'molecules'):
                    current_molecule[current_section] = []
            elif current_section == 'moleculetype':
                if current_molecule:
                    molecules[current_molecule['molname']] = current_molecule
                molname, nrexcl = line.split()
                current_molecule = {'molname': molname}
            elif current_section in ('atoms'):
                words = int(line.split()[0])
                current_molecule[current_section].append(words)
            elif current_section in ('constraints'):
                words = line.split()
                current_molecule[current_section].append([int(words[0]), int(words[1]), float(words[3]) * 10.0]);
            elif current_section in ('settles'):
                words = line.split()
                words[0] = int(words[0])
                current_molecule[current_section].append([words[0], words[0] + 1, float(words[2]) *10.0]);
                current_molecule[current_section].append([words[0], words[0] + 2, float(words[2]) *10.0]);
                current_molecule[current_section].append([words[0] + 1, words[0] + 2, float(words[3]) *10.0]);
            elif current_section in ('virtual_sites2', 'dummies2'):
                words = line.split()
                words[0] = int(words[0])
                words[1] = int(words[1])
                words[2] = int(words[2])
                words[4] = float(words[4])
                current_molecule[current_section].append([1, words[0], words[1], words[2], words[4]])
            elif current_section in ('virtual_sites3', 'dummies3'):
                words = line.split()
                words[0] = int(words[0])
                words[1] = int(words[1])
                words[2] = int(words[2])
                words[3] = int(words[3])
                words[4] = int(words[4])
                words[5] = float(words[5])
                words[6] = float(words[6])
                temp = []
                if words[4] == 1:
                    temp.append(2)
                elif words[4] == 2:
                    temp.append(3)
                    words[5], words[6] = words[6] * 10, words[5]
                else:
                    print("Ignored %s. Becuase not support for funct %d of virtual_sites3\n"%(line, words[4]))
                    continue
                temp.append(words[0])
                temp.append(words[1])
                temp.append(words[2])
                temp.append(words[3])
                temp.append(words[5])
                temp.append(words[6])
                current_molecule[current_section].append(temp)
            elif current_section in ('system'):
                if current_molecule:
                    molecules[current_molecule['molname']] = current_molecule
            elif current_section in ('molecules'):
                words = line.split()
                current_molecule = molecules[words[0]]
                current_atoms = -1
                molecule_numbers = int(words[1])
                for i in range(molecule_numbers):
                    #constraints
                    if "constraints" in current_molecule.keys():
                        for j in current_molecule["constraints"]:
                            constraints.append([j[0] + current_atoms, j[1] + current_atoms, j[2]])
                    if "settles" in current_molecule.keys():
                        for j in current_molecule["settles"]:
                            constraints.append([j[0] + current_atoms, j[1] + current_atoms, j[2]])
                    #virtual_atoms
                    if "virtual_sites3" in current_molecule.keys():
                        for j in current_molecule["virtual_sites3"]:
                            virtual_atoms.append("%d %d %d %d %d %f %f"%(j[0], j[1] + current_atoms, 
                            j[2] + current_atoms, j[3] + current_atoms, j[4] + current_atoms,
                            j[5], j[6]))
                    if "dummies3" in current_molecule.keys():
                        for j in current_molecule["dummies3"]:
                            virtual_atoms.append("%d %d %d %d %d %f %f"%(j[0], j[1] + current_atoms, 
                            j[2] + current_atoms, j[3] + current_atoms, j[4] + current_atoms,
                            j[5], j[6]))
                    if "virtual_sites2" in current_molecule.keys():
                        for j in current_molecule["virtual_sites2"]:
                            virtual_atoms.append("%d %d %d %d %f"%(j[0], j[1] + current_atoms, 
                            j[2] + current_atoms, j[3] + current_atoms, j[4]))
                    if "dummies2" in current_molecule.keys():
                        for j in current_molecule["dummies2"]:
                            virtual_atoms.append("%d %d %d %d %f"%(j[0], j[1] + current_atoms, 
                            j[2] + current_atoms, j[3] + current_atoms, j[4]))
                    
                    current_atoms += len(current_molecule["atoms"])
        #体系处理
        if constraints:
            constrain_file_name = os.path.join(folder,"%s_constrain.txt"%args.prefix)
            constrain_file = open(constrain_file_name,"w")
            towrite = ""
            for constraint in constraints:
                towrite += "{} {} {}\n".format(*constraint)
            constrain_file.write(towrite)
            constrain_file.close()
            print("constrain_in_file = %s"%(constrain_file_name))
        if virtual_atoms:
            virtual_atom_file_name = os.path.join(folder,"%s_vatom.txt"%args.prefix)
            virtual_atom_file = open(virtual_atom_file_name,"w")
            towrite = ""
            for virtual_atom in virtual_atoms:
                temp = "{} " * len(virtual_atom)
                temp = temp[:-1] + "\n"
                towrite += temp.format(*virtual_atom)
            virtual_atom_file.write(towrite)
            virtual_atom_file.close()
            print("virtual_atom_in_file = %s"%(virtual_atom_file_name))



