PyRun_SimpleString(R"XYJ(    #如果直接使用python跑，注释这行与最后一行。"
import numpy as np
import sponge
from subprocess import run, PIPE

def Get_Results_From_Output(FILE_NAME, atom_numbers, debug = 0):
    fp = open(FILE_NAME)
    contents = fp.read()
    fp.close()
    
    ene = float(contents.split("Final scf result")[1].split("E_tot =")[1].split()[0])
    if debug:
        print(ene)

    forces = np.zeros((atom_numbers,3))
    temp = contents.split("Gradient contribution from Tot-egrad\n")[1].split("\n")[:atom_numbers]
    for i, tempi in enumerate(temp):
        tempi = list(map(float, tempi.split()))
        forces[i][0] = tempi[1]
        forces[i][1] = tempi[2]
        forces[i][2] = tempi[3]
    if debug:
        print(forces)
    
    charges = np.zeros(atom_numbers)
    temp = contents.split("[Mulliken Population Analysis]\n")[1].split("\n")[1:1+atom_numbers]
    for i, tempi in enumerate(temp):
        charges[i] = float(tempi.split()[1])
    if debug:
        print(charges)
    return ene, forces, charges


print("Initializing Quantum Mechanics Module")


if sponge.command_exist("qm_command"):
    qm_command = sponge.command("qm_command")
else:
    raise Exception("qm_command needed")
    exit(1)

qm_temp_prefix = "qm_tmp_file"
if sponge.command_exist("qm_temp_prefix"):
    qm_atoms = sponge.command("qm_temp_prefix")

if sponge.command_exist("qm_atoms"):
    qm_atoms = int(sponge.command("qm_atoms"))
else:
    raise Exception("qm_atoms needed")
    exit(1)


labels = sponge.get_label(length = qm_atoms)

print("\tqm atom number is %d"%qm_atoms)

if sponge.command_exist("qm_charge"):
    qm_charge = int(sponge.command("qm_charge"))
else:
    qm_charge = int(np.sum(sponge.get_charge(length=qm_atoms)) / 18.2223)
print("\tqm charge is %d"%qm_charge)
    

qm_spin = 1
if sponge.command_exist("qm_spin"):
    qm_spin = int(sponge.command("qm_spin"))
print("\tqm spin is %d"%qm_spin)


qm_basis='6-31g'
if sponge.command_exist("qm_basis"):
    qm_basis = sponge.command("qm_basis")
print("\tqm basis is %s"%qm_basis)


if sponge.command_exist("qm_method"):
    qm_method = sponge.command("qm_method").lower()
else:
    raise Exception("qm_method needed.")
    exit(1)
    
if qm_method not in ("uhf", "rhf", "rohf"):
    qm_method = "dft\n %s"%qm_method
    if sponge.command_exist("qm_ks_method"):
        qm_method = sponge.command("qm_ks_method").lower() + "\n" + qm_method
    else:
        raise Exception("qm_ks_method needed for Kohn-Sham Calculation. RKS/UKS/ROKS.")
        exit(1)
print("\tqm method is %s"%qm_method.replace("\n", " "))


sponge.add_print_head("_QM_ENE_")
print("Initializing Quantum Mechanics Module Finished")

qm_inp_pre = """$COMPASS
Title
 sponge
Basis
 %s
Geometry
{}
End geometry
skeleton
nosymm
$end

$XUANYUAN
Direct
Schwarz
$END

$SCF
%s
spin
 %d
charge
 %d
#coulpot+cosx
$END
 

# GS-Gradient
$resp
geom
norder
 1
method
 1
#coulpot+cosx
grid
 fine
numinttype
 0
$end"""%(qm_basis, qm_method, qm_spin, qm_charge)


ene = 0
old_crd = None

@sponge.register("Calculate_Force")
def cal_force():
    global ene, old_crd
    crd = sponge.get_coordinate(length = qm_atoms)
    old_crd = crd
    temp = ""
    for i in zip(labels,crd.tolist()):
        temp += "{} {} {} {}\n".format(i[0],*i[1])
    qm_inp = qm_inp_pre.format(temp)
    ftemp = open(qm_temp_prefix + ".inp", "w")
    ftemp.write(qm_inp)
    ftemp.close()
    run(qm_command + " " + qm_temp_prefix, stdout = PIPE, shell = True)
    ene, force, charge = Get_Results_From_Output(qm_temp_prefix + ".out", qm_atoms)
    ene *= 627.51
    if sponge.need_atom_energy():
        atom_ene = sponge.get_atom_energy(start = 0, length = 1)
        atom_ene += ene
        sponge.set_atom_energy(atom_ene, start = 0, length = 1)

    force_old = sponge.get_force(length = qm_atoms)
    force *=  -1186.2
    force = force.astype(np.float32)
    if sponge.need_virial:
        virial = sponge.get_virial()
        virial += np.sum(force * crd)
        sponge.set_virial(virial)
    force += force_old
    sponge.set_force(force, length = qm_atoms)
   
     
@sponge.register("Calculate_Energy")
def cal_ene():
    global ene, old_crd
    crd = sponge.get_coordinate(length = qm_atoms)
    if (crd==old_crd).all():
      pass
    else:
        temp = ""
        for i in zip(labels,crd.tolist()):
            temp += "{} {} {} {}\n".format(i[0],*i[1])
        qm_inp = qm_inp_pre.format(temp)
        ftemp = open(qm_temp_prefix + ".inp", "w")
        ftemp.write(qm_inp)
        ftemp.close()
        run(qm_command + " " + qm_temp_prefix, stdout = PIPE, shell = True)
        ene, force, charge = Get_Results_From_Output(qm_temp_prefix + ".out", qm_atoms)
        ene *= 627.51

@sponge.register("After_Calculate_Energy")
def after_cal_ene():
    global ene
    ene += sponge.get_potential_energy()
    sponge.set_potential_energy(ene)

@sponge.register("Print")
def main_print():
    global ene
    sponge.add_print("%8.2f"%ene)

)XYJ");
