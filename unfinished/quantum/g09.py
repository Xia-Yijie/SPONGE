PyRun_SimpleString(R"XYJ( #直接python跑注释掉第一行和最后一行"
import numpy as np
import sponge
import re
from subprocess import run, PIPE

def Get_Results_From_Output(outcome, atom_numbers, debug = 0):
  energy = re.findall("SCF Done:\s*E\(.*\)\s*=[\s]*[+-]?\d+.\d*\s",outcome,re.I)
  if energy:
    energy = float(energy[0].split("=")[1].split()[0])
    if debug:
      print(energy)
  else:
    raise Exception("energy not found.")
  line = "Center\s*Atomic\s*Forces" + ".*\n"*(3 + atom_numbers)
  forces = re.findall(line,outcome,re.I)[0]
  frc = np.zeros((atom_numbers, 3),dtype=np.float32)
  if forces:
    forces = forces.split("\n")[3:3 + atom_numbers]
    for fi,force in enumerate(forces):
      temp = force.split()
      frc[fi][0] =float(temp[2])
      frc[fi][1] =float(temp[3])
      frc[fi][2] =float(temp[4])
  else:
    raise Exception("forces not found.")
  if debug:
    print(frc)
  line =  "Mulliken charges:.*Sum of Mulliken charges"
  charges = re.findall(line,outcome,re.I | re.DOTALL)

  if charges:
    line = "\d+\s*[A-Z]+\s*[+-]?\d+.\d*"
    charges = re.findall(line,charges[0],re.I)
    cg = np.zeros(atom_numbers, dtype=np.float32)

    for ci, charge in enumerate(charges):
      cg[ci] =float(charge.split()[2])
  if debug:
    print(cg)
  return energy, frc, cg

qm_command = "g09"
if sponge.command_exist("qm_command"):
    qm_command = sponge.command("qm_command")

qm_temp_file_name = "qm_tmp_file.gjf"
if sponge.command_exist("qm_temp_file_name"):
    qm_atoms = sponge.command("qm_temp_file_name")

if sponge.command_exist("qm_atoms"):
    qm_atoms = int(sponge.command("qm_atoms"))
else:
    raise Exception("qm_atoms needed")
    exit(1)

print("Initializing Quantum Mechanics Module")
labels = sponge.get_label(length = qm_atoms)

print("\tqm atom number is %d"%qm_atoms)

g09_keywords = ["force", "NoSymm"]

qm_max_memory = 1000
if sponge.command_exist("qm_max_memory"):
    qm_max_memory = int(sponge.command("qm_max_memory"))
print("\tqm max memory is %d Mb"%(qm_max_memory))

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
    qm_method = sponge.command("qm_method")
else:
    raise Exception("qm_method needed.")
    exit(1)
print("\tqm method is %s"%qm_method)

g09_keywords.append(qm_method+"/"+qm_basis)

if sponge.command_exist("qm_extra_key_words"):
    g09_keywords.append(sponge.original_command("qm_extra_key_words"))



sponge.add_print_head("_QM_ENE_")
print("Initializing Quantum Mechanics Module Finished")

qm_gjf_pre = "%%mem=%dMB\n#"%qm_max_memory + " ".join(g09_keywords) + "\n\n\n"
qm_gjf_pre += "%d %d\n"%(qm_charge, qm_spin)
qm_gjf_pre += "{}"+"\n"*4


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
   qm_gjf = qm_gjf_pre.format(temp)
   ftemp = open(qm_temp_file_name,"w")
   ftemp.write(qm_gjf)
   ftemp.close()
   ftemp = open(qm_temp_file_name,"r")
   out = run(qm_command, stdin=ftemp, stdout=PIPE, shell = True).stdout.decode(encoding="utf-8")
   ftemp.close()
   ene, force, charge = Get_Results_From_Output(str(out), qm_atoms)
   ene *= 627.51
   if sponge.need_atom_energy():
      atom_ene = sponge.get_atom_energy(start = 0, length = 1)
      atom_ene += ene
      sponge.set_atom_energy(atom_ene, start = 0, length = 1)
   
   force_old = sponge.get_force(length = qm_atoms)
   force *=  1186.2
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
     qm_gjf = qm_gjf_pre.format(temp)
     ftemp = open(qm_temp_file_name,"w")
     ftemp.write(qm_gjf)
     ftemp.close()
     ftemp = open(qm_temp_file_name,"r")
     out = run(qm_command, stdin=ftemp, stdout=PIPE, shell = True).stdout.decode(encoding="utf-8")
     ftemp.close()
     ene, force, charge = Get_Results_From_Output(str(out), qm_atoms)
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
