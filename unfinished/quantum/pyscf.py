PyRun_SimpleString(R"XYJ( #直接python跑注释掉第一行和最后一行"
from pyscf import gto
import numpy as np
import sponge

if sponge.command_exist("qm_atoms"):
    qm_atoms = int(sponge.command("qm_atoms"))
else:
    raise Exception("qm_atoms needed")
    exit(1)

print("Initializing Quantum Mechanics Module")
labels = sponge.get_label(length = qm_atoms)

print("\tqm atom number is %d"%qm_atoms)

mol = gto.Mole()
mol.verbose= 0
mol.symmetry = False

mol.max_memory = 1000
if sponge.command_exist("qm_max_memory"):
    mol.max_memory = float(sponge.command("qm_max_memory"))
print("\tqm max memory is %.1f Mb"%(mol.max_memory))

if sponge.command_exist("qm_charge"):
    mol.charge = int(sponge.command("qm_charge"))
else:
    mol.charge = int(np.sum(sponge.get_charge(length=qm_atoms)) / 18.2223)
print("\tqm charge is %d"%mol.charge)
    

mol.spin = 0
if sponge.command_exist("qm_spin"):
    mol.spin = int(sponge.command("qm_spin"))
print("\tqm spin is %d"%mol.spin)


mol.basis='6-31g'
if sponge.command_exist("qm_basis"):
    mol.basis = sponge.command("qm_basis")
print("\tqm basis is %s"%mol.basis)

mff = None
if sponge.command_exist("qm_method"):
    mff = sponge.command("qm_method").lower()
else:
    raise Exception("qm_method needed.")
    exit(1)

xc_func = None

if mff == "rhf":
    from pyscf import scf
    mff = scf.RHF
elif mff == "uhf":
    from pyscf import scf
    mff = scf.UHF
elif mff == "rohf":
    from pyscf import scf
    mff = scf.ROHF
else:
    from pyscf import dft
    xc_func = mff
    if sponge.command_exist("qm_ks_method"):
        temp =  sponge.command("qm_ks_method").lower()
    else:
        raise Exception("qm_ks_method needed for Kohn-Sham Calculation. RKS/UKS/ROKS.")
        exit(1)
    if temp == "rks":
        mff = dft.RKS
    elif temp == "uks":
        mff = dft.UKS
    elif temp == "roks":
        mff = dft.ROKS
    else:
        raise Exception("qm_ks_method can only be RKS/UKS/ROKS.")
        exit(1)
    print("\tqm method is %s dft  %s"%(temp, xc_func))

sponge.add_print_head("_QM_ENE_")
print("Initializing Quantum Mechanics Module Finished")

ene = 0
old_crd = None

@sponge.register("Calculate_Force")
def cal_force():
   global mol, ene, old_crd, mff, qm_atoms
   crd = sponge.get_coordinate(length = qm_atoms)
   old_crd = crd
   
   mol.atom = list(zip(labels,crd))
   mol.build()
   
   mf = mff(mol)
   if xc_func:
       mf.xc = xc_func
   ene = mf.kernel()
   ene *= 627.51
   if sponge.need_atom_energy():
      atom_ene = sponge.get_atom_energy(start = 0, length = 1)
      atom_ene += ene
      sponge.set_atom_energy(atom_ene, start = 0, length = 1)

   g = mf.nuc_grad_method()
   force_old = sponge.get_force(length = qm_atoms)
   force = g.kernel()
   force *= -1186.2
   force = force.astype(np.float32)
   if sponge.need_virial:
      virial = sponge.get_virial()
      virial += np.sum(force * crd)
      sponge.set_virial(virial)
      
   force += force_old
   sponge.set_force(force, length = qm_atoms)
     
@sponge.register("Calculate_Energy")
def cal_ene():
   global mol, ene, old_crd, mff, qm_atoms
   crd = sponge.get_coordinate(length = qm_atoms)
   if (crd==old_crd).all():
      pass
   else:
      mol.atom = list(zip(labels,crd))
      mol.build()
   
      mf = mff(mol)
      if xc_func:
          mf.xc = xc_func
      ene = mf.kernel()
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
