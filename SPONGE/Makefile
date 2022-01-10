

SUBDIRS=$(shell ls -l | grep ^d | awk '{print $$9}')
compiler=nvcc
COMMON_COMMAND=-arch=sm_50 -rdc=true -lcudadevrt -lcufft --use_fast_math -O4 -std=c++11
BIN_NAME=SPONGE

install: main.o common.o control.o MD_core/MD_core.o MD_core/FEP_core.o MD_core/TI_core.o bond/bond.o bond/bond_soft.o angle/angle.o angle/Urey_Bradley_force.o dihedral/dihedral.o dihedral/improper_dihedral.o cmap/cmap.o nb14/nb14.o neighbor_list/neighbor_list.o Lennard_Jones_force/Lennard_Jones_force.o Lennard_Jones_force/LJ_soft_core.o PME_force/PME_force.o  thermostats/Langevin_MD.o thermostats/Middle_Langevin_MD.o thermostats/Andersen_thermostat.o thermostats/Berendsen_thermostat.o thermostats/nose_hoover_chain.o barostats/MC_barostat.o barostats/Berendsen_barostat.o barostats/andersen_barostat.o restrain/restrain.o constrain/constrain.o constrain/simple_constrain.o constrain/SETTLE.o constrain/SHAKE.o virtual_atoms/virtual_atoms.o crd_molecular_map/crd_molecular_map.o
	$(compiler) -o $(BIN_NAME)  $^ $(COMMON_COMMAND)


all: main.o FEP.o TI.o common.o control.o MD_core/MD_core.o MD_core/FEP_core.o MD_core/TI_core.o bond/bond.o bond/bond_soft.o angle/angle.o angle/Urey_Bradley_force.o dihedral/dihedral.o dihedral/improper_dihedral.o cmap/cmap.o nb14/nb14.o neighbor_list/neighbor_list.o Lennard_Jones_force/Lennard_Jones_force.o Lennard_Jones_force/LJ_soft_core.o PME_force/PME_force.o  thermostats/Langevin_MD.o thermostats/Middle_Langevin_MD.o thermostats/Andersen_thermostat.o thermostats/Berendsen_thermostat.o thermostats/nose_hoover_chain.o barostats/MC_barostat.o barostats/Berendsen_barostat.o barostats/andersen_barostat.o restrain/restrain.o constrain/constrain.o constrain/simple_constrain.o constrain/SETTLE.o constrain/SHAKE.o virtual_atoms/virtual_atoms.o crd_molecular_map/crd_molecular_map.o
	
FEP: FEP.o common.o control.o MD_core/FEP_core.o bond/bond.o bond/bond_soft.o angle/angle.o dihedral/dihedral.o Lennard_Jones_force/LJ_soft_core.o nb14/nb14.o neighbor_list/neighbor_list.o PME_force/PME_force.o
	$(compiler) -o $(BIN_NAME)_FEP  $^ $(COMMON_COMMAND)

TI: TI.o common.o control.o MD_core/FEP_core.o bond/bond.o bond/bond_soft.o angle/angle.o dihedral/dihedral.o Lennard_Jones_force/LJ_soft_core.o nb14/nb14.o neighbor_list/neighbor_list.o PME_force/PME_force.o
	$(compiler) -o $(BIN_NAME)_TI  $^ $(COMMON_COMMAND)

main.o: main.cu main.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

FEP.o: FEP.cu FEP.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

TI.o: TI.cu TI.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

common.o: common.cu common.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


control.o: control.cu common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


MD_core/MD_core.o: MD_core/MD_core.cu MD_core/MD_core.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

MD_core/FEP_core.o: MD_core/FEP_core.cu MD_core/FEP_core.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

MD_core/TI_core.o: MD_core/TI_core.cu MD_core/TI_core.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

bond/bond.o: bond/bond.cu bond/bond.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

bond/bond_soft.o: bond/bond_soft.cu bond/bond_soft.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


angle/angle.o: angle/angle.cu angle/angle.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

angle/Urey_Bradley_force.o: angle/Urey_Bradley_force.cu angle/Urey_Bradley_force.cuh bond/bond.cuh angle/angle.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

dihedral/dihedral.o: dihedral/dihedral.cu dihedral/dihedral.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

dihedral/improper_dihedral.o: dihedral/improper_dihedral.cu dihedral/improper_dihedral.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

cmap/cmap.o: cmap/cmap.cu cmap/cmap.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

nb14/nb14.o: nb14/nb14.cu nb14/nb14.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


neighbor_list/neighbor_list.o: neighbor_list/neighbor_list.cu neighbor_list/neighbor_list.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


Lennard_Jones_force/Lennard_Jones_force.o: Lennard_Jones_force/Lennard_Jones_force.cu Lennard_Jones_force/Lennard_Jones_force.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

Lennard_Jones_force/LJ_soft_core.o: Lennard_Jones_force/LJ_soft_core.cu Lennard_Jones_force/LJ_soft_core.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

PME_force/PME_force.o: PME_force/PME_force.cu PME_force/PME_force.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

thermostats/Langevin_MD.o: thermostats/Langevin_MD.cu thermostats/Langevin_MD.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

thermostats/Middle_Langevin_MD.o: thermostats/Middle_Langevin_MD.cu thermostats/Middle_Langevin_MD.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


thermostats/Andersen_thermostat.o: thermostats/Andersen_thermostat.cu thermostats/Andersen_thermostat.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


thermostats/Berendsen_thermostat.o: thermostats/Berendsen_thermostat.cu thermostats/Berendsen_thermostat.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


thermostats/nose_hoover_chain.o: thermostats/nose_hoover_chain.cu thermostats/nose_hoover_chain.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


barostats/MC_barostat.o: barostats/MC_barostat.cu barostats/MC_barostat.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


barostats/Berendsen_barostat.o: barostats/Berendsen_barostat.cu barostats/Berendsen_barostat.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

barostats/andersen_barostat.o: barostats/andersen_barostat.cu barostats/andersen_barostat.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

restrain/restrain.o: restrain/restrain.cu restrain/restrain.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

constrain/constrain.o: constrain/constrain.cu constrain/constrain.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

constrain/SETTLE.o: constrain/SETTLE.cu constrain/SETTLE.cuh constrain/constrain.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

constrain/SHAKE.o: constrain/SHAKE.cu constrain/SHAKE.cuh constrain/constrain.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

constrain/simple_constrain.o: constrain/simple_constrain.cu constrain/simple_constrain.cuh constrain/constrain.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

virtual_atoms/virtual_atoms.o: virtual_atoms/virtual_atoms.cu virtual_atoms/virtual_atoms.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

crd_molecular_map/crd_molecular_map.o: crd_molecular_map/crd_molecular_map.cu crd_molecular_map/crd_molecular_map.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

clean:
	rm -f *.o
	rm -f $(foreach i, $(SUBDIRS), $(i)/*.o)
