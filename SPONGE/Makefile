

SUBDIRS=$(shell ls -l | grep ^d | awk '{print $$9}')
compiler=nvcc
COMMON_COMMAND=-arch=sm_50 -rdc=true -lcudadevrt -lcufft --use_fast_math -O4 -std=c++11
BIN_NAME=SPONGE

install: main.o common.o control.o MD_core/MD_core.o bond/bond.o angle/angle.o dihedral/dihedral.o nb14/nb14.o neighbor_list/neighbor_list.o Lennard_Jones_force/Lennard_Jones_force.o PME_force/PME_force.o Langevin_MD/Langevin_MD.o Langevin_MD/Middle_Langevin_MD.o MC_barostat/MC_barostat.o Berendsen_barostat/Berendsen_barostat.o restrain/restrain.o simple_constrain/simple_constrain.o virtual_atoms/virtual_atoms.o crd_molecular_map/crd_molecular_map.o
	$(compiler) -o $(BIN_NAME)  $^ $(COMMON_COMMAND)


all: main.o common.o control.o MD_core/MD_core.o bond/bond.o angle/angle.o dihedral/dihedral.o nb14/nb14.o neighbor_list/neighbor_list.o Lennard_Jones_force/Lennard_Jones_force.o PME_force/PME_force.o Langevin_MD/Langevin_MD.o Langevin_MD/Middle_Langevin_MD.o MC_barostat/MC_barostat.o Berendsen_barostat/Berendsen_barostat.o restrain/restrain.o simple_constrain/simple_constrain.o virtual_atoms/virtual_atoms.o crd_molecular_map/crd_molecular_map.o
	


main.o: main.cu main.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


common.o: common.cu common.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


control.o: control.cu common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


MD_core/MD_core.o: MD_core/MD_core.cu MD_core/MD_core.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


bond/bond.o: bond/bond.cu bond/bond.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


angle/angle.o: angle/angle.cu angle/angle.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


dihedral/dihedral.o: dihedral/dihedral.cu dihedral/dihedral.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


nb14/nb14.o: nb14/nb14.cu nb14/nb14.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


neighbor_list/neighbor_list.o: neighbor_list/neighbor_list.cu neighbor_list/neighbor_list.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


Lennard_Jones_force/Lennard_Jones_force.o: Lennard_Jones_force/Lennard_Jones_force.cu Lennard_Jones_force/Lennard_Jones_force.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


PME_force/PME_force.o: PME_force/PME_force.cu PME_force/PME_force.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


Langevin_MD/Langevin_MD.o: Langevin_MD/Langevin_MD.cu Langevin_MD/Langevin_MD.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


Langevin_MD/Middle_Langevin_MD.o: Langevin_MD/Middle_Langevin_MD.cu Langevin_MD/Middle_Langevin_MD.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


MC_barostat/MC_barostat.o: MC_barostat/MC_barostat.cu MC_barostat/MC_barostat.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


Berendsen_barostat/Berendsen_barostat.o: Berendsen_barostat/Berendsen_barostat.cu Berendsen_barostat/Berendsen_barostat.cuh common.cuh control.cuh main.cuh MD_core/MD_core.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


restrain/restrain.o: restrain/restrain.cu restrain/restrain.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


simple_constrain/simple_constrain.o: simple_constrain/simple_constrain.cu simple_constrain/simple_constrain.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)


virtual_atoms/virtual_atoms.o: virtual_atoms/virtual_atoms.cu virtual_atoms/virtual_atoms.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

crd_molecular_map/crd_molecular_map.o: crd_molecular_map/crd_molecular_map.cu crd_molecular_map/crd_molecular_map.cuh common.cuh control.cuh
	$(compiler) -o $@ -c $< $(COMMON_COMMAND)

clean:
	rm -f *.o
	rm -f $(foreach i, $(SUBDIRS), $(i)/*.o)
