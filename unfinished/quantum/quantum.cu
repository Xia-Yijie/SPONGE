#include "quantum.cuh"

static void init_pyscf()
{
        printf("\tqm program: pyscf\n");
        #include "pyscf.py"	
}

static void init_gaussian()
{
	printf("\tqm program: Gaussian\n");
	#include "g09.py"
}

static void init_bdf()
{
	printf("\tqm program: Beijing Density Functional\n");
	#include "bdf.py"
}

void QUANTUM_INFORMATION::Initial(CONTROLLER *controller)
{
	if (controller[0].Command_Exist("qm_program"))
	{
		program = PROGRAM(atoi(controller[0].Command("qm_program")));
		if (program == NONE)
			printf("qm program: not set.\n");
		else if (program == PYSCF)
			init_pyscf();
	        else if (program == GAUSSIAN)
			init_gaussian();
		else if (program == BDF)
			init_bdf();
	}
}
