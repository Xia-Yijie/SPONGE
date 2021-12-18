#ifndef QUANTUM_CUH
#define QUANTUM_CUH
#include <Python.h>
#include "../control.cuh"
#include "../common.cuh"

struct QUANTUM_INFORMATION
{
    int is_init = 0;
    enum PROGRAM
    {
       NONE = 0,
       PYSCF = 1,
       GAUSSIAN = 2,
       BDF = 3,
    } program;
    void Initial(CONTROLLER *controller);
};

#endif
