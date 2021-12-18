#ifndef PYPLUGIN_CUH
#define PYPLUGIN_CUH
#include <Python.h>
#include <numpy/arrayobject.h>
#include "../common.cuh"
#include "../control.cuh"
#include "../MD_core/MD_core.cuh"
#include "../main.cuh"

struct PYTHON_INFORMATION
{
public:
	npy_intp vector_dim[2];
	npy_intp float_dim[1];
        npy_intp temp_dim[1];
	PyArrayObject* pArrayTemp = NULL;

    void Initial();
    void Before_Calculate_Force();
    void Calculate_Force();
    void Calculate_Energy();
    void After_Calculate_Energy();
    void After_Calculate_Force();
    void Iteration_1();
    void Iteration_2();
    void Before_Iteration();
    void After_Iteration();
    void Volume_Change(double factor);
    void Print();
    void Destroy();
};

/*
关于编译变量在windows上的设置：下面是夏义杰的例子
PYTHON=python37
PYTHON_INCLUDE=C:\Users\xyj\AppData\Local\Programs\Python\Python37\include
NUMPY_INCLUDE=C:\Users\xyj\AppData\Local\Programs\Python\Python37\Lib\site-packages\numpy\core\include
PYTHON_LIBRARY=C:\Users\xyj\AppData\Local\Programs\Python\Python37\libs
*/

#endif
