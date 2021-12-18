#include "pyplugin.cuh"


extern PYTHON_INFORMATION py_info;
extern MD_INFORMATION md_info;
extern CONTROLLER controller;

//the construction of Python module "SPONGE" 
//sponge.get_atom_numbers()
//return int atom_numbers
static PyObject* Sponge_Get_Atom_Numbers(PyObject* self, PyObject*args)
{
	if (!PyArg_ParseTuple(args,""))
		return NULL;
	return Py_BuildValue("l", md_info.atom_numbers);
}


//sponge.get_step()
//return int step
static PyObject* Sponge_Get_Step(PyObject* self, PyObject*args)
{
	if (!PyArg_ParseTuple(args,""))
		return NULL;
	return Py_BuildValue("l", md_info.steps);
}

//sponge.get_potential_energy()
//return float potential_energy
static PyObject* Sponge_Get_Potential_Energy(PyObject* self, PyObject*args)
{
	if (!PyArg_ParseTuple(args, ""))
		return NULL;
	return Py_BuildValue("f", md_info.total_potential_energy);
}

//sponge.set_potential_energy(float energy)
//return None
static PyObject* Sponge_Set_Potential_Energy(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"energy", NULL};
	if (!PyArg_ParseTupleAndKeywords(args, kw, "f", kwlist, &md_info.total_potential_energy))
		return NULL;
	return Py_BuildValue("");
}
//sponge.get_box_length()
//return numpy.ndarray box_length
static PyObject* Sponge_Get_Box_Length(PyObject* self, PyObject*args)
{
	if (!PyArg_ParseTuple(args,""))
		return NULL;
	return PyArray_SimpleNewFromData(1, py_info.temp_dim, NPY_FLOAT32, &md_info.box_length);
}

//sponge.get_mass(int start = 0, int length = atom_numbers, bool cuda_memcpy = False)
//return numpy.ndarray mass
static PyObject* Sponge_Get_Mass(PyObject* self, PyObject* args, PyObject* kw)
{
	static char *kwlist[] = {"start","length","cuda_memcpy", NULL};
	int length = md_info.atom_numbers;
	int start = 0;
	int cuda_memcpy = 0;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "|iii",kwlist, &start, &length, &cuda_memcpy))
		return NULL;
	py_info.float_dim[0] =length;
	if (cuda_memcpy)
		cudaMemcpy(md_info.h_mass + start, md_info.d_mass + start, sizeof(float) * length, cudaMemcpyDeviceToHost);
	return PyArray_SimpleNewFromData(1, py_info.float_dim, NPY_FLOAT32, md_info.h_mass + start);
}



//sponge.get_charge(int start = 0, int length = atom_numbers, bool cuda_memcpy = False)
//return numpy.ndarray charge
static PyObject* Sponge_Get_Charge(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"start","length","cuda_memcpy", NULL};
	int length = md_info.atom_numbers;
	int start = 0;
	int cuda_memcpy = 0;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "|iii",kwlist, &start, &length, &cuda_memcpy))
		return NULL;
	py_info.float_dim[0] = length;
	if (cuda_memcpy)
		cudaMemcpy(md_info.h_charge + start, md_info.d_charge + start, sizeof(float)*length, cudaMemcpyDeviceToHost);
	return PyArray_SimpleNewFromData(1, py_info.float_dim, NPY_FLOAT32, md_info.h_charge + start);
}

//sponge.set_charge(numpy.ndarray charge, int start = 0, int length = atom_numbers, bool cuda_memcpy = True)
//return None
static PyObject* Sponge_Set_Charge(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"charge","start","length","cuda_memcpy", NULL};
	int length = md_info.atom_numbers;
	int start = 0;
	int cuda_memcpy = 1;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "O!|iii", kwlist, &PyArray_Type, &py_info.pArrayTemp, &start, &length, &cuda_memcpy))
		return NULL;
	if (py_info.pArrayTemp->strides[1] != 4 || !PyArray_ISFLOAT(py_info.pArrayTemp))
	{
		PyErr_SetString(PyExc_Exception, "float32 needed");		
		return NULL;
	}
	if (py_info.pArrayTemp->nd != 1 || py_info.pArrayTemp->dimensions[0] != length)
	{
		PyErr_SetString(PyExc_Exception, "array shape is not right");
		return NULL;
	}
	cudaMemcpy(md_info.h_charge + start, py_info.pArrayTemp->data, sizeof(float)*length, cudaMemcpyHostToHost);
	if (cuda_memcpy)
		cudaMemcpy(md_info.d_charge + start, md_info.h_charge + start, sizeof(float)*length, cudaMemcpyHostToDevice);
	return Py_BuildValue("");
}

//sponge.get_coordinate(int start = 0, int length = atom_numbers, bool cuda_memcpy = True)
//return numpy.ndarray coordinate
static PyObject* Sponge_Get_Coordinate(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"start","length","cuda_memcpy", NULL};
	int length = md_info.atom_numbers;
	int start = 0;
	int cuda_memcpy = 1;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "|iii",kwlist, &start, &length, &cuda_memcpy))
		return NULL;
	py_info.vector_dim[0] = length;
	if (cuda_memcpy)
		cudaMemcpy(md_info.coordinate + start, md_info.crd + start, sizeof(VECTOR)*length, cudaMemcpyDeviceToHost);
	return PyArray_SimpleNewFromData(2, py_info.vector_dim, NPY_FLOAT32, md_info.coordinate + start);
}

//sponge.set_coordinate(numpy.array coordinate, int start = 0, int length = atom_numbers, bool cuda_memcpy = True)
//return None
static PyObject* Sponge_Set_Coordinate(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"coordinate","start","length","cuda_memcpy", NULL};
	int length = md_info.atom_numbers;
	int start = 0;
	int cuda_memcpy = 1;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "O!|iii", kwlist, &PyArray_Type, &py_info.pArrayTemp, &start, &length, &cuda_memcpy))
		return NULL;
	py_info.vector_dim[0] = length;
	if (py_info.pArrayTemp->strides[1] != 4 || !PyArray_ISFLOAT(py_info.pArrayTemp))
	{
		PyErr_SetString(PyExc_Exception, "float32 needed");		
		return NULL;
	}
	if (py_info.pArrayTemp->nd != 2 || py_info.pArrayTemp->dimensions[0] != length || py_info.pArrayTemp->dimensions[1] != 3)
	{
		PyErr_SetString(PyExc_Exception, "array shape is not right");
		return NULL;
	}
	cudaMemcpy(md_info.coordinate + start, py_info.pArrayTemp->data, sizeof(VECTOR)*length, cudaMemcpyHostToHost);
	if (cuda_memcpy)
		cudaMemcpy(md_info.crd + start, md_info.coordinate + start, sizeof(VECTOR)* length, cudaMemcpyHostToDevice);
	return Py_BuildValue("");
}

//sponge.get_force(int start = 0, int length = atom_numbers, bool cuda_memcpy = True)
//return numpy.ndarray force
static PyObject* Sponge_Get_Force(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"start","length","cuda_memcpy", NULL};
	int length = md_info.atom_numbers;
	int start = 0;
	int cuda_memcpy = 1;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "|iii",kwlist, &start, &length, &cuda_memcpy))
		return NULL;
	py_info.vector_dim[0] = length;
	if (cuda_memcpy)
		cudaMemcpy(md_info.force + start, md_info.frc + start, sizeof(VECTOR)*length, cudaMemcpyDeviceToHost);
	return PyArray_SimpleNewFromData(2, py_info.vector_dim, NPY_FLOAT32, md_info.force + start);
}

//sponge.set_force(numpy.ndarray force, int start = 0, int length = atom_numbers, bool cuda_memcpy = True)
//return None
static PyObject* Sponge_Set_Force(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"force","start","length","cuda_memcpy", NULL};
	int length = md_info.atom_numbers;
	int start = 0;
	int cuda_memcpy = 1;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "O!|iii", kwlist, &PyArray_Type, &py_info.pArrayTemp, &start, &length, &cuda_memcpy))
		return NULL;
	py_info.vector_dim[0] = length;
	if (py_info.pArrayTemp->strides[1] != 4 || !PyArray_ISFLOAT(py_info.pArrayTemp))
	{
		PyErr_SetString(PyExc_Exception, "float32 needed");		
		return NULL;
	}
	if (py_info.pArrayTemp->nd != 2 || py_info.pArrayTemp->dimensions[0] != length || py_info.pArrayTemp->dimensions[1] != 3)
	{
		PyErr_SetString(PyExc_Exception, "array shape is not right");
		return NULL;
	}
	cudaMemcpy(md_info.force + start, py_info.pArrayTemp->data, sizeof(VECTOR)*length, cudaMemcpyHostToHost);
	if (cuda_memcpy)
		cudaMemcpy(md_info.frc + start, md_info.force + start, sizeof(VECTOR)* length, cudaMemcpyHostToDevice);
	return Py_BuildValue("");
}

//sponge.get_label(int start = 0, int length = atom_numbers)
//return list label
static PyObject* Sponge_Get_Label(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"start", "length", NULL};
	int length = md_info.atom_numbers;
	int start = 0;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "|ii",kwlist, &start, &length))
		return NULL;
	PyObject* label = PyList_New((Py_ssize_t)length);
	for (Py_ssize_t i = 0; i < length; i++)
		PyList_SetItem(label, i, Py_BuildValue("s", md_info.atom_labels[start + i]));
	return label;
}

//sponge.set_label(list label, int start = 0, int length = atom_numbers)
//return None
static PyObject* Sponge_Set_Label(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"label", "start", "length", NULL};
	PyObject* label;
	int length = md_info.atom_numbers;
	int start = 0;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "O!|ii",kwlist, &PyList_Type, &label, &start, &length))
		return NULL;
	PyObject* temp = PyTuple_New(1);
	char *buffer;
	char *temp2;
	for (Py_ssize_t i = 0; i < length; i++)
	{
		PyTuple_SetItem(temp, 0, PyList_GetItem(label, i));
		PyArg_ParseTuple(temp, "s", &buffer);
		temp2 = md_info.atom_labels[start+i];
		strcpy(temp2, buffer);
		
	}	
	return Py_BuildValue("");
}

//sponge.command_exist(str command)
//return bool existance
static PyObject* Sponge_Command_Exist(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"command", NULL};
	char *buffer;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "s",kwlist, &buffer))
		return NULL;
	if (controller.Command_Exist(buffer))
		Py_RETURN_TRUE;
	else
		Py_RETURN_FALSE;
}

//sponge.command(str command)
//return str value
static PyObject* Sponge_Command(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"command", NULL};
	char *buffer;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "s",kwlist, &buffer))
		return NULL;
	return Py_BuildValue("s", controller.Command(buffer));
}

//sponge.original_command(str command)
//return str value
static PyObject* Sponge_Original_Command(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"command", NULL};
	char *buffer;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "s",kwlist, &buffer))
		return NULL;
	return Py_BuildValue("s", controller.Original_Command(buffer));
}

//sponge.need_atom_energy()
//return bool need
static PyObject* Sponge_Need_Atom_Energy(PyObject* self, PyObject* args)
{
	if (!PyArg_ParseTuple(args, ""))
		return NULL;
	if (md_info.need_atom_energy)
		Py_RETURN_TRUE;
	else
		Py_RETURN_FALSE;
}

//sponge.need_virial()
//return bool need
static PyObject* Sponge_Need_Virial(PyObject* self, PyObject* args)
{
	if (!PyArg_ParseTuple(args, ""))
		return NULL;
	if (md_info.need_virial)
		Py_RETURN_TRUE;
	else
		Py_RETURN_FALSE;
}

//sponge.get_atom_energy(int start = 0, int length = atom_numbers, bool cuda_memcpy = True)
//return numpy.ndarray atom_energy
static PyObject* Sponge_Get_Atom_Energy(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"start","length","cuda_memcpy", NULL};
	int length = md_info.atom_numbers;
	int start = 0;
	int cuda_memcpy = 1;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "|iii",kwlist, &start, &length, &cuda_memcpy))
		return NULL;
	py_info.float_dim[0] = length;
	if (cuda_memcpy)
		cudaMemcpy(md_info.h_atom_energy + start, md_info.d_atom_energy + start, sizeof(float)*length, cudaMemcpyDeviceToHost);
	return PyArray_SimpleNewFromData(1, py_info.float_dim, NPY_FLOAT32, md_info.h_atom_energy + start);
}

//sponge.set_atom_energy(numpy.ndarray atom_energy, int start = 0, int length = atom_numbers, bool cuda_memcpy = True)
//return None
static PyObject* Sponge_Set_Atom_Energy(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"atom_energy","start","length","cuda_memcpy", NULL};
	int length = md_info.atom_numbers;
	int start = 0;
	int cuda_memcpy = 1;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "O!|iii", kwlist, &PyArray_Type, &py_info.pArrayTemp, &start, &length, &cuda_memcpy))
		return NULL;
	if (py_info.pArrayTemp->strides[1] != 4 || !PyArray_ISFLOAT(py_info.pArrayTemp))
	{
		PyErr_SetString(PyExc_Exception, "float32 needed");		
		return NULL;
	}
	if (py_info.pArrayTemp->nd != 1 || py_info.pArrayTemp->dimensions[0] != length)
	{
		PyErr_SetString(PyExc_Exception, "array shape is not right");
		return NULL;
	}
	cudaMemcpy(md_info.h_atom_energy + start, py_info.pArrayTemp->data, sizeof(float)*length, cudaMemcpyHostToHost);
	if (cuda_memcpy)
		cudaMemcpy(md_info.d_atom_energy + start, md_info.h_atom_energy + start, sizeof(float)* length, cudaMemcpyHostToDevice);
	return Py_BuildValue("");
}


//sponge.get_atom_virial(int start = 0, int length = atom_numbers, bool cuda_memcpy = True)
//return numpy.ndarray atom_virial
static PyObject* Sponge_Get_Atom_Virial(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"start","length","cuda_memcpy", NULL};
	int length = md_info.atom_numbers;
	int start = 0;
	int cuda_memcpy = 1;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "|iii",kwlist, &start, &length, &cuda_memcpy))
		return NULL;
	py_info.float_dim[0] = length;
	if (cuda_memcpy)
		cudaMemcpy(md_info.h_atom_virial + start, md_info.d_atom_virial + start, sizeof(float)*length, cudaMemcpyDeviceToHost);
	return PyArray_SimpleNewFromData(1, py_info.float_dim, NPY_FLOAT32, md_info.h_atom_virial + start);
}

//sponge.set_atom_virial(numpy.ndarray atom_virial, int start = 0, int length = atom_numbers, bool cuda_memcpy = True)
//return None
static PyObject* Sponge_Set_Atom_Virial(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"atom_virial","start","length","cuda_memcpy", NULL};
	int length = md_info.atom_numbers;
	int start = 0;
	int cuda_memcpy = 1;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "O!|iii", kwlist, &PyArray_Type, &py_info.pArrayTemp, &start, &length, &cuda_memcpy))
		return NULL;
	if (py_info.pArrayTemp->strides[1] != 4 || !PyArray_ISFLOAT(py_info.pArrayTemp))
	{
		PyErr_SetString(PyExc_Exception, "float32 needed");		
		return NULL;
	}
	if (py_info.pArrayTemp->nd != 1 || py_info.pArrayTemp->dimensions[0] != length)
	{
		PyErr_SetString(PyExc_Exception, "array shape is not right");
		return NULL;
	}
	cudaMemcpy(md_info.h_atom_virial + start, py_info.pArrayTemp->data, sizeof(float)*length, cudaMemcpyHostToHost);
	if (cuda_memcpy)
		cudaMemcpy(md_info.d_atom_virial + start, md_info.h_atom_virial + start, sizeof(float)* length, cudaMemcpyHostToDevice);
	return Py_BuildValue("");
}

//sponge.get_virial(cuda_memcpy = True)
//return float virial
static PyObject* Sponge_Get_Virial(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"cuda_memcpy", NULL};
	int cuda_memcpy = 1;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "|i", kwlist, &cuda_memcpy))
		return NULL;
	if (cuda_memcpy)
		cudaMemcpy(&md_info.h_virial, md_info.d_virial, sizeof(float), cudaMemcpyDeviceToHost);
	return Py_BuildValue("f", md_info.h_virial);
}

//sponge.set_virial(float virial, cuda_memcpy = True)
//return None
static PyObject* Sponge_Set_Virial(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"virial", "cuda_memcpy", NULL};
	int cuda_memcpy = 1;
	if (!PyArg_ParseTupleAndKeywords(args, kw, "f|i", kwlist, &md_info.h_virial, &cuda_memcpy))
		return NULL;
	if (cuda_memcpy)
		cudaMemcpy(md_info.d_virial, &md_info.h_virial, sizeof(float), cudaMemcpyHostToDevice);

	return Py_BuildValue("");
}

//sponge.add_print_head(str head)
//return None
static PyObject* Sponge_Add_Print_Head(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"head", NULL};
	char *head = "";
	if (!PyArg_ParseTupleAndKeywords(args, kw, "s", kwlist, &head))
		return NULL;
	if (strlen(head) != 0)
	{
		strcat(controller.out_line_format, " ");
		strcat(controller.out_line_format, head);
	}
	return Py_BuildValue("");
}

//sponge.add_print(str content)
//return None
static PyObject* Sponge_Add_Print(PyObject* self, PyObject* args, PyObject* kw)
{
	static char* kwlist[] = {"content", NULL};
	char *content = "";
	if (!PyArg_ParseTupleAndKeywords(args, kw, "s", kwlist, &content))
		return NULL;
	if (strlen(content) != 0)
	{
		strcat(controller.out_line, " ");
		strcat(controller.out_line, content);
	}
	return Py_BuildValue("");
}

static PyMethodDef SpongeMethod[] =
{
  {"get_step", (PyCFunction)Sponge_Get_Step, METH_VARARGS, "get the md step"},
  {"get_atom_numbers", (PyCFunction)Sponge_Get_Atom_Numbers, METH_VARARGS, "get the total atom number"},
  {"get_potential_energy", (PyCFunction)Sponge_Get_Potential_Energy, METH_VARARGS, "get the total potential energy"},
  {"set_potential_energy", (PyCFunction)Sponge_Set_Potential_Energy, METH_VARARGS|METH_KEYWORDS, "set the total potential energy" },
  {"get_box_length", (PyCFunction)Sponge_Get_Box_Length, METH_VARARGS, "get the box length"},
  {"get_mass", (PyCFunction)Sponge_Get_Mass, METH_VARARGS|METH_KEYWORDS, "get the md mass" },
  {"get_charge", (PyCFunction)Sponge_Get_Charge, METH_VARARGS|METH_KEYWORDS, "get the md charge" },
  {"set_charge", (PyCFunction)Sponge_Set_Charge, METH_VARARGS|METH_KEYWORDS, "set the md charge" },
  {"get_coordinate", (PyCFunction)Sponge_Get_Coordinate, METH_VARARGS|METH_KEYWORDS, "get the md coordinate" },
  {"set_coordinate", (PyCFunction)Sponge_Set_Coordinate, METH_VARARGS|METH_KEYWORDS, "set the md coordinate" },
  {"get_force", (PyCFunction)Sponge_Get_Force, METH_VARARGS|METH_KEYWORDS, "get the md force" },
  {"set_force", (PyCFunction)Sponge_Set_Force, METH_VARARGS|METH_KEYWORDS, "set the md force" },
  {"get_label", (PyCFunction)Sponge_Get_Label, METH_VARARGS|METH_KEYWORDS, "get the atom labels" },
  {"set_label", (PyCFunction)Sponge_Set_Label, METH_VARARGS|METH_KEYWORDS, "set the atom labels" },
  {"command_exist", (PyCFunction)Sponge_Command_Exist, METH_VARARGS|METH_KEYWORDS, "check if the command exist" },
  {"command", (PyCFunction)Sponge_Command, METH_VARARGS|METH_KEYWORDS, "get the value of the command"},
  {"original_command", (PyCFunction)Sponge_Original_Command, METH_VARARGS|METH_KEYWORDS, "get the original value of the command"},
  {"need_atom_energy", (PyCFunction)Sponge_Need_Atom_Energy, METH_VARARGS, "check if need atom_energy"},
  {"need_virial", (PyCFunction)Sponge_Need_Virial, METH_VARARGS, "check if need virial"},
  {"get_atom_energy", (PyCFunction)Sponge_Get_Atom_Energy, METH_VARARGS|METH_KEYWORDS, "get the atom energy" },
  {"set_atom_energy", (PyCFunction)Sponge_Set_Atom_Energy, METH_VARARGS|METH_KEYWORDS, "set the atom energy" },
  {"get_atom_virial", (PyCFunction)Sponge_Get_Atom_Virial, METH_VARARGS|METH_KEYWORDS, "get the atom virial" },
  {"set_atom_virial", (PyCFunction)Sponge_Set_Atom_Virial, METH_VARARGS|METH_KEYWORDS, "set the atom virial" },
  {"get_virial", (PyCFunction)Sponge_Get_Virial, METH_VARARGS|METH_KEYWORDS, "get the scaler total virial" },
  {"set_virial", (PyCFunction)Sponge_Set_Virial, METH_VARARGS|METH_KEYWORDS, "set the scaler total virial" },
  {"add_print_head", (PyCFunction)Sponge_Add_Print_Head, METH_VARARGS|METH_KEYWORDS, "add some contents to the print head" },
  {"add_print", (PyCFunction)Sponge_Add_Print, METH_VARARGS|METH_KEYWORDS, "add some contents to the every-step print" },
  {NULL,NULL,0,NULL}
};

static PyModuleDef SpongeModule = 
{
  PyModuleDef_HEAD_INIT, "sponge",NULL,-1,SpongeMethod,
  NULL,NULL,NULL,NULL
};

static PyObject* PyInit_sponge(void)
{
  return PyModule_Create(&SpongeModule);
}
///////////////////////////////

///////numpy init////
int init_numpy()
{
    import_array();
    return 0;
}


//class method
void PYTHON_INFORMATION::Initial()
{
    printf("INITIALIZING PYTHON\n");
    char buffer[256];
    vector_dim[0] = md_info.atom_numbers;
    vector_dim[1] = 3;
    temp_dim[0] = 3;
    float_dim[0] = md_info.atom_numbers;
    
    PyImport_AppendInittab("sponge",&PyInit_sponge);
    Py_Initialize();
    if (!Py_IsInitialized())
    {
        printf("\tPython Init Failed.\n");
    }
	else
	{
		PyRun_SimpleString("print('\tPython Initialized')");
	}
    wchar_t *temp_args[1] = {L"SPONGE"};
    PySys_SetArgv(1,temp_args); 
    init_numpy();
    //初始化python系统文件路径，保证可以访问到 .py文件  
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.dont_write_bytecode = True");
    PyRun_SimpleString("import os");
    PyRun_SimpleString("import sponge");
 

    PyRun_SimpleString("sponge.plugin_funcs={}");
    PyRun_SimpleString("sponge.plugin_funcs['Before_Calculate_Force']=[]");
    PyRun_SimpleString("sponge.plugin_funcs['Calculate_Force']=[]");
    PyRun_SimpleString("sponge.plugin_funcs['Calculate_Energy']=[]");
    PyRun_SimpleString("sponge.plugin_funcs['After_Calculate_Energy']=[]");
    PyRun_SimpleString("sponge.plugin_funcs['After_Calculate_Force']=[]");
    PyRun_SimpleString("sponge.plugin_funcs['Before_Iteration']=[]");
    PyRun_SimpleString("sponge.plugin_funcs['After_Iteration']=[]");
    PyRun_SimpleString("sponge.plugin_funcs['Iteration_1']=[]");
    PyRun_SimpleString("sponge.plugin_funcs['Iteration_2']=[]");
    PyRun_SimpleString("sponge.plugin_funcs['Volume_Change']=[]");
    PyRun_SimpleString("sponge.plugin_funcs['Print']=[]");
    PyRun_SimpleString("sponge.plugin_funcs['Destroy']=[]");
    PyRun_SimpleString("def register(place):\n"\
                       "  def decorator(func):\n"\
                       "     sponge.plugin_funcs[place].append(func)\n"\
                       "     return func\n"\
                       "  return decorator\n");
    PyRun_SimpleString("sponge.register = register");

    if (controller.Command_Exist("py"))
    {
        sprintf(buffer, "sponge.fname = '%s'", controller.Command("py")); 
        PyRun_SimpleString(buffer);  
        PyRun_SimpleString( "sys.path.append( os.path.dirname( os.path.abspath( sponge.fname ) ) )" );
        PyRun_SimpleString("sponge_pyplugin = __import__(os.path.splitext(os.path.basename(sponge.fname))[0])");
        PyRun_SimpleString("if 'sponge_pyplugin' not in dir():\n"\
		               "    print('\tImport Module Failed.')\n"\
					   "    input()\n"\
					   "    exit(1)\n"\
					   "else:\n"\
					   "    print('\tModule Imported.')\n"
					   );
    }
    printf("END (INITIALIZING PYTHON)\n");

}

void PYTHON_INFORMATION::Before_Calculate_Force()
{
    PyRun_SimpleString(R"XYJ(
try:
  for func in sponge.plugin_funcs['Before_Calculate_Force']:
    func()
except:
  import traceback
  traceback.print_exc()
  exit(1)
)XYJ");
}


void PYTHON_INFORMATION::Calculate_Force()
{
    PyRun_SimpleString(R"XYJ(
try:
  for func in sponge.plugin_funcs['Calculate_Force']:
    func()
except:
  import traceback
  traceback.print_exc()
  exit(1)
)XYJ");
}

void PYTHON_INFORMATION::Calculate_Energy()
{
    PyRun_SimpleString(R"XYJ(
try:
  for func in sponge.plugin_funcs['Calculate_Energy']:
    func()
except:
  import traceback
  traceback.print_exc()
  exit(1)
)XYJ");
}

void PYTHON_INFORMATION::After_Calculate_Energy()
{
    PyRun_SimpleString(R"XYJ(
try:
  for func in sponge.plugin_funcs['After_Calculate_Energy']:
    func()
except:
  import traceback
  traceback.print_exc()
  exit(1)
)XYJ");

}

void PYTHON_INFORMATION::After_Calculate_Force()
{
    PyRun_SimpleString(R"XYJ(
try:
  for func in sponge.plugin_funcs['After_Calculate_Force']:
    func()
except:
  import traceback
  traceback.print_exc()
  exit(1)
)XYJ");
}

void PYTHON_INFORMATION::Iteration_1()
{
    PyRun_SimpleString(R"XYJ(
try:
  for func in sponge.plugin_funcs['Iteration_1']:
    func()
except:
  import traceback
  traceback.print_exc()
  exit(1)
)XYJ");
}

void PYTHON_INFORMATION::Iteration_2()
{
    PyRun_SimpleString(R"XYJ(
try:
  for func in sponge.plugin_funcs['Iteration_2']:
    func()
except:
  import traceback
  traceback.print_exc()
  exit(1)
)XYJ");
}

void PYTHON_INFORMATION::Before_Iteration()
{
    PyRun_SimpleString(R"XYJ(
try:
  for func in sponge.plugin_funcs['Before_Iteration']:
    func()
except:
  import traceback
  traceback.print_exc()
  exit(1)
)XYJ");
}

void PYTHON_INFORMATION::After_Iteration()
{
    PyRun_SimpleString(R"XYJ(
try:
  for func in sponge.plugin_funcs['After_Iteration']:
    func()
except:
  import traceback
  traceback.print_exc()
  exit(1)
)XYJ");
}

void PYTHON_INFORMATION::Volume_Change(double factor)
{
    char buffer[256];
    //sprintf(buffer, "try:\n  for func in sponge.plugin_funcs['Volume_Change']:\n    func(%lf)\nexcept:\n  exit(1)", factor);
	sprintf(buffer, R"XYJ(
try:
  for func in sponge.plugin_funcs['Volume_Change']:
    func(%lf)
except:
  import traceback
  traceback.print_exc()
  exit(1))XYJ", factor);
    PyRun_SimpleString(buffer);
}

void PYTHON_INFORMATION::Print()
{
    PyRun_SimpleString(R"XYJ(
try:
  for func in sponge.plugin_funcs['Print']:
    func()
except:
  import traceback
  traceback.print_exc()
  exit(1)
)XYJ");
}

void PYTHON_INFORMATION::Destroy()
{
    PyRun_SimpleString(R"XYJ(
try:
  for func in sponge.plugin_funcs['Destroy']:
    func()
except:
  import traceback
  traceback.print_exc()
  exit(1)
)XYJ");
    Py_Finalize();
}
