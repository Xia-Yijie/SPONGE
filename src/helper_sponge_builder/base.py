import sys
sys.dont_write_bytecode = True 

#À©³äÄ£¿é
from helper_sponge_builder import html


class Module:
    def __init__(self, dic):
        for i, j in dic.items():
            setattr(self, '_'.join(i.split()), j)
        if not (hasattr(self, "Name") and hasattr(self, "Author") and hasattr(self, "Category")):
            raise Exception("Attribute: Name, Author and Category needed")
        if hasattr(self,"HTML"):
            self.HTML = html.HTML_INFO(self.HTML)
    def __repr__(self):
        return "<Module '{}' written by {}>".format(self.Name, self.Author)
    def __eq__(self, other):
        if (type(other) == type("")):
            return hasattr(self, "Name") and self.Name == other
        if type(self) != type(other):
            return False
        if not (hasattr(self, "Name") and hasattr(other, "Name")):
            return False
        return self.Name == other.Name
        
class Package:
    def __init__(self, dic, modules):
        for i, j in dic.items():
            setattr(self, '_'.join(i.split()), j)
            if i == "Modules":
                module_names = self.Modules.split()
                self.Modules = []
                for module_name in module_names:
                    notfound = 1
                    for module in modules:
                        if module.Name == module_name:
                            self.Modules.append(module)
                            notfound = 0
                            break
                    if notfound:
                        raise Exception("Module Name Not Found", module_name)
                            
    def __repr__(self):
        towrite = "<<Package '{0}'\n Including Modules:\n".format(self.Name)
        for module in self.Modules:
            towrite +="{}\n".format(repr(module))
        return towrite+"Combined by {0}.>>".format(self.Author)
    
    def __eq__(self, other):
        if (type(other) == type("")):
            return hasattr(self, "Name") and self.Name == other
        if type(self) != type(other):
            return False
        if not (hasattr(self, "Name") and hasattr(other, "Name")):
            return False
        return self.Name == other.Name
    
    def Generate_Main_Head_File(self):
        towrite = """#ifndef SMALL_MAIN_CUH
#define SMALL_MAIN_CUH
\n\n"""
        for module in self.Modules:
            if hasattr(module, "Main_Include"):
                towrite += module.Main_Include
                towrite += "\n"
        towrite += """\n\n
void Main_Initial(int argc, char *argv[]);
void Main_Destroy();

void Main_Iteration_1();
void Main_Iteration_2();
void Main_Before_Iteration();
void Main_After_Iteration();

void Main_Before_Calculate_Force();
void Main_Calculate_Force();
void Main_After_Calculate_Force();

void Main_Calculate_Energy();
void Main_After_Calculate_Energy();
void Main_Print();

void Main_Volume_Change(double factor);
void Main_Volume_Change_Largely();

#endif //SMALL_MAIN_CUH"""
        return towrite
    def _Generate_Main_Declaration(self):
        towrite = """#include "main.cuh"\n"""
        for module in self.Modules:
            if hasattr(module, "Main_Declaration"):
                towrite += module.Main_Declaration
                towrite += "\n"
        return towrite
    def _Generate_Main_MD(self):
        return """
int main(int argc, char * argv[])
{
    Main_Initial(argc, argv);
    for (md_info.steps = 1; md_info.steps <= md_info.step_limit; md_info.steps++)
    {
        Main_Iteration_1();
        
        Main_Before_Calculate_Force();
        Main_Calculate_Force();
        Main_After_Calculate_Force();
        
        if ((md_info.steps % md_info.ntwx) == 0)
        {
            Main_Calculate_Energy();
            Main_After_Calculate_Energy();
            Main_Print();
        }
        
        Main_Iteration_2();
    }
    Main_Destroy();
    return 0;
}\n
"""
    def _Generate_Some_Main_Part(self, func_name, attribute_name=None, begin="", end=""):
        towrite = "void %s\n{\n"%(func_name)
        towrite += begin
        towrite += "\n"
        if not attribute_name:
            attribute_name = func_name.split("(")[0]
        priority_name = attribute_name + "_Priority"
        contents_and_priority = []
        for module in self.Modules:
            if hasattr(module, attribute_name) and getattr(module, attribute_name).strip():
                content = getattr(module, attribute_name).rstrip().split("\n")
                
                if hasattr(module, priority_name) and getattr(module, priority_name).strip():
                    priority = list(map(float, getattr(module, priority_name).rstrip().split("\n")))
                else:
                    priority = [0.0]
                
                if len(priority) == 1:
                    priority = [ priority[0] for i in range(len(content))]
                else:
                    priority = [ priority[i] for i in range(len(priority))] + [priority[-1] for i in range(len(priority), len(content))]    
                contents_and_priority.extend(list(zip(content, priority)))
        
        contents_and_priority.sort(key = lambda x: -x[1])
        if contents_and_priority:
            towrite += "\n".join(list(zip(*contents_and_priority))[0])
            towrite += "\n"
        if end:
            towrite += end + "\n"
        towrite += "\n}\n"
        return towrite
        
    def _Generate_Main_Initial(self):
        func_name = "Main_Initial(int argc, char *argv[])"
        attribute_name = "Main_Initial"
        return self._Generate_Some_Main_Part(func_name, attribute_name)
        
    def _Generate_Main_Before_Calculate_Force(self):
        return self._Generate_Some_Main_Part("Main_Before_Calculate_Force()")
        
    def _Generate_Main_Calculate_Force(self):
        return self._Generate_Some_Main_Part("Main_Calculate_Force()")

    def _Generate_Main_Calculate_Energy(self):
        return self._Generate_Some_Main_Part("Main_Calculate_Energy()")

    def _Generate_Main_After_Calculate_Energy(self):
        return self._Generate_Some_Main_Part("Main_After_Calculate_Energy()")

    def _Generate_Main_After_Calculate_Force(self):
        return self._Generate_Some_Main_Part("Main_After_Calculate_Force()")
    def _Generate_Main_Iteration_1(self):
        return self._Generate_Some_Main_Part("Main_Iteration_1()")        
    def _Generate_Main_Iteration_2(self):
        return self._Generate_Some_Main_Part("Main_Iteration_2()")
    def _Generate_Main_Before_Iteration(self):
        return self._Generate_Some_Main_Part("Main_Before_Iteration()")
    def _Generate_Main_After_Iteration(self):
        return self._Generate_Some_Main_Part("Main_After_Iteration()")
    def _Generate_Main_Volume_Change(self):
        return self._Generate_Some_Main_Part("Main_Volume_Change(double factor)")
    def _Generate_Main_Volume_Change_Largely(self):
        end = """	printf("Volume changed too largely. Please restart.\\n");
	Main_Destroy();"""
        return self._Generate_Some_Main_Part("Main_Volume_Change_Largely()", end = end)
    def _Generate_Main_Destroy(self):
        return self._Generate_Some_Main_Part("Main_Destroy()")
    def _Generate_Main_Print(self):
        return self._Generate_Some_Main_Part("Main_Print()")
    
    def Generate_Main_File(self):
        towrite = ""
        towrite += self._Generate_Main_Declaration()
        towrite += self._Generate_Main_MD()
        towrite += self._Generate_Main_Initial()
        towrite += self._Generate_Main_Before_Calculate_Force()
        towrite += self._Generate_Main_Calculate_Force()
        towrite += self._Generate_Main_After_Calculate_Force()
        towrite += self._Generate_Main_Iteration_1()
        towrite += self._Generate_Main_Iteration_2()
        towrite += self._Generate_Main_Before_Iteration()
        towrite += self._Generate_Main_After_Iteration()
        towrite += self._Generate_Main_Volume_Change()
        towrite += self._Generate_Main_Volume_Change_Largely()
        towrite += self._Generate_Main_Calculate_Energy()
        towrite += self._Generate_Main_After_Calculate_Energy()
        towrite += self._Generate_Main_Print()
        towrite += self._Generate_Main_Destroy()
        return towrite
    
    def Generate_Makefile(self):
        towrite = "\n"
        for module in self.Modules:
            if hasattr(module, "Make_Variables") and getattr(module, "Make_Variables").strip():
                lines = getattr(module, "Make_Variables").split("\n")
                for line in lines:
                    line = line.strip()
                    if line:
                        towrite += line
                        towrite += "\n"
        towrite += "\nSUBDIRS=$(shell ls -l | grep ^d | awk '{print $$9}')\ncompiler=nvcc\nCOMMON_COMMAND=-arch=sm_50 -rdc=true -lcudadevrt -lcufft --use_fast_math -O4 -std=c++11"
        for module in self.Modules:
            if hasattr(module, "Make_Command") and getattr(module, "Make_Command").strip():
                line = getattr(module, "Make_Command").strip()
                towrite += " " + line
        towrite += "\nBIN_NAME=SPONGE\n"
        
        towrite += "\ninstall: main.o"
        for module in self.Modules:
            towrite += " " + module.Module_File.replace(".cu", ".o")
        towrite += "\n\t$(compiler) -o $(BIN_NAME)  $^ $(COMMON_COMMAND)\n\tmv $(BIN_NAME) ../bin/$(BIN_NAME)\n\n"
        
        towrite += "\nall: main.o"
        for module in self.Modules:
            towrite += " " + module.Module_File.replace(".cu", ".o")
        towrite += "\n\t\n\n"
        
        towrite += "\nmain.o: main.cu main.cuh"
        towrite += "\n\t$(compiler) -o $@ -c $< $(COMMON_COMMAND)\n\n"
        
        for module in self.Modules:
            towrite += "\n%s: %s"%(module.Module_File.replace(".cu", ".o"), module.Module_File)
            if hasattr(module, "Module_Include") and getattr(module, "Module_Include").strip():
                lines = getattr(module, "Module_Include").split("\n")
                for line in lines:
                    line = line.strip()
                    if line:
                        towrite += " " + line
            towrite += "\n\t$(compiler) -o $@ -c $< $(COMMON_COMMAND)\n\n"
        
        towrite += """clean:
	rm -f *.o
	rm -f $(foreach i, $(SUBDIRS), $(i)/*.o)"""
        return towrite
    
    def Generate_HTML(self):
        return html.Generate_Html(self.Modules)
    

                        
