import sys
sys.dont_write_bytecode = True

import os
from helper_sponge_builder import base
import PySimpleGUI as sg
import subprocess
import threading
import time
import traceback
import shutil

def parse_smf(string):
    dic = {}

    temp = string.split(":::")
    if temp[0] != 'SPONGE Module File\n':
        return
    
    for i in range(1,len(temp),2):
        if temp[i+1]:
            dic[temp[i]] = temp[i+1][1:]
        else:
            dic[temp[i]] = temp[i+1]
    
    for i in ["Description"]:
        if i in dic:
            dic[i] = dic[i].rstrip()
            
    for i in ["Name", "Author", "Progress", "Module File", "Main Include", "Main Declaration", "Category"]:
        if i in dic:
            dic[i] = dic[i].strip()
    return dic

def parse_spf(string):
    dic = {}

    temp = string.split(":::")
    if temp[0] != 'SPONGE Package File\n':
        return
    
    for i in range(1,len(temp),2):
        if temp[i+1]:
            dic[temp[i]] = temp[i+1][1:]
        else:
            dic[temp[i]] = temp[i+1]
    for i in ["Description"]:
        if i in dic:
            dic[i] = dic[i].rstrip()
    for i in ["Name", "Author", "Progress"]:
        dic[i] = dic[i].strip()
    
    return dic


def Search_Modules(path):
    Modules = []
    for root, dirs, files in os.walk(path):
        for f in files:
            if f.endswith(".smf"):
                #print(f)
                fr = open(os.path.join(root,f), encoding="utf-8")
                
                temp = fr.read()
                fr.close()
                temp = parse_smf(temp)
                if temp:
                    if temp["Name"] not in Modules:
                        Modules.append(base.Module(temp))
                    else:
                        raise Exception("Modules have the same name: ", temp["Name"])
    return Modules

def Search_Packages(path, modules):
    Packages = []
    for root, dirs, files in os.walk(path):
        for f in files:
            if f.endswith(".spf"):
                fr = open(os.path.join(root,f), encoding="utf-8")
                temp = fr.read()
                fr.close()
                temp = parse_spf(temp)
                if temp:
                    Packages.append(base.Package(temp, modules))
    return Packages

def Update_Graph(gui,content):
    if content:
        gui.Graph.update(value = "\n".join(content)+"\n")
    else:
        gui.Graph.update(value = "")
     

def Check_GUI_Graph(value, gui):
    content = value.split("\n")
    while '' in content:
        content.remove('')
    while "md_core" in content:
        content.remove("md_core")
    content.insert(0,"md_core")
    while "control" in content:
        content.remove("control")
    content.insert(0,"control")
    while "common" in content:
        content.remove("common")
    content.insert(0,"common")
    for i in content:
        if i not in gui.Modules:
            gui.stderr.print("警告："+i+"不在支持的Modules之中，已自动剔除")
            content.remove(i)
    for mi, Module in enumerate(gui.Modules):
        if Module in content:
            gui.Module_Buttons[mi].update(value = True)
        else:
            gui.Module_Buttons[mi].update(value = False)
    Update_Graph(gui,content)
        


threadWorking = []
hint = ""
def Module_compile(module, command):
    global threadWorking, hint
    threadWorking.append(module.Name)
    try:
        obj = module.Module_File.replace(".cu",".o")
        need_compile = 0
        if os.path.exists(obj):
            selftime = os.path.getmtime(obj)
            if os.path.getmtime(obj.replace(".o",".cu")) > selftime:
                need_compile = 1
            else:
                if hasattr(module, "Module_Include") and getattr(module, "Module_Include").strip():
                    lines = getattr(module, "Module_Include").split("\n")
                    for line in lines:
                        line = line.strip()
                        if line and os.path.getmtime(line) > selftime:
                            need_compile = 1
                            break
        else:
            need_compile = 1
        if need_compile:
            run_command = "nvcc -o {} -c {} {}".format(obj, module.Module_File, command)
            hint += run_command + "\n"
            out = subprocess.run(run_command, universal_newlines=True, capture_output=True)
            hint += out.stdout + "\n"
            hint += out.stderr + "\n"
    except:
        hint += run_command + "\n"
    finally:
        threadWorking.remove(module.Name)

def Obj_compile(obj, command, gui, *includes):
    global threadWorking, hint
    threadkey = obj[:-2]
    threadWorking.append(threadkey)
    try:
        need_compile = 0
        if os.path.exists(obj):
            selftime = os.path.getmtime(obj)
            if os.path.getmtime(obj.replace(".o",".cu")) > selftime:
                need_compile = 1
            else:
                for line in includes:
                    if line and os.path.getmtime(line) > selftime:
                        need_compile = 1
                        break
        else:
            need_compile = 1
        if need_compile:
            run_command = "nvcc -o {} -c {} {}".format(obj, obj.replace(".o", ".cu"), command)
            hint += run_command + "\n"
            out = subprocess.run(run_command, shell = True, universal_newlines=True, capture_output=True)
            hint += out.stdout + "\n"
            hint += out.stderr + "\n"
    except:
        hint += run_command + "\n"
    finally:
        threadWorking.remove(threadkey)
 

def Link_compile(objs, command):
    global threadWorking, hint
    threadWorking.append("Linking")
    try:
        run_command = "nvcc -o {} {} {}".format("SPONGE", objs, command)
        hint += run_command + "\n"
        out = subprocess.run(run_command, shell = True, universal_newlines=True, capture_output=True)
        
        hint += out.stdout + "\n"
        hint += out.stderr + "\n"
    except:
        hint += run_command + "\n"
    finally:
        threadWorking.remove("Linking")
    
def Ask_For_Complete(**kw):
    layout = []
    default_text = ""
    for key, value in kw.items():
        default_text += key + "=" + value + "\n"
    layout.append([sg.Multiline(default_text, key = "0", size=(80,40))])
    layout.append([sg.Button("确定")])
    window = sg.Window("编译参数确定", layout)
    event, values = window.read(close=True)
    value = values["0"].split("\n")
    helper = set([])
    for line in value:
        if "=" not in line:
            continue
        k, v = line.split("=")
        k = k.strip()
        v = v.strip()
        if k not in kw.keys():
            return None
        helper.add(k)
        kw[k] = v
    for k in kw.keys():
        if k not in helper:
            return None
    return kw

def Compile_Directly(package, gui, objs, **kw):
    command = "-arch=sm_50 -rdc=true -lcudadevrt -lcufft --use_fast_math -O4 -std=c++11"
    for module in package.Modules:
        if hasattr(module, "Make_Command") and getattr(module, "Make_Command").strip():
            line = getattr(module, "Make_Command").strip()
            command += " " + line.replace("$(","{").replace(")","}")

    command = command.format(**kw)

    for module in package.Modules:
        thread = threading.Thread(target = Module_compile, args = (module, command, ))
        thread.start()
  
    
    thread = threading.Thread(target = Obj_compile, args = ("main.o", command, "main.cuh"))
    thread.start()
    
    while threadWorking:
        gui.stdout.update(value = hint)
        gui.stdout.print("\nNow working: %s"%(" ".join(threadWorking)))
        event, values = gui.window.read(timeout = 1e-2)
        if event == None:
            break
    
    thread = threading.Thread(target = Link_compile, args = (objs, command))
    thread.start()
    while threadWorking:
        gui.stdout.update(value = hint)
        gui.stdout.print("\nNow working: %s"%(" ".join(threadWorking)))
        #gui.stdout.print(threadWorking)
        event, values = gui.window.read(timeout = 1e-2)
        if event == None:
            break
    gui.stdout.update(value = hint)
    if os.path.exists("SPONGE.exe"):
        gui.stderr.print("编译成功")
        gui.stderr.print("SPONGE.exe移动至bin文件夹")
        shutil.move("./SPONGE.exe", "../bin/SPONGE.exe")
    else:
        gui.stderr.print("编译未成功")
    


def Event_Module(event, values, gui):
    if values[event]:
        gui.Graph.print(event)
    else:
        content = gui.Graph.get().split("\n")
        while event in content:
            content.remove(event)
        while '' in content:
            content.remove('')
        Update_Graph(gui,content)

def Event_Package(event, values, gui):
    for i in gui.Packages:
        if event == i:
            event = i
            break
    content = []
    for i in event.Modules:
        content.append(i.Name)
    Update_Graph(gui,content)
    for j, module in enumerate(gui.Modules):
        if module.Name in event.Modules:
            gui.Module_Buttons[j].update(value = True)
        else:
            gui.Module_Buttons[j].update(value = False)
    
def Event_Generate_Main(event, values, gui):
    Check_GUI_Graph(values["GUI_Graph"], gui)
    temp = base.Package({"Modules": gui.Graph.get()},gui.Modules)
    print("已生成包含了%d个模块的main.cu"%len(temp.Modules), file = sys.stderr)
    f = open("main.cu","w", encoding="utf-8")
    f.write(temp.Generate_Main_File())
    f.close()
    
def Event_Generate_Main_Head(event, values, gui):
    Check_GUI_Graph(values["GUI_Graph"], gui)
    temp = base.Package({"Modules": gui.Graph.get()},gui.Modules)
    print("已生成包含了%d个模块的main.cuh"%len(temp.Modules), file = sys.stderr)
    f = open("main.cuh","w", encoding="utf-8")
    f.write(temp.Generate_Main_Head_File())
    f.close()

def Event_Generate_Makefile(event, values, gui):
    Check_GUI_Graph(values["GUI_Graph"], gui)
    temp = base.Package({"Modules": gui.Graph.get()},gui.Modules)
    print("已生成包含了%d个模块的Makefile"%len(temp.Modules), file = sys.stderr)
    f = open("Makefile","w", encoding="utf-8")
    f.write(temp.Generate_Makefile())
    f.close()

def Event_Generate_HTML_GUI(event, values, gui):
    Check_GUI_Graph(values["GUI_Graph"], gui)
    temp = base.Package({"Modules": gui.Graph.get()},gui.Modules)
    print("已生成包含了%d个模块的图形化界面"%len(temp.Modules), file = sys.stderr)
    f = open("SPONGE.html","w", encoding="utf-8")
    f.write(temp.Generate_HTML())
    f.close()
    shutil.move("./SPONGE.html", "../bin/SPONGE.html")
    

def Event_Compile_Directly(event, values, gui):
    global hint
    hint = gui.stdout.get()
    Check_GUI_Graph(values["GUI_Graph"], gui)
    temp = base.Package({"Modules": gui.Graph.get()},gui.Modules)
    print("已生成包含了%d个模块的main.cuh"%len(temp.Modules), file = sys.stderr)
    f = open("main.cuh","w", encoding="utf-8")
    f.write(temp.Generate_Main_Head_File())
    f.close()
    print("已生成包含了%d个模块的main.cu"%len(temp.Modules), file = sys.stderr)
    f = open("main.cu","w", encoding="utf-8")
    f.write(temp.Generate_Main_File())
    f.close()
    
    Make_Variables = {}
    objs = "main.o"
    for module in temp.Modules:
        if hasattr(module, "Make_Variables") and getattr(module, "Make_Variables").strip():
            lines = getattr(module, "Make_Variables").split("\n")
            for line in lines:
                line = line.strip()
                if line:
                    key, value = line.split("=")
                    Make_Variables[key] = value
        if hasattr(module, "Module_File") and getattr(module, "Module_File").strip():
            objs += " " + module.Module_File.replace(".cu", ".o")
    
    if Make_Variables:
        Make_Variables = Ask_For_Complete(**Make_Variables)
        if not Make_Variables:
             print("变量编辑有误", file = sys.stderr)
             return None
    print("开始编译，在输出“编译结束”前请勿乱点", file = sys.stderr)
    Compile_Directly(temp, gui, objs, **Make_Variables)
    gui.stdout.update(value = hint)
    print("编译结束", file = sys.stderr)


def Event_All(event, values, gui):
    if event in gui.Modules:
        Event_Module(event, values, gui)
    elif event in gui.Packages:
        Event_Package(event, values, gui)
    elif event == "生成主文件":
        Event_Generate_Main(event, values, gui)
    elif event == "生成主头文件":
        Event_Generate_Main_Head(event, values, gui)
    elif event == "生成Makefile":
        Event_Generate_Makefile(event, values, gui)
    elif event == "直接编译":
        Event_Compile_Directly(event, values, gui)
    elif event == "生成图形化界面":
        Event_Generate_HTML_GUI(event, values, gui)
