# -*- coding:utf-8 -*-
import sys
sys.dont_write_bytecode = True 

from helper_sponge_builder import funcs

import argparse

parser = argparse.ArgumentParser(description='Build Sponge')

parser.add_argument('-search_path', metavar = "path", default = "..", help = "the path to search SPONGE packages and modules.")
parser.add_argument('-write_path', metavar = "path", default = ".", help = "the path to write main files and Makefile.")
parser.add_argument('-bin_path', metavar = "path", default = "../bin", help = "the path to the SPONGE bin folder")

subparsers = parser.add_subparsers(dest = "gui", help='subprogram to use gui or not')
parser_a = subparsers.add_parser('nogui', help='do not use gui. Use "sponge_builder.py nogui -h" to see more detail helps.')
parser_b = subparsers.add_parser('gui', help='use gui. This is the default choice.')

parser_a.add_argument('-package', type = str, metavar = "package", help = "add the package name you want use")
parser_a.add_argument('-modules', type = str, nargs = "+", metavar = "module", help = "add the module names you want use")
parser_a.add_argument('-show_all', action='store_true', help = "show all modules and packages searched and exit")

args = parser.parse_args()

if not args.gui or args.gui == "gui":
    from helper_sponge_builder import gui_main
    from traceback import print_exc
    gui = gui_main.Build_GUI(args.search_path)
    if not gui:
        raise Exception("gui build failed")    
    while 1:
        event, values = gui.window.read()
        if event == None:
            break
        else:
            try:
                funcs.Event_All(event, values, gui)
            except:
                gui.stderr.print(print_exc())
    gui.window.close()
else:
    from helper_sponge_builder import base
    from os.path import join
    modules = funcs.Search_Modules(args.search_path)
    packages = funcs.Search_Packages(args.search_path, modules)
    if args.show_all:
        print("Packages:")
        for p in packages:
            print("\t"+p.Name)
        print("Modules:")
        for m in modules:
            print("\t"+m.Name)
        exit(0)
        
    if not args.package:
        args.package = base.Package({"Modules": "common control md_core"}, modules)
    else:
        for p in packages:
            if args.package == p:
                args.package = p
                break
    if args.modules:
        for m1 in args.modules:
            finded = False
            for m2 in modules:
                if m1 == m2:
                    finded = True
                    if m2 in args.package.Modules:
                        raise Exception("Module has already in the package", m1)
                    else:
                        args.package.Modules.append(m2)
                    break
            if not finded:
                raise Exception("Module not found", m1)
    
    print("已生成包含了%d个模块的main.cuh"%len(args.package.Modules))
    f = open(join(args.write_path, "main.cuh"),"w", encoding="utf-8")
    f.write(args.package.Generate_Main_Head_File())
    f.close()
    print("已生成包含了%d个模块的main.cu"%len(args.package.Modules))
    f = open(join(args.write_path, "main.cu"),"w", encoding="utf-8")
    f.write(args.package.Generate_Main_File())
    f.close()
    print("已生成包含了%d个模块的Makefile"%len(args.package.Modules))
    f = open(join(args.write_path, "Makefile"),"w", encoding="utf-8")
    f.write(args.package.Generate_Makefile())
    f.close()
    print("已生成包含了%d个模块的图形化界面"%len(args.package.Modules))
    f = open(join(args.bin_path, "SPONGE.html"),"w", encoding="utf-8")
    f.write(args.package.Generate_HTML())
    f.close()
        
